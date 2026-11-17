# non-centred update of the persistence rho and the stationary variance w of AR(1) coordinate paths.
# with the whitened innovations e held fixed (u_1 = sqrt(w) e_1, u_t = rho u_t-1 + sqrt(w (1 - rho^2)) e_t),
# the paths follow (rho, w), and the density of e does not depend on them, so a random-walk Metropolis step on
# (atanh(rho / 0.99), log w) with the Gaussian residual likelihood of the position kernels targets the joint
# posterior. it crosses the ridge between the paths and (rho, w) that the conditional draws cross slowly when
# the paths are persistent and the data weak. the residual model is the one rUV_dynamic_fc_cpp() and
# rUV_dynamic_bip_fc_cpp() draw the positions under: every cell of E_t - U_t V_t' (U_t G V_t' in a bipartite
# fit, missing cells at zero) with variance s2, the two cells of a directed dyad with correlation rho_dyad

# with `carry`, binary and poisson fits take a second pair of proposals that carry the latent array Z with
# the product, so the residual and its likelihood are unchanged and the ratio is that of the data given Z:
# a binary carry keeps every observed cell on the side of zero its tie sets, a poisson one has the poisson
# log-likelihood ratio

# `E` the residual array the positions were drawn under, `s2` its variance, `G` the bipartite interaction
# matrix (NULL for a directed fit); `shape`, `scale` and `stationary` the variance prior (R/uv_var_prior.R);
# `stat_var_cap`, `sigma_max`, `max_abs` the support bounds; `carry` NULL or list(family, Z, Y, log_exposure)
# with the per-period offset of a poisson fit. returns the paths, rho, sigma, Z and the accepts
.lame_uv_nc_hyper_step <- function(
	U_cube,
	V_cube,
	E,
	s2,
	rho,
	sigma,
	prior_mean,
	prior_sd,
	shape,
	scale,
	stationary = TRUE,
	G = NULL,
	rho_dyad = 0,
	stat_var_cap = Inf,
	sigma_max = Inf,
	max_abs = Inf,
	tries = 4L,
	step_rho = 0.15,
	step_w = 0.1,
	carry = NULL
) {
	limit <- 0.99
	still <- list(U = U_cube, V = V_cube, rho = rho, sigma = sigma, Z = carry$Z, accepts = 0L, attempts = 0L)
	N <- dim(U_cube)[3L]
	if (is.null(N) || N < 2L || !is.finite(rho) || abs(rho) >= limit || !is.finite(sigma) || sigma <= 0 ||
			!is.finite(s2) || s2 <= 0 || !is.finite(rho_dyad) || abs(rho_dyad) >= 1) return(still)
	bip <- !is.null(G)
	if (bip) E[!is.finite(E)] <- 0
	if (!all(is.finite(E))) return(still)
	dU <- dim(U_cube)
	dV <- dim(V_cube)
	nu <- prod(dU[1:2])
	X <- rbind(matrix(U_cube, nu, N), matrix(V_cube, prod(dV[1:2]), N))
	s2t <- s2 * (1 - rho_dyad) * (1 + rho_dyad)
	loglik <- function(X) {
		q <- 0
		for (t in seq_len(N)) {
			Ut <- matrix(X[seq_len(nu), t], dU[1L], dU[2L])
			Vt <- matrix(X[-seq_len(nu), t], dV[1L], dV[2L])
			Rt <- E[, , t] - if (bip) Ut %*% G %*% t(Vt) else tcrossprod(Ut, Vt)
			q <- q + sum(Rt * Rt)
			if (rho_dyad != 0) q <- q - rho_dyad * sum(Rt * t(Rt))
		}
		-q / (2 * s2t)
	}
	product <- function(X) {
		P <- array(0, dim(E))
		for (t in seq_len(N)) {
			Ut <- matrix(X[seq_len(nu), t], dU[1L], dU[2L])
			Vt <- matrix(X[-seq_len(nu), t], dV[1L], dV[2L])
			P[, , t] <- if (bip) Ut %*% G %*% t(Vt) else tcrossprod(Ut, Vt)
		}
		P
	}
	# log-likelihood of the data given the latent array, up to a constant
	if (!is.null(carry)) {
		obs <- is.finite(carry$Y)
		Zc <- carry$Z
		if (!identical(dim(Zc), dim(E)) || !all(is.finite(Zc[obs]))) carry <- NULL
	}
	if (!is.null(carry)) {
		le <- array(rep(carry$log_exposure %||% 0, each = prod(dim(E)[1:2]), length.out = length(E)), dim(E))
		data_ll <- function(Zn) {
			if (identical(carry$family, "binary")) {
				if (all((Zn[obs] > 0) == (carry$Y[obs] > 0))) 0 else -Inf
			} else sum(carry$Y[obs] * Zn[obs] - exp(Zn[obs] + le[obs]))
		}
		dll <- data_ll(Zc)
		# a latent array the data do not support is left to its own update
		if (!is.finite(dll)) carry <- NULL
	}
	# log prior density of (rho, w): on w itself, or the inverse gamma on sigma^2 = w (1 - rho^2)
	log_prior <- function(r, w) {
		om <- 1 - r * r
		lv <- if (stationary) {
			-(shape + 1) * log(w) - scale / w
		} else -(shape + 1) * log(w * om) - scale / (w * om) + log(om)
		lv + stats::dnorm(r, prior_mean, prior_sd, log = TRUE)
	}
	paths <- function(eps, r, w) {
		out <- eps
		out[, 1L] <- sqrt(w) * eps[, 1L]
		s <- sqrt(w * (1 - r * r))
		for (t in 2:N) out[, t] <- r * out[, t - 1L] + s * eps[, t]
		out
	}
	ok <- function(r, w, Xn) {
		abs(r) < limit && is.finite(w) && w > 0 && w <= stat_var_cap &&
			sqrt(w * (1 - r * r)) <= sigma_max && (!is.finite(max_abs) || max(abs(Xn)) <= max_abs)
	}
	w <- sigma * sigma / (1 - rho * rho)
	eps <- X
	eps[, 1L] <- X[, 1L] / sqrt(w)
	eps[, -1L] <- (X[, -1L, drop = FALSE] - rho * X[, -N, drop = FALSE]) / sqrt(w * (1 - rho * rho))
	ll <- loglik(X)
	accepts <- 0L
	for (k in seq_len(tries)) {
		# persistence at a fixed stationary variance, on the atanh scale with its Jacobian
		a <- atanh(rho / limit)
		rp <- limit * tanh(a + stats::rnorm(1, 0, step_rho))
		if (abs(rp) < limit) {
			Xn <- paths(eps, rp, w)
			if (ok(rp, w, Xn)) {
				lln <- loglik(Xn)
				ratio <- lln - ll + log_prior(rp, w) - log_prior(rho, w) +
					log1p(-(rp / limit)^2) - log1p(-(rho / limit)^2)
				if (is.finite(ratio) && log(stats::runif(1)) < ratio) {
					rho <- rp
					X <- Xn
					ll <- lln
					accepts <- accepts + 1L
				}
			}
		}
		# stationary variance at a fixed persistence, on the log scale with its Jacobian
		wp <- w * exp(stats::rnorm(1, 0, step_w))
		Xn <- paths(eps, rho, wp)
		if (ok(rho, wp, Xn)) {
			lln <- loglik(Xn)
			ratio <- lln - ll + log_prior(rho, wp) - log_prior(rho, w) + log(wp / w)
			if (is.finite(ratio) && log(stats::runif(1)) < ratio) {
				w <- wp
				X <- Xn
				ll <- lln
				accepts <- accepts + 1L
			}
		}
	}
	# the same proposals with the latent array carried along
	if (!is.null(carry)) P <- product(X)
	if (!is.null(carry)) for (k in seq_len(tries)) {
		a <- atanh(rho / limit)
		rp <- limit * tanh(a + stats::rnorm(1, 0, step_rho))
		wp <- w * exp(stats::rnorm(1, 0, step_w))
		for (move in c("rho", "w")) {
			r_new <- if (move == "rho") rp else rho
			w_new <- if (move == "rho") w else wp
			if (abs(r_new) >= limit) next
			Xn <- paths(eps, r_new, w_new)
			if (!ok(r_new, w_new, Xn)) next
			Pn <- product(Xn)
			Zn <- Zc + (Pn - P)
			dlln <- data_ll(Zn)
			jac <- if (move == "rho") log1p(-(r_new / limit)^2) - log1p(-(rho / limit)^2) else log(w_new / w)
			ratio <- dlln - dll + log_prior(r_new, w_new) - log_prior(rho, w) + jac
			if (is.finite(ratio) && log(stats::runif(1)) < ratio) {
				rho <- r_new
				w <- w_new
				X <- Xn
				E <- E + (Pn - P)
				Zc <- Zn
				P <- Pn
				dll <- dlln
				accepts <- accepts + 1L
			}
		}
	}
	n_prop <- if (is.null(carry)) 2L * tries else 4L * tries
	if (!accepts) {
		still$attempts <- n_prop
		return(still)
	}
	U_cube[] <- X[seq_len(nu), ]
	V_cube[] <- X[-seq_len(nu), ]
	list(
		U = U_cube,
		V = V_cube,
		rho = rho,
		sigma = sqrt(w * (1 - rho * rho)),
		Z = if (is.null(carry)) still$Z else Zc,
		accepts = accepts,
		attempts = n_prop
	)
}
