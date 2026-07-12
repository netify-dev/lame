# ffbs sampler for the time-varying bipartite interaction matrix g_t.
# state: g_t = vec(g_t) of length p = ra * rb.
# ar(1) state-space prior:
#   g_t = rho_g * g_{t-1} + eta_t,  eta_t ~ n(0, sigma_g^2 * i_p)
#   g_1 ~ n(0, (sigma_g^2 / (1 - rho_g^2)) * i_p)   (stationary)
# observation per period:
#   vec(e_t) = (v_t kron u_t) g_t + epsilon_t,  epsilon_t ~ n(0, s2 * i)
# carter-kohn / frühwirth-schnatter forward filter + backward sample.
#
# references:
#   carter & kohn (1994) biometrika; frühwirth-schnatter (1994) jtsa.

#' Forward-filter / backward-sample for vec(G_t) under AR(1) state prior
#'
#' One Carter-Kohn sweep over the per-period observations of the form
#' \code{vec(E_t) = H_t g_t + eps_t}, with state transition
#' \code{g_t = rho * g_{t-1} + eta_t}.
#'
#' @param E_cube nA x nB x T residual cube (Z minus base - a - b - UV')
#' @param U_cube nA x RA x T latent row factor cube
#' @param V_cube nB x RB x T latent column factor cube
#' @param s2 scalar observation variance
#' @param rho_G AR(1) coefficient in \code{(-1, 1)} (use \code{rho = 1} for
#'   the RW1 limit; the forward variance is clamped via a tiny floor).
#' @param sigma_G2 state innovation variance
#' @return list with \code{G_cube} (RA x RB x T) and \code{vecG_path}
#'   (p x T draws of the vectorised state).
#' @keywords internal
ffbs_vecG <- function(E_cube, U_cube, V_cube, s2, rho_G, sigma_G2) {
	if (length(dim(E_cube)) != 3L) {
		cli::cli_abort("ffbs_vecG: {.arg E_cube} must be a 3-D array.")
	}
	N <- dim(E_cube)[3L]
	RA <- dim(U_cube)[2L]
	RB <- dim(V_cube)[2L]
	p  <- RA * RB
	if (p == 0L) return(list(G_cube = array(0, dim = c(RA, RB, N)),
	                          vecG_path = matrix(0, 0L, N)))
	# stationary prior variance per state coord (clamped for rho near 1)
	denom <- max(1 - rho_G^2, 1e-6)
	v_stat <- sigma_G2 / denom
	# storage for forward filter
	m_filt <- matrix(0, p, N)
	P_filt <- array(0, dim = c(p, p, N))
	# state prior at t = 1: stationary
	m_prior <- rep(0, p)
	P_prior <- diag(v_stat, p)
	for (t in seq_len(N)) {
		U_t <- U_cube[, , t]
		V_t <- V_cube[, , t]
		# observation matrix h_t = v_t kron u_t (rows correspond to vec(e_t))
		H_t <- kronecker(V_t, U_t)
		y_t <- as.numeric(E_cube[, , t])
		ok <- is.finite(y_t)
		if (any(ok)) {
			Hk <- H_t[ok, , drop = FALSE]
			yk <- y_t[ok]
			# innovation: y - h m_prior; innovation covariance s = h p h' + s2 i
			# we use information-form update for numerical stability when n_obs
			# is large compared to p
			HtH <- crossprod(Hk)
			Hty <- crossprod(Hk, yk)
			# posterior precision = prior precision + hth / s2
			Prec_post <- safe_sympd_inverse(P_prior) + HtH / s2
			# (information-form) posterior mean: prec_post m = prec_prior m_prior + hty/s2
			rhs <- safe_sympd_inverse(P_prior) %*% m_prior + Hty / s2
			P_t <- safe_sympd_inverse(Prec_post)
			m_t <- as.numeric(P_t %*% rhs)
		} else {
			# no obs at t (e.g. all-na panel slice): keep prior as filter posterior
			m_t <- m_prior
			P_t <- P_prior
		}
		m_filt[, t] <- m_t
		P_filt[, , t] <- P_t
		# predict next-step prior
		if (t < N) {
			m_prior <- rho_G * m_t
			P_prior <- (rho_G^2) * P_t + diag(sigma_G2, p)
		}
	}
	# backward sample
	vecG_path <- matrix(0, p, N)
	# draw g_n from n(m_filt_n, p_filt_n)
	L_N <- safe_chol(P_filt[, , N])
	vecG_path[, N] <- m_filt[, N] + as.numeric(L_N %*% stats::rnorm(p))
	for (t in (N - 1L):1L) {
		if (t < 1L) break
		# smoothing: g_t | g_{t+1}, y_{1:n}
		# predictive cov of g_{t+1} | y_{1:t}: p_pred = rho^2 p_t + sigma^2 i
		P_t <- P_filt[, , t]
		m_t <- m_filt[, t]
		P_pred <- (rho_G^2) * P_t + diag(sigma_G2, p)
		# kalman smoother gain: j = rho * p_t * p_pred^-1
		J <- rho_G * P_t %*% safe_sympd_inverse(P_pred)
		m_smooth <- m_t + J %*% (vecG_path[, t + 1L] - rho_G * m_t)
		P_smooth <- P_t - J %*% (rho_G * P_t)
		# symmetric clean-up (numerical)
		P_smooth <- (P_smooth + t(P_smooth)) / 2
		L_s <- safe_chol(P_smooth)
		vecG_path[, t] <- as.numeric(m_smooth) +
			as.numeric(L_s %*% stats::rnorm(p))
	}
	G_cube <- array(0, dim = c(RA, RB, N))
	for (t in seq_len(N)) G_cube[, , t] <- matrix(vecG_path[, t], RA, RB)
	list(G_cube = G_cube, vecG_path = vecG_path)
}

# safe symmetric-pd inverse with jitter ladder (mirrors safe_sympd_inverse
# in dynamic_beta.cpp / r helpers; this is the r fallback used by the ffbs).
safe_sympd_inverse <- function(M, jitters = c(0, 1e-8, 1e-6, 1e-4)) {
	for (j in jitters) {
		Mj <- if (j > 0) M + diag(j, nrow(M)) else M
		out <- tryCatch(chol2inv(chol(Mj)), error = function(e) NULL)
		if (!is.null(out)) return(out)
	}
	# final fallback: pseudo-inverse via svd (rare path)
	sv <- svd(M)
	tol <- max(sv$d) * .Machine$double.eps * length(sv$d)
	d_inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
	sv$v %*% diag(d_inv, length(d_inv)) %*% t(sv$u)
}

safe_chol <- function(M, jitters = c(0, 1e-8, 1e-6, 1e-4)) {
	for (j in jitters) {
		Mj <- if (j > 0) M + diag(j, nrow(M)) else M
		out <- tryCatch(chol(Mj), error = function(e) NULL)
		if (!is.null(out)) return(t(out))  # lower-triangular for l %*% l'
	}
	# final fallback: svd-based factor
	sv <- svd(M)
	d_pos <- pmax(sv$d, 0)
	sv$u %*% diag(sqrt(d_pos), length(d_pos))
}

#' Logit-scale MH update on the AR(1) coefficient rho_G
#'
#' Proposes \code{rho* = tanh(atanh(rho) + tau * N(0,1))}; under a uniform
#' \code{rho ~ U(-1, 1)} prior the proposal is symmetric in the Fisher-z
#' scale and the Jacobian is \code{log(1 - rho^2)}. The MH ratio
#' compares the AR(1) prior density of \code{vec(G)_path} at the
#' current vs proposed \code{rho}.
#'
#' @param rho_G current AR(1) coefficient
#' @param vecG_path p x N matrix of FFBS-sampled states
#' @param sigma_G2 state innovation variance
#' @param tau RW proposal SD on the Fisher-z scale
#' @return list with \code{rho} (updated), \code{accept} (logical)
#' @keywords internal
sample_rho_G_mh <- function(rho_G, vecG_path, sigma_G2, tau = 0.3) {
	N <- ncol(vecG_path)
	if (N < 2L) return(list(rho = rho_G, accept = FALSE))
	z_curr <- atanh(rho_G)
	z_prop <- z_curr + tau * stats::rnorm(1)
	# clamp the persistence strictly below 1. the absolute scale of vec(g_t)
	# is only weakly identified (the model identifies u g v', not g alone),
	# so an unclamped rho_g drifts to the unit root and the state-space prior
	# degenerates into an unbounded random walk -- vec(g_t) then inflates by
	# orders of magnitude with no change to the fitted linear predictor.
	# holding rho_g <= 0.98 keeps the prior stationary with finite marginal
	# variance sigma_g2 / (1 - rho_g^2), bounding the g_t scale.
	rho_prop <- max(min(tanh(z_prop), 0.98), -0.98)
	# log target on rho-scale = ar(1) log-prior on g + log|j|;
	# the fisher-z proposal is symmetric on z, so the jacobian from
	# d rho / d z = 1 - rho^2 appears once.
	loglik <- function(rho) {
		# stationary marginal for t = 1
		v_stat <- sigma_G2 / max(1 - rho^2, 1e-6)
		ll <- sum(stats::dnorm(vecG_path[, 1L], 0, sqrt(v_stat), log = TRUE))
		if (N >= 2L) {
			for (t in 2L:N) {
				ll <- ll + sum(stats::dnorm(vecG_path[, t],
					rho * vecG_path[, t - 1L], sqrt(sigma_G2), log = TRUE))
			}
		}
		ll
	}
	la_curr <- loglik(rho_G) + log(max(1 - rho_G^2, 1e-12))
	la_prop <- loglik(rho_prop) + log(max(1 - rho_prop^2, 1e-12))
	log_a <- la_prop - la_curr
	if (!is.finite(log_a)) return(list(rho = rho_G, accept = FALSE))
	if (log(stats::runif(1)) < log_a) {
		list(rho = rho_prop, accept = TRUE)
	} else {
		list(rho = rho_G, accept = FALSE)
	}
}

#' Inverse-gamma posterior draw for sigma_G^2 given the FFBS path
#'
#' Conjugate IG update on \code{sigma_G^2} given the AR(1) innovations
#' \code{eta_t = g_t - rho * g_{t-1}}, t = 2..N, plus the stationary
#' density at t = 1. Default prior is IG(2, 1) on \code{sigma_G^2}.
#'
#' Scale identification. The model identifies only the product
#' \eqn{U_t G_t V_t'}, not \eqn{G_t} alone, so the overall scale of
#' \code{vec(G_t)} is free: left unchecked the chain finds a degenerate
#' mode where \code{sigma_G^2} (and hence \code{G_t}) inflates by orders
#' of magnitude while \code{U,V} collapse to compensate, leaving the
#' linear predictor unchanged but the reported \code{G_cube} meaningless.
#' Clamping \code{rho_G} alone does not bound this because the stationary
#' state variance is \code{sigma_G^2 / (1 - rho_G^2)}. We therefore cap
#' the implied stationary variance at \code{v_cap_mult * s2_obs} (a few
#' observation-noise units), which forces the multiplicative scale onto
#' \code{U,V} -- where the regularising N(0, s2) prior pins it -- and
#' keeps \code{G_t} on a scale comparable to a static-G fit. The cap is
#' loose enough never to bind on a genuinely small-variation \code{G_t}.
#'
#' @param s2_obs observation-noise variance (1 for probit/binary).
#' @param v_cap_mult cap on the stationary G-state variance in units of
#'   \code{s2_obs} (default 4).
#' @keywords internal
sample_sigma_G2 <- function(vecG_path, rho_G, prior_shape = 2, prior_rate = 1,
                            s2_obs = 1, v_cap_mult = 4) {
	N <- ncol(vecG_path)
	p <- nrow(vecG_path)
	if (p == 0L || N < 2L) return(min(1.0, v_cap_mult * max(s2_obs, 1e-8)))
	# stationary at t = 1: g_1 ~ n(0, sigma^2 / (1 - rho^2)) ->
	# (1 - rho^2) g_1^2 ~ sigma^2 * chi^2_1 in expectation; we treat the
	# stationary contribution as an effective sum of squares too for
	# stability.
	rho2 <- max(1 - rho_G^2, 1e-6)
	ss_stat <- sum(vecG_path[, 1L]^2) * rho2
	# innovations 2..n
	if (N >= 2L) {
		eta <- vecG_path[, -1L, drop = FALSE] -
			rho_G * vecG_path[, -N, drop = FALSE]
		ss_innov <- sum(eta^2)
	} else {
		ss_innov <- 0
	}
	shape_post <- prior_shape + (p * N) / 2
	rate_post  <- prior_rate  + (ss_stat + ss_innov) / 2
	sig2 <- 1 / stats::rgamma(1, shape = shape_post, rate = rate_post)
	# scale-identification cap: stationary variance sig2/(1-rho^2) <= cap
	sig2_cap <- v_cap_mult * max(s2_obs, 1e-8) * rho2
	min(sig2, sig2_cap)
}

#' Rotation-drift diagnostic for the canonical (U_t, G_t, V_t) trio
#'
#' Reports a scalar ratio comparing the variance of the *raw* G_t
#' entries (per element across t) to the variance of the *canonical*
#' G_t entries (after per-period SVD). When the raw / canonical
#' variance ratio is large (>= 5), the apparent G_t time-variation is
#' dominated by rotation drift in U_t, V_t rather than real temporal
#' change. The canonicalisation removes that rotation and the
#' canonical entries should be the user-facing summary.
#'
#' @param G_cube_raw RA x RB x T raw (FFBS-sampled) G cube
#' @param U_cube nA x RA x T raw U cube
#' @param V_cube nB x RB x T raw V cube
#' @return list with \code{ratio} (numeric), \code{var_raw} (RA x RB),
#'   \code{var_canonical} (RA x RB), \code{flag} (logical: ratio >= 5)
#' @keywords internal
.rotation_drift_diagnostic <- function(G_cube_raw, U_cube, V_cube) {
	canon <- .canonicalize_UGV(U_cube, G_cube_raw, V_cube)
	var_raw <- apply(G_cube_raw, c(1, 2), stats::var)
	var_can <- apply(canon$G, c(1, 2), stats::var)
	# scalar summary: geometric mean of element-wise ratios on nonzero cells
	cells <- which(var_can > 1e-12, arr.ind = TRUE)
	if (nrow(cells) == 0L) {
		return(list(ratio = NA_real_, var_raw = var_raw,
		            var_canonical = var_can, flag = FALSE))
	}
	r_vec <- numeric(nrow(cells))
	for (k in seq_len(nrow(cells))) {
		i <- cells[k, 1L]; j <- cells[k, 2L]
		r_vec[k] <- var_raw[i, j] / var_can[i, j]
	}
	ratio <- exp(mean(log(pmax(r_vec, 1e-12))))
	list(ratio = ratio, var_raw = var_raw, var_canonical = var_can,
	     flag = isTRUE(ratio >= 5))
}
