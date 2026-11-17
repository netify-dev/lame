# exact draws of the size of the dynamic position product of binary and Poisson fits (?lame,
# prior$uv_size_move). the gibbs kernel moves a period's product U_t G V_t' only as far as the latent
# array lets it, so its size and the latent variance trade slowly. two one-parameter groups move it with
# the latent array integrated out (binary: the probit factor of each observed cell, with a directed
# network's dyadic correlation kept through the shared dyad factor; the array is redrawn after the step)
# or carried (Poisson: the cells move with the predictor, so every residual and the gaussian density of
# the array stay, and the Poisson factor of the observed cells is the likelihood):
#   every period: the product times c = e^u, through G -> c G (bipartite; prior N(0, g_var) cut at
#     g_max, jacobian c^dim(G)) or U, V and sigma_uv times c^(1/2) (directed), which keeps every
#     standardised innovation, so the paths see only the reset positions of a snap fit with a fixed
#     snap_kappa and sigma_uv its inverse-gamma prior (jacobian c^(1/2) for sigma_uv)
#   one period: U_t and V_t times s = e^(u/2), the transitions into and out of t (w, k as in
#     R/uv_split_dynamic.R) giving the prior -(A s^2 - 2 B s) / (2 sigma^2) and the jacobian s^(dim U_t + dim V_t);
#     under kind t the weights of those transitions (and of the first period) are integrated out
#     (multivariate t), the split U_t c, V_t / c of the period is drawn under the same marginal, and the
#     weights are drawn again
# each log factor is drawn by a slice step from 0 within the bounds of G, sigma_uv, its stationary
# variance and max_abs

# observed cells of the panel and their data; `offset` the log exposure of each period (Poisson)
.lame_uvps_context <- function(Y, family, mode, offset = NULL) {
	d <- dim(Y)
	N <- d[3L]
	obs <- is.finite(Y)
	dg <- array(FALSE, d)
	if (!identical(mode, "bip")) dg[cbind(rep(seq_len(d[1L]), N), rep(seq_len(d[1L]), N), rep(seq_len(N), each = d[1L]))] <- TRUE
	obs <- obs & !dg
	idx <- which(obs)
	per <- arrayInd(idx, d)[, 3L]
	cx <- list(family = family, mode = mode, idx = idx, y = Y[idx], per = lapply(seq_len(N), function(t) which(per == t)))
	if (identical(family, "normal")) {
		cx$obs <- obs
		if (!identical(mode, "bip")) cx$upper <- upper.tri(matrix(0, d[1L], d[2L]))
	}
	if (identical(family, "poisson")) cx$off <- if (is.null(offset)) numeric(length(idx)) else as.numeric(offset)[per]
	if (identical(family, "binary")) {
		cx$sg <- 2 * Y[idx] - 1
		if (!identical(mode, "bip")) {
			cx$i1 <- which(Y == 1 & !dg)
			cx$i0 <- which(Y == 0 & !dg)
			cx$ina <- which(is.na(Y) & !dg)
			cx$idg <- which(dg)
			cx$upper <- upper.tri(matrix(0, d[1L], d[2L]))
		}
	}
	cx
}

# likelihood of a change `d` of the predictor on the observed cells (`k` a subset of them) given the
# latent array `Z` and its mean `EZ`, and the latent array of the state after a change `D` of the whole
# predictor; `rho` the dyadic correlation of a directed fit, `s2` the latent variance of a normal fit, whose
# unobserved cells are carried and whose pairs keep their correlation (`ll_full` takes the whole change)
.lame_uvps_like <- function(cx, Z, EZ, Y, rho = 0, s2 = 1) {
	if (identical(cx$family, "normal")) {
		e0 <- Z - EZ
		e0[!is.finite(e0)] <- 0
		dir <- !identical(cx$mode, "bip")
		q_of <- function(e) {
			if (!dir) return(sum(e[cx$obs]^2))
			q <- 0
			for (t in seq_len(dim(e)[3L])) {
				m <- e[, , t]
				a <- m[cx$upper]
				b <- t(m)[cx$upper]
				q <- q + sum(a * a - 2 * rho * a * b + b * b)
			}
			q / (1 - rho * rho)
		}
		q0 <- q_of(e0)
		ll_full <- function(D) -(q_of(e0 - D * cx$obs) - q0) / (2 * s2)
		return(list(ll_full = ll_full, finish = function(D) Z + D * !cx$obs))
	}
	if (identical(cx$family, "poisson")) {
		z0 <- Z[cx$idx] + cx$off
		ll <- function(d, k = NULL) {
			if (is.null(k)) return(sum(cx$y * (z0 + d) - exp(z0 + d)))
			z <- z0[k] + d
			sum(cx$y[k] * z - exp(z))
		}
		return(list(ll = ll, ll_full = function(D) ll(D[cx$idx]), finish = function(D) Z + D))
	}
	dir <- !identical(cx$mode, "bip")
	off <- if (dir) .lame_dyad_factor_draw(Z, EZ, rho, cx$upper) else NULL
	sd_cell <- if (dir) sqrt(1 - abs(rho)) else 1
	b0 <- EZ[cx$idx] + (if (dir) off[cx$idx] else 0)
	ll <- function(d, k = NULL) {
		if (is.null(k)) return(sum(stats::pnorm(cx$sg * (b0 + d) / sd_cell, log.p = TRUE)))
		sum(stats::pnorm(cx$sg[k] * (b0[k] + d) / sd_cell, log.p = TRUE))
	}
	finish <- function(D) {
		if (dir) .lame_binary_z_given_factor(Z, EZ + D, off, rho, cx$i1, cx$i0, cx$ina, cx$idg)
		else rZ_bin_bip_batch_cpp(Z, EZ + D, Y)
	}
	list(ll = ll, ll_full = function(D) ll(D[cx$idx]), finish = finish)
}

# the product of each period (rows x columns x N)
.lame_uvps_product <- function(U_cube, V_cube, G = NULL) {
	N <- dim(U_cube)[3L]
	P <- array(0, c(dim(U_cube)[1L], dim(V_cube)[1L], N))
	for (t in seq_len(N)) {
		Ut <- matrix(U_cube[, , t], dim(U_cube)[1L])
		if (!is.null(G)) Ut <- Ut %*% G
		P[, , t] <- tcrossprod(Ut, matrix(V_cube[, , t], dim(V_cube)[1L]))
	}
	P
}

# log density of the transitions into and out of period t of one side scaled by s, the weights of kind t
# integrated out (`nu` their degrees of freedom; the first period a t of scale sigma / sqrt(1 - rho^2));
# `X` actors x R x N
.lame_uvps_t_terms <- function(X, t, rho, s2, nu) {
	n <- dim(X)[1L]
	R <- dim(X)[2L]
	N <- dim(X)[3L]
	xt <- matrix(X[, , t], n, R)
	sq <- rowSums(xt * xt)
	if (t == 1L) {
		into <- function(s) -(nu + R) / 2 * sum(log1p((1 - rho * rho) * s * s * sq / (nu * s2)))
	} else {
		xp <- matrix(X[, , t - 1L], n, R)
		cr <- rowSums(xt * xp)
		pp <- rowSums(xp * xp)
		into <- function(s) -(nu + R) / 2 * sum(log1p((s * s * sq - 2 * s * rho * cr + rho * rho * pp) / (nu * s2)))
	}
	if (t == N) return(into)
	xn <- matrix(X[, , t + 1L], n, R)
	cn <- rowSums(xn * xt)
	nn <- rowSums(xn * xn)
	function(s) into(s) - (nu + R) / 2 * sum(log1p((nn - 2 * s * rho * cn + rho * rho * s * s * sq) / (nu * s2)))
}

# fresh weights of kind t for the transitions into and out of period t of one side, the first period's
# weight at t = 1 (`lam` actors x N)
.lame_uvps_t_weights <- function(X, lam, t, rho, s2, nu) {
	n <- dim(X)[1L]
	R <- dim(X)[2L]
	for (s in intersect(c(t, t + 1L), seq_len(dim(X)[3L]))) {
		e <- matrix(X[, , s], n, R)
		q <- if (s == 1L) (1 - rho * rho) * rowSums(e * e) / s2 else {
			e <- e - rho * matrix(X[, , s - 1L], n, R)
			rowSums(e * e) / s2
		}
		lam[, s] <- stats::rgamma(n, (nu + R) / 2, rate = (nu + q) / 2)
	}
	lam
}

# one visit: the size of every period's product, then each period alone. `cx` from .lame_uvps_context(),
# `Z` the latent array, `EZ` its mean with the current product, `rho` the dyadic correlation (directed
# binary), `G` the interaction (bipartite), `rho_uv`, `sigma_uv` the path parameters with the inverse-gamma
# prior of sigma_uv^2 given rho_uv (`shape`, `scale`, R/uv_var_prior.R) and bounds `sigma_max`, `w_cap` of
# sigma_uv; `kind` the path kind with the weights `lambda_u`, `lambda_v` and their degrees of freedom `nu`
# (kind t) or the reset scale `kappa`, fixed or NULL when linked to sigma_uv, and the resets `delta_u`,
# `delta_v` (snap). returns the state and the number of steps that moved
.lame_uvps_visit <- function(
	cx,
	Z,
	EZ,
	Y,
	U_cube,
	V_cube,
	rho_uv,
	sigma_uv,
	shape,
	scale,
	rho = 0,
	G = NULL,
	g_var = 1,
	g_max = Inf,
	sigma_max = Inf,
	w_cap = Inf,
	max_abs = Inf,
	kind = "ar1",
	lambda_u = NULL,
	lambda_v = NULL,
	nu = 4,
	kappa = NULL,
	delta_u = NULL,
	delta_v = NULL,
	width = 0.25
) {
	out <- list(Z = Z, U = U_cube, V = V_cube, G = G, sigma = sigma_uv, lambda_u = lambda_u, lambda_v = lambda_v, moved = 0L)
	if (!is.finite(rho_uv) || abs(rho_uv) >= 1 || !is.finite(sigma_uv) || sigma_uv <= 0) return(out)
	if (!all(is.finite(U_cube)) || !all(is.finite(V_cube)) || !all(is.finite(EZ[cx$idx]))) return(out)
	bip <- identical(cx$mode, "bip")
	if (bip && (is.null(G) || !all(is.finite(G)))) return(out)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	if (!is.finite(g_max) || g_max <= 0) g_max <- Inf
	N <- dim(U_cube)[3L]
	P <- .lame_uvps_product(U_cube, V_cube, G)
	like <- .lame_uvps_like(cx, Z, EZ, Y, rho)
	po <- P[cx$idx]
	D <- array(0, dim(P))
	moved <- 0L
	# every period: the product times e^u
	if (bip) {
		g2 <- sum(G * G)
		dG <- length(G)
		hi <- if (is.finite(g_max)) log(g_max / max(abs(G))) else Inf
		lg <- function(u) like$ll(expm1(u) * po) - exp(2 * u) * g2 / (2 * g_var) + dG * u
	} else {
		# reset positions of a snap fit drawn at a fixed scale do not scale with sigma_uv
		sr <- dr <- 0
		if (identical(kind, "snap") && !is.null(kappa)) {
			st <- .lame_snap_reset_ss(U_cube, V_cube, FALSE, delta_u, delta_v)
			sr <- st$ss / (kappa * kappa)
			dr <- st$cnt
		}
		lsig <- function(sg) -(2 * shape + 1) * log(sg) - scale / (sg * sg)
		cap <- min(sigma_max, sqrt(w_cap * (1 - rho_uv * rho_uv)))
		hi <- min(2 * log(cap / sigma_uv), 2 * log(max_abs / max(abs(U_cube), abs(V_cube))))
		lg <- function(u) {
			like$ll(expm1(u) * po) + lsig(sigma_uv * exp(u / 2)) + u / 2 - expm1(u) * sr / 2 + dr * u / 2
		}
	}
	if (hi > 0) {
		u <- .lame_nc_slice(lg, 0, -Inf, hi, w = width)$x
		if (is.finite(u) && u != 0) {
			D <- D + expm1(u) * P
			P <- exp(u) * P
			po <- exp(u) * po
			if (bip) {
				out[["G"]] <- exp(u) * G
			} else {
				U_cube <- U_cube * exp(u / 2)
				V_cube <- V_cube * exp(u / 2)
				out$sigma <- sigma_uv * exp(u / 2)
			}
			moved <- moved + 1L
		}
	}
	sigma_uv <- out$sigma
	s2 <- sigma_uv * sigma_uv
	p1 <- 1 - rho_uv * rho_uv
	# one period at a time: U_t, V_t times e^(u/2), the reset weights at the sigma_uv of the step above
	if (!identical(kind, "t")) {
		kw <- if (identical(kind, "snap")) .lame_uv_dyn_kind_weights(
			"snap",
			NULL,
			NULL,
			delta_u,
			delta_v,
			s2 / .lame_snap_kappa(kappa, rho_uv, sigma_uv)^2
		)
		one <- function(x, n) if (is.null(x)) matrix(1, n, N) else matrix(as.numeric(x), n, N)
		w_u <- one(kw$w_u, dim(U_cube)[1L]); w_v <- one(kw$w_v, dim(V_cube)[1L])
		k_u <- one(kw$k_u, dim(U_cube)[1L]); k_v <- one(kw$k_v, dim(V_cube)[1L])
		su <- .lame_uvsd_sums(U_cube)
		sv <- .lame_uvsd_sums(V_cube)
	}
	dt <- prod(dim(U_cube)[1:2]) + prod(dim(V_cube)[1:2])
	for (t in seq_len(N)) {
		k <- cx$per[[t]]
		if (length(k) == 0L) next
		pt <- P[cx$idx[k]]
		# the change of the observed cells: the base likelihood holds every earlier change of the visit
		dk <- D[cx$idx[k]]
		if (identical(kind, "t")) {
			fu <- .lame_uvps_t_terms(U_cube, t, rho_uv, s2, nu)
			fv <- .lame_uvps_t_terms(V_cube, t, rho_uv, s2, nu)
			lp <- function(s) fu(s) + fv(s)
		} else {
			pu <- .lame_uvsd_period(su, w_u, k_u, rho_uv, p1, t)
			pv <- .lame_uvsd_period(sv, w_v, k_v, rho_uv, p1, t)
			A <- pu[["a"]] + pv[["a"]]
			B <- pu[["b"]] + pv[["b"]]
			lp <- function(s) -(A * s * s - 2 * B * s) / (2 * s2)
		}
		hi <- 2 * log(max_abs / max(abs(U_cube[, , t]), abs(V_cube[, , t])))
		if (!(hi > 0)) next
		lg <- function(u) like$ll(dk + expm1(u) * pt, k) + lp(exp(u / 2)) + dt * u / 2
		u <- .lame_nc_slice(lg, 0, -Inf, hi, w = width)$x
		if (is.finite(u) && u != 0) {
			D[, , t] <- D[, , t] + expm1(u) * P[, , t]
			P[, , t] <- exp(u) * P[, , t]
			U_cube[, , t] <- U_cube[, , t] * exp(u / 2)
			V_cube[, , t] <- V_cube[, , t] * exp(u / 2)
			moved <- moved + 1L
			if (!identical(kind, "t")) {
				# rescale the sums that touch period t
				rescale <- function(sm) {
					sm$sq[, t] <- sm$sq[, t] * exp(u)
					sm$cr[, t] <- sm$cr[, t] * exp(u / 2)
					if (t < N) sm$cr[, t + 1L] <- sm$cr[, t + 1L] * exp(u / 2)
					sm
				}
				su <- rescale(su)
				sv <- rescale(sv)
			}
		}
		if (identical(kind, "t")) {
			# the split of the period, U_t c and V_t / c, with the same weights integrated out
			fu <- .lame_uvps_t_terms(U_cube, t, rho_uv, s2, nu)
			fv <- .lame_uvps_t_terms(V_cube, t, rho_uv, s2, nu)
			lj <- (prod(dim(U_cube)[1:2]) - prod(dim(V_cube)[1:2]))
			lm <- c(log(max(abs(V_cube[, , t])) / max_abs), log(max_abs / max(abs(U_cube[, , t]))))
			if (lm[1L] <= 0 && lm[2L] >= 0) {
				x <- .lame_nc_slice(function(v) fu(exp(v)) + fv(exp(-v)) + lj * v, 0, lm[1L], lm[2L], w = width)$x
				if (is.finite(x) && x != 0) {
					U_cube[, , t] <- U_cube[, , t] * exp(x)
					V_cube[, , t] <- V_cube[, , t] * exp(-x)
				}
			}
			out$lambda_u <- .lame_uvps_t_weights(U_cube, out$lambda_u, t, rho_uv, s2, nu)
			out$lambda_v <- .lame_uvps_t_weights(V_cube, out$lambda_v, t, rho_uv, s2, nu)
		}
	}
	if (moved > 0L) out$Z <- like$finish(D)
	out$U <- U_cube
	out[["V"]] <- V_cube
	out$moved <- moved
	out
}

# non-centred draw of the persistence of directed and bipartite paths (?lame, prior$uv_rho_nc_move). with the
# stationary variance w = sigma^2 / (1 - rho^2) and the standardised innovations held,
#     x_1 = (w / lambda_1)^(1/2) e_1,   x_t = k_t rho x_t-1 + sd_t e_t,   sd_t = (w (1 - rho^2) / lambda_t)^(1/2) (a drift),
# kappa or w^(1/2) (a snap reset at a fixed or linked scale), the paths follow rho, the innovations keep their
# standard normal law and sigma^2 = w (1 - rho^2) its inverse-gamma prior, so rho has the density
#     pi(rho) (1 - rho^2)^(-a) exp(-b / (w (1 - rho^2))) L(rho)
# with L the likelihood of the product the paths give (latent array integrated out, carried or, normal, with
# its unobserved cells carried). under the inverse gamma on w itself (`stationary`, R/uv_var_prior.R) the prior
# of (rho, w) is pi(rho) IG(w), and with w held the density is pi(rho) L(rho). rho is drawn by a slice step within
# +-0.99, sigma_uv_max and max_abs

# the standardised innovations of one side (`lam` the weights of kind t, `del` the snap resets, `kap` the
# reset scale or NULL when linked) and the paths they give at persistence `r`
.lame_uvnc_eps <- function(X, rho, w, lam = NULL, del = NULL, kap = NULL) {
	N <- dim(X)[3L]
	E <- X
	E[, , 1L] <- X[, , 1L] / sqrt(w / (if (is.null(lam)) 1 else lam[, 1L]))
	for (t in seq_len(N)[-1L]) {
		lt <- if (is.null(lam)) 1 else lam[, t]
		E[, , t] <- (X[, , t] - rho * X[, , t - 1L]) / sqrt(w * (1 - rho * rho) / lt)
		if (!is.null(del) && any(del[, t] > 0.5)) {
			hit <- del[, t] > 0.5
			E[hit, , t] <- X[hit, , t] / (if (is.null(kap)) sqrt(w) else kap)
		}
	}
	E
}
.lame_uvnc_paths <- function(E, r, w, lam = NULL, del = NULL, kap = NULL) {
	N <- dim(E)[3L]
	X <- E
	X[, , 1L] <- E[, , 1L] * sqrt(w / (if (is.null(lam)) 1 else lam[, 1L]))
	for (t in seq_len(N)[-1L]) {
		lt <- if (is.null(lam)) 1 else lam[, t]
		X[, , t] <- r * X[, , t - 1L] + E[, , t] * sqrt(w * (1 - r * r) / lt)
		if (!is.null(del) && any(del[, t] > 0.5)) {
			hit <- del[, t] > 0.5
			X[hit, , t] <- E[hit, , t] * (if (is.null(kap)) sqrt(w) else kap)
		}
	}
	X
}

# arguments as for .lame_uvps_visit, with the normal prior (`prior_mean`, `prior_sd`) of rho, `s2` the latent
# variance of a normal fit and `stationary` TRUE when `shape`, `scale` are those of the inverse gamma on w;
# returns the latent array, the paths, rho_uv, sigma_uv and whether rho moved
.lame_uv_nc_rho_draw <- function(
	cx,
	Z,
	EZ,
	Y,
	U_cube,
	V_cube,
	rho_uv,
	sigma_uv,
	prior_mean,
	prior_sd,
	shape,
	scale,
	rho = 0,
	s2 = 1,
	G = NULL,
	sigma_max = Inf,
	max_abs = Inf,
	kind = "ar1",
	lambda_u = NULL,
	lambda_v = NULL,
	kappa = NULL,
	delta_u = NULL,
	delta_v = NULL,
	bound = 0.99,
	width = 0.05,
	stationary = FALSE
) {
	out <- list(Z = Z, U = U_cube, V = V_cube, rho = rho_uv, sigma = sigma_uv, moved = FALSE)
	if (!is.finite(rho_uv) || abs(rho_uv) >= bound || !is.finite(sigma_uv) || sigma_uv <= 0) return(out)
	if (dim(U_cube)[3L] < 2L || !all(is.finite(U_cube)) || !all(is.finite(V_cube))) return(out)
	if (!all(is.finite(EZ[cx$idx]))) return(out)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	snap <- identical(kind, "snap")
	lu <- if (identical(kind, "t")) lambda_u
	lv <- if (identical(kind, "t")) lambda_v
	du <- if (snap) delta_u
	dv <- if (snap) delta_v
	kap <- if (snap) kappa
	w <- sigma_uv^2 / (1 - rho_uv^2)
	Eu <- .lame_uvnc_eps(U_cube, rho_uv, w, lu, du, kap)
	Ev <- .lame_uvnc_eps(V_cube, rho_uv, w, lv, dv, kap)
	P0 <- .lame_uvps_product(U_cube, V_cube, G)
	like <- .lame_uvps_like(cx, Z, EZ, Y, rho, s2)
	lg <- function(r) {
		if (!is.finite(r) || abs(r) >= bound) return(-Inf)
		om <- 1 - r * r
		if (w * om > sigma_max^2) return(-Inf)
		Un <- .lame_uvnc_paths(Eu, r, w, lu, du, kap)
		Vn <- .lame_uvnc_paths(Ev, r, w, lv, dv, kap)
		if (is.finite(max_abs) && max(abs(Un), abs(Vn)) > max_abs) return(-Inf)
		lp <- stats::dnorm(r, prior_mean, prior_sd, log = TRUE)
		if (!stationary) lp <- lp - shape * log(om) - scale / (w * om)
		lp + like$ll_full(.lame_uvps_product(Un, Vn, G) - P0)
	}
	r <- .lame_nc_slice(lg, rho_uv, -bound, bound, w = width)$x
	if (!is.finite(r) || r == rho_uv) return(out)
	out$U <- .lame_uvnc_paths(Eu, r, w, lu, du, kap)
	out[["V"]] <- .lame_uvnc_paths(Ev, r, w, lv, dv, kap)
	out$Z <- like$finish(.lame_uvps_product(out$U, out[["V"]], G) - P0)
	out$rho <- r
	out$sigma <- sqrt(w * (1 - r * r))
	out$moved <- TRUE
	out
}
