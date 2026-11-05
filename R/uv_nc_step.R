# non-centred draws of the persistence and scale of dynamic_uv coordinate
# paths, with the paths carried along (given the paths, rho_uv and sigma_uv are
# pinned far more tightly than their posteriors).
#   symmetric fits (unit stationary variance): the standardised innovations
#     eps_1 = u_1, eps_t = (u_t - rho u_{t-1}) / sqrt(1 - rho^2) are N(0, 1) for
#     any rho, so given eps rho is drawn by a slice step (neal 2003) from
#         pi(rho) proportional to  prior(rho) * exp(-sum_t sum_{i<j} (e_tij - sum_r L_r u_tir u_tjr)^2 / (2 s2))
#     with the paths u(rho, eps) following rho
#   bipartite fits (U G V'): scaling both sides' paths and sigma by c and G by
#     1/c^2 leaves the likelihood unchanged; log c is drawn by a slice step from
#     the priors of sigma and G with Jacobian c^(1 - 2 d_G) (group move of liu and
#     wu 1999), cut at the support bounds

# slice step on (lo, hi) for the log density `lg`, started at x0 (stepping out
# then shrinkage)
.lame_nc_slice <- function(lg, x0, lo, hi, w = 0.05, m = 40L) {
	lg0 <- lg(x0)
	if (!is.finite(lg0)) return(list(x = x0, evals = 1L))
	thr <- lg0 + log(stats::runif(1L))
	left <- max(lo, x0 - w * stats::runif(1L))
	right <- min(hi, left + w)
	j <- floor(m * stats::runif(1L))
	k <- (m - 1L) - j
	evals <- 1L
	while (j > 0L && left > lo) {
		evals <- evals + 1L
		if (!(lg(left) > thr)) break
		left <- max(lo, left - w)
		j <- j - 1L
	}
	while (k > 0L && right < hi) {
		evals <- evals + 1L
		if (!(lg(right) > thr)) break
		right <- min(hi, right + w)
		k <- k - 1L
	}
	repeat {
		x <- stats::runif(1L, left, right)
		val <- lg(x)
		evals <- evals + 1L
		if (is.finite(val) && val > thr) return(list(x = x, evals = evals))
		if (x < x0) left <- x else right <- x
		if (right - left < 1e-12) return(list(x = x0, evals = evals))
	}
}

# paths of the unit-variance gauge with persistence rho and standardised innovations `eps` (n x R x N)
.lame_uv_unit_paths <- function(eps, rho) {
	N <- dim(eps)[3L]
	U <- eps
	s <- sqrt(1 - rho * rho)
	if (N > 1L) for (t in 2:N) U[, , t] <- rho * U[, , t - 1L] + s * eps[, , t]
	U
}

# non-centred draw of the persistence of a symmetric fit. `U_cube` the paths
# (n x R x N), `L` the direction weights, `E` the latent residual without the
# product of the positions, `s2` its variance; normal prior on rho cut at
# +- bound, `max_abs` the bound on the coordinates
.lame_uv_sym_nc_rho <- function(U_cube, L, E, s2, rho, prior_mean = 0.9, prior_sd = 0.1, max_abs = Inf, bound = 0.99) {
	d <- dim(U_cube)
	n <- d[1L]; R <- d[2L]; N <- d[3L]
	still <- list(rho = rho, U = U_cube, evals = 0L)
	if (is.null(N) || N < 2L || !is.finite(rho) || abs(rho) >= bound || !is.finite(s2) || s2 <= 0) return(still)
	s <- sqrt(1 - rho * rho)
	eps <- U_cube
	for (t in 2:N) eps[, , t] <- (U_cube[, , t] - rho * U_cube[, , t - 1L]) / s
	ut <- upper.tri(matrix(0, n, n))
	ok <- array(FALSE, c(n, n, N))
	for (t in seq_len(N)) ok[, , t] <- ut & is.finite(E[, , t])
	E0 <- E
	E0[!ok] <- 0
	Ld <- diag(L, R)
	ss <- function(U) {
		tot <- 0
		for (t in seq_len(N)) {
			Ut <- matrix(U[, , t], n, R)
			res <- E0[, , t] - Ut %*% Ld %*% t(Ut)
			tot <- tot + sum(res[ok[, , t]]^2)
		}
		tot
	}
	lg <- function(r) {
		if (!is.finite(r) || abs(r) >= bound) return(-Inf)
		U <- .lame_uv_unit_paths(eps, r)
		if (is.finite(max_abs) && max(abs(U)) > max_abs) return(-Inf)
		-ss(U) / (2 * s2) + stats::dnorm(r, prior_mean, prior_sd, log = TRUE)
	}
	sl <- .lame_nc_slice(lg, rho, -bound, bound, w = 0.03)
	list(rho = sl$x, U = .lame_uv_unit_paths(eps, sl$x), evals = sl$evals)
}

# scale ridge of a bipartite fit: both sides' paths and sigma times c, G divided
# by c^2. `shape`, `scale` the inverse-gamma prior of sigma^2, `sigma_max`,
# `w_cap` (stationary variance), `g_cap` and `max_abs` the support bounds,
# `g_var` the prior variance of an element of G
.lame_uv_bip_scale_step <- function(
	U_cube,
	V_cube,
	G,
	rho,
	sigma,
	shape,
	scale,
	sigma_max = Inf,
	w_cap = Inf,
	g_var = 1,
	g_cap = Inf,
	max_abs = Inf
) {
	still <- list(U = U_cube, V = V_cube, G = G, sigma = sigma, c = 1, evals = 0L)
	if (!is.finite(sigma) || sigma <= 0 || !is.finite(rho) || abs(rho) >= 1 || length(G) == 0L || !all(is.finite(G))) return(still)
	dG <- length(G)
	umax <- max(abs(c(U_cube, V_cube)))
	g2 <- sum(G^2)
	gmax <- max(abs(G))
	om <- 1 - rho * rho
	lg <- function(x) {
		sg <- sigma * exp(x)
		if (sg > sigma_max || sg * sg / om > w_cap) return(-Inf)
		if (is.finite(g_cap) && gmax * exp(-2 * x) > g_cap) return(-Inf)
		if (is.finite(max_abs) && umax * exp(x) > max_abs) return(-Inf)
		-(2 * shape + 1) * log(sg) - scale / (sg * sg) + (1 - 2 * dG) * x - g2 * exp(-4 * x) / (2 * g_var)
	}
	sl <- .lame_nc_slice(lg, 0, -Inf, Inf, w = 0.25)
	cc <- exp(sl$x)
	if (!is.finite(cc) || cc == 1) { still$evals <- sl$evals; return(still) }
	list(U = U_cube * cc, V = V_cube * cc, G = G / (cc * cc), sigma = sigma * cc, c = cc, evals = sl$evals)
}

# the same draw when the latent array enters through a pointwise factor
# (Poisson): all latent values are carried by the change of the product, so no
# Gaussian residual changes and the conditional of rho is its prior times the
# Poisson factor of the observed cells. `Z` the latent array (mirrored when
# symmetric), `cx` the context of .lame_ltm_context()
.lame_uv_sym_nc_rho_carry <- function(
	U_cube,
	L,
	Z,
	cx,
	rho,
	prior_mean = 0.9,
	prior_sd = 0.1,
	max_abs = Inf,
	bound = 0.99
) {
	d <- dim(U_cube)
	n <- d[1L]; R <- d[2L]; N <- d[3L]
	still <- list(rho = rho, U = U_cube, Z = Z, evals = 0L)
	if (is.null(N) || N < 2L || !is.finite(rho) || abs(rho) >= bound || !all(is.finite(Z))) return(still)
	s <- sqrt(1 - rho * rho)
	eps <- U_cube
	for (t in 2:N) eps[, , t] <- (U_cube[, , t] - rho * U_cube[, , t - 1L]) / s
	Ld <- diag(L, R)
	prod_of <- function(U) {
		P <- array(0, c(n, n, N))
		for (t in seq_len(N)) { Ut <- matrix(U[, , t], n, R); P[, , t] <- Ut %*% Ld %*% t(Ut) }
		P
	}
	P0 <- prod_of(U_cube)
	lg <- function(r) {
		if (!is.finite(r) || abs(r) >= bound) return(-Inf)
		U <- .lame_uv_unit_paths(eps, r)
		if (is.finite(max_abs) && max(abs(U)) > max_abs) return(-Inf)
		z <- Z[cx$ok_idx] + (prod_of(U) - P0)[cx$ok_idx]
		sum(cx$y_ok * z - exp(z + cx$off_ok)) + stats::dnorm(r, prior_mean, prior_sd, log = TRUE)
	}
	sl <- .lame_nc_slice(lg, rho, -bound, bound, w = 0.03)
	if (sl$x == rho) { still$evals <- sl$evals; return(still) }
	U <- .lame_uv_unit_paths(eps, sl$x)
	list(rho = sl$x, U = U, Z = Z + (prod_of(U) - P0), evals = sl$evals)
}
