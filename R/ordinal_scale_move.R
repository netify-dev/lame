# exact rank-preserving moves along the scale of the latent array of a data-induced
# ordinal fit. the rank windows of Z move by about 1 / cells per sweep while the
# posterior of the scale has width 1 / sqrt(cells), so the parameters tied to that
# scale mix slowly. each sweep rescales the whole state by c > 0 (a group move,
# Liu and Wu 1999) and shifts and rescales each period's latent array, Z_t -> c Z_t + d.
# the metric is that of the sampler's Gaussian conditionals: dyads (e_ij, e_ji) with
# variance s2 and correlation rho when directed, one value per unordered pair when
# symmetric, independent cells when bipartite

# the number of cells of a period's grid
.lame_osm_ncell <- function(nr, nc, mode) {
	switch(mode, dir = nr * nr, sym = nr * (nr - 1) / 2, bip = nr * nc)
}

# inner product of two matrices in the metric of the covariance of the latent array;
# `up` the upper-triangle mask (square networks)
.lame_osm_dot <- function(x, y, rho, s2, mode, up = NULL) {
	if (identical(mode, "bip")) return(sum(x * y) / s2)
	if (is.null(up)) up <- upper.tri(x)
	if (identical(mode, "sym")) {
		xs <- (x + t(x)) / 2
		ys <- (y + t(y)) / 2
		return(sum(xs[up] * ys[up]) / s2)
	}
	xt <- t(x)
	yt <- t(y)
	sum(((x + xt) * (y + yt))[up]) / (2 * (1 + rho) * s2) +
		sum(((x - xt) * (y - yt))[up]) / (2 * (1 - rho) * s2) +
		sum(diag(x) * diag(y)) / ((1 + rho) * s2)
}

# a draw from the density c^(D - 1) exp(-A c^2 / 2 + B c) on c > 0: exact (gamma) when B
# is zero, otherwise an independence Metropolis-Hastings step from the Laplace
# approximation of the density, started at c = 1 (the current state)
.lame_osm_draw_c <- function(D, A, B) {
	if (!is.finite(D) || !is.finite(A) || !is.finite(B) || A <= 0 || D <= 1) return(1)
	if (B == 0) return(sqrt(stats::rgamma(1L, D / 2, A / 2)))
	cm <- (B + sqrt(B^2 + 4 * A * (D - 1))) / (2 * A)
	curv <- (D - 1) / cm^2 + A
	cp <- stats::rnorm(1L, cm, 1 / sqrt(curv))
	if (!(cp > 0)) return(1)
	lf <- function(x) (D - 1) * log(x) - A * x^2 / 2 + B * x
	lq <- function(x) -0.5 * curv * (x - cm)^2
	if (log(stats::runif(1L)) < lf(cp) - lf(1) + lq(1) - lq(cp)) cp else 1
}

# the shift and the scale of the latent array of one period given its predictor
.lame_osm_period <- function(Zt, Et, rho, s2, mode) {
	nr <- nrow(Zt)
	nc <- ncol(Zt)
	up <- if (identical(mode, "bip")) NULL else upper.tri(Zt)
	ncell <- .lame_osm_ncell(nr, nc, mode)
	if (identical(mode, "sym")) {
		Zt <- (Zt + t(Zt)) / 2
		Et <- (Et + t(Et)) / 2
	}
	# shift: the density of Z + d around E is Gaussian in d
	res <- Zt - Et
	if (identical(mode, "sym")) {
		prec <- ncell / s2
		mu <- -sum(res[up]) / ncell
	} else if (identical(mode, "bip")) {
		prec <- ncell / s2
		mu <- -sum(res) / ncell
	} else {
		prec <- ncell / ((1 + rho) * s2)
		mu <- -sum(res) / ncell
	}
	Zt <- Zt + stats::rnorm(1L, mu, 1 / sqrt(prec))
	# scale
	A <- .lame_osm_dot(Zt, Zt, rho, s2, mode, up)
	B <- .lame_osm_dot(Zt, Et, rho, s2, mode, up)
	cc <- .lame_osm_draw_c(ncell, A, B)
	if (cc != 1) Zt <- cc * Zt
	Zt
}

# reference visit written directly on the latent array (the definition that
# .lame_osm_visit computes from one-pass sums); arguments as for .lame_osm_visit
.lame_osm_visit_ref <- function(Z, E, m, mode, rho, s2, A_pri, B_pri, d_th, affine = TRUE, orbit = TRUE) {
	if (!all(is.finite(Z)) || !all(is.finite(E))) return(list(Z = Z, c = 1, E = E))
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	nt <- dim(Z)[3L]
	rho <- rep_len(rho, nt)
	up <- if (identical(mode, "bip")) NULL else upper.tri(matrix(0, nr, nc))
	ncell <- .lame_osm_ncell(nr, nc, mode)
	if (identical(mode, "sym")) {
		sym_t <- function(a) {
			for (tt in seq_len(dim(a)[3L])) a[, , tt] <- (a[, , tt] + t(a[, , tt])) / 2
			a
		}
		E <- sym_t(E)
		if (!is.null(m)) m <- sym_t(m)
	}
	# the orbit of the whole state: u = Z - E_scaled, the unscaled part m
	cc <- 1
	if (orbit) {
		Es <- if (is.null(m)) E else E - m
		A <- A_pri
		B <- B_pri
		for (tt in seq_len(nt)) {
			u <- Z[, , tt] - Es[, , tt]
			A <- A + .lame_osm_dot(u, u, rho[tt], s2, mode, up)
			if (!is.null(m)) B <- B + .lame_osm_dot(u, m[, , tt], rho[tt], s2, mode, up)
		}
		cc <- .lame_osm_draw_c(nt * ncell + d_th, A, B)
		if (cc != 1) {
			Z <- cc * Z
			E <- if (is.null(m)) cc * E else cc * Es + m
		}
	}
	if (affine) {
		for (tt in seq_len(nt)) Z[, , tt] <- .lame_osm_period(Z[, , tt], E[, , tt], rho[tt], s2, mode)
	}
	list(Z = Z, c = cc, E = E)
}


# one visit. `Z` the latent array (nr x nc x T); `E` the predictor of the state and `m` the part the
# scale does not carry (an array or NULL), both averaged over the pair when symmetric; `mode` "dir", "sym"
# or "bip"; `rho` the dyadic correlation per period (length T or one value); `s2` the dyadic variance;
# `A_pri`, `B_pri`, `d_th` the quadratic term, linear term and number of the scaled parameters under
# their priors. `affine` adds each period's shift and scale; `orbit = FALSE` skips the scale of the state.
# returns the updated latent array and the factor c for the scaled parameters. every quantity is a
# polynomial in c and d of the Gram sums of Z, E and m, taken in one pass (osm_gram_cpp)
.lame_osm_visit <- function(Z, E, m, mode, rho, s2, A_pri, B_pri, d_th, affine = TRUE, orbit = TRUE) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	nt <- dim(Z)[3L]
	rho <- rep_len(rho, nt)
	G <- osm_gram_cpp(
		Z,
		E,
		if (is.null(m)) numeric(0) else m,
		c(nr, nc, nt),
		switch(mode, dir = 1L, sym = 2L, bip = 3L)
	)
	if (!all(is.finite(G))) return(list(Z = Z, c = 1))
	ncell <- .lame_osm_ncell(nr, nc, mode)
	dir <- identical(mode, "dir")
	# the metric: <x, y> = sum over the three forms of (form sum x coefficient); <1, x> = (sum of x) ko
	kf <- if (dir) cbind(
		1 / (2 * (1 + rho) * s2),
		1 / (2 * (1 - rho) * s2),
		1 / ((1 + rho) * s2)
	) else cbind(rep(1 / s2, nt), 0, 0)
	ko <- if (dir) 1 / ((1 + rho) * s2) else rep(1 / s2, nt)
	dotm <- function(p) rowSums(G[, 3L * (p - 1L) + 1:3, drop = FALSE] * kf)
	ZZ <- dotm(1L)
	ZE <- dotm(2L)
	ZM <- dotm(3L)
	EE <- dotm(4L)
	EM <- dotm(5L)
	MM <- dotm(6L)
	sZ <- G[, 19L]
	sE <- G[, 20L]
	sM <- G[, 21L]
	# the orbit of the whole state: u = Z - E + m, the unscaled part m
	cc <- 1
	if (orbit) {
		A <- A_pri + sum(ZZ + EE + MM - 2 * ZE + 2 * ZM - 2 * EM)
		B <- B_pri + sum(ZM - EM + MM)
		cc <- .lame_osm_draw_c(nt * ncell + d_th, A, B)
	}
	a_t <- rep(cc, nt)
	b_t <- numeric(nt)
	if (affine) {
		for (tt in seq_len(nt)) {
			# after the state scale the latent array is c Z and the predictor c E + (1 - c) m; the shift is a normal draw
			mu <- -(cc * (sZ[tt] - sE[tt]) - (1 - cc) * sM[tt]) / ncell
			d <- stats::rnorm(1L, mu, 1 / sqrt(ncell * ko[tt]))
			A2 <- cc^2 * ZZ[tt] + 2 * cc * d * sZ[tt] * ko[tt] + d^2 * ncell * ko[tt]
			B2 <- cc^2 * ZE[tt] + cc * (1 - cc) * ZM[tt] + d * cc * sE[tt] * ko[tt] + d * (1 - cc) * sM[tt] * ko[tt]
			ct <- .lame_osm_draw_c(ncell, A2, B2)
			a_t[tt] <- ct * cc
			b_t[tt] <- ct * d
		}
	}
	if (any(a_t != 1) || any(b_t != 0)) {
		cells <- nr * nc
		Z <- Z * array(rep(a_t, each = cells), c(nr, nc, nt)) + array(rep(b_t, each = cells), c(nr, nc, nt))
	}
	list(Z = Z, c = cc)
}

# quadratic form and number of coordinates of the additive effects under N(0, Sab): actor pairs when
# directed, one effect per actor when symmetric, separate rows and columns when bipartite. a margin
# not in the model is held at zero; NULL when the covariance cannot be inverted
.lame_osm_prior_ab <- function(a, b, Sab, rvar, cvar, symmetric = FALSE, bip = FALSE) {
	if (isTRUE(symmetric)) {
		if (!isTRUE(rvar) || !(Sab[1L, 1L] > 0)) return(c(0, 0))
		return(c(sum(a^2) / Sab[1L, 1L], length(a)))
	}
	if (isTRUE(bip) || !(isTRUE(rvar) && isTRUE(cvar))) {
		A <- 0
		d <- 0
		if (isTRUE(rvar)) {
			if (!(Sab[1L, 1L] > 0)) return(NULL)
			A <- A + sum(a^2) / Sab[1L, 1L]
			d <- d + length(a)
		}
		if (isTRUE(cvar)) {
			if (!(Sab[2L, 2L] > 0)) return(NULL)
			A <- A + sum(b^2) / Sab[2L, 2L]
			d <- d + length(b)
		}
		return(c(A, d))
	}
	Si <- tryCatch(solve(Sab), error = function(e) NULL)
	if (is.null(Si) || any(!is.finite(Si))) return(NULL)
	ab <- cbind(a, b)
	c(sum((ab %*% Si) * ab), 2 * length(a))
}

# priors of the scaled parameters of lame(): quadratic term A, linear term B and count d, or NULL when
# Sab cannot be inverted. static coefficients take the sampler's g-prior (`P0` unipartite, I / (g s2)
# bipartite), dynamic coefficients their path prior (as in R/poisson_period_move.R), static additive
# effects N(0, Sab); dynamic additive effects stay in the unscaled part of the predictor
.lame_osm_prior_lame <- function(
	beta_dyn,
	beta,
	beta_dyn_path,
	beta_static,
	P0,
	g,
	s2,
	a,
	b,
	Sab,
	rvar,
	cvar,
	symmetric,
	bip,
	dynamic_ab,
	kind,
	rho_beta,
	sigma_beta,
	lam,
	m0,
	v0,
	tpos,
	length_scale,
	nt
) {
	A <- 0
	B <- 0
	d <- 0
	if (isTRUE(beta_dyn$any)) {
		for (j in seq_along(beta_dyn$dynamic_idx)) {
			pri <- .lame_ppm_path_prior(kind, nt, rho_beta[j], sigma_beta[j], lam[j], tpos, length_scale, m0[j], v0[j])
			x <- beta_dyn_path[, j]
			A <- A + sum(x * drop(pri$P %*% x))
			B <- B + sum(pri$b * x)
			d <- d + nt
		}
		if (length(beta_static)) {
			A <- A + sum(beta_static^2 / (.lame_static_g_scale(g, beta_dyn$static_idx) * s2))
			d <- d + length(beta_static)
		}
	} else if (length(beta)) {
		A <- A + if (isTRUE(bip)) sum(beta^2) / (as.numeric(g[1L]) * s2) else drop(crossprod(beta, P0 %*% beta))
		d <- d + length(beta)
	}
	if (!isTRUE(dynamic_ab)) {
		ab <- .lame_osm_prior_ab(a, b, Sab, rvar, cvar, symmetric = symmetric, bip = bip)
		if (is.null(ab)) return(NULL)
		A <- A + ab[1L]
		d <- d + ab[2L]
	}
	if (!all(is.finite(c(A, B, d)))) return(NULL)
	list(A = A, B = B, d = d)
}

# whether every observed latent value lies strictly inside its explicit-cutpoint interval, so a
# rescaled state is accepted only if rounding left no value on a boundary
.lame_osm_ex_ok <- function(Z, Y_int, alpha) {
	obs <- !is.na(Y_int)
	if (!any(obs)) return(TRUE)
	bounds <- c(-Inf, alpha, Inf)
	code <- Y_int[obs]
	all(Z[obs] > bounds[code] & Z[obs] < bounds[code + 1L])
}
