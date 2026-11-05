# AR(1) conditional of a dynamic effect at period t given its neighbours.
# the effect is a per-actor level plus a zero-mean AR(1) deviation, so the conditional is formed on the deviations
.lame_ab_ar1_cond <- function(mat, t, N, rho, sigma, level = NULL) {
	v_innov <- sigma * sigma
	if (!is.finite(v_innov) || v_innov <= 0 || !is.finite(rho) || abs(rho) >= 1) stop("Invalid additive AR parameters.")
	n <- nrow(mat)
	lev <- if (is.null(level)) rep(0, n) else level
	lev[!is.finite(lev)] <- 0
	if (N == 1L) {
		return(list(m = lev, v = v_innov / (1 - rho^2)))
	}
	prev <- if (t > 1L) mat[, t - 1L] - lev else rep(0, n)
	nxt <- if (t < N) mat[, t + 1L] - lev else rep(0, n)
	prev[!is.finite(prev)] <- 0
	nxt[!is.finite(nxt)] <- 0
	if (t == 1L) {
		list(m = lev + rho * nxt, v = v_innov)
	} else if (t == N) {
		list(m = lev + rho * prev, v = v_innov)
	} else {
		list(
			m = lev + rho * (prev + nxt) / (1 + rho * rho),
			v = v_innov / (1 + rho * rho)
		)
	}
}
# sufficient statistics of one margin's per-actor levels from its AR(1) paths;
# the caller adds the prior jointly over each actor's pair
.lame_ab_level_path_stats <- function(mat, rho, sigma) {
	n <- nrow(mat); N <- ncol(mat)
	if (!is.finite(sigma) || sigma <= 0 || !is.finite(rho) || abs(rho) >= 1) stop("Invalid additive AR parameters.")
	s2i <- 1 / (sigma * sigma)
	one_m_rho2 <- 1 - rho * rho
	prec <- numeric(n)
	num <- numeric(n)
	x1 <- mat[, 1L]
	ok1 <- is.finite(x1)
	prec[ok1] <- prec[ok1] + one_m_rho2 * s2i
	num[ok1] <- num[ok1] + one_m_rho2 * s2i * x1[ok1]
	if (N >= 2L) {
		innov <- mat[, -1L, drop = FALSE] - rho * mat[, -N, drop = FALSE]
		okt <- is.finite(innov)
		prec <- prec + rowSums(okt) * (1 - rho)^2 * s2i
		innov[!okt] <- 0
		num <- num + (1 - rho) * s2i * rowSums(innov)
	}
	list(num = num, prec = prec)
}
# draws x_i ~ N(m_i, v_i) conditioned on sum(x) = 0
.lame_ab_center_weighted <- function(x, v) {
	ok <- is.finite(x) & is.finite(v) & v > 0
	if (!any(ok)) return(x)
	x[ok] <- x[ok] - v[ok] * (sum(x[ok]) / sum(v[ok]))
	x
}
# joint draw of the row and column levels, each actor's pair a priori N(0, Sab).
# with `center` the draw is conditioned on each margin summing to zero; a symmetric network has one level per actor
.lame_ab_level_joint <- function(
	num_a,
	prec_a,
	num_b,
	prec_b,
	Sab,
	symmetric = FALSE,
	center = TRUE
) {
	n <- length(num_a)
	if (n == 0L) return(list(a = numeric(0), b = numeric(0)))
	if (isTRUE(symmetric)) {
		d <- prec_a + 1 / Sab[1, 1]
		if (any(!is.finite(d) | d <= 0)) stop("Invalid additive level precision.")
		v <- 1 / d
		x <- v * num_a + sqrt(v) * rnorm(n)
		if (center) x <- .lame_ab_center_weighted(x, v)
		x[!is.finite(x)] <- 0
		return(list(a = x, b = x))
	}
	if (length(num_a) != length(num_b)) {
		# bipartite: row and column actors are different sets, so each margin is drawn on its own
		draw_one <- function(num, prec, v_prior) {
			d <- prec + 1 / v_prior
			if (any(!is.finite(d) | d <= 0)) stop("Invalid additive level precision.")
			v <- 1 / d
			x <- v * num + sqrt(v) * rnorm(length(num))
			if (center) x <- .lame_ab_center_weighted(x, v)
			x[!is.finite(x)] <- 0
			x
		}
		return(list(
			a = draw_one(num_a, prec_a, Sab[1, 1]),
			b = draw_one(num_b, prec_b, Sab[2, 2])
		))
	}
	iS <- tryCatch(solve(Sab), error = function(e) NULL)
	if (is.null(iS) || !all(is.finite(iS))) {
		stop("Unable to resolve the additive-effect covariance precision.")
	}
	q11 <- prec_a + iS[1, 1]
	q22 <- prec_b + iS[2, 2]
	q12 <- iS[1, 2]
	det <- q11 * q22 - q12 * q12
	bad <- !is.finite(det) | det <= 0
	if (any(bad)) {
		stop("Unable to resolve the coupled additive level precision.")
	}
	idet <- 1 / det
	v11 <- q22 * idet; v22 <- q11 * idet; v12 <- (-q12) * idet
	# lower Cholesky factor of each actor's 2 x 2 covariance block
	l11 <- sqrt(v11)
	l21 <- v12 / l11
	l22_sq <- v22 - l21 * l21
	# Schur complement form of the variance avoids cancellation at the boundary
	if (any(l22_sq <= 0)) l22_sq[l22_sq <= 0] <- 1 / q22[l22_sq <= 0]
	l22 <- sqrt(l22_sq)
	z1 <- rnorm(n); z2 <- rnorm(n)
	xa <- v11 * num_a + v12 * num_b + l11 * z1
	xb <- v12 * num_a + v22 * num_b + l21 * z1 + l22 * z2
	if (center) {
		# condition on both margins summing to zero
		s11 <- sum(v11); s12 <- sum(v12); s22 <- sum(v22)
		dS <- s11 * s22 - s12 * s12
		if (is.finite(dS) && abs(dS) > 1e-12) {
			ra <- sum(xa); rb <- sum(xb)
			t1 <- (s22 * ra - s12 * rb) / dS
			t2 <- (s11 * rb - s12 * ra) / dS
			xa <- xa - (v11 * t1 + v12 * t2)
			xb <- xb - (v12 * t1 + v22 * t2)
		} else {
			constraint <- matrix(c(s11, s12, s12, s22), 2L)
			t <- solve(constraint, c(sum(xa), sum(xb)))
			xa <- xa - v11 * t[1L] - v12 * t[2L]
			xb <- xb - v12 * t[1L] - v22 * t[2L]
		}
	}
	if (any(!is.finite(c(xa, xb)))) stop("Nonfinite joint additive level draw.")
	list(a = xa, b = xb)
}
# per-actor levels of the dynamic additive effects given their paths
.lame_ab_level_draw <- function(
	a_mat,
	b_mat,
	rho,
	sigma,
	Sab,
	symmetric = FALSE,
	center = TRUE,
	rvar = TRUE,
	cvar = TRUE
) {
	sa <- .lame_ab_level_path_stats(a_mat, rho, sigma)
	sb <- if (isTRUE(symmetric)) sa else .lame_ab_level_path_stats(b_mat, rho, sigma)
	if (!isTRUE(symmetric) && (!isTRUE(rvar) || !isTRUE(cvar))) {
		draw_one <- function(st, variance, active) {
			if (!isTRUE(active)) return(numeric(length(st$num)))
			v <- 1 / (st$prec + 1 / variance)
			x <- v * st$num + sqrt(v) * rnorm(length(v))
			if (isTRUE(center)) x <- .lame_ab_center_weighted(x, v)
			x
		}
		return(list(a = draw_one(sa, Sab[1, 1], rvar), b = draw_one(sb, Sab[2, 2], cvar)))
	}
	if (isTRUE(symmetric) && !isTRUE(rvar)) return(list(a = numeric(nrow(a_mat)), b = numeric(nrow(b_mat))))
	.lame_ab_level_joint(
		sa$num,
		sa$prec,
		sb$num,
		sb$prec,
		Sab,
		symmetric = symmetric,
		center = center
	)
}
# conditional prior of one margin's levels given the other margin's under (a_i, b_i) ~ N(0, Sab);
# the marginal is used when Sab carries no usable covariance
.lame_ab_cond_prior <- function(Sab, margin, other_lev) {
	k <- margin; o <- 3L - margin
	v_k <- Sab[k, k]; v_o <- Sab[o, o]; c_ko <- Sab[k, o]
	if (any(!is.finite(c(v_k, v_o, c_ko))) || v_k <= 0 || v_o < 0)
		stop("Invalid additive-effect covariance.")
	# a disabled margin is a point mass at zero and leaves the active prior marginal
	if (v_o == 0 && (c_ko != 0 || any(!is.finite(other_lev)) || any(other_lev != 0)))
		stop("Invalid disabled additive-effect margin.")
	if (c_ko == 0 || length(other_lev) == 0L) return(list(mean = NULL, var = v_k))
	variance <- v_k - c_ko * c_ko / v_o
	if (!is.finite(variance) || variance <= 0 || any(!is.finite(other_lev)))
		stop("Invalid conditional additive-effect covariance or levels.")
	list(mean = (c_ko / v_o) * other_lev, var = variance)
}

# sufficient statistics for the per-actor levels of a dynamic additive effect given the deviations and the data.
# with a dyadic correlation the two cells of a dyad are read together; the unipartite diagonal is skipped
.lame_ab_level_stats <- function(
	resid,
	dev,
	other,
	s2,
	margin = 1L,
	rho = 0,
	skip_diag = FALSE
) {
	dims <- dim(resid)
	N <- dims[3L]
	if (length(rho) > 1L) return(.lame_ab_level_stats_period(resid, dev, other, s2, margin, rho, skip_diag))
	n_this <- if (margin == 1L) dims[1L] else dims[2L]
	square <- dims[1L] == dims[2L]
	rho <- if (square) rho else 0
	if (!is.finite(rho) || abs(rho) >= 1) stop("Residual correlation must lie strictly between -1 and 1.")
	if (!is.finite(s2) || s2 <= 0) stop("Residual variance must be finite and positive.")
	dev[!is.finite(dev)] <- NA_real_
	other[!is.finite(other)] <- 0
	if (rho == 0) {
		# uncorrelated dyads: every observed cell is one observation of its level
		m <- resid -
			aperm(
				array(dev, c(nrow(dev), N, if (margin == 1L) dims[2L] else dims[1L])),
				if (margin == 1L) c(1L, 3L, 2L) else c(3L, 1L, 2L)
			) -
			aperm(
				array(other, c(nrow(other), N, if (margin == 1L) dims[1L] else dims[2L])),
				if (margin == 1L) c(3L, 1L, 2L) else c(1L, 3L, 2L)
			)
		if (skip_diag && square) {
			m[cbind(rep(seq_len(dims[1L]), N), rep(seq_len(dims[1L]), N), rep(seq_len(N), each = dims[1L]))] <- NA_real_
		}
		ok <- is.finite(m)
		m[!ok] <- 0
		if (margin == 1L) {
			num <- rowSums(m, dims = 1L) / s2
			prec_d <- rowSums(ok, dims = 1L) / s2
		} else {
			num <- rowSums(colSums(m)) / s2
			prec_d <- rowSums(colSums(ok)) / s2
		}
		return(list(num = num, prec_d = prec_d, prec_off = NULL))
	}
	# correlated dyads: fully observed pairs couple the two levels, half-observed pairs enter with 1 / s2
	q <- 1 / (s2 * (1 - rho * rho))
	p <- -rho / (s2 * (1 - rho * rho))
	num <- numeric(n_this)
	prec_d <- numeric(n_this)
	prec_off <- matrix(0, n_this, n_this)
	for (tt in seq_len(N)) {
		m <- resid[, , tt]
		d_t <- dev[, tt]
		o_t <- other[, tt]
		m <- if (margin == 1L) {
			m - d_t - matrix(o_t, dims[1L], dims[2L], byrow = TRUE)
		} else {
			m - matrix(o_t, dims[1L], dims[2L]) - matrix(d_t, dims[1L], dims[2L], byrow = TRUE)
		}
		if (skip_diag && square) diag(m) <- NA_real_
		ok <- is.finite(m)
		m[!ok] <- 0
		both <- ok & t(ok)
		w <- q * both + (ok & !both) / s2
		if (margin == 1L) {
			num <- num + rowSums(m * w) + p * rowSums(t(m) * both)
			prec_d <- prec_d + rowSums(w)
		} else {
			num <- num + colSums(m * w) + p * colSums(t(m) * both)
			prec_d <- prec_d + colSums(w)
		}
		prec_off <- prec_off + p * both
	}
	list(num = num, prec_d = prec_d, prec_off = prec_off)
}
# conjugate draw of the per-actor levels from the statistics above and their conditional prior:
# independent normals for uncorrelated dyads, one joint draw from the pair-coupled precision otherwise
.lame_ab_level_draw_lik <- function(
	resid,
	dev,
	other,
	var_level,
	s2,
	margin = 1L,
	rho = 0,
	skip_diag = FALSE,
	mean_level = NULL,
	center = FALSE
) {
	st <- .lame_ab_level_stats(resid, dev, other, s2, margin, rho, skip_diag)
	vl <- var_level
	if (!is.finite(vl) || vl <= 0) stop("Additive prior variance must be finite and positive.")
	prec_d <- st$prec_d + 1 / vl
	num <- st$num
	if (!is.null(mean_level)) {
		ml <- mean_level; ml[!is.finite(ml)] <- 0
		num <- num + ml / vl
	}
	n_this <- length(prec_d)
	lev <- NULL
	if (!is.null(st$prec_off)) {
		Q <- st$prec_off
		diag(Q) <- prec_d
		L <- chol(Q)
		if (!is.null(L)) {
			mu <- backsolve(L, forwardsolve(t(L), num))
			lev <- as.numeric(mu + backsolve(L, rnorm(n_this)))
			if (center) {
				# condition on sum(lev) = 0
				u <- backsolve(L, forwardsolve(t(L), rep(1, n_this)))
				su <- sum(u)
				if (is.finite(su) && su > 0) lev <- lev - u * (sum(lev) / su)
			}
		}
	}
	if (is.null(lev)) {
		lev <- rnorm(n_this, num / prec_d, 1 / sqrt(prec_d))
		if (center) lev <- .lame_ab_center_weighted(lev, 1 / prec_d)
	}
	if (any(!is.finite(lev))) stop("Nonfinite additive level draw.")
	lev
}
# centre the dynamic additive effects in every period and return the removed constants,
# so the caller can hand them to the intercept (identification against the intercept)
.lame_ab_center_dynamic <- function(a_mat, b_mat) {
	ma <- colMeans(a_mat, na.rm = TRUE)
	mb <- colMeans(b_mat, na.rm = TRUE)
	ma[!is.finite(ma)] <- 0
	mb[!is.finite(mb)] <- 0
	a_mat <- a_mat - rep(ma, each = nrow(a_mat))
	b_mat <- b_mat - rep(mb, each = nrow(b_mat))
	list(a_mat = a_mat, b_mat = b_mat, shift = ma + mb)
}
# prior scale of the deviation innovation variance: the prior mean of the level variance times (1 - rho^2),
# so a deviation path's stationary variance matches that of the levels. it is fixed before sampling
.lame_sigma_ab_prior_scale <- function(Sab0, eta0, rho, rvar, cvar) {
	v <- c(if (isTRUE(rvar)) Sab0[1, 1], if (isTRUE(cvar)) Sab0[2, 2])
	v <- v[is.finite(v) & v > 0]
	if (!length(v)) return(1)
	# prior mean of the inverse Wishart; it exists only for eta0 > 3
	mult <- if (is.finite(eta0) && eta0 > 3) eta0 / (eta0 - 3) else 1
	# floor relative to the scale of the effects
	max(mult * mean(v) * (1 - min(rho * rho, 0.9801)), 1e-6 * mean(v))
}
# residual variance of the two-way-centred response on the model's latent scale.
# probit and rank families fix that scale at one; poisson reads it on the log scale
.lame_response_resid_var <- function(Y, family) {
	if (!(family %in% c("normal", "poisson"))) return(1)
	A <- if (identical(family, "poisson")) log1p(pmax(as.array(Y), 0)) else as.array(Y)
	d <- dim(A)
	if (is.null(d) || length(d) < 2L) return(1)
	slices <- if (length(d) >= 3L) seq_len(d[3L]) else 1L
	res <- unlist(lapply(slices, function(t) {
		M <- if (length(d) >= 3L) A[, , t] else A
		M <- as.matrix(M)
		r <- rowMeans(M, na.rm = TRUE); r[!is.finite(r)] <- 0
		cc <- colMeans(M, na.rm = TRUE); cc[!is.finite(cc)] <- 0
		g <- mean(M, na.rm = TRUE); if (!is.finite(g)) g <- 0
		M - outer(r, cc, "+") + g
	}), use.names = FALSE)
	v <- stats::var(res, na.rm = TRUE)
	if (!is.finite(v) || v <= 0) 1 else v
}
# prior scale of the dynamic-coefficient innovation variance, on the response scale.
# the default puts the prior mean of the innovation sd at a fifth of a residual sd per period
.lame_sigma_beta_prior_scale <- function(Y, family, shape = 2, frac = 0.2) {
	unit <- .lame_response_resid_var(Y, family)
	mult <- if (is.finite(shape) && shape > 1) shape - 1 else 1
	# floor relative to the residual scale
	max(mult * frac * frac * unit, 1e-8 * unit)
}
# independence MH step for rho_uv under the stationary AR(1) initial condition of the dynamic UV kernels.
# the proposal is the truncated conjugate normal; the ratio supplies the t = 1 stationary factor it omits
.lame_rho_uv_mh_stationary <- function(
	rho_current,
	mean_post,
	var_post,
	S1,
	n1,
	sigma_uv,
	stat_var_cap = Inf
) {
	lo <- -0.99
	hi <- 0.99
	if (!is.finite(mean_post) || !is.finite(var_post) || var_post <= 0) {
		return(rho_current)
	}
	sd_post <- sqrt(var_post)
	prop <- .lame_ar_truncnorm_draw(mean_post, sd_post, lo, hi)
	rc <- min(max(rho_current, lo), hi)
	s2i <- 1 / (sigma_uv * sigma_uv)
	log_f1 <- function(r) {
		0.5 * n1 * log(1 - r * r) - 0.5 * (1 - r * r) * S1 * s2i
	}
	la <- log_f1(prop) - log_f1(rc)
	candidate <- if (is.finite(la) && log(runif(1)) < la) prop else rho_current
	.lame_ar_rho_support(candidate, rho_current, sigma_uv, stat_var_cap)
}
# scale of the position a snap reset lands on: prior$snap_kappa when set,
# otherwise the stationary sd of the position paths, sigma_uv / sqrt(1 - rho_uv^2)
.lame_snap_kappa <- function(kappa_fixed, rho_uv, sigma_uv) {
	if (!is.null(kappa_fixed)) return(kappa_fixed)
	sigma_uv / sqrt(1 - min(rho_uv * rho_uv, 0.9801))
}
# sum of squares and coordinate count of the positions drawn at a reset (delta = 1, periods 2 onward)
.lame_snap_reset_ss <- function(U_cube, V_cube, symmetric, delta_u, delta_v = NULL) {
	n_time <- dim(U_cube)[3]
	R <- dim(U_cube)[2]
	ss <- 0
	cnt <- 0
	if (is.null(n_time) || n_time < 2 || is.null(delta_u)) return(list(ss = 0, cnt = 0))
	for (t in 2:n_time) {
		hit <- delta_u[, t] > 0.5
		if (any(hit)) {
			ss <- ss + sum(U_cube[hit, , t]^2)
			cnt <- cnt + sum(hit) * R
		}
		if (!symmetric && !is.null(delta_v)) {
			hit <- delta_v[, t] > 0.5
			if (any(hit)) {
				ss <- ss + sum(V_cube[hit, , t]^2)
				cnt <- cnt + sum(hit) * dim(V_cube)[2]
			}
		}
	}
	list(ss = ss, cnt = cnt)
}
.lame_sample_rho_uv_drift_only <- function(
	U_cube,
	V_cube,
	sigma_uv,
	rho_current,
	symmetric,
	delta_u,
	delta_v = NULL,
	prior_mean = 0.9,
	prior_sd = 0.1,
	reset_stationary = FALSE,
	stat_var_cap = Inf
) {
	n_time <- dim(U_cube)[3]
	if (is.null(n_time) || n_time < 2 || is.null(delta_u)) return(rho_current)
	sigma2_inv <- 1 / (sigma_uv * sigma_uv)
	prior_prec <- 1 / (prior_sd * prior_sd)
	sum_yt_yt1 <- 0
	sum_yt1_yt1 <- 0
	for (t in 2:n_time) {
		keep_u <- delta_u[, t] == 0
		if (any(keep_u)) {
			U_t <- U_cube[keep_u, , t, drop = FALSE]
			U_t1 <- U_cube[keep_u, , t - 1, drop = FALSE]
			sum_yt_yt1 <- sum_yt_yt1 + sum(U_t * U_t1)
			sum_yt1_yt1 <- sum_yt1_yt1 + sum(U_t1 * U_t1)
		}
		if (!symmetric && !is.null(delta_v)) {
			keep_v <- delta_v[, t] == 0
			if (any(keep_v)) {
				V_t <- V_cube[keep_v, , t, drop = FALSE]
				V_t1 <- V_cube[keep_v, , t - 1, drop = FALSE]
				sum_yt_yt1 <- sum_yt_yt1 + sum(V_t * V_t1)
				sum_yt1_yt1 <- sum_yt1_yt1 + sum(V_t1 * V_t1)
			}
		}
	}
	if (!is.finite(sum_yt1_yt1) || sum_yt1_yt1 <= 0) return(rho_current)
	var_post <- 1 / (sum_yt1_yt1 * sigma2_inv + prior_prec)
	mean_post <- var_post * (sum_yt_yt1 * sigma2_inv + prior_mean * prior_prec)
	# stationary initial condition at t = 1, which always contributes
	S1 <- sum(U_cube[, , 1]^2)
	n1 <- dim(U_cube)[1] * dim(U_cube)[2]
	if (!symmetric && !is.null(delta_v)) {
		S1 <- S1 + sum(V_cube[, , 1]^2)
		n1 <- n1 + dim(V_cube)[1] * dim(V_cube)[2]
	}
	# under the linked snap scale every reset position is a stationary draw too
	if (isTRUE(reset_stationary)) {
		st <- .lame_snap_reset_ss(U_cube, V_cube, symmetric, delta_u, delta_v)
		S1 <- S1 + st$ss
		n1 <- n1 + st$cnt
	}
	.lame_rho_uv_mh_stationary(
		rho_current,
		mean_post,
		var_post,
		S1,
		n1,
		sigma_uv,
		stat_var_cap = stat_var_cap
	)
}
# counter for draws pushed to the bound of a stationary-variance cap, reported by the fit
.lame_cap_env <- new.env(parent = emptyenv())
# the stationary-variance caps are on the latent scale: relative to s2 for the normal family, absolute otherwise
.lame_cap_unit <- function(family, s2) {
	if (identical(family, "normal")) s2 else 1
}
# innovation scale of the dynamic additive effects, drawn exactly from the inverse gamma truncated so that
# sigma^2 / (1 - rho^2) stays at or below `stat_var_cap`, which keeps a sparse panel from separating.
# a non-finite cap disables it
.lame_sample_sigma_ab <- function(
	a_mat,
	b_mat,
	rho_ab,
	sigma_current,
	symmetric,
	prior_shape = 2,
	prior_scale = 1,
	stat_var_cap = Inf,
	center = FALSE,
	rvar = TRUE,
	cvar = TRUE
) {
	N <- ncol(a_mat)
	if (is.null(N) || N < 2L) return(sigma_current)
	one_m_rho2 <- 1 - min(rho_ab * rho_ab, 0.9801)
	paths <- list()
	if (isTRUE(rvar)) paths[[length(paths) + 1L]] <- a_mat
	if (!isTRUE(symmetric) && isTRUE(cvar)) paths[[length(paths) + 1L]] <- b_mat
	ss <- 0
	cnt <- 0L
	for (M in paths) {
		x1 <- M[, 1L]
		x1 <- x1[is.finite(x1)]
		ss <- ss + one_m_rho2 * sum(x1 * x1)
		cnt <- cnt + length(x1)
		innov <- M[, -1L, drop = FALSE] - rho_ab * M[, -N, drop = FALSE]
		innov <- innov[is.finite(innov)]
		ss <- ss + sum(innov * innov)
		cnt <- cnt + length(innov)
		if (isTRUE(center)) {
			if (any(!is.finite(M))) stop("Centered actor variance updates require finite latent paths.", call. = FALSE)
			cnt <- cnt - if (nrow(M) > 0L) N else 0L
		}
	}
	if (!is.finite(ss) || cnt < 0L) return(sigma_current)
	shape <- prior_shape + cnt / 2
	rate <- prior_scale + ss / 2
	sigma_max <- if (is.finite(stat_var_cap) && stat_var_cap > 0) {
		sqrt(stat_var_cap * one_m_rho2)
	} else Inf
	if (is.finite(sigma_max)) {
		return(.lame_rinvgamma_sd_trunc(shape, rate, sigma_max, sigma_current))
	}
	g <- rgamma(1, shape = shape, rate = rate)
	if (!is.finite(g) || g <= 0) stop("Invalid additive innovation variance draw.")
	sqrt(1 / g)
}
# truncated inverse-gamma draw in log survival space, which keeps tail resolution when the Gamma survival underflows
.lame_rinvgamma_sd_trunc <- function(shape, rate, sigma_max, current) {
	glo <- 1 / (sigma_max * sigma_max)
	log_hi <- pgamma(glo, shape = shape, rate = rate, lower.tail = FALSE, log.p = TRUE)
	.lame_cap_env$draws <- (.lame_cap_env$draws %||% 0L) + 1L
	if (!is.finite(log_hi)) return(current)
	if (log_hi < log(.5)) .lame_cap_env$hits <- (.lame_cap_env$hits %||% 0L) + 1L
	log_u <- log(runif(1)) + log_hi
	g <- qgamma(log_u, shape = shape, rate = rate, lower.tail = FALSE, log.p = TRUE)
	if (!is.finite(g) || g < glo) return(current)
	sqrt(1 / g)
}
# compare predictive-draw averages with the observed level for dynamic binary or poisson fits.
# returns TRUE when the level is off (and warns), FALSE when it holds, NA when the check does not apply
.lame_predictive_level_check <- function(family, dynamic, YPS, Y, n_stored) {
	# fewer than 100 stored draws do not pin down the mean of the replicated counts
	if (!family %in% c("binary", "poisson") || !isTRUE(dynamic) || is.null(YPS) ||
			!is.numeric(n_stored) || length(n_stored) != 1L || !is.finite(n_stored) || n_stored < 100) {
		return(NA)
	}
	obs_cells <- is.finite(Y)
	if (!any(obs_cells) || !identical(dim(YPS), dim(Y))) return(NA)
	obs_level <- mean(Y[obs_cells])
	pred_level <- mean(YPS[obs_cells]) / n_stored
	if (!is.finite(obs_level) || obs_level <= 0 || !is.finite(pred_level)) return(NA)
	level_ratio <- pred_level / obs_level
	off <- level_ratio > 1.5 || level_ratio < 1 / 1.5
	if (off) {
		level_word <- if (identical(family, "poisson")) "mean count" else "tie rate"
		cli::cli_abort(
			c(
				"The stored predictive draws have {level_word} {round(pred_level, 3)} against an observed {round(obs_level, 3)}.",
				"i" = "The chain has not found the level of the panel, and no fit is returned."
			),
			class = "lame_predictive_level"
		)
	}
	off
}
# centre one period's latent positions and return the additive parts the centring removes from the product,
# so the linear predictor is unchanged. returns NULL when the positions are already centred
.lame_center_uv_period <- function(Ut, Vt, Gt) {
	ubar <- colMeans(Ut, na.rm = TRUE); vbar <- colMeans(Vt, na.rm = TRUE)
	ubar[!is.finite(ubar)] <- 0; vbar[!is.finite(vbar)] <- 0
	if (all(abs(ubar) < 1e-10) && all(abs(vbar) < 1e-10)) return(NULL)
	Uc <- sweep(Ut, 2L, ubar); Vc <- sweep(Vt, 2L, vbar)
	# a_row and b_col go to the additive effects, const (the grand mean of the product) to the intercept;
	# a_add carries the row part and the constant for callers with no intercept
	a_row <- as.numeric(Uc %*% Gt %*% vbar)
	b_col <- as.numeric(Vc %*% t(Gt) %*% ubar)
	const <- as.numeric(ubar %*% Gt %*% vbar)
	list(
		U = Uc,
		V = Vc,
		a_add = a_row + const,
		b_add = b_col,
		a_row = a_row,
		b_col = b_col,
		const = const
	)
}
.lame_sample_sigma_uv_drift_only <- function(
	U_cube,
	V_cube,
	rho_uv,
	sigma_current,
	symmetric,
	delta_u,
	delta_v = NULL,
	prior_shape = 2,
	prior_scale = 1,
	sigma_max = Inf,
	stat_var_cap = Inf,
	reset_stationary = FALSE
) {
	n_time <- dim(U_cube)[3]
	n_rank <- dim(U_cube)[2]
	if (is.null(n_time) || n_time < 2 || is.null(delta_u)) return(sigma_current)
	ss <- 0
	n_trans <- 0
	# stationary initial condition at t = 1, u_1 ~ N(0, sigma^2 / (1 - rho^2))
	one_m_rho2 <- 1 - min(rho_uv * rho_uv, 0.9801)
	ss <- ss + one_m_rho2 * sum(U_cube[, , 1]^2)
	count_init <- dim(U_cube)[1] * n_rank
	if (!symmetric && !is.null(delta_v)) {
		ss <- ss + one_m_rho2 * sum(V_cube[, , 1]^2)
		count_init <- count_init + dim(V_cube)[1] * dim(V_cube)[2]
	}
	# reset positions under the linked snap scale are stationary draws too
	if (isTRUE(reset_stationary)) {
		st <- .lame_snap_reset_ss(U_cube, V_cube, symmetric, delta_u, delta_v)
		ss <- ss + one_m_rho2 * st$ss
		count_init <- count_init + st$cnt
	}
	for (t in 2:n_time) {
		keep_u <- delta_u[, t] == 0
		if (any(keep_u)) {
			innov_u <- U_cube[keep_u, , t, drop = FALSE] -
				rho_uv * U_cube[keep_u, , t - 1, drop = FALSE]
			ss <- ss + sum(innov_u * innov_u)
			n_trans <- n_trans + sum(keep_u) * dim(U_cube)[2L]
		}
		if (!symmetric && !is.null(delta_v)) {
			keep_v <- delta_v[, t] == 0
			if (any(keep_v)) {
				innov_v <- V_cube[keep_v, , t, drop = FALSE] -
					rho_uv * V_cube[keep_v, , t - 1, drop = FALSE]
				ss <- ss + sum(innov_v * innov_v)
				n_trans <- n_trans + sum(keep_v) * dim(V_cube)[2L]
			}
		}
	}
	count <- n_trans + count_init
	if (!is.finite(ss) || count <= 0) return(sigma_current)
	shape <- count / 2 + prior_shape
	scale <- ss / 2 + prior_scale
	# bound sigma absolutely and through a cap on the stationary variance of a coordinate path
	if (is.finite(stat_var_cap) && stat_var_cap > 0) {
		sigma_max <- min(sigma_max, sqrt(stat_var_cap * one_m_rho2))
	}
	if (is.finite(sigma_max) && sigma_max > 0) {
		# exact truncated draw; projecting onto the cap would put an atom there
		return(.lame_rinvgamma_sd_trunc(shape, scale, sigma_max, sigma_current))
	}
	# X ~ IG(shape, scale) <=> 1/X ~ Gamma(shape, rate = scale)
	sqrt(1 / rgamma(1, shape = shape, rate = scale))
}
.lame_safe_inv_gamma <- function(shape, rate) {
	if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0)
		stop("Inverse-gamma shape and rate must be finite and positive.")
	gdraw <- rgamma(1, shape = shape, rate = rate)
	if (!is.finite(gdraw) || gdraw <= 0)
		stop("Unable to represent the inverse-gamma conditional draw.")
	out <- 1 / gdraw
	if (!is.finite(out) || out <= 0)
		stop("Unable to represent the inverse-gamma conditional draw.")
	out
}
# Gibbs draw of the diagonal L of a symmetric dynamic product U_t L U_t'.
# the coordinate paths have unit stationary variance, so L carries scale and sign; each L_r is normal under its prior
.lame_sym_dyn_L_draw <- function(E, U_cube, L, s2, l_var = 1) {
	R <- length(L)
	N <- dim(U_cube)[3L]
	n <- dim(U_cube)[1L]
	ut <- upper.tri(matrix(0, n, n))
	if (!is.finite(s2) || s2 <= 0) stop("Invalid symmetric factor observation variance.")
	for (r in seq_len(R)) {
		num <- 0
		den <- 0
		for (t in seq_len(N)) {
			Ut <- matrix(U_cube[, , t], n, R)
			others <- if (R > 1L) {
				Ut[, -r, drop = FALSE] %*% diag(L[-r], R - 1L) %*% t(Ut[, -r, drop = FALSE])
			} else 0
			Er <- E[, , t] - others
			uu <- tcrossprod(Ut[, r])
			ok <- ut & is.finite(Er)
			num <- num + sum(Er[ok] * uu[ok])
			den <- den + sum(uu[ok]^2)
		}
		# the weights are N(0, l_var): 1 for a unit response, the response variance for the normal family
		prec <- 1 / l_var + den / s2
		L[r] <- stats::rnorm(1, (num / s2) / prec, sqrt(1 / prec))
	}
	L
}
# Metropolis draw of rho_uv in the symmetric gauge (unit stationary variance, innovation variance 1 - rho^2),
# with the normal prior on rho truncated to (-0.99, 0.99)
.lame_rho_uv_unit_mh <- function(
	U_cube,
	rho_current,
	prior_mean = 0.9,
	prior_sd = 0.1
) {
	N <- dim(U_cube)[3L]
	if (is.null(N) || N < 2L) return(rho_current)
	n <- dim(U_cube)[1L]
	R <- dim(U_cube)[2L]
	s11 <- 0; s12 <- 0; s22 <- 0
	for (t in 2:N) {
		U1 <- matrix(U_cube[, , t - 1L], n, R)
		U2 <- matrix(U_cube[, , t], n, R)
		s11 <- s11 + sum(U1 * U1)
		s12 <- s12 + sum(U1 * U2)
		s22 <- s22 + sum(U2 * U2)
	}
	m <- n * R * (N - 1L)
	if (!is.finite(s11) || s11 <= 0) return(rho_current)
	log_post <- function(r) {
		v <- 1 - r * r
		-0.5 * m * log(v) - (s22 - 2 * r * s12 + r * r * s11) / (2 * v) +
			stats::dnorm(r, prior_mean, prior_sd, log = TRUE)
	}
	lo <- -0.99
	hi <- 0.99
	rc <- min(max(rho_current, lo), hi)
	step <- min(max(2 / sqrt(m), 0.01), 0.2)
	prop <- rc + stats::rnorm(1, 0, step)
	if (prop <= lo || prop >= hi) return(rc)
	la <- log_post(prop) - log_post(rc)
	if (is.finite(la) && log(stats::runif(1)) < la) prop else rc
}
# innovation sd of the symmetric gauge, whose coordinate paths have unit stationary variance
.lame_sym_uv_innov_sd <- function(rho_uv) {
	sqrt(max(1 - rho_uv^2, 1e-6))
}
# additive-effects prior for the directed coefficient kernel: a symmetric fit ties the row and column
# effects (a_i = b_i) through a rank-one covariance
.lame_sab_gls <- function(Sab, symmetric) {
	if (isTRUE(symmetric)) matrix(Sab[1L, 1L], 2L, 2L) else Sab
}
# data-adaptive scale of the additive-effects prior for the ordinal, cbin and frn families, from the
# probit-moment linearisation the binary family uses. returns the pooled scale, a per-side scale and the df multiplier.
# for cbin and frn with a cap, rows at their cap are censored, so the moments are read off the rows below it
.lame_probit_ab_scale <- function(Y, family, odmax = NULL) {
	fallback <- list(vscale = 1, vdfmlt = 1, vscale_row = 1, vscale_col = 1)
	slices <- if (is.list(Y)) Y
	else if (length(dim(Y)) == 3L) lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t])
	else list(Y)
	slices <- Filter(function(m) !is.null(m) && length(m) > 0L, slices)
	if (!length(slices)) return(fallback)
	slices <- lapply(slices, function(m) if (is.matrix(m)) m else as.matrix(m))
	first <- slices[[1L]]
	n_a <- nrow(first); n_b <- ncol(first)
	if (is.null(n_a) || is.null(n_b) || n_a < 1L || n_b < 1L) return(fallback)
	dich <- family %in% c("cbin", "frn")
	# cells of rows at their cap are left out; an actor seen only in such rows is left out of its margin's spread
	drop_a <- rep(FALSE, n_a); drop_b <- rep(FALSE, n_b)
	slices_all <- slices
	below <- NULL
	# a single slice reads every row
	if (dich && !is.null(odmax) && length(slices) >= 2L) {
		below <- .lame_prior_rows_below_cap(slices, odmax)
		if (!is.null(below)) {
			had_a <- numeric(n_a); had_b <- numeric(n_b)
			kept_a <- numeric(n_a); kept_b <- numeric(n_b)
			for (t in seq_along(slices)) {
				m <- slices[[t]]
				if (nrow(m) != n_a || ncol(m) != n_b) next
				ob <- !is.na(m)
				if (n_a == n_b) diag(ob) <- FALSE
				kp <- ob & below[[t]]
				had_a <- had_a + rowSums(ob); had_b <- had_b + colSums(ob)
				kept_a <- kept_a + rowSums(kp); kept_b <- kept_b + colSums(kp)
				m[!below[[t]], ] <- NA
				slices[[t]] <- m
			}
			drop_a <- had_a > 0 & kept_a == 0
			drop_b <- had_b > 0 & kept_b == 0
		}
	}
	yv <- unlist(lapply(slices, as.numeric)); yv <- yv[is.finite(yv)]
	if (length(yv) < 2L) return(fallback)
	if (dich) {
		ybar <- mean(yv > 0)
		if (!(ybar > 0 && ybar < 1)) return(fallback)
		slope <- stats::dnorm(stats::qnorm(ybar))
		score <- function(m) (1 * (m > 0) - ybar) / slope
	} else {
		lev <- sort(unique(yv))
		if (length(lev) < 2L) return(fallback)
		idx <- match(yv, lev)
		shares <- cumsum(tabulate(idx, length(lev)) / length(idx))[-length(lev)]
		slope <- sum(stats::dnorm(stats::qnorm(shares)))
		mean_idx <- mean(idx)
		score <- function(m) (match(m, lev) - mean_idx) / slope
	}
	if (!is.finite(slope) || slope <= 0) return(fallback)
	num_a <- numeric(n_a); cnt_a <- numeric(n_a)
	num_b <- numeric(n_b); cnt_b <- numeric(n_b)
	for (m in slices) {
		if (nrow(m) != n_a || ncol(m) != n_b) next
		E <- matrix(score(m), n_a, n_b)
		E[!is.finite(m)] <- NA_real_
		if (n_a == n_b) diag(E) <- NA_real_
		ok <- is.finite(E)
		E[!ok] <- 0
		num_a <- num_a + rowSums(E); cnt_a <- cnt_a + rowSums(ok)
		num_b <- num_b + colSums(E); cnt_b <- cnt_b + colSums(ok)
	}
	a_tmp <- ifelse(cnt_a > 0, num_a / pmax(cnt_a, 1), 0)[!drop_a]
	b_tmp <- ifelse(cnt_b > 0, num_b / pmax(cnt_b, 1), 0)[!drop_b]
	vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
	if (!is.finite(vscale) || vscale <= 0) vscale <- 1
	vdfmlt <- 1
	if (dich) {
		PHAT <- stats::pnorm(stats::qnorm(ybar) + outer(a_tmp, b_tmp, "+"))
		vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
		if (!is.finite(vdfmlt) || vdfmlt < 1) vdfmlt <- 1
	}
	out <- list(vscale = vscale, vdfmlt = vdfmlt, vscale_row = vscale, vscale_col = vscale)
	if (!dich) return(out)
	# single-period panels, and panels whose rows at the cap reach .lame_odmax_sat_threshold
	# (or cannot be ruled out to), keep the pooled scale: too few rows sit below the cap to read a sender spread
	if (length(slices_all) < 2L) return(out)
	if (!is.null(odmax)) {
		sat <- .lame_odmax_saturated_share(slices_all, odmax, family)
		if (is.finite(sat) && (sat >= .lame_odmax_sat_threshold || isTRUE(attr(sat, "not_ruled_out")))) {
			return(out)
		}
	}
	# cbin and frn take one scale per side; the column side is the spread of the column means over every row
	vrow <- stats::var(a_tmp)
	vcol <- stats::var(b_tmp)
	if (!is.null(below)) {
		vcol <- .lame_probit_col_scale(slices_all)
		# the rows below the cap are a truncated sample, so the rows at the cap enter as censored observations
		vr <- .lame_censored_row_scale(
			slices_all,
			odmax,
			center = if (is.finite(vcol) && vcol > 0) vcol / (1 + vcol) else NULL
		)
		if (is.finite(vr) && vr > 0 && is.finite(vcol) && vcol > 0) {
			# the censored fit's probit link absorbs the receivers' spread, hence the factor 1 + vcol; the result is held
			# at or above the spread of the rows below the cap and within a factor of four of the column side
			vrow <- min(max(vr * (1 + vcol), vrow, vcol / 4), 4 * vcol)
		} else if (is.finite(vcol) && vcol > 0) {
			vrow <- min(max(vrow, vcol / 4), 4 * vcol)
		}
	}
	if (is.finite(vrow) && vrow > 0) out$vscale_row <- vrow
	if (is.finite(vcol) && vcol > 0) out$vscale_col <- vcol
	out
}
# column side of the cbin / frn prior scale: spread of the column means of the probit scores over every row
.lame_probit_col_scale <- function(slices) {
	yv <- unlist(lapply(slices, as.numeric)); yv <- yv[is.finite(yv)]
	if (length(yv) < 2L) return(NA_real_)
	ybar <- mean(yv > 0)
	if (!(ybar > 0 && ybar < 1)) return(NA_real_)
	slope <- stats::dnorm(stats::qnorm(ybar))
	n_b <- ncol(slices[[1L]])
	num_b <- numeric(n_b); cnt_b <- numeric(n_b)
	for (m in slices) {
		if (ncol(m) != n_b) next
		E <- (1 * (m > 0) - ybar) / slope
		E[!is.finite(m)] <- NA_real_
		if (nrow(m) == ncol(m)) diag(E) <- NA_real_
		ok <- is.finite(E)
		E[!ok] <- 0
		num_b <- num_b + colSums(E); cnt_b <- cnt_b + colSums(ok)
	}
	if (sum(cnt_b > 0) < 2L) return(NA_real_)
	stats::var((num_b / pmax(cnt_b, 1))[cnt_b > 0])
}
# Gauss-Hermite nodes and weights against a standard normal (Golub-Welsch)
.lame_gh_normal <- function(k = 20L) {
	J <- matrix(0, k, k)
	for (i in seq_len(k - 1L)) J[i, i + 1L] <- J[i + 1L, i] <- sqrt(i)
	e <- eigen(J, symmetric = TRUE)
	list(x = e$values, w = e$vectors[1L, ]^2)
}
# maximum-likelihood variance of the sender effects on the probit scale from the row outdegrees,
# with rows at their cap as censored observations and the random effect integrated by Gauss-Hermite quadrature.
# center, when given, is the column side the fit is drawn toward on panels with hidden cells
.lame_censored_row_scale <- function(slices, odmax, center = NULL) {
	cap0 <- suppressWarnings(as.numeric(odmax))
	if (!length(cap0) || !any(is.finite(cap0))) return(NA_real_)
	od <- mm <- np <- ff <- cp <- id <- numeric(0)
	cen <- logical(0)
	for (m in slices) {
		cap <- if (length(cap0) == 1L) rep(cap0, nrow(m)) else cap0
		if (length(cap) != nrow(m)) next
		sq <- nrow(m) == ncol(m) && nrow(m) > 1L
		ob <- !is.na(m)
		if (sq) diag(ob) <- FALSE
		n_obs <- rowSums(ob)
		present <- colSums(ob) > 0L
		n_poss <- rep(sum(present), nrow(m))
		if (sq) n_poss <- n_poss - present
		o <- rowSums(m > 0 & ob, na.rm = TRUE)
		f <- pmin(n_obs / pmax(n_poss, 1), 1)
		capr <- pmin(round(cap), pmax(n_poss, 1))
		keep <- n_obs > 0L & is.finite(capr)
		od <- c(od, o[keep]); mm <- c(mm, n_obs[keep]); np <- c(np, n_poss[keep])
		ff <- c(ff, f[keep]); cp <- c(cp, capr[keep])
		cen <- c(cen, !(o[keep] < capr[keep] * f[keep]))
		id <- c(id, which(keep))
	}
	if (length(od) < 10L || !any(cen) || all(cen)) return(NA_real_)
	gh <- .lame_gh_normal(20L)
	ids <- factor(id)
	thin <- ff < 1
	nll <- function(par) {
		s <- exp(par[2L])
		P <- stats::pnorm(outer(rep(par[1L], length(od)), s * gh$x, "+"))
		# one row per row-period, one column per node
		L <- matrix(stats::dbinom(od, mm, P, log = TRUE), length(od), length(gh$x))
		Lc <- matrix(
			stats::pbinom(od - 1, mm, P, lower.tail = FALSE, log.p = TRUE),
			length(od),
			length(gh$x)
		)
		L[cen, ] <- Lc[cen, ]
		if (any(thin)) {
			# rows with hidden cells: nominations thinned at the row's observed share
			Pt <- P[thin, , drop = FALSE]
			nt <- np[thin]; ct <- cp[thin]; ot <- od[thin]; ft <- ff[thin]
			acc <- matrix(0, nrow(Pt), ncol(Pt))
			for (k in 0:max(ct)) {
				pk <- matrix(0, nrow(Pt), ncol(Pt))
				lt <- k < ct; eq <- k == ct
				if (any(lt)) pk[lt, ] <- matrix(stats::dbinom(k, nt, Pt), nrow(Pt))[lt, ]
				if (any(eq)) pk[eq, ] <- matrix(stats::pbinom(k - 1, nt, Pt, lower.tail = FALSE), nrow(Pt))[eq, ]
				acc <- acc + pk * stats::dbinom(ot, k, ft)
			}
			L[thin, ] <- log(pmax(acc, 1e-300))
		}
		# per actor: sum over its row-periods, then integrate over the node
		La <- rowsum(L, ids)
		mx <- apply(La, 1L, max)
		ll <- mx + log(as.numeric(exp(La - mx) %*% gh$w))
		v <- -sum(ll)
		# with hidden cells the likelihood is flat in the spread, so a normal prior on the log variance
		# draws the fit toward the column side
		if (any(thin) && !is.null(center)) v <- v + 0.5 * (2 * par[2L] - log(center))^2
		if (is.finite(v)) v else 1e100
	}
	p0 <- min(max(mean(od / mm), 1e-3), 1 - 1e-3)
	fit <- tryCatch(
		stats::optim(c(stats::qnorm(p0), log(0.5)), nll, method = "BFGS"),
		error = function(e) NULL
	)
	if (is.null(fit) || fit$convergence != 0L) return(NA_real_)
	s2 <- exp(2 * fit$par[2L])
	if (!is.finite(s2) || s2 > 25) NA_real_ else s2
}
# which rows of each slice sit below their nomination cap (TRUE = use the row), with hidden cells
# thinning the cap at the row's observed share. returns NULL when no row sits at its cap or fewer
# than five row actors sit below it, in which case every row is used
.lame_prior_rows_below_cap <- function(slices, odmax) {
	cap0 <- suppressWarnings(as.numeric(odmax))
	if (!length(cap0) || !any(is.finite(cap0))) return(NULL)
	n_a <- nrow(slices[[1L]])
	any_capped <- FALSE
	below_actor <- rep(FALSE, n_a)
	keep <- lapply(slices, function(m) {
		k <- rep(TRUE, nrow(m))
		cap <- if (length(cap0) == 1L) rep(cap0, nrow(m)) else cap0
		if (length(cap) != nrow(m)) return(k)
		sq <- nrow(m) == ncol(m) && nrow(m) > 1L
		ob <- !is.na(m)
		if (sq) diag(ob) <- FALSE
		n_obs <- rowSums(ob)
		# the cells a row could have nominated: columns observed anywhere in the slice, less its own
		present <- colSums(ob) > 0L
		n_poss <- rep(sum(present), nrow(m))
		if (sq) n_poss <- n_poss - present
		od <- rowSums(m > 0 & ob, na.rm = TRUE)
		f <- pmin(n_obs / pmax(n_poss, 1), 1)
		capr <- pmin(round(cap), pmax(n_poss, 1))
		at_cap <- n_obs > 0L & is.finite(capr) & capr >= 1 & !(od < capr * f)
		k[at_cap] <- FALSE
		k
	})
	for (t in seq_along(slices)) {
		m <- slices[[t]]
		obs_row <- rowSums(!is.na(m)) > 0L
		any_capped <- any_capped || any(!keep[[t]])
		if (length(keep[[t]]) == n_a) below_actor <- below_actor | (keep[[t]] & obs_row)
	}
	if (!any_capped || sum(below_actor) < 5L) return(NULL)
	keep
}
# share of the observed rows at their cap above which the row side of the prior scale is pooled
# with the column side (a sampler setting, not a reporting rule)
.lame_odmax_sat_threshold <- 0.8
# the same setting for a single period
.lame_odmax_sat_threshold_single <- 0.6
.lame_odmax_threshold_for <- function(n_periods) {
	if (isTRUE(n_periods <= 1L)) .lame_odmax_sat_threshold_single else .lame_odmax_sat_threshold
}
# share of observed rows whose observed outdegree sits at their nomination cap, the quantity
# behind the cap rules in ?lame. NA for the families without a cap; a missing cap resolves as in lame()
.lame_odmax_saturated_share <- function(Y, odmax, family) {
	if (length(family) != 1L || !is.character(family) ||
			!family %in% c("cbin", "frn")) return(NA_real_)
	slices <- if (is.list(Y)) Y
	else if (length(dim(Y)) == 3L) lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t])
	else list(Y)
	slices <- lapply(slices, function(m) {
		m <- suppressWarnings(as.matrix(m))
		if (!is.numeric(m) && !is.logical(m)) return(NULL)
		# the diagonal is never a nomination
		if (nrow(m) == ncol(m) && nrow(m) > 1L) diag(m) <- NA
		m
	})
	slices <- slices[!vapply(slices, is.null, logical(1))]
	if (!length(slices)) return(NA_real_)
	cap0 <- suppressWarnings(as.numeric(odmax))
	if (!length(cap0) || !any(is.finite(cap0))) {
		cap0 <- suppressWarnings(max(
			vapply(
				slices,
				function(m)
					max(rowSums(m > 0, na.rm = TRUE)),
				numeric(1)
			),
			na.rm = TRUE
		))
	}
	if (!any(is.finite(cap0))) return(NA_real_)
	od_all <- f_all <- cap_all <- numeric(0)
	for (m in slices) {
		# a per-row cap applies only to a slice it matches
		cap <- if (length(cap0) == 1L) rep(cap0, nrow(m)) else cap0
		if (length(cap) != nrow(m)) next
		# a row with no observed cell, or no meaningful cap, says nothing
		n_obs <- rowSums(!is.na(m))
		keep <- (n_obs > 0L) & is.finite(cap) & cap >= 1
		if (!any(keep)) next
		# the cells a row could have nominated: columns observed anywhere in the slice, less its own
		present <- colSums(!is.na(m)) > 0L
		n_poss <- rep(sum(present), nrow(m))
		if (nrow(m) == ncol(m)) n_poss <- n_poss - present
		od_all <- c(od_all, rowSums(m > 0, na.rm = TRUE)[keep])
		f_all <- c(f_all, pmin(n_obs[keep] / pmax(n_poss[keep], 1), 1))
		cap_all <- c(cap_all, pmin(round(cap[keep]), pmax(n_poss[keep], 1)))
	}
	if (!length(od_all)) return(NA_real_)
	complete <- f_all >= 1
	if (all(complete)) return(mean(od_all >= cap_all))
	# with hidden cells the share at the cap is estimated: EM for the distribution of the gap between
	# the cap and the full outdegree, the share being the mass at a gap of zero
	gmax <- max(cap_all)
	gaps <- 0:gmax
	lik <- matrix(0, length(od_all), length(gaps))
	for (g in gaps) {
		ok <- cap_all - g >= 0
		lik[ok, g + 1L] <- stats::dbinom(od_all[ok], cap_all[ok] - g, f_all[ok])
	}
	# a row showing more nominations than its cap allows sits at its cap
	over <- rowSums(lik) <= 0
	if (any(over)) lik[over, 1L] <- 1
	# no row could sit at its cap
	if (max(lik[, 1L]) < 1e-12) return(0)
	# gaps no row gives any likelihood to carry no mass
	g_keep <- sort(union(1L, which(apply(lik, 2L, max) >= 1e-12)))
	lik <- lik[, g_keep, drop = FALSE]
	gaps <- gaps[g_keep]
	# EM for the gap distribution, optionally with the mass at zero held fixed
	em_gap <- function(p0_fixed = NULL) {
		p_gap <- rep(1 / length(gaps), length(gaps))
		if (!is.null(p0_fixed)) {
			p_gap <- c(p0_fixed, rep((1 - p0_fixed) / max(length(gaps) - 1L, 1), length(gaps) - 1L))
		}
		for (it in seq_len(1000L)) {
			w <- sweep(lik, 2L, p_gap, "*")
			w <- w / pmax(rowSums(w), 1e-300)
			p_new <- colMeans(w)
			if (!is.null(p0_fixed) && length(gaps) > 1L) {
				rest <- p_new[-1L]
				p_new <- c(p0_fixed, (1 - p0_fixed) * rest / max(sum(rest), 1e-300))
			}
			done <- max(abs(p_new - p_gap)) < 1e-9
			p_gap <- p_new
			if (done) break
		}
		list(p = p_gap, ll = sum(log(pmax(lik %*% p_gap, 1e-300))))
	}
	fit_free <- em_gap()
	share <- fit_free$p[1L]
	thr <- .lame_odmax_threshold_for(length(slices))
	# with hidden cells, a share at the pooling threshold that a likelihood ratio test (5 percent)
	# does not reject marks the panel saturated; the test runs only within 0.08 of the threshold
	if (share < thr && share >= thr - 0.08 &&
			gmax >= 1L) {
		fit_thr <- em_gap(thr)
		if (2 * (fit_free$ll - fit_thr$ll) < stats::qchisq(0.95, 1)) {
			attr(share, "not_ruled_out") <- TRUE
		}
	}
	share
}
# events per actor-period below which a directed binary panel does not report va and vb,
# and its intercept, without and with a multiplicative term (see ?lame)
.lame_sparse_period_rate <- 1.2
.lame_sparse_period_rate_intercept <- 0.6
.lame_sparse_period_rate_intercept_rank <- 1.6
.lame_ab_sparse <- function(per_period) {
	pp <- suppressWarnings(as.numeric(per_period))
	length(pp) == 1L && is.finite(pp) && pp < .lame_sparse_period_rate
}
# the rate is higher with a multiplicative term, whose mean the intercept trades off against
.lame_sparse_intercept_rate <- function(fit) {
	R <- suppressWarnings(as.numeric(fit[["R"]] %||% fit[["rank"]]))
	if (length(R) == 1L && is.finite(R) && R > 0) .lame_sparse_period_rate_intercept_rank else .lame_sparse_period_rate_intercept
}
# whether the user supplied the additive-effects prior (Sab0, eta0 or etaab)
.lame_ab_prior_given <- function(prior) {
	if (!is.list(prior) || !length(prior)) return(FALSE)
	any(!vapply(prior[c("Sab0", "eta0", "etaab")], is.null, logical(1)))
}
# number of periods a fit or its summary carries, NA when it does not say
.lame_fit_n_periods <- function(x) {
	n <- x[["n.periods"]] %||% x[["n_time"]]
	if (is.numeric(n) && length(n) == 1L && is.finite(n)) return(n)
	if (inherits(x, c("lame", "summary.lame"))) {
		return(if (is.list(x[["YPM"]])) length(x[["YPM"]]) else NA_real_)
	}
	if (inherits(x, c("ame", "summary.ame"))) return(1)
	NA_real_
}
# shares of rows at the cap from which a one-period cbin or frn design leaves its row level
# unreported: every row, or 0.5 with hidden cells (see odmax in ?lame)
.lame_cap_level_share <- 0.99
.lame_cap_level_share_hidden <- 0.5
.lame_cap_level_periods <- 1L
# whether a cbin or frn fit leaves the intercept and row-covariate coefficients unreported
.lame_cap_level_unreported <- function(fit) {
	if (!is.list(fit)) return(FALSE)
	fam <- fit[["family"]] %||% ""
	if (length(fam) != 1L || !is.character(fam) || !fam %in% c("cbin", "frn")) return(FALSE)
	sat <- suppressWarnings(as.numeric(fit[["odmax_saturated_share"]]))
	n <- .lame_fit_n_periods(fit)
	if (length(sat) != 1L || !is.finite(sat) || !is.finite(n) || n > .lame_cap_level_periods) return(FALSE)
	sat >= .lame_cap_level_share ||
		(isTRUE(fit[["cells_hidden"]]) && sat >= .lame_cap_level_share_hidden)
}
# a cbin or frn panel of two periods or more with at least .lame_cap_shrunk_share of rows at the cap
# and posterior mean of va in .lame_cap_shrunk_va withholds its row level: va stays near its prior there
.lame_cap_shrunk_share <- 0.5
.lame_cap_shrunk_va <- c(0.5, 1.5)
.lame_cap_shrunk_unreported <- function(fit) {
	if (!is.list(fit)) return(FALSE)
	fam <- fit[["family"]] %||% ""
	if (length(fam) != 1L || !is.character(fam) || !fam %in% c("cbin", "frn")) return(FALSE)
	if (!identical(fit[["mode"]] %||% "unipartite", "unipartite")) return(FALSE)
	sat <- suppressWarnings(as.numeric(fit[["odmax_saturated_share"]]))
	n <- .lame_fit_n_periods(fit)
	if (length(sat) != 1L || !is.finite(sat) || !is.finite(n) || n <= .lame_cap_level_periods) return(FALSE)
	VC <- fit[["VC"]]
	if (!is.matrix(VC) || !"va" %in% colnames(VC)) return(FALSE)
	va <- mean(VC[, "va"], na.rm = TRUE)
	is.finite(va) && sat >= .lame_cap_shrunk_share && va >= .lame_cap_shrunk_va[1L] && va < .lame_cap_shrunk_va[2L]
}
# whether a sparse directed binary panel withholds its effect variances or intercept (see ?lame)
.lame_ab_sparse_unreported <- function(fit, variances = TRUE) {
	if (!is.list(fit)) return(FALSE)
	if (!identical(fit[["family"]], "binary") ||
			!identical(fit[["mode"]] %||% "unipartite", "unipartite")) return(FALSE)
	ev <- suppressWarnings(as.numeric(fit[["events_per_actor_period"]]))
	if (length(ev) != 1L || !is.finite(ev)) return(FALSE)
	ev < (if (variances) .lame_sparse_period_rate else .lame_sparse_intercept_rate(fit))
}
# why va and vb are not reported, or NA when they are
.lame_ab_data_reason <- function(fit) {
	if (.lame_ab_sparse_unreported(fit, variances = TRUE)) "sparse" else NA_character_
}
# why the intercept and row-covariate coefficients are not reported, or NA
.lame_ab_coef_reason <- function(fit) {
	if (.lame_cap_level_unreported(fit)) return("saturated")
	if (.lame_cap_shrunk_unreported(fit)) return("cap_shrunk")
	if (.lame_ab_sparse_unreported(fit, variances = FALSE)) return("sparse")
	NA_character_
}
# why the additive-effect variances are not reported, or NA: the data's reason, unless the fit
# ran under a supplied effects prior, whose posterior is reported with that prior named
.lame_ab_var_unident_reason <- function(fit) {
	reason <- .lame_ab_data_reason(fit)
	if (!is.na(reason) && isTRUE(fit[["prior_supplied"]])) return(NA_character_)
	reason
}
# whether a fit carries too little information to report the additive-effect variances
.lame_ab_var_unidentified <- function(fit) {
	!is.na(.lame_ab_var_unident_reason(fit))
}
# how a nomination-cap share reads in the messages
.lame_sat_phrase <- function(share, hidden = FALSE) {
	pct <- round(100 * as.numeric(share)[1L])
	if (isTRUE(hidden)) {
		paste0("an estimated ", pct, " percent of the observed rows nominate the full odmax with cells hidden")
	} else {
		paste0(pct, " percent of the observed rows nominate the full odmax")
	}
}
# variance rows a reason hides: va and vb of a sparse panel
.lame_ab_var_hidden_rows <- function(reason) {
	if (identical(reason, "sparse")) return(c("va", "vb"))
	character(0)
}
# coefficient rows a reason hides: the intercept of a sparse panel; the intercept and
# row-covariate coefficients when the cap leaves the level of the rows unfixed
.lame_ab_coef_hidden_rows <- function(reason, coef_names) {
	if (is.na(reason) || !length(coef_names)) return(character(0))
	if (reason %in% c("saturated", "cap_shrunk")) return(
		coef_names[coef_names == "intercept" | grepl("_row$", coef_names)]
	)
	if (identical(reason, "sparse")) return(coef_names[coef_names == "intercept"])
	character(0)
}
# the effects prior as the reporting methods name it: diag(Sab0) and the degrees of freedom
.lame_ab_prior_text <- function(fit_or_summary) {
	txt <- fit_or_summary[["ab_prior_text"]]
	if (is.character(txt) && length(txt) == 1L) return(txt)
	pr <- fit_or_summary[["prior"]]
	if (!is.list(pr)) return("")
	fmt <- function(v) format(signif(as.numeric(v), 3), trim = TRUE)
	parts <- character(0)
	S0 <- pr[["Sab0"]]
	if (is.matrix(S0)) {
		parts <- c(parts, paste0("Sab0 diag = ", paste(fmt(diag(S0)), collapse = ", ")))
	} else if (is.numeric(S0) && length(S0)) {
		parts <- c(parts, paste0("Sab0 = ", paste(fmt(S0), collapse = ", ")))
	}
	eta0 <- pr[["eta0"]]
	if (is.numeric(eta0) && length(eta0)) parts <- c(parts, paste0("eta0 = ", fmt(eta0[1L])))
	etaab <- pr[["etaab"]]
	if (is.numeric(etaab) && length(etaab) &&
			!(is.numeric(eta0) && length(eta0) && isTRUE(all.equal(etaab[1L], eta0[1L])))) {
		parts <- c(parts, paste0("etaab = ", fmt(etaab[1L])))
	}
	paste(parts, collapse = ", ")
}
# what a fit or its summary does not report, as plain sentences (none when it reports everything)
.lame_ab_statements <- function(x) {
	pick <- function(field, fun) if (is.character(x[[field]]) && length(x[[field]]) == 1L) x[[field]] else fun(x)
	coef_reason <- pick("ab_coef_reason", .lame_ab_coef_reason)
	data_reason <- pick("ab_data_reason", .lame_ab_data_reason)
	out <- character(0)
	ev <- round(suppressWarnings(as.numeric(x[["events_per_actor_period"]]))[1L], 2)
	if (identical(data_reason, "sparse")) {
		thr <- .lame_sparse_period_rate
		if (isTRUE(x[["prior_supplied"]])) {
			out <- c(
				out,
				paste0(
					"va and vb are reported under the supplied prior (",
					.lame_ab_prior_text(x),
					"): the panel carries ",
					ev,
					" events per actor-period, below ",
					thr,
					", so the data alone do not estimate them."
				)
			)
		} else {
			out <- c(
				out,
				paste0(
					"va and vb are not estimated: the panel carries ",
					ev,
					" events per actor-period, below ",
					thr,
					", too few for their intervals to cover the generating variances. ",
					"The draws are in fit$VC; supplying prior$Sab0 or prior$eta0 reports them under that prior."
				)
			)
		}
	}
	if (identical(coef_reason, "cap_shrunk")) {
		share <- round(100 * as.numeric(x[["odmax_saturated_share"]])[1L])
		out <- c(
			out,
			paste0(
				"The intercept and any row-covariate coefficient are not estimated: ",
				share,
				" percent of the observed rows nominate the full odmax and the posterior mean of va lies between ",
				.lame_cap_shrunk_va[1L],
				" and ",
				.lame_cap_shrunk_va[2L],
				", where the cap hides the level of the rows. The draws are in fit$BETA. va is estimated."
			)
		)
	}
	if (identical(coef_reason, "sparse")) {
		thr <- .lame_sparse_intercept_rate(x)
		out <- c(
			out,
			paste0(
				"The intercept is not estimated: the panel carries ",
				ev,
				" events per actor-period, below ",
				thr,
				if (thr > .lame_sparse_period_rate_intercept) " for a fit with a multiplicative term" else "",
				", too few for its interval to cover the generating value. The draws are in fit$BETA."
			)
		)
	}
	if (identical(coef_reason, "saturated")) {
		share <- x[["odmax_saturated_share"]]
		hidden <- isTRUE(x[["cells_hidden"]]) && as.numeric(share)[1L] < .lame_cap_level_share
		out <- c(
			out,
			paste0(
				"The intercept and any row-covariate coefficient are not estimated: ",
				.lame_sat_phrase(share, hidden),
				" in a one-period design, where a row at its cap shows only that its level fills the cap and nothing bounds the level of a row above. ",
				"va, vb and cab are estimated. The draws are in fit$BETA."
			)
		)
	}
	out
}

# whether the snap and t position hyperparameters are withheld: at rank two or more their chains
# mix poorly (see ?lame); a fit that does not record its rank counts as rank two or more
.lame_uv_rank_withheld <- function(fit) {
	R <- fit[["R"]]
	if (is.null(R) && length(dim(fit[["U"]])) >= 2L) R <- dim(fit[["U"]])[2L]
	!(is.numeric(R) && length(R) == 1L && is.finite(R) && R < 2)
}
.lame_snap_rate_unidentified <- function(fit) {
	if (!is.list(fit)) return(FALSE)
	if (!isTRUE(fit[["dynamic_uv"]])) return(FALSE)
	kind <- fit[["dynamic_uv_kind"]] %||% ((fit[["dynamic"]] %||% list())[["uv_kind"]] %||% "ar1")
	identical(kind, "snap") && .lame_uv_rank_withheld(fit)
}
.lame_uv_hyper_unidentified <- function(fit) {
	if (!is.list(fit)) return(FALSE)
	if (!isTRUE(fit[["dynamic_uv"]])) return(FALSE)
	kind <- fit[["dynamic_uv_kind"]] %||% ((fit[["dynamic"]] %||% list())[["uv_kind"]] %||% "ar1")
	if (!.lame_uv_rank_withheld(fit)) return(FALSE)
	identical(kind, "snap") || (identical(kind, "t") && !isTRUE(fit[["als_start"]]))
}
.lame_expand_rho <- function(rho, n_time) {
	if (!is.numeric(rho) || !length(rho) || any(!is.finite(rho)) ||
			any(abs(rho) >= 1)) {
		cli::cli_abort("Dyadic correlations must be finite and strictly between -1 and 1.")
	}
	if (length(n_time) != 1L || !is.finite(n_time) || n_time < 1 ||
			n_time != as.integer(n_time) || !length(rho) %in% c(1L, n_time)) {
		cli::cli_abort("Dyadic correlations must have length one or match the number of periods.")
	}
	if (length(rho) == n_time) as.numeric(rho) else rep(as.numeric(rho), n_time)
}
.lame_rho_at <- function(rho, t) {
	if (length(rho) <= 1L) as.numeric(rho)[1L] else as.numeric(rho)[t]
}
.lame_rho_scalar <- function(rho) {
	mean(.lame_expand_rho(rho, length(rho)))
}
.lame_sample_s2_dynamic_rho <- function(E, rho, current = NULL, s20 = 1, beta_df = 0, beta_ss = 0) {
	s20 <- .lame_validate_s20(s20)
	if (!is.numeric(E) || length(dim(E)) != 3L || nrow(E) != ncol(E)) {
		cli::cli_abort("Residuals must be a square numeric array with one slice per period.")
	}
	n_time <- dim(E)[3L]
	rho <- .lame_expand_rho(rho, n_time)
	if (length(beta_df) != 1L || !is.finite(beta_df) || beta_df < 0 ||
			length(beta_ss) != 1L || !is.finite(beta_ss) || beta_ss < 0) {
		cli::cli_abort("Coefficient prior contributions must be finite nonnegative scalars.")
	}
	rate <- s20 / 2 + beta_ss / 2
	count <- 0
	for (t in seq_len(n_time)) {
		rt <- E[, , t, drop = FALSE][, , 1L]
		ut <- upper.tri(rt)
		e1 <- rt[ut]
		e2 <- t(rt)[ut]
		if (any(is.infinite(e1)) || any(is.infinite(e2))) {
			cli::cli_abort("Residuals contain non-finite observed values.")
		}
		ok <- is.finite(e1) & is.finite(e2)
		# orthogonal dyad coordinates avoid cancellation near either boundary
		plus <- (e1[ok] / 2 + e2[ok] / 2) / sqrt(1 + rho[t])
		minus <- (e1[ok] / 2 - e2[ok] / 2) / sqrt(1 - rho[t])
		rate <- rate + sum(plus^2 + minus^2)
		count <- count + 2L * sum(ok)
	}
	if (!is.finite(rate) || rate <= 0) {
		cli::cli_abort("Residual variance conditional has a non-finite or nonpositive rate.")
	}
	out <- unname(rate / rgamma(1, shape = (count + 1 + beta_df) / 2, rate = 1))
	if (!is.finite(out) || out <= 0) {
		cli::cli_abort("Residual variance draw is not finite and positive.")
	}
	out
}
.lame_center_ab_state <- function(a, b, beta, intercept = TRUE, symmetric = FALSE) {
	# without an intercept the level of the predictor lives in a and b
	if (!isTRUE(intercept) || isTRUE(symmetric) || length(a) == 0L || length(b) == 0L) {
		return(list(a = a, b = b, beta = beta))
	}
	af <- is.finite(a)
	bf <- is.finite(b)
	if (!any(af) || !any(bf)) {
		return(list(a = a, b = b, beta = beta))
	}
	ma <- mean(a[af])
	mb <- mean(b[bf])
	if (!is.finite(ma)) ma <- 0
	if (!is.finite(mb)) mb <- 0
	a[af] <- a[af] - ma
	b[bf] <- b[bf] - mb
	if (isTRUE(intercept) && length(beta) > 0L && is.finite(beta[1L])) {
		beta[1L] <- beta[1L] + ma + mb
	}
	list(a = a, b = b, beta = beta)
}
.lame_uv_ok <- function(U, V, max_abs = Inf) {
	if (is.null(U) || is.null(V)) return(FALSE)
	if (!all(is.finite(U)) || !all(is.finite(V))) return(FALSE)
	if (is.finite(max_abs) && max(abs(c(U, V)), na.rm = TRUE) > max_abs) {
		return(FALSE)
	}
	TRUE
}
.lame_shrink_sparse_ab <- function(
	a,
	b,
	row_obs,
	col_obs,
	min_obs = 0,
	symmetric = FALSE
) {
	if (isTRUE(symmetric) || !is.finite(min_obs) || min_obs <= 1) {
		return(list(a = a, b = b))
	}
	if (length(a) == length(row_obs)) {
		wr <- pmin(1, pmax(0, row_obs / min_obs))
		a <- a * wr
	}
	if (length(b) == length(col_obs)) {
		wc <- pmin(1, pmax(0, col_obs / min_obs))
		b <- b * wc
	}
	list(a = a, b = b)
}
# joint scale move for the binary sampler: the predictor is multiplied by c (beta and the additive
# effects by c, the positions by sqrt(c)) with z integrated out, a long step along the weakly
# identified ridge between level and spread. the dyadic correlation is kept through a shared dyad
# factor, which makes the sweep partially collapsed and exact.
#
# shared dyad factors of the equicorrelated pair error, as the per-cell mean offset sqrt(|rho|) * s_ij * f_d
.lame_dyad_factor_draw <- function(Z, EZ, rho, upper_idx) {
	out <- array(0, dim = dim(Z))
	if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
		stop("Binary auxiliary correlation must lie strictly between -1 and 1.", call. = FALSE)
	ar <- abs(rho)
	if (ar <= 0) return(out)
	n <- dim(Z)[1L]
	sgn <- if (rho < 0) -1 else 1
	k <- sqrt(ar)
	sd_f <- sqrt((1 - ar) / (1 + ar))
	for (t in seq_len(dim(Z)[3L])) {
		E <- Z[, , t] - EZ[, , t]
		E[!is.finite(E)] <- 0
		u <- E + sgn * t(E)
		fm <- matrix(0, n, n)
		m <- k * u[upper_idx] / (1 + ar)
		fm[upper_idx] <- k * (m + sd_f * rnorm(length(m)))
		out[, , t] <- fm + sgn * t(fm)
	}
	out
}
# Gaussian sum of squares and term count of the additive deviations under their AR(1) prior,
# matching .lame_sample_sigma_ab()
.lame_ab_dev_ss <- function(a_dev, b_dev, rho_ab, symmetric = FALSE) {
	N <- ncol(a_dev)
	if (is.null(N) || N < 1L) return(list(ss = 0, cnt = 0L))
	one_m_rho2 <- 1 - min(rho_ab * rho_ab, 0.9801)
	paths <- if (isTRUE(symmetric)) list(a_dev) else list(a_dev, b_dev)
	ss <- 0; cnt <- 0L
	for (M in paths) {
		x1 <- M[, 1L]; x1 <- x1[is.finite(x1)]
		ss <- ss + one_m_rho2 * sum(x1 * x1); cnt <- cnt + length(x1)
		if (N >= 2L) {
			innov <- M[, -1L, drop = FALSE] - rho_ab * M[, -N, drop = FALSE]
			innov <- innov[is.finite(innov)]
			ss <- ss + sum(innov * innov); cnt <- cnt + length(innov)
		}
	}
	list(ss = ss, cnt = cnt)
}
# same for the dynamic latent positions, matching .lame_sample_sigma_uv_drift_only()
.lame_uv_ss <- function(U_cube, V_cube, rho_uv, symmetric = FALSE) {
	N <- dim(U_cube)[3L]
	if (is.null(N) || N < 1L) return(list(ss = 0, cnt = 0L))
	one_m_rho2 <- 1 - min(rho_uv * rho_uv, 0.9801)
	ss <- one_m_rho2 * sum(U_cube[, , 1L]^2)
	cnt <- prod(dim(U_cube)[1:2])
	use_v <- !isTRUE(symmetric) && !is.null(V_cube)
	if (use_v) {
		ss <- ss + one_m_rho2 * sum(V_cube[, , 1L]^2)
		cnt <- cnt + prod(dim(V_cube)[1:2])
	}
	if (N >= 2L) {
		iu <- U_cube[, , -1L, drop = FALSE] - rho_uv * U_cube[, , -N, drop = FALSE]
		ss <- ss + sum(iu * iu); cnt <- cnt + length(iu)
		if (use_v) {
			iv <- V_cube[, , -1L, drop = FALSE] - rho_uv * V_cube[, , -N, drop = FALSE]
			ss <- ss + sum(iv * iv); cnt <- cnt + length(iv)
		}
	}
	list(ss = ss, cnt = cnt)
}
# log prior ratio, Jacobian included, of the deviation innovation variance under the scale move;
# v_max is the bound prior$ab_var_cap puts on the stationary variance
.lame_scale_ig_logratio <- function(v, v_new, shape, scale, v_max = Inf) {
	if (!is.finite(v) || v <= 0 || !is.finite(v_new) || v_new <= 0) return(0)
	if (!is.finite(scale) || scale <= 0) return(0)
	if (is.finite(v_max) && v_new > v_max) return(-Inf)
	-shape * log(v_new / v) - scale / v_new + scale / v
}
# log prior ratio plus Jacobian of the joint scale move, with the variance parameters held fixed.
# a Gaussian block of dimension d and quadratic form q moved by f contributes d log(f) - (f^2 - 1) q / 2
.lame_scale_prior_logratio <- function(
	cc,
	beta,
	iV0,
	q_lev,
	n_lev,
	q_dev = 0,
	n_dev = 0L,
	q_uv = 0,
	n_uv = 0L
) {
	lg <- log(cc)
	c2m1 <- cc * cc - 1
	out <- length(beta) * lg - 0.5 * c2m1 * sum(beta * as.numeric(iV0 %*% beta))
	out <- out + n_lev * lg - 0.5 * c2m1 * q_lev
	out <- out + n_dev * lg - 0.5 * c2m1 * q_dev
	# the positions move by sqrt(c), so their quadratic form moves by c
	out <- out + 0.5 * n_uv * lg - 0.5 * (cc - 1) * q_uv
	out
}
# redraw the binary latent array given the shared dyad factors, within the sign region Y implies
.lame_binary_z_given_factor <- function(
	Z,
	EZ,
	off,
	rho,
	idx1,
	idx0,
	idx_na,
	idx_diag
) {
	if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
		stop("Binary auxiliary correlation must lie strictly between -1 and 1.", call. = FALSE)
	ar <- abs(rho)
	sd_cell <- sqrt(1 - ar)
	Znew <- Z
	if (length(idx1) > 0L) {
		mu <- EZ[idx1] + off[idx1]
		lt <- pnorm(-mu / sd_cell, lower.tail = FALSE, log.p = TRUE)
		Znew[idx1] <- mu + sd_cell * qnorm(
			log(runif(length(mu))) + lt,
			lower.tail = FALSE,
			log.p = TRUE
		)
	}
	if (length(idx0) > 0L) {
		mu <- EZ[idx0] + off[idx0]
		lt <- pnorm(-mu / sd_cell, log.p = TRUE)
		Znew[idx0] <- mu + sd_cell * qnorm(
			log(runif(length(mu))) + lt,
			log.p = TRUE
		)
	}
	if (length(idx_na) > 0L) {
		mu <- EZ[idx_na] + off[idx_na]
		Znew[idx_na] <- mu + sd_cell * rnorm(length(mu))
	}
	if (length(idx_diag) > 0L) {
		Znew[idx_diag] <- rnorm(length(idx_diag), EZ[idx_diag], sqrt(1 + rho))
	}
	bad <- !is.finite(Znew)
	if (any(bad)) stop("Non-finite binary latent-response refresh.", call. = FALSE)
	Znew
}
#' AME model fitting routine for longitudinal relational data
#'
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to longitudinal (time-series) relational data of
#' various types. Supports both unipartite (square) and bipartite (rectangular)
#' network structures. For cross-sectional (single time point) networks,
#' use the \code{ame} function.
#'
#' @details
#' This command provides posterior inference for parameters in AME models of
#' longitudinal relational data, assuming one of six possible data
#' types/models. The model supports both unipartite networks (square adjacency
#' matrices) and bipartite networks (rectangular adjacency matrices with
#' distinct row and column node sets) across multiple time points.
#'
#' \strong{Tempered symmetric binary sampling:}
#'
#' Set \code{sampler = "tempered"} for a common actor grid with
#' \code{family = "binary"}, \code{symmetric = TRUE}, \code{dynamic_ab = TRUE},
#' positive static rank, a static intercept and optional symmetric dyadic
#' covariates. At least four actors and two equally spaced periods are required.
#' The observed likelihood contains each unordered dyad once:
#' \deqn{P(Y_{ijt}=1)=\Phi(x_{ijt}^{T}\beta+a_{it}+a_{jt}+u_i^T\Lambda u_j).}
#' Shared actor effects comprise a zero-sum level and a zero-sum stationary
#' Gaussian AR(1) deviation. Their level variance has an inverse-gamma prior
#' with shape \code{etaab/2} and scale \code{etaab * Sab0[1,1]/2};
#' the innovation variance has shape \code{sigma_ab_shape} and scale
#' \code{sigma_ab_scale}. Persistence has its specified normal prior truncated
#' to \eqn{(-0.99,0.99)}. Coefficients retain the panel-design Gaussian prior.
#'
#' \code{prior$uv_prior = "hierarchical"} retains the static symmetric factor
#' prior: signed weights are independent standard normals, and each raw
#' coordinate column has variance \eqn{v_k\sim IG(1,0.5)}. This is the default
#' for the tempered sampler. The explicitly selected alternative
#' \code{prior$uv_prior = "unit"} fixes every coordinate variance at one;
#' it changes the model's prior, not just its sampling algorithm. Raw factors
#' are centered for reporting, with compensating actor-level and intercept
#' shifts that preserve every joint predictor.
#'
#' The sampler uses bounded Langevin proposals with their complete forward
#' and reverse Gaussian densities, a fixed ladder of 21 likelihood powers,
#' replica exchanges, prior refreshes at power zero, and signed-weight
#' reflections. Step sizes adapt only during burn-in. \code{fit$tempering}
#' records transition counts, step sizes and replica sign traces;
#' convergence diagnostics are computed from the stored draws of
#' independently seeded fits.
#'
#' This route accepts only its documented factor and additive prior settings;
#' \code{Sab0} must be a positive multiple of the two by two identity and
#' \code{ab_var_cap} must be infinite. Nodal covariates, changing actor sets,
#' time-varying coefficients/factors, exposure, and chunked likelihood output
#' are not supported. \code{save_log_lik = TRUE} stores one column per observed
#' unordered dyad-period, also the unit counted by \code{nobs()}.
#' Posterior storage uses the same \code{posterior_opts} interface, and compact
#' joint terminal states are always retained for forecasts. LFO training
#' windows need at least two periods. Tempered checkpoints restore every
#' replica, the frozen model/prior and both random-number streams; only
#' sampling length, output and time-limit settings can change on continuation.
#' Built-in RNGs with serializable state are supported; Box-Muller and
#' user-supplied RNGs are refused with an error.
#'
#' \strong{Dynamic Effects Implementation:}
#'
#' The dynamic_uv and dynamic_ab parameters enable time-varying latent representations
#' through autoregressive processes. These extensions are particularly useful for
#' understanding how network structure evolves over time.
#'
#' \emph{Dynamic Multiplicative Effects (dynamic_uv=TRUE):}
#'
#' The latent factors U and V evolve according to the transition selected by
#' \code{dynamic_uv_kind}. The default \code{"ar1"} model uses Gaussian AR(1)
#' drift:
#' \deqn{U_{i,k,t} = \rho_{uv} U_{i,k,t-1} + \epsilon_{i,k,t}}
#' where \eqn{\epsilon_{i,k,t} \sim N(0, \sigma_{uv}^2)}, i indexes actors, k indexes
#' latent dimensions, and t indexes time. The parameter \eqn{\rho_{uv}} controls
#' temporal persistence (values near 1 indicate slow evolution). This captures
#' time-varying homophily, latent community structure, and transitivity dynamics.
#' The \code{"snap"} transition uses a mixture of AR(1) drift and discontinuous
#' reset transitions, while \code{"t"} uses heavy-tailed Student-t innovations
#' represented by local transition scales. On the MCMC path, snap and t
#' transitions are supported for directed unipartite models. With
#' \code{method = "als"}, Student-t dynamic UV is supported for directed and
#' symmetric normal, binary, and Poisson panels, including named panels where
#' actors enter or exit; the snap kind has no ALS route.
#'
#' Key references:
#' \itemize{
#'   \item Sewell & Chen (2015): Introduced dynamic latent space models with
#'         actor-specific evolution rates
#'   \item Durante & Dunson (2014): Nonparametric Bayesian approach allowing
#'         flexible evolution of network structure
#'   \item Hoff (2011): Hierarchical multilinear models providing theoretical
#'         foundation for temporal dependencies
#' }
#'
#' \emph{Dynamic Additive Effects (dynamic_ab=TRUE):}
#' Each sender (a) and receiver (b) effect is a per-actor level plus a
#' zero-mean AR(1) deviation:
#' \deqn{a_{i,t} = \mu^a_i + d_{i,t}, \quad d_{i,t} = \rho_{ab} d_{i,t-1} + \epsilon_{i,t}}
#' \deqn{b_{j,t} = \mu^b_j + e_{j,t}, \quad e_{j,t} = \rho_{ab} e_{j,t-1} + \eta_{j,t}}
#' where \eqn{\epsilon_{i,t}, \eta_{j,t} \sim N(0, \sigma_{ab}^2)}. The variance
#' of the levels across actors is the reported \code{va} / \code{vb};
#' \eqn{\rho_{ab}} and \eqn{\sigma_{ab}} describe the deviations. This models
#' time-varying individual activity levels (outdegree) and popularity (indegree).
#'
#' \emph{Identified reporting:}
#' Stored static factor draws are centered only when both resolved actor
#' margins are enabled. The removed row, column and grand-mean components of
#' \eqn{U G V'} are transferred to the actor effects and intercept; for a
#' symmetric shared actor without an intercept, the grand constant is divided
#' equally between its two contributions. Otherwise, the saved draws retain
#' their sampled coordinates because disabled effects cannot absorb the offsets.
#' The reported map preserves the linear predictor and leaves the sampler state
#' unchanged. Dynamic-coordinate reporting follows the supported fitted route;
#' inspect fitted products and predictions when comparing latent geometries.
#'
#' Dynamic additive effects are useful when sender activity and receiver
#' popularity change over time rather than staying fixed across the whole panel.
#' A \code{dynamic_ab} fit returns the posterior-mean paths on
#' \code{fit$a_dynamic} / \code{fit$b_dynamic} (actor x period) and their
#' per-period posterior standard deviations on \code{fit$a_dynamic_sd} /
#' \code{fit$b_dynamic_sd}; \code{ab_plot(fit, plot_type = "ribbon")} uses the
#' latter to draw a credible band around each actor's path.
#'
#' \emph{Prior Specification for Dynamic Parameters:}
#' \itemize{
#'   \item \eqn{\rho_{uv}, \rho_{ab} \sim TruncNormal(mean, sd, -0.99, 0.99)}:
#'         Ensures stationarity of AR(1) process
#'   \item \eqn{\sigma_{uv}^2, \sigma_{ab}^2 \sim InverseGamma(shape, scale)}:
#'         Controls innovation variance
#'   \item Default priors (\eqn{\rho_{uv}} mean=0.9, \eqn{\rho_{ab}} mean=0.8)
#'         favor smooth evolution
#'   \item Adjust rho_*_mean closer to 1 for slower evolution, closer to 0 for
#'         more rapid changes
#' }
#'
#' \emph{Computational Considerations:}
#' \itemize{
#'   \item Dynamic effects add path updates; runtime depends on actor count, rank and number of periods
#'   \item Memory usage scales as O(n*R*T) for dynamic_uv, O(n*T) for dynamic_ab
#'   \item C++ implementation is substantially faster than a pure R loop
#'   \item Convergence diagnostics: Monitor rho and sigma parameters carefully
#'   \item Effective sample sizes typically lower due to temporal correlation
#'   \item A chain that is still moving when the run ends is continued, and refused if it keeps moving (see the convergence policy below)
#' }
#'
#' \emph{Chains and convergence policy:}
#' \code{lame()} runs \code{chains} chains (default 2; \code{options(lame.chains = k)}
#' sets another default). Chain \eqn{k} is seeded with \code{seed + 1000 * (k - 1)},
#' so \code{chains = 1} is the first chain of a two-chain fit. On Unix the chains
#' run forked, on at most \code{cores} workers (all cores but one by default); on
#' Windows they run one after the other. The draws do not depend on the number of
#' workers. The result is the pooled fit of \code{\link{combine_ame_chains}}, with
#' \code{fit$chains}, \code{fit$chain_seeds}, \code{fit$rhat} (the R-hat of the
#' monitored set), the state each chain stopped in, and the \code{chain_indicator}
#' that \code{as_draws()} and \code{\link{trace_plot}} use. Two chains
#' cost twice the CPU time and about the same wall time on a multicore machine.
#'
#' At the end of the post-burn run every chain is checked, and the chains are
#' compared with each other. Within a chain, the first and last thirds of the
#' stored draws of the intercept (or first coefficient), \code{va}, \code{vb},
#' \code{ve}, \code{rho}, \code{rho_ab}, \code{sigma_ab}, \code{rho_uv},
#' \code{sigma_uv} and \code{pi_snap} are compared with batch-means standard errors
#' (ten batches per third); a quantity is flagged when the two means differ by more
#' than 3 standard errors and by more than half a posterior standard deviation of the
#' whole series. Between chains, the split rank-normalised R-hat (Vehtari et al.,
#' 2021) of the monitored set is compared with 1.05: every regression coefficient
#' (each period's under \code{dynamic_beta}), \code{va}, \code{vb}, \code{cab},
#' \code{ve}, \code{rho}, the dynamic hyperparameters, and summaries of the
#' additive and latent structure recorded at every stored draw (the rms of
#' \code{a} and \code{b}, the largest actor effect, the rms of the latent product
#' \eqn{UGV'}, the norm of \code{G} of a bipartite fit and the split of scale
#' between \code{U} and \code{V} of a directed fit, per period for dynamic fits).
#' The dyadic products are not monitored. Quantities a fit reports as not estimated
#' are left out, and a run with fewer than 60 stored draws per chain is not checked.
#' The gate is computed on the draws the returned fit holds (the last block of
#' every chain). \code{fit$diagnostics$rhat}, from
#' \code{\link{compute_mcmc_diagnostics}}, is the same statistic of the same
#' draws and equals \code{fit$rhat} for every monitored coefficient and variance
#' component.
#'
#' When a chain drifts or an R-hat exceeds 1.05, every chain is extended by a block
#' of \code{nscan} sweeps from its own final state, as \code{\link{lame_resume}} does
#' (the seeds advance by the sweeps already run); the flagged block counts as extra
#' burn-in, its draws are discarded, and both checks run again. Up to
#' \code{prior$chain_extensions} extensions are made (default 3, so at most four
#' times \code{nscan} post-burn sweeps per chain). The returned fit has the same
#' number of stored draws as a fit that needed no extension, and
#' \code{fit$auto_extend} records the extensions, the sweeps discarded and the
#' quantities that triggered each. If the last block is still flagged, \code{lame()}
#' stops with an error of class \code{lame_chain_unsettled}: its \code{fits} field
#' holds the chains, \code{fit} the pooled fit of the last block, \code{rhat} the
#' monitored table and \code{moving} the offending quantities, and
#' \code{lame_resume(e$fit, nscan = nscan)} continues every chain
#' (\code{lame_resume(e$fits[[k]], nscan = nscan)} one). \code{prior = list(chain_extensions
#' = 0)} refuses at the first flag; \code{chains = 1} runs one chain with the drift
#' check only. There is no mode that returns a flagged fit.
#'
#' Measured on chains that were already stationary (112 blocks of 400 draws from 28
#' simulated configurations, two chains each), the largest monitored R-hat had a
#' median of 1.007 and a maximum of 1.040, and no block was flagged by either
#' check (the 95 percent upper bound on the false-extension rate is 3.2 percent;
#' three consecutive false flags have probability below \eqn{4 \times 10^{-5}}).
#' Measured on the 13 simulated configurations of the replication study that
#' failed a two-chain R-hat gate in earlier builds (four pairs of seeds each, 52
#' two-chain fits at the defaults, \code{burn = 500}, \code{nscan = 10000},
#' \code{odens = 25}), the policy settled 47 fits as run, extended 5 (four once,
#' one twice) and refused none. 47 of the 52 pass the book's gate, an R-hat below
#' 1.05 on every identified parameter including every dyadic product; the other five
#' sit at R-hat 1.04 to 1.08 on dyadic products or one actor effect, which the
#' monitored set leaves out.
#'
#' The dyadic products are left out because their largest R-hat crosses 1.05 on
#' settled chains. A fit has about \eqn{10^4} of them, and on stationary chains of
#' 17 simulated configurations with a latent product (68 two-chain blocks) the
#' largest product R-hat of a block reached 1.055; settled ordinal directed fits
#' with \code{dynamic_uv} crossed 1.05 in 2 of their 4 blocks, so a product gate at
#' 1.05 would extend and refuse settled fits. In the three acceptance configurations
#' whose products failed the book's gate at the defaults (a Poisson undirected static
#' rank-2 fit, and frn and ordinal directed fits with \code{dynamic_uv}), the
#' smallest product effective sample size of a failing run was 19 to 87 of 800
#' draws. Run at four times the default \code{nscan}, all three passed, and the
#' posterior-mean product \code{fit$UVPM} moved by at most 0.27 posterior standard
#' deviations on any dyad. A study that reports individual dyadic products, rather
#' than \code{fit$UVPM} or the monitored summaries, keeps the position draws
#' (\code{posterior_opts = list(save_UV = TRUE)}) and continues the chains with
#' \code{lame_resume()} until the R-hat of every reported product is below 1.05.
#'
#'
#' \emph{Model Selection Guidelines:}
#' Use both dynamic_uv and dynamic_ab when:
#' \itemize{
#'   \item Networks show clear temporal trends in density or clustering
#'   \item Individual node behavior changes systematically over time
#'   \item Community structure evolves (merging, splitting, drift)
#' }
#'
#' Use only dynamic_uv when:
#' \itemize{
#'   \item Latent structure/communities change but individual effects are stable
#'   \item Focus is on evolving homophily or clustering patterns
#'   \item Network shows structural reconfiguration over time
#' }
#'
#' Use only dynamic_ab when:
#' \itemize{
#'   \item Individual heterogeneity varies but overall structure is stable
#'   \item Actors' activity/popularity changes over observation period
#'   \item Focus is on individual-level temporal dynamics
#' }
#'
#' \strong{Bipartite Network Models:}
#'
#' When mode="bipartite", the model handles rectangular adjacency matrices Y with
#' dimensions n_A x n_B, where n_A and n_B represent the number of row and column
#' nodes respectively.
#'
#' \emph{Static Bipartite Case:}
#' The model uses separate latent factor matrices:
#' \itemize{
#'   \item U: n_A x R_row matrix of row node latent positions
#'   \item V: n_B x R_col matrix of column node latent positions
#'   \item G: R_row x R_col interaction matrix mapping between latent spaces
#'   \item Multiplicative term: U G V' captures bipartite community structure
#' }
#'
#' \emph{Dynamic Bipartite Case:}
#' When dynamic_uv=TRUE for bipartite networks:
#' \deqn{U_{i,k,t} = \rho_{uv} U_{i,k,t-1} + \epsilon_{i,k,t}}
#' \deqn{V_{j,k,t} = \rho_{uv} V_{j,k,t-1} + \eta_{j,k,t}}
#' where i indexes row nodes, j indexes column nodes, k indexes latent dimensions.
#'
#' A time-varying bipartite interaction matrix is not supported:
#' \code{dynamic_G = TRUE} is refused on both the MCMC and the ALS route. With
#' moving positions a per-period interaction matrix is not identified apart
#' from them -- \eqn{U_t G_t V_t'} and \eqn{(U_t G_t) V_t'} describe the same
#' products, so a fit cannot separate movement in \eqn{G_t} from movement in
#' \eqn{U_t} and \eqn{V_t}. With static positions, \eqn{U G_t V'} is identified
#' as a static \eqn{G} is, but it is not implemented. Use
#' \code{dynamic_uv = TRUE} for a time-varying multiplicative
#' term and read the per-period product \code{fit$UVPM}, which is the identified
#' quantity; every bipartite fit still estimates a single interaction matrix and
#' returns it as \code{fit$G}.
#'
#' \emph{Key Differences from Unipartite Models:}
#' \itemize{
#'   \item No dyadic correlation (rho): Bipartite edges are inherently directed
#'   \item Separate dimensions: R_row and R_col can differ for row/column spaces
#'   \item Rectangular structure: Network density patterns differ from square matrices
#'   \item Community interpretation: Captures affiliation patterns between node types
#' }
#'
#' \strong{Missing covariates:}
#'
#' A cell whose dyadic covariate, sender (row) covariate or receiver (column)
#' covariate is missing in its period is an unobserved tie, for unipartite and
#' bipartite panels alike: it is set to \code{NA} in \code{Y} before the design
#' is built, and the covariate coefficients are estimated from the complete
#' cells. A message gives the number of such cells. A missing value is never
#' modelled as 0, which would bias the coefficient of a covariate that is not
#' centred. In a symmetric panel a dyad has one covariate value: a hole in one
#' triangle of a dyadic covariate is filled from the other triangle (a message
#' gives the number of values), so the slice stays symmetric and the dyad stays
#' observed, and a dyad missing in both triangles is unobserved. A missing node
#' covariate removes the actor's row and column in that period. A covariate
#' that is missing in every period stops the fit. \code{\link{ame}} and the ALS
#' routes apply the same rule and give the same messages.
#'
#' \strong{Collinear covariates:}
#'
#' A regression design whose columns are exactly aliased over the observed
#' cells of the whole panel (a covariate that duplicates another, or one that
#' is a linear combination of others, counting the intercept and the sender and
#' receiver covariates), static or with \code{dynamic_beta}, is an error of
#' class \code{lame_design_rank}. It names the aliased columns and the
#' combination each equals: the data cannot separate their coefficients, so the
#' posterior along that direction would be the prior's alone. The periods are
#' pooled, so a covariate that duplicates another in some periods only is
#' fitted, and period indicators that span the intercept still drop it, with the
#' warning. A column is aliased when it is a combination of the columns before it
#' to a relative 1e-8 (as \code{lm()} reports them); a design that is only
#' nearly collinear is fitted. \code{method = "als"} refuses the same designs
#' with the same condition (fields \code{aliased}, \code{rank} and
#' \code{n_col}) for every estimator and family.
#'
#' \strong{Standard AME Model Types:}
#'
#' The following describes the six standard data types/models available:
#'
#' "normal": A normal AME model.
#'
#' "binary": A binary probit AME model.
#'
#' "ordinal": An ordinal probit AME model. An intercept is not identifiable with
#' data-induced cutpoints. Explicit cutpoints anchored at zero support an intercept.
#'
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
#'
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#'
#' For both families a missing cell (\code{NA}) of an actor present in the
#' period is an unobserved nomination or non-nomination: it counts against
#' its row's cap while it holds a nomination, and the sampler draws it from
#' its conditional given the row's observed cells. A row at its cap is bounded
#' below only, so the sweep moves its level slowly; each sender's level is
#' therefore also moved together with its latent row once per sweep, drawn
#' from its conditional. On simulated designs (30 actors and one period, 24
#' actors and two to six, cap 10, more than 3,000 fits) the 95 percent intervals of
#' \code{vb} and \code{cab} covered 0.97 to 1.00 of the panels at every share
#' of rows at the cap, with and without a fifth of the cells hidden, and those
#' of \code{va} 0.90 to 1.00 outside the region \code{odmax} describes and
#' 0.91 inside it. Two designs withhold the level of the rows (see
#' \code{odmax}): one period, where a row at its cap shows only that
#' its level fills the cap, with every row at the cap, or with cells hidden
#' and a share at the cap from 0.5 (the intercept interval covered 0.64 and
#' 0.84 of the panels; \code{ame()} and one-period \code{lame()} report the
#' intercept and any row-covariate coefficient as not estimated then); and two
#' periods or more with a share from 0.5 and a posterior-mean \code{va} from
#' 0.5 to 1.5. At a generating sender sd of 0.25 the posterior mean of
#' \code{va} runs up to 2.3 times the generating variance while the interval
#' still covers. \code{fit$hidden_cell_ess}
#' records the effective sample size of the intercept and \code{va} of a fit
#' with hidden cells. The posterior predictive
#' (\code{YPM}, \code{fitted()}, the goodness-of-fit replicates and
#' \code{simulate()}) lets a missing cell hold a nomination, as the sampler
#' does, so \code{YPM} also predicts those cells; the diagonal and the cells
#' of an actor absent from the period are left out.
#'
#' "poisson": An overdispersed Poisson-lognormal model for counts (negative
#' or non-integer observed values are an error):
#' \eqn{Y \sim \mathrm{Poisson}(e \exp(z))} with
#' \eqn{z \sim N(\eta, \sigma^2)}. Given \eqn{\eta} and \eqn{\sigma^2}, a
#' cell's mean count is \eqn{e \exp(\eta + \sigma^2/2)}. \code{YPM},
#' \code{fitted()} and \code{predict(type = "response")} report it, for the
#' exposure \eqn{e_t} of the period (\code{period_exposure}; one without), at
#' the posterior mean \eqn{\eta} of the linear predictor (\code{EZ}) and the
#' posterior mean \eqn{\sigma^2} of the residual variance, the form the fast
#' estimators report. Predictive quantiles and event probabilities come from
#' posterior predictive draws (\code{simulate()}, \code{predict(fit, h = 1, interval = "prediction")}).
#' Zero exposure gives a structural zero count.
#'
#' \strong{Joint scale move (binary, unipartite and bipartite):}
#'
#' A probit panel identifies the scale of the linear predictor only weakly:
#' a higher intercept with less spread in the additive effects and smaller
#' latent positions fits a sparse panel nearly as well, so the level of the
#' predictor and the spread of the effects trade off along a ridge. The
#' conditional sweep cannot step along it, because conditioning on the latent
#' \eqn{Z} pins the scale of every block, and a cold start can take many
#' thousands of sweeps to arrive. Binary unipartite fits therefore carry a
#' Metropolis step that multiplies the coefficients and additive effects by a
#' factor \eqn{c} and the latent positions by \eqn{\sqrt{c}}, with the
#' variance parameters held fixed over the step and following it through their
#' own draws. Under \code{dynamic_ab} the same visit walks a second direction:
#' the per-period deviations of the additive effects together with
#' \eqn{\sigma_{ab}}, whose amplitude relative to that innovation scale is a
#' ridge of its own (the deviation draw pins the deviations given
#' \eqn{\sigma_{ab}} and its draw pins \eqn{\sigma_{ab}} given them). Moving
#' the two together leaves their AR(1) prior invariant, so only the likelihood
#' and the inverse-gamma prior on \eqn{\sigma_{ab}^2} resist; set
#' \code{prior$dev_move = FALSE} to leave that direction out. The step
#' integrates \eqn{Z} out (the dyad errors are split into a shared factor and
#' an independent part; the factor is drawn from its conditional and the
#' independent part integrated out), so it is accepted against the exact
#' conditional posterior and leaves it unchanged; \eqn{Z} is redrawn from its
#' conditional when the step is accepted. It is visited every \code{prior$level_move_every} sweeps
#' (default 30) with \code{prior$level_move_tries} proposals per visit
#' (default 8), adapts its step size during burn-in only, and reports its
#' acceptance rate in \code{fit$level_move}. It is off for other families,
#' symmetric fits and \code{dynamic_rho}.
#' Changing coefficients with static factors are supported for ranks 1 and 2
#' when both actor margins and \code{dynamic_ab} are enabled,
#' \code{dcor = FALSE}, and the coefficient law is AR(1) or RW1 on the
#' regular grid, with zero initial means and diagonal innovation scaling.
#' The static factor covariance is integrated in the acceptance density.
#'
#' That configuration also receives a joint Hamiltonian update of coefficients,
#' actor effects and factors with the binary latent responses integrated out.
#' It requires at least three actors and one static coefficient. The update
#' uses exact actor contrasts, keeps hyperparameters fixed, and adapts its
#' step size during burn-in only. Bounded momentum kicks preserve reversal and
#' phase-space volume; the Metropolis correction retains the original target.
#' It runs every 20 sweeps; acceptance and
#' numerical rejections are recorded in \code{fit$marginal_move}. Set
#' \code{prior$marginal_move = FALSE} to omit this update.
#' With static additive effects, an intercept, both actor margins and the
#' same rank-one or rank-two changing-coefficient configuration, a joint
#' observed-data update moves static actor contrasts, coefficients and factors.
#' It integrates out the static factor covariance under its original prior.
#' Bounded momentum updates and a Metropolis correction preserve the original
#' target. The move begins after half of burn-in and runs every 20 sweeps;
#' its step size adapts only during burn-in. Acceptance and numerical
#' rejections are recorded in \code{fit$marginal_move}. Set
#' \code{prior$marginal_move = FALSE} to omit this static-actor update.
#'
#' Other changing-coefficient and static-factor configurations retain their
#' existing updates.
#'
#' Binary bipartite fits with an intercept, both margins, static coefficients
#' and, when they have a multiplicative term, a static interaction matrix carry
#' the same step on the whole predictor: the coefficients, the additive effects
#' (and their deviations under \code{dynamic_ab}) and the interaction matrix
#' \eqn{G} are multiplied by \eqn{c} with the latent positions held. The cells
#' of a bipartite panel have no dyadic correlation, so no shared factor is
#' drawn, and the prior terms of every block the factor touches (the
#' coefficient prior, the effects, their AR(1) deviations, the \eqn{N(0, 1)}
#' prior of the entries of \eqn{G} cut at \code{prior$G_max}) are included in
#' the acceptance ratio.
#'
#' Binary unipartite fits (directed or symmetric, with any additive margin in
#' the model, and with static coefficients or \code{dynamic_beta} of kind
#' \code{"ar1"} or \code{"rw1"}) also receive a collapsed per-actor
#' update. The data-augmentation sweep draws an actor's
#' parameters from a normal whose precision is the number of cells the actor
#' appears in, so an actor whose cells hold no events (a flat likelihood on a
#' half line) moves by about \eqn{1/\sqrt{m}} per sweep however wide its
#' posterior is, and relaxes in about \eqn{v m} sweeps: thousands of sweeps
#' for the additive level of an inactive actor, and \eqn{10^5} for static
#' positions that have spread out, because the positions also feed the factor
#' covariance they are drawn under. The update integrates the latent array out
#' (as the scale move does) and draws, one actor at a time, an exact elliptical
#' slice step of the actor's sender and receiver level translations and, for
#' static positions, of its sender and receiver positions. The level
#' translation adds a constant to one actor's predictor while the other levels
#' and the intercept absorb the re-centring, so it keeps the sum-to-zero
#' constraint and touches one row (or column) of the likelihood. The latent
#' array is redrawn afterwards, so the posterior is unchanged; there is
#' nothing to tune. It is available for \code{dynamic_ab} and static effects,
#' ranks 0 to any with static positions (positions that move over time keep
#' their own updates; only the levels are moved), and is visited every
#' \code{prior$actor_move_every} sweeps (default 10). The positions are left to
#' the sampler during the first half of burn-in, when it runs under a fixed
#' diffuse prior that lets the positions of actors without information take any
#' scale. In a symmetric fit an actor has one level, which enters both its row
#' and its column, and the positions are those of the eigenmodel (\eqn{u_i}
#' with the diagonal weights held fixed during the visit, drawn under the
#' sampler's own precision prior); there is no dyadic correlation. Panels with
#' absent actors (changing composition) are covered: only observed cells enter
#' a likelihood, an actor
#' with no observed cell in a block has the prior as its conditional there and
#' the slice step moves it across that prior, the translation moves the level
#' of absent periods with the rest, and the missing cells are redrawn from the
#' model afterwards, as the sampler does.
#' Set
#' \code{prior$actor_move = FALSE} to omit it; its work is recorded in
#' \code{fit$actor_move}. For directed fits with \code{dcor = TRUE} the same
#' visit also updates the dyadic correlation with both latent responses of each
#' dyad integrated out (an independent proposal from its \eqn{1 - \rho^2}
#' prior, accepted on the exact pair likelihood, followed by an exact redraw of
#' the latent array), because given the augmented array the correlation moves
#' by about \eqn{1/\sqrt{\mathrm{dyads}}} and mixes slowly when few dyads are
#' mutual; \code{prior$actor_move_rho = FALSE} omits that part.
#'
#' The Poisson and ordinal families have the same slow direction, and the
#' latent array cannot be integrated out, so their update carries it along:
#' an actor's level is translated together with the latent values of its row
#' (its column, or the pairs of a symmetric network), which leaves every
#' residual, hence the Gaussian density of the latent array, unchanged. The
#' shift is then an exact draw from the level prior (with the intercept that
#' absorbs the re-centring) times the data factor of the row: the log-concave
#' Poisson factor of the observed cells (an elliptical slice step) or, for the
#' ordinal family, the indicator that the shifted values keep their ranks (or
#' stay inside their cutpoint intervals), which cuts the prior at both ends
#' and is drawn as a truncated normal. The intercept, the spread of the
#' levels and the levels of the inactive actors, which the sweep alone crosses
#' over thousands of sweeps on a sparse panel, then mix within a few visits.
#' With static coefficients the Poisson family also translates each
#' coefficient with the latent array, and moves the dyadic variance and the
#' dyadic correlation together with the latent residuals (a Metropolis step
#' whose ratio is the change of the Poisson factor and the priors), because
#' given the array the residuals are prior draws wherever the counts say
#' little, and, with static positions, redraws each actor's position with the
#' latent values of its row or column carried along (an elliptical slice step
#' under the factor prior the sampler draws under; the positions are left to
#' the sampler during the first half of burn-in, when it runs under a fixed
#' diffuse prior); the ordinal family draws the location of the latent array
#' (the means of the levels, or the intercept), which the rank likelihood
#' leaves free and explicit cutpoints confine to a strip. Visits are every
#' \code{prior$actor_move_every} sweeps (default 10).
#'
#' With moving positions (\code{dynamic_uv = TRUE} with a positive rank, of the
#' \code{"ar1"} kind, or the \code{"t"} kind of a directed network) the sweep
#' moves the position path of an actor whose cells say little (no events, or
#' absent from most periods) by a step of about \eqn{\sqrt{s_2 / m}} per sweep,
#' \eqn{m} the number of its cells, however wide the posterior of the path is
#' (\eqn{s_2} is the dyadic variance, near 0.05 for a sparse Poisson panel and 1
#' for the probit families). The paths of such actors then wander at no likelihood
#' cost, and the squared product of the positions, the persistence and
#' innovation scale of the paths and the intercept follow them with
#' autocorrelation times of hundreds to thousands of sweeps. Every
#' \code{prior$actor_move_every} sweeps the update therefore redraws each actor's
#' whole path, the senders' and then the receivers' (a symmetric network's
#' once), with one elliptical slice step under the AR(1) prior of the path
#' (stationary start; the local innovation scales of the \code{"t"} kind), the
#' persistence and the innovation scale being held. For the Poisson family the
#' latent values of the actor's row or column are moved with the predictor, which
#' changes no residual, so the conditional of the path is its prior times the
#' Poisson factor of the cells. For the binary family the latent array is
#' integrated out as in the collapsed update above (the shared dyad factor of a
#' directed network is drawn once per visit, each observed cell contributes a
#' probit factor, and the latent array is redrawn at the end). Bipartite fits use
#' their interaction matrix and symmetric fits the signed weights of their
#' unit-variance coordinate paths. A directed binary fit with
#' \code{dynamic_ab} at rank two already has the paired-actor slice step
#' described below and is left to it; the \code{"snap"} kind and ordinal fits
#' keep the sampler's own updates. A symmetric Poisson or
#' binary fit also draws, every sweep, the scale ridge between its coordinate
#' paths and their weights (scaling a coordinate of every path by \eqn{c} and its
#' weight by \eqn{1/c^2} changes no cell of the predictor, so the draw is a
#' Metropolis step on \eqn{\log c} with the ratio of the priors; it is
#' recorded as \code{dpm_scale_attempts} and \code{dpm_scale_accepts}) and draws
#' the persistence \code{rho_uv} with a proposal scale taken from the
#' transitions of the paths, where the fixed scale of the sampler's own draw is
#' 10 to 90 times the width of the conditional at persistences from 0.9 to 0.985
#' and moves \code{rho_uv} in one call of every ten to ninety.
#' \code{prior$actor_move_paths = FALSE} omits this part only.
#'
#' With \code{dynamic_beta} the prior of the intercept enters the density of a
#' level translation through its path: the translation moves the whole path by
#' the same constant, which changes the prior of the first state and, unless
#' the path is a random walk, of the transitions. A static intercept next to
#' moving coefficients enters through its own prior, and moving coefficients
#' other than the intercept are not touched. A fit with no additive margin
#' (\code{rvar = cvar = FALSE}, or a symmetric fit with \code{nvar = FALSE}) and
#' static positions runs the position steps alone, and a fit with one margin
#' translates the levels of that margin. Coefficient paths of the other kinds
#' (\code{"rw2"}, \code{"matern32"}) or on a non-unit period grid keep their
#' own updates.
#'
#' The same update runs for bipartite binary fits (and in \code{\link{ame}}
#' for a bipartite cross-section): each row and column actor's level
#' translation, which moves the intercept with it and touches one row or
#' column of the likelihood, and, for static positions, its position, as exact
#' elliptical slice draws with the latent array integrated out; under
#' \code{dynamic_ab} it also draws the deviation paths of random pairs of actors
#' of a margin with their sum held. A bipartite
#' cell has no partner cell, so there is no reciprocity update, and the
#' positions have a fixed unit-variance prior, so there is no factor covariance
#' to draw. Missing cells and the cells of an actor absent from a period carry
#' no likelihood, and the update is never switched off by missingness. With
#' \code{dynamic_beta} of kind \code{"ar1"} or \code{"rw1"} on the unit period
#' grid the translation moves the whole intercept path by the same constant
#' (a static intercept next to moving coefficients by the same amount), and the
#' prior of the intercept enters the density of the translation through its
#' path, as for unipartite fits.
#'
#' A \code{"poisson"} fit translates the coefficients that move as well. The
#' counts of a period fix the log of its mean to about
#' \eqn{1/\sqrt{\mathrm{events}}}, while given the latent array the intercept of
#' the period is known to \eqn{\sqrt{s^2/\mathrm{cells}}}, 10 to 20 times less on a
#' sparse panel (0.015 to 0.022 against 0.0012 at 170 actors); the sweep moved the
#' intercepts of the periods, and with them the coefficient of a dyadic covariate
#' and the dyadic variance and correlation, by random walk. The gate panel of 170
#' actors and four periods had 3 to 14 effective draws of the period intercepts in
#' 1,600 (split R-hat 1.5), 6 of the covariate coefficient, 6 of the correlation and 16 of
#' the variance. Raising coefficient \eqn{k} of one period by \eqn{\delta} and the
#' latent values of the period by \eqn{\delta x_{kt}} leaves every residual unchanged,
#' so the shift is an exact draw (an elliptical slice step) from the prior of the
#' shifted path times the log-concave Poisson factor of the period's cells. The prior
#' is the one of the path's kind (\code{"ar1"}, \code{"rw1"}, \code{"rw2"},
#' \code{"matern32"}, on any time grid), read from the precision the path sampler
#' uses. The whole path of a coefficient is translated the same way, and so are the
#' static coefficients next to moving ones, with the latent array of every period, at
#' the visits of the translations of the actor levels
#' (\code{prior$poisson_period_move = FALSE} and
#' \code{prior$poisson_static_move = FALSE} omit them). Those fits, and every other
#' \code{"poisson"} fit without \code{dynamic_rho}, then move the dyadic variance and
#' correlation at every sweep with the latent residuals carried along their
#' conditional law: given the counts the residuals of a dyad are close to a Gaussian
#' with mean \eqn{m} and covariance \eqn{V} (the Gaussian prior times the first-order
#' terms of the Poisson factors), the move maps \eqn{e' = m' + L'L^{-1}(e - m)} for a
#' proposed variance and correlation, and accepts with the Metropolis ratio of the
#' joint density (the change of the Poisson factors, of the Gaussian density and
#' the Jacobian of the map). The residuals of the cells with a count are tilted by about
#' \eqn{s^2 (y - \lambda)}, which follows the variance and not its square root, so scaling
#' them as the translations of the actor levels do accepts steps of a few hundredths in
#' the log variance, and this map accepts steps of 0.1 to 0.8. The visits of the translations of
#' the actor levels are every second sweep for the Poisson family (every tenth left the
#' least active actors with 30 to 80 effective draws of 2,000); \code{fit$poisson_period_move}
#' records the moves. On the gate panel the period intercepts then have 1,370 to 1,450 effective draws of
#' 1,600, the covariate coefficient 1,320, the correlation 650, the variance 850 and the
#' slowest actor effect 400, and the split R-hat is 1.006; a sweep takes about a third longer.
#'
#' An \code{"ordinal"} fit draws the scale of its latent array along the transformations that
#' keep every rank. With the default data-induced cutpoints a cell of the latent array is drawn between
#' the extreme values of the neighbouring categories, and those extremes are order statistics of
#' thousands of cells: they move by steps of \eqn{1/\mathrm{cells}} per sweep, while the scale of the
#' whole configuration (the spread of the implied cutpoints against the unit residual variance) has a
#' posterior of width \eqn{1/\sqrt{\mathrm{cells}}}, and given the latent array the coefficients, the
#' additive effects and the dyadic correlation are tied to it. On 150 actors the coefficient of a dyadic
#' covariate had 8 effective draws of 600 and the correlation 95; on 100 actors in four periods, whose
#' rank windows are separate, 3 of 400 and 3. Multiplying the latent array, the coefficients and the
#' static additive effects by a factor \eqn{c} keeps every rank, and \eqn{c} has an exact
#' conditional (a group move whose invariant measure is \eqn{dc/c}: a gamma draw when the priors
#' have mean zero and the predictor is fully scaled, an independence Metropolis step from its Laplace
#' approximation otherwise). The shift and the scale of the latent array of each period alone keep that
#' period's ranks and have exact conditionals too. Both are drawn at every sweep, and the same panels then
#' have 600 and 400 effective draws of the coefficient and 400 and 270 of the correlation. The
#' multiplicative effects, the dynamic additive effects and the dynamic positions stay unscaled
#' parts of the predictor: the moves stay exact, and their efficiency is the share of the signal they
#' scale. With explicit cutpoints the scale of the whole state carries the free cutpoints with it; the
#' cutpoints of a directed fit with \code{dcor = TRUE} are drawn given the latent array, and on 50 actors
#' the coefficient had 56 effective draws of 500 and the correlation 47, now 440 and 290
#' (\code{prior$ordinal_scale_move = FALSE} omits the moves; \code{fit$ordinal_scale_move} records
#' them).
#'
#' A \code{"binary"} fit draws the scale of its linear predictor with the latent array integrated
#' out. When a covariate separates the ties (every tie lies above every non-tie on it) the
#' likelihood increases without a bound along the ray that multiplies the coefficients and the
#' additive effects by a factor \eqn{c > 1}, and the posterior of the coefficient is its prior
#' times a likelihood that saturates only where the dyads nearest the boundary are classified.
#' The data-augmentation sweep draws every block given the latent array and crosses that
#' direction by steps of \eqn{1/\sqrt{\mathrm{cells}}}: on a 30-actor directed panel the
#' separating coefficient had 12 effective draws of 2,000 (steps of 0.03 against a posterior
#' sd of 0.6), and on a bipartite panel, whose coefficients have the independent prior
#' \eqn{N(0, g)} with \eqn{g} the number of cells, it was still climbing after 8,000 sweeps
#' (14.0 and 15.6 for two seeds, against a posterior mean of 59). Every
#' \code{prior$binary_scale_every}-th sweep (default 5) the coefficients and the static additive
#' effects are therefore multiplied by a factor drawn exactly from its conditional density, a
#' group move whose invariant measure is \eqn{dc/c}: the dyad errors of a directed network are
#' split into a shared factor, drawn from its conditional given the latent array, and independent
#' parts that are integrated out, so that each cell contributes one probit factor, and the latent
#' array is redrawn from its conditional afterwards (the device of the collapsed actor updates).
#' The multiplicative effects, the dynamic additive effects and the actor-specific slopes stay where
#' they are and the draw stays exact; its efficiency is the share of the signal it scales. The same
#' panels then have 460 and 140 effective draws of 2,000 after 8,000 sweeps, two seeds agree, and the
#' posterior mean is the same at 2,000, 8,000 and 32,000 sweeps. The draw runs in every binary
#' configuration of \code{ame()} and \code{lame()} except \code{dynamic_rho}, next to the joint
#' scale move and the Hamiltonian updates described above where those run: in the Hamiltonian
#' configuration (changing coefficients and effects, rank one) the separating coefficient had 16 to
#' 20 effective draws of 2,000 after 8,000 sweeps (the two seeds agreed on its mean to within the
#' Monte Carlo error) and now has 1,000 to 1,200, and in the static directed configuration 140 to
#' 200 and now 630 to 910.
#' \code{prior$binary_scale_move = FALSE} omits it; \code{fit$binary_scale_move} records it.
#'
#' The actor-period effects of a \code{"poisson"} fit with \code{dynamic_ab} are moved together with the
#' latent values of their rows and columns. Given the latent array an effect has the precision of the
#' cells it touches over the residual variance, so a Poisson panel of sparse counts moved the effects,
#' and with them the innovation scale \code{sigma_ab}, the persistence and the levels, by small steps
#' (on 50 actors in five periods \code{sigma_ab} had 6 to 34 effective draws of 800 and the effects of
#' the first and last period, whose paths have one neighbour, R-hat up to 1.34). Each sweep moves an effect
#' and the latent values of its row by one amount, the cell residuals staying where they are, so that the
#' conditional is the prior of the effect given its neighbours and level times a Poisson factor that is
#' two sums per row; the draw is an elliptical slice step, and with an intercept the effects of a period
#' are moved in pairs, which keeps their sum. \code{sigma_ab} then has 216 to 311 effective draws
#' (\code{prior$ab_period_move = FALSE} omits the move; \code{fit$ab_period_move} records it).
#'
#' Under \code{dynamic_ab} (every family; bipartite, directed and symmetric
#' unipartite networks) the effect of an actor in a period is a level plus a
#' stationary AR(1) deviation, and the likelihood sees only their sum, so a
#' short or sparse panel splits the common variance of the effects between the
#' level variance and a persistent deviation almost freely, and the sweep,
#' which draws the levels given the hyperparameters and the hyperparameters
#' given the levels, crosses that ridge only over thousands of sweeps. Each
#' sweep therefore draws the level variances (for a directed network their
#' covariance, the pair of an actor's levels being correlated), the
#' innovation scale and the persistence from their density given the totals,
#' with the levels integrated out, before the levels are redrawn from their
#' conditional. The update uses only the totals, so it needs no latent array
#' and does not depend on the family, leaves the posterior unchanged and is
#' recorded in \code{fit$ab_hyper_move}; \code{prior$ab_hyper_move = FALSE}
#' omits it and \code{prior$ab_hyper_move_tries} sets the proposals per sweep
#' (default 8). For bipartite fits with static positions the scale of the
#' positions trades against the scale of the interaction matrix, which the
#' unit priors anchor weakly. Each sweep also draws, exactly from its
#' conditional, the factors by which the columns of the positions and the rows
#' and columns of the interaction matrix are rescaled without changing their
#' product (recorded in \code{fit$uv_g_scale_move};
#' \code{prior$uv_g_scale_move = FALSE} omits it).
#'
#' Directed \code{"cbin"} and \code{"frn"} fits receive a receiver-level
#' counterpart of the update their sender levels already have, which
#' translates a sender's level together with its whole latent row. The
#' data-augmentation sweep alone moves a receiver's level and its latent
#' column by about \eqn{1/\sqrt{m}} per sweep, \eqn{m} the number of cells in
#' the column, however wide the posterior is, and for a receiver nobody
#' nominates (a column bounded above only) the effective sample size of the
#' level falls roughly as \eqn{1/n}. One sweep of the update visits each
#' receiver once and translates its level and its latent column by the same
#' amount. Every residual is unchanged, so the amount enters the density only
#' through the prior of the shifted level and the order constraints of the
#' translated cells (a nominee stays above zero and above the non-nominees of
#' its row, a non-nominee stays below the weakest nominee and, while the row
#' has room, below zero, and the ranks of \code{"frn"} rows stay in order). It
#' is drawn exactly from the prior of the level cut to the interval those
#' constraints leave: there is no proposal and nothing to tune. With an
#' intercept the levels are centred and the translation is compensated by the
#' other levels and the intercept, with the prior of the intercept entering:
#' the coefficient prior of a static intercept or, in a \code{dynamic_beta} fit
#' of kind \code{"ar1"} or \code{"rw1"} on the unit period grid, the prior of
#' the intercept's path, which the translation moves by one constant (a static
#' intercept in a \code{dynamic_beta} fit enters through its own prior). Hidden
#' cells are labelled by their latent values, as in the sampler, and keep their
#' labels during the move. Panels with absent actors and bipartite fits are
#' covered, and so are \code{dynamic_ab} fits, whose sender and receiver
#' translations move the level together with its whole path (the deviations
#' from the level are untouched), and \code{dynamic_beta_per_actor} fits, whose
#' actor-specific slopes are a function of their own parameters and the design
#' alone and so stay where they are while the levels, the intercept and the
#' latent array move. The sender translation of a \code{dynamic_ab} fit, which
#' a static fit has always made every sweep, runs under the same switch and
#' schedule as the others. Symmetric \code{"cbin"} and \code{"frn"}
#' fits are not available. The update runs every sweep (the compiled sweep adds
#' roughly a tenth to the time of a sweep); \code{prior$actor_move_every}
#' changes that, \code{prior$actor_move = FALSE} omits the update, and
#' \code{fit$actor_move} records its work.
#'
#' Two more updates serve the same fits. A row at its cap is bounded below
#' only, but a row with room (fewer nominations than its cap) is pinned from
#' both sides, its non-nominations below zero and its nominations above, and a
#' sparse panel has a few of them (3 to 15 percent of the rows): translating
#' the whole latent array together with the intercept has room for a step of
#' the order of 0.01, so the intercept moved by about \eqn{1/\sqrt{m}} per
#' sweep however wide its posterior was (4 to 12 effective draws of 1,200 at 170
#' actors, and a level of the predictor that failed the convergence gate). The
#' level update translates the rows that are free above (the rows at their cap
#' in every period) by one amount and leaves the pinned rows where they are: the
#' intercept takes the compensation (with centred levels it and the pinned rows'
#' levels take it between them), every residual is unchanged, and the amount is
#' drawn exactly from the prior of the levels and the intercept cut at the weakest
#' nomination of the moving rows. It needs an intercept (a static coefficient,
#' or the path of a \code{dynamic_beta} fit of kind \code{"ar1"} or
#' \code{"rw1"} on the unit period grid) and the row effects, since the row
#' effects are what keeps the rows with room in place while the other rows rise
#' (with \code{rvar = FALSE} the receiver and reciprocity updates still run);
#' a bipartite fit of \code{ame()} (levels
#' free of centring) also draws, exactly, along the ridge the intercept shares
#' with the means of its levels. The dyadic correlation of a directed fit with
#' \code{dcor = TRUE} has the same trouble on a larger scale: the non-nominations
#' of rows at their cap only have to stay below their row's weakest nomination,
#' so their pair correlation is whatever \eqn{\rho} and the latent draws make
#' it, and given the latent array \eqn{\rho} moves by about
#' \eqn{1/\sqrt{m}} per sweep, \eqn{m} the number of dyads, against a posterior
#' a few tenths wide (4 to 8 effective draws of 1,200 at 170 actors). Every tenth
#' visit the reciprocity update integrates the non-nominations and the cells no
#' nomination can sit in out of the density of \eqn{\rho}, draws it by
#' Metropolis steps given the nominations, and redraws the cells it integrated
#' out from their conditionals given the nominations and the new value
#' (\code{prior$actor_move_rho = FALSE} omits it; \code{fit$actor_move} records its
#' proposals and accepts).
#'
#' In a \code{dynamic_beta} fit of kind \code{"ar1"} or \code{"rw1"} on the unit
#' period grid with a dynamic intercept, two more updates move the path of the
#' intercept period by period. A period in which almost every row is at its cap
#' (0 and 1 rows with room of 170 in two of four periods of the simulated
#' panels) has a level the data hardly inform, and the level update above, which
#' moves one amount for every period, is stopped by the rows with room of the
#' other periods: the intercepts of those periods had 4 to 7 effective draws of
#' 1,600 (two chains of 800), the least mobile rows with room 13 to 24, and the
#' innovation scale of the path, which follows its shape, 6 to 16; the two-seed
#' convergence gate failed on it (split R-hat 1.10 for \code{"cbin"}, 1.27 for
#' \code{"frn"}). The period update translates the periods that no row with room
#' ties to another one (the connected components of the graph of periods and
#' rows with room) by one amount each, the rows with room compensated by their
#' own levels and the intercepts of the other periods, and draws the amounts
#' exactly. The room update moves the intercept of a period that has rows with
#' room by itself. A row with room is a probit row (its nominations above zero,
#' its non-nominations below), so the Gaussian conditional of the intercept
#' given the latent cells of those rows is narrower than its posterior by the
#' information a latent cell loses when it is integrated out (0.13 to 0.37 per
#' non-nomination). Given a shift of the intercept of a period, the latent slice
#' of the period translates with it, the non-nominations of its rows with room
#' are integrated out of the density of the shift and redrawn from their
#' conditional afterwards, the dyadic correlation entering through the partner
#' cell; the shift is a Metropolis step on that density, so that the intercept
#' moves by its posterior width, not by the width of the Gaussian conditional.
#' A period with a row whose hidden cells leave its room open (fewer observed
#' nominations than the cap and an observed non-nomination) is left to the
#' exact translation. \code{prior$actor_move_period = FALSE} and
#' \code{prior$actor_move_room = FALSE} omit the two updates
#' (\code{prior$actor_move = FALSE} omits every nomination update);
#' \code{fit$actor_move} records their work.
#'
#' Directed binary panels with rank-two Gaussian AR(1) position paths,
#' dynamic additive effects, static coefficients, an intercept, both actor
#' margins and the default unit period grid also receive joint Hamiltonian
#' updates. These integrate both binary latent responses out and retain the
#' full reciprocal-pair likelihood. Factor persistence and stationary scale
#' move in the same Hamiltonian trajectory as their whitened paths, with
#' their original priors and transformation Jacobians retained. The metric
#' uses fixed prior references rather than either proposed endpoint.
#' Additional updates move persistence with
#' path scales, translate factors while integrating and redrawing actor levels,
#' and update dyadic correlation with the latent responses integrated out.
#' Hamiltonian step sizes adapt during burn-in only and survive checkpoints.
#' Set \code{prior$marginal_move = FALSE} to disable this group.
#' Other configurations retain their existing updates. These moves preserve
#' the posterior target.
#'
#' A directed unipartite \code{dynamic_uv} fit of the \code{"ar1"} or
#' \code{"t"} kind with a normal, Poisson or binary outcome starts its chain
#' from a dynamic AR(1) ALS fit of the same model, mapped through
#' \code{\link{als_start_vals}}, unless
#' \code{start_vals} is supplied or \code{prior$als_start = FALSE}. Bipartite
#' fits keep the cold start: a bipartite ALS solution parks its scale in the
#' interaction matrix. An ALS start places the initial path on the observed
#' temporal structure; the chain then samples the paths and hyperparameters
#' from there.
#' The ALS fit is deterministic, so the start does not depend on \code{seed};
#' it is the penalised point estimate (\code{bias_correct = "none"}), not the
#' empirical-Bayes fit. A start fit that is refused as unconverged (class
#' \code{"lame_als_convergence"}) contributes its last iterate; a start with
#' non-finite values, or positions beyond \code{prior$uv_max_abs}, is not used.
#' \code{fit$als_start} records whether it was used.
#'
#' @param Y a T length list of relational matrices, or a 3D array of dimensions
#' \code{[n_row, n_col, T]}, where T is the number of time periods. Named lists
#' may have changing actor composition; actors are aligned by row and column
#' names. A longitudinal \pkg{netify} object is also accepted and converted with
#' \code{netify::to_lame(lame = TRUE)}. When \code{family} or \code{mode} is
#' omitted, the value inferred by \pkg{netify} is used. See family below for
#' data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#'   With \code{method = "als"} and a normal, binary or Poisson dynamic fit by
#'   the empirical-Bayes estimator, \code{rvar}, \code{cvar}, \code{nvar} and
#'   \code{dcor} are honoured (a disabled margin is a point mass at zero, as in
#'   the sampler; see \code{dcor} for the correlation); the penalised routes
#'   name them in the warning about ignored arguments.
#' @param cvar logical: fit column random effects (asymmetric case)? See
#'   \code{rvar} for \code{method = "als"}.
#' @param dcor logical: fit a dyadic correlation (asymmetric case)? The
#'   correlation \eqn{\rho} of the two residuals of a dyad carries the prior
#'   \eqn{p(\rho) \propto 1 - \rho^2} on \eqn{(-1, 1)}, which vanishes at the
#'   boundaries. Information about it depends on joint presence and absence
#'   of ties and the fitted means. With fewer than one mutual dyad per
#'   period, or fewer than 30 in all, the data inform it mainly through the
#'   absence of mutual ties, and a message gives the count. Its posterior can
#'   then stay near the prior or move well away from it: toward -1 when the
#'   fitted means predict more mutual ties than the panel has (on the
#'   40-state MID panel, -0.58 with 95 percent interval -0.92 to -0.13). The
#'   additive effects, coefficients and latent positions all retain the same
#'   sampled correlation, including on panels with few mutual dyads.
#'   \code{dynamic_uv = TRUE} keeps the correlation: the per-period position
#'   kernel whitens each dyad's pair of residuals by the current \eqn{\rho},
#'   exactly as the static-position kernel does, so the reciprocity is carried
#'   by \eqn{\rho} instead of being absorbed into the coordinates. A
#'   period-varying \eqn{\rho} (\code{dynamic_rho}) needs static positions. A
#'   Poisson residual variance below 0.1 leaves \eqn{\rho} to its prior in
#'   the same way, and a message says so.
#'   An exactly symmetric \code{Y} fitted with \code{symmetric = FALSE} has
#'   equal reciprocal residuals, so the posterior of \eqn{\rho} concentrates
#'   near one (with \code{dcor = FALSE} the two cells of a dyad count as two
#'   observations); a message says so and points at \code{symmetric = TRUE},
#'   which fits the undirected model. A normal model with \code{dcor = TRUE}
#'   has no proper posterior once the reciprocal pairs of the panel outnumber
#'   what its mean parameters absorb by four (the bound \code{\link{ame}}
#'   applies, pooled over the periods: they share \eqn{\rho}, so the excess of
#'   every period adds, and a parameter they share absorbs its contrasts
#'   once): \eqn{\rho} runs to one, with static or dynamic positions or none,
#'   and the variance draw breaks down. Such a fit is an error of class
#'   \code{lame_symmetric_y}.
#'   With \code{method = "als"} a directed binary fit (the probit link of the
#'   empirical-Bayes default) models the correlation: the two cells of a dyad
#'   are one bivariate probit observation, the correlation is one number for
#'   all periods (\code{fit$VC["rho"]}; the record, with its standard deviation,
#'   in \code{fit$rho_estimation} and \code{fit$eb$pair}), and it is estimated
#'   with the other hyperparameters under the prior above. A directed Poisson
#'   fit takes the two counts of a dyad as one observation at the fitted state
#'   of independent cells: the correlation is the posterior mean, under the same
#'   prior, of the profile of the pair probability, and the information of the
#'   pairs enters the standard errors. A directed normal fit takes the two cells
#'   of a dyad as one bivariate normal observation, in closed form, and
#'   estimates the correlation under the same prior; a directed network whose
#'   reciprocal pairs are all exactly equal is an error of class
#'   \code{"lame_improper_normal_reciprocity"}, as in \code{\link{ame}} (fit it
#'   with \code{symmetric = TRUE} or \code{dcor = FALSE}).
#'   \code{dcor = FALSE} fits independent
#'   cells. Undirected and bipartite networks have no pair of cells. On the
#'   dynamic route the pair of cells of a binary fit is integrated over the
#'   uncertainty of its latent factors or moving positions, and the fit takes
#'   1.4 to 2.7 times as long as the independent-cell fit (1.4 with additive
#'   effects only, 1.9 to 2.7 with static factors or moving positions);
#'   \code{options(lame.als_threads)} splits the work over threads and changes
#'   nothing else.
#' @param nvar logical: fit nodal random effects (symmetric case)? See
#'   \code{rvar} for \code{method = "als"}.
#' @param R integer: dimension of the multiplicative effects (can be zero).
#' A binary fit whose multiplicative term reproduces the observed ties
#' exactly is separated: the size of the term is then set by its prior, the
#' fit says so, and the signs and the fitted probabilities are what the data
#' determine; a smaller \code{R} (down to zero for a panel whose ties barely
#' change over time) reproduces them without saturating.
#' With \code{method = "als"} and \code{dynamic_uv = FALSE}, the latent factors
#' of a normal, binary or Poisson dynamic fit are static and have the variational
#' posterior of the static estimator (\code{\link{ame_als}}): their posterior
#' means are reported as \code{U} and \code{V}, the scale of their prior in
#' \code{fit$eb$factor_prior}, and the standard errors are marginal over the
#' factors and carry the uncertainty of that scale
#' (\code{fit$eb$factor_scale$sd_log};
#' \code{options(lame.eb_scale_uncertainty = FALSE)} leaves it out).
#' @param R_row integer: for bipartite networks, dimension of row node multiplicative effects (defaults to R)
#' @param R_col integer: for bipartite networks, dimension of column node multiplicative effects (defaults to R)
#' @param mode character: either "unipartite" (default) for square networks or "bipartite" for rectangular networks
#' @param dynamic_uv logical: fit dynamic multiplicative effects (latent
#' factors). The transition model is selected by \code{dynamic_uv_kind}.
#' A directed fit with dynamic positions carries the dyadic correlation of
#' \code{dcor} through the position kernel (see \code{dcor}).
#' \code{fit$U} / \code{fit$V} hold the positions factored from each period's
#' posterior mean product (see \code{\link{latent_positions}}); the per-iteration
#' cubes are kept as \code{fit$U_draws} / \code{fit$V_draws} by the default
#' \code{posterior_opts} (or with \code{posterior_opts = list(save_UV_draws =
#' TRUE)}), and give the posterior standard deviations. Default \code{FALSE}.
#' With \code{method = "als"} and a normal, binary or Poisson outcome the positions are
#' estimated by the dynamic empirical-Bayes estimator (see
#' \code{bias_correct}): \code{fit$U} and \code{fit$V} are the posterior means
#' (n x R x T) with posterior standard deviations \code{fit$U_sd} and
#' \code{fit$V_sd}, \code{fit$Oarr} (alias \code{fit$UVPM}) is the
#' posterior-mean product, \code{fit$rho_uv} and \code{fit$sigma_uv} are the
#' persistence and innovation scale, and \code{fit$eb$uv_scale} and
#' \code{fit$eb$uv_search} are the stationary variance and the search that chose
#' it. Directed, bipartite (\code{R_row} and \code{R_col} may differ: the
#' interaction has rank \code{min(R_row, R_col)}) and undirected networks are
#' covered.
#' @param dynamic_ab logical: fit dynamic additive effects (sender/receiver effects) that
#' evolve over time. When TRUE, the row effects (a) and column effects (b)
#' become time-varying as a per-actor level plus an AR(1) deviation (see the
#' model description in the details).
#' This captures temporal heterogeneity in actors' baseline propensities to send and
#' receive ties, allowing for smooth changes in activity levels and popularity over time.
#' For example, an actor's tendency to form outgoing ties might gradually increase or
#' decrease across observation periods. The AR(1) specification ensures temporal smoothness
#' while allowing for actor-specific evolution patterns. Implementation uses conjugate
#' updates where possible and C++ for computational efficiency. Default FALSE.
#'   Each dynamic effect is a per-actor level plus a zero-mean AR(1)
#'   deviation: the level is pooled over every period (its variance is the
#'   reported \code{va} / \code{vb}) and the deviation carries the
#'   persistence \code{rho_ab} and innovation scale \code{sigma_ab}. The
#'   marginal variance of a per-period effect is therefore
#'   \code{va + sigma_ab^2 / (1 - rho_ab^2)}, and \code{rho_ab} describes the
#'   deviations, so on short panels it reads lower than the persistence of
#'   the raw effect path. Both parts are drawn from the likelihood: the
#'   levels given the deviations and the deviations given the levels, so
#'   persistent actor heterogeneity lives in the effects. The regression
#'   coefficients are drawn conditional on the current effects, and the
#'   effects are drawn conditional on summing to zero in every period when
#'   an intercept is included. Under that constraint an individual actor's
#'   marginal variance has the additional factor \code{1 - 1/n}, where
#'   \code{n} is the number of actors in that margin.
#'
#'   The reported intercept and effects are the identified quantities when
#'   the fit has a multiplicative term; see the details.
#'
#'   A directed binary panel carrying fewer than
#'   `r .lame_sparse_period_rate` events per actor-period (two periods or
#'   more) reports \code{va} and \code{vb} as not estimated, and fewer than
#'   `r .lame_sparse_period_rate_intercept` (`r .lame_sparse_period_rate_intercept_rank` with a
#'   multiplicative term) its intercept. On simulated panels of 30 actors and 5
#'   periods (static and dynamic effects, ranks 0 and 1, 64 to 252 panels a
#'   rate bin) the 95 percent intervals of \code{va} and \code{vb} covered
#'   0.31 and 0.28 of the panels below 0.6 events per actor-period, 0.73 and
#'   0.68 from 0.6 to 0.9, 0.91 and 0.81 from 0.9 to 1.2, and 0.91 to 1.00
#'   from 1.2; the intercept covered 0.57 to 0.87 of the panels below its rate
#'   and 0.92 or more from it. The persistence \code{rho_ab} (coverage 1.00 at every rate) and the
#'   innovation scale \code{sigma_ab} (0.84 to 1.00) are reported at every
#'   rate, and so is \code{cab} (0.97 to 1.00). The draws remain in
#'   \code{fit$VC} and \code{fit$BETA}; \code{fit$events_per_actor_period} and
#'   \code{fit$events_per_actor} record the event rates. Supplying an effects
#'   prior (\code{prior$Sab0}, \code{eta0} or \code{etaab}) reports \code{va} and
#'   \code{vb} under that explicit prior, which the summary names.
#'   Count families keep reporting

#'   their variance components. The stored paths \code{a_dynamic} and
#'   \code{b_dynamic} contain the full level-plus-deviation effects.
#'   Static latent positions are centered for reporting, with the removed
#'   product components carried by the reported intercept and effects.
#'   Dynamic latent coordinates retain the coordinates drawn under their
#'   proper path prior. Compare fitted predictors and per-period products
#'   across chains; raw coordinates can differ through latent rotations.
#' @param dynamic_G logical. Not supported: \code{TRUE} is refused with an
#' error on both the MCMC and the ALS route, and the argument is kept only so
#' that existing calls and saved call objects still parse. With moving
#' positions a per-period interaction matrix is not identified apart from
#' them, because \eqn{U_t G_t V_t'} and \eqn{(U_t G_t) V_t'} describe the same
#' products; with static positions it is identified but not implemented. The
#' interaction matrix exists only in bipartite fits. Use \code{dynamic_uv = TRUE} and
#' read the per-period product \code{fit$UVPM}. Default FALSE.
#' @param dynamic_rho logical. For directed unipartite normal models, allow the
#' residual dyadic reciprocity parameter \eqn{\rho_t} to vary by period instead
#' of using one pooled dyadic correlation across the whole panel. It requires
#' static latent positions (\code{dynamic_uv = FALSE}). MCMC fits store
#' per-draw paths in \code{fit$RHO} and the posterior mean path in
#' \code{fit$rho_path}. This option is not supported with \code{method = "als"};
#' use a supported MCMC configuration for period-specific reciprocity.
#' Dynamic ALS fits use one pooled correlation, \code{fit$VC["rho"]}, and
#' report period-specific residual diagnostics in \code{fit$residual_rho_path}.
#' The legacy \code{fit$rho_path} field remains a diagnostic alias, labelled
#' by \code{fit$rho_path_role = "residual_diagnostic"}; it is not a fitted
#' covariance path. \code{als_start_vals()} uses the pooled ALS correlation,
#' including for legacy ALS objects. Default \code{FALSE}.
#' @param dynamic_uv_kind Character string selecting the transition model for
#' dynamic multiplicative latent positions when \code{dynamic_uv = TRUE}.
#' \code{"ar1"} uses Gaussian AR(1) drift, \code{"snap"} uses a mixture of
#' AR(1) drift and discontinuous reset transitions, and \code{"t"} uses
#' heavy-tailed Student-t innovations represented by local transition scales.
#' On the MCMC path, \code{"snap"} and \code{"t"} are available for directed
#' unipartite fits; symmetric fits use \code{"ar1"} (the symmetric product's
#' unit-variance coordinates do not hold their scale under the reset and
#' scale-mixture transitions). Under \code{"snap"}, a reset lands on the
#' stationary distribution of the position paths unless
#' \code{prior$snap_kappa} fixes its scale, and \code{rho_uv} and
#' \code{sigma_uv} are updated jointly with the latent reset indicators.
#' At rank one, \code{summary()} and \code{print()} report \code{rho_uv},
#' \code{sigma_uv} and (for \code{"snap"}) \code{pi_snap}. At rank two or
#' more, simulation studies found the \code{"snap"} chains settle where extra
#' resets stand in for drift, so those three are not reported for it (the
#' draws stay on the fit). A \code{"t"} chain starts from a dynamic ALS fit
#' and then recovers \code{rho_uv} and \code{sigma_uv} at any rank; if that
#' start is unavailable they are not reported at rank two or more.
#' Both kinds receive conditional sign updates that preserve the fitted
#' products and help chains explore alternative trajectory orientations.
#' With \code{method = "als"}, \code{"snap"} is refused. The ALS
#' \code{"t"} kind supports directed and symmetric normal, binary and
#' Poisson panels; bipartite Student ALS is refused.
#' A normal, binary or Poisson fit with \code{bias_correct = NULL} is fitted by the
#' dynamic empirical-Bayes estimator for \code{"ar1"} (Gaussian) and \code{"t"}
#' transitions (directed and undirected networks): the weights of the scale
#' mixture have a variational posterior, and the degrees of freedom are
#' \code{prior$t_nu} (default 4).
#' ALS Student-t dynamic UV fits attach the final local transition-weight
#' matrices as \code{fit$lambda_u} and \code{fit$lambda_v} (n x T, one in the
#' first period).
#' @param dynamic_beta logical, character, integer, or logical-vector flag
#' selecting which regression coefficients evolve over time via independent
#' AR(1) processes. Default FALSE keeps every coefficient static (the
#' historical behaviour). Accepted forms:
#' \itemize{
#'   \item \code{FALSE} / \code{NULL}: no coefficient is dynamic.
#'   \item \code{TRUE}: every coefficient is dynamic.
#'   \item character vector of block shortcuts (\code{"intercept"},
#'         \code{"dyad"}, \code{"row"}, \code{"col"}) or specific coefficient
#'         names from \code{colnames(BETA)} -- those become dynamic.
#'   \item integer vector: 1-based column indices into \code{colnames(BETA)}.
#'   \item logical vector of length \code{p}: per-coefficient mask.
#' }
#' Each dynamic-coefficient block (one per distinct intercept / dyad / row /
#' col label that has at least one dynamic coefficient) gets its own AR(1)
#' parameters \eqn{\rho_\beta} (truncated-Normal prior, default mean 0.8,
#' bounded between 0 and 0.999) and \eqn{\sigma_\beta^2} (inverse-Gamma
#' prior, default shape 2 with a scale set from the response scale, see
#' \code{prior}). The joint posterior over the time path \eqn{\beta_t}
#' is drawn by Forward-Filter / Backward-Sample (FFBS) inside the MCMC
#' loop, conditional on the current static-block beta and on (a, b, U, V).
#' When the intercept or a nodal coefficient is dynamic, a sum-to-zero
#' contrast basis is applied to the additive-effect sampler to keep
#' \eqn{a_i + intercept_t} identified. Requires at least 2 time periods.
#' \code{dynamic_beta} composes with multiplicative latent factors
#' (\code{R} > 0) and additive random effects for every family, across
#' unipartite/bipartite and directed/symmetric panels and in combination
#' with \code{dynamic_ab} / \code{dynamic_uv}.
#' The state equation is a zero-mean AR(1) on the coefficient itself,
#' \eqn{\beta_t = \rho_\beta \beta_{t-1} + e_t}, with no level term, so
#' \eqn{\rho_\beta} is the persistence of the coefficient's level toward zero
#' rather than the persistence of its movement around a level. A coefficient
#' whose level is far from zero therefore has a persistence close to 1 whatever
#' its movement: for a path with mean \eqn{m}, deviation sd \eqn{s} and
#' deviation persistence \eqn{r} the value this state equation implies is
#' \eqn{(m^2 + r s^2) / (m^2 + s^2)}, which is 0.95 for a coefficient of 0.7
#' whose deviation has sd 0.3 and persistence 0.7. Read \eqn{\rho_\beta} that
#' way, and read the path \code{BETA[, k, ]} for how much the coefficient
#' moves.
#' With a single dynamic coefficient block and small \eqn{T}, only
#' \eqn{T - 1} transitions inform its AR(1) persistence, so the data inform
#' \eqn{\rho_\beta} weakly and its posterior stays near the prior; the path
#' \code{BETA[, k, ]} is read separately from the hyperparameters, and
#' \code{\link{dynamic_beta_prior_summary}} shows the paths the prior implies.
#' The hyperparameters of each block (\eqn{\rho_\beta} and \eqn{\sigma_\beta},
#' and the pooled scale and persistence of \code{dynamic_beta_pool}) are drawn
#' from their exact conditionals given the path at every sweep (a draw costs
#' well under a millisecond). Earlier versions drew them every 20th sweep; on normal and
#' binary panels that left the scale and the persistence with 4 to 20 times
#' fewer effective draws (panels of 5 to 170 actors; at 170 actors 357 of 1,600
#' draws of the scale on a normal panel and 304 on a binary one, where every
#' sweep gives 1,599 and 1,335).
#' An interweaving step also moves the hyperparameters together with the path,
#' before the path is redrawn given them. It proposes new values of
#' \eqn{(\sigma_\beta, \rho_\beta)} by a symmetric walk on
#' \eqn{(\log \sigma_\beta, \mathrm{logit}\,\rho_\beta)}, rebuilds the path
#' from the same standardised innovations (the state equation divided by the
#' innovation scale; \code{"rw2"}: the part of the path that its second
#' differences leave flat stays fixed; \code{"matern32"}: the path scales with
#' \eqn{\sigma_\beta}, the length scale being fixed) and accepts with the
#' Metropolis ratio of the likelihood of the rebuilt path and the
#' hyperprior (the pooled hierarchies included). It targets the funnel of
#' hierarchical models: where the data say little about the path (few actors,
#' many periods) and the prior of the scale is diffuse, the centred draws alone
#' move the scale slowly, because the scale given the path is pinned at the
#' spread of the path in hand. On a panel of six actors and 30 periods with
#' \code{prior$sigma_beta_shape = 0.5} and \code{prior$sigma_beta_scale = 1e-4} the
#' step raised the effective draws of the scale 2.8 times and those of the persistence 2.0 times over the draws at every
#' sweep; at the default prior it raised them 1.0 to 1.5 times on the panels
#' checked (5 to 30 actors, 12 to 60 periods), which about repays the time it
#' takes there (10 to 60 percent of a sweep). Where the data pin the path (a
#' dense panel) about nine proposals in ten are rejected and the step costs 3
#' percent of a sweep at 170 actors. \code{prior$beta_hyper_nc = FALSE} omits the
#' step, and \code{fit$beta_hyper_move} records its attempts and accepts.
#' With \code{method = "als"}, selected intercept, dyadic, row-node, and
#' column-node coefficients are supported on normal, binary, and Poisson panels,
#' including named panels where actors enter or exit. Dynamic ALS uses
#' period-specific \code{Xrow} or
#' \code{Xcol} values for selected dynamic node coefficients; static node
#' coefficients still use per-actor means. Dynamic intercepts are
#' unavailable when \code{family = "ordinal"}.
#' Default FALSE.
#' @param dynamic_beta_kind character: state-space prior on the dynamic
#'   coefficient block(s). \code{"ar1"} (default) gives mean-reverting
#'   AR(1) with a truncated-Normal prior on \eqn{\rho_\beta}. \code{"rw1"}
#'   gives a random walk (\eqn{\rho_\beta} pinned at 1, no truncation),
#'   appropriate when you expect permanent drift with no mean reversion --
#'   e.g. trade-gravity coefficients in a permanently changing world economy.
#'   The alias \code{"random_walk"} is also accepted for \code{"rw1"}.
#'   When the stationarity warning in \code{summary(fit)} fires under the
#'   default AR(1), refit with \code{dynamic_beta_kind = "rw1"} for cleaner
#'   inference. Decision tree: \strong{mean-reverting?} AR(1).
#'   \strong{Permanent / unit-root drift?} RW1.
#'   \strong{Smooth with curvature?} \code{"rw2"} (second-order random walk).
#'   \strong{Smooth with a known length-scale?} \code{"matern32"} (Matern 3/2).
#'   Note: \code{"rw2"} and \code{"matern32"} use an R-level joint Gaussian
#'   sampler and run substantially slower per iteration than the C++ FFBS
#'   used for \code{"ar1"} / \code{"rw1"} (roughly 2-4x on n = 100, T = 10).
#'   On \code{method = "als"} the coefficient paths are \code{"ar1"},
#'   \code{"rw1"} or \code{"rw2"}, including the empirical-Bayes estimator;
#'   \code{"matern32"} is refused there and needs \code{method = "mcmc"}.
#' @param family character: one of "normal","binary","ordinal","cbin","frn","poisson" - see
#' the details below
#' @param intercept logical: fit model with an intercept?
#' @param symmetric logical: is the sociomatrix symmetric? Not available for
#'   the \code{"cbin"} and \code{"frn"} families, whose per-actor ranked
#'   nominations are directed by construction. A symmetric fit needs a
#'   symmetric \code{Y} in every period and symmetric \code{Xdyad} slices, as
#'   in \code{\link{ame}}: the regression term of a symmetric fit is
#'   symmetric, so an asymmetric slice is an error (class
#'   \code{lame_symmetric_xdyad}) that names the slice and the first period;
#'   so is it with \code{method = "als"}.
#'   Replace it by \code{(X + t(X)) / 2}. A value missing in one triangle of a
#'   slice is not an asymmetry: it is filled from the other triangle.
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' families. A cap below a row's observed out-degree contradicts \code{Y} and
#' is an error.
#'   A cap that most rows' observed outdegrees reach erases the cross-row
#'   spread in outdegree that the sender effects rest on: a row at its cap
#'   shows only that its level is high enough to fill the cap. The sampler
#'   pools the prior scale of the row side with the column side once
#'   `r round(100 * .lame_odmax_sat_threshold)` percent of the observed rows
#'   (`r round(100 * .lame_odmax_sat_threshold_single)` percent in a single
#'   period) nominate the full \code{odmax}, counting a row with hidden cells
#'   by what a full cap thins to at its observed share;
#'   \code{fit$odmax_saturated_share} records the share and
#'   \code{fit$odmax_sat_threshold} the threshold. \code{va}, \code{vb} and
#'   \code{cab} are reported at every share: on simulated panels of 24 to 30
#'   actors their 95 percent intervals covered the generating values in 0.88
#'   to 1.00 of the panels at every share, with and without a fifth of the
#'   cells hidden and at generating sender sds of 0.25, 0.5 and 1 (at 0.25
#'   the posterior mean of \code{va} runs up to 2.3 times the generating
#'   variance while the interval still covers). \code{va} and the level of the rows follow the design. A design of two periods
#'   or more with a cap share of 0.5 or more and a posterior-mean \code{va} from
#'   0.5 to 1.5 reports the intercept and the row-covariate coefficients as not
#'   estimated: the rows below the cap carry too little of the spread of the
#'   senders and the posterior of \code{va} sits near its prior. In simulated
#'   panels (24 actors, 4 and 6 periods, one chain at the package defaults,
#'   1440 fits at generating sender sds of 0.5, 1 and 1.5) the intercept was
#'   covered in 0.89 of the 312 panels in that region (0.76 at sd 0.5 and 1.5,
#'   0.97 at sd 1) and in 0.96 of the 1128 outside it. \code{va} is reported in
#'   the region: its interval covered 0.91 of the 312 panels (0.79 at sd 0.5,
#'   1.00 at sd 1, 0.72 at sd 1.5) and 0.99 outside it. The posterior mean
#'   of \code{va} ran 0.18 and 0.23 above the generating variance at sd 0.5
#'   (shares 0.3 and 0.7), 0.08 and 0.00 above it at sd 1 and 0.02 and 0.18
#'   below it at sd 1.5. In a single period (\code{ame()} and
#'   one-period \code{lame()} fits) a row at its cap leaves the level of the
#'   rows unbounded above: with every row at the cap the intercept interval
#'   covered 0.64 of the panels, and with cells hidden and a share from 0.5 it
#'   covered 0.84, against 0.98 or more elsewhere. The intercept and any
#'   row-covariate coefficient are reported as not estimated there, with the
#'   draws kept in \code{fit$BETA}. On 200 simulated 170 by 100 bipartite panels
#'   (cap 10, shares 0.3 and 0.7) the intercept was covered in 1.00 of the
#'   panels, \code{vb} in 0.98 to 1.00 and \code{va} in 0.90 to 1.00, its
#'   posterior mean 0.09 to 0.15 below the generating variance; nothing is withheld there.
#'   When the effects prior is supplied (\code{prior$Sab0}, \code{eta0} or
#'   \code{etaab}), the fit names the prior in the summary
#'   (\code{fit$prior_supplied} records it); a prior on the effects does not
#'   fix the level of a row, so the intercept stays withheld in those designs, and
#'   \code{va} is reported under that prior.
#'   A cap the design itself imposes cannot be relaxed: read \code{APM} /
#'   \code{BPM} as relative sender and receiver levels there.
#' @param prior a list containing hyperparameters for the prior distributions.
#' Available options and their defaults:
#' \describe{
#'   \item{s20}{Fixed positive prior scale for the dyadic residual variance
#'         \eqn{v_e}. The prior has one degree of freedom. By default, normal
#'         fits use the variance of the observed responses, computed once before
#'         sampling; constant responses use their mean square, or 1 if all zero.
#'         Poisson fits use 1 on the latent log-rate scale. The resolved scale is
#'         reported as \code{fit$s2_prior_scale}. It is never recomputed from the
#'         current residuals. A supplied value overrides this default. The coefficient prior is controlled separately by \code{g}.}
#'   \item{uv_prior}{For \code{sampler = "tempered"}, \code{"hierarchical"}
#'   (default) or the explicitly different \code{"unit"} coordinate prior.
#'   This option is refused on the Gibbs route. See Details.}
#'   \item{Sab0}{Prior scale matrix for the additive-effects covariance. A 2x2
#'         matrix where Sab0\\[1,1\\] is the prior variance for row effects,
#'         Sab0\\[2,2\\] is the prior variance for column effects, and off-diagonals
#'         control correlation between row and column effects. For the
#'         \code{"normal"} and \code{"poisson"} families, and for unipartite
#'         \code{"binary"}, \code{"ordinal"}, \code{"cbin"} and \code{"frn"},
#'         this defaults to \code{diag(2)} scaled to the observed
#'         sender/receiver heterogeneity in \code{Y} (for the ordinal, cbin
#'         and frn families, of its probit-moment scores) rather than to a
#'         fixed \code{diag(2)}; see \code{\link{ame}} for the rationale. For
#'         \code{"cbin"} and \code{"frn"} panels of at least two periods the
#'         two sides get their own scale, the row side read with the rows at
#'         their nomination cap as censored observations. This keeps the
#'         receiver spread from setting the sender's prior scale when capped
#'         outdegrees contain little sender information. Where the rows
#'         at the cap reach the reporting threshold of \code{odmax} the pooled
#'         scale is kept. Bipartite
#'         fits of the probit-scale families use \code{diag(2)}. Pass
#'         \code{diag(2)} explicitly for a fixed unit-scale prior.}
#'   \item{eta0}{Prior degrees of freedom for the additive-effects covariance
#'         \eqn{\Sigma_{ab}} (default: round(4 \\* vdfmlt) for unipartite
#'         \code{"binary"}, \code{"ordinal"}, \code{"cbin"} and \code{"frn"}
#'         fits, where vdfmlt is a probit-moment variance multiplier estimated
#'         from \code{Y}, and round(4 + 3 \\* n/100) otherwise, where n is the
#'         number of actors). Higher values impose stronger shrinkage of the
#'         row/column effects toward the prior scale. The multiplicative-effects
#'         static unipartite prior is the one of \code{\link{ame}}: an
#'         inverse-Wishart covariance with the degrees of freedom and the
#'         data-scaled \code{Suv0} of \code{kappa0} and \code{Suv0} in
#'         \code{\link{ame}}, read off every period of the panel (the static
#'         factors are shared by the periods). \code{Suv0} and
#'         \code{kappa0} apply to \code{\link{ame}} of a unipartite network only; a bipartite network has
#'         no such covariance, and supplying either is an error of class \code{"lame_prior_unused"}.}
#'   \item{etaab}{Prior degrees of freedom for covariance of additive effects
#'         (default: 4 + 3 \\* n/100). Controls shrinkage of row/column random effects.}
#'   \item{ab_min_observed}{Observation count used to shrink sparse actor
#'         initial values (default: max(10, T)). Subsequent actor updates use
#'         the specified Gaussian prior without deterministic shrinkage.}
#'   \item{rho_uv_mean}{For dynamic_uv=TRUE: Prior mean for UV AR(1) persistence.
#'         Gaussian bipartite AR(1) MCMC fits default to 0.5; other routes
#'         default to 0.9. Values close to 1 indicate high persistence.}
#'   \item{rho_uv_sd}{For dynamic_uv=TRUE: Prior SD for UV AR(1) persistence.
#'         Gaussian bipartite AR(1) MCMC fits default to 0.2; other routes
#'         default to 0.1. Explicit settings override these defaults.}
#'   \item{sigma_uv_shape}{For dynamic_uv=TRUE: Shape parameter for inverse-gamma
#'         prior on UV innovation variance (default: 2). Symmetric fits do not
#'         draw the innovation variance: their coordinate paths have unit
#'         stationary variance, so the innovation SD is
#'         \eqn{\sqrt{1 - \rho_{uv}^2}} and the diagonal \code{L} of the
#'         product carries the scale.}
#'   \item{sigma_uv_scale}{For dynamic_uv=TRUE: Scale parameter for inverse-gamma
#'         prior on UV innovation variance (default: 1; for the normal family of
#'         a directed unipartite network the typical size of the additive effects,
#'         the scale of the factor prior above, so that the prior is the one of a
#'         response measured in other units. The positions of a bipartite or
#'         undirected network are unit free: the interaction matrix and the
#'         weights carry the unit).}
#'   \item{uv_max_abs}{For dynamic_uv=TRUE: absolute bound on each free
#'         latent coordinate. Each actor-vector conditional rejects proposals
#'         outside this joint prior support. Default 8 for binary, ordinal,
#'         cbin and frn families and 50 for normal and Poisson families (the
#'         bound of a unipartite normal fit is multiplied by the square root of
#'         \code{sd(Y)}: its coordinates carry the square root of the unit of the
#'         response).
#'         Set to \code{Inf} to remove the restriction. A finite bound changes
#'         the prior; no deterministic coordinate normalization is applied.
#'         The box is part of the prior on the scale of the positions, which
#'         the data do not fix in a bipartite fit (the interaction matrix
#'         carries it) and which \code{uv_var_cap} already limits: on the
#'         validation panels (50 actors, 5 periods; bipartite 40 by 30) the
#'         largest coordinate drawn under the default box was 7.8 (the box is
#'         8 on the probit scale) and 6.9 with the box removed, a proposal
#'         fell outside it in at most 1 of 10000 sweeps, and the coefficients,
#'         variance components, fitted values and posterior mean of the product
#'         \eqn{U G V'} with and without the box differed by about as much as
#'         two chains of the same fit do (at most 0.0014 root-mean-square for
#'         the product, and none for the cbin and Poisson panels, where no
#'         proposal came near the box). It matters when a fit separates the data:
#'         without it the product runs away.}
#'   \item{G_max}{For bipartite fits with a static interaction matrix: absolute
#'         bound on entries of \eqn{G}. A joint Gaussian conditional proposal
#'         outside the bound is rejected. Default 20 for binary, ordinal,
#'         cbin and frn families with dynamic_uv=TRUE and 100 otherwise
#'         (multiplied by \code{sd(Y)} for the normal family: the interaction
#'         matrix carries the unit of the response).
#'         Set to \code{Inf} to remove this prior restriction. On the
#'         validation panels the largest entry drawn was 2.7, a tenth of the
#'         default bound, and no proposal fell outside it.}
#'   \item{sigma_uv_max}{For dynamic_uv=TRUE: rail on the innovation SD of
#'         the latent coordinate paths. Default 10 for bipartite fits, where
#'         \eqn{G} absorbs the scale, and none for unipartite fits. On the
#'         validation panels the posterior of the innovation SD stayed below
#'         1.4, a seventh of the rail.}
#'   \item{uv_var_cap}{For dynamic_uv=TRUE: bound on the stationary variance
#'   of a latent coordinate path. This is joint prior support for
#'   \code{sigma_uv} and \code{rho_uv}; both updates respect it.
#'   A finite cap also constrains the residual-variance update for normal
#'   responses. The bound is a multiple of
#'   the observation variance for the normal family and absolute on the
#'   latent scale otherwise: the binary, ordinal, cbin and frn families fix
#'   that variance at one, and for the Poisson family the latent scale is
#'   the log rate, whose overdispersion variance does not enter.
#'   Default 4 for the binary,
#'   Poisson, ordinal, cbin and frn families and \code{Inf} for the normal
#'   family. In a bipartite fit the interaction matrix carries the scale of the
#'   paths and the data do not identify it, so the posterior of the stationary
#'   variance lies against this bound (within 10 percent of the default in 29
#'   percent of the sweeps of the validation panel, against none in the
#'   unipartite ones; without the bound its 1st to 99th percentiles run from 3
#'   to 76). Between a cap of 4, a cap of 16 and none the coefficients and the
#'   variance components of the validation fit moved by at most 1.2 percent,
#'   while \code{rho_uv}, \code{sigma_uv} and the norms of \code{U} and
#'   \code{V} follow the cap: read a bipartite fit through the product
#'   \eqn{U G V'}, not through its factors. For the binary, ordinal, cbin and frn families \code{uv_max_abs}
#'   also defaults to 8 and \code{G_max} to 20 rather than 50 and 100, so that
#'   once a fit separates the data the product \eqn{U G V'} is bounded on the
#'   latent scale instead of running away.}
#'   \item{uv_temper_move}{Logical, default \code{TRUE}. For symmetric Poisson fits with \code{dynamic_uv} of the
#'   \code{"ar1"} kind: every \code{uv_temper_every} sweeps (a positive whole number, default 5; a fit with the
#'   transition and the \code{uv_nc_move} step costs 1.4 times one without them) make a tempered transition (Neal 1996) of one randomly chosen
#'   eigen-direction (its weight and coordinate paths, the latent array carried with the predictor), which crosses
#'   between the signatures of the weights that the other draws cross only in thousands of sweeps. Recorded in
#'   \code{fit$uv_temper_move}.}
#'   \item{uv_nc_move}{Logical, default \code{TRUE}. For \code{dynamic_uv} of the \code{"ar1"} kind: in
#'   symmetric fits, draw the persistence \code{rho_uv} given the standardised innovations of the paths
#'   (the paths then follow the new persistence), and in bipartite fits draw the factor by which the paths, the
#'   innovation scale \code{sigma_uv} and the interaction matrix trade scale at no cost to the predictor. The
#'   conditional draws cross these directions slowly (30 to 70 effective draws of 400 for the
#'   persistence of a 50-actor symmetric panel). \code{FALSE} omits them.}
#'   \item{cov_ridge_move}{Logical, default \code{TRUE}. For directed fits with a row or column
#'         covariate (a design slice constant along the columns or the rows) and static
#'         coefficients: each sweep also draws, exactly, the amount by which the coefficient
#'         and the additive effects of that margin (the levels and, with \code{dynamic_ab},
#'         the deviation paths; the intercept takes the constant of the centring) trade
#'         against one another without changing any cell's predictor, up to a
#'         constant of each period that the latent array absorbs. The conditional draws cross
#'         this ridge slowly (the coefficient of a period-varying nodal covariate had
#'         38 to 99 effective draws of 800 on a 170-actor ordinal panel). \code{FALSE} omits it.}
#'   \item{sparse_position_move}{For directed \code{"cbin"},
#'         \code{"frn"} and \code{"ordinal"} fits with latent factors (\code{R > 0}, not \code{dynamic_uv}): each
#'         sweep also draws the position of every sender and receiver together with the
#'         latent values of its row or column, exactly. Moving \eqn{u_i} changes the
#'         predictor of row \eqn{i} by a different amount in each column; carrying the
#'         latent values along changes no residual, so the conditional of the position is
#'         its prior times the indicator that the shifted row is still a latent row the
#'         data allow, drawn by an elliptical slice step. Actors the data bound on one
#'         side only (no nominations, no nominators, only the lowest category) move across
#'         their wide posterior in a few sweeps instead of by random walk. Not applied to
#'         \code{"cbin"} / \code{"frn"} panels with hidden cells. Recorded in
#'         \code{fit$sparse_position_move}.
#'         The coordinate paths of \code{"cbin"} and \code{"frn"} fits with \code{dynamic_uv} of the \code{"ar1"}
#'         kind (bipartite or directed) move in the same way, each path as a whole, every second sweep unless
#'         \code{sparse_position_every} says otherwise.
#'         An ordinal fit that runs it visits the level translations of
#'         \code{actor_move} every sweep unless \code{actor_move_every} says otherwise.}
#'   \item{sparse_position_every}{Sweeps between visits of \code{sparse_position_move}
#'         (default 1: every sweep).}
#'   \item{uv_split_move}{For directed fits with static positions (\code{R > 0}, not
#'   \code{dynamic_uv}): draw the split of the scale between the two factors exactly
#'   (default \code{TRUE}). Multiplying \eqn{U} by \eqn{c} and \eqn{V} by \eqn{1/c}
#'   leaves the product, so the likelihood, unchanged; with the covariance of the rows
#'   integrated out the factors have a matrix-t prior, and \eqn{\log c} is drawn from
#'   its density along the orbit by a slice step. The conditional draws of the sweep
#'   moved the split by a few percent per sweep against a posterior sd of 0.2 to 0.5
#'   (24 to 54 effective draws of 1,600 in the static directed cells of the acceptance
#'   sweep, 4 on a 170-actor ordinal panel); the coordinates of the actors the data say
#'   little about, and their products, followed it. \code{FALSE} omits the draw.}
#'   \item{symmetric_node_tie}{For symmetric binary and ordinal fits with a nodal
#'   covariate: draw the row and the column coefficient of the covariate as one
#'   (default \code{TRUE}). The design carries a nodal covariate twice, and the
#'   coefficients of its two slices enter the predictor of a symmetric network
#'   only through their mean, so the data leave their difference to its prior.
#'   The exact symmetric draws of the binary and ordinal families return a
#'   difference of the width of that prior, and the draws of the additive effects,
#'   of the positions and of the dynamic coefficients that read the predictor
#'   carried it as noise: a fit with \code{dynamic_ab} and a time-varying nodal
#'   covariate let the node coefficient and the effect paths walk off together
#'   (the innovation sd of the paths ran to 20 and more against 0.3), and the
#'   persistence of dynamic positions was pulled down. Each draw of the
#'   coefficients is now followed by replacing the pair by its mean, which keeps
#'   the predictor symmetric. Fits without dynamic effects, positions or
#'   coefficients, and the normal and Poisson families, are not affected.
#'   \code{FALSE} restores the earlier draws.}
#'   \item{rho_ab_mean}{For dynamic_ab=TRUE: Prior mean for additive effects AR(1)
#'         parameter (default: 0.8). Controls temporal smoothness of sender/receiver effects.}
#'   \item{rho_ab_sd}{For dynamic_ab=TRUE: Prior SD for additive effects AR(1)
#'         parameter (default: 0.15).}
#'   \item{sigma_ab_shape}{For dynamic_ab=TRUE: Shape parameter for inverse-gamma
#'         prior on additive effects innovation variance (default: 2).}
#'   \item{sigma_ab_scale}{For dynamic_ab=TRUE: scale parameter of the
#'   inverse-gamma prior on the innovation variance of the effect deviations.
#'   It is set once before sampling to the prior mean of
#'   \code{Sab} (that is \code{etaab * Sab0 / (etaab - 3)}, the mean of the
#'   inverse Wishart the covariance is drawn from) times
#'   \code{1 - rho_ab_mean^2}, so that a priori the stationary variance of an
#'   actor's deviation path matches the variance of the levels and per-period
#'   deviations are still pulled toward the actor's level unless the data
#'   insist; it is reported as \code{fit$sigma_ab_prior_scale}. This is a
#'   prior-scale convention. Because
#'   it is read off the
#'   prior of the additive-effect covariance rather than its current
#'   draw, it still follows the scale of the outcome while leaving the
#'   \code{sigma_ab} and \code{Sab} draws as conditionals of one joint
#'   density. Set a value to override it.}
#'   \item{ab_var_cap}{For dynamic_ab=TRUE: optional bound on the stationary
#'   variance of an additive-effect path, in the same units as
#'   \code{uv_var_cap}; the innovation scale is then drawn so that
#'   \code{sigma_ab^2 / (1 - rho_ab^2)} stays at or below it. Default
#'   \code{Inf} (no bound). Set it, for example to 4, only when the additive
#'   effects themselves are running away; on a near-separable panel a bound
#'   here pushes the separation into the interaction matrix instead.}
#'   \item{actor_factor_slice}{For supported directed binary rank-two AR(1)
#'     factor panels with dynamic actor effects and static coefficients, use
#'     paired actor/factor elliptical slice updates of the observed likelihood.
#'     Default \code{TRUE}; \code{FALSE} disables this additional transition.
#'     Actor centering and factor bounds are preserved.}
#'   \item{bipartite_scale_move}{For Gaussian, binary, Poisson or ordinal bipartite panels
#'   (the last two with \code{actor_move}) with moving AR(1)
#'   factors and a static interaction matrix: jointly rescale both factor paths,
#'   their innovation scale and the interaction matrix (default \code{TRUE}).
#'   This preserves fitted values and the original posterior target, including
#'   all prior bounds. Set \code{FALSE} to retain the conditional updates alone.}
#'   \item{marginal_move}{Run the joint observed-data updates and associated
#'   persistence, scale and translation moves for supported binary
#'   configurations (default \code{TRUE}); see Details.}
#'   \item{binary_scale_move}{For \code{"binary"} fits (every mode, not with
#'   \code{dynamic_rho}): draw the factor by which the coefficients and the static
#'   additive effects are multiplied, exactly and with the latent array
#'   integrated out, every \code{binary_scale_every} sweeps (default
#'   \code{TRUE}); \code{FALSE} omits it. See Details.}
#'   \item{binary_scale_every}{For \code{"binary"} fits: number of sweeps between
#'   draws of that factor (default 5).}
#'   \item{level_move_every}{For binary unipartite and bipartite fits: number of
#'   sweeps between visits of the scale move (default 30).}
#'   \item{level_move_tries}{For binary unipartite and bipartite fits: proposals per visit
#'   of the scale move (default 8). Raising either this or the frequency buys
#'   a longer effective sample along the level / spread ridge at a cost per
#'   sweep; the visit pays for the factor draw and the latent redraw once, so
#'   more proposals per visit are cheaper than more visits.}
#'   \item{dev_move}{For binary unipartite \code{dynamic_ab} fits: walk the
#'   second direction of the scale move, the per-period deviations together
#'   with \code{sigma_ab} (default \code{TRUE}). See Details.}
#'   \item{dev_move_tries}{Proposals per visit along that direction
#'   (default 8).}
#'   \item{actor_move}{For \code{"poisson"} and \code{"ordinal"} fits (unipartite,
#'   directed or symmetric, or bipartite; static or \code{dynamic_ab} effects;
#'   static coefficients or, for the level translations, \code{dynamic_beta} of
#'   kind \code{"ar1"} or \code{"rw1"}): run the translations of the levels with
#'   their latent rows and columns, and the further updates listed in Details
#'   (default \code{TRUE}). For binary fits (unipartite: directed or symmetric, with
#'   any additive margin and static coefficients or \code{dynamic_beta} of kind
#'   \code{"ar1"} or \code{"rw1"}, incomplete panels included; or bipartite with
#'   an additive margin or static positions, with static coefficients or
#'   \code{dynamic_beta} of kind \code{"ar1"} or \code{"rw1"}): run the
#'   collapsed per-actor update (default \code{TRUE}). For \code{"cbin"} and
#'   \code{"frn"} fits (unipartite or bipartite, static or \code{dynamic_ab}
#'   levels, static coefficients or \code{dynamic_beta} of kind \code{"ar1"} or
#'   \code{"rw1"} on the unit period grid, with or without
#'   \code{dynamic_beta_per_actor}): run the receiver-level translation, with an
#'   intercept the level update, and for directed fits with \code{dcor = TRUE}
#'   the reciprocity update; the sender translation of a \code{dynamic_ab} fit
#'   runs under the same switch (default \code{TRUE}). With
#'   \code{dynamic_uv}, Poisson and binary fits also redraw the position paths
#'   (\code{actor_move_paths}). See Details.}
#'   \item{chain_extensions}{Number of times a chain that is still moving at
#'   the end of the run is continued before the call is refused (default 3;
#'   0 refuses at the first flag). See the convergence policy.}
#'   \item{actor_move_every}{Number of sweeps between visits of the collapsed
#'   per-actor update (default 10), of the translations of the Poisson (default
#'   2) and ordinal (default 10) families or of the receiver-level translation and the
#'   level, period and room updates of \code{"cbin"} and \code{"frn"} fits, and
#'   of the sender translation of their \code{dynamic_ab} fits (default 1; their
#'   reciprocity update runs on every tenth visit).}
#'   \item{actor_move_rho}{For directed fits with \code{dcor = TRUE}, also
#'   update the dyadic correlation with the latent responses integrated out at
#'   each visit of the collapsed per-actor update (default \code{TRUE}); for
#'   \code{"cbin"} and \code{"frn"} fits it switches the reciprocity update.}
#'   \item{actor_move_period}{For \code{"cbin"} and \code{"frn"} fits with a
#'   \code{dynamic_beta} intercept: also translate the intercept of the periods
#'   that the rows with room do not tie together, period by period (default
#'   \code{TRUE}; \code{FALSE} keeps the common translation of the level update
#'   only). See Details.}
#'   \item{actor_move_room}{For \code{"cbin"} and \code{"frn"} fits with a
#'   \code{dynamic_beta} intercept: also move the intercept of each period that
#'   has rows with room, with the non-nominations of those rows integrated out
#'   (default \code{TRUE}). See Details.}
#'   \item{beta_hyper_nc}{For \code{dynamic_beta} fits: run the interweaving
#'   step that moves the hyperparameters of the coefficient paths with the path
#'   (default \code{TRUE}); \code{FALSE} omits it and keeps the draws at every
#'   sweep. See the \code{dynamic_beta} argument.}
#'   \item{poisson_period_move}{For \code{"poisson"} fits with \code{dynamic_beta}:
#'   translate every period, and the whole path, of each dynamic coefficient with
#'   its latent slice at every sweep (default \code{TRUE}); \code{FALSE} omits
#'   them. See Details.}
#'   \item{poisson_static_move}{For \code{"poisson"} fits with \code{dynamic_beta}
#'   and static coefficients: translate each static coefficient with the whole
#'   latent array every \code{actor_move_every} sweeps (default \code{TRUE}).}
#'   \item{ordinal_scale_move}{For \code{"ordinal"} fits (unipartite and
#'   bipartite, static or dynamic coefficients and effects): draw, at every
#'   sweep, the scale of the whole state (the coefficients, the static additive
#'   effects and the latent array together, with the free cutpoints of
#'   \code{ordinal_cutpoints = "explicit"}) and, with the default data-induced
#'   cutpoints, the shift and the scale of the latent array of each period alone,
#'   each exactly along a transformation that keeps every rank (default
#'   \code{TRUE}); \code{FALSE} omits them. See Details.}
#'   \item{actor_hyper_nc}{For \code{dynamic_beta_per_actor} fits: run the
#'   interweaving step that moves the innovation variance and the persistence of
#'   the actor-specific paths with the paths (default \code{TRUE}); \code{FALSE}
#'   omits it and keeps the draws at every sweep.}
#'   \item{actor_move_paths}{For Poisson and binary fits with moving positions
#'   (\code{dynamic_uv = TRUE}): redraw the whole position path of each actor
#'   at each visit (default \code{TRUE}); \code{FALSE} omits it and leaves the
#'   other updates of \code{actor_move} as they are. See Details.}
#'   \item{ab_hyper_move}{For \code{dynamic_ab} fits with an additive margin
#'   (every family and mode): draw the level variances (their covariance for
#'   directed unipartite fits), the innovation scale and the persistence of
#'   the additive effects given their totals with the levels integrated out
#'   (default \code{TRUE}). See Details.}
#'   \item{ab_hyper_move_tries}{Proposals per sweep of that update (default 8).}
#'   \item{ab_period_move}{For \code{"poisson"} fits with \code{dynamic_ab}: move each actor-period
#'   effect together with the latent values of its row (column), in pairs when the effects of a period
#'   sum to zero, exactly and at every sweep (default \code{TRUE}). See Details.}
#'   \item{uv_g_scale_move}{For bipartite fits with static positions: draw the
#'   scale split between the positions and the interaction matrix exactly
#'   (default \code{TRUE}). See Details.}
#'   \item{als_start}{For directed unipartite \code{dynamic_uv} fits of the
#'   \code{"ar1"} or \code{"t"} kind (normal, Poisson or binary outcome):
#'   start the chain from a dynamic AR(1) ALS fit of the same model (default
#'   \code{TRUE}; see Details).
#'   \code{FALSE} keeps the cold start.}
#'   \item{snap_kappa}{For \code{dynamic_uv_kind = "snap"}: the sd of the
#'   position a reset lands on. Unset (the default), a reset is a fresh draw
#'   from the stationary distribution of the position paths,
#'   \eqn{\kappa^2 = \sigma_{uv}^2 / (1 - \rho_{uv}^2)}, recomputed every sweep,
#'   and the reset positions enter the \code{rho_uv} and \code{sigma_uv}
#'   draws as stationary draws. A positive number fixes the scale. A fixed
#'   scale well above the positions' spread makes a reset implausible for any
#'   ordinary jump; choose this scale in the units of the latent coordinates.}
#'   \item{rho_beta_mean}{For dynamic_beta: Prior mean for the per-block AR(1)
#'         parameter on time-varying regression coefficients
#'         (default: 0.8). Closer to 1 = smoother evolution.}
#'   \item{rho_beta_sd}{For dynamic_beta: Prior SD for the per-block AR(1)
#'         parameter (default: 0.15).}
#'   \item{rho_beta_lower}{For dynamic_beta: Lower truncation bound on the
#'         AR(1) parameter (default: 0). Pass a negative value to allow
#'         negative autoregression.}
#'   \item{rho_beta_upper}{For dynamic_beta: Upper truncation bound (default:
#'         0.999). Closer to 1 admits near-unit-root behaviour.}
#'   \item{sigma_beta_shape}{For dynamic_beta: shape parameter of the
#'         inverse-Gamma prior on the per-block innovation variance
#'         (default: 2). The conditional combines this proper prior with the
#'         coefficient path's quadratic form and effective dimension under
#'         its selected temporal law. With few transitions informing a
#'         coefficient block the posterior of the innovation variance stays
#'         near this prior.}
#'   \item{sigma_beta_scale}{For dynamic_beta: scale parameter of the
#'         inverse-Gamma prior on the per-block innovation variance. It is set
#'         once before sampling from the response scale: the innovation
#'         covariance of the coefficient path is
#'         \code{sigma_beta^2 Lambda} with \code{Lambda} the g-prior scale
#'         matrix, so \code{sigma_beta} is the per-period innovation sd of a
#'         covariate's contribution to the linear predictor in the units of
#'         the response. The default prior mean innovation variance is
#'         0.04 times the residual variance (the variance of the two-way
#'         centred outcome for the normal and poisson families, 1 for the
#'         probit and rank families, whose latent scale is fixed). It is
#'         reported as \code{fit$sigma_beta_prior_scale}. Set a value to
#'         override it.}
#'   \item{sigma_beta_init}{For dynamic_beta: Initial value of the per-block
#'         innovation standard deviation (default: 0.25; \code{0.25 * sd(Y)} for
#'         the normal family). Affects mixing, not
#'         the stationary distribution.}
#'   \item{beta0_mean}{For dynamic_beta: Mean of the Gaussian prior on the
#'         first observed coefficient state under AR(1) or RW1 (default: 0).}
#'   \item{beta0_var}{For dynamic_beta: Variance of the Gaussian prior on the
#'         first observed coefficient state under AR(1) or RW1 (default: 10;
#'         \code{10 * var(Y)} for the normal family, in the unit of the response).
#'         This prior is independent of the sampled evolution parameters.
#'         RW2 uses its proper difference-precision prior and Matern 3/2 its
#'         covariance kernel over all observed periods instead.}
#' }
#' Common usage: prior = list(Sab0 = diag(c(1, 1)), eta0 = 10) for stronger
#' shrinkage, or prior = list(rho_uv_mean = 0.95) for higher temporal persistence,
#' or prior = list(rho_beta_mean = 0.95, sigma_beta_scale = 0.01) for very smooth
#' time-varying coefficients with tight innovations.
#' @param g optional scalar or vector for g-prior on regression coefficients.
#' Default is \code{n^2}, where \code{n} is the number of actors (\code{nA} for a
#' bipartite fit). The prior of the normal family scales with the response, for
#' every unit of \code{Y}: a static unipartite fit has the coefficient precision
#' \code{sum_t(crossprod(X_t)) / (T g)}, whatever the residual variance, and takes
#' \code{n^2 * var(Y)}; the coefficients of a bipartite fit, and the static
#' coefficients next to dynamic ones, have the prior \code{N(0, g s2)}, whose
#' \code{s2} already carries the unit, and take \code{n^2}. The g-prior
#' controls the variance of regression coefficients: larger values allow for
#' larger coefficient values. Can be a vector of length p for parameter-specific
#' control.
#' @param seed random seed for the MCMC sampler (default 6886). The sampler
#'   is seeded internally with this value, so results are reproducible by
#'   default and an external \code{set.seed()} call has \emph{no} effect on
#'   the chain -- pass a different \code{seed} here to vary the draws (e.g.
#'   when running multiple chains). The caller's \code{.Random.seed} is
#'   restored on exit, so fitting never perturbs your RNG stream.
#'   Changing \code{nscan} or \code{burn} changes the whole chain, not only
#'   its length, because the sweep permutations are drawn up front; to
#'   lengthen a finished chain use \code{\link{lame_resume}}.
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param chains integer: the number of independently seeded chains (default 2,
#'   or \code{getOption("lame.chains")}). Chain \eqn{k} runs with seed
#'   \code{seed + 1000 * (k - 1)}. With \code{chains > 1} the chains are checked,
#'   compared and pooled (see the convergence policy); \code{chains = 1} runs one
#'   chain and the drift check only. ALS fits and \code{sampler = "tempered"} run one
#'   chain.
#' @param cores integer or \code{NULL}: workers for the chains. The chains run
#'   forked on Unix, at most one worker per chain; \code{NULL} uses
#'   \code{getOption("lame.cores")} if set and otherwise all cores but one.
#'   Windows runs the chains one after the other.
#' @param plot logical: plot results while running?
#' @param verbose logical: print progress while running? Default FALSE.
#' @param method character: \code{"mcmc"} (default, the Bayesian MCMC fit) or
#'   \code{"als"} (the fast, MCMC-free point-estimation path). With no dynamic
#'   arguments, \code{method = "als"} uses the pooled-static
#'   \code{\link{lame_als}} estimator. For normal, binary,
#'   and Poisson panels, \code{method = "als"} also supports a dynamic
#'   estimator (an empirical-Bayes fit for normal, binary and Poisson panels,
#'   or with \code{bias_correct = "none"} the penalized point estimate; see
#'   \code{bias_correct}) for \code{dynamic_ab}, selected
#'   \code{dynamic_beta}, smooth AR(1) \code{dynamic_uv} in directed,
#'   symmetric, and bipartite panels, and
#'   Student-t \code{dynamic_uv} in directed and symmetric panels (bipartite
#'   t requests are refused).
#'   \code{dynamic_uv_kind = "snap"} is refused on this route (use
#'   \code{method = "mcmc"}). For a normal, binary or Poisson panel the dynamic
#'   route is by default an empirical-Bayes fit (see \code{bias_correct}):
#'   \code{rho_ab}, \code{sigma_ab}, \code{rho_uv}, \code{sigma_uv},
#'   \code{rho_beta} and \code{sigma_beta} are then estimates under the
#'   hyperpriors of \code{lame()} (the variance components and persistences by
#'   EM, the stationary variance of moving positions by a profile of the
#'   marginal objective), not the prior means and penalty-implied scales of the
#'   penalised fit. The uncorrected penalised fit (\code{bias_correct = "none"}
#'   or \code{"jackknife"}) reports \code{rho_uv},
#'   \code{rho_ab} and \code{rho_beta} as the AR(1) persistences its smoothness
#'   penalties used, which are the prior means unless \code{prior} sets them,
#'   not as estimates, and \code{sigma_uv}, \code{sigma_ab} and
#'   \code{sigma_beta} as the innovation scales implied by the fitted paths,
#'   which the penalty shrinks below the posterior scale. Either can seed a
#'   warm-started MCMC run via \code{\link{als_start_vals}}. Under
#'   \code{dynamic_ab} the reported \code{VC} entries \code{va} and \code{vb}
#'   are the variances of the actor \emph{levels}, as in an MCMC fit: the
#'   empirical-Bayes fit estimates them by EM (\code{fit$eb$sab}), and the
#'   penalised fits take them from the
#'   cross-period covariance profile of the fitted paths, whose lag-\eqn{k}
#'   means follow \eqn{v_a + \tau^2 \rho^k}, so the AR(1) deviation and the
#'   per-period estimation noise both decay out of the reported level variance.
#'   The mean per-period variance of the paths -- level plus deviation -- is
#'   kept as \code{fit$vc_period}, and the uncorrected variance of the
#'   time-averaged paths as \code{fit$vc_raw}. For a penalised normal-family
#'   fit (\code{bias_correct = "none"}), \code{\link{ame_als_bootstrap}}
#'   redraws the actor paths (level plus AR(1) deviation) in each replicate,
#'   so its intervals for \code{va} and \code{vb} carry the uncertainty of
#'   that split.
#'   Node-covariate coefficients use the same orthogonal additive-effect
#'   decomposition as \code{\link{lame_als}}; dynamic node coefficients use
#'   period-specific node values when selected by \code{dynamic_beta}, while
#'   static node coefficients use per-actor means. Named changing-composition
#'   panels are aligned to the union actor set; actor-entry gaps break the
#'   smoothing penalties of the penalised fits, and the empirical-Bayes fit
#'   treats an absent actor-period as a free AR state. Posterior draws and
#'   rank/censored families
#'   are not handled by the ALS routes. Dynamic ALS prediction evaluates
#'   the fitted periods, optionally with replacement dyadic covariates;
#'   \code{predict(fit, h = ...)} forecasting is not implemented and a
#'   positive horizon is refused.
#' @param als_stability character (only used when \code{method = "als"} and a
#'   dynamic ALS path is selected): \code{"none"} (default), \code{"quick"}, or
#'   \code{"validation"}. The non-\code{"none"} presets rerun the dynamic ALS
#'   fit from jittered starts and attach \code{fit$stability}, which reports
#'   fitted-surface and coefficient-path differences across starts.
#' @param als_max_iter positive integer (only used when \code{method = "als"}):
#'   maximum iterations for the selected ALS solver. Default \code{1000}.
#'   Dynamic binary and count fits count IRLS updates, each with at most 40
#'   block-coordinate sweeps. Stalled Gaussian bipartite and directed Poisson
#'   fits with AR position paths and static coefficients also allow this many
#'   joint gradient-refinement iterations (with static interactions in the
#'   bipartite case). They retain that result only when the same objective
#'   decreases and the likelihood-score convergence criterion passes.
#'   The dynamic empirical-Bayes estimator needs up to a few hundred passes on
#'   a weakly identified panel; the Gaussian-transition fit of a 40 x 10 panel
#'   takes about 150 passes at the scale its search settles on.
#'   A dynamic fit (empirical-Bayes or penalised) that has not met its
#'   convergence criterion in this many iterations is an error of class
#'   \code{"lame_als_convergence"}; its field \code{fit} holds the last
#'   iterate. The sweep limit of a penalised normal fit is raised once, to four
#'   times this value, before the fit is refused. A static fit that does not
#'   converge is refused in the same way, without a \code{fit} field (see
#'   \code{max_iter} in \code{\link{ame_als}}).
#' @param als_tol optional positive scalar (only used when
#'   \code{method = "als"}): convergence tolerance for the ALS objective and
#'   fitted values. \code{NULL} keeps each ALS estimator's built-in default.
#'   The dynamic empirical-Bayes estimator stops when every block's scaled
#'   score is below \code{sqrt(als_tol)} and, for moving positions, the
#'   objective has settled to within \code{als_tol} deviance units of the fixed
#'   point; for a Poisson fit the estimated distance of the noise variance
#'   \code{ve} to the limit of its iteration must be below the threshold as
#'   well.
#' @param als_lowrank_method optional character (only used when
#'   \code{method = "als"}): inner solver for the low-rank factor block,
#'   \code{"mm"}, \code{"als"}, or \code{"hybrid"}. \code{NULL} (default)
#'   keeps each route's own default -- \code{"mm"} on the static
#'   \code{\link{lame_als}} route, \code{"hybrid"} on the dynamic route. The
#'   dynamic empirical-Bayes estimator (the default for normal, binary and
#'   Poisson panels) solves the factors jointly with the other effects by Gauss-Newton
#'   and warns that the argument does not apply; \code{bias_correct = "none"}
#'   uses the alternating low-rank solvers.
#' @param als_linear_solver optional character (only used when
#'   \code{method = "als"} with no dynamic arguments): linear solver for the
#'   static block updates, \code{"eigen"}, \code{"qr"}, or \code{"auto"}.
#'   \code{NULL} (default) keeps \code{\link{lame_als}}'s default
#'   (\code{"eigen"}). The dynamic ALS route uses its own path solvers and
#'   aborts if a value is given.
#' @param als_multistart optional character (only used when
#'   \code{method = "als"} with no dynamic arguments, and \code{R > 0}):
#'   \code{"none"}, \code{"cheap"}, or \code{"full"}, forwarded to
#'   \code{\link{lame_als}}. \code{NULL} (default) defers to
#'   \code{multistart}. Supplying both uses \code{als_multistart}.
#' @param multistart character (only used when \code{method = "als"} with no
#'   dynamic arguments, and \code{R > 0}): \code{"none"} (default),
#'   \code{"cheap"} (4 starts), or \code{"full"} (8 starts). Forwarded to
#'   \code{\link{lame_als}}, which reruns the block-coordinate fit from
#'   several starts and keeps the one with the lowest objective. The dynamic
#'   ALS route has no multi-start and ignores this argument with a warning;
#'   use \code{als_stability} there instead.
#' @param bias_correct character or \code{NULL} (only used when
#'   \code{method = "als"}): how the ALS point estimate handles the
#'   finite-sample (incidental-parameter) bias of a fixed-effects fit.
#'   For the static ALS point estimate it is forwarded to
#'   \code{\link{lame_als}}. \code{NULL} (default) fits a normal, binary or count
#'   static model by the empirical-Bayes estimator (\code{"eb"}: a Gaussian
#'   prior on the actor effects with its covariance estimated from the data,
#'   Laplace standard errors from \code{vcov()}, the exact Gaussian cell for a normal response and for counts the
#'   Poisson-lognormal model that \code{lame()} samples; see
#'   \code{bias_correct} in \code{\link{ame_als}}), for directed, bipartite and
#'   undirected networks, with or without a latent rank. \code{"eb"} asks for
#'   the empirical-Bayes estimator explicitly.
#'
#'   A normal, binary or Poisson fit with \code{dynamic_ab},
#'   \code{dynamic_beta} or \code{dynamic_uv} (of the \code{"ar1"} or
#'   \code{"t"} kind, with static or moving latent positions) and
#'   \code{bias_correct = NULL} is fitted by the dynamic empirical-Bayes
#'   estimator of the model \code{lame()} samples: the actor levels and
#'   zero-mean AR(1) deviation paths are random effects and the coefficients
#'   have a Gaussian prior; the variance of the cell-level noise (normal and
#'   Poisson), the covariance of the actor levels, the persistence and
#'   innovation scale of the additive paths (and of moving coefficients and
#'   positions) are estimated by EM under the \code{lame()} hyperpriors; the
#'   coefficients have Laplace standard errors (\code{vcov()},
#'   \code{confint()}, \code{fit$eb$se}; exact for a normal response); static
#'   factors are the variational posterior of \code{ame_als()}
#'   (\code{options(lame.eb_scale)} selects the rule that sets the scale of
#'   their prior, \code{"cavity"} by default) and moving positions have a
#'   mean-field Gaussian posterior whose product uncertainty is carried
#'   through the cell likelihood; and a directed panel gets the dyadic
#'   correlation (\code{dcor}).
#'
#'   A normal response has an exact Gaussian cell: the posterior of the
#'   coefficients, the levels and the paths at the hyperparameters is exact,
#'   the standard errors are the marginal ones given the hyperparameters, the
#'   dyadic correlation is an exact bivariate normal pair, and the reported
#'   intercept is the population intercept. A Poisson response is the
#'   Poisson-lognormal model: the correlation of a directed panel is taken at
#'   the fitted state, with the information of the pairs in the standard
#'   errors; the paths and levels are reported at their first-order posterior
#'   mean, the linear predictor (\code{EZ}) and the intercept are on the
#'   latent scale, and \code{fitted()} and \code{predict()} return
#'   \code{exp(EZ + ve / 2)}.
#'
#'   A dynamic fit that did not meet its convergence criterion, or whose
#'   moving positions have no identified scale, is not returned: the call
#'   stops with an error of class \code{"lame_als_convergence"},
#'   \code{"lame_dynamic_uv_scale_edge"} or
#'   \code{"lame_dynamic_uv_scale_unsettled"} whose field \code{fit} holds the
#'   last iterate.
#'
#'   \code{"none"} returns the fixed-effects fit of the earlier releases at any
#'   size (on the dynamic route the penalised fit with fixed smoothness
#'   penalties). The penalised dynamic fit descends at every sweep and is
#'   returned only if it met its first-order criterion; the static factors of
#'   a binary or Poisson fit carry the level prior of the dynamic factors
#'   (precision 1). \code{"jackknife"} requires the split-sample correction of a
#'   static fixed-effects fit. On the dynamic route it is available for
#'   directed Poisson panels with static coefficients, at least four actors,
#'   and rank zero or AR(1) positions, and corrects the penalised fit for
#'   latent count variance and reciprocal covariance using two complementary
#'   dyad histories. The correction combines natural variances and covariances
#'   as twice the full estimate minus the observation-fraction-weighted
#'   halves. It preserves response means by shifting the latent intercept by
#'   half the variance change. Unrestricted moments, boundary projections and
#'   every half-fit attempt are retained in \code{fit$bias_correct};
#'   \code{fit$raw} holds the original summaries. Each half uses
#'   \code{als_max_iter}, with one warm retry at five times that budget if
#'   necessary. All final halves must pass the same likelihood-score criterion.
#'   The correction removes the first-order term of the finite-sample bias of
#'   the latent variance and reciprocal covariance, the term that falls as the
#'   inverse of the number of observed cells per actor. Other
#'   dynamic specifications refuse an explicit \code{"jackknife"} request
#'   (class \code{"lame_jk_split"}). \code{"bootstrap"}, the
#'   parametric-bootstrap correction, is for static binary fits.
#'
#'   On the dynamic route the empirical-Bayes estimator
#'   estimates the persistences and innovation scales, so it refuses the
#'   smoothness penalties of the penalised fit (\code{prior$lambda_ab_als},
#'   \code{prior$lambda_beta_als}, \code{prior$lambda_uv_als}) with an error of
#'   class \code{"lame_dynamic_eb_unsupported"} (which also has class
#'   \code{"lame_als_eb"}), and a bootstrap (\code{bootstrap > 0}) with an
#'   error of class \code{"lame_dynamic_eb_unsupported"}; use
#'   \code{bias_correct = "none"} to fix the penalties. A bootstrap of a
#'   penalised dynamic fit exists for the normal family only.
#' @param bootstrap integer (only used when \code{method = "als"}): number of
#'   bootstrap replicates. \code{0} (default) skips the bootstrap; \code{N >= 30}
#'   runs \code{N} replicates and attaches the result so that \code{\link{confint}}
#'   returns bootstrap intervals (a smaller positive \code{N} is an error, on
#'   the static and the dynamic route alike, raised before anything is
#'   fitted). The dynamic route bootstraps the penalised normal-family fit
#'   (\code{bias_correct = "none"}) only. Each replicate is drawn as
#'   \code{\link{ame_als_bootstrap}} draws it (the additive paths redrawn, then
#'   the residual) and refit from the fitted trajectories; a replicate that has
#'   not converged within \code{als_max_iter} iterations is retried at four
#'   times that budget and otherwise dropped, never kept; more than 5 percent
#'   of the replicates dropped is an error of class
#'   \code{"lame_als_bootstrap_convergence"}, as is a bootstrap of a fit that
#'   has not converged. A binary static network fitted with
#'   \code{bias_correct = "none"} has no bootstrap (an error of class
#'   \code{"lame_als_bootstrap_bias"}), nor has a penalised binary or Poisson
#'   dynamic fit (class \code{"lame_als_bootstrap_unsupported"}).
#' @param bootstrap_type character (only used when \code{method = "als"}):
#'   \code{"parametric"} (default) or \code{"block"}. The block scheme needs
#'   at least ten resampling blocks
#'   (\code{ceiling(T / bootstrap_block_length)}); on a shorter panel, and on
#'   any dynamic ALS fit, it warns and uses the parametric scheme instead.
#'   The block scheme returns no standard error or interval for the parameters
#'   that resampling time slices cannot vary (see
#'   \code{\link{ame_als_bootstrap}}); the parametric scheme gives them.
#'   A dynamic ALS fit has a bootstrap only when it is a penalised normal fit
#'   (\code{bias_correct = "none"}). The default dynamic empirical-Bayes
#'   estimator of every family refuses it with an error of class
#'   \code{"lame_dynamic_eb_unsupported"}, and \code{vcov()}, \code{confint()}
#'   and \code{summary()} carry its Laplace standard errors; a penalised binary
#'   or Poisson fit refuses it too.
#' @param bootstrap_block_length integer: block length for the block bootstrap.
#' @param bootstrap_seed optional integer seed for the bootstrap.
#' @param print Deprecated. Use \code{verbose} instead.
#' @param gof logical: calculate goodness of fit statistics?
#' @param start_vals list of parameter starting values for the MCMC chain.
#'   Tempered fits accept their own \code{fit$start_vals} replica matrix as a
#'   warm start; use \code{\link{lame_resume}} for exact continuation.
#'   \code{dynamic_uv} fits of the \code{"ar1"} kind start from a dynamic ALS
#'   fit by default when this is \code{NULL} (see Details).
#'   Use \code{\link{als_start_vals}} to start MCMC from an ALS fit, or the
#'   \code{start_vals} element of a previous fit to continue from its final
#'   state. Recognised slots are \code{Z}, \code{beta}, \code{a}, \code{b},
#'   \code{U}, \code{V}, \code{rho}, \code{s2}, \code{Sab} and, for dynamic
#'   fits, the scalars \code{rho_uv}, \code{sigma_uv}, \code{rho_ab} and
#'   \code{sigma_ab}; bipartite fits also take \code{G}, use a slot only
#'   when its shape matches the panel and warn about the ones they skip.
#'   The \code{rho_uv} and \code{rho_ab} values that
#'   \code{\link{als_start_vals}} emits are the penalty settings of the ALS
#'   fit rather than estimates (for an empirical-Bayes fit, the default for
#'   normal, binary and Poisson panels, they are estimates); they only seed the chain.
#' @param periodic_save logical: indicating whether to periodically save MCMC results
#' @param out_file character vector indicating name and path in which file should be stored if periodic_save is selected. For example, on an Apple OS out_file="~/Desktop/ameFit.rda".
#' @param save_interval quantile interval indicating when to save during the post-burn-in period.
#' @param log_lik_method one of \code{"observed_exact"} (default),
#'   \code{"observed_ghk"}, or \code{"augmented"}. Observed marginal scores
#'   are supported for normal, binary, Poisson-lognormal and explicit-cutpoint
#'   ordinal fits. Poisson integration uses an adaptive log-integral with
#'   concavity bounds on omitted tails. Capped nominations, fixed ranks and data-induced ordinal
#'   models do not have the cellwise observed likelihood claimed by
#'   \code{"observed_exact"}; such requests fail explicitly. The existing
#'   \code{"observed_ghk"} rank estimators are approximate diagnostics,
#'   and \code{"augmented"} stores Gaussian-on-Z latent diagnostics;
#'   \code{loo()} and \code{waic()} refuse both with an error of class
#'   \code{"lame_loo_factorization_unsupported"}, and \code{read_log_lik()}
#'   returns the scores. Even correct cell
#'   marginal scores require independent residual cells, or unique unordered
#'   dyads in symmetric fits, before LOO or WAIC can use them. Directed
#'   correlated cell marginals require joint dyad factors, which are not
#'   supplied by this interface. Inspect \code{fit$log_lik_contract}. \code{"augmented"} uses the augmented-data Gaussian-on-\eqn{Z}
#'   contribution. Note the GHK cost grows exponentially with the
#'   per-row dimension \eqn{D} (the number of observed dyads in a
#'   sender row): the per-row Monte Carlo budget escalates to
#'   \code{min(2048, max(64, 8 * 2^D))} draws across a bank of 8
#'   randomized Halton sequences, once per stored draw, so
#'   \code{"observed_ghk"} on a unipartite \code{"frn"} fit with more
#'   than roughly 12 actors (where \eqn{D = n - 1}) is expensive; rows
#'   with \eqn{D > 15} fall back to a fast pairwise-normal
#'   approximation. Power users can fix the per-row budget exactly via
#'   \code{prior$ghk_n_mc} (a single positive integer), which both
#'   raises the 2048 ceiling and lowers the \code{8 * 2^D} accuracy
#'   floor; lower budgets increase the downward Jensen bias of the
#'   log-lik estimate.
#' @param keep_snap_draws controls storage of post-burn-in MCMC snap-shift
#'   indicators when \code{dynamic_uv = TRUE} and
#'   \code{dynamic_uv_kind = "snap"}. \code{"none"} keeps only the posterior
#'   mean snap probabilities. \code{"summary"} stores cell-level count, mean,
#'   and variance summaries. \code{"draws"} stores retained binary indicator
#'   arrays as \code{fit$snap_draws} and, for directed fits,
#'   \code{fit$snap_draws_v}, with dimensions draw x actor x time and first
#'   period set to \code{NA}. \code{"chunked"} currently stores the same arrays
#'   in memory and labels the storage mode for downstream tooling.
#' @param dynamic_beta_per_actor optional, one of \code{NULL} (default),
#'   \code{"row"}, or \code{"col"}. When non-\code{NULL}, fits
#'   actor-specific time-varying deviations for one dyadic covariate on the
#'   selected margin. The default \code{per_actor_identifiability = "center"}
#'   keeps the population coefficient in \code{coef(fit)} and enforces a
#'   per-period sum-to-zero constraint on the actor deviations.
#'   Their innovation variance defaults to an inverse-Gamma(2, 1) prior
#'   (scale \code{var(Y)} for the normal family: the unit of the response),
#'   configurable through \code{prior$sigma_actor_shape} and
#'   \code{prior$sigma_actor_scale}. The transformed
#'   persistence \eqn{(1 + \rho)/2} defaults to a Beta(2, 2) prior, with
#'   shapes \code{prior$rho_actor_shape1} and \code{prior$rho_actor_shape2}, restricted to
#'   \eqn{-0.99 < \rho < 0.99}, with a stationary initial state. Priors act on
#'   the centered subspace with one fewer dimension than the actor count.
#'   All four hyperparameters must be finite and positive. Choose the innovation
#'   scale for the outcome and covariate units.
#'   Stored fits record these choices in \code{per_actor_prior}. Draws are
#'   stored in \code{fit$THETA_ACTOR} when retained; summary mode stores
#'   \code{fit$theta_actor_mean} and \code{fit$theta_actor_sd}.
#'   The innovation variance and the persistence are drawn from their exact
#'   conditionals at every sweep from the first one, and an interweaving step
#'   moves them together with the paths (\code{prior$actor_hyper_nc}). Where each
#'   actor's data say little about its path, which is the usual case, the
#'   draws every 20th sweep of earlier versions left them with 18 to 75 effective draws in 2,000
#'   and the draws at every sweep with 220 to 900; the step adds a factor of
#'   1.1 to 2.4.
#' @param per_actor_covariate_idx positive integer; index into the
#'   dyadic covariate cube to slope on for the per-actor extension.
#'   Default \code{1L}.
#' @param per_actor_identifiability one of \code{"center"} (default),
#'   \code{"exact_center"} or \code{"drop_population"}. \code{"center"}
#'   preserves the population coefficient and constrains the per-actor
#'   deviations to sum to zero per period via pairwise contrast FFBS.
#'   \code{"exact_center"} draws each actor's path from its own conditional
#'   AR(1) posterior and conditions on all period-wise sum-to-zero constraints
#'   using the complete temporal covariance. Dyadic reciprocity couples actors;
#'   that case uses a joint Gaussian draw on the constrained subspace, at a
#'   higher cost per sweep. Actor-specific slopes require directed or bipartite
#'   data and unit-spaced periods.
#'   \code{"drop_population"} is reserved and currently falls back to
#'   \code{"center"} with a warning.
#' @param keep_per_actor one of \code{"auto"} (default; full draws
#'   when memory cost < 250 MB else streaming summary), \code{"draws"}
#'   (full \code{[n_iter, n_actors, T]} cube), \code{"summary"}
#'   (streaming posterior mean + sd only), or \code{"none"} (only
#'   hyperparameter chains).
#' @param dynamic_beta_pool one of \code{"none"} (default), \code{"rho"},
#'   \code{"sigma"}, or \code{"both"}. When two or more coefficient blocks are
#'   dynamic, non-\code{"none"} values update shared hyperparameters across
#'   blocks for the selected persistence and / or innovation-scale component.
#'   Persistence pooling uses a Gaussian hierarchy on the logit of
#'   \code{(rho - rho_beta_lower) / (rho_beta_upper - rho_beta_lower)},
#'   with a normal mean prior and inverse-Gamma variance prior. Its current
#'   mean and scale are retained in checkpoint state.
#'   With a single dynamic block the value is recorded but there is no other
#'   block to pool against.
#'   On \code{method = "als"} the empirical-Bayes estimator
#'   honours it: \code{"sigma"} shares the innovation scales across the dynamic
#'   blocks, \code{"rho"} pools the first-order persistences through the same
#'   logit-normal hierarchy as the sampler (a random-walk block has none to
#'   pool), \code{"both"} does both; the penalised routes name it in the
#'   warning about ignored arguments.
#' @param time_index optional numeric vector of length \eqn{T} giving
#'   observation times. Default \code{NULL} treats periods as equally
#'   spaced. Strictly-increasing values are required when supplied. Unequal
#'   gaps scale the AR(1) / RW1 transition variance and are also used by the
#'   RW2 and Matern 3/2 dynamic-beta precision builders.
#' @param period_exposure optional non-negative numeric vector of length
#'   \eqn{T} giving period-level exposure offsets. **Wired for Poisson:**
#'   when supplied with any value != 1, the Poisson observation
#'   likelihood becomes
#'   \eqn{Y_{ij,t} \sim \text{Poisson}(e_t \cdot \exp(Z_{ij,t}))},
#'   while the latent \eqn{Z} retains its existing semantic of
#'   "unexposed log-rate". For other families, non-trivial
#'   \code{period_exposure} is an error (rescale \code{Y} or use a
#'   covariate offset). \code{NULL} (default) or all-ones uses the
#'   unscaled likelihood path. Zero exposure is supported when all observed
#'   counts in that period are zero; its observed log-likelihood is zero.
#' @param max_seconds optional positive scalar; if the MCMC wall-clock
#'   time exceeds this many seconds, the chain terminates cleanly and
#'   \code{fit$terminated_early} is set to \code{TRUE}. The draw arrays hold
#'   the draws stored until then. When the time runs out before the first
#'   stored draw, \code{fit$no_draws_stored} is \code{TRUE}, every draw array
#'   (\code{BETA}, \code{VC} and the others) has no draws along its draw
#'   dimension, the posterior means are \code{NA}, \code{print()} says so, and
#'   \code{coef()}, \code{summary()} and the other methods that summarise
#'   draws stop with an error of class \code{lame_no_draws};
#'   \code{lame_resume(fit, nscan_more = N, max_seconds = Inf)} continues the
#'   chains.
#' @param checkpoint_path optional file path. When set, the chain
#'   periodically writes \code{BETA}, \code{VC}, the full sampler state,
#'   seed, completed iteration count, remaining burn-in and original call
#'   to this file. Use \code{\link{lame_resume}}
#'   or pass the same path back as \code{resume_from = path} to
#'   continue.
#' @param checkpoint_every positive integer; iterations between
#'   checkpoint writes. Default \code{100L}.
#' @param sampler MCMC transition method. Default \code{"gibbs"} retains the
#'   established sampler. \code{"tempered"} uses exact observed-probit replica
#'   exchange for symmetric binary networks with positive static rank, a static
#'   intercept and optional symmetric dyadic covariates, and shared dynamic
#'   additive effects on a common actor grid. Its default burn-in is 15,000
#'   sweeps and default post-burn length is 200,000 sweeps. See Details for its
#'   explicit latent-coordinate prior and supported options.
#' @param resume_from optional path to a checkpoint file produced by a
#'   previous \code{lame(..., checkpoint_path = path)} call. When
#'   non-\code{NULL}, \code{lame()} short-circuits all other input
#'   parsing and delegates to \code{\link{lame_resume}} with the
#'   user-supplied overrides forwarded. Pass \code{nscan = K} to
#'   request \code{K} additional post-burn sweeps on the continuation.
#'   Supply \code{checkpoint_path} and \code{max_seconds} again to
#'   checkpoint or limit the continuation. Unfinished burn-in is completed
#'   before any new posterior draws are stored.
#'   Default \code{NULL}.
#' @param freeze_call logical: if \code{TRUE}, store a snapshot of the
#'   evaluated \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol} on
#'   \code{fit$data_snapshot} so that a later \code{update(fit, ...)} refits
#'   against the same data even if the caller has mutated those objects in
#'   their workspace. Memory cost equals the size of the data; default
#'   \code{FALSE}.
#' @param save_log_lik one of \code{FALSE} (default), \code{TRUE}, or
#'   \code{"chunked"}. When \code{TRUE}, stores the per-iteration pointwise
#'   log-likelihood matrix on \code{fit$log_lik} (an
#'   \code{[n_stored, n_obs]} double matrix). When \code{"chunked"},
#'   streams the log-lik values to per-column-chunk binary files under
#'   \code{log_lik_path} so the in-memory cost during MCMC is just one
#'   chunk's width; the chunks are recovered later via
#'   \code{read_log_lik(fit)}. Required for \code{loo::loo(fit)} and
#'   \code{loo::waic(fit)}.
#' @param posterior_opts optional list of posterior draw-storage options,
#'   either a plain \code{list()} or one built with
#'   \code{\link{posterior_options}} (which spells the position option
#'   \code{save_UV}). Recognised names (unknown names trigger a warning):
#'   \code{save_UV}, \code{save_UV_draws}, \code{save_ab}, \code{thin_UV},
#'   \code{thin_ab}. \code{save_UV = TRUE} and \code{save_UV_draws = TRUE} are
#'   aliases on this path. With
#'   \code{dynamic_uv = TRUE} they store the per-iteration dynamic
#'   latent-position cubes at every stored (\code{odens}) iteration on
#'   \code{fit$U_draws} / \code{fit$V_draws}, each an
#'   \code{[actor, dim, period, draw]} array. On a static fit with a
#'   positive multiplicative rank they store the per-iteration
#'   latent-position matrices on \code{fit$U_samples} / \code{fit$V_samples}
#'   (\code{[actor, dim, draw]} arrays, as \code{\link{ame}} does). Bipartite
#'   fits with positive latent ranks additionally store the
#'   interaction-matrix draws on \code{fit$G_samples}
#'   (\code{[R_row, R_col, draw]}).
#'   \code{save_ab = TRUE} stores identified additive-effect draws in
#'   \code{a_samples} / \code{b_samples} (actor by draw) for static effects,
#'   or \code{a_draws} / \code{b_draws} (actor by period by draw) for
#'   dynamic effects, with actor-level draws in \code{a_level_samples} /
#'   \code{b_level_samples}. \code{thin_ab} thins the retained coefficient
#'   draws further; \code{ab_draw_index} records their matching row indices.
#'   Compact joint terminal states are always retained on
#'   \code{fit$forecast_state} for posterior forecasts, even when full
#'   actor and position trajectories are not requested.
#'   \code{thin_UV} remains accepted for compatibility but is not applied
#'   here (\code{odens} already thins the stored chain). Requesting latent draw
#'   storage with a zero multiplicative rank warns and stores nothing.
#'   Default \code{NULL}.
#' @param ordinal_cutpoints character: cutpoint convention for
#'   \code{family = "ordinal"}. \code{"data_induced"} (default) uses the
#'   data-induced cutpoints; \code{"explicit"} samples ordered cutpoints with the first fixed at zero and unit latent residual
#'   variance. Independent directed and symmetric models use the Cowles (1996)
#'   marginal Metropolis-Hastings update. Directed models with \code{dcor = TRUE}
#'   use cutpoint Gibbs updates conditional on the correlated latent observations.
#'   Explicit cutpoints are supported for unipartite models; predictions retain
#'   the original numeric category labels and use each saved cutpoint draw.
#'   Ignored for other families.
#' @param log_lik_path directory to write log-lik chunks to when
#'   \code{save_log_lik = "chunked"}. Default \code{NULL} creates a
#'   session-scoped subdirectory under \code{tempdir()}.
#' @param log_lik_chunk_size column-width of each on-disk chunk when
#'   \code{save_log_lik = "chunked"}. Larger chunks mean fewer files (and
#'   one \code{readBin} per chunk on read) but more memory during MCMC.
#'   Default \code{10000L}.
#' @param model.name optional string for model selection output
#' @return \item{auto_extend}{the record of the convergence policy: \code{max_extensions},
#' \code{extensions} (made), \code{discarded_sweeps}, \code{triggers} (a list with
#' the flagged quantities of each extension) and the thresholds used.}
#' @return \item{BETA}{posterior samples of regression coefficients. A
#'   2-dimensional matrix \code{[n_stored, p]} when all coefficients are
#'   static. When \code{dynamic_beta} flags any coefficient as time-varying,
#'   \code{BETA} is a 3-dimensional array \code{[n_stored, p, T]} whose third
#'   dimension is the time period. Static coefficients are still present with
#'   their values replicated across the third dimension. \code{coef(fit)}
#'   collapses this to a \code{[p, T]} posterior-mean matrix.
#'
#'   \strong{Migration from amen.} Under \code{amen::ame()} the
#'   \code{BETA} slot was always 2-D. Scripts that compute
#'   \code{apply(fit$BETA, 2, mean)} on an amen fit will silently
#'   aggregate across periods when run against a \code{lame()} fit
#'   with \code{dynamic_beta} active (the second margin is the
#'   coefficient index in both shapes; the new third margin is time).
#'   Use \code{length(dim(fit$BETA))} to detect the shape, or call
#'   \code{coef(fit)} which returns a \code{[p, T]} matrix in either
#'   case.}
#' \item{VC}{posterior samples of the variance parameters. On a unipartite
#'   directed binary panel of two periods or more carrying fewer than 1.2
#'   events per actor-period, the draws of \code{va} and \code{vb} are kept
#'   here but \code{summary()}, \code{print()} and \code{confint()} report
#'   them as not estimated, with the measured coverage. Under a supplied
#'   effects prior the variance draws are reported under it (see
#'   \code{prior_supplied})}
#' \item{events_per_actor_period}{for the families with an event rate
#'   (binary, cbin, frn, poisson) on a panel of at least two periods, the mean
#'   over periods of the events per row actor; the rate the reporting methods
#'   read}
#' \item{events_per_actor}{on the same panels, the events per row actor over
#'   the whole panel (every period's events over the actors present in any
#'   period), quoted with the rate}
#' \item{prior_supplied}{\code{TRUE} when \code{prior} carried
#'   \code{Sab0}, \code{eta0} or \code{etaab}: where the data alone leave
#'   the effect variances to the prior (the sparse rule, a saturated cap) the
#'   reporting methods then show them under the supplied prior; the
#'   intercept and row-covariate coefficients a saturated cap withholds stay
#'   withheld}
#' \item{odmax_saturated_share}{for the cbin and frn families, the share of
#'   observed rows whose observed outdegree sits at their \code{odmax}; the
#'   reporting methods read it alongside \code{events_per_actor_period} and
#'   stop reporting \code{va}, \code{cab}, the intercept and any
#'   row-covariate coefficient at
#'   `r round(100 * .lame_odmax_sat_threshold)` percent or above
#'   (`r round(100 * .lame_odmax_sat_threshold_single)` percent for a panel of
#'   one period; \code{odmax_sat_threshold} records the one applied)
#'   (\code{vb} stays reported; see \code{odmax})}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b}
#' \item{APM_sd, BPM_sd}{posterior standard deviation of each additive effect
#'   (static effects), used by \code{simulate()} to draw the effects}
#' \item{U}{multiplicative row effects u. For \code{dynamic_uv = TRUE}, a 3D
#'   array (n x R x T) whose periods factor the posterior mean product of each
#'   period at rank R and are aligned across periods; \code{NA} for an actor
#'   absent from a period (see \code{\link{latent_positions}})}
#' \item{V}{multiplicative column effects v (asymmetric case), with the same
#'   layout as \code{U}}
#' \item{U_sd, V_sd}{for \code{dynamic_uv = TRUE} with saved latent draws (the
#'   default), the posterior standard deviation of each reported coordinate,
#'   same shape as \code{U} and \code{V}; \code{latent_positions()} reports
#'   them as \code{posterior_sd}}
#' \item{UVPM}{on the union of actors; under changing composition the dynamic array has NA rows and columns for actors absent from a period, and a static fit gives one matrix for every pair. posterior mean of UV}
#' \item{ULUPM}{posterior mean of ULU (symmetric case)}
#' \item{L}{for symmetric fits, the diagonal of the product \eqn{U L U'}: the
#'   eigenvalues of \code{ULUPM} when the positions are static; when they are
#'   dynamic, one diagonal for every period that carries the scale and sign of
#'   each direction, so that \eqn{U_t L U_t'} is the best rank-R part of each
#'   period's \code{ULUPM} with those signs}
#'  \item{EZ}{posterior mean of the linear predictor Z, on the link scale
#' for every family (the log rate for Poisson, the probit scale for the
#' binary, ordinal, cbin and frn families); \code{YPM} is the response-scale
#' counterpart. For \code{mode = "bipartite"}, \code{EZ} is a list of
#' per-period matrices with the row and column actors as dimnames, in the
#' lexicographic sort of the input actor names (so \code{"r10"} sorts
#' before \code{"r2"} unless you zero-pad).}
#' \item{YPM}{Stored response-scale fitted value, one matrix per period: the
#' average of the outcomes replicated at each stored draw, except for Poisson
#' fits, where it is the mean count \eqn{e_t \exp(\eta + \sigma^2 / 2)} at the
#' posterior means of the linear predictor and of the residual variance; see the
#' Poisson model details.}
#' \item{coefficient_prior}{Recorded specification of the coefficient prior
#'   used by the fitting build, as displayed by \code{prior_summary()}.}
#' \item{static_uv_prior}{For static asymmetric unipartite factors, the
#'   inverse-Wishart covariance prior's degrees of freedom, scale and dimension,
#'   the data scale \eqn{s} of \code{Suv0 = s I} (\code{data_scale}) and the
#'   multiplier of \eqn{2R + 2} in the degrees of freedom
#'   (\code{df_multiplier}; see \code{kappa0} in \code{\link{ame}}).
#'   Older fit objects do not carry this record.}
#' \item{sampler_health}{A data frame separating numerical failures from
#'   valid prior-bound rejections for each conditional sampler step. The
#'   compatibility counts in \code{tryErrorChecks} and \code{mh_counters}
#'   contain both. A numerical failure, in a conditional step or in a joint
#'   Hamiltonian, actor-factor slice or bipartite scale move, aborts the fit
#'   with an error of class \code{lame_sampler_failure} (the counts are in the
#'   \code{health} field of the condition), so a returned fit has none. A
#'   proposal outside the prior's box (\code{prior$uv_max_abs},
#'   \code{prior$G_max}) is a Metropolis rejection of the model, counted here;
#'   \code{summary()} states it when it is on one percent or more of the
#'   sweeps. Joint Hamiltonian acceptance and rejection
#'   records are kept in their corresponding move fields.
#'   \code{sampler_health_iterations} gives the number of sweeps counted,
#'   including burn-in.}
#' \item{var_cap_hits}{number of sweeps in which a stationary-variance cap
#'   (\code{prior$uv_var_cap}, \code{prior$ab_var_cap}) removed most of a
#'   conditional draw}
#' \item{level_move}{for binary unipartite and bipartite fits, the attempts,
#'   accepts, acceptance rate and final step SD of the level direction of the
#'   joint scale move (absent when the move did not run)}
#' \item{actor_move}{for Poisson and ordinal fits, the work of the translations:
#'   \code{visits}, the level translations drawn (\code{level_calls}) with the
#'   likelihood evaluations of the Poisson ones (\code{level_evals}) and the
#'   steps that failed to terminate (\code{stuck}) or that were left to the
#'   sampler because the current state lay outside the interval the translation
#'   needs (\code{skipped}), the draws of the location of the latent array
#'   (\code{location_calls}, ordinal), the coefficient translations
#'   (\code{coef_calls}, \code{coef_evals}, \code{coef_stuck}; Poisson), the
#'   position steps (\code{uv_calls}, \code{uv_evals}, \code{uv_stuck};
#'   Poisson with static positions) and the attempts and accepts of the scale
#'   move of the dyadic variance and correlation (\code{scale_attempts},
#'   \code{scale_accepts}, and \code{rho_attempts}, \code{rho_accepts} of its
#'   correlation-only proposals; Poisson: the carried scale move, with
#'   \code{dynamic_rho} the one of the translations); for binary fits, the work of the collapsed per-actor
#'   update: \code{visits}, the slice steps of the level translations
#'   (\code{level_calls}) and of the static positions (\code{uv_calls}), the
#'   likelihood evaluations they used (\code{level_evals}, \code{uv_evals}), the
#'   steps that failed to terminate (\code{stuck}) and the attempts and accepts of
#'   the reciprocity update (\code{rho_attempts}, \code{rho_accepts}; zero for
#'   bipartite fits, which add the slice steps of the deviation paths of actor
#'   pairs under \code{dynamic_ab}, \code{dev_calls} and \code{dev_evals});
#'   for \code{"cbin"} and \code{"frn"} fits, the work of the receiver-level
#'   translation: \code{visits}, the receivers moved (\code{col_calls}) and the
#'   receivers left to the sampler because the current state lay outside the
#'   interval the translation needs (\code{col_skips}), the sweeps in which the
#'   level update moved a set of rows (\code{level_calls}) and the rows it
#'   moved (\code{level_rows}), the sweeps in which the period update
#'   translated a set of periods (\code{period_calls}) and the sets it moved
#'   (\code{period_comps}), the proposals and accepts of the collapsed
#'   translation of the periods that have rows with room (\code{room_attempts},
#'   \code{room_accepts}) and the moves refused because a redraw did not
#'   terminate (\code{room_stuck}), and the proposals and accepts of the
#'   reciprocity update (\code{rho_attempts}, \code{rho_accepts}); with moving
#'   positions, Poisson and binary fits add the path steps drawn
#'   (\code{dpm_calls}), the likelihood evaluations they used
#'   (\code{dpm_evals}) and those that failed to terminate (\code{dpm_stuck}),
#'   and symmetric fits the proposals and accepts of the scale ridge draw
#'   (\code{dpm_scale_attempts}, \code{dpm_scale_accepts});
#'   absent when the updates did not run}
#' \item{beta_hyper_move}{for \code{dynamic_beta} fits, the sweeps, proposals and
#'   accepts of the interweaving step of the path hyperparameters and its
#'   acceptance rate; absent when the step did not run}
#' \item{poisson_period_move}{for \code{"poisson"} fits, the visits, the translations
#'   of the periods (\code{period_calls}) and of the whole paths (\code{path_calls}) and
#'   of the static coefficients (\code{static_calls}) drawn, the evaluations their
#'   slice steps used, and the proposals and accepts of the carried scale move
#'   (\code{scale_attempts}, \code{scale_accepts}, and for the correlation alone
#'   \code{scale_rho_attempts}, \code{scale_rho_accepts}); absent when the moves did not run}
#' \item{binary_scale_move}{for \code{"binary"} fits, the draws of the scale of the linear
#' predictor (\code{visits}), the evaluations of its density they used (\code{evals}) and the
#' draws that moved the state (\code{moved}); absent when the draw did not run}
#' \item{ordinal_scale_move}{for \code{"ordinal"} fits, the sweeps the scale moves visited
#'   (\code{visits}) and the sweeps at which the scale of the whole state moved
#'   (\code{orbit_moved}); absent when the moves did not run}
#' \item{pa_hyper_move}{for \code{dynamic_beta_per_actor} fits, the sweeps, proposals and
#'   accepts of the interweaving step of the actor-specific slope hyperparameters and its
#'   acceptance rate; absent when the step did not run}
#' \item{uv_nc_move}{for \code{dynamic_uv} fits of the \code{"ar1"} kind, the calls and the draws that
#'   changed the value (\code{sym_calls}, \code{sym_moved}, \code{bip_calls}, \code{bip_moved}) of the non-centred
#'   persistence step (symmetric) and scale step (bipartite); absent when the move did not run}
#' \item{uv_temper_move}{for symmetric Poisson \code{dynamic_uv} fits, the tempered transitions made, accepted
#'   and those that changed the sign of the weight; absent when the move did not run}
#' \item{ab_period_move}{for \code{"poisson"} fits with \code{dynamic_ab}, the sweeps, the pairs (or single
#'   effects) moved and the slice evaluations of the translations of the actor-period effects; absent
#'   when the move did not run}
#' \item{ab_hyper_move}{for \code{dynamic_ab} fits, the proposals and accepts of
#'   the update of the additive-effect hyperparameters given the totals and the
#'   acceptance rate; absent when the update did not run}
#' \item{uv_g_scale_move}{for bipartite fits with static positions, the number of
#'   sweeps and of rescaled coordinates of the draw along the scale ridge of
#'   the positions and the interaction matrix; absent when the update did not run}
#' \item{marginal_move}{For supported binary fits, joint-update attempts,
#'   accepts, acceptance rate, final step size and numerical rejections.}
#' \item{bipartite_scale_move}{For supported Gaussian, binary, Poisson or ordinal bipartite fits, joint
#' scale attempts, accepts, acceptance rate, prior-bound rejections and
#' numerical rejections. The fixed proposal requires no adaptive tuning.}
#' \item{actor_factor_slice}{For supported binary fits, completed paired
#'   actor/factor rounds, block and likelihood evaluation counts, prior-bound
#'   rejections, and separate numerical failures.}
#' \item{dynamic_marginal_move}{For supported binary fits with moving factors,
#'   Hamiltonian attempts, accepts, acceptance rate, tuned step size,
#'   numerical rejections and prior-bound rejections.}
#' \item{nc_ar_move}{Attempts and accepts for the joint persistence and
#'   path-scale updates of the factors and actor effects.}
#' \item{translation_move}{Number of visits to the factor/actor translation
#'   update, counting each margin separately.}
#' \item{rho_collapsed_move}{Attempts, accepts and acceptance rate for dyadic
#'   correlation updates with binary latent responses integrated out.}
#' \item{uv_scale_move}{Attempts, accepts and acceptance rate for joint
#'   latent-factor and innovation-scale updates.}
#' \item{dev_move}{the same for the deviation direction, present for binary
#'   unipartite \code{dynamic_ab} fits}
#' \item{als_start}{for \code{dynamic_uv} fits, whether the chain started from
#'   a dynamic ALS fit (\code{FALSE} when start values were supplied, when
#'   \code{prior$als_start = FALSE}, for the snap kind, symmetric and
#'   bipartite fits and the rank families, or when the ALS fit failed or fell
#'   outside the range the sampler holds)}
#' \item{sigma_ab_prior_scale}{for \code{dynamic_ab} fits, the inverse-gamma
#'   scale of the deviation innovation variance, fixed before sampling}
#' \item{sigma_beta_prior_scale}{for \code{dynamic_beta} fits, the
#'   inverse-gamma scale of the coefficient innovation variance, fixed before
#'   sampling from the response scale}
#' \item{G}{for bipartite fits, the interaction matrix of \eqn{U G V'}: the
#'   posterior mean of the gauged draws when the positions are static; when they
#'   are dynamic, a diagonal matrix fixed across periods that carries the scale,
#'   so that \eqn{U_t G V_t'} is the best rank-R part of each period's
#'   \code{UVPM}}
#' \item{GOF}{a named list with one element per goodness-of-fit statistic
#' (\code{sd.rowmean}, \code{sd.colmean}, \code{dyad.dep}, \code{cycle.dep},
#' \code{trans.dep}; for \code{mode = "bipartite"}, \code{sd.rowmean},
#' \code{sd.colmean}, \code{four.cycles}, with \code{four.cycle.dep}
#' replacing the count for Gaussian outcomes). Weighted rectangle dependence
#' is \code{NA} when no complete rectangle is observed or the values have
#' zero variance; see \code{\link{gof_stats_bipartite}}. Each element is a
#' \code{[T, 1 + n_stored]} matrix: one row per period, the first column the
#' observed value and the rest the posterior predictive draws. For the
#' data-induced ordinal family the predictive draws are quantile-matched to
#' observed category frequencies and assess dependence conditional on those
#' frequencies. Explicit ordinal fits instead use each draw's sampled
#' cutpoints, allowing checks of category frequencies as well.
#' See \code{\link{gof}} for post-hoc computation and \code{\link{gof_plot}} for visualization.}
#' \item{a_level, b_level}{for \code{dynamic_ab = TRUE}, the posterior means of
#' the per-actor levels behind the effect paths; \code{predict()} and
#' \code{simulate()} propagate the deviation around them}
#' \item{a_dynamic, b_dynamic}{for \code{dynamic_ab = TRUE}, the posterior mean
#' effect paths (actors by periods), with \code{a_dynamic_sd} /
#' \code{b_dynamic_sd} their posterior sds}
#' \item{rho_ab, sigma_ab, rho_uv, sigma_uv}{draws of the AR(1) persistence and
#' innovation scale of the dynamic effects and positions, when fitted}
#' \item{start_vals}{Final parameter values from MCMC, can be used as the input
#' for a future model run.}
#' \item{model.name}{Name of the model (if provided)}
#' @seealso \code{\link{ame}} for cross-sectional models,
#'   \code{\link{lame_als}} for the fast MCMC-free point estimator,
#'   \code{\link{als_dynamic_beta}} for a regression-only penalised
#'   smoother on the time-varying coefficient path,
#'   \code{\link{lame_resume}} for resuming saved checkpoints,
#'   \code{\link{gof}} for post-hoc goodness-of-fit computation,
#'   \code{\link{gof_plot}} for visualizing GOF results,
#'   \code{\link{latent_positions}} for extracting latent positions as a tidy data frame,
#'   \code{\link{procrustes_align}} for Procrustes alignment of latent positions,
#'   \code{\link{summary.lame}} for model summaries,
#'   \code{\link{coef.lame}} for coefficient extraction.
#' Long-format edgelists and covariates can be prepared with
#' \code{netify::netify()} and passed directly to \code{lame()}, or converted
#' explicitly with \code{netify::to_lame(lame = TRUE)}.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#'
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit<-lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X,
#'           burn=5, nscan=5, odens=1, family="binary", chains = 1)
#' # you should run the Markov chain much longer than this
#'
#' \donttest{
#' ## Time-varying regression coefficients (dynamic_beta).
#' ## Make every dyadic coefficient evolve as an AR(1):
#' fit_dyn <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X,
#'                 family = "binary", R = 0,
#'                 nscan = 60, burn = 15, odens = 5,
#'                 dynamic_beta = "dyad", chains = 1)
#' dim(fit_dyn$BETA)        # [n_stored, p, T] -- 3-D when dynamic
#' coef(fit_dyn)            # [p, T] posterior-mean coefficient paths
#' confint(fit_dyn)         # per-period 95% credible intervals
#' summary(fit_dyn)         # prints a "Dynamic coefficients per period" block
#' }
#'
#' @export lame
lame <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	rvar = TRUE,
	cvar = TRUE,
	dcor = !symmetric,
	nvar = TRUE,
	R = 0,
	R_row = NULL,
	R_col = NULL, # separate ranks for bipartite
	mode = c("unipartite", "bipartite"), # network mode
	dynamic_uv = FALSE,
	dynamic_ab = FALSE,
	dynamic_G = FALSE,
	dynamic_beta = FALSE,
	dynamic_rho = FALSE,
	dynamic_beta_kind = c("ar1", "rw1", "rw2", "matern32"),
	dynamic_uv_kind = c("ar1", "snap", "t"),
	family = "normal",
	intercept = !(family == "ordinal"),
	symmetric = FALSE,
	odmax = NULL,
	prior = list(),
	g = NA,
	seed = 6886,
	nscan = 10000,
	burn = 500,
	odens = 25,
	plot = FALSE,
	verbose = FALSE,
	gof = TRUE,
	start_vals = NULL,
	periodic_save = FALSE,
	out_file = NULL,
	save_interval = 0.25,
	model.name = NULL,
	save_log_lik = FALSE,
	posterior_opts = NULL,
	log_lik_path = NULL,
	log_lik_chunk_size = 10000L,
	keep_snap_draws = c("none", "summary", "draws", "chunked"),
	freeze_call = FALSE,
	dynamic_beta_pool = c("none", "rho", "sigma", "both"),
	dynamic_beta_per_actor = NULL,
	per_actor_covariate_idx = 1L,
	per_actor_identifiability = c("center", "exact_center", "drop_population"),
	keep_per_actor = c("auto", "draws", "summary", "none"),
	time_index = NULL,
	period_exposure = NULL,
	max_seconds = Inf,
	checkpoint_path = NULL,
	checkpoint_every = 100L,
	log_lik_method = c("observed_exact", "observed_ghk", "augmented"),
	ordinal_cutpoints = c("data_induced", "explicit"),
	method = c("mcmc", "als"),
	als_stability = c("none", "quick", "validation"),
	als_max_iter = 1000L,
	als_tol = NULL,
	als_lowrank_method = NULL,
	als_linear_solver = NULL,
	als_multistart = NULL,
	multistart = c("none", "cheap", "full"),
	bias_correct = NULL,
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	chains = .lame_chains_default(),
	cores = NULL,
	resume_from = NULL,
	print,
	sampler = c("gibbs", "tempered")
) {
	`%||%` <- function(x, y) if (is.null(x)) y else x
	# catch abbreviations of model.name: `model = "bin"` would partial-match it and silently fit
	# the default family. exact matches of other formals win over partial matching
	sc_nms <- names(sys.call()) %||% character()
	part <- sc_nms[nzchar(sc_nms) & sc_nms != "model.name" &
			startsWith("model.name", sc_nms) &
			!(sc_nms %in% names(formals(sys.function())))]
	if (length(part) > 0L) {
		cli::cli_abort(c(
			"{.arg {part[1]}} partially matched the {.arg model.name} argument.",
			"i" = "If you meant to choose the distribution, use {.arg family} (e.g. {.code family = \"binary\"}).",
			"i" = "{.arg model.name} is only a label; spell it out in full if intended."
		))
	}
	# resume from a checkpoint through lame_resume(), forwarding user overrides; the caller's frame
	# goes along as `.envir` so the saved call's references resolve when y is not supplied here
	if (!is.null(resume_from)) {
		caller_env <- parent.frame()
		mc0 <- match.call()
		# user-supplied named args except resume_from, skipping any still-missing promises
		user_named <- setdiff(names(as.list(mc0))[-1L], "resume_from")
		usr <- list()
		for (nm in user_named) {
			val <- tryCatch(
				get(nm, envir = environment(), inherits = FALSE),
				error = function(e) NULL
			)
			if (!is.null(val)) usr[[nm]] <- val
		}
		nscan_more <- usr$nscan
		usr$nscan <- NULL
		return(do.call(
			lame_resume,
			c(
				list(path = resume_from, nscan_more = nscan_more, .envir = caller_env),
				usr
			)
		))
	}
	sampler <- match.arg(sampler)
	if (identical(sampler, "tempered")) {
		if (missing(burn)) burn <- 15000L
		if (missing(nscan)) nscan <- 200000L
	} else if (is.list(prior) && !is.null(prior[["uv_prior"]])) {
		cli::cli_abort("{.code prior$uv_prior} is supported only with {.code sampler = \"tempered\"}.")
	}
	# a time-varying interaction matrix is refused on every route (see .lame_refuse_dynamic_G())
	.lame_refuse_dynamic_G(dynamic_G)
	family_missing <- missing(family)
	mode_missing <- missing(mode)
	symmetric_missing <- missing(symmetric)
	netify_input <- .lame_prepare_netify_input(
		Y = Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		mode = mode,
		family = family,
		mode_missing = mode_missing,
		family_missing = family_missing,
		symmetric = symmetric,
		symmetric_missing = symmetric_missing,
		longitudinal = TRUE,
		caller = "lame",
		fit_method = if (identical(method[1L], "als")) "als" else "gibbs"
	)
	if (isTRUE(netify_input$converted)) {
		Y <- netify_input$Y
		Xdyad <- netify_input$Xdyad
		Xrow <- netify_input$Xrow
		Xcol <- netify_input$Xcol
		mode <- netify_input$mode
		family <- netify_input$family
		symmetric <- netify_input$symmetric
	}
	# deprecated print argument
	if (!missing(print)) {
		# warn only when verbose was not also passed; multi-chain wrappers pass both
		if (missing(verbose)) {
			cli::cli_warn(c(
				"The {.arg print} argument is deprecated.",
				"i" = "Use {.arg verbose} instead."
			))
			verbose <- print
		}
	}
	# capture the call so summary.lame can echo it
	mc <- match.call()
	# the chain policy: a continuation block of an extension, or the user's call
	.chain_ctx <- .lame_chain_enter(environment(), mc, parent.frame())
	# the call as evaluated (the ALS start re-evaluates it); a chain's block shows the user's call
	.mc_eval <- mc
	if (isTRUE(.chain_ctx$block)) mc <- .chain_ctx$mc
	# several chains: each is a call of this function with chains = 1; tempered and ALS fits have none
	if (!isTRUE(.chain_ctx$block)) {
		if (missing(chains) && identical(sampler[1L], "tempered")) chains <- 1L
		chains <- .lame_chains_check(chains)
		if (chains > 1L && identical(method[1L], "mcmc") && identical(sampler[1L], "gibbs")) {
			return(.lame_chains_from_lame(environment(), .chain_ctx, chains))
		}
		if (chains > 1L && identical(sampler[1L], "tempered")) {
			cli::cli_abort(
				"{.code sampler = \"tempered\"} runs one chain; set {.code chains = 1}.", class = "lame_chains_arg"
			)
		}
	}
	if (!is.null(attr(prior, "lame_chain_block", exact = TRUE))) attr(prior, "lame_chain_block") <- NULL
	.lame_input_prior <- attr(prior, "lame_checkpoint_input_prior", exact = TRUE)
	if (!is.null(.lame_input_prior)) {
		attr(prior, "lame_checkpoint_input_prior") <- NULL
	} else .lame_input_prior <- prior
	# keep the user's priors and design for training-window refits, so defaults resolved from the
	# complete panel do not leak future data
	.lame_refit_options <- mget(
		c(
			"rvar", "cvar", "dcor", "nvar", "R", "R_row", "R_col",
			"mode", "dynamic_uv", "dynamic_ab", "dynamic_G", "dynamic_beta", "dynamic_rho",
			"dynamic_beta_kind", "dynamic_uv_kind", "family", "intercept", "symmetric",
			"odmax", "prior", "g", "seed", "nscan", "burn", "odens", "gof",
			"dynamic_beta_pool", "dynamic_beta_per_actor", "per_actor_covariate_idx",
			"per_actor_identifiability", "keep_per_actor", "time_index", "period_exposure",
			"ordinal_cutpoints", "sampler"
		),
		envir = environment(),
		inherits = FALSE
	)
	.lame_refit_options$prior <- .lame_input_prior
	.lame_refit_data <- list(Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)
	# snapshot the user's data before internal conversion so update.lame can refit the original input
	# (only with freeze_call = TRUE)
	if (isTRUE(freeze_call)) {
		.lame_freeze_snapshot <- list(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol
		)
	} else {
		.lame_freeze_snapshot <- NULL
	}
	mode <- match.arg(mode)
	bip <- identical(mode, "bipartite")
	# validate the family early so every route, ALS included, gets the same amen-era aliases
	valid_families <- c(
		"normal", "binary", "ordinal",
		"cbin", "frn", "poisson"
	)
	amen_aliases <- c(
		nrm = "normal", bin = "binary", ord = "ordinal",
		pois = "poisson"
	)
	if (length(family) == 1L && is.character(family) &&
			family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		# announce the alias once per session
		if (!isTRUE(.lame_state$family_alias_msg_shown)) {
			cli::cli_inform(c(
				"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
				"i" = "Use {.code family = \"{new_family}\"} to avoid the alias message.",
				"i" = "This note is shown once per session."
			))
			.lame_state$family_alias_msg_shown <- TRUE
		}
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
			!family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."
		))
	}
	# cbin and frn model each actor's ranked nominations row by row, so they are directed by construction
	if (isTRUE(symmetric) && !bip && family %in% c("cbin", "frn")) {
		cli::cli_abort(c(
			"{.arg symmetric} = TRUE is not supported for {.code family = \"{family}\"}.",
			"i" = "The {.val cbin} and {.val frn} families model each actor's own ranked nominations, so they are directed by construction.",
			"i" = "Fit with {.code symmetric = FALSE}, or recode the ties and use {.val binary} / {.val ordinal} for a symmetric fit."
		))
	}
	# dynamic_beta_kind, with "random_walk" as an alias for "rw1"
	if (length(dynamic_beta_kind) > 1L) dynamic_beta_kind <- dynamic_beta_kind[[1L]]
	if (identical(dynamic_beta_kind, "random_walk")) dynamic_beta_kind <- "rw1"
	if (!dynamic_beta_kind %in% c("ar1", "rw1", "rw2", "matern32")) {
		cli::cli_abort(c(
			"{.arg dynamic_beta_kind} must be one of {.val ar1}, {.val rw1}, {.val rw2}, {.val matern32}.",
			"i" = "Got {.val {dynamic_beta_kind}}.",
			"i" = "{.val random_walk} is accepted as an alias for {.val rw1}."
		))
	}
	dynamic_beta_kind_requested <- dynamic_beta_kind
	# dynamic_uv_kind: "ar1" smooth drift, or the "snap" and "t" transition models
	dynamic_uv_kind <- match.arg(dynamic_uv_kind)
	if (!isTRUE(dynamic_uv) && !identical(dynamic_uv_kind, "ar1")) {
		cli::cli_abort(c(
			"{.arg dynamic_uv_kind = \"{dynamic_uv_kind}\"} requires {.arg dynamic_uv = TRUE}.",
			"i" = "Set {.arg dynamic_uv = TRUE} or leave {.arg dynamic_uv_kind} at {.val ar1}."
		))
	}
	dynamic_uv_kind_requested <- dynamic_uv_kind
	if (!is.logical(dynamic_rho) || length(dynamic_rho) != 1L || is.na(dynamic_rho)) {
		cli::cli_abort("{.arg dynamic_rho} must be TRUE or FALSE.")
	}
	dynamic_beta_requested_global <- !is.null(dynamic_beta) &&
		!(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta))
	# choose the static, dynamic, or snap point estimator
	method <- match.arg(method)
	if (identical(sampler, "tempered")) return(.lame_tempered_fit(environment()))
	als_stability <- match.arg(als_stability)
	if (!is.numeric(als_max_iter) || length(als_max_iter) != 1L ||
			anyNA(als_max_iter) || als_max_iter <= 0 ||
			!isTRUE(all.equal(als_max_iter, round(als_max_iter)))) {
		cli::cli_abort("{.arg als_max_iter} must be a positive integer.")
	}
	als_max_iter <- as.integer(als_max_iter)
	if (!is.null(als_tol) &&
			(!is.numeric(als_tol) || length(als_tol) != 1L ||
					!is.finite(als_tol) || als_tol <= 0)) {
		cli::cli_abort("{.arg als_tol} must be NULL or a positive finite number.")
	}
	# ALS tuning knobs default to NULL so each route keeps its own default and refuses values it cannot honour
	als_opt <- function(value, arg, allowed) {
		if (is.null(value)) return(NULL)
		if (!is.character(value) || length(value) != 1L || !value %in% allowed) {
			cli::cli_abort("{.arg {arg}} must be NULL or one of {.val {allowed}}.")
		}
		value
	}
	als_lowrank_method <- als_opt(
		als_lowrank_method,
		"als_lowrank_method",
		c("mm", "als", "hybrid")
	)
	als_linear_solver <- als_opt(
		als_linear_solver,
		"als_linear_solver",
		c("eigen", "qr", "auto")
	)
	als_multistart <- als_opt(
		als_multistart,
		"als_multistart",
		c("none", "cheap", "full")
	)
	if (identical(method, "als")) {
		# check family support before the intercept guard, since the ordinal default sets intercept = FALSE itself
		if (!family %in% .ae_supported_families) {
			supported <- .ae_supported_families
			cli::cli_abort(c(
				"family = {.val {family}} is not supported by the ALS estimator.",
				"i" = "Supported families: {.val {supported}}.",
				"i" = "Use {.code method = \"mcmc\"} to fit {.val {family}} models longitudinally."
			))
		}
		# poisson responses are counts on the ALS route too
		if (identical(family, "poisson")) .lame_check_poisson_counts(Y)
		if (isTRUE(dynamic_rho)) {
			cli::cli_abort(
				c(
					"{.arg dynamic_rho = TRUE} is not supported with {.code method = \"als\"}.",
					"i" = "ALS uses one pooled residual dyadic correlation; its per-period correlations are residual diagnostics.",
					"i" = "Use {.code method = \"mcmc\"} with static latent positions ({.code dynamic_uv = FALSE}) for period-specific reciprocity."
				),
				class = "lame_als_dynamic_rho"
			)
		}
		if (!isTRUE(intercept)) {
			cli::cli_abort(c(
				"{.fn lame} {.arg method} = {.val als} does not currently support {.code intercept = FALSE}.",
				"i" = "Use {.code method = \"mcmc\"} for a no-intercept longitudinal fit, or keep the ALS intercept and center covariates as needed."
			))
		}
		bootstrap_type <- match.arg(bootstrap_type)
		if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
				!isTRUE(all.equal(bootstrap, round(bootstrap)))) {
			cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
		}
		bootstrap <- as.integer(bootstrap)
		multistart <- match.arg(multistart)
		bias_correct <- .ae_check_bias_correct(bias_correct)
		# only the static lame_als() route reruns from several starts
		warn_als_no_multistart <- function(route) {
			if (multistart != "none") {
				cli::cli_warn(c(
					"{.arg multistart} is not available for {route} and is ignored.",
					"i" = "{.arg als_stability} is still available for start-sensitivity checks."
				))
			}
		}
		dynamic_beta_requested <- !is.null(dynamic_beta) &&
			!(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta))
		als_user_args <- names(as.list(mc))[-1L]
		warn_als_mcmc_controls <- function(route) {
			ignored <- intersect(
				als_user_args,
				c(
					"nscan", "burn", "odens", "g", "keep_snap_draws",
					"save_log_lik", "posterior_opts",
					"log_lik_path", "log_lik_chunk_size",
					"log_lik_method", "start_vals", "periodic_save", "out_file",
					"save_interval", "freeze_call",
					"dynamic_beta_per_actor", "per_actor_covariate_idx",
					"per_actor_identifiability", "keep_per_actor", "time_index",
					"period_exposure", "max_seconds", "checkpoint_path",
					"checkpoint_every", "ordinal_cutpoints"
				)
			)
			if (length(ignored) > 0L) {
				ignored_text <- paste(ignored, collapse = ", ")
				cli::cli_warn(c(
					"{.fn lame} {.arg method} = {.val als}: ignored MCMC-only arguments on the {route} route: {ignored_text}.",
					"i" = "Use {.code method = \"mcmc\"} if you need those sampler controls."
				))
			}
		}
		# no ALS route for the snap kind (.lame_snap_als_refuse), nor for the t kind on bipartite panels
		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap")) {
			.lame_snap_als_refuse(if (bip) "bipartite" else if (isTRUE(symmetric)) "symmetric" else "directed")
		}
		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") && bip) {
			.lda_refuse_bipartite_t()
		}
		dyn_requested <- character()
		if (isTRUE(dynamic_uv)) dyn_requested <- c(dyn_requested, "dynamic_uv")
		if (isTRUE(dynamic_ab)) dyn_requested <- c(dyn_requested, "dynamic_ab")
		if (isTRUE(dynamic_rho)) dyn_requested <- c(dyn_requested, "dynamic_rho")
		if (dynamic_beta_requested) {
			dyn_requested <- c(dyn_requested, "dynamic_beta")
		}
		# static lame_als() uses one latent rank
		R_als <- R
		R_row_als <- R_row
		R_col_als <- R_col
		if (!is.null(R_row) || !is.null(R_col)) {
			rr <- R_row %||% R; rc <- R_col %||% R
			if (rr != rc && !(length(dyn_requested) > 0L &&
					identical(mode, "bipartite"))) {
				cli::cli_warn(c(
					"{.arg R_row} ({rr}) != {.arg R_col} ({rc}), but the ALS estimator uses a single latent dimension.",
					"i" = "Using {.code R = {max(rr, rc)}} for the ALS fit (the larger of the two)."
				))
			}
			R_als <- max(rr, rc, R)
			R_row_als <- rr
			R_col_als <- rc
		}
		if (length(dyn_requested) > 0L) {
			warn_als_mcmc_controls("dynamic ALS")
			warn_als_no_multistart("dynamic ALS")
			if (!is.null(als_linear_solver)) {
				cli::cli_abort(c(
					"{.arg als_linear_solver} is not available for dynamic ALS fits.",
					"i" = "The dynamic route solves whole-path blocks with its own eigen solver.",
					"i" = "Drop the argument, or use the static ALS route (no dynamic arguments) where it applies."
				))
			}
			if (!is.null(als_multistart)) {
				cli::cli_abort(c(
					"{.arg als_multistart} is not available for dynamic ALS fits.",
					"i" = "Use {.arg als_stability} for start-sensitivity checks on this route."
				))
			}
			dyn_args <- list(
				Y = Y,
				Xdyad = Xdyad,
				Xrow = Xrow,
				Xcol = Xcol,
				R = R_als,
				R_row = R_row_als,
				R_col = R_col_als,
				family = family,
				mode = mode,
				symmetric = symmetric,
				dynamic_ab = dynamic_ab,
				dynamic_uv = dynamic_uv,
				dynamic_uv_kind = dynamic_uv_kind,
				dynamic_beta = dynamic_beta,
				dynamic_beta_kind = dynamic_beta_kind,
				prior = prior,
				bootstrap = bootstrap,
				bootstrap_type = bootstrap_type,
				bootstrap_block_length = bootstrap_block_length,
				bootstrap_seed = bootstrap_seed,
				stability = als_stability,
				bias_correct = bias_correct,
				max_iter = als_max_iter,
				tol = if (is.null(als_tol)) 5e-5 else als_tol,
				verbose = verbose,
				seed = seed,
				call = mc
			)
			# pass the factor solver only when one was asked for, so the fitter can tell it from its own default
			if (!is.null(als_lowrank_method)) {
				dyn_args$lowrank_method <- als_lowrank_method
			}
			# pass the sampler switches only when the caller set them; lame_dynamic_als() names those its route cannot honour
			for (nm in intersect(als_user_args, c("rvar", "cvar", "nvar", "dcor", "dynamic_beta_pool"))) {
				dyn_args[[nm]] <- get(nm, inherits = FALSE)
			}
			# quote = TRUE: the list already holds values, and `call` would otherwise be re-evaluated
			fit <- do.call(lame_dynamic_als, dyn_args, quote = TRUE)
			fit$call <- mc
			return(fit)
		}
		user_args <- names(as.list(mc))[-1L]
		forwarded <- c(
			"Y", "Xdyad", "Xrow", "Xcol", "R", "R_row", "R_col",
			"family", "mode", "symmetric",
			"bootstrap", "bootstrap_type",
			"bootstrap_block_length", "bootstrap_seed",
			"als_stability", "als_max_iter", "als_tol",
			"als_lowrank_method", "als_linear_solver",
			"als_multistart", "multistart", "bias_correct", "dcor",
			"verbose", "seed", "method", "intercept", "odmax", "prior",
			# the dynamic flags choose between the static and dynamic ALS estimators
			"dynamic_uv", "dynamic_ab", "dynamic_G", "dynamic_rho",
			"dynamic_beta", "dynamic_beta_kind", "dynamic_uv_kind",
			# cosmetic MCMC args with no ALS equivalent
			"plot", "model.name"
		)
		dropped <- setdiff(user_args, forwarded)
		if (length(dropped) > 0L) {
			cli::cli_warn(c(
				"{.fn lame} {.arg method} = {.val als}: ignored {length(dropped)} arguments that apply only to the MCMC path: {.arg {dropped}}.",
				"i" = "Use {.code method = \"mcmc\"} (the default) if you need custom MCMC-only features."
			))
		}
		fit <- lame_als(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			R = R_als,
			family = family,
			mode = mode,
			symmetric = symmetric,
			bootstrap = bootstrap,
			bootstrap_type = bootstrap_type,
			bootstrap_block_length = bootstrap_block_length,
			bootstrap_seed = bootstrap_seed,
			multistart = if (is.null(als_multistart)) multistart else als_multistart,
			bias_correct = bias_correct,
			lowrank_method = if (is.null(als_lowrank_method)) "mm" else
				als_lowrank_method,
			linear_solver = if (is.null(als_linear_solver)) "eigen" else
				als_linear_solver,
			max_iter = als_max_iter,
			tol = if (is.null(als_tol)) 1e-6 else als_tol,
			verbose = verbose,
			seed = seed,
			dcor = dcor,
			prior = prior
		)
		fit$call <- mc
		return(fit)
	}
	# dynamic "ar1" positions start from the deterministic dynamic ALS solution unless start values are given,
	# since the cold start can settle on too little persistence; prior$als_start = FALSE keeps the cold start.
	# a start that is not finite or leaves prior$uv_max_abs is dropped
	als_start_used <- FALSE
	dyn_beta_req <- !is.null(dynamic_beta) &&
		!(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta))
	if (is.null(start_vals) && isTRUE(dynamic_uv) && dynamic_uv_kind %in% c("ar1", "t") &&
			!isTRUE(symmetric) && !identical(mode, "bipartite") &&
			family %in% c("normal", "poisson", "binary") &&
			isTRUE(intercept) && !isTRUE(dynamic_rho) &&
			is.null(dynamic_beta_per_actor) && !dyn_beta_req &&
			!isFALSE(prior[["als_start"]])) {
		als_call <- .mc_eval
		als_call$method <- "als"
		# a "t" chain starts from the AR(1) solution
		als_call$dynamic_uv_kind <- "ar1"
		als_call$verbose <- FALSE
		als_call$plot <- FALSE
		als_call$gof <- FALSE
		for (nm in intersect(
			names(als_call),
			c(
				"nscan", "burn", "odens", "seed", "start_vals", "periodic_save",
				"out_file", "save_interval", "model.name", "save_log_lik",
				"posterior_opts", "log_lik_method", "log_lik_path",
				"log_lik_chunk_size", "keep_snap_draws", "g", "bootstrap",
				"bootstrap_type", "multistart", "bias_correct", "prior"
			)
		)) {
			als_call[[nm]] <- NULL
		}
		# the start needs a point estimate only, so skip the bias correction (which refits half samples)
		als_call$bias_correct <- "none"
		t_als <- proc.time()[["elapsed"]]
		# a fit that stops short of its criterion still gives a usable start; the checks below drop unusable ones
		als_fit <- tryCatch(
			suppressWarnings(suppressMessages(eval(als_call, parent.frame()))),
			lame_als_convergence = function(e) e$fit,
			error = function(e) NULL
		)
		if (!is.null(als_fit)) {
			sv <- tryCatch(als_start_vals(als_fit), error = function(e) NULL)
			sv_num <- if (is.list(sv)) unlist(
				lapply(sv, function(z) if (is.numeric(z)) as.numeric(z) else NULL)
			) else NA_real_
			uv_cap <- as.numeric(prior$uv_max_abs %||% Inf)
			sv_ok <- is.list(sv) && !is.null(sv$U) && length(sv_num) > 0L &&
				all(is.finite(sv_num)) &&
				max(abs(sv$U), abs(sv[["V"]] %||% 0)) <= uv_cap
			if (isTRUE(sv_ok)) {
				start_vals <- sv
				als_start_used <- TRUE
				if (isTRUE(verbose)) {
					cli::cli_inform("Dynamic positions start from a dynamic ALS fit ({round(proc.time()[['elapsed']] - t_als, 1)} s).")
				}
			}
		}
		rm(als_call, als_fit)
	}
	log_lik_method <- match.arg(log_lik_method)
	keep_snap_draws <- match.arg(keep_snap_draws)
	if (identical(log_lik_method, "observed_ghk")) {
		# the GHK marginal log-likelihood is a Monte Carlo estimate whose log is biased downward; say so once per fit
		cli::cli_inform(c(
			"i" = "{.arg log_lik_method} = {.val observed_ghk}: GHK Monte Carlo marginal for rank/ordinal families.",
			"!" = "The log of an unbiased probability estimate is biased downward; the stored scores are Monte Carlo estimates."
		))
	}
	# ordinal cutpoints use either data-induced thresholds or explicit MH updates
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	if (identical(ordinal_cutpoints, "explicit") && family != "ordinal") {
		cli::cli_warn(c(
			"{.arg ordinal_cutpoints} = {.val explicit} only takes effect when {.code family = \"ordinal\"}.",
			"i" = "Argument ignored for {.code family = \"{family}\"}."
		))
		ordinal_cutpoints <- "data_induced"
	}
	.lame_explicit_ordinal_support(
		family,
		ordinal_cutpoints,
		bip,
		symmetric,
		dcor || isTRUE(dynamic_rho)
	)
	.lame_validate_score_request(
		family,
		ordinal_cutpoints,
		log_lik_method,
		isTRUE(save_log_lik) || identical(save_log_lik, "chunked")
	)
	dynamic_beta_pool <- match.arg(dynamic_beta_pool)
	# dynamic_beta_pool "none" runs per-block updates; "rho", "sigma" and "both" add a hierarchical pool across blocks
	if (!is.null(time_index)) {
		if (!is.numeric(time_index)) {
			cli::cli_abort("{.arg time_index} must be numeric.")
		}
		if (any(diff(time_index) <= 0)) {
			cli::cli_abort("{.arg time_index} must be strictly increasing.")
		}
		# unequal gaps use the joint Gaussian sampler
	}
	if (!is.null(period_exposure)) {
		if (!is.numeric(period_exposure) || any(!is.finite(period_exposure)) ||
				any(period_exposure < 0)) {
			cli::cli_abort("{.arg period_exposure} must be a non-negative finite numeric vector.")
		}
		if (family != "poisson" && any(period_exposure != 1)) {
			cli::cli_abort(c(
				"{.arg period_exposure} != 1 is currently meaningful only for {.code family = \"poisson\"}.",
				"i" = "Got {.val {family}} family with non-trivial exposure.",
				"i" = "For other families, scale {.arg Y} manually or use a covariate offset."
			))
		}
	}
	# poisson exposures are active only when some value differs from one
	poisson_exposure_active <- !is.null(period_exposure) &&
		identical(family, "poisson") &&
		any(period_exposure != 1)
	# dynamic_beta_per_actor: pairwise contrast FFBS with sum-to-zero per period
	per_actor_identifiability <- match.arg(per_actor_identifiability)
	keep_per_actor <- match.arg(keep_per_actor)
	per_actor_active <- !is.null(dynamic_beta_per_actor)
	if (per_actor_active) {
		if (!is.character(dynamic_beta_per_actor) ||
				!dynamic_beta_per_actor %in% c("row", "col")) {
			cli::cli_abort("{.arg dynamic_beta_per_actor} must be one of {.val NULL}, {.val row}, or {.val col}.")
		}
		# drop_population is reserved; the fitted path stays centered
		if (identical(per_actor_identifiability, "drop_population")) {
			cli::cli_warn(c(
				"!" = "{.arg per_actor_identifiability} = {.val drop_population} is reserved but not available yet.",
				"i" = "Using {.val center}: zero-sum per-actor deviations on top of the population coefficient.",
				"i" = "Pass {.val center} to silence this warning."
			))
		}
		if (isTRUE(symmetric)) cli::cli_abort("Actor-specific row or column coefficient paths require a directed or bipartite model.")
		if (!is.null(time_index) && any(abs(diff(as.numeric(time_index)) - 1) > 1e-10))
			cli::cli_abort("Actor-specific coefficient paths currently require unit-spaced periods.")
		if (!is.numeric(per_actor_covariate_idx) || length(per_actor_covariate_idx) != 1L ||
				!is.finite(per_actor_covariate_idx) || per_actor_covariate_idx != floor(per_actor_covariate_idx) ||
				per_actor_covariate_idx < 1L) {
			cli::cli_abort("{.arg per_actor_covariate_idx} must be a positive integer.")
		}
	}
	# max_seconds and checkpointing
	if (!is.null(max_seconds) && (!is.numeric(max_seconds) || max_seconds <= 0)) {
		cli::cli_abort("{.arg max_seconds} must be a positive number or {.val Inf}.")
	}
	if (!is.null(checkpoint_path)) {
		ck_dir <- dirname(checkpoint_path)
		if (!dir.exists(ck_dir)) dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
	}
	if (!is.numeric(checkpoint_every) || checkpoint_every <= 0) {
		cli::cli_abort("{.arg checkpoint_every} must be a positive integer.")
	}
	.lame_start_time <- Sys.time()
	.lame_terminated_early <- FALSE
	# convert a 3d Y array to list format
	if (is.array(Y) && length(dim(Y)) == 3L) {
		dn <- dimnames(Y)
		Y <- lapply(seq_len(dim(Y)[3]), function(t) {
			mat <- Y[, , t]
			if (!is.null(dn)) {
				rownames(mat) <- dn[[1]]
				colnames(mat) <- dn[[2]]
			}
			mat
		})
		if (!is.null(dn) && !is.null(dn[[3]])) {
			names(Y) <- dn[[3]]
		}
	}
	# split a 4-D Xdyad ([n, n, p, T], the amen layout) into a list of [n, n, p] arrays, keeping slice names
	if (!is.null(Xdyad) && is.array(Xdyad) && length(dim(Xdyad)) == 4L) {
		dn <- dimnames(Xdyad)
		Xdyad <- lapply(seq_len(dim(Xdyad)[4L]), function(tt) {
			sl <- array(Xdyad[, , , tt], dim = dim(Xdyad)[1:3])
			if (!is.null(dn)) dimnames(sl) <- dn[1:3]
			sl
		})
		if (!is.null(dn) && !is.null(dn[[4L]])) names(Xdyad) <- dn[[4L]]
	}
	# validate Y and the key MCMC arguments up front
	if (!is.list(Y) || length(Y) == 0L) {
		cli::cli_abort("{.arg Y} must be a non-empty list of relational matrices (or a 3D array).")
	}
	obs_y <- unlist(lapply(Y, function(y) y[is.finite(y)]), use.names = FALSE)
	if (length(obs_y) == 0L) {
		cli::cli_abort("{.arg Y} has no observed (non-missing, finite) values.")
	}
	if (length(nscan) != 1L || anyNA(nscan) || nscan < 1L) {
		cli::cli_abort("{.arg nscan} must be a positive integer (post-burn-in iterations).")
	}
	if (length(burn) != 1L || anyNA(burn) || burn < 0L) {
		cli::cli_abort("{.arg burn} must be a non-negative integer.")
	}
	if (length(odens) != 1L || anyNA(odens) || odens < 1L || odens > nscan) {
		cli::cli_abort(c(
			"{.arg odens} must be a positive integer no greater than {.arg nscan} = {nscan}.",
			"i" = "Otherwise no posterior draws would be stored."
		))
	}
	# poisson responses are non-negative integer counts
	if (family == "poisson") .lame_check_poisson_counts(obs_y)
	# an entirely missing row or column in a slice would abort the bipartite binary z sampler; refuse it cleanly
	.lame_empty_periods <- list()
	.lame_check_all_na <- function(Y_per, t_label) {
		if (is.null(Y_per)) return(invisible())
		drow <- apply(Y_per, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y_per, 2, function(c) all(is.na(c) | !is.finite(c)))
		if (bip) {
			if (any(drow) || any(dcol)) {
				if (family %in% c("binary")) {
					cli::cli_abort(c(
						"{.arg Y[[{t_label}]]} (bipartite, {.val {family}}) has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
						"x" = "The bipartite {.val {family}} sampler cannot proceed with rows/columns that have no observed data.",
						"i" = "Drop the empty actor(s) from {.arg Y} (and matching covariates) before fitting, or use {.val normal} family."
					))
				}
				.lame_empty_periods[[length(.lame_empty_periods) + 1L]] <<-
					data.frame(
						period = t_label,
						rows = sum(drow),
						cols = sum(dcol),
						stringsAsFactors = FALSE
					)
			}
		} else {
			# unipartite: the diagonal does not count toward an all-missing row
			drow_ex <- vapply(
				seq_len(nrow(Y_per)),
				function(i)
					all(is.na(Y_per[i, -i]) | !is.finite(Y_per[i, -i])),
				logical(1)
			)
			dcol_ex <- vapply(
				seq_len(ncol(Y_per)),
				function(j)
					all(is.na(Y_per[-j, j]) | !is.finite(Y_per[-j, j])),
				logical(1)
			)
			if (any(drow_ex) || any(dcol_ex)) {
				.lame_empty_periods[[length(.lame_empty_periods) + 1L]] <<-
					data.frame(
						period = t_label,
						rows = sum(drow_ex),
						cols = sum(dcol_ex),
						stringsAsFactors = FALSE
					)
			}
		}
	}
	for (.t in seq_along(Y)) {
		.lame_check_all_na(Y[[.t]], if (!is.null(names(Y))) names(Y)[.t] else as.character(.t))
	}
	if (length(.lame_empty_periods) > 0L) {
		.empty_summary <- do.call(rbind, .lame_empty_periods)
		first_periods <- paste(utils::head(.empty_summary$period, 5L), collapse = ", ")
		n_empty <- nrow(.empty_summary)
		row_max <- max(.empty_summary$rows)
		col_max <- max(.empty_summary$cols)
		period_word <- if (n_empty == 1L) "period" else "periods"
		row_word <- if (row_max == 1L) "row" else "rows"
		col_word <- if (col_max == 1L) "column" else "columns"
		period_label <- if (n_empty == 1L) "First affected period" else "First affected periods"
		period_suffix <- if (n_empty > 5L) " ..." else ""
		cli::cli_inform(c(
			"i" = "{.arg Y} has fully-unobserved rows or columns in {n_empty} {period_word}; largest count is {row_max} {row_word} and {col_max} {col_word}.",
			"i" = "{period_label}: {first_periods}{period_suffix}. Their additive effects are pooled with the prior in those periods."
		))
	}
	# a constant response cannot be fit under any family (as in ame())
	{
		eff_y <- if (family %in% c("binary", "cbin")) 1 * (obs_y > 0) else obs_y
		if (length(unique(eff_y)) == 1L) {
			cli::cli_abort(c(
				"{.arg Y} is constant; there is no relational variation to model.",
				"i" = "Every observed cell equals {.val {unique(eff_y)}}."
			))
		}
	}
	# binary/cbin threshold via 1*(y>0); warn when counts are collapsed
	if (family %in% c("binary", "cbin")) {
		uy <- unique(obs_y)
		if (any(!uy %in% c(0, 1))) {
			cli::cli_warn(c(
				"{.arg family} = {.val {family}} but {.arg Y} contains values other than 0/1.",
				"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val ordinal}/{.val normal} as appropriate."
			))
		}
	}
	if (family == "normal") {
		uy <- unique(obs_y)
		if (length(uy) <= 2L && all(uy %in% c(0, 1))) {
			cli::cli_inform(c(
				"i" = "{.arg Y} looks binary (values in {.val {{0, 1}}}) but {.arg family} = {.val normal}.",
				"i" = "If this is a 0/1 network, pass {.code family = \"binary\"}."
			))
		}
	}
	if (burn > nscan) {
		cli::cli_warn(c(
			"{.arg burn} = {burn} is greater than {.arg nscan} = {nscan}.",
			"i" = "Most iterations will be discarded; consider {.arg burn} <= {.arg nscan}."
		))
	}
	# AR(1) dynamics need at least two time slices
	if ((isTRUE(dynamic_uv) || isTRUE(dynamic_ab) || isTRUE(dynamic_rho)) && length(Y) < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_uv}, {.arg dynamic_ab}, and {.arg dynamic_rho} require at least 2 time periods.",
			"i" = "Got {length(Y)} slice{?s}; use {.fn ame} for a single cross-section, or supply more periods."
		))
	}
	# events per actor-period, the information that per-period effects, the level of a saturating
	# panel and the dyadic correlation rest on; warn up front when per-period effects are prior-driven
	.per_actor_period <- NA_real_
	.per_actor_total <- NA_real_
	if (family %in% c("binary", "cbin", "frn", "poisson") && is.list(Y) && length(Y) >= 2L) {
		.events_per_period <- suppressWarnings(vapply(Y, function(m) {
			sum(as.matrix(m) > 0, na.rm = TRUE)
		}, numeric(1)))
		.rows_per_period <- vapply(Y, NROW, numeric(1))
		.events_per_actor_period <- ifelse(
			.rows_per_period > 0,
			.events_per_period / pmax(.rows_per_period, 1),
			NA_real_
		)
		.per_actor_period <- mean(.events_per_actor_period, na.rm = TRUE)
		# events per row actor over the panel, an actor present in fewer periods counting at its own share
		.row_names <- lapply(Y, rownames)
		.n_row_actors <- if (all(!vapply(.row_names, is.null, logical(1)))) {
			length(unique(unlist(.row_names)))
		} else max(.rows_per_period)
		if (.n_row_actors > 0) .per_actor_total <- sum(.events_per_period, na.rm = TRUE) / .n_row_actors
	}
	# a supplied effects prior (Sab0, eta0 or etaab): the effect variances are reported under it
	.ab_prior_supplied <- .lame_ab_prior_given(.lame_input_prior)
	# share of rows at the nomination cap and whether cells are hidden, read by the reporting rules
	.odmax_sat_share <- .lame_odmax_saturated_share(Y, odmax, family)
	.odmax_n_per <- if (is.list(Y)) length(Y) else if (length(dim(Y)) == 3L) dim(Y)[3L] else 1L
	.cells_hidden <- family %in% c("cbin", "frn") && any(vapply(
		if (is.list(Y)) Y else lapply(seq_len(dim(Y)[3L]), function(tt) Y[, , tt]),
		function(m) any(.lame_nom_hidden_cells(m, bip)),
		logical(1)
	))
	# what the fit will not report, said before it runs
	for (st_ in .lame_ab_statements(list(
		family = family,
		mode = if (bip) "bipartite" else "unipartite",
		dynamic_ab = isTRUE(dynamic_ab),
		events_per_actor_period = .per_actor_period,
		odmax_saturated_share = .odmax_sat_share,
		cells_hidden = .cells_hidden,
		n_time = .odmax_n_per,
		prior_supplied = .ab_prior_supplied,
		R = R,
		ab_prior_text = .lame_ab_prior_text(list(prior = .lame_input_prior))
	))) {
		cli::cli_inform(c("i" = "{st_}"))
	}
	# the dynamic position kernel whitens each dyad's residual pair by one rho, so it cannot carry a period-varying rho
	if (isTRUE(dcor) && !isTRUE(bip) && !isTRUE(symmetric) && isTRUE(dynamic_uv) &&
			isTRUE(dynamic_rho)) {
		cli::cli_abort(c(
			"{.arg dynamic_rho} is not available with {.arg dynamic_uv = TRUE}.",
			"i" = "The dynamic position kernel whitens each dyad pair by a single {.code rho}; fit static positions ({.code dynamic_uv = FALSE}) for a period-varying {.code rho}."
		))
	}
	# note a panel with few mutual dyads; this does not change the likelihood
	if (isTRUE(dcor) && !isTRUE(bip) && !isTRUE(symmetric) &&
			family %in% c("binary", "cbin", "frn", "poisson") && is.list(Y) && length(Y) >= 1L) {
		.mutual_per_period <- suppressWarnings(vapply(Y, function(m) {
			m <- as.matrix(m)
			if (nrow(m) != ncol(m)) return(NA_real_)
			sum(m > 0 & t(m) > 0, na.rm = TRUE) / 2
		}, numeric(1)))
		.mutual <- mean(.mutual_per_period, na.rm = TRUE)
		if (is.finite(.mutual) && .mutual < 1) {
			cli::cli_inform(c(
				"i" = "This panel has about {round(.mutual, 1)} mutual dyads per period, so the data inform the dyadic correlation {.code rho} mainly through the absence of mutual ties.",
				"i" = .lame_rho_few_mutual_text
			))
		}
		# fewer than 30 mutual dyads in total, with at least one per period on average
		.mutual_total <- sum(.mutual_per_period, na.rm = TRUE)
		if (is.finite(.mutual) && .mutual >= 1 && .mutual_total < 30) {
			cli::cli_inform(c(
				"i" = "This panel has {round(.mutual_total)} mutual dyads in total, so the data inform the dyadic correlation {.code rho} mainly through the absence of mutual ties.",
				"i" = .lame_rho_few_mutual_text
			))
		}
	}
	if (isTRUE(dynamic_rho)) {
		if (isTRUE(bip) || isTRUE(symmetric) || !isTRUE(dcor)) {
			cli::cli_abort("{.arg dynamic_rho} is for directed unipartite models with {.arg dcor = TRUE}.")
		}
		if (!identical(family, "normal")) {
			cli::cli_abort("{.arg dynamic_rho} currently supports {.code family = \"normal\"}.")
		}
		if (isTRUE(dynamic_beta_requested_global)) {
			cli::cli_abort(c(
				"{.arg dynamic_rho} is not yet combined with {.arg dynamic_beta}.",
				"i" = "Fit time-varying reciprocity with static coefficients, or keep {.arg dynamic_beta} and scalar {.code rho}."
			))
		}
	}
	# dynamic_beta needs at least two periods, like the other dynamics
	if (!is.null(dynamic_beta) &&
			!(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta)) &&
			length(Y) < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_beta} requires at least 2 time periods.",
			"i" = "Got {length(Y)} slice{?s}; use {.fn ame} for a single cross-section."
		))
	}
	if (isTRUE(bip)) .lame_refuse_bipartite_suv0(prior)
	# unknown prior-list keys usually mean a typo
	if (length(prior) > 0L && !is.null(names(prior))) {
		known_prior <- c(
			"Sab0", "eta0", "etaab", "Suv0", "kappa0",
			"s20", "s2u0",
			"ab_min_observed",
			"rho_uv_mean", "rho_uv_sd",
			"rho_ab_mean", "rho_ab_sd",
			"sigma_uv_shape", "sigma_uv_scale",
			"uv_max_abs", "G_max", "sigma_uv_max",
			"sigma_ab_shape", "sigma_ab_scale",
			"sigma_actor_shape", "sigma_actor_scale",
			"rho_actor_shape1", "rho_actor_shape2",
			"ab_var_cap",
			"uv_var_cap",
			"rho_beta_mean", "rho_beta_sd",
			"sigma_beta_shape", "sigma_beta_scale",
			"sigma_beta_init",
			"rho_beta_lower", "rho_beta_upper",
			"beta0_mean", "beta0_var",
			"matern32_length_scale",
			"snap_kappa", "snap_pi_a", "snap_pi_b", "t_nu",
			"ghk_n_mc",
			"symmetric_node_tie",
			"uv_split_move",
			"uv_nc_move",
			"uv_temper_move", "uv_temper_every",
			"cov_ridge_move",
			"sparse_position_move",
			"sparse_position_every",
			"level_move_every",
			"level_move_tries", "marginal_move", "bipartite_scale_move", "actor_factor_slice",
			"actor_move", "actor_move_every", "actor_move_rho", "actor_move_paths", "actor_move_period", "actor_move_room",
			"beta_hyper_nc", "poisson_period_move", "poisson_static_move", "actor_hyper_nc", "ordinal_scale_move",
			"binary_scale_move", "binary_scale_every",
			"dev_move", "dev_move_tries", "als_start",
			"ab_hyper_move", "ab_hyper_move_tries", "uv_g_scale_move", "ab_period_move",
			"chain_extensions"
		)
		bad <- setdiff(names(prior), known_prior)
		if (length(bad) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg prior} list element{?s}: {.val {bad}}.",
				"i" = "Recognised names: {.val {known_prior}}.",
				"i" = "Note: {.code g} is a top-level argument to {.fn lame}, not a {.arg prior} entry."
			))
		}
	}
	.lame_chain_check_prior(prior)
	if (!is.null(prior[["cov_ridge_move"]]) &&
			(!is.logical(prior[["cov_ridge_move"]]) || length(prior[["cov_ridge_move"]]) != 1L ||
					is.na(prior[["cov_ridge_move"]]))) cli::cli_abort("prior$cov_ridge_move must be TRUE or FALSE.")
	if (!is.null(prior[["uv_temper_move"]]) &&
			(!is.logical(prior[["uv_temper_move"]]) || length(prior[["uv_temper_move"]]) != 1L ||
					is.na(prior[["uv_temper_move"]]))) cli::cli_abort("prior$uv_temper_move must be TRUE or FALSE.")
	if (!is.null(prior[["uv_temper_every"]]) &&
			(!is.numeric(prior[["uv_temper_every"]]) || length(prior[["uv_temper_every"]]) != 1L ||
					!is.finite(prior[["uv_temper_every"]]) || prior[["uv_temper_every"]] < 1)) cli::cli_abort("prior$uv_temper_every must be a positive integer.")
	if (!is.null(prior[["uv_nc_move"]]) &&
			(!is.logical(prior[["uv_nc_move"]]) || length(prior[["uv_nc_move"]]) != 1L ||
					is.na(prior[["uv_nc_move"]]))) cli::cli_abort("prior$uv_nc_move must be TRUE or FALSE.")
	if (!is.null(prior[["uv_split_move"]]) &&
			(!is.logical(prior[["uv_split_move"]]) || length(prior[["uv_split_move"]]) != 1L ||
					is.na(prior[["uv_split_move"]]))) cli::cli_abort("prior$uv_split_move must be TRUE or FALSE.")
	if (!is.null(prior[["sparse_position_every"]]) &&
			(!is.numeric(prior[["sparse_position_every"]]) || length(prior[["sparse_position_every"]]) != 1L ||
					!is.finite(prior[["sparse_position_every"]]) || prior[["sparse_position_every"]] < 1))
		cli::cli_abort("prior$sparse_position_every must be a positive integer.")
	if (!is.null(prior[["sparse_position_move"]]) &&
			(!is.logical(prior[["sparse_position_move"]]) || length(prior[["sparse_position_move"]]) != 1L ||
					is.na(
						prior[["sparse_position_move"]]
					))) cli::cli_abort("prior$sparse_position_move must be TRUE or FALSE.")
	if (!is.null(prior[["symmetric_node_tie"]]) &&
			(!is.logical(prior[["symmetric_node_tie"]]) || length(prior[["symmetric_node_tie"]]) != 1L ||
					is.na(prior[["symmetric_node_tie"]]))) cli::cli_abort("prior$symmetric_node_tie must be TRUE or FALSE.")
	if (!is.null(prior[["marginal_move"]]) &&
			(!is.logical(prior[["marginal_move"]]) || length(prior[["marginal_move"]]) != 1L ||
					is.na(prior[["marginal_move"]]))) cli::cli_abort("prior$marginal_move must be TRUE or FALSE.")
	if (!is.null(prior[["actor_move"]]) &&
			(!is.logical(prior[["actor_move"]]) || length(prior[["actor_move"]]) != 1L ||
					is.na(prior[["actor_move"]]))) cli::cli_abort("prior$actor_move must be TRUE or FALSE.")
	for (nm_ in c("ab_hyper_move", "uv_g_scale_move", "ab_period_move")) {
		if (!is.null(prior[[nm_]]) && (!is.logical(prior[[nm_]]) || length(prior[[nm_]]) != 1L || is.na(prior[[nm_]])))
			cli::cli_abort(paste0("prior$", nm_, " must be TRUE or FALSE."))
	}
	if (!is.null(prior[["actor_move_rho"]]) &&
			(!is.logical(prior[["actor_move_rho"]]) || length(prior[["actor_move_rho"]]) != 1L ||
					is.na(prior[["actor_move_rho"]]))) cli::cli_abort("prior$actor_move_rho must be TRUE or FALSE.")
	if (!is.null(prior[["actor_move_paths"]]) &&
			(!is.logical(prior[["actor_move_paths"]]) || length(prior[["actor_move_paths"]]) != 1L ||
					is.na(prior[["actor_move_paths"]]))) cli::cli_abort("prior$actor_move_paths must be TRUE or FALSE.")
	if (!is.null(prior[["actor_move_period"]]) &&
			(!is.logical(prior[["actor_move_period"]]) || length(prior[["actor_move_period"]]) != 1L ||
					is.na(prior[["actor_move_period"]]))) cli::cli_abort("prior$actor_move_period must be TRUE or FALSE.")
	if (!is.null(prior[["actor_move_room"]]) &&
			(!is.logical(prior[["actor_move_room"]]) || length(prior[["actor_move_room"]]) != 1L ||
					is.na(prior[["actor_move_room"]]))) cli::cli_abort("prior$actor_move_room must be TRUE or FALSE.")
	for (nm_ in c(
		"beta_hyper_nc", "poisson_period_move", "poisson_static_move", "actor_hyper_nc", "ordinal_scale_move", "binary_scale_move"
	)) {
		if (!is.null(prior[[nm_]]) && (!is.logical(prior[[nm_]]) || length(prior[[nm_]]) != 1L || is.na(prior[[nm_]])))
			cli::cli_abort(paste0("prior$", nm_, " must be TRUE or FALSE."))
	}
	if (!is.null(prior[["actor_move_every"]]) &&
			(!is.numeric(prior[["actor_move_every"]]) || length(prior[["actor_move_every"]]) != 1L ||
					!is.finite(prior[["actor_move_every"]]) || prior[["actor_move_every"]] < 1 ||
					prior[["actor_move_every"]] != floor(prior[["actor_move_every"]])))
		cli::cli_abort("prior$actor_move_every must be a positive whole number.")
	if (!is.null(prior[["binary_scale_every"]]) &&
			(!is.numeric(prior[["binary_scale_every"]]) || length(prior[["binary_scale_every"]]) != 1L ||
					!is.finite(prior[["binary_scale_every"]]) || prior[["binary_scale_every"]] < 1 ||
					prior[["binary_scale_every"]] != floor(prior[["binary_scale_every"]])))
		cli::cli_abort("prior$binary_scale_every must be a positive whole number.")
	if (!is.null(prior[["actor_factor_slice"]]) &&
			(!is.logical(prior[["actor_factor_slice"]]) || length(prior[["actor_factor_slice"]]) != 1L ||
					is.na(prior[["actor_factor_slice"]]))) cli::cli_abort("prior$actor_factor_slice must be TRUE or FALSE.")
	if (!is.null(prior[["bipartite_scale_move"]]) &&
			(!is.logical(prior[["bipartite_scale_move"]]) || length(prior[["bipartite_scale_move"]]) != 1L ||
					is.na(
						prior[["bipartite_scale_move"]]
					))) cli::cli_abort("prior$bipartite_scale_move must be TRUE or FALSE.")
	# posterior_opts: warn on unknown keys and accept save_UV (as in ame()) as an alias for save_UV_draws
	if (!is.null(posterior_opts)) {
		if (!is.list(posterior_opts)) {
			cli::cli_abort(c(
				"{.arg posterior_opts} must be a list.",
				"i" = "Build it with {.fn posterior_options}, e.g. {.code posterior_opts = posterior_options(save_UV = TRUE)}."
			))
		}
		known_posterior_opts <- c(
			"save_UV", "save_UV_draws", "save_ab",
			"thin_UV", "thin_ab"
		)
		bad_opts <- setdiff(names(posterior_opts), known_posterior_opts)
		if (length(bad_opts) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg posterior_opts} list element{?s}: {.val {bad_opts}}.",
				"i" = "Recognised names: {.val {known_posterior_opts}}."
			))
		}
		if (isTRUE(posterior_opts$save_UV)) {
			posterior_opts$save_UV_draws <- TRUE
		}
		if (isTRUE(posterior_opts$save_ab)) {
			thin_ab <- posterior_opts[["thin_ab"]] %||% 1L
			if (!is.numeric(thin_ab) || length(thin_ab) != 1L ||
					!is.finite(thin_ab) || thin_ab < 1 || thin_ab != floor(thin_ab)) {
				cli::cli_abort("{.code posterior_opts$thin_ab} must be a positive integer.")
			}
			posterior_opts$thin_ab <- as.integer(thin_ab)
		}
	}
	# post-hoc gof(), posterior simulation and the pointwise log-likelihood need the saved additive and
	# latent draws, so keep them by default; a latent cube too large to keep is dropped further down
	posterior_opts_default <- is.null(posterior_opts)
	if (posterior_opts_default) {
		posterior_opts <- list(save_ab = TRUE, thin_ab = 1L, save_UV_draws = TRUE)
	}
	if (isTRUE(periodic_save) && is.null(out_file)) {
		cli::cli_abort(c(
			"{.arg periodic_save} = TRUE requires an {.arg out_file} path.",
			"i" = "Supply {.arg out_file} or set {.arg periodic_save} = FALSE."
		))
	}
	if (nscan %% odens != 0) { cli::cli_abort("{.arg nscan} must be a multiple of {.arg odens}.")}
	# set the seed locally and restore the global RNG stream on exit
	.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (.had_seed) {
			assign(".Random.seed", .old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)
	set.seed(seed)
	# actors are identified by the row/column names of Y; unnamed matrices get default names,
	# assuming a fixed actor composition
	y_named <- vapply(Y, function(y) {
		!is.null(rownames(y)) || !is.null(colnames(y))
	}, logical(1))
	if (!any(y_named)) {
		# give Y and every covariate array the same default actor names
		nr0 <- nrow(Y[[1]]); nc0 <- ncol(Y[[1]])
		row_names_def <- if (bip) paste0("row", seq_len(nr0)) else paste0("a", seq_len(nr0))
		col_names_def <- if (bip) paste0("col", seq_len(nc0)) else row_names_def
		Y <- lapply(Y, function(y) {
			rownames(y) <- row_names_def; colnames(y) <- col_names_def; y
		})
		.apply_default_dimnames <- function(arr, row_nm, col_nm) {
			if (is.null(arr)) return(arr)
			d <- dim(arr)
			if (length(d) == 2L) {
				if (d[1] == length(row_nm)) rownames(arr) <- row_nm
				if (d[2] == length(col_nm)) colnames(arr) <- col_nm
			} else if (length(d) == 3L) {
				dn <- dimnames(arr); if (is.null(dn)) dn <- list(NULL, NULL, NULL)
				if (d[1] == length(row_nm) && is.null(dn[[1]])) dn[[1]] <- row_nm
				if (d[2] == length(col_nm) && is.null(dn[[2]])) dn[[2]] <- col_nm
				dimnames(arr) <- dn
			}
			arr
		}
		if (!is.null(Xdyad)) {
			Xdyad <- lapply(Xdyad, .apply_default_dimnames, row_names_def, col_names_def)
		}
		if (!is.null(Xrow)) {
			Xrow <- lapply(
				Xrow,
				.apply_default_dimnames,
				row_names_def,
				colnames(Xrow[[1]]) %||% character(0)
			)
		}
		if (!is.null(Xcol)) {
			Xcol <- lapply(
				Xcol,
				.apply_default_dimnames,
				col_names_def,
				colnames(Xcol[[1]]) %||% character(0)
			)
		}
	} else if (!all(y_named)) {
		cli::cli_abort(c(
			"Some {.arg Y} slices carry actor names and others do not.",
			"i" = "Give every slice row/column names, or leave them all unnamed (names are then generated automatically)."
		))
	} else {
		# covariate arrays without dimnames take Y's per-slice dimnames, so list_to_array sorts X and Y the same way.
		# Xcol inherits colnames(Y); Xrow and the first dimension of Xdyad inherit rownames(Y)
		.inherit_dimnames <- function(X, row_nm, col_nm) {
			if (is.null(X)) return(X)
			lapply(seq_along(X), function(t) {
				arr <- X[[t]]
				if (is.null(arr)) return(arr)
				rn <- row_nm[[t]] ; cn <- col_nm[[t]]
				d <- dim(arr)
				if (length(d) == 2L) {
					if (is.null(rownames(arr)) && !is.null(rn) && d[1] == length(rn))
						rownames(arr) <- rn
					if (is.null(colnames(arr)) && !is.null(cn) && d[2] == length(cn))
						colnames(arr) <- cn
				} else if (length(d) == 3L) {
					dn <- dimnames(arr); if (is.null(dn)) dn <- list(NULL, NULL, NULL)
					if (is.null(dn[[1L]]) && !is.null(rn) && d[1] == length(rn))
						dn[[1L]] <- rn
					if (is.null(dn[[2L]]) && !is.null(cn) && d[2] == length(cn))
						dn[[2L]] <- cn
					dimnames(arr) <- dn
				}
				arr
			})
		}
		y_rn <- lapply(Y, rownames)
		y_cn <- lapply(Y, colnames)
		y_no <- vector("list", length(Y))
		if (!is.null(Xdyad)) Xdyad <- .inherit_dimnames(Xdyad, y_rn, y_cn)
		if (!is.null(Xrow)) Xrow <- .inherit_dimnames(Xrow, y_rn, y_no)
		if (!is.null(Xcol)) Xcol <- .inherit_dimnames(Xcol, y_cn, y_no)
	}
	# a covariate that is entirely NA in every period can never be identified
	.lame_cov_na_audit <- function(lst, nm, dyadic) {
		if (is.null(lst) || !is.list(lst) || length(lst) == 0L) return(FALSE)
		has_na <- any(vapply(lst, anyNA, logical(1L)))
		p <- max(vapply(lst, function(x) {
			d <- dim(x)
			if (dyadic && length(d) == 3L) d[3L] else if (!dyadic && length(d) == 2L) d[2L] else 1L
		}, numeric(1L)))
		for (j in seq_len(p)) {
			# count only periods whose shape carries slice j (a 2-D dyadic element is a single covariate)
			slice_na <- vapply(lst, function(x) {
				d <- dim(x)
				if (dyadic && length(d) == 3L && d[3L] >= j) return(all(is.na(x[, , j])))
				if (dyadic && length(d) == 2L && j == 1L) return(all(is.na(x)))
				if (!dyadic && length(d) == 2L && d[2L] >= j) return(all(is.na(x[, j])))
				NA
			}, logical(1L))
			slice_na <- slice_na[!is.na(slice_na)]
			if (length(slice_na) > 0L && all(slice_na)) {
				cli::cli_abort(c(
					"{.arg {nm}} covariate {.val {j}} is entirely NA in every period.",
					"i" = "An all-NA covariate cannot be identified; drop it or supply values."
				))
			}
		}
		has_na
	}
	.lame_cov_na_audit(Xdyad, "Xdyad", dyadic = TRUE)
	.lame_cov_na_audit(Xrow, "Xrow", dyadic = FALSE)
	.lame_cov_na_audit(Xcol, "Xcol", dyadic = FALSE)
	# per-period checks on Y: rectangular slices need mode = "bipartite", and duplicated actor
	# names would scramble the name-based alignment below
	for (tt in seq_along(Y)) {
		if (!bip && nrow(Y[[tt]]) != ncol(Y[[tt]])) {
			cli::cli_abort(c(
				"{.arg Y} period {tt} is rectangular ({nrow(Y[[tt]])} x {ncol(Y[[tt]])}) but {.arg mode} = {.val unipartite}.",
				"i" = "Pass {.code mode = \"bipartite\"} for two-mode (rectangular) networks."
			))
		}
		dup <- unique(c(
			rownames(Y[[tt]])[duplicated(rownames(Y[[tt]]))],
			colnames(Y[[tt]])[duplicated(colnames(Y[[tt]]))]
		))
		if (length(dup) > 0L) {
			cli::cli_abort(c(
				"{.arg Y} period {tt} has duplicated actor name{?s}: {.val {head(dup, 5)}}.",
				"i" = "Actors are aligned by name across periods; every actor needs a unique name."
			))
		}
	}
	# actor sets
	if (bip) {
		row_actor_by_yr <- lapply(Y, rownames)
		colActorByYr <- lapply(Y, colnames)
		rowActorSet <- sort(unique(unlist(row_actor_by_yr)))
		colActorSet <- sort(unique(unlist(colActorByYr)))
		nA <- length(rowActorSet)
		nB <- length(colActorSet)
		# overlapping bipartite row and column names would label two different actors with the same key
		overlap <- intersect(rowActorSet, colActorSet)
		if (length(overlap) > 0L) {
			cli::cli_warn(c(
				"{.arg mode} = {.val bipartite} but {length(overlap)} actor label{?s} appear{?s/} on both the row and column sides: {.val {head(overlap, 5)}}{.if {length(overlap) > 5} (and {length(overlap) - 5} more)}.",
				"i" = "In bipartite networks the row and column actor sets should be disjoint; APM and BPM will carry the same label but refer to different actors.",
				"i" = "Rename to disambiguate (e.g. prefix row actors with {.val r_} and column actors with {.val c_}) or refit with {.code mode = \"unipartite\"} if the data is square."
			))
		}
		n <- nA
		actorSet <- rowActorSet
		actorByYr <- row_actor_by_yr
	} else {
		actorByYr <- lapply(Y, rownames)
		actorSet <- sort(unique(unlist(actorByYr)))
		n <- length(actorSet)
		nA <- nB <- n
	}
	.lame_presence_matrix <- function(actor_set, actors_by_period) {
		out <- matrix(
			FALSE,
			length(actor_set),
			length(actors_by_period),
			dimnames = list(actor_set, names(Y))
		)
		for (tt in seq_along(actors_by_period)) {
			out[actors_by_period[[tt]], tt] <- TRUE
		}
		out
	}
	.lame_transition_eligible <- function(presence) {
		out <- matrix(
			FALSE,
			nrow(presence),
			ncol(presence),
			dimnames = dimnames(presence)
		)
		if (ncol(presence) > 1L) {
			out[, -1L] <- presence[, -1L, drop = FALSE] &
				presence[, -ncol(presence), drop = FALSE]
		}
		out
	}
	if (bip) {
		snap_eligible_u <- .lame_transition_eligible(
			.lame_presence_matrix(rowActorSet, row_actor_by_yr)
		)
		snap_eligible_v <- .lame_transition_eligible(
			.lame_presence_matrix(colActorSet, colActorByYr)
		)
	} else {
		snap_presence <- .lame_presence_matrix(actorSet, actorByYr)
		snap_eligible_u <- .lame_transition_eligible(snap_presence)
		snap_eligible_v <- snap_eligible_u
	}

	# odmax defaults to the largest observed outdegree; a scalar is expanded per actor
	if (is.null(odmax)) {
		odmax <- rep(max(unlist(lapply(Y, function(y) { apply(y > 0, 1, sum, na.rm = TRUE) }))), n)
	} else if (length(odmax) == 1L) {
		odmax <- rep(odmax, n)
	}

	# save points
	save_points <- (burn:(nscan + burn))[(burn:(nscan + burn)) %% odens == 0]
	save_points <- save_points[round(quantile(1:length(save_points), probs = seq(save_interval, 1, save_interval)))]

	check_format(Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)

	N <- length(Y) ; pd_labs <- names(Y)
	# time_index and period_exposure need one entry per period
	if (!is.null(time_index) && length(time_index) != N) {
		cli::cli_abort(c(
			"{.arg time_index} must have one entry per period.",
			"i" = "{.arg Y} has {N} period{?s}, but {.arg time_index} has {length(time_index)}."
		))
	}
	if (!is.null(period_exposure) && length(period_exposure) != N) {
		cli::cli_abort(c(
			"{.arg period_exposure} must have one entry per period.",
			"i" = "{.arg Y} has {N} period{?s}, but {.arg period_exposure} has {length(period_exposure)}."
		))
	}
	if (!is.null(period_exposure)) for (t in which(period_exposure == 0)) {
		observed <- Y[[t]]
		if (!bip) diag(observed) <- NA_real_
		if (any(observed != 0, na.rm = TRUE))
			cli::cli_abort(c(
				"An observed positive count requires strictly positive {.arg period_exposure}.",
				"i" = "Zero exposure is valid only when every observed count in that period is zero."
			))
	}
	.lame_sym_y <- FALSE
	if (!bip) {
		Y <- lapply(Y, function(y) {diag(y) <- NA; return(y)})
		# symmetric data under the directed model: with a dyadic correlation its posterior concentrates near one
		# (refused for the normal family once the design is built)
		if (!isTRUE(symmetric)) {
			# per period: NA when no dyad has both cells observed, else whether every observed pair agrees
			sym_stat <- vapply(Y, function(y) {
				off <- y - t(y)
				f <- is.finite(off)
				if (!any(f)) NA else max(abs(off[f])) <= 1e-8
			}, logical(1))
			if (any(!is.na(sym_stat)) && all(sym_stat[!is.na(sym_stat)])) {
				.lame_sym_y <- TRUE
			}
		}
	}

	# convert into large array format
	if (bip) {
		arrayObj <- list_to_array_bipartite(rowActorSet, colActorSet, Y, Xdyad, Xrow, Xcol)
	} else {
		arrayObj <- list_to_array(actorSet, Y, Xdyad, Xrow, Xcol)
	}
	Y <- arrayObj$Y ; Xdyad <- arrayObj$Xdyad ; Xrow <- arrayObj$Xrow
	Xcol <- arrayObj$Xcol ; rm(arrayObj)
	if (bip && is.list(Xdyad) && !is.null(Xdyad) && length(Xdyad) > 0) {
		pd_tmp <- dim(Xdyad[[1]])[3]
		# keep the per-slice covariate names of the first period
		slice_names <- dimnames(Xdyad[[1]])[[3]]
		Xdyad_4d <- array(
			NA,
			dim = c(nA, nB, pd_tmp, length(Xdyad)),
			dimnames = list(NULL, NULL, slice_names, NULL)
		)
		for (tt in seq_along(Xdyad)) {
			Xdyad_4d[, , , tt] <- Xdyad[[tt]]
		}
		Xdyad <- Xdyad_4d
		rm(Xdyad_4d)
	}
	if (bip) {
		# a cell whose covariates are missing (an actor absent from a period) is treated as unobserved
		# with its covariates set to zero, as in the unipartite design builder
		n_na_before <- sum(is.na(Y))
		if (is.array(Xdyad) && length(dim(Xdyad)) == 4L && anyNA(Xdyad)) {
			for (tt in seq_len(dim(Xdyad)[4L])) {
				Xt <- Xdyad[, , , tt, drop = FALSE]
				na_cov <- apply(is.na(Xt), c(1L, 2L), any)
				if (any(na_cov)) {
					Yt <- Y[, , tt]
					Yt[na_cov] <- NA
					Y[, , tt] <- Yt
					Xt[is.na(Xt)] <- 0
					Xdyad[, , , tt] <- Xt
				}
			}
		}
		if (is.array(Xrow) && length(dim(Xrow)) == 3L && anyNA(Xrow)) {
			for (tt in seq_len(dim(Xrow)[3L])) {
				na_row <- apply(is.na(Xrow[, , tt, drop = FALSE]), 1L, any)
				if (any(na_row)) Y[na_row, , tt] <- NA
			}
			Xrow[is.na(Xrow)] <- 0
		}
		if (is.array(Xcol) && length(dim(Xcol)) == 3L && anyNA(Xcol)) {
			for (tt in seq_len(dim(Xcol)[3L])) {
				na_col <- apply(is.na(Xcol[, , tt, drop = FALSE]), 1L, any)
				if (any(na_col)) Y[, na_col, tt] <- NA
			}
			Xcol[is.na(Xcol)] <- 0
		}
		.lame_inform_cov_na(sum(is.na(Y)) - n_na_before)
	}
	# g-prior parameter when not provided; `g` may be a scalar or a per-coefficient vector
	if (length(g) == 1L && is.na(g)) {
		# Zellner's g-prior default: for the static unipartite prior g carries the unit of the response
		# (n^2 times its variance, as in ame()); bipartite and dynamic fits use g = n^2 since their prior is N(0, g s2)
		g <- n^2 * (if (isTRUE(bip) || dynamic_beta_requested_global || !is.null(dynamic_beta_per_actor)) 1 else .lame_unit_var(Y, family))
	}
	# the unit of the response: 1 where the link fixes the latent scale, the response variance for the normal family.
	# starting values, several prior variances and the rails of the dynamic priors are written in it
	unit_var <- .lame_unit_var(Y, family)
	unit_sd <- sqrt(unit_var)
	weight_prior_var <- .lame_weight_prior_var(Y, family)
	rail_units <- .lame_rail_units(Y, family, bip, symmetric)
	# the interaction matrix is cut at 100 units of the response unless prior$G_max is set
	g_rail_default <- 100 * rail_units$G
	# prior variance of the initial dynamic coefficients: 10 in the unit of the response
	beta0_var_default <- .lame_beta0_var_default(Y, family)

	if (is.element(family, c("binary", "cbin"))) { Y <- 1 * (Y > 0) }

	ab_obs_mask <- is.finite(Y)
	.lame_row_obs_ab <- apply(ab_obs_mask, 1L, sum)
	.lame_col_obs_ab <- apply(ab_obs_mask, 2L, sum)
	# observed and max outdegrees
	if (is.element(family, c("cbin", "frn"))) {
		odobs <- apply(Y > 0, c(1, 3), sum, na.rm = TRUE)
		# the posterior-predictive template per period, built on first use
		yo_pred <- NULL
		# hidden cells hold a nomination when their latent value is positive, so a row's nomination count is
		# recomputed every sweep; the diagonal and an absent actor's cells are never nominations
		nom_hid_cells <- array(FALSE, dim(Y))
		for (tt in seq_len(dim(Y)[3L])) {
			nom_hid_cells[, , tt] <- .lame_nom_hidden_cells(Y[, , tt], bip)
		}
		.odobs_room <- function(t) {
			.lame_nom_room(odobs[, t], nom_hid_cells[, , t], Z[, , t])
		}
		# reject odmax = 0 with observed nominations, and odmax = NA
		if (!is.null(odmax)) {
			if (anyNA(odmax)) {
				cli::cli_abort(c(
					"{.arg odmax} contains {.val NA}.",
					"i" = "Supply a positive integer (scalar or per-row vector). Pass {.code odmax = NULL} to let {.fn lame} infer it from {.arg Y}."
				))
			}
			if (any(odmax <= 0L) && any(Y > 0, na.rm = TRUE)) {
				cli::cli_abort(c(
					"{.arg odmax} <= 0 but {.arg Y} contains nominations.",
					"i" = "{.code odmax = 0} prohibits any nominations; this is inconsistent with the observed data."
				))
			}
		}
		if (!is.null(odmax) && !(length(odmax) %in% c(1L, nrow(Y[, , 1])))) {
			cli::cli_abort(c(
				"{.arg odmax} must be a single value or one value per row.",
				"i" = "Got length {length(odmax)} for {nrow(Y[,,1])} rows."
			))
		}
		if (length(odmax) == 1) { odmax <- rep(odmax, nrow(Y[, , 1])) }
		# a cap below a row's observed outdegree in any period contradicts Y
		.lame_check_odmax_outdegree(odobs, odmax, dimnames(Y)[[1L]])
	}

	# symmetric = TRUE needs a symmetric Y: the sampler averages the two triangles
	if (symmetric && !bip) {
		off_max <- 0
		one_sided <- FALSE
		for (tt in seq_len(dim(Y)[3])) {
			Yt <- Y[, , tt]
			off <- Yt - t(Yt)
			if (any(is.finite(off))) off_max <- max(off_max, max(abs(off), na.rm = TRUE))
			ut <- upper.tri(Yt) ; lt <- lower.tri(Yt)
			if ((sum(is.finite(Yt[ut])) == 0L) != (sum(is.finite(Yt[lt])) == 0L)) {
				one_sided <- TRUE
			}
		}
		if (off_max > 1e-8) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric (max |Y - t(Y)| = {signif(off_max, 4)}).",
				"i" = "Symmetrize every slice (e.g. {.code lapply(Y, function(y) (y + t(y))/2)}) or set {.arg symmetric} = FALSE."
			))
		}
		if (one_sided) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but only one triangle of {.arg Y} is observed (the other is entirely {.val NA}).",
				"i" = "The symmetric sampler would silently mirror the populated triangle.",
				"i" = "Symmetrize explicitly with {.code y[lower.tri(y)] <- t(y)[lower.tri(y)]} (or vice versa) before fitting."
			))
		}
		# a symmetric dyad has one covariate value: holes in one triangle are filled from the other,
		# and the regression term must be symmetric (as in ame())
		if (!is.null(Xdyad)) {
			Xdyad <- .lame_fill_symmetric_dyadic(Xdyad)
			.lame_check_symmetric_xdyad(Xdyad)
		}
		Xcol <- Xrow ; rvar <- cvar <- nvar
	}
	if (bip && symmetric) {
		cli::cli_warn("`symmetric` ignored for bipartite networks")
		symmetric <- FALSE
	}
	# validate latent ranks as in ame(); bipartite ranks are each bounded by their own mode,
	# and a shared R by the smaller one
	rmax_of <- c(
		R = if (bip) min(nA, nB) else n,
		R_row = if (bip) nA else n,
		R_col = if (bip) nB else n
	)
	for (rr_nm in c("R", "R_row", "R_col")) {
		rv <- get(rr_nm)
		if (!is.null(rv)) {
			if (length(rv) != 1L || anyNA(rv) || rv < 0 ||
					!isTRUE(all.equal(as.numeric(rv), round(as.numeric(rv))))) {
				cli::cli_abort("{.arg {rr_nm}} must be a single non-negative integer.")
			}
			rmax <- rmax_of[[rr_nm]]
			if (rv >= rmax) {
				cli::cli_abort("{.arg {rr_nm}} = {rv} must be smaller than the network dimension ({rmax}).")
			}
		}
	}
	# rank dimensions
	if (bip) {
		RA <- R_row %||% R %||% 0L
		RB <- R_col %||% R %||% 0L
		dcor <- FALSE # no dyadic correlation in bipartite
		# the bilinear term U G V' needs both ranks positive
		if (xor(RA > 0L, RB > 0L)) {
			cli::cli_warn(c(
				"!" = "Only one of {.arg R_row} ({RA}) / {.arg R_col} ({RB}) is positive.",
				"i" = "The multiplicative term U G V' needs both ranks > 0; it is dropped from this fit.",
				"i" = "Set both {.arg R_row} and {.arg R_col} (or pass {.arg R}) to include latent factors."
			))
		}
	} else {
		if (!is.null(R_row) || !is.null(R_col)) {
			cli::cli_warn("R_row/R_col ignored for unipartite networks")
		}
		RA <- RB <- as.integer(R %||% 0L)
	}
	.lame_has_uv_rank <- function(bip, R, RA, RB) {
		if (isTRUE(bip)) return(isTRUE(RA > 0L) && isTRUE(RB > 0L))
		isTRUE(R > 0L)
	}
	if (isTRUE(dynamic_uv) && !.lame_has_uv_rank(bip, R, RA, RB)) {
		cli::cli_abort(c(
			"{.arg dynamic_uv = TRUE} requires a positive multiplicative latent rank.",
			"i" = if (bip) {
				"Set positive {.arg R_row} and {.arg R_col} values, or set {.arg R} > 0 for both."
			} else {
				"Set {.arg R} > 0 or disable {.arg dynamic_uv}."
			}
		))
	}
	if (isTRUE(dynamic_uv) && isTRUE(bip) &&
			dynamic_uv_kind %in% c("snap", "t")) {
		cli::cli_abort(c(
			"{.arg dynamic_uv_kind = \"{dynamic_uv_kind}\"} is not available for this bipartite MCMC path.",
			"i" = "Use {.arg dynamic_uv_kind = \"ar1\"}, here or with {.code method = \"als\"}; the ALS route refuses bipartite snap and t as well."
		))
	}
	# the symmetric product's unit-variance coordinates do not hold their scale under the reset
	# and scale-mixture transitions
	if (isTRUE(dynamic_uv) && isTRUE(symmetric) && !isTRUE(bip) &&
			dynamic_uv_kind %in% c("snap", "t")) {
		cli::cli_abort(c(
			"{.arg dynamic_uv_kind = \"{dynamic_uv_kind}\"} is not available for symmetric fits on the MCMC path.",
			"i" = "Use {.arg dynamic_uv_kind = \"ar1\"}, which also follows paths with resets or heavy-tailed steps."
		))
	}
	if (dynamic_uv && verbose) {
		msg <- switch(
			dynamic_uv_kind,
			ar1 = c(
				"i" = "Fitting AR(1) time-varying multiplicative latent effects.",
				"i" = "Latent coordinates are gauge-dependent; interpret aligned trajectories, bilinear contributions, and predictive diagnostics."
			),
			snap = c(
				"i" = "Fitting snap-shift dynamic multiplicative latent effects.",
				"i" = "First-period snap probabilities are transition-ineligible and reported as NA."
			),
			t = c(
				"i" = "Fitting heavy-tailed Student-t dynamic multiplicative latent effects.",
				"i" = "Local transition scales diagnose unusually large movement; they are not discrete snap probabilities."
			)
		)
		cli::cli_inform(msg)
	}
	# design matrix
	pr <- length(Xrow[, , 1]) / n
	pc <- length(Xcol[, , 1]) / n

	# default Xlist for bipartite before building the design
	if (bip) {
		pd <- if (!is.null(Xdyad)) length(Xdyad[, , , 1]) / (nA * nB) else 0
		if (!is.null(Xdyad) && length(Xdyad) > 0) {
			# convert Xdyad (nA x nB x pd x N) to a list of 3d arrays, keeping the dyadic slice names
			Xlist <- vector("list", N)
			n_cov_d <- dim(Xdyad)[3]
			n_cov_total <- n_cov_d
			if (intercept) n_cov_total <- n_cov_total + 1L
			dyad_names <- dimnames(Xdyad)[[3]]
			if (is.null(dyad_names) || length(dyad_names) == 0L) {
				dyad_names <- paste0("x", seq_len(n_cov_d))
			}
			slice_names <- if (intercept) c("intercept", dyad_names) else dyad_names
			for (t in 1:N) {
				Xt <- array(
					0,
					dim = c(nA, nB, n_cov_total),
					dimnames = list(NULL, NULL, slice_names)
				)
				k_start <- 1L
				if (intercept) {
					Xt[, , 1] <- 1
					k_start <- 2L
				}
				for (k in 1:n_cov_d) {
					Xt[, , k_start + k - 1] <- Xdyad[, , k, t]
				}
				Xlist[[t]] <- Xt
			}
		} else {
			# intercept-only Xlist
			if (intercept) {
				Xlist <- replicate(
					N,
					array(
						1,
						dim = c(nA, nB, 1L),
						dimnames = list(NULL, NULL, "intercept")
					),
					simplify = FALSE
				)
			} else {
				Xlist <- replicate(N, array(0, dim = c(nA, nB, 0L)), simplify = FALSE)
			}
		}
		# add row/column covariates to Xlist
		if (!is.null(Xrow) && length(Xrow) > 0) {
			n_xr <- if (length(dim(Xrow)) >= 2) dim(Xrow)[2] else 1L
			# carry the user's column names and append `_row` only when needed
			xr_names <- if (length(dim(Xrow)) >= 2) dimnames(Xrow)[[2]] else NULL
			if (is.null(xr_names) || length(xr_names) == 0L)
				xr_names <- paste0("xrow", seq_len(n_xr))
			xr_names <- .lame_apply_suffix(xr_names, "row")
			for (t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				old_nm <- dimnames(Xlist[[t]])[[3]]
				new_Xt <- array(
					0,
					dim = c(nA, nB, old_p + n_xr),
					dimnames = list(NULL, NULL, c(old_nm, xr_names))
				)
				new_Xt[, , 1:old_p] <- Xlist[[t]]
				for (k in 1:n_xr) {
					xr_vals <- if (length(dim(Xrow)) == 3) Xrow[, k, t] else Xrow[, t]
					new_Xt[, , old_p + k] <- matrix(xr_vals, nA, nB)
				}
				Xlist[[t]] <- new_Xt
			}
		}
		if (!is.null(Xcol) && length(Xcol) > 0) {
			n_xc <- if (length(dim(Xcol)) >= 2) dim(Xcol)[2] else 1L
			xc_names <- if (length(dim(Xcol)) >= 2) dimnames(Xcol)[[2]] else NULL
			if (is.null(xc_names) || length(xc_names) == 0L)
				xc_names <- paste0("xcol", seq_len(n_xc))
			xc_names <- .lame_apply_suffix(xc_names, "col")
			for (t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				old_nm <- dimnames(Xlist[[t]])[[3]]
				new_Xt <- array(
					0,
					dim = c(nA, nB, old_p + n_xc),
					dimnames = list(NULL, NULL, c(old_nm, xc_names))
				)
				new_Xt[, , 1:old_p] <- Xlist[[t]]
				for (k in 1:n_xc) {
					xc_vals <- if (length(dim(Xcol)) == 3) Xcol[, k, t] else Xcol[, t]
					new_Xt[, , old_p + k] <- t(matrix(xc_vals, nB, nA))
				}
				Xlist[[t]] <- new_Xt
			}
		}
	} else {
		# unipartite networks
		pd <- length(Xdyad[, , , 1]) / n^2
	}
	# use the whole observed panel to decide whether the intercept is redundant
	if (bip && isTRUE(intercept) && length(Xlist) > 0L && dim(Xlist[[1L]])[3L] > 1L) {
		covariates <- lapply(Xlist, function(Xt) Xt[, , -1L, drop = FALSE])
		if (.lame_panel_intercept_redundant(covariates, Y)) {
			cli::cli_warn(c(
				"The covariates are collinear with the intercept across the observed panel.",
				"i" = "The intercept is dropped; the covariates carry the level."
			))
			Xlist <- covariates
			intercept <- FALSE
		}
	}

	if (!bip) {
		# get_design_rep is unipartite only; a cell with a missing covariate comes back unobserved in Y
		n_na_before <- sum(is.na(Y))
		design_obj <- get_design_rep(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			actorSet = actorSet,
			intercept = intercept,
			n = n,
			N = N,
			pr = pr,
			pc = pc,
			pd = pd,
			symmetric = symmetric
		)
		Y <- design_obj$Y ; X <- design_obj$X ; Xlist <- design_obj$Xlist
		.lame_inform_cov_na(sum(is.na(Y)) - n_na_before)
		XrLong <- design_obj$XrLong ; XcLong <- design_obj$XcLong
		mXLong <- design_obj$mXLong ; mXtLong <- design_obj$mXtLong
		xxLong <- design_obj$xxLong ; xxTLong <- design_obj$xxTLong
		# adopt the panel-wide flag so a/b centring shifts the intercept only when the design contains it
		intercept <- design_obj$intercept ; rm(design_obj)
		.x_slice_sd <- .lame_check_covariate_scale(X)
	} else {
		# bipartite: Xlist was built above
		X <- Xdyad
		.lame_check_covariate_scale(Xlist)
		XrLong <- NULL ; XcLong <- NULL
		mXLong <- NULL ; mXtLong <- NULL
		xxLong <- NULL ; xxTLong <- NULL
	}
	# a design aliased over the observed panel has coefficients the data cannot separate
	.lame_check_design_rank(Xlist, Y)
	# a directed normal model with a dyadic correlation has an improper posterior when every reciprocal
	# pair is exactly equal and the pairs outnumber what the mean parameters absorb: check the pooled
	# bound (.lame_panel_reciprocity_bound) and the per-period ame() check
	if (.lame_sym_y && identical(family, "normal") && !bip && !isTRUE(symmetric) && isTRUE(dcor)) {
		for (tt in seq_len(N)) {
			tryCatch(
				.ame_check_normal_reciprocity(Y[, , tt], Xlist[[tt]], rvar, cvar, dcor, symmetric, RA, family),
				lame_improper_normal_reciprocity = function(e) {
					cli::cli_abort(
						c(
							"Every period of {.arg Y} is exactly symmetric, but {.arg symmetric} = FALSE with {.arg dcor}.",
							"i" = "The two cells of each dyad are equal, and the reciprocal pairs of period {tt} outnumber what the mean parameters absorb, so the posterior of the directed normal model is improper."
						),
						class = c("lame_symmetric_y", "lame_improper_normal_reciprocity"),
						parent = e
					)
				}
			)
		}
		if (!isTRUE(dynamic_rho) && is.null(dynamic_beta_per_actor)) {
			bound <- .lame_panel_reciprocity_bound(
				Y,
				Xlist,
				rvar,
				cvar,
				RA,
				dynamic_ab = isTRUE(dynamic_ab),
				dynamic_uv = isTRUE(dynamic_uv),
				dynamic_beta = dynamic_beta_requested_global
			)
			if (bound$improper) {
				need <- .lame_rho_improper_dyads()
				cli::cli_abort(
					c(
						"Every period of {.arg Y} is exactly symmetric, but {.arg symmetric} = FALSE with {.arg dcor}.",
						"i" = "The two cells of each dyad are equal: all {bound$pairs} reciprocal pairs of the panel are exactly equal, and the mean parameters absorb at most {bound$rank} of their contrasts ({bound$actor} for the actor effects, {bound$coef} for the coefficients, {bound$factor} for the latent factors). The posterior of the directed normal model is improper once the pairs exceed that bound by {need}.",
						"i" = "Use {.code symmetric = TRUE} for an undirected network, or {.code dcor = FALSE} for fixed zero reciprocity."
					),
					class = c("lame_symmetric_y", "lame_improper_normal_reciprocity")
				)
			}
		}
	}
	# an exactly symmetric Y under the directed model that passes the bound: say what it is
	if (.lame_sym_y) .lame_inform_symmetric_y(dcor, panel = TRUE, family = family)

	# warn about a constant covariate under data-induced ordinal cutpoints
	if (!bip && is.element(family, c("ordinal")) && identical(ordinal_cutpoints, "data_induced") &&
			!is.null(X) && length(dim(X)) == 4 &&
			any(apply(X, c(3), function(x) {var(c(x))}) == 0)) {
		cli::cli_warn("An intercept is not estimable using this procedure")
	}

	# ranked nominations for frn: the sampler renumbers ranks to consecutive integers within each row;
	# the fit stores the user's ranks
	if (is.element(family, c("frn"))) {
		Y_user <- Y
		ymx <- max(apply(1 * (Y > 0), c(1, 3), sum, na.rm = TRUE))
		YL <- list()
		for (t in 1:N) {
			YL.t <- NULL
			warn <- FALSE
			for (i in 1:nrow(Y[, , 1])) {
				yi <- Y[i, , t] ; rnkd <- which(!is.na(yi) & yi > 0)
				if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {warn <- TRUE}
				yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
				Y[i, , t] <- yi
				YL.t <- rbind(YL.t, match(1:ymx, yi))
			}
			YL[[t]] <- YL.t
			if (warn) {cli::cli_warn("Random reordering used to break ties in ranks")}
		}
	}

	# the assembled bipartite design includes the intercept
	if (bip && !is.null(Xlist) && length(Xlist) > 0) {
		p <- dim(Xlist[[1]])[3]
	} else if (!is.null(X) && length(dim(X)) >= 3) {
		p <- dim(X)[3]
	} else {
		p <- if (intercept) 1L else 0L
	}
	# number of observed dyads, as in amen; `g` may also be a per-coefficient vector
	if (length(g) == 1L && is.na(g)) {
		n_dyads <- sum(!is.na(Y))
		g <- max(n_dyads, n^2) # use number of observed dyads or n^2, whichever is larger
	}
	if (length(g) == 1) { g <- rep(g, p) }

	# prior list
	if (!is.list(prior)) {
		if (!is.null(prior)) {
			cli::cli_warn(c(
				"{.arg prior} must be a named list; got {.cls {class(prior)[1]}}.",
				"i" = "Ignoring it and fitting with the default priors."
			))
		}
		prior <- list()
	}
	s2_prior_scale <- .lame_residual_prior_scale(Y, family, prior$s20)
	if (is.null(prior$Sab0)) {
		if (family %in% c("normal", "poisson")) {
			# scale the additive-effects prior to the data, as in ame_unipartite.R
			if (family == "poisson") {
				# poisson uses log-rate moments (.ame_pois_effect_scale)
				vscale_ab <- .ame_pois_effect_scale(Y)
			} else {
				Ec <- Y - mean(Y, na.rm = TRUE)
				a_tmp <- apply(Ec, 1, mean, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
				b_tmp <- apply(Ec, 2, mean, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
				# variance per side (a and b can differ in length in bipartite)
				vscale_ab <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
				if (!is.finite(vscale_ab) || vscale_ab <= 0) {
					vscale_ab <- max(stats::var(c(Y), na.rm = TRUE), 1e-12)
				}
			}
			prior$Sab0 <- diag(2) * vscale_ab
		} else if (family == "binary") {
			# the probit-moment adaptive scale of ame_unipartite, so ame() and lame() share the prior on binary data
			Y1 <- if (is.array(Y)) Y[, , 1] else Y
			ydist <- table(Y1)
			ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
			YB <- 1 * (Y1 != ymode)
			ybar <- mean(YB, na.rm = TRUE)
			if (is.finite(ybar) && ybar > 0 && ybar < 1) {
				E <- (YB - ybar) / dnorm(qnorm(ybar))
				if (nrow(E) == ncol(E)) diag(E) <- 0
				a_tmp <- rowMeans(E, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
				b_tmp <- colMeans(E, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
				vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)))
				if (!is.finite(vscale) || vscale <= 0) vscale <- 1
				prior$Sab0 <- diag(2) * vscale
			} else {
				prior$Sab0 <- diag(2)
			}
		} else if (family %in% c("ordinal", "cbin", "frn")) {
			# rank-based families take the probit-moment scale (.lame_probit_ab_scale); cbin and frn
			# read it off the rows below their cap, one scale per side
			sc_ab <- .lame_probit_ab_scale(Y, family, odmax)
			prior$Sab0 <- if (family %in% c("cbin", "frn")) {
				diag(c(sc_ab$vscale_row, sc_ab$vscale_col))
			} else diag(2) * sc_ab$vscale
		} else {
			prior$Sab0 <- diag(2)
		}
	}
	# binary eta0
	if (family == "binary" && is.null(prior$eta0)) {
		Y1 <- if (is.array(Y)) Y[, , 1] else Y
		ydist <- table(Y1)
		ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
		YB <- 1 * (Y1 != ymode)
		ybar <- mean(YB, na.rm = TRUE)
		if (ybar > 0 && ybar < 1) {
			mu_bin <- qnorm(ybar)
			E <- (YB - ybar) / dnorm(qnorm(ybar))
			if (!bip && nrow(E) == ncol(E)) diag(E) <- 0
			a_tmp <- rowMeans(E, na.rm = TRUE)
			b_tmp <- colMeans(E, na.rm = TRUE)
			a_tmp[is.na(a_tmp)] <- 0
			b_tmp[is.na(b_tmp)] <- 0
			if (bip) {
				vscale <- (var(a_tmp) + var(b_tmp)) / 2
			} else {
				vscale <- mean(diag(cov(cbind(a_tmp, b_tmp))))
			}
			PHAT <- pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
			vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
			prior$eta0 <- round(4 * vdfmlt)
		} else {
			prior$eta0 <- round(4 + 3 * n / 100)
		}
	} else if (is.null(prior$eta0) && family %in% c("ordinal", "cbin", "frn")) {
		# prior degrees of freedom inflated by the inverse Bernoulli information, as for the binary family
		prior$eta0 <- round(4 * .lame_probit_ab_scale(Y, family, odmax)$vdfmlt)
	} else if (is.null(prior$eta0)) {
		prior$eta0 <- round(4 + 3 * n / 100)
	}
	if (is.null(prior$etaab)) { prior$etaab <- prior$eta0 }
	if (is.null(prior$ab_min_observed)) { prior$ab_min_observed <- max(10, N) }
	# prior of the static multiplicative effects of a unipartite network, shared with ame() (R/uv_prior.R)
	uv_prior_static <- if (!bip && !symmetric && R > 0L) .lame_uv_prior(Y, family, R, odmax)
	# the same prior for the position moves, which apply it from the second half of burn-in
	uv_prior_move <- if (!is.null(uv_prior_static)) .lame_uv_prior_move(uv_prior_static, R)

	# dynamic UV needs the relevant multiplicative ranks
	.dyn_uv_has_rank <- isTRUE(dynamic_uv) && .lame_has_uv_rank(bip, R, RA, RB)
	# at rank two or more the snap sampler, and a t sampler without an ALS start, do not recover their
	# position hyperparameters (see .lame_uv_hyper_unidentified)
	if (.dyn_uv_has_rank && isTRUE(R >= 2) && (identical(dynamic_uv_kind, "snap") ||
			(identical(dynamic_uv_kind, "t") && !isTRUE(als_start_used)))) {
		cli::cli_inform(c(
			"i" = if (identical(dynamic_uv_kind, "snap")) {
				"With {.code dynamic_uv_kind = \"snap\"} at rank two or more, {.fn summary} and {.fn print} do not report {.code rho_uv}, {.code sigma_uv} or the reset rate {.code pi_snap}: extra resets stand in for drift, so the data do not separate the reset rate from the drift."
			} else {
				"This {.code dynamic_uv_kind = \"t\"} chain could not start from a dynamic ALS fit, so at rank two or more {.fn summary} and {.fn print} do not report {.code rho_uv} or {.code sigma_uv}: started cold, the chain lets heavy-tailed steps stand in for drift."
			},
			"i" = "The positions and the stationary spread of their paths were recovered.",
			"i" = if (identical(dynamic_uv_kind, "snap")) {
				"Read {.code fit$UVPM} and the positions; the draws stay in {.code fit$rho_uv}, {.code fit$sigma_uv} and {.code fit$pi_snap}."
			} else {
				"Read {.code fit$UVPM} and the positions; the draws stay in {.code fit$rho_uv} and {.code fit$sigma_uv}."
			}
		))
	}
	# dynamic UV priors
	if (.dyn_uv_has_rank) {
		gaussian_bipartite_ar <- isTRUE(bip) && identical(family, "normal") &&
			identical(dynamic_uv_kind, "ar1")
		if (is.null(prior$rho_uv_mean)) {
			prior$rho_uv_mean <- if (gaussian_bipartite_ar) .5 else .9
		}
		if (is.null(prior$rho_uv_sd)) {
			prior$rho_uv_sd <- if (gaussian_bipartite_ar) .2 else .1
		}
		if (is.null(prior$sigma_uv_shape)) { prior$sigma_uv_shape <- 2 }
		# scale of the inverse-gamma prior on the position innovation variance (see .lame_sigma_uv_scale_default)
		if (is.null(prior$sigma_uv_scale)) { prior$sigma_uv_scale <- .lame_sigma_uv_scale_default(Y, family, bip, symmetric) }
		# a saturating likelihood can push the product scale up without limit once the fit separates the data,
		# so those families cap the stationary variance of the coordinate paths; the probit-scale families also
		# take a tighter coordinate rail and a bounded interaction matrix
		.saturating <- !identical(family, "normal")
		.probit_scale <- family %in% c("binary", "ordinal", "cbin", "frn")
		if (is.null(prior$uv_var_cap)) { prior$uv_var_cap <- if (.saturating) 4 else Inf }
		# the rails are in the unit of the response: unipartite coordinates carry its square root, G the unit
		if (is.null(prior$uv_max_abs)) { prior$uv_max_abs <- (if (.probit_scale) 8 else 50) * rail_units$uv }
		if (is.null(prior$G_max)) { prior$G_max <- (if (.probit_scale) 20 else 100) * rail_units$G }
		if (is.null(prior$sigma_uv_max)) { prior$sigma_uv_max <- if (bip) 10 else Inf }
		# snap jump scale and beta prior on the snap rate
		if (identical(dynamic_uv_kind, "snap")) {
			# a NULL snap_kappa is the linked scale (.lame_snap_kappa)
			if (!is.null(prior$snap_kappa) &&
					(!is.numeric(prior$snap_kappa) || length(prior$snap_kappa) != 1L ||
							!is.finite(prior$snap_kappa) || prior$snap_kappa <= 0)) {
				cli::cli_abort(
					c(
						"{.code prior$snap_kappa} must be a single positive number.",
						"i" = "Leave it unset to draw a reset from the stationary distribution of the position paths."
					),
					class = "lame_snap_kappa"
				)
			}
			if (is.null(prior$snap_pi_a)) { prior$snap_pi_a <- 1 }
			if (is.null(prior$snap_pi_b)) { prior$snap_pi_b <- 9 }
		}
		if (identical(dynamic_uv_kind, "t")) {
			if (is.null(prior$t_nu)) { prior$t_nu <- 4 }
		}
	}

	# dynamic ab priors
	if (dynamic_ab) {
		if (is.null(prior$rho_ab_mean)) { prior$rho_ab_mean <- 0.8 }
		if (is.null(prior$rho_ab_sd)) { prior$rho_ab_sd <- 0.15 }
		if (is.null(prior$sigma_ab_shape)) { prior$sigma_ab_shape <- 2 }
		# fixed once from the prior inputs so the sigma_ab and Sab draws share one joint density
		# (.lame_sigma_ab_prior_scale); reported as fit$sigma_ab_prior_scale
		if (is.null(prior$sigma_ab_scale)) {
			prior$sigma_ab_scale <- .lame_sigma_ab_prior_scale(
				prior$Sab0,
				prior$etaab,
				prior$rho_ab_mean,
				rvar,
				cvar
			)
		}
		# optional bound on the stationary variance of an additive-effect path, as a multiple of the
		# observation variance; off by default, since the latent scale is held by the coordinate cap and G rail
		if (is.null(prior$ab_var_cap)) { prior$ab_var_cap <- Inf }
	}
	# dynamic beta priors, mirroring the dynamic_ab block
	if (is.null(prior$rho_beta_mean)) { prior$rho_beta_mean <- 0.8 }
	if (is.null(prior$rho_beta_sd)) { prior$rho_beta_sd <- 0.15 }
	if (is.null(prior$sigma_beta_shape)) { prior$sigma_beta_shape <- 2 }
	# fixed once from the response scale (.lame_sigma_beta_prior_scale); reported as fit$sigma_beta_prior_scale
	if (is.null(prior$sigma_beta_scale)) {
		prior$sigma_beta_scale <- if (isTRUE(dynamic_beta_requested_global)) {
			.lame_sigma_beta_prior_scale(Y, family, prior$sigma_beta_shape)
		} else 1
	}
	if (is.null(prior$sigma_beta_init)) { prior$sigma_beta_init <- 0.25 * unit_sd }
	if (is.null(prior$rho_beta_lower)) { prior$rho_beta_lower <- 0 }
	if (is.null(prior$rho_beta_upper)) { prior$rho_beta_upper <- 0.999}

	if (!bip && (is.null(Xlist) || length(Xlist) == 0L)) {
		# unipartite: n x n intercept matrices if still needed
		Xlist <- replicate(N, array(1, dim = c(n, n, 1L)), simplify = FALSE)
	}

	# starting values for MCMC (without the labels of their source, see .start_vals_plain)
	start_vals <- .start_vals_plain(start_vals)
	if (bip) {
		# bipartite initialization
		start_vals_obj <- init_bipartite_startvals(
			Y,
			family,
			nA,
			nB,
			RA,
			RB,
			N,
			Xlist,
			odmax,
			unit_var = unit_var
		)
		# user-supplied start values replace the initialiser's defaults slot by slot
		start_vals_obj <- .lame_bip_overlay_start(
			start_vals_obj,
			start_vals,
			nA,
			nB,
			RA,
			RB,
			N,
			dynamic_ab
		)
		Z <- start_vals_obj$Z
		beta <- start_vals_obj$beta
		a <- start_vals_obj$a
		b <- start_vals_obj$b
		U_full <- start_vals_obj$U
		V_full <- start_vals_obj[["V"]]
		# static case: first time point
		if (!dynamic_uv) {
			if (RA > 0) {
				U <- U_full[, , 1, drop = FALSE] ; dim(U) <- c(nA, RA)
			} else {
				U <- matrix(0, nA, 1) # minimal dimension
			}
			if (RB > 0) {
				V <- V_full[, , 1, drop = FALSE] ; dim(V) <- c(nB, RB)
			} else {
				V <- matrix(0, nB, 1) # minimal dimension
			}
		} else {
			U <- U_full
			V <- V_full
		}
		G <- start_vals_obj$G # rectangular interaction matrix
		G_expected_dim <- c(max(1L, RA), max(1L, RB))
		if (is.null(G) || length(G) == 0L ||
				is.null(dim(G)) ||
				!all(dim(G) == G_expected_dim) || !all(is.finite(G))) {
			G <- if (RA > 0 && RB > 0) {
				diag(min(RA, RB), RA, RB)
			} else {
				matrix(0, max(1L, RA), max(1L, RB))
			}
		}
		sigma2_a <- start_vals_obj$sigma2_a
		sigma2_b <- start_vals_obj$sigma2_b
		s2 <- start_vals_obj$s2
		rho <- 0 # no dyadic correlation in bipartite (0, not NULL)
		Sab <- diag(c(sigma2_a, sigma2_b)) # independent variances
	} else {
		# unipartite initialization
		start_vals_obj <- get_start_vals(
			start_vals,
			Y,
			family,
			xP = dim(X)[3],
			rvar,
			cvar,
			R,
			odmax = odmax,
			bip = bip,
			unit_var = unit_var
		)
		Z <- start_vals_obj$Z ; beta <- start_vals_obj$beta ; a <- start_vals_obj$a
		b <- start_vals_obj$b ; U <- start_vals_obj$U ; V <- start_vals_obj[["V"]]
		rho <- start_vals_obj$rho ; s2 <- start_vals_obj$s2 ; Sab <- start_vals_obj$Sab
		G <- NULL # no G matrix for unipartite
		# slopes start on the latent scale, so a covariate in large units does not start the predictor
		# far from the data
		if (is.null(start_vals[["beta"]]) && length(beta) == length(.x_slice_sd)) {
			.x_big <- is.finite(.x_slice_sd) & .x_slice_sd > 0
			beta[.x_big] <- beta[.x_big] / .x_slice_sd[.x_big]
		}
		# static symmetric factors start from the signed eigenvectors of the starting residual, as in ame();
		# from U = V = 0 a negative weight is reached only by passing through zero
		if (symmetric && R > 0 && !isTRUE(dynamic_uv) &&
				is.null(start_vals[["U"]]) && is.null(start_vals[["V"]])) {
			UV0 <- .lame_signed_sym_start(Z, R, unit_sd = unit_sd, canonical_sign = identical(family, "normal"))
			if (!is.null(UV0)) { U <- UV0$U ; V <- UV0$V }
		}

		# symmetric fits and fits without a dyadic correlation keep rho at zero
		if (symmetric || !isTRUE(dcor)) {
			rho <- 0
			use_rho <- FALSE
		} else {
			use_rho <- TRUE
		}
		if (isTRUE(dynamic_rho)) {
			rho <- .lame_expand_rho(rho, N)
		} else {
			rho <- mean(.lame_expand_rho(rho, N), na.rm = TRUE)
		}
	}
	# explicit ordinal cutpoints: Y_int is the 1..K recoded array, alpha the free cutpoints anchored at
	# alpha[1] = 0, tau_prop the random-walk proposal sd in delta space (Robbins-Monro, burn-in only)
	use_explicit_cutpoints <- (family == "ordinal" &&
			identical(ordinal_cutpoints, "explicit"))
	if (use_explicit_cutpoints) {
		# explicit cutpoints have no bipartite branch
		if (bip) {
			cli::cli_warn(c(
				"{.arg ordinal_cutpoints} = {.val explicit} is not available for bipartite ordinal fits.",
				"i" = "Using {.val data_induced} instead."
			))
			use_explicit_cutpoints <- FALSE
		}
	}
	if (use_explicit_cutpoints) {
		ord_lvls <- sort(unique(c(Y[is.finite(Y)])))
		K_ord <- length(ord_lvls)
		if (K_ord < 2L) {
			cli::cli_abort(c(
				"{.arg ordinal_cutpoints} = {.val explicit} requires at least 2 distinct observed categories.",
				"i" = "Found {K_ord}; check {.arg Y} or use a different family."
			))
		}
		Y_int <- array(match(Y, ord_lvls), dim = dim(Y), dimnames = dimnames(Y))
		alpha <- .init_alpha_from_data(Y)
		# random-walk proposal sd in delta = log(diff(alpha)) space, adapted in burn-in and capped to [0.01, 2]
		tau_prop <- 0.5
		tau_target <- 0.30 # MH acceptance target for the Cowles update
		alpha_accept_n <- 0L
		alpha_accept_total <- 0L
		# alpha storage is allocated after beta storage below
	} else {
		alpha <- numeric(0)
		Y_int <- NULL
		tau_prop <- NA_real_
		K_ord <- 0L
	}
	# dynamic ab parameters
	if (dynamic_ab) {
		n_a <- if (bip) nA else n
		n_b <- if (bip) nB else n
		# a and b as n_a x N and n_b x N matrices; a margin given as a path is kept, one given as levels
		# (or not at all) is spread over the periods
		a_path <- is.matrix(a) && ncol(a) == N
		b_path <- is.matrix(b) && ncol(b) == N
		if (length(a) == n_a && !is.matrix(a)) {
			a_mat <- matrix(rep(a, N), nrow = n_a, ncol = N)
			b_mat <- if (b_path) b else matrix(rep(b, N), nrow = n_b, ncol = N)
			# small random variation across time
			for (t in 2:N) {
				a_mat[, t] <- prior$rho_ab_mean * a_mat[, t - 1] + rnorm(n_a, 0, 0.1 * unit_sd)
				if (!b_path) b_mat[, t] <- prior$rho_ab_mean * b_mat[, t - 1] + rnorm(n_b, 0, 0.1 * unit_sd)
			}
		} else if (a_path) {
			a_mat <- a
			if (b_path) {
				b_mat <- b
			} else if (symmetric) {
				b_mat <- a
			} else {
				b_mat <- matrix(rep(if (length(b) == n_b) as.numeric(b) else numeric(n_b), N), nrow = n_b, ncol = N)
				for (t in 2:N) b_mat[, t] <- prior$rho_ab_mean * b_mat[, t - 1] + rnorm(n_b, 0, 0.1 * unit_sd)
			}
		} else {
			# initialize from scratch
			a_init <- init_dynamic_ab_cpp(n_a, N, prior$rho_ab_mean, 0.1 * unit_sd)
			b_init <- init_dynamic_ab_cpp(n_b, N, prior$rho_ab_mean, 0.1 * unit_sd)
			a_mat <- a_init$a
			b_mat <- b_init$b
		}
		# disabled margins are fixed at zero, including their initial paths
		if (!symmetric) {
			if (!isTRUE(rvar)) { a_mat[] <- 0; a[] <- 0 }
			if (!isTRUE(cvar)) { b_mat[] <- 0; b[] <- 0 }
		}
		# AR(1) parameters
		rho_ab <- if (!is.null(start_vals_obj$rho_ab)) start_vals_obj$rho_ab
		else if (!is.null(prior$rho_ab_mean)) prior$rho_ab_mean
		else 0.8
		sigma_ab <- if (!is.null(start_vals_obj$sigma_ab)) start_vals_obj$sigma_ab
		else 0.1 * unit_sd
		if (!symmetric) {
			ab_shift <- numeric(N)
			for (tt in seq_len(N)) {
				af <- is.finite(a_mat[, tt])
				bf <- is.finite(b_mat[, tt])
				ma <- if (any(af)) mean(a_mat[af, tt]) else 0
				mb <- if (any(bf)) mean(b_mat[bf, tt]) else 0
				if (!is.finite(ma)) ma <- 0
				if (!is.finite(mb)) mb <- 0
				a_mat[af, tt] <- a_mat[af, tt] - ma
				b_mat[bf, tt] <- b_mat[bf, tt] - mb
				ab_shift[tt] <- ma + mb
			}
			if (isTRUE(intercept) && length(beta) > 0L && is.finite(beta[1L])) {
				beta[1L] <- beta[1L] + mean(ab_shift, na.rm = TRUE)
			}
		}
		# per-actor levels of the dynamic effects; the paths carry the deviations around them
		a_lev <- rowMeans(a_mat, na.rm = TRUE); a_lev[!is.finite(a_lev)] <- 0
		b_lev <- rowMeans(b_mat, na.rm = TRUE); b_lev[!is.finite(b_lev)] <- 0
	} else {
		ab_centered <- .lame_center_ab_state(a, b, beta, intercept, symmetric)
		a <- ab_centered$a
		b <- ab_centered$b
		beta <- ab_centered$beta
		ab_sparse <- .lame_shrink_sparse_ab(
			a,
			b,
			.lame_row_obs_ab,
			.lame_col_obs_ab,
			prior$ab_min_observed,
			symmetric
		)
		a <- ab_sparse$a
		b <- ab_sparse$b
		a_mat <- NULL
		b_mat <- NULL
		rho_ab <- NULL
		sigma_ab <- NULL
	}

	# dynamic UV parameters
	if (.dyn_uv_has_rank) {
		n_u <- if (bip) nA else n
		n_v <- if (bip) nB else n
		R_u <- if (bip) RA else R
		R_v <- if (bip) RB else R
		# U and V as 3d arrays
		if (length(dim(U)) == 2) {
			U_cube <- array(0, dim = c(n_u, R_u, N))
			V_cube <- array(0, dim = c(n_v, R_v, N))
			for (t in 1:N) {
				U_cube[, , t] <- U
				V_cube[, , t] <- V
			}
		} else {
			U_cube <- U
			V_cube <- V
		}
		U <- matrix(U_cube[, , 1L], nrow = n_u, ncol = R_u)
		V <- matrix(V_cube[, , 1L], nrow = n_v, ncol = R_v)
		# a symmetric fit carries one set of unit-variance coordinates and a diagonal L holding the scale and sign
		# of every direction: the product is U_t L U_t', and V_t = U_t L is what the directed kernel sees
		L_dyn <- NULL
		L_SUM <- NULL
		L_N <- 0L
		if (symmetric && !bip && R_u > 0L) {
			L_dyn <- vapply(seq_len(R_u), function(r) {
				den <- sum(U[, r]^2)
				if (is.finite(den) && den > 1e-12) sum(V[, r] * U[, r]) / den else 1
			}, numeric(1))
			L_dyn[!is.finite(L_dyn) | L_dyn == 0] <- 1
			if (all(U == 0)) {
				# no start positions: directions, signs and scale from the eigen-decomposition of the double-centred
				# mean latent matrix, eigenvectors scaled to unit variance
				Zbar <- apply(Z, c(1L, 2L), mean, na.rm = TRUE)
				Zbar[!is.finite(Zbar)] <- 0
				Zbar <- (Zbar + t(Zbar)) / 2
				diag(Zbar) <- 0
				Zc <- Zbar - outer(rowMeans(Zbar), rep(1, n_u)) -
					outer(rep(1, n_u), colMeans(Zbar)) + mean(Zbar)
				ev <- eigen(Zc, symmetric = TRUE)
				ord <- order(abs(ev$values), decreasing = TRUE)[seq_len(R_u)]
				# a normal response takes the sign that makes the largest coordinate positive, so y and c y start alike
				if (identical(family, "normal")) {
					for (j in ord) {
						if (ev$vectors[which.max(abs(ev$vectors[, j])), j] < 0) ev$vectors[, j] <- -ev$vectors[, j]
					}
				}
				L_dyn <- ev$values[ord] / n_u
				L_dyn[!is.finite(L_dyn)] <- 0
				U0 <- ev$vectors[, ord, drop = FALSE] * sqrt(n_u)
				for (t in seq_len(N)) U_cube[, , t] <- U0
				U <- U0
			}
			# put each start direction at unit variance and fold its scale into L, leaving the product unchanged
			for (r in seq_len(R_u)) {
				sdr <- stats::sd(U_cube[, r, ])
				if (is.finite(sdr) && sdr > 1e-8) {
					U_cube[, r, ] <- U_cube[, r, ] / sdr
					L_dyn[r] <- L_dyn[r] * sdr^2
				}
			}
			U <- matrix(U_cube[, , 1L], nrow = n_u, ncol = R_u)
			for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n_u, R_u) %*% diag(L_dyn, R_u)
			V <- matrix(V_cube[, , 1L], nrow = n_v, ncol = R_v)
		}
		# persistence from the saved state or the resolved prior; never left unset
		rho_uv <- if (!is.null(start_vals_obj$rho_uv)) start_vals_obj$rho_uv
		else if (!is.null(prior$rho_uv_mean)) prior$rho_uv_mean
		else 0.9
		# the coordinate innovation sd carries the unit of the coordinates (rail_units$uv)
		sigma_uv <- if (!is.null(start_vals_obj$sigma_uv)) start_vals_obj$sigma_uv
		else 0.1 * rail_units$uv
		# the symmetric gauge ties the innovation sd to rho_uv
		if (!is.null(L_dyn)) sigma_uv <- .lame_sym_uv_innov_sd(rho_uv)
		# snap rate starts at the beta prior mean
		pi_snap <- if (identical(dynamic_uv_kind, "snap")) {
			prior$snap_pi_a / (prior$snap_pi_a + prior$snap_pi_b)
		} else { NULL }
	} else {
		U_cube <- NULL
		V_cube <- NULL
		rho_uv <- NULL
		sigma_uv <- NULL
		pi_snap <- NULL
	}

	# keep the adapted step sizes of the scale move that ride in on start_vals
	lev_sd_start <- start_vals_obj[["level_move_sd"]]
	dev_sd_start <- start_vals_obj[["dev_move_sd"]]
	marginal_sd_start <- start_vals_obj[["marginal_move_sd"]]
	dhmc_sd_start <- start_vals[["dynamic_marginal_sd"]]
	# sweeps the symmetric reciprocity ramp has already run, so a resumed chain continues it
	sym_ramp_done <- start_vals[["sym_ramp_sweeps"]]
	if (is.null(sym_ramp_done)) sym_ramp_done <- 0
	if (!is.numeric(sym_ramp_done) || length(sym_ramp_done) != 1L ||
			!is.finite(sym_ramp_done) || sym_ramp_done < 0) {
		cli::cli_abort("{.arg start_vals$sym_ramp_sweeps} must be a single non-negative number.")
	}
	checkpoint_start <- start_vals[["checkpoint_state"]]
	rm(list = c('start_vals_obj', 'start_vals'))

	sym_loop_ids <- lapply(1:(nscan + burn), function(x) { rep(sample(1:nrow(U)), 4) })
	asym_loop_ids <- lapply(1:(nscan + burn), function(x) { sample(1:R) })
	# rail hits and per-parameter sampler failures, reported after MCMC
	uv_rail_count <- 0L
	g_rail_count <- 0L
	.lame_cap_env$hits <- 0L
	.lame_cap_env$draws <- 0L
	tryErrorChecks <- list(
		s2 = 0, betaAB = 0, rho = 0, UV = 0, Z = 0, beta = 0,
		rho_sigma_beta = 0, # dynamic_beta hyperparameter updates
		beta_static = 0, # static-block solve in dynamic_ab bipartite path
		ab = 0, # dynamic additive finite-state repairs
		pool = 0, # hierarchical pooling MH step
		G = 0 # G update
	)
	iter <- 1
	.chain_acc <- NULL
	# per-actor slope state
	if (per_actor_active) {
		# the actor-slope innovation variance has default scale unit_var (squared response unit per squared covariate unit)
		pa_defaults <- list(
			sigma_actor_shape = 2,
			sigma_actor_scale = unit_var,
			rho_actor_shape1 = 2,
			rho_actor_shape2 = 2
		)
		for (nm in names(pa_defaults)) {
			if (is.null(prior[[nm]])) prior[[nm]] <- pa_defaults[[nm]]
			value <- prior[[nm]]
			if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
				cli::cli_abort("{.code prior${nm}} must be a finite positive number.")
		}
		if (per_actor_covariate_idx > pd) cli::cli_abort("{.arg per_actor_covariate_idx} exceeds the number of dyadic covariates.")
		per_actor_design_idx <- as.integer(intercept) + (if (bip) 0L else pr + pc) + as.integer(per_actor_covariate_idx)
		per_actor_X <- lapply(Xlist, function(x) matrix(x[, , per_actor_design_idx], dim(x)[1L], dim(x)[2L]))
		per_actor_offset_sum <- array(0, dim(Y))
		n_actors_init <- if (identical(dynamic_beta_per_actor, "row")) {
			if (bip) nA else n
		} else {
			if (bip) nB else n
		}
		theta_actor <- matrix(0, n_actors_init, N)
		theta_actor_prev <- theta_actor
		rho_actor <- 0.8
		sigma_actor_sq <- 0.1
		# storage for per-actor draws / summary
		n_stored_iters_pa <- nscan / odens
		if (keep_per_actor %in% c("draws", "auto")) {
			# memory check for auto mode
			pa_bytes <- 8 * n_stored_iters_pa * n_actors_init * N
			if (keep_per_actor == "auto" && pa_bytes > 250 * 1024^2) {
				cli::cli_warn(c(
					"!" = "{.code keep_per_actor = \"auto\"}: per-actor draws would require {.val {round(pa_bytes/1024^2, 1)}} MB.",
					"i" = "Falling back to summary mode. Pass {.code keep_per_actor = \"draws\"} to force full storage."
				))
				keep_per_actor <- "summary"
			} else {
				THETA_ACTOR <- array(NA_real_, dim = c(n_stored_iters_pa, n_actors_init, N))
			}
		}
		if (keep_per_actor == "summary") {
			theta_actor_mean <- matrix(0, n_actors_init, N)
			theta_actor_m2 <- matrix(0, n_actors_init, N) # streaming second moment
			theta_actor_n_obs <- 0L
		}
		RHO_ACTOR <- numeric(n_stored_iters_pa)
		SIGMA_ACTOR <- numeric(n_stored_iters_pa)
	}
	# actor-slope hyperparameters (R/per_actor_hyper_move.R): exact conditionals every sweep plus an
	# interweaving step with the standardised innovations held fixed (prior$actor_hyper_nc = FALSE omits it)
	pa_hyper_nc <- per_actor_active && !identical(prior[["actor_hyper_nc"]], FALSE)
	pa_nc_stats <- c(sweeps = 0, attempts = 0, accepts = 0)
	# output items
	if (bip) {
		n_beta <- if (!is.null(Xlist) && length(Xlist) > 0) dim(Xlist[[1]])[3] else 1
		BETA <- matrix(nrow = nscan / odens, ncol = n_beta)
	} else {
		BETA <- matrix(nrow = nscan / odens, ncol = dim(X)[3] - pr * symmetric)
	}
	VC <- matrix(nrow = nscan / odens, ncol = 5 - 3 * symmetric)
	# a static bipartite G is summarised by the posterior mean of its gauged draws
	G_sum <- if (isTRUE(bip) && !is.null(G) && is.matrix(G)) array(0, dim = dim(G)) else NULL
	G_n <- 0L
	# storage for the Cowles explicit cutpoints: alpha[1] = 0 is the anchor, so K - 2 free cutpoints are
	# stored; none for K = 2
	if (use_explicit_cutpoints && K_ord >= 3L) {
		ALPHA <- matrix(
			NA_real_,
			nrow = nscan / odens,
			ncol = K_ord - 2L,
			dimnames = list(NULL, paste0("alpha_", seq.int(2L, K_ord - 1L)))
		)
	} else {
		ALPHA <- NULL
	}
	# storage for the dynamic parameters: persistence and innovation-sd chains (fit$sigma_ab, fit$sigma_uv)
	if (dynamic_ab) {
		RHO_AB <- numeric(nscan / odens)
		SIGMA_AB <- numeric(nscan / odens)
	} else {
		RHO_AB <- NULL
		SIGMA_AB <- NULL
	}
	if (dynamic_uv) {
		RHO_UV <- numeric(nscan / odens)
		SIGMA_UV <- numeric(nscan / odens)
	} else {
		RHO_UV <- NULL
		SIGMA_UV <- NULL
	}
	# optional draw-level latent-position storage: dynamic fits keep U_cube / V_cube per stored draw in
	# fit$U_draws / fit$V_draws, static fits keep fit$U_samples / fit$V_samples, bipartite fits also fit$G_samples.
	# additive draws use the same identified parameterisation as the summaries
	A_SAMPLES <- B_SAMPLES <- A_LEVEL_SAMPLES <- B_LEVEL_SAMPLES <- NULL
	ab_draw_index <- integer(0)
	if (isTRUE(posterior_opts$save_ab)) {
		thin_ab <- posterior_opts$thin_ab
		n_ab_draws <- floor((nscan / odens) / thin_ab)
		ab_draw_index <- seq_len(n_ab_draws) * thin_ab
		ab_dims <- if (dynamic_ab) c(N, n_ab_draws) else n_ab_draws
		A_SAMPLES <- array(NA_real_, c(if (bip) nA else n, ab_dims))
		B_SAMPLES <- array(NA_real_, c(if (bip) nB else n, ab_dims))
		if (dynamic_ab) {
			A_LEVEL_SAMPLES <- matrix(NA_real_, if (bip) nA else n, n_ab_draws)
			B_LEVEL_SAMPLES <- matrix(NA_real_, if (bip) nB else n, n_ab_draws)
		}
	}
	# one terminal state per retained draw supplies the joint forecast posterior
	FC_STATE <- list(
		a = matrix(NA_real_, if (bip) nA else n, nscan / odens),
		b = matrix(NA_real_, if (bip) nB else n, nscan / odens)
	)
	if (dynamic_ab) {
		FC_STATE$a_level <- FC_STATE$a
		FC_STATE$b_level <- FC_STATE$b
	}
	if (.lame_has_uv_rank(bip, R, RA, RB)) {
		FC_STATE$U <- array(NA_real_, c(if (bip) nA else n, if (bip) RA else R, nscan / odens))
		FC_STATE$V <- array(NA_real_, c(if (bip) nB else n, if (bip) RB else R, nscan / odens))
		if (bip) FC_STATE$G <- array(NA_real_, c(RA, RB, nscan / odens))
	}
	U_DRAWS <- NULL
	V_DRAWS <- NULL
	U_SAMPLES <- NULL
	V_SAMPLES <- NULL
	G_SAMPLES <- NULL
	if (posterior_opts_default && isTRUE(posterior_opts$save_UV_draws)) {
		if (.dyn_uv_has_rank) {
			est_bytes <- 8 * (nscan / odens) * N * (n_u * R_u + n_v * R_v)
			if (est_bytes > 100 * 1024^2) {
				posterior_opts$save_UV_draws <- FALSE
				est_size <- format(structure(est_bytes, class = "object_size"), units = "auto")
				cli::cli_inform(c(
					"i" = "Latent trajectory draws ({est_size}) are not kept by default, so post-hoc {.fn gof} and {.code simulate(type = \"posterior\")} are unavailable for this fit.",
					"i" = "Set {.code posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE, thin_UV = 1, thin_ab = 1)} to keep them."
				))
			}
		} else if (!.lame_has_uv_rank(bip, R, RA, RB)) {
			posterior_opts$save_UV_draws <- FALSE
		}
	}
	if (isTRUE(posterior_opts$save_UV_draws)) {
		n_uv_draws <- nscan / odens
		if (.dyn_uv_has_rank) {
			U_DRAWS <- array(NA_real_, dim = c(n_u, R_u, N, n_uv_draws))
			V_DRAWS <- array(NA_real_, dim = c(n_v, R_v, N, n_uv_draws))
			est_bytes <- 8 * (prod(dim(U_DRAWS)) + prod(dim(V_DRAWS)))
			if (est_bytes > 100 * 1024^2) {
				est_size <- format(
					structure(est_bytes, class = "object_size"),
					units = "auto"
				)
				cli::cli_inform(c(
					"i" = "{.arg posterior_opts$save_UV_draws} will allocate approximately {est_size} for the latent trajectory draws ({n_uv_draws} draws x {N} periods).",
					"i" = "Increase {.arg odens} or decrease {.arg nscan} to bound this."
				))
			}
		} else if (.lame_has_uv_rank(bip, R, RA, RB)) {
			# static multiplicative effects: per-draw latent-position matrices
			U_SAMPLES <- array(
				NA_real_,
				dim = c(if (bip) nA else n, if (bip) RA else R, n_uv_draws)
			)
			# symmetric fits store V = U L per draw so U L U' can be rebuilt draw by draw
			V_SAMPLES <- array(
				NA_real_,
				dim = c(if (bip) nB else n, if (bip) RB else R, n_uv_draws)
			)
		} else {
			cli::cli_warn(c(
				"{.arg posterior_opts} latent draw storage ({.code save_UV} / {.code save_UV_draws}) requires a positive multiplicative rank; ignoring it.",
				"i" = "Refit with {.code R > 0} (or positive {.arg R_row} / {.arg R_col}) to store latent-position draws."
			))
		}
		if (bip && RA > 0 && RB > 0) {
			G_SAMPLES <- array(NA_real_, dim = c(RA, RB, n_uv_draws))
		}
	}
	if (dynamic_uv && identical(dynamic_uv_kind, "snap")) {
		PI_SNAP <- numeric(nscan / odens)
	} else {
		PI_SNAP <- NULL
	}

	# UV posterior sums, dynamic or static
	R_eff <- if (bip) max(RA, RB) else R
	# snap-indicator accumulators (sender and receiver sides)
	DELTA_SUM <- NULL
	DELTA_SUM_V <- NULL
	DELTA_COUNT <- NULL
	DELTA_COUNT_V <- NULL
	SNAP_DRAWS <- NULL
	SNAP_DRAWS_V <- NULL
	SNAP_DRAWS_MODE <- keep_snap_draws
	if (dynamic_uv && identical(dynamic_uv_kind, "snap") && R_eff > 0) {
		DELTA_SUM <- matrix(0, if (bip) nA else n, N)
		DELTA_SUM_V <- matrix(0, if (bip) nB else n, N)
		DELTA_COUNT <- matrix(0L, if (bip) nA else n, N)
		DELTA_COUNT_V <- matrix(0L, if (bip) nB else n, N)
		if (SNAP_DRAWS_MODE %in% c("draws", "chunked")) {
			snap_bytes <- 8 * (nscan / odens) * N *
				((if (bip) nA else n) + if (!symmetric) (if (bip) nB else n) else 0)
			if (identical(SNAP_DRAWS_MODE, "chunked")) {
				cli::cli_warn(c(
					"{.code keep_snap_draws = \"chunked\"} currently stores snap indicators as in-memory arrays.",
					"i" = "Estimated storage for this fit is {.val {round(snap_bytes / 1024^2, 1)}} MB.",
					"i" = "Use {.code keep_snap_draws = \"summary\"} for cell-level summaries without retained draw arrays."
				))
			}
		}
		if (SNAP_DRAWS_MODE %in% c("draws", "chunked")) {
			SNAP_DRAWS <- array(NA_real_, dim = c(nscan / odens, if (bip) nA else n, N))
			if (!symmetric) {
				SNAP_DRAWS_V <- array(NA_real_, dim = c(nscan / odens, if (bip) nB else n, N))
			}
		}
	} else if (!identical(SNAP_DRAWS_MODE, "none")) {
		cli::cli_warn(c(
			"{.arg keep_snap_draws} only applies to MCMC fits with {.code dynamic_uv = TRUE} and {.code dynamic_uv_kind = \"snap\"}.",
			"i" = "No snap-indicator draws will be stored for this fit."
		))
	}
	if (dynamic_uv && R_eff > 0) {
		if (bip) {
			UVPS <- array(0, dim = c(nA, nB, N))
			U_SUM <- array(0, dim = c(nA, RA, N))
			V_SUM <- array(0, dim = c(nB, RB, N))
			U_SUM_SQ <- array(0, dim = c(nA, RA, N))
			V_SUM_SQ <- array(0, dim = c(nB, RB, N))
		} else {
			UVPS <- array(0, dim = c(n, n, N))
			U_SUM <- array(0, dim = c(n, R_eff, N))
			V_SUM <- array(0, dim = c(n, R_eff, N))
			U_SUM_SQ <- array(0, dim = c(n, R_eff, N))
			V_SUM_SQ <- array(0, dim = c(n, R_eff, N))
		}
	} else if (R_eff > 0) {
		if (bip) {
			UVPS <- matrix(0, nA, nB)
		} else {
			UVPS <- U %*% t(V) * 0
		}
	} else {
		if (bip) {
			UVPS <- matrix(0, nA, nB)
		} else {
			UVPS <- array(0, dim = c(n, n))
		}
	}

	# APS and BPS sized for bipartite
	if (bip) {
		APS <- rep(0, nA) # row/sender effects
		names(APS) <- rowActorSet
		BPS <- rep(0, nB) # column/receiver effects
		names(BPS) <- colActorSet
		APS2 <- rep(0, nA)
		BPS2 <- rep(0, nB)
	} else {
		APS2 <- BPS2 <- rep(0, n)
		APS <- BPS <- rep(0, nrow(Y[, , 1]))
	}
	YPS <- array(0, dim = dim(Y), dimnames = dimnames(Y))
	# GOF storage
	if (!bip) {
		GOF <- array(NA, dim = c(5, N, (nscan / odens) + 1))
		dimnames(GOF)[[1]] <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		if (!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
		dimnames(GOF)[[3]] <- c('obs', 1:(nscan / odens))
		GOF[, , 1] <- apply(Y, 3, gof_stats, mode = "unipartite")
	} else {
		# bipartite: actor heterogeneity and family-appropriate rectangle dependence
		GOF <- array(NA, dim = c(3, N, (nscan / odens) + 1))
		dimnames(GOF)[[1]] <- .lame_bipartite_gof_names(family)
		if (!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
		dimnames(GOF)[[3]] <- c('obs', 1:(nscan / odens))
		GOF[, , 1] <- apply(Y, 3, gof_stats_bipartite, warn_square = FALSE, family = family)
	}
	if (bip) {
		names(APS) <- rownames(Y[, , 1]) # row actors
		names(BPS) <- colnames(Y[, , 1]) # column actors
		rownames(U) <- rownames(Y[, , 1])
		rownames(V) <- colnames(Y[, , 1])
	} else {
		names(APS) <- names(BPS) <- rownames(U) <- rownames(V) <- rownames(Y[, , 1])
	}

	# parameter names, asymmetric case
	if (!symmetric) {
		colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
		if (bip) {
			# bipartite: names from Xlist
			if (!is.null(Xlist) && length(Xlist) > 0) {
				p_names <- dimnames(Xlist[[1]])[[3]]
				if (is.null(p_names) || length(p_names) != ncol(BETA)) {
					nc <- ncol(BETA)
					n_extra <- nc - as.integer(intercept)
					p_names <- if (nc == 0L) {
						character(0)
					} else if (intercept && n_extra > 0) {
						c("intercept", paste0("x", seq_len(n_extra)))
					} else if (intercept) {
						"intercept"
					} else {
						paste0("x", seq_len(nc))
					}
				}
				colnames(BETA) <- p_names
			}
		} else {
			colnames(BETA) <- dimnames(X)[[3]]
		}
	}

	# parameter names, symmetric case
	if (symmetric) {
		colnames(VC) <- c("va", "ve")
		rb <- intercept + seq(1, pr, length = pr) ; cb <- intercept + pr + seq(1, pr, length = pr)
		bnames <- dimnames(X)[[3]]
		# index only when there are nodal covariates to drop: bnames[-numeric(0)] is character(0)
		# and would clear the names
		excl <- c(rb, cb)
		if (intercept) {
			bni <- bnames[1]
			bnn <- gsub("row", bnames[rb], replacement = "node")
			bnd <- bnames[-c(1, excl)]
			colnames(BETA) <- c(bni, bnn, bnd)
		} else {
			bnn <- gsub("row", bnames[rb], replacement = "node")
			bnd <- if (length(excl) > 0) bnames[-excl] else bnames
			colnames(BETA) <- c(bnn, bnd)
		}
	}

	# dynamic_beta: parse, initialise state and build the design caches the FFBS sampler needs;
	# with dynamic_beta = FALSE every branch falls through to the static path
	beta_dyn_labels <- .lame_coef_blocks(colnames(BETA), Xrow, Xcol, Xdyad, intercept, symmetric, bip)
	beta_dyn <- parse_dynamic_beta(
		dynamic_beta = dynamic_beta,
		coef_names = colnames(BETA),
		coef_block = beta_dyn_labels,
		intercept = intercept,
		Tn = N,
		family = family,
		mode = mode
	)
	# the fit reports the mask and block labels on the stored layout
	beta_dyn_fit <- beta_dyn
	# a symmetric fit stores one coefficient per nodal covariate but samples its row and column slices
	# separately, so the sampler's mask, indices and labels are expanded to the design layout
	if (symmetric && beta_dyn$any && dim(X)[3] > ncol(BETA)) {
		to_store <- integer(dim(X)[3])
		if (intercept) to_store[1L] <- 1L
		to_store[rb] <- intercept + seq_len(pr)
		to_store[cb] <- intercept + seq_len(pr)
		rest <- setdiff(seq_len(dim(X)[3]), c(if (intercept) 1L else integer(0), rb, cb))
		to_store[rest] <- intercept + pr + seq_along(rest)
		beta_dyn$mask <- beta_dyn$mask[to_store]
		beta_dyn$groups <- beta_dyn$groups[to_store]
		beta_dyn$dynamic_idx <- which(beta_dyn$mask)
		beta_dyn$static_idx <- which(!beta_dyn$mask)
		beta_dyn_labels <- beta_dyn_labels[to_store]
	}
	# per-period long-format design matrices and the XtX / n_dyads caches for the FFBS, built only when
	# dynamic_beta is active
	if (beta_dyn$any) {
		p_beta <- ncol(BETA)
		Xdyn_list_per_t <- vector("list", N)
		Xstat_list_per_t <- vector("list", N)
		XtX_per_t <- vector("list", N)
		n_dyads_per_t <- vector("list", N)
		# unipartite X[,,,t] -> n*n x p; bipartite Xlist[[t]] -> nA*nB x p, in the loop's actor order
		for (t in seq_len(N)) {
			if (bip) {
				Xt <- Xlist[[t]] # nA x nB x p
				n_long <- nA * nB
			} else {
				Xt <- X[, , , t] # n x n x p
				if (length(dim(Xt)) == 2L) Xt <- array(Xt, c(n, n, 1))
				n_long <- n * n
			}
			p_t <- dim(Xt)[3]
			Xmat <- matrix(0, n_long, p_t)
			for (k in seq_len(p_t)) Xmat[, k] <- as.numeric(Xt[, , k])
			if (length(beta_dyn$dynamic_idx) > 0L) {
				Xdyn_list_per_t[[t]] <- Xmat[, beta_dyn$dynamic_idx, drop = FALSE]
			} else {
				Xdyn_list_per_t[[t]] <- matrix(0, n_long, 0)
			}
			if (length(beta_dyn$static_idx) > 0L) {
				Xstat_list_per_t[[t]] <- Xmat[, beta_dyn$static_idx, drop = FALSE]
			} else {
				Xstat_list_per_t[[t]] <- matrix(0, n_long, 0)
			}
			XtX_per_t[[t]] <- crossprod(Xmat)
			n_dyads_per_t[[t]] <- n_long
		}
		# innovation scales: one diag(p_dyn) matrix per period
		Lambda_list <- build_beta_state_scales(
			beta_dyn,
			XtX_per_t,
			n_dyads_per_t
		)
		Lambda_full <- Lambda_list[[1]]
		# per-coefficient rho / sigma, broadcast from the per-block parameters
		rho_beta_by_group <- init_rho_beta(beta_dyn$group_names, prior)
		# rw1 pins rho to 1 for every block and skips its sampling
		if (identical(dynamic_beta_kind, "rw1") && length(rho_beta_by_group) > 0L) {
			rho_beta_by_group[] <- 1
		}
		sigma_beta_by_group <- init_sigma_beta(beta_dyn$group_names, prior)
		if (dynamic_beta_pool %in% c("rho", "both") &&
				identical(dynamic_beta_kind, "ar1") && beta_dyn$n_groups >= 2L) {
			pool_rho_mu <- 0
			pool_rho_tau <- 1
		}
		# integer group_id for the C++ FFBS / rho / sigma helpers
		dyn_groups <- beta_dyn$groups[beta_dyn$dynamic_idx]
		group_id <- as.integer(match(dyn_groups, beta_dyn$group_names))
		rho_by_coef <- rho_beta_by_group[group_id]
		sigma_by_coef <- sigma_beta_by_group[group_id]
		# AR(1) state prior at t = 0
		p_dyn <- length(beta_dyn$dynamic_idx)
		beta0_mean <- rep(prior$beta0_mean %||% 0, p_dyn)
		beta0_var <- prior$beta0_var %||% beta0_var_default
		# the default prior variance is per unit of the covariate: the variance over the mean square of the
		# design column (one for the intercept); a variance the caller sets is taken as given
		beta0_ms <- if (is.null(prior$beta0_var)) {
			vapply(beta_dyn$dynamic_idx, function(k) {
				num <- sum(vapply(XtX_per_t, function(m) m[k, k], numeric(1)))
				ms <- num / max(sum(unlist(n_dyads_per_t)), 1)
				if (is.finite(ms) && ms > 0) ms else 1
			}, numeric(1))
		} else rep(1, p_dyn)
		beta0_cov <- diag(beta0_var / beta0_ms, p_dyn)
		# initial beta path: the static draw replicated across periods
		beta_dyn_path <- matrix(0, N, p_dyn)
		for (t in seq_len(N)) beta_dyn_path[t, ] <- beta[beta_dyn$dynamic_idx]
		beta_static <- beta[beta_dyn$static_idx]
		# 3d BETA storage for dynamic beta, [iter, p, t]
		BETA <- array(
			NA_real_,
			dim = c(nscan / odens, p_beta, N),
			dimnames = list(
				NULL,
				colnames(BETA),
				if (!is.null(pd_labs)) pd_labs else paste0("t", seq_len(N))
			)
		)
		# storage for the per-block AR(1) hyperparameters
		RHO_BETA <- matrix(
			NA_real_,
			nrow = nscan / odens,
			ncol = beta_dyn$n_groups,
			dimnames = list(NULL, beta_dyn$group_names)
		)
		SIGMA_BETA <- matrix(
			NA_real_,
			nrow = nscan / odens,
			ncol = beta_dyn$n_groups,
			dimnames = list(NULL, beta_dyn$group_names)
		)
		# memory warning for large 3d BETA storage
		est_bytes <- prod(dim(BETA)) * 8
		if (est_bytes > 1024^3) {
			cli::cli_warn(c(
				"{.arg dynamic_beta} will allocate {.val {round(est_bytes/1024^3, 2)}} GB for 3-D BETA storage ({nscan/odens} draws x {p_beta} coefs x {N} periods).",
				"i" = "Increase {.arg odens}, decrease {.arg nscan}, or pass {.arg dynamic_beta} as a subset of coefficients to bound this."
			))
		}
	} else {
		# defaults for variables downstream branches may reference when beta_dyn$any is FALSE
		Xdyn_list_per_t <- NULL; Xstat_list_per_t <- NULL
		XtX_per_t <- NULL; n_dyads_per_t <- NULL
		Lambda_full <- NULL; Lambda_list <- NULL
		rho_beta_by_group <- numeric(0)
		sigma_beta_by_group <- numeric(0)
		rho_by_coef <- numeric(0); sigma_by_coef <- numeric(0)
		group_id <- integer(0)
		beta_dyn_path <- NULL; beta_static <- NULL
		beta0_mean <- NULL; beta0_cov <- NULL
		RHO_BETA <- NULL; SIGMA_BETA <- NULL
	}
	if (isTRUE(dynamic_rho)) {
		RHO <- matrix(
			NA_real_,
			nrow = nscan / odens,
			ncol = N,
			dimnames = list(NULL, pd_labs %||% paste0("t", seq_len(N)))
		)
	} else {
		RHO <- NULL
	}
	# pointwise log-likelihood storage for loo() / waic(): an in-memory [nscan/odens x n_obs] matrix, or with
	# save_log_lik = "chunked" per-column-chunk binary files that read_log_lik(fit) stitches back together
	log_lik_mode <- if (identical(save_log_lik, "chunked")) "chunked"
	else if (isTRUE(save_log_lik)) "memory"
	else "off"
	if (log_lik_mode != "off") {
		obs_idx <- .lame_score_index(Y, bipartite = bip, symmetric = symmetric)
		n_obs_per_t <- tabulate(obs_idx[, 3L], nbins = N)
		n_obs_total <- sum(n_obs_per_t)
		n_stored <- nscan / odens
		est_bytes <- as.numeric(n_stored) * n_obs_total * 8
		if (log_lik_mode == "memory" && est_bytes > 1024^3) {
			cli::cli_warn(c(
				"{.arg save_log_lik} = TRUE will allocate {.val {round(est_bytes/1024^3, 2)}} GB for the log-lik matrix ({n_stored} draws x {n_obs_total} obs).",
				"i" = "Increase {.arg odens}, switch to {.code save_log_lik = \"chunked\"}, or set {.arg save_log_lik} = FALSE."
			))
		}
		if (log_lik_mode == "memory") {
			LOG_LIK <- matrix(NA_real_, nrow = n_stored, ncol = n_obs_total)
			LOG_LIK_CONNS <- NULL
			LOG_LIK_META <- NULL
		} else {
			# chunked: one binary connection per chunk
			LOG_LIK <- NULL
			ll_dir <- if (is.null(log_lik_path)) {
				tempfile(pattern = "lame_loglik_", tmpdir = tempdir())
			} else log_lik_path
			dir.create(ll_dir, showWarnings = FALSE, recursive = TRUE)
			chunk_size <- max(1L, as.integer(log_lik_chunk_size))
			n_chunks <- as.integer(ceiling(n_obs_total / chunk_size))
			chunk_starts <- seq.int(1L, n_obs_total, by = chunk_size)
			chunk_ends <- pmin(chunk_starts + chunk_size - 1L, n_obs_total)
			chunk_files <- file.path(
				ll_dir,
				sprintf("loglik_chunk_%03d.bin", seq_len(n_chunks))
			)
			LOG_LIK_CONNS <- lapply(
				chunk_files,
				function(p)
					file(p, open = "wb")
			)
			LOG_LIK_META <- list(
				path = ll_dir,
				files = chunk_files,
				chunk_size = chunk_size,
				n_obs_total = n_obs_total,
				n_stored = n_stored,
				chunk_starts = chunk_starts,
				chunk_ends = chunk_ends,
				n_chunks = n_chunks
			)
		}
	} else {
		LOG_LIK <- NULL
		LOG_LIK_CONNS <- NULL
		LOG_LIK_META <- NULL
		obs_idx <- NULL
	}
	# randomised Halton shift for the GHK pointwise log-likelihood, deterministic given the fit `seed`
	.lame_log_lik_halton_shift <- if (identical(log_lik_method, "observed_ghk")) {
		# a small RNG kick that leaves the main MCMC stream untouched (.Random.seed is restored)
		old_seed <- if (exists(
			".Random.seed",
			envir = globalenv(),
			inherits = FALSE
		)) {
			get(".Random.seed", envir = globalenv())
		} else NULL
		set.seed(seed + 7919L) # offset so it is not collinear with main seed
		shift <- stats::runif(1)
		if (!is.null(old_seed)) {
			assign(".Random.seed", old_seed, envir = globalenv())
		}
		shift
	} else 0
	# optional GHK Monte Carlo budget: prior$ghk_n_mc sets the exact per-row draw count,
	# overriding the adaptive 8 * 2^D floor
	.lame_ghk_n_mc <- if (is.list(prior)) prior$ghk_n_mc else NULL
	if (!is.null(.lame_ghk_n_mc)) {
		if (!is.numeric(.lame_ghk_n_mc) || length(.lame_ghk_n_mc) != 1L ||
				!is.finite(.lame_ghk_n_mc) || .lame_ghk_n_mc < 1) {
			cli::cli_abort(c(
				"{.code prior$ghk_n_mc} must be a single positive integer.",
				"i" = "It sets the exact per-row Monte Carlo budget for {.code log_lik_method = \"observed_ghk\"}."
			))
		}
		.lame_ghk_n_mc <- as.integer(.lame_ghk_n_mc)
	}
	# MCMC
	have_coda <- suppressWarnings(
		try(requireNamespace("coda", quietly = TRUE), silent = TRUE)
	)

	# numerical stability message, shown once
	if (verbose && symmetric) {
		cli::cli_div(theme = list(span.note = list(color = "grey60", "font-style" = "italic")))
		cli::cli_inform(c(
			"i" = "Covariance matrices are adjusted to positive definiteness when a draw needs it.",
			" " = cli::col_grey("The corrections are small."),
			" " = cli::col_grey("This is normal and does not affect model validity.")
		))
		cli::cli_end()
	}

	# residual array, allocated once
	E.nrm <- array(dim = dim(Z))
	if (burn != 0 && verbose) {
		cli::cli_h3("Starting burn-in period...")
		cli::cli_progress_bar("Burn-in", total = burn, .envir = environment())
	}
	# with dynamic_beta the per-period X beta cube is added to the EZ of the static helpers,
	# which are called with beta = 0
	BETA_static_view <- function() {
		matrix(0, 1, length(beta))
	}
	.lame_beta_full_path <- function() {
		# N x p matrix of per-period beta when dynamic_beta is on, NULL otherwise
		if (!beta_dyn$any) return(NULL)
		out <- matrix(0, N, length(beta))
		for (t_ in seq_len(N)) {
			out[t_, ] <- assemble_beta_path(
				beta_dyn_path[t_, ],
				beta_static,
				beta_dyn$dynamic_idx,
				beta_dyn$static_idx,
				length(beta)
			)
		}
		out
	}
	.lame_xbeta_cube <- function(beta_full_path) {
		if (is.null(beta_full_path)) return(NULL)
		n_a_ <- if (bip) nA else n
		n_b_ <- if (bip) nB else n
		out <- array(0, dim = c(n_a_, n_b_, N))
		for (t_ in seq_len(N)) {
			if (bip) Xt_ <- Xlist[[t_]] else Xt_ <- X[, , , t_]
			if (length(dim(Xt_)) == 2L) Xt_ <- array(Xt_, c(n_a_, n_b_, 1))
			p_t_ <- dim(Xt_)[3]
			b_t_ <- beta_full_path[t_, ]
			if (p_t_ > 0 && length(b_t_) >= p_t_) {
				for (k_ in seq_len(p_t_)) {
					if (b_t_[k_] != 0) out[, , t_] <- out[, , t_] + b_t_[k_] * Xt_[, , k_]
				}
			}
		}
		out
	}
	# snap indicators carried across scans for the storage gate
	delta_uv <- NULL
	delta_uv_v <- NULL
	# local scales carried across scans for the heavy-tailed dynamic UV kind
	if (dynamic_uv && identical(dynamic_uv_kind, "t")) {
		lambda_u <- matrix(1, n, N)
		lambda_v <- matrix(1, n, N)
	}
	checkpoint_layout <- list(
		family = family,
		mode = mode,
		symmetric = symmetric,
		dim_y = dim(Y),
		names_y = dimnames(Y),
		rank = R,
		intercept = intercept,
		rvar = rvar,
		cvar = cvar,
		dynamic_ab = dynamic_ab,
		dynamic_uv = dynamic_uv,
		dynamic_uv_kind = dynamic_uv_kind,
		dynamic_rho = dynamic_rho,
		dynamic_beta = beta_dyn_fit,
		dynamic_beta_kind = dynamic_beta_kind,
		dynamic_beta_pool = dynamic_beta_pool,
		per_actor = dynamic_beta_per_actor,
		ordinal_cutpoints = ordinal_cutpoints
	)
	if (per_actor_active) checkpoint_layout$per_actor_model_version <- 2L
	# the current sampler state as a start_vals list, written into every checkpoint and returned as fit$start_vals
	.lame_snapshot_state <- function() {
		st <- list(
			Z = Z,
			beta = beta,
			a = if (dynamic_ab) a_mat else a,
			b = if (dynamic_ab) b_mat else b,
			U = if (.dyn_uv_has_rank) U_cube else U,
			V = if (.dyn_uv_has_rank) V_cube else V,
			rho = rho,
			s2 = s2,
			Sab = Sab
		)
		if (.dyn_uv_has_rank) {
			st$rho_uv <- rho_uv
			st$sigma_uv <- sigma_uv
		}
		if (dynamic_ab) {
			st$rho_ab <- rho_ab
			st$sigma_ab <- sigma_ab
		}
		# the adapted step sizes of the scale move, so a continued run keeps its tuned step
		if (isTRUE(level_move)) {
			st$level_move_sd <- lev_sd
			if (isTRUE(dev_move)) st$dev_move_sd <- dev_sd
		}
		if (isTRUE(bip_level_move)) st$level_move_sd <- lev_sd
		if (isTRUE(marginal_move)) st$marginal_move_sd <- marginal_sd
		if (isTRUE(dhmc_enabled)) st$dynamic_marginal_sd <- dhmc_sd
		# the cumulative sweep count of the symmetric reciprocity ramp
		if (symmetric && !bip) st$sym_ramp_sweeps <- sym_ramp_done
		if (bip && !is.null(G)) {
			st$G <- G
		}
		st$checkpoint_state <- .lame_checkpoint_capture(
			parent.env(environment()),
			checkpoint_layout
		)
		st
	}
	# joint move of a sender's level and its latent row for cbin / frn panels (.lame_nom_row_shift): every sweep
	# for a static fit, on the schedule of prior$actor_move / prior$actor_move_every under dynamic_ab
	nom_shift <- family %in% c("cbin", "frn") && !symmetric && isTRUE(rvar) &&
		(if (dynamic_ab) length(a_lev) == dim(Y)[1L] else (is.null(dim(a)) && length(a) == dim(Y)[1L]))
	nom_shift_every <- 1L
	if (nom_shift && dynamic_ab) {
		if (!is.null(prior[["actor_move"]])) nom_shift <- isTRUE(prior[["actor_move"]])
		nom_shift_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
	}
	if (nom_shift) {
		nom_shift_masks <- .lame_nom_shift_masks(Y, bip)
		# rows the sparse-actor shrink leaves alone
		nom_shift_rows <- .lame_row_obs_ab >= prior$ab_min_observed
	}
	# receiver counterpart (R/nomination_actor_move.R): a receiver's level moves with its latent column, compensated
	# by the other levels and the intercept; prior$actor_move and prior$actor_move_every control it
	nom_gaps <- !is.null(time_index) && any(abs(diff(as.numeric(time_index)) - 1) > 1e-10)
	nom_beta_ok <- !isTRUE(beta_dyn$any) || (dynamic_beta_kind %in% c("ar1", "rw1") && !nom_gaps)
	nom_col <- family %in% c("cbin", "frn") && !symmetric && isTRUE(cvar) &&
		(if (dynamic_ab) length(b_lev) == dim(Y)[2L] else (is.null(dim(b)) && length(b) == dim(Y)[2L])) &&
		min(dim(Y)[1:2]) >= 3L &&
		(!isTRUE(intercept) || (length(beta) > 0L && nom_beta_ok))
	if (!is.null(prior[["actor_move"]])) nom_col <- nom_col && isTRUE(prior[["actor_move"]])
	nom_col_stats <- c(visits = 0, col_calls = 0, col_skips = 0)
	if (nom_col) {
		nom_col_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
		nom_col_ok <- .lame_col_obs_ab >= prior$ab_min_observed
		# precision of the coefficient prior (the g-prior of the replicated draw)
		nom_col_iV0 <- if (isTRUE(intercept) && !bip) apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1])) else NULL
	}
	# level move (.lame_nom_level_shift): rows free above move with the intercept against the rows with room;
	# it needs the row effects, which is why nom_shift requires rvar
	nom_lvl <- nom_shift && isTRUE(intercept) && length(beta) > 0L && nom_beta_ok
	if (!is.null(prior[["actor_move"]])) nom_lvl <- nom_lvl && isTRUE(prior[["actor_move"]])
	nom_lvl_stats <- c(visits = 0, calls = 0, free_rows = 0)
	if (nom_lvl) {
		nom_lvl_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
		nom_lvl_iV0 <- if (!bip) apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1])) else NULL
	}
	# period-level counterpart (R/nomination_period_move.R): a dynamic intercept moves period by period;
	# prior$actor_move_period = FALSE omits it
	nom_per <- nom_lvl && isTRUE(beta_dyn$any) && N >= 2L && !identical(prior[["actor_move_period"]], FALSE)
	nom_per_stats <- c(visits = 0, calls = 0, comps = 0)
	# collapsed translation of the periods that have rows with room (R/nomination_room_move.R);
	# prior$actor_move_room = FALSE omits it
	nom_room <- nom_per && !identical(prior[["actor_move_room"]], FALSE) && !is.na(match(1L, beta_dyn$dynamic_idx))
	nom_room_stats <- c(visits = 0, attempts = 0, accepts = 0, stuck = 0)
	if (nom_room) {
		nom_room_ctx <- .lame_nom_room_context(Y, odmax, bip)
		# nothing to integrate where no period has a row with room
		nom_room <- any(vapply(nom_room_ctx, function(ct) any(ct$pinned), NA))
		nom_room_k <- match(1L, beta_dyn$dynamic_idx)
	}
	# collapsed reciprocity update (.lame_nom_rho_collapsed): rho with the non-nominations integrated out,
	# every tenth visit of the translations
	nom_rho <- family %in% c("cbin", "frn") && !bip && !symmetric && isTRUE(dcor) && !isTRUE(dynamic_rho)
	if (!is.null(prior[["actor_move"]])) nom_rho <- nom_rho && isTRUE(prior[["actor_move"]])
	if (identical(prior[["actor_move_rho"]], FALSE)) nom_rho <- FALSE
	nom_rho_stats <- c(calls = 0, attempts = 0, accepts = 0)
	if (nom_rho) nom_rho_every <- 10L * max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
	# static positions of directed cbin, frn and ordinal fits move with the latent rows and columns they bound
	# (R/sparse_position_move.R); prior$sparse_position_move = FALSE omits it
	spm_on <- family %in% c("cbin", "frn", "ordinal") && !bip && !symmetric && R > 0L && !isTRUE(dynamic_uv) &&
		!isTRUE(.dyn_uv_has_rank) && !identical(prior[["sparse_position_move"]], FALSE)
	spm_stats <- c(visits = 0, spm_calls = 0, spm_evals = 0, spm_stuck = 0)
	spm_every <- max(1L, as.integer(prior[["sparse_position_every"]] %||% 1L))
	# coordinate paths of cbin and frn fits with dynamic_uv move likewise (R/nomination_path_move.R);
	# prior$sparse_position_every sets the schedule
	npm_on <- family %in% c("cbin", "frn") && isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "ar1") && !symmetric &&
		(if (bip) RA > 0L && RB > 0L else R > 0L) && !identical(prior[["sparse_position_move"]], FALSE) &&
		!any(vapply(seq_len(dim(Y)[3L]), function(tt) any(.lame_nom_hidden_cells(Y[, , tt], bip)), logical(1)))
	npm_stats <- c(visits = 0, npm_calls = 0, npm_evals = 0, npm_stuck = 0)
	npm_ctx <- if (npm_on) .lame_npm_context(
		Y,
		family,
		if (identical(family, "frn")) YL else NULL,
		odmax,
		bip
	) else NULL
	npm_every <- max(1L, as.integer(prior[["sparse_position_every"]] %||% 2L))
	# ridge of a row or column covariate with the additive effects of its margin (R/covariate_ridge_move.R);
	# prior$cov_ridge_move = FALSE omits it
	cr_on <- !bip && !symmetric && (isTRUE(rvar) || isTRUE(cvar)) && !isTRUE(beta_dyn$any) && length(beta) > 0L &&
		is.null(dynamic_beta_per_actor) && !is.null(xxLong) && !identical(prior[["cov_ridge_move"]], FALSE)
	cr_stats <- c(visits = 0, ridge_calls = 0, ridge_moved = 0)
	if (cr_on) {
		cr_P0 <- apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1]))
		cr_P0 <- (cr_P0 + t(cr_P0)) / 2
		# ordinal under the rank likelihood moves every latent value of a period with its predictor; other families
		# keep the latent values of observed cells (cbin / frn panels with hidden cells are not served)
		cr_M <- if (identical(family, "ordinal") && !isTRUE(use_explicit_cutpoints)) NULL else array(is.finite(Y), dim(Y))
		if (family %in% c("cbin", "frn") && any(
			vapply(seq_len(N), function(tt) any(.lame_nom_hidden_cells(Y[, , tt], bip)), logical(1))
		)) cr_on <- FALSE
	}
	if (spm_on) {
		spm_ctx <- .lame_spm_context(
			Y,
			family,
			if (identical(family, "frn")) YL else NULL,
			if (family %in% c("cbin", "frn")) odmax else NULL,
			Y_int = if (use_explicit_cutpoints) Y_int else NULL
		)
		spm_on <- isTRUE(spm_ctx$ok)
	}
	# dynamic_beta path hyperparameters (R/dynamic_beta_hyper_move.R): exact conditionals every sweep plus an
	# interweaving step (prior$beta_hyper_nc = FALSE omits it)
	beta_hyper_nc <- isTRUE(beta_dyn$any) && !identical(prior[["beta_hyper_nc"]], FALSE)
	beta_nc_stats <- c(sweeps = 0, attempts = 0, accepts = 0)
	# joint scale move along the level / spread ridge of the binary unipartite sampler, where a probit panel
	# identifies the predictor scale weakly; prior$level_move_every sets how often it runs
	static_path_move <- identical(family, "binary") && !bip && !symmetric &&
		isTRUE(intercept) && isTRUE(rvar) && isTRUE(cvar) && !isTRUE(dcor) &&
		isTRUE(dynamic_ab) && !isTRUE(.dyn_uv_has_rank) && R %in% c(1L, 2L) &&
		isTRUE(beta_dyn$any) && dynamic_beta_kind %in% c("ar1", "rw1") &&
		is.null(dynamic_beta_per_actor) && is.null(time_index) &&
		all(beta0_mean == 0) && all(Lambda_full[row(Lambda_full) != col(Lambda_full)] == 0)
	level_move <- identical(family, "binary") && !bip && !symmetric &&
		isTRUE(intercept) && length(beta) > 0L &&
		isTRUE(rvar) && isTRUE(cvar) &&
		!isTRUE(beta_dyn$any) && is.null(dynamic_beta_per_actor) &&
		!isTRUE(dynamic_rho) &&
		(R == 0L || !isTRUE(dynamic_uv) ||
				(isTRUE(.dyn_uv_has_rank) && identical(dynamic_uv_kind, "ar1")))
	# static positions scale with the predictor under their integrated inverse-Wishart (matrix-t) prior
	static_uv_move <- level_move && R > 0L && !isTRUE(dynamic_uv)
	level_move <- level_move || static_path_move
	static_uv_move <- static_uv_move && level_move
	lev_uv <- level_move && R > 0L && (isTRUE(.dyn_uv_has_rank) || static_path_move || static_uv_move)
	lev_att <- 0L
	translation_att <- 0L
	# the bipartite scale step leaves the predictor and so the likelihood unchanged; poisson and ordinal
	# take it with the actor moves
	bip_scale_enabled <- bip && isTRUE(.dyn_uv_has_rank) &&
		(family %in% c("normal", "binary") || (family %in% c("poisson", "ordinal") && !identical(prior[["actor_move"]], FALSE))) &&
		identical(dynamic_uv_kind, "ar1") && RA > 0L && RB > 0L &&
		!identical(prior[["bipartite_scale_move"]], FALSE)
	bip_scale_att <- bip_scale_acc <- bip_scale_bounds <- bip_scale_numerical <- 0L
	uv_scale_att <- uv_scale_acc <- 0L
	rho_collapsed_att <- rho_collapsed_acc <- 0L
	lev_acc <- 0L
	# a resumed run carries the step it had adapted (see .lame_snapshot_state)
	lev_sd <- if (is.numeric(lev_sd_start)) {
		min(max(lev_sd_start, 1e-5), 1)
	} else 0.02
	if (level_move) {
		lev_every <- max(1L, as.integer(prior[["level_move_every"]] %||% 30L))
		# several proposals per visit: the factor draw and z redraw are paid once, so extra proposals are cheap
		lev_tries <- max(1L, as.integer(prior[["level_move_tries"]] %||% 8L))
		lev_uv_max <- as.numeric(prior$uv_max_abs %||% Inf)
		if (!is.finite(lev_uv_max) || lev_uv_max <= 0) lev_uv_max <- Inf
		# actor levels and deviations live on centred subspaces; latent coordinates keep all n * R * T dimensions
		lev_cen_lev <- 2L
		lev_cen_dev <- if (dynamic_ab) 2L * N else 0L
		lev_cen_uv <- 0L # dynamic coordinates are not centered by the state kernel
		lev_upper <- upper.tri(matrix(0, n, n))
		lev_dg <- array(FALSE, dim = dim(Y))
		lev_dg[cbind(
			rep(seq_len(n), N),
			rep(seq_len(n), N),
			rep(seq_len(N), each = n)
		)] <- TRUE
		lev_i1 <- which(Y == 1 & !lev_dg)
		lev_i0 <- which(Y == 0 & !lev_dg)
		lev_ina <- which(is.na(Y) & !lev_dg)
		lev_idg <- which(lev_dg)
		lev_obs <- c(lev_i1, lev_i0)
		lev_sY <- c(rep(1, length(lev_i1)), rep(-1, length(lev_i0)))
		rm(lev_dg)
		# prior precision of the coefficients: the kernel adds sum_t XX_t / (g N)
		lev_iV0 <- apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1]))
		# a panel with no events either way leaves the step no signal
		if (length(lev_i1) == 0L || length(lev_i0) == 0L) level_move <- FALSE
	}
	# non-centred draws of the persistence of a symmetric fit and of the scale of a bipartite fit
	# (R/uv_nc_step.R); prior$uv_nc_move = FALSE omits them
	uv_nc_move <- !identical(prior[["uv_nc_move"]], FALSE) && isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "ar1")
	uv_nc_stats <- c(sym_calls = 0, sym_moved = 0, bip_calls = 0, bip_moved = 0)
	# a count panel carries its latent array along with the new persistence
	uv_nc_cx <- NULL
	if (uv_nc_move && symmetric && identical(family, "poisson"))
		uv_nc_cx <- .lame_ltm_context(
			Y,
			"poisson",
			"sym",
			offset = if (poisson_exposure_active) log(period_exposure) else NULL
		)
	# tempered transition of one eigen-direction of a symmetric count fit, crossing between weight signatures
	# (R/uv_tempered_move.R); prior$uv_temper_move and prior$uv_temper_every control it
	uv_temper_dt <- NULL
	uv_temper_every <- max(1L, as.integer(prior[["uv_temper_every"]] %||% 5L))
	uv_temper_bs <- seq(1, 0.7, length.out = 33L)
	uv_temper_stats <- c(calls = 0, accepted = 0, flipped = 0)
	if (!is.null(uv_nc_cx) && isTRUE(
		prior[["uv_temper_move"]] %||% TRUE
	)) uv_temper_dt <- .lame_temper_data(uv_nc_cx, dim(Y)[1L], dim(Y)[3L])
	# exact draw of the split of scale between the static factors U and V (R/uv_scale_split.R);
	# prior$uv_split_move = FALSE switches it off
	uv_split_on <- !bip && !symmetric && R > 0L && !isTRUE(dynamic_uv) && !isTRUE(.dyn_uv_has_rank) && !is.null(uv_prior_static) &&
		!identical(prior[["uv_split_move"]], FALSE)
	uv_split_stats <- c(visits = 0, moved = 0)
	# tie the row and column coefficients of a nodal covariate of a symmetric binary or ordinal fit after
	# every draw (R/symmetric_node_tie.R); prior$symmetric_node_tie = FALSE omits it
	sym_tie <- symmetric && !bip && length(rb) > 0L && family %in% c("binary", "ordinal") &&
		!identical(prior[["symmetric_node_tie"]], FALSE)
	# bipartite counterpart of the level direction (R/bipartite_predictor_scale.R) for binary panels with an
	# intercept, both margins and static coefficients; schedule and proposals as in the unipartite move
	bip_level_move <- bip && identical(family, "binary") && isTRUE(intercept) && length(beta) > 0L &&
		isTRUE(rvar) && isTRUE(cvar) && !isTRUE(beta_dyn$any) && is.null(dynamic_beta_per_actor) &&
		nA >= 3L && nB >= 3L
	if (bip_level_move) {
		bpl_cx <- .lame_bpl_context(Y)
		# a panel with no events either way leaves the step no signal
		if (bpl_cx$n1 == 0L || bpl_cx$n0 == 0L) bip_level_move <- FALSE
		bpl_every <- max(1L, as.integer(prior[["level_move_every"]] %||% 30L))
		bpl_tries <- max(1L, as.integer(prior[["level_move_tries"]] %||% 8L))
		# the coefficients' prior precision, N(0, g s2) each
		bpl_iV0 <- diag(length(beta)) / (as.numeric(g[1]) * s2)
	}
	# whitened dynamic observed-data HMC and its companion moves (translation, factor scale, non-centred AR,
	# collapsed rho, actor slice); valid at any rank, enabled at rank two for efficiency
	dynamic_joint_move <- level_move && isTRUE(.dyn_uv_has_rank) &&
		dynamic_ab && R == 2L && !isTRUE(beta_dyn$any) && is.null(time_index) &&
		!identical(prior[["marginal_move"]], FALSE) && rvar && cvar &&
		is.null(dynamic_beta_per_actor)
	actor_slice_enabled <- dynamic_joint_move && n >= 3L &&
		!identical(prior[["actor_factor_slice"]], FALSE)
	actor_slice_rounds <- actor_slice_blocks <- actor_slice_evaluations <-
		actor_slice_bounds <- actor_slice_bad <- 0L
	dhmc_enabled <- dynamic_joint_move
	rho_collapsed_context <- if (dynamic_joint_move && isTRUE(dcor)) .lame_binary_dynamic_collapsed_rho_context(Y) else NULL
	dhmc_context <- if (dhmc_enabled) .lame_binary_dynamic_pair_hmc_context(Y, Xlist, R) else NULL
	dhmc_sd <- if (is.null(dhmc_sd_start)) .12 else dhmc_sd_start
	if (!is.numeric(dhmc_sd) || length(dhmc_sd) != 1L || !is.finite(dhmc_sd) ||
			dhmc_sd < .001 || dhmc_sd > .5) cli::cli_abort("Invalid dynamic marginal move step size.")
	ncar_enabled <- dynamic_joint_move
	ncar_uv_att <- ncar_uv_acc <- ncar_ab_att <- ncar_ab_acc <- 0L
	dhmc_att <- dhmc_acc <- dhmc_nonfinite <- dhmc_bounds <- 0L
	# second direction of the visit: the additive deviations together with their innovation scale, a ridge
	# the separate deviation and sigma_ab draws cross slowly; the joint move leaves their Gaussian prior invariant
	dev_move <- level_move && isTRUE(dynamic_ab)
	if (!is.null(prior[["dev_move"]])) {
		dev_move <- dev_move && isTRUE(prior[["dev_move"]])
	}
	dev_att <- 0L
	dev_acc <- 0L
	dev_sd <- if (is.numeric(dev_sd_start)) {
		min(max(dev_sd_start, 1e-5), 1)
	} else 0.05
	if (dev_move) {
		dev_tries <- max(1L, as.integer(prior[["dev_move_tries"]] %||% 8L))
		# index the deviation paths of the observed cells once, so the visit only gathers
		lev_row <- ((lev_obs - 1L) %% n) + 1L
		lev_col <- ((lev_obs - 1L) %/% n) %% n + 1L
		lev_per <- ((lev_obs - 1L) %/% (n * n)) + 1L
		dev_ia <- lev_row + (lev_per - 1L) * n
		dev_ib <- lev_col + (lev_per - 1L) * n
		rm(lev_row, lev_col, lev_per)
	}
	# collapsed per-actor update (R/binary_actor_collapsed_move.R): with the augmentation integrated out, each
	# actor's level translation and static position is an exact elliptical slice draw; data augmentation alone
	# moves an actor without events only slowly. prior$actor_move and prior$actor_move_every control it
	am_lev <- if (symmetric) rep(isTRUE(nvar), 2L) else c(isTRUE(rvar), isTRUE(cvar))
	am_centred <- isTRUE(intercept) && length(beta) > 0L
	am_uv_ok <- R > 0L && !isTRUE(dynamic_uv) && !isTRUE(.dyn_uv_has_rank)
	am_gaps <- !is.null(time_index) && any(abs(diff(as.numeric(time_index)) - 1) > 1e-10)
	actor_move <- identical(family, "binary") && !bip && (any(am_lev) || am_uv_ok) &&
		is.null(dynamic_beta_per_actor) && !isTRUE(dynamic_rho) && n >= 3L &&
		(!isTRUE(beta_dyn$any) || (dynamic_beta_kind %in% c("ar1", "rw1") && !am_gaps))
	if (!is.null(prior[["actor_move"]])) actor_move <- actor_move && isTRUE(prior[["actor_move"]])
	actor_move_uv <- actor_move && am_uv_ok
	actor_move_stats <- c(visits = 0, level_calls = 0, level_evals = 0, uv_calls = 0, uv_evals = 0, stuck = 0)
	if (actor_move) {
		am_cx <- if (symmetric) .lame_actor_move_context_sym(Y) else .lame_actor_move_context(Y)
		am_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 10L))
		am_iV0 <- if (am_centred && !isTRUE(beta_dyn$any)) apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1]))
	}
	# exact reciprocity update with both latent responses of every dyad integrated out, visited with the
	# actor update, since given the augmented array rho moves by only about 1 / sqrt(dyads)
	actor_move_rho <- actor_move && !symmetric && isTRUE(dcor) && !identical(prior[["actor_move_rho"]], FALSE)
	if (actor_move_rho) am_rho_cx <- .lame_binary_dynamic_collapsed_rho_context(Y)
	actor_move_stats <- c(actor_move_stats, rho_attempts = 0, rho_accepts = 0)
	# collapsed per-actor update for bipartite networks (R/binary_actor_collapsed_move_bip.R): row and column
	# level translations and static positions as exact elliptical slice draws; same switches as the unipartite update
	actor_move_bip <- identical(family, "binary") && bip &&
		(isTRUE(rvar) || isTRUE(cvar) || (RA > 0L && RB > 0L)) &&
		(!isTRUE(beta_dyn$any) || (dynamic_beta_kind %in% c("ar1", "rw1") && !am_gaps)) &&
		is.null(dynamic_beta_per_actor) && nA >= 2L && nB >= 2L &&
		isTRUE(s2 == 1)
	if (!is.null(prior[["actor_move"]])) actor_move_bip <- actor_move_bip && isTRUE(prior[["actor_move"]])
	actor_move_bip_uv <- actor_move_bip && RA > 0L && RB > 0L && !isTRUE(dynamic_uv) && !isTRUE(.dyn_uv_has_rank)
	if (actor_move_bip) {
		am_cx_bip <- .lame_actor_move_context_bip(Y)
		am_every_bip <- max(1L, as.integer(prior[["actor_move_every"]] %||% 10L))
		# coefficient prior of rbeta_ab_bip_gibbs_cpp, N(0, g s2) each with the intercept first
		am_iV0_bip <- diag(length(beta)) / (as.numeric(g[1]) * s2)
		am_comp_bip <- isTRUE(intercept) && length(beta) > 0L
		# dynamic effects: the deviation paths of pairs of actors are drawn too, their sum held
		am_dev_bip <- isTRUE(dynamic_ab) && N >= 2L
		if (am_dev_bip) actor_move_stats <- c(actor_move_stats, dev_calls = 0, dev_evals = 0)
	}
	# position paths of a binary fit with moving positions: each actor's whole path moves with the latent array
	# integrated out (R/dynamic_position_move.R). prior$actor_move / prior$actor_move_paths switch it off,
	# prior$actor_move_every sets the sweeps between visits
	dpm_bin <- identical(family, "binary") && isTRUE(dynamic_uv) && isTRUE(.dyn_uv_has_rank) &&
		!identical(prior[["actor_move"]], FALSE) && !identical(prior[["actor_move_paths"]], FALSE) &&
		is.null(dynamic_beta_per_actor) &&
		dynamic_uv_kind %in% c("ar1", "t") && (identical(dynamic_uv_kind, "ar1") || (!bip && !symmetric)) &&
		(!symmetric || !is.null(L_dyn)) &&
		(if (bip) nA >= 2L && nB >= 2L && isTRUE(s2 == 1) else n >= 3L && !isTRUE(dynamic_rho) && !actor_slice_enabled)
	dpm_stats <- c(dpm_visits = 0, dpm_calls = 0, dpm_evals = 0, dpm_stuck = 0)
	if (dpm_bin) {
		dpm_cx <- .lame_dpm_binary_context(
			Y,
			if (bip) "bip" else if (symmetric) "sym" else "dir",
			if (bip && actor_move_bip) am_cx_bip else if (!bip && actor_move) am_cx
		)
		dpm_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 10L))
		dpm_uv_max <- as.numeric(prior$uv_max_abs %||% Inf)
		if (!is.finite(dpm_uv_max) || dpm_uv_max <= 0) dpm_uv_max <- Inf
	}
	# unit-variance gauge of a symmetric dynamic_uv fit (poisson and binary): exact draws along the scale ridge
	# between the paths and their weights, and of the persistence (R/dynamic_position_move.R)
	dpm_sym <- family %in% c("poisson", "binary") && symmetric && isTRUE(dynamic_uv) && isTRUE(.dyn_uv_has_rank) &&
		identical(dynamic_uv_kind, "ar1") && !is.null(L_dyn) &&
		!identical(prior[["actor_move"]], FALSE) && !identical(prior[["actor_move_paths"]], FALSE)
	dpm_sym_stats <- c(dpm_scale_attempts = 0, dpm_scale_accepts = 0)
	# translations of the additive levels with their latent rows and columns for the poisson and ordinal
	# families (R/latent_translation_move.R), which leave the Gaussian density of the latent array unchanged;
	# ordinal also draws the location of the latent array. prior$actor_move and prior$actor_move_every control them
	ltm_on <- family %in% c("poisson", "ordinal") && is.null(dynamic_beta_per_actor) &&
		!identical(prior[["actor_move"]], FALSE) &&
		(!isTRUE(beta_dyn$any) || (dynamic_beta_kind %in% c("ar1", "rw1") && !am_gaps)) &&
		(if (bip) nA >= 2L && nB >= 2L else n >= 3L) &&
		(identical(family, "poisson") ||
				(if (bip) (isTRUE(rvar) || isTRUE(cvar)) else (if (symmetric) isTRUE(nvar) else isTRUE(rvar) || isTRUE(cvar))))
	ltm_stats <- c(visits = 0, level_calls = 0, level_evals = 0, stuck = 0, skipped = 0, location_calls = 0)
	if (ltm_on) {
		ltm_cx <- .lame_ltm_context(
			Y,
			family,
			if (bip) "bip" else if (symmetric) "sym" else "dir",
			offset = if (poisson_exposure_active) log(period_exposure) else NULL,
			Y_int = if (use_explicit_cutpoints) Y_int else NULL
		)
		ltm_on <- isTRUE(ltm_cx$ok)
	}
	if (ltm_on) {
		# the poisson default visits every 2nd sweep, which the least active actors' levels need for mixing
		ltm_every <- max(
			1L,
			as.integer(
				prior[["actor_move_every"]] %||% (if (identical(family, "poisson")) 2L else if (spm_on) 1L else 10L)
			)
		)
		ltm_lev <- if (symmetric) rep(isTRUE(nvar), 2L) else c(isTRUE(rvar), isTRUE(cvar))
		# the coefficient prior precision (bipartite panels: 1 / (g s2), set at the visit)
		ltm_iV0 <- if (!bip && length(beta) > 0L && !isTRUE(beta_dyn$any)) apply(xxLong, c(1, 2), sum) / (N * as.numeric(g[1]))
		# location of the latent array: free under the rank likelihood, a strip under explicit cutpoints
		ltm_loc <- identical(family, "ordinal")
		# poisson: static coefficients, the dyadic variance and the correlation move with the latent array as well
		ltm_coef <- identical(family, "poisson") && length(beta) > 0L && !isTRUE(beta_dyn$any) &&
			all(vapply(Xlist, function(x) length(dim(x)) == 3L && dim(x)[3L] == length(beta), logical(1)))
		if (ltm_coef) {
			ltm_Xk <- lapply(
				seq_along(beta),
				function(k) array(
					vapply(
						seq_len(N),
						function(t) as.numeric(Xlist[[t]][, , k]),
						numeric(ltm_cx$nr * ltm_cx$nc)
					),
					c(ltm_cx$nr, ltm_cx$nc, N)
				)
			)
			ltm_coef <- all(vapply(ltm_Xk, function(x) all(is.finite(x)), logical(1)))
			ltm_Xz <- if (symmetric) lapply(ltm_Xk, .lame_ppm_sym_array) else ltm_Xk
		}
		ltm_scale <- identical(family, "poisson") && !isTRUE(dynamic_rho)
		# poisson: static positions move with the latent array as well
		ltm_uv <- identical(family, "poisson") && (if (bip) RA > 0L && RB > 0L else R > 0L) &&
			!isTRUE(dynamic_uv) && !isTRUE(.dyn_uv_has_rank)
		ltm_stats <- c(
			ltm_stats,
			coef_calls = 0,
			coef_evals = 0,
			coef_stuck = 0,
			scale_attempts = 0,
			scale_accepts = 0,
			rho_attempts = 0,
			rho_accepts = 0,
			uv_calls = 0,
			uv_evals = 0,
			uv_stuck = 0
		)
		# poisson with moving positions: each actor's whole position path moves with the latent array
		# (R/dynamic_position_move.R); a symmetric fit needs the weights of its product and has no t kind
		ltm_duv <- identical(family, "poisson") && isTRUE(dynamic_uv) && isTRUE(.dyn_uv_has_rank) &&
			!identical(prior[["actor_move_paths"]], FALSE) &&
			dynamic_uv_kind %in% c("ar1", "t") && (!symmetric || (identical(dynamic_uv_kind, "ar1") && !is.null(L_dyn)))
		if (ltm_duv) {
			ltm_uv_max <- as.numeric(prior$uv_max_abs %||% Inf)
			if (!is.finite(ltm_uv_max) || ltm_uv_max <= 0) ltm_uv_max <- Inf
			ltm_stats <- c(ltm_stats, dpm_calls = 0, dpm_evals = 0, dpm_stuck = 0)
		}
	}
	# translations of the actor-period effects of a poisson dynamic_ab fit with the latent array
	# (R/poisson_ab_period_move.R); prior$ab_period_move = FALSE omits them
	pab_on <- identical(family, "poisson") && isTRUE(dynamic_ab) && N >= 2L &&
		!identical(prior[["ab_period_move"]], FALSE) && (if (symmetric) isTRUE(nvar) else isTRUE(rvar) || isTRUE(cvar))
	pab_stats <- c(visits = 0, pairs = 0, evals = 0)
	if (pab_on) pab_cx <- .lame_pab_context(
		Y,
		if (bip) "bip" else if (symmetric) "sym" else "dir",
		offset = if (poisson_exposure_active) log(period_exposure) else NULL
	)
	ppm_on <- identical(family, "poisson") && isTRUE(beta_dyn$any) && N >= 2L &&
		!identical(prior[["poisson_period_move"]], FALSE)
	ppm_scale_on <- identical(family, "poisson") && !isTRUE(dynamic_rho)
	# the carried scale move takes the place of the plain residual scaling of the actor translations
	if (ltm_on && ppm_scale_on) ltm_scale <- FALSE
	ppm_static <- FALSE
	ppm_static_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 2L))
	ppm_stats <- c(
		visits = 0, period_calls = 0, period_evals = 0, path_calls = 0, path_evals = 0, stuck = 0,
		static_calls = 0, static_evals = 0, static_stuck = 0, scale_attempts = 0, scale_accepts = 0,
		scale_rho_attempts = 0, scale_rho_accepts = 0
	)
	# the context also serves the static coefficient translations of fits without dynamic_beta
	if (ppm_on || ppm_scale_on || (ltm_on && ltm_coef)) {
		ppm_cx <- .lame_ppm_context(
			Y,
			if (bip) "bip" else if (symmetric) "sym" else "dir",
			offset = if (poisson_exposure_active) log(period_exposure) else NULL,
			ltm = if (ltm_on) ltm_cx else NULL
		)
	}
	if (ppm_on) {
		ppm_tpos <- if (!is.null(time_index)) as.numeric(time_index) else seq_len(N)
		ppm_static <- length(beta_dyn$static_idx) > 0L && !identical(prior[["poisson_static_move"]], FALSE)
		if (ppm_static) {
			ppm_Xs <- lapply(
				seq_along(beta_dyn$static_idx),
				function(j)
					array(
						vapply(seq_len(N), function(t) Xstat_list_per_t[[t]][, j], numeric(ppm_cx$nr * ppm_cx$nc)),
						c(ppm_cx$nr, ppm_cx$nc, N)
					)
			)
			ppm_static <- all(vapply(ppm_Xs, function(x) all(is.finite(x)), logical(1)))
			ppm_Xzs <- if (symmetric) lapply(ppm_Xs, .lame_ppm_sym_array) else ppm_Xs
		}
		ppm_Xz <- if (symmetric) lapply(Xdyn_list_per_t, .lame_ppm_sym_cols, nr = n) else Xdyn_list_per_t
	}
	# exact draws of the scale of the rank-likelihood latent state, whole and per period (R/ordinal_scale_move.R);
	# prior$ordinal_scale_move = FALSE switches them off
	osm_on <- identical(family, "ordinal") && !identical(prior[["ordinal_scale_move"]], FALSE)
	osm_has_m <- R_eff > 0 || isTRUE(dynamic_ab) || per_actor_active
	# the centred additive effects live on a subspace with one dimension fewer per margin
	scale_cen <- if (isTRUE(intercept) && length(beta) > 0L && !isTRUE(dynamic_ab)) {
		if (symmetric) as.integer(isTRUE(nvar)) else as.integer(isTRUE(rvar)) + as.integer(isTRUE(cvar))
	} else 0L
	osm_stats <- c(visits = 0, orbit_moved = 0)
	# exact draw of the scale of a binary probit predictor with the latent array integrated out (R/binary_scale_orbit.R),
	# a direction the augmented sweep crosses slowly. prior$binary_scale_move and prior$binary_scale_every control it
	bso_on <- identical(family, "binary") && !isTRUE(dynamic_rho) && !identical(prior[["binary_scale_move"]], FALSE)
	bso_stats <- c(visits = 0, evals = 0, moved = 0)
	if (bso_on) {
		bso_every <- max(1L, as.integer(prior[["binary_scale_every"]] %||% 5L))
		bso_cx <- .lame_bso_context(Y, if (bip) "bip" else if (symmetric) "sym" else "dir")
		bso_P0 <- if (!bip && !isTRUE(beta_dyn$any) && length(beta)) .lame_replicated_prior(
			mXLong,
			xxLong,
			as.numeric(g[1L]),
			length(beta),
			N
		) else NULL
		bso_tpos <- if (!is.null(time_index)) as.numeric(time_index) else seq_len(N)
	}
	if (osm_on) {
		osm_mode <- if (bip) "bip" else if (symmetric) "sym" else "dir"
		osm_P0 <- if (!bip && !isTRUE(beta_dyn$any) && length(beta)) .lame_replicated_prior(
			mXLong,
			xxLong,
			as.numeric(g[1L]),
			length(beta),
			N
		) else NULL
		osm_tpos <- if (!is.null(time_index)) as.numeric(time_index) else seq_len(N)
	}
	# ridge moves of the dynamic-effect and bipartite samplers (R/dynamic_ab_hyper_move.R,
	# R/binary_bipartite_scale_moves.R): dynamic_ab hyperparameters with the levels integrated out, and the scale
	# split between bipartite positions and G. prior$ab_hyper_move and prior$uv_g_scale_move switch them off
	ab_hyper_move <- isTRUE(dynamic_ab) && N >= 2L && !identical(prior[["ab_hyper_move"]], FALSE) &&
		(if (bip) (isTRUE(rvar) || isTRUE(cvar)) && nA >= 2L && nB >= 2L
			else n >= 2L && (if (symmetric) isTRUE(nvar) else isTRUE(rvar) || isTRUE(cvar)))
	ab_hyper_par <- NULL
	if (ab_hyper_move) {
		ab_hyper_cx <- .lame_ab_hyper_context(
			a_mat,
			b_mat,
			if (symmetric) isTRUE(nvar) else isTRUE(rvar),
			isTRUE(cvar),
			symmetric,
			bip,
			isTRUE(intercept),
			prior,
			(prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
		)
		ab_hyper_tries <- max(1L, as.integer(prior[["ab_hyper_move_tries"]] %||% 8L))
	}
	uv_g_scale_bip <- bip && RA > 0L && RB > 0L && !isTRUE(dynamic_uv) &&
		!identical(prior[["uv_g_scale_move"]], FALSE)
	uv_g_scale_stats <- c(visits = 0, drawn = 0)
	# collapsed joint update for sparse panels with static factors and changing coefficients,
	# with the covariance prior integrated and orthonormal actor contrasts
	static_actor_marginal <- identical(family, "binary") && !bip && !symmetric &&
		isTRUE(intercept) && isTRUE(rvar) && isTRUE(cvar) && !isTRUE(dcor) &&
		!isTRUE(dynamic_ab) && !isTRUE(.dyn_uv_has_rank) && R %in% c(1L, 2L) &&
		isTRUE(beta_dyn$any) && dynamic_beta_kind %in% c("ar1", "rw1") &&
		is.null(dynamic_beta_per_actor) && is.null(time_index) && n >= 3L && N >= 2L &&
		length(beta_dyn$static_idx) > 0L && all(beta0_mean == 0) &&
		all(Lambda_full[row(Lambda_full) != col(Lambda_full)] == 0)
	marginal_move <- (static_path_move && level_move && n >= 3L && N >= 2L &&
			length(beta_dyn$static_idx) > 0L) || static_actor_marginal
	if (!is.null(prior[["marginal_move"]])) marginal_move <- marginal_move && prior[["marginal_move"]]
	marginal_att <- 0L
	marginal_acc <- 0L
	marginal_bad <- 0L
	if (marginal_move) {
		marginal_sd <- if (is.null(marginal_sd_start)) 0.12 else marginal_sd_start
		if (!is.numeric(marginal_sd) || length(marginal_sd) != 1L ||
				!is.finite(marginal_sd) || marginal_sd <= 0 || marginal_sd > .5)
			cli::cli_abort("The marginal move step size must be finite and in (0, 0.5].")
		marginal_context <- if (static_actor_marginal)
			.lame_binary_static_hmc_context(Y, Xlist, beta_dyn$static_idx, beta_dyn$dynamic_idx, R, uv_prior_static)
		else .lame_binary_hmc_context(Y, Xlist, beta_dyn$static_idx, beta_dyn$dynamic_idx, R, uv_prior_static)
		if (static_actor_marginal) {
			mdg <- array(rep(diag(TRUE, n), N), dim(Y))
			lev_i1 <- which(Y == 1 & !mdg); lev_i0 <- which(Y == 0 & !mdg)
			lev_ina <- which(is.na(Y) & !mdg); lev_idg <- which(mdg)
		}
	}
	# joint Gaussian block removing the coupling between coefficients and actor levels
	joint_levels_enabled <- identical(family, "binary") && !bip && !symmetric &&
		isTRUE(dynamic_ab) && !isTRUE(beta_dyn$any) && R == 0L &&
		isTRUE(intercept) && isTRUE(rvar) && isTRUE(cvar) &&
		!isTRUE(dynamic_rho) && is.null(dynamic_beta_per_actor) && is.null(time_index)
	joint_levels_context <- if (joint_levels_enabled)
		.lame_joint_levels_context(Xlist, apply(xxLong, c(1L, 2L), sum) / (N * g[1L])) else NULL
	joint_levels_enabled <- joint_levels_enabled && !is.null(joint_levels_context)
	joint_levels_attempts <- 0L
	if (joint_levels_enabled) tryErrorChecks$joint_levels <- 0L
	# binary symmetric dyads have one unit-variance latent observation; use that exact Gaussian law
	sym_gaussian_context <- if (family %in% c("binary", "ordinal") && !bip &&
			isTRUE(symmetric) && !isTRUE(beta_dyn$any)) {
		q <- if (length(beta)) apply(xxLong, c(1L, 2L), sum) / (N * g[1L]) else matrix(0, 0L, 0L)
		.lame_symmetric_gaussian_context(
			Xlist,
			q,
			actor = !isTRUE(dynamic_ab) && isTRUE(nvar),
			center = isTRUE(intercept)
		)
	} else NULL
	# restore after all initializers so the saved state enters the first transition exactly
	if (!is.null(checkpoint_start)) {
		.lame_checkpoint_restore(checkpoint_start, environment(), checkpoint_layout)
		if (dhmc_enabled && (!is.numeric(dhmc_sd) || length(dhmc_sd) != 1L || !is.finite(dhmc_sd) || dhmc_sd < .001 || dhmc_sd > .5)) cli::cli_abort("Invalid checkpoint dynamic marginal move step size.")
		if (marginal_move && (!is.numeric(marginal_sd) || length(marginal_sd) != 1L ||
				!is.finite(marginal_sd) || marginal_sd <= 0 || marginal_sd > .5))
			cli::cli_abort("The checkpoint marginal move step size must be finite and in (0, 0.5].")
	}
	.lame_explicit_ordinal_support(
		family,
		ordinal_cutpoints,
		bip,
		symmetric,
		dcor || isTRUE(dynamic_rho),
		rho
	)
	# stationary-variance limits are joint prior support; start inside it, and a resumed state must already satisfy it
	cap_unit_start <- .lame_cap_unit(family, s2)
	if (dynamic_ab) {
		ab_start <- .lame_ar_start_support(
			rho_ab,
			sigma_ab,
			(prior$ab_var_cap %||% Inf) * cap_unit_start,
			resume = !is.null(checkpoint_start)
		)
		rho_ab <- ab_start$rho
		sigma_ab <- ab_start$sigma
	}
	if (.dyn_uv_has_rank && is.null(L_dyn)) {
		uv_start <- .lame_ar_start_support(
			rho_uv,
			sigma_uv,
			(prior$uv_var_cap %||% Inf) * cap_unit_start,
			sigma_max = prior$sigma_uv_max %||% Inf,
			resume = !is.null(checkpoint_start)
		)
		rho_uv <- uv_start$rho
		sigma_uv <- uv_start$sigma
	}
	if (.dyn_uv_has_rank && !is.null(L_dyn) &&
			is.finite(prior$uv_var_cap %||% Inf) && prior$uv_var_cap * cap_unit_start < 1) {
		cli::cli_abort(
			"The stationary-variance cap excludes the symmetric unit-variance position model at its starting residual scale.", call = NULL
		)
	}
	if (.dyn_uv_has_rank) {
		uv_start_bound <- as.numeric(prior$uv_max_abs %||% Inf)
		if (is.finite(uv_start_bound) && uv_start_bound > 0) {
			free_max <- if (symmetric) max(abs(U_cube)) else max(abs(U_cube), abs(V_cube))
			if (free_max > uv_start_bound) {
				if (!is.null(checkpoint_start)) cli::cli_abort(
					"The checkpoint coordinates exceed prior$uv_max_abs.", call = NULL
				)
				# initial values are arbitrary; put them strictly inside support
				U_cube <- U_cube * (0.5 * uv_start_bound / free_max)
				if (symmetric && !is.null(L_dyn)) {
					for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
				} else V_cube <- V_cube * (0.5 * uv_start_bound / free_max)
				U <- apply(U_cube, c(1, 2), mean)
				V <- apply(V_cube, c(1, 2), mean)
			}
		}
	}
	for (s in 1:(nscan + burn)) {
		# stop cleanly once past max_seconds, writing a final checkpoint so a resumer has a state
		if (is.finite(max_seconds)) {
			elapsed <- as.numeric(difftime(Sys.time(), .lame_start_time, units = "secs"))
			if (elapsed > max_seconds) {
				if (verbose) cli::cli_inform(c(
					"i" = "Reached {.arg max_seconds} = {.val {max_seconds}}s at iteration {s} (of {nscan + burn}); terminating early."
				))
				.lame_terminated_early <- TRUE
				if (!is.null(checkpoint_path)) {
					.ckpt <- list(
						iter = s - 1L,
						burn_remaining = max(0L, burn - s + 1L),
						BETA = BETA,
						VC = VC,
						seed = seed,
						start_vals = .lame_snapshot_state(),
						call_args = mc,
						resolved_prior = prior,
						requested_prior = .lame_input_prior,
						timestamp = Sys.time()
					)
					try(saveRDS(.ckpt, file = checkpoint_path), silent = TRUE)
				}
				break
			}
		}
		# periodic checkpoint for lame_resume(path): stored draws, seed and full sampler state
		if (!is.null(checkpoint_path) && s > 1L && (s %% checkpoint_every == 0L)) {
			.ckpt <- list(
				iter = s - 1L,
				burn_remaining = max(0L, burn - s + 1L),
				BETA = BETA,
				VC = VC,
				seed = seed,
				start_vals = .lame_snapshot_state(),
				call_args = mc,
				resolved_prior = prior,
				requested_prior = .lame_input_prior,
				timestamp = Sys.time()
			)
			# best-effort write; on failure the MCMC continues
			try(saveRDS(.ckpt, file = checkpoint_path), silent = TRUE)
		}
		# symmetric dyads are counted once: as in ame_unipartite, rho ramps to min(.9999, 1 - 1/sqrt(s)) so the
		# pair-whitened GLS collapses each mirrored pair into one observation. the ramp runs on the cumulative
		# sweep count, carried by start_vals and checkpoints
		if (symmetric && !bip) {
			rho <- min(.9999, 1 - 1 / sqrt(sym_ramp_done + 1))
		}
		# update z
		Xlist_tmp <- Xlist
		# beta_full_path and the X beta cube offset for the dynamic branch (NULL when dynamic_beta is off)
		beta_full_path_iter <- .lame_beta_full_path()
		Xbeta_cube_iter <- .lame_xbeta_cube(beta_full_path_iter)
		# with dynamic_beta the EZ helpers see a zero beta and Xbeta_cube_iter is added instead
		beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta
		if (dynamic_ab) {
			# EZ with time-varying additive effects
			EZ <- get_EZ_dynamic_ab(
				Xlist_tmp,
				beta_for_EZ,
				a_mat,
				b_mat,
				if (.dyn_uv_has_rank) U_cube else U,
				if (.dyn_uv_has_rank) V_cube else V,
				N,
				bip = bip,
				G = G,
				nA = nA,
				nB = nB
			)
		} else {
			if (bip) {
				if (length(a) == 1) a <- rep(a, nA)
				if (length(b) == 1) b <- rep(b, nB)
				a_mat <- matrix(a[1:nA], nA, N)
				b_mat <- matrix(b[1:nB], nB, N)
				if (is.null(dim(U)) || length(dim(U)) == 2) {
					U_3d <- array(0, dim = c(nA, max(1, ncol(U)), N))
					for (t in 1:N) {
						U_3d[, , t] <- U
					}
					U <- U_3d
				}

				if (is.null(dim(V)) || length(dim(V)) == 2) {
					V_3d <- array(0, dim = c(nB, max(1, ncol(V)), N))
					for (t in 1:N) {
						V_3d[, , t] <- V
					}
					V <- V_3d
				}

				base_cube <- array(0, dim = c(nA, nB, N))
				if (!beta_dyn$any) {
					for (t in 1:N) {
						if (length(beta) > 0) {
							X_dims <- dim(Xlist_tmp[[t]])
							if (X_dims[1] != nA || X_dims[2] != nB) {
								cli::cli_abort(paste(
									"Bipartite X dimension mismatch at time",
									t,
									": expected",
									nA,
									"x",
									nB,
									"but got",
									X_dims[1],
									"x",
									X_dims[2]
								))
							}
							# sum X beta over the covariates
							base_mat <- matrix(0, nA, nB)
							n_covs <- min(length(beta), X_dims[3])
							for (i in 1:n_covs) {
								base_mat <- base_mat + beta[i] * Xlist_tmp[[t]][, , i]
							}
							base_cube[, , t] <- base_mat
						}
					}
				}
				# with dynamic_beta, X beta is in Xbeta_cube_iter and added later; the predictor carries each period's
				# own positions
				EZ <- get_EZ_bip_cpp(
					base_cube,
					a_mat,
					b_mat,
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					G
				)
			} else {
				ab_mat <- outer(a, b, "+")
				storage.mode(ab_mat) <- "double"
				EZ <- get_EZ_unip_time(
					Xlist_tmp,
					beta_for_EZ,
					ab_mat,
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N
				)
			}
		}
		# add the per-period X beta contribution when dynamic_beta is on
		if (beta_dyn$any && !is.null(Xbeta_cube_iter)) {
			EZ <- EZ + Xbeta_cube_iter
		}
		if (per_actor_active) {
			per_actor_offset <- .lame_per_actor_offset(theta_actor, per_actor_X, dynamic_beta_per_actor)
			EZ <- EZ + per_actor_offset
		}
		# batch C++ z sampling where possible
		if (family == "normal") {
			if (isTRUE(dynamic_rho)) {
				for (t in 1:N) {
					rho_t <- .lame_rho_at(rho, t)
					Z[, , t] <- rZ_nrm_fc(Z[, , t], EZ[, , t], rho_t, s2, Y[, , t])
					E.nrm[, , t] <- Z[, , t] - EZ[, , t]
				}
			} else {
				# single C++ call for all periods
				Z_batch <- try(rZ_nrm_batch_cpp(Z, EZ, rho, s2, Y), silent = TRUE)
				if (!inherits(Z_batch, 'try-error')) {
					Z <- Z_batch$Z
					E.nrm <- Z_batch$E_nrm
				} else {
					# fallback: per-period R sampling
					for (t in 1:N) { Z[, , t] <- rZ_nrm_fc(
						Z[, , t],
						EZ[, , t],
						rho,
						s2,
						Y[, , t]
					); E.nrm[, , t] <- Z[, , t] - EZ[, , t] }
					tryErrorChecks$Z <- tryErrorChecks$Z + 1
				}
			}
		} else if (family == "binary" && bip && rho == 0) {
			# bipartite binary z via the batch C++ path when available
			Z_new <- tryCatch(
				rZ_bin_bip_batch_cpp(Z, EZ, Y),
				error = function(e) { tryErrorChecks$Z <<- tryErrorChecks$Z + 1; NULL }
			)
			if (!is.null(Z_new)) {
				Z <- Z_new
			} else {
				for (t in 1:N) Z[, , t] <- rZ_bin_fc(Z[, , t], EZ[, , t], 0, Y[, , t])
			}
			E.nrm <- Z - EZ
		} else {
			# per-period loop for the remaining families; bipartite slices use the rectangular `_bip_fc` samplers
			# (R/rZ_bipartite.R), with rho = 0 since a bipartite cell has no reciprocal partner
			for (t in 1:N) {
				if (family == "binary") {
					if (symmetric) {
						# one latent per unordered pair
						Z[, , t] <- rZ_bin_sym_fc(Z[, , t], EZ[, , t], Y[, , t], .lame_rho_scalar(rho))
					} else {
						# fused C++ kernel (same conditional and sweep order as rZ_bin_fc), with the R sampler as fallback
						z_new <- tryCatch(
							rZ_bin_fused_cpp(Z[, , t], EZ[, , t], rho, Y[, , t]),
							error = function(e) { tryErrorChecks$Z <<- tryErrorChecks$Z + 1; NULL }
						)
						Z[, , t] <- if (!is.null(z_new)) z_new else rZ_bin_fc(Z[, , t], EZ[, , t], rho, Y[, , t])
					}
				}
				if (family == "ordinal") {
					if (bip) {
						Z[, , t] <- rZ_ord_bip_fc(Z[, , t], EZ[, , t], Y[, , t])
					} else if (use_explicit_cutpoints) {
						# explicit cutpoints (Cowles 1996): truncated-normal z given the current alpha
						if (symmetric) {
							Z[, , t] <- rZ_ord_sym_explicit_fc(
								Z[, , t],
								EZ[, , t],
								Y_int[, , t],
								alpha
							)
						} else {
							Z[, , t] <- rZ_ord_explicit_fc(
								Z[, , t],
								EZ[, , t],
								Y_int[, , t],
								alpha,
								rho = if (dcor) .lame_rho_scalar(rho) else 0
							)
						}
					} else if (symmetric) {
						# one latent per unordered pair, as in ame()
						Z[, , t] <- rZ_ord_sym_fc(Z[, , t], EZ[, , t], Y[, , t], .lame_rho_scalar(rho))
					} else {
						Z[, , t] <- rZ_ord_fc(Z[, , t], EZ[, , t], rho, Y[, , t])
					}
				}
				if (family == "cbin") {
					if (bip) {
						# per-period row count of nominations made (see .odobs_room)
						Z[, , t] <- rZ_cbin_bip_fc(Z[, , t], EZ[, , t], Y[, , t], odmax, .odobs_room(t))
					} else {
						Z[, , t] <- rZ_cbin_fc(Z[, , t], EZ[, , t], rho, Y[, , t], odmax, .odobs_room(t))
					}
				}
				if (family == "frn") {
					if (bip) {
						Z[, , t] <- rZ_frn_bip_fc(Z[, , t], EZ[, , t], Y[, , t], YL[[t]], odmax, .odobs_room(t))
					} else {
						Z[, , t] <- rZ_frn_fc(Z[, , t], EZ[, , t], rho, Y[, , t], YL[[t]], odmax, .odobs_room(t))
					}
				}
				if (family == "poisson") {
					if (bip) {
						# rectangular MH on z, rho = 0, with the per-period exposure offset
						log_exp <- if (poisson_exposure_active) log(period_exposure[t]) else 0
						Z[, , t] <- rZ_pois_bip_fc(
							Z[, , t],
							EZ[, , t],
							s2,
							Y[, , t],
							log_exposure = log_exp
						)
					} else if (symmetric) {
						# one latent per unordered pair
						Z[, , t] <- rZ_pois_sym_fc(
							Z[, , t],
							EZ[, , t],
							s2,
							Y[, , t],
							log_exposure = if (poisson_exposure_active) log(period_exposure[t]) else 0,
							rho = .lame_rho_scalar(rho)
						)
					} else if (poisson_exposure_active) {
						# unipartite poisson with a per-period exposure offset
						Z[, , t] <- suppressWarnings(rZ_pois_fc_exposure_cpp(
							Z[, , t],
							EZ[, , t],
							rho,
							s2,
							Y[, , t],
							log_exposure = log(period_exposure[t])
						))
					} else {
						Z[, , t] <- rZ_pois_fc(Z[, , t], EZ[, , t], rho, s2, Y[, , t])
					}
					E.nrm[, , t] <- Z[, , t] - EZ[, , t]
				}
			}
		}

		# the state EZ describes, for the collapsed translations below
		if (nom_room) nom_top <- list(
			a = if (dynamic_ab) a_mat else a,
			b = if (dynamic_ab) b_mat else b,
			path = beta_dyn_path[, nom_room_k]
		)
		# cbin / frn: each sender's level moves with its latent row (.lame_nom_row_shift); under dynamic_ab the
		# level and its whole path move, the deviations untouched
		if (nom_shift && s %% nom_shift_every == 0L) {
			nom_sh <- .lame_nom_row_shift(
				Z,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				Sab,
				odmax,
				nom_shift_masks,
				bip,
				nom_shift_rows,
				center = isTRUE(intercept)
			)
			Z <- nom_sh$Z
			if (dynamic_ab) {
				a_mat <- a_mat + (nom_sh$a - a_lev)
				a_lev <- nom_sh$a
			}
			a <- nom_sh$a
		}
		# the intercept of a dynamic_beta fit moves with the levels as a path: the two moves below read its prior
		# from nom_ic, and it is written back to the coefficient state once
		nom_ic <- NULL
		nom_ic_moved <- FALSE
		if (isTRUE(intercept) && isTRUE(beta_dyn$any) &&
				((nom_lvl && (s %% nom_lvl_every == 0L)) || (nom_col && (s %% nom_col_every == 0L)))) {
			nom_ic <- .lame_nom_icpt_state(
				beta_dyn,
				beta_dyn_path,
				beta_static,
				rho_beta_by_group,
				sigma_beta_by_group,
				Lambda_full,
				group_id,
				beta0_mean,
				beta0_cov,
				g,
				s2
			)
		}
		# the level of the rows free above moves with the intercept against the rows with room (.lame_nom_level_shift)
		if (nom_lvl && (s %% nom_lvl_every == 0L)) {
			nom_lv <- .lame_nom_level_shift(
				Z,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				Sab,
				odmax,
				nom_shift_masks,
				bip,
				nom_shift_rows,
				beta,
				if (bip) diag(1 / (as.numeric(g[1]) * s2), length(beta)) else nom_lvl_iV0,
				ic = nom_ic
			)
			Z <- nom_lv$Z
			if (dynamic_ab) {
				a_mat <- a_mat + (nom_lv$a - a_lev)
				a_lev <- nom_lv$a
			}
			a <- nom_lv$a
			beta <- nom_lv$beta
			nom_ic <- nom_lv$ic
			nom_ic_moved <- nom_ic_moved || nom_lv$calls > 0
			nom_lvl_stats <- nom_lvl_stats + c(1, nom_lv$calls, nom_lv$free)
		}
		# the intercept of each period moves with the rows it does not tie to the others (.lame_nom_period_shift)
		if (nom_per && (s %% nom_lvl_every == 0L) && !is.null(nom_ic) && identical(nom_ic$kind, "path")) {
			nom_pm <- .lame_nom_period_shift(
				Z,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				Sab,
				odmax,
				nom_shift_masks,
				bip,
				nom_shift_rows,
				nom_ic
			)
			Z <- nom_pm$Z
			if (dynamic_ab) {
				a_mat <- a_mat + (nom_pm$a - a_lev)
				a_lev <- nom_pm$a
			}
			a <- nom_pm$a
			nom_ic <- nom_pm$ic
			nom_ic_moved <- nom_ic_moved || nom_pm$calls > 0
			nom_per_stats <- nom_per_stats + c(1, nom_pm$calls, nom_pm$comps)
		}
		# each receiver's level moves with its latent column (.lame_nom_col_shift)
		if (nom_col && (s %% nom_col_every == 0L)) {
			nom_cs <- .lame_nom_col_shift(
				Z,
				Y,
				family,
				if (identical(family, "frn")) YL else NULL,
				odmax,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				Sab,
				bip,
				nom_col_ok,
				center = isTRUE(intercept),
				beta = beta,
				iV0 = if (isTRUE(intercept) && bip) diag(1 / (as.numeric(g[1]) * s2), length(beta)) else nom_col_iV0,
				ic = nom_ic
			)
			Z <- nom_cs$Z
			if (dynamic_ab) {
				b_mat <- b_mat + (nom_cs$b - b_lev)
				b_lev <- nom_cs$b
			}
			b <- nom_cs$b
			if (isTRUE(intercept)) beta <- nom_cs$beta
			nom_ic <- nom_cs$ic
			nom_ic_moved <- nom_ic_moved || nom_cs$calls > 0
			nom_col_stats <- nom_col_stats + c(1, nom_cs$calls, nom_cs$skips)
		}
		if (nom_ic_moved && !is.null(nom_ic)) {
			if (identical(nom_ic$kind, "path")) {
				beta_dyn_path[, nom_ic$k] <- nom_ic$path
				beta[1L] <- mean(nom_ic$path)
			} else {
				beta_static[nom_ic$k] <- nom_ic$value
				beta[1L] <- nom_ic$value
			}
		}
		# collapsed translation of the periods with rows with room (.lame_nom_room_period), starting from the predictor
		# of the latent sweep plus what the translations since moved
		if (nom_room && (s %% nom_lvl_every == 0L)) {
			nom_pic <- .lame_nom_icpt_state(
				beta_dyn,
				beta_dyn_path,
				beta_static,
				rho_beta_by_group,
				sigma_beta_by_group,
				Lambda_full,
				group_id,
				beta0_mean,
				beta0_cov,
				g,
				s2
			)
			if (identical(nom_pic$kind, "path")) {
				nom_dpath <- nom_pic$path - nom_top$path
				nom_ez <- EZ
				if (dynamic_ab) {
					nom_da <- a_mat - nom_top$a
					nom_db <- b_mat - nom_top$b
					for (t in seq_len(N)) nom_ez[, , t] <- EZ[, , t] + outer(nom_da[, t], nom_db[, t], "+") + nom_dpath[t]
				} else {
					nom_dab <- outer(a - nom_top$a, b - nom_top$b, "+")
					for (t in seq_len(N)) nom_ez[, , t] <- EZ[, , t] + nom_dab + nom_dpath[t]
				}
				nom_rho_t <- if (!bip && !symmetric && isTRUE(dcor)) vapply(
					seq_len(N),
					function(t) .lame_rho_at(rho, t),
					0
				) else rep(0, N)
				nom_rm <- .lame_nom_room_period(Z, nom_ez, nom_room_ctx, bip, nom_rho_t, nom_pic)
				Z <- nom_rm$Z
				beta_dyn_path[, nom_pic$k] <- nom_rm$ic$path
				beta[1L] <- mean(nom_rm$ic$path)
				nom_room_stats <- nom_room_stats + c(1, nom_rm$stats[c("attempts", "accepts", "stuck")])
			}
		}

		# joint scale move: one visit draws the shared dyad factors, walks two ridges and redraws z once.
		# level: the predictor times c (coefficients and additive effects by c, positions by sqrt(c));
		# deviation: the additive deviations and their innovation scale times d. other variance parameters are held fixed
		if (level_move && (s %% lev_every == 0L)) {
			rho_lev <- .lame_rho_scalar(rho)
			lev_off <- .lame_dyad_factor_draw(Z, EZ, rho_lev, lev_upper)
			sd_cell <- sqrt(1 - abs(rho_lev))
			ez_o <- EZ[lev_obs] / sd_cell
			off_o <- lev_off[lev_obs] / sd_cell
			ll_cur <- sum(pnorm(lev_sY * (ez_o + off_o), log.p = TRUE))
			# quadratic forms of the blocks the level direction moves, at the current state
			if (dynamic_ab) {
				a_dev_now <- a_mat - a_lev
				b_dev_now <- b_mat - b_lev
				# the levels are jointly normal with covariance Sab, so the quadratic form is the full bilinear one
				iSab_lev <- try(solve(Sab), silent = TRUE)
				if (inherits(iSab_lev, "try-error") || !all(is.finite(iSab_lev))) {
					stop("Unable to resolve the joint level-move covariance precision.")
				}
				AB_lev <- cbind(a_lev, b_lev)
				q_lev <- sum(AB_lev * (AB_lev %*% iSab_lev))
				n_lev <- length(a_lev) + length(b_lev)
				rm(iSab_lev, AB_lev)
				dss <- .lame_ab_dev_ss(a_dev_now, b_dev_now, rho_ab)
				q_dev <- dss$ss / sigma_ab^2
				n_dev <- dss$cnt
			} else {
				iSab <- try(solve(Sab), silent = TRUE)
				if (inherits(iSab, "try-error") || !all(is.finite(iSab))) {
					stop("Unable to resolve the joint level-move covariance precision.")
				}
				AB <- cbind(a, b)
				q_lev <- sum(AB * (AB %*% iSab))
				n_lev <- length(a) + length(b)
				q_dev <- 0; n_dev <- 0L
			}
			if (static_path_move || static_uv_move) {
				uv_eigen <- eigen(crossprod(cbind(U, V)), symmetric = TRUE, only.values = TRUE)$values
				uv_eigen <- pmax(uv_eigen, 0)
				q_uv <- sum(uv_eigen) / n; n_uv <- length(U) + length(V)
			}
			if (static_path_move) {
				lev_beta_terms <- .lame_path_scale_terms(
					beta_static,
					beta_dyn_path,
					rho_beta_by_group[group_id],
					sigma_beta_by_group[group_id],
					Lambda_full,
					beta0_cov,
					.lame_static_g_scale(g, beta_dyn$static_idx),
					s2
				)
			} else if (static_uv_move) {
				# the static positions enter through their matrix-t term below
			} else if (lev_uv) {
				uss <- .lame_uv_ss(U_cube, V_cube, rho_uv)
				q_uv <- uss$ss / sigma_uv^2
				n_uv <- uss$cnt
			} else {
				q_uv <- 0; n_uv <- 0L
			}
			# keep the step inside the coordinate rail prior$uv_max_abs, which the UV draw enforces
			uv_room <- if (lev_uv && !static_path_move && !static_uv_move && is.finite(lev_uv_max)) {
				(lev_uv_max / max(abs(c(U_cube, V_cube)), 1e-12))^2
			} else Inf
			ctot <- 1
			for (k_lev in seq_len(lev_tries)) {
				cc <- ctot * exp(rnorm(1, 0, lev_sd))
				lev_att <- lev_att + 1L
				lev_ok <- FALSE
				if (cc < uv_room && is.finite(cc) && cc > 0) {
					ll_new <- sum(pnorm(lev_sY * (cc * ez_o + off_o), log.p = TRUE))
					lr <- ll_new - ll_cur + if (static_path_move) {
						.lame_static_scale_logratio(
							cc / ctot,
							ctot^2 * (lev_beta_terms$q + q_lev + q_dev),
							lev_beta_terms$dimension + n_lev - lev_cen_lev + max(n_dev - lev_cen_dev, 0L),
							ctot * uv_eigen,
							n,
							R,
							s > .5 * burn,
							uv_prior_static$df,
							uv_prior_static$psi
						)
					} else if (static_uv_move) {
						.lame_scale_prior_logratio(
							cc / ctot,
							ctot * beta,
							lev_iV0,
							ctot * ctot * q_lev,
							n_lev - lev_cen_lev,
							ctot * ctot * q_dev,
							max(n_dev - lev_cen_dev, 0L)
						) +
							.lame_static_scale_logratio(
								cc / ctot,
								0,
								0,
								ctot * uv_eigen,
								n,
								R,
								s > .5 * burn,
								uv_prior_static$df,
								uv_prior_static$psi
							)
					} else .lame_scale_prior_logratio(
						cc / ctot,
						ctot * beta,
						lev_iV0,
						ctot * ctot * q_lev,
						n_lev - lev_cen_lev,
						ctot * ctot * q_dev,
						max(n_dev - lev_cen_dev, 0L),
						ctot * q_uv,
						max(n_uv - lev_cen_uv, 0L)
					)
					lev_ok <- is.finite(lr) && log(runif(1)) < lr
					if (lev_ok) {
						ctot <- cc
						ll_cur <- ll_new
						lev_acc <- lev_acc + 1L
					}
				}
				# Robbins-Monro step adaptation in burn-in only, so the kernel is fixed afterwards
				if (s <= burn) {
					lev_sd <- min(
						max(
							exp(log(lev_sd) +
									(as.numeric(lev_ok) - 0.35) / sqrt(lev_att)),
							1e-5
						),
						1
					)
				}
			}
			# deviation direction: scaling the deviations with sigma_ab leaves their AR(1) prior invariant, so only the
			# likelihood and the inverse-gamma prior on sigma_ab^2 remain
			dtot <- 1
			if (dev_move && n_dev > 0L) {
				dev_o <- (a_dev_now[dev_ia] + b_dev_now[dev_ib]) / sd_cell
				w_o <- ctot * ez_o
				dev_cur <- ctot * dev_o
				v_ab <- sigma_ab * sigma_ab
				ab_v_max <- (prior$ab_var_cap %||% Inf) *
					.lame_cap_unit(family, s2) * (1 - min(rho_ab * rho_ab, 0.9801))
				for (k_dev in seq_len(dev_tries)) {
					dd <- dtot * exp(rnorm(1, 0, dev_sd))
					dev_att <- dev_att + 1L
					dev_ok <- FALSE
					if (is.finite(dd) && dd > 0) {
						ll_new <- sum(pnorm(
							lev_sY *
								(w_o + (dd / dtot - 1) * dev_cur + off_o),
							log.p = TRUE
						))
						# the centred Gaussian density and the state Jacobian cancel, leaving the variance prior, its Jacobian
						# and the likelihood ratio
						lr <- ll_new - ll_cur + .lame_scale_ig_logratio(
							dtot * dtot * v_ab,
							dd * dd * v_ab,
							prior$sigma_ab_shape %||% 2,
							prior$sigma_ab_scale,
							ab_v_max
						)
						dev_ok <- is.finite(lr) && log(runif(1)) < lr
						if (dev_ok) {
							w_o <- w_o + (dd / dtot - 1) * dev_cur
							dev_cur <- (dd / dtot) * dev_cur
							dtot <- dd
							ll_cur <- ll_new
							dev_acc <- dev_acc + 1L
						}
					}
					if (s <= burn) {
						dev_sd <- min(
							max(
								exp(log(dev_sd) +
										(as.numeric(dev_ok) - 0.35) / sqrt(dev_att)),
								1e-5
							),
							1
						)
					}
				}
				rm(dev_o, w_o, dev_cur)
			}
			if (ctot != 1 || dtot != 1) {
				beta <- ctot * beta
				if (static_path_move) {
					beta_dyn_path <- ctot * beta_dyn_path
					beta_static <- ctot * beta_static
				}
				a <- ctot * a
				b <- ctot * b
				if (dynamic_ab) {
					a_lev <- ctot * a_lev
					b_lev <- ctot * b_lev
					a_mat <- a_lev + (ctot * dtot) * a_dev_now
					b_mat <- b_lev + (ctot * dtot) * b_dev_now
					if (dtot != 1 && !is.null(sigma_ab) && is.finite(sigma_ab)) {
						sigma_ab <- dtot * sigma_ab
					}
				}
				if (lev_uv && ctot != 1) {
					rc <- sqrt(ctot)
					if (!static_path_move && !static_uv_move) {
						U_cube <- rc * U_cube
						V_cube <- rc * V_cube
					}
					U <- rc * U
					V <- rc * V
				}
				# EZ is linear in the effect paths: c EZ moves the whole predictor and c (d - 1) the rest of the deviation part
				EZ <- ctot * EZ
				if (dtot != 1) {
					for (t_dev in seq_len(N)) {
						EZ[, , t_dev] <- EZ[, , t_dev] + (ctot * (dtot - 1)) *
							outer(a_dev_now[, t_dev], b_dev_now[, t_dev], "+")
					}
				}
				Z <- .lame_binary_z_given_factor(
					Z,
					EZ,
					lev_off,
					rho_lev,
					lev_i1,
					lev_i0,
					lev_ina,
					lev_idg
				)
			}
			if (dynamic_ab) rm(a_dev_now, b_dev_now)
			rm(lev_off, ez_o, off_o)
		}
		# bipartite counterpart: the whole predictor times a factor with the latent array integrated out,
		# then the latent array is redrawn
		if (bip_level_move && (s %% bpl_every == 0L)) {
			bpl_G <- !is.null(G) && RA > 0L && RB > 0L
			bpl_ez <- EZ[bpl_cx$idx]
			if (dynamic_ab) {
				bpl_adev <- a_mat - a_lev
				bpl_bdev <- b_mat - b_lev
				bpl_dss <- .lame_ab_dev_ss(bpl_adev, bpl_bdev, rho_ab)
				bpl_ql <- sum(a_lev^2) / Sab[1L, 1L] + sum(b_lev^2) / Sab[2L, 2L]
				bpl_qd <- bpl_dss$ss / sigma_ab^2
				bpl_nd <- max(bpl_dss$cnt - 2L * N, 0L)
				rm(bpl_adev, bpl_bdev, bpl_dss)
			} else {
				bpl_ql <- sum(a^2) / Sab[1L, 1L] + sum(b^2) / Sab[2L, 2L]
				bpl_qd <- 0
				bpl_nd <- 0L
			}
			bpl_o <- .lame_bpl_scale(
				bpl_ez,
				bpl_cx$sY,
				beta,
				bpl_iV0,
				bpl_ql,
				nA + nB - 2L,
				bpl_qd,
				bpl_nd,
				if (bpl_G) sum(G^2) else 0,
				if (bpl_G) length(G) else 0L,
				room = if (bpl_G && is.finite(
					prior$G_max %||% g_rail_default
				) && max(abs(G)) > 0) (prior$G_max %||% g_rail_default) / max(abs(G)) else Inf,
				tries = bpl_tries,
				sd = lev_sd,
				adapt = s <= burn,
				att = lev_att
			)
			lev_att <- bpl_o$att
			lev_acc <- lev_acc + bpl_o$acc
			lev_sd <- bpl_o$sd
			if (bpl_o$c != 1) {
				beta <- bpl_o$c * beta
				a <- bpl_o$c * a
				b <- bpl_o$c * b
				if (dynamic_ab) {
					a_lev <- bpl_o$c * a_lev
					b_lev <- bpl_o$c * b_lev
					a_mat <- bpl_o$c * a_mat
					b_mat <- bpl_o$c * b_mat
				}
				if (bpl_G) G <- bpl_o$c * G
				EZ <- bpl_o$c * EZ
				Z <- rZ_bin_bip_batch_cpp(Z, EZ, Y)
				E.nrm <- Z - EZ
			}
			rm(bpl_ez, bpl_ql, bpl_qd, bpl_nd, bpl_o, bpl_G)
		}
		# collapsed per-actor update: exact slice draws of level translations and static positions with the
		# latent array integrated out (with no additive margin, positions only after the first half of burn-in)
		if (actor_move && (s %% am_every == 0L) && (any(am_lev) || s > .5 * burn)) {
			# prior of the intercept of a dynamic_beta fit: its path, or N(0, g s2) when static
			am_ic <- NULL
			if (am_centred && isTRUE(beta_dyn$any)) {
				am_kd <- match(1L, beta_dyn$dynamic_idx)
				if (!is.na(am_kd)) {
					am_ic <- list(
						kind = "path",
						path = beta_dyn_path[, am_kd],
						rho = rho_beta_by_group[group_id[am_kd]],
						sd2 = sigma_beta_by_group[group_id[am_kd]]^2 * Lambda_full[am_kd, am_kd],
						m0 = beta0_mean[am_kd],
						v0 = beta0_cov[am_kd, am_kd]
					)
				} else {
					am_ks <- match(1L, beta_dyn$static_idx)
					am_ic <- list(
						kind = "static",
						value = beta_static[am_ks],
						prec = 1 / (.lame_static_g_scale(g, beta_dyn$static_idx)[am_ks] * s2)
					)
				}
			}
			if (symmetric) {
				# symmetric positions are V = U diag(L), the weights read back from the pivot coordinate;
				# a state not of that form leaves the positions to the host sampler
				am_Ld <- NULL
				if (actor_move_uv) {
					am_Ld <- vapply(seq_len(R), function(r) {
						pv <- which.max(abs(U[, r]))
						if (U[pv, r] == 0) 1 else V[pv, r] / U[pv, r]
					}, numeric(1))
					if (!all(is.finite(am_Ld)) ||
							max(abs(V - U %*% diag(am_Ld, R))) > 1e-8 * max(1, max(abs(V))))
						am_Ld <- NULL
				}
				am_out <- .lame_actor_move_visit_sym(
					Z,
					EZ,
					Y,
					.lame_rho_scalar(rho),
					beta,
					if (dynamic_ab) a_lev else a,
					if (dynamic_ab) a_mat else NULL,
					if (is.null(am_Ld)) NULL else U,
					am_Ld,
					Sab[1L, 1L],
					am_iV0,
					s > .5 * burn,
					am_cx,
					centred = am_centred,
					ic = am_ic,
					lev = am_lev[1L]
				)
				Z <- am_out$Z; EZ <- am_out$EZ; beta <- am_out$beta
				if (dynamic_ab) {
					a_lev <- am_out$la; b_lev <- a_lev
					a_mat <- am_out$a_mat; b_mat <- a_mat
					a <- a_lev; b <- b_lev
				} else {
					a <- am_out$la; b <- a
				}
				if (!is.null(am_Ld)) { U <- am_out$U; V <- U %*% diag(am_Ld, R) }
			} else {
				am_out <- .lame_actor_move_visit(
					Z,
					EZ,
					.lame_rho_scalar(rho),
					beta,
					if (dynamic_ab) a_lev else a,
					if (dynamic_ab) b_lev else b,
					if (dynamic_ab) a_mat else NULL,
					if (dynamic_ab) b_mat else NULL,
					if (actor_move_uv) U else NULL,
					if (actor_move_uv) V else NULL,
					Sab,
					am_iV0,
					s > .5 * burn,
					am_cx,
					centred = am_centred,
					uv_prior = if (s > .5 * burn && actor_move_uv) uv_prior_move,
					ic = am_ic,
					lev = am_lev
				)
				Z <- am_out$Z; EZ <- am_out$EZ; beta <- am_out$beta
				if (dynamic_ab) {
					a_lev <- am_out$la; b_lev <- am_out$lb
					a_mat <- am_out$a_mat; b_mat <- am_out$b_mat
					a <- a_lev; b <- b_lev
				} else {
					a <- am_out$la; b <- am_out$lb
				}
				if (actor_move_uv) { U <- am_out$U; V <- am_out$V }
			}
			# the intercept moved with the levels: write it back to the coefficient state
			if (!is.null(am_ic)) {
				if (identical(am_ic$kind, "path")) {
					beta_dyn_path[, am_kd] <- am_out$ic$path
					beta[1L] <- mean(am_out$ic$path)
				} else {
					beta_static[am_ks] <- am_out$ic$value
					beta[1L] <- am_out$ic$value
				}
			}
			actor_move_stats[c("visits", names(am_out$stats))] <-
				actor_move_stats[c("visits", names(am_out$stats))] + c(1, am_out$stats)
			rm(am_out)
			if (actor_move_rho) {
				am_r <- .lame_binary_dynamic_collapsed_rho_mh(.lame_rho_scalar(rho), EZ, am_rho_cx)
				actor_move_stats[c("rho_attempts", "rho_accepts")] <-
					actor_move_stats[c("rho_attempts", "rho_accepts")] + c(am_r$attempts, am_r$accepts)
				if (am_r$rho != rho) {
					rho <- am_r$rho
					Z <- .lame_binary_dynamic_collapsed_rho_z(EZ, Y, rho, am_rho_cx)
				}
				rm(am_r)
			}
		}
		# bipartite counterpart: exact slice draws for row and column actors, then the latent array is redrawn
		if (actor_move_bip && (s %% am_every_bip == 0L)) {
			# prior of the intercept of a dynamic_beta fit (as for the unipartite update)
			am_ic <- NULL
			if (am_comp_bip && isTRUE(beta_dyn$any)) {
				am_kd <- match(1L, beta_dyn$dynamic_idx)
				if (!is.na(am_kd)) {
					am_ic <- list(
						kind = "path",
						path = beta_dyn_path[, am_kd],
						rho = rho_beta_by_group[group_id[am_kd]],
						sd2 = sigma_beta_by_group[group_id[am_kd]]^2 * Lambda_full[am_kd, am_kd],
						m0 = beta0_mean[am_kd],
						v0 = beta0_cov[am_kd, am_kd]
					)
				} else {
					am_ks <- match(1L, beta_dyn$static_idx)
					am_ic <- list(
						kind = "static",
						value = beta_static[am_ks],
						prec = 1 / (.lame_static_g_scale(g, beta_dyn$static_idx)[am_ks] * s2)
					)
				}
			}
			am_out <- .lame_actor_move_visit_bip(
				EZ,
				beta,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				if (dynamic_ab) a_mat else NULL,
				if (dynamic_ab) b_mat else NULL,
				if (actor_move_bip_uv) matrix(if (length(dim(U)) == 3L) U[, , 1L] else U, nA, RA) else NULL,
				if (actor_move_bip_uv) matrix(if (length(dim(V)) == 3L) V[, , 1L] else V, nB, RB) else NULL,
				G,
				Sab[1L, 1L],
				Sab[2L, 2L],
				am_iV0_bip,
				am_comp_bip,
				am_cx_bip,
				isTRUE(rvar),
				isTRUE(cvar),
				if (am_dev_bip) list(rho = rho_ab, sigma = sigma_ab) else NULL,
				ic = am_ic
			)
			EZ <- am_out$EZ; beta <- am_out$beta
			# the intercept moved with the levels: write it back to the coefficient state
			if (!is.null(am_ic)) {
				if (identical(am_ic$kind, "path")) {
					beta_dyn_path[, am_kd] <- am_out$ic$path
					beta[1L] <- mean(am_out$ic$path)
				} else {
					beta_static[am_ks] <- am_out$ic$value
					beta[1L] <- am_out$ic$value
				}
			}
			if (dynamic_ab) {
				a_lev <- am_out$la; b_lev <- am_out$lb
				a_mat <- am_out$a_mat; b_mat <- am_out$b_mat
				a <- a_lev; b <- b_lev
			} else {
				a <- am_out$la; b <- am_out$lb
			}
			if (actor_move_bip_uv) {
				U <- array(am_out$U, dim = dim(U))
				V <- array(am_out$V, dim = dim(V))
			}
			actor_move_stats[c("visits", names(am_out$stats))] <-
				actor_move_stats[c("visits", names(am_out$stats))] + c(1, am_out$stats)
			rm(am_out)
			Z <- rZ_bin_bip_batch_cpp(Z, EZ, Y)
			E.nrm <- Z - EZ
		}
		# position paths of a binary fit with moving positions, with the latent array integrated out and redrawn
		# (persistence and innovation scale held)
		if (dpm_bin && (s %% dpm_every == 0L)) {
			dpm_o <- .lame_dpm_visit_binary(
				dpm_cx,
				Z,
				EZ,
				Y,
				if (bip) 0 else .lame_rho_scalar(rho),
				U_cube,
				if (symmetric) NULL else V_cube,
				if (bip) G else NULL,
				if (symmetric) L_dyn else NULL,
				rho_uv,
				sigma_uv,
				max_abs = dpm_uv_max,
				lambda_u = if (identical(dynamic_uv_kind, "t")) lambda_u else NULL,
				lambda_v = if (identical(dynamic_uv_kind, "t")) lambda_v else NULL
			)
			Z <- dpm_o$Z; EZ <- dpm_o$EZ
			U_cube <- dpm_o$U
			if (symmetric) {
				for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
			} else V_cube <- dpm_o$V
			U <- apply(U_cube, c(1, 2), mean)
			V <- apply(V_cube, c(1, 2), mean)
			E.nrm <- Z - EZ
			dpm_stats["dpm_visits"] <- dpm_stats["dpm_visits"] + 1
			dpm_stats[names(dpm_o$stats)] <- dpm_stats[names(dpm_o$stats)] + dpm_o$stats
			rm(dpm_o)
		}
		# ridge of a row or column covariate with the additive effects of its margin (R/covariate_ridge_move.R)
		if (cr_on) {
			cr_dyn <- isTRUE(dynamic_ab)
			cr_o <- .lame_cov_ridge_visit(
				Xlist,
				beta,
				cr_P0,
				if (cr_dyn) a_lev else a,
				if (cr_dyn) b_lev else b,
				if (cr_dyn) a_mat else NULL,
				if (cr_dyn) b_mat else NULL,
				Sab,
				rho_ab,
				sigma_ab,
				rvar,
				cvar,
				center = isTRUE(intercept),
				E = Z - EZ,
				M = cr_M,
				rho_t = vapply(seq_len(N), function(tt) if (isTRUE(dcor)) .lame_rho_at(rho, tt) else 0, numeric(1)),
				s2 = s2,
				skip = if (isTRUE(intercept)) 1L else integer(0)
			)
			beta <- cr_o$beta
			if (cr_dyn) {
				a_lev <- cr_o$a_lev; b_lev <- cr_o$b_lev; a_mat <- cr_o$a_mat; b_mat <- cr_o$b_mat
				a <- a_lev; b <- b_lev
			} else {
				a <- cr_o$a_lev; b <- cr_o$b_lev
			}
			for (tt in which(cr_o$shift != 0)) {
				EZ[, , tt] <- EZ[, , tt] + cr_o$shift[tt]
				free <- if (is.null(cr_M)) TRUE else !cr_M[, , tt]
				Zt <- Z[, , tt]
				Zt[free] <- Zt[free] + cr_o$shift[tt]
				Z[, , tt] <- Zt
			}
			cr_stats[c("visits", names(cr_o$stats))] <- cr_stats[c("visits", names(cr_o$stats))] + c(1, cr_o$stats)
			rm(cr_o)
		}
		# static positions with the latent rows and columns of the actors they bound (R/sparse_position_move.R)
		if (spm_on && (s %% spm_every == 0L)) {
			spm_o <- .lame_spm_visit(
				spm_ctx,
				Z,
				U,
				V,
				shrink = s > .5 * burn,
				uv_prior = if (s > .5 * burn) uv_prior_move,
				alpha = if (use_explicit_cutpoints) alpha else NULL
			)
			Z <- array(spm_o$Z, dim(Z)); EZ <- EZ + spm_o$shift; U <- spm_o$U; V <- spm_o$V
			spm_stats[c("visits", names(spm_o$stats))] <- spm_stats[c("visits", names(spm_o$stats))] + c(1, spm_o$stats)
			rm(spm_o)
		}
		# translations of the levels with their latent rows and columns (poisson, ordinal)
		if (ltm_on && (s %% ltm_every == 0L)) {
			# prior of the intercept of a dynamic_beta fit (as for the binary update)
			ltm_ic <- NULL
			if (am_centred && isTRUE(beta_dyn$any)) {
				am_kd <- match(1L, beta_dyn$dynamic_idx)
				if (!is.na(am_kd)) {
					ltm_ic <- list(
						kind = "path",
						path = beta_dyn_path[, am_kd],
						rho = rho_beta_by_group[group_id[am_kd]],
						sd2 = sigma_beta_by_group[group_id[am_kd]]^2 * Lambda_full[am_kd, am_kd],
						m0 = beta0_mean[am_kd],
						v0 = beta0_cov[am_kd, am_kd]
					)
				} else {
					am_ks <- match(1L, beta_dyn$static_idx)
					ltm_ic <- list(
						kind = "static",
						value = beta_static[am_ks],
						prec = 1 / (.lame_static_g_scale(g, beta_dyn$static_idx)[am_ks] * s2)
					)
				}
			}
			ltm_iV0_now <- if (bip) diag(length(beta)) / (as.numeric(g[1]) * s2) else ltm_iV0
			ltm_st <- .lame_ltm_visit(
				ltm_cx,
				Z,
				EZ,
				if (dynamic_ab) a_lev else a,
				if (dynamic_ab) b_lev else b,
				if (dynamic_ab) a_mat else NULL,
				if (dynamic_ab) b_mat else NULL,
				beta,
				Sab,
				ltm_iV0_now,
				centred = am_centred,
				ic = ltm_ic,
				lev = ltm_lev,
				alpha = if (use_explicit_cutpoints) alpha else NULL
			)
			ltm_stats[c("visits", names(ltm_st$stats))] <- ltm_stats[c("visits", names(ltm_st$stats))] + c(1, ltm_st$stats)
			if (ltm_loc) {
				ltm_st <- .lame_ltm_location(
					ltm_cx,
					ltm_st$Z,
					ltm_st$EZ,
					ltm_st$la,
					ltm_st$lb,
					ltm_st$a_mat,
					ltm_st$b_mat,
					ltm_st$beta,
					Sab,
					ltm_iV0_now,
					centred = am_centred,
					ic = ltm_st$ic,
					lev = ltm_lev,
					alpha = if (use_explicit_cutpoints) alpha else NULL
				)
				ltm_stats["location_calls"] <- ltm_stats["location_calls"] + 1
			}
			Z <- ltm_st$Z; EZ <- ltm_st$EZ; beta <- ltm_st$beta
			if (symmetric) { ltm_st$lb <- ltm_st$la; ltm_st$b_mat <- ltm_st$a_mat }
			if (dynamic_ab) {
				a_lev <- ltm_st$la; b_lev <- ltm_st$lb
				a_mat <- ltm_st$a_mat; b_mat <- ltm_st$b_mat
				a <- a_lev; b <- b_lev
			} else {
				a <- ltm_st$la; b <- ltm_st$lb
			}
			# the intercept moved with the levels: write it back to the coefficient state
			if (!is.null(ltm_ic)) {
				if (identical(ltm_ic$kind, "path")) {
					beta_dyn_path[, am_kd] <- ltm_st$ic$path
					beta[1L] <- mean(ltm_st$ic$path)
				} else {
					beta_static[am_ks] <- ltm_st$ic$value
					beta[1L] <- ltm_st$ic$value
				}
			}
			# poisson: the coefficients move with the latent array, then the dyadic variance and correlation
			if (ltm_coef) {
				ltm_cf <- .lame_ppm_static(ppm_cx, Z, EZ, beta, ltm_Xk, ltm_iV0_now, Xz = ltm_Xz)
				Z <- ltm_cf$Z; EZ <- ltm_cf$EZ; beta <- ltm_cf$beta
				ltm_stats[names(ltm_cf$stats)] <- ltm_stats[names(ltm_cf$stats)] + ltm_cf$stats
				rm(ltm_cf)
			}
			if (ltm_uv) {
				ltm_uo <- NULL
				if (bip) {
					ltm_uo <- .lame_ltm_uv_visit(
						ltm_cx,
						Z,
						EZ,
						matrix(if (length(dim(U)) == 3L) U[, , 1L] else U, nA, RA),
						matrix(if (length(dim(V)) == 3L) V[, , 1L] else V, nB, RB),
						G
					)
					U <- array(ltm_uo$U, dim = dim(U)); V <- array(ltm_uo$V, dim = dim(V))
				} else if (symmetric) {
					# symmetric positions are V = U diag(L), the weights read back from the pivot coordinate;
					# a state not of that form leaves the positions to the host sampler
					ltm_Ld <- vapply(seq_len(R), function(r) {
						pv <- which.max(abs(U[, r]))
						if (U[pv, r] == 0) 1 else V[pv, r] / U[pv, r]
					}, numeric(1))
					if (all(is.finite(ltm_Ld)) && max(abs(V - U %*% diag(ltm_Ld, R))) <= 1e-8 * max(1, max(abs(V)))) {
						ltm_uo <- .lame_ltm_uv_visit(ltm_cx, Z, EZ, U, NULL, NULL, ltm_Ld, shrink = s > .5 * burn)
						U <- ltm_uo$U; V <- U %*% diag(ltm_Ld, R)
					}
				} else {
					ltm_uo <- .lame_ltm_uv_visit(
						ltm_cx,
						Z,
						EZ,
						U,
						V,
						shrink = s > .5 * burn,
						uv_prior = if (s > .5 * burn) uv_prior_move
					)
					U <- ltm_uo$U; V <- ltm_uo$V
				}
				if (!is.null(ltm_uo)) {
					Z <- ltm_uo$Z; EZ <- ltm_uo$EZ
					ltm_stats[names(ltm_uo$stats)] <- ltm_stats[names(ltm_uo$stats)] + ltm_uo$stats
				}
				rm(ltm_uo)
			}
			if (ltm_duv) {
				# the position paths move with the latent array, persistence and innovation scale held
				ltm_dp <- .lame_dpm_visit(
					ltm_cx,
					Z,
					EZ,
					U_cube,
					if (symmetric) NULL else V_cube,
					if (bip) G else NULL,
					if (symmetric) L_dyn else NULL,
					rho_uv,
					sigma_uv,
					max_abs = ltm_uv_max,
					lambda_u = if (identical(dynamic_uv_kind, "t")) lambda_u else NULL,
					lambda_v = if (identical(dynamic_uv_kind, "t")) lambda_v else NULL
				)
				Z <- ltm_dp$Z; EZ <- ltm_dp$EZ
				U_cube <- ltm_dp$U
				if (symmetric) {
					for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
				} else V_cube <- ltm_dp$V
				U <- apply(U_cube, c(1, 2), mean)
				V <- apply(V_cube, c(1, 2), mean)
				ltm_stats[names(ltm_dp$stats)] <- ltm_stats[names(ltm_dp$stats)] + ltm_dp$stats
				rm(ltm_dp)
			}
			if (ltm_scale) {
				ltm_s2b <- if (beta_dyn$any) .lame_s2_beta_terms(
					beta_static,
					.lame_static_g_scale(g, beta_dyn$static_idx)
				)
				else if (bip) .lame_s2_beta_terms(beta, g[1L]) else c(df = 0, ss = 0)
				ltm_two <- !bip && !symmetric && isTRUE(dcor)
				ltm_sc <- .lame_ltm_scale(
					ltm_cx,
					Z,
					EZ,
					s2,
					if (ltm_two) .lame_rho_scalar(rho) else 0,
					c(shape0 = (1 + ltm_s2b[["df"]]) / 2, rate0 = (s2_prior_scale + ltm_s2b[["ss"]]) / 2),
					two = ltm_two,
					rho_logprior = .lame_rho_prior_log
				)
				Z <- ltm_sc$Z; s2 <- ltm_sc$s2
				if (ltm_two) rho <- ltm_sc$rho
				ltm_stats[names(ltm_sc$stats)] <- ltm_stats[names(ltm_sc$stats)] + ltm_sc$stats
				rm(ltm_sc)
			}
			if (identical(family, "poisson")) E.nrm <- Z - EZ
			rm(ltm_st)
		}
		# actor-period effects of a poisson dynamic_ab fit move with the latent array
		if (pab_on) {
			pab_o <- .lame_pab_visit(
				pab_cx,
				Z,
				EZ,
				a_mat,
				b_mat,
				a_lev,
				b_lev,
				rho_ab,
				sigma_ab,
				rv = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
				cv = isTRUE(cvar),
				centred = isTRUE(intercept) && length(beta) > 0L
			)
			Z <- pab_o$Z; EZ <- pab_o$EZ; a_mat <- pab_o$a_mat; b_mat <- pab_o$b_mat
			pab_stats <- pab_stats + c(1, pab_o$stats)
			rm(pab_o)
		}
		# translations of the dynamic_beta coefficients with their latent slices (poisson), then the static
		# coefficients (R/poisson_period_move.R)
		if (ppm_on) {
			ppm_out <- .lame_ppm_visit(
				ppm_cx,
				Z,
				EZ,
				beta_dyn_path,
				Xdyn_list_per_t,
				list(
					kind = dynamic_beta_kind,
					rho = rho_beta_by_group[group_id],
					sigma = sigma_beta_by_group[group_id],
					lam = diag(Lambda_full),
					m0 = beta0_mean,
					v0 = diag(beta0_cov),
					time_positions = ppm_tpos,
					length_scale = prior$matern32_length_scale %||% NULL
				),
				Xz = ppm_Xz
			)
			Z <- ppm_out$Z; EZ <- ppm_out$EZ; beta_dyn_path <- ppm_out$path
			beta[beta_dyn$dynamic_idx] <- colMeans(beta_dyn_path)
			ppm_stats[c(
				"visits",
				names(ppm_out$stats)
			)] <- ppm_stats[c("visits", names(ppm_out$stats))] + c(1, ppm_out$stats)
			if (ppm_static && s %% ppm_static_every == 0L) {
				ppm_sc <- .lame_ppm_static(
					ppm_cx,
					Z,
					EZ,
					beta_static,
					ppm_Xs,
					diag(1 / (.lame_static_g_scale(g, beta_dyn$static_idx) * s2), length(beta_static)),
					Xz = ppm_Xzs
				)
				Z <- ppm_sc$Z; EZ <- ppm_sc$EZ; beta_static <- ppm_sc$beta
				beta[beta_dyn$static_idx] <- beta_static
				ppm_stats[c(
					"static_calls", "static_evals", "static_stuck"
				)] <- ppm_stats[c("static_calls", "static_evals", "static_stuck")] +
					ppm_sc$stats[c("coef_calls", "coef_evals", "coef_stuck")]
				rm(ppm_sc)
			}
			rm(ppm_out)
		}
		# dyadic variance and correlation with the latent residuals carried along their conditional law
		if (ppm_scale_on) {
			ppm_s2b <- if (beta_dyn$any) .lame_s2_beta_terms(beta_static, .lame_static_g_scale(g, beta_dyn$static_idx))
			else if (bip) .lame_s2_beta_terms(beta, g[1L]) else c(df = 0, ss = 0)
			ppm_two <- !bip && !symmetric && isTRUE(dcor)
			ppm_sc <- .lame_ppm_scale(
				ppm_cx,
				Z,
				EZ,
				s2,
				if (ppm_two) .lame_rho_scalar(rho) else 0,
				c(shape0 = (1 + ppm_s2b[["df"]]) / 2, rate0 = (s2_prior_scale + ppm_s2b[["ss"]]) / 2),
				ppm_two,
				.lame_rho_prior_log
			)
			Z <- ppm_sc$Z; s2 <- ppm_sc$s2
			if (ppm_two) rho <- ppm_sc$rho
			# the residuals the sweep's own s2 and rho draws read
			E.nrm <- Z - EZ
			ppm_stats[names(ppm_sc$stats)] <- ppm_stats[names(ppm_sc$stats)] + ppm_sc$stats
			# this scale move stands in for that of the actor translations, so it is recorded under their stats
			if (ltm_on) {
				ltm_stats[c("scale_attempts", "scale_accepts", "rho_attempts", "rho_accepts")] <-
					ltm_stats[c("scale_attempts", "scale_accepts", "rho_attempts", "rho_accepts")] +
					ppm_sc$stats[c("scale_attempts", "scale_accepts", "scale_rho_attempts", "scale_rho_accepts")]
			}
			if (!ppm_on) ppm_stats["visits"] <- ppm_stats["visits"] + 1
			rm(ppm_sc)
		}
		# joint update with the probit auxiliary responses integrated out; fixed hierarchy, exact actor contrasts
		if (marginal_move && s > .5 * burn && s %% 20L == 0L) {
			cx <- marginal_context
			hmc_state <- if (static_actor_marginal) list(
				static = beta_static,
				path = beta_dyn_path,
				level = crossprod(cx$H, cbind(a, b)),
				U = U,
				V = V
			)
			else list(
				static = beta_static,
				path = beta_dyn_path,
				level = crossprod(cx$H, cbind(a_lev, b_lev)),
				adev = t(crossprod(cx$H, a_mat - a_lev)),
				bdev = t(crossprod(cx$H, b_mat - b_lev)),
				U = U,
				V = V
			)
			hmc_h <- list(
				g = .lame_static_g_scale(g, beta_dyn$static_idx) * s2,
				rho_beta = rho_beta_by_group[group_id],
				var_beta = sigma_beta_by_group[group_id]^2 * diag(Lambda_full),
				iC0 = solve(beta0_cov),
				iSab = solve(Sab),
				rho_ab = rho_ab,
				var_ab = sigma_ab^2
			)
			hmc_step <- if (static_actor_marginal) .lame_binary_static_hmc_step else .lame_binary_hmc_step
			hmc_unpack <- if (static_actor_marginal) .lame_binary_static_hmc_unpack else .lame_binary_hmc_unpack
			hmc_evaluate <- if (static_actor_marginal) .lame_binary_static_hmc_evaluate else .lame_binary_hmc_evaluate
			hmc_out <- hmc_step(.lame_binary_hmc_pack(hmc_state), cx, hmc_h, marginal_sd)
			marginal_att <- marginal_att + 1L
			marginal_acc <- marginal_acc + as.integer(hmc_out$accepted)
			marginal_bad <- marginal_bad + as.integer(hmc_out$nonfinite)
			if (s <= burn) marginal_sd <- min(
				.5,
				max(1e-5, exp(log(marginal_sd) + (hmc_out$prob - .65) / sqrt(marginal_att)))
			)
			if (hmc_out$accepted) {
				z_hmc <- hmc_unpack(hmc_out$q, cx)
				beta_static <- z_hmc$static; beta_dyn_path <- z_hmc$path
				beta[beta_dyn$static_idx] <- beta_static; beta[beta_dyn$dynamic_idx] <- colMeans(beta_dyn_path)
				lv <- cx$H %*% z_hmc$level; a_lev <- lv[, 1L]; b_lev <- lv[, 2L]
				a <- a_lev; b <- b_lev
				if (!static_actor_marginal) {
					a_mat <- a_lev + cx$H %*% t(z_hmc$adev); b_mat <- b_lev + cx$H %*% t(z_hmc$bdev)
				}
				U <- z_hmc$U; V <- z_hmc$V
				EZ <- hmc_evaluate(hmc_out$q, cx, hmc_h)$eta
				Z <- .lame_binary_z_given_factor(
					Z,
					EZ,
					array(0, dim(Z)),
					0,
					lev_i1,
					lev_i0,
					lev_ina,
					lev_idg
				)
			}
		}
		# exact translations: move the factor means with compensating centred actor deviations and intercept,
		# preserving the full predictor
		if (dynamic_joint_move && s > .5 * burn &&
				s %% 10L == 0L &&
				is.null(time_index)) {
			tr <- .lame_binary_dynamic_translation_collapsed_step(
				U_cube,
				V_cube,
				b_mat,
				b_lev,
				a_lev,
				beta,
				lev_iV0,
				rho_uv,
				sigma_uv,
				rho_ab,
				sigma_ab,
				Sab,
				2L,
				lev_uv_max
			)
			U_cube <- tr$U; b_mat <- tr$margin; b_lev <- tr$level; beta <- tr$beta
			tr <- .lame_binary_dynamic_translation_collapsed_step(
				V_cube,
				U_cube,
				a_mat,
				a_lev,
				b_lev,
				beta,
				lev_iV0,
				rho_uv,
				sigma_uv,
				rho_ab,
				sigma_ab,
				Sab,
				1L,
				lev_uv_max
			)
			V_cube <- tr$U; a_mat <- tr$margin; a_lev <- tr$level; beta <- tr$beta
			a <- a_lev; b <- b_lev
			U <- apply(U_cube, c(1, 2), mean)
			V <- apply(V_cube, c(1, 2), mean)
			translation_att <- translation_att + 2L
		}
		# non-centred factor scale with binary Z integrated out, given the shared reciprocity factor
		if (dynamic_joint_move && s > .5 * burn && s %% 30L == 0L) {
			rho_us <- .lame_rho_scalar(rho)
			off_us <- .lame_dyad_factor_draw(Z, EZ, rho_us, lev_upper)
			sd_us <- sqrt(1 - abs(rho_us))
			uv_us <- array(0, dim(EZ))
			for (t_us in seq_len(N)) uv_us[, , t_us] <-
				tcrossprod(matrix(U_cube[, , t_us], n, R), matrix(V_cube[, , t_us], n, R))
			q_us <- .lame_binary_dynamic_uv_scale_step(
				EZ[lev_obs] / sd_us,
				uv_us[lev_obs] / sd_us,
				off_us[lev_obs] / sd_us,
				lev_sY,
				sigma_uv,
				prior$sigma_uv_shape,
				prior$sigma_uv_scale,
				variance_max = min(
					(prior$sigma_uv_max %||% Inf)^2,
					(prior$uv_var_cap %||% Inf) * (1 - rho_uv^2)
				),
				factor_max = lev_uv_max / max(abs(U_cube), abs(V_cube))
			)
			uv_scale_att <- uv_scale_att + q_us$attempts
			uv_scale_acc <- uv_scale_acc + q_us$accepts
			if (q_us$factor != 1) {
				U_cube <- U_cube * q_us$factor; V_cube <- V_cube * q_us$factor
				U <- U * q_us$factor; V <- V * q_us$factor
				sigma_uv <- q_us$sigma
				EZ <- EZ + (q_us$factor^2 - 1) * uv_us
				Z <- .lame_binary_z_given_factor(
					Z,
					EZ,
					off_us,
					rho_us,
					lev_i1,
					lev_i0,
					lev_ina,
					lev_idg
				)
			}
		}
		# reciprocity update under the observed binary-pair likelihood with both auxiliaries integrated out;
		# Z is refreshed exactly
		if (dynamic_joint_move && isTRUE(dcor) && s > .5 * burn && s %% 30L == 0L) {
			mr <- .lame_binary_dynamic_collapsed_rho_mh(rho, EZ, rho_collapsed_context)
			rho_collapsed_att <- rho_collapsed_att + mr$attempts
			rho_collapsed_acc <- rho_collapsed_acc + mr$accepts
			if (mr$rho != rho) {
				rho <- mr$rho
				Z <- .lame_binary_dynamic_collapsed_rho_z(EZ, Y, rho, rho_collapsed_context)
			}
		}
		# persistence / innovation-variance interweaving at fixed stationary variance and standard innovations;
		# the density includes both Jacobians, and bounds reject rather than project
		if (ncar_enabled && s > .5 * burn && s %% 20L == 0L) {
			nc_raw <- .lame_binary_dynamic_rho_factor_draw(Z, EZ, rho)
			nc_offset <- .lame_binary_dynamic_rho_offset(nc_raw, rho)
			nc_sd <- sqrt(1 - abs(rho))
			nc_uv <- array(0, dim(EZ))
			for (nt in seq_len(N)) nc_uv[, , nt] <- tcrossprod(
				matrix(U_cube[, , nt], n, R),
				matrix(V_cube[, , nt], n, R)
			)
			nc_base <- EZ - nc_uv
			nc_like_uv <- function(paths) {
				pu <- array(paths[seq_len(n * R), , drop = FALSE], c(n, R, N))
				pv <- array(paths[n * R + seq_len(n * R), , drop = FALSE], c(n, R, N))
				eta <- nc_base
				for (tt in seq_len(N)) eta[, , tt] <- eta[, , tt] +
					tcrossprod(matrix(pu[, , tt], n, R), matrix(pv[, , tt], n, R))
				list(
					loglik = sum(pnorm(
						lev_sY * (eta[lev_obs] + nc_offset[lev_obs]) / nc_sd,
						log.p = TRUE
					)),
					eta = eta
				)
			}
			ncu <- .lame_binary_dynamic_nc_ar_step(
				rbind(matrix(U_cube, n * R, N), matrix(V_cube, n * R, N)),
				rho_uv,
				sigma_uv,
				nc_like_uv,
				prior$rho_uv_mean,
				prior$rho_uv_sd,
				prior$sigma_uv_shape,
				prior$sigma_uv_scale,
				stationary_max = prior$uv_var_cap %||% Inf,
				sigma_max = prior$sigma_uv_max %||% Inf,
				max_abs = lev_uv_max
			)
			ncar_uv_att <- ncar_uv_att + ncu$attempts; ncar_uv_acc <- ncar_uv_acc + ncu$accepts
			if (ncu$accepts > 0L) {
				U_cube <- array(ncu$paths[seq_len(n * R), , drop = FALSE], c(n, R, N))
				V_cube <- array(ncu$paths[n * R + seq_len(n * R), , drop = FALSE], c(n, R, N))
				U <- apply(U_cube, c(1, 2), mean);V <- apply(V_cube, c(1, 2), mean)
				rho_uv <- ncu$rho;sigma_uv <- ncu$sigma;EZ <- ncu$evaluation$eta
			}
			nc_ad <- a_mat - a_lev;nc_bd <- b_mat - b_lev
			nc_base <- EZ
			for (nt in seq_len(N)) nc_base[, , nt] <- nc_base[, , nt] - outer(nc_ad[, nt], nc_bd[, nt], '+')
			nc_like_ab <- function(paths) {
				eta <- nc_base
				for (tt in seq_len(N)) eta[, , tt] <- eta[, , tt] + outer(
					paths[seq_len(n), tt],
					paths[n + seq_len(n), tt],
					'+'
				)
				list(
					loglik = sum(pnorm(
						lev_sY * (eta[lev_obs] + nc_offset[lev_obs]) / nc_sd,
						log.p = TRUE
					)),
					eta = eta
				)
			}
			nca <- .lame_binary_dynamic_nc_ar_step(
				rbind(nc_ad, nc_bd),
				rho_ab,
				sigma_ab,
				nc_like_ab,
				prior$rho_ab_mean,
				prior$rho_ab_sd,
				prior$sigma_ab_shape,
				prior$sigma_ab_scale,
				stationary_max = (prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
			)
			ncar_ab_att <- ncar_ab_att + nca$attempts;ncar_ab_acc <- ncar_ab_acc + nca$accepts
			if (nca$accepts > 0L) {
				a_mat <- a_lev + nca$paths[seq_len(n), , drop = FALSE]
				b_mat <- b_lev + nca$paths[n + seq_len(n), , drop = FALSE]
				rho_ab <- nca$rho;sigma_ab <- nca$sigma;EZ <- nca$evaluation$eta
			}
			if (ncu$accepts + nca$accepts > 0L) Z <- .lame_binary_dynamic_rho_z_draw(
				Z,
				EZ,
				nc_offset,
				rho,
				lev_i1,
				lev_i0,
				lev_ina,
				lev_idg
			)
		}
		# skipped when rho_uv sits at the AR endpoint, which cannot be Fisher transformed
		if (dhmc_enabled && s > .5 * burn && s %% 20L == 0L && abs(rho_uv) < .99) {
			dh_offset <- NULL # both Gaussian auxiliaries integrated out
			dh_hyper <- .lame_binary_dynamic_ar_hmc_hyper(
				dhmc_context,
				lev_iV0,
				Sab,
				rho_ab,
				sigma_ab,
				rho,
				prior$sigma_uv_shape,
				prior$sigma_uv_scale,
				prior$rho_uv_mean,
				prior$rho_uv_sd,
				variance_max = (prior$sigma_uv_max %||% Inf)^2,
				stationary_max = prior$uv_var_cap %||% Inf,
				max_abs = lev_uv_max
			)
			dh_q <- .lame_binary_dynamic_ar_hmc_encode(
				beta,
				cbind(a_lev, b_lev),
				a_mat - a_lev,
				b_mat - b_lev,
				U_cube,
				V_cube,
				sigma_uv,
				rho_uv,
				dhmc_context,
				dh_hyper
			)
			dh_proposal <- .lame_binary_dynamic_ar_hmc_step(dh_q, dhmc_context, dh_hyper, dhmc_sd)
			dhmc_att <- dhmc_att + 1L
			dhmc_nonfinite <- dhmc_nonfinite + as.integer(dh_proposal$nonfinite)
			dhmc_bounds <- dhmc_bounds + as.integer(dh_proposal$bound)
			if (dh_proposal$accepted) {
				dh_state <- .lame_binary_dynamic_ar_hmc_decode(dh_proposal$q, dhmc_context, dh_hyper)
				sigma_uv <- dh_state$sigma; rho_uv <- dh_state$rho_uv
				beta <- dh_state$beta
				a_lev <- dh_state$level[, 1]; b_lev <- dh_state$level[, 2]
				a_mat <- a_lev + dh_state$adev; b_mat <- b_lev + dh_state$bdev
				a <- a_lev; b <- b_lev
				U_cube <- dh_state$U; V_cube <- dh_state$V
				U <- apply(U_cube, c(1, 2), mean); V <- apply(V_cube, c(1, 2), mean)
				EZ <- dh_proposal$eta
				Z <- .lame_binary_dynamic_collapsed_rho_z(EZ, Y, rho, dhmc_context$pairs)
				dhmc_acc <- dhmc_acc + 1L
			}
			if (s <= burn) dhmc_sd <- min(
				.5,
				max(
					.001,
					exp(log(dhmc_sd) +
							(dh_proposal$prob - .65) / sqrt(dhmc_att))
				)
			)
		}
		# paired actor and factor paths under the observed binary likelihood, keeping centred pair sums; Z is refreshed
		if (actor_slice_enabled && s > .5 * burn && s %% 20L == 0L) {
			as_state <- tryCatch({
				z <- .lame_binary_actor_slice_round(
					U_cube,
					V_cube,
					a_mat,
					b_mat,
					a_lev,
					b_lev,
					EZ,
					Y,
					Sab,
					rho_ab,
					sigma_ab,
					rho_uv,
					sigma_uv,
					rho,
					lev_uv_max
				)
				z$Z <- .lame_binary_dynamic_collapsed_rho_z(z$eta, Y, rho, dhmc_context$pairs)
				z
			}, error = function(e) NULL)
			actor_slice_rounds <- actor_slice_rounds + 1L
			if (is.null(as_state)) actor_slice_bad <- actor_slice_bad + 1L else {
				U_cube <- as_state$U; V_cube <- as_state$V
				a_mat <- as_state$A; b_mat <- as_state$B
				a_lev <- as_state$alevel; b_lev <- as_state$blevel
				a <- a_lev; b <- b_lev
				U <- apply(U_cube, c(1, 2), mean); V <- apply(V_cube, c(1, 2), mean)
				EZ <- as_state$eta; Z <- as_state$Z
				actor_slice_blocks <- actor_slice_blocks + as_state$blocks
				actor_slice_evaluations <- actor_slice_evaluations + as_state$evaluations
				actor_slice_bounds <- actor_slice_bounds + as_state$bound_rejections
			}
		}
		# explicit cutpoints after the z sweep: correlated directed models draw alpha given Z, independent and
		# symmetric ones use the Z-marginalised Cowles update
		if (use_explicit_cutpoints) {
			alpha_step <- if (!symmetric && dcor) .sample_alpha_given_z(alpha, Y_int, Z)
			else sample_alpha_cowles(
				alpha = alpha,
				Y_int = Y_int,
				EZ = EZ,
				tau_prop = tau_prop,
				symmetric = symmetric
			)
			alpha <- alpha_step$alpha
			alpha_accept_total <- alpha_accept_total + isTRUE(alpha_step$accept)
			alpha_accept_n <- alpha_accept_n + 1L
			# Robbins-Monro drift of log(tau) toward the target acceptance, in burn-in only
			if ((symmetric || !dcor) && s <= burn && alpha_accept_n >= 25L) {
				acc_rate <- alpha_accept_total / alpha_accept_n
				log_tau <- log(tau_prop) + (acc_rate - tau_target) / sqrt(s)
				tau_prop <- min(max(exp(log_tau), 0.01), 2.0)
				alpha_accept_n <- 0L
				alpha_accept_total <- 0L
			}
			# redraw Z given the updated alpha before the Gaussian parameter moves
			for (t in seq_len(N)) {
				Z[, , t] <- if (symmetric)
					rZ_ord_sym_explicit_fc(Z[, , t], EZ[, , t], Y_int[, , t], alpha)
				else rZ_ord_explicit_fc(
					Z[, , t],
					EZ[, , t],
					Y_int[, , t],
					alpha,
					rho = if (dcor) .lame_rho_scalar(rho) else 0
				)
			}
		}
		# the other blocks work on Z minus the actor-slope offset; subtracting it from both Z and EZ keeps
		# the residuals for s2 and rho
		if (per_actor_active) {
			Z <- Z - per_actor_offset
			EZ <- EZ - per_actor_offset
		}
		# update s2
		s2_before_cap <- s2
		if (is.element(family, c("normal", "poisson"))) {
			# include every N(0, g_j * s2) coefficient prior; dynamic paths have their own innovation variance
			s2_beta <- if (beta_dyn$any) {
				.lame_s2_beta_terms(
					beta_static,
					.lame_static_g_scale(g, beta_dyn$static_idx)
				)
			} else if (bip) {
				.lame_s2_beta_terms(beta, g[1L])
			} else c(df = 0, ss = 0)
			if (bip) {
				# bipartite: direct s2 draw (no dyadic correlation; rs2_rep_fc_cpp assumes square matrices)
				resid_all <- Z - EZ
				n_obs <- sum(!is.na(resid_all))
				ss <- sum(resid_all^2, na.rm = TRUE)
				s2 <- 1 / rgamma(
					1,
					(1 + n_obs + s2_beta["df"]) / 2,
					(s2_prior_scale + ss + s2_beta["ss"]) / 2
				)
			} else if (symmetric) {
				# symmetric: draw s2 from the upper-triangle residuals pooled across periods, each dyad once, since the
				# whitened update at the ramped rho near one would halve it
				e <- unlist(lapply(seq_len(N), function(t) {
					r <- E.nrm[, , t]; r[upper.tri(r)]
				}))
				e <- e[is.finite(e)]
				n_e <- length(e)
				if (n_e > 1) {
					ssr <- sum(e^2)
					s2 <- 1 / rgamma(
						1,
						shape = (n_e + 1 + s2_beta["df"]) / 2,
						rate = (ssr + s2_prior_scale + s2_beta["ss"]) / 2
					)
				}
			} else {
				s2New <- try({
					if (isTRUE(dynamic_rho)) {
						.lame_sample_s2_dynamic_rho(
							E.nrm,
							rho,
							current = s2,
							s20 = s2_prior_scale,
							beta_df = s2_beta["df"],
							beta_ss = s2_beta["ss"]
						)
					} else {
						rs2_rep_fc_cpp(
							E.nrm,
							solve(matrix(c(1, rho, rho, 1), 2, 2)),
							s20 = s2_prior_scale,
							beta_df = s2_beta["df"],
							beta_ss = s2_beta["ss"]
						)
					}
				}, silent = TRUE)
				if (!inherits(s2New, 'try-error')) { s2 <- s2New } else { tryErrorChecks$s2 <- tryErrorChecks$s2 + 1 }
			}
		}

		if (identical(family, "normal")) {
			s2 <- .lame_ar_s2_support(
				s2,
				s2_before_cap,
				rho_ab = if (dynamic_ab) rho_ab else NULL,
				sigma_ab = if (dynamic_ab) sigma_ab else NULL,
				ab_cap = prior$ab_var_cap %||% Inf,
				rho_uv = if (.dyn_uv_has_rank) rho_uv else NULL,
				sigma_uv = if (.dyn_uv_has_rank) sigma_uv else NULL,
				uv_cap = prior$uv_var_cap %||% Inf
			)
		}

		# the unipartite variance draw excludes the diagonal, so redraw it under the new variance
		if (!bip && family %in% c("normal", "poisson")) {
			Z <- .lame_refresh_diagonal(Z, .lame_diagonal(EZ), rho, s2)
		}

		# update beta, a and b with the g-prior
		if (beta_dyn$any) {
			# dynamic_beta: FFBS draw of beta_dyn_path given a, b, U, V, a conjugate draw of beta_static, then
			# a and b given the per-period beta. first the per-period offset for the FFBS
			n_a_ <- if (bip) nA else n
			n_b_ <- if (bip) nB else n
			offset_list_ffbs <- vector("list", N)
			Z_list_ffbs <- vector("list", N)
			# per-period a, b, U, V from the current state
			a_mat_eff <- if (dynamic_ab) a_mat else matrix(rep(a, N), n_a_, N)
			b_mat_eff <- if (dynamic_ab) b_mat else matrix(rep(b, N), n_b_, N)
			# dynamic UV uses its cubes; bipartite stores U/V as 3-d cubes even when static; only the unipartite
			# static case is replicated across periods
			U_cube_eff <- if (.dyn_uv_has_rank) U_cube
			else if (length(dim(U)) == 3L) U
			else {
				if (R > 0 || (bip && RA > 0)) {
					Uc <- array(0, dim = c(nrow(U), max(1L, ncol(U)), N))
					for (t in seq_len(N)) Uc[, , t] <- U
					Uc
				} else array(0, dim = c(n_a_, 1L, N))
			}
			V_cube_eff <- if (.dyn_uv_has_rank) V_cube
			else if (length(dim(V)) == 3L) V
			else {
				if (R > 0 || (bip && RB > 0)) {
					Vc <- array(0, dim = c(nrow(V), max(1L, ncol(V)), N))
					for (t in seq_len(N)) Vc[, , t] <- V
					Vc
				} else array(0, dim = c(n_b_, 1L, N))
			}
			G_eff <- if (bip) G else diag(1, max(1L, dim(U_cube_eff)[2]))
			# offset = a + b + UV (no X beta)
			for (t in seq_len(N)) {
				# bipartite ranks live in RA/RB, and positions are taken per period from the cubes
				if (bip) {
					off_t <- outer(a_mat_eff[, t], b_mat_eff[, t], "+")
					if (RA > 0 && RB > 0 && nrow(G_eff) > 0 && ncol(G_eff) > 0) {
						Ut <- U_cube_eff[, , t]; Vt <- V_cube_eff[, , t]
						if (!is.matrix(Ut)) Ut <- matrix(Ut, n_a_, 1)
						if (!is.matrix(Vt)) Vt <- matrix(Vt, n_b_, 1)
						off_t <- off_t + Ut %*% G_eff %*% t(Vt)
					}
				} else {
					off_t <- outer(a_mat_eff[, t], b_mat_eff[, t], "+")
					if (R > 0) {
						Ut <- U_cube_eff[, , t]; Vt <- V_cube_eff[, , t]
						if (!is.matrix(Ut)) Ut <- matrix(Ut, n_a_, 1)
						if (!is.matrix(Vt)) Vt <- matrix(Vt, n_b_, 1)
						off_t <- off_t + Ut %*% t(Vt)
					}
				}
				offset_list_ffbs[[t]] <- off_t
				Z_list_ffbs[[t]] <- Z[, , t]
			}
			use_dyad_rho <- (!bip) && (!symmetric) && isTRUE(dcor) && (.lame_rho_scalar(rho) != 0)
			# interweaving step: the hyperparameters and path move together with the standardised innovations held
			# fixed, before the path is redrawn; it hands on the per-period statistics the filter needs
			nc_st <- NULL
			if (beta_hyper_nc && N >= 2L && length(beta_dyn$dynamic_idx) > 0L) {
				nc_st <- tryCatch(
					beta_dynamic_stats_cpp(
						Xdyn_list = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_static = if (length(beta_static) > 0L) beta_static else numeric(0),
						s2 = s2,
						dyad_rho = if (is.null(rho)) 0 else rho,
						bipartite = bip,
						symmetric = symmetric,
						use_dyad_rho = use_dyad_rho
					),
					error = function(e) NULL
				)
				if (!is.null(nc_st)) {
					nc_pool_scale <- prior$sigma_beta_scale
					if (dynamic_beta_pool %in% c("sigma", "both") && beta_dyn$n_groups >= 2L) {
						nc_pool_scale <- .pool_update_sigma_scale(
							sigma_beta_by_group,
							alpha = prior$sigma_beta_shape,
							hyper_a = prior$sigma_beta_shape,
							hyper_b = prior$sigma_beta_scale
						)
					}
					nc_gaps <- if (!is.null(time_index)) diff(as.numeric(time_index)) else rep(1, N - 1L)
					nc_out <- .lame_beta_nc_step(
						beta_dyn_path,
						nc_st,
						group_id,
						beta_dyn$n_groups,
						dynamic_beta_kind,
						diag(Lambda_full),
						sigma_beta_by_group,
						rho_beta_by_group,
						prior$sigma_beta_shape,
						nc_pool_scale,
						prior$rho_beta_mean,
						prior$rho_beta_sd,
						prior$rho_beta_lower,
						prior$rho_beta_upper,
						rho_pool = if (dynamic_beta_pool %in% c("rho", "both") && identical(dynamic_beta_kind, "ar1") &&
								beta_dyn$n_groups >= 2L) list(mu = pool_rho_mu, tau = pool_rho_tau) else NULL,
						movable_rho = identical(dynamic_beta_kind, "ar1"),
						gaps = nc_gaps
					)
					beta_dyn_path <- nc_out$path
					sigma_beta_by_group[] <- nc_out$sigma
					rho_beta_by_group[] <- nc_out$rho
					beta_nc_stats <- beta_nc_stats + c(1, nc_out$attempts, nc_out$accepts)
				}
			}
			# per-coefficient rho and sigma from the per-block values
			rho_by_coef <- rho_beta_by_group[group_id]
			sigma_by_coef <- sigma_beta_by_group[group_id]
			# draw the dynamic coefficient block
			.nonunit_gaps <- !is.null(time_index) &&
				any(abs(diff(as.numeric(time_index)) - 1) > 1e-10)
			ffbs_out <- if (dynamic_beta_kind %in% c("rw2", "matern32") || .nonunit_gaps) {
				tryCatch(
					.sample_beta_dynamic_alt(
						Xdyn_list = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_static = if (length(beta_static) > 0L) beta_static else numeric(0),
						sigma_by_coef = sigma_by_coef,
						Lambda = Lambda_full,
						s2 = s2,
						kind = dynamic_beta_kind,
						length_scale = prior$matern32_length_scale %||% NULL,
						time_positions = if (!is.null(time_index)) as.numeric(time_index) else NULL,
						rho_by_coef = rho_by_coef,
						beta0_mean = beta0_mean,
						beta0_cov = beta0_cov,
						dyad_rho = if (is.null(rho)) 0 else rho,
						use_dyad_rho = use_dyad_rho,
						bipartite = bip,
						symmetric = symmetric
					),
					error = function(e) NULL
				)
			} else {
				tryCatch(
					sample_beta_dynamic_cpp(
						Xdyn_list = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_static = if (length(beta_static) > 0L) beta_static else numeric(0),
						rho_by_coef = rho_by_coef,
						sigma_by_coef = sigma_by_coef,
						Lambda = Lambda_full,
						beta0_mean = beta0_mean,
						beta0_cov = beta0_cov,
						s2 = s2,
						dyad_rho = if (is.null(rho)) 0 else rho,
						bipartite = bip,
						symmetric = symmetric,
						use_dyad_rho = use_dyad_rho,
						Omega_in = if (!is.null(nc_st)) nc_st$Omega,
						h_in = if (!is.null(nc_st)) nc_st$h
					),
					error = function(e) NULL
				)
			}
			if (!is.null(ffbs_out) && all(is.finite(ffbs_out$path))) {
				beta_dyn_path <- ffbs_out$path
				if (sym_tie) beta_dyn_path <- .lame_sym_tie_state(
					beta_dyn_path,
					beta_static,
					beta_dyn$dynamic_idx,
					beta_dyn$static_idx,
					rb,
					cb
				)$path
				tryErrorChecks$beta <- tryErrorChecks$beta + (ffbs_out$chol_fail %||% 0L)
			} else {
				# on failure keep the current draw and count it
				tryErrorChecks$beta <- tryErrorChecks$beta + 1L
			}
			# conjugate update of the static coefficients
			if (length(beta_dyn$static_idx) > 0L) {
				# prior precision diag(1 / (g * s2)), the g-prior form
				p_static <- length(beta_dyn$static_idx)
				g_diag <- .lame_static_g_scale(g, beta_dyn$static_idx)
				prior_prec <- diag(1 / (g_diag * s2), p_static)
				ffbs_stat <- tryCatch(
					sample_beta_static_cpp(
						Xdyn_list = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_dyn_path = beta_dyn_path,
						prior_mean = rep(0, p_static),
						prior_prec = prior_prec,
						s2 = s2,
						dyad_rho = if (is.null(rho)) 0 else rho,
						bipartite = bip,
						symmetric = symmetric,
						use_dyad_rho = use_dyad_rho
					),
					error = function(e) NULL
				)
				if (!is.null(ffbs_stat) && all(is.finite(ffbs_stat$beta))) {
					beta_static <- as.numeric(ffbs_stat$beta)
					if (sym_tie) beta_static <- .lame_sym_tie_state(
						beta_dyn_path,
						beta_static,
						beta_dyn$dynamic_idx,
						beta_dyn$static_idx,
						rb,
						cb
					)$stat
					tryErrorChecks$beta <- tryErrorChecks$beta + (ffbs_stat$chol_fail %||% 0L)
				} else {
					tryErrorChecks$beta <- tryErrorChecks$beta + 1L
				}
			}
			# the p-length beta (time mean of the path) for summaries; beta_full_path feeds EZ next iteration
			beta_path_now <- matrix(0, N, length(beta))
			for (t in seq_len(N)) {
				beta_path_now[t, ] <- assemble_beta_path(
					beta_dyn_path[t, ],
					beta_static,
					beta_dyn$dynamic_idx,
					beta_dyn$static_idx,
					length(beta)
				)
			}
			beta <- colMeans(beta_path_now)
			# redraw a and b given the per-period beta path, from the residual cube Z - X beta - UV
			n_a__ <- if (bip) nA else n
			n_b__ <- if (bip) nB else n
			resid_cube <- array(0, dim = c(n_a__, n_b__, N))
			for (t in seq_len(N)) {
				if (bip) Xt_ <- Xlist[[t]] else Xt_ <- X[, , , t]
				if (length(dim(Xt_)) == 2L) Xt_ <- array(Xt_, c(n_a__, n_b__, 1))
				p_t_ <- dim(Xt_)[3]
				xb <- matrix(0, n_a__, n_b__)
				if (p_t_ > 0) {
					b_t_ <- beta_path_now[t, ]
					for (k in seq_len(p_t_)) if (b_t_[k] != 0) xb <- xb + b_t_[k] * Xt_[, , k]
				}
				uv_ <- if ((bip && RA > 0 && RB > 0) || (!bip && R > 0)) {
					Ut <- U_cube_eff[, , t]; Vt <- V_cube_eff[, , t]
					if (!is.matrix(Ut)) Ut <- matrix(Ut, n_a__, 1)
					if (!is.matrix(Vt)) Vt <- matrix(Vt, n_b__, 1)
					if (bip) Ut %*% G_eff %*% t(Vt) else Ut %*% t(Vt)
				} else 0
				resid_cube[, , t] <- Z[, , t] - xb - uv_
			}
			# draw a and b from the residual conditional; the block labels are passed on the sampler's layout
			# so the helper knows which block is the intercept
			ab_bases <- build_a_b_constraint_bases(
				beta_dyn,
				beta_dyn_labels,
				n_a__,
				n_b__,
				bip,
				symmetric,
				intercept = intercept
			)
			s2_floor <- s2
			sab_a_floor <- Sab[1, 1]
			# projection of a/b onto the constraint space
			.apply_H <- function(v, H) if (is.null(H)) v else as.numeric(H %*% (t(H) %*% v))
			if (!bip) {
				if (dynamic_ab) {
					# levels from the likelihood given the deviations, then the whole (a, b) panel from the AR(1) kernel with
					# the levels folded into the known predictor; the dyadic correlation enters when dcor is on
					dyad_rho_ab <- if (!symmetric && isTRUE(dcor)) .lame_rho_scalar(rho) else 0
					a_dev <- a_mat - a_lev
					b_dev <- b_mat - b_lev
					if (isTRUE(rvar) || (symmetric && isTRUE(nvar))) {
						cp_a <- if (symmetric) list(
							mean = NULL,
							var = sab_a_floor
						) else .lame_ab_cond_prior(Sab, 1L, b_lev)
						a_lev <- if (symmetric) {
							.lame_symmetric_level_draw(
								resid_cube,
								a_dev,
								Sab[1, 1],
								s2_floor,
								center = isTRUE(intercept)
							)
						} else .lame_ab_level_draw_lik(
							resid_cube,
							a_dev,
							b_mat,
							cp_a$var,
							s2_floor,
							1L,
							mean_level = cp_a$mean,
							center = isTRUE(intercept),
							rho = dyad_rho_ab,
							skip_diag = TRUE
						)
						a_mat <- a_lev + a_dev
					}
					if (symmetric) {
						b_lev <- a_lev
						b_mat <- a_mat
					} else if (isTRUE(cvar)) {
						cp_b <- .lame_ab_cond_prior(Sab, 2L, a_lev)
						b_lev <- .lame_ab_level_draw_lik(
							resid_cube,
							b_dev,
							a_mat,
							cp_b$var,
							s2_floor,
							2L,
							mean_level = cp_b$mean,
							center = isTRUE(intercept),
							rho = dyad_rho_ab,
							skip_diag = TRUE
						)
						b_mat <- b_lev + b_dev
					}
					lev_field <- outer(a_lev, b_lev, "+")
					lev_cube <- array(rep(as.numeric(lev_field), N), dim = c(n, n, N))
					ab_update <- tryCatch(
						sample_dynamic_ab_cpp(
							a_mat - a_lev,
							b_mat - b_lev,
							resid_cube,
							lev_cube,
							rho_ab,
							sigma_ab,
							s2_floor,
							symmetric,
							dyad_rho_ab,
							center = isTRUE(intercept),
							rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
							cvar = isTRUE(cvar)
						),
						error = function(e) {
							tryErrorChecks$betaAB <<- tryErrorChecks$betaAB + 1L
							NULL
						}
					)
					if (!is.null(ab_update)) {
						a_mat <- ab_update$a + a_lev
						b_mat <- if (symmetric) a_mat else ab_update$b + b_lev
						# disabled margins are exactly zero
						if (symmetric) {
							if (!isTRUE(nvar)) a_mat[is.finite(a_mat)] <- 0
							b_mat <- a_mat
						} else {
							if (!rvar) a_mat[is.finite(a_mat)] <- 0
							if (!cvar) b_mat[is.finite(b_mat)] <- 0
						}
						a_mat[!is.finite(a_mat)] <- 0
						b_mat[!is.finite(b_mat)] <- 0
						# identification against the intercept: the effects are centred in every period; a dynamic intercept
						# takes each period's constant, a static one their average
						if (isTRUE(intercept) && length(beta) > 0L) {
							jd <- match(1L, beta_dyn$dynamic_idx)
							js <- match(1L, beta_dyn$static_idx)
							cen_ab <- .lame_ab_center_dynamic(a_mat, b_mat)
							a_mat <- cen_ab$a_mat
							b_mat <- cen_ab$b_mat
							if (!is.na(jd)) {
								beta_dyn_path[, jd] <- beta_dyn_path[, jd] + cen_ab$shift
							} else if (!is.na(js)) {
								beta_static[js] <- beta_static[js] + mean(cen_ab$shift)
							}
							for (t in seq_len(N)) beta_path_now[t, ] <- assemble_beta_path(
								beta_dyn_path[t, ],
								beta_static,
								beta_dyn$dynamic_idx,
								beta_dyn$static_idx,
								length(beta)
							)
							beta <- colMeans(beta_path_now)
						} else if (!is.null(ab_bases$H_a) || !is.null(ab_bases$H_b)) {
							# without an intercept the effects still carry the contrast a dynamic nodal coefficient needs
							for (t in seq_len(N)) {
								a_mat[, t] <- .apply_H(a_mat[, t], ab_bases$H_a)
								b_mat[, t] <- if (symmetric) a_mat[, t] else .apply_H(b_mat[, t], ab_bases$H_b)
							}
						}
					}
				} else {
					# static a, b pooled over periods by the joint additive-effects sampler with an empty design; resid_cube
					# holds Z - X beta - U V'. the diagonal is refreshed first, since this Gaussian block reads it
					actor_diagonal <- matrix(a + b, n, N)
					old_diagonal <- .lame_diagonal(resid_cube)
					resid_cube <- .lame_refresh_diagonal(resid_cube, actor_diagonal, rho, s2)
					for (tt in seq_len(N)) {
						z_slice <- Z[, , tt]
						diag(z_slice) <- diag(z_slice) + diag(resid_cube[, , tt]) - old_diagonal[, tt]
						Z[, , tt] <- z_slice
					}
					Zr_ab <- resid_cube
					Zr_ab[!is.finite(Zr_ab)] <- 0
					# the same whitening rho as the standard rbeta_ab path; for symmetric fits it is ramped near one
					# so each unordered dyad counts once
					dyad_rho_eff <- if (symmetric || isTRUE(dcor)) .lame_rho_scalar(rho) else 0
					abd <- tryCatch(
						.lame_static_actor_draw(
							Zr_ab,
							Sab,
							dyad_rho_eff,
							s2_floor,
							symmetric = symmetric,
							rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
							cvar = if (symmetric) isTRUE(nvar) else isTRUE(cvar),
							center = isTRUE(intercept)
						),
						error = function(e) {
							tryErrorChecks$betaAB <<- tryErrorChecks$betaAB + 1L
							NULL
						}
					)
					if (!is.null(abd) && all(is.finite(abd$a)) && all(is.finite(abd$b))) {
						if (symmetric) {
							a <- ((as.numeric(abd$a) + as.numeric(abd$b)) / 2) * nvar
							a <- .apply_H(a, ab_bases$H_a)
							b <- a
						} else {
							a <- as.numeric(abd$a) * rvar
							b <- as.numeric(abd$b) * cvar
							a <- .apply_H(a, ab_bases$H_a)
							b <- .apply_H(b, ab_bases$H_b)
						}
					}
				}
			} else {
				# bipartite case (asymmetric by construction)
				if (dynamic_ab) {
					# per-period joint (a_t, b_t) draw with the bipartite rbeta_ab kernel (empty design), under the AR(1)
					# conditional prior, so the effects pool across time
					Xempty_bip1 <- list(array(0, dim = c(nA, nB, 0L)))
					for (t in seq_len(N)) {
						cond_a <- .lame_ab_ar1_cond(a_mat, t, N, rho_ab, sigma_ab, level = a_lev)
						cond_b <- .lame_ab_ar1_cond(b_mat, t, N, rho_ab, sigma_ab, level = b_lev)
						Zb1 <- resid_cube[, , t, drop = FALSE]
						Zb1[, , 1] <- Zb1[, , 1] - outer(cond_a$m, rep(1, nB)) -
							outer(rep(1, nA), cond_b$m)
						Zb1[!is.finite(Zb1)] <- 0
						abd_t <- tryCatch(
							rbeta_ab_bip_gibbs_cpp(
								Zb1,
								Xempty_bip1,
								matrix(0, nA, nB),
								a_mat[, t] - cond_a$m,
								b_mat[, t] - cond_b$m,
								s2,
								as.numeric(g[1]),
								cond_a$v,
								cond_b$v,
								rvar,
								cvar,
								center = isTRUE(intercept)
							),
							error = function(e) {
								tryErrorChecks$betaAB <<- tryErrorChecks$betaAB + 1L
								NULL
							}
						)
						if (!is.null(abd_t) && all(is.finite(abd_t$a)) &&
								all(is.finite(abd_t$b))) {
							a_mat[, t] <- .apply_H(cond_a$m + as.numeric(abd_t$a), ab_bases$H_a)
							b_mat[, t] <- .apply_H(cond_b$m + as.numeric(abd_t$b), ab_bases$H_b)
						}
					}
				} else {
					# static a, b (bipartite): the bipartite joint Gibbs sampler with an empty design and a zero
					# multiplicative offset, since resid_cube already holds Z - X beta - U G V'
					Zbip <- resid_cube
					Zbip[!is.finite(Zbip)] <- 0
					Xempty_bip <- lapply(
						seq_len(N),
						function(t) array(0, dim = c(nA, nB, 0L))
					)
					abd <- tryCatch(
						rbeta_ab_bip_gibbs_cpp(
							Zbip,
							Xempty_bip,
							matrix(0, nA, nB),
							if (length(a) == nA) a else rep(0, nA),
							if (length(b) == nB) b else rep(0, nB),
							s2,
							as.numeric(g[1]),
							Sab[1, 1],
							Sab[2, 2],
							rvar,
							cvar,
							center = isTRUE(intercept)
						),
						error = function(e) {
							tryErrorChecks$betaAB <<- tryErrorChecks$betaAB + 1L
							NULL
						}
					)
					if (!is.null(abd) && all(is.finite(abd$a)) && all(is.finite(abd$b))) {
						a <- as.numeric(abd$a)
						b <- as.numeric(abd$b)
						if (!is.null(ab_bases$H_a)) a <- as.numeric(ab_bases$H_a %*% (t(ab_bases$H_a) %*% a))
						if (!is.null(ab_bases$H_b)) b <- as.numeric(ab_bases$H_b %*% (t(ab_bases$H_b) %*% b))
					}
				}
			}
			if (dynamic_ab) {
				# hyperparameters given the totals, the levels integrated out
				if (ab_hyper_move) {
					hm <- .lame_ab_hyper_step(
						a_mat,
						b_mat,
						Sab,
						sigma_ab,
						rho_ab,
						ab_hyper_cx,
						ab_hyper_par,
						ab_hyper_tries,
						adapt = s <= burn,
						stat_var_cap = (prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
					)
					Sab <- hm$Sab
					sigma_ab <- hm$sigma
					rho_ab <- hm$rho
					ab_hyper_par <- hm$state
				}
				# levels given the paths; they stand in for the static effects wherever the sweep needs one
				ab_lev <- .lame_ab_level_draw(
					a_mat,
					b_mat,
					rho_ab,
					sigma_ab,
					Sab,
					symmetric = symmetric,
					center = isTRUE(intercept),
					rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
					cvar = isTRUE(cvar)
				)
				a_lev <- ab_lev$a
				b_lev <- ab_lev$b
				a <- a_lev; b <- b_lev
			}
			# additive-effect AR(1) hyperparameters every sweep
			if (dynamic_ab && N >= 2L &&
					!is.null(rho_ab) && !is.null(sigma_ab)) {
				hh <- .lame_ab_hyper_update(
					a_mat - a_lev,
					b_mat - b_lev,
					rho_ab,
					sigma_ab,
					prior,
					symmetric = symmetric,
					center = isTRUE(intercept),
					rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
					cvar = isTRUE(cvar),
					stat_var_cap = (prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
				)
				rho_ab <- hh$rho_ab
				sigma_ab <- hh$sigma_ab
			}
			# per-block rho_beta and sigma_beta; rw1 keeps rho at 1 and updates sigma_beta only
			is_rw1 <- identical(dynamic_beta_kind, "rw1")
			is_alt_kind <- dynamic_beta_kind %in% c("rw2", "matern32")
			if (N >= 2 && length(beta_dyn$dynamic_idx) > 0L) {
				# Lambda inverse for the C++ helpers
				Lambda_inv_full <- .lame_beta_scale_precision(Lambda_full)
				# hierarchical pooling updates the shared per-block prior
				pool_sigma_scale <- prior$sigma_beta_scale
				pool_rho_mean_vec <- rep(prior$rho_beta_mean, beta_dyn$n_groups)
				pool_rho_sd_vec <- rep(prior$rho_beta_sd, beta_dyn$n_groups)
				if (!identical(dynamic_beta_pool, "none") &&
						beta_dyn$n_groups >= 2L) {
					if (dynamic_beta_pool %in% c("sigma", "both")) {
						pool_sigma_scale <- .pool_update_sigma_scale(
							sigma_beta_by_group,
							alpha = prior$sigma_beta_shape,
							hyper_a = prior$sigma_beta_shape,
							hyper_b = prior$sigma_beta_scale
						)
					}
					if (dynamic_beta_pool %in% c("rho", "both") && !is_rw1 && !is_alt_kind) {
						hyper <- tryCatch(
							.pool_update_rho_hyper(
								rho_beta_by_group,
								prior_mu_mean = 0,
								prior_mu_sd = 10,
								current_tau = pool_rho_tau,
								rho_lower = prior$rho_beta_lower,
								rho_upper = prior$rho_beta_upper
							),
							error = function(e) {
								tryErrorChecks$pool <<- tryErrorChecks$pool + 1L
								NULL
							}
						)
						if (!is.null(hyper)) {
							pool_rho_mu <- hyper$mu
							pool_rho_tau <- hyper$tau
						}
					}
				}
				rho_new <- if (is_rw1 || is_alt_kind) {
					# rho is not meaningful for rw1/rw2/matern32; report it pinned at 1
					rep(1, beta_dyn$n_groups)
				} else if (dynamic_beta_pool %in% c("rho", "both") && beta_dyn$n_groups >= 2L) {
					.pool_sample_rho(
						beta_dyn_path,
						group_id,
						beta_dyn$n_groups,
						Lambda_full,
						sigma_by_coef,
						rho_beta_by_group,
						pool_rho_mu,
						pool_rho_tau,
						prior$rho_beta_lower,
						prior$rho_beta_upper,
						if (!is.null(time_index)) as.numeric(time_index) else seq_len(N)
					)
				} else tryCatch(
					(if (.nonunit_gaps) function(...) .lame_sample_rho_beta_gaps(
						...,
						time_positions = as.numeric(time_index)
					) else sample_rho_beta_cpp)(
						beta_path = beta_dyn_path,
						group_id = as.integer(group_id),
						n_groups = as.integer(beta_dyn$n_groups),
						Lambda_inv = Lambda_inv_full,
						sigma_by_coef = sigma_by_coef,
						rho_current = rho_beta_by_group,
						rho_prior_mean = pool_rho_mean_vec,
						rho_prior_sd = pool_rho_sd_vec,
						rho_lower = prior$rho_beta_lower,
						rho_upper = prior$rho_beta_upper),
					error = function(e) {
						tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
						NULL
					}
				)
				if (!is.null(rho_new) && all(is.finite(rho_new))) {
					rho_beta_by_group <- as.numeric(rho_new)
					names(rho_beta_by_group) <- beta_dyn$group_names
				}
				sigma_new <- if (is_alt_kind || .nonunit_gaps) {
					tryCatch(
						.sample_sigma_beta_alt(
							beta_path = beta_dyn_path,
							kind = dynamic_beta_kind,
							length_scale = prior$matern32_length_scale %||% NULL,
							prior_shape = prior$sigma_beta_shape,
							prior_scale = pool_sigma_scale,
							Lambda = Lambda_full,
							group_id = group_id,
							n_groups = beta_dyn$n_groups,
							time_positions = if (!is.null(time_index)) as.numeric(time_index) else NULL,
							rho_by_group = rho_beta_by_group
						),
						error = function(e) {
							tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
							NULL
						}
					)
				} else tryCatch(
					sample_sigma_beta_cpp(
						beta_path = beta_dyn_path,
						group_id = as.integer(group_id),
						n_groups = as.integer(beta_dyn$n_groups),
						Lambda_inv = Lambda_inv_full,
						rho_by_group = rho_beta_by_group,
						prior_shape = rep(prior$sigma_beta_shape, beta_dyn$n_groups),
						prior_scale = rep(pool_sigma_scale, beta_dyn$n_groups)
					),
					error = function(e) {
						tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
						NULL
					}
				)
				if (!is.null(sigma_new) && all(is.finite(sigma_new))) {
					sigma_beta_by_group <- as.numeric(sigma_new)
					names(sigma_beta_by_group) <- beta_dyn$group_names
				}
			}
		} else if (dynamic_ab) {
			if (bip) {
				# bipartite dynamic_ab: the beta draw is the full coefficient vector, so its residual keeps X beta in place
				# (the base is built with beta zeroed and only the UV term removed)
				zero_a <- matrix(0, nA, N)
				zero_b <- matrix(0, nB, N)
				EZ_uv_only <- get_EZ_dynamic_ab(
					Xlist,
					rep(0, length(beta)),
					zero_a,
					zero_b,
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N,
					bip = TRUE,
					G = G,
					nA = nA,
					nB = nB
				)
				# bipartite conjugate beta update
				p_bip <- length(beta)
				if (p_bip > 0) {
					resid_b <- Z - EZ_uv_only
					for (t in 1:N) {
						resid_b[, , t] <- resid_b[, , t] - a_mat[, t]
						for (j in 1:nB) resid_b[, j, t] <- resid_b[, j, t] - b_mat[j, t]
					}
					# actor paths and the multiplicative term are already removed: the exact coefficient conditional of the
					# static bipartite route, with both actor updates disabled
					beta <- as.numeric(
						rbeta_ab_bip_gibbs_cpp(
							resid_b,
							Xlist,
							matrix(0, nA, nB),
							numeric(nA),
							numeric(nB),
							s2,
							as.numeric(g[1]),
							1,
							1,
							FALSE,
							FALSE
						)$beta
					)
				}
				# recompute EZ_no_ab after beta changes
				EZ_no_ab <- get_EZ_dynamic_ab(
					Xlist,
					beta,
					zero_a,
					zero_b,
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N,
					bip = TRUE,
					G = G,
					nA = nA,
					nB = nB
				)
			} else if ((pr + pc + pd + intercept) > 0) {
				# unipartite dynamic_ab: coefficient draw given the actor paths
				EZ_no_ab <- get_EZ_unip_time(
					Xlist,
					beta,
					outer(rep(0, n), rep(0, n), "+"),
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N
				)
				# the coefficients are drawn conditional on the current additive effects, whose centring already pins them
				# to the intercept; the correlation enters only when the panel identifies it
				rho_eff <- if (symmetric || isTRUE(dcor)) .lame_rho_scalar(rho) else 0
				# a zero actor covariance keeps this coefficient update conditional
				a_fin <- a_mat
				b_fin <- b_mat
				a_fin[!is.finite(a_fin)] <- 0
				b_fin[!is.finite(b_fin)] <- 0
				Z_temp <- Z -
					aperm(array(a_fin, c(n, N, n)), c(1L, 3L, 2L)) -
					aperm(array(b_fin, c(n, N, n)), c(3L, 1L, 2L))
				if (.dyn_uv_has_rank) {
					for (tt in seq_len(N)) Z_temp[, , tt] <- Z_temp[, , tt] - U_cube[, , tt] %*% t(V_cube[, , tt])
				} else {
					Z_temp <- Z_temp - as.vector(U %*% t(V))
				}
				na_mask <- is.na(Z_temp)
				Z_temp[na_mask] <- 0
				g_num <- as.numeric(g[1])
				beta_ab_calc <- if (isTRUE(dynamic_rho)) {
					.lame_dynamic_rho_beta_ab_draw(
						Z_temp,
						mXLong,
						mXtLong,
						xxLong,
						matrix(0, 2L, 2L),
						rho,
						s2,
						g_num,
						center = isTRUE(intercept)
					)
				} else try(
					if (!is.null(sym_gaussian_context))
						.lame_symmetric_gaussian_draw(Z_temp, Sab[1L, 1L], sym_gaussian_context)
					else .lame_replicated_gaussian_draw(
						Z_temp,
						mXLong,
						matrix(0, 2L, 2L),
						as.numeric(rho_eff),
						s2,
						g_num,
						FALSE,
						XX = xxLong
					),
					silent = FALSE
				)
				if (!inherits(beta_ab_calc, 'try-error')) {
					beta <- c(beta_ab_calc$beta)
					if (sym_tie) beta <- .lame_sym_tie(beta, rb, cb)
				} else {
					tryErrorChecks$betaAB <- tryErrorChecks$betaAB + 1
				}
				EZ_no_ab <- get_EZ_unip_time(
					Xlist,
					beta,
					outer(rep(0, n), rep(0, n), "+"),
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N
				)
			} else {
				# no covariates: EZ_no_ab is just UV
				if (bip) {
					zero_a <- matrix(0, nA, N)
					zero_b <- matrix(0, nB, N)
					EZ_no_ab <- get_EZ_dynamic_ab(
						Xlist,
						beta,
						zero_a,
						zero_b,
						if (.dyn_uv_has_rank) U_cube else U,
						if (.dyn_uv_has_rank) V_cube else V,
						N,
						bip = TRUE,
						G = G,
						nA = nA,
						nB = nB
					)
				} else {
					EZ_no_ab <- get_EZ_unip_time(
						Xlist,
						beta,
						outer(rep(0, n), rep(0, n), "+"),
						if (.dyn_uv_has_rank) U_cube else U,
						if (.dyn_uv_has_rank) V_cube else V,
						N
					)
				}
			}
			# update dynamic a and b
			if (bip) {
				# bipartite: a (nA x N) and b (nB x N) are updated separately by conjugate normal Gibbs steps
				var_innov <- sigma_ab^2
				s2_floor_ab <- s2
				# levels from the likelihood given the deviations (see the unipartite branch)
				resid_lev <- Z - EZ_no_ab
				a_dev <- a_mat - a_lev
				b_dev <- b_mat - b_lev
				if (isTRUE(rvar)) {
					a_lev <- .lame_ab_level_draw_lik(
						resid_lev,
						a_dev,
						b_mat,
						Sab[1, 1],
						s2,
						1L,
						center = isTRUE(intercept)
					)
					a_mat <- a_lev + a_dev
				}
				if (isTRUE(cvar)) {
					b_lev <- .lame_ab_level_draw_lik(
						resid_lev,
						b_dev,
						a_mat,
						Sab[2, 2],
						s2,
						2L,
						center = isTRUE(intercept)
					)
					b_mat <- b_lev + b_dev
				}
				for (t in 1:N) {
					resid_t <- Z[, , t] - EZ_no_ab[, , t] - outer(a_lev, b_lev, "+")
					resid_t[!is.finite(resid_t)] <- NA_real_
					# AR(1) prior
					v_a <- numeric(nA)
					for (i in 1:nA) {
						if (!isTRUE(rvar)) { a_mat[i, t] <- 0; next }
						prev_a <- if (t > 1 && is.finite(a_mat[i, t - 1])) a_mat[i, t - 1] - a_lev[i] else 0
						next_a <- if (t < N && is.finite(a_mat[i, t + 1])) a_mat[i, t + 1] - a_lev[i] else 0
						if (N == 1L) {
							pprec <- 1 / var_innov
							pm <- 0
						} else if (t == 1L) {
							pprec <- 1 / var_innov
							pm <- rho_ab * next_a
						} else if (t == N) {
							pprec <- 1 / var_innov
							pm <- rho_ab * prev_a
						} else {
							pprec <- (1 + rho_ab^2) / var_innov
							pm <- rho_ab * (prev_a + next_a) / (1 + rho_ab^2)
						}
						row_residual <- resid_t[i, ] - (b_mat[, t] - b_lev)
						row_obs <- sum(is.finite(row_residual))
						row_sum <- sum(row_residual, na.rm = TRUE)
						dprec <- row_obs / s2_floor_ab
						post_prec <- pprec + dprec
						post_mean <- (pprec * pm + row_sum / s2_floor_ab) / post_prec
						if (is.finite(post_prec) && post_prec > 0 &&
								is.finite(post_mean)) {
							a_mat[i, t] <- a_lev[i] + rnorm(1, post_mean, 1 / sqrt(post_prec))
							v_a[i] <- 1 / post_prec
						} else {
							a_mat[i, t] <- a_lev[i] + pm
							tryErrorChecks$ab <- tryErrorChecks$ab + 1L
						}
					}
					# condition on the row effects summing to zero
					if (isTRUE(intercept)) a_mat[, t] <- a_lev + .lame_ab_center_weighted(a_mat[, t] - a_lev, v_a)
					# residual for the b step
					for (i in 1:nA) resid_t[i, ] <- resid_t[i, ] - (a_mat[i, t] - a_lev[i])
					v_b <- numeric(nB)
					for (j in 1:nB) {
						if (!isTRUE(cvar)) { b_mat[j, t] <- 0; next }
						prev_b <- if (t > 1 && is.finite(b_mat[j, t - 1])) b_mat[j, t - 1] - b_lev[j] else 0
						next_b <- if (t < N && is.finite(b_mat[j, t + 1])) b_mat[j, t + 1] - b_lev[j] else 0
						if (N == 1L) {
							pprec <- 1 / var_innov
							pm <- 0
						} else if (t == 1L) {
							pprec <- 1 / var_innov
							pm <- rho_ab * next_b
						} else if (t == N) {
							pprec <- 1 / var_innov
							pm <- rho_ab * prev_b
						} else {
							pprec <- (1 + rho_ab^2) / var_innov
							pm <- rho_ab * (prev_b + next_b) / (1 + rho_ab^2)
						}
						col_obs <- sum(is.finite(resid_t[, j]))
						col_sum <- sum(resid_t[, j], na.rm = TRUE)
						dprec <- col_obs / s2_floor_ab
						post_prec <- pprec + dprec
						post_mean <- (pprec * pm + col_sum / s2_floor_ab) / post_prec
						if (is.finite(post_prec) && post_prec > 0 &&
								is.finite(post_mean)) {
							b_mat[j, t] <- b_lev[j] + rnorm(1, post_mean, 1 / sqrt(post_prec))
							v_b[j] <- 1 / post_prec
						} else {
							b_mat[j, t] <- b_lev[j] + pm
							tryErrorChecks$ab <- tryErrorChecks$ab + 1L
						}
					}
					if (isTRUE(intercept)) b_mat[, t] <- b_lev + .lame_ab_center_weighted(b_mat[, t] - b_lev, v_b)
				}
			} else {
				# levels from the likelihood given the deviations, so persistent actor heterogeneity sits in the
				# additive effects rather than the multiplicative term
				resid_lev <- Z - EZ_no_ab
				# the dyadic correlation enters when dcor is on, as in the z draw
				dyad_rho_ab <- if (!symmetric && isTRUE(dcor)) rho else 0
				a_dev <- a_mat - a_lev
				b_dev <- b_mat - b_lev
				if (isTRUE(rvar) || (symmetric && isTRUE(nvar))) {
					cp_a <- if (symmetric) list(mean = NULL, var = Sab[1, 1]) else .lame_ab_cond_prior(Sab, 1L, b_lev)
					a_lev <- if (symmetric) {
						.lame_symmetric_level_draw(
							resid_lev,
							a_dev,
							Sab[1, 1],
							s2,
							center = isTRUE(intercept)
						)
					} else .lame_ab_level_draw_lik(
						resid_lev,
						a_dev,
						b_mat,
						cp_a$var,
						s2,
						1L,
						mean_level = cp_a$mean,
						center = isTRUE(intercept),
						rho = dyad_rho_ab,
						skip_diag = TRUE
					)
					a_mat <- a_lev + a_dev
				}
				if (symmetric) {
					b_lev <- a_lev
					b_mat <- a_mat
				} else if (isTRUE(cvar)) {
					cp_b <- .lame_ab_cond_prior(Sab, 2L, a_lev)
					b_lev <- .lame_ab_level_draw_lik(
						resid_lev,
						b_dev,
						a_mat,
						cp_b$var,
						s2,
						2L,
						mean_level = cp_b$mean,
						center = isTRUE(intercept),
						rho = dyad_rho_ab,
						skip_diag = TRUE
					)
					b_mat <- b_lev + b_dev
				}
				# the kernel draws a zero-mean AR(1) state, so the levels are folded into the known predictor
				lev_field <- outer(a_lev, b_lev, "+")
				EZ_lev <- EZ_no_ab
				for (tt in seq_len(N)) EZ_lev[, , tt] <- EZ_no_ab[, , tt] + lev_field
				ab_update <- .lame_sample_dynamic_ab_period(
					a_mat - a_lev,
					b_mat - b_lev,
					Z,
					EZ_lev,
					rho_ab,
					sigma_ab,
					s2,
					symmetric,
					dyad_rho_ab,
					center = isTRUE(intercept),
					rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
					cvar = isTRUE(cvar)
				)
				a_mat <- ab_update$a + a_lev
				b_mat <- ab_update$b + b_lev
				# the kernel has no rvar / cvar switch, so a disabled margin is zeroed here
				if (!rvar) a_mat[is.finite(a_mat)] <- 0
				if (!cvar || (symmetric && !rvar)) b_mat[is.finite(b_mat)] <- 0
			}
			if (any(!is.finite(a_mat))) {
				a_mat[!is.finite(a_mat)] <- 0
				tryErrorChecks$ab <- tryErrorChecks$ab + 1L
			}
			if (any(!is.finite(b_mat))) {
				b_mat[!is.finite(b_mat)] <- 0
				tryErrorChecks$ab <- tryErrorChecks$ab + 1L
			}
			# centre the effects in every period; the static intercept takes the average of the removed constants
			if (isTRUE(intercept) && length(beta) > 0L && is.finite(beta[1L])) {
				cen_ab <- .lame_ab_center_dynamic(a_mat, b_mat)
				a_mat <- cen_ab$a_mat
				b_mat <- cen_ab$b_mat
				beta[1L] <- beta[1L] + mean(cen_ab$shift)
			}
			# hyperparameters given the totals, the levels integrated out
			if (ab_hyper_move) {
				hm <- .lame_ab_hyper_step(
					a_mat,
					b_mat,
					Sab,
					sigma_ab,
					rho_ab,
					ab_hyper_cx,
					ab_hyper_par,
					ab_hyper_tries,
					adapt = s <= burn,
					stat_var_cap = (prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
				)
				Sab <- hm$Sab
				sigma_ab <- hm$sigma
				rho_ab <- hm$rho
				ab_hyper_par <- hm$state
			}
			# levels given the paths; the deviations feed the AR(1) hyperparameters below
			ab_lev <- .lame_ab_level_draw(
				a_mat,
				b_mat,
				rho_ab,
				sigma_ab,
				Sab,
				symmetric = symmetric,
				center = isTRUE(intercept),
				rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
				cvar = isTRUE(cvar)
			)
			a_lev <- ab_lev$a
			b_lev <- ab_lev$b
			a <- a_lev
			b <- b_lev
			# AR(1) parameters every sweep, burn-in included
			if (N >= 2L) {
				hh <- .lame_ab_hyper_update(
					a_mat - a_lev,
					b_mat - b_lev,
					rho_ab,
					sigma_ab,
					prior,
					symmetric = symmetric,
					center = isTRUE(intercept),
					rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
					cvar = isTRUE(cvar),
					stat_var_cap = (prior$ab_var_cap %||% Inf) * .lame_cap_unit(family, s2)
				)
				rho_ab <- hh$rho_ab
				sigma_ab <- hh$sigma_ab
			}
		} else {
			# static update: subtract each period's current multiplicative offset before the draw
			Z_static_ab <- Z
			if (.dyn_uv_has_rank) {
				for (tt in seq_len(N)) {
					Ut_ab <- matrix(U_cube[, , tt], nrow = dim(U_cube)[1L])
					Vt_ab <- matrix(V_cube[, , tt], nrow = dim(V_cube)[1L])
					Z_static_ab[, , tt] <- Z_static_ab[, , tt] -
						if (bip) Ut_ab %*% G %*% t(Vt_ab) else tcrossprod(Ut_ab, Vt_ab)
				}
			} else if (bip) {
				Ut_ab <- if (length(dim(U)) == 3L) matrix(U[, , 1L], nrow = dim(U)[1L]) else U
				Vt_ab <- if (length(dim(V)) == 3L) matrix(V[, , 1L], nrow = dim(V)[1L]) else V
				if (RA > 0 && RB > 0) for (tt in seq_len(N))
					Z_static_ab[, , tt] <- Z_static_ab[, , tt] - Ut_ab %*% G %*% t(Vt_ab)
			} else Z_static_ab <- Z - as.vector(tcrossprod(U, V))
			if (!is.null(sym_gaussian_context)) {
				beta_ab_calc <- try(
					.lame_symmetric_gaussian_draw(
						Z_static_ab,
						Sab[1L, 1L],
						sym_gaussian_context
					),
					silent = FALSE
				)
			} else if (bip) {
				# bipartite Gibbs update in C++
				beta_ab_calc <- tryCatch(
					rbeta_ab_bip_gibbs_cpp(
						Z_static_ab,
						Xlist_tmp,
						matrix(0, nA, nB),
						if (length(a) == nA) a else rep(0, nA),
						if (length(b) == nB) b else rep(0, nB),
						s2,
						as.numeric(g[1]),
						Sab[1, 1],
						Sab[2, 2],
						rvar,
						cvar,
						center = isTRUE(intercept)
					),
					error = function(e) {
						list(
							beta = beta,
							a = if (length(a) == nA) a else rep(0, nA),
							b = if (length(b) == nB) b else rep(0, nB)
						)
					}
				)
			} else if ((pr + pc + pd + intercept) > 0) {
				# the sampled correlation enters the effects draw whenever dcor is on
				rho_eff <- if (symmetric || isTRUE(dcor)) .lame_rho_scalar(rho) else 0
				Sab_active <- .lame_sab_gls(Sab, symmetric)
				if (!isTRUE(rvar)) Sab_active[1L, ] <- Sab_active[, 1L] <- 0
				if (!isTRUE(cvar)) Sab_active[2L, ] <- Sab_active[, 2L] <- 0
				# zero the missing cells of Z before passing it to C++
				Z_temp <- Z_static_ab
				na_mask <- is.na(Z_temp)
				Z_temp[na_mask] <- 0
				g_num <- as.numeric(g[1])
				beta_ab_calc <- if (isTRUE(dynamic_rho)) {
					.lame_dynamic_rho_beta_ab_draw(
						Z_temp,
						mXLong,
						mXtLong,
						xxLong,
						Sab_active,
						rho,
						s2,
						g_num,
						center = isTRUE(intercept)
					)
				} else try(
					.lame_replicated_gaussian_draw(
						Z_temp,
						mXLong,
						Sab_active,
						as.numeric(rho_eff),
						s2,
						g_num,
						isTRUE(intercept),
						XX = xxLong
					),
					silent = FALSE
				)
			} else {
				# no coefficients to draw: the effects alone
				beta_ab_calc <- if (isTRUE(dynamic_rho)) {
					Sab_active <- Sab
					if (!isTRUE(rvar)) Sab_active[1L, ] <- Sab_active[, 1L] <- 0
					if (!isTRUE(cvar)) Sab_active[2L, ] <- Sab_active[, 2L] <- 0
					.lame_dynamic_rho_beta_ab_draw(
						Z_static_ab,
						array(numeric(0), c(n * n, 0L, N)),
						array(numeric(0), c(n * n, 0L, N)),
						array(numeric(0), c(0L, 0L, N)),
						Sab_active,
						rho,
						s2,
						1,
						center = isTRUE(intercept)
					)
				} else try(
					.lame_static_actor_draw(
						Z_static_ab,
						Sab,
						if (symmetric || isTRUE(dcor)) .lame_rho_scalar(rho) else 0,
						s2,
						symmetric = symmetric,
						rvar = if (symmetric) isTRUE(nvar) else isTRUE(rvar),
						cvar = if (symmetric) isTRUE(nvar) else isTRUE(cvar),
						center = isTRUE(intercept)
					),
					silent = FALSE
				)
			}
			if (!inherits(beta_ab_calc, 'try-error')) {
				beta <- c(beta_ab_calc$beta)
				if (sym_tie) beta <- .lame_sym_tie(beta, rb, cb)
				a <- c(beta_ab_calc$a) * rvar
				b <- c(beta_ab_calc$b) * cvar
				if (symmetric) { a <- b <- (a + b) / 2 }
			} else {
				# count a failed betaAB draw
				tryErrorChecks$betaAB <- tryErrorChecks$betaAB + 1
			}
		}

		# update Sab
		if (bip) {
			# under dynamic effects the pooled variance is that of the actor levels
			a_sab <- if (dynamic_ab && exists("a_lev")) a_lev else a
			b_sab <- if (dynamic_ab && exists("b_lev")) b_lev else b
			a_sab <- a_sab[is.finite(a_sab)]
			b_sab <- b_sab[is.finite(b_sab)]
			if (rvar) {
				Sab[1, 1] <- .lame_safe_inv_gamma(
					shape = (prior$etaab + length(a_sab) - as.integer(isTRUE(intercept) && length(a_sab) > 0L)) / 2,
					rate = (prior$etaab * prior$Sab0[1, 1] + sum(a_sab^2)) / 2
				)
			}
			if (cvar) {
				Sab[2, 2] <- .lame_safe_inv_gamma(
					shape = (prior$etaab + length(b_sab) - as.integer(isTRUE(intercept) && length(b_sab) > 0L)) / 2,
					rate = (prior$etaab * prior$Sab0[2, 2] + sum(b_sab^2)) / 2
				)
			}
			Sab[1, 2] <- Sab[2, 1] <- 0 # bipartite row and column effects are separate
		} else {
			# Sab0 is passed as is: rSab_fc scales the inverse Wishart by eta0 * Sab0 itself (the amen convention)
			Sab <- rSab_fc(
				a,
				b,
				Sab0 = prior$Sab0,
				eta0 = prior$etaab,
				tiny = 1e-8 * unit_var,
				rvar = rvar,
				cvar = cvar,
				symmetric = symmetric,
				center = isTRUE(intercept)
			)
		}
		if (joint_levels_enabled && s %% 10L == 0L) {
			ade <- a_mat - a_lev; bde <- b_mat - b_lev
			joint <- tryCatch(
				.lame_joint_levels_draw(
					Z,
					ade,
					bde,
					Sab,
					if (isTRUE(dcor)) .lame_rho_scalar(rho) else 0,
					joint_levels_context
				),
				error = function(e) {
					tryErrorChecks$joint_levels <<- tryErrorChecks$joint_levels + 1L
					NULL
				}
			)
			joint_levels_attempts <- joint_levels_attempts + 1L
			if (!is.null(joint)) {
				beta <- joint$beta; a <- a_lev <- joint$a; b <- b_lev <- joint$b
				a_mat <- a_lev + ade; b_mat <- b_lev + bde
			}
		}
		# update rho (asymmetric models only)
		if (dcor && !symmetric) {
			# refresh the per-period X beta, since beta may have moved this iteration
			if (beta_dyn$any) {
				beta_full_path_iter <- .lame_beta_full_path()
				Xbeta_cube_iter <- .lame_xbeta_cube(beta_full_path_iter)
			}
			beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta
			if (dynamic_ab) {
				E.T <- Z - get_EZ_dynamic_ab(
					Xlist,
					beta_for_EZ,
					a_mat,
					b_mat,
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N,
					bip = bip,
					G = G,
					nA = nA,
					nB = nB
				)
			} else {
				E.T <- Z - get_EZ_unip_time(
					Xlist,
					beta_for_EZ,
					outer(a, b, "+"),
					if (.dyn_uv_has_rank) U_cube else U,
					if (.dyn_uv_has_rank) V_cube else V,
					N
				)
			}
			if (beta_dyn$any && !is.null(Xbeta_cube_iter)) E.T <- E.T - Xbeta_cube_iter
			if (nom_rho && s %% nom_rho_every == 0L && isTRUE(s2 == 1)) {
				# collapsed rho update with the non-nominations integrated out, on the complete latent array including
				# the actor-slope offset; the cells come back without the offset
				nom_ez <- Z - E.T
				nom_pa <- if (per_actor_active) per_actor_offset else NULL
				nom_r <- .lame_nom_rho_collapsed(
					if (is.null(nom_pa)) Z else Z + nom_pa,
					if (is.null(nom_pa)) nom_ez else nom_ez + nom_pa,
					.lame_rho_scalar(rho),
					Y,
					family,
					if (identical(family, "frn")) YL else NULL,
					odmax
				)
				nom_rho_stats <- nom_rho_stats + c(1, 4, nom_r$accepts)
				if (nom_r$accepts > 0L) {
					Z <- if (is.null(nom_pa)) nom_r$Z else nom_r$Z - nom_pa
					E.T <- Z - nom_ez
					rho <- nom_r$rho
				}
			}
			if (isTRUE(dynamic_rho)) {
				rho_draw <- try({
					rho_path_new <- .lame_expand_rho(rho, N)
					for (t in seq_len(N)) {
						E_slice <- array(E.T[, , t], dim = c(dim(E.T)[1:2], 1L))
						rho_path_new[t] <- rrho_mh_rep_cpp(
							E_slice,
							rho_path_new[t],
							s2,
							prior_power = .lame_rho_prior_power()
						)
					}
					.lame_expand_rho(rho_path_new, N)
				}, silent = TRUE)
			} else {
				# several MH sweeps per scan: rho mixes slowly and the kernel is cheap
				rho_draw <- try({
					rho_tmp <- rho
					for (k_rho in 1:4) rho_tmp <- rrho_mh_rep_cpp(
						E.T,
						rho_tmp,
						s2,
						prior_power = .lame_rho_prior_power()
					)
					rho_tmp
				}, silent = TRUE)
			}
			if (!inherits(rho_draw, 'try-error')) { rho <- rho_draw } else { tryErrorChecks$rho <- tryErrorChecks$rho + 1 }
			# the rho update also refreshes the diagonal, using the current mean in E.T
			Z <- .lame_refresh_diagonal(
				Z,
				.lame_diagonal(Z) - .lame_diagonal(E.T),
				rho,
				s2
			)
		}

		# update U, V (bipartite uses RA / RB)
		if ((!bip && R > 0) || (bip && (RA > 0 || RB > 0))) {
			# refresh the per-period X beta cube when dynamic_beta is on
			if (beta_dyn$any) {
				beta_full_path_iter <- .lame_beta_full_path()
				Xbeta_cube_iter <- .lame_xbeta_cube(beta_full_path_iter)
			}
			beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta
			if (dynamic_ab) {
				# dynamic ab: time-varying effects with UV zeroed
				U_zero <- if (bip) array(0, dim = c(nA, max(1, RA), N)) else U * 0
				V_zero <- if (bip) array(0, dim = c(nB, max(1, RB), N)) else V * 0
				E <- Z - get_EZ_dynamic_ab(
					Xlist,
					beta_for_EZ,
					a_mat,
					b_mat,
					U_zero,
					V_zero,
					N,
					bip = bip,
					G = G,
					nA = nA,
					nB = nB
				)
			} else if (bip) {
				U_zero <- array(0, dim = c(nA, max(1, ncol(U)), N))
				V_zero <- array(0, dim = c(nB, max(1, ncol(V)), N))
				# base cube with X beta and a/b effects (X beta added per period under dynamic_beta)
				base_cube <- array(0, dim = c(nA, nB, N))
				if (!beta_dyn$any) {
					for (t in 1:N) {
						if (length(beta) > 0) {
							base_mat <- matrix(0, nA, nB)
							n_covs <- min(length(beta), dim(Xlist_tmp[[t]])[3])
							for (i in 1:n_covs) {
								base_mat <- base_mat + beta[i] * Xlist_tmp[[t]][, , i]
							}
							base_cube[, , t] <- base_mat
						}
					}
				}
				# time-expanded a and b matrices
				a_mat_temp <- matrix(a[1:nA], nA, N)
				b_mat_temp <- matrix(b[1:nB], nB, N)
				E <- Z - get_EZ_bip_cpp(base_cube, a_mat_temp, b_mat_temp, U_zero, V_zero, G)
			} else {
				E <- Z - get_EZ_unip_time(
					Xlist,
					beta_for_EZ,
					outer(a, b, "+"),
					U * 0,
					V * 0,
					N
				)
			}
			if (beta_dyn$any && !is.null(Xbeta_cube_iter)) E <- E - Xbeta_cube_iter
			if (!bip && !isTRUE(dcor) && isTRUE(dynamic_ab)) {
				# redraw the diagonal under the new actor effects before the UV kernel reads it (the rho update does
				# this when enabled)
				old_diagonal <- .lame_diagonal(E)
				uv_diagonal <- matrix(0, n, N)
				for (tt in seq_len(N)) {
					ut <- if (.dyn_uv_has_rank) matrix(U_cube[, , tt], n, R) else U
					vt <- if (.dyn_uv_has_rank) matrix(V_cube[, , tt], n, R) else V
					uv_diagonal[, tt] <- rowSums(ut * vt)
				}
				E <- .lame_refresh_diagonal(E, uv_diagonal, rho, s2)
				new_diagonal <- .lame_diagonal(E)
				for (tt in seq_len(N)) {
					z_slice <- Z[, , tt]
					diag(z_slice) <- diag(z_slice) + new_diagonal[, tt] - old_diagonal[, tt]
					Z[, , tt] <- z_slice
				}
			}
			shrink <- (s > .5 * burn)

			if (dynamic_uv) {
				uv_max_abs_iter <- as.numeric(prior$uv_max_abs %||% Inf)
				if (!is.finite(uv_max_abs_iter) || uv_max_abs_iter <= 0) {
					uv_max_abs_iter <- Inf
				}
				if (bip) {
					if (!identical(dynamic_uv_kind, "ar1")) {
						cli::cli_abort(c(
							"Internal dynamic-UV dispatch error.",
							"x" = "Bipartite dynamic UV can only use {.val ar1}; got {.val {dynamic_uv_kind}}."
						))
					}
					U_prev <- U_cube
					V_prev <- V_cube
					UV_try <- tryCatch(
						.lame_uv_sign_ar1(
							.lame_uv_svd_sign_ar1(
								rUV_dynamic_bip_fc_cpp(
									U_cube,
									V_cube,
									E,
									G,
									rho_uv,
									sigma_uv,
									s2,
									max_abs = uv_max_abs_iter
								),
								G,
								rho_uv,
								sigma_uv,
								max_abs = uv_max_abs_iter
							),
							rho_uv,
							sigma_uv
						),
						error = function(e) NULL
					)
					if (!is.null(UV_try) &&
							.lame_uv_ok(UV_try$U, UV_try[["V"]], uv_max_abs_iter)) {
						if ((UV_try$bound_rejections %||% 0L) > 0L) {
							uv_rail_count <- uv_rail_count + 1L
							tryErrorChecks$UV <- tryErrorChecks$UV + 1L
						}
						U_cube <- UV_try$U; V_cube <- UV_try[["V"]]
						U <- apply(U_cube, c(1, 2), mean)
						V <- apply(V_cube, c(1, 2), mean)
					} else {
						tryErrorChecks$UV <- tryErrorChecks$UV + 1
						if (!is.null(UV_try) && all(is.finite(UV_try$U)) &&
								all(is.finite(UV_try[["V"]]))) {
							stop(
								"Internal bipartite UV update left its coordinate support.", call. = FALSE
							)
						}
						U_cube <- U_prev; V_cube <- V_prev
						UV_try <- NULL # the AR(1) guard treats this as a failure
					}
					# after an accepted position draw: path moves, then the AR(1) hyperparameters every sweep
					if (!is.null(UV_try)) {
						if (npm_on && s %% npm_every == 0L && !is.null(G) && all(dim(G) == c(RA, RB))) {
							# the paths with the latent rows and columns they bound (R/nomination_path_move.R)
							npm_o <- .lame_npm_visit(
								Z,
								Y,
								family,
								if (identical(family, "frn")) YL else NULL,
								odmax,
								U_cube,
								V_cube,
								G,
								rho_uv,
								sigma_uv,
								max_abs = uv_max_abs_iter,
								ctx = npm_ctx
							)
							Z <- array(npm_o$Z, dim(Z)); U_cube <- npm_o$U; V_cube <- npm_o$V
							U <- apply(U_cube, c(1, 2), mean)
							V <- apply(V_cube, c(1, 2), mean)
							npm_stats[c("visits", names(npm_o$stats))] <- npm_stats[c("visits", names(npm_o$stats))] + c(1, npm_o$stats)
							rm(npm_o)
						}
						if (uv_nc_move && !is.null(G) && all(dim(G) == c(RA, RB))) {
							# scale ridge of the paths and the interaction matrix (R/uv_nc_step.R)
							nc_o <- .lame_uv_bip_scale_step(
								U_cube,
								V_cube,
								G,
								rho_uv,
								sigma_uv,
								prior$sigma_uv_shape,
								prior$sigma_uv_scale,
								sigma_max = prior$sigma_uv_max %||% 10,
								w_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2),
								g_var = weight_prior_var,
								g_cap = prior$G_max %||% g_rail_default,
								max_abs = uv_max_abs_iter
							)
							uv_nc_stats[c("bip_calls", "bip_moved")] <- uv_nc_stats[c("bip_calls", "bip_moved")] + c(1, nc_o$c != 1)
							if (nc_o$c != 1) {
								U_cube <- nc_o$U; V_cube <- nc_o$V; G <- nc_o$G; sigma_uv <- nc_o$sigma
								U <- apply(U_cube, c(1, 2), mean)
								V <- apply(V_cube, c(1, 2), mean)
							}
							rm(nc_o)
						}
						pmean_uv <- prior$rho_uv_mean %||% 0.9
						psd_uv <- prior$rho_uv_sd %||% 0.1
						# persistence and stationary variance as a block (R/uv_rho_w_step.R)
						uvh <- .lame_uv_rho_w_step(
							U_cube,
							V_cube,
							rho_uv,
							sigma_uv,
							pmean_uv,
							psd_uv,
							prior$sigma_uv_shape,
							prior$sigma_uv_scale,
							sigma_max = prior$sigma_uv_max %||% 10,
							stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2)
						)
						rho_uv <- uvh$rho
						sigma_uv <- uvh$sigma
					}
				} else {
					U_prev <- U_cube
					V_prev <- V_cube
					uv_update_ok <- FALSE
					# positions are drawn under the pair-whitened observation model at the sampled rho, so reciprocity
					# is carried by rho rather than absorbed into the coordinates
					rho_dyad_uv <- if (!symmetric && isTRUE(dcor))
						.lame_rho_scalar(rho) else 0
					if (identical(dynamic_uv_kind, "snap")) {
						if (is.null(delta_uv)) { delta_uv <- matrix(0, n, N) }
						if (is.null(delta_uv_v)) { delta_uv_v <- matrix(0, n, N) }
						UV <- try(
							rUV_dynamic_snap_fc_cpp(
								U_cube,
								V_cube,
								E,
								rho_uv,
								sigma_uv,
								s2,
								.lame_snap_kappa(prior$snap_kappa, rho_uv, sigma_uv),
								pi_snap,
								delta_uv,
								delta_uv_v,
								shrink,
								symmetric,
								rho_dyad_uv,
								max_abs = uv_max_abs_iter
							),
							silent = FALSE
						)
					} else if (identical(dynamic_uv_kind, "t")) {
						UV <- try(
							rUV_dynamic_t_fc_cpp(
								U_cube,
								V_cube,
								E,
								rho_uv,
								sigma_uv,
								s2,
								prior$t_nu,
								lambda_u,
								lambda_v,
								shrink,
								symmetric,
								rho_dyad_uv,
								max_abs = uv_max_abs_iter
							),
							silent = FALSE
						)
					} else {
						UV <- try(
							.lame_uv_orientation_ar1(
								rUV_dynamic_fc_cpp(
									U_cube,
									V_cube,
									E,
									rho_uv,
									sigma_uv,
									s2,
									shrink,
									symmetric,
									rho_dyad_uv,
									max_abs = uv_max_abs_iter
								),
								rho_uv,
								sigma_uv,
								enabled = (identical(family, "normal") || dynamic_joint_move) && !symmetric && R == 2L,
								align = shrink,
								max_abs = uv_max_abs_iter
							),
							silent = FALSE
						)
					}
					if (!inherits(UV, "try-error") && identical(dynamic_uv_kind, "snap")) {
						UV <- .lame_uv_sign_ar1(
							UV,
							rho_uv,
							sigma_uv,
							1 - UV$delta_u,
							1 - UV$delta_v
						)
					} else if (!inherits(UV, "try-error") && identical(dynamic_uv_kind, "t")) {
						UV <- .lame_uv_sign_ar1(
							UV,
							rho_uv,
							sigma_uv,
							UV$lambda_u,
							UV$lambda_v
						)
					}

					if (inherits(UV, 'try-error')) {
						UV <- list(U = U_cube, V = V_cube)
						tryErrorChecks$UV <- tryErrorChecks$UV + 1
					} else if (!.lame_uv_ok(UV$U, UV[["V"]], uv_max_abs_iter)) {
						if (all(is.finite(UV$U)) && all(is.finite(UV[["V"]]))) {
							stop("Internal dynamic UV update left its coordinate support.", call. = FALSE)
						}
						U_cube <- U_prev
						V_cube <- V_prev
						UV <- list(U = U_prev, V = V_prev)
						tryErrorChecks$UV <- tryErrorChecks$UV + 1
					} else {
						if ((UV$bound_rejections %||% 0L) > 0L) {
							uv_rail_count <- uv_rail_count + 1L
							tryErrorChecks$UV <- tryErrorChecks$UV + 1L
						}
						uv_update_ok <- TRUE
						U_cube <- UV$U
						V_cube <- UV[["V"]]
						if (symmetric && !is.null(L_dyn)) {
							# diagonal of the symmetric product, then the mirrored coordinates the directed kernel expects
							L_dyn <- .lame_sym_dyn_L_draw(E, U_cube, L_dyn, s2, weight_prior_var)
							for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
						}
						U <- apply(U_cube, c(1, 2), mean)
						V <- apply(V_cube, c(1, 2), mean)
						if (identical(dynamic_uv_kind, "snap")) {
							delta_uv <- UV$delta_u
							delta_uv_v <- UV$delta_v
							n_snap <- sum(delta_uv[, -1]) +
								if (symmetric) 0 else sum(delta_uv_v[, -1])
							n_elig <- (if (symmetric) 1 else 2) * nrow(delta_uv) * (N - 1)
							pi_snap <- stats::rbeta(
								1,
								prior$snap_pi_a + n_snap,
								prior$snap_pi_b + (n_elig - n_snap)
							)
						}
						if (identical(dynamic_uv_kind, "t")) {
							lambda_u <- UV$lambda_u
							lambda_v <- UV$lambda_v
						}
					}
					# after an accepted position draw: path moves, then the AR(1) hyperparameters every sweep, with the
					# sigma_uv draw capped by prior$uv_var_cap as on the bipartite path
					if (isTRUE(uv_update_ok)) {
						if (npm_on && s %% npm_every == 0L) {
							# the paths with the latent rows and columns they bound (R/nomination_path_move.R)
							npm_o <- .lame_npm_visit(
								Z,
								Y,
								family,
								if (identical(family, "frn")) YL else NULL,
								odmax,
								U_cube,
								V_cube,
								NULL,
								rho_uv,
								sigma_uv,
								max_abs = uv_max_abs_iter,
								ctx = npm_ctx
							)
							Z <- array(npm_o$Z, dim(Z)); U_cube <- npm_o$U; V_cube <- npm_o$V
							U <- apply(U_cube, c(1, 2), mean)
							V <- apply(V_cube, c(1, 2), mean)
							npm_stats[c("visits", names(npm_o$stats))] <- npm_stats[c("visits", names(npm_o$stats))] + c(1, npm_o$stats)
							rm(npm_o)
						}
						pmean_uv <- prior$rho_uv_mean %||% 0.9
						psd_uv <- prior$rho_uv_sd %||% 0.1
						if (symmetric && !is.null(L_dyn)) {
							# symmetric gauge: rho_uv from the transitions of the unit-variance paths; the innovation sd follows
							rho_uv <- if (dpm_sym) .lame_dpm_rho_unit_mh(U_cube, rho_uv, pmean_uv, psd_uv)
							else .lame_rho_uv_unit_mh(U_cube, rho_uv, pmean_uv, psd_uv)
							sigma_uv <- .lame_sym_uv_innov_sd(rho_uv)
							if (uv_nc_move) {
								# persistence given the standardised innovations, the paths following it (R/uv_nc_step.R)
								nc_o <- if (is.null(uv_nc_cx)) .lame_uv_sym_nc_rho(
									U_cube,
									L_dyn,
									E,
									s2,
									rho_uv,
									pmean_uv,
									psd_uv,
									max_abs = uv_max_abs_iter
								)
								else .lame_uv_sym_nc_rho_carry(
									U_cube,
									L_dyn,
									Z,
									uv_nc_cx,
									rho_uv,
									pmean_uv,
									psd_uv,
									max_abs = uv_max_abs_iter
								)
								uv_nc_stats[c("sym_calls", "sym_moved")] <- uv_nc_stats[c("sym_calls", "sym_moved")] + c(1, nc_o$rho != rho_uv)
								if (nc_o$rho != rho_uv) {
									rho_uv <- nc_o$rho
									U_cube <- nc_o$U
									if (!is.null(uv_nc_cx)) Z <- array(nc_o$Z, dim(Z))
									sigma_uv <- .lame_sym_uv_innov_sd(rho_uv)
									for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
									U <- apply(U_cube, c(1, 2), mean)
									V <- apply(V_cube, c(1, 2), mean)
								}
								rm(nc_o)
							}
							if (dpm_sym) {
								# scale ridge: the product, hence the predictor and the latent array, is unchanged
								dpm_rg <- .lame_dpm_sym_scale(
									U_cube,
									L_dyn,
									rho_uv,
									max_abs = uv_max_abs_iter,
									l_var = weight_prior_var
								)
								U_cube <- dpm_rg$U; L_dyn <- dpm_rg$Ld
								for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
								U <- apply(U_cube, c(1, 2), mean)
								V <- apply(V_cube, c(1, 2), mean)
								dpm_sym_stats <- dpm_sym_stats + dpm_rg$stats
								rm(dpm_rg)
							}
							if (!is.null(uv_temper_dt) && s %% uv_temper_every == 0L) {
								# tempered transition of one direction, the latent array carried (R/uv_tempered_move.R)
								tm_o <- .lame_uv_sym_temper(
									U_cube,
									L_dyn,
									Z,
									uv_temper_dt,
									sample.int(R, 1L),
									rho_uv,
									weight_prior_var,
									uv_temper_bs,
									cycles = 1L,
									max_abs = uv_max_abs_iter
								)
								uv_temper_stats <- uv_temper_stats + c(1, tm_o$accepted, isTRUE(tm_o$flipped))
								if (tm_o$accepted) {
									U_cube <- tm_o$U; L_dyn <- tm_o$L; Z <- array(tm_o$Z, dim(Z))
									for (t in seq_len(N)) V_cube[, , t] <- matrix(U_cube[, , t], n, R) %*% diag(L_dyn, R)
									U <- apply(U_cube, c(1, 2), mean)
									V <- apply(V_cube, c(1, 2), mean)
								}
								rm(tm_o)
							}
						} else if (identical(dynamic_uv_kind, "snap") && !is.null(delta_uv)) {
							snap_linked <- is.null(prior$snap_kappa)
							rho_uv <- .lame_sample_rho_uv_drift_only(
								U_cube,
								V_cube,
								sigma_uv,
								rho_uv,
								symmetric,
								delta_uv,
								delta_uv_v,
								prior_mean = pmean_uv,
								prior_sd = psd_uv,
								reset_stationary = snap_linked,
								stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2)
							)
							sigma_uv <- .lame_sample_sigma_uv_drift_only(
								U_cube,
								V_cube,
								rho_uv,
								sigma_uv,
								symmetric,
								delta_uv,
								delta_uv_v,
								prior$sigma_uv_shape,
								prior$sigma_uv_scale,
								sigma_max = prior$sigma_uv_max %||% Inf,
								stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2),
								reset_stationary = snap_linked
							)
						} else if (identical(dynamic_uv_kind, "t") && !is.null(lambda_u)) {
							# persistence and stationary variance as a block, the transitions weighted (R/uv_rho_w_step.R)
							uvh <- .lame_uv_rho_w_step(
								U_cube,
								V_cube,
								rho_uv,
								sigma_uv,
								pmean_uv,
								psd_uv,
								prior$sigma_uv_shape,
								prior$sigma_uv_scale,
								sigma_max = prior$sigma_uv_max %||% Inf,
								stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2),
								symmetric = symmetric,
								lambda_u = lambda_u,
								lambda_v = lambda_v
							)
							rho_uv <- uvh$rho
							sigma_uv <- uvh$sigma
						} else {
							# persistence and stationary variance as a block (R/uv_rho_w_step.R)
							uvh <- .lame_uv_rho_w_step(
								U_cube,
								V_cube,
								rho_uv,
								sigma_uv,
								pmean_uv,
								psd_uv,
								prior$sigma_uv_shape,
								prior$sigma_uv_scale,
								sigma_max = prior$sigma_uv_max %||% Inf,
								stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2),
								symmetric = symmetric
							)
							rho_uv <- uvh$rho
							sigma_uv <- uvh$sigma
						}
					}
				}
			} else {
				# static UV update
				if (bip) {
					if (RA > 0 && RB > 0) {
						# drop = FALSE keeps U[, , 1] a matrix when RA = 1
						U_2d <- if (length(dim(U)) == 3)
							U[, , 1, drop = FALSE] else U
						V_2d <- if (length(dim(V)) == 3)
							V[, , 1, drop = FALSE] else V
						if (length(dim(U_2d)) == 3) dim(U_2d) <- dim(U_2d)[1:2]
						if (length(dim(V_2d)) == 3) dim(V_2d) <- dim(V_2d)[1:2]
						# exact row draws given the other factor and G, over every period (R/bipartite_uv_gibbs.R)
						U_2d[] <- .lame_bip_factor_draw(E, V_2d %*% t(G), s2)
						V_2d[] <- .lame_bip_factor_draw(aperm(E, c(2L, 1L, 3L)), U_2d %*% G, s2)
						# back to 3d for get_EZ_bip_cpp
						if (length(dim(U)) == 3) {
							for (t in 1:N) { U[, , t] <- U_2d; V[, , t] <- V_2d }
						} else {
							U <- U_2d; V <- V_2d
						}
					}
				} else if (symmetric) {
					EA <- apply(E, c(1, 2), mean) ; EA <- .5 * (EA + t(EA))
					UV <- try(
						rUV_sym_fc_cpp(
							EA,
							U,
							V,
							s2 / dim(E)[3],
							shrink,
							sym_loop_ids[[s]] - 1,
							weight_prior_var
						),
						silent = TRUE
					)
					if (inherits(UV, 'try-error')) { UV <- list(U = U, V = V) ; tryErrorChecks$UV <- tryErrorChecks$UV + 1 }
					U <- UV$U ; V <- UV[["V"]]
				} else if (isTRUE(dynamic_rho)) {
					UV <- rUV_rep_fc(E, U, V, rho, s2, shrink, uv_prior_static$df, uv_prior_static$psi)
					U <- UV$U; V <- UV$V
				} else if (!symmetric) {
					rho_eff <- .lame_rho_scalar(rho)
					UV <- try(
						rUV_rep_fc_cpp(
							E,
							U,
							V,
							rho_eff,
							s2,
							diag(2), # unused whitening slot; rho and s2 set the law
							maxmargin = 1e-6,
							shrink,
							asym_loop_ids[[s]] - 1,
							uv_prior_static$df,
							uv_prior_static$psi
						),
						silent = TRUE
					)
					if (inherits(UV, 'try-error')) {
						UV <- try(
							rUV_rep_fc(E, U, V, rho_eff, s2, shrink, uv_prior_static$df, uv_prior_static$psi),
							silent = TRUE
						)
						if (inherits(UV, 'try-error')) { UV <- list(U = U, V = V) ; tryErrorChecks$UV <- tryErrorChecks$UV + 1 }
					}
					U <- UV$U ; V <- UV[["V"]]
				}
			}
			# scale split of the static factors of a directed network, the covariance integrated out (R/uv_scale_split.R)
			if (uv_split_on) {
				uv_sp <- .lame_uv_split_draw(U, V, uv_prior_static$df, uv_prior_static$psi, hier = shrink)
				U <- uv_sp$U; V <- uv_sp$V
				uv_split_stats <- uv_split_stats + c(1, uv_sp$c != 1)
				rm(uv_sp)
			}
			# update G for bipartite (interaction between the row and column latent spaces)
			if (bip && RA > 0 && RB > 0 && !is.null(G)) {
				U_3d <- if (.dyn_uv_has_rank) U_cube else if (length(dim(U)) == 3) U else {
					U_tmp <- array(0, dim = c(nA, RA, N))
					for (t in 1:N) U_tmp[, , t] <- U
					U_tmp
				}
				V_3d <- if (.dyn_uv_has_rank) V_cube else if (length(dim(V)) == 3) V else {
					V_tmp <- array(0, dim = c(nB, RB, N))
					for (t in 1:N) V_tmp[, , t] <- V
					V_tmp
				}
				G_prev <- G
				if (is.null(G_prev) || length(G_prev) == 0L ||
						is.null(dim(G_prev)) ||
						!all(dim(G_prev) == c(RA, RB)) || !all(is.finite(G_prev))) {
					G_prev <- if (RA > 0 && RB > 0) {
						diag(min(RA, RB), RA, RB)
					} else {
						matrix(0, max(1L, RA), max(1L, RB))
					}
				}
				G_new <- tryCatch(
					sample_G_bip_cpp(E, U_3d, V_3d, lambdaG = 1 / weight_prior_var, s2 = rep(s2, N)),
					error = function(e) {
						.lame_state$last_G_error <- conditionMessage(e)
						tryErrorChecks$G <<- tryErrorChecks$G + 1L
						G_prev
					}
				)
				# truncated prior on the interaction scale: once a binary likelihood saturates |G| is unidentified upward,
				# so a draw beyond the cap (prior$G_max) keeps the current G, an exact MH step for the truncated target
				g_cap <- prior$G_max %||% g_rail_default
				if (!is.null(G_new) && all(dim(G_new) == c(RA, RB)) &&
						all(is.finite(G_new)) && max(abs(G_new)) <= g_cap) {
					G <- G_new
				} else {
					# count draws beyond the cap apart from numerical failures
					if (!is.null(G_new) && all(is.finite(G_new)) &&
							max(abs(G_new)) > g_cap) {
						g_rail_count <- g_rail_count + 1L
					}
					G <- G_prev
					tryErrorChecks$G <- tryErrorChecks$G + 1L
				}
			}
		}

		# exact draw along the scale ridge of the static positions and G (U D1, V D2, D1^-1 G D2^-1 give the same product)
		if (uv_g_scale_bip && !dynamic_uv && !is.null(G) && all(dim(G) == c(RA, RB))) {
			rg <- .lame_bip_uv_g_scale_step(
				matrix(if (length(dim(U)) == 3L) U[, , 1L] else U, nA, RA),
				matrix(if (length(dim(V)) == 3L) V[, , 1L] else V, nB, RB),
				G,
				prior$G_max %||% g_rail_default,
				g_var = weight_prior_var
			)
			U <- array(rg$U, dim = dim(U))
			V <- array(rg$V, dim = dim(V))
			G <- rg$G
			uv_g_scale_stats <- uv_g_scale_stats + c(1, rg$drawn)
			rm(rg)
		}
		# rescale the factor paths and innovation scale with a compensating interaction; the likelihood is unchanged
		# and the prior ratio includes the Jacobian
		if (bip_scale_enabled && s > .5 * burn) {
			bip_scale_att <- bip_scale_att + 1L
			bs <- tryCatch(
				.lame_bipartite_scale_step(
					U_cube,
					V_cube,
					G,
					sigma_uv,
					rho_uv,
					prior$sigma_uv_shape,
					prior$sigma_uv_scale,
					max_abs = uv_max_abs_iter,
					max_G = prior$G_max %||% g_rail_default,
					sigma_max = prior$sigma_uv_max %||% 10,
					stat_var_cap = (prior$uv_var_cap %||% Inf) * .lame_cap_unit(family, s2),
					g_var = weight_prior_var
				),
				error = function(e) NULL
			)
			if (is.null(bs)) {
				bip_scale_numerical <- bip_scale_numerical + 1L
			} else {
				bip_scale_acc <- bip_scale_acc + bs$accepts
				bip_scale_bounds <- bip_scale_bounds + bs$bound_rejections
				bip_scale_numerical <- bip_scale_numerical + bs$numerical_rejections
				U_cube <- bs$U; V_cube <- bs$V; G <- bs$G; sigma_uv <- bs$sigma
				U_3d <- U_cube; V_3d <- V_cube
				U <- apply(U_cube, c(1L, 2L), mean)
				V <- apply(V_cube, c(1L, 2L), mean)
			}
		}
		# actor-specific slopes given the current complete base model
		if (per_actor_active) {
			Z <- Z + per_actor_offset
			n_a_pa <- dim(Y)[1L]; n_b_pa <- dim(Y)[2L]
			pa_beta <- if (beta_dyn$any) .lame_beta_full_path() else
				matrix(rep(beta, each = N), N, length(beta))
			pa_a <- if (dynamic_ab) a_mat else matrix(a, n_a_pa, N)
			pa_b <- if (dynamic_ab) b_mat else matrix(b, n_b_pa, N)
			pa_base <- .lame_per_actor_base_predictor(
				Xlist,
				pa_beta,
				pa_a,
				pa_b,
				if (.dyn_uv_has_rank) U_cube else U,
				if (.dyn_uv_has_rank) V_cube else V,
				if (bip) G else NULL,
				rank_active = R_eff > 0
			)
			if (!bip) Z <- .lame_refresh_diagonal(
				Z,
				.lame_diagonal(pa_base + per_actor_offset),
				rho,
				s2
			)
			n_actors_pa <- nrow(theta_actor)
			H_mat <- h_mat <- matrix(0, n_actors_pa, N)
			pa_precision <- vector("list", N); pa_coupled <- FALSE
			for (t_pa in seq_len(N)) {
				suff <- .lame_per_actor_suffstats(
					Z[, , t_pa] - pa_base[, , t_pa],
					per_actor_X[[t_pa]],
					dynamic_beta_per_actor,
					s2,
					rho = if (bip) 0 else .lame_rho_at(rho, t_pa),
					bipartite = bip
				)
				H_mat[, t_pa] <- suff$H; h_mat[, t_pa] <- suff$h
				pa_precision[[t_pa]] <- if (is.null(suff$precision)) diag(suff$H, n_actors_pa) else suff$precision
				pa_coupled <- pa_coupled || !is.null(suff$precision)
			}
			if (!pa_coupled) pa_precision <- NULL
			if (identical(per_actor_identifiability, "exact_center")) {
				theta_actor <- .sweep_per_actor_exact(
					H_mat,
					h_mat,
					rho_actor,
					sigma_actor_sq,
					precision_list = pa_precision
				)
			} else {
				Q_AR1 <- .ar1_path_precision(N, rho_actor, sigma_actor_sq)
				theta_actor <- .sweep_per_actor(
					theta_actor,
					H_mat,
					h_mat,
					Q_AR1,
					precision_list = pa_precision
				)
			}
			# interweaving step: (sigma, rho) and the paths move together with the standardised innovations held
			# fixed, before the hyperparameters are drawn given the paths
			if (pa_hyper_nc && nrow(theta_actor) >= 2L) {
				pa_nc <- .lame_pa_nc_step(
					theta_actor,
					H_mat,
					h_mat,
					pa_precision,
					rho_actor,
					sigma_actor_sq,
					prior$sigma_actor_shape,
					prior$sigma_actor_scale,
					prior$rho_actor_shape1,
					prior$rho_actor_shape2
				)
				theta_actor <- pa_nc$theta
				sigma_actor_sq <- pa_nc$sigma2
				rho_actor <- pa_nc$rho
				pa_nc_stats <- pa_nc_stats + c(1, pa_nc$attempts, pa_nc$accepts)
				rm(pa_nc)
			}
			sigma_actor_sq <- .sample_sigma_actor(
				theta_actor,
				rho_actor,
				prior$sigma_actor_shape,
				prior$sigma_actor_scale
			)
			rho_actor <- .sample_rho_actor(
				theta_actor,
				sigma_actor_sq,
				rho_actor,
				prior$rho_actor_shape1,
				prior$rho_actor_shape2
			)
			per_actor_offset <- .lame_per_actor_offset(theta_actor, per_actor_X, dynamic_beta_per_actor)
			EZ <- pa_base + per_actor_offset
		}

		# scale of the rank-likelihood latent state, whole and per period
		if (osm_on) {
			osm_bp <- if (beta_dyn$any) .lame_beta_full_path() else matrix(rep(beta, each = N), N, length(beta))
			osm_ap <- if (dynamic_ab) a_mat else matrix(a, dim(Y)[1L], N)
			osm_bq <- if (dynamic_ab) b_mat else matrix(b, dim(Y)[2L], N)
			osm_E <- .lame_saved_predictor(
				Xlist,
				osm_bp,
				osm_ap,
				osm_bq,
				if (.dyn_uv_has_rank) U_cube else U,
				if (.dyn_uv_has_rank) V_cube else V,
				G = if (bip) G else NULL,
				rank_active = R_eff > 0,
				offset = if (per_actor_active) per_actor_offset else NULL,
				symmetric = symmetric
			)
			# the part of the predictor the scale does not carry (positions, dynamic effects, actor-specific slopes)
			osm_Es <- if (osm_has_m) .lame_saved_predictor(
				Xlist,
				osm_bp,
				if (dynamic_ab) 0 * osm_ap else osm_ap,
				if (dynamic_ab) 0 * osm_bq else osm_bq,
				NULL,
				NULL,
				G = NULL,
				rank_active = FALSE,
				offset = NULL,
				symmetric = symmetric
			) else NULL
			osm_pr <- .lame_osm_prior_lame(
				beta_dyn,
				beta,
				if (beta_dyn$any) beta_dyn_path else NULL,
				if (beta_dyn$any) beta_static else NULL,
				osm_P0,
				g,
				s2,
				a,
				b,
				Sab,
				if (symmetric) isTRUE(nvar) else isTRUE(rvar),
				isTRUE(cvar),
				symmetric,
				bip,
				dynamic_ab,
				dynamic_beta_kind,
				rho_beta_by_group[group_id],
				sigma_beta_by_group[group_id],
				diag(Lambda_full),
				beta0_mean,
				diag(beta0_cov),
				osm_tpos,
				prior$matern32_length_scale %||% NULL,
				N
			)
			osm_o <- .lame_osm_visit(
				Z,
				osm_E,
				if (osm_has_m) osm_E - osm_Es else NULL,
				osm_mode,
				vapply(seq_len(N), function(t_) .lame_rho_at(rho, t_), 0),
				s2,
				A_pri = if (is.null(osm_pr)) 0 else osm_pr$A,
				B_pri = if (is.null(osm_pr)) 0 else osm_pr$B,
				d_th = if (is.null(osm_pr)) 0 else osm_pr$d - scale_cen - (if (sym_tie) .lame_sym_tie_dims(
					rb,
					cb,
					beta_dyn$dynamic_idx,
					beta_dyn$static_idx,
					N,
					isTRUE(beta_dyn$any)
				) else 0L) + (if (use_explicit_cutpoints) length(alpha) - 1L else 0L),
				affine = !use_explicit_cutpoints,
				orbit = !is.null(osm_pr)
			)
			# explicit cutpoints scale with the state (the first is the anchor at zero), which must stay inside its intervals
			if (use_explicit_cutpoints && osm_o$c != 1) {
				osm_alpha <- c(alpha[1L], osm_o$c * alpha[-1L])
				if (.lame_osm_ex_ok(osm_o$Z, Y_int, osm_alpha)) alpha <- osm_alpha else osm_o$c <- 1
			}
			if (osm_o$c != 1 || !use_explicit_cutpoints) Z <- osm_o$Z
			osm_stats <- osm_stats + c(1, osm_o$c != 1)
			if (osm_o$c != 1) {
				if (beta_dyn$any) {
					beta_dyn_path <- osm_o$c * beta_dyn_path
					beta_static <- osm_o$c * beta_static
					beta[beta_dyn$dynamic_idx] <- colMeans(beta_dyn_path)
					beta[beta_dyn$static_idx] <- beta_static
				} else beta <- osm_o$c * beta
				if (!dynamic_ab) {
					a <- osm_o$c * a
					b <- osm_o$c * b
				}
			}
			rm(osm_bp, osm_ap, osm_bq, osm_E, osm_Es, osm_pr, osm_o)
		}

		# scale of the linear predictor of a binary fit, with the latent array integrated out
		if (bso_on && (s %% bso_every == 0L)) {
			bso_bp <- if (beta_dyn$any) .lame_beta_full_path() else matrix(rep(beta, each = N), N, length(beta))
			bso_ap <- if (dynamic_ab) a_mat else matrix(a, dim(Y)[1L], N)
			bso_bq <- if (dynamic_ab) b_mat else matrix(b, dim(Y)[2L], N)
			bso_E <- .lame_saved_predictor(
				Xlist,
				bso_bp,
				bso_ap,
				bso_bq,
				if (.dyn_uv_has_rank) U_cube else U,
				if (.dyn_uv_has_rank) V_cube else V,
				G = if (bip) G else NULL,
				rank_active = R_eff > 0,
				offset = if (per_actor_active) per_actor_offset else NULL,
				symmetric = symmetric
			)
			bso_Es <- if (osm_has_m) .lame_saved_predictor(
				Xlist,
				bso_bp,
				if (dynamic_ab) 0 * bso_ap else bso_ap,
				if (dynamic_ab) 0 * bso_bq else bso_bq,
				NULL,
				NULL,
				G = NULL,
				rank_active = FALSE,
				offset = NULL,
				symmetric = symmetric
			) else NULL
			bso_pr <- .lame_osm_prior_lame(
				beta_dyn,
				beta,
				if (beta_dyn$any) beta_dyn_path else NULL,
				if (beta_dyn$any) beta_static else NULL,
				bso_P0,
				g,
				s2,
				a,
				b,
				Sab,
				if (symmetric) isTRUE(nvar) else isTRUE(rvar),
				isTRUE(cvar),
				symmetric,
				bip,
				dynamic_ab,
				dynamic_beta_kind,
				rho_beta_by_group[group_id],
				sigma_beta_by_group[group_id],
				diag(Lambda_full),
				beta0_mean,
				diag(beta0_cov),
				bso_tpos,
				prior$matern32_length_scale %||% NULL,
				N
			)
			if (!is.null(bso_pr)) {
				bso_o <- .lame_bso_visit(
					Z,
					bso_E,
					if (osm_has_m) bso_E - bso_Es else NULL,
					bso_cx,
					if (bip) 0 else .lame_rho_scalar(rho),
					A_pri = bso_pr$A,
					B_pri = bso_pr$B,
					d_th = bso_pr$d - scale_cen - (if (sym_tie) .lame_sym_tie_dims(
						rb,
						cb,
						beta_dyn$dynamic_idx,
						beta_dyn$static_idx,
						N,
						isTRUE(beta_dyn$any)
					) else 0L)
				)
				bso_stats <- bso_stats + c(1, bso_o$evals, bso_o$c != 1)
				if (bso_o$c != 1) {
					Z <- bso_o$Z
					if (beta_dyn$any) {
						beta_dyn_path <- bso_o$c * beta_dyn_path
						beta_static <- bso_o$c * beta_static
						beta[beta_dyn$dynamic_idx] <- colMeans(beta_dyn_path)
						beta[beta_dyn$static_idx] <- beta_static
					} else beta <- bso_o$c * beta
					if (!dynamic_ab) {
						a <- bso_o$c * a
						b <- bso_o$c * b
					}
				}
				rm(bso_o)
			}
			rm(bso_bp, bso_ap, bso_bq, bso_E, bso_Es, bso_pr)
		}

		# burn-in countdown
		if (burn != 0 && s <= burn && verbose) {
			cli::cli_progress_update()
		}

		# store parameter values and monitor the chain
		if (s == burn + 1 && verbose) {
			# with burn = 0 (the usual resume path) there is no burn-in bar to close; the sampling bar still opens
			if (burn != 0) {
				cli::cli_progress_done()
				cli::cli_alert_success("Burn-in period complete")
			}
			cli::cli_progress_bar("Sampling", total = nscan, .envir = environment())
		}
		if (s %% odens == 0 & s > burn) {

			# per-period full beta path when dynamic_beta is on, for 3-d BETA storage and summaries
			beta_path_now_store <- if (beta_dyn$any) {
				bp <- matrix(0, N, length(beta))
				for (t in seq_len(N)) {
					bp[t, ] <- assemble_beta_path(
						beta_dyn_path[t, ],
						beta_static,
						beta_dyn$dynamic_idx,
						beta_dyn$static_idx,
						length(beta)
					)
				}
				bp
			} else NULL
			# saved scores and predictive replicates use the same current state as the stored parameters
			EZ <- .lame_saved_predictor(
				Xlist,
				if (beta_dyn$any) beta_path_now_store else matrix(rep(beta, each = N), N, length(beta)),
				if (dynamic_ab) a_mat else matrix(a, dim(Y)[1L], N),
				if (dynamic_ab) b_mat else matrix(b, dim(Y)[2L], N),
				if (.dyn_uv_has_rank) U_cube else U,
				if (.dyn_uv_has_rank) V_cube else V,
				G = if (bip) G else NULL,
				rank_active = R_eff > 0,
				offset = if (per_actor_active) .lame_per_actor_offset(
					theta_actor,
					per_actor_X,
					dynamic_beta_per_actor
				) else NULL,
				symmetric = symmetric
			)
			# reported parameterisation: with static positions the draw is stored in identified form, positions
			# centred over actors with the removed row / column parts added to the effects and the constant to the
			# intercept, leaving the predictor unchanged. dynamic positions and the sampler state are left as they are
			beta_rep <- beta
			beta_path_rep <- beta_path_now_store
			a_rep <- a
			b_rep <- b
			a_mat_rep <- if (dynamic_ab) a_mat else NULL
			b_mat_rep <- if (dynamic_ab) b_mat else NULL
			a_lev_rep <- if (dynamic_ab) a_lev else NULL
			b_lev_rep <- if (dynamic_ab) b_lev else NULL
			U_rep <- U
			V_rep <- V
			if (isTRUE(rvar) && isTRUE(cvar) && !.dyn_uv_has_rank &&
					R_eff > 0 && !is.null(U) && !is.null(V) &&
					NCOL(U) > 0L && NCOL(V) > 0L) {
				U_st <- if (length(dim(U)) == 3L) matrix(U[, , 1L], nrow = dim(U)[1L]) else as.matrix(U)
				V_st <- if (length(dim(V)) == 3L) matrix(V[, , 1L], nrow = dim(V)[1L]) else as.matrix(V)
				G_st <- if (bip && !is.null(G)) G else diag(1, ncol(U_st))
				cen_rep <- if (all(is.finite(U_st)) && all(is.finite(V_st)) && all(is.finite(G_st)) &&
						ncol(U_st) == nrow(G_st) && ncol(V_st) == ncol(G_st) &&
						nrow(U_st) == length(a_rep) && nrow(V_st) == length(b_rep)) {
					.lame_center_uv_period(U_st, V_st, G_st)
				} else NULL
				if (!is.null(cen_rep)) {
					a_row <- cen_rep$a_row
					b_col <- cen_rep$b_col
					cst <- cen_rep$const
					if (symmetric) {
						# one effect per actor: under V = U L the row and column halves are the same vector; averaging keeps
						# that under rounding
						a_row <- (a_row + b_col) / 2
						b_col <- a_row
					}
					if (isTRUE(intercept) && length(beta_rep) > 0L && is.finite(beta_rep[1L])) {
						beta_rep[1L] <- beta_rep[1L] + cst
						if (!is.null(beta_path_rep)) beta_path_rep[, 1L] <- beta_path_rep[, 1L] + cst
					} else if (symmetric) {
						# a shared actor effect appears twice in each predictor
						a_row <- a_row + cst / 2
						b_col <- b_col + cst / 2
					} else {
						a_row <- a_row + cst
					}
					a_rep <- a_rep + a_row
					b_rep <- b_rep + b_col
					if (dynamic_ab) {
						a_mat_rep <- a_mat_rep + a_row
						b_mat_rep <- b_mat_rep + b_col
						a_lev_rep <- a_lev_rep + a_row
						b_lev_rep <- b_lev_rep + b_col
					}
					if (length(dim(U)) == 3L) {
						U_rep <- U
						for (tt in seq_len(dim(U)[3L])) U_rep[, , tt] <- cen_rep$U
					} else U_rep <- cen_rep$U
					if (length(dim(V)) == 3L) {
						V_rep <- V
						for (tt in seq_len(dim(V)[3L])) V_rep[, , tt] <- cen_rep[["V"]]
					} else V_rep <- cen_rep[["V"]]
				}
			}
			# the summaries the chain policy compares across the run
			.chain_acc <- .lame_chain_acc_update(
				.chain_acc,
				iter,
				nscan / odens,
				if (dynamic_ab) a_mat_rep else a_rep,
				if (dynamic_ab) b_mat_rep else b_rep,
				if (R_eff > 0) { if (.dyn_uv_has_rank) U_cube else U_rep } else NULL,
				if (R_eff > 0) { if (.dyn_uv_has_rank) V_cube else V_rep } else NULL,
				if (bip) G else NULL,
				symmetric = symmetric,
				actors = list(dimnames(Y)[[1L]], dimnames(Y)[[2L]])
			)
			FC_STATE$a[, iter] <- if (dynamic_ab) a_mat_rep[, N] else a_rep
			FC_STATE$b[, iter] <- if (dynamic_ab) b_mat_rep[, N] else b_rep
			if (dynamic_ab) {
				FC_STATE$a_level[, iter] <- a_lev_rep
				FC_STATE$b_level[, iter] <- b_lev_rep
			}
			if (!is.null(FC_STATE$U)) {
				FC_STATE$U[, , iter] <- if (.dyn_uv_has_rank) U_cube[, , N] else
				if (length(dim(U_rep)) == 3L) U_rep[, , 1L] else U_rep
				FC_STATE$V[, , iter] <- if (.dyn_uv_has_rank) V_cube[, , N] else
				if (length(dim(V_rep)) == 3L) V_rep[, , 1L] else V_rep
				if (bip) FC_STATE$G[, , iter] <- G
			}
			# store beta and VC, symmetric case
			if (!is.null(A_SAMPLES) && iter %% thin_ab == 0L) {
				ab_idx <- iter %/% thin_ab
				if (dynamic_ab) {
					A_SAMPLES[, , ab_idx] <- a_mat_rep
					B_SAMPLES[, , ab_idx] <- b_mat_rep
					A_LEVEL_SAMPLES[, ab_idx] <- a_lev_rep
					B_LEVEL_SAMPLES[, ab_idx] <- b_lev_rep
				} else {
					A_SAMPLES[, ab_idx] <- a_rep
					B_SAMPLES[, ab_idx] <- b_rep
				}
			}
			if (symmetric) {
				br <- beta_rep[rb] ; bc <- beta_rep[cb] ; bn <- (br + bc) / 2
				# build the exclude index explicitly, since beta[-integer(0)] would be empty when there are no nodal covariates
				excl_idx <- c(if (intercept) 1L else integer(0), rb, cb)
				bint <- if (intercept) beta_rep[1L] else numeric(0)
				bdyad <- if (length(excl_idx) > 0L) beta_rep[-excl_idx] else beta_rep
				sbeta <- c(bint, bn, bdyad)
				if (beta_dyn$any) {
					# 3-d BETA storage: the same rb/cb merge applied to each period's beta
					for (t in seq_len(N)) {
						b_t <- beta_path_rep[t, ]
						brt <- b_t[rb]; bct <- b_t[cb]; bnt <- (brt + bct) / 2
						excl_t <- c(if (intercept) 1L else integer(0), rb, cb)
						bint_t <- if (intercept) b_t[1L] else numeric(0)
						bdyad_t <- if (length(excl_t) > 0L) b_t[-excl_t] else b_t
						st <- c(bint_t, bnt, bdyad_t)
						BETA[iter, , t] <- st
					}
				} else {
					BETA[iter, ] <- sbeta
				}
				VC[iter, ] <- c(Sab[1, 1], s2)
			}
			# store beta and VC, asymmetric case
			if (!symmetric) {
				if (beta_dyn$any) {
					BETA[iter, , ] <- t(beta_path_rep) # p x t slot
				} else {
					BETA[iter, ] <- beta_rep
				}
				rho_store <- if (isTRUE(dynamic_rho)) mean(rho, na.rm = TRUE) else rho
				VC[iter, ] <- c(Sab[upper.tri(Sab, diag = TRUE)], rho_store, s2)
			}
			# store dynamic parameters
			if (dynamic_ab && !is.null(rho_ab)) {
				RHO_AB[iter] <- rho_ab
				if (!is.null(sigma_ab)) SIGMA_AB[iter] <- sigma_ab
			}
			if (dynamic_uv && !is.null(rho_uv)) {
				RHO_UV[iter] <- rho_uv
				if (!is.null(sigma_uv)) SIGMA_UV[iter] <- sigma_uv
			}
			# optional draw-level storage of the full dynamic cubes (skipped on non-dynamic branches)
			if (!is.null(U_DRAWS) &&
					exists("U_cube", inherits = FALSE) && !is.null(U_cube) &&
					length(dim(U_cube)) == 3L &&
					exists("V_cube", inherits = FALSE) && !is.null(V_cube) &&
					length(dim(V_cube)) == 3L) {
				U_DRAWS[, , , iter] <- U_cube
				V_DRAWS[, , , iter] <- V_cube
			}
			# static latent-position draw storage; a static bipartite U / V held as identical 3-d slices stores slice 1
			if (!is.null(U_SAMPLES)) {
				U_SAMPLES[, , iter] <- if (length(dim(U_rep)) == 3L)
					matrix(U_rep[, , 1L], nrow = dim(U_rep)[1L]) else U_rep
				if (!is.null(V_SAMPLES)) {
					V_SAMPLES[, , iter] <- if (length(dim(V_rep)) == 3L)
						matrix(V_rep[, , 1L], nrow = dim(V_rep)[1L]) else V_rep
				}
			}
			# bipartite G draws accompany the latent draws so U G V' can be rebuilt per draw (as in ame())
			if (!is.null(G_SAMPLES) && !is.null(G) &&
					identical(dim(G), dim(G_SAMPLES)[1:2])) {
				G_SAMPLES[, , iter] <- G
			}
			if (isTRUE(dynamic_rho) && !is.null(RHO)) {
				RHO[iter, ] <- .lame_expand_rho(rho, N)
			}
			if (dynamic_uv && identical(dynamic_uv_kind, "snap") && !is.null(pi_snap)) {
				PI_SNAP[iter] <- pi_snap
			}
			# store dynamic_beta hyperparameters
			if (beta_dyn$any) {
				RHO_BETA[iter, ] <- rho_beta_by_group
				SIGMA_BETA[iter, ] <- sigma_beta_by_group
			}
			# store explicit ordinal cutpoints alpha_2..alpha_{K-1} (none for K = 2)
			if (use_explicit_cutpoints && !is.null(ALPHA)) {
				ALPHA[iter, ] <- alpha[-1L]
			}
			# accumulate the static bipartite interaction matrix
			if (!is.null(G_sum) && !is.null(G) && all(is.finite(G))) {
				G_sum <- G_sum + G
				G_n <- G_n + 1L
			}
			# pointwise log-likelihood at this draw, per observed dyad, from the refreshed EZ
			if (!is.null(LOG_LIK) || !is.null(LOG_LIK_CONNS)) {
				n_obs_iter <- if (!is.null(LOG_LIK)) ncol(LOG_LIK) else LOG_LIK_META$n_obs_total
				ll_row <- numeric(n_obs_iter)
				# observed y and matching EZ in the same order
				y_obs <- Y[obs_idx]
				ez_obs <- EZ[obs_idx]
				# explicit ordinal and poisson observations use their fitted marginal laws
				if (use_explicit_cutpoints && !identical(log_lik_method, "augmented")) {
					ll_row <- .lame_explicit_ordinal_loglik(y_obs, ez_obs, alpha, ord_lvls)
				} else if (family == "poisson" && !identical(log_lik_method, "augmented")) {
					exposure_obs <- if (poisson_exposure_active) period_exposure[obs_idx[, 3L]] else rep(1, length(y_obs))
					positive_exposure <- exposure_obs > 0
					ll_row <- numeric(length(y_obs))
					ll_row[positive_exposure] <- .lame_pois_marginal_loglik(
						round(y_obs[positive_exposure]),
						ez_obs[positive_exposure] + log(exposure_obs[positive_exposure]),
						s2
					)
					ll_row[!positive_exposure & y_obs != 0] <- -Inf
				} else if (identical(log_lik_method, "observed_ghk")) {
					# GHK + Gauss-Hermite + dimension-capped pipeline; halton_shift comes from fit$log_lik_meta when present
					ll_row <- .pointwise_loglik_observed_ghk_rank(
						y_obs = y_obs,
						ez_obs = ez_obs,
						z_obs = Z[obs_idx],
						obs_idx = obs_idx,
						family = family,
						s2 = s2,
						dyad_rho = if (is.null(rho)) 0 else .lame_rho_scalar(rho),
						n_mc = .lame_ghk_n_mc %||% 64L,
						n_mc_fixed = !is.null(.lame_ghk_n_mc),
						halton_shift = .lame_log_lik_halton_shift %||% 0
					)
				} else if (identical(log_lik_method, "augmented")) {
					# augmented-data contribution: Gaussian density of the latent z around its mean
					ll_row <- dnorm(
						Z[obs_idx],
						mean = ez_obs,
						sd = sqrt(s2),
						log = TRUE
					)
				} else if (family == "normal") {
					ll_row <- dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
				} else if (family == "binary" || family == "cbin") {
					# signed log tails keep extreme probit scores
					ll_row <- .binary_pointwise_loglik(y_obs, ez_obs)
				} else if (family == "ordinal") {
					# closed-form cumulative-probit log-likelihood with empirical cutpoints
					ll_row <- .ordinal_pointwise_loglik(y_obs, ez_obs)
				} else {
					# frn ranked likelihoods use the augmented-data normal contribution unless observed_ghk is requested
					if (iter == 1L) {
						cli::cli_warn(c(
							"{.arg save_log_lik} = TRUE: family {.val {family}} pointwise log-lik uses the augmented-data normal approximation.",
							"i" = "{.fn loo::loo} and {.fn waic} refuse these scores (error class {.code lame_loo_factorization_unsupported}); read them with {.fn read_log_lik}."
						))
					}
					ll_row <- dnorm(Z[obs_idx], mean = ez_obs, sd = sqrt(s2), log = TRUE)
				}
				if (!is.null(LOG_LIK)) {
					LOG_LIK[iter, ] <- ll_row
				} else {
					# chunked: append this iteration's values to each chunk file, row-major within its column range
					for (cc in seq_len(LOG_LIK_META$n_chunks)) {
						writeBin(
							ll_row[LOG_LIK_META$chunk_starts[cc]:LOG_LIK_META$chunk_ends[cc]],
							LOG_LIK_CONNS[[cc]],
							size = 8L
						)
					}
				}
			}
			# update posterior sums of random effects
			if (.dyn_uv_has_rank && !is.null(U_cube) && !is.null(V_cube)) {
				if (!.lame_uv_ok(U_cube, V_cube, Inf)) {
					U_cube[!is.finite(U_cube)] <- 0
					V_cube[!is.finite(V_cube)] <- 0
					tryErrorChecks$UV <- tryErrorChecks$UV + 1L
				}
				U_SUM <- U_SUM + U_cube
				V_SUM <- V_SUM + V_cube
				if (exists("U_SUM_SQ")) { U_SUM_SQ <- U_SUM_SQ + U_cube^2 }
				if (exists("V_SUM_SQ")) { V_SUM_SQ <- V_SUM_SQ + V_cube^2 }
				if (exists("L_dyn", inherits = FALSE) && !is.null(L_dyn)) {
					L_SUM <- if (is.null(L_SUM)) L_dyn else L_SUM + L_dyn
					L_N <- L_N + 1L
				}
				if (identical(dynamic_uv_kind, "snap") && !is.null(DELTA_SUM) && !is.null(delta_uv)) {
					du_store <- delta_uv
					du_store[!snap_eligible_u] <- NA_real_
					DELTA_SUM <- DELTA_SUM + ifelse(is.finite(du_store), du_store, 0)
					DELTA_COUNT <- DELTA_COUNT + ifelse(is.finite(du_store), 1L, 0L)
					if (!is.null(DELTA_SUM_V) && !is.null(delta_uv_v)) {
						dv_store <- delta_uv_v
						dv_store[!snap_eligible_v] <- NA_real_
						DELTA_SUM_V <- DELTA_SUM_V + ifelse(is.finite(dv_store), dv_store, 0)
						DELTA_COUNT_V <- DELTA_COUNT_V + ifelse(is.finite(dv_store), 1L, 0L)
					}
					if (!is.null(SNAP_DRAWS)) {
						SNAP_DRAWS[iter, , ] <- du_store
						if (!symmetric && !is.null(SNAP_DRAWS_V) && !is.null(delta_uv_v)) {
							SNAP_DRAWS_V[iter, , ] <- dv_store
						}
					}
				}
				# UV products for each time point
				for (t in 1:N) {
					U_t <- U_cube[, , t, drop = FALSE]
					V_t <- V_cube[, , t, drop = FALSE]
					if (bip) {
						U_t <- matrix(U_t, nrow = nA, ncol = RA)
						V_t <- matrix(V_t, nrow = nB, ncol = RB)
						UVPS[, , t] <- UVPS[, , t] + U_t %*% G %*% t(V_t)
					} else {
						U_t <- matrix(U_t, nrow = n, ncol = R)
						V_t <- matrix(V_t, nrow = n, ncol = R)
						UVPS[, , t] <- UVPS[, , t] + U_t %*% t(V_t)
					}
				}
			} else if (R_eff > 0) {
				if (bip) {
					# bipartite: 2d slices with the G interaction matrix
					U_2d <- if (length(dim(U_rep)) == 3) U_rep[, , 1] else U_rep
					V_2d <- if (length(dim(V_rep)) == 3) V_rep[, , 1] else V_rep
					UVPS <- UVPS + U_2d %*% G %*% t(V_2d)
				} else {
					UVPS <- UVPS + U_rep %*% t(V_rep)
				}
			}
			if (dynamic_ab) {
				# dynamic ab: accumulate time-specific effects
				if (iter == 1) {
					APS_dyn <- a_mat_rep
					BPS_dyn <- b_mat_rep
					# running sums of the actor levels behind the paths
					A_LEV_SUM <- a_lev_rep
					B_LEV_SUM <- b_lev_rep
					# Welford streaming mean and second moment for the per-period posterior sd (fit$a_dynamic_sd, fit$b_dynamic_sd)
					ab_dyn_n <- 0L
					a_dyn_mean <- matrix(0, nrow(a_mat), ncol(a_mat))
					a_dyn_m2 <- matrix(0, nrow(a_mat), ncol(a_mat))
					b_dyn_mean <- matrix(0, nrow(b_mat), ncol(b_mat))
					b_dyn_m2 <- matrix(0, nrow(b_mat), ncol(b_mat))
				} else {
					APS_dyn <- APS_dyn + a_mat_rep
					BPS_dyn <- BPS_dyn + b_mat_rep
					A_LEV_SUM <- A_LEV_SUM + a_lev_rep
					B_LEV_SUM <- B_LEV_SUM + b_lev_rep
				}
				# streaming variance update at every stored draw
				ab_dyn_n <- ab_dyn_n + 1L
				d_a <- a_mat_rep - a_dyn_mean
				a_dyn_mean <- a_dyn_mean + d_a / ab_dyn_n
				a_dyn_m2 <- a_dyn_m2 + d_a * (a_mat_rep - a_dyn_mean)
				d_b <- b_mat_rep - b_dyn_mean
				b_dyn_mean <- b_dyn_mean + d_b / ab_dyn_n
				b_dyn_m2 <- b_dyn_m2 + d_b * (b_mat_rep - b_dyn_mean)
				APS <- APS + rowMeans(a_mat_rep)
				BPS <- BPS + rowMeans(b_mat_rep)
			} else {
				APS <- APS + a_rep
				BPS <- BPS + b_rep
				APS2 <- APS2 + a_rep^2
				BPS2 <- BPS2 + b_rep^2
			}

			# simulate from the posterior predictive, using the complete current predictor
			if (!bip) {
				dimnames(EZ) <- dimnames(Y)
				Ys <- EZ * 0
				# cbin / frn: a hidden cell may hold a nomination in the replicate, as in the sampler; the diagonal and
				# absent actors stay out (.lame_nom_predictive_yo)
				if (family %in% c("cbin", "frn") && is.null(yo_pred)) {
					yo_pred <- lapply(1:N, function(tt) .lame_nom_predictive_yo(Y[, , tt], bip))
				}
				for (t in 1:N) {
					if (symmetric) { EZ[, , t] <- (EZ[, , t] + t(EZ[, , t])) / 2 }
					rho_t <- if (isTRUE(dynamic_rho)) .lame_rho_at(rho, t) else rho

					if (family == "binary") { Ys[, , t] <- simY_bin(EZ[, , t], rho_t) }
					if (family == "cbin") { Ys[, , t] <- 1 * (simY_frn(EZ[, , t], rho_t, odmax, YO = yo_pred[[t]]) > 0)}
					if (family == "frn") { Ys[, , t] <- simY_frn(EZ[, , t], rho_t, odmax, YO = yo_pred[[t]]) }
					if (family == "normal") { Ys[, , t] <- simY_nrm(EZ[, , t], rho_t, s2) }
					if (family == "ordinal") {
						Ys[, , t] <- if (use_explicit_cutpoints)
							.simY_ord_explicit(EZ[, , t], rho_t, alpha, ord_lvls, symmetric)
						else simY_ord(EZ[, , t], rho_t, Y[, , t])
					}
					if (family == "poisson") {
						# exposure-aware simulation when exposures are non-trivial
						if (poisson_exposure_active) {
							# with the latent overdispersion layer: y ~ Pois(exp(z)), z ~ N(EZ, s2 Sigma_rho)
							Ys[, , t] <- simY_pois_exposure(
								simZ(EZ[, , t], rho_t, s2),
								period_exposure[t]
							)
						} else {
							Ys[, , t] <- simY_pois(simZ(EZ[, , t], rho_t, s2))
						}
					}

					if (symmetric) {
						Yst <- Ys[, , t] ; Yst[lower.tri(Yst)] <- 0 ; Ys[, , t] <- Yst + t(Yst)
					}
				}
			} else {
				# bipartite posterior predictive simulation (no dyadic correlation)
				Ys <- array(0, dim = c(nA, nB, N))
				# cbin / frn: hidden cells may hold a nomination, as in the sampler; absent actors' cells may not
				if (family %in% c("cbin", "frn") && is.null(yo_pred)) {
					yo_pred <- lapply(1:N, function(tt) .lame_nom_predictive_yo(Y[, , tt], bip))
				}
				for (t in 1:N) {
					EZ_t <- EZ[, , t]
					if (family == "binary") {
						noise <- matrix(rnorm(nA * nB), nA, nB)
						Ys[, , t] <- 1 * (EZ_t + noise > 0)
					} else if (family == "normal") {
						Ys[, , t] <- EZ_t + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
					} else if (family == "poisson") {
						# overdispersed predictive as on the unipartite path: y ~ Pois(exposure * e^z), z ~ N(EZ, s2)
						Zs <- EZ_t + matrix(rnorm(nA * nB, 0, sqrt(max(s2, 0))), nA, nB)
						exposure_t <- if (poisson_exposure_active) period_exposure[t] else 1
						Ys[, , t] <- simY_pois_exposure(Zs, exposure_t, bipartite = TRUE)
					} else if (family == "ordinal") {
						# quantile-based discretization
						uY <- sort(unique(c(Y[, , t])))
						uY <- uY[!is.na(uY)]
						FY <- table(c(Y[, , t]))
						FY <- FY / sum(FY)
						FY <- cumsum(FY)
						ZS <- EZ_t + matrix(rnorm(nA * nB), nA, nB)
						qZ <- quantile(ZS, FY[-length(FY)], na.rm = TRUE)
						YS <- ZS * 0 + max(uY)
						for (w in 1:(length(uY) - 1)) {
							YS[(YS == max(uY)) & ZS <= qZ[w]] <- uY[w]
						}
						Ys[, , t] <- YS
					} else if (family == "cbin") {
						# censored binary: the predictive draw honours the per-row cap odmax, as in the unipartite branch
						ZS <- EZ_t + matrix(rnorm(nA * nB), nA, nB)
						ZS[!is.finite(ZS)] <- 0
						out_t <- is.na(yo_pred[[t]])
						if (any(out_t)) ZS[out_t] <- -Inf
						YS <- ZS * 0
						for (i in 1:nA) {
							rs <- rank(ZS[i, ]) - (nB - odmax[i])
							YS[i, ] <- 1 * ((rs > 0) & (ZS[i, ] > 0))
						}
						if (any(out_t)) YS[out_t] <- NA
						Ys[, , t] <- YS
					} else if (family == "frn") {
						# fixed rank nomination; non-finite ZS is zeroed so rank() gives no NA
						ZS <- EZ_t + matrix(rnorm(nA * nB), nA, nB)
						ZS[!is.finite(ZS)] <- 0
						out_t <- is.na(yo_pred[[t]])
						if (any(out_t)) ZS[out_t] <- -Inf
						YS <- ZS * 0
						for (i in 1:nA) {
							rs <- rank(ZS[i, ]) - (nB - odmax[i])
							YS[i, ] <- rs * (rs > 0) * (ZS[i, ] > 0)
							pos <- YS[i, ] > 0
							if (any(pos)) {
								YS[i, pos] <- match(YS[i, pos], sort(unique(YS[i, pos])))
							}
						}
						if (any(out_t)) YS[out_t] <- NA
						Ys[, , t] <- YS
					}
				}
			}

			# update posterior sum
			YPS <- YPS + Ys

			# posterior predictive GOF stats
			if (gof) {
				Ys[is.na(Y)] <- NA
				if (bip) {
					GOF[, , (iter) + 1] <- apply(Ys, 3, gof_stats_bipartite, warn_square = FALSE, family = family)
				} else {
					GOF[, , (iter) + 1] <- apply(Ys, 3, gof_stats, mode = "unipartite")
				}
			}

			# keep the time dimension for temporary fits; flatten it only for progress summaries
			if (verbose || (plot && iter %% 10 == 0)) {
				current_BETA <- if (beta_dyn$any) BETA[seq_len(iter), , , drop = FALSE] else
					BETA[seq_len(iter), , drop = FALSE]
			}
			# progress report
			if (verbose) {
				beta_progress <- if (beta_dyn$any) matrix(current_BETA, nrow = iter) else current_BETA
				beta_means <- round(colMeans(beta_progress), 2)
				vc_means <- round(apply(VC[1:iter, , drop = FALSE], 2, mean), 2)
				cli::cli_text("Iteration {.val {s}}: beta = [{.field {paste(beta_means, collapse=', ')}}], VC = [{.field {paste(vc_means, collapse=', ')}}]")
				if (have_coda & nrow(VC[1:iter, , drop = FALSE]) > 3 & length(beta) > 0) {
					eff_sizes <- round(coda::effectiveSize(beta_progress))
					cli::cli_text("  Effective sizes: [{.emph {paste(eff_sizes, collapse=', ')}}]")
				}
			}

			# periodic save
			if (periodic_save & s %in% save_points & !is.null(out_file)) {
				# start_vals for future runs
				start_vals <- .lame_snapshot_state()
				fit <- get_fit_object(
					APS = APS,
					BPS = BPS,
					UVPS = UVPS,
					YPS = YPS,
					BETA = BETA,
					VC = VC,
					GOF = GOF,
					Xlist = Xlist,
					actorByYr = actorByYr,
					colActorByYr = if (bip) colActorByYr else NULL,
					start_vals = start_vals,
					symmetric = symmetric,
					tryErrorChecks = tryErrorChecks,
					model.name = model.name,
					family = family,
					odmax = odmax,
					nA = if (bip) nA else NULL,
					nB = if (bip) nB else NULL,
					n_time = N,
					dynamic_uv = dynamic_uv,
					dynamic_ab = dynamic_ab,
					bip = bip,
					G = if (bip) G else NULL,
					dynamic_beta = beta_dyn$any,
					beta_dynamic_mask = beta_dyn_fit$mask,
					beta_dynamic_groups = beta_dyn_fit$groups,
					rho_beta = rho_beta_by_group,
					sigma_beta = sigma_beta_by_group,
					RHO_BETA = RHO_BETA,
					SIGMA_BETA = SIGMA_BETA,
					dynamic_rho = dynamic_rho,
					RHO = if (isTRUE(dynamic_rho)) RHO[seq_len(iter), , drop = FALSE] else NULL,
					period_exposure = period_exposure
				)
				save(fit, file = out_file) ; rm(list = c('fit', 'start_vals'))
			}

			# live plotting
			if (plot && iter %% 10 == 0) {
				if (iter == 10) {
					cli::cli_alert_info("Real-time plotting enabled - updating every 10 iterations")
				}
				# temporary fit object for plotting
				current_APS <- if (dynamic_ab) APS_dyn / iter else APS
				current_BPS <- if (dynamic_ab) BPS_dyn / iter else BPS
				current_UVPS <- UVPS / iter
				current_U <- if (.dyn_uv_has_rank) U_SUM / iter else U
				current_V <- if (.dyn_uv_has_rank) V_SUM / iter else V

				plot_result <- try({
					temp_fit <- get_fit_object(
						APS = current_APS,
						BPS = current_BPS,
						UVPS = current_UVPS,
						YPS = YPS,
						BETA = current_BETA,
						VC = VC[1:iter, , drop = FALSE],
						GOF = GOF,
						Xlist = Xlist,
						actorByYr = actorByYr,
						colActorByYr = if (bip) colActorByYr else NULL,
						start_vals = NULL,
						symmetric = symmetric,
						tryErrorChecks = tryErrorChecks,
						model.name = model.name,
						U = current_U,
						V = current_V,
						dynamic_uv = dynamic_uv,
						dynamic_ab = dynamic_ab,
						bip = bip,
						rho_ab = if (dynamic_ab) RHO_AB[1:iter] else NULL,
						rho_uv = if (dynamic_uv) RHO_UV[1:iter] else NULL,
						family = family,
						odmax = odmax,
						nA = if (bip) nA else NULL,
						nB = if (bip) nB else NULL,
						n_time = N,
						G = if (bip) G else NULL,
						dynamic_beta = beta_dyn$any,
						beta_dynamic_mask = beta_dyn_fit$mask,
						beta_dynamic_groups = beta_dyn_fit$groups,
						dynamic_rho = dynamic_rho,
						RHO = if (isTRUE(dynamic_rho)) RHO[1:iter, , drop = FALSE] else NULL,
						period_exposure = period_exposure
					)
					class(temp_fit) <- "lame"
					suppressMessages(plot(temp_fit, which = c(1, 2), pages = "single"))
				}, silent = TRUE)

				if (inherits(plot_result, "try-error") && iter == 10) {
					cli::cli_alert_warning("Real-time plotting failed - continuing without plots")
				}
			}
			# store per-actor draws / streaming summaries
			if (per_actor_active) {
				per_actor_offset_sum <- per_actor_offset_sum + per_actor_offset
				if (keep_per_actor == "draws" || keep_per_actor == "auto") {
					if (exists("THETA_ACTOR", inherits = FALSE)) {
						THETA_ACTOR[iter, , ] <- theta_actor
					}
				}
				if (keep_per_actor == "summary") {
					theta_actor_n_obs <- theta_actor_n_obs + 1L
					# Welford streaming mean and second moment
					delta_run <- theta_actor - theta_actor_mean
					theta_actor_mean <- theta_actor_mean + delta_run / theta_actor_n_obs
					theta_actor_m2 <- theta_actor_m2 + delta_run *
						(theta_actor - theta_actor_mean)
				}
				RHO_ACTOR[iter] <- rho_actor
				SIGMA_ACTOR[iter] <- sqrt(sigma_actor_sq)
			}
			iter <- iter + 1
		} # post burn-in
		if (verbose && s > burn) {
			cli::cli_progress_update()
		}
		sym_ramp_done <- sym_ramp_done + 1
	} # end MCMC
	if (verbose) {
		cli::cli_progress_done()
		cli::cli_alert_success("MCMC sampling complete")
	}
	total_iters <- nscan + burn
	# a time budget that ran out before the first stored draw leaves nothing to summarise; the fit is
	# still returned so its checkpoint can be resumed
	no_draws_stored <- isTRUE(.lame_terminated_early) &&
		(!exists("iter", inherits = FALSE) || !is.numeric(iter) || iter < 2L)
	if (no_draws_stored) {
		cli::cli_warn(c(
			"{.arg max_seconds} elapsed before any draw was stored; the returned fit carries no posterior draws.",
			"i" = "Increase {.arg max_seconds} or reduce {.arg burn}; a checkpoint written with {.arg checkpoint_path} can be resumed with {.fn lame_resume}."
		))
	}
	# after an early exit under max_seconds, trim every per-draw store to the draws that exist
	n_planned_draws <- nscan / odens
	n_stored_draws <- n_planned_draws
	if (isTRUE(.lame_terminated_early) && exists("iter", inherits = FALSE) &&
			is.numeric(iter) && iter >= 2L && (iter - 1L) < n_planned_draws) {
		n_stored_draws <- iter - 1L
		keep <- seq_len(n_stored_draws)
		trim_draws <- function(x) {
			if (is.null(x)) return(x)
			nd <- length(dim(x))
			if (nd == 0L && length(x) == n_planned_draws) return(x[keep])
			if (nd == 2L && dim(x)[1L] == n_planned_draws) return(x[keep, , drop = FALSE])
			if (nd == 3L && dim(x)[1L] == n_planned_draws) return(x[keep, , , drop = FALSE])
			x
		}
		for (nm in c(
			"BETA", "VC", "ALPHA", "RHO_AB", "SIGMA_AB",
			"RHO_UV", "SIGMA_UV", "PI_SNAP", "SNAP_DRAWS", "SNAP_DRAWS_V",
			"RHO_BETA", "SIGMA_BETA", "RHO", "LOG_LIK", "THETA_ACTOR",
			"RHO_ACTOR", "SIGMA_ACTOR"
		)) {
			if (exists(nm, inherits = FALSE)) assign(nm, trim_draws(get(nm, inherits = FALSE)))
		}
		if (exists("GOF", inherits = FALSE) && !is.null(GOF) && length(dim(GOF)) == 3L &&
				dim(GOF)[3L] == n_planned_draws + 1L) {
			GOF <- GOF[, , c(1L, keep + 1L), drop = FALSE]
		}
		if (exists("U_DRAWS", inherits = FALSE) && !is.null(U_DRAWS) && length(dim(U_DRAWS)) == 4L &&
				dim(U_DRAWS)[4L] == n_planned_draws) {
			U_DRAWS <- U_DRAWS[, , , keep, drop = FALSE]
			V_DRAWS <- V_DRAWS[, , , keep, drop = FALSE]
		}
		# static position draws and bipartite G draws keep the draw in their last margin
		for (nm in c("U_SAMPLES", "V_SAMPLES", "G_SAMPLES")) {
			x <- if (exists(nm, inherits = FALSE)) get(nm, inherits = FALSE) else NULL
			if (!is.null(x) && length(dim(x)) == 3L && dim(x)[3L] == n_planned_draws) {
				assign(nm, x[, , keep, drop = FALSE])
			}
		}
	}
	# compare the stored predictive-draw average with the observed level (.lame_predictive_level_check)
	.lame_predictive_level_check(
		family,
		(isTRUE(dynamic_ab) || isTRUE(dynamic_uv)) && !isTRUE(.lame_terminated_early),
		if (exists("YPS", inherits = FALSE)) YPS else NULL,
		Y,
		nrow(VC)
	)
	# the stationary-variance caps bind without rejecting a draw, so they are reported here rather than in tryErrorChecks
	var_cap_hits <- .lame_cap_env$hits %||% 0L
	var_cap_draws <- .lame_cap_env$draws %||% 0L
	if (var_cap_hits > max(1, round(0.05 * total_iters))) {
		cli::cli_inform(c(
			"i" = "A variance cap removed most of the probability from {var_cap_hits} of {var_cap_draws} conditional updates ({round(100 * var_cap_hits / max(var_cap_draws, 1L))}%).",
			"i" = "The limits ({.arg prior$uv_var_cap}, {.arg prior$ab_var_cap}, {.arg prior$sigma_uv_max}) are part of the prior: the draws are from the posterior under them."
		))
	}
	# a poisson residual variance below 0.1 leaves the dyadic correlation to its prior; say so
	if (identical(family, "poisson") && isTRUE(dcor) && !bip && !symmetric &&
			is.matrix(VC) && "ve" %in% colnames(VC) && nrow(VC) > 0L) {
		.ve_mean <- mean(VC[, "ve"], na.rm = TRUE)
		if (is.finite(.ve_mean) && .ve_mean < 0.1) {
			cli::cli_inform(c(
				"i" = "The residual variance of this Poisson fit is about {signif(.ve_mean, 2)}; with so little latent variation the dyadic correlation {.code rho} rests largely on its prior."
			))
		}
	}
	if (.dyn_uv_has_rank && !.lame_uv_ok(U_cube, V_cube, Inf)) {
		U_cube[!is.finite(U_cube)] <- 0
		V_cube[!is.finite(V_cube)] <- 0
		tryErrorChecks$UV <- tryErrorChecks$UV + 1L
	}
	if (dynamic_ab) {
		if (exists("a_mat") && any(!is.finite(a_mat))) {
			a_mat[!is.finite(a_mat)] <- 0
			tryErrorChecks$ab <- tryErrorChecks$ab + 1L
		}
		if (exists("b_mat") && any(!is.finite(b_mat))) {
			b_mat[!is.finite(b_mat)] <- 0
			tryErrorChecks$ab <- tryErrorChecks$ab + 1L
		}
	}
	sampler_health <- .lame_sampler_health(tryErrorChecks, uv_rail_count, g_rail_count)
	.lame_report_sampler_health(sampler_health, total_iters)
	# start_vals for future runs
	start_vals <- .lame_snapshot_state()

	# output
	if (.dyn_uv_has_rank) {
		# 3d arrays for dynamic UV
		n_draws_uv <- n_stored_draws
		U_final <- U_SUM / n_draws_uv
		V_final <- V_SUM / n_draws_uv
		UVPS_final <- UVPS / n_draws_uv
		# posterior sd of latent positions from the running sum of squares
		if (exists("U_SUM_SQ")) {
			U_var <- U_SUM_SQ / n_draws_uv - U_final^2
			U_var[U_var < 0] <- 0
			U_SD <- sqrt(U_var)
		} else { U_SD <- NULL }
		if (exists("V_SUM_SQ")) {
			V_var <- V_SUM_SQ / n_draws_uv - V_final^2
			V_var[V_var < 0] <- 0
			V_SD <- sqrt(V_var)
		} else { V_SD <- NULL }
	} else {
		U_final <- U
		V_final <- V
		# the last draw is centred the same way as the reported positions, so the raw matrices carry the split
		# the effects and intercept assume
		if (R_eff > 0 && !is.null(U) && !is.null(V) && NCOL(U) > 0L && NCOL(V) > 0L) {
			U_st <- if (length(dim(U)) == 3L) matrix(U[, , 1L], nrow = dim(U)[1L]) else as.matrix(U)
			V_st <- if (length(dim(V)) == 3L) matrix(V[, , 1L], nrow = dim(V)[1L]) else as.matrix(V)
			G_st <- if (bip && !is.null(G)) G else diag(1, ncol(U_st))
			cen_fin <- if (all(is.finite(U_st)) && all(is.finite(V_st)) && all(is.finite(G_st)) &&
					ncol(U_st) == nrow(G_st) && ncol(V_st) == ncol(G_st)) {
				.lame_center_uv_period(U_st, V_st, G_st)
			} else NULL
			if (!is.null(cen_fin)) {
				if (length(dim(U)) == 3L) { for (tt in seq_len(dim(U)[3L])) U_final[, , tt] <- cen_fin$U } else U_final <- cen_fin$U
				if (length(dim(V)) == 3L) { for (tt in seq_len(dim(V)[3L])) V_final[, , tt] <- cen_fin[["V"]] } else V_final <- cen_fin[["V"]]
			}
		}
		UVPS_final <- UVPS / n_stored_draws
	}

	# dynamic additive effects output
	if (dynamic_ab && exists("APS_dyn", inherits = FALSE)) {
		APS_final <- APS_dyn / n_stored_draws
		BPS_final <- BPS_dyn / n_stored_draws
	} else if (dynamic_ab) {
		# no draw was stored: the accumulators were never created
		APS_final <- matrix(0, nrow(a_mat), ncol(a_mat))
		BPS_final <- matrix(0, nrow(b_mat), ncol(b_mat))
	} else {
		# raw sums; get_fit_object divides by nrow(VC)
		APS_final <- APS
		BPS_final <- BPS
	}

	fit <- get_fit_object(
		APS = APS_final,
		BPS = BPS_final,
		UVPS = UVPS_final,
		YPS = YPS,
		BETA = BETA,
		VC = VC,
		GOF = GOF,
		Xlist = Xlist,
		actorByYr = actorByYr,
		colActorByYr = if (bip) colActorByYr else NULL,
		start_vals = start_vals,
		symmetric = symmetric,
		tryErrorChecks = tryErrorChecks,
		model.name = model.name,
		U = U_final,
		V = V_final,
		dynamic_uv = dynamic_uv,
		dynamic_ab = dynamic_ab,
		bip = bip,
		rho_ab = RHO_AB,
		rho_uv = RHO_UV,
		family = family,
		odmax = odmax,
		nA = if (bip) nA else NULL,
		nB = if (bip) nB else NULL,
		n_time = N,
		Y_obs = if (exists("Y_user", inherits = FALSE)) Y_user else Y,
		G = if (bip) G else NULL,
		L_dyn = if (exists("L_SUM", inherits = FALSE) && !is.null(L_SUM)) L_SUM / max(L_N, 1L) else NULL,
		dynamic_beta = beta_dyn$any,
		beta_dynamic_mask = beta_dyn_fit$mask,
		beta_dynamic_groups = beta_dyn_fit$groups,
		rho_beta = rho_beta_by_group,
		sigma_beta = sigma_beta_by_group,
		RHO_BETA = RHO_BETA,
		SIGMA_BETA = SIGMA_BETA,
		dynamic_rho = dynamic_rho,
		RHO = RHO,
		period_exposure = period_exposure
	)
	if (!is.null(A_SAMPLES)) {
		actual_ab <- ab_draw_index <= if (no_draws_stored) 0L else n_stored_draws
		ab_names <- if (any(actual_ab)) paste0("draw", ab_draw_index[actual_ab]) else character(0)
		if (dynamic_ab) {
			fit$a_draws <- A_SAMPLES[, , actual_ab, drop = FALSE]
			fit$b_draws <- B_SAMPLES[, , actual_ab, drop = FALSE]
			dimnames(fit$a_draws) <- c(dimnames(fit$a_dynamic), list(ab_names))
			b_names <- dimnames(fit[["b_dynamic"]] %||% fit[["a_dynamic"]])
			dimnames(fit$b_draws) <- c(b_names, list(ab_names))
			fit$a_level_samples <- A_LEVEL_SAMPLES[, actual_ab, drop = FALSE]
			fit$b_level_samples <- B_LEVEL_SAMPLES[, actual_ab, drop = FALSE]
			dimnames(fit$a_level_samples) <- list(rownames(fit$a_dynamic), ab_names)
			dimnames(fit$b_level_samples) <- list(b_names[[1L]], ab_names)
		} else {
			fit$a_samples <- A_SAMPLES[, actual_ab, drop = FALSE]
			fit$b_samples <- B_SAMPLES[, actual_ab, drop = FALSE]
			dimnames(fit$a_samples) <- list(names(fit$APM), ab_names)
			dimnames(fit$b_samples) <- list(names(fit[["BPM"]] %||% fit[["APM"]]), ab_names)
		}
		fit$ab_draw_index <- ab_draw_index[actual_ab]
	}
	fit$forecast_state <- .lame_forecast_trim(
		FC_STATE,
		if (no_draws_stored) 0L else n_stored_draws
	)
	# attach the log-likelihood matrix so loo.lame() can find it
	if (!is.null(LOG_LIK)) {
		fit$log_lik <- LOG_LIK
		fit$log_lik_index <- obs_idx
	}
	# chunked mode: close all connections and attach metadata
	if (!is.null(LOG_LIK_CONNS)) {
		for (cc in LOG_LIK_CONNS) close(cc)
		fit$log_lik_index <- obs_idx
		fit$log_lik_chunks <- LOG_LIK_META
	}
	# per-period posterior sd of the dynamic additive effects from the streaming second moment;
	# symmetric fits carry only a_dynamic_sd
	if (dynamic_ab && exists("a_dyn_m2", inherits = FALSE) && ab_dyn_n > 1L) {
		.ab_denom <- ab_dyn_n - 1L
		if (!is.null(fit$a_dynamic)) {
			.a_sd <- sqrt(a_dyn_m2 / .ab_denom)
			if (identical(dim(.a_sd), dim(fit$a_dynamic)))
				dimnames(.a_sd) <- dimnames(fit$a_dynamic)
			fit$a_dynamic_sd <- .a_sd
		}
		if (!is.null(fit$b_dynamic)) {
			.b_sd <- sqrt(b_dyn_m2 / .ab_denom)
			if (identical(dim(.b_sd), dim(fit$b_dynamic)))
				dimnames(.b_sd) <- dimnames(fit$b_dynamic)
			fit$b_dynamic_sd <- .b_sd
		}
	}
	# record dynamic argument values on the fit for downstream tools
	fit$dynamic_beta_pool <- dynamic_beta_pool
	# frn has no closed-form observed marginal, so "observed_exact" is recorded as "augmented" (as in ame())
	fit$log_lik_method <- if (family == "frn" &&
			identical(log_lik_method, "observed_exact")) {
		"augmented"
	} else log_lik_method
	fit$dynamic <- list(
		uv = isTRUE(dynamic_uv),
		uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
		ab = isTRUE(dynamic_ab),
		beta = isTRUE(beta_dyn$any),
		beta_kind = if (isTRUE(beta_dyn$any)) dynamic_beta_kind else NULL,
		snap_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap") &&
			isTRUE(.dyn_uv_has_rank),
		t_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
			isTRUE(.dyn_uv_has_rank),
		mode = mode
	)
	# innovation-sd chains of the dynamic blocks (fit$sigma_ab, fit$sigma_uv), which get_fit_object() does not take
	if (dynamic_ab && !is.null(SIGMA_AB)) fit$sigma_ab <- SIGMA_AB
	# posterior means of the actor levels; forecasts and simulations put the AR(1) on the deviations around them
	if (dynamic_ab && exists("A_LEV_SUM", inherits = FALSE)) {
		fit$a_level <- as.numeric(A_LEV_SUM / n_stored_draws)
		fit$b_level <- as.numeric(B_LEV_SUM / n_stored_draws)
		nm_a <- if (is.matrix(fit$a_dynamic)) rownames(fit$a_dynamic) else names(fit$APM)
		nm_b <- if (is.matrix(fit$b_dynamic)) rownames(fit$b_dynamic) else names(fit$BPM)
		if (!is.null(nm_a) && length(nm_a) == length(fit$a_level)) names(fit$a_level) <- nm_a
		if (!is.null(nm_b) && length(nm_b) == length(fit$b_level)) names(fit$b_level) <- nm_b
	}
	if (dynamic_uv && !is.null(SIGMA_UV)) fit$sigma_uv <- SIGMA_UV
	# posterior sd of the dynamic latent positions (same dims as fit$U)
	if (.dyn_uv_has_rank && exists("U_SD") && !is.null(U_SD)) {
		dimnames(U_SD) <- dimnames(fit$U)
		fit$U_sd <- U_SD
	}
	# same for the receiver positions, so latent_positions() reports an sd for the V rows too
	if (.dyn_uv_has_rank && exists("V_SD") && !is.null(V_SD) &&
			!is.null(fit[["V"]]) && identical(dim(V_SD), dim(fit[["V"]]))) {
		dimnames(V_SD) <- dimnames(fit[["V"]])
		fit$V_sd <- V_SD
	}
	# optional draw-level latent trajectories: [actor, dim, period, draw] arrays matching fit$U / fit$V
	if (!is.null(U_DRAWS) && !is.null(V_DRAWS)) {
		.uv_draw_dimnames <- function(mean_cube, draws) {
			dn <- if (!is.null(mean_cube) &&
					identical(dim(mean_cube), dim(draws)[1:3]) &&
					!is.null(dimnames(mean_cube))) {
				dimnames(mean_cube)
			} else {
				list(NULL, NULL, NULL)
			}
			# generic labels where the mean cube carries none
			if (is.null(dn[[2L]])) dn[[2L]] <- paste0("dim", seq_len(dim(draws)[2L]))
			if (is.null(dn[[3L]])) {
				dn[[3L]] <- pd_labs %||% paste0("t", seq_len(dim(draws)[3L]))
			}
			c(dn, list(paste0("draw", seq_len(dim(draws)[4L]))))
		}
		dimnames(U_DRAWS) <- .uv_draw_dimnames(fit$U, U_DRAWS)
		# symmetric signed receiver draws share the sender actor grid
		v_names_src <- if (symmetric) fit[["U"]] else fit[["V"]]
		dimnames(V_DRAWS) <- .uv_draw_dimnames(v_names_src, V_DRAWS)
		fit$U_draws <- U_DRAWS
		fit$V_draws <- V_DRAWS
	}
	# static-fit latent-position draws as [actor, dim, draw] arrays under ame()'s names, fit$U_samples / fit$V_samples
	if (!is.null(U_SAMPLES)) {
		# fit$U is a cube of identical slices on bipartite static fits, so match on the first two margins only
		.uv_sample_dimnames <- function(mean_arr, draws) {
			dn <- list(NULL, NULL)
			if (!is.null(mean_arr) && length(dim(mean_arr)) >= 2L &&
					identical(dim(mean_arr)[1:2], dim(draws)[1:2]) &&
					!is.null(dimnames(mean_arr))) {
				dn <- dimnames(mean_arr)[1:2]
			}
			if (is.null(dn[[2L]])) dn[[2L]] <- paste0("dim", seq_len(dim(draws)[2L]))
			c(dn, list(paste0("draw", seq_len(dim(draws)[3L]))))
		}
		dimnames(U_SAMPLES) <- .uv_sample_dimnames(fit[["U"]], U_SAMPLES)
		fit$U_samples <- U_SAMPLES
		if (!is.null(V_SAMPLES)) {
			# symmetric fits have no V ($V would partial-match $VC): use exact indexing and fall back to U
			v_exact <- fit[["V"]]
			v_names_src <- if (!is.null(v_exact) && length(dim(v_exact)) >= 2L &&
					identical(dim(v_exact)[1:2], dim(V_SAMPLES)[1:2]) &&
					!is.null(dimnames(v_exact))) v_exact else fit[["U"]]
			dimnames(V_SAMPLES) <- .uv_sample_dimnames(v_names_src, V_SAMPLES)
			fit$V_samples <- V_SAMPLES
		}
	}
	# bipartite interaction-matrix draws (as fit$G_samples in ame())
	if (!is.null(G_SAMPLES)) {
		dimnames(G_SAMPLES) <- list(
			NULL,
			NULL,
			paste0("draw", seq_len(dim(G_SAMPLES)[3L]))
		)
		fit$G_samples <- G_SAMPLES
	}
	# snap outputs: posterior snap probabilities, optional indicator draws and the snap-rate chain
	if (isTRUE(fit$dynamic$snap_model_estimated)) {
		fit$pi_snap <- PI_SNAP
		if (!is.null(DELTA_SUM)) {
			.snap_prob_from_counts <- function(sum_mat, count_mat) {
				out <- sum_mat / count_mat
				out[!is.finite(out) | count_mat <= 0L] <- NA_real_
				out
			}
			fit$snap_prob <- .snap_prob_from_counts(DELTA_SUM, DELTA_COUNT)
			if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
				dimnames(fit$snap_prob) <- list(dimnames(fit$U)[[1]], dimnames(fit$U)[[3]])
				names(dimnames(fit$snap_prob)) <- c("actor", "time")
				dimnames(DELTA_COUNT) <- dimnames(fit$snap_prob)
			}
			if (!symmetric) {
				fit$snap_prob_v <- .snap_prob_from_counts(DELTA_SUM_V, DELTA_COUNT_V)
				if (!is.null(fit[["V"]]) && length(dim(fit[["V"]])) == 3) {
					dimnames(fit$snap_prob_v) <- list(dimnames(fit[["V"]])[[1]], dimnames(fit[["V"]])[[3]])
					names(dimnames(fit$snap_prob_v)) <- c("actor", "time")
					dimnames(DELTA_COUNT_V) <- dimnames(fit$snap_prob_v)
				}
			}
			.snap_summary_from_prob <- function(prob, count_mat) {
				n_mat <- count_mat
				dimnames(n_mat) <- dimnames(prob)
				n_mat[!is.finite(prob)] <- 0
				var_mat <- matrix(NA_real_, nrow(prob), ncol(prob), dimnames = dimnames(prob))
				ok <- is.finite(prob) & n_mat > 1L
				var_mat[ok] <- prob[ok] * (1 - prob[ok]) * n_mat[ok] / (n_mat[ok] - 1L)
				list(n = n_mat, mean = prob, variance = var_mat)
			}
			if (SNAP_DRAWS_MODE %in% c("summary", "draws", "chunked")) {
				fit$snap_draws_summary <- .snap_summary_from_prob(fit$snap_prob, DELTA_COUNT)
				if (!symmetric && !is.null(fit$snap_prob_v)) {
					fit$snap_draws_summary_v <- .snap_summary_from_prob(fit$snap_prob_v, DELTA_COUNT_V)
				}
			}
			if (!is.null(SNAP_DRAWS)) {
				dimnames(SNAP_DRAWS) <- list(
					draw = paste0("draw", seq_len(dim(SNAP_DRAWS)[1L])),
					actor = rownames(fit$snap_prob),
					time = colnames(fit$snap_prob)
				)
				fit$snap_draws <- SNAP_DRAWS
				if (!symmetric && !is.null(SNAP_DRAWS_V) && !is.null(fit$snap_prob_v)) {
					dimnames(SNAP_DRAWS_V) <- list(
						draw = paste0("draw", seq_len(dim(SNAP_DRAWS_V)[1L])),
						actor = rownames(fit$snap_prob_v),
						time = colnames(fit$snap_prob_v)
					)
					fit$snap_draws_v <- SNAP_DRAWS_V
				}
				fit$snap_draws_meta <- list(
					mode = SNAP_DRAWS_MODE,
					storage = "in_memory_array",
					dim_order = c("draw", "actor", "time"),
					n_draws = dim(SNAP_DRAWS)[1L],
					first_period = "ineligible_na",
					transition_eligibility = "actor must be present in both adjacent periods"
				)
			}
		}
	}
	# randomised Halton metadata for reproducible GHK scores: the fit seed, the [0, 1) shift and the base
	# chain length n_mc (exact when prior$ghk_n_mc was set)
	if (identical(log_lik_method, "observed_ghk")) {
		fit$log_lik_meta <- list(
			seed = seed,
			halton_shift = .lame_log_lik_halton_shift,
			n_mc = .lame_ghk_n_mc %||% 64L,
			n_mc_fixed = !is.null(.lame_ghk_n_mc),
			method = "GHK + Halton (randomised) + Gauss-Hermite (cbin)",
			dim_cap = 15L
		)
	}
	# explicit ordinal cutpoint posterior on alpha, alpha_1 = 0 not stored; `ordinal_levels` maps
	# categories 1..K back to the observed levels
	fit$ordinal_cutpoints <- ordinal_cutpoints
	if (use_explicit_cutpoints) fit$ordinal_cutpoint_update <-
		if (!symmetric && dcor) "conditional_gibbs" else "marginal_cowles"
	if (use_explicit_cutpoints) {
		fit$ordinal_levels <- ord_lvls
		if (!is.null(ALPHA)) {
			fit$ALPHA <- ALPHA
			fit$alpha_post_mean <- c(0, colMeans(ALPHA, na.rm = TRUE))
		} else {
			# for K = 2 the posterior is a point mass at alpha_1 = 0
			fit$alpha_post_mean <- 0
		}
	}
	# a static bipartite G is reported as the posterior mean of its gauged draws, like the positions
	if (exists("G_sum", inherits = FALSE) && !is.null(G_sum) && G_n > 0L) {
		fit[["G"]] <- G_sum / G_n
	}
	# per-actor slopes (draws or streaming summary)
	if (isTRUE(per_actor_active)) {
		fit$dynamic_beta_per_actor <- dynamic_beta_per_actor
		fit$per_actor_model_version <- 2L
		fit$per_actor_prior <- list(
			rho = list(
				distribution = "scaled_beta",
				shape1 = prior$rho_actor_shape1,
				shape2 = prior$rho_actor_shape2,
				lower = -.99,
				upper = .99
			),
			innovation_variance = list(
				distribution = "inverse_gamma",
				shape = prior$sigma_actor_shape,
				scale = prior$sigma_actor_scale
			),
			initial_state = "stationary",
			centered_dimensions = nrow(theta_actor) - 1L
		)
		fit$per_actor_design_idx <- per_actor_design_idx
		fit$per_actor_covariate_idx <- per_actor_covariate_idx
		fit$per_actor_identifiability <- per_actor_identifiability
		fit$RHO_ACTOR <- RHO_ACTOR
		fit$SIGMA_ACTOR <- SIGMA_ACTOR
		if (keep_per_actor %in% c("draws", "auto") &&
				exists("THETA_ACTOR", inherits = FALSE)) {
			fit$THETA_ACTOR <- THETA_ACTOR
			# posterior-mean summary
			fit$theta_actor_mean <- apply(THETA_ACTOR, c(2, 3), mean)
			fit$theta_actor_sd <- apply(THETA_ACTOR, c(2, 3), stats::sd)
		} else if (keep_per_actor == "summary") {
			fit$theta_actor_mean <- theta_actor_mean
			fit$theta_actor_sd <- sqrt(
				theta_actor_m2 /
					max(1L, theta_actor_n_obs - 1L)
			)
		}
		# actor / period labels for the per-actor summaries
		.pa_margin <- if (identical(dynamic_beta_per_actor, "col")) 2L else 1L
		.pa_names <- dimnames(Y)[[.pa_margin]]
		if (is.null(.pa_names)) .pa_names <- tryCatch({
			if (.pa_margin == 2L) colnames(fit$YPM[[1]]) else rownames(fit$YPM[[1]])
		}, error = function(e) NULL)
		.pa_periods <- tryCatch(names(fit$YPM), error = function(e) NULL)
		if (!is.null(fit$THETA_ACTOR) && length(.pa_names) == dim(fit$THETA_ACTOR)[2L])
			dimnames(fit$THETA_ACTOR) <- list(NULL, .pa_names, .pa_periods)
		for (.fld in c("theta_actor_mean", "theta_actor_sd")) {
			.m <- fit[[.fld]]
			if (is.matrix(.m) && !is.null(.pa_names) &&
					length(.pa_names) == nrow(.m)) {
				dimnames(.m) <- list(
					.pa_names,
					if (!is.null(.pa_periods) && length(.pa_periods) == ncol(.m))
						.pa_periods else NULL
				)
				fit[[.fld]] <- .m
			}
		}
	}
	if (!is.null(time_index)) fit$time_index <- time_index
	if (!is.null(period_exposure)) fit$period_exposure <- period_exposure
	if (is.finite(max_seconds)) fit$max_seconds <- max_seconds
	fit$terminated_early <- isTRUE(.lame_terminated_early)
	fit$var_cap_hits <- var_cap_hits
	fit$var_cap_draws <- var_cap_draws
	# events per actor-period and per actor, the rates the reporting methods read (.lame_ab_sparse)
	if (is.finite(.per_actor_period)) fit$events_per_actor_period <- .per_actor_period
	if (is.finite(.per_actor_total)) fit$events_per_actor <- .per_actor_total
	# whether the effects prior was supplied (.lame_ab_var_unident_reason)
	fit$prior_supplied <- .ab_prior_supplied
	# share of observed rows at their nomination cap (cbin / frn only)
	if (is.finite(.odmax_sat_share)) {
		fit$odmax_saturated_share <- .odmax_sat_share
		fit$odmax_sat_threshold <- .lame_odmax_threshold_for(.odmax_n_per)
	}
	fit$cells_hidden <- .cells_hidden
	# whether the dynamic positions started from a dynamic ALS fit
	if (isTRUE(dynamic_uv)) fit$als_start <- als_start_used
	# acceptance record of the joint scale move, when it ran
	if (marginal_move) fit$marginal_move <- c(
		attempts = marginal_att,
		accepts = marginal_acc,
		rate = marginal_acc / max(marginal_att, 1L),
		step_sd = marginal_sd,
		numerical_rejections = marginal_bad
	)
	if ((isTRUE(level_move) || isTRUE(bip_level_move)) && lev_att > 0L) {
		fit$level_move <- c(
			attempts = lev_att,
			accepts = lev_acc,
			rate = lev_acc / lev_att,
			step_sd = lev_sd
		)
	}
	if (isTRUE(dev_move) && dev_att > 0L) {
		fit$dev_move <- c(
			attempts = dev_att,
			accepts = dev_acc,
			rate = dev_acc / dev_att,
			step_sd = dev_sd
		)
	}
	# acceptance record of the interweaving step of the dynamic_beta hyperparameters
	if (isTRUE(beta_hyper_nc) && beta_nc_stats[["attempts"]] > 0) {
		fit$beta_hyper_move <- c(
			sweeps = beta_nc_stats[["sweeps"]],
			attempts = beta_nc_stats[["attempts"]],
			accepts = beta_nc_stats[["accepts"]],
			rate = beta_nc_stats[["accepts"]] / beta_nc_stats[["attempts"]]
		)
	}
	# the collapsed per-actor update records its visits and the likelihood evaluations of its slice steps
	if ((actor_move || actor_move_bip) && actor_move_stats["visits"] > 0) fit$actor_move <- actor_move_stats
	if (ltm_on && ltm_stats["visits"] > 0) fit$actor_move <- ltm_stats
	if (spm_on && spm_stats["visits"] > 0) fit$sparse_position_move <- spm_stats
	if (npm_on && npm_stats["visits"] > 0) fit$sparse_position_move <- npm_stats
	if (cr_on && cr_stats["ridge_moved"] > 0) fit$cov_ridge_move <- cr_stats
	# ordinal scale move: visits and the visits at which the scale of the whole state moved
	if (osm_on && osm_stats[["visits"]] > 0) fit$ordinal_scale_move <- osm_stats
	if (pab_on && pab_stats[["visits"]] > 0) fit$ab_period_move <- pab_stats
	if (uv_split_on && uv_split_stats[["visits"]] > 0) fit$uv_split_move <- uv_split_stats
	if (uv_nc_move && (uv_nc_stats[["sym_calls"]] + uv_nc_stats[["bip_calls"]]) > 0) fit$uv_nc_move <- uv_nc_stats
	if (!is.null(uv_temper_dt) && uv_temper_stats[["calls"]] > 0) fit$uv_temper_move <- uv_temper_stats
	if (bso_on && bso_stats[["visits"]] > 0) fit$binary_scale_move <- bso_stats
	if (dpm_sym && dpm_sym_stats[["dpm_scale_attempts"]] > 0) {
		fit$actor_move <- if (is.null(fit$actor_move)) dpm_sym_stats else c(fit$actor_move, dpm_sym_stats)
	}
	# the position-path move of a binary fit adds its counts to the collapsed update's record
	if (dpm_bin && dpm_stats[["dpm_visits"]] > 0) {
		fit$actor_move <- if (is.null(fit$actor_move)) c(
			visits = dpm_stats[["dpm_visits"]],
			dpm_stats[-1L]
		) else c(fit$actor_move, dpm_stats[-1L])
	}
	# translations of the poisson coefficients: visits, period and path translations, static draws and slice evaluations
	if ((ppm_on || ppm_scale_on) && ppm_stats[["visits"]] > 0) fit$poisson_period_move <- ppm_stats
	# interweaving step of the actor-specific slope hyperparameters: sweeps, proposals and accepts
	if (isTRUE(pa_hyper_nc) && pa_nc_stats[["attempts"]] > 0) {
		fit$pa_hyper_move <- c(
			sweeps = pa_nc_stats[["sweeps"]],
			attempts = pa_nc_stats[["attempts"]],
			accepts = pa_nc_stats[["accepts"]],
			rate = pa_nc_stats[["accepts"]] / pa_nc_stats[["attempts"]]
		)
	}
	# counts of the nomination-family receiver and level moves
	if ((nom_col && nom_col_stats["visits"] > 0) || (nom_lvl && nom_lvl_stats["visits"] > 0) || (nom_rho && nom_rho_stats[["calls"]] > 0)) {
		nom_rec <- if (nom_col) nom_col_stats else if (nom_lvl) c(
			visits = nom_lvl_stats[["visits"]]
		) else c(visits = nom_rho_stats[["calls"]])
		if (nom_lvl) nom_rec <- c(
			nom_rec,
			level_calls = nom_lvl_stats[["calls"]],
			level_rows = nom_lvl_stats[["free_rows"]]
		)
		if (nom_per && nom_per_stats[["visits"]] > 0) nom_rec <- c(
			nom_rec,
			period_calls = nom_per_stats[["calls"]],
			period_comps = nom_per_stats[["comps"]]
		)
		if (nom_room && nom_room_stats[["visits"]] > 0) nom_rec <- c(
			nom_rec,
			room_attempts = nom_room_stats[["attempts"]],
			room_accepts = nom_room_stats[["accepts"]],
			room_stuck = nom_room_stats[["stuck"]]
		)
		if (nom_rho && nom_rho_stats[["calls"]] > 0) nom_rec <- c(
			nom_rec,
			rho_attempts = nom_rho_stats[["attempts"]],
			rho_accepts = nom_rho_stats[["accepts"]]
		)
		fit$actor_move <- nom_rec
	}
	if (ab_hyper_move && !is.null(ab_hyper_par) && ab_hyper_par[3L] > 0) {
		fit$ab_hyper_move <- c(
			attempts = ab_hyper_par[3L],
			accepts = ab_hyper_par[4L],
			rate = ab_hyper_par[4L] / ab_hyper_par[3L]
		)
	}
	if (uv_g_scale_stats[["visits"]] > 0) fit$uv_g_scale_move <- uv_g_scale_stats
	if (bip_scale_att > 0L) fit$bipartite_scale_move <- c(
		attempts = bip_scale_att,
		accepts = bip_scale_acc,
		rate = bip_scale_acc / bip_scale_att,
		bound_rejections = bip_scale_bounds,
		numerical_rejections = bip_scale_numerical
	)
	if (ncar_uv_att > 0L) fit$nc_ar_move <- c(
		uv_attempts = ncar_uv_att,
		uv_accepts = ncar_uv_acc,
		ab_attempts = ncar_ab_att,
		ab_accepts = ncar_ab_acc
	)
	if (actor_slice_rounds > 0L) fit$actor_factor_slice <- c(
		rounds = actor_slice_rounds,
		blocks = actor_slice_blocks,
		evaluations = actor_slice_evaluations,
		bound_rejections = actor_slice_bounds,
		numerical_rejections = actor_slice_bad
	)
	if (dhmc_att > 0L) fit$dynamic_marginal_move <- c(
		attempts = dhmc_att,
		accepts = dhmc_acc,
		rate = dhmc_acc / dhmc_att,
		step_sd = dhmc_sd,
		numerical_rejections = dhmc_nonfinite,
		bound_rejections = dhmc_bounds
	)
	if (joint_levels_enabled) fit$joint_levels_move <- c(
		attempts = joint_levels_attempts,
		numerical_failures = tryErrorChecks$joint_levels
	)
	if (rho_collapsed_att > 0L) fit$rho_collapsed_move <- c(
		attempts = rho_collapsed_att,
		accepts = rho_collapsed_acc,
		rate = rho_collapsed_acc / rho_collapsed_att
	)
	if (uv_scale_att > 0L) fit$uv_scale_move <- c(
		attempts = uv_scale_att,
		accepts = uv_scale_acc,
		rate = uv_scale_acc / uv_scale_att
	)
	if (translation_att > 0L) fit$translation_move <- c(attempts = translation_att)
	# the fixed inverse-gamma scale of the deviation innovation variance
	if (isTRUE(dynamic_ab) && !is.null(prior$sigma_ab_scale)) {
		fit$sigma_ab_prior_scale <- as.numeric(prior$sigma_ab_scale)
	}
	# the same, for the dynamic-coefficient innovation variance
	if (isTRUE(beta_dyn$any) && !is.null(prior$sigma_beta_scale)) {
		fit$sigma_beta_prior_scale <- as.numeric(prior$sigma_beta_scale)
	}
	if (!is.null(checkpoint_path)) fit$checkpoint_path <- checkpoint_path
	# scalar typing
	if (!is.null(fit$s2)) fit$s2 <- as.numeric(fit$s2)
	# normalize EZ to a list of matrices
	.as_ez_list <- function(EZ) {
		if (is.null(EZ)) return(NULL)
		if (is.list(EZ)) return(EZ)
		if (length(dim(EZ)) == 3L) {
			Tn <- dim(EZ)[3]
			out <- vector("list", Tn)
			for (t in seq_len(Tn)) out[[t]] <- EZ[, , t]
			return(out)
		}
		list(EZ)
	}
	fit$EZ <- .as_ez_list(fit$EZ)
	if (per_actor_active && !no_draws_stored) {
		for (tt in seq_len(N)) {
			offset <- matrix(
				per_actor_offset_sum[, , tt] / n_stored_draws,
				dim(Y)[1L],
				dim(Y)[2L],
				dimnames = dimnames(Y)[1:2]
			)
			rn <- rownames(fit$EZ[[tt]]); cn <- colnames(fit$EZ[[tt]])
			if (!is.null(rn) && !is.null(rownames(offset))) offset <- offset[rn, , drop = FALSE]
			if (!is.null(cn) && !is.null(colnames(offset))) offset <- offset[, cn, drop = FALSE]
			fit$EZ[[tt]] <- fit$EZ[[tt]] + offset
		}
	}
	# latent positions from the posterior mean products, aligned over time
	if (isTRUE(dynamic_uv)) {
		fit <- .lame_dynamic_positions(fit)
		rho_hat <- estimate_rho_from_U(fit$U)
		if (is.finite(rho_hat)) fit$rho_uv_aligned <- as.numeric(rho_hat)
	}
	# posterior sd of the static additive effects, for simulate(); dynamic effects carry a_dynamic_sd / b_dynamic_sd
	if (!dynamic_ab && exists("APS2", inherits = FALSE) && length(APS2) == length(fit$APM)) {
		fit$APM_sd <- stats::setNames(sqrt(pmax(APS2 / nrow(fit$VC) - fit$APM^2, 0)), names(fit$APM))
		if (!is.null(fit$BPM) && length(BPS2) == length(fit$BPM)) {
			fit$BPM_sd <- stats::setNames(sqrt(pmax(BPS2 / nrow(fit$VC) - fit$BPM^2, 0)), names(fit$BPM))
		}
	}
	fit$call <- mc
	fit$refit_options <- .lame_refit_options
	fit$refit_data <- .lame_refit_data
	# record the dynamic_beta_kind used, and the one requested, for summaries and forecasting
	if (isTRUE(beta_dyn$any)) {
		fit$dynamic_beta_kind <- dynamic_beta_kind
		if (!identical(dynamic_beta_kind, dynamic_beta_kind_requested)) {
			fit$dynamic_beta_kind_requested <- dynamic_beta_kind_requested
		}
	}
	if (isTRUE(dynamic_uv)) {
		fit$dynamic_uv_kind <- dynamic_uv_kind
		if (!identical(dynamic_uv_kind, dynamic_uv_kind_requested)) {
			fit$dynamic_uv_kind_requested <- dynamic_uv_kind_requested
		}
	}
	# the pre-conversion data snapshot, so update(fit, ...) can refit the original input
	if (isTRUE(freeze_call)) {
		fit$data_snapshot <- .lame_freeze_snapshot
		fit$freeze_call <- TRUE
	}
	# the resolved g and the prior list with sampler defaults filled in, for prior_summary()
	fit$g <- g
	fit$prior <- prior
	# the coefficient kernel this fit used
	fit$coefficient_prior <- if (!is.null(dynamic_beta_per_actor)) {
		list(
			kind = "actor_dynamic", description = "Actor-specific coefficient paths use the fitted initialization and transition priors."
		)
	} else if (isTRUE(beta_dyn$any)) {
		list(
			kind = "dynamic", description = "Dynamic coefficients use the fitted initialization and transition priors; static coefficients have independent N(0, g_j * s2) priors."
		)
	} else if (bip) {
		list(
			kind = "independent_scaled", description = "Regression coefficients have independent N(0, g[1] * s2) priors."
		)
	} else {
		list(
			kind = "panel_design", description = "Regression coefficients have zero-mean Gaussian prior precision sum_t(X_t'X_t) / (T * g[1]), independent of residual s2."
		)
	}
	if (!bip && !symmetric && R > 0L && !isTRUE(.dyn_uv_has_rank)) {
		fit$static_uv_prior <- list(
			distribution = "inverse_wishart",
			df = uv_prior_static$df,
			scale = diag(uv_prior_static$psi, 2L * R),
			dimension = 2L * R,
			data_scale = uv_prior_static$scale,
			df_multiplier = uv_prior_static$df_multiplier
		)
	}
	if (!is.null(s2_prior_scale)) fit$s2_prior_scale <- s2_prior_scale
	# raw sampling-failure counts for auditing; `fit$mh_counters` is an alias
	fit$tryErrorChecks <- tryErrorChecks
	fit$sampler_health <- sampler_health
	fit$sampler_health_iterations <- total_iters
	# fraction of sweeps with a rejected conditional or a numerical failure
	fit$refusal_share <- c(
		UV = (tryErrorChecks$UV %||% 0L) / max(total_iters, 1L),
		G = (tryErrorChecks$G %||% 0L) / max(total_iters, 1L)
	)
	# cbin / frn with hidden cells: the intercept's effective sample size (.lame_nom_hidden_ess)
	fit$hidden_cell_ess <- .lame_nom_hidden_ess(fit, .cells_hidden)
	fit$mh_counters <- tryErrorChecks
	# structural flags and latent dims so format_ame_call() can render the model formula
	fit$intercept <- isTRUE(intercept)
	fit$rvar <- isTRUE(rvar)
	fit$cvar <- isTRUE(cvar)
	fit$dcor <- isTRUE(dcor)
	if (!is.null(LOG_LIK) || !is.null(LOG_LIK_CONNS)) {
		fit$log_lik_contract <- .lame_score_contract(
			family,
			fit$log_lik_method,
			ordinal_cutpoints,
			bip,
			symmetric,
			dcor || isTRUE(dynamic_rho),
			obs_idx,
			if (bip || symmetric) NULL else rho
		)
	}
	fit$R <- if (bip) max(RA, RB) else R
	if (bip) {
		fit$R_row <- RA
		fit$R_col <- RB
	}
	.lame_report_joint_sampler_health(fit)
	.lame_report_bipartite_scale_health(fit)
	.lame_report_actor_slice_health(fit)
	# the chain policy: extend a chain that is still moving, or refuse it
	if (isTRUE(.chain_ctx$record)) fit$chain_summary <- .lame_chain_acc_finish(.chain_acc)
	.chain_out <- .lame_chain_policy(fit, .chain_ctx, environment())
	if (isTRUE(.chain_out$done)) return(.chain_out$fit)
	fit$auto_extend <- .chain_out$record
	.lame_separation_check(fit)
	# a run stopped before its first stored draw returns empty draw arrays
	if (no_draws_stored) fit <- .lame_empty_draws(fit)
	class(fit) <- c("lame", "ame")
	return(fit)

}
