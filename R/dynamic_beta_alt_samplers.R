# joint Gaussian samplers for dynamic_beta paths of kind rw2 or matern32, and
# for AR(1) paths over unequal time gaps; each kind supplies a T x T prior
# precision K (.lame_beta_change_precision), and the path is drawn from its
# joint Gaussian conditional

#' @noRd
# tridiagonal AR(1) precision K over time gaps, with transition factor
# rho^gap and variance (1 - rho^(2 gap)) / (1 - rho^2), so beta ~ N(0, sigma^2 K^{-1})
.ar1_precision <- function(T_per, gaps = NULL, rho = 0, jitter = 1e-10) {
	if (T_per < 2L) return(matrix(1, 1L, 1L))
	if (is.null(gaps)) gaps <- rep(1, T_per - 1L)
	rho <- max(min(rho, 0.999999), -0.999999)
	v_stat <- 1 / (1 - rho^2)
	K <- matrix(0, T_per, T_per)
	for (t in seq_len(T_per - 1L)) {
		dt <- as.numeric(gaps[t])
		Ft <- rho^dt
		qt <- (1 - rho^(2 * dt)) / (1 - rho^2)
		# terms from -0.5 * (beta_{t+1} - Ft * beta_t)^2 / qt
		K[t, t] <- K[t, t] + Ft^2 / qt
		K[t + 1, t + 1] <- K[t + 1, t + 1] + 1 / qt
		K[t, t + 1] <- K[t, t + 1] - Ft / qt
		K[t + 1, t] <- K[t + 1, t] - Ft / qt
	}
	# stationary initial state beta_1 ~ N(0, v_stat)
	K[1, 1] <- K[1, 1] + 1 / v_stat
	K + diag(jitter, T_per)
}

# prior precision, in units of sigma^2 lam, of the change of an rw2 or matern32 path since its first
# period: rw2 a first increment and second differences (the first increment over gap g has variance
# sigma^2 lam g, as under rw1), matern32 the change since period 1 of a Matern 3/2 process (length scale
# as in .matern32_precision()). its null space is the constant, the level, which has its own prior
.lame_beta_change_precision <- function(kind, T_per, length_scale = NULL, time_positions = seq_len(T_per)) {
	tt <- as.numeric(time_positions)
	gaps <- diff(tt)
	if (identical(kind, "rw2")) {
		D <- matrix(0, T_per - 1L, T_per)
		D[1L, 1:2] <- c(-1, 1) / sqrt(gaps[1L])
		for (i in seq_len(T_per - 2L)) {
			weight <- sqrt(mean(gaps[i + 0:1]))
			D[i + 1L, i] <- 1 / gaps[i] / weight
			D[i + 1L, i + 1L] <- -(1 / gaps[i] + 1 / gaps[i + 1L]) / weight
			D[i + 1L, i + 2L] <- 1 / gaps[i + 1L] / weight
		}
		return(crossprod(D))
	}
	if (is.null(length_scale) || !is.finite(length_scale) || length_scale <= 0) length_scale <- max(2, T_per / 4)
	a <- sqrt(3) * abs(outer(tt, tt, "-")) / length_scale
	M <- cbind(-1, diag(T_per - 1L))
	crossprod(M, solve(M %*% ((1 + a) * exp(-a)) %*% t(M), M))
}

#' @noRd
.dynamic_beta_kind_to_prec <- function(
	kind,
	T_per,
	length_scale = NULL,
	time_positions = NULL,
	gaps = NULL,
	rho = NULL
) {
	# integer positions when none are supplied
	if (is.null(time_positions)) time_positions <- seq_len(T_per)
	if (is.null(gaps)) gaps <- diff(time_positions)
	switch(
		kind,
		ar1 = .ar1_precision(T_per, gaps = gaps, rho = if (is.null(rho)) 0 else rho),
		rw1 = .ar1_precision(T_per, gaps = gaps, rho = 0.999999),
		rw2 = .lame_beta_change_precision("rw2", T_per, time_positions = time_positions),
		matern32 = .lame_beta_change_precision("matern32", T_per, length_scale, time_positions),
		cli::cli_abort(sprintf("kind '%s' not handled by alt sampler", kind))
	)
}

# Gaussian likelihood sufficient statistics with the observed dyads and
# reciprocal covariance of the compiled coefficient sampler
.lame_beta_observation_stats <- function(
	X,
	residual,
	s2,
	rho = 0,
	bipartite = FALSE,
	symmetric = FALSE
) {
	y <- as.numeric(residual)
	keep <- is.finite(y)
	if (!bipartite) {
		keep <- keep & as.vector(if (symmetric) upper.tri(residual) else row(residual) != col(residual))
	}
	if (!any(keep)) return(list(Q = matrix(0, ncol(X), ncol(X)), h = numeric(ncol(X))))
	Xk <- X[keep, , drop = FALSE]
	yk <- y[keep]
	QX <- Xk / s2
	Qy <- yk / s2
	if (!bipartite && !symmetric && rho != 0) {
		partner <- as.vector(t(matrix(seq_along(y), nrow(residual))))[keep]
		paired <- is.finite(y[partner])
		QX[paired, ] <- (Xk[paired, , drop = FALSE] - rho * X[partner[paired], , drop = FALSE]) /
			(s2 * (1 - rho^2))
		Qy[paired] <- (yk[paired] - rho * y[partner[paired]]) / (s2 * (1 - rho^2))
	}
	list(Q = crossprod(Xk, QX), h = drop(crossprod(Xk, Qy)))
}

# unit-innovation AR(1) transition factors and variances over the time gaps;
# the first state has its own prior, independent of rho and sigma
.lame_beta_gap_transition <- function(rho, gaps) {
	if (any(!is.finite(gaps) | gaps <= 0)) stop("Coefficient time gaps must be positive.")
	if (rho < 0 && any(abs(gaps - round(gaps)) > 1e-10)) {
		stop("Negative coefficient persistence requires integer time gaps.")
	}
	f <- rho^gaps
	q <- if (abs(abs(rho) - 1) < 1e-12) gaps else (1 - rho^(2 * gaps)) / (1 - rho^2)
	if (any(!is.finite(q) | q <= 0)) stop("Invalid coefficient transition variance.")
	list(f = f, q = q)
}

.lame_beta_transition_precision <- function(rho, gaps) {
	nt <- length(gaps) + 1L
	z <- .lame_beta_gap_transition(rho, gaps)
	D <- matrix(0, nt - 1L, nt)
	for (t in seq_along(gaps)) {
		D[t, t] <- -z$f[t] / sqrt(z$q[t])
		D[t, t + 1L] <- 1 / sqrt(z$q[t])
	}
	crossprod(D)
}

# draw the whole coefficient path from its joint Gaussian conditional
.sample_beta_dynamic_alt <- function(
	Xdyn_list,
	Xstat_list,
	Z_list,
	offset_list,
	beta_static,
	sigma_by_coef,
	Lambda,
	s2,
	kind,
	length_scale = NULL,
	time_positions = NULL,
	rho_by_coef = NULL,
	dyad_rho = 0,
	use_dyad_rho = FALSE,
	bipartite = FALSE,
	symmetric = FALSE,
	beta0_mean = NULL,
	beta0_cov = NULL
) {
	nt <- length(Z_list)
	p <- ncol(Xdyn_list[[1L]])
	if (p == 0L) return(list(path = matrix(0, nt, 0L), chol_fail = 0L))
	if (is.null(time_positions)) time_positions <- seq_len(nt)
	if (is.null(beta0_mean)) beta0_mean <- numeric(p)
	if (is.null(beta0_cov)) beta0_cov <- diag(10, p)
	if (is.null(rho_by_coef)) rho_by_coef <- rep(0, p)
	sigma_by_coef <- rep_len(sigma_by_coef, p)
	Lambda <- as.matrix(Lambda)
	if (any(Lambda != diag(diag(Lambda), p))) stop("Alternative coefficient evolution requires diagonal Lambda.")
	if (any(!is.finite(diag(Lambda)) | diag(Lambda) <= 0)) stop("Invalid coefficient scale.")
	Q <- matrix(0, nt * p, nt * p)
	h <- numeric(nt * p)
	for (k in seq_len(p)) {
		idx <- seq_len(nt) + (k - 1L) * nt
		K <- if (kind %in% c("ar1", "rw1")) {
			.lame_beta_transition_precision(if (kind == "rw1") 1 else rho_by_coef[k], diff(time_positions))
		} else .dynamic_beta_kind_to_prec(kind, nt, length_scale, time_positions)
		Q[idx, idx] <- K / (sigma_by_coef[k]^2 * Lambda[k, k])
	}
	# the first state: its prior under ar1 and rw1, the level of an rw2 or matern32 path
	idx <- 1L + (seq_len(p) - 1L) * nt
	Q0 <- solve(beta0_cov)
	Q[idx, idx] <- Q[idx, idx] + Q0
	h[idx] <- Q0 %*% beta0_mean
	for (t in seq_len(nt)) {
		r <- Z_list[[t]] - offset_list[[t]]
		if (length(beta_static)) r <- r - matrix(Xstat_list[[t]] %*% beta_static, nrow(r), ncol(r))
		zz <- .lame_beta_observation_stats(
			Xdyn_list[[t]],
			r,
			s2,
			if (use_dyad_rho) rep_len(dyad_rho, nt)[t] else 0,
			bipartite,
			symmetric
		)
		idx <- t + (seq_len(p) - 1L) * nt
		Q[idx, idx] <- Q[idx, idx] + zz$Q
		h[idx] <- h[idx] + zz$h
	}
	L <- chol(Q)
	mu <- backsolve(L, forwardsolve(t(L), h))
	list(path = matrix(mu + backsolve(L, rnorm(nt * p)), nt, p), chol_fail = 0L)
}

# draw one innovation sd per coefficient group, pooling quadratic forms and
# degrees of freedom across the group; a `stationary` ar1 group also counts the
# first state of each path, N(0, sigma^2 lam / (1 - rho^2))
.sample_sigma_beta_alt <- function(
	beta_path,
	kind,
	length_scale = NULL,
	prior_shape,
	prior_scale,
	Lambda = diag(ncol(beta_path)),
	group_id = seq_len(ncol(beta_path)),
	n_groups = max(group_id),
	time_positions = NULL,
	rho_by_group = rep(0, n_groups),
	stationary = rep(FALSE, n_groups)
) {
	nt <- nrow(beta_path)
	if (is.null(time_positions)) time_positions <- seq_len(nt)
	shape <- rep_len(prior_shape, n_groups)
	rate <- rep_len(prior_scale, n_groups)
	for (k in seq_len(ncol(beta_path))) {
		g <- group_id[k]
		K <- if (kind %in% c("ar1", "rw1")) {
			.lame_beta_transition_precision(if (kind == "rw1") 1 else rho_by_group[g], diff(time_positions))
		} else .dynamic_beta_kind_to_prec(kind, nt, length_scale, time_positions)
		# every kind informs sigma through nt - 1 changes; the first state has its own prior
		df <- nt - 1L
		if (identical(kind, "ar1") && isTRUE(stationary[g])) {
			K[1L, 1L] <- K[1L, 1L] + 1 - rho_by_group[g]^2
			df <- nt
		}
		shape[g] <- shape[g] + df / 2
		rate[g] <- rate[g] + drop(crossprod(beta_path[, k], K %*% beta_path[, k])) / (2 * Lambda[k, k])
	}
	sqrt(1 / rgamma(n_groups, shape = shape, rate = rate))
}

# Metropolis-Hastings update of the persistence rho over unequal time gaps; the
# paths of a `stationary` group are deviations whose first state has variance
# sigma^2 lam / (1 - rho^2)
.lame_sample_rho_beta_gaps <- function(
	beta_path,
	group_id,
	n_groups,
	Lambda_inv,
	sigma_by_coef,
	rho_current,
	rho_prior_mean,
	rho_prior_sd,
	rho_lower,
	rho_upper,
	time_positions,
	stationary = rep(FALSE, n_groups)
) {
	gaps <- diff(time_positions)
	log_target <- function(rho, idx, g) {
		if (rho <= rho_lower || rho >= rho_upper ||
				(rho < 0 && any(abs(gaps - round(gaps)) > 1e-10))) return(-Inf)
		z <- .lame_beta_gap_transition(rho, gaps)
		x <- beta_path[, idx, drop = FALSE]
		err <- x[-1L, , drop = FALSE] - z$f * x[-nrow(x), , drop = FALSE]
		ss <- colSums(err^2 / z$q) * diag(Lambda_inv)[idx] / sigma_by_coef[idx]^2
		lp <- dnorm(rho, rho_prior_mean[g], rho_prior_sd[g], log = TRUE) -
			.5 * (sum(ss) + length(idx) * sum(log(z$q)))
		if (isTRUE(stationary[g])) {
			lp <- lp + .5 * length(idx) * log1p(-rho^2) -
				.5 * (1 - rho^2) * sum(x[1L, ]^2 * diag(Lambda_inv)[idx] / sigma_by_coef[idx]^2)
		}
		lp
	}
	for (g in seq_len(n_groups)) {
		idx <- which(group_id == g)
		for (step in seq_len(4L)) {
			proposal <- rho_current[g] + rnorm(1, sd = .08)
			lr <- log_target(proposal, idx, g) - log_target(rho_current[g], idx, g)
			if (is.finite(lr) && log(runif(1)) < lr) rho_current[g] <- proposal
		}
	}
	rho_current
}
