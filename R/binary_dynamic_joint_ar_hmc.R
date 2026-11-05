# joint HMC update of whitened paths, stationary scale and transformed persistence
# the mass matrix uses a fixed prior reference, independent of proposed endpoints

.lame_binary_dynamic_ar_hmc_hyper <- function(
	cx,
	P,
	Sab,
	rho_ab,
	sigma_ab,
	rho,
	shape,
	scale,
	prior_mean,
	prior_sd,
	variance_max = Inf,
	stationary_max = Inf,
	max_abs = Inf
) {
	if (any(!is.finite(c(shape, scale, prior_mean, prior_sd))) || shape <= 0 || scale <=
			0 || prior_sd <= 0)
		stop(
			"Joint dynamic factor HMC requires positive variance-prior shape/scale and persistence-prior SD, with a finite persistence-prior mean.",
			call. = FALSE
		)
	ref_rho <- max(-0.98, min(0.98, prior_mean))
	ref_sigma <- sqrt(scale / (shape + 1))
	h <- .lame_binary_dynamic_pair_hmc_hyper(
		cx,
		P,
		Sab,
		rho_ab,
		sigma_ab,
		ref_rho,
		ref_sigma,
		rho,
		max_abs
	)
	h$shape <- shape
	h$scale <- scale
	h$prior_mean <- prior_mean
	h$prior_sd <- prior_sd
	h$variance_max <- variance_max
	h$stationary_max <- stationary_max
	h$limit <- 0.99
	h$reference_rho <- ref_rho
	h
}

.lame_binary_dynamic_ar_hmc_actual <- function(h, l, alpha, cx) {
	r <- h$limit * tanh(alpha)
	w <- exp(2 * l)
	v <- w * (1 - r * r)
	h$uv_root <- chol(w * r^abs(outer(seq_len(cx$nt), seq_len(cx$nt), "-")))
	h$uv_variance <- w
	h$sigma <- sqrt(v)
	h$rho_uv <- r
	h
}

.lame_binary_dynamic_ar_hmc_encode <- function(
	beta,
	level,
	adev,
	bdev,
	U,
	V,
	sigma,
	rho_uv,
	cx,
	h
) {
	stopifnot(abs(rho_uv) < h$limit)
	l <- log(sigma) - 0.5 * log1p(-rho_uv * rho_uv)
	a <- atanh(rho_uv / h$limit)
	ah <- .lame_binary_dynamic_ar_hmc_actual(h, l, a, cx)
	c(.lame_binary_dynamic_pair_hmc_encode(beta, level, adev, bdev, U, V, cx, ah), l, a)
}

.lame_binary_dynamic_ar_hmc_decode <- function(q, cx, h) {
	ah <- .lame_binary_dynamic_ar_hmc_actual(h, q[cx$size + 1L], q[cx$size + 2L], cx)
	z <- .lame_binary_dynamic_pair_hmc_decode(q[seq_len(cx$size)], cx, ah)
	z$sigma <- ah$sigma
	z$rho_uv <- ah$rho_uv
	z
}

.lame_binary_dynamic_ar_hmc_path_derivative <- function(path, r) {
	d <- matrix(0, nrow(path), ncol(path))
	if (ncol(path) > 1L)
		for (t in 2:ncol(path)) d[, t] <- path[, t - 1L] + r * d[, t - 1L] - r / (1 - r * r) *
		(path[, t] - r * path[, t - 1L])
	d
}

.lame_binary_dynamic_ar_hmc_evaluate <- function(q, cx, h) {
	x <- q[seq_len(cx$size)]
	l <- q[cx$size + 1L]
	alpha <- q[cx$size + 2L]
	r <- h$limit * tanh(alpha)
	w <- exp(2 * l)
	v <- w * (1 - r * r)
	# zero prior density is an ordinary rejection, not a numerical failure
	if (!is.finite(v) || v <= 0 || !is.finite(alpha))
		return(list(logp = -Inf, gradient = NULL, bound = TRUE))
	if (abs(r) >= h$limit || v > h$variance_max || w > h$stationary_max)
		return(list(logp = -Inf, gradient = NULL, bound = TRUE))
	ah <- .lame_binary_dynamic_ar_hmc_actual(h, l, alpha, cx)
	z <- .lame_binary_dynamic_pair_hmc_evaluate(x, cx, ah, NULL)
	if (!is.finite(z$logp))
		return(z)
	uv <- c(cx$indices$U, cx$indices$V)
	gp <- -2 * h$shape + 2 * h$scale / v
	gl <- sum(x[uv] * (z$gradient[uv] + x[uv])) + gp
	gr <- 0
	for (ix in list(cx$indices$U, cx$indices$V)) {
		standard <- matrix(x[ix], ncol = cx$nt)
		path <- standard %*% ah$uv_root
		gq <- matrix(z$gradient[ix] + x[ix], ncol = cx$nt)
		physical_gradient <- t(backsolve(ah$uv_root, t(gq)))
		gr <- gr + sum(
			physical_gradient * .lame_binary_dynamic_ar_hmc_path_derivative(
				path,
				r
			)
		)
	}
	dr <- h$limit * (1 - (r / h$limit)^2)
	ga <- (gr - gp * r / (1 - r * r) - (r - h$prior_mean) / (h$prior_sd^2)) * dr - 2 * r / h$limit
	# log-scale Jacobian and transformed persistence; whitening cancels the path determinants
	z$logp <- z$logp - h$shape * log(v) - h$scale / v + dnorm(
		r,
		h$prior_mean,
		h$prior_sd,
		log = TRUE
	) + log(dr)
	z$gradient <- c(z$gradient, gl, ga)
	z
}

.lame_binary_dynamic_ar_hmc_mass <- function(cx, h) {
	C <- h$uv_root
	r <- h$reference_rho
	dC <- matrix(0, nrow(C), ncol(C))
	if (cx$nt > 1L)
		for (t in 2:cx$nt) {
			dC[, t] <- C[, t - 1L] + r * dC[, t - 1L]
			dC[t, t] <- dC[t, t] - r / (1 - r * r) * C[t, t]
		}
	w <- h$uv_variance
	dr <- h$limit * (1 - (r / h$limit)^2)
	derivative_variance <- 2 * cx$rank * (w * colSums(dC * dC) + colSums(C * dC)^2)
	ma <- 1 + (dr / h$prior_sd)^2 + 4 * h$shape * (r * dr / (1 - r * r))^2 + cx$weight * cx$n *
		(cx$n - 1) * sum(derivative_variance) * dr^2
	ml <- 4 * h$shape + 4 * cx$weight * cx$n * (cx$n - 1) * cx$nt * cx$rank * w * w
	c(.lame_binary_dynamic_pair_hmc_mass(cx, h), ml, ma)
}

.lame_binary_dynamic_ar_hmc_integrate <- function(q, p, cx, h, mass, step, L) {
	z <- .lame_binary_dynamic_ar_hmc_evaluate(q, cx, h)
	p <- p + step / 2 * z$gradient
	for (k in seq_len(L)) {
		q <- q + step * p / mass
		z <- tryCatch(.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h), error = function(e) NULL)
		# zero density rejects like a bound; NaN or an error is a numerical failure
		if (!is.null(z) && identical(z$logp, -Inf))
			return(list(failed = TRUE, bound = TRUE))
		if (is.null(z) || !is.finite(z$logp) || any(!is.finite(z$gradient)))
			return(list(failed = TRUE, bound = !is.null(z) && isTRUE(z$bound)))
		p <- p + step * (if (k == L)
				0.5
			else 1) * z$gradient
	}
	list(q = q, p = p, evaluation = z, failed = FALSE)
}

.lame_binary_dynamic_ar_hmc_step <- function(q, cx, h, step, L = sample(12:20, 1)) {
	mass <- .lame_binary_dynamic_ar_hmc_mass(cx, h)
	p <- rnorm(length(q)) * sqrt(mass)
	old <- tryCatch(.lame_binary_dynamic_ar_hmc_evaluate(q, cx, h), error = function(e) NULL)
	# the current state should always evaluate, so a failure here is numerical
	if (is.null(old) || !is.finite(old$logp))
		return(list(q = q, accepted = FALSE, prob = 0, bound = FALSE, nonfinite = TRUE))
	new <- .lame_binary_dynamic_ar_hmc_integrate(q, p, cx, h, mass, step, L)
	if (new$failed)
		return(list(q = q, accepted = FALSE, prob = 0, bound = isTRUE(new$bound), nonfinite = !isTRUE(new$bound)))
	lr <- .lame_hmc_log_acceptance(new$evaluation$logp, old$logp, new$p, p, mass)
	if (!is.finite(lr))
		return(list(q = q, accepted = FALSE, prob = 0, bound = FALSE, nonfinite = TRUE))
	prob <- min(1, exp(lr))
	accept <- is.finite(prob) && runif(1) < prob
	list(
		q = if (accept) new$q else q,
		accepted = accept,
		prob = if (is.finite(prob)) prob else 0,
		bound = FALSE,
		nonfinite = FALSE,
		eta = if (accept) new$evaluation$eta else NULL
	)
}
