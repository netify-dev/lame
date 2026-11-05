# a factor column has one likelihood eigenvector along its partner and a repeated perpendicular
# eigenvalue; symmetric and antisymmetric residual scores stay separate so the latter never leaks
# into the parallel score
.lame_factor_common_contrast <- function(E, rho, s2 = 1) {
	dims <- dim(E)
	if (length(dims) != 3L || dims[1L] != dims[2L] || any(dims <= 0L) ||
			any(!is.finite(E)) || !(length(rho) %in% c(1L, dims[3L])) ||
			any(!is.finite(rho)) || any(abs(rho) >= 1) ||
			length(s2) != 1L || !is.finite(s2) || s2 <= 0)
		stop('Factor Gaussian residuals, correlations and variance must be valid.')
	rho <- rep_len(rho, dims[3L])
	plus <- (1 / s2) / (1 + rho); minus <- (1 / s2) / (1 - rho)
	total_plus <- sum(plus); total_minus <- sum(minus)
	if (any(!is.finite(c(total_plus, total_minus))) || total_plus <= 0 || total_minus <= 0)
		stop('Factor Gaussian likelihood weights cannot be represented.')
	symmetric <- antisymmetric <- matrix(0, dims[1L], dims[1L])
	for (tt in seq_len(dims[3L])) {
		E_t <- E[, , tt]
		symmetric <- symmetric + (plus[tt] / total_plus) * (E_t / 2 + t(E_t) / 2)
		antisymmetric <- antisymmetric + (minus[tt] / total_minus) * (E_t / 2 - t(E_t) / 2)
	}
	list(
		symmetric = symmetric,
		antisymmetric = antisymmetric,
		plus = total_plus,
		minus = total_minus
	)
}

.lame_factor_reflect <- function(x, reflector) {
	if (is.null(reflector)) return(x)
	x - 2 * reflector * sum(reflector * x)
}

.lame_factor_gaussian_system <- function(
	partner,
	signal,
	prior_mean,
	prior_precision,
	transpose = FALSE
) {
	n <- length(partner)
	if (!n || length(prior_mean) != n || any(!is.finite(c(partner, prior_mean))) ||
			length(prior_precision) != 1L || !is.finite(prior_precision) || prior_precision <= 0 ||
			!identical(dim(signal$symmetric), c(n, n)) || !identical(dim(signal$antisymmetric), c(n, n)) ||
			any(!is.finite(c(signal$symmetric, signal$antisymmetric, signal$plus, signal$minus))) ||
			signal$plus <= 0 || signal$minus <= 0)
		stop('Factor Gaussian signal, partner and conditional prior must be valid.')
	scale <- max(abs(partner))
	if (scale == 0) return(list(
		mean = as.numeric(prior_mean),
		reflector = NULL,
		inverse_roots = rep(1 / sqrt(prior_precision), n)
	))
	scaled <- partner / scale; norm_scaled <- sqrt(sum(scaled^2))
	direction <- scaled / norm_scaled
	reflector <- direction
	reflector[1L] <- reflector[1L] + if (direction[1L] >= 0) 1 else -1
	reflector <- reflector / sqrt(sum(reflector^2))
	plus <- signal$plus; minus <- signal$minus
	self <- plus / 2 + minus / 2
	weighted_parallel <- (sqrt(plus) * scale) * norm_scaled
	weighted_perpendicular <- (sqrt(self) * scale) * norm_scaled
	parallel <- prior_precision + weighted_parallel^2
	perpendicular <- prior_precision + weighted_perpendicular^2
	if (any(!is.finite(c(parallel, perpendicular))) || min(parallel, perpendicular) <= 0)
		stop('Factor Gaussian eigenprecision cannot be represented.')
	prior_coordinates <- .lame_factor_reflect(prior_mean, reflector)
	symmetric_coordinates <- .lame_factor_reflect(drop(signal$symmetric %*% direction), reflector)
	antisymmetric_coordinates <- .lame_factor_reflect(drop(signal$antisymmetric %*% direction), reflector)
	if (isTRUE(transpose)) antisymmetric_coordinates <- -antisymmetric_coordinates
	# ||partner|| * weight / precision, computed without forming a possibly enormous score
	coefficient <- function(weight, precision)
		(((sqrt(weight) * scale) * norm_scaled) / sqrt(precision)) *
		(sqrt(weight) / sqrt(precision))
	coordinates <- numeric(n)
	coordinates[1L] <- coefficient(plus, parallel) * symmetric_coordinates[1L] +
		(prior_precision / parallel) * prior_coordinates[1L]
	if (n > 1L) {
		at <- seq.int(2L, n)
		coordinates[at] <- coefficient(plus, perpendicular) * symmetric_coordinates[at] +
			coefficient(minus, perpendicular) * antisymmetric_coordinates[at] +
			(prior_precision / perpendicular) * prior_coordinates[at]
	}
	mean <- .lame_factor_reflect(coordinates, reflector)
	inverse_roots <- c(1 / sqrt(parallel), rep(1 / sqrt(perpendicular), n - 1L))
	if (any(!is.finite(c(mean, inverse_roots))) || any(inverse_roots <= 0))
		stop('Factor Gaussian moments cannot be represented.')
	list(mean = mean, reflector = reflector, inverse_roots = inverse_roots)
}

.lame_factor_gaussian_innovation <- function(system, e) {
	if (length(e) != length(system$mean) || any(!is.finite(e)))
		stop('Factor Gaussian innovation must have the expected finite dimension.')
	.lame_factor_reflect(system$inverse_roots * e, system$reflector)
}

.lame_factor_gaussian_draw <- function(
	partner,
	signal,
	prior_mean,
	prior_precision,
	transpose = FALSE
) {
	system <- .lame_factor_gaussian_system(partner, signal, prior_mean, prior_precision, transpose)
	draw <- system$mean + .lame_factor_gaussian_innovation(system, stats::rnorm(length(partner)))
	if (any(!is.finite(draw)))
		stop('Factor Gaussian draw is not finite.')
	draw
}

# precision of the stacked factor prior covariance: equilibrated Cholesky inverse, with an
# eigenvalue-floored inverse when the covariance is not numerically positive definite
.lame_factor_prior_precision <- function(covariance) {
	d <- nrow(covariance)
	if (!is.matrix(covariance) || d < 1L || ncol(covariance) != d ||
			any(!is.finite(covariance)) || any(diag(covariance) <= 0))
		cli::cli_abort("Factor prior covariance must be finite with positive diagonal.")
	scales <- sqrt(diag(covariance))
	balanced <- covariance / outer(scales, scales)
	balanced <- balanced / 2 + t(balanced) / 2
	root <- tryCatch(chol(balanced), error = function(e) NULL)
	precision <- if (!is.null(root)) chol2inv(root) else {
		e <- eigen(balanced, symmetric = TRUE)
		ev <- pmax(e$values, max(e$values, 1) * sqrt(.Machine$double.eps))
		e$vectors %*% (t(e$vectors) / ev)
	}
	precision <- (precision + t(precision)) / 2 / outer(scales, scales)
	if (any(!is.finite(precision)))
		cli::cli_abort("Factor prior precision is not finite.")
	precision
}

.lame_factor_validate_coordinates <- function(E, U, V, rho, s2) {
	dims <- dim(E)
	if (length(dims) != 3L || dims[1L] != dims[2L] || any(dims <= 0L) ||
			!is.matrix(U) || !is.matrix(V) || nrow(U) != dims[1L] ||
			!identical(dim(U), dim(V)) || any(!is.finite(c(E, U, V))))
		cli::cli_abort("Factor conditionals require finite complete residuals and matching coordinates.")
	.lame_factor_common_contrast(E, rho, s2)
	invisible(TRUE)
}
