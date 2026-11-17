# covariance of the coefficients of an ALS fit: the Laplace covariance, or that of an
# attached bootstrap

#' Covariance of the coefficients of a fast AME fit
#'
#' @description
#' Returns the covariance matrix of every regression coefficient of an
#' \code{\link{ame_als}} or \code{\link{lame_als}} fit: the intercept, the
#' dyadic coefficients and the node-covariate coefficients.
#'
#' @details
#' Every fit the fast route returns is an empirical-Bayes fit (see
#' \code{\link{ame_als}}), and
#' its covariance is the Laplace covariance: the inverse of the penalised
#' Fisher information of all the parameters (the coefficients, the actor
#' effects, the node-covariate coefficients and the latent factors) at the
#' final weights (for a normal fit the exact information of the Gaussian cell,
#' whose Laplace approximation is exact for the coefficients and the effects),
#' with the actor effects' covariance \code{fit$eb$sab} and the factor scales
#' held at their estimates. It propagates the uncertainty of the actor effects
#' and the factors. The covariance conditional on the factors is
#' \code{fit$eb$vcov_cond}. A fit an earlier release made with a setting
#' this release no longer has stops with an error of class
#' \code{lame_als_removed_setting}; refit it with \code{\link{ame_als}} or
#' \code{\link{lame_als}}.
#'
#' A bootstrap attached to the fit (\code{bootstrap = N}) takes precedence: its
#' covariance is that of the replicates, and a coefficient the bootstrap gives
#' no standard error for (see \code{\link{ame_als_bootstrap}}) keeps its
#' analytic variance, with zero covariance with the others
#' (\code{fit$bootstrap$se_source} names it).
#'
#' @param object an \code{ame_als} fit.
#' @param ... ignored.
#' @return A covariance matrix with matching row/column names, covering every
#'   coefficient of \code{coef()}.
#' @method vcov ame_als
#' @seealso \code{\link{ame_als_bootstrap}} for bootstrap uncertainty.
#' @export
vcov.ame_als <- function(object, ...) {
	# an attached bootstrap supplies the covariance, so vcov(), summary() and confint()
	# agree on one fit (cov() of the replicates reproduces $bootstrap$se)
	boot_coefs <- object$bootstrap[["coefs", exact = TRUE]]
	if (is.matrix(boot_coefs)) {
		# failed replicates are NA rows, dropped as the se and interval code drops them
		boot_coefs <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]
	}
	if (is.matrix(boot_coefs) && nrow(boot_coefs) > 1L) {
		# a parameter the bootstrap returns no standard error for takes its analytic covariance
		return(vcov.boot_ame(object$bootstrap))
	}
	# empirical-Bayes fit: the Laplace covariance of the penalised Hessian
	.ae_fit_estimator(object)
	object$eb$vcov
}

#' Confidence intervals for a fast AME fit
#'
#' @description
#' Wald confidence intervals for every regression coefficient of an
#' \code{\link{ame_als}} fit, from its covariance (\code{\link{vcov.ame_als}}).
#'
#' @details
#' The intervals use the Laplace covariance of the empirical-Bayes fit, which
#' propagates the uncertainty of the actor effects and the latent factors and
#' covers the node-covariate coefficients too. A fit carrying a bootstrap
#' returns the bootstrap intervals of \code{\link{confint.boot_ame}} instead,
#' with analytic Wald intervals for the quantities the bootstrap gives none.
#'
#' @param object an \code{ame_als} fit.
#' @param parm character vector of parameter names, or integer indices; if
#'   \code{NULL} (default) all available coefficients are returned.
#' @param level confidence level (default 0.95).
#' @param ... passed to \code{\link{vcov.ame_als}}.
#' @return A matrix with one row per coefficient and lower/upper bound columns.
#' @method confint ame_als
#' @seealso \code{\link{ame_als_bootstrap}} for bootstrap intervals on all
#'   parameters.
#' @importFrom stats vcov
#' @export
confint.ame_als <- function(object, parm = NULL, level = 0.95, ...) {
	if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number in (0, 1).")
	}
	# a fit with a bootstrap returns its fully propagated intervals
	if (!is.null(object$bootstrap)) {
		return(stats::confint(object$bootstrap, parm = parm, level = level, ...))
	}
	V <- stats::vcov(object, ...)
	est <- object$coefficients[rownames(V)]
	se <- sqrt(diag(V))
	a <- (1 - level) / 2
	z <- stats::qnorm(c(a, 1 - a))
	ci <- cbind(est + z[1] * se, est + z[2] * se)
	rownames(ci) <- rownames(V)
	colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
	if (!is.null(parm)) {
		if (is.character(parm)) {
			bad <- setdiff(parm, rownames(ci))
			if (length(bad)) {
				cli::cli_abort(c(
					"Unknown parameter name in {.arg parm}: {.val {bad}}.",
					"i" = "Available: {.val {rownames(ci)}}."
				))
			}
		} else if (is.numeric(parm) &&
				(any(parm < 1) || any(parm > nrow(ci)))) {
			cli::cli_abort("{.arg parm} index is out of range (1:{nrow(ci)}).")
		}
		ci <- ci[parm, , drop = FALSE]
	}
	ci
}
