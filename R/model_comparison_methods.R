#' @importFrom stats nobs logLik formula
NULL

####
#' Number of observed dyads in an AME / LAME fit
#'
#' Number of finite (non-missing) cells in \code{Y}, summed over time slices
#' for longitudinal fits. Useful as a denominator for sample-size reporting
#' and as the basis for AIC/BIC \emph{if} you also have a log-likelihood
#' (which \code{ame()} / \code{lame()} do not directly expose; see
#' \code{\link{logLik.ame}}).
#'
#' @param object an \code{ame} or \code{lame} fit.
#' @param ... ignored.
#' @return Integer: number of observed dyads in \code{Y}.
#' @method nobs ame
#' @export
nobs.ame <- function(object, ...) {
	Y <- object$Y
	if (is.null(Y)) Y <- object$Y_obs
	if (is.null(Y)) {
		cli::cli_abort("Cannot compute {.fn nobs}: {.field Y} not stored on the fit.")
	}
	if (is.list(Y)) {
		sum(vapply(Y, function(y) sum(is.finite(y)), integer(1)))
	} else {
		sum(is.finite(Y))
	}
}

#' @rdname nobs.ame
#' @method nobs lame
#' @export
nobs.lame <- nobs.ame

####
#' Log-likelihood is not directly exposed for ame() / lame() fits
#'
#' \code{ame()} and \code{lame()} produce a posterior sample, not a maximum-
#' likelihood point. A pointwise log-likelihood is computable from the
#' posterior draws but is not stored on the fit object, so \code{logLik()}
#' (and the AIC / BIC generics that dispatch through it) error out
#' informatively rather than return a misleading number.
#'
#' For Bayesian model comparison use posterior-predictive checks via
#' \code{\link{gof}} / \code{\link{gof_plot}}, or compute WAIC / LOO yourself
#' from the per-draw log-likelihoods (e.g. via the \pkg{loo} package on the
#' BETA / VC chains).
#'
#' @param object an \code{ame} or \code{lame} fit.
#' @param ... ignored.
#' @return Never returns; raises an error.
#' @method logLik ame
#' @export
logLik.ame <- function(object, ...) {
	cli::cli_abort(c(
		"A point log-likelihood is not stored on an {.fn ame} / {.fn lame} fit.",
		"i" = "These are Bayesian MCMC fits, so {.fn AIC} / {.fn BIC} dispatch through {.fn logLik} would be misleading.",
		"i" = "For model comparison use {.fn gof} / {.fn gof_plot} (posterior predictive checks), or compute WAIC / LOO from the posterior draws with the {.pkg loo} package."))
}

#' @rdname logLik.ame
#' @method logLik lame
#' @export
logLik.lame <- logLik.ame

####
#' formula() is not defined for an ame() / lame() fit
#'
#' \code{ame()} and \code{lame()} use named arguments (\code{Xrow}, \code{Xcol},
#' \code{Xdyad}, \code{symmetric}, etc.) rather than a formula interface. The
#' \code{summary()} output renders a pseudo-formula for diagnostic display, but
#' it is not a valid R formula and cannot be passed back to a fitting function.
#'
#' @param x an \code{ame} or \code{lame} fit.
#' @param ... ignored.
#' @return Never returns; raises an error.
#' @method formula ame
#' @export
formula.ame <- function(x, ...) {
	cli::cli_abort(c(
		"{.fn ame} / {.fn lame} do not use a formula interface.",
		"i" = "Covariates are passed via {.arg Xrow}, {.arg Xcol}, and {.arg Xdyad}.",
		"i" = "The {.field Call} line printed by {.fn summary} is a diagnostic, not a re-usable formula."))
}

#' @rdname formula.ame
#' @method formula lame
#' @export
formula.lame <- formula.ame
