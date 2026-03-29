#' Extract model coefficients from AME model
#'
#' Returns posterior means of regression coefficients from a fitted AME or
#' LAME model. The coefficients correspond to the columns of the \code{BETA}
#' matrix stored in the model object, which records one draw per post-burn-in
#' MCMC iteration. The posterior mean is computed as \code{colMeans(BETA)}.
#'
#' For binary models, these are on the probit (latent) scale. Use
#' \code{\link{predict.ame}} with \code{type = "response"} to get predicted
#' probabilities.
#'
#' @param object Fitted AME model (class \code{"ame"} or \code{"lame"}).
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector of posterior mean coefficients.
#'
#' @seealso \code{\link{vcov.ame}} for the posterior covariance matrix,
#'   \code{\link{confint.ame}} for credible intervals,
#'   \code{\link{summary.ame}} for a full summary table
#'
#' @method coef ame
#' @export
coef.ame <- function(object, ...) {
	if(is.null(object$BETA) || nrow(object$BETA) == 0) {
		return(numeric(0))
	}
	out <- colMeans(object$BETA)
	nms <- colnames(object$BETA)
	if(!is.null(nms)) names(out) <- nms
	out
}

#' @rdname coef.ame
#' @method coef lame
#' @export
coef.lame <- function(object, ...) {
	coef.ame(object, ...)
}

#' Posterior covariance of AME model coefficients
#'
#' Returns the posterior covariance matrix of regression coefficients.
#'
#' @param object fitted AME model (class "ame")
#' @param ... additional arguments (ignored)
#'
#' @return p x p covariance matrix of posterior BETA draws
#'
#' @method vcov ame
#' @export
vcov.ame <- function(object, ...) {
	if(is.null(object$BETA) || nrow(object$BETA) < 2) {
		p <- if(!is.null(object$BETA)) ncol(object$BETA) else 0
		return(matrix(NA, p, p))
	}
	out <- cov(object$BETA)
	nms <- colnames(object$BETA)
	if(!is.null(nms)) {
		rownames(out) <- colnames(out) <- nms
	}
	out
}

#' @rdname vcov.ame
#' @method vcov lame
#' @export
vcov.lame <- function(object, ...) {
	vcov.ame(object, ...)
}

#' Bayesian credible intervals for AME model coefficients
#'
#' Returns posterior quantile-based credible intervals for regression
#' coefficients. These are Bayesian credible intervals, not frequentist
#' confidence intervals.
#'
#' @param object fitted AME model (class "ame")
#' @param parm character vector of parameter names, or numeric indices.
#'   If NULL (default), intervals for all parameters are returned.
#' @param level credible level (default 0.95)
#' @param ... additional arguments (ignored)
#'
#' @return matrix with columns for lower and upper bounds
#'
#' @method confint ame
#' @export
confint.ame <- function(object, parm = NULL, level = 0.95, ...) {
	if(is.null(object$BETA) || nrow(object$BETA) == 0) {
		return(matrix(nrow = 0, ncol = 2,
									dimnames = list(NULL, c("lower", "upper"))))
	}

	alpha <- (1 - level) / 2
	probs <- c(alpha, 1 - alpha)

	ci <- apply(object$BETA, 2, quantile, probs = probs)
	ci <- t(ci)

	nms <- colnames(object$BETA)
	if(!is.null(nms)) rownames(ci) <- nms
	pct <- paste0(round(probs * 100, 1), "%")
	colnames(ci) <- pct

	# subset if parm is specified

	if(!is.null(parm)) {
		if(is.character(parm)) {
			ci <- ci[parm, , drop = FALSE]
		} else if(is.numeric(parm)) {
			ci <- ci[parm, , drop = FALSE]
		}
	}

	ci
}

#' @rdname confint.ame
#' @method confint lame
#' @export
confint.lame <- function(object, parm = NULL, level = 0.95, ...) {
	confint.ame(object, parm = parm, level = level, ...)
}
