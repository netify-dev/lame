#' Summary of a LAME object
#' 
#' Provides a comprehensive summary of a fitted LAME (Longitudinal Additive and
#' Multiplicative Effects) model, including parameter estimates, standard errors,
#' credible intervals, and model diagnostics.
#' 
#' @details
#' The summary includes:
#' \describe{
#'   \item{Regression coefficients}{Posterior means, posterior standard deviations,
#'         z-values, approximate p-values, and 95% credible intervals for dyadic,
#'         sender, and receiver covariates. Note: the z-values are computed as
#'         posterior mean / posterior SD, and the p-values are derived from a
#'         normal approximation. These are convenient screening statistics but
#'         are not formal frequentist test statistics. For rigorous inference,
#'         use the credible intervals or examine the full posterior via the
#'         BETA matrix directly.}
#'   \item{Variance components}{Estimates and standard errors for:
#'     \describe{
#'       \item{va}{Variance of additive sender/row effects}
#'       \item{cab}{Covariance between sender and receiver effects}
#'       \item{vb}{Variance of additive receiver/column effects}
#'       \item{rho}{Dyadic correlation (reciprocity)}
#'       \item{ve}{Residual variance}
#'     }}
#' }
#' 
#' @param object an object of class "lame", typically the result of fitting a 
#'        longitudinal AME model using the \code{lame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.lame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates}
#'   \item{n.periods}{Number of time periods in the longitudinal data}
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
#' @seealso \code{\link{lame}}, \code{\link{print.summary.lame}}
#' @method summary lame
#' @export
summary.lame <- function(object, ...) {
	fit <- object
	
	# regression coefficient summaries
	beta_mean <- apply(fit$BETA, 2, mean, na.rm = TRUE)
	beta_sd <- apply(fit$BETA, 2, sd, na.rm = TRUE)
	beta_z <- ifelse(beta_sd > 0, beta_mean / beta_sd, NA_real_)
	beta_p <- 2 * (1 - pnorm(abs(beta_z)))
	
	# 95% quantile-based credible interval
	beta_lower <- apply(fit$BETA, 2, quantile, probs = 0.025, na.rm = TRUE)
	beta_upper <- apply(fit$BETA, 2, quantile, probs = 0.975, na.rm = TRUE)

	# coefficient table
	beta_table <- cbind(
		Estimate = beta_mean,
		StdError = beta_sd,
		z_value = beta_z,
		p_value = beta_p,
		CI_lower = beta_lower,
		CI_upper = beta_upper
	)
	
	# variance component summaries
	vc_mean <- apply(fit$VC, 2, mean)
	vc_sd <- apply(fit$VC, 2, sd)
	
	vc_table <- cbind(
		Estimate = vc_mean,
		StdError = vc_sd
	)
	
	n_periods <- fit$n_time

	# assemble summary object
	sum_obj <- list(
		call = fit$call,
		beta = beta_table,
		variance = vc_table,
		n.periods = n_periods,
		symmetric = isTRUE(fit$symmetric),
		mode = fit$mode %||% "unipartite",
		family = fit$family,
		dynamic_uv = isTRUE(fit$dynamic_uv),
		dynamic_ab = isTRUE(fit$dynamic_ab),
		rho_ab = if(isTRUE(fit$dynamic_ab) && !is.null(fit$rho_ab)) mean(fit$rho_ab) else NULL,
		rho_uv = if(isTRUE(fit$dynamic_uv) && !is.null(fit$rho_uv)) mean(fit$rho_uv) else NULL
	)

	class(sum_obj) <- "summary.lame"
	return(sum_obj)
}

#' Print method for summary.lame objects
#' 
#' Prints a formatted summary of a LAME model fit
#' 
#' @param x a summary.lame object
#' @param digits number of digits to display (default: 3)
#' @param ... additional arguments (not used)
#' @return the summary.lame object invisibly
#' @export
print.summary.lame <- function(x, digits = 3, ...) {
	cat("\n=== Longitudinal AME Model Summary ===\n")
	cat("\nCall:\n")
	print(x$call)
	
	if (!is.null(x$n.periods)) {
		cat("\nTime periods:", x$n.periods, "\n")
	}
	if (!is.null(x$family)) {
		cat("Family:", x$family, "\n")
	}
	if (!is.null(x$mode)) {
		cat("Mode:", x$mode, "\n")
	}
	if (isTRUE(x$dynamic_uv)) {
		cat("Dynamic latent positions: enabled")
		if (!is.null(x$rho_uv)) cat(" (rho_uv =", round(x$rho_uv, 3), ")")
		cat("\n")
	}
	if (isTRUE(x$dynamic_ab)) {
		cat("Dynamic additive effects: enabled")
		if (!is.null(x$rho_ab)) cat(" (rho_ab =", round(x$rho_ab, 3), ")")
		cat("\n")
	}
	
	cat("\nRegression coefficients:\n")
	cat("------------------------\n")
	coef_table <- x$beta
	coef_table <- round(coef_table, digits)
	stars <- symnum(x$beta[, "p_value"], 
									corr = FALSE,
									na = FALSE,
									cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
									symbols = c("***", "**", "*", ".", " "))
	coef_table <- cbind(coef_table, " " = stars)
	print(coef_table, quote = FALSE, right = TRUE)
	cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
	cat("Note: p-values are approximate (posterior mean / SD); use credible intervals for inference.\n")

	cat("\nVariance components:\n")
	cat("-------------------\n")
	var_table <- round(x$variance, digits)
	print(var_table, quote = FALSE, right = TRUE)

	if(isTRUE(x$symmetric)) {
		cat("  (va = nodal variance, ve = residual variance)\n")
	} else {
		cat("  (va = sender, cab = sender-receiver covariance, vb = receiver,\n")
		cat("   rho = dyadic correlation, ve = residual variance)\n")
	}
	if(!is.null(x$mode) && x$mode == "bipartite") {
		cat("  Note: bipartite model (rho fixed to 0, cab fixed to 0)\n")
	}

	invisible(x)
}