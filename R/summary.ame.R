#' Summary of an AME object
#' 
#' Provides a comprehensive summary of a fitted AME (Additive and Multiplicative 
#' Effects) model, including parameter estimates, standard errors, confidence 
#' intervals, and model diagnostics.
#' 
#' @details
#' The summary includes:
#' \describe{
#'   \item{Regression coefficients}{Point estimates, standard errors, z-values,
#'         p-values, and 95% confidence intervals for dyadic, sender, and 
#'         receiver covariates}
#'   \item{Variance components}{Estimates and standard errors for:
#'     \describe{
#'       \item{va}{Variance of additive sender/row effects (asymmetric networks)}
#'       \item{cab}{Covariance between sender and receiver effects}
#'       \item{vb}{Variance of additive receiver/column effects (asymmetric networks)}
#'       \item{rho}{Dyadic correlation (reciprocity in directed networks)}
#'       \item{ve}{Residual variance}
#'     }
#'     For symmetric networks, only va and ve are estimated.}
#' }
#' 
#' @param object an object of class "ame", typically the result of fitting an 
#'        AME model using the \code{ame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.ame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{ame}}, \code{\link{print.summary.ame}}
#' @method summary ame
#' @export
summary.ame <- function(object, ...) {
	fit <- object
	
	# regression coefficient summaries
	beta_mean <- apply(fit$BETA, 2, mean)
	beta_sd <- apply(fit$BETA, 2, sd)
	beta_z <- beta_mean / beta_sd
	beta_p <- 2 * (1 - pnorm(abs(beta_z)))
	
	# 95% credible interval
	z_val <- 1.96
	beta_lower <- beta_mean - z_val * beta_sd
	beta_upper <- beta_mean + z_val * beta_sd
	
	# coefficient names
	if(!is.null(fit$X_names)) {
		beta_names <- fit$X_names
	} else if(!is.null(dimnames(fit$X)[[3]])) {
		beta_names <- dimnames(fit$X)[[3]]
	} else {
		beta_names <- paste0("beta", 0:(ncol(fit$BETA)-1))
	}
	
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
	
	# label rows
	if(length(beta_names) == nrow(beta_table)) {
		rownames(beta_table) <- beta_names
	}
	
	# format call
	display_call <- tryCatch({
		format_ame_call(fit)
	}, error = function(e) {
		fit$call
	})
	
	# assemble summary object
	sum_obj <- list(
		call = display_call,
		beta = beta_table,
		variance = vc_table,
		symmetric = isTRUE(fit$symmetric),
		mode = fit$mode %||% "unipartite",
		family = fit$family,
		dynamic_uv = isTRUE(fit$dynamic_uv),
		dynamic_ab = isTRUE(fit$dynamic_ab)
	)

	class(sum_obj) <- "summary.ame"
	return(sum_obj)
}

#' Print method for summary.ame objects
#' 
#' Prints a formatted summary of an AME model fit
#' 
#' @param x a summary.ame object
#' @param digits number of digits to display (default: 3)
#' @param ... additional arguments (not used)
#' @return the summary.ame object invisibly
#' @export
print.summary.ame <- function(x, digits = 3, ...) {
	cat("\n=== AME Model Summary ===\n")
	cat("\nCall:\n")
	print(x$call)
	
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