#' Summary of an AME object
#' 
#' Provides a comprehensive summary of a fitted AME (Additive and Multiplicative
#' Effects) model, including parameter estimates, standard errors, credible
#' intervals, and model diagnostics.
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
	beta_z <- ifelse(beta_sd > 0, beta_mean / beta_sd, NA_real_)
	beta_p <- 2 * (1 - pnorm(abs(beta_z)))
	
	# 95% quantile-based credible interval
	beta_lower <- apply(fit$BETA, 2, quantile, probs = 0.025)
	beta_upper <- apply(fit$BETA, 2, quantile, probs = 0.975)

	# coefficient names
	# the bipartite path leaves fit$X_names and dimnames(fit$X)[[3]] NULL even
	# when fit$BETA has proper colnames, so fall back to colnames(fit$BETA)
	# before the generic beta0..betaN labels.
	if(!is.null(fit$X_names)) {
		beta_names <- fit$X_names
	} else if(!is.null(dimnames(fit$X)[[3]])) {
		beta_names <- dimnames(fit$X)[[3]]
	} else if(!is.null(colnames(fit$BETA))) {
		beta_names <- colnames(fit$BETA)
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
	
	# surface multi-chain diagnostics if the fit carries them so a brms-style
	# user can read Rhat/ESS right out of summary()
	diag_block <- NULL
	if (!is.null(fit$diagnostics) &&
	    (!is.null(fit$diagnostics$rhat) || !is.null(fit$diagnostics$ess))) {
		n_chains_val <- fit$diagnostics$n_chains %||% NA_integer_
		rhat_all <- fit$diagnostics$rhat
		ess_all  <- fit$diagnostics$ess
		# beta_table rownames are the BETA columns; fish out the matching subset
		if (!is.null(rhat_all)) {
			r_beta <- rhat_all[intersect(names(rhat_all), rownames(beta_table))]
			# pad to match rows so cbind aligns
			beta_table <- cbind(beta_table,
			                    Rhat = unname(rhat_all[match(rownames(beta_table),
			                                                 names(rhat_all))]))
		}
		if (!is.null(ess_all)) {
			beta_table <- cbind(beta_table,
			                    ESS = unname(ess_all[match(rownames(beta_table),
			                                               names(ess_all))]))
		}
		# vc table
		if (!is.null(rhat_all)) {
			vc_table <- cbind(vc_table,
			                  Rhat = unname(rhat_all[match(rownames(vc_table),
			                                               names(rhat_all))]))
		}
		if (!is.null(ess_all)) {
			vc_table <- cbind(vc_table,
			                  ESS = unname(ess_all[match(rownames(vc_table),
			                                             names(ess_all))]))
		}
		max_rhat <- if (length(rhat_all) > 0) max(rhat_all, na.rm = TRUE) else NA_real_
		min_ess  <- if (length(ess_all)  > 0) min(ess_all,  na.rm = TRUE) else NA_real_
		diag_block <- list(n_chains = n_chains_val,
		                   max_rhat = max_rhat,
		                   min_ess  = min_ess)
	}

	# assemble summary object. `coefficients` is an alias for `beta` to match
	# the `summary(lm())$coefficients` convention that base-R users reach for
	sum_obj <- list(
		call = display_call,
		beta = beta_table,
		coefficients = beta_table,
		variance = vc_table,
		symmetric = isTRUE(fit$symmetric),
		mode = fit$mode %||% "unipartite",
		family = fit$family,
		dynamic_uv = isTRUE(fit$dynamic_uv),
		dynamic_ab = isTRUE(fit$dynamic_ab),
		diagnostics = diag_block
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
	cat("Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.\n")

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