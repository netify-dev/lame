#' Summary of a LAME object
#' 
#' Summarizes a fitted LAME (Longitudinal Additive and Multiplicative Effects)
#' model, including parameter estimates, standard errors, credible intervals,
#' and model diagnostics.
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
#'   \item{Dynamic coefficients per period}{Only printed when the fit was
#'         produced with \code{dynamic_beta} on at least one coefficient. The
#'         table has one row per coefficient with columns:
#'         \describe{
#'           \item{Mean}{average of the per-period posterior means across t}
#'           \item{Min, Max}{smallest and largest per-period posterior mean}
#'           \item{Drift_pct}{\code{100 * (Max - Min) / |Mean|}; a coarse
#'                 measure of how much the coefficient moves across periods
#'                 relative to its average level}
#'           \item{Dynamic}{\code{"Y"} if the coefficient was flagged as
#'                 dynamic, \code{"N"} if it was held static}
#'         }
#'         The block also prints the per-block AR(1) hyperparameters
#'         (\code{rho_beta = ...}). For per-period credible intervals use
#'         \code{\link{confint.lame}}.}
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

	# rendered model formula including multiplicative term and family.
	# falls back to the raw match.call if format_ame_call fails.
	display_call <- tryCatch(format_ame_call(fit),
	                         error = function(e) fit$call)

	# regression coefficient summaries. for 3-d dynamic_beta beta, the main
	# table reports the overall (time-collapsed) posterior summary; the
	# per-period table is attached as `beta_dynamic` for the print method.
	beta_is_dyn <- length(dim(fit$BETA)) == 3L
	if (beta_is_dyn) {
		# overall: collapse across iterations and periods
		BETA_overall <- matrix(fit$BETA, nrow = dim(fit$BETA)[1],
		                       ncol = dim(fit$BETA)[2] * dim(fit$BETA)[3])
		# but the overall mean is meaningless for time-varying coefficients
		# (the time-average could even pass through zero while every period is
		# nonzero). instead, report the overall time-average per-coefficient
		# computed as mean over iterations and periods, plus its drift range.
		beta_mean_overall <- apply(fit$BETA, 2, mean, na.rm = TRUE)
		beta_sd_overall   <- apply(apply(fit$BETA, c(1, 2), mean), 2, sd, na.rm = TRUE)
		beta_min   <- apply(apply(fit$BETA, c(2, 3), mean), 1, min)
		beta_max   <- apply(apply(fit$BETA, c(2, 3), mean), 1, max)
		drift_pct  <- ifelse(abs(beta_mean_overall) > 1e-8,
		                     100 * (beta_max - beta_min) / abs(beta_mean_overall),
		                     NA_real_)
		beta_z <- ifelse(beta_sd_overall > 0, beta_mean_overall / beta_sd_overall, NA_real_)
		beta_p <- 2 * (1 - pnorm(abs(beta_z)))
		beta_lower <- apply(apply(fit$BETA, c(1, 2), mean), 2, quantile,
		                    probs = 0.025, na.rm = TRUE)
		beta_upper <- apply(apply(fit$BETA, c(1, 2), mean), 2, quantile,
		                    probs = 0.975, na.rm = TRUE)
		beta_table <- cbind(
			Estimate = beta_mean_overall,
			StdError = beta_sd_overall,
			z_value  = beta_z,
			p_value  = beta_p,
			CI_lower = beta_lower,
			CI_upper = beta_upper
		)
		# per-period table (mean | sd | min | max | %drift)
		beta_per_t_mean <- apply(fit$BETA, c(2, 3), mean)
		beta_per_t_sd   <- apply(fit$BETA, c(2, 3), sd)
	} else {
		beta_mean <- apply(fit$BETA, 2, mean, na.rm = TRUE)
		beta_sd <- apply(fit$BETA, 2, sd, na.rm = TRUE)
		beta_z <- ifelse(beta_sd > 0, beta_mean / beta_sd, NA_real_)
		beta_p <- 2 * (1 - pnorm(abs(beta_z)))
		# 95% quantile-based credible interval
		beta_lower <- apply(fit$BETA, 2, quantile, probs = 0.025, na.rm = TRUE)
		beta_upper <- apply(fit$BETA, 2, quantile, probs = 0.975, na.rm = TRUE)
		beta_table <- cbind(
			Estimate = beta_mean,
			StdError = beta_sd,
			z_value = beta_z,
			p_value = beta_p,
			CI_lower = beta_lower,
			CI_upper = beta_upper
		)
		beta_per_t_mean <- NULL
		beta_per_t_sd   <- NULL
		beta_min <- NULL; beta_max <- NULL; drift_pct <- NULL
	}
	
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
		call = display_call,
		beta = beta_table,
		variance = vc_table,
		n.periods = n_periods,
		symmetric = isTRUE(fit$symmetric),
		mode = fit$mode %||% "unipartite",
		family = fit$family,
		dynamic_uv = isTRUE(fit$dynamic_uv),
		dynamic_ab = isTRUE(fit$dynamic_ab),
		dynamic_beta = isTRUE(fit$dynamic_beta),
		dynamic_beta_kind = fit$dynamic_beta_kind,
		beta_dynamic_mask = fit$beta_dynamic_mask,
		beta_dynamic_per_t = beta_per_t_mean,
		beta_dynamic_per_t_sd = beta_per_t_sd,
		beta_dynamic_min = beta_min,
		beta_dynamic_max = beta_max,
		beta_dynamic_drift_pct = drift_pct,
		# report the posterior mean of the rho_beta draws (fit$rho_beta holds
		# the last iteration); fall back to fit$rho_beta only when the draws
		# are unavailable.
		rho_beta = if (isTRUE(fit$dynamic_beta)) {
			RB <- fit$RHO_BETA
			if (is.matrix(RB) && nrow(RB) > 0L) colMeans(RB, na.rm = TRUE)
			else if (!is.null(RB)) mean(RB, na.rm = TRUE)
			else fit$rho_beta
		} else NULL,
		sigma_beta = if (isTRUE(fit$dynamic_beta)) fit$sigma_beta else NULL,
		# per-block stationarity diagnostic: q05 of the rho_beta posterior
		# and its iqr. the print method warns when q05 >= 0.97 and iqr < 0.1.
		rho_beta_stationarity = if (isTRUE(fit$dynamic_beta) && !is.null(fit$RHO_BETA)) {
			RB <- fit$RHO_BETA
			if (is.matrix(RB) && nrow(RB) > 0L) {
				q05 <- apply(RB, 2, quantile, probs = 0.05, na.rm = TRUE)
				q25 <- apply(RB, 2, quantile, probs = 0.25, na.rm = TRUE)
				q75 <- apply(RB, 2, quantile, probs = 0.75, na.rm = TRUE)
				iqr <- q75 - q25
				warn <- (q05 >= 0.97) & (iqr < 0.1)
				data.frame(block = colnames(RB) %||% paste0("b", seq_along(q05)),
				           q05 = q05, iqr = iqr, warn = warn,
				           stringsAsFactors = FALSE, row.names = NULL)
			} else NULL
		} else NULL,
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
	if (isTRUE(x$dynamic_beta)) {
		kind_str <- if (!is.null(x$dynamic_beta_kind))
			paste0(", kind = '", x$dynamic_beta_kind, "'") else ""
		cat("Dynamic regression coefficients: enabled", kind_str, sep = "")
		if (!is.null(x$rho_beta) && length(x$rho_beta) > 0) {
			cat(" (rho_beta = ",
			    paste(names(x$rho_beta), ":", round(x$rho_beta, 2), collapse = ", "),
			    ")", sep = "")
		}
		cat("\n")
		# stationarity warning: when the rho_beta posterior 5% quantile
		# is >= 0.97 and the iqr is < 0.1, the chain is sitting near a
		# unit root. recommend refitting with dynamic_beta_kind = "rw1".
		stat <- x$rho_beta_stationarity
		if (!is.null(stat)) {
			bad <- stat[stat$warn, , drop = FALSE]
			if (nrow(bad) > 0L) {
				cat("\n  Note: rho_beta posterior is concentrated near 1 for block(s): ",
				    paste(bad$block, collapse = ", "), ".\n", sep = "")
				cat("  q05 = ", paste(sprintf("%.3f", bad$q05), collapse = ", "),
				    "; IQR = ", paste(sprintf("%.3f", bad$iqr), collapse = ", "), ".\n",
				    sep = "")
				cat("  The data are consistent with a near-random-walk;\n")
				cat("  consider refitting with dynamic_beta_kind = \"rw1\".\n")
			}
		}
	}
		# note when a longitudinal fit is static
	if (!is.null(x$n.periods) && x$n.periods > 1L &&
	    !isTRUE(x$dynamic_uv) && !isTRUE(x$dynamic_ab) && !isTRUE(x$dynamic_beta)) {
		cat("\nNote: STATIC fit pooled across", x$n.periods,
		    "time periods --\n")
		cat("  U, V, a, b are time-invariant; per-period predictions vary\n")
		cat("  only through per-period covariates. For time-varying effects,\n")
		cat("  refit with dynamic_uv = TRUE and/or dynamic_ab = TRUE.\n")
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
	cat("Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.\n")

	# per-period dynamic coefficients (when dynamic_beta is on)
	if (isTRUE(x$dynamic_beta) && !is.null(x$beta_dynamic_per_t)) {
		cat("\nDynamic coefficients per period:\n")
		cat("-------------------------------\n")
		num_tbl <- data.frame(
			Mean      = rowMeans(x$beta_dynamic_per_t),
			Min       = x$beta_dynamic_min,
			Max       = x$beta_dynamic_max,
			Drift_pct = x$beta_dynamic_drift_pct,
			row.names = rownames(x$beta_dynamic_per_t)
		)
		num_tbl <- round(num_tbl, digits)
			# attach dynamic-mask flag after rounding so numeric columns stay numeric
		if (!is.null(x$beta_dynamic_mask) && length(x$beta_dynamic_mask) == nrow(num_tbl)) {
			num_tbl$Dynamic <- ifelse(x$beta_dynamic_mask, "Y", "N")
		}
		print(num_tbl)
	}

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
