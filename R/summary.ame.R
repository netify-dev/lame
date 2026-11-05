#' Summary of an AME object
#'
#' Summarizes a fitted AME (Additive and Multiplicative Effects) model,
#' including parameter estimates, standard errors, credible intervals, and
#' model diagnostics.
#'
#' @details
#' The summary includes:
#' \describe{
#'   \item{Regression coefficients}{Posterior means, posterior standard deviations,
#'         z-values, normal-tail p-values, and 95% credible intervals for dyadic,
#'         sender, and receiver covariates. The z-value is the posterior mean
#'         divided by the posterior SD, and the p-value is the two-sided
#'         normal tail probability of that z-value: a posterior summary, not a
#'         frequentist test. The credible intervals are posterior quantiles of
#'         the draws, and the full posterior is in the BETA matrix.}
#'   \item{Variance components}{Estimates and standard errors for:
#'     \describe{
#'       \item{va}{Variance of additive sender/row effects (asymmetric networks)}
#'       \item{cab}{Covariance between sender and receiver effects}
#'       \item{vb}{Variance of additive receiver/column effects (asymmetric networks)}
#'       \item{rho}{Dyadic correlation (reciprocity in directed networks)}
#'       \item{ve}{Residual variance}
#'     }
#'     For symmetric networks, only va and ve are estimated.}
#'   \item{Rhat and ESS}{A fit of two or more chains adds the columns
#'     \code{Rhat} and \code{ESS} to both tables. \code{Rhat} is the split
#'     rank-normalised R-hat of Vehtari et al. (2021) that the between-chain
#'     gate compares with 1.05, on the draws the fit holds (the last block of
#'     every chain); for a monitored quantity it is the value in
#'     \code{fit$rhat}. \code{ESS} is the effective sample size of the pooled
#'     draws. A one-chain fit has neither column: R-hat compares chains.}
#' }
#'
#' @param object an object of class "ame", typically the result of fitting an
#'        AME model using the \code{ame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.ame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates; every row is
#'     reported (a single period carries no event-rate rule), with a note on
#'     the intercept where a saturated cap withholds it}
#'   \item{ab_var_unidentified}{\code{TRUE} when variance rows are \code{NA}}
#'   \item{ab_coef_reason}{\code{"saturated"} when a cbin or frn fit whose
#'     rows are all at the cap, or half of them with cells hidden, withholds
#'     the intercept and row-covariate coefficients (\code{NA} otherwise); a
#'     row at its cap leaves the level of the rows unbounded above, and in
#'     simulated one-period designs the intercept interval covered 0.64 of the
#'     panels with every row at the cap and 0.84 of those with cells hidden
#'     and a share from 0.5, against 0.98 or more elsewhere}
#'   \item{prior_supplied}{\code{TRUE} when the fit's \code{prior} carried
#'     \code{Sab0}, \code{eta0} or \code{etaab}; the print method names the
#'     prior}
#'   \item{odmax_saturated_share}{share of observed rows at their nomination
#'     cap, which sets the sampler's prior pooling and the intercept rule (\code{NULL} off cbin / frn)}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{ame}}, \code{\link{print.summary.ame}}
#' @method summary ame
#' @export
summary.ame <- function(object, ...) {
	.lame_require_draws(object, "summary")
	fit <- object

	# regression coefficient summaries
	beta_mean <- apply(fit$BETA, 2, mean)
	beta_sd <- apply(fit$BETA, 2, sd)
	beta_z <- ifelse(beta_sd > 0, beta_mean / beta_sd, NA_real_)
	beta_p <- 2 * (1 - pnorm(abs(beta_z)))

	# 95% quantile interval, computed as in confint() and tidy() so the three agree
	a <- (1 - 0.95) / 2
	beta_lower <- apply(fit$BETA, 2, quantile, probs = a)
	beta_upper <- apply(fit$BETA, 2, quantile, probs = 1 - a)

	# coefficient names, falling back to colnames(fit$BETA) (bipartite fits) and then beta0..betan
	if (!is.null(fit$X_names)) {
		beta_names <- fit$X_names
	} else if (!is.null(dimnames(fit[["X"]])[[3]])) {
		beta_names <- dimnames(fit[["X"]])[[3]]
	} else if (!is.null(colnames(fit$BETA))) {
		beta_names <- colnames(fit$BETA)
	} else {
		beta_names <- paste0("beta", 0:(ncol(fit$BETA) - 1))
	}

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

	if (length(beta_names) == nrow(beta_table)) {
		rownames(beta_table) <- beta_names
	}

	display_call <- tryCatch({
		format_ame_call(fit)
	}, error = function(e) {
		fit$call
	})

	# multi-chain diagnostics (R-hat, ESS) when the fit carries them
	diag_block <- NULL
	if (!is.null(fit$diagnostics) &&
			(!is.null(fit$diagnostics$rhat) || !is.null(fit$diagnostics$ess))) {
		n_chains_val <- fit$diagnostics$n_chains %||% NA_integer_
		rhat_all <- fit$diagnostics$rhat
		# monitored quantities take the gate's split rank-normalised R-hat, so the column matches fit$rhat
		gate <- fit$rhat
		if (!is.null(rhat_all) && is.data.frame(gate) && nrow(gate) > 0L) {
			hit <- intersect(names(rhat_all), gate$parameter)
			rhat_all[hit] <- gate$rhat[match(hit, gate$parameter)]
		}
		ess_all <- fit$diagnostics$ess
		# match the diagnostics to the table rows by name
		if (!is.null(rhat_all)) {
			beta_table <- cbind(
				beta_table,
				Rhat = unname(
					rhat_all[match(
						rownames(beta_table),
						names(rhat_all)
					)]
				)
			)
		}
		if (!is.null(ess_all)) {
			beta_table <- cbind(
				beta_table,
				ESS = unname(
					ess_all[match(
						rownames(beta_table),
						names(ess_all)
					)]
				)
			)
		}
		# same for the variance components
		if (!is.null(rhat_all)) {
			vc_table <- cbind(
				vc_table,
				Rhat = unname(
					rhat_all[match(
						rownames(vc_table),
						names(rhat_all)
					)]
				)
			)
		}
		if (!is.null(ess_all)) {
			vc_table <- cbind(
				vc_table,
				ESS = unname(
					ess_all[match(
						rownames(vc_table),
						names(ess_all)
					)]
				)
			)
		}
		max_rhat <- if (length(rhat_all) > 0) max(rhat_all, na.rm = TRUE) else NA_real_
		min_ess <- if (length(ess_all) > 0) min(ess_all, na.rm = TRUE) else NA_real_
		diag_block <- list(
			n_chains = n_chains_val,
			max_rhat = max_rhat,
			min_ess = min_ess
		)
	}

	# a sparse panel reports its variances as not estimated (the draws stay in
	# fit$VC) unless the effects prior was supplied
	ab_reason <- .lame_ab_var_unident_reason(fit)
	ab_data_reason <- .lame_ab_data_reason(fit)
	ab_prior_supplied <- isTRUE(fit$prior_supplied)
	ab_unident <- !is.na(ab_reason)
	if (ab_unident) {
		hide <- intersect(rownames(vc_table), .lame_ab_var_hidden_rows(ab_reason))
		if (length(hide)) vc_table[hide, ] <- NA_real_
	}
	# the intercept and row-covariate coefficients go with the sender effects;
	# an effects prior does not identify the level of a row
	coef_hide <- .lame_ab_coef_hidden_rows(.lame_ab_coef_reason(fit), rownames(beta_table))
	if (length(coef_hide)) beta_table[coef_hide, ] <- NA_real_

	# `coefficients` aliases `beta`, as in summary(lm())$coefficients
	sum_obj <- list(
		call = display_call,
		beta = beta_table,
		coefficients = beta_table,
		variance = vc_table,
		ab_var_unidentified = ab_unident,
		ab_var_unidentified_reason = ab_reason,
		ab_data_reason = ab_data_reason,
		prior_supplied = ab_prior_supplied,
		ab_prior_text = if (ab_prior_supplied) .lame_ab_prior_text(fit) else NULL,
		events_per_actor_period = fit$events_per_actor_period,
		events_per_actor = fit$events_per_actor,
		odmax_saturated_share = fit$odmax_saturated_share,
		ab_coef_reason = .lame_ab_coef_reason(fit),
		R = fit$R,
		cells_hidden = fit$cells_hidden,
		dynamic_ab = fit$dynamic_ab,
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
#' @method print summary.ame
#' @export
print.summary.ame <- function(x, digits = 3, ...) {
	.lame_cli_stdout({
		cli::cli_h1("AME model summary")
		cli::cli_h3("Call")
		# a formatted call is one string; anything else prints as R would show it
		call_lines <- if (is.character(x$call)) x$call else utils::capture.output(print(x$call))
		cli::cli_verbatim(call_lines)

		cli::cli_h3("Regression coefficients")
		coef_table <- x$beta
		coef_table <- round(coef_table, digits)
		stars <- symnum(
			x$beta[, "p_value"],
			corr = FALSE,
			na = FALSE,
			cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
			symbols = c("***", "**", "*", ".", " ")
		)
		coef_table <- cbind(coef_table, " " = stars)
		cli::cli_verbatim(utils::capture.output(print(coef_table, quote = FALSE, right = TRUE)))
		cli::cli_verbatim(
			"---",
			"Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"
		)
		cli::cli_bullets(c(
			"i" = "Stars are a visual hint from the posterior mean / SD only; for inference use the credible intervals."
		))

		cli::cli_h3("Variance components")
		var_table <- round(x$variance, digits)
		cli::cli_verbatim(utils::capture.output(print(var_table, quote = FALSE, right = TRUE)))

		if (isTRUE(x$symmetric)) {
			cli::cli_bullets(c(" " = "va = nodal variance, ve = residual variance"))
		} else {
			cli::cli_bullets(c(
				" " = "va = sender variance, cab = sender-receiver covariance, vb = receiver variance, rho = dyadic correlation, ve = residual variance"
			))
		}
		if (!is.null(x$mode) && x$mode == "bipartite") {
			cli::cli_bullets(c("i" = "Bipartite model: rho and cab are fixed to 0."))
		}
		# why va and cab are not estimated, or the prior they are reported under
		for (st_ in .lame_ab_statements(x)) cli::cli_alert_info("{st_}", wrap = TRUE)
	})
	invisible(x)
}
