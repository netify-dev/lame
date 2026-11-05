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
#'         z-values, normal-tail p-values, and 95% credible intervals for dyadic,
#'         sender, and receiver covariates. The z-value is the posterior mean
#'         divided by the posterior SD, and the p-value is the two-sided
#'         normal tail probability of that z-value: a posterior summary, not a
#'         frequentist test. The credible intervals are posterior quantiles of
#'         the draws, and the full posterior is in the BETA matrix.}
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
#'           \item{Drift}{\code{Max - Min}: the absolute range of the
#'                 per-period posterior means, in coefficient units}
#'           \item{Drift_pct}{\code{100 * Drift / |Mean|}, the drift as a
#'                 percentage of the average level. Reported as \code{NA} when
#'                 the time-average is near zero (within 5\% of the drift
#'                 range), because a trajectory that crosses zero makes the
#'                 percentage explode; read \code{Drift} in that case}
#'           \item{Dynamic}{\code{"Y"} if the coefficient was flagged as
#'                 dynamic, \code{"N"} if it was held static}
#'         }
#'         The block also prints the per-block AR(1) hyperparameters
#'         (\code{rho_beta = ...}). For per-period credible intervals use
#'         \code{\link{confint.lame}}.}
#' }
#'
#'
#' For dynamic fits the summary reports the posterior means of the
#' persistence parameters \code{rho_uv} (latent positions) and
#' \code{rho_ab} (additive effects). These are learned from period-to-period
#' transitions, so on short panels, and on sparse bipartite panels where each
#' period's positions are themselves uncertain, the data inform them weakly.
#' The summary measures how weakly. For each reported persistence it divides
#' the posterior standard deviation by the standard deviation of the prior the
#' fit used (the normal prior of \code{prior$rho_uv_mean} and
#' \code{prior$rho_uv_sd}, or of \code{prior$rho_ab_mean} and
#' \code{prior$rho_ab_sd}, truncated to the support of the sampler,
#' \eqn{|\rho| < 0.99}; see \code{prior_summary()}) and returns the ratio in
#' \code{persistence_prior_check}. A ratio near one says that the posterior is
#' as wide as the prior, so the data add little to that persistence. When the
#' ratio is above 0.6 the printed summary states it, for example
#' \code{rho_uv: posterior sd 0.073 is 0.92 of the prior sd 0.079: the data
#' barely inform it}. The cutoff, a posterior sd above 0.6 of the prior sd, is
#' the one \code{rho_beta_prior_check} applies to the coefficient persistences.
#'
#' @param object an object of class "lame", typically the result of fitting a
#'        longitudinal AME model using the \code{lame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.lame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates. On a directed
#'     binary panel of two periods or more with fewer than 1.2 events per
#'     actor-period the \code{va} and \code{vb} rows are \code{NA}: in simulated
#'     panels their 95 percent intervals covered the generating variances in
#'     0.28 to 0.91 of the panels below that rate. The draws stay in
#'     \code{fit$VC}. Under a supplied effects prior the rows are reported
#'     under it.}
#'   \item{ab_var_unidentified}{\code{TRUE} when those rows are \code{NA}}
#'   \item{ab_var_unidentified_reason}{\code{"sparse"} when the event rate put
#'     them there (\code{NA} when they are reported)}
#'   \item{ab_data_reason}{the same reason read off the data alone: it differs
#'     from \code{ab_var_unidentified_reason} when the effects prior was
#'     supplied, in which case the rows are reported under that prior}
#'   \item{ab_coef_reason}{\code{"sparse"}, \code{"saturated"} or
#'     \code{"cap_shrunk"} when the
#'     intercept and row-covariate coefficients are \code{NA} (an event rate
#'     below 0.6, or 1.6 with a multiplicative term; a one-period cbin or frn
#'     design whose rows are at the cap; a cbin or frn design of two periods or
#'     more with a cap share of 0.5 or more and a posterior-mean \code{va} from
#'     0.5 to 1.5), \code{NA} otherwise}
#'   \item{prior_supplied}{\code{TRUE} when the fit's \code{prior} carried
#'     \code{Sab0}, \code{eta0} or \code{etaab}; the print method then names
#'     the prior the variances are reported under}
#'   \item{events_per_actor_period}{events per actor per period (\code{NULL}
#'     on families without an event rate), the rate the sparse-panel rules read}
#'   \item{events_per_actor}{events per actor over the panel (\code{NULL} on
#'     families without an event rate and on fits saved before it was recorded)}
#'   \item{odmax_saturated_share}{share of observed rows at their nomination
#'     cap, which sets the sampler's prior pooling and the one-period intercept rule (\code{NULL} off cbin / frn)}
#'   \item{persistence_prior_check}{for dynamic fits, a data frame with one row
#'     per reported persistence (\code{rho_uv}, \code{rho_ab}): the posterior
#'     standard deviation \code{post_sd}, the standard deviation of the
#'     truncated prior it started from \code{prior_sd}, their ratio
#'     \code{sd_ratio}, and \code{weak}, \code{TRUE} when the ratio is above
#'     0.6 (\code{NULL} when no persistence is reported)}
#'   \item{n.periods}{Number of time periods in the longitudinal data}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{lame}}, \code{\link{print.summary.lame}}
#' @method summary lame
#' @export
summary.lame <- function(object, ...) {
	.lame_require_draws(object, "summary")
	fit <- object

	# formatted model call, or the raw call if formatting fails
	display_call <- tryCatch(
		format_ame_call(fit),
		error = function(e) fit$call
	)

	# coefficient summaries; under dynamic_beta the main table is time-collapsed
	# and the per-period means go to the print method
	beta_is_dyn <- length(dim(fit$BETA)) == 3L
	if (beta_is_dyn) {
		# time-averaged coefficient with the range of its per-period means
		beta_mean_overall <- apply(fit$BETA, 2, mean, na.rm = TRUE)
		beta_sd_overall <- apply(apply(fit$BETA, c(1, 2), mean), 2, sd, na.rm = TRUE)
		beta_min <- apply(apply(fit$BETA, c(2, 3), mean), 1, min)
		beta_max <- apply(apply(fit$BETA, c(2, 3), mean), 1, max)
		# drift in coefficient units
		beta_drift <- beta_max - beta_min
		# percentage drift relative to |time-average|; NA when |average| is below
		# 5% of the range, where the ratio is unstable
		drift_pct <- ifelse(
			beta_drift <= 0,
			0,
			ifelse(
				abs(beta_mean_overall) >= 0.05 * beta_drift,
				100 * beta_drift / abs(beta_mean_overall),
				NA_real_
			)
		)
		beta_z <- ifelse(beta_sd_overall > 0, beta_mean_overall / beta_sd_overall, NA_real_)
		beta_p <- 2 * (1 - pnorm(abs(beta_z)))
		beta_lower <- apply(
			apply(fit$BETA, c(1, 2), mean),
			2,
			quantile,
			probs = 0.025,
			na.rm = TRUE
		)
		beta_upper <- apply(
			apply(fit$BETA, c(1, 2), mean),
			2,
			quantile,
			probs = 0.975,
			na.rm = TRUE
		)
		beta_table <- cbind(
			Estimate = beta_mean_overall,
			StdError = beta_sd_overall,
			z_value = beta_z,
			p_value = beta_p,
			CI_lower = beta_lower,
			CI_upper = beta_upper
		)
		# per-period posterior means and sds
		beta_per_t_mean <- apply(fit$BETA, c(2, 3), mean)
		beta_per_t_sd <- apply(fit$BETA, c(2, 3), sd)
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
		beta_per_t_sd <- NULL
		beta_min <- NULL; beta_max <- NULL; beta_drift <- NULL; drift_pct <- NULL
	}

	# a saturated nomination cap leaves row levels unfixed: withhold the
	# intercept and row-covariate coefficients (draws stay in fit$BETA)
	coef_hide <- .lame_ab_coef_hidden_rows(
		.lame_ab_coef_reason(fit),
		rownames(beta_table)
	)
	if (length(coef_hide)) beta_table[coef_hide, ] <- NA_real_

	vc_mean <- apply(fit$VC, 2, mean)
	vc_sd <- apply(fit$VC, 2, sd)

	vc_table <- cbind(
		Estimate = vc_mean,
		StdError = vc_sd
	)
	# position hyperparameters reported only where identified (see .lame_uv_hyper_unidentified)
	uv_unident <- .lame_uv_hyper_unidentified(fit)
	pi_unident <- .lame_snap_rate_unidentified(fit)
	# effect variances of a sparse panel are withheld (draws stay in fit$VC)
	# unless the effects prior was supplied
	ab_reason <- .lame_ab_var_unident_reason(fit)
	ab_prior_supplied <- isTRUE(fit$prior_supplied)
	ab_unident <- !is.na(ab_reason)
	if (ab_unident) {
		hide <- intersect(rownames(vc_table), .lame_ab_var_hidden_rows(ab_reason))
		if (length(hide)) vc_table[hide, ] <- NA_real_
	}

	n_periods <- fit$n_time

	# `coefficients` aliases `beta`, as in summary(lm())
	sum_obj <- list(
		call = display_call,
		beta = beta_table,
		coefficients = beta_table,
		variance = vc_table,
		ab_var_unidentified = ab_unident,
		ab_var_unidentified_reason = ab_reason,
		ab_data_reason = .lame_ab_data_reason(fit),
		prior_supplied = ab_prior_supplied,
		ab_prior_text = if (ab_prior_supplied) .lame_ab_prior_text(fit) else NULL,
		events_per_actor_period = fit$events_per_actor_period,
		events_per_actor = fit$events_per_actor,
		odmax_saturated_share = fit$odmax_saturated_share,
		ab_coef_reason = .lame_ab_coef_reason(fit),
		R = fit$R,
		cells_hidden = fit$cells_hidden,
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
		beta_dynamic_drift = beta_drift,
		beta_dynamic_drift_pct = drift_pct,
		# posterior mean of the rho_beta draws; fit$rho_beta is the last draw
		rho_beta = if (isTRUE(fit$dynamic_beta)) {
			RB <- fit$RHO_BETA
			if (is.matrix(RB) && nrow(RB) > 0L) colMeans(RB, na.rm = TRUE)
			else if (!is.null(RB)) mean(RB, na.rm = TRUE)
			else fit$rho_beta
		} else NULL,
		sigma_beta = if (isTRUE(fit$dynamic_beta)) fit$sigma_beta else NULL,
		# flag blocks whose rho_beta posterior sd exceeds 0.6 of the prior sd
		rho_beta_prior_check = if (isTRUE(fit$dynamic_beta) &&
				is.matrix(fit$RHO_BETA) &&
				nrow(fit$RHO_BETA) > 1L) {
			prior_sd <- fit$prior$rho_beta_sd %||% 0.15
			post_sd <- apply(fit$RHO_BETA, 2, stats::sd, na.rm = TRUE)
			data.frame(
				block = colnames(fit$RHO_BETA) %||%
					paste0("b", seq_along(post_sd)),
				post_sd = post_sd,
				prior_sd = prior_sd,
				dominated = is.finite(post_sd) & post_sd > 0.6 * prior_sd,
				stringsAsFactors = FALSE,
				row.names = NULL
			)
		} else NULL,
		# near-unit-root check per block: q05 >= 0.97 and IQR < 0.1
		rho_beta_stationarity = if (isTRUE(fit$dynamic_beta) && !is.null(fit$RHO_BETA)) {
			RB <- fit$RHO_BETA
			if (is.matrix(RB) && nrow(RB) > 0L) {
				q05 <- apply(RB, 2, quantile, probs = 0.05, na.rm = TRUE)
				q25 <- apply(RB, 2, quantile, probs = 0.25, na.rm = TRUE)
				q75 <- apply(RB, 2, quantile, probs = 0.75, na.rm = TRUE)
				iqr <- q75 - q25
				warn <- (q05 >= 0.97) & (iqr < 0.1)
				data.frame(
					block = colnames(RB) %||% paste0("b", seq_along(q05)),
					q05 = q05,
					iqr = iqr,
					warn = warn,
					stringsAsFactors = FALSE,
					row.names = NULL
				)
			} else NULL
		} else NULL,
		rho_ab = if (isTRUE(fit$dynamic_ab) && !is.null(fit$rho_ab)) mean(fit$rho_ab) else NULL,
		sigma_ab = if (isTRUE(fit$dynamic_ab) && !is.null(fit$sigma_ab)) mean(fit$sigma_ab) else NULL,
		# posterior against prior sd of rho_uv and rho_ab
		persistence_prior_check = .lame_persistence_prior_check(fit, uv_unident),
		# withheld where not identified (see .lame_uv_hyper_unidentified)
		rho_uv = if (isTRUE(fit$dynamic_uv) && !is.null(fit$rho_uv) && !uv_unident)
			mean(fit$rho_uv) else NULL,
		sigma_uv = if (isTRUE(fit$dynamic_uv) && !is.null(fit$sigma_uv) && !uv_unident)
			mean(fit$sigma_uv) else NULL,
		uv_hyper_unidentified = uv_unident,
		# snap reset rate under the same rule (see .lame_snap_rate_unidentified)
		pi_snap = if (isTRUE(fit$dynamic_uv) && is.numeric(fit$pi_snap) &&
				length(fit$pi_snap) > 0L && !pi_unident) mean(fit$pi_snap) else NULL,
		pi_snap_unidentified = pi_unident,
		dynamic_uv_kind = fit$dynamic_uv_kind,
		refusal_share = fit$refusal_share,
		bound_line = .lame_bound_share_line(fit)
	)

	class(sum_obj) <- "summary.lame"
	return(sum_obj)
}

# sd of a normal(mean, sd^2) truncated to (lo, hi)
.lame_truncnorm_sd <- function(mean, sd, lo, hi) {
	a <- (lo - mean) / sd
	b <- (hi - mean) / sd
	z <- stats::pnorm(b) - stats::pnorm(a)
	if (!is.finite(z) || z <= 0) return(NA_real_)
	v <- 1 + (a * stats::dnorm(a) - b * stats::dnorm(b)) / z -
		((stats::dnorm(a) - stats::dnorm(b)) / z)^2
	sd * sqrt(max(v, 0))
}

# posterior against prior sd of the persistences rho_uv and rho_ab, the prior
# truncated to |rho| < 0.99; weak means a ratio above 0.6
.lame_persistence_prior_check <- function(fit, uv_unident = FALSE) {
	# exact lookup, since $prior would partial-match prior_supplied
	prior <- fit[["prior", exact = TRUE]]
	if (!is.list(prior)) return(NULL)
	spec <- list(
		rho_uv = list(
			on = isTRUE(fit$dynamic_uv) && !isTRUE(uv_unident),
			draws = fit$rho_uv,
			mean = prior$rho_uv_mean,
			sd = prior$rho_uv_sd
		),
		rho_ab = list(
			on = isTRUE(fit$dynamic_ab),
			draws = fit$rho_ab,
			mean = prior$rho_ab_mean,
			sd = prior$rho_ab_sd
		)
	)
	rows <- lapply(names(spec), function(nm) {
		x <- spec[[nm]]
		if (!x$on || !is.numeric(x$draws) || length(x$draws) < 2L ||
				length(x$mean) != 1L || length(x$sd) != 1L ||
				!is.finite(x$mean) || !is.finite(x$sd) || x$sd <= 0) return(NULL)
		post_sd <- stats::sd(x$draws, na.rm = TRUE)
		prior_sd <- .lame_truncnorm_sd(x$mean, x$sd, -0.99, 0.99)
		# a chain that never moved has no posterior sd to compare
		if (!is.finite(post_sd) || post_sd <= 0 || !is.finite(prior_sd) || prior_sd <= 0)
			return(NULL)
		data.frame(
			parameter = nm,
			post_sd = post_sd,
			prior_sd = prior_sd,
			sd_ratio = post_sd / prior_sd,
			weak = post_sd > 0.6 * prior_sd,
			stringsAsFactors = FALSE
		)
	})
	rows <- rows[!vapply(rows, is.null, logical(1))]
	if (!length(rows)) return(NULL)
	out <- do.call(rbind, rows)
	rownames(out) <- NULL
	out
}

# summary line for a persistence whose posterior is nearly as wide as its prior
.lame_persistence_prior_line <- function(check, parameter) {
	if (is.null(check)) return(NULL)
	r <- check[check$parameter == parameter, , drop = FALSE]
	if (nrow(r) != 1L || !isTRUE(r$weak)) return(NULL)
	sprintf(
		"%s: posterior sd %.3f is %.2f of the prior sd %.3f: the data %s.",
		parameter,
		r$post_sd,
		r$sd_ratio,
		r$prior_sd,
		if (r$sd_ratio >= 0.8) "barely inform it" else "inform it only weakly"
	)
}

#' Print method for summary.lame objects
#'
#' Prints a formatted summary of a LAME model fit
#'
#' @param x a summary.lame object
#' @param digits number of digits to display (default: 3)
#' @param ... additional arguments (not used)
#' @return the summary.lame object invisibly
#' @method print summary.lame
#' @export
print.summary.lame <- function(x, digits = 3, ...) {
	.lame_cli_stdout({
		cli::cli_h1("Longitudinal AME model summary")
		cli::cli_h3("Call")
		# a formatted call is one string; anything else prints as R would show it
		call_lines <- if (is.character(x$call)) x$call else utils::capture.output(print(x$call))
		cli::cli_verbatim(call_lines)

		cli::cli_h3("Model")
		model_info <- character()
		if (!is.null(x$n.periods)) {
			model_info <- c(model_info, "*" = "Time periods: {.val {x$n.periods}}")
		}
		if (!is.null(x$family)) {
			model_info <- c(model_info, "*" = "Family: {.field {x$family}}")
		}
		if (!is.null(x$mode)) {
			model_info <- c(model_info, "*" = "Mode: {.field {x$mode}}")
		}
		cli::cli_bullets(model_info)
		if (!is.null(x$bound_line)) cli::cli_alert_info("{x$bound_line}", wrap = TRUE)
		if (isTRUE(x$dynamic_uv)) {
			uv_parts <- c(
				if (!is.null(x$dynamic_uv_kind)) paste0("kind = ", x$dynamic_uv_kind),
				if (!is.null(x$rho_uv)) paste0("rho_uv = ", round(x$rho_uv, 3)),
				if (!is.null(x$pi_snap)) paste0("pi_snap = ", round(x$pi_snap, 3))
			)
			uv_txt <- if (length(uv_parts)) paste0(" (", paste(uv_parts, collapse = ", "), ")") else ""
			cli::cli_bullets(c("*" = "Dynamic latent positions: enabled{uv_txt}"))
			if (isTRUE(x$uv_hyper_unidentified)) {
				cli::cli_bullets(c(
					" " = paste0(
						"rho_uv and sigma_uv are not reported for this transition kind at rank two or more ",
						"(the draws are in {.code fit$rho_uv} and {.code fit$sigma_uv})."
					)
				))
			}
			if (isTRUE(x$pi_snap_unidentified)) {
				cli::cli_bullets(c(
					" " = paste0(
						"The reset rate pi_snap is not reported at rank two or more ",
						"(the draws are in {.code fit$pi_snap})."
					)
				))
			}
			pl_uv <- .lame_persistence_prior_line(x$persistence_prior_check, "rho_uv")
			if (!is.null(pl_uv)) cli::cli_bullets(c(" " = "{pl_uv}"))
		}
		if (isTRUE(x$dynamic_ab)) {
			ab_txt <- if (!is.null(x$rho_ab)) paste0(" (rho_ab = ", round(x$rho_ab, 3), ")") else ""
			cli::cli_bullets(c("*" = "Dynamic additive effects: enabled{ab_txt}"))
			pl_ab <- .lame_persistence_prior_line(x$persistence_prior_check, "rho_ab")
			if (!is.null(pl_ab)) cli::cli_bullets(c(" " = "{pl_ab}"))
		}
		if (isTRUE(x$dynamic_beta)) {
			beta_parts <- c(
				if (!is.null(x$dynamic_beta_kind)) paste0("kind = ", x$dynamic_beta_kind),
				if (!is.null(x$rho_beta) && length(x$rho_beta) > 0) {
					paste0(
						"rho_beta = ",
						paste(names(x$rho_beta), ":", round(x$rho_beta, 2), collapse = ", ")
					)
				}
			)
			beta_txt <- if (length(beta_parts)) paste0(" (", paste(beta_parts, collapse = "; "), ")") else ""
			cli::cli_bullets(c("*" = "Dynamic regression coefficients: enabled{beta_txt}"))
			# rho_beta weakly informed by the data
			pc <- x$rho_beta_prior_check
			if (!is.null(pc) && any(pc$dominated)) {
				dom <- pc[pc$dominated, , drop = FALSE]
				dom_blocks <- paste(dom$block, collapse = ", ")
				dom_post <- paste(sprintf("%.2f", dom$post_sd), collapse = ", ")
				dom_prior <- sprintf("%.2f", dom$prior_sd[1])
				cli::cli_bullets(c(
					"i" = "rho_beta is weakly informed by the data for block(s): {dom_blocks} (posterior sd {dom_post} vs prior sd {dom_prior}).",
					" " = "The posterior of rho_beta stays near its prior; {.fn dynamic_beta_prior_summary} shows that prior, and {.code prior = list(rho_beta_mean = , rho_beta_sd = )} sets it."
				))
			}
			# rho_beta near a unit root: suggest rw1
			stat <- x$rho_beta_stationarity
			if (!is.null(stat)) {
				bad <- stat[stat$warn, , drop = FALSE]
				if (nrow(bad) > 0L) {
					bad_blocks <- paste(bad$block, collapse = ", ")
					bad_q05 <- paste(sprintf("%.3f", bad$q05), collapse = ", ")
					bad_iqr <- paste(sprintf("%.3f", bad$iqr), collapse = ", ")
					cli::cli_bullets(c(
						"!" = "The rho_beta posterior is concentrated near 1 for block(s): {bad_blocks}.",
						" " = "q05 = {bad_q05}; IQR = {bad_iqr}.",
						" " = "The data are consistent with a near-random-walk; consider refitting with {.code dynamic_beta_kind = \"rw1\"}."
					))
				}
			}
		}
		# note a static fit to a longitudinal panel
		if (!is.null(x$n.periods) && x$n.periods > 1L &&
				!isTRUE(x$dynamic_uv) && !isTRUE(x$dynamic_ab) && !isTRUE(x$dynamic_beta)) {
			cli::cli_alert_warning(
				paste0(
					"This is a {.strong STATIC} fit pooled across {x$n.periods} time periods: ",
					"U, V, a, b are time-invariant, and per-period predictions vary only through per-period covariates."
				),
				wrap = TRUE
			)
			cli::cli_alert_info(
				"For time-varying effects, refit with {.code dynamic_uv = TRUE} and/or {.code dynamic_ab = TRUE}.",
				wrap = TRUE
			)
		}

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

		# per-period dynamic coefficients
		if (isTRUE(x$dynamic_beta) && !is.null(x$beta_dynamic_per_t)) {
			cli::cli_h3("Dynamic coefficients per period")
			num_tbl <- data.frame(
				Mean = rowMeans(x$beta_dynamic_per_t),
				Min = x$beta_dynamic_min,
				Max = x$beta_dynamic_max,
				Drift = x$beta_dynamic_drift %||%
					(x$beta_dynamic_max - x$beta_dynamic_min),
				Drift_pct = x$beta_dynamic_drift_pct,
				row.names = rownames(x$beta_dynamic_per_t)
			)
			num_tbl <- round(num_tbl, digits)
			# add the dynamic flag after rounding so the numeric columns stay numeric
			if (!is.null(x$beta_dynamic_mask) && length(x$beta_dynamic_mask) == nrow(num_tbl)) {
				num_tbl$Dynamic <- ifelse(x$beta_dynamic_mask, "Y", "N")
			}
			cli::cli_verbatim(utils::capture.output(print(num_tbl)))
		}

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
		if (isTRUE(x$dynamic_ab)) {
			cli::cli_bullets(c(
				" " = "Under {.code dynamic_ab}, va and vb are the variances of the actor levels; rho_ab and sigma_ab describe the per-period deviations around them."
			))
		}
		if (!is.null(x$mode) && x$mode == "bipartite") {
			cli::cli_bullets(c("i" = "Bipartite model: rho and cab are fixed to 0."))
		}
		# why va, vb and cab are withheld, or the prior they are reported under
		for (st_ in .lame_ab_statements(x)) cli::cli_alert_info("{st_}", wrap = TRUE)
	})
	invisible(x)
}
