# global variables for R CMD check
if (getRversion() >= "2.15.1") {
	utils::globalVariables("pit")
}

# probability-integral-transform (PIT) calibration check of h-step forecasts from a
# dynamic lame() fit; discrete families use the randomised PIT

# deterministic Poisson-lognormal CDF and PMF with bounded omitted tails; for y >= 1,
# ppois(y - 1, lambda) = P(Gamma(y, 1) > lambda), whose quantiles isolate sharp
# transitions before integration
.lame_pois_cdf_mass <- function(y, eta, ve, tolerance = 1e-9) {
	stopifnot(length(y) == length(eta), length(ve) == length(eta))
	if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 || tolerance >= .01 ||
			any(!is.finite(y) | y < 0 | y != floor(y)) || any(!is.finite(eta)) ||
			any(!is.finite(ve) | ve < 0))
		cli::cli_abort("Poisson predictive scoring requires finite counts, predictors, nonnegative variances and a valid accuracy tolerance.")
	answer <- matrix(NA_real_, length(y), 2L, dimnames = list(NULL, c("lower", "mass")))
	for (i in seq_along(y)) {
		if (ve[i] == 0) {
			lambda <- exp(eta[i])
			answer[i, ] <- c(stats::ppois(y[i] - 1, lambda), stats::dpois(y[i], lambda))
			next
		}
		sigma <- sqrt(ve[i]); epsilon <- tolerance / 128
		bound <- stats::qnorm(1 - epsilon)
		integrate_one <- function(fn, lo, hi) {
			if (lo >= hi) return(0)
			res <- stats::integrate(
				fn,
				lo,
				hi,
				subdivisions = 1000L,
				abs.tol = tolerance / 64,
				rel.tol = tolerance / 64,
				stop.on.error = TRUE
			)
			if (!is.finite(res$value) || !is.finite(res$abs.error) || res$abs.error > tolerance / 32)
				cli::cli_abort("Poisson predictive integration did not meet its accuracy tolerance.")
			res$value
		}
		transition <- (log(max(y[i], 1)) - eta[i]) / sigma
		gamma_bounds <- (log(stats::qgamma(c(epsilon, 1 - epsilon), shape = max(y[i], 1))) - eta[i]) / sigma
		if (y[i] == 0) {
			lower <- 0
		} else {
			# the truncated tails cost at most 3 epsilon of CDF error; splitting at zero and at
			# the transition avoids quadrature-node aliasing
			lo <- max(gamma_bounds[1L], -bound); hi <- min(gamma_bounds[2L], bound)
			cdf_fn <- function(z) stats::ppois(y[i] - 1, exp(eta[i] + sigma * z)) * stats::dnorm(z)
			lower <- stats::pnorm(gamma_bounds[1L])
			if (lo < hi) {
				cuts <- sort(unique(c(lo, hi, c(0, transition)[c(0, transition) > lo & c(0, transition) < hi])))
				lower <- lower + sum(vapply(
					seq_len(length(cuts) - 1L),
					function(j) integrate_one(cdf_fn, cuts[j], cuts[j + 1L]),
					numeric(1)
				))
			}
		}
		# omitted normal tails cost at most 2 epsilon of PMF; the transition cuts keep a narrow peak
		bound <- stats::qnorm(1 - epsilon)
		width <- 12 / (sigma * sqrt(max(y[i], 1)))
		cuts <- sort(unique(pmax(
			-bound,
			pmin(
				bound,
				c(
					-bound,
					bound,
					gamma_bounds,
					transition,
					transition - width,
					transition + width
				)
			)
		)))
		mass_fn <- function(z) stats::dpois(y[i], exp(eta[i] + sigma * z)) * stats::dnorm(z)
		mass <- sum(vapply(
			seq_len(length(cuts) - 1L),
			function(j) integrate_one(mass_fn, cuts[j], cuts[j + 1L]),
			numeric(1)
		))
		answer[i, ] <- c(lower, mass)
	}
	if (any(!is.finite(answer)) || any(answer < 0) || any(rowSums(answer) > 1 + tolerance))
		cli::cli_abort("Invalid integrated Poisson predictive CDF.")
	answer
}

.forecast_pit_robust <- function(fit, y_future, h, n_draws, newdata, newdesign, newexposure) {
	if (!is.null(newdesign) || !is.null(newexposure))
		cli::cli_abort("Robust PIT uses newdata dyadic covariates and has no exposure or prepared-design argument.")
	z <- fit$draws; D <- fit$data
	required <- c("mu", "A", "B", "a0", "b0", "U", "V", "beta", "lambda")
	hyper_names <- c("ve", "rho_mu", "rho_ab", "rho_uv", "var_mu", "var_ab", "var_uv")
	if (!is.list(D) || !is.list(z) || !is.matrix(z$hyper) || nrow(z$hyper) < 1L ||
			!all(hyper_names %in% colnames(z$hyper)) || any(
				vapply(required, function(nm) is.null(z[[nm]]), logical(1))
			))
		cli::cli_abort(
			c(
				"Robust PIT requires saved joint future-state draws and persistent pair precisions.",
				"i" = "Refit with the current lame_robust() implementation retaining its complete draws."
			),
			class = "lame_forecast_robust_draws"
		)
	fc <- predict.lame_robust(
		fit,
		h = h,
		newdata = newdata,
		type = "link",
		by_draw = TRUE,
		nsim = n_draws
	)
	idx <- attr(fc, "draw_index"); dm <- dim(fc)
	z <- fit$draws; df <- fit$observation_df %||% Inf
	if (length(df) != 1L || is.na(df) || df <= 2 ||
			length(dm) != 4L || dm[3L] != h || length(idx) != dm[4L] ||
			!is.matrix(z$hyper) || !"ve" %in% colnames(z$hyper) ||
			any(!is.finite(idx) | idx != floor(idx) | idx < 1 | idx > nrow(z$hyper)) ||
			length(dim(z$lambda)) != 3L || !identical(dim(z$lambda)[1:2], dm[1:2]) ||
			max(idx) > dim(z$lambda)[3L])
		cli::cli_abort("Robust PIT requires matched future states, observation variances and persistent pair precisions.")
	ve <- z$hyper[idx, "ve"]
	if (any(!is.finite(ve) | ve <= 0)) cli::cli_abort("Robust observation variances must be finite and positive.")
	pit <- numeric(0); ids <- list()
	for (k in seq_len(h)) {
		Y <- .lfo_align_heldout(y_future[[k]], dimnames(fc)[[1L]], dimnames(fc)[[2L]], dm[1L], dm[2L], k)
		ok <- is.finite(Y); diag(ok) <- FALSE
		ix <- which(ok); if (!length(ix)) next
		eta <- matrix(fc[, , k, , drop = FALSE], prod(dm[1:2]), dm[4L])[ix, , drop = FALSE]
		lambda <- matrix(z$lambda[, , idx, drop = FALSE], prod(dm[1:2]), dm[4L])[ix, , drop = FALSE]
		if (any(!is.finite(eta)) || any(!is.finite(lambda) | lambda <= 0))
			cli::cli_abort("Robust held-out cells require finite predictors and positive pair precisions.")
		scale <- sqrt(matrix(rep(ve, each = length(ix)), length(ix)) / lambda)
		if (any(
			!is.finite(scale) | scale <= 0
		)) cli::cli_abort("Robust predictive scales cannot be resolved numerically.")
		standardized <- (Y[ix] - eta) / scale
		cdf <- if (is.infinite(df)) stats::pnorm(standardized) else stats::pt(standardized, df)
		pit <- c(pit, rowMeans(cdf))
		ij <- arrayInd(ix, dim(Y))
		ids[[length(ids) + 1L]] <- data.frame(
			period = k,
			row_actor = (rownames(Y) %||% seq_len(nrow(Y)))[ij[, 1L]],
			col_actor = (colnames(Y) %||% seq_len(ncol(Y)))[ij[, 2L]]
		)
	}
	.lame_pit_result(pit, pit, pit, idx, ids, fit$family, h)
}
.lame_pit_result <- function(pit, lower, upper, idx, ids, family, h) {
	ks <- if (length(pit)) suppressWarnings(stats::ks.test(pit, "punif")) else NULL
	out <- list(
		pit = pit,
		ks_stat = if (is.null(ks)) NA_real_ else unname(ks$statistic),
		ks_p = NA_real_,
		cover_95 = if (length(pit)) mean(pit >= .025 & pit <= .975) else NA_real_,
		family = family,
		h = h,
		n_obs = length(pit),
		ks_inference = "Global p-value unavailable: network dependence is not calibrated",
		cdf_lower = lower,
		cdf_upper = upper,
		draw_index = idx,
		cells = if (length(ids)) do.call(rbind, ids) else data.frame()
	)
	class(out) <- c("forecast_pit", "list"); out
}

#' Probability-integral-transform calibration check for h-step forecasts
#'
#' For a fit with at least one dynamic component, evaluates how well the
#' \code{h}-step posterior-predictive distribution covers the actually
#' observed outcomes at the corresponding period(s). A well-calibrated
#' forecast produces approximately Uniform(0, 1) PIT values.
#'
#' Workflow: fit on periods \code{1:(T - h)}, forecast \code{h} steps
#' ahead, score against the held-out \code{T - h + 1 : T} observations.
#' Returns a numeric vector of PIT values (one per observed dyad in the
#' held-out periods), plus a one-row summary with the Kolmogorov-Smirnov
#' statistic against Uniform(0, 1) and a fraction-of-observed-dyads
#' coverage diagnostic.
#'
#' \strong{Families.} For \code{family = "normal"}, the
#' PIT averages the Gaussian CDF over matched predictor and residual-variance
#' draws. For
#' discrete families (\code{"binary"}, \code{"poisson"}, \code{"ordinal"})
#' the randomised PIT of Czado, Gneiting & Held (2009) is used. Poisson forecasts
#' integrate the Gaussian observation effect and future exposure. Ordinal scoring
#' requires saved explicit posterior cutpoints; data-induced joint rank simulation
#' does not define a fixed future threshold law. For rank
#' families (\code{"cbin"}, \code{"frn"}) PIT is not
#' implemented; the function returns NULL with an informational note.
#' Robust fits from \code{lame_robust()} use their genuine future latent states,
#' matched persistent pair precisions and Gaussian or transient Student residual CDF.
#' Fitted temporary observation precisions are never reused for future outcomes.
#'
#' \strong{Plotting.} The companion \code{plot()} method on the returned
#' object renders a histogram with a uniform reference line. Use
#' \code{ggplot2::ggplot(pit$pit) + geom_histogram()} for custom plots.
#'
#' @param fit A fitted \code{lame} or \code{lame_robust} object trained on the first
#'   \code{T - h} periods (the function does not refit; the user is
#'   responsible for the train/forecast split).
#' @param y_future A list of length \code{h} of held-out observed
#'   matrices, one per forecast period. The same shape as
#'   the fitted actor grid. Named periods are aligned to shared actors; unshared
#'   actors and missing outcomes are excluded. Unipartite self cells are excluded
#'   and symmetric networks score each unordered dyad once. Bipartite grids retain
#'   cells whose row and column positions coincide.
#' @param h The forecast horizon (length of \code{y_future}). Inferred
#'   from \code{length(y_future)} when not supplied.
#' @param n_draws Number of posterior draws to use for the forecast.
#'   Default uses all stored draws.
#' @param newdata Future dyadic covariate arrays passed to the forecaster. Required
#'   for robust fits that contain dyadic covariates, as in \code{predict.lame_robust()}.
#' @param newdesign Full prepared future designs, including nodal roles. Supported
#'   only by ordinary lame fits.
#' @param newexposure Future Poisson exposures, one per horizon. Defaults to the
#'   last fitted exposure. Zero exposure defines a degenerate zero-count law.
#' @return A list of class \code{forecast_pit} with components
#'   \code{pit} (numeric vector), \code{ks_stat} (KS statistic against
#'   Uniform), \code{ks_p} (\code{NA}; a dependent-network null reference is
#'   unavailable), \code{cover_95} (fraction of
#'   PIT values in \code{[0.025, 0.975]}), matched \code{draw_index}, predictive
#'   \code{cdf_lower}/\code{cdf_upper}, and scored actor/period \code{cells}.
#'   Network cells need not be independent. The KS statistic, histogram and coverage
#'   are descriptive. No global model-adequacy p-value is reported without a
#'   calibrated joint network reference. \code{ks_inference} records this restriction.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' Y_train <- lapply(YX_bin_list$Y, function(y) 1 * (y > 0))[1:3]
#' Y_test <- lapply(YX_bin_list$Y, function(y) 1 * (y > 0))[4]
#' X_train <- YX_bin_list$X[1:3]
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(Y_train, Xdyad = X_train,
#'             family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE, chains = 1)
#' pit <- forecast_pit(fit, y_future = Y_test)
#' pit$ks_stat
#' plot(pit)
#' }
#' @references Czado, C., Gneiting, T., & Held, L. (2009). Predictive
#'   model assessment for count data. \emph{Biometrics}, 65(4), 1254-1261.
#' @export
forecast_pit <- function(
	fit,
	y_future,
	h = NULL,
	n_draws = NULL,
	newdata = NULL,
	newdesign = NULL,
	newexposure = NULL
) {
	if (!inherits(fit, "lame") && !inherits(fit, "lame_robust"))
		cli::cli_abort("PIT requires a fitted lame or lame_robust object.")
	if (!is.list(y_future)) y_future <- list(y_future)
	if (is.null(h)) h <- length(y_future)
	if (length(h) != 1L || !is.finite(h) || h < 1 || h != floor(h) ||
			length(y_future) != h)
		cli::cli_abort("h must be a positive integer matching the held-out period list.")
	if (!all(vapply(
		y_future,
		function(y) is.matrix(y) && is.numeric(y) &&
			all(is.na(y) | is.finite(y)),
		logical(1)
	)))
		cli::cli_abort("Held-out outcomes must be numeric matrices with finite values or NA.")
	for (y in y_future) for (nm in list(rownames(y), colnames(y)))
		if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)))
		cli::cli_abort("Held-out actor names must be nonmissing, nonempty and unique.")
	if (inherits(fit, "lame_robust"))
		return(.forecast_pit_robust(fit, y_future, h, n_draws, newdata, newdesign, newexposure))
	family <- fit$family
	if (length(family) != 1L || !family %in% c("normal", "binary", "poisson", "ordinal", "cbin", "frn"))
		cli::cli_abort("Unsupported PIT observation family.")
	if (family %in% c("cbin", "frn")) {
		cli::cli_inform("PIT is not implemented for rank-based observation families.")
		return(invisible(NULL))
	}
	if (family == "ordinal" && !identical(fit$ordinal_cutpoints, "explicit"))
		cli::cli_abort(
			c(
				"Ordinal PIT requires saved explicit posterior cutpoints.",
				"i" = "Refit with ordinal_cutpoints = 'explicit'. Data-induced ordinal simulation uses joint random network quantiles, not a saved future threshold law."
			),
			class = "lame_forecast_ordinal_model"
		)
	if (family == "poisson") {
		newexposure <- newexposure %||% rep(utils::tail(fit$period_exposure %||% 1, 1L), h)
		if (length(newexposure) != h || any(!is.finite(newexposure) | newexposure < 0))
			cli::cli_abort("Future exposures must be finite nonnegative values, one per horizon.")
	} else if (!is.null(newexposure)) cli::cli_abort("newexposure is used only by Poisson forecasts.")
	# forecast the latent mean states once; their posterior indices select observation
	# variances and cutpoints from the same joint draws
	fc <- .forecast_lame(
		fit,
		h = h,
		type = "link",
		by_draw = TRUE,
		n_draws = n_draws,
		newdata = newdata,
		newdesign = newdesign
	)
	dm <- dim(fc); idx <- attr(fc, "draw_index")
	if (length(dm) != 4L || dm[3L] != h || dm[4L] < 1L ||
			length(idx) != dm[4L] || any(
				!is.finite(idx) | idx != floor(idx) |
					idx < 1 | idx > nrow(fit$BETA)
			))
		cli::cli_abort("Forecast scoring requires matched joint posterior draw indices.")
	S <- dm[4L]
	ve <- rep(1, S)
	if (family %in% c("normal", "poisson")) {
		if (!is.matrix(fit$VC) || !"ve" %in% colnames(fit$VC) || max(idx) > nrow(fit$VC))
			cli::cli_abort("Forecast scoring requires matched named observation variance draws.")
		ve <- fit$VC[idx, "ve"]
		if (any(!is.finite(ve) | if (family == "normal") ve <= 0 else ve < 0))
			cli::cli_abort("Forecast observation variance draws are invalid.")
	}
	cuts <- NULL; levels <- NULL
	if (family == "ordinal") {
		levels <- fit$ordinal_levels
		if (!is.numeric(levels) || length(levels) < 2L || any(!is.finite(levels)) ||
				any(diff(levels) <= 0)) cli::cli_abort("Ordinal category levels are missing or invalid.")
		K <- length(levels)
		if (K > 2L && (!is.matrix(fit$ALPHA) || ncol(fit$ALPHA) != K - 2L ||
				max(idx) > nrow(fit$ALPHA))) cli::cli_abort("Ordinal PIT requires matched posterior cutpoint draws.")
		cuts <- cbind(-Inf, 0, if (K > 2L) fit$ALPHA[idx, , drop = FALSE] else NULL, Inf)
		if (nrow(cuts) != S) cuts <- cuts[rep(1L, S), , drop = FALSE]
		if (any(!is.finite(cuts[, 2:K, drop = FALSE])) ||
				any(apply(cuts, 1L, function(x) any(diff(x) <= 0))))
			cli::cli_abort("Ordinal posterior cutpoints must be finite and strictly increasing.")
	}
	pit <- lower_all <- upper_all <- numeric(0); ids <- list()
	for (k in seq_len(h)) {
		Y <- .lfo_align_heldout(y_future[[k]], dimnames(fc)[[1L]], dimnames(fc)[[2L]], dm[1L], dm[2L], k)
		ok <- is.finite(Y)
		if (!identical(fit$mode, "bipartite")) {
			if (nrow(Y) != ncol(Y)) cli::cli_abort("Unipartite forecast grids must be square.")
			diag(ok) <- FALSE
			if (isTRUE(fit$symmetric)) ok[lower.tri(ok)] <- FALSE
		}
		ix <- which(ok); if (!length(ix)) next
		y <- Y[ix]
		eta <- matrix(fc[, , k, , drop = FALSE], prod(dm[1:2]), S)[ix, , drop = FALSE]
		if (any(!is.finite(eta))) cli::cli_abort("Observed held-out cells require finite forecast predictors in every selected draw.")
		if (family == "normal") {
			lower <- upper <- rowMeans(matrix(
				stats::pnorm(
					y - eta,
					sd = matrix(rep(sqrt(ve), each = length(ix)), length(ix))
				),
				length(ix),
				S
			))
		} else if (family == "binary") {
			if (any(!y %in% c(0, 1))) cli::cli_abort("Binary held-out outcomes must be zero or one.")
			p <- rowMeans(stats::pnorm(eta))
			lower <- ifelse(y == 0, 0, 1 - p); upper <- ifelse(y == 0, 1 - p, 1)
		} else if (family == "poisson") {
			if (any(y < 0 | y != floor(y))) cli::cli_abort("Poisson held-out outcomes must be nonnegative integer counts.")
			if (newexposure[k] == 0) {
				lower <- ifelse(y == 0, 0, 1); upper <- rep(1, length(y))
			} else {
				cdf <- .lame_pois_cdf_mass(
					rep(y, S),
					c(eta) + log(newexposure[k]),
					rep(ve, each = length(ix))
				)
				lower <- rowMeans(matrix(cdf[, "lower"], length(ix), S))
				upper <- lower + rowMeans(matrix(cdf[, "mass"], length(ix), S))
			}
		} else {
			category <- match(y, levels)
			if (anyNA(category)) cli::cli_abort("Held-out ordinal outcomes contain an unknown category label.")
			L <- U <- matrix(NA_real_, length(ix), S)
			for (s in seq_len(S)) {
				L[, s] <- stats::pnorm(cuts[s, category] - eta[, s])
				U[, s] <- stats::pnorm(cuts[s, category + 1L] - eta[, s])
			}
			lower <- rowMeans(L); upper <- rowMeans(U)
		}
		if (any(!is.finite(lower) | !is.finite(upper) | lower < 0 | upper > 1 + 1e-9 | lower > upper))
			cli::cli_abort("Invalid posterior predictive CDF bounds.")
		lower <- pmin(lower, 1); upper <- pmin(upper, 1) # bounded quadrature roundoff only
		value <- if (family == "normal") lower else lower + stats::runif(length(ix)) * (upper - lower)
		pit <- c(pit, value); lower_all <- c(lower_all, lower); upper_all <- c(upper_all, upper)
		ij <- arrayInd(ix, dim(Y))
		ids[[length(ids) + 1L]] <- data.frame(
			period = k,
			row_actor = (rownames(Y) %||% seq_len(nrow(Y)))[ij[, 1L]],
			col_actor = (colnames(Y) %||% seq_len(ncol(Y)))[ij[, 2L]]
		)
	}
	.lame_pit_result(pit, lower_all, upper_all, idx, ids, family, h)
}
#' @method print forecast_pit
#' @export
print.forecast_pit <- function(x, digits = 3, ...) {
	.lame_cli_stdout({
		cli::cli_h2("Forecast PIT calibration check")
		cli::cli_ul(c(
			"Family: {.val {x$family}}",
			"Forecast horizon: {.val {x$h}}",
			"Held-out dyads: {.val {x$n_obs}}",
			"Kolmogorov-Smirnov statistic against Uniform(0, 1): {.val {signif(x$ks_stat, digits)}}",
			"Global calibration p-value: unavailable for dependent network cells",
			"Fraction in [0.025, 0.975]: {.val {signif(x$cover_95, digits)}} (expect {.val {0.95}})"
		))
		cli::cli_text("These are descriptive diagnostics. Network dependence requires a calibrated joint reference for global inference.")
	})
	invisible(x)
}

#' @method plot forecast_pit
#' @export
plot.forecast_pit <- function(x, ...) {
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg ggplot2} is required for the PIT plot.")
	}
	df <- data.frame(pit = x$pit)
	ggplot2::ggplot(df, ggplot2::aes(x = pit)) +
		ggplot2::geom_histogram(bins = 20) +
		ggplot2::geom_hline(
			yintercept = nrow(df) / 20,
			linetype = "dashed"
		) +
		ggplot2::labs(
			x = "PIT value",
			y = "Count",
			title = "Forecast calibration check",
			subtitle = sprintf(
				"Horizon %d, %s family; descriptive KS statistic = %.3f",
				x$h,
				x$family,
				x$ks_stat
			)
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			plot.caption = ggplot2::element_text(hjust = 0),
			axis.ticks = ggplot2::element_blank()
		)
}
