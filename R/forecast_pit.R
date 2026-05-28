# global variables for R CMD check
if(getRversion() >= "2.15.1") {
	utils::globalVariables("pit")
}

# probability-integral-transform (PIT) calibration check for h-step
# forecasts from a dynamic lame() fit. for any continuous family the PIT
# values P(Y_obs <= y | forecast) should be Uniform(0, 1) when the forecast
# is calibrated. for discrete families we use the randomised PIT.

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
#' \strong{Families.} For \code{family = "normal"} / \code{"tobit"}, the
#' PIT is computed analytically from the per-draw forecast variance. For
#' discrete families (\code{"binary"}, \code{"poisson"}, \code{"ordinal"})
#' the randomised PIT of Czado, Gneiting & Held (2009) is used. For rank
#' families (\code{"cbin"}, \code{"frn"}, \code{"rrl"}) PIT is not
#' implemented; the function returns NULL with an informational note.
#'
#' \strong{Plotting.} The companion \code{plot()} method on the returned
#' object renders a histogram with a uniform reference line. Use
#' \code{ggplot2::ggplot(pit$pit) + geom_histogram()} for custom plots.
#'
#' @param fit A fitted \code{lame} object trained on the first
#'   \code{T - h} periods (the function does not refit; the user is
#'   responsible for the train/forecast split).
#' @param y_future A list of length \code{h} of held-out observed
#'   matrices, one per forecast period. The same shape as
#'   \code{fit$YPM[[t]]} elements.
#' @param h The forecast horizon (length of \code{y_future}). Inferred
#'   from \code{length(y_future)} when not supplied.
#' @param n_draws Number of posterior draws to use for the forecast.
#'   Default uses all stored draws.
#' @return A list of class \code{forecast_pit} with components
#'   \code{pit} (numeric vector), \code{ks_stat} (KS statistic against
#'   Uniform), \code{ks_p} (KS p-value), \code{cover_95} (fraction of
#'   PIT values in \code{[0.025, 0.975]}).
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' Y_train <- YX_bin_list$Y[1:3]
#' Y_test <- YX_bin_list$Y[4]
#' X_train <- YX_bin_list$X[1:3]
#' fit <- lame(Y_train, Xdyad = X_train,
#'             family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 200, burn = 50, odens = 5, verbose = FALSE)
#' pit <- forecast_pit(fit, y_future = Y_test)
#' pit$ks_p
#' plot(pit)
#' }
#' @references Czado, C., Gneiting, T., & Held, L. (2009). Predictive
#'   model assessment for count data. \emph{Biometrics}, 65(4), 1254-1261.
#' @export
forecast_pit <- function(fit, y_future, h = NULL, n_draws = NULL) {
	if (!is.list(y_future)) y_future <- list(y_future)
	if (is.null(h)) h <- length(y_future)
	family <- fit$family %||% "normal"

	if (family %in% c("cbin", "frn", "rrl")) {
		cli::cli_inform(c(
			"i" = "PIT is not implemented for {.val {family}} (rank-based likelihood)."))
		return(invisible(NULL))
	}

	# build h-step draw-level forecasts on the response scale
	fc_full <- predict(fit, h = h, type = "response", by_draw = TRUE,
	                   n_draws = n_draws)
	# fc_full is a 4-D [n, n, h, S] array
	S <- dim(fc_full)[4L]

	pit_vec <- numeric(0)
	for (k in seq_len(h)) {
		Y_k <- y_future[[k]]
		ok_idx <- which(is.finite(Y_k))
		if (length(ok_idx) == 0L) next
		# extract per-draw forecast probabilities / means for the
		# observed cells at horizon k
		fc_k <- matrix(fc_full[, , k, ], nrow = length(Y_k))
		fc_obs <- fc_k[ok_idx, , drop = FALSE]
		y_obs  <- as.vector(Y_k)[ok_idx]

		pit_k <- switch(family,
			normal = {
				# gaussian PIT: Phi((y - mu) / sd), with sd combining draw spread and residual variance
				mu  <- rowMeans(fc_obs)
				sdr <- apply(fc_obs, 1, stats::sd)
				ve  <- if (!is.null(fit$VC)) mean(fit$VC[, ncol(fit$VC)],
				                                  na.rm = TRUE) else 1
				sdy <- sqrt(sdr^2 + ve)
				stats::pnorm(y_obs, mean = mu, sd = pmax(sdy, 1e-8))
			},
			tobit = {
				# tobit treated as gaussian on the latent scale
				mu  <- rowMeans(fc_obs)
				sdr <- apply(fc_obs, 1, stats::sd)
				ve  <- if (!is.null(fit$VC)) mean(fit$VC[, ncol(fit$VC)],
				                                  na.rm = TRUE) else 1
				sdy <- sqrt(sdr^2 + ve)
				stats::pnorm(y_obs, mean = mu, sd = pmax(sdy, 1e-8))
			},
			binary = {
				# randomised bernoulli PIT
				p_mean <- rowMeans(fc_obs)
				p_mean <- pmin(pmax(p_mean, 1e-12), 1 - 1e-12)
				u <- stats::runif(length(y_obs))
				ifelse(y_obs == 0, u * (1 - p_mean),
				       (1 - p_mean) + u * p_mean)
			},
			poisson = {
				# randomised poisson PIT (Czado, Gneiting & Held 2009)
				lam_mean <- rowMeans(fc_obs)
				lam_mean <- pmax(lam_mean, 1e-12)
				u <- stats::runif(length(y_obs))
				stats::ppois(y_obs - 1, lambda = lam_mean) +
				 u * stats::dpois(y_obs, lambda = lam_mean)
			},
			ordinal = {
				# ordinal scored like binary on the predicted probability
				p_mean <- rowMeans(fc_obs)
				p_mean <- pmin(pmax(p_mean, 1e-12), 1 - 1e-12)
				u <- stats::runif(length(y_obs))
				ifelse(y_obs == 0, u * (1 - p_mean),
				       (1 - p_mean) + u * p_mean)
			},
			# gaussian fallback for any unrecognised family
			{
				mu <- rowMeans(fc_obs)
				sdy <- apply(fc_obs, 1, stats::sd)
				stats::pnorm(y_obs, mean = mu, sd = pmax(sdy, 1e-8))
			}
		)
		pit_vec <- c(pit_vec, pit_k)
	}

	# KS test against Uniform(0, 1)
	ks <- suppressWarnings(stats::ks.test(pit_vec, "punif"))
	cover_95 <- mean(pit_vec >= 0.025 & pit_vec <= 0.975, na.rm = TRUE)

	out <- list(
		pit       = pit_vec,
		ks_stat   = unname(ks$statistic),
		ks_p      = unname(ks$p.value),
		cover_95  = cover_95,
		family    = family,
		h         = h,
		n_obs     = length(pit_vec)
	)
	class(out) <- c("forecast_pit", "list")
	out
}

#' @method print forecast_pit
#' @export
print.forecast_pit <- function(x, digits = 3, ...) {
	cli::cli_h2("Forecast PIT calibration check")
	cli::cli_ul(c(
		"Family: {.val {x$family}}",
		"Forecast horizon: {.val {x$h}}",
		"Held-out dyads: {.val {x$n_obs}}",
		"KS stat vs Uniform(0,1): {.val {signif(x$ks_stat, digits)}}",
		"KS p-value: {.val {signif(x$ks_p, digits)}}",
		"Fraction in [0.025, 0.975]: {.val {signif(x$cover_95, digits)}} (expect {.val {0.95}})"))
	if (is.finite(x$ks_p) && x$ks_p < 0.05) {
		cli::cli_text("KS rejects uniformity at 5%: the forecast distribution is mis-calibrated.")
	} else {
		cli::cli_text("PIT values are consistent with Uniform(0, 1) at the 5% level.")
	}
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
		ggplot2::geom_hline(yintercept = nrow(df) / 20,
		                    linetype = "dashed") +
		ggplot2::labs(
			x = "PIT Value",
			y = "Count",
			title = "Forecast Calibration Check",
			subtitle = sprintf("h = %d, family = %s, KS p = %.3f",
			                    x$h, x$family, x$ks_p)
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank()
		)
}
