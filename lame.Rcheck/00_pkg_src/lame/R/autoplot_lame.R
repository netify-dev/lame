# autoplot.lame() ribbon plot for dynamic regression coefficients.

# re-export ggplot2's autoplot generic so plain autoplot(fit) dispatches.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# silence r cmd check globals for aes() bare-name references
utils::globalVariables(c("period", "med", "lo", "hi", "coef", "term",
                         "estimate", "conf.low", "conf.high"))

#' Ribbon plot of time-varying coefficients (or coefplot for static fits)
#'
#' For a \code{lame} fit with \code{dynamic_beta} on, returns a faceted
#' ggplot of the posterior mean coefficient path per period with a 95
#' percent credible-interval ribbon. For a static fit (no
#' \code{dynamic_beta}), falls back to a \code{tidy()}-driven horizontal
#' coefplot with posterior-mean point estimate and credible-interval
#' bars so that \code{autoplot(fit)} returns a ggplot regardless of fit
#' type.
#'
#' @param object A fitted \code{ame} / \code{lame} object.
#' @param which One of \code{"beta"} (default; coefficient plot --
#'   ribbon when dynamic, coefplot when static), \code{"ab"} (sender /
#'   receiver effects when \code{dynamic_ab}), \code{"uv"} (latent
#'   positions when \code{dynamic_uv}).
#' @param probs Length-3 vector of quantiles to plot. Default
#'   \code{c(0.025, 0.5, 0.975)} for 95 percent intervals.
#' @param coefs Optional character vector of coefficient names to subset.
#' @param ... Ignored.
#'
#' @return A \code{ggplot2} object that can be further customised.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   autoplot(fit)
#' }
#' }
#'
#' @method autoplot lame
#' @export
autoplot.lame <- function(object,
                          which = c("beta", "ab", "uv"),
                          probs = c(0.025, 0.5, 0.975),
                          coefs = NULL,
                          ...) {
	which <- match.arg(which)
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg ggplot2} is required for {.fn autoplot.lame}.")
	}

	if (which == "beta") {
		# static fit: fall back to a tidy()-driven coefplot rather than
		# erroring. dynamic_beta off is the common case for ame() fits and
		# for static lame() fits, so autoplot() should still return a ggplot.
		if (!isTRUE(object$dynamic_beta) ||
		    is.null(object$BETA) ||
		    length(dim(object$BETA)) != 3L) {
			return(.autoplot_static_coefplot(object, probs = probs, coefs = coefs))
		}
		# fit$beta is [iter x p x t] under dynamic_beta
		B <- object$BETA
		if (length(dim(B)) != 3L) {
			cli::cli_abort("{.code fit$BETA} must be a 3-D array for {.code which = \"beta\"}.")
		}
		dn <- dimnames(B)
		coef_nms <- if (!is.null(dn)) dn[[2L]] else paste0("v", seq_len(dim(B)[2L]))
		per_nms  <- if (!is.null(dn) && !is.null(dn[[3L]])) dn[[3L]] else paste0("t", seq_len(dim(B)[3L]))
		if (is.null(coefs)) coefs <- coef_nms

		dfs <- list()
		for (cf in coefs) {
			k <- match(cf, coef_nms)
			if (is.na(k)) next
			mat <- B[, k, , drop = FALSE]
			mat <- matrix(mat, nrow = dim(mat)[1L], ncol = dim(mat)[3L])
			q <- apply(mat, 2, stats::quantile, probs = probs, na.rm = TRUE)
			dfs[[cf]] <- data.frame(
				coef   = cf,
				period = seq_along(per_nms),
				period_label = per_nms,
				lo     = q[1L, ],
				med    = q[2L, ],
				hi     = q[3L, ],
				stringsAsFactors = FALSE
			)
		}
		df <- do.call(rbind, dfs)

		ggplot2::ggplot(df, ggplot2::aes(x = period, y = med)) +
			ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.25) +
			ggplot2::geom_line() +
			ggplot2::geom_point(size = 1.5) +
			ggplot2::facet_wrap(~ coef, scales = "free_y") +
			ggplot2::labs(
				x = "Period",
				y = "Coefficient",
				title = "Time-Varying Regression Coefficients",
				subtitle = sprintf("Posterior Median and %d%% Credible Interval",
				                   round(100 * (probs[3L] - probs[1L])))
			) +
			ggplot2::theme_bw() +
			ggplot2::theme(
				panel.border    = ggplot2::element_blank(),
				axis.ticks      = ggplot2::element_blank(),
				legend.position = "top",
				strip.background = ggplot2::element_rect(fill = "black", color = "black"),
				strip.text.x    = ggplot2::element_text(color = "white", hjust = 0),
				strip.text.y    = ggplot2::element_text(color = "white", hjust = 0)
			)
	} else if (which == "uv") {
		# delegate to uv_plot; it returns a ggplot and works for both static
		# and dynamic_uv fits. trajectory mode for longitudinal/dynamic
		# fits, snapshot for static (single time point).
		plot_type <- if (length(dim(object$U)) == 3L) "trajectory" else "snapshot"
		uv_plot(object, plot_type = plot_type)
	} else if (which == "ab") {
		# delegate to ab_plot; sender effect (a) by default. user can call
		# ab_plot() directly for the receiver-side plot. dynamic_ab fits
		# get a trajectory plot; static fits get the lollipop.
		is_dyn_ab <- isTRUE(object$dynamic_ab)
		plot_type <- if (is_dyn_ab) "trajectory" else "snapshot"
		ab_plot(object, effect = "sender", plot_type = plot_type)
	} else {
		cli::cli_abort("which = {.val {which}} is not implemented. Supported: {.val beta}, {.val uv}, {.val ab}.")
	}
}

# coefplot fallback used when which = "beta" is requested on a static
# (non-dynamic_beta) fit. builds a horizontal point-and-interval plot
# from tidy(fit) so that autoplot(fit) always returns a ggplot.
#' @noRd
.autoplot_static_coefplot <- function(object,
                                      probs = c(0.025, 0.5, 0.975),
                                      coefs = NULL) {
	conf_level <- max(probs) - min(probs)
	td <- tidy.ame(object, conf.int = TRUE, conf.level = conf_level)
	if (is.null(td) || nrow(td) == 0L) {
		cli::cli_abort("No coefficients to plot from {.fn tidy}.")
	}
	# tidy() on a dynamic fit yields one row per coef x period: collapse to
	# the first period for the static fallback case (we should not reach
	# here on a dynamic fit, but guard anyway)
	if ("period" %in% names(td)) {
		td <- td[td$period == td$period[1L], , drop = FALSE]
	}
	if (!is.null(coefs)) td <- td[td$term %in% coefs, , drop = FALSE]
	td$term <- factor(td$term, levels = rev(td$term))
	# pointrange renders horizontal error bars and the point in one geom.
	ggplot2::ggplot(td, ggplot2::aes(x = estimate, y = term)) +
		ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
		ggplot2::geom_pointrange(
			ggplot2::aes(xmin = conf.low, xmax = conf.high),
			size = 0.4) +
		ggplot2::labs(
			x = "Coefficient",
			y = NULL,
			title = "Posterior coefficient estimates",
			subtitle = sprintf("Posterior mean and %d%% credible interval",
			                   round(100 * conf_level))
		) +
		ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border    = ggplot2::element_blank(),
			axis.ticks      = ggplot2::element_blank(),
			legend.position = "top"
		)
}

# autoplot.ame: dispatch to autoplot.lame so that cross-sectional ame()
# fits also get a ggplot from the same generic. sharing the body means
# autoplot(fit) returns a coefplot regardless of which constructor
# (ame or lame) produced the fit.
#' @rdname autoplot.lame
#' @method autoplot ame
#' @export
autoplot.ame <- function(object,
                         which = c("beta", "ab", "uv"),
                         probs = c(0.025, 0.5, 0.975),
                         coefs = NULL,
                         ...) {
	autoplot.lame(object, which = which, probs = probs, coefs = coefs, ...)
}

#' autoplot method for ALS fits
#'
#' Coefficient point-and-interval plot for an \code{\link{ame_als}} or
#' \code{\link{lame_als}} fit. The body builds the same frame
#' \code{\link{autoplot.lame}} uses, sourcing the point estimate from
#' \code{coef(fit)} and the interval from \code{confint(fit)} (sandwich
#' or bootstrap, depending on which is available). \code{which =
#' "beta"} is the only mode supported; \code{which = "uv"} /
#' \code{"ab"} return an informative error.
#'
#' @param object A fitted \code{ame_als} / \code{lame_als} object.
#' @param which One of \code{"beta"} (currently the only supported
#'   value).
#' @param conf.level Confidence level for the interval. Default
#'   \code{0.95}.
#' @param ... Passed to \code{\link{vcov.ame_als}} / \code{\link{confint.ame_als}}
#'   (e.g. \code{cluster}).
#' @return A \code{ggplot} object.
#' @method autoplot ame_als
#' @export
autoplot.ame_als <- function(object,
                              which = c("beta", "uv", "ab"),
                              conf.level = 0.95, ...) {
	which <- match.arg(which)
	if (!identical(which, "beta")) {
		cli::cli_abort(c(
			"{.code which = \"{which}\"} is not supported on {.cls ame_als} fits.",
			"i" = "ALS produces a point estimate without posterior draws for {.code u} / {.code v} / {.code a} / {.code b}.",
			"i" = "Use the MCMC path ({.fn ame} / {.fn lame}) or {.code which = \"beta\"}."))
	}
	tdf <- tidy(object, conf.level = conf.level, ...)
	# require ggplot2 for the plot itself; mirrors autoplot.lame
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
		cli::cli_abort(c(
			"{.pkg ggplot2} is required for {.fn autoplot.ame_als}.",
			"i" = "Install it with {.code install.packages(\"ggplot2\")}."))
	}
	# colour encodes significance: an interval excluding 0 is significant,
	# signed by the point estimate; intervals covering 0 are not significant.
	sig <- !(tdf$conf.low <= 0 & tdf$conf.high >= 0)
	tdf$significance <- factor(
		ifelse(sig & tdf$estimate > 0, "Significant positive",
			ifelse(sig, "Significant negative", "Not significant")),
		levels = c("Significant positive", "Significant negative", "Not significant"))
	ggplot2::ggplot(tdf, ggplot2::aes(x = stats::reorder(.data$term, .data$estimate),
	                                    y = .data$estimate)) +
		ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$conf.low,
		                                       ymax = .data$conf.high,
		                                       colour = .data$significance)) +
		ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
		                     colour = "grey40") +
		ggplot2::scale_colour_manual(
			name = NULL,
			values = c("Significant positive" = "#0072B2",
			           "Significant negative" = "#D55E00",
			           "Not significant"      = "grey60")) +
		ggplot2::coord_flip() +
		ggplot2::labs(x = NULL, y = "Coefficient (sandwich/bootstrap CI)",
		               title = "ALS coefficient estimate") +
		ggplot2::theme_bw() +
		ggplot2::theme(panel.border = ggplot2::element_blank(),
		                axis.ticks      = ggplot2::element_blank(),
		                legend.position = "top")
}

#' @rdname autoplot.ame_als
#' @method autoplot lame_als
#' @export
autoplot.lame_als <- function(object,
                               which = c("beta", "uv", "ab"),
                               conf.level = 0.95, ...) {
	autoplot.ame_als(object, which = which, conf.level = conf.level, ...)
}

#' @noRd
.register_autoplot_method <- function() {
	if (requireNamespace("ggplot2", quietly = TRUE)) {
		# best-effort registration; ggplot2 may not have an `autoplot` generic
		# in extremely old installations, in which case our exported autoplot.lame
		# is still callable as `autoplot.lame(fit)`.
		try(registerS3method("autoplot", "lame", autoplot.lame,
		                     envir = asNamespace("ggplot2")),
		    silent = TRUE)
		try(registerS3method("autoplot", "ame", autoplot.ame,
		                     envir = asNamespace("ggplot2")),
		    silent = TRUE)
		try(registerS3method("autoplot", "ame_als", autoplot.ame_als,
		                     envir = asNamespace("ggplot2")),
		    silent = TRUE)
		try(registerS3method("autoplot", "lame_als", autoplot.lame_als,
		                     envir = asNamespace("ggplot2")),
		    silent = TRUE)
	}
}
