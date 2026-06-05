# posterior-predictive test for temporal trends. fits a linear time
# trend to a chosen network statistic and compares its slope to slopes
# from simulate(fit) replicates. p_pp is the two-sided gelman-style
# posterior-predictive p-value: 2 * min(p_up, 1 - p_up) where
# p_up = mean(slope_rep >= slope_obs). large p_pp -> observed slope is
# central in the predictive distribution (well covered); small p_pp ->
# observed slope is in either tail (incompatible).

#' Posterior-predictive temporal-trend test
#'
#' For a fitted \code{lame} object, computes a network statistic
#' (\code{density}, \code{reciprocity}, or \code{transitivity}) at each
#' observed period, fits a least-squares linear trend on period index,
#' and compares the observed slope to slopes from posterior-predictive
#' replicates. The two-sided value is \code{p_pp = 2 * min(p_up, 1 - p_up)},
#' which runs from 0 (observed slope in the extreme tail) to 1 (observed
#' slope dead-centre). A static fit on truly trending data yields
#' \code{p_pp} near 0; a dynamic fit that captures the trend yields
#' \code{p_pp} near 1.
#'
#' @param fit A fitted \code{lame} object.
#' @param stat One of \code{"auto"} (default; picks
#'   \code{"reciprocity"} for unipartite directed fits with \eqn{T \ge 3}
#'   and family \code{normal} / \code{poisson} / \code{tobit}, where
#'   \code{"density"} is constant and uninformative; otherwise
#'   \code{"density"}), \code{"density"}, \code{"mean"} (mean of
#'   off-diagonal Y; the right "density" analogue for continuous
#'   outcomes), \code{"reciprocity"}, \code{"transitivity"}.
#' @param n_rep Number of posterior-predictive replicates to draw (each
#'   replicate is a full \eqn{T}-period network from \code{simulate(fit)}).
#' @param seed Optional RNG seed.
#'
#' @return A list with
#' \itemize{
#'   \item \code{stat} (chosen statistic name),
#'   \item \code{slope_obs} (observed slope of stat on period index),
#'   \item \code{slope_rep} (length-\code{n_rep} vector of replicate slopes),
#'   \item \code{p_pp} (two-sided posterior-predictive p-value),
#'   \item \code{stat_obs_by_t}, \code{stat_rep_by_t} (per-period statistics)
#' }
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 200, burn = 50, odens = 5, verbose = FALSE)
#' gof_temporal(fit, stat = "density", n_rep = 100)
#' }
#'
#' @export
gof_temporal <- function(fit,
                         stat = c("auto", "density", "mean",
                                  "reciprocity", "transitivity"),
                         n_rep = 500,
                         seed = NULL) {
	stat <- match.arg(stat)
	if (!inherits(fit, "lame")) {
		cli::cli_abort("{.arg fit} must be a fitted {.cls lame} object.")
	}
	if (!is.null(seed)) set.seed(seed)

	# observed y: get a list of per-period matrices
	Y_obs <- .get_Y_list_for_gof(fit)
	T_periods <- length(Y_obs)
	if (T_periods < 3L) {
		cli::cli_abort("{.fn gof_temporal} requires at least 3 time periods.")
	}

	# resolve "auto" stat: pick something non-degenerate for the fit's
	# family. "density" measures the fraction of nonzero entries, which is
	# constant ~ 1 for continuous-outcome families and gives a spuriously
	# significant p_pp = 1. for those families default to "mean" (mean
	# off-diagonal value), which actually moves with the temporal trend.
	if (identical(stat, "auto")) {
		fam <- tolower(fit$family %||% "normal")
		stat <- switch(fam,
		               binary  = "density",
		               bin     = "density",
		               cbin    = "density",
		               ordinal = "density",
		               "mean")
	}

	stat_fun <- .pick_temporal_stat(stat)

	# observed slope
	stat_obs_by_t <- vapply(Y_obs, stat_fun, numeric(1))
	period <- seq_len(T_periods)

	# guard: when the chosen statistic is constant or near-constant across
	# observed periods, the slope is 0 and the p_pp value below is
	# meaningless. surface this with an informative warning + bail with a
	# clearly degenerate flag so the print method does not report a
	# false-positive "incompatible" alarm.
	if (!all(is.finite(stat_obs_by_t)) ||
	    stats::sd(stat_obs_by_t, na.rm = TRUE) < 1e-12) {
		cli::cli_warn(c(
			"Statistic {.val {stat}} is constant or near-constant across observed periods.",
			"i" = "The temporal-trend test is uninformative here.",
			"i" = "Try {.code stat = \"mean\"} (continuous outcomes) or {.code stat = \"reciprocity\"} (directed networks)."))
		res <- list(
			stat            = stat,
			slope_obs       = 0,
			slope_rep       = rep(0, n_rep),
			p_pp            = NA_real_,
			stat_obs_by_t   = stat_obs_by_t,
			stat_rep_by_t   = matrix(NA_real_, n_rep, T_periods),
			n_rep           = n_rep,
			degenerate      = TRUE
		)
		class(res) <- c("gof_temporal", "list")
		return(res)
	}
	slope_obs <- as.numeric(stats::coef(stats::lm(stat_obs_by_t ~ period))[2L])

	# posterior-predictive: draw n_rep replicate networks via simulate(fit).
	# simulate.lame returns an `lame.sim` object whose $y is a list of
	# length nsim, each element itself a list of t per-period matrices.
	sim_obj <- stats::simulate(fit, nsim = n_rep, seed = seed)
	sim_list <- if (is.list(sim_obj) && !is.null(sim_obj$Y)) sim_obj$Y else sim_obj
	if (!is.list(sim_list) || length(sim_list) < 1L) {
		cli::cli_abort("{.fn simulate.lame} did not return the expected list.")
	}
	# normalise: each element should be a list-of-t matrices
	get_t_slice <- function(rep, t) {
		if (is.list(rep)) return(rep[[t]])
		if (length(dim(rep)) == 3L) return(rep[, , t])
		cli::cli_abort("Unexpected {.fn simulate.lame} replicate format.")
	}

	stat_rep_by_t <- matrix(NA_real_, nrow = n_rep, ncol = T_periods)
	slope_rep <- numeric(n_rep)
	for (r in seq_len(n_rep)) {
		for (t in seq_len(T_periods)) {
			stat_rep_by_t[r, t] <- stat_fun(get_t_slice(sim_list[[r]], t))
		}
		slope_rep[r] <- as.numeric(
			stats::coef(stats::lm(stat_rep_by_t[r, ] ~ period))[2L])
	}

	# two-sided posterior-predictive p-value (gelman et al. style)
	# two-sided gelman-style posterior-predictive p-value: large p_pp
	# means well covered, small p_pp means in either tail.
	p_up <- mean(slope_rep >= slope_obs, na.rm = TRUE)
	p_pp <- 2 * min(p_up, 1 - p_up)

	res <- list(
		stat            = stat,
		slope_obs       = slope_obs,
		slope_rep       = slope_rep,
		p_pp            = p_pp,
		stat_obs_by_t   = stat_obs_by_t,
		stat_rep_by_t   = stat_rep_by_t,
		n_rep           = n_rep,
		degenerate      = FALSE
	)
	class(res) <- c("gof_temporal", "list")
	res
}

#' @noRd
.get_Y_list_for_gof <- function(fit) {
	# prefer y_obs (3-d array) over y (sometimes list, sometimes array)
	Y <- fit$Y_obs %||% fit$Y
	if (is.list(Y)) return(Y)
	if (is.array(Y) && length(dim(Y)) == 3L) {
		return(lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t]))
	}
	cli::cli_abort("Could not extract per-period Y matrices from the fit.")
}

#' @noRd
.pick_temporal_stat <- function(stat) {
	switch(stat,
		density = function(Y) {
			# fraction of nonzero off-diagonal entries.
			# meaningful only for binary / count-like networks; for
			# continuous outcomes mean(y != 0) is ~1 by construction.
			if (nrow(Y) == ncol(Y)) diag(Y) <- NA
			mean(Y != 0 & !is.na(Y), na.rm = TRUE)
		},
		mean = function(Y) {
			# mean off-diagonal value: the "density" analogue for
			# continuous outcomes.
			if (nrow(Y) == ncol(Y)) diag(Y) <- NA
			mean(Y, na.rm = TRUE)
		},
		reciprocity = function(Y) {
			# pearson correlation of y with t(y), off-diagonal only
			if (nrow(Y) != ncol(Y)) return(NA_real_)
			ut <- Y[upper.tri(Y)]
			lt <- t(Y)[upper.tri(Y)]
			ok <- is.finite(ut) & is.finite(lt)
			if (sum(ok) < 3L) return(NA_real_)
			suppressWarnings(stats::cor(ut[ok], lt[ok]))
		},
		transitivity = function(Y) {
			# fraction of length-2 paths that complete a triangle
			if (nrow(Y) != ncol(Y)) return(NA_real_)
			A <- (!is.na(Y) & Y != 0) * 1
			diag(A) <- 0
			paths <- A %*% A
			triangles <- paths * A
			denom <- sum(paths) - sum(diag(paths))
			if (denom <= 0) return(NA_real_)
			sum(triangles) / denom
		})
}

#' Print method for gof_temporal output
#' @param x A \code{gof_temporal} result.
#' @param ... Ignored.
#' @export
#' @method print gof_temporal
print.gof_temporal <- function(x, ...) {
	cli::cli_h2("Temporal-trend posterior-predictive check")
	if (isTRUE(x$degenerate)) {
		cli::cli_ul(c(
			"Statistic: {.val {x$stat}}",
			"Status: {.val degenerate} (statistic constant across periods)",
			"Replicates: {.val {x$n_rep}}",
			"Posterior-predictive p-value: {.val NA}"
		))
		cli::cli_text("Test is uninformative for this fit / statistic combination.")
		return(invisible(x))
	}
	cli::cli_ul(c(
		"Statistic: {.val {x$stat}}",
		"Observed slope (per-period): {.val {signif(x$slope_obs, 3)}}",
		"Replicates: {.val {x$n_rep}}",
		"Posterior-predictive p-value (two-sided): {.val {signif(x$p_pp, 3)}}"
	))
	# the two-sided p-value is small in the tails and large in the centre,
	# so "well covered" is high p_pp, "incompatible" is low p_pp.
	interp <- if (!is.finite(x$p_pp)) {
		"Posterior-predictive p-value is NA; the temporal-trend test is uninformative."
	} else if (x$p_pp < 0.05) {
		"Observed temporal trend is incompatible with the fitted model."
	} else if (x$p_pp < 0.2) {
		"Observed temporal trend is in the tails of the predictive distribution."
	} else {
		"Observed temporal trend is well covered by the fitted model."
	}
	cli::cli_text(interp)
	invisible(x)
}
