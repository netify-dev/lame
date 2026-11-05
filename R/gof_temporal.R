#' Posterior-predictive temporal-trend test
#'
#' For a fitted \code{lame} object, computes a network statistic
#' (\code{density}, \code{reciprocity}, or \code{transitivity}) at each
#' observed period, fits a least-squares linear trend on period index,
#' and compares the observed slope to slopes from posterior-predictive
#' replicates. The two-sided value is \code{min(1, 2 * min(p_up, p_down))},
#' where both empirical tails include ties. The value runs from 0 (observed slope in the extreme tail) to 1 (observed
#' slope central in the predictive distribution).
#'
#' @param fit A fitted \code{lame} object.
#' @param stat One of \code{"auto"} (default: density for binary, cbin and
#'   ordinal fits, otherwise mean), \code{"density"}, \code{"mean"},
#'   \code{"reciprocity"}, or \code{"transitivity"}. Density is the fraction
#'   of nonzero observed cells. Mean is the mean observed value. Unipartite
#'   checks exclude self ties; square bipartite checks retain diagonal cells.
#'   Reciprocity requires directed unipartite data; transitivity requires
#'   unipartite data. Gaussian density and nonzero transitivity are refused;
#'   use mean or weighted reciprocity for continuous outcomes.
#' @param n_rep Number of posterior-predictive replicates to draw (each
#'   replicate is a full \eqn{T}-period network). Must be an integer at least 2.
#' @param seed Optional RNG seed.
#'
#' @return A list with
#' \itemize{
#'   \item \code{stat} (chosen statistic name),
#'   \item \code{slope_obs} (observed slope of stat on period index),
#'   \item \code{slope_rep} (length-\code{n_rep} vector of replicate slopes),
#'   \item \code{p_pp} (two-sided posterior-predictive p-value),
#'   \item \code{stat_obs_by_t}, \code{stat_rep_by_t} (per-period statistics),
#'   \item \code{n_rep} (number of replicate networks drawn),
#'   \item \code{degenerate} (\code{TRUE} when the replicate slopes have no
#'   spread or a statistic is undefined, in which case \code{p_pp} is \code{NA}),
#'   \item \code{reason} and \code{n_undefined} describe an unestimated check,
#'   \item \code{replication_method} and \code{draw_index} identify joint draws
#' }
#'
#' @details
#' Observed and replicated statistics use the same actor labels and observed
#' cells in every period. Transitivity counts only fully observed triples.
#' Joint posterior draws are required; missing storage is refused. Fit with
#' \code{posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1,
#' thin_UV = 1)}, omitting latent storage at rank zero. A constant observed
#' series is still tested when the predictive slopes vary.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE, chains = 1)
#' gof_temporal(fit, stat = "density", n_rep = 50)
#' }
#'
#' @export
gof_temporal <- function(
	fit,
	stat = c("auto", "density", "mean", "reciprocity", "transitivity"),
	n_rep = 500,
	seed = NULL
) {
	stat <- match.arg(stat)
	if (!inherits(fit, "lame"))
		cli::cli_abort("{.arg fit} must be a fitted {.cls lame} object.")
	if (length(n_rep) != 1L || !is.finite(n_rep) || n_rep < 2L ||
			n_rep > .Machine$integer.max || n_rep != as.integer(n_rep))
		cli::cli_abort("n_rep must be an integer at least 2.")
	Y_obs <- .get_Y_list_for_gof(fit)
	n_time <- length(Y_obs)
	if (n_time < 3L) cli::cli_abort("{.fn gof_temporal} requires at least 3 time periods.")
	family <- tolower(fit$family %||% "normal")
	bip <- identical(fit$mode, "bipartite")
	if (stat == "auto") stat <- if (family %in% c("binary", "bin", "cbin", "ordinal"))
		"density" else "mean"
	if ((bip && stat %in% c("reciprocity", "transitivity")) ||
			(isTRUE(fit$symmetric) && stat == "reciprocity"))
		cli::cli_abort(
			"This network mode does not support the requested temporal statistic; use density or mean.",
			class = "lame_gof_uninformative"
		)
	if (family == "normal" && stat %in% c("density", "transitivity"))
		cli::cli_abort(
			"Nonzero-edge density and transitivity cannot check Gaussian predictive fit; use mean or directed reciprocity.",
			class = "lame_gof_uninformative"
		)
	stat_fun <- .pick_temporal_stat(stat, if (bip) "bipartite" else "unipartite")
	stat_obs_by_t <- vapply(Y_obs, stat_fun, numeric(1))
	period <- seq_len(n_time) - mean(seq_len(n_time))
	slope <- function(x) if (all(is.finite(x))) sum(period * x) / sum(period^2) else NA_real_
	slope_obs <- slope(stat_obs_by_t)

	sim_obj <- stats::simulate(fit, nsim = n_rep, seed = seed, type = "posterior")
	if (!is.list(sim_obj$Y) || length(sim_obj$Y) != n_rep)
		cli::cli_abort("Posterior simulation returned an unexpected number of networks.")
	stat_rep_by_t <- matrix(NA_real_, n_rep, n_time)
	for (r in seq_len(n_rep)) {
		ys <- .lame_gof_periods(sim_obj$Y[[r]])
		if (length(ys) != n_time) cli::cli_abort(
			"Predictive and observed GOF periods differ.",
			class = "lame_gof_alignment"
		)
		for (t in seq_len(n_time)) {
			y <- .lame_gof_align_replicate(ys[[t]], Y_obs[[t]])
			stat_rep_by_t[r, t] <- stat_fun(y)
		}
	}
	slope_rep <- apply(stat_rep_by_t, 1L, slope)
	n_undefined <- sum(!is.finite(stat_obs_by_t)) + sum(!is.finite(stat_rep_by_t))
	reason <- if (n_undefined) "The statistic is undefined in observed or replicated periods." else
	if (length(unique(slope_rep)) == 1L) "All predictive slopes are identical." else NULL
	degenerate <- !is.null(reason)
	p_pp <- if (degenerate) NA_real_ else
		min(1, 2 * min(mean(slope_rep >= slope_obs), mean(slope_rep <= slope_obs)))
	if (degenerate) cli::cli_warn(
		c(
			"Temporal-trend check is not estimated: {reason}",
			"i" = "Use a statistic defined on the observed cells, such as density or mean."
		),
		class = "lame_gof_uninformative"
	)
	res <- list(
		stat = stat,
		slope_obs = slope_obs,
		slope_rep = slope_rep,
		p_pp = p_pp,
		stat_obs_by_t = stat_obs_by_t,
		stat_rep_by_t = stat_rep_by_t,
		n_rep = n_rep,
		degenerate = degenerate,
		reason = reason,
		n_undefined = n_undefined,
		replication_method = sim_obj$method,
		draw_index = sim_obj$draw_index
	)
	class(res) <- c("gof_temporal", "list")
	res
}

# observed period matrices, keeping singleton axes and actor labels
.get_Y_list_for_gof <- function(fit) {
	.lame_gof_periods(fit$Y_obs %||% fit$Y)
}

# temporal statistics on observed cells and complete observed triples
.pick_temporal_stat <- function(stat, mode = "unipartite") {
	bip <- identical(mode, "bipartite")
	prepare <- function(Y) {
		Y[!is.finite(Y)] <- NA_real_
		if (!bip) diag(Y) <- NA_real_
		Y
	}
	switch(
		stat,
		density = function(Y) {
			Y <- prepare(Y)
			mean(Y[is.finite(Y)] != 0)
		},
		mean = function(Y) mean(prepare(Y), na.rm = TRUE),
		reciprocity = function(Y) {
			if (bip || nrow(Y) != ncol(Y)) return(NA_real_)
			Y <- prepare(Y)
			ut <- Y[upper.tri(Y)]; lt <- t(Y)[upper.tri(Y)]
			ok <- is.finite(ut) & is.finite(lt)
			if (sum(ok) < 3L || stats::sd(ut[ok]) == 0 || stats::sd(lt[ok]) == 0)
				return(NA_real_)
			stats::cor(ut[ok], lt[ok])
		},
		transitivity = function(Y) {
			if (bip || nrow(Y) != ncol(Y)) return(NA_real_)
			Y <- prepare(Y)
			observed <- is.finite(Y)
			A <- (observed & Y != 0) * 1
			paths <- A %*% A
			denom <- sum(paths * observed)
			if (denom <= 0) return(NA_real_)
			sum(paths * A) / denom
		}
	)
}

#' Print method for gof_temporal output
#' @param x A \code{gof_temporal} result.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   summary of the temporal-trend posterior-predictive check.
#' @export
#' @method print gof_temporal
print.gof_temporal <- function(x, ...) {
	.lame_cli_stdout(.lame_print_gof_temporal_body(x))
	invisible(x)
}

# cli body of print.gof_temporal, written to standard output
.lame_print_gof_temporal_body <- function(x) {
	cli::cli_h2("Temporal-trend posterior-predictive check")
	if (isTRUE(x$degenerate)) {
		cli::cli_ul(c(
			"Statistic: {.val {x$stat}}",
			"Status: {.val not estimated}",
			"Reason: {x$reason %||% 'Degenerate predictive statistic.'}",
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
	# two-sided p-value: high means well covered, low means incompatible
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
