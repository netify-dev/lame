# a fit that stores no posterior draws (max_seconds elapsed before the first stored draw): its draw
# arrays are empty along the draw axis, its means are NA, methods that summarise draws refuse
# with class lame_no_draws, and its start values are kept for lame_resume()

# the per-draw fields of a fit, by the position of the draw axis
.lame_draw_fields_first <- c(
	"BETA", "VC", "RHO_BETA", "SIGMA_BETA", "RHO", "log_lik",
	"snap_draws", "snap_draws_v", "THETA_ACTOR", "RHO_ACTOR", "SIGMA_ACTOR", "chain_summary",
	"rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap", "ab_draw_index", "uv_draw_index",
	"chain_indicator"
)
.lame_draw_fields_last <- c(
	"a_samples", "b_samples", "U_samples", "V_samples", "G_samples",
	"U_draws", "V_draws", "a_draws", "b_draws", "a_level_samples", "b_level_samples"
)
# the posterior means and standard deviations computed from the draws
.lame_draw_summary_fields <- c(
	"APM", "BPM", "APM_sd", "BPM_sd", "U", "V", "G", "UVPM",
	"EZ", "YPM", "a_level", "b_level", "theta_actor_mean", "theta_actor_sd"
)

# `x` with no entries along its draw axis ("first" or "last"); the other
# dimensions and their names are kept
.lame_zero_draws <- function(x, axis) {
	if (is.null(x) || !(is.numeric(x) || is.logical(x) || is.integer(x))) return(x)
	d <- dim(x)
	if (is.null(d)) return(x[0L])
	k <- if (identical(axis, "last")) length(d) else 1L
	idx <- lapply(d, seq_len)
	idx[[k]] <- integer(0)
	do.call(`[`, c(list(x), idx, list(drop = FALSE)))
}

# `x` (a vector, an array or a list of them) with every value NA
.lame_na_like <- function(x) {
	if (is.list(x) && !is.data.frame(x)) return(lapply(x, .lame_na_like))
	if (is.numeric(x) || is.logical(x)) x[] <- NA_real_
	x
}

# the fit of a run that stored no draws
.lame_empty_draws <- function(fit) {
	for (nm in intersect(.lame_draw_fields_first, names(fit))) fit[[nm]] <- .lame_zero_draws(fit[[nm]], "first")
	for (nm in intersect(.lame_draw_fields_last, names(fit))) fit[[nm]] <- .lame_zero_draws(fit[[nm]], "last")
	for (nm in intersect(.lame_draw_summary_fields, names(fit))) {
		if (!is.null(fit[[nm]])) fit[[nm]] <- .lame_na_like(fit[[nm]])
	}
	fs <- fit$forecast_state
	if (is.list(fs)) {
		for (nm in intersect(
			c("a", "b", "a_level", "b_level", "U", "V", "G"),
			names(fs)
		)) fs[[nm]] <- .lame_zero_draws(fs[[nm]], "last")
		if (!is.null(fs$draw_index)) fs$draw_index <- .lame_zero_draws(fs$draw_index, "first")
		fit$forecast_state <- fs
	}
	# the goodness-of-fit series keep the observed statistic, their first column
	if (is.list(fit$GOF)) {
		fit$GOF <- lapply(fit$GOF, function(g) if (is.matrix(g) && ncol(g) >= 1L) g[, 1L, drop = FALSE] else g)
	}
	fit$no_draws_stored <- TRUE
	fit
}

# whether a fit stores no posterior draws
.lame_has_no_draws <- function(fit) is.list(fit) && isTRUE(fit[["no_draws_stored"]])

# the refusal of a method that summarises draws
.lame_require_draws <- function(object, what) {
	if (!.lame_has_no_draws(object)) return(invisible(TRUE))
	cli::cli_abort(
		c(
			"The fit stores no posterior draws, so {.fn {what}} has nothing to summarise.",
			"i" = "{.arg max_seconds} elapsed before the first stored draw.",
			"i" = "Refit with a larger {.arg max_seconds} or a shorter {.arg burn}, or continue the chains with {.code lame_resume(fit, nscan_more = N, max_seconds = Inf)}."
		),
		class = "lame_no_draws",
		call = NULL
	)
}
