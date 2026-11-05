# the snap kind of dynamic_uv has no ALS route: its drift scale would be read
# from paths the penalty has already shrunk. called by lame() for snap with
# method = "als" and by the exported lame_snap_als(), which only refuses
.lame_snap_als_refuse <- function(
	shape = c("directed", "bipartite", "symmetric"),
	caller = c("lame", "lame_snap_als")
) {
	shape <- match.arg(shape)
	caller <- match.arg(caller)
	reason <- switch(
		shape,
		directed = "On directed snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it either tracks the latent product less closely than the {.val ar1} kind or misses resets.",
		bipartite = "On bipartite snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it tracks the latent product less closely than the {.val ar1} kind.",
		symmetric = "On symmetric snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it tracks the latent product less closely than the {.val ar1} kind."
	)
	use <- if (identical(shape, "directed")) {
		"Use {.code dynamic_uv_kind = \"snap\"} with {.code method = \"mcmc\"}, or {.code dynamic_uv_kind = \"ar1\"} with {.code method = \"als\"}."
	} else {
		"Use {.code lame(..., dynamic_uv = TRUE, dynamic_uv_kind = \"ar1\")} with {.code method = \"als\"} or {.code method = \"mcmc\"}."
	}
	head <- if (identical(caller, "lame")) {
		"{.code dynamic_uv_kind = \"snap\"} with {.code method = \"als\"} is not available."
	} else {
		"{.fn lame_snap_als} was removed in lame 1.4.0; the snap kind has no ALS estimator."
	}
	cli::cli_abort(
		c(
			head,
			"i" = reason,
			"i" = use
		),
		call = parent.frame()
	)
}

#' Removed: fast snap-shift ALS estimator
#'
#' @description
#' \code{lame_snap_als()} fitted a fast point estimate of the snap-shift
#' dynamic latent-factor model. It was removed in lame 1.4.0 and now stops with
#' a message: on simulated snap panels its per-period latent product tracked
#' the generating one less closely than the \code{"ar1"} kind of dynamic ALS
#' (correlation 0.72 against 0.80 on directed panels, 0.24 against 0.67 on
#' bipartite and 0.44 against 0.72 on symmetric ones), and holding its drift at
#' the \code{"ar1"} kind's smoothing, which closed that gap, missed planted
#' resets (both ranked first in 15 of 20 panels against 20 of 20).
#'
#' Fit snap-shift positions with \code{\link{lame}} and
#' \code{dynamic_uv = TRUE, dynamic_uv_kind = "snap"} (the MCMC, directed
#' panels). For bipartite and symmetric panels, or for a point estimate, use
#' \code{dynamic_uv_kind = "ar1"} with \code{method = "als"} or the MCMC.
#'
#' @param ... ignored: every call stops with the message above. A call with
#'   \code{mode = "bipartite"}, \code{R_row} / \code{R_col}, or
#'   \code{symmetric = TRUE} gets the message for that kind of panel.
#' @return Does not return.
#' @seealso \code{\link{lame}}
#' @export
lame_snap_als <- function(...) {
	args <- list(...)
	shape <- if (identical(args[["mode"]], "bipartite") ||
			!is.null(args[["R_row"]]) || !is.null(args[["R_col"]])) {
		"bipartite"
	} else if (isTRUE(args[["symmetric"]])) {
		"symmetric"
	} else {
		"directed"
	}
	.lame_snap_als_refuse(shape, caller = "lame_snap_als")
}
