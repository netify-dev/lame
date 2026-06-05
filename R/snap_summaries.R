#' Extract posterior draws of snap indices
#'
#' Computes draw-level averages of MCMC snap indicators over selected actors and
#' years. The result is useful for system-level rupture indices and focal-actor
#' rupture indices. Fits should be run with \code{keep_snap_draws = "draws"} for
#' posterior uncertainty. If only \code{snap_prob} is present, the function
#' returns a single score row and warns that intervals cannot be read as
#' posterior intervals.
#'
#' @param fit a \code{\link{lame}} fit with snap-shift MCMC output, or a
#'   \code{\link{lame_snap_als}} fit. MCMC fits should be run with retained
#'   snap draws for posterior uncertainty; ALS fits contribute their
#'   \code{snap_prob} scores as one score row.
#' @param actors optional actor names or indices. Default uses every actor on
#'   the selected side.
#' @param years optional year/period names or indices. Default uses every
#'   period.
#' @param side \code{"u"} for sender/row-side snap indicators or \code{"v"} for
#'   receiver/column-side indicators.
#' @return A data frame with columns \code{draw}, \code{year}, \code{side},
#'   \code{n_actors}, and \code{snap_index}.
#' @export
snap_index_draws <- function(fit, actors = NULL, years = NULL,
                             side = c("u", "v")) {
	side <- match.arg(side)
	obj <- .lame_snap_draw_array(fit, side)
	arr <- obj$array
	actor_idx <- .lame_snap_select_dim(actors, dimnames(arr)[[2L]],
	                                   dim(arr)[2L], "actors")
	year_idx <- .lame_snap_select_dim(years, dimnames(arr)[[3L]],
	                                  dim(arr)[3L], "years")
	mat <- .lame_snap_index_matrix(arr, actor_idx, year_idx)
	draw_names <- dimnames(arr)[[1L]]
	year_names <- dimnames(arr)[[3L]][year_idx]
	out <- data.frame(
		draw = rep(draw_names, times = length(year_idx)),
		year = rep(year_names, each = length(draw_names)),
		side = side,
		n_actors = length(actor_idx),
		snap_index = as.vector(mat),
		stringsAsFactors = FALSE)
	attr(out, "snap_draw_source") <- obj$source
	out
}

#' Summarize posterior snap indices
#'
#' @inheritParams snap_index_draws
#' @param probs numeric quantiles to report.
#' @return A data frame with one row per selected year.
#' @export
snap_index_summary <- function(fit, actors = NULL, years = NULL,
                               side = c("u", "v"),
                               probs = c(.025, .1, .5, .9, .975)) {
	side <- match.arg(side)
	d <- snap_index_draws(fit, actors = actors, years = years, side = side)
	years_out <- unique(d$year)
	rows <- lapply(years_out, function(yy) {
		x <- d$snap_index[d$year == yy]
		q <- stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
		if (all(!is.finite(x))) q[] <- NA_real_
		data.frame(
			year = yy,
			side = side,
			n_actors = d$n_actors[d$year == yy][1L],
			mean = mean(x, na.rm = TRUE),
			sd = stats::sd(x, na.rm = TRUE),
			t(q),
			check.names = FALSE,
			stringsAsFactors = FALSE)
	})
	out <- do.call(rbind, rows)
	q_names <- .lame_snap_prob_names(probs)
	names(out)[seq.int(ncol(out) - length(q_names) + 1L, ncol(out))] <- q_names
	out$mean[!is.finite(out$mean)] <- NA_real_
	out$sd[!is.finite(out$sd)] <- NA_real_
	attr(out, "snap_draw_source") <- attr(d, "snap_draw_source")
	out
}

#' Summarize snap indices by actor category
#'
#' @inheritParams snap_index_draws
#' @param groups a named list mapping group names to actor names/indices, or a
#'   named vector whose values are group labels and whose names are actors.
#' @return A data frame with one row per group-year.
#' @export
snap_category_summary <- function(fit, groups, years = NULL,
                                  side = c("u", "v")) {
	side <- match.arg(side)
	groups <- .lame_snap_normalize_groups(groups)
	rows <- lapply(names(groups), function(g) {
		out <- snap_index_summary(fit, actors = groups[[g]], years = years,
		                          side = side)
		out$group <- g
		out
	})
	out <- do.call(rbind, rows)
	out <- out[, c("group", setdiff(names(out), "group")), drop = FALSE]
	rownames(out) <- NULL
	out
}

#' Summarize posterior rank uncertainty for snap years
#'
#' Ranks year-level snap indices within each retained draw, with rank 1 assigned
#' to the highest snap index.
#'
#' @inheritParams snap_index_draws
#' @return A data frame with year-level index means and rank probabilities.
#' @export
snap_rank_summary <- function(fit, actors = NULL, years = NULL,
                              side = c("u", "v")) {
	side <- match.arg(side)
	obj <- .lame_snap_draw_array(fit, side)
	arr <- obj$array
	actor_idx <- .lame_snap_select_dim(actors, dimnames(arr)[[2L]],
	                                   dim(arr)[2L], "actors")
	year_idx <- .lame_snap_select_dim(years, dimnames(arr)[[3L]],
	                                  dim(arr)[3L], "years")
	mat <- .lame_snap_index_matrix(arr, actor_idx, year_idx)
	rank_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat),
	                   dimnames = dimnames(mat))
	for (ii in seq_len(nrow(mat))) {
		x <- mat[ii, ]
		if (!all(!is.finite(x))) {
			rank_mat[ii, ] <- rank(-x, ties.method = "min", na.last = "keep")
		}
	}
	year_names <- dimnames(arr)[[3L]][year_idx]
	top_k <- min(3L, length(year_idx))
	out <- data.frame(
		year = year_names,
		side = side,
		n_actors = length(actor_idx),
		index_mean = colMeans(mat, na.rm = TRUE),
		prob_rank_1 = colMeans(rank_mat == 1, na.rm = TRUE),
		prob_top_3 = colMeans(rank_mat <= top_k, na.rm = TRUE),
		rank_mean = colMeans(rank_mat, na.rm = TRUE),
		rank_median = apply(rank_mat, 2L, stats::median, na.rm = TRUE),
		stringsAsFactors = FALSE)
	out[] <- lapply(out, function(x) {
		if (is.numeric(x)) {
			x[!is.finite(x)] <- NA_real_
		}
		x
	})
	attr(out, "snap_draw_source") <- obj$source
	out
}

.lame_snap_draw_array <- function(fit, side) {
	draw_name <- if (identical(side, "u")) "snap_draws" else "snap_draws_v"
	prob_name <- if (identical(side, "u")) "snap_prob" else "snap_prob_v"
	arr <- fit[[draw_name]]
	if (!is.null(arr)) {
		if (length(dim(arr)) != 3L) {
			cli::cli_abort("{.code fit${draw_name}} must have dimensions draw x actor x time.")
		}
		arr <- .lame_snap_fill_dimnames(arr, draw_prefix = "draw")
		return(list(array = arr, source = "draws"))
	}
	prob <- fit[[prob_name]]
	if (is.null(prob)) {
		cli::cli_abort(c(
			"No snap output found for side {.val {side}}.",
			"i" = "For MCMC snap uncertainty, refit with {.code keep_snap_draws = \"draws\"}."))
	}
	is_als <- inherits(fit, "lame_snap_als")
	source <- if (is_als) "als_score" else "posterior_mean"
	draw_label <- if (is_als) "als_score" else "posterior_mean"
	cli::cli_warn(c(
		"No retained snap-indicator draws found on {.code fit${draw_name}}.",
		"i" = if (is_als) {
			"Using {.code fit${prob_name}} as one ALS score row; intervals and rank probabilities are not posterior uncertainty."
		} else {
			"Using {.code fit${prob_name}} as one posterior-mean row; intervals and rank probabilities are not posterior uncertainty."
		}))
	actor_names <- rownames(prob)
	if (is.null(actor_names)) actor_names <- paste0("actor", seq_len(nrow(prob)))
	time_names <- colnames(prob)
	if (is.null(time_names)) time_names <- paste0("t", seq_len(ncol(prob)))
	arr <- array(prob, dim = c(1L, nrow(prob), ncol(prob)))
	dimnames(arr) <- list(draw = draw_label,
	                      actor = actor_names,
	                      time = time_names)
	list(array = arr, source = source)
}

.lame_snap_fill_dimnames <- function(arr, draw_prefix = "draw") {
	dn <- dimnames(arr)
	if (is.null(dn)) dn <- vector("list", 3L)
	if (is.null(dn[[1L]])) dn[[1L]] <- paste0(draw_prefix, seq_len(dim(arr)[1L]))
	if (is.null(dn[[2L]])) dn[[2L]] <- paste0("actor", seq_len(dim(arr)[2L]))
	if (is.null(dn[[3L]])) dn[[3L]] <- paste0("t", seq_len(dim(arr)[3L]))
	dimnames(arr) <- dn
	arr
}

.lame_snap_select_dim <- function(x, dim_names, n, arg) {
	if (is.null(dim_names)) dim_names <- as.character(seq_len(n))
	if (is.null(x)) return(seq_len(n))
	if (is.character(x)) {
		if (is.null(dim_names)) {
			cli::cli_abort("{.arg {arg}} cannot be character because the snap array has no names.")
		}
		m <- match(x, dim_names)
		if (anyNA(m)) {
			cli::cli_abort("{.arg {arg}} not found: {.val {x[is.na(m)]}}.")
		}
		return(m)
	}
	if (!is.numeric(x) || anyNA(x) || any(x < 1L) || any(x > n)) {
		cli::cli_abort("{.arg {arg}} must be valid names or 1-based indices.")
	}
	as.integer(x)
}

.lame_snap_index_matrix <- function(arr, actor_idx, year_idx) {
	out <- matrix(NA_real_, nrow = dim(arr)[1L], ncol = length(year_idx))
	for (j in seq_along(year_idx)) {
		slice <- arr[, actor_idx, year_idx[j], drop = FALSE]
		dim(slice) <- c(dim(arr)[1L], length(actor_idx))
		out[, j] <- rowMeans(slice, na.rm = TRUE)
	}
	out[!is.finite(out)] <- NA_real_
	colnames(out) <- dimnames(arr)[[3L]][year_idx]
	rownames(out) <- dimnames(arr)[[1L]]
	out
}

.lame_snap_prob_names <- function(probs) {
	vals <- trimws(formatC(100 * probs, format = "fg", digits = 6))
	paste0("q", gsub("\\.", "_", vals))
}

.lame_snap_normalize_groups <- function(groups) {
	if (is.list(groups)) {
		if (is.null(names(groups)) || any(!nzchar(names(groups)))) {
			names(groups) <- paste0("group", seq_along(groups))
		}
		return(groups)
	}
	if (!is.null(names(groups))) {
		return(split(names(groups), groups))
	}
	cli::cli_abort("{.arg groups} must be a named list or a named vector of group labels.")
}
