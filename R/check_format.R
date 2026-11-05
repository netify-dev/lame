#' Validate input data format for lame function
#'
#' Internal validation function that checks the format and consistency of input
#' data for longitudinal AME models. Ensures that network data and covariates
#' are properly formatted as lists with consistent dimensions across time periods.
#'
#' @usage check_format(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL)
#' @param Y a list of T network matrices, where T is the number of time periods.
#' Each element should be an n x m matrix representing the network at time t.
#' @param Xdyad an optional list of T dyadic covariates. Each element can be
#' either an n x m matrix (single covariate) or an n x m x p array (p covariates).
#' Must have the same length as Y if provided.
#' @param Xrow an optional list of T matrices of row/sender covariates.
#' Each element should be an n x pr matrix where pr is the number of row covariates.
#' Must have the same length as Y if provided.
#' @param Xcol an optional list of T matrices of column/receiver covariates.
#' Each element should be an m x pc matrix where pc is the number of column covariates.
#' Must have the same length as Y if provided.
#'
#' @details
#' Validates input data for longitudinal network analysis:
#' \itemize{
#'   \item Verifying Y is a non-empty list of matrices
#'   \item Checking that all covariates (if provided) are lists of appropriate length
#'   \item Validating data types for all elements
#'   \item Warning about dimension inconsistencies across time periods
#' }
#'
#' The function uses informative error messages via the cli package to help
#' users identify and correct data formatting issues.
#'
#' @return Invisible TRUE if all checks pass. Throws an error with an informative
#' message if any validation fails.
#'
#' @note This is an internal function primarily used by \code{lame()} but exported
#' for advanced users who want to validate their data before model fitting.
#'
#' @examples
#' n <- 12
#' Y <- lapply(1:3, function(t) matrix(rnorm(n * n), n, n))
#' check_format(Y)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#'
#' @export check_format

check_format <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL) {

	# validate Y
	if (!is.list(Y)) {
		cli::cli_abort(c(
			"Y must be a list of matrices.",
			"i" = "Received object of class {class(Y)[1]}."
		))
	}

	if (length(Y) == 0) {
		cli::cli_abort("Y must contain at least one matrix.")
	}

	for (t in seq_along(Y)) {
		if (!is.matrix(Y[[t]]) && !is.array(Y[[t]])) {
			cli::cli_abort(c(
				"All elements of Y must be matrices.",
				"i" = "Y[[{t}]] is of class {class(Y[[t]])[1]}."
			))
		}
	}

	N <- length(Y)

	# validate Xdyad
	if (!is.null(Xdyad)) {
		if (!is.list(Xdyad)) {
			cli::cli_abort(c(
				"Xdyad must be NULL or a list of matrices/arrays.",
				"i" = "Received object of class {class(Xdyad)[1]}."
			))
		}
		if (length(Xdyad) != N) {
			cli::cli_abort(c(
				"Xdyad must have the same length as Y.",
				"i" = "Y has {N} time periods, but Xdyad has {length(Xdyad)}."
			))
		}
		for (t in seq_along(Xdyad)) {
			# every period carries the same dyadic covariates
			if (is.null(Xdyad[[t]])) {
				cli::cli_abort(c(
					"Xdyad[[{t}]] is NULL; every period needs its dyadic covariate array.",
					"i" = "Supply an array for that period; a cell with an NA covariate is treated as an unobserved tie."
				))
			}
			if (!is.matrix(Xdyad[[t]]) && !is.array(Xdyad[[t]])) {
				cli::cli_abort(c(
					"All elements of Xdyad must be matrices or arrays.",
					"i" = "Xdyad[[{t}]] is of class {class(Xdyad[[t]])[1]}."
				))
			}
		}
	}

	# validate Xrow
	if (!is.null(Xrow)) {
		if (!is.list(Xrow)) {
			cli::cli_abort(c(
				"Xrow must be NULL or a list of matrices.",
				"i" = "Received object of class {class(Xrow)[1]}."
			))
		}
		if (length(Xrow) != N) {
			cli::cli_abort(c(
				"Xrow must have the same length as Y.",
				"i" = "Y has {N} time periods, but Xrow has {length(Xrow)}."
			))
		}
		for (t in seq_along(Xrow)) {
			if (!is.null(Xrow[[t]]) && !is.matrix(Xrow[[t]]) && !is.data.frame(Xrow[[t]])) {
				cli::cli_abort(c(
					"All non-NULL elements of Xrow must be matrices or data frames.",
					"i" = "Xrow[[{t}]] is of class {class(Xrow[[t]])[1]}."
				))
			}
		}
	}

	# validate Xcol
	if (!is.null(Xcol)) {
		if (!is.list(Xcol)) {
			cli::cli_abort(c(
				"Xcol must be NULL or a list of matrices.",
				"i" = "Received object of class {class(Xcol)[1]}."
			))
		}
		if (length(Xcol) != N) {
			cli::cli_abort(c(
				"Xcol must have the same length as Y.",
				"i" = "Y has {N} time periods, but Xcol has {length(Xcol)}."
			))
		}
		for (t in seq_along(Xcol)) {
			if (!is.null(Xcol[[t]]) && !is.matrix(Xcol[[t]]) && !is.data.frame(Xcol[[t]])) {
				cli::cli_abort(c(
					"All non-NULL elements of Xcol must be matrices or data frames.",
					"i" = "Xcol[[{t}]] is of class {class(Xcol[[t]])[1]}."
				))
			}
		}
	}

	# covariate dimensions must match Y
	for (t in seq_along(Y)) {
		nr <- nrow(Y[[t]])
		nc <- ncol(Y[[t]])
		if (!is.null(Xdyad) && !is.null(Xdyad[[t]])) {
			xd <- Xdyad[[t]]
			if (nrow(xd) != nr || ncol(xd) != nc) {
				cli::cli_abort(c(
					"Xdyad dimensions must match Y at each time period.",
					"i" = "Y[[{t}]] is {nr} x {nc}, but Xdyad[[{t}]] is {nrow(xd)} x {ncol(xd)}."
				))
			}
		}
		if (!is.null(Xrow) && !is.null(Xrow[[t]])) {
			xr <- Xrow[[t]]
			if (is.matrix(xr) || is.data.frame(xr)) {
				if (nrow(xr) != nr) {
					cli::cli_abort(c(
						"Xrow row count must match Y row count at each time period.",
						"i" = "Y[[{t}]] has {nr} rows, but Xrow[[{t}]] has {nrow(xr)} rows."
					))
				}
			}
		}
		if (!is.null(Xcol) && !is.null(Xcol[[t]])) {
			xc <- Xcol[[t]]
			if (is.matrix(xc) || is.data.frame(xc)) {
				if (nrow(xc) != nc) {
					cli::cli_abort(c(
						"Xcol row count must match Y column count at each time period.",
						"i" = "Y[[{t}]] has {nc} columns, but Xcol[[{t}]] has {nrow(xc)} rows."
					))
				}
			}
		}
	}

	# covariate counts must agree across periods: list_to_array() sizes the array
	# from period 1 and fills every period positionally
	.cov_count <- function(x, dyadic) {
		d <- dim(x)
		if (dyadic) {
			if (length(d) == 3L) d[3L] else 1L
		} else {
			if (length(d) >= 2L) d[2L] else 1L
		}
	}
	.cov_names <- function(x, dyadic) {
		d <- dim(x) ; dn <- dimnames(x)
		if (dyadic) {
			if (length(d) == 3L && !is.null(dn)) dn[[3L]] else NULL
		} else {
			if (!is.null(dn)) dn[[2L]] else NULL
		}
	}
	for (nm in c("Xdyad", "Xrow", "Xcol")) {
		lst <- switch(nm, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)
		if (is.null(lst)) next
		idx <- which(!vapply(lst, is.null, logical(1)))
		if (length(idx) < 2L) next
		dyadic <- identical(nm, "Xdyad")
		cnt <- vapply(lst[idx], .cov_count, numeric(1), dyadic = dyadic)
		bad <- which(cnt != cnt[1L])
		if (length(bad) > 0L) {
			cli::cli_abort(c(
				"{.arg {nm}} must carry the same number of covariates in every period.",
				"i" = "Period {idx[1]} has {cnt[1]}, but period {idx[bad[1]]} has {cnt[bad[1]]}.",
				"i" = "Covariates are matched positionally across periods; pad an unobserved covariate with {.val NA} rather than dropping it."
			))
		}
		nms <- lapply(lst[idx], .cov_names, dyadic = dyadic)
		ref <- nms[[1L]]
		if (!is.null(ref)) {
			bad <- which(
				!vapply(
					nms,
					function(z) is.null(z) || identical(z, ref),
					logical(1)
				)
			)
			if (length(bad) > 0L) {
				cli::cli_warn(c(
					"{.arg {nm}} covariate names differ across periods.",
					"i" = "Period {idx[1]}: {.val {ref}}; period {idx[bad[1]]}: {.val {nms[[bad[1]]]}}.",
					"i" = "Covariates are matched positionally, so the fit uses the period-{idx[1]} names throughout."
				))
			}
		}
	}

	# dimensions must agree across the periods of Y
	if (N > 1) {
		first_nrow <- nrow(Y[[1]])
		first_ncol <- ncol(Y[[1]])
		for (t in 2:N) {
			if (nrow(Y[[t]]) != first_nrow || ncol(Y[[t]]) != first_ncol) {
				cli::cli_inform(c(
					"Network dimensions vary across time periods.",
					"i" = "Y[[1]] has dimensions {first_nrow} x {first_ncol}.",
					"i" = "Y[[{t}]] has dimensions {nrow(Y[[t]])} x {ncol(Y[[t]])}.",
					"i" = "This is allowed but actors will be aligned by names."
				))
				break
			}
		}
	}

	invisible(TRUE)
}

# covariates far from unit scale only prompt a suggestion to rescale; returns the
# slice standard deviations. `X` is an [n, m, p] or [n, m, p, T] array or a list of [n, m, p] arrays
.lame_check_covariate_scale <- function(X, inform_sd = 1e2) {
	if (is.null(X)) return(invisible(NULL))
	slices <- if (is.list(X)) {
		nm <- dimnames(X[[1]])[[3]]
		lapply(seq_len(dim(X[[1]])[3]), function(k) unlist(lapply(X, function(a) a[, , k])))
	} else if (length(dim(X)) == 4L) {
		nm <- dimnames(X)[[3]]
		lapply(seq_len(dim(X)[3]), function(k) as.numeric(X[, , k, ]))
	} else if (length(dim(X)) == 3L) {
		nm <- dimnames(X)[[3]]
		lapply(seq_len(dim(X)[3]), function(k) as.numeric(X[, , k]))
	} else return(invisible(NULL))
	if (is.null(nm)) nm <- paste0("slice ", seq_along(slices))
	sds <- vapply(slices, function(v) stats::sd(v, na.rm = TRUE), numeric(1))
	sds[!is.finite(sds)] <- 0
	large <- which(sds > inform_sd)
	if (length(large) > 0L) {
		large_txt <- paste0(nm[large], " (sd ", signif(sds[large], 2), ")", collapse = ", ")
		cli::cli_inform(c(
			"i" = "{cli::qty(length(large))}Covariate{?s} with a large standard deviation: {large_txt}; the fit does not depend on the units, but rescaling ({.fn scale} or {.fn log}) makes the coefficients easier to read."
		))
	}
	invisible(sds)
}

# the ALS route refuses covariates with standard deviation above `abort_sd`, since
# the solvers work on the covariates' own scale (binary IRLS stops converging near 1e5)
.lame_als_check_covariate_scale <- function(Xdyad = NULL, Xrow = NULL, Xcol = NULL, abort_sd = 1e4) {
	slice_sd <- function(x) {
		if (is.null(x)) return(numeric(0))
		if (is.list(x)) {
			parts <- lapply(x, slice_sd)
			k <- max(lengths(parts), 0L)
			if (k == 0L) return(numeric(0))
			return(vapply(
				seq_len(k),
				function(j) max(vapply(
					parts,
					function(p)
						if (length(p) >= j) p[j] else 0,
					numeric(1)
				)),
				numeric(1)
			))
		}
		d <- dim(x)
		if (is.null(d)) return(stats::sd(as.numeric(x), na.rm = TRUE))
		k_dim <- if (length(d) >= 3L) 3L else 2L
		vapply(seq_len(d[k_dim]), function(j) {
			v <- if (length(d) >= 3L) as.numeric(x[slice.index(x, 3L) == j]) else as.numeric(x[, j])
			stats::sd(v, na.rm = TRUE)
		}, numeric(1))
	}
	sds <- c(
		Xdyad = max(c(slice_sd(Xdyad), 0), na.rm = TRUE),
		Xrow = max(c(slice_sd(Xrow), 0), na.rm = TRUE),
		Xcol = max(c(slice_sd(Xcol), 0), na.rm = TRUE)
	)
	big <- names(sds)[is.finite(sds) & sds > abort_sd]
	if (length(big)) {
		cli::cli_abort(
			c(
				"{.arg {big}} contain{?s/} a covariate with a standard deviation above {abort_sd}.",
				"i" = "The ALS route does not converge reliably at that scale; rescale the covariate (for example with {.fn scale} or {.fn log}), or use {.code method = \"mcmc\"}, which is exact at any scale."
			),
			class = "lame_als_covariate_scale"
		)
	}
	invisible(sds)
}

# note for a directed fit whose observed dyads are all exactly symmetric
# (symmetric = TRUE fits the undirected model)
.lame_inform_symmetric_y <- function(dcor, panel = FALSE, family = "binary") {
	property <- if (isTRUE(dcor)) {
		# the nominations of a rank family only draw rho toward one
		paste0(
			"The residuals of reciprocal cells are equal, so the posterior of the dyadic correlation {.code rho} ",
			if (identical(family, "frn")) "is drawn toward one." else "concentrates near one."
		)
	} else {
		"The two equal cells of each dyad count as two observations."
	}
	what <- if (isTRUE(panel)) "Every period of {.arg Y} is exactly symmetric" else "{.arg Y} is exactly symmetric"
	cli::cli_inform(
		c(
			"i" = paste0(what, " and {.arg symmetric} = FALSE."),
			"i" = property,
			"i" = "{.code symmetric = TRUE} fits the undirected model."
		),
		class = "lame_symmetric_y_note"
	)
	invisible(NULL)
}
