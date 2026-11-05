# canonical regression design for newdata and counterfactual prediction
#
# coefficient layout: [intercept?, <name>_row..., <name>_col..., <name>_dyad...].
# newdata supplies the dyadic slices; intercept and nodal slices are held at
# their fitted values from the stored design

#' Map a coefficient vector onto a design array's slice layout
#'
#' symmetric fits store one merged nodal coefficient per covariate (the
#' sampler averages the row and column coefficients at storage time), so
#' \code{colnames(fit$BETA)} is shorter than \code{dim(fit$X)[3]} and a
#' positional contraction of beta against the design pairs coefficients with
#' the wrong slices. this returns a coefficient per design slice, matched by
#' name, with the nodal coefficient \code{<v>_node} used for both
#' \code{<v>_row} and \code{<v>_col}. a vector that already matches the
#' design slice for slice is returned unchanged.
#'
#' @param beta named numeric coefficient vector.
#' @param design_names \code{dimnames(X)[[3]]} of the design array.
#' @return a numeric vector of length \code{length(design_names)}.
#' @keywords internal
.beta_on_design <- function(beta, design_names) {
	if (is.null(design_names) || is.null(names(beta))) return(beta)
	if (length(beta) == length(design_names)) return(beta)
	out <- stats::setNames(numeric(length(design_names)), design_names)
	# a slice takes the coefficient of its own name, or the merged <v>_node
	# coefficient when row and column coefficients are pooled
	hit <- match(design_names, names(beta))
	key <- sub("_(row|col)$", "_node", design_names)
	miss <- is.na(hit)
	hit[miss] <- match(key[miss], names(beta))
	ok <- !is.na(hit)
	if (!any(ok)) {
		cli::cli_abort(c(
			"None of the design slices match the coefficient names.",
			"i" = "Design slices: {.val {design_names}}; coefficients: {.val {names(beta)}}."
		))
	}
	out[ok] <- beta[hit[ok]]
	out
}

#' Build the full canonical design for a newdata prediction
#'
#' @param newdata an \code{n x m x n_dyad} array (or \code{n x m} matrix) of
#'   new dyadic covariates. If its 3rd-dim names match the model's dyadic
#'   coefficient names they are mapped by name; otherwise the slices are
#'   taken in order as the model's dyadic covariates.
#' @param model_names the model's coefficient names
#'   (\code{colnames(fit$BETA)} / \code{dimnames(fit$BETA)[[2]]}), in
#'   canonical layout.
#' @param fitted_design the model's fitted design array (\code{fit$X} or
#'   \code{Xlist[[t]]}); used to hold intercept + nodal slices fixed. May be
#'   \code{NULL} only when the model has no nodal covariates.
#' @param n,m output dimensions.
#' @return an \code{n x m x length(model_names)} array whose slices are in
#'   \code{model_names} order: dyadic slices from \code{newdata}, all other
#'   slices held at their fitted values (intercept = 1).
#' @keywords internal
.build_full_design <- function(newdata, model_names, fitted_design, n, m) {
	if (is.null(model_names)) {
		cli::cli_abort(c(
			"Cannot reconcile {.arg newdata} against an unnamed coefficient vector.",
			"i" = "The fit must carry {.code colnames(BETA)} so dyadic vs nodal covariates can be told apart."
		))
	}
	if (length(dim(newdata)) == 2L) newdata <- array(newdata, c(dim(newdata), 1L))
	p <- length(model_names)
	is_dyad <- grepl("_dyad$|\\.dyad$", model_names)
	dyad_pos <- which(is_dyad)
	n_dyad <- length(dyad_pos)

	if (dim(newdata)[3] != n_dyad) {
		cli::cli_abort(c(
			"{.arg newdata} has {dim(newdata)[3]} dyadic slice{?s} but the model has {n_dyad} dyadic coefficient{?s}{if (n_dyad > 0) paste0(' (', paste(model_names[dyad_pos], collapse = ', '), ')') else ''}.",
			"i" = "Supply one slice per dyadic covariate; the intercept and any row/col effects are held at their fitted values."
		))
	}

	out <- array(0, dim = c(n, m, p), dimnames = list(NULL, NULL, model_names))

	# intercept and nodal slices are held at fitted values
	nondyad_pos <- which(!is_dyad)
	# match slices to coefficients by name (own name first, then the merged
	# <v>_node key), since a symmetric fit pools <v>_row and <v>_col
	design_names <- if (!is.null(fitted_design)) dimnames(fitted_design)[[3]] else NULL
	has_design <- !is.null(fitted_design) && length(dim(fitted_design)) == 3L &&
		!is.null(design_names)
	nodal_pos <- nondyad_pos[model_names[nondyad_pos] != "intercept"]
	if (length(nodal_pos) > 0L && !has_design) {
		cli::cli_abort(c(
			"Cannot apply {.arg newdata} to a fit with nodal (row/col) covariates without the stored design.",
			"i" = "The intercept + nodal contributions must be held at their fitted values, which requires the named design array stored on the fit.",
			"i" = "Refit without dropping the design (avoid compaction / {.code use_sparse_matrices}), or predict without {.arg newdata}."
		))
	}
	design_key <- if (has_design) sub("_(row|col)$", "_node", design_names) else character(0)
	for (k in nondyad_pos) {
		if (identical(model_names[k], "intercept")) {
			out[, , k] <- 1
		} else {
			src <- which(design_names == model_names[k])
			if (length(src) == 0L) src <- which(design_key == model_names[k])
			if (length(src) == 0L) {
				cli::cli_abort(c(
					"No design slice matches the coefficient {.val {model_names[k]}}.",
					"i" = "The stored design's slice names are {.val {design_names}}."
				))
			}
			out[, , k] <- Reduce(`+`, lapply(src, function(j) fitted_design[, , j]))
		}
	}

	# dyadic slices: from newdata, mapped by name when possible
	if (n_dyad > 0L) {
		nd_names <- dimnames(newdata)[[3]]
		dyad_names <- model_names[dyad_pos]
		if (!is.null(nd_names) && all(dyad_names %in% nd_names)) {
			for (k in dyad_pos) out[, , k] <- newdata[, , match(model_names[k], nd_names)]
		} else {
			# positional: newdata dyad slices taken in the model's dyad order
			for (jj in seq_along(dyad_pos)) out[, , dyad_pos[jj]] <- newdata[, , jj]
		}
	}
	out
}

#' Name-safe linear predictor from a newdata prediction
#'
#' Convenience wrapper: builds the full canonical design via
#' \code{\link{.build_full_design}} and contracts it against a named
#' coefficient vector. Because the design is rebuilt in the model's slice
#' order, the positional contraction is correct even with nodal covariates.
#'
#' @param beta named numeric coefficient vector (a single draw or the
#'   posterior mean).
#' @param newdata dyadic covariate array (see \code{.build_full_design}).
#' @param fitted_design fitted design array (\code{fit$X}); held-fixed source
#'   for intercept + nodal slices.
#' @param n,m output dimensions.
#' @return an \code{n x m} matrix of the linear predictor.
#' @keywords internal
.xbeta_newdata <- function(beta, newdata, fitted_design, n, m) {
	Xfull <- .build_full_design(newdata, names(beta), fitted_design, n, m)
	XB <- matrix(0, n, m)
	for (k in seq_len(dim(Xfull)[3])) XB <- XB + beta[k] * Xfull[, , k]
	XB
}

#' The fit's internal (canonical) actor order
#'
#' \code{lame()} sorts actors alphabetically when it ingests a list of \code{Y}
#' matrices (via \code{list_to_array}), so the stored additive / multiplicative
#' effects and design slices are in sorted actor order, which need not match the
#' order a user later supplies in \code{newdata}. This returns the row and column
#' actor names the fit was estimated in, tried from the most authoritative source
#' down, so \code{newdata} can be realigned by name before prediction.
#'
#' @param object a fitted \code{lame} (or \code{ame}) object.
#' @return a list with \code{rows} and \code{cols} character vectors (or
#'   \code{NULL} when the order cannot be recovered).
#' @keywords internal
.fit_actor_order <- function(object) {
	# lame fits store the per-period designs as Xlist, ame fits as X
	fdx <- object[["Xlist", exact = TRUE]] %||% object[["X", exact = TRUE]]
	fd1 <- if (is.list(fdx)) fdx[[1L]] else fdx
	rows <- NULL; cols <- NULL
	if (!is.null(fd1) && length(dim(fd1)) >= 2L) {
		dn <- dimnames(fd1); rows <- dn[[1L]]; cols <- dn[[2L]]
	}
	if (is.null(rows)) {
		rows <- names(object$APM) %||% rownames(object$U) %||%
			rownames(object$a_dynamic)
	}
	if (is.null(cols)) {
		cols <- names(object$BPM) %||% rownames(object[["V"]]) %||%
			rownames(object$b_dynamic)
	}
	if (isTRUE(object$symmetric) && !identical(object$mode, "bipartite")) cols <- cols %||% rows
	list(rows = rows, cols = cols)
}

#' Realign a newdata covariate array to the fit's internal actor order
#'
#' Reorders the rows and columns of each \code{n x m x p} slice by their
#' dimnames so they line up with the fit's stored (sorted) additive /
#' multiplicative effects. Only acts when \code{newdata} carries dimnames whose
#' name set exactly matches the fit's; otherwise it falls back to the supplied
#' positional order unchanged (the historical behaviour, correct when the user
#' already passes actors in the fit's order or supplies no names).
#'
#' @param X an \code{n x m x p} array (already promoted from a matrix).
#' @param row_order,col_order canonical actor names from
#'   \code{\link{.fit_actor_order}}.
#' @return \code{X} with rows / columns permuted into canonical order, or
#'   unchanged.
#' @keywords internal
.reorder_newdata_actors <- function(X, row_order, col_order) {
	if (is.null(row_order) && is.null(col_order)) return(X)
	dn <- dimnames(X)
	rn <- if (length(dn) >= 1L) dn[[1L]] else NULL
	cn <- if (length(dn) >= 2L) dn[[2L]] else NULL
	ri <- seq_len(dim(X)[1L]); ci <- seq_len(dim(X)[2L])
	if (!is.null(row_order) && !is.null(rn) &&
			length(rn) == length(row_order) && setequal(rn, row_order)) {
		ri <- match(row_order, rn)
	}
	if (!is.null(col_order) && !is.null(cn) &&
			length(cn) == length(col_order) && setequal(cn, col_order)) {
		ci <- match(col_order, cn)
	}
	if (identical(ri, seq_len(dim(X)[1L])) &&
			identical(ci, seq_len(dim(X)[2L]))) return(X)
	X[ri, ci, , drop = FALSE]
}
