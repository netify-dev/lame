# canonical regression design for newdata and counterfactual prediction
#
# the coefficient vector has a structured layout:
# [intercept?, <name>_row..., <name>_col..., <name>_dyad...].
# newdata supplies dyadic covariate slices; intercept and nodal slices are
# held at fitted values from the stored design.

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
			"i" = "The fit must carry {.code colnames(BETA)} so dyadic vs nodal covariates can be told apart."))
	}
	if (length(dim(newdata)) == 2L) newdata <- array(newdata, c(dim(newdata), 1L))
	p <- length(model_names)
	is_dyad <- grepl("_dyad$|\\.dyad$", model_names)
	dyad_pos <- which(is_dyad)
	n_dyad <- length(dyad_pos)

	if (dim(newdata)[3] != n_dyad) {
		cli::cli_abort(c(
			"{.arg newdata} has {dim(newdata)[3]} dyadic slice{?s} but the model has {n_dyad} dyadic coefficient{?s}{if (n_dyad > 0) paste0(' (', paste(model_names[dyad_pos], collapse = ', '), ')') else ''}.",
			"i" = "Supply one slice per dyadic covariate; the intercept and any row/col effects are held at their fitted values."))
	}

	out <- array(0, dim = c(n, m, p), dimnames = list(NULL, NULL, model_names))

	# intercept + nodal slices: held at fitted values
	nondyad_pos <- which(!is_dyad)
	has_design <- !is.null(fitted_design) && length(dim(fitted_design)) == 3L &&
		dim(fitted_design)[3] == p
	nodal_pos <- nondyad_pos[model_names[nondyad_pos] != "intercept"]
	if (length(nodal_pos) > 0L && !has_design) {
		cli::cli_abort(c(
			"Cannot apply {.arg newdata} to a fit with nodal (row/col) covariates without the stored design.",
			"i" = "The intercept + nodal contributions must be held at their fitted values, which requires {.code fit$X}.",
			"i" = "Refit without dropping the design (avoid compaction / {.code use_sparse_matrices}), or predict without {.arg newdata}."))
	}
	for (k in nondyad_pos) {
		if (identical(model_names[k], "intercept")) {
			out[, , k] <- 1
		} else {
			out[, , k] <- fitted_design[, , k]
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
	fdx <- object$X
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
		cols <- names(object$BPM) %||% rownames(object$V) %||%
			rownames(object$b_dynamic)
	}
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
