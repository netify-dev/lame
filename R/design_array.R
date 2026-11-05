#' Assemble the dyadic design socioarray for an AME model
#'
#' Stacks the row, column and dyadic covariates into a single \eqn{n \times n
#' \times p} array in which slice \eqn{k} holds the dyad-level values of the
#' \eqn{k}th predictor. A sender attribute is broadcast down the rows, a
#' receiver attribute is broadcast across the columns, and a dyadic covariate
#' is copied verbatim. An intercept slice of ones is prepended unless it is
#' switched off or would duplicate an already-constant predictor.
#'
#' @usage design_array(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n,warn=TRUE)
#' @param Xrow an n x pr matrix of row (sender) covariates
#' @param Xcol an n x pc matrix of column (receiver) covariates
#' @param Xdyad an n x n x pd array of dyadic covariates
#' @param intercept logical; prepend a constant slice of ones
#' @param n number of rows/columns
#' @param warn logical; warn when missing covariate values are zero-filled
#'   (default \code{TRUE}). Set \code{FALSE} when the caller has already
#'   propagated the missingness into the response.
#' @return an n x n x (pr+pc+pd+intercept) 3-way array
#' @author lame authors
#' @keywords internal
#' @export design_array
design_array <- function(Xrow = NULL, Xcol = NULL, Xdyad = NULL, intercept = TRUE, n, warn = TRUE) {

	# coerce data frames first, since length() of a data frame counts columns
	if (!is.null(Xrow) && is.data.frame(Xrow)) Xrow <- as.matrix(Xrow)
	if (!is.null(Xcol) && is.data.frame(Xcol)) Xcol <- as.matrix(Xcol)

	# collect n x n covariate layers and their labels, then build the array once
	layers <- vector("list", 0L)
	labels <- character(0L)

	# row covariates: entry (i, j) carries actor i's value
	if (!is.null(Xrow)) {
		Xrow <- as.matrix(Xrow)
		pr <- ncol(Xrow)
		nms <- colnames(Xrow)
		if (is.null(nms)) nms <- paste0("Xrow", seq_len(pr))
		for (k in seq_len(pr)) {
			layers <- c(layers, list(matrix(Xrow[, k], n, n)))
		}
		labels <- c(labels, .lame_apply_suffix(nms, "row"))
	}

	# column covariates: entry (i, j) carries actor j's value
	if (!is.null(Xcol)) {
		Xcol <- as.matrix(Xcol)
		pc <- ncol(Xcol)
		nms <- colnames(Xcol)
		if (is.null(nms)) nms <- paste0("Xcol", seq_len(pc))
		for (k in seq_len(pc)) {
			layers <- c(layers, list(matrix(Xcol[, k], n, n, byrow = TRUE)))
		}
		labels <- c(labels, .lame_apply_suffix(nms, "col"))
	}

	# dyadic covariates enter as they are
	if (!is.null(Xdyad)) {
		dnms <- dimnames(Xdyad)[[3L]]
		Xd <- if (length(dim(Xdyad)) == 3L) Xdyad else array(Xdyad, dim = c(n, n, 1L))
		pd <- dim(Xd)[3L]
		if (is.null(dnms)) dnms <- paste0("Xdyad", seq_len(pd))
		for (k in seq_len(pd)) {
			layers <- c(layers, list(Xd[, , k]))
		}
		labels <- c(labels, .lame_apply_suffix(dnms, "dyad"))
	}

	p <- length(layers)
	if (p > 0L) {
		X <- array(unlist(layers, use.names = FALSE), dim = c(n, n, p))
	} else {
		X <- array(numeric(0L), dim = c(n, n, 0L))
	}

	# drop the intercept when a supplied covariate is already constant
	use_intercept <- isTRUE(intercept)
	if (use_intercept && p > 0L) {
		layer_var <- apply(X, 3L, function(s) var(as.vector(s), na.rm = TRUE))
		if (any(layer_var == 0, na.rm = TRUE)) use_intercept <- FALSE
	}

	if (use_intercept) {
		ones <- array(1, dim = c(n, n, 1L))
		if (p > 0L) {
			X <- array(c(as.vector(ones), as.vector(X)), dim = c(n, n, p + 1L))
		} else {
			X <- ones
		}
		labels <- c("intercept", labels)
	}

	# predictor names go on the third dimension only
	np <- dim(X)[3L]
	if (np >= 1L) dimnames(X)[[3L]] <- labels

	# zero-fill gaps, warning only for missing off-diagonal values
	if (anyNA(X)) {
		if (warn) {
			n_missing <- sum(is.na(X))
			n_diag <- sum(is.na(apply(X, 3L, diag)))
			if (n_missing > n_diag && getOption("lame.warn.na", TRUE)) {
				cli::cli_warn("Replacing NAs in design matrix with zeros")
			}
		}
		X[is.na(X)] <- 0
	}

	X
}

#' Propagate covariate missingness into the response
#'
#' For an AME fit, a dyad whose covariate value is missing cannot contribute a
#' covariate-coefficient observation. Rather than silently imputing the missing
#' covariate as 0 (which biases the coefficient when the covariate is not
#' mean-centred), the cell is treated as an unobserved tie and handled by data
#' augmentation. This helper sets the affected \code{Y} cells to \code{NA}: the
#' row of a sender whose covariate is missing, the column of a receiver whose
#' covariate is missing, and each cell with a missing dyadic covariate.
#' \code{\link{ame}} uses it for unipartite and bipartite networks;
#' \code{\link{lame}} does the same per period (see \code{\link{get_design_rep}}).
#'
#' @param Y an n_row x n_col response matrix.
#' @param Xrow,Xcol,Xdyad row, column and dyadic covariates (or \code{NULL}).
#'   A 2-D \code{Xdyad} is one n_row x n_col covariate.
#' @param symmetric logical; in a symmetric network a dyad has one covariate
#'   value, so a dyad is unobserved when the dyadic covariate of either of its
#'   two cells is missing (the node covariates already act on a row and its
#'   column together). \code{\link{ame}} and \code{\link{lame}} first fill a
#'   hole in one triangle from the other, so only a dyad missing in both
#'   triangles reaches this rule.
#' @return \code{Y} with covariate-missing cells set to \code{NA}.
#' @keywords internal
.ame_propagate_cov_na <- function(
	Y,
	Xrow = NULL,
	Xcol = NULL,
	Xdyad = NULL,
	symmetric = FALSE
) {
	if (!is.null(Xrow)) {
		Xr <- as.matrix(Xrow)
		bad <- which(!stats::complete.cases(Xr))
		if (length(bad)) Y[bad, ] <- NA
	}
	if (!is.null(Xcol)) {
		Xc <- as.matrix(Xcol)
		bad <- which(!stats::complete.cases(Xc))
		if (length(bad)) Y[, bad] <- NA
	}
	if (!is.null(Xdyad)) {
		Xd <- if (length(dim(Xdyad)) == 3L) Xdyad else
			array(Xdyad, dim = c(nrow(Y), ncol(Y), 1L))
		miss <- rowSums(is.na(Xd), dims = 2L) > 0
		if (isTRUE(symmetric) && nrow(Y) == ncol(Y)) miss <- miss | t(miss)
		Y[miss] <- NA
	}
	Y
}

# shared note for cells made unobserved by missing covariates
.lame_inform_cov_na <- function(n_cov_na) {
	if (!is.null(n_cov_na) && n_cov_na > 0L && getOption("lame.warn.na", TRUE)) {
		cli::cli_inform(c(
			"i" = "{n_cov_na} {.arg Y} cell{?s} had a missing covariate and {?is/are} treated as unobserved.",
			"i" = "The covariate coefficients are estimated from the complete dyads only."
		))
	}
	invisible(NULL)
}

# fill a dyadic covariate missing in one triangle from the other triangle;
# x is n x n x p, or n x n x p x t for a panel
.lame_fill_symmetric_dyadic <- function(x) {
	d <- dim(x)
	if (length(d) < 3L || d[1L] != d[2L]) return(x)
	flip <- aperm(x, c(2L, 1L, seq_along(d)[-(1:2)]))
	hole <- is.na(x) & !is.na(flip)
	n_fill <- sum(hole)
	if (n_fill > 0L) {
		x[hole] <- flip[hole]
		if (getOption("lame.warn.na", TRUE)) {
			cli::cli_inform(
				c(
					"i" = "{n_fill} {.arg Xdyad} value{?s} missing in one triangle {?was/were} filled from the other triangle.",
					"i" = "A symmetric dyad has one covariate value."
				),
				class = "lame_symmetric_xdyad_fill"
			)
		}
	}
	x
}

# refuse a dyadic covariate slice that is not symmetric under a symmetric fit;
# x is n x n x p, or n x n x p x t for a panel
.lame_check_symmetric_xdyad <- function(x) {
	d <- dim(x)
	if (length(d) < 3L || d[1L] != d[2L]) return(invisible(NULL))
	off <- abs(x - aperm(x, c(2L, 1L, seq_along(d)[-(1:2)])))
	bad <- !is.na(off) & off > 1e-8
	if (!any(bad)) return(invisible(NULL))
	asym <- which(apply(bad, 3L, any))
	slice_nms <- dimnames(x)[[3L]] %||% paste0("k=", seq_len(d[3L]))
	panel <- length(d) == 4L
	where <- ""
	if (panel) {
		first <- which(apply(bad, 4L, any))[1L]
		where <- paste0(" (first in period ", dimnames(x)[[4L]][first] %||% first, ")")
	}
	fix <- if (panel) {
		"e.g. {.code Xdyad[[t]][, , k] <- (Xdyad[[t]][, , k] + t(Xdyad[[t]][, , k])) / 2} in every period"
	} else {
		"e.g. {.code Xdyad[, , k] <- (Xdyad[, , k] + t(Xdyad[, , k])) / 2}"
	}
	cli::cli_abort(
		c(
			"{.arg symmetric} = TRUE, but {length(asym)} {.arg Xdyad} slice{?s} {.val {slice_nms[asym]}} {?is/are} not symmetric{where}.",
			"i" = paste0(
				"A symmetric fit needs a symmetric regression term. Replace each such slice by {.code (X + t(X)) / 2}, ",
				fix,
				", or set {.arg symmetric} = FALSE."
			)
		),
		class = "lame_symmetric_xdyad"
	)
}
