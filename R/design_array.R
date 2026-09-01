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
design_array <- function(Xrow=NULL, Xcol=NULL, Xdyad=NULL, intercept=TRUE, n, warn=TRUE) {

	# data.frames must be coerced before anything measures their size:
	# length(df) reports the column count, not the cell count, which would
	# corrupt every downstream dimension.
	if (!is.null(Xrow) && is.data.frame(Xrow)) Xrow <- as.matrix(Xrow)
	if (!is.null(Xcol) && is.data.frame(Xcol)) Xcol <- as.matrix(Xcol)

	# collect the covariate layers one at a time; each layer is a full n x n
	# matrix and gets a matching label. the array is materialised only once,
	# after every layer is known.
	layers <- vector("list", 0L)
	labels <- character(0L)

	# --- sender (row) effects: entry (i, j) carries actor i's attribute ---
	if (!is.null(Xrow)) {
		Xrow <- as.matrix(Xrow)
		pr   <- ncol(Xrow)
		nms  <- colnames(Xrow)
		if (is.null(nms)) nms <- paste0("Xrow", seq_len(pr))
		for (k in seq_len(pr)) {
			# recycling a length-n column down a matrix fills entry (i, j)
			# with Xrow[i, k] for every j
			layers <- c(layers, list(matrix(Xrow[, k], n, n)))
		}
		labels <- c(labels, .lame_apply_suffix(nms, "row"))
	}

	# --- receiver (col) effects: entry (i, j) carries actor j's attribute ---
	if (!is.null(Xcol)) {
		Xcol <- as.matrix(Xcol)
		pc   <- ncol(Xcol)
		nms  <- colnames(Xcol)
		if (is.null(nms)) nms <- paste0("Xcol", seq_len(pc))
		for (k in seq_len(pc)) {
			# filling by row puts Xcol[j, k] in entry (i, j) for every i
			layers <- c(layers, list(matrix(Xcol[, k], n, n, byrow = TRUE)))
		}
		labels <- c(labels, .lame_apply_suffix(nms, "col"))
	}

	# --- dyadic effects: the covariate surface is used as-is ---
	if (!is.null(Xdyad)) {
		dnms <- dimnames(Xdyad)[[3L]]
		Xd   <- if (length(dim(Xdyad)) == 3L) Xdyad else array(Xdyad, dim = c(n, n, 1L))
		pd   <- dim(Xd)[3L]
		if (is.null(dnms)) dnms <- paste0("Xdyad", seq_len(pd))
		for (k in seq_len(pd)) {
			layers <- c(layers, list(Xd[, , k]))
		}
		labels <- c(labels, .lame_apply_suffix(dnms, "dyad"))
	}

	# stack the collected layers into the socioarray
	p <- length(layers)
	if (p > 0L) {
		X <- array(unlist(layers, use.names = FALSE), dim = c(n, n, p))
	} else {
		X <- array(numeric(0L), dim = c(n, n, 0L))
	}

	# decide on the intercept: keep it unless the user declined or one of the
	# supplied covariates is already constant (which would make it redundant
	# and the design rank-deficient).
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

	# attach predictor names on the third dimension only (the first two
	# dimensions stay anonymous, matching the estimation routines' expectations)
	np <- dim(X)[3L]
	if (np >= 1L) dimnames(X)[[3L]] <- labels

	# zero-fill any remaining gaps. self-ties on the diagonal are structurally
	# absent, so only genuinely missing off-diagonal covariate values are worth
	# a warning.
	if (anyNA(X)) {
		if (warn) {
			n_missing <- sum(is.na(X))
			n_diag    <- sum(is.na(apply(X, 3L, diag)))
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
#' mean-centred), the dyad is treated as an unobserved tie and handled by data
#' augmentation. This helper sets the affected \code{Y} cells to \code{NA};
#' \code{\link{lame}} already does this internally, and this brings
#' \code{\link{ame}} into line.
#'
#' @param Y an n x n response matrix.
#' @param Xrow,Xcol,Xdyad row, column and dyadic covariates (or \code{NULL}).
#' @return \code{Y} with covariate-missing cells set to \code{NA}.
#' @keywords internal
.ame_propagate_cov_na <- function(Y, Xrow = NULL, Xcol = NULL, Xdyad = NULL) {
	n <- nrow(Y)
	if (!is.null(Xrow)) {
		Xr  <- as.matrix(Xrow)
		bad <- which(!stats::complete.cases(Xr))
		if (length(bad)) Y[bad, ] <- NA
	}
	if (!is.null(Xcol)) {
		Xc  <- as.matrix(Xcol)
		bad <- which(!stats::complete.cases(Xc))
		if (length(bad)) Y[, bad] <- NA
	}
	if (!is.null(Xdyad)) {
		Xd <- if (length(dim(Xdyad)) == 3L) Xdyad else
			array(Xdyad, dim = c(n, n, 1L))
		miss <- apply(Xd, c(1, 2), function(z) any(is.na(z)))
		Y[miss] <- NA
	}
	Y
}
