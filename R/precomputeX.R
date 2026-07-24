#' Precompute design-array cross-product summaries
#'
#' Attaches to a dyadic covariate array the row/column margins and the
#' Gram matrices of its vectorised (and dyad-transposed) design, so that
#' the regression update inside the MCMC need not recompute them.
#'
#' @param X n x n x p array of dyadic covariates.
#' @return \code{X} unchanged in value but carrying the extra attributes
#'   \code{Xr} (n x p row margins), \code{Xc} (n x p column margins),
#'   \code{mX} (n^2 x p vectorised design), \code{mXt} (n^2 x p
#'   dyad-transposed vectorised design), \code{XX} (p x p Gram matrix
#'   \code{crossprod(mX)}) and \code{XXt} (p x p cross Gram matrix
#'   \code{crossprod(mX, mXt)}).
#' @keywords internal
#' @author lame authors
#' @export precomputeX
precomputeX <- function(X) {
	dims <- dim(X)
	n <- dims[1]
	p <- dims[3]
	slice_names <- dimnames(X)[[3]]

	# per-slice margins: sum a dyadic matrix over its columns (row margin)
	# and over its rows (column margin)
	row_margin <- apply(X, 3, rowSums)
	col_margin <- apply(X, 3, colSums)

	# stack each n-by-n slice column-major into one long design column
	design <- matrix(X, nrow = n * n, ncol = p)
	colnames(design) <- slice_names

	# the dyad transpose reorders the vectorised entries by the
	# commutation permutation vec(t(M)) = vec(M)[perm]
	perm <- as.vector(t(matrix(seq_len(n * n), n, n)))
	design_t <- design[perm, , drop = FALSE]

	gram <- crossprod(design)
	gram_cross <- crossprod(design, design_t)

	attr(X, "Xr") <- row_margin
	attr(X, "Xc") <- col_margin
	attr(X, "mX") <- design
	attr(X, "mXt") <- design_t
	attr(X, "XX") <- gram
	attr(X, "XXt") <- gram_cross

	X
}