#' Simulate a binary relational matrix from latent means
#'
#' Draws a latent Gaussian array with the given dyadic correlation and
#' thresholds it at zero to produce a binary (0/1) sociomatrix. Self-ties on
#' the diagonal of a unipartite matrix are set to \code{NA};
#' rectangular matrices and explicit \code{bipartite = TRUE} inputs keep all cells.
#'
#' @usage simY_bin(EZ, rho, bipartite = FALSE)
#' @param EZ matrix giving the expected value of the latent \code{Z} matrix
#' @param rho dyadic correlation
#' @param bipartite whether row and column actors form disjoint sets; keeps
#'   diagonal cells even when the matrix is square.
#' @return a binary matrix matching the dimensions of \code{EZ}
#' @author lame authors
#' @keywords internal
#' @export simY_bin
simY_bin <- function(EZ, rho, bipartite = FALSE) {
	z_draw <- simZ(EZ, rho, bipartite = bipartite)
	y_out <- 1 * (z_draw > 0)
	if (!isTRUE(bipartite) && nrow(y_out) == ncol(y_out)) {
		diag(y_out) <- NA
	}
	y_out
}
