#' Simulate a binary relational matrix from latent means
#'
#' Draws a latent Gaussian array with the given dyadic correlation and
#' thresholds it at zero to produce a binary (0/1) sociomatrix. Self-ties on
#' the diagonal of a square (unipartite) matrix are set to \code{NA};
#' rectangular (bipartite) matrices are returned with every entry populated.
#'
#' @usage simY_bin(EZ, rho)
#' @param EZ matrix giving the expected value of the latent \code{Z} matrix
#' @param rho dyadic correlation
#' @return a binary matrix matching the dimensions of \code{EZ}
#' @author lame authors
#' @keywords internal
#' @export simY_bin
simY_bin <- function(EZ, rho) {
	z_draw <- simZ(EZ, rho)
	y_out <- matrix(0, nrow(z_draw), ncol(z_draw), dimnames = dimnames(z_draw))
	y_out[z_draw > 0] <- 1
	if (nrow(y_out) == ncol(y_out)) {
		diag(y_out) <- NA
	}
	y_out
}
