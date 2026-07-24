#' Simulate a normal relational matrix
#'
#' Draws a Gaussian sociomatrix whose entries have expectation \code{EY},
#' dyadic variance \code{s2} and within-dyad correlation \code{rho}.
#'
#' @usage simY_nrm(EY, rho, s2)
#' @param EY square matrix giving the expected value of the relational matrix
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a square matrix
#' @author lame authors
#' @keywords internal
#' @export simY_nrm
simY_nrm <-
function(EY, rho, s2) {
	# For a normal outcome the observed relation coincides with its own
	# latent variable, so a draw is exactly a draw of Z about mean EY.
	y_draw <- simZ(EY, rho, s2)
	# One-mode (square) networks carry no self-ties; blank the diagonal.
	if (nrow(y_draw) == ncol(y_draw)) {
		diag(y_draw) <- NA
	}
	y_draw
}
