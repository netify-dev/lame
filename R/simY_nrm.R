#' Simulate a normal relational matrix
#'
#' Draws a Gaussian sociomatrix whose entries have expectation \code{EY},
#' dyadic variance \code{s2} and within-dyad correlation \code{rho}.
#'
#' @usage simY_nrm(EY, rho, s2, bipartite = FALSE)
#' @param EY matrix giving the expected value of the relational matrix
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param bipartite whether row and column actors form disjoint sets; keeps
#'   diagonal cells even when the matrix is square.
#' @return a matrix with the dimensions of \code{EY}
#' @author lame authors
#' @keywords internal
#' @export simY_nrm
simY_nrm <-
	function(EY, rho, s2, bipartite = FALSE) {
		# a normal outcome is its own latent variable, so draw Z about EY
		y_draw <- simZ(EY, rho, s2, bipartite = bipartite)
		# one-mode networks have no self-ties
		if (!isTRUE(bipartite) && nrow(y_draw) == ncol(y_draw)) {
			diag(y_draw) <- NA
		}
		y_draw
	}
