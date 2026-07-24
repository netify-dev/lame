#' Simulate a latent relational Gaussian array
#'
#' Draws a latent normal array from the social relations model: each
#' directed cell has mean \code{EZ} and marginal variance \code{s2}, and
#' the two cells of a dyad share within-dyad correlation \code{rho}.
#' Rectangular (bipartite) inputs carry no within-dyad correlation and are
#' drawn independently.
#'
#' @usage simZ(EZ, rho, s2 = 1)
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a simulated value of Z
#' @author lame authors
#' @keywords internal
#' @export simZ
simZ <- function(EZ, rho, s2 = 1) {
	nr <- nrow(EZ)
	nc <- ncol(EZ)
	sd <- sqrt(s2)

	if (nr != nc) {
		# rectangular: independent cells, no reciprocity structure
		return(EZ + matrix(rnorm(nr * nc, 0, sd), nr, nc))
	}

	# square case. build the correlated noise from the symmetric and
	# antisymmetric parts of an i.i.d. Gaussian matrix: the symmetric part
	# is perfectly correlated across a dyad, the antisymmetric part is
	# perfectly anti-correlated, and they are mutually independent. mixing
	# them with weights sqrt((1+rho)/2) and sqrt((1-rho)/2) gives unit
	# marginal variance and within-dyad correlation exactly rho.
	g    <- matrix(rnorm(nr * nr), nr, nr)
	sym  <- (g + t(g)) / sqrt(2)
	asym <- (g - t(g)) / sqrt(2)
	noise <- sqrt((1 + rho) / 2) * sym + sqrt((1 - rho) / 2) * asym

	EZ + sd * noise
}
