#' Simulation from a Wishart distribution
#'
#' Simulates a random Wishart-distributed matrix
#'
#'
#' @usage rwish(S0, nu = dim(S0)[1] + 2)
#' @param S0 a positive definite matrix
#' @param nu a positive integer
#' @return a positive definite matrix
#' @author lame authors
#' @keywords internal
#' @examples
#'
#' ## The expectation is S0*nu
#'
#' S0<-rwish(diag(3))
#'
#' SS<-matrix(0,3,3)
#' for(s in 1:1000) { SS<-SS+rwish(S0,5) }
#'
#' SS/s
#'
#' S0*5
#'
#'
#' @export rwish
rwish <- function(S0, nu = dim(S0)[1] + 2) {
	# Bartlett decomposition: draw a lower-triangular factor A whose
	# squared diagonal entries are chi-square with decreasing degrees
	# of freedom and whose strictly-lower entries are standard normal,
	# then rotate by the Cholesky factor of the scale matrix.
	p <- nrow(S0)
	lower_factor <- t(chol(S0))            # S0 = lower_factor %*% t(lower_factor)

	bartlett <- matrix(0, p, p)
	for (i in seq_len(p)) {
		bartlett[i, i] <- sqrt(rchisq(1, df = nu - i + 1))
	}
	if (p > 1) {
		below <- lower.tri(bartlett)
		bartlett[below] <- rnorm(sum(below))
	}

	root <- lower_factor %*% bartlett      # W = root %*% t(root) ~ Wishart(S0, nu)
	tcrossprod(root)
}
