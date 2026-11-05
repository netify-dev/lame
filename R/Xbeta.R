#' Linear combinations of submatrices of an array
#'
#' Computes a matrix of expected values based on an array X of predictors and a
#' vector beta of regression coefficients.
#'
#'
#' @usage Xbeta(X, beta)
#' @param X an n by n by p array
#' @param beta a p by 1 vector
#' @return An n by n matrix
#' @examples
#' set.seed(1)
#' n <- 12
#' X <- array(rnorm(n * n * 2), dim = c(n, n, 2))
#' beta <- c(0.5, -0.3)
#' EZ <- Xbeta(X, beta)
#' dim(EZ)
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export Xbeta
Xbeta <-
	function(X, beta) {
		d <- dim(X)
		if (length(d) != 3L) {
			cli::cli_abort("{.arg X} must be a 3D array.")
		}
		if (length(beta) == 0) {
			return(matrix(0, nrow = d[1], ncol = d[2]))
		}
		if (length(beta) != d[3L]) {
			cli::cli_abort(c(
				"{.arg beta} length must match the third dimension of {.arg X}.",
				"i" = "Got length(beta) = {length(beta)} and dim(X)[3] = {d[3L]}."
			))
		}
		Xbeta_cpp(X, beta)
	}
