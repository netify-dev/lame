#' Simulation from a multivariate normal distribution
#'
#' Simulates a matrix where the rows are i.i.d. samples from a multivariate
#' normal distribution. Positive semidefinite covariances are supported,
#' including degenerate distributions with zero-variance directions.
#'
#'
#' @usage rmvnorm(n, mu, Sigma, Sigma.chol = NULL)
#' @param n nonnegative integer sample size
#' @param mu finite numeric multivariate mean vector
#' @param Sigma finite symmetric positive semidefinite covariance matrix
#' @param Sigma.chol optional finite square factor whose crossproduct is
#'   \code{Sigma}; when supplied, this factor determines the covariance
#' @return a matrix with \code{n} rows
#' @examples
#' set.seed(1)
#' mu <- c(0, 1)
#' Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
#' draws <- rmvnorm(100, mu, Sigma)
#' colMeans(draws)
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rmvnorm
rmvnorm <- function(n, mu, Sigma, Sigma.chol = NULL) {
	p <- length(mu)
	if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
			n < 0 || n != floor(n))
		cli::cli_abort("n must be a nonnegative integer.", call = NULL)
	if (!is.numeric(mu) || any(!is.finite(mu)))
		cli::cli_abort("mu must be finite and numeric.", call = NULL)
	if (is.null(Sigma.chol)) {
		if (!is.matrix(Sigma) || !is.numeric(Sigma) ||
				!identical(dim(Sigma), c(p, p)) || any(!is.finite(Sigma)))
			cli::cli_abort("Sigma must be a finite covariance matrix matching mu.", call = NULL)
		if (p == 0L) return(matrix(numeric(), nrow = n, ncol = 0L))
		tolerance <- 8 * p * .Machine$double.eps * max(abs(Sigma))
		if (any(abs(Sigma - t(Sigma)) > tolerance))
			cli::cli_abort("Sigma must be symmetric.", call = NULL)
		Sigma.chol <- tryCatch(chol(Sigma), error = function(e) NULL)
		if (is.null(Sigma.chol)) {
			# singular covariance: factor through the eigendecomposition, zeroing only negative roundoff
			eS <- eigen(Sigma / 2 + t(Sigma) / 2, symmetric = TRUE)
			if (any(!is.finite(eS$values)))
				cli::cli_abort("Sigma could not be factored at its numerical scale.", call = NULL)
			tolerance <- 8 * p * .Machine$double.eps * max(abs(eS$values))
			if (any(eS$values < -tolerance))
				cli::cli_abort("Sigma must be positive semidefinite.", call = NULL)
			Sigma.chol <- t(sweep(eS$vectors, 2L, sqrt(pmax(eS$values, 0)), "*"))
		}
	} else if (!is.matrix(Sigma.chol) || !is.numeric(Sigma.chol) ||
			!identical(dim(Sigma.chol), c(p, p)) || any(!is.finite(Sigma.chol))) {
		cli::cli_abort("Sigma.chol must be a finite matrix matching mu.", call = NULL)
	}
	E <- matrix(rnorm(n * p), n, p)
	t(t(E %*% Sigma.chol) + c(mu))
}
