#' Precompute design matrix statistics
#'
#' Precomputes summary statistics of the design array X to speed up
#' MCMC calculations in the AME model
#'
#' @param X n x n x p covariate array
#' @return The input array with precomputed statistics as attributes:
#' \item{Xr}{row sums}
#' \item{Xc}{column sums}
#' \item{mX}{matricized version}
#' \item{mXt}{dyad-transposed matricized version}
#' \item{XX}{regression sums of squares}
#' \item{XXt}{crossproduct sums of squares}
#' @author Peter Hoff
#' @export precomputeX
precomputeX <- function(X) {
	Xr <- apply(X, c(1,3), sum)
	Xc <- apply(X, c(2,3), sum)
	mX <- apply(X, 3, c)
	mXt <- apply(aperm(X, c(2,1,3)), 3, c)
	XX <- t(mX) %*% mX
	XXt <- t(mX) %*% mXt

	attributes(X)$Xr <- Xr
	attributes(X)$Xc <- Xc
	attributes(X)$mX <- mX
	attributes(X)$mXt <- mXt
	attributes(X)$XX <- XX
	attributes(X)$XXt <- XXt

	X
}
