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
  # Try C++ version first if available
  if(exists("precomputeX_cpp", mode = "function")) {
    tryCatch({
      result <- precomputeX_cpp(X)
      attributes(X)$Xr <- result$Xr
      attributes(X)$Xc <- result$Xc
      attributes(X)$mX <- result$mX
      attributes(X)$mXt <- result$mXt
      attributes(X)$XX <- result$XX
      attributes(X)$XXt <- result$XXt
      return(X)
    }, error = function(e) {
      # Fall back to R implementation
    })
  }
  
  # R implementation
  Xr <- apply(X, c(1,3), sum)           # row sum
  Xc <- apply(X, c(2,3), sum)           # col sum
  mX <- apply(X, 3, c)                  # design matrix
  mXt <- apply(aperm(X, c(2,1,3)), 3, c) # dyad-transposed design matrix
  XX <- t(mX) %*% mX                    # regression sums of squares
  XXt <- t(mX) %*% mXt                   # crossproduct sums of squares
  
  attributes(X)$Xr <- Xr
  attributes(X)$Xc <- Xc
  attributes(X)$mX <- mX
  attributes(X)$mXt <- mXt
  attributes(X)$XX <- XX
  attributes(X)$XXt <- XXt
  
  X
}