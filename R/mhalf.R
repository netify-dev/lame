#' Symmetric square root of a matrix
#' 
#' Computes the symmetric square root of a positive definite matrix
#' 
#' 
#' @usage mhalf(M)
#' @param M a positive definite matrix
#' @return a matrix \code{H} such that \code{H^2} equals \code{M}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export mhalf
mhalf <-
  function(M) 
  { 
    tryCatch({
      return(mhalf_cpp(M))
    }, error = function(e) {
      # Ensure symmetry first
      M <- (M + t(M))/2
      tmp <- eigen(M, symmetric = TRUE)
      # Floor eigenvalues at small positive value for stability
      D <- pmax(tmp$values, 1e-12)
      tmp$vectors %*% diag(sqrt(D), nrow(M)) %*% t(tmp$vectors)
    })
  }