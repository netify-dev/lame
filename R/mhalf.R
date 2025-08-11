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
      tmp<-eigen(M)
      tmp$vec%*%sqrt(diag(tmp$val,nrow=nrow(M)))%*%t(tmp$vec)
    })
  }