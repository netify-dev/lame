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
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export Xbeta
Xbeta <-
  function(X,beta)
  {
    if(length(beta) == 0) {
      return(matrix(0, nrow=dim(X)[1], ncol=dim(X)[2]))
    }
    Xbeta_cpp(X, beta)
  }