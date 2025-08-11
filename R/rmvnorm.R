#' Simulation from a multivariate normal distribution
#' 
#' Simulates a matrix where the rows are i.i.d. samples from a multivariate
#' normal distribution
#' 
#' 
#' @usage rmvnorm(n, mu, Sigma, Sigma.chol = NULL)
#' @param n sample size
#' @param mu multivariate mean vector
#' @param Sigma covariance matrix
#' @param Sigma.chol Cholesky factorization of \code{Sigma}
#' @return a matrix with \code{n} rows
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rmvnorm
rmvnorm <-
  function(n,mu,Sigma,Sigma.chol=NULL)
  {
    tryCatch({
      return(rmvnorm_cpp(n, mu, Sigma))
    }, error = function(e) {
      p <- length(mu)
      
      if(any(!is.finite(Sigma))) {
        Sigma.chol <- diag(p)
      } else if(is.null(Sigma.chol)) {
        tryCatch({
          Sigma.chol <- chol(Sigma)
        }, error = function(e) {
          tryCatch({
            eS <- eigen(Sigma, symmetric=TRUE)
            eS$values[eS$values < 1e-10] <- 1e-10
            Sigma.chol <- t(eS$vectors %*% diag(sqrt(eS$values)))
          }, error = function(e2) {
            dS <- diag(Sigma)
            if(any(!is.finite(dS)) || any(dS <= 0)) {
              Sigma.chol <- diag(p) * 0.001
            } else {
              Sigma.chol <- diag(sqrt(pmax(dS, 1e-10)))
            }
          })
        })
      }
      E<-matrix(rnorm(n*p),n,p)
      t(  t(E%*%Sigma.chol) +c(mu))
    })
  }