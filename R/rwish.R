#' Simulation from a Wishart distribution
#' 
#' Simulates a random Wishart-distributed matrix
#' 
#' 
#' @usage rwish(S0, nu = dim(S0)[1] + 2)
#' @param S0 a positive definite matrix
#' @param nu a positive integer
#' @return a positive definite matrix
#' @author Peter Hoff
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
rwish <-
function(S0,nu=dim(S0)[1]+2)
  {
    # Use C++ version only
    rwish_cpp(S0, nu)
  }