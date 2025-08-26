#' Simulate Z given its expectation and covariance
#' 
#' Simulate Z given its expectation and covariance
#' 
#' 
#' @usage simZ(EZ, rho, s2 = 1)
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a simulated value of Z
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export simZ
simZ <-
function(EZ,rho,s2=1)
{ 
  tryCatch({
    return(simZ_cpp(EZ, rho, s2))
  }, error = function(e) {
    # Handle both square and rectangular matrices
    nr <- nrow(EZ)
    nc <- ncol(EZ)
    
    if(nr == nc) {
      # Square matrix - original code for unipartite
      c<-(sqrt(1+rho) + sqrt(1-rho))/2
      d<-(sqrt(1+rho) - sqrt(1-rho))/2
      EC<-matrix(rnorm(length(EZ)), nr, nc)
      EC<- sqrt(s2)*( c*EC + d*t(EC) )
      EZ+EC
    } else {
      # Rectangular matrix - bipartite, no dyadic correlation
      # For bipartite, dyadic correlation doesn't apply since i->j and j->i don't exist
      EC <- matrix(rnorm(length(EZ), 0, sqrt(s2)), nr, nc)
      EZ + EC
    }
  })
}