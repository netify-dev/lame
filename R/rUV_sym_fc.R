#' Gibbs sampling of U and V
#' 
#' A Gibbs sampler for updating the multiplicative effect matrices U and V
#' in the symmetric case. In this case \code{U\%*\%t(V)} is symmetric, so
#' this is parameterized as \code{V=U\%*\%L} where \code{L} is the 
#' diagonal matrix of eigenvalues of \code{U\%*\%t(V)}. 
#' 
#' @usage rUV_sym_fc(E, U, V, s2 = 1, shrink=TRUE)
#' @param E square residual relational matrix
#' @param U current value of U
#' @param V current value of V
#' @param s2 dyadic variance
#' @param shrink adaptively shrink the factors with a hierarchical prior
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @author Peter Hoff
#' @examples
#' 
#' U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
#' E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
#  rUV_sym_fc(E,U0,V0) 
#' rUV_sym_fc 
#' 
#' @export rUV_sym_fc
rUV_sym_fc<-function(E,U,V,s2=1,shrink=TRUE)
{
  n<-nrow(U)
  # Use C++ version only
  # Generate random loop IDs for the C++ version (C++ uses 0-based indexing)
  uLoopIDs <- as.integer(rep(sample(1:n),4) - 1)  # Convert to 0-based indexing
  rUV_sym_fc_cpp(E, U, V, s2, shrink, uLoopIDs)
}