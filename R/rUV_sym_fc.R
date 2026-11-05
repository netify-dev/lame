#' Gibbs sampling of U and V
#'
#' A Gibbs sampler for updating the multiplicative effect matrices U and V
#' in the symmetric case. In this case \code{U\%*\%t(V)} is symmetric, so
#' this is parameterized as \code{V=U\%*\%L} where \code{L} is the
#' diagonal matrix of eigenvalues of \code{U\%*\%t(V)}.
#'
#' @usage rUV_sym_fc(E, U, V, s2 = 1, shrink=TRUE, l_var = 1)
#' @param E square residual relational matrix
#' @param U current value of U
#' @param V current value of V
#' @param s2 dyadic variance
#' @param shrink adaptively shrink the factors with a hierarchical prior
#' @param l_var prior variance of an eigen-weight: 1 for a unit response (the
#'   default), the variance of the response for a normal response of free scale
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#'
#' U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2))
#' E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30)
#' rUV_sym_fc(E, U0, V0)
#'
#' @export rUV_sym_fc
rUV_sym_fc <- function(E, U, V, s2 = 1, shrink = TRUE, l_var = 1) {
	n <- nrow(U)
	# random update order, 0-based for C++
	uLoopIDs <- as.integer(rep(sample(1:n), 4) - 1)
	rUV_sym_fc_cpp(E, U, V, s2, shrink, uLoopIDs, l_var)
}
