#' Gibbs update for multiplicative effects covariance
#' 
#' Gibbs sampling for the covariance matrix of multiplicative effects U and V
#' in the AME model. This function implements the inverse-Wishart posterior 
#' update for the covariance matrix of the stacked UV effects.
#' 
#' @usage rSuv_fc(U, V, Suv0=NULL, kappa0=NULL)
#' @param U matrix of multiplicative row effects (n x R matrix where n is the
#'   number of nodes and R is the dimension of multiplicative effects)
#' @param V matrix of multiplicative column effects (n x R matrix)
#' @param Suv0 prior (inverse) scale matrix for the prior distribution. 
#'   Default is identity matrix of dimension 2R x 2R, providing a weakly 
#'   informative prior.
#' @param kappa0 prior degrees of freedom for the prior distribution.
#'   Default is 2 + 2R, which is the minimum for a proper prior with a 
#'   2R x 2R covariance matrix.
#' @return Updated covariance matrix Suv (2R x 2R matrix) for the stacked
#'   effects \\[U, V\\]. The first R x R block contains covariances for U,
#'   the last R x R block contains covariances for V, and the off-diagonal
#'   blocks contain cross-covariances between U and V.
#' @details
#' The function updates the full covariance matrix for multiplicative effects
#' using an inverse-Wishart distribution. The posterior distribution is:
#' \deqn{Suv ~ IW(kappa0 * Suv0 + t(UV) \%*\% UV, n + kappa0)}
#' where UV = cbind(U, V) is the stacked matrix of effects.
#' 
#' This hierarchical prior allows for adaptive shrinkage of the multiplicative
#' effects, with the amount of shrinkage determined by the data through the
#' posterior update.
#' @author Peter Hoff, Shahryar Minhas
#' @export rSuv_fc
rSuv_fc <- function(U, V, Suv0=NULL, kappa0=NULL) 
{
  # Get dimensions
  n <- nrow(U)
  R <- ncol(U)
  
  # Set defaults if not provided
  if(is.null(Suv0)){ Suv0 <- diag(2*R) }
  if(is.null(kappa0)){ kappa0 <- 2 + 2*R }
  
  # Stack U and V
  UV <- cbind(U, V)
  
  # Update covariance using inverse-Wishart
  Suv <- solve(rwish(solve(kappa0*Suv0 + t(UV)%*%UV), n + kappa0))
  
  return(Suv)
}