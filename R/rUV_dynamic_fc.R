#' Gibbs sampling of dynamic U and V with AR(1) evolution
#'
#' Updates latent factor positions U and V that evolve over time according to
#' an AR(1) process: u_\{i,t\} = rho * u_\{i,t-1\} + epsilon_\{i,t\}
#'
#' @usage rUV_dynamic_fc(U, V, ET, rho_uv, sigma_uv, s2, shrink=TRUE, symmetric=FALSE)
#' @param U 3D array of current U positions (n x R x T)
#' @param V 3D array of current V positions (n x R x T)
#' @param ET 3D array of residuals (n x n x T)
#' @param rho_uv AR(1) autoregressive parameter for latent positions
#' @param sigma_uv Innovation standard deviation for latent positions
#' @param s2 dyadic variance
#' @param shrink whether to apply shrinkage (default TRUE)
#' @param symmetric whether the network is symmetric (default FALSE)
#' @return list with updated U and V arrays
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rUV_dynamic_fc
rUV_dynamic_fc <- function(U, V, ET, rho_uv, sigma_uv, s2, shrink=TRUE, symmetric=FALSE) {
  # Use C++ implementation
  result <- rUV_dynamic_fc_cpp(U, V, ET, rho_uv, sigma_uv, s2, shrink, symmetric)
  return(result)
}