#' Gibbs update for latent variable in a Poisson AME model
#' 
#' Updates the latent variable Z in a Poisson AME model using a
#' Metropolis-Hastings step. The model assumes y_\{i,j\} ~ Poisson(exp(z_\{i,j\}))
#' where z_\{i,j\} is the latent variable representing the log mean.
#' 
#' @usage rZ_pois_fc(Z, EZ, rho, s2, Y)
#' @param Z a square matrix, the current value of the latent variable
#' @param EZ expected value of Z (regression effects + random effects)
#' @param rho dyadic correlation
#' @param s2 dyadic variance (overdispersion parameter)
#' @param Y square relational matrix of observed counts
#' @return updated value of Z
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats dpois
#' @export rZ_pois_fc
rZ_pois_fc <- function(Z, EZ, rho, s2, Y) {
  # Use C++ version only
  suppressWarnings(rZ_pois_fc_cpp(Z, EZ, rho, s2, Y))
}