#' Simulate Z based on a probit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix Y
#' 
#' 
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix , the new value of Z
#' @author Peter Hoff
#' @export rZ_bin_fc
rZ_bin_fc <-
function(Z,EZ,rho,Y)
{ 
  # Use C++ version only
  rZ_bin_fc_single(Z, EZ, rho, Y)
}