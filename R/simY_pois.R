#' Simulate a Poisson relational matrix
#' 
#' Simulates a relational matrix from a Poisson distribution given
#' the log mean values
#' 
#' @usage simY_pois(EZ)
#' @param EZ square matrix giving the log expected values of the relational matrix
#' @return a square matrix of counts
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats rpois
#' @export simY_pois
simY_pois <- function(EZ) {
  # Use C++ version only
  suppressWarnings(.Call(`_lame_simY_pois`, EZ))
}