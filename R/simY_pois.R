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
  tryCatch({
    return(.Call(`_lame_simY_pois`, EZ))
  }, error = function(e) {
    n <- nrow(EZ)
    Y <- matrix(NA, n, n)
    for(i in 1:n) {
      for(j in 1:n) {
        if(i != j) {
          lambda <- exp(EZ[i,j])
          # Cap lambda to avoid numerical issues
          lambda <- min(lambda, 1e6)
          Y[i,j] <- rpois(1, lambda)
        }
      }
    }
    Y
  })
}