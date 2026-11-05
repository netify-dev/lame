#' Simulate a Poisson relational matrix
#'
#' Simulates a relational matrix from a Poisson distribution given
#' the log mean values
#'
#' @usage simY_pois(EZ, bipartite = FALSE)
#' @param EZ matrix giving the log expected values of the relational matrix
#' @param bipartite whether row and column actors form disjoint sets; keeps
#'   diagonal cells even when the matrix is square.
#' @return a matrix with the dimensions of \code{EZ} of counts
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats rpois
#' @keywords internal
#' @export simY_pois
simY_pois <- function(EZ, bipartite = FALSE) {
	suppressWarnings(.Call(`_lame_simY_pois`, EZ, isTRUE(bipartite)))
}
