#' Goodness of fit statistics for bipartite networks
#' 
#' Calculates goodness of fit statistics specifically designed for bipartite 
#' networks, evaluating degree heterogeneity and higher-order dependencies.
#' 
#' @param Y a bipartite relational data matrix (nA x nB rectangular matrix) 
#'   where Y\\[i,j\\] represents the relationship from node i in set A to 
#'   node j in set B. Missing values (NA) are allowed and will be handled 
#'   appropriately.
#' @return A named numeric vector containing bipartite-specific goodness-of-fit 
#'   statistics:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means. Measures the 
#'       heterogeneity in out-degree from set A nodes (sender effects). 
#'       Higher values indicate more variation in how active A nodes are.}
#'     \item{sd.colmean}{Standard deviation of column means. Measures the
#'       heterogeneity in in-degree to set B nodes (receiver effects). 
#'       Higher values indicate more variation in how popular B nodes are.}
#'     \item{four.cycles}{Count of four-cycles (also called 4-paths or squares) in 
#'       the bipartite network. A four-cycle occurs when two nodes from set A 
#'       (e.g., i and k) both connect to the same two nodes in set B (e.g., j and l),
#'       forming a closed path: i->j->k->l->i. This measures the tendency for pairs 
#'       of A-nodes to share multiple common B-node connections, capturing a form of 
#'       clustering specific to bipartite networks. High four-cycle counts indicate 
#'       that connections are not random but show patterns of shared preferences or 
#'       co-occurrence. For example, in a user-item network, many four-cycles suggest 
#'       that users who like one item tend to also like other items that co-occur 
#'       with it.}
#'   }
#' @details
#' For bipartite networks, reciprocity and triadic closure are not meaningful
#' concepts since edges only exist between the two node sets. Instead, this
#' function focuses on:
#' \itemize{
#'   \item Degree heterogeneity in both node sets
#'   \item Four-cycles as the simplest higher-order dependence pattern
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
#' @examples
#' \dontrun{
#' # Create a random bipartite network
#' Y <- matrix(rnorm(10*12), 10, 12)
#' 
#' # Calculate GOF statistics
#' gof_stats_bipartite(Y)
#' }
#' @export
gof_stats_bipartite <- function(Y) {
  
  # Check that matrix is rectangular (bipartite)
  if(nrow(Y) == ncol(Y)) {
    warning("Y appears to be square. Use gof_stats() for unipartite networks.")
  }
  
  # Core bipartite statistics
  gof <- numeric(3)
  
  # 1. Sender heterogeneity (std dev of row means)
  gof[1] <- sd(rowMeans(Y, na.rm = TRUE), na.rm = TRUE)
  
  # 2. Receiver heterogeneity (std dev of column means)  
  gof[2] <- sd(colMeans(Y, na.rm = TRUE), na.rm = TRUE)
  
  # 3. Four-cycles (requires C++ function)
  # Convert to 3D array for C++ function (nA x nB x 1)
  Y_arr <- array(Y, dim = c(nrow(Y), ncol(Y), 1))
  gof[3] <- as.numeric(lame:::count_four_cycles_bip_cpp(Y_arr))
  
  names(gof) <- c("sd.rowmean", "sd.colmean", "four.cycles")
  
  return(gof)
}