#' Goodness of fit statistics for unipartite networks
#' 
#' Calculates goodness of fit statistics specifically for unipartite (square) networks,
#' evaluating second-order (dyadic) and third-order (triadic) dependence patterns.
#' 
#' @usage gof_stats_unipartite(Y)
#' @param Y a square n x n relational data matrix where Y[i,j] represents the 
#'   relationship from node i to node j. Missing values (NA) are allowed and will 
#'   be handled appropriately. Diagonal values are typically NA for non-self-loop networks.
#' @return A named numeric vector containing five goodness-of-fit statistics:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means. Measures the 
#'       heterogeneity in out-degree centrality (sender effects).}
#'     \item{sd.colmean}{Standard deviation of column means. Measures the
#'       heterogeneity in in-degree centrality (receiver effects).}
#'     \item{dyad.dep}{Dyadic dependence/reciprocity correlation. Measures the
#'       correlation between Y[i,j] and Y[j,i], capturing reciprocity patterns.}
#'     \item{cycle.dep}{Cyclic triadic dependence. Measures the tendency for
#'       directed cycles (i->j->k->i) in the network.}
#'     \item{trans.dep}{Transitive triadic dependence. Measures the tendency for
#'       transitivity (if i->j and j->k, then i->k) in the network.}
#'   }
#' @details
#' This function computes network statistics that capture different aspects of
#' network structure beyond simple density. These statistics are particularly
#' useful for evaluating how well a model captures the observed network patterns.
#' 
#' The dyadic dependence statistic captures reciprocity - the tendency for 
#' relationships to be mutual. The triadic statistics capture different forms
#' of triadic closure that are common in social networks.
#' 
#' Missing values in Y are handled by pairwise deletion for correlations and
#' are excluded from matrix products in triadic calculations.
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' 
#' # Create a random unipartite network
#' Y <- matrix(rnorm(100), 10, 10)
#' diag(Y) <- NA
#' gof_stats_unipartite(Y)
#' 
#' @export
gof_stats_unipartite <- function(Y) {
  # Check that Y is square
  if(nrow(Y) != ncol(Y)) {
    stop("Y must be a square matrix for unipartite networks. Use gof_stats_bipartite() for rectangular matrices.")
  }
  
  # Try to use C++ implementation first for speed
  tryCatch({
    gof <- as.vector(gof_stats_cpp(as.matrix(Y)))
    names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    gof[is.na(gof)] <- 0
    return(gof)
  }, error = function(e) {
    # Fallback to R implementation if C++ fails
    
    # Degree heterogeneity statistics
    sd.rowmean <- sd(rowMeans(Y, na.rm=TRUE), na.rm=TRUE) 
    sd.colmean <- sd(colMeans(Y, na.rm=TRUE), na.rm=TRUE)
    
    # Dyadic dependence (reciprocity)
    dyad.dep <- tryCatch({
      suppressWarnings(cor(c(Y), c(t(Y)), use="complete.obs"))
    }, error = function(e2) {
      tryCatch({
        suppressWarnings(cor(c(Y), c(t(Y)), use="pairwise.complete.obs"))
      }, error = function(e3) {
        0  # No dyadic dependence detected
      })
    })
    
    # Prepare for triadic statistics
    Y.mean <- mean(Y, na.rm=TRUE)
    E <- Y - Y.mean  # Centered matrix
    D <- 1 * (!is.na(E))  # Indicator matrix for non-missing values
    E[is.na(E)] <- 0  # Replace NA with 0 for matrix multiplication
    Y.sd <- sd(c(Y), na.rm=TRUE)
    
    # Calculate triadic statistics
    if(is.na(Y.sd) || Y.sd == 0) {
      # No variation in Y, so triadic stats are 0
      cycle.dep <- 0
      trans.dep <- 0
    } else {
      Y.sd3 <- Y.sd^3
      
      # Cyclic triadic dependence: i->j->k->i
      EEE <- E %*% E %*% E
      DDD <- D %*% D %*% D
      ddd_diag_sum <- sum(diag(DDD))
      cycle.dep <- if(ddd_diag_sum > 0) sum(diag(EEE)) / (ddd_diag_sum * Y.sd3) else 0
      
      # Transitive triadic dependence: i->j, j->k, i->k
      EtE <- E %*% t(E) %*% E
      DtD <- D %*% t(D) %*% D
      dtd_diag_sum <- sum(diag(DtD))
      trans.dep <- if(dtd_diag_sum > 0) sum(diag(EtE)) / (dtd_diag_sum * Y.sd3) else 0
    }
    
    # Combine and return
    gof <- c(sd.rowmean, sd.colmean, dyad.dep, cycle.dep, trans.dep)
    gof[is.na(gof)] <- 0 
    names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    return(gof)
  })
}