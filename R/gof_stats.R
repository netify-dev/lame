#' Goodness of fit statistics
#' 
#' Calculates goodness of fit statistics for relational data matrices,
#' evaluating second-order (dyadic) and third-order (triadic) dependence
#' patterns. These statistics are useful for assessing model fit in
#' network analysis and relational data modeling.
#' 
#' @usage gof_stats(Y)
#' @param Y a relational data matrix (n x n square matrix) where Y\\[i,j\\]
#'   represents the relationship from node i to node j. Missing values (NA)
#'   are allowed and will be handled appropriately.
#' @return A named numeric vector containing five goodness-of-fit statistics:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means. Measures the 
#'       heterogeneity in out-degree centrality (sender effects). Higher values
#'       indicate more variation in how active nodes are as senders.}
#'     \item{sd.colmean}{Standard deviation of column means. Measures the
#'       heterogeneity in in-degree centrality (receiver effects). Higher values
#'       indicate more variation in how popular nodes are as receivers.}
#'     \item{dyad.dep}{Dyadic dependence/reciprocity correlation. Pearson correlation
#'       between Y\\[i,j\\] and Y\\[j,i\\] across all dyads. Positive values indicate
#'       reciprocity (mutual relationships), negative values indicate
#'       anti-reciprocity. Range: \\[-1, 1\\].}
#'     \item{cycle.dep}{Cyclic/transitive triadic dependence. Normalized sum of
#'       products along three-cycles (i to j to k to i). Positive values indicate
#'       transitivity clustering, where 'a friend of a friend is a friend'.
#'       Based on the trace of the cubed centered matrix, normalized by the 
#'       trace of the cubed data availability matrix and the cubed standard deviation.}
#'     \item{trans.dep}{Transitive triadic dependence. Normalized sum of products
#'       along two-paths that close into triangles (i to j to k with k to i). Measures
#'       the tendency for open triads to close. Based on the trace of the product
#'       E*E'*E where E is the centered matrix, normalized appropriately.}
#'   }
#' @details
#' The function computes network statistics that capture different aspects of
#' network structure beyond simple density. These statistics are particularly
#' useful for:
#' \itemize{
#'   \item Model checking: comparing observed statistics to those from simulated networks
#'   \item Model selection: choosing between models that better capture network dependencies
#'   \item Descriptive analysis: summarizing key structural features of the network
#' }
#' 
#' Missing values in Y are handled by pairwise deletion for correlations and
#' are excluded from matrix products in triadic calculations.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' 
#' data(YX_nrm) 
#' 
#' gof_stats(YX_nrm$Y) 
#' 
#' 
#' @export gof_stats
gof_stats<-function(Y)
{
  # Check if matrix is square - GOF stats require square matrices
  if(nrow(Y) != ncol(Y)) {
    # For non-square (bipartite) matrices, return limited stats
    gof <- c(
      sd(rowMeans(Y, na.rm=TRUE), na.rm=TRUE),  # sd.rowmean
      sd(colMeans(Y, na.rm=TRUE), na.rm=TRUE),  # sd.colmean
      NA,  # dyad.dep (not applicable for bipartite)
      NA,  # cycle.dep (not applicable for bipartite)
      NA   # trans.dep (not applicable for bipartite)
    )
    names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    return(gof)
  }
  
  tryCatch({
    gof <- as.vector(gof_stats_cpp(as.matrix(Y)))
    names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    gof[is.na(gof)] <- 0
    return(gof)
  }, error = function(e) {
    sd.rowmean <- sd(rowMeans(Y, na.rm=TRUE), na.rm=TRUE) 
    sd.colmean <- sd(colMeans(Y, na.rm=TRUE), na.rm=TRUE)
    
    # Handle case where there are no complete pairs for correlation
    dyad.dep <- tryCatch({
      suppressWarnings(cor(c(Y), c(t(Y)), use="complete.obs"))
    }, error = function(e2) {
      # If no complete pairs, try pairwise.complete.obs
      tryCatch({
        suppressWarnings(cor(c(Y), c(t(Y)), use="pairwise.complete.obs"))
      }, error = function(e3) {
        # If still fails, return 0 (no dyadic dependence detected)
        0
      })
    })
    
    Y.mean <- mean(Y, na.rm=TRUE)
    E <- Y - Y.mean
    D <- 1 * (!is.na(E))
    E[is.na(E)] <- 0
    Y.sd <- sd(c(Y), na.rm=TRUE)
    
    # Handle case where Y.sd is 0 or NA
    if(is.na(Y.sd) || Y.sd == 0) {
      Y.sd3 <- 1  # Avoid division by zero
      cycle.dep <- 0
      trans.dep <- 0
    } else {
      Y.sd3 <- Y.sd^3
      EEE <- E %*% E %*% E
      DDD <- D %*% D %*% D
      EtE <- E %*% t(E) %*% E
      DtD <- D %*% t(D) %*% D
      
      ddd_diag_sum <- sum(diag(DDD))
      dtd_diag_sum <- sum(diag(DtD))
      
      cycle.dep <- if(ddd_diag_sum > 0) sum(diag(EEE)) / (ddd_diag_sum * Y.sd3) else 0
      trans.dep <- if(dtd_diag_sum > 0) sum(diag(EtE)) / (dtd_diag_sum * Y.sd3) else 0
    }
    
    gof <- c(sd.rowmean, sd.colmean, dyad.dep, cycle.dep, trans.dep)
    gof[is.na(gof)] <- 0 
    names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    gof
  })
}
