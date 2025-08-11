#' Goodness of fit statistics
#' 
#' Calculates goodness of fit statistics for relational data matrices,
#' evaluating second-order (dyadic) and third-order (triadic) dependence
#' patterns. These statistics are useful for assessing model fit in
#' network analysis and relational data modeling.
#' 
#' @usage gofstats(Y)
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
#' @author Peter Hoff
#' @examples
#' 
#' data(YX_nrm) 
#' 
#' gofstats(YX_nrm$Y) 
#' 
#' 
#' @export gofstats
gofstats<-function(Y)
{
  # Calculate row and column standard deviations
  sd.rowmean <- sd(rowMeans(Y, na.rm=TRUE), na.rm=TRUE) 
  sd.colmean <- sd(colMeans(Y, na.rm=TRUE), na.rm=TRUE)
  
  # Calculate dyadic dependence (reciprocity)
  dyad.dep <- suppressWarnings(cor(c(Y), c(t(Y)), use="complete.obs")) 
  
  # Center the matrix and create data availability matrix
  Y.mean <- mean(Y, na.rm=TRUE)
  E <- Y - Y.mean
  D <- 1 * (!is.na(E))
  E[is.na(E)] <- 0
  
  # Pre-compute common values to avoid repetition
  Y.sd <- sd(c(Y), na.rm=TRUE)
  Y.sd3 <- Y.sd^3
  
  # Optimize triadic calculations by computing matrix products once
  # For cycle dependence: E*E*E (all same matrix)
  EE <- E %*% E  # Compute once, reuse
  EEE <- EE %*% E
  DD <- D %*% D  # Compute once, reuse  
  DDD <- DD %*% D
  
  # For transitive dependence: E*t(E)*E
  Et <- t(E)  # Transpose once
  Dt <- t(D)  # Transpose once
  EtE <- E %*% Et %*% E
  DtD <- D %*% Dt %*% D
  
  # Calculate triadic dependencies with pre-computed values
  cycle.dep <- sum(diag(EEE)) / (sum(diag(DDD)) * Y.sd3)
  trans.dep <- sum(diag(EtE)) / (sum(diag(DtD)) * Y.sd3)
  
  # Compile results
  gof <- c(sd.rowmean, sd.colmean, dyad.dep, cycle.dep, trans.dep)
  
  # Handle any NaN/NA values (e.g., from empty networks)
  gof[is.na(gof)] <- 0 
  
  names(gof) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
  gof
}
