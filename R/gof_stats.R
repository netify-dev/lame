#' Goodness of fit statistics
#' 
#' Calculates goodness of fit statistics for relational data matrices,
#' evaluating second-order (dyadic) and third-order (triadic) dependence
#' patterns. Can handle both unipartite and bipartite networks.
#' 
#' @usage gof_stats(Y, mode = NULL, custom_gof = NULL)
#' @param Y a relational data matrix. For unipartite networks, a square n x n 
#'   matrix where Y\\[i,j\\] represents the relationship from node i to node j.
#'   For bipartite networks, an nA x nB matrix where Y\\[i,j\\] 
#'   represents the relationship from node i in set A to node j in set B.
#'   Missing values (NA) are allowed and will be handled appropriately.
#' @param mode character string specifying the network type: "unipartite" or "bipartite".
#'   If NULL (default), attempts to infer from matrix dimensions (rectangular = bipartite,
#'   square = unipartite). Note: square bipartite networks must specify mode="bipartite".
#' @param custom_gof optional function or list of functions to compute custom GOF statistics.
#'   Each function should take Y as input and return a named numeric value or vector.
#'   Custom statistics will be added to the standard statistics.
#' @return A named numeric vector containing goodness-of-fit statistics.
#'   For unipartite networks:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means. Measures the 
#'       heterogeneity in out-degree centrality (sender effects).}
#'     \item{sd.colmean}{Standard deviation of column means. Measures the
#'       heterogeneity in in-degree centrality (receiver effects).}
#'     \item{dyad.dep}{Dyadic dependence/reciprocity correlation.}
#'     \item{cycle.dep}{Cyclic/transitive triadic dependence.}
#'     \item{trans.dep}{Transitive triadic dependence.}
#'   }
#'   For bipartite networks:
#'   \describe{
#'     \item{sd.rowmean}{Standard deviation of row means (sender heterogeneity).}
#'     \item{sd.colmean}{Standard deviation of column means (receiver heterogeneity).}
#'     \item{four.cycles}{Count of four-cycles in the bipartite network.}
#'   }
#'   
#'   If custom_gof is provided, additional statistics will be included with their
#'   user-specified names.
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
#' For bipartite networks with square dimensions (nA = nB), you must explicitly 
#' specify mode="bipartite" to ensure correct statistics are calculated.
#' 
#' Missing values in Y are handled by pairwise deletion for correlations and
#' are excluded from matrix products in triadic calculations.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' 
#' data(YX_nrm) 
#' 
#' # Auto-detect unipartite
#' gof_stats(YX_nrm$Y) 
#' 
#' # Explicitly specify mode for square bipartite
#' # Y_bip <- matrix(rnorm(100), 10, 10) # 10x10 bipartite
#' # gof_stats(Y_bip, mode = "bipartite")
#' 
#' # Custom GOF function
#' # my_stat <- function(Y) { c(my_measure = sum(Y > 0, na.rm=TRUE)) }
#' # gof_stats(YX_nrm$Y, custom_gof = my_stat)
#' 
#' @export gof_stats
gof_stats <- function(Y, mode = NULL, custom_gof = NULL) {
  
  # Determine network mode
  if(is.null(mode)) {
    # Auto-detect based on dimensions
    if(nrow(Y) != ncol(Y)) {
      mode <- "bipartite"
    } else {
      mode <- "unipartite"
      # Only show message once per session using an option
      if(is.null(getOption("lame.gof_stats.msg_shown"))) {
        message("Note: Square matrix assumed to be unipartite. ",
                "Use mode='bipartite' for square bipartite networks.")
        options(lame.gof_stats.msg_shown = TRUE)
      }
    }
  } else {
    mode <- match.arg(mode, c("unipartite", "bipartite"))
  }
  
  # Calculate standard GOF statistics based on mode
  if(mode == "bipartite") {
    standard_gof <- gof_stats_bipartite(Y)
  } else {
    standard_gof <- gof_stats_unipartite(Y)
  }
  
  # Add custom GOF statistics if provided
  if(!is.null(custom_gof)) {
    custom_stats <- NULL
    
    if(is.function(custom_gof)) {
      # Single custom function
      custom_stats <- custom_gof(Y)
    } else if(is.list(custom_gof)) {
      # List of custom functions
      custom_stats <- numeric()
      for(i in seq_along(custom_gof)) {
        if(is.function(custom_gof[[i]])) {
          stat_result <- custom_gof[[i]](Y)
          # Get function name or use index
          stat_name <- names(custom_gof)[i]
          if(is.null(stat_name) || stat_name == "") {
            stat_name <- paste0("custom", i)
          }
          names(stat_result) <- stat_name
          custom_stats <- c(custom_stats, stat_result)
        }
      }
    }
    
    # Combine standard and custom statistics
    if(!is.null(custom_stats)) {
      gof_result <- c(standard_gof, custom_stats)
    } else {
      gof_result <- standard_gof
    }
  } else {
    gof_result <- standard_gof
  }
  
  return(gof_result)
}