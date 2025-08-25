#' AME model fitting for longitudinal bipartite networks
#' 
#' Internal function that handles AME model fitting for longitudinal bipartite networks.
#' This function would contain the core logic for time-series of rectangular adjacency 
#' matrices with distinct row and column node sets.
#' 
#' @inheritParams lame
#' @param R_row integer: dimension of row node multiplicative effects
#' @param R_col integer: dimension of column node multiplicative effects
#' @param dynamic_G logical: fit dynamic interaction matrix G
#' @return A lame object with posterior samples and estimates
#' @keywords internal
#' @noRd
lame_bipartite <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
                          rvar = TRUE, cvar = TRUE, 
                          R_row = 0, R_col = 0,
                          dynamic_uv = FALSE, dynamic_ab = FALSE, dynamic_G = FALSE,
                          family="normal", intercept=TRUE,
                          odmax=NULL, prior=list(), g=NA,
                          seed = 6886, nscan = 10000, burn = 500, odens = 25,
                          plot=FALSE, print = FALSE, gof=TRUE, 
                          start_vals=NULL, periodic_save=FALSE, out_file=NULL,
                          save_interval=0.25, model.name=NULL) {
  
  # Bipartite longitudinal networks are not yet fully implemented
  cli::cli_abort(c(
    "Longitudinal bipartite network models are not yet fully implemented.",
    "i" = "Please use unipartite networks for longitudinal analysis.",
    "i" = "Full bipartite support is coming in a future release."
  ))
  
  # TODO: Implement bipartite-specific logic here
  # Key differences from unipartite:
  # - No diagonal to set to NA
  # - Separate row (nA) and column (nB) dimensions across time
  # - No dyadic correlation (rho)
  # - No symmetric option
  # - Separate U (nA x R_row x T) and V (nB x R_col x T) arrays
  # - Interaction matrix G (R_row x R_col x T) for UV interaction
  # - Dynamic G option for time-varying interaction matrix
  # - Different GOF statistics (no transitivity/cycles)
}