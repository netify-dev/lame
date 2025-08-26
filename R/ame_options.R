#' AME model fitting options
#' 
#' @description
#' Configure options for fitting AME models. Memory efficiency is now
#' handled automatically based on network size.
#' 
#' @param parallel_chains Number of parallel chains to run (default: 1)
#' @param verbose Logical; print progress information (default: TRUE)
#' @param odens Output density - save every odens iterations (default: 25)
#' @param use_sparse_matrices Logical; use sparse matrices for storing results (default: FALSE).
#'   Set to TRUE if your network is actually sparse (many zero/NA entries) and memory is a concern.
#' 
#' @return List of options to pass to ame()
#' 
#' @details
#' Memory optimization features:
#' \itemize{
#'   \item Redundant matrices (EZ, UVPM) are never stored - they can be reconstructed if needed
#'   \item use_sparse_matrices = TRUE: Converts large matrices to sparse format
#'   \item Posterior samples are always thinned appropriately
#' }
#' 
#' When to use sparse matrices:
#' \itemize{
#'   \item Your network has < 10\% non-zero entries (truly sparse)
#'   \item Memory usage is a critical concern
#'   \item You're willing to trade computational speed for memory efficiency
#' }
#' 
#' Note: For dense networks (most edges observed), sparse matrices will be 
#' slower and may use MORE memory than dense storage.
#' 
#' @examples
#' \dontrun{
#' # Standard usage (memory optimization is automatic)
#' fit <- ame(Y, X)
#' 
#' # Run parallel chains for better convergence
#' fit <- ame(Y, X, parallel_chains = 4)
#' 
#' # Reduce output frequency for very long chains
#' fit <- ame(Y, X, nscan = 100000, odens = 100)
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
ame_options <- function(
  parallel_chains = 1,
  verbose = TRUE,
  odens = 25,
  use_sparse_matrices = FALSE
  ) {
  
  list(
    parallel_chains = parallel_chains,
    verbose = verbose,
    odens = odens,
    use_sparse_matrices = use_sparse_matrices
  )
}
