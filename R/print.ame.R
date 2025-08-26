#' Print method for AME model objects
#' 
#' @description
#' Displays a formatted summary of a fitted AME (Additive and Multiplicative Effects) model.
#' This method provides a concise overview of model structure, parameter estimates,
#' and goodness-of-fit statistics without generating new data.
#' 
#' @param x an object of class "ame" from fitting an AME model
#' @param ... additional arguments (not currently used)
#' 
#' @return Invisibly returns the input object (for method chaining)
#' 
#' @details
#' The print method displays:
#' \itemize{
#'   \item Model type (unipartite/bipartite, symmetric/asymmetric)
#'   \item Network dimensions and observation count
#'   \item Family and link function used
#'   \item Number of MCMC iterations
#'   \item Parameter counts for regression coefficients and latent factors
#'   \item Basic convergence diagnostics if available
#' }
#' 
#' Unlike \code{simulate}, this method only formats existing results
#' for display and does not perform any new computations or data generation.
#' 
#' @examples
#' \dontrun{
#' # Fit model
#' fit <- ame(Y, X, R = 2)
#' 
#' # Display summary
#' print(fit)
#' # or simply:
#' fit
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @seealso 
#' \code{\link{summary.ame}} for detailed summaries,
#' \code{\link{simulate.ame}} for generating new networks,
#' \code{\link{predict.ame}} for predictions
#' @method print ame
#' @export
print.ame <- function(x, ...) {
  cat("\nAdditive and Multiplicative Effects (AME) Model\n")
  cat("================================================\n")
  
  # Basic model info
  n <- nrow(x$APM)
  nscan <- nrow(x$BETA)
  
  cat("\nNetwork dimensions: ", n, "x", n, "\n")
  cat("MCMC iterations: ", nscan, "\n")
  
  # Model type info
  if (!is.null(x$model.name)) {
    cat("Model type: ", x$model.name, "\n")
  }
  
  # Number of parameters
  cat("\nNumber of parameters:\n")
  cat("  Regression coefficients: ", ncol(x$BETA), "\n")
  if (!is.null(x$U)) {
    cat("  Multiplicative effects dimension: ", ncol(x$U), "\n")
  }
  
  
  cat("\nUse summary(object) for detailed results\n")
  
  invisible(x)
}