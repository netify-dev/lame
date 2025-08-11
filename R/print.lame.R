#' Print method for LAME objects
#' 
#' Provides a concise print output for fitted LAME models
#' 
#' @param x an object of class "lame"
#' @param ... additional arguments (not used)
#' @return the lame object invisibly
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
#' @method print lame
#' @export
print.lame <- function(x, ...) {
  cat("\nLongitudinal Additive and Multiplicative Effects (LAME) Model\n")
  cat("==============================================================\n")
  
  # Basic model info
  n_periods <- length(x$Y_T)
  if (n_periods > 0) {
    n_nodes <- nrow(x$Y_T[[1]])
    cat("\nNumber of time periods: ", n_periods, "\n")
    cat("Network dimensions: ", n_nodes, "x", n_nodes, " (per period)\n")
  }
  
  nscan <- nrow(x$BETA)
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
  
  # Composition changes if applicable
  if (!is.null(x$nodeID)) {
    total_nodes <- length(unique(unlist(x$nodeID)))
    cat("  Total unique nodes across periods: ", total_nodes, "\n")
  }
  
  # Model fit if available
  if (!is.null(x$AIC) || !is.null(x$BIC)) {
    cat("\nModel fit:\n")
    if (!is.null(x$AIC)) cat("  AIC: ", round(x$AIC, 2), "\n")
    if (!is.null(x$BIC)) cat("  BIC: ", round(x$BIC, 2), "\n")
  }
  
  cat("\nUse summary(object) for detailed results\n")
  
  invisible(x)
}