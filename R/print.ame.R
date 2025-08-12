#' Print method for AME objects
#' 
#' Provides a concise print output for fitted AME models
#' 
#' @param x an object of class "ame"
#' @param ... additional arguments (not used)
#' @return the ame object invisibly
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
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