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
#' @import cli
print.lame <- function(x, ...) {
  cli::cli_h1("Longitudinal Additive and Multiplicative Effects (LAME) Model")
  
  # Basic model info
  n_periods <- length(x$Y_T)
  if (n_periods > 0) {
    n_nodes <- nrow(x$Y_T[[1]])
    cli::cli_h3("Network Structure")
    cli::cli_bullets(c(
      "*" = "Time periods: {.val {n_periods}}",
      "*" = "Network dimensions: {.val {n_nodes}} x {.val {n_nodes}} per period"
    ))
  }
  
  nscan <- nrow(x$BETA)
  cli::cli_bullets(c(
    "*" = "MCMC iterations: {.val {nscan}}"
  ))
  
  # Model type info
  if (!is.null(x$model.name)) {
    cli::cli_bullets(c(
      "*" = "Model type: {.field {x$model.name}}"
    ))
  }
  
  # Check for dynamic UV
  is_dynamic_uv <- !is.null(x$U) && length(dim(x$U)) == 3
  
  # Number of parameters
  cli::cli_h3("Model Parameters")
  param_info <- character()
  param_info <- c(param_info, "*" = "Regression coefficients: {.val {ncol(x$BETA)}}")
  
  if (!is.null(x$U)) {
    if (is_dynamic_uv) {
      R <- dim(x$U)[2]
      T <- dim(x$U)[3]
      param_info <- c(param_info, 
        "*" = "Multiplicative effects: {.val {R}}-dimensional {.emph (dynamic over {T} periods)}")
    } else {
      param_info <- c(param_info, 
        "*" = "Multiplicative effects: {.val {ncol(x$U)}}-dimensional")
    }
  }
  cli::cli_bullets(param_info)
  
  # Composition changes if applicable
  if (!is.null(x$nodeID)) {
    total_nodes <- length(unique(unlist(x$nodeID)))
    cli::cli_bullets(c(
      "i" = "Total unique nodes across periods: {.val {total_nodes}}"
    ))
  }
  
  
  cli::cli_alert_info("Use {.code summary(object)} for detailed results")
  
  invisible(x)
}