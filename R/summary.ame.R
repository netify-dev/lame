#' Summary of an AME object
#' 
#' Summary method for an AME object
#' 
#' @param object the result of fitting an AME model
#' @param ... additional parameters (not used)
#' @return a summary of parameter estimates and confidence intervals for an AME
#' fit
#' @author Peter Hoff, Shahryar Minhas, Tosin Salau
#' @method summary ame
#' @export
summary.ame <- function(object, ...) {
  fit <- object
  require(amen)
  
  #coefficient statistics
  beta_mean <- apply(fit$BETA, 2, mean)
  beta_sd <- apply(fit$BETA, 2, sd)
  beta_z <- beta_mean / beta_sd
  beta_p <- 2 * (1 - pnorm(abs(beta_z)))
  
  #CI bounds
  z_val <- 1.96 
  beta_lower <- beta_mean - z_val * beta_sd
  beta_upper <- beta_mean + z_val * beta_sd
  
  #table
  beta_table <- cbind(
    Estimate = beta_mean,
    StdError = beta_sd,
    z_value = beta_z,
    p_value = beta_p,
    CI_lower = beta_lower,
    CI_upper = beta_upper
  )
  
  # Variance components
  vc_mean <- apply(fit$VC, 2, mean)
  vc_sd <- apply(fit$VC, 2, sd)
  
  vc_table <- cbind(
    Estimate = vc_mean,
    StdError = vc_sd
  )
  
  
  # Create a summary object
  sum_obj <- list(
    call = fit$call,
    beta = beta_table,
    variance = vc_table
  )
  
  class(sum_obj) <- "summary.ame"
  return(sum_obj)
}

#' Print method for summary.ame objects
#' 
#' @param x a summary.ame object
#' @param ... additional arguments (not used)
#' @return the summary.ame object invisibly
#' @export
print.summary.ame <- function(x, ...) {
  # Print model call
  cat("\nCall:\n")
  print(x$call)
  
  # Print regression coefficients
  cat("\nRegression coefficients:\n")
  coef_table <- x$beta
  coef_table <- round(coef_table, 3)
  # Add significance stars
  stars <- symnum(x$beta[, "p_value"], 
                  corr = FALSE,
                  na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
  coef_table <- cbind(coef_table, " " = stars)
  print(coef_table, quote = FALSE, right = TRUE)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  # Print variance components
  cat("\nVariance components:\n")
  var_table <- round(x$variance, 3)
  print(var_table, quote = FALSE, right = TRUE)
  
  invisible(x)
}