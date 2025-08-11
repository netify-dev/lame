#' Summary of an AME object
#' 
#' Provides a comprehensive summary of a fitted AME (Additive and Multiplicative 
#' Effects) model, including parameter estimates, standard errors, confidence 
#' intervals, and model diagnostics.
#' 
#' @details
#' The summary includes:
#' \describe{
#'   \item{Regression coefficients}{Point estimates, standard errors, z-values,
#'         p-values, and 95% confidence intervals for dyadic, sender, and 
#'         receiver covariates}
#'   \item{Variance components}{Estimates and standard errors for:
#'     \describe{
#'       \item{va}{Variance of additive sender/row effects (asymmetric networks)}
#'       \item{cab}{Covariance between sender and receiver effects}
#'       \item{vb}{Variance of additive receiver/column effects (asymmetric networks)}
#'       \item{rho}{Dyadic correlation (reciprocity in directed networks)}
#'       \item{ve}{Residual variance}
#'     }
#'     For symmetric networks, only va and ve are estimated.}
#'   \item{Model fit statistics}{If available, AIC and BIC for model comparison}
#' }
#' 
#' @param object an object of class "ame", typically the result of fitting an 
#'        AME model using the \code{ame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.ame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates}
#'   \item{model.fit}{Model fit statistics (AIC, BIC) if available}
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
#' @seealso \code{\link{ame}}, \code{\link{print.summary.ame}}
#' @method summary ame
#' @export
summary.ame <- function(object, ...) {
  fit <- object
  
  # Coefficient statistics
  beta_mean <- apply(fit$BETA, 2, mean)
  beta_sd <- apply(fit$BETA, 2, sd)
  beta_z <- beta_mean / beta_sd
  beta_p <- 2 * (1 - pnorm(abs(beta_z)))
  
  # CI bounds
  z_val <- 1.96 
  beta_lower <- beta_mean - z_val * beta_sd
  beta_upper <- beta_mean + z_val * beta_sd
  
  # Table
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
  
  # Model fit statistics if available
  model_fit <- NULL
  if (!is.null(fit$AIC) || !is.null(fit$BIC)) {
    model_fit <- c(AIC = fit$AIC, BIC = fit$BIC)
  }
  
  # Create summary object
  sum_obj <- list(
    call = fit$call,
    beta = beta_table,
    variance = vc_table,
    model.fit = model_fit
  )
  
  class(sum_obj) <- "summary.ame"
  return(sum_obj)
}

#' Print method for summary.ame objects
#' 
#' Prints a formatted summary of an AME model fit
#' 
#' @param x a summary.ame object
#' @param digits number of digits to display (default: 3)
#' @param ... additional arguments (not used)
#' @return the summary.ame object invisibly
#' @export
print.summary.ame <- function(x, digits = 3, ...) {
  # Print model call
  cat("\n=== AME Model Summary ===\n")
  cat("\nCall:\n")
  print(x$call)
  
  # Print regression coefficients
  cat("\nRegression coefficients:\n")
  cat("------------------------\n")
  coef_table <- x$beta
  coef_table <- round(coef_table, digits)
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
  cat("-------------------\n")
  var_table <- round(x$variance, digits)
  print(var_table, quote = FALSE, right = TRUE)
  
  # Print model fit statistics if available
  if (!is.null(x$model.fit)) {
    cat("\nModel fit statistics:\n")
    cat("--------------------\n")
    print(round(x$model.fit, digits))
  }
  
  invisible(x)
}