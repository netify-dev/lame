#' Summary of a LAME object
#' 
#' Provides a comprehensive summary of a fitted LAME (Longitudinal Additive and
#' Multiplicative Effects) model, including parameter estimates, standard errors,
#' confidence intervals, and model diagnostics.
#' 
#' @details
#' The summary includes:
#' \describe{
#'   \item{Regression coefficients}{Point estimates, standard errors, z-values,
#'         p-values, and 95% confidence intervals for dyadic, sender, and 
#'         receiver covariates}
#'   \item{Variance components}{Estimates and standard errors for:
#'     \describe{
#'       \item{va}{Variance of additive sender/row effects}
#'       \item{cab}{Covariance between sender and receiver effects}
#'       \item{vb}{Variance of additive receiver/column effects}
#'       \item{rho}{Dyadic correlation (reciprocity)}
#'       \item{ve}{Residual variance}
#'     }}
#' }
#' 
#' @param object an object of class "lame", typically the result of fitting a 
#'        longitudinal AME model using the \code{lame} function
#' @param ... additional parameters (currently not used)
#' @return A list of class "summary.lame" containing:
#'   \item{call}{The original function call}
#'   \item{beta}{Matrix of regression coefficient estimates and statistics}
#'   \item{variance}{Matrix of variance component estimates}
#'   \item{n.periods}{Number of time periods in the longitudinal data}
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
#' @seealso \code{\link{lame}}, \code{\link{print.summary.lame}}
#' @method summary lame
#' @export
summary.lame <- function(object, ...) {
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
  
  # Model fit statistics removed
  
  # Number of time periods if available
  n_periods <- if (!is.null(fit$n.periods)) fit$n.periods else NULL
  
  # Create summary object
  sum_obj <- list(
    call = fit$call,
    beta = beta_table,
    variance = vc_table,
    n.periods = n_periods
  )
  
  class(sum_obj) <- "summary.lame"
  return(sum_obj)
}

#' Print method for summary.lame objects
#' 
#' Prints a formatted summary of a LAME model fit
#' 
#' @param x a summary.lame object
#' @param digits number of digits to display (default: 3)
#' @param ... additional arguments (not used)
#' @return the summary.lame object invisibly
#' @export
print.summary.lame <- function(x, digits = 3, ...) {
  # Print model call
  cat("\n=== Longitudinal AME Model Summary ===\n")
  cat("\nCall:\n")
  print(x$call)
  
  # Print number of periods if available
  if (!is.null(x$n.periods)) {
    cat("\nNumber of time periods:", x$n.periods, "\n")
  }
  
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
  
  
  invisible(x)
}