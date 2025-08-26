# Investigation of coverage issues in lame package
# Focus on why additive effects model has 80% coverage instead of 95%

library(lame)
library(testthat)

# Run a more focused simulation to diagnose the coverage issue
diagnose_coverage <- function(n_sims = 50, n = 40, burn = 400, nscan = 1500) {
  set.seed(6886)
  
  mu_true <- 0
  beta_true <- 1
  
  results <- list()
  
  for(i in 1:n_sims) {
    set.seed(i)
    
    # Generate data with additive effects
    xw <- matrix(rnorm(n*2), n, 2)
    X <- tcrossprod(xw[,1])
    diag(X) <- NA
    
    # True additive effects
    a_true <- rnorm(n, 0, 0.5)
    b_true <- rnorm(n, 0, 0.5)
    
    eta <- mu_true + beta_true * X + 
           outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
    
    sigma2_true <- 1
    Y <- eta + matrix(rnorm(n*n, 0, sqrt(sigma2_true)), n, n)
    diag(Y) <- NA
    
    # Fit model
    fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
              rvar = TRUE, cvar = TRUE, dcor = FALSE,
              burn = burn, nscan = nscan,
              print = FALSE)
    
    # Extract estimates
    beta_samples <- fit$BETA[,2]
    beta_hat <- median(beta_samples)
    beta_se <- sd(beta_samples)
    beta_ci <- quantile(beta_samples, c(0.025, 0.975))
    
    # Effective sample size
    ess <- if(requireNamespace("coda", quietly = TRUE)) {
      coda::effectiveSize(beta_samples)
    } else {
      NA
    }
    
    results[[i]] <- list(
      beta_hat = beta_hat,
      beta_se = beta_se,
      beta_ci_lower = beta_ci[1],
      beta_ci_upper = beta_ci[2],
      beta_covered = (beta_true >= beta_ci[1]) & (beta_true <= beta_ci[2]),
      ess = ess,
      ci_width = beta_ci[2] - beta_ci[1]
    )
  }
  
  return(results)
}

test_that("Diagnose coverage for additive effects model", {
  results <- diagnose_coverage(n_sims = 30, n = 40, burn = 400, nscan = 1500)
  
  # Calculate summary statistics
  coverage <- mean(sapply(results, function(x) x$beta_covered))
  mean_estimate <- mean(sapply(results, function(x) x$beta_hat))
  mean_se <- mean(sapply(results, function(x) x$beta_se))
  mean_ci_width <- mean(sapply(results, function(x) x$ci_width))
  
  # Check which simulations failed to cover
  not_covered <- which(!sapply(results, function(x) x$beta_covered))
  
  cat("\n=== Coverage Diagnosis ===\n")
  cat("Coverage rate:", coverage, "\n")
  cat("Mean beta estimate:", mean_estimate, "(true = 1)\n")
  cat("Mean SE:", mean_se, "\n")
  cat("Mean CI width:", mean_ci_width, "\n")
  
  if(length(not_covered) > 0) {
    cat("\nSimulations that failed to cover true value:\n")
    for(i in not_covered[1:min(5, length(not_covered))]) {
      cat(sprintf("Sim %d: beta_hat=%.3f, CI=[%.3f, %.3f]\n",
                 i, results[[i]]$beta_hat, 
                 results[[i]]$beta_ci_lower,
                 results[[i]]$beta_ci_upper))
    }
  }
  
  # Check if CIs are systematically too narrow
  theoretical_se <- mean_se
  empirical_se <- sd(sapply(results, function(x) x$beta_hat))
  
  cat("\nSE comparison:\n")
  cat("Empirical SE of estimates:", empirical_se, "\n")
  cat("Mean estimated SE:", theoretical_se, "\n")
  cat("Ratio (should be ~1):", empirical_se/theoretical_se, "\n")
  
  # If ratio > 1, the model is underestimating uncertainty
  if(empirical_se/theoretical_se > 1.1) {
    warning("Model appears to be underestimating uncertainty")
  }
})

# Test different prior specifications
test_that("Test effect of prior specifications", {
  set.seed(6886)
  n <- 40
  
  # Generate one dataset
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  eta <- 0 + 1 * X + 
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  
  Y <- eta + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Fit with default priors
  fit_default <- ame(Y, Xdyad = X, R = 0, family = "normal",
                    rvar = TRUE, cvar = TRUE, dcor = FALSE,
                    burn = 500, nscan = 2000,
                    print = FALSE)
  
  # Fit with specified g-prior
  prior_g <- list(g = n * var(c(Y), na.rm = TRUE))
  fit_gprior <- ame(Y, Xdyad = X, R = 0, family = "normal",
                   rvar = TRUE, cvar = TRUE, dcor = FALSE,
                   burn = 500, nscan = 2000,
                   print = FALSE,
                   prior = prior_g)
  
  # Fit with stronger prior on variance components
  prior_strong <- list(
    Sab0 = diag(2) * 0.25,  # Stronger prior on additive effects
    eta0 = 10  # More degrees of freedom
  )
  fit_strong <- ame(Y, Xdyad = X, R = 0, family = "normal",
                   rvar = TRUE, cvar = TRUE, dcor = FALSE,
                   burn = 500, nscan = 2000,
                   print = FALSE,
                   prior = prior_strong)
  
  cat("\n=== Effect of Priors ===\n")
  cat("Default prior:\n")
  cat("  Beta estimate:", median(fit_default$BETA[,2]), "\n")
  cat("  Beta SE:", sd(fit_default$BETA[,2]), "\n")
  
  cat("G-prior specified:\n")
  cat("  Beta estimate:", median(fit_gprior$BETA[,2]), "\n")
  cat("  Beta SE:", sd(fit_gprior$BETA[,2]), "\n")
  
  cat("Strong variance prior:\n")
  cat("  Beta estimate:", median(fit_strong$BETA[,2]), "\n")
  cat("  Beta SE:", sd(fit_strong$BETA[,2]), "\n")
})

# Check MCMC convergence
test_that("Check MCMC convergence", {
  set.seed(6886)
  n <- 30
  
  # Simple data
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  Y <- 1 * X + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Run with longer chain
  fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
            rvar = TRUE, cvar = TRUE, dcor = FALSE,
            burn = 1000, nscan = 5000,
            print = FALSE)
  
  if(requireNamespace("coda", quietly = TRUE)) {
    beta_chain <- coda::mcmc(fit$BETA[,2])
    ess <- coda::effectiveSize(beta_chain)
    
    cat("\n=== MCMC Diagnostics ===\n")
    cat("Effective sample size for beta:", ess, "\n")
    cat("Efficiency:", ess/length(fit$BETA[,2]), "\n")
    
    # Geweke diagnostic
    geweke <- coda::geweke.diag(beta_chain)
    cat("Geweke z-score:", geweke$z, "\n")
    cat("Geweke p-value:", 2*pnorm(-abs(geweke$z)), "\n")
    
    if(abs(geweke$z) > 2) {
      warning("Potential convergence issue detected")
    }
  }
})