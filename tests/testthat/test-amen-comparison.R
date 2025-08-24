# Comparison test between original amen package and lame package
# Testing parameter recovery and confidence interval calibration

library(testthat)

# Try to install amen if not already installed
if (!requireNamespace("amen", quietly = TRUE)) {
  tryCatch({
    install.packages("amen", repos = "https://cloud.r-project.org")
  }, error = function(e) {
    message("Could not install amen package: ", e$message)
  })
}

# Load packages
library(lame)
if (requireNamespace("amen", quietly = TRUE)) {
  library(amen)
}

# Helper function to run simulation with both packages
compare_amen_lame <- function(seed, n = 40, mu = 0, beta = 1, 
                              rvar = FALSE, cvar = FALSE, R = 0,
                              burn = 500, nscan = 2000) {
  set.seed(seed)
  
  # Generate data
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Build linear predictor
  eta <- mu + beta * X
  
  # Add additive effects if requested
  if(rvar || cvar) {
    a_true <- if(rvar) rnorm(n, 0, 0.5) else rep(0, n)
    b_true <- if(cvar) rnorm(n, 0, 0.5) else rep(0, n)
    eta <- eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  }
  
  # Add multiplicative effects if R > 0
  if(R > 0) {
    U_true <- matrix(rnorm(n*R, 0, 1), n, R)
    V_true <- matrix(rnorm(n*R, 0, 1), n, R)
    eta <- eta + tcrossprod(U_true, V_true)
  }
  
  # Generate Gaussian outcome with noise
  sigma2 <- 1
  Y <- eta + matrix(rnorm(n*n, 0, sqrt(sigma2)), n, n)
  diag(Y) <- NA
  
  results <- list()
  
  # Fit with amen package if available
  if (requireNamespace("amen", quietly = TRUE)) {
    fit_amen <- amen::ame(Y, Xdyad = X, R = R, family = "nrm",
                          rvar = rvar, cvar = cvar, dcor = FALSE,
                          burn = burn, nscan = nscan, 
                          print = FALSE, plot = FALSE)
    
    beta_hat_amen <- median(fit_amen$BETA[,2])
    beta_ci_lower_amen <- quantile(fit_amen$BETA[,2], 0.025)
    beta_ci_upper_amen <- quantile(fit_amen$BETA[,2], 0.975)
    beta_covered_amen <- (beta >= beta_ci_lower_amen) & (beta <= beta_ci_upper_amen)
    
    mu_hat_amen <- median(fit_amen$BETA[,1])
    mu_ci_lower_amen <- quantile(fit_amen$BETA[,1], 0.025)
    mu_ci_upper_amen <- quantile(fit_amen$BETA[,1], 0.975)
    mu_covered_amen <- (mu >= mu_ci_lower_amen) & (mu <= mu_ci_upper_amen)
    
    results$amen <- list(
      beta_hat = beta_hat_amen,
      beta_covered = beta_covered_amen,
      mu_hat = mu_hat_amen,
      mu_covered = mu_covered_amen,
      sigma2_hat = if(!is.null(fit_amen$VC)) median(fit_amen$VC[,"ve"]) else NA
    )
  }
  
  # Fit with lame package (ame function)
  fit_lame <- lame::ame(Y, Xdyad = X, R = R, family = "normal",
                        rvar = rvar, cvar = cvar, dcor = FALSE,
                        burn = burn, nscan = nscan,
                        print = FALSE, plot = FALSE)
  
  beta_hat_lame <- median(fit_lame$BETA[,2])
  beta_ci_lower_lame <- quantile(fit_lame$BETA[,2], 0.025)
  beta_ci_upper_lame <- quantile(fit_lame$BETA[,2], 0.975)
  beta_covered_lame <- (beta >= beta_ci_lower_lame) & (beta <= beta_ci_upper_lame)
  
  mu_hat_lame <- median(fit_lame$BETA[,1])
  mu_ci_lower_lame <- quantile(fit_lame$BETA[,1], 0.025)
  mu_ci_upper_lame <- quantile(fit_lame$BETA[,1], 0.975)
  mu_covered_lame <- (mu >= mu_ci_lower_lame) & (mu <= mu_ci_upper_lame)
  
  # Extract sigma2 from lame (might be in different locations)
  sigma2_lame <- NA
  if(!is.null(fit_lame$s2)) {
    sigma2_lame <- median(fit_lame$s2)
  } else if(!is.null(fit_lame$VC) && "va" %in% colnames(fit_lame$VC)) {
    sigma2_lame <- median(fit_lame$VC[,"va"])
  } else if(!is.null(fit_lame$VC) && "ve" %in% colnames(fit_lame$VC)) {
    sigma2_lame <- median(fit_lame$VC[,"ve"])
  }
  
  results$lame <- list(
    beta_hat = beta_hat_lame,
    beta_covered = beta_covered_lame,
    mu_hat = mu_hat_lame,
    mu_covered = mu_covered_lame,
    sigma2_hat = sigma2_lame
  )
  
  return(results)
}

# Run comparison tests
test_that("Compare amen and lame packages - covariates only", {
  skip_if_not_installed("amen")
  
  set.seed(6886)
  n_sims <- 30
  
  results <- lapply(1:n_sims, function(i) {
    compare_amen_lame(seed = i, n = 40, mu = 0, beta = 1,
                     rvar = FALSE, cvar = FALSE, R = 0,
                     burn = 300, nscan = 1000)
  })
  
  # Extract coverage for amen
  beta_coverage_amen <- mean(sapply(results, function(x) x$amen$beta_covered))
  mu_coverage_amen <- mean(sapply(results, function(x) x$amen$mu_covered))
  
  # Extract coverage for lame
  beta_coverage_lame <- mean(sapply(results, function(x) x$lame$beta_covered))
  mu_coverage_lame <- mean(sapply(results, function(x) x$lame$mu_covered))
  
  cat("\n=== Covariates Only Model ===\n")
  cat("AMEN Beta coverage:", beta_coverage_amen, "\n")
  cat("LAME Beta coverage:", beta_coverage_lame, "\n")
  cat("AMEN Mu coverage:", mu_coverage_amen, "\n")
  cat("LAME Mu coverage:", mu_coverage_lame, "\n")
  
  # Both should have good coverage
  expect_gt(beta_coverage_amen, 0.85)
  expect_gt(beta_coverage_lame, 0.85)
})

test_that("Compare amen and lame packages - with additive effects", {
  skip_if_not_installed("amen")
  
  set.seed(6886)
  n_sims <- 30
  
  results <- lapply(1:n_sims, function(i) {
    compare_amen_lame(seed = i, n = 40, mu = 0, beta = 1,
                     rvar = TRUE, cvar = TRUE, R = 0,
                     burn = 300, nscan = 1000)
  })
  
  # Extract coverage for amen
  beta_coverage_amen <- mean(sapply(results, function(x) x$amen$beta_covered))
  mu_coverage_amen <- mean(sapply(results, function(x) x$amen$mu_covered))
  
  # Extract coverage for lame
  beta_coverage_lame <- mean(sapply(results, function(x) x$lame$beta_covered))
  mu_coverage_lame <- mean(sapply(results, function(x) x$lame$mu_covered))
  
  # Extract bias
  beta_bias_amen <- mean(sapply(results, function(x) x$amen$beta_hat)) - 1
  beta_bias_lame <- mean(sapply(results, function(x) x$lame$beta_hat)) - 1
  
  cat("\n=== Additive Effects Model ===\n")
  cat("AMEN Beta coverage:", beta_coverage_amen, "\n")
  cat("LAME Beta coverage:", beta_coverage_lame, "\n")
  cat("AMEN Beta bias:", beta_bias_amen, "\n")
  cat("LAME Beta bias:", beta_bias_lame, "\n")
  cat("AMEN Mu coverage:", mu_coverage_amen, "\n")
  cat("LAME Mu coverage:", mu_coverage_lame, "\n")
  
  # Check if amen has better coverage
  if (beta_coverage_amen > beta_coverage_lame + 0.05) {
    warning("AMEN has notably better beta coverage than LAME for additive effects model")
  }
})

test_that("Compare amen and lame packages - full model", {
  skip_if_not_installed("amen")
  
  set.seed(6886)
  n_sims <- 30
  
  results <- lapply(1:n_sims, function(i) {
    compare_amen_lame(seed = i, n = 40, mu = 0, beta = 1,
                     rvar = TRUE, cvar = TRUE, R = 2,
                     burn = 300, nscan = 1000)
  })
  
  # Extract coverage for amen
  beta_coverage_amen <- mean(sapply(results, function(x) x$amen$beta_covered))
  mu_coverage_amen <- mean(sapply(results, function(x) x$amen$mu_covered))
  
  # Extract coverage for lame
  beta_coverage_lame <- mean(sapply(results, function(x) x$lame$beta_covered))
  mu_coverage_lame <- mean(sapply(results, function(x) x$lame$mu_covered))
  
  cat("\n=== Full Model (Additive + Multiplicative) ===\n")
  cat("AMEN Beta coverage:", beta_coverage_amen, "\n")
  cat("LAME Beta coverage:", beta_coverage_lame, "\n")
  cat("AMEN Mu coverage:", mu_coverage_amen, "\n")
  cat("LAME Mu coverage:", mu_coverage_lame, "\n")
  
  # Check relative performance
  if (beta_coverage_amen > beta_coverage_lame + 0.05) {
    warning("AMEN has notably better beta coverage than LAME for full model")
  }
})

# Single detailed comparison test
test_that("Detailed single run comparison", {
  skip_if_not_installed("amen")
  
  set.seed(6886)
  n <- 50
  mu_true <- 0
  beta_true <- 1
  
  # Generate data
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Add additive effects
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  eta <- mu_true + beta_true * X + 
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  
  Y <- eta + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Fit both models
  fit_amen <- amen::ame(Y, Xdyad = X, R = 0, family = "nrm",
                       rvar = TRUE, cvar = TRUE, dcor = FALSE,
                       burn = 500, nscan = 2000,
                       print = FALSE, plot = FALSE)
  
  fit_lame <- lame::ame(Y, Xdyad = X, R = 0, family = "normal",
                       rvar = TRUE, cvar = TRUE, dcor = FALSE,
                       burn = 500, nscan = 2000,
                       print = FALSE, plot = FALSE)
  
  # Compare estimates
  beta_amen <- median(fit_amen$BETA[,2])
  beta_lame <- median(fit_lame$BETA[,2])
  
  cat("\n=== Detailed Comparison ===\n")
  cat("True beta:", beta_true, "\n")
  cat("AMEN beta estimate:", beta_amen, "\n")
  cat("LAME beta estimate:", beta_lame, "\n")
  cat("AMEN beta SD:", sd(fit_amen$BETA[,2]), "\n")
  cat("LAME beta SD:", sd(fit_lame$BETA[,2]), "\n")
  
  # Check if estimates are similar
  expect_lt(abs(beta_amen - beta_lame), 0.2)
})