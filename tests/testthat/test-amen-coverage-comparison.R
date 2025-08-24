# Direct comparison of coverage between amen and lame packages
library(testthat)

# Check if amen is available
if (!requireNamespace("amen", quietly = TRUE)) {
  message("Installing amen package for comparison...")
  install.packages("amen", repos = "https://cloud.r-project.org")
}

library(amen)
library(lame)

test_that("Compare coverage rates between amen and lame - additive effects", {
  set.seed(6886)
  n_sims <- 50  # More simulations for better estimate
  n <- 40
  mu_true <- 0
  beta_true <- 1
  
  coverage_amen <- numeric(n_sims)
  coverage_lame <- numeric(n_sims)
  se_amen <- numeric(n_sims)
  se_lame <- numeric(n_sims)
  beta_amen <- numeric(n_sims)
  beta_lame <- numeric(n_sims)
  
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
    
    Y <- eta + matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    # Fit with original amen
    fit_amen <- suppressWarnings(
      amen::ame(Y, Xdyad = X, R = 0, family = "nrm",
                rvar = TRUE, cvar = TRUE, dcor = FALSE,
                burn = 400, nscan = 1500,
                print = FALSE, plot = FALSE)
    )
    
    # Fit with lame
    fit_lame <- lame::ame(Y, Xdyad = X, R = 0, family = "normal",
                         rvar = TRUE, cvar = TRUE, dcor = FALSE,
                         burn = 400, nscan = 1500,
                         print = FALSE, plot = FALSE)
    
    # Extract results
    beta_ci_amen <- quantile(fit_amen$BETA[,2], c(0.025, 0.975))
    beta_ci_lame <- quantile(fit_lame$BETA[,2], c(0.025, 0.975))
    
    coverage_amen[i] <- (beta_true >= beta_ci_amen[1]) & (beta_true <= beta_ci_amen[2])
    coverage_lame[i] <- (beta_true >= beta_ci_lame[1]) & (beta_true <= beta_ci_lame[2])
    
    se_amen[i] <- sd(fit_amen$BETA[,2])
    se_lame[i] <- sd(fit_lame$BETA[,2])
    
    beta_amen[i] <- median(fit_amen$BETA[,2])
    beta_lame[i] <- median(fit_lame$BETA[,2])
    
    if(i %% 10 == 0) cat("Completed", i, "simulations\n")
  }
  
  # Calculate summary statistics
  cov_rate_amen <- mean(coverage_amen)
  cov_rate_lame <- mean(coverage_lame)
  
  mean_se_amen <- mean(se_amen)
  mean_se_lame <- mean(se_lame)
  
  empirical_se_amen <- sd(beta_amen)
  empirical_se_lame <- sd(beta_lame)
  
  cat("\n=== COVERAGE COMPARISON - ADDITIVE EFFECTS MODEL ===\n")
  cat("Sample size: n =", n, ", Simulations:", n_sims, "\n\n")
  
  cat("Coverage rates:\n")
  cat("  AMEN:", sprintf("%.1f%%", cov_rate_amen * 100), "\n")
  cat("  LAME:", sprintf("%.1f%%", cov_rate_lame * 100), "\n")
  cat("  Difference:", sprintf("%.1f%%", (cov_rate_amen - cov_rate_lame) * 100), "\n\n")
  
  cat("Standard error comparison:\n")
  cat("  Mean estimated SE:\n")
  cat("    AMEN:", sprintf("%.4f", mean_se_amen), "\n")
  cat("    LAME:", sprintf("%.4f", mean_se_lame), "\n")
  cat("  Empirical SE of estimates:\n")
  cat("    AMEN:", sprintf("%.4f", empirical_se_amen), "\n")
  cat("    LAME:", sprintf("%.4f", empirical_se_lame), "\n")
  cat("  SE ratio (empirical/estimated):\n")
  cat("    AMEN:", sprintf("%.3f", empirical_se_amen/mean_se_amen), "\n")
  cat("    LAME:", sprintf("%.3f", empirical_se_lame/mean_se_lame), "\n\n")
  
  # Check if both have the same issue
  if(cov_rate_amen < 0.9 && cov_rate_lame < 0.9) {
    cat("⚠️ BOTH packages show undercoverage with additive effects\n")
  } else if(cov_rate_amen > cov_rate_lame + 0.05) {
    cat("⚠️ LAME has worse coverage than AMEN\n")
  } else if(abs(cov_rate_amen - cov_rate_lame) < 0.05) {
    cat("✓ Both packages have similar coverage\n")
  }
})

# Test without additive effects for comparison
test_that("Compare coverage rates between amen and lame - simple model", {
  set.seed(6886)
  n_sims <- 30
  n <- 40
  mu_true <- 0
  beta_true <- 1
  
  coverage_amen <- numeric(n_sims)
  coverage_lame <- numeric(n_sims)
  
  for(i in 1:n_sims) {
    set.seed(i)
    
    # Generate simple data (no additive effects)
    xw <- matrix(rnorm(n*2), n, 2)
    X <- tcrossprod(xw[,1])
    diag(X) <- NA
    
    eta <- mu_true + beta_true * X
    Y <- eta + matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    # Fit with original amen
    fit_amen <- suppressWarnings(
      amen::ame(Y, Xdyad = X, R = 0, family = "nrm",
                rvar = FALSE, cvar = FALSE, dcor = FALSE,
                burn = 300, nscan = 1000,
                print = FALSE, plot = FALSE)
    )
    
    # Fit with lame
    fit_lame <- lame::ame(Y, Xdyad = X, R = 0, family = "normal",
                         rvar = FALSE, cvar = FALSE, dcor = FALSE,
                         burn = 300, nscan = 1000,
                         print = FALSE, plot = FALSE)
    
    # Extract results
    beta_ci_amen <- quantile(fit_amen$BETA[,2], c(0.025, 0.975))
    beta_ci_lame <- quantile(fit_lame$BETA[,2], c(0.025, 0.975))
    
    coverage_amen[i] <- (beta_true >= beta_ci_amen[1]) & (beta_true <= beta_ci_amen[2])
    coverage_lame[i] <- (beta_true >= beta_ci_lame[1]) & (beta_true <= beta_ci_lame[2])
  }
  
  cov_rate_amen <- mean(coverage_amen)
  cov_rate_lame <- mean(coverage_lame)
  
  cat("\n=== COVERAGE COMPARISON - SIMPLE MODEL (NO EFFECTS) ===\n")
  cat("Coverage rates:\n")
  cat("  AMEN:", sprintf("%.1f%%", cov_rate_amen * 100), "\n")
  cat("  LAME:", sprintf("%.1f%%", cov_rate_lame * 100), "\n")
  cat("  Difference:", sprintf("%.1f%%", (cov_rate_amen - cov_rate_lame) * 100), "\n\n")
  
  if(cov_rate_amen > 0.85 && cov_rate_lame > 0.85) {
    cat("✓ Both packages have good coverage for simple models\n")
  }
})