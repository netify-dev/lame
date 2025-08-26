# Quick comparison test between amen and lame
library(testthat)

# Check if amen is installed
amen_available <- requireNamespace("amen", quietly = TRUE)

if (amen_available) {
  library(amen)
  library(lame)
  
  test_that("Quick comparison of amen vs lame", {
    set.seed(6886)
    n <- 30
    
    # Generate simple test data
    X <- matrix(rnorm(n*n), n, n)
    diag(X) <- NA
    
    Y <- 0.5 + 1.2 * X + matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    # Run with fewer iterations for speed
    burn <- 100
    nscan <- 500
    
    # Fit with amen
    fit_amen <- amen::ame(Y, Xdyad = X, R = 0, family = "nrm",
                          rvar = FALSE, cvar = FALSE, dcor = FALSE,
                          burn = burn, nscan = nscan,
                          print = FALSE)
    
    # Fit with lame
    fit_lame <- lame::ame(Y, Xdyad = X, R = 0, family = "normal",
                         rvar = FALSE, cvar = FALSE, dcor = FALSE,
                         burn = burn, nscan = nscan,
                         print = FALSE)
    
    # Compare beta estimates
    beta_amen <- median(fit_amen$BETA[,2])
    beta_lame <- median(fit_lame$BETA[,2])
    
    cat("\n=== Quick Comparison Results ===\n")
    cat("True beta: 1.2\n")
    cat("AMEN beta:", beta_amen, "\n")
    cat("LAME beta:", beta_lame, "\n")
    cat("Difference:", abs(beta_amen - beta_lame), "\n")
    
    # Check 95% CI coverage
    ci_amen <- quantile(fit_amen$BETA[,2], c(0.025, 0.975))
    ci_lame <- quantile(fit_lame$BETA[,2], c(0.025, 0.975))
    
    cat("\nAMEN 95% CI:", ci_amen, "\n")
    cat("LAME 95% CI:", ci_lame, "\n")
    
    covered_amen <- (1.2 >= ci_amen[1]) & (1.2 <= ci_amen[2])
    covered_lame <- (1.2 >= ci_lame[1]) & (1.2 <= ci_lame[2])
    
    cat("AMEN covers true value:", covered_amen, "\n")
    cat("LAME covers true value:", covered_lame, "\n")
    
    # Compare variance estimates
    if (!is.null(fit_amen$VC) && !is.null(fit_lame$VC)) {
      ve_amen <- median(fit_amen$VC[,"ve"])
      ve_lame <- if("ve" %in% colnames(fit_lame$VC)) {
        median(fit_lame$VC[,"ve"])
      } else if("va" %in% colnames(fit_lame$VC)) {
        median(fit_lame$VC[,"va"])
      } else {
        NA
      }
      
      cat("\nVariance estimates:\n")
      cat("AMEN ve:", ve_amen, "\n")
      cat("LAME ve:", ve_lame, "\n")
    }
  })
  
  # Test with additive effects
  test_that("Comparison with additive effects", {
    set.seed(6886)
    n <- 30
    
    # Generate data with additive effects
    X <- matrix(rnorm(n*n), n, n)
    diag(X) <- NA
    
    a_true <- rnorm(n, 0, 0.5)
    b_true <- rnorm(n, 0, 0.5)
    
    Y <- 0.5 + 1.0 * X + 
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true) +
         matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    burn <- 200
    nscan <- 800
    
    # Fit with amen
    fit_amen <- amen::ame(Y, Xdyad = X, R = 0, family = "nrm",
                          rvar = TRUE, cvar = TRUE, dcor = FALSE,
                          burn = burn, nscan = nscan,
                          print = FALSE)
    
    # Fit with lame
    fit_lame <- lame::ame(Y, Xdyad = X, R = 0, family = "normal",
                         rvar = TRUE, cvar = TRUE, dcor = FALSE,
                         burn = burn, nscan = nscan,
                         print = FALSE)
    
    beta_amen <- median(fit_amen$BETA[,2])
    beta_lame <- median(fit_lame$BETA[,2])
    
    cat("\n=== Additive Effects Comparison ===\n")
    cat("True beta: 1.0\n")
    cat("AMEN beta:", beta_amen, "\n")
    cat("LAME beta:", beta_lame, "\n")
    
    # Check CI width (as proxy for efficiency)
    ci_width_amen <- diff(quantile(fit_amen$BETA[,2], c(0.025, 0.975)))
    ci_width_lame <- diff(quantile(fit_lame$BETA[,2], c(0.025, 0.975)))
    
    cat("\nCI widths:\n")
    cat("AMEN:", ci_width_amen, "\n")
    cat("LAME:", ci_width_lame, "\n")
    
    # Check if lame CI is much wider (indicating less efficiency)
    if (ci_width_lame > ci_width_amen * 1.2) {
      warning("LAME has notably wider CI than AMEN")
    }
  })
  
} else {
  test_that("amen package availability", {
    skip("amen package not available for comparison")
  })
}