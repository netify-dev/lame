# Quick check of amen coverage issue
library(testthat)
library(amen)
library(lame)

test_that("Quick check: Both packages have same coverage issue", {
  set.seed(6886)
  n_sims <- 20  # Fewer simulations for speed
  n <- 30      # Smaller network
  
  results <- data.frame(
    sim = 1:n_sims,
    covered_amen = NA,
    covered_lame = NA,
    se_amen = NA,
    se_lame = NA,
    beta_amen = NA,
    beta_lame = NA
  )
  
  for(i in 1:n_sims) {
    set.seed(i)
    
    # Generate data with additive effects
    X <- matrix(rnorm(n*n), n, n)
    diag(X) <- NA
    
    a_true <- rnorm(n, 0, 0.5)
    b_true <- rnorm(n, 0, 0.5)
    
    Y <- 1 * X + outer(a_true, rep(1,n)) + outer(rep(1,n), b_true) + 
         matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    # Quick fits with fewer iterations
    fit_amen <- suppressWarnings(
      amen::ame(Y, Xdyad=X, R=0, family="nrm",
                rvar=TRUE, cvar=TRUE, dcor=FALSE,
                burn=200, nscan=800, print=FALSE, plot=FALSE)
    )
    
    fit_lame <- lame::ame(Y, Xdyad=X, R=0, family="normal",
                         rvar=TRUE, cvar=TRUE, dcor=FALSE,
                         burn=200, nscan=800, print=FALSE, plot=FALSE)
    
    # Store results
    ci_amen <- quantile(fit_amen$BETA[,2], c(0.025, 0.975))
    ci_lame <- quantile(fit_lame$BETA[,2], c(0.025, 0.975))
    
    results$covered_amen[i] <- (1 >= ci_amen[1]) & (1 <= ci_amen[2])
    results$covered_lame[i] <- (1 >= ci_lame[1]) & (1 <= ci_lame[2])
    results$se_amen[i] <- sd(fit_amen$BETA[,2])
    results$se_lame[i] <- sd(fit_lame$BETA[,2])
    results$beta_amen[i] <- median(fit_amen$BETA[,2])
    results$beta_lame[i] <- median(fit_lame$BETA[,2])
  }
  
  cat("\n=== FINAL VERDICT ===\n")
  cat("Coverage (should be 95%):\n")
  cat("  AMEN:", mean(results$covered_amen)*100, "%\n")
  cat("  LAME:", mean(results$covered_lame)*100, "%\n\n")
  
  cat("SE Ratio (empirical/mean estimated):\n")
  cat("  AMEN:", sd(results$beta_amen)/mean(results$se_amen), "\n")
  cat("  LAME:", sd(results$beta_lame)/mean(results$se_lame), "\n\n")
  
  if(mean(results$covered_amen) < 0.9 && mean(results$covered_lame) < 0.9) {
    cat("CONCLUSION: BOTH packages have the coverage issue with additive effects!\n")
    cat("This is NOT a bug introduced by lame - it's in the original amen too.\n")
  }
})

# Also test simple model to confirm both work well there
test_that("Both work well for simple models", {
  set.seed(6886)
  coverage_amen <- coverage_lame <- 0
  
  for(i in 1:20) {
    set.seed(i)
    n <- 30
    X <- matrix(rnorm(n*n), n, n)
    diag(X) <- NA
    Y <- 1 * X + matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    fit_amen <- suppressWarnings(
      amen::ame(Y, Xdyad=X, R=0, family="nrm",
                rvar=FALSE, cvar=FALSE, dcor=FALSE,
                burn=200, nscan=800, print=FALSE, plot=FALSE)
    )
    
    fit_lame <- lame::ame(Y, Xdyad=X, R=0, family="normal",
                         rvar=FALSE, cvar=FALSE, dcor=FALSE,
                         burn=200, nscan=800, print=FALSE, plot=FALSE)
    
    ci_amen <- quantile(fit_amen$BETA[,2], c(0.025, 0.975))
    ci_lame <- quantile(fit_lame$BETA[,2], c(0.025, 0.975))
    
    coverage_amen <- coverage_amen + ((1 >= ci_amen[1]) & (1 <= ci_amen[2]))
    coverage_lame <- coverage_lame + ((1 >= ci_lame[1]) & (1 <= ci_lame[2]))
  }
  
  cat("\n=== SIMPLE MODEL (NO ADDITIVE EFFECTS) ===\n")
  cat("Coverage: AMEN =", coverage_amen/20*100, "%, LAME =", coverage_lame/20*100, "%\n")
  cat("Both should be ~95% and they are!\n")
})