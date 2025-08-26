# Diagnostic test to understand negative predictions in count models
library(lame)
library(testthat)

test_that("Diagnose negative predictions in Poisson models", {
  set.seed(6886)
  n <- 25
  
  # Simple case - just intercept and one covariate
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  
  # Generate count data
  eta <- 1.0 + 0.5 * X  # Log scale
  lambda <- exp(eta)
  lambda[lambda > 30] <- 30
  Y <- matrix(rpois(n*n, as.vector(lambda)), n, n)
  diag(Y) <- NA
  
  cat("\n=== Input Data Diagnostics ===\n")
  cat("Y range:", range(Y, na.rm=TRUE), "\n")
  cat("Y mean:", mean(Y, na.rm=TRUE), "\n")
  cat("Y variance:", var(c(Y), na.rm=TRUE), "\n")
  cat("Proportion of zeros:", mean(Y == 0, na.rm=TRUE), "\n")
  
  # Fit simple model
  fit <- ame(Y, Xdyad=X, R=0, family="poisson",
            rvar=FALSE, cvar=FALSE, dcor=FALSE,
            burn=200, nscan=800, print=FALSE)
  
  cat("\n=== Model Output Diagnostics ===\n")
  cat("Beta estimate:", median(fit$BETA[,2]), "\n")
  cat("Intercept estimate:", median(fit$BETA[,1]), "\n")
  
  # Check EZ (expected values)
  cat("\nEZ (predictions) diagnostics:\n")
  EZ <- reconstruct_EZ(fit)
  cat("EZ range:", range(EZ, na.rm=TRUE), "\n")
  cat("Number of negative values:", sum(EZ < 0, na.rm=TRUE), "\n")
  cat("Proportion negative:", mean(EZ < 0, na.rm=TRUE), "\n")
  
  if(any(EZ < 0, na.rm=TRUE)) {
    neg_vals <- EZ[EZ < 0 & !is.na(EZ)]
    cat("Negative values:", head(neg_vals, 10), "\n")
    cat("Most negative:", min(EZ, na.rm=TRUE), "\n")
  }
  
  # Check if it's a transformation issue
  # For Poisson, EZ should be on the count scale (not log scale)
  cat("\nChecking scale of predictions:\n")
  cat("Mean of Y:", mean(Y, na.rm=TRUE), "\n")
  cat("Mean of EZ:", mean(EZ, na.rm=TRUE), "\n")
  cat("Median of Y:", median(Y, na.rm=TRUE), "\n")
  cat("Median of EZ:", median(EZ, na.rm=TRUE), "\n")
})

test_that("Check internal Z matrix for Poisson", {
  set.seed(6886)
  n <- 20
  
  # Very simple network
  Y <- matrix(rpois(n*n, 3), n, n)
  diag(Y) <- NA
  
  # Fit and look at internals
  fit <- ame(Y, R=0, family="poisson",
            burn=100, nscan=400, print=FALSE)
  
  # The issue might be in how Z (latent variable) is handled
  # For Poisson, there might be a log-link involved
  
  cat("\n=== Checking for log-scale vs count-scale issue ===\n")
  cat("If EZ is on log scale, exp(EZ) should match Y better:\n")
  EZ <- reconstruct_EZ(fit)
  cat("Correlation Y vs EZ:", cor(c(Y), c(EZ), use="complete.obs"), "\n")
  cat("Correlation Y vs exp(EZ):", cor(c(Y), exp(c(EZ)), use="complete.obs"), "\n")
  
  # Check if negative values become reasonable after exp()
  if(any(EZ < 0, na.rm=TRUE)) {
    cat("\nTransforming negative values:\n")
    neg_idx <- which(EZ < 0 & !is.na(EZ))
    cat("Example negative EZ values:", EZ[neg_idx[1:min(5, length(neg_idx))]], "\n")
    cat("After exp():", exp(EZ[neg_idx[1:min(5, length(neg_idx))]]), "\n")
  }
})

test_that("Compare with sparse vs dense count networks", {
  set.seed(6886)
  n <- 25
  
  # Sparse network
  Y_sparse <- matrix(0, n, n)
  nonzero <- sample(1:(n*n), size=floor(0.2*n*n))
  Y_sparse[nonzero] <- rpois(length(nonzero), 2)
  diag(Y_sparse) <- NA
  
  fit_sparse <- ame(Y_sparse, R=0, family="poisson",
                   burn=200, nscan=800, print=FALSE)
  
  cat("\n=== Sparse Network ===\n")
  cat("Proportion zeros in Y:", mean(Y_sparse == 0, na.rm=TRUE), "\n")
  EZ_sparse <- reconstruct_EZ(fit_sparse)
  cat("Proportion negative in EZ:", mean(EZ_sparse < 0, na.rm=TRUE), "\n")
  cat("Min EZ:", min(EZ_sparse, na.rm=TRUE), "\n")
  
  # Dense network
  Y_dense <- matrix(rpois(n*n, 5), n, n)
  diag(Y_dense) <- NA
  
  fit_dense <- ame(Y_dense, R=0, family="poisson",
                  burn=200, nscan=800, print=FALSE)
  
  cat("\n=== Dense Network ===\n")
  cat("Proportion zeros in Y:", mean(Y_dense == 0, na.rm=TRUE), "\n")
  EZ_dense <- reconstruct_EZ(fit_dense)
  cat("Proportion negative in EZ:", mean(EZ_dense < 0, na.rm=TRUE), "\n")
  cat("Min EZ:", min(EZ_dense, na.rm=TRUE), "\n")
  
  # Pattern: sparse networks likely have more negative predictions
  # EZ removed - check using reconstruct_EZ instead
  EZ_sparse <- reconstruct_EZ(fit_sparse)
  EZ_dense <- reconstruct_EZ(fit_dense)
  expect_true(mean(EZ_sparse < 0, na.rm=TRUE) >= 
              mean(EZ_dense < 0, na.rm=TRUE))
})