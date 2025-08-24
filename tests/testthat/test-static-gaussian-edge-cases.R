# Edge cases and additional tests for static Gaussian models
library(lame)
library(testthat)

# Test 1: Missing data handling
test_that("Gaussian AME handles missing data correctly", {
  set.seed(6886)
  n <- 30
  
  # Generate complete data
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  Y_complete <- 0.5 + 1.2 * X + matrix(rnorm(n*n), n, n)
  
  # Add 20% missing values
  Y_missing <- Y_complete
  n_missing <- floor(0.2 * n * (n-1))
  missing_idx <- sample(which(!is.na(Y_complete)), n_missing)
  Y_missing[missing_idx] <- NA
  diag(Y_missing) <- NA
  
  # Fit model
  fit <- ame(Y_missing, Xdyad = X, R = 1, family = "normal",
            rvar = TRUE, cvar = TRUE,
            burn = 200, nscan = 800, print = FALSE, plot = FALSE)
  
  # Check that predictions exist for missing values
  expect_true(!is.null(fit$EZ))
  expect_true(all(!is.na(fit$EZ[missing_idx])))
  
  # Check that predictions are reasonable
  pred_corr <- cor(c(Y_complete[missing_idx]), c(fit$EZ[missing_idx]))
  expect_gt(pred_corr, 0.3)  # Should have some predictive power
})

# Test 2: Different network sizes
test_that("Gaussian AME works across different network sizes", {
  set.seed(6886)
  
  for(n in c(10, 25, 50)) {
    X <- matrix(rnorm(n*n), n, n)
    diag(X) <- NA
    Y <- 1 * X + matrix(rnorm(n*n), n, n)
    diag(Y) <- NA
    
    fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
              burn = 100, nscan = 400, print = FALSE, plot = FALSE)
    
    beta_est <- median(fit$BETA[,2])
    expect_lt(abs(beta_est - 1), 0.5)
  }
})

# Test 3: Multiple covariates
test_that("Gaussian AME handles multiple covariates correctly", {
  set.seed(6886)
  n <- 30
  p <- 3  # Number of covariates
  
  # Generate multiple covariates
  X_array <- array(rnorm(n*n*p), dim = c(n, n, p))
  for(k in 1:p) {
    diag(X_array[,,k]) <- NA
  }
  
  # True coefficients
  beta_true <- c(0.5, -0.3, 0.8)
  
  # Generate outcome
  Y <- matrix(0, n, n)
  for(k in 1:p) {
    Y <- Y + beta_true[k] * X_array[,,k]
  }
  Y <- Y + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad = X_array, R = 0, family = "normal",
            burn = 300, nscan = 1200, print = FALSE, plot = FALSE)
  
  # Check all coefficients are recovered (skip intercept)
  for(k in 1:p) {
    beta_est <- median(fit$BETA[, k+1])
    expect_lt(abs(beta_est - beta_true[k]), 0.3)
  }
})

# Test 4: Symmetric networks
test_that("Gaussian AME handles symmetric networks", {
  set.seed(6886)
  n <- 25
  
  # Generate symmetric network
  X <- matrix(rnorm(n*n), n, n)
  X <- (X + t(X))/2  # Make symmetric
  diag(X) <- NA
  
  # Generate symmetric outcome with single variance component
  a_true <- rnorm(n, 0, 0.5)
  Y <- 0.5 + 0.8 * X + outer(a_true, rep(1,n)) + outer(rep(1,n), a_true)
  Y <- (Y + t(Y))/2  # Ensure symmetry
  Y <- Y + matrix(rnorm(n*n, 0, 0.5), n, n)
  Y <- (Y + t(Y))/2
  diag(Y) <- NA
  
  # Fit symmetric model
  fit <- ame(Y, Xdyad = X, R = 1, family = "normal",
            symmetric = TRUE, nvar = TRUE,
            burn = 200, nscan = 800, print = FALSE, plot = FALSE)
  
  # Check symmetry is preserved in predictions
  expect_true(isSymmetric(fit$EZ, check.attributes = FALSE, tol = 0.1))
  
  # For symmetric model, U and L should exist (not V)
  expect_true(!is.null(fit$U))
  if(!is.null(fit$L)) {
    expect_true(is.numeric(fit$L))
  }
})

# Test 5: High R (many latent dimensions)
test_that("Gaussian AME handles higher dimensional latent space", {
  set.seed(6886)
  n <- 40
  R_true <- 3
  
  # Generate data with R=3 latent dimensions
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  U_true <- matrix(rnorm(n*R_true), n, R_true)
  V_true <- matrix(rnorm(n*R_true), n, R_true)
  
  Y <- 1 * X + tcrossprod(U_true, V_true) + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Fit with R=3
  fit <- ame(Y, Xdyad = X, R = 3, family = "normal",
            burn = 300, nscan = 1200, print = FALSE, plot = FALSE)
  
  # Check dimensions
  expect_equal(ncol(fit$U), 3)
  expect_equal(ncol(fit$V), 3)
  
  # Check that UVPM captures variance
  var_explained <- var(c(fit$UVPM)) / var(c(Y), na.rm = TRUE)
  expect_gt(var_explained, 0.1)  # Should explain some variance
})

# Test 6: Variance components estimation
test_that("Gaussian AME estimates variance components correctly", {
  set.seed(6886)
  n <- 30
  
  # Known variance components
  var_row <- 0.5
  var_col <- 0.3
  var_dyad <- 1.0
  
  # Generate data
  a_true <- rnorm(n, 0, sqrt(var_row))
  b_true <- rnorm(n, 0, sqrt(var_col))
  
  Y <- outer(a_true, rep(1,n)) + outer(rep(1,n), b_true) + 
       matrix(rnorm(n*n, 0, sqrt(var_dyad)), n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, R = 0, family = "normal",
            rvar = TRUE, cvar = TRUE,
            burn = 500, nscan = 2000, print = FALSE, plot = FALSE)
  
  # Extract variance estimates
  if(!is.null(fit$VC) && ncol(fit$VC) >= 3) {
    var_row_est <- median(fit$VC[,1])
    var_col_est <- median(fit$VC[,3])
    
    # Check order of magnitude is correct
    expect_lt(var_row_est, var_row * 3)
    expect_gt(var_row_est, var_row / 3)
    expect_lt(var_col_est, var_col * 3) 
    expect_gt(var_col_est, var_col / 3)
  }
})

# Test 7: Extreme values handling
test_that("Gaussian AME handles extreme values gracefully", {
  set.seed(6886)
  n <- 25
  
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  # Add some extreme values
  Y <- 1 * X + matrix(rnorm(n*n), n, n)
  Y[1,2] <- 100  # Extreme outlier
  Y[3,4] <- -100
  diag(Y) <- NA
  
  # Should not crash
  expect_error(
    fit <- ame(Y, Xdyad = X, R = 1, family = "normal",
              burn = 100, nscan = 400, print = FALSE, plot = FALSE),
    NA  # Expect no error
  )
  
  # Should still recover reasonable beta
  if(exists("fit")) {
    beta_est <- median(fit$BETA[,2])
    # More tolerant due to outliers
    expect_lt(abs(beta_est - 1), 1.0)
  }
})

# Test 8: Zero variance networks
test_that("Gaussian AME handles low variance networks", {
  set.seed(6886)
  n <- 20
  
  # Very low variance network
  Y <- matrix(0.1, n, n) + matrix(rnorm(n*n, 0, 0.01), n, n)
  diag(Y) <- NA
  
  # Should not crash
  expect_error(
    fit <- ame(Y, R = 0, family = "normal",
              burn = 100, nscan = 400, print = FALSE, plot = FALSE),
    NA  # Expect no error
  )
})

# Test 9: Convergence diagnostics
test_that("Gaussian AME provides convergence diagnostics", {
  skip_if_not_installed("coda")
  
  set.seed(6886)
  n <- 20
  
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  Y <- X + matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
            burn = 200, nscan = 1000, print = FALSE, plot = FALSE)
  
  # Check MCMC output
  expect_true(nrow(fit$BETA) > 0)
  
  # Calculate effective sample size
  library(coda)
  ess <- effectiveSize(fit$BETA[,2])
  expect_gt(ess, 50)  # Should have reasonable ESS
})

# Test 10: Prediction accuracy
test_that("Gaussian AME predictions are accurate", {
  set.seed(6886)
  n <- 30
  
  # Generate data
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  eta_true <- 0.5 + 1.2 * X + 
              outer(a_true, rep(1,n)) + outer(rep(1,n), b_true)
  Y <- eta_true + matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
            rvar = TRUE, cvar = TRUE,
            burn = 300, nscan = 1200, print = FALSE, plot = FALSE)
  
  # Check prediction accuracy
  # EZ should be close to true linear predictor
  pred_corr <- cor(c(eta_true[!is.na(Y)]), c(fit$EZ[!is.na(Y)]))
  expect_gt(pred_corr, 0.8)  # Should have high correlation
  
  # Check residuals
  residuals <- Y - fit$EZ
  expect_lt(abs(mean(residuals, na.rm = TRUE)), 0.1)  # Centered
  expect_lt(sd(residuals, na.rm = TRUE), 1.0)  # Reasonable variance
})