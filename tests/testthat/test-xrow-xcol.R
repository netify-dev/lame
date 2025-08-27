# Test xrow and xcol covariates for static Gaussian models
library(lame)
library(testthat)

test_that("xrow/xcol work correctly for static unipartite Gaussian models", {
  set.seed(123)
  n <- 30
  
  # Generate row and column covariates
  Xrow <- matrix(rnorm(n * 2), n, 2)
  colnames(Xrow) <- c("row_cov1", "row_cov2")
  Xcol <- matrix(rnorm(n * 2), n, 2)
  colnames(Xcol) <- c("col_cov1", "col_cov2")
  
  # True coefficients
  beta_row <- c(0.5, -0.3)
  beta_col <- c(0.8, 0.4)
  intercept <- 0.2
  
  # Generate Y based on row and column effects
  Y <- matrix(intercept, n, n)
  
  # Add row effects (constant across columns)
  for(i in 1:n) {
    Y[i, ] <- Y[i, ] + Xrow[i, 1] * beta_row[1] + Xrow[i, 2] * beta_row[2]
  }
  
  # Add column effects (constant across rows)
  for(j in 1:n) {
    Y[, j] <- Y[, j] + Xcol[j, 1] * beta_col[1] + Xcol[j, 2] * beta_col[2]
  }
  
  # Add noise
  Y <- Y + matrix(rnorm(n * n, 0, 0.5), n, n)
  diag(Y) <- NA
  
  # Fit model with row and column covariates
  fit <- ame(Y, Xrow = Xrow, Xcol = Xcol, 
            R = 0, family = "normal",
            burn = 200, nscan = 800, odens = 25,
            print = FALSE)
  
  # Check dimensions of BETA
  # Should have: intercept + 2 row + 2 col = 5 coefficients
  expect_equal(ncol(fit$BETA), 5)
  
  # Extract coefficient estimates
  beta_est <- colMeans(fit$BETA)
  
  # Check coefficient names if available
  if(!is.null(colnames(fit$BETA))) {
    expect_true("intercept" %in% colnames(fit$BETA) || 
                colnames(fit$BETA)[1] == "(Intercept)")
    expect_true(any(grepl("row_cov1", colnames(fit$BETA))))
    expect_true(any(grepl("col_cov1", colnames(fit$BETA))))
  }
  
  # Coefficients should be reasonably close to truth
  # Note: order is typically intercept, row covs, col covs
  expect_equal(unname(beta_est[1]), intercept, tolerance = 0.5, 
              label = "Intercept estimate")
  expect_equal(unname(beta_est[2]), beta_row[1], tolerance = 0.5,
              label = "Row covariate 1 estimate")
  expect_equal(unname(beta_est[3]), beta_row[2], tolerance = 0.5,
              label = "Row covariate 2 estimate")
  expect_equal(unname(beta_est[4]), beta_col[1], tolerance = 0.5,
              label = "Column covariate 1 estimate")
  expect_equal(unname(beta_est[5]), beta_col[2], tolerance = 0.5,
              label = "Column covariate 2 estimate")
  
  # Check that X array is properly constructed
  expect_false(is.null(fit$X))
  expect_equal(dim(fit$X), c(n, n, 5))
  
  # First layer should be intercept (all 1s)
  expect_true(all(fit$X[,,1] == 1, na.rm = TRUE))
  
  # Row covariate layers should be constant across columns
  for(k in 2:3) {
    for(i in 1:n) {
      row_vals <- fit$X[i, , k]
      expect_true(length(unique(row_vals[!is.na(row_vals)])) == 1,
                 info = paste("Row covariate", k-1, "not constant across columns"))
    }
  }
  
  # Column covariate layers should be constant across rows
  for(k in 4:5) {
    for(j in 1:n) {
      col_vals <- fit$X[, j, k]
      expect_true(length(unique(col_vals[!is.na(col_vals)])) == 1,
                 info = paste("Column covariate", k-3, "not constant across rows"))
    }
  }
  
  # Test prediction with row/col covariates
  pred <- predict(fit, type = "link")
  expect_equal(dim(pred), c(n, n))
  
  # Predictions should correlate with true linear predictor
  true_eta <- Y - matrix(rnorm(n * n, 0, 0.5), n, n)  # Approximate
  diag(true_eta) <- NA
  corr <- cor(c(true_eta), c(pred), use = "pairwise.complete.obs")
  expect_gt(corr, 0.7, label = "Predictions should correlate with truth")
})

test_that("xrow/xcol work correctly for static bipartite Gaussian models", {
  set.seed(456)
  nA <- 20  # Row nodes
  nB <- 15  # Column nodes
  
  # Generate row covariates (for set A)
  Xrow <- matrix(rnorm(nA * 2), nA, 2)
  colnames(Xrow) <- c("rowA_cov1", "rowA_cov2")
  
  # Generate column covariates (for set B)
  Xcol <- matrix(rnorm(nB * 3), nB, 3)
  colnames(Xcol) <- c("colB_cov1", "colB_cov2", "colB_cov3")
  
  # True coefficients
  beta_row <- c(0.6, -0.4)
  beta_col <- c(0.3, 0.7, -0.2)
  intercept <- -0.1
  
  # Generate Y based on effects
  Y <- matrix(intercept, nA, nB)
  
  # Add row effects (from set A)
  for(i in 1:nA) {
    Y[i, ] <- Y[i, ] + Xrow[i, 1] * beta_row[1] + Xrow[i, 2] * beta_row[2]
  }
  
  # Add column effects (from set B)
  for(j in 1:nB) {
    Y[, j] <- Y[, j] + Xcol[j, 1] * beta_col[1] + 
               Xcol[j, 2] * beta_col[2] + Xcol[j, 3] * beta_col[3]
  }
  
  # Add noise
  Y <- Y + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
  
  # Fit bipartite model
  fit <- ame(Y, Xrow = Xrow, Xcol = Xcol,
            mode = "bipartite",
            R_row = 0, R_col = 0, family = "normal",
            burn = 200, nscan = 800, odens = 25,
            print = FALSE)
  
  # Check dimensions of BETA
  # Should have: intercept + 2 row + 3 col = 6 coefficients
  expect_equal(ncol(fit$BETA), 6)
  
  # Extract estimates
  beta_est <- colMeans(fit$BETA)
  
  # Check recovery (wider tolerance for bipartite)
  expect_equal(beta_est[1], intercept, tolerance = 0.5,
              label = "Bipartite intercept")
  expect_equal(beta_est[2], beta_row[1], tolerance = 0.5,
              label = "Bipartite row cov 1")
  expect_equal(beta_est[3], beta_row[2], tolerance = 0.5,
              label = "Bipartite row cov 2")
  expect_equal(beta_est[4], beta_col[1], tolerance = 0.5,
              label = "Bipartite col cov 1")
  expect_equal(beta_est[5], beta_col[2], tolerance = 0.5,
              label = "Bipartite col cov 2")
  expect_equal(beta_est[6], beta_col[3], tolerance = 0.5,
              label = "Bipartite col cov 3")
  
  # Check X array structure
  expect_false(is.null(fit$X))
  expect_equal(dim(fit$X), c(nA, nB, 6))
  
  # Verify bipartite-specific attributes
  expect_equal(fit$nA, nA)
  expect_equal(fit$nB, nB)
  expect_equal(fit$mode, "bipartite")
  
  # Test reconstruction
  EZ <- reconstruct_EZ(fit)
  expect_equal(dim(EZ), c(nA, nB))
  
  # Should correlate with truth
  true_eta <- Y - matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
  corr <- cor(c(true_eta), c(EZ), use = "pairwise.complete.obs")
  expect_gt(corr, 0.7)
})

test_that("Mixed dyadic and nodal covariates work together", {
  set.seed(789)
  n <- 25
  
  # Generate all three types of covariates
  Xdyad <- matrix(rnorm(n * n), n, n)
  diag(Xdyad) <- NA
  Xrow <- matrix(rnorm(n), n, 1)
  Xcol <- matrix(rnorm(n), n, 1)
  
  # True coefficients
  beta_dyad <- 1.0
  beta_row <- 0.5
  beta_col <- -0.5
  intercept <- 0.0
  
  # Generate Y
  Y <- intercept + beta_dyad * Xdyad
  
  # Add row effect
  for(i in 1:n) {
    Y[i, ] <- Y[i, ] + Xrow[i] * beta_row
  }
  
  # Add column effect
  for(j in 1:n) {
    Y[, j] <- Y[, j] + Xcol[j] * beta_col
  }
  
  # Add noise
  Y <- Y + matrix(rnorm(n * n, 0, 0.3), n, n)
  diag(Y) <- NA
  
  # Fit model with all covariate types
  fit <- ame(Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
            R = 0, family = "normal",
            burn = 200, nscan = 800, odens = 25,
            print = FALSE)
  
  # Should have 4 coefficients
  expect_equal(ncol(fit$BETA), 4)
  
  beta_est <- colMeans(fit$BETA)
  
  # All coefficients should be recovered  
  # Note: order might be intercept, row, col, dyad based on design_array
  expect_equal(unname(beta_est[1]), intercept, tolerance = 0.5)
  # The remaining coefficients may be in different order
  # Check that all true values are recovered somewhere
  all_true <- c(beta_dyad, beta_row, beta_col)
  all_est <- unname(beta_est[2:4])
  
  # Each true coefficient should match one estimated coefficient
  expect_true(any(abs(all_est - beta_dyad) < 0.3), 
              info = "Dyadic coefficient not recovered")
  expect_true(any(abs(all_est - beta_row) < 0.5),
              info = "Row coefficient not recovered")
  expect_true(any(abs(all_est - beta_col) < 0.5),
              info = "Column coefficient not recovered")
  
  # Verify design array construction
  expect_equal(dim(fit$X), c(n, n, 4))
})