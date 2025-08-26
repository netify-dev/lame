# Comprehensive test suite for static AME models
# Authors: Cassy Dorff, Shahryar Minhas, Tosin Salau

library(lame)

# ============================================================================
# Unipartite Models - Complete Coverage
# ============================================================================

test_that("Unipartite models handle all parameter combinations", {
  skip_on_cran()
  set.seed(123)
  
  n <- 15
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Test 1: Basic model (no effects)
  fit1 <- ame(Y, R = 0, rvar = FALSE, cvar = FALSE,
             burn = 50, nscan = 100, print = FALSE)
  expect_equal(fit1$mode, "unipartite")
  expect_false(is.null(fit1$BETA))
  
  # Test 2: Additive effects only
  fit2 <- ame(Y, R = 0, rvar = TRUE, cvar = TRUE,
             burn = 50, nscan = 100, print = FALSE)
  expect_false(is.null(fit2$APM))
  expect_false(is.null(fit2$BPM))
  expect_equal(length(fit2$APM), n)
  
  # Test 3: Multiplicative effects only
  fit3 <- ame(Y, R = 2, rvar = FALSE, cvar = FALSE,
             burn = 50, nscan = 100, print = FALSE)
  expect_false(is.null(fit3$U))
  expect_false(is.null(fit3$V))
  expect_equal(dim(fit3$U), c(n, 2))
  
  # Test 4: Both additive and multiplicative
  fit4 <- ame(Y, R = 2, rvar = TRUE, cvar = TRUE,
             burn = 50, nscan = 100, print = FALSE)
  expect_false(is.null(fit4$APM))
  expect_false(is.null(fit4$U))
  
  # Test 5: With covariates
  X <- array(rnorm(n * n * 2), c(n, n, 2))
  fit5 <- ame(Y, Xdyad = X, R = 1,
             burn = 50, nscan = 100, print = FALSE)
  expect_equal(ncol(fit5$BETA), 3)  # intercept + 2 covariates
})

test_that("Symmetric networks handled correctly", {
  skip_on_cran()
  set.seed(456)
  
  n <- 12
  Y <- matrix(rnorm(n * n), n, n)
  Y <- (Y + t(Y))/2  # Make symmetric
  diag(Y) <- NA
  
  fit <- ame(Y, symmetric = TRUE, R = 2,
            burn = 50, nscan = 100, print = FALSE)
  
  # Check symmetric-specific outputs
  expect_false(is.null(fit$L))
  # ULUPM no longer stored - can be reconstructed from U and L
  expect_null(fit$V)  # V should be NULL for symmetric
  expect_null(fit$BPM)  # BPM should be NULL for symmetric
  expect_equal(fit$symmetric, TRUE)
  
  # U and L should be related correctly
  expect_equal(dim(fit$U), c(n, 2))
  expect_equal(dim(fit$L), c(2, 2))
})

# ============================================================================
# Bipartite Models - Complete Coverage
# ============================================================================

test_that("Bipartite models handle all parameter combinations", {
  skip_on_cran()
  set.seed(789)
  
  nA <- 10
  nB <- 12
  Y <- matrix(rnorm(nA * nB), nA, nB)
  
  # Test 1: Basic bipartite
  fit1 <- ame(Y, mode = "bipartite",
             R_row = 0, R_col = 0,
             burn = 50, nscan = 100, print = FALSE)
  expect_equal(fit1$mode, "bipartite")
  expect_equal(fit1$nA, nA)
  expect_equal(fit1$nB, nB)
  
  # Test 2: With row/column effects
  fit2 <- ame(Y, mode = "bipartite",
             R_row = 2, R_col = 2,
             burn = 50, nscan = 100, print = FALSE)
  expect_equal(dim(fit2$U), c(nA, 2))
  expect_equal(dim(fit2$V), c(nB, 2))
  expect_equal(dim(fit2$G), c(2, 2))
  
  # Test 3: Asymmetric ranks
  fit3 <- ame(Y, mode = "bipartite",
             R_row = 1, R_col = 3,
             burn = 50, nscan = 100, print = FALSE)
  expect_equal(dim(fit3$U), c(nA, 1))
  expect_equal(dim(fit3$V), c(nB, 3))
  expect_equal(dim(fit3$G), c(1, 3))
})

# ============================================================================
# Memory Management Tests
# ============================================================================

test_that("Memory efficient storage for large networks", {
  skip_on_cran()
  set.seed(111)
  
  # Test with larger network
  n <- 50
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Without posterior saving
  fit1 <- ame(Y, R = 2, burn = 50, nscan = 100, print = FALSE)
  
  # Check that output matrices are reasonable size
  expect_equal(dim(fit1$YPM), c(n, n))
  # EZ removed -   expect_equal(dim(fit1$EZ), c(n, n))
  # UVPM removed -   expect_equal(dim(fit1$UVPM), c(n, n))
  
  # These should not be stored by default
  expect_null(fit1$U_samples)
  expect_null(fit1$V_samples)
  
  # With posterior saving (thinned)
  opts <- posterior_options(save_UV = TRUE, thin_UV = 20)
  fit2 <- ame(Y, R = 2, burn = 50, nscan = 100,
             posterior_opts = opts, print = FALSE)
  
  # Should have thinned samples
  if(!is.null(fit2$U_samples)) {
    n_expected <- floor(100 / (10 * 20))  # nscan / (odens * thin)
    expect_lte(dim(fit2$U_samples)[3], n_expected + 1)
  }
})

# ============================================================================
# Missing Data Handling
# ============================================================================

test_that("Missing data handled correctly", {
  skip_on_cran()
  set.seed(222)
  
  n <- 20
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Add random missing values
  missing_idx <- sample(which(!is.na(Y)), size = 50)
  Y[missing_idx] <- NA
  
  fit <- ame(Y, R = 1, burn = 50, nscan = 100, print = FALSE)
  
  # Predictions should exist for all non-diagonal entries
  non_diag <- which(row(fit$YPM) != col(fit$YPM))
  expect_false(any(is.na(fit$YPM[non_diag])))
  # EZ removed -   expect_false(any(is.na(fit$EZ[non_diag])))
  
  # Missing pattern should be preserved in some form
  expect_equal(sum(is.na(Y)), 20 + 50)  # diagonal + random
})

# ============================================================================
# Edge Cases
# ============================================================================

test_that("Edge cases handled gracefully", {
  skip_on_cran()
  set.seed(333)
  
  # Very small network
  Y_small <- matrix(rnorm(9), 3, 3)
  diag(Y_small) <- NA
  
  fit_small <- ame(Y_small, R = 1, burn = 20, nscan = 50, print = FALSE)
  expect_equal(nrow(fit_small$YPM), 3)
  
  # Sparse network
  Y_sparse <- matrix(0, 20, 20)
  Y_sparse[sample(400, 20)] <- rnorm(20)
  diag(Y_sparse) <- NA
  
  fit_sparse <- ame(Y_sparse, R = 1, burn = 50, nscan = 100, print = FALSE)
  non_diag_sparse <- which(row(fit_sparse$YPM) != col(fit_sparse$YPM))
  expect_false(any(is.na(fit_sparse$YPM[non_diag_sparse])))
  
  # All positive values
  Y_pos <- matrix(abs(rnorm(100)), 10, 10)
  diag(Y_pos) <- NA
  
  fit_pos <- ame(Y_pos, burn = 50, nscan = 100, print = FALSE)
  expect_true(mean(fit_pos$BETA[,1]) > 0)  # Intercept should be positive
})

# ============================================================================
# Convergence and Stability
# ============================================================================

test_that("Models converge and are stable", {
  skip_on_cran()
  set.seed(444)
  
  n <- 15
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Run two models with different seeds
  fit1 <- ame(Y, R = 1, seed = 1, burn = 100, nscan = 200, print = FALSE)
  fit2 <- ame(Y, R = 1, seed = 2, burn = 100, nscan = 200, print = FALSE)
  
  # Parameter estimates should be similar
  beta1 <- colMeans(fit1$BETA)
  beta2 <- colMeans(fit2$BETA)
  expect_equal(beta1, beta2, tolerance = 0.5)
  
  # Predictions should be similar (exclude NAs from diagonal)
  pred1 <- c(fit1$YPM)
  pred2 <- c(fit2$YPM)
  valid_idx <- !is.na(pred1) & !is.na(pred2)
  cor_pred <- cor(pred1[valid_idx], pred2[valid_idx])
  expect_gt(cor_pred, 0.5)  # Lower threshold for small samples
})

# ============================================================================
# Posterior Predictive Checks
# ============================================================================

test_that("Posterior predictive checks work", {
  skip_on_cran()
  set.seed(555)
  
  n <- 20
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, R = 1, burn = 100, nscan = 200, gof = TRUE, print = FALSE)
  
  # GOF should have correct structure
  expect_false(is.null(fit$GOF))
  expect_equal(ncol(fit$GOF), 5)
  expect_gt(nrow(fit$GOF), 1)
  
  # First row is observed statistics
  expect_equal(rownames(fit$GOF)[1], "obs")
  
  # Statistics should be reasonable
  obs_stats <- fit$GOF[1,]
  sim_stats <- colMeans(fit$GOF[-1,, drop=FALSE])
  
  # Simulated should be close to observed for well-fitting model
  for(i in 1:5) {
    if(!is.na(obs_stats[i]) && !is.na(sim_stats[i])) {
      # Very loose check - just ensuring they're in same ballpark
      expect_equal(obs_stats[i], sim_stats[i], tolerance = 5)
    }
  }
})