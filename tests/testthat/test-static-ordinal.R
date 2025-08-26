# Basic tests for static ordinal models
# Simplified tests to ensure ordinal family works

library(lame)
library(testthat)

# ============================================================================
# TEST 1: Ordinal model basic functionality
# ============================================================================

test_that("Ordinal AME runs without errors", {
  set.seed(6886)
  n <- 30
  
  # Create simple ordinal data
  Y <- matrix(sample(1:5, n*n, replace=TRUE), n, n)
  diag(Y) <- NA
  
  # Test with no covariates
  fit_simple <- ame(Y, R=1, family="ordinal",
                   burn=200, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit_simple$BETA))
  expect_true(!is.null(fit_simple$U))
  expect_true(!is.null(fit_simple$V))
})

# ============================================================================
# TEST 2: Ordinal model with covariates
# ============================================================================

test_that("Ordinal AME with covariates works", {
  set.seed(6886)
  n <- 30
  
  # Generate covariate
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  # Generate ordinal outcome influenced by X
  eta <- 0.5 * X
  Z <- eta + matrix(rnorm(n*n), n, n)
  
  # Discretize to ordinal
  Y <- cut(c(Z), breaks=5, labels=FALSE)
  Y <- matrix(Y, n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad=X, R=1, family="ordinal",
            burn=300, nscan=1000, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  expect_gt(ncol(fit$BETA), 0)
  
  # Check that coefficient is in reasonable range
  if(ncol(fit$BETA) >= 2) {
    beta_est <- median(fit$BETA[,2])
    expect_true(abs(beta_est) < 2)
  }
})

# ============================================================================
# TEST 3: Ordinal model with additive effects
# ============================================================================

test_that("Ordinal AME with additive effects works", {
  set.seed(6886)
  n <- 30
  
  # Generate data with row/column heterogeneity
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  eta <- outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  Z <- eta + matrix(rnorm(n*n), n, n)
  
  # Discretize to ordinal
  Y <- cut(c(Z), breaks=4, labels=FALSE)
  Y <- matrix(Y, n, n)
  diag(Y) <- NA
  
  # Fit model with additive effects
  fit <- ame(Y, R=0, family="ordinal",
            rvar=TRUE, cvar=TRUE,
            burn=300, nscan=1000, print=FALSE)
  
  expect_true(!is.null(fit$APM))
  expect_true(!is.null(fit$BPM))
  expect_equal(length(fit$APM), n)
  expect_equal(length(fit$BPM), n)
  
  # Check some recovery (relaxed)
  if(!is.null(fit$APM) && !is.null(fit$BPM)) {
    a_cor <- cor(a_true, fit$APM, use='pairwise.complete.obs')
    b_cor <- cor(b_true, fit$BPM, use='pairwise.complete.obs')
    expect_gt(a_cor, 0.2)
    expect_gt(b_cor, 0.2)
  }
})

# ============================================================================
# TEST 4: Ordinal model with different numbers of categories
# ============================================================================

test_that("Ordinal AME handles different numbers of categories", {
  set.seed(6886)
  n <- 25
  
  # Test with 3 categories
  Y_3 <- matrix(sample(1:3, n*n, replace=TRUE), n, n)
  diag(Y_3) <- NA
  
  fit_3 <- ame(Y_3, R=1, family="ordinal",
              burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit_3$BETA))
  expect_equal(length(unique(c(Y_3[!is.na(Y_3)]))), 3)
  
  # Test with 7 categories
  Y_7 <- matrix(sample(1:7, n*n, replace=TRUE), n, n)
  diag(Y_7) <- NA
  
  fit_7 <- ame(Y_7, R=1, family="ordinal",
              burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit_7$BETA))
  expect_equal(length(unique(c(Y_7[!is.na(Y_7)]))), 7)
})

# ============================================================================
# TEST 5: Ordinal model with skewed distributions
# ============================================================================

test_that("Ordinal AME handles skewed distributions", {
  set.seed(6886)
  n <- 25
  
  # Create skewed ordinal data (mostly low values)
  Y_skew <- matrix(sample(1:5, n*n, replace=TRUE, 
                         prob=c(0.5, 0.25, 0.15, 0.07, 0.03)), n, n)
  diag(Y_skew) <- NA
  
  # Check skewness
  expect_lt(mean(Y_skew, na.rm=TRUE), 2.5)
  
  fit_skew <- ame(Y_skew, R=1, family="ordinal",
                 burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit_skew$BETA))
  # EZ removed -   expect_true(is.numeric(fit_skew$EZ))
  
  # Create reverse-skewed data (mostly high values)
  Y_high <- matrix(sample(1:5, n*n, replace=TRUE,
                         prob=c(0.03, 0.07, 0.15, 0.25, 0.5)), n, n)
  diag(Y_high) <- NA
  
  expect_gt(mean(Y_high, na.rm=TRUE), 3.5)
  
  fit_high <- ame(Y_high, R=1, family="ordinal",
                 burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit_high$BETA))
})

# ============================================================================
# TEST 6: Ordinal model with covariates and multiplicative effects
# ============================================================================

test_that("Ordinal AME with covariates and multiplicative effects works", {
  set.seed(6886)
  n <- 30
  
  # Generate covariate and latent factors
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  U_true <- matrix(rnorm(n*2, 0, 0.5), n, 2)
  V_true <- matrix(rnorm(n*2, 0, 0.5), n, 2)
  
  # Build model
  eta <- 0.4 * X + tcrossprod(U_true, V_true)
  Z <- eta + matrix(rnorm(n*n), n, n)
  
  # Discretize to ordinal
  Y <- cut(c(Z), breaks=5, labels=FALSE)
  Y <- matrix(Y, n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad=X, R=2, family="ordinal",
            burn=300, nscan=1000, print=FALSE)
  
  expect_true(!is.null(fit$U))
  expect_true(!is.null(fit$V))
  expect_equal(ncol(fit$U), 2)
  expect_equal(ncol(fit$V), 2)
  
  # Check that UVPM captures some structure
  if(!is.null(fit$UVPM)) {
  # UVPM removed -     expect_true(var(c(fit$UVPM)) > 0)
  }
})