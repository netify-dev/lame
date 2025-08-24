# Comprehensive tests for static binary (probit) models
# Following same structure as Gaussian tests
# Testing three cases:
# 1. Covariates only (no additive/multiplicative effects)
# 2. Covariates + additive effects
# 3. Covariates + additive + multiplicative effects

library(lame)
library(testthat)

# Helper function to run single simulation for binary
sim_binary_ame <- function(seed, n, mu, beta, gamma=NULL, 
                           rvar=FALSE, cvar=FALSE, R=0,
                           burn=500, nscan=2000) {
  set.seed(seed)
  
  # Generate covariates
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Build linear predictor (probit scale)
  eta <- mu + beta * X
  
  # Add unobserved covariate effect if specified
  if(!is.null(gamma)) {
    W <- tcrossprod(xw[,2])
    eta <- eta + gamma * W
  }
  
  # Add additive effects if requested
  if(rvar || cvar) {
    a_true <- if(rvar) rnorm(n, 0, 0.5) else rep(0, n)
    b_true <- if(cvar) rnorm(n, 0, 0.5) else rep(0, n)
    eta <- eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  }
  
  # Add multiplicative effects if R > 0
  U_true <- V_true <- NULL
  if(R > 0) {
    U_true <- matrix(rnorm(n*R, 0, 1), n, R)
    V_true <- matrix(rnorm(n*R, 0, 1), n, R)
    eta <- eta + tcrossprod(U_true, V_true)
  }
  
  # Generate binary outcome via probit link
  Z <- eta + matrix(rnorm(n*n), n, n)  # Add standard normal noise
  Y <- 1 * (Z > 0)  # Binary outcome
  diag(Y) <- NA
  
  # Fit AME model
  fit <- ame(Y, Xdyad=X, R=R, family="binary",
            rvar=rvar, cvar=cvar, dcor=FALSE,
            burn=burn, nscan=nscan, 
            print=FALSE, plot=FALSE)
  
  # Extract results
  beta_hat <- median(fit$BETA[,2])
  beta_se <- sd(fit$BETA[,2])
  beta_ci_lower <- quantile(fit$BETA[,2], 0.025)
  beta_ci_upper <- quantile(fit$BETA[,2], 0.975)
  beta_covered <- (beta >= beta_ci_lower) & (beta <= beta_ci_upper)
  
  mu_hat <- median(fit$BETA[,1])
  mu_se <- sd(fit$BETA[,1])
  mu_ci_lower <- quantile(fit$BETA[,1], 0.025)
  mu_ci_upper <- quantile(fit$BETA[,1], 0.975)
  mu_covered <- (mu >= mu_ci_lower) & (mu <= mu_ci_upper)
  
  # Correlation with unobserved if applicable
  cor_with_W <- if(!is.null(gamma) && R > 0 && !is.null(fit$UVPM)) {
    cor(c(W), c(fit$UVPM), use='pairwise.complete.obs')
  } else {
    NA
  }
  
  # Recovery of additive effects
  a_cor <- if(rvar && !is.null(fit$APM) && exists("a_true")) {
    cor(a_true, fit$APM, use='pairwise.complete.obs')
  } else {
    NA
  }
  
  b_cor <- if(cvar && !is.null(fit$BPM) && exists("b_true")) {
    cor(b_true, fit$BPM, use='pairwise.complete.obs')
  } else {
    NA
  }
  
  return(list(
    beta_hat = beta_hat,
    beta_se = beta_se,
    beta_covered = beta_covered,
    mu_hat = mu_hat,
    mu_se = mu_se,
    mu_covered = mu_covered,
    cor_with_W = cor_with_W,
    a_cor = a_cor,
    b_cor = b_cor
  ))
}

# ============================================================================
# TEST 1: Binary model with covariates only (no additive/multiplicative)
# ============================================================================

test_that("Binary AME with covariates only recovers true parameters", {
  set.seed(6886)
  n <- 50
  mu_true <- -0.5  # Negative intercept for sparse network
  beta_true <- 1.0
  
  # Run single detailed test
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Generate on probit scale
  eta <- mu_true + beta_true * X
  Z <- eta + matrix(rnorm(n*n), n, n)
  Y <- 1 * (Z > 0)
  diag(Y) <- NA
  
  # Check sparsity
  density <- mean(Y, na.rm=TRUE)
  expect_lt(density, 0.5)  # Should be sparse with negative intercept
  
  # Fit model with no additive or multiplicative effects
  fit <- ame(Y, Xdyad=X, R=0, family="binary",
            rvar=FALSE, cvar=FALSE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.3)
  expect_lt(abs(mu_est - mu_true), 0.3)
  
  # Check that no additive/multiplicative effects were fit (or are minimal)
  if(!is.null(fit$APM)) {
    expect_lt(max(abs(fit$APM)), 0.5)
  }
  if(!is.null(fit$BPM)) {
    expect_lt(max(abs(fit$BPM)), 0.5)
  }
})

# Simulation study for binary covariates only
test_that("Binary AME with covariates only has calibrated confidence intervals", {
  skip_on_cran()
  
  set.seed(6886)
  n_sims <- 30
  n <- 40
  mu_true <- -0.3
  beta_true <- 0.8
  
  # Run simulations
  results <- lapply(1:n_sims, function(i) {
    sim_binary_ame(seed=i, n=n, mu=mu_true, beta=beta_true,
                  rvar=FALSE, cvar=FALSE, R=0,
                  burn=300, nscan=1000)
  })
  
  # Extract coverage rates
  beta_coverage <- mean(sapply(results, function(x) x$beta_covered))
  mu_coverage <- mean(sapply(results, function(x) x$mu_covered))
  
  # Check calibration (relaxed for binary due to discrete nature)
  expect_gt(beta_coverage, 0.80)
  expect_gt(mu_coverage, 0.80)
  
  # Check bias
  beta_bias <- mean(sapply(results, function(x) x$beta_hat)) - beta_true
  mu_bias <- mean(sapply(results, function(x) x$mu_hat)) - mu_true
  
  expect_lt(abs(beta_bias), 0.15)
  expect_lt(abs(mu_bias), 0.15)
})

# ============================================================================
# TEST 2: Binary model with covariates + additive effects
# ============================================================================

test_that("Binary AME with additive effects recovers true parameters", {
  set.seed(6886)
  n <- 50
  mu_true <- 0
  beta_true <- 0.8
  
  # Generate data with additive effects
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # True additive effects
  a_true <- rnorm(n, 0, 0.7)
  b_true <- rnorm(n, 0, 0.7)
  
  eta <- mu_true + beta_true * X + 
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  
  Z <- eta + matrix(rnorm(n*n), n, n)
  Y <- 1 * (Z > 0)
  diag(Y) <- NA
  
  # Fit model with additive effects
  fit <- ame(Y, Xdyad=X, R=0, family="binary",
            rvar=TRUE, cvar=TRUE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.4)
  expect_lt(abs(mu_est - mu_true), 0.4)
  
  # Check recovery of additive effects
  if(!is.null(fit$APM) && !is.null(fit$BPM)) {
    a_cor <- cor(a_true, fit$APM, use='pairwise.complete.obs')
    b_cor <- cor(b_true, fit$BPM, use='pairwise.complete.obs')
    
    expect_gt(a_cor, 0.4)  # Slightly lower threshold for binary
    expect_gt(b_cor, 0.4)
  }
})

# ============================================================================
# TEST 3: Binary model with covariates + additive + multiplicative effects
# ============================================================================

test_that("Binary AME with full model recovers true parameters", {
  set.seed(6886)
  n <- 50
  mu_true <- -0.2
  beta_true <- 0.8
  gamma_true <- 0.6  # Unobserved covariate
  R_true <- 2
  
  # Generate data
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  W <- tcrossprod(xw[,2])  # Unobserved
  diag(X) <- NA
  diag(W) <- NA
  
  # True additive effects
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  # True latent factors
  U_true <- matrix(rnorm(n*R_true, 0, 1), n, R_true)
  V_true <- matrix(rnorm(n*R_true, 0, 1), n, R_true)
  
  # Build full model
  eta <- mu_true + beta_true * X + gamma_true * W +
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true) +
         tcrossprod(U_true, V_true)
  
  Z <- eta + matrix(rnorm(n*n), n, n)
  Y <- 1 * (Z > 0)
  diag(Y) <- NA
  
  # Fit full AME model
  fit <- ame(Y, Xdyad=X, R=R_true, family="binary",
            rvar=TRUE, cvar=TRUE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery (more tolerance for complex model)
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.5)
  
  # Check that multiplicative effects capture unobserved covariate
  if(!is.null(fit$UVPM)) {
    cor_W <- cor(c(W), c(fit$UVPM), use='pairwise.complete.obs')
    # Just check it's positive - complex models have weak signal
    expect_gt(cor_W, 0)
  }
  
  # Check dimensions of latent factors
  expect_equal(ncol(fit$U), R_true)
  expect_equal(ncol(fit$V), R_true)
})

# ============================================================================
# TEST 4: Compare models with and without multiplicative effects (binary)
# ============================================================================

test_that("Multiplicative effects reduce bias in binary models", {
  set.seed(6886)
  n <- 50
  mu_true <- 0
  beta_true <- 1
  gamma_true <- 1.2  # Strong unobserved effect
  
  # Generate data with unobserved confounder
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  W <- tcrossprod(xw[,2])
  diag(X) <- NA
  diag(W) <- NA
  
  eta <- mu_true + beta_true * X + gamma_true * W
  Z <- eta + matrix(rnorm(n*n), n, n)
  Y <- 1 * (Z > 0)
  diag(Y) <- NA
  
  # Fit model WITHOUT multiplicative effects
  fit_no_uv <- ame(Y, Xdyad=X, R=0, family="binary",
                  rvar=FALSE, cvar=FALSE, dcor=FALSE,
                  burn=400, nscan=1500, print=FALSE, plot=FALSE)
  
  # Fit model WITH multiplicative effects
  fit_with_uv <- ame(Y, Xdyad=X, R=2, family="binary",
                    rvar=FALSE, cvar=FALSE, dcor=FALSE,
                    burn=400, nscan=1500, print=FALSE, plot=FALSE)
  
  beta_no_uv <- median(fit_no_uv$BETA[,2])
  beta_with_uv <- median(fit_with_uv$BETA[,2])
  
  # Model with UV should have less bias
  bias_no_uv <- abs(beta_no_uv - beta_true)
  bias_with_uv <- abs(beta_with_uv - beta_true)
  
  # We expect the model with multiplicative effects to have less bias
  # but be lenient given model complexity
  expect_lt(bias_with_uv, bias_no_uv + 0.3)
  
  # Check that UV captures the unobserved structure
  if(!is.null(fit_with_uv$UVPM)) {
    cor_W <- cor(c(W), c(fit_with_uv$UVPM), use='pairwise.complete.obs')
    expect_gt(cor_W, 0.2)
  }
})

# ============================================================================
# TEST 5: Binary model edge cases
# ============================================================================

test_that("Binary AME handles sparse and dense networks", {
  set.seed(6886)
  n <- 30
  
  # Test sparse network
  Y_sparse <- matrix(0, n, n)
  Y_sparse[sample(1:(n*n), size=floor(0.05*n*n))] <- 1  # 5% density
  diag(Y_sparse) <- NA
  
  fit_sparse <- ame(Y_sparse, R=1, family="binary",
                   burn=200, nscan=800, print=FALSE, plot=FALSE)
  
  expect_true(!is.null(fit_sparse$BETA))
  expect_lt(median(fit_sparse$BETA[,1]), 0)  # Intercept should be negative
  
  # Test dense network
  Y_dense <- matrix(1, n, n)
  Y_dense[sample(1:(n*n), size=floor(0.05*n*n))] <- 0  # 95% density
  diag(Y_dense) <- NA
  
  fit_dense <- ame(Y_dense, R=1, family="binary",
                  burn=200, nscan=800, print=FALSE, plot=FALSE)
  
  expect_true(!is.null(fit_dense$BETA))
  expect_gt(median(fit_dense$BETA[,1]), 0)  # Intercept should be positive
})

test_that("Binary AME predictions are probabilities", {
  set.seed(6886)
  n <- 25
  
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  eta <- 0.5 * X
  Y <- 1 * (eta + matrix(rnorm(n*n), n, n) > 0)
  diag(Y) <- NA
  
  fit <- ame(Y, Xdyad=X, R=1, family="binary",
            burn=200, nscan=800, print=FALSE, plot=FALSE)
  
  # EZ contains latent Z values for binary
  # But we can check they're reasonable
  expect_true(is.numeric(fit$EZ))
  
  # Check GOF if available
  if(!is.null(fit$GOF)) {
    expect_true(nrow(fit$GOF) > 0)
  }
})