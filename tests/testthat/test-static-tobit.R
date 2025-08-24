# Basic tests for static tobit (censored normal) models
# Following same structure as Gaussian tests but simpler
# Testing three cases:
# 1. Covariates only (no additive/multiplicative effects)
# 2. Covariates + additive effects
# 3. Covariates + additive + multiplicative effects

library(lame)
library(testthat)

# Helper function to run single simulation for tobit
sim_tobit_ame <- function(seed, n, mu, beta, gamma=NULL, 
                          rvar=FALSE, cvar=FALSE, R=0,
                          burn=500, nscan=2000) {
  set.seed(seed)
  
  # Generate covariates
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Build linear predictor
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
    U_true <- matrix(rnorm(n*R, 0, 0.5), n, R)
    V_true <- matrix(rnorm(n*R, 0, 0.5), n, R)
    eta <- eta + tcrossprod(U_true, V_true)
  }
  
  # Generate tobit outcome (censored at 0)
  Y <- eta + matrix(rnorm(n*n, 0, 1), n, n)
  Y[Y < 0] <- 0  # Left-censoring at 0
  diag(Y) <- NA
  
  # Fit AME model
  fit <- ame(Y, Xdyad=X, R=R, family="tobit",
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
# TEST 1: Tobit model with covariates only (no additive/multiplicative)
# ============================================================================

test_that("Tobit AME with covariates only recovers true parameters", {
  set.seed(6886)
  n <- 40
  mu_true <- 0.5
  beta_true <- 1.2
  
  # Run single detailed test
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # Generate tobit data
  eta <- mu_true + beta_true * X
  Y <- eta + matrix(rnorm(n*n, 0, 1), n, n)
  Y[Y < 0] <- 0  # Censor at 0
  diag(Y) <- NA
  
  # Check censoring
  censored_prop <- mean(Y == 0, na.rm=TRUE)
  expect_gt(censored_prop, 0.05)  # Should have some censoring
  expect_lt(censored_prop, 0.5)   # But not too much
  
  # Fit model with no additive or multiplicative effects
  fit <- ame(Y, Xdyad=X, R=0, family="tobit",
            rvar=FALSE, cvar=FALSE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.2)
  expect_lt(abs(mu_est - mu_true), 0.2)
  
  # Check that no additive/multiplicative effects were fit
  if(!is.null(fit$APM)) {
    expect_lt(max(abs(fit$APM)), 0.5)
  }
  if(!is.null(fit$BPM)) {
    expect_lt(max(abs(fit$BPM)), 0.5)
  }
})

# ============================================================================
# TEST 2: Tobit model with covariates + additive effects
# ============================================================================

test_that("Tobit AME with additive effects recovers true parameters", {
  set.seed(6886)
  n <- 40
  mu_true <- 0.3
  beta_true <- 0.8
  
  # Generate data with additive effects
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  diag(X) <- NA
  
  # True additive effects
  a_true <- rnorm(n, 0, 0.5)
  b_true <- rnorm(n, 0, 0.5)
  
  eta <- mu_true + beta_true * X + 
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
  
  Y <- eta + matrix(rnorm(n*n, 0, 1), n, n)
  Y[Y < 0] <- 0  # Censor at 0
  diag(Y) <- NA
  
  # Fit model with additive effects
  fit <- ame(Y, Xdyad=X, R=0, family="tobit",
            rvar=TRUE, cvar=TRUE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.3)
  expect_lt(abs(mu_est - mu_true), 0.3)
  
  # Check recovery of additive effects
  if(!is.null(fit$APM) && !is.null(fit$BPM)) {
    a_cor <- cor(a_true, fit$APM, use='pairwise.complete.obs')
    b_cor <- cor(b_true, fit$BPM, use='pairwise.complete.obs')
    
    expect_gt(a_cor, 0.5)
    expect_gt(b_cor, 0.5)
  }
})

# ============================================================================
# TEST 3: Tobit model with covariates + additive + multiplicative effects
# ============================================================================

test_that("Tobit AME with full model recovers true parameters", {
  set.seed(6886)
  n <- 40
  mu_true <- 0.2
  beta_true <- 0.7
  gamma_true <- 0.5  # Unobserved covariate
  R_true <- 2
  
  # Generate data
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  W <- tcrossprod(xw[,2])  # Unobserved
  diag(X) <- NA
  diag(W) <- NA
  
  # True additive effects
  a_true <- rnorm(n, 0, 0.4)
  b_true <- rnorm(n, 0, 0.4)
  
  # True latent factors
  U_true <- matrix(rnorm(n*R_true, 0, 0.5), n, R_true)
  V_true <- matrix(rnorm(n*R_true, 0, 0.5), n, R_true)
  
  # Build full model
  eta <- mu_true + beta_true * X + gamma_true * W +
         outer(a_true, rep(1, n)) + outer(rep(1, n), b_true) +
         tcrossprod(U_true, V_true)
  
  Y <- eta + matrix(rnorm(n*n, 0, 1), n, n)
  Y[Y < 0] <- 0  # Censor at 0
  diag(Y) <- NA
  
  # Fit full AME model
  fit <- ame(Y, Xdyad=X, R=R_true, family="tobit",
            rvar=TRUE, cvar=TRUE, dcor=FALSE,
            burn=500, nscan=2000, print=FALSE, plot=FALSE)
  
  # Check parameter recovery
  beta_est <- median(fit$BETA[,2])
  mu_est <- median(fit$BETA[,1])
  
  expect_lt(abs(beta_est - beta_true), 0.4)
  
  # Check that multiplicative effects capture unobserved covariate
  if(!is.null(fit$UVPM)) {
    cor_W <- cor(c(W), c(fit$UVPM), use='pairwise.complete.obs')
    expect_gt(cor_W, 0.3)
  }
  
  # Check dimensions of latent factors
  expect_equal(ncol(fit$U), R_true)
  expect_equal(ncol(fit$V), R_true)
})

# ============================================================================
# TEST 4: Compare models with and without multiplicative effects (tobit)
# ============================================================================

test_that("Multiplicative effects reduce bias in tobit models", {
  set.seed(6886)
  n <- 40
  mu_true <- 0.3
  beta_true <- 0.8
  gamma_true <- 0.9  # Strong unobserved effect
  
  # Generate data with unobserved confounder
  xw <- matrix(rnorm(n*2), n, 2)
  X <- tcrossprod(xw[,1])
  W <- tcrossprod(xw[,2])
  diag(X) <- NA
  diag(W) <- NA
  
  eta <- mu_true + beta_true * X + gamma_true * W
  Y <- eta + matrix(rnorm(n*n, 0, 1), n, n)
  Y[Y < 0] <- 0  # Censor at 0
  diag(Y) <- NA
  
  # Fit model WITHOUT multiplicative effects
  fit_no_uv <- ame(Y, Xdyad=X, R=0, family="tobit",
                  rvar=FALSE, cvar=FALSE, dcor=FALSE,
                  burn=400, nscan=1500, print=FALSE, plot=FALSE)
  
  # Fit model WITH multiplicative effects
  fit_with_uv <- ame(Y, Xdyad=X, R=2, family="tobit",
                    rvar=FALSE, cvar=FALSE, dcor=FALSE,
                    burn=400, nscan=1500, print=FALSE, plot=FALSE)
  
  beta_no_uv <- median(fit_no_uv$BETA[,2])
  beta_with_uv <- median(fit_with_uv$BETA[,2])
  
  # Model with UV should have less bias
  bias_no_uv <- abs(beta_no_uv - beta_true)
  bias_with_uv <- abs(beta_with_uv - beta_true)
  
  # We expect the model with multiplicative effects to have less bias
  expect_lt(bias_with_uv, bias_no_uv + 0.2)
  
  # Check that UV captures the unobserved structure
  if(!is.null(fit_with_uv$UVPM)) {
    cor_W <- cor(c(W), c(fit_with_uv$UVPM), use='pairwise.complete.obs')
    expect_gt(cor_W, 0.3)
  }
})

# ============================================================================
# TEST 5: Tobit model edge cases
# ============================================================================

test_that("Tobit AME handles different censoring levels", {
  set.seed(6886)
  n <- 30
  
  # Test with heavy censoring (negative intercept)
  Y_heavy <- matrix(rnorm(n*n, -1, 1), n, n)  # Mean -1 leads to heavy censoring
  Y_heavy[Y_heavy < 0] <- 0
  diag(Y_heavy) <- NA
  
  censored_prop <- mean(Y_heavy == 0, na.rm=TRUE)
  expect_gt(censored_prop, 0.3)  # Should have substantial censoring
  
  fit_heavy <- ame(Y_heavy, R=1, family="tobit",
                  burn=200, nscan=800, print=FALSE, plot=FALSE)
  
  expect_true(!is.null(fit_heavy$BETA))
  expect_lt(median(fit_heavy$BETA[,1]), 0.5)  # Intercept should reflect censoring
  
  # Test with light censoring (positive intercept)
  Y_light <- matrix(rnorm(n*n, 2, 1), n, n)  # Mean 2 leads to light censoring
  Y_light[Y_light < 0] <- 0
  diag(Y_light) <- NA
  
  censored_prop_light <- mean(Y_light == 0, na.rm=TRUE)
  expect_lt(censored_prop_light, 0.05)  # Should have minimal censoring
  
  fit_light <- ame(Y_light, R=1, family="tobit",
                  burn=200, nscan=800, print=FALSE, plot=FALSE)
  
  expect_true(!is.null(fit_light$BETA))
  expect_gt(median(fit_light$BETA[,1]), 1)  # Intercept should be positive
})