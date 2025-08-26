# Tests for posterior utility functions
library(lame)

test_that("posterior_options creates correct structure", {
  # Test default options
  opts1 <- posterior_options()
  expect_equal(opts1$save_UV, FALSE)
  expect_equal(opts1$save_ab, FALSE)
  expect_equal(opts1$thin_UV, 10)
  expect_equal(opts1$thin_ab, 10)
  
  # Test custom options
  opts2 <- posterior_options(save_UV = TRUE, save_ab = TRUE, 
                            thin_UV = 5, thin_ab = 20)
  expect_equal(opts2$save_UV, TRUE)
  expect_equal(opts2$save_ab, TRUE)
  expect_equal(opts2$thin_UV, 5)
  expect_equal(opts2$thin_ab, 20)
})

test_that("simulate_posterior works for beta coefficients", {
  set.seed(123)
  n <- 20
  
  # Generate simple network
  X <- array(rnorm(n * n * 2), c(n, n, 2))
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad = X, R = 0, 
            burn = 100, nscan = 500, print = FALSE)
  
  # Simulate posterior for beta
  beta_post <- simulate_posterior(fit, "beta", n_samples = 100)
  
  # Check dimensions (may have fewer samples than requested)
  expect_lte(nrow(beta_post), 100)
  expect_gt(nrow(beta_post), 0)
  expect_equal(ncol(beta_post), ncol(fit$BETA))
  
  # Check values are numeric
  expect_true(all(is.finite(beta_post)))
  
  # Check mean is close to MCMC mean
  beta_mcmc_mean <- colMeans(fit$BETA)
  beta_sim_mean <- colMeans(beta_post)
  expect_equal(beta_sim_mean, beta_mcmc_mean, tolerance = 0.5)
})

test_that("simulate_posterior works for additive effects", {
  set.seed(456)
  n <- 20
  
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model with additive effects
  fit <- ame(Y, R = 0, rvar = TRUE, cvar = TRUE,
            burn = 100, nscan = 500, print = FALSE)
  
  # Simulate posterior for additive effects
  ab_post <- simulate_posterior(fit, "ab", n_samples = 50)
  
  # Check structure
  expect_true(is.list(ab_post))
  expect_equal(names(ab_post), c("a", "b"))
  
  # Check dimensions
  expect_equal(dim(ab_post$a), c(n, 50))
  expect_equal(dim(ab_post$b), c(n, 50))
  
  # Check values are numeric
  expect_true(all(is.finite(ab_post$a)))
  expect_true(all(is.finite(ab_post$b)))
  
  # Check means are reasonable
  a_mean <- rowMeans(ab_post$a)
  expect_equal(mean(a_mean), 0, tolerance = 1)
})

test_that("simulate_posterior works for multiplicative effects", {
  set.seed(789)
  n <- 20
  R <- 2
  
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model with multiplicative effects
  fit <- ame(Y, R = R, burn = 100, nscan = 500, print = FALSE)
  
  # Simulate posterior for UV
  UV_post <- simulate_posterior(fit, "UV", n_samples = 10)
  
  # Check dimensions
  expect_lte(dim(UV_post)[3], 10)
  expect_gt(dim(UV_post)[3], 0)
  expect_equal(dim(UV_post)[1:2], c(n, n))
  
  # Check values are numeric
  expect_true(all(is.finite(UV_post)))
  
  # Check that UVPM exists
  # UVPM removed -   expect_false(is.null(fit$UVPM))
})

test_that("posterior_quantiles computes correct intervals", {
  set.seed(321)
  n <- 20
  
  X <- array(rnorm(n * n * 2), c(n, n, 2))
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad = X, R = 0,
            burn = 100, nscan = 500, print = FALSE)
  
  # Get quantiles for beta
  beta_quants <- posterior_quantiles(fit, "beta", probs = c(0.025, 0.5, 0.975))
  
  # Check dimensions
  expect_equal(nrow(beta_quants), 3)  # 3 quantiles
  expect_equal(ncol(beta_quants), ncol(fit$BETA))  # number of parameters
  
  # Check ordering (lower < median < upper)
  expect_true(all(beta_quants[1,] <= beta_quants[2,]))
  expect_true(all(beta_quants[2,] <= beta_quants[3,]))
  
  # Check median is close to mean for approximately normal posteriors
  beta_means <- colMeans(fit$BETA)
  beta_medians <- beta_quants[2,]
  expect_equal(beta_medians, beta_means, tolerance = 0.5)
})

test_that("posterior samples are saved when requested", {
  set.seed(654)
  n <- 15
  R <- 2
  
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model with posterior saving
  opts <- posterior_options(save_UV = TRUE, save_ab = TRUE, 
                           thin_UV = 5, thin_ab = 5)
  
  fit <- ame(Y, R = R, rvar = TRUE, cvar = TRUE,
            burn = 50, nscan = 200, odens = 5,
            posterior_opts = opts, print = FALSE)
  
  # Check that samples were saved
  if(!is.null(opts) && opts$save_UV) {
    expect_false(is.null(fit$U_samples))
    expect_false(is.null(fit$V_samples))
    
    # Check dimensions
    n_expected <- floor(200 / (5 * 5))  # nscan / (odens * thin)
    expect_equal(dim(fit$U_samples)[1], n)
    expect_equal(dim(fit$U_samples)[2], R)
    expect_lte(dim(fit$U_samples)[3], n_expected + 1)
  }
  
  if(!is.null(opts) && opts$save_ab) {
    expect_false(is.null(fit$a_samples))
    expect_false(is.null(fit$b_samples))
    
    # Check dimensions
    expect_equal(nrow(fit$a_samples), n)
    expect_equal(nrow(fit$b_samples), n)
  }
})

test_that("simulate_posterior works with saved samples", {
  set.seed(987)
  n <- 15
  R <- 2
  
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  # Fit model with posterior saving
  opts <- posterior_options(save_UV = TRUE, thin_UV = 10)
  
  fit <- ame(Y, R = R, burn = 50, nscan = 200, odens = 5,
            posterior_opts = opts, print = FALSE)
  
  # Simulate from saved samples
  if(!is.null(fit$U_samples) && !is.null(fit$V_samples)) {
    UV_post <- simulate_posterior(fit, "UV", n_samples = 5)
    
    # Should use saved samples
    expect_equal(dim(UV_post)[1:2], c(n, n))
    expect_lte(dim(UV_post)[3], 5)
  }
})

test_that("posterior functions work for bipartite networks", {
  set.seed(1234)
  nA <- 10
  nB <- 12
  
  Y <- matrix(rnorm(nA * nB), nA, nB)
  
  # Fit bipartite model
  fit <- ame(Y, mode = "bipartite", R_row = 1, R_col = 1,
            burn = 50, nscan = 200, print = FALSE)
  
  # Check we can simulate posteriors
  beta_post <- simulate_posterior(fit, "beta", n_samples = 20)
  expect_lte(nrow(beta_post), 20)  # May have fewer if not enough MCMC samples
  expect_gt(nrow(beta_post), 0)    # But should have some
  
  # Check quantiles work
  beta_quants <- posterior_quantiles(fit, "beta")
  expect_equal(nrow(beta_quants), 3)  # Default 3 quantiles
})

test_that("simulate_Y_posterior generates reasonable predictions", {
  set.seed(5678)
  n <- 20
  
  # Simple model
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, R = 0, burn = 100, nscan = 300, print = FALSE)
  
  # Simulate posterior predictive
  Y_post <- simulate_posterior(fit, "Y", n_samples = 10)
  
  # Check dimensions
  expect_equal(dim(Y_post), c(n, n, 10))
  
  # Check missing values are preserved
  for(i in 1:10) {
    expect_true(all(is.na(diag(Y_post[,,i]))))
  }
  
  # Check values are reasonable (similar scale to original)
  Y_mean <- mean(Y, na.rm = TRUE)
  Y_sd <- sd(Y, na.rm = TRUE)
  Y_post_mean <- mean(Y_post, na.rm = TRUE)
  Y_post_sd <- sd(Y_post, na.rm = TRUE)
  
  expect_equal(Y_post_mean, Y_mean, tolerance = 1)
  expect_equal(Y_post_sd, Y_sd, tolerance = 1)
})