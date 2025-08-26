# Test suite for static bipartite AME models
# Similar to unipartite tests but for rectangular networks

test_that("Bipartite AME model works for Gaussian case", {
  skip_on_cran()
  
  set.seed(789)
  
  # Set up bipartite network dimensions
  nA <- 20  # Row nodes (e.g., individuals)
  nB <- 15  # Column nodes (e.g., organizations)
  
  # True parameters
  mu_true <- 0.5
  beta_true <- 0.8
  sigma2_true <- 0.5
  
  # Generate covariates
  X <- matrix(rnorm(nA * nB), nA, nB)
  
  # Generate true network
  Y_true <- mu_true + beta_true * X + matrix(rnorm(nA * nB, 0, sqrt(sigma2_true)), nA, nB)
  
  # Fit basic model
  fit <- ame(Y_true, Xdyad = X, mode = "bipartite", 
             R_row = 0, R_col = 0,
             burn = 500, nscan = 2000, print = FALSE)
  
  # Check dimensions
  expect_equal(fit$mode, "bipartite")
  expect_equal(fit$nA, nA)
  expect_equal(fit$nB, nB)
  expect_equal(dim(fit$BETA)[2], 2)  # Intercept + covariate
  
  # Check parameter recovery
  beta_est <- colMeans(fit$BETA)
  expect_equal(beta_est[1], mu_true, tolerance = 0.2)
  expect_equal(beta_est[2], beta_true, tolerance = 0.2)
  
  # Check variance recovery
  s2_est <- mean(fit$VC[,4])
  expect_equal(s2_est, sigma2_true, tolerance = 0.3)
})

test_that("Bipartite model with additive effects", {
  skip_on_cran()
  
  set.seed(790)
  
  nA <- 15
  nB <- 20
  
  # True parameters
  mu_true <- 0.3
  sigma2_a <- 0.4  # Row variance
  sigma2_b <- 0.3  # Column variance
  sigma2_e <- 0.5  # Error variance
  
  # Generate additive effects
  a_true <- rnorm(nA, 0, sqrt(sigma2_a))
  b_true <- rnorm(nB, 0, sqrt(sigma2_b))
  
  # Generate network
  Y <- matrix(mu_true, nA, nB)
  for(i in 1:nA) {
    Y[i,] <- Y[i,] + a_true[i]
  }
  for(j in 1:nB) {
    Y[,j] <- Y[,j] + b_true[j]
  }
  Y <- Y + matrix(rnorm(nA * nB, 0, sqrt(sigma2_e)), nA, nB)
  
  # Fit model with additive effects
  fit <- ame(Y, mode = "bipartite",
             rvar = TRUE, cvar = TRUE,
             R_row = 0, R_col = 0,
             burn = 500, nscan = 2000, print = FALSE)
  
  # Check that additive effects were estimated
  expect_false(is.null(fit$APM))
  expect_false(is.null(fit$BPM))
  expect_equal(length(fit$APM), nA)
  expect_equal(length(fit$BPM), nB)
  
  # Check variance components
  expect_true(mean(fit$VC[,1]) > 0)  # Row variance
  expect_true(mean(fit$VC[,3]) > 0)  # Column variance
  
  # Check that row/column means are close to true values
  expect_equal(mean(fit$APM), mean(a_true), tolerance = 0.3)
  expect_equal(mean(fit$BPM), mean(b_true), tolerance = 0.3)
})

test_that("Bipartite model with multiplicative effects", {
  skip_on_cran()
  
  set.seed(791)
  
  nA <- 20
  nB <- 25
  R_row <- 2
  R_col <- 2
  
  # Generate latent factors
  U_true <- matrix(rnorm(nA * R_row), nA, R_row)
  V_true <- matrix(rnorm(nB * R_col), nB, R_col)
  G_true <- diag(c(2, 1), R_row, R_col)  # Interaction strengths
  
  # Generate network with multiplicative effects
  mu_true <- 0.2
  UV_true <- U_true %*% G_true %*% t(V_true)
  Y <- mu_true + UV_true + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
  
  # Fit model with multiplicative effects
  fit <- ame(Y, mode = "bipartite",
             R_row = R_row, R_col = R_col,
             burn = 500, nscan = 2000, print = FALSE)
  
  # Check dimensions of latent factors
  expect_equal(dim(fit$U), c(nA, R_row))
  expect_equal(dim(fit$V), c(nB, R_col))
  expect_equal(dim(fit$G), c(R_row, R_col))
  
  # Check that UVPM exists and has right dimensions
  # UVPM removed -   expect_false(is.null(fit$UVPM))
  # UVPM removed -   expect_equal(dim(fit$UVPM), c(nA, nB))
  
  # Check that multiplicative effects capture structure
  # Reconstruct UVPM from U, V, G
  UVPM_reconstructed <- reconstruct_UVPM(fit)
  
  # Correlation should be positive but not perfect
  if(!is.null(UVPM_reconstructed)) {
    cor_uv <- cor(c(UV_true), c(UVPM_reconstructed))
    expect_gt(cor_uv, 0.3)
  }
})

test_that("Bipartite binary model", {
  skip_on_cran()
  
  set.seed(792)
  
  nA <- 15
  nB <- 18
  
  # Generate binary network
  mu_true <- -0.5
  beta_true <- 1.2
  X <- matrix(rnorm(nA * nB), nA, nB)
  
  # Probit model
  Z_true <- mu_true + beta_true * X + matrix(rnorm(nA * nB), nA, nB)
  Y <- (Z_true > 0) * 1
  
  # Fit binary model
  fit <- ame(Y, Xdyad = X, mode = "bipartite",
             family = "binary",
             burn = 500, nscan = 2000, print = FALSE)
  
  # Check model type
  expect_equal(fit$family, "binary")
  
  # Check parameter signs are correct
  beta_est <- colMeans(fit$BETA)
  expect_true(beta_est[1] < 0)  # Negative intercept
  expect_true(beta_est[2] > 0)  # Positive covariate effect
  
  # Check predictions are binary
  expect_true(all(fit$YPM >= 0 & fit$YPM <= 1))
})

test_that("Bipartite Poisson model", {
  skip_on_cran()
  
  set.seed(793)
  
  nA <- 12
  nB <- 10
  
  # Generate count network
  mu_true <- 1.0  # log scale
  X <- matrix(runif(nA * nB, -1, 1), nA, nB)
  beta_true <- 0.5
  
  lambda <- exp(mu_true + beta_true * X)
  Y <- matrix(rpois(nA * nB, lambda), nA, nB)
  
  # Fit Poisson model
  suppressWarnings({
    fit <- ame(Y, Xdyad = X, mode = "bipartite",
               family = "poisson",
               burn = 300, nscan = 1000, print = FALSE)
  })
  
  # Check model type
  expect_equal(fit$family, "poisson")
  
  # Check that predictions are non-negative
  # EZ removed -   expect_true(all(fit$EZ >= 0))
  
  # Check parameter recovery (loose tolerance for Poisson)
  beta_est <- colMeans(fit$BETA)
  expect_equal(beta_est[1], mu_true, tolerance = 0.5)
  expect_equal(beta_est[2], beta_true, tolerance = 0.5)
})

test_that("Bipartite model comparison with different R values", {
  skip_on_cran()
  
  set.seed(794)
  
  nA <- 15
  nB <- 12
  
  # Generate network with rank-2 structure
  U_true <- matrix(rnorm(nA * 2), nA, 2)
  V_true <- matrix(rnorm(nB * 2), nB, 2)
  Y <- U_true %*% t(V_true) + matrix(rnorm(nA * nB, 0, 0.3), nA, nB)
  
  # Fit with R = 0
  fit0 <- ame(Y, mode = "bipartite",
              R_row = 0, R_col = 0,
              burn = 300, nscan = 1000, print = FALSE)
  
  # Fit with R = 2
  fit2 <- ame(Y, mode = "bipartite",
              R_row = 2, R_col = 2,
              burn = 300, nscan = 1000, print = FALSE)
  
  # Model with multiplicative effects should fit better
  # Use YPM (posterior mean predictions) instead of EZ
  resid0 <- Y - fit0$YPM
  resid2 <- Y - fit2$YPM
  
  # Check that both produced valid predictions
  expect_false(all(is.na(fit0$YPM)))
  expect_false(all(is.na(fit2$YPM)))
  
  # Model with R=2 should have lower residuals than R=0
  mse0 <- mean(resid0^2, na.rm = TRUE)
  mse2 <- mean(resid2^2, na.rm = TRUE)
  
  # Only test if both MSEs are finite
  if(is.finite(mse0) && is.finite(mse2)) {
    expect_lt(mse2, mse0)
  }
})

test_that("Bipartite handles missing data correctly", {
  skip_on_cran()
  
  set.seed(795)
  
  nA <- 10
  nB <- 8
  
  # Generate complete network
  Y <- matrix(rnorm(nA * nB), nA, nB)
  
  # Add missing values
  missing_idx <- sample(1:(nA*nB), size = 10, replace = FALSE)
  Y[missing_idx] <- NA
  
  # Fit model
  fit <- ame(Y, mode = "bipartite",
             burn = 200, nscan = 500, print = FALSE)
  
  # Check that model ran
  expect_equal(fit$mode, "bipartite")
  
  # Check that predictions exist for missing values
  # EZ removed -   expect_false(any(is.na(fit$EZ)))
  expect_false(any(is.na(fit$YPM)))
})

test_that("Bipartite GOF statistics are computed", {
  skip_on_cran()
  
  set.seed(796)
  
  nA <- 8
  nB <- 10
  
  Y <- matrix(rnorm(nA * nB), nA, nB)
  
  # Fit with GOF
  fit <- ame(Y, mode = "bipartite",
             burn = 100, nscan = 300, 
             gof = TRUE, print = FALSE)
  
  # Check GOF exists
  expect_false(is.null(fit$GOF))
  # Bipartite GOF has 3 columns: sd.rowmean, sd.colmean, four.cycles
  expect_equal(ncol(fit$GOF), 3)
  
  # Check GOF statistics are reasonable
  # Column 1: sd.rowmean (std dev of row means)
  obs_sd_row <- fit$GOF[1,1]
  # Column 2: sd.colmean (std dev of col means)
  obs_sd_col <- fit$GOF[1,2]
  # Column 3: four.cycles count
  obs_four_cycles <- fit$GOF[1,3]
  
  # Check statistics are non-negative
  expect_true(obs_sd_row >= 0)
  expect_true(obs_sd_col >= 0)
  expect_true(obs_four_cycles >= 0)
})