# Comprehensive tests for covariates across all families
# This test file was created after discovering that reconstruct_EZ
# was not including covariate effects due to missing X in fit object

test_that("Binary AME correctly handles covariates in all functions", {
  set.seed(123)
  n <- 15
  
  # Generate covariate
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  # True model: logit(p) = -0.5 + 0.8*X
  eta_true <- -0.5 + 0.8 * X
  p_true <- pnorm(eta_true)  # probit link
  # Handle NAs when generating Y
  Y <- matrix(NA, n, n)
  non_na_idx <- which(!is.na(p_true))
  Y[non_na_idx] <- rbinom(length(non_na_idx), 1, p_true[non_na_idx])
  diag(Y) <- NA
  
  # Fit model with covariates
  fit <- ame(Y, Xdyad = X, R = 0, family = "binary",
            burn = 200, nscan = 500, print = FALSE, seed = 456)
  
  # Test 1: X is stored in fit object
  expect_false(is.null(fit$X), label = "X should be stored in fit object")
  expect_equal(dim(fit$X), c(n, n, 2))  # intercept + covariate
  
  # Test 2: reconstruct_EZ includes covariate effects
  EZ <- reconstruct_EZ(fit)
  expect_false(is.null(EZ))
  
  # Should have reasonable correlation with true eta
  corr <- cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
  expect_gt(corr, 0.5)
  
  # Test 3: predict(type="link") matches reconstruct_EZ
  pred_link <- predict(fit, type = "link")
  expect_equal(pred_link, EZ, tolerance = 1e-10,
              info = "predict(type='link') should match reconstruct_EZ")
  
  # Test 4: predict(type="response") gives probabilities
  pred_resp <- predict(fit, type = "response")
  expect_true(all(pred_resp >= 0 & pred_resp <= 1, na.rm = TRUE),
              info = "Binary predictions should be probabilities")
  
  # Test 5: YPM approximates Φ(EZ)
  expected_probs <- pnorm(EZ)
  mean_diff <- mean(abs(fit$YPM - expected_probs), na.rm = TRUE)
  expect_lt(mean_diff, 0.1, label = "YPM should approximate Φ(EZ)")
})

test_that("Poisson AME correctly handles covariates in all functions", {
  set.seed(234)
  n <- 15
  
  # Generate covariate
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  
  # True model: log(λ) = 0.5 + 1.2*X
  eta_true <- 0.5 + 1.2 * X
  lambda_true <- exp(eta_true)
  # Handle NAs when generating Y
  Y <- matrix(NA, n, n)
  non_na_idx <- which(!is.na(lambda_true))
  Y[non_na_idx] <- rpois(length(non_na_idx), lambda_true[non_na_idx])
  diag(Y) <- NA
  
  # Fit model with covariates
  fit <- ame(Y, Xdyad = X, R = 0, family = "poisson",
            burn = 200, nscan = 500, odens = 1, print = FALSE, seed = 567)
  
  # Test 1: X is stored in fit object
  expect_false(is.null(fit$X), label = "X should be stored in fit object")
  
  # Test 2: reconstruct_EZ includes covariate effects
  EZ <- reconstruct_EZ(fit)
  expect_false(is.null(EZ))
  
  # Should have reasonable correlation with true eta
  corr <- cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
  expect_gt(corr, 0.5, label = "reconstruct_EZ should correlate with true log(lambda)")
  
  # Test 3: predict(type="link") matches reconstruct_EZ
  pred_link <- predict(fit, type = "link")
  expect_equal(pred_link, EZ, tolerance = 1e-10)
  
  # Test 4: YPM approximates exp(EZ)
  expected_lambda <- exp(EZ)
  ratio <- mean(fit$YPM / expected_lambda, na.rm = TRUE)
  expect_true(abs(ratio - 1) < 0.15, 
              info = paste("YPM/exp(EZ) ratio should be near 1, got", round(ratio, 3)))
  
  # Test 5: YPM is non-negative
  expect_true(all(fit$YPM >= 0, na.rm = TRUE),
              info = "Poisson YPM should be non-negative")
})

test_that("Gaussian AME correctly handles covariates in all functions", {
  set.seed(345)
  n <- 15
  
  # Generate covariates (2 predictors)
  X1 <- matrix(rnorm(n*n), n, n)
  X2 <- matrix(rnorm(n*n), n, n)
  diag(X1) <- diag(X2) <- NA
  X_array <- array(c(X1, X2), dim = c(n, n, 2))
  
  # True model: Y = 1.0 + 0.5*X1 - 0.3*X2 + noise
  eta_true <- 1.0 + 0.5 * X1 - 0.3 * X2
  Y <- eta_true + matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(Y) <- NA
  
  # Fit model with covariates
  fit <- ame(Y, Xdyad = X_array, R = 0, family = "normal",
            burn = 200, nscan = 500, print = FALSE, seed = 678)
  
  # Test 1: X is stored in fit object
  expect_false(is.null(fit$X), label = "X should be stored in fit object")
  expect_equal(dim(fit$X)[3], 3)  # intercept + 2 covariates
  
  # Test 2: reconstruct_EZ includes covariate effects
  EZ <- reconstruct_EZ(fit)
  expect_false(is.null(EZ))
  
  # Should have high correlation with true eta
  corr <- cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
  expect_gt(corr, 0.8, label = "reconstruct_EZ should highly correlate with true mean")
  
  # Test 3: For normal family, YPM = EZ (identity link)
  mean_diff <- mean(abs(fit$YPM - EZ), na.rm = TRUE)
  expect_lt(mean_diff, 0.1, label = "For normal family, YPM should equal EZ")
  
  # Test 4: Residuals should be centered
  residuals <- Y - EZ
  expect_lt(abs(mean(residuals, na.rm = TRUE)), 0.1)
  
  # Test 5: Beta estimates should be reasonable
  beta_est <- colMeans(fit$BETA)
  expect_equal(length(beta_est), 3)  # intercept + 2 covariates
  # Check if estimates are in right ballpark
  expect_true(abs(beta_est[1] - 1.0) < 0.5, label = "Intercept estimate")
  expect_true(abs(beta_est[2] - 0.5) < 0.5, label = "Beta1 estimate")
  expect_true(abs(beta_est[3] - (-0.3)) < 0.5, label = "Beta2 estimate")
})

test_that("Tobit AME correctly handles covariates in all functions", {
  set.seed(456)
  n <- 12
  
  # Generate covariate
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  
  # True model: latent Y* = 0.5 + 0.8*X + noise
  eta_true <- 0.5 + 0.8 * X
  Y_latent <- eta_true + matrix(rnorm(n*n, 0, 0.5), n, n)
  Y <- matrix(pmax(0, Y_latent), n, n)  # Censor at 0
  diag(Y) <- NA
  
  # Fit model with covariates
  fit <- ame(Y, Xdyad = X, R = 0, family = "tobit",
            burn = 200, nscan = 500, print = FALSE, seed = 789)
  
  # Test 1: X is stored in fit object
  expect_false(is.null(fit$X), label = "X should be stored in fit object for tobit")
  
  # Test 2: reconstruct_EZ includes covariate effects
  EZ <- reconstruct_EZ(fit)
  expect_false(is.null(EZ))
  
  # Test 3: YPM is non-negative (censored)
  expect_true(all(fit$YPM >= 0, na.rm = TRUE),
              info = "Tobit YPM should be non-negative")
  
  # Test 4: EZ can be negative (latent scale)
  expect_true(any(EZ < 0, na.rm = TRUE) || min(EZ, na.rm = TRUE) < 0.1,
              info = "Tobit EZ (latent) can be negative")
})

test_that("Ordinal AME correctly handles covariates in all functions", {
  skip("Ordinal with covariates - implement after verifying ordinal model structure")
})

test_that("Bipartite models correctly handle covariates across families", {
  set.seed(567)
  nA <- 8
  nB <- 10
  
  # Generate covariates
  Xdyad <- array(rnorm(nA * nB * 2), c(nA, nB, 2))
  
  # Test binary bipartite
  eta_true <- 0.5 + 0.7 * Xdyad[,,1] - 0.4 * Xdyad[,,2]
  p_true <- pnorm(eta_true)
  Y_bin <- matrix(rbinom(nA * nB, 1, c(p_true)), nA, nB)
  
  fit_bin <- ame(Y_bin, Xdyad = Xdyad, mode = "bipartite", family = "binary",
                burn = 200, nscan = 400, print = FALSE, seed = 890)
  
  # Test 1: X is stored (bipartite already stores X, but verify)
  expect_false(is.null(fit_bin$X), label = "Bipartite should store X")
  
  # Test 2: reconstruct_EZ works
  EZ_bin <- reconstruct_EZ(fit_bin)
  expect_false(is.null(EZ_bin))
  expect_equal(dim(EZ_bin), c(nA, nB))
  
  # Test 3: Correlation with truth
  corr <- cor(c(eta_true), c(EZ_bin))
  expect_gt(corr, 0.4, label = "Bipartite binary EZ should correlate with truth")
  
  # Test Gaussian bipartite
  Y_norm <- eta_true + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
  
  fit_norm <- ame(Y_norm, Xdyad = Xdyad, mode = "bipartite", family = "normal",
                 burn = 200, nscan = 400, print = FALSE, seed = 901)
  
  EZ_norm <- reconstruct_EZ(fit_norm)
  corr_norm <- cor(c(eta_true), c(EZ_norm))
  expect_gt(corr_norm, 0.7, label = "Bipartite normal EZ should highly correlate with truth")
})

test_that("predict() function works correctly with covariates", {
  set.seed(678)
  n <- 12
  
  # Fit a model with covariates
  X <- matrix(rnorm(n*n), n, n)
  diag(X) <- NA
  Y <- 1 + 0.5 * X + matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, Xdyad = X, R = 1, family = "normal",
            burn = 200, nscan = 400, print = FALSE)
  
  # Test 1: predict with original data (no newdata)
  pred_link <- predict(fit, type = "link")
  pred_resp <- predict(fit, type = "response")
  
  expect_equal(pred_link, reconstruct_EZ(fit), tolerance = 1e-10)
  expect_equal(pred_resp, fit$YPM, tolerance = 1e-10)
  
  # Test 2: predict with new covariates
  X_new <- array(matrix(rnorm(n*n), n, n), dim = c(n, n, 1))
  diag(X_new[,,1]) <- NA
  
  pred_new <- predict(fit, newdata = X_new, type = "link")
  expect_false(is.null(pred_new))
  expect_equal(dim(pred_new), c(n, n))
  
  # Should be different from original predictions
  expect_true(mean(abs(pred_new - pred_link), na.rm = TRUE) > 0.1,
              info = "New predictions should differ from original")
})

test_that("simulate() function uses stored X correctly", {
  set.seed(789)
  n <- 10
  
  # Fit a model with covariates
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  Y <- 2 + 1.5 * X + matrix(rnorm(n*n, 0, 0.3), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, Xdyad = X, R = 0, family = "normal",
            burn = 200, nscan = 400, print = FALSE)
  
  # Test 1: simulate without warning about missing covariates
  warnings_caught <- character()
  withCallingHandlers({
    sims <- simulate(fit, nsim = 10, seed = 123)
  }, warning = function(w) {
    warnings_caught <<- c(warnings_caught, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  
  # Should not warn about missing covariates if X is stored
  expect_false(any(grepl("covariates not stored", warnings_caught)),
               label = "Should not warn about missing X when X is stored")
  
  # Test 2: simulated values should be reasonable
  sim_mean <- Reduce(`+`, sims$Y) / length(sims$Y)
  
  # Mean should be similar to original data
  expect_true(abs(mean(sim_mean, na.rm = TRUE) - mean(Y, na.rm = TRUE)) < 1,
              info = "Simulated mean should be close to original")
})

test_that("Y is available in fit object for gof()", {
  set.seed(890)
  n <- 10
  Y <- matrix(rnorm(n*n), n, n)
  diag(Y) <- NA
  
  # Test unipartite
  fit_uni <- ame(Y, R = 0, family = "normal",
                burn = 100, nscan = 200, print = FALSE)
  
  expect_false(is.null(fit_uni$Y), label = "Unipartite should store Y")
  
  # Test bipartite
  Y_bip <- matrix(rnorm(8*10), 8, 10)
  fit_bip <- ame(Y_bip, mode = "bipartite", family = "normal",
                burn = 100, nscan = 200, print = FALSE)
  
  expect_false(is.null(fit_bip$Y), label = "Bipartite should store Y")
  
  # Test that gof() works without providing Y
  gof_result <- gof(fit_uni, nsim = 5, verbose = FALSE)
  expect_false(is.null(gof_result))
  expect_true(nrow(gof_result) == 6)  # 1 observed + 5 simulated
})

