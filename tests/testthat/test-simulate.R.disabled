# Test simulate.ame and simulate.lame functions
# Basic tests to ensure simulation functions work correctly

library(testthat)
library(lame)

test_that("simulate.ame works for binary networks", {
  # Load data and fit model
  data(YX_bin)
  set.seed(6886)
  
  # Fit a simple binary AME model
  fit <- ame(YX_bin$Y, YX_bin$X[,,1:2], 
             burn = 20, nscan = 30, odens = 10,
             family = "binary", print = FALSE, plot = FALSE)
  
  # Test 1: Basic simulation works
  sims <- simulate(fit, nsim = 5, seed = 456, 
                   newdata = list(Xdyad = YX_bin$X[,,1:2]))
  
  expect_s3_class(sims, "ame.sim")
  expect_equal(length(sims$Y), 5)
  expect_equal(sims$family, "binary")
  expect_equal(sims$mode, "unipartite")
  
  # Test 2: Simulated networks have correct dimensions
  expect_equal(dim(sims$Y[[1]]), dim(YX_bin$Y))
  
  # Test 3: Binary networks contain only 0s, 1s, and NAs
  vals <- unique(as.vector(sims$Y[[1]]))
  vals <- vals[!is.na(vals)]
  expect_true(all(vals %in% c(0, 1)))
  
  # Test 4: Diagonal is NA (for unipartite)
  expect_true(all(is.na(diag(sims$Y[[1]]))))
})

test_that("simulate.ame captures uncertainty", {
  # Load data and fit model
  data(YX_bin)
  set.seed(6886)
  
  fit <- ame(YX_bin$Y, YX_bin$X[,,1:2],
             burn = 20, nscan = 30, odens = 10,
             family = "binary", print = FALSE, plot = FALSE)
  
  # Simulate multiple networks
  sims <- simulate(fit, nsim = 20, seed = 789,
                   newdata = list(Xdyad = YX_bin$X[,,1:2]))
  
  # Calculate densities
  densities <- sapply(sims$Y, function(y) mean(y == 1, na.rm = TRUE))
  
  # Test 1: Different simulations produce different networks
  expect_true(sd(densities) > 0)
  
  # Test 2: Mean density is reasonable
  orig_density <- mean(YX_bin$Y == 1, na.rm = TRUE)
  expect_true(abs(mean(densities) - orig_density) < 0.1)
  
  # Test 3: Can return latent Z matrices
  sims_z <- simulate(fit, nsim = 3, return_latent = TRUE, seed = 111,
                     newdata = list(Xdyad = YX_bin$X[,,1:2]))
  expect_false(is.null(sims_z$Z))
  expect_equal(length(sims_z$Z), 3)
})

test_that("simulate.lame works for longitudinal binary networks", {
  # Load data
  data(YX_bin_list)
  set.seed(6886)
  
  # Use subset for faster testing
  Y_sub <- YX_bin_list$Y[1:3]
  X_sub <- YX_bin_list$X[1:3]
  
  # Fit longitudinal model
  fit <- lame(Y_sub, X_sub,
              burn = 20, nscan = 30, odens = 10,
              family = "binary", print = FALSE, plot = FALSE)
  
  # Test 1: Basic simulation works
  sims <- simulate(fit, nsim = 5, n_time = 4, seed = 456,
                   newdata = list(Xdyad = X_sub))
  
  expect_s3_class(sims, "lame.sim")
  expect_equal(length(sims$Y), 5)  # 5 trajectories
  expect_equal(sims$n_time, 4)     # 4 time points each
  expect_equal(sims$family, "binary")
  
  # Test 2: Each trajectory has correct structure
  expect_equal(length(sims$Y[[1]]), 4)  # 4 time points
  
  # Test 3: Networks have correct dimensions at each time
  for (t in 1:4) {
    expect_equal(dim(sims$Y[[1]][[t]]), dim(Y_sub[[1]]))
  }
  
  # Test 4: Can specify different starting points
  sims_random <- simulate(fit, nsim = 3, n_time = 3, 
                          start_from = "random", seed = 222)
  expect_equal(length(sims_random$Y), 3)
})

test_that("simulate.lame captures temporal dynamics", {
  # Load data
  data(YX_bin_list)
  set.seed(6886)
  
  # Use subset
  Y_sub <- YX_bin_list$Y[1:3]
  X_sub <- YX_bin_list$X[1:3]
  
  # Fit model with dynamic effects (if supported)
  fit <- lame(Y_sub, X_sub,
              burn = 20, nscan = 30, odens = 10,
              family = "binary", print = FALSE, plot = FALSE)
  
  # Simulate longer time series
  sims <- simulate(fit, nsim = 10, n_time = 10, seed = 333)
  
  # Test 1: Time series have requested length
  expect_equal(length(sims$Y[[1]]), 10)
  
  # Test 2: Calculate densities over time for each trajectory
  density_matrix <- matrix(NA, 10, 10)  # 10 sims x 10 time points
  
  for (i in 1:10) {
    for (t in 1:10) {
      density_matrix[i, t] <- mean(sims$Y[[i]][[t]] == 1, na.rm = TRUE)
    }
  }
  
  # Test 3: Densities vary over time and across simulations
  # (shows both temporal dynamics and uncertainty)
  time_variation <- apply(density_matrix, 1, sd)  # SD across time for each sim
  sim_variation <- apply(density_matrix, 2, sd)   # SD across sims for each time
  
  expect_true(mean(time_variation) > 0)  # Variation over time
  expect_true(mean(sim_variation) > 0)   # Variation across simulations
  
  # Test 4: Print and summary methods work
  expect_no_error(print(sims))
  expect_no_error(summary(sims))
})