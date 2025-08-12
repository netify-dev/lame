################################################################################
# Basic Test Script for simulate.ame and simulate.lame Functions
# 
# Four essential tests to verify simulation functionality
################################################################################

library(lame)

# Set up test tracking
tests_passed <- 0
tests_total <- 4

cat("Running basic simulation tests...\n")
cat("================================\n\n")

################################################################################
# TEST 1: simulate.ame works for binary networks
################################################################################

cat("TEST 1: simulate.ame for binary networks\n")
tryCatch({
  # Load data and fit model
  data(YX_bin)
  set.seed(6886)
  
  fit <- ame(YX_bin$Y, YX_bin$X[,,1:2], 
             burn = 20, nscan = 30, odens = 10,
             family = "binary", print = FALSE, plot = FALSE)
  
  # Simulate networks
  sims <- simulate(fit, nsim = 5, seed = 456, 
                   newdata = list(Xdyad = YX_bin$X[,,1:2]))
  
  # Check results
  stopifnot(class(sims)[1] == "ame.sim")
  stopifnot(length(sims$Y) == 5)
  stopifnot(all(dim(sims$Y[[1]]) == dim(YX_bin$Y)))
  
  cat("  ✓ Passed: Generated 5 binary networks with correct dimensions\n")
  tests_passed <- tests_passed + 1
  
}, error = function(e) {
  cat("  ✗ Failed:", e$message, "\n")
})

################################################################################
# TEST 2: simulate.ame captures uncertainty
################################################################################

cat("\nTEST 2: Uncertainty quantification in simulate.ame\n")
tryCatch({
  # Use previous fit
  set.seed(789)
  
  # Simulate multiple networks
  sims <- simulate(fit, nsim = 20, seed = 789,
                   newdata = list(Xdyad = YX_bin$X[,,1:2]))
  
  # Calculate densities
  densities <- sapply(sims$Y, function(y) mean(y == 1, na.rm = TRUE))
  
  # Check variation exists (uncertainty)
  stopifnot(sd(densities) > 0)
  stopifnot(length(unique(densities)) > 1)
  
  cat(sprintf("  ✓ Passed: Density varies across simulations (SD = %.4f)\n", 
              sd(densities)))
  tests_passed <- tests_passed + 1
  
}, error = function(e) {
  cat("  ✗ Failed:", e$message, "\n")
})

################################################################################
# TEST 3: simulate.lame works for longitudinal networks
################################################################################

cat("\nTEST 3: simulate.lame for longitudinal binary networks\n")
tryCatch({
  # Load data
  data(YX_bin_list)
  set.seed(6886)
  
  # Use subset for faster testing
  Y_sub <- YX_bin_list$Y[1:3]
  X_sub <- YX_bin_list$X[1:3]
  
  # Fit model
  fit_lame <- lame(Y_sub, X_sub,
                   burn = 20, nscan = 30, odens = 10,
                   family = "binary", print = FALSE, plot = FALSE)
  
  # Simulate trajectories - use default n_time from model
  sims_lame <- simulate(fit_lame, nsim = 5, seed = 456)
  
  # Check results
  stopifnot(class(sims_lame)[1] == "lame.sim")
  stopifnot(length(sims_lame$Y) == 5)  # 5 trajectories
  stopifnot(length(sims_lame$Y[[1]]) > 0)  # Has time points
  
  n_time_sim <- length(sims_lame$Y[[1]])
  cat(sprintf("  ✓ Passed: Generated 5 trajectories with %d time points each\n", 
              n_time_sim))
  tests_passed <- tests_passed + 1
  
}, error = function(e) {
  cat("  ✗ Failed:", e$message, "\n")
})

################################################################################
# TEST 4: simulate with return_latent option
################################################################################

cat("\nTEST 4: Return latent Z matrices\n")
tryCatch({
  # Use the AME fit from test 1
  data(YX_bin)
  set.seed(6886)
  
  fit <- ame(YX_bin$Y, YX_bin$X[,,1:2], 
             burn = 20, nscan = 30, odens = 10,
             family = "binary", print = FALSE, plot = FALSE)
  
  # Simulate with latent Z
  sims_z <- simulate(fit, nsim = 3, return_latent = TRUE, seed = 111,
                     newdata = list(Xdyad = YX_bin$X[,,1:2]))
  
  # Check Z matrices are included
  stopifnot(!is.null(sims_z$Z))
  stopifnot(length(sims_z$Z) == 3)
  stopifnot(all(dim(sims_z$Z[[1]]) == dim(YX_bin$Y)))
  
  # Z should be continuous, Y should be binary
  z_vals <- sims_z$Z[[1]][!is.na(sims_z$Z[[1]])]
  y_vals <- sims_z$Y[[1]][!is.na(sims_z$Y[[1]])]
  
  stopifnot(length(unique(z_vals)) > 2)  # Z is continuous
  stopifnot(all(y_vals %in% c(0, 1)))    # Y is binary
  
  cat("  ✓ Passed: Latent Z matrices returned correctly\n")
  tests_passed <- tests_passed + 1
  
}, error = function(e) {
  cat("  ✗ Failed:", e$message, "\n")
})

################################################################################
# Summary
################################################################################

cat("\n================================\n")
cat(sprintf("Tests Passed: %d/%d\n", tests_passed, tests_total))

if (tests_passed == tests_total) {
  cat("\n✓ All tests passed successfully!\n")
  cat("\nThe simulate functions are working correctly for:\n")
  cat("• Cross-sectional networks (simulate.ame)\n")
  cat("• Longitudinal networks (simulate.lame)\n")
  cat("• Uncertainty quantification\n")
  cat("• Latent variable extraction\n")
} else {
  cat("\n⚠ Some tests failed. Please check the error messages above.\n")
}