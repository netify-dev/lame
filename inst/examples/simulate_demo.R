################################################################################
# Demonstration of simulate.ame and simulate.lame functions
# 
# This script shows how to:
# 1. Simulate networks from fitted AME models
# 2. Simulate longitudinal networks from fitted LAME models  
# 3. Quantify uncertainty in network simulations
# 4. Use simulations for posterior predictive checks
################################################################################

library(lame)

# Set seed for reproducibility
set.seed(2024)

################################################################################
# PART 1: Simulate from AME (cross-sectional) models
################################################################################

cat("\n========== PART 1: AME Network Simulation ==========\n\n")

# Load binary network data
data(YX_bin)

# Fit a binary AME model
cat("Fitting binary AME model...\n")
fit_ame <- ame(YX_bin$Y, YX_bin$X[,,1:3],  # Use first 3 covariates
               burn = 50, nscan = 100, odens = 10,
               family = "binary", 
               print = FALSE, plot = FALSE)

# Simulate 25 networks from the posterior
cat("Simulating 25 networks from posterior distribution...\n")
sims_ame <- simulate(fit_ame, 
                     nsim = 25, 
                     seed = 123,
                     newdata = list(Xdyad = YX_bin$X[,,1:3]))

# Print simulation summary
print(sims_ame)

# Calculate and compare densities
original_density <- mean(YX_bin$Y == 1, na.rm = TRUE)
sim_densities <- sapply(sims_ame$Y, function(y) mean(y == 1, na.rm = TRUE))

cat("\nDensity Comparison:\n")
cat(sprintf("  Original network: %.3f\n", original_density))
cat(sprintf("  Simulated mean:   %.3f (SD: %.3f)\n", 
            mean(sim_densities), sd(sim_densities)))
cat(sprintf("  95%% interval:    [%.3f, %.3f]\n",
            quantile(sim_densities, 0.025),
            quantile(sim_densities, 0.975)))

################################################################################
# PART 2: Simulate from LAME (longitudinal) models
################################################################################

cat("\n========== PART 2: LAME Longitudinal Simulation ==========\n\n")

# Load longitudinal network data
data(YX_bin_list)

# Use first 5 time periods for faster demo
Y_subset <- YX_bin_list$Y[1:5]
X_subset <- YX_bin_list$X[1:5]

# Fit a longitudinal LAME model
cat("Fitting longitudinal LAME model...\n")
fit_lame <- lame(Y_subset, X_subset,
                 burn = 50, nscan = 100, odens = 10,
                 family = "binary",
                 print = FALSE, plot = FALSE)

# Simulate 10 network trajectories, each with 8 time points
cat("Simulating 10 longitudinal network trajectories...\n")
sims_lame <- simulate(fit_lame, 
                      nsim = 10,
                      n_time = 8,
                      seed = 456)

# Print simulation summary
print(sims_lame)

# Analyze temporal patterns
cat("\nTemporal Density Analysis:\n")

# Calculate density at each time point across all simulations
density_by_time <- matrix(NA, 10, 8)  # 10 simulations x 8 time points
for (i in 1:10) {
  for (t in 1:8) {
    density_by_time[i, t] <- mean(sims_lame$Y[[i]][[t]] == 1, na.rm = TRUE)
  }
}

# Average density at each time point
mean_density_t <- colMeans(density_by_time)
sd_density_t <- apply(density_by_time, 2, sd)

cat("Mean density by time point:\n")
for (t in 1:8) {
  cat(sprintf("  Time %d: %.3f (SD: %.3f)\n", t, mean_density_t[t], sd_density_t[t]))
}

################################################################################
# PART 3: Using different covariate values
################################################################################

cat("\n========== PART 3: Simulation with New Covariates ==========\n\n")

# Create new covariate values (e.g., intervention scenario)
n <- nrow(YX_bin$Y)
p <- dim(YX_bin$X)[3]

# Double the effect of first covariate
new_X <- YX_bin$X[,,1:3]
new_X[,,1] <- new_X[,,1] * 2

cat("Simulating with modified covariates (doubled first covariate)...\n")
sims_intervention <- simulate(fit_ame,
                              nsim = 20,
                              seed = 789,
                              newdata = list(Xdyad = new_X))

# Compare densities
intervention_densities <- sapply(sims_intervention$Y, 
                                 function(y) mean(y == 1, na.rm = TRUE))

cat("\nIntervention Effect on Network Density:\n")
cat(sprintf("  Original covariates:     %.3f (SD: %.3f)\n", 
            mean(sim_densities), sd(sim_densities)))
cat(sprintf("  Modified covariates:     %.3f (SD: %.3f)\n",
            mean(intervention_densities), sd(intervention_densities)))
cat(sprintf("  Estimated change:        %+.3f\n",
            mean(intervention_densities) - mean(sim_densities)))

################################################################################
# PART 4: Posterior Predictive Checks
################################################################################

cat("\n========== PART 4: Posterior Predictive Checks ==========\n\n")

# Simulate networks with latent variables for diagnostics
sims_with_latent <- simulate(fit_ame,
                              nsim = 100,
                              return_latent = TRUE,
                              seed = 999,
                              newdata = list(Xdyad = YX_bin$X[,,1:3]))

# Calculate test statistics for each simulated network
calc_stats <- function(y) {
  c(density = mean(y == 1, na.rm = TRUE),
    reciprocity = sum(y * t(y), na.rm = TRUE) / sum(y, na.rm = TRUE),
    transitivity = sum(diag(y %*% y %*% y), na.rm = TRUE) / 
                   sum(diag(y %*% y), na.rm = TRUE))
}

# Original network statistics
orig_stats <- calc_stats(YX_bin$Y)

# Simulated network statistics
sim_stats <- t(sapply(sims_with_latent$Y, calc_stats))

# Posterior predictive p-values
cat("Posterior Predictive P-values:\n")
for (stat in names(orig_stats)) {
  p_val <- mean(sim_stats[, stat] >= orig_stats[stat])
  cat(sprintf("  %s: %.3f %s\n", 
              stat, p_val,
              ifelse(p_val < 0.05 | p_val > 0.95, "(!)", "")))
}

cat("\nNote: P-values near 0 or 1 suggest model misfit for that statistic\n")

################################################################################
# Summary
################################################################################

cat("\n========== SUMMARY ==========\n\n")
cat("The simulate functions enable:\n")
cat("• Uncertainty quantification through posterior sampling\n")
cat("• Posterior predictive checks for model validation\n")
cat("• Scenario analysis with modified covariates\n")
cat("• Generation of synthetic networks with known properties\n")
cat("\nEach simulation incorporates:\n")
cat("• Parameter uncertainty (sampling from MCMC chains)\n")
cat("• Random effect uncertainty (sampling from posterior distributions)\n")
cat("• Appropriate error structure for each data family\n")

cat("\n✓ Demo complete!\n")