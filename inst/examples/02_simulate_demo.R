####
# simulate from fitted ame and lame models
library(lame)
set.seed(2024)
####

####
# part 1: simulate from ame (cross-sectional)
cat("\n========== PART 1: AME Network Simulation ==========\n\n")

# load data
data(YX_bin)

# fit binary ame
cat("Fitting binary AME model...\n")
fit_ame = ame(YX_bin$Y, YX_bin$X[,,1:3],
	burn = 50, nscan = 100, odens = 10,
	family = "binary",
	print = FALSE, plot = FALSE)

# simulate 25 networks
cat("Simulating 25 networks from posterior distribution...\n")
sims_ame = simulate(fit_ame,
	nsim = 25,
	seed = 123,
	newdata = list(Xdyad = YX_bin$X[,,1:3]))
print(sims_ame)

# density comparison
original_density = mean(YX_bin$Y == 1, na.rm = TRUE)
sim_densities = sapply(sims_ame$Y, function(y) mean(y == 1, na.rm = TRUE))

cat("\nDensity Comparison:\n")
cat(sprintf("  Original network: %.3f\n", original_density))
cat(sprintf("  Simulated mean:   %.3f (SD: %.3f)\n",
	mean(sim_densities), sd(sim_densities)))
cat(sprintf("  95%% interval:    [%.3f, %.3f]\n",
	quantile(sim_densities, 0.025),
	quantile(sim_densities, 0.975)))
####

####
# part 2: simulate from lame (longitudinal)
cat("\n========== PART 2: LAME Longitudinal Simulation ==========\n\n")

# load data
data(YX_bin_list)
Y_subset = YX_bin_list$Y[1:5]
X_subset = YX_bin_list$X[1:5]

# fit longitudinal model
cat("Fitting longitudinal LAME model...\n")
fit_lame = lame(Y_subset, X_subset,
	burn = 50, nscan = 100, odens = 10,
	family = "binary",
	print = FALSE, plot = FALSE)

# simulate 10 trajectories
cat("Simulating 10 longitudinal network trajectories...\n")
sims_lame = simulate(fit_lame,
	nsim = 10,
	n_time = 8,
	seed = 456)
print(sims_lame)

# temporal density
cat("\nTemporal Density Analysis:\n")
density_by_time = matrix(NA, 10, 8)
for (i in 1:10) {
	for (t in 1:8) {
		density_by_time[i, t] = mean(sims_lame$Y[[i]][[t]] == 1, na.rm = TRUE)
	}
}

mean_density_t = colMeans(density_by_time)
sd_density_t = apply(density_by_time, 2, sd)

cat("Mean density by time point:\n")
for (t in 1:8) {
	cat(sprintf("  Time %d: %.3f (SD: %.3f)\n", t, mean_density_t[t], sd_density_t[t]))
}
####

####
# part 3: simulation with modified covariates
cat("\n========== PART 3: Simulation with New Covariates ==========\n\n")

n = nrow(YX_bin$Y)
p = dim(YX_bin$X)[3]

# double the first covariate
new_X = YX_bin$X[,,1:3]
new_X[,,1] = new_X[,,1] * 2

cat("Simulating with modified covariates (doubled first covariate)...\n")
sims_intervention = simulate(fit_ame,
	nsim = 20,
	seed = 789,
	newdata = list(Xdyad = new_X))

intervention_densities = sapply(sims_intervention$Y,
	function(y) mean(y == 1, na.rm = TRUE))

cat("\nIntervention Effect on Network Density:\n")
cat(sprintf("  Original covariates:     %.3f (SD: %.3f)\n",
	mean(sim_densities), sd(sim_densities)))
cat(sprintf("  Modified covariates:     %.3f (SD: %.3f)\n",
	mean(intervention_densities), sd(intervention_densities)))
cat(sprintf("  Estimated change:        %+.3f\n",
	mean(intervention_densities) - mean(sim_densities)))
####

####
# part 4: posterior predictive checks
cat("\n========== PART 4: Posterior Predictive Checks ==========\n\n")

sims_with_latent = simulate(fit_ame,
	nsim = 100,
	return_latent = TRUE,
	seed = 999,
	newdata = list(Xdyad = YX_bin$X[,,1:3]))

# test statistics
calc_stats = function(y) {
	c(density = mean(y == 1, na.rm = TRUE),
		reciprocity = sum(y * t(y), na.rm = TRUE) / sum(y, na.rm = TRUE),
		transitivity = sum(diag(y %*% y %*% y), na.rm = TRUE) /
			sum(diag(y %*% y), na.rm = TRUE))
}

orig_stats = calc_stats(YX_bin$Y)
sim_stats = t(sapply(sims_with_latent$Y, calc_stats))

# posterior predictive p-values
cat("Posterior Predictive P-values:\n")
for (stat in names(orig_stats)) {
	p_val = mean(sim_stats[, stat] >= orig_stats[stat])
	cat(sprintf("  %s: %.3f %s\n",
		stat, p_val,
		ifelse(p_val < 0.05 | p_val > 0.95, "(!)", "")))
}
cat("\nNote: P-values near 0 or 1 suggest model misfit for that statistic\n")
####

####
# summary
cat("\n========== SUMMARY ==========\n\n")
cat("The simulate functions enable:\n")
cat("- Uncertainty quantification through posterior sampling\n")
cat("- Posterior predictive checks for model validation\n")
cat("- Scenario analysis with modified covariates\n")
cat("- Generation of synthetic networks with known properties\n")
cat("\nDemo complete.\n")
####
