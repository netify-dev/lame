# Test simulate.ame and simulate.lame functions
library(lame)

cat("Testing simulation functions for AME and LAME models\n")
cat("====================================================\n\n")

# Test 1: simulate.ame with binary data
cat("Test 1: Binary AME model simulation\n")
data(YX_bin)
set.seed(6886)
fit_bin <- ame(YX_bin$Y, YX_bin$X, burn=50, nscan=100, odens=10, 
               family="binary", print=FALSE, plot=FALSE)

# Simulate from fitted model
sims_bin <- simulate(fit_bin, nsim=10, seed=456)
print(sims_bin)
cat("  Simulated", length(sims_bin$Y), "binary networks\n")
cat("  First network density:", mean(sims_bin$Y[[1]], na.rm=TRUE), "\n\n")

# Test 2: simulate.lame with binary longitudinal data
cat("Test 2: Binary LAME model simulation\n")
data(YX_bin_list)
# Take subset for faster testing
Y_sub <- YX_bin_list$Y[1:3]
X_sub <- YX_bin_list$X[1:3]
set.seed(6886)
fit_lame_bin <- lame(Y_sub, X_sub, burn=50, nscan=100, odens=10,
                     family="binary", print=FALSE, plot=FALSE)

# Simulate from fitted model
sims_lame <- simulate(fit_lame_bin, nsim=5, n_time=5, seed=789)
print(sims_lame)
cat("  Simulated", length(sims_lame$Y), "network trajectories\n")
cat("  Each trajectory has", length(sims_lame$Y[[1]]), "time points\n")
cat("  First network at t=1 density:", 
    mean(sims_lame$Y[[1]][[1]], na.rm=TRUE), "\n\n")

# Test 3: simulate.ame with normal data
cat("Test 3: Normal AME model simulation\n")
data(YX_nrm)
set.seed(6886)
fit_nrm <- ame(YX_nrm$Y, YX_nrm$X, burn=50, nscan=100, odens=10,
               family="normal", print=FALSE, plot=FALSE)

sims_nrm <- simulate(fit_nrm, nsim=10, seed=111)
print(sims_nrm)
cat("  Mean of first simulated network:", 
    mean(sims_nrm$Y[[1]], na.rm=TRUE), "\n\n")

# Test 4: Test with new covariates
cat("Test 4: Simulation with new covariates\n")
# Create new covariates
n <- nrow(YX_bin$Y)
new_X <- array(rnorm(n * n * dim(YX_bin$X)[3]), 
               dim=dim(YX_bin$X))

sims_new <- simulate(fit_bin, nsim=5, newdata=list(Xdyad=new_X), seed=222)
cat("  Simulated", length(sims_new$Y), "networks with new covariates\n")
cat("  Density with new covariates:", 
    mean(sims_new$Y[[1]], na.rm=TRUE), "\n\n")

# Test 5: Return latent Z matrices
cat("Test 5: Simulation with latent Z matrices\n")
sims_with_z <- simulate(fit_bin, nsim=3, return_latent=TRUE, seed=333)
cat("  Returned", length(sims_with_z$Z), "latent Z matrices\n")
cat("  First Z matrix mean:", mean(sims_with_z$Z[[1]], na.rm=TRUE), "\n\n")

# Summary
cat("Summary of simulations\n")
summary(sims_bin)
summary(sims_lame)

cat("\n✓ All simulation tests completed successfully!\n")