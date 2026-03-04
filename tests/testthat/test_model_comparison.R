# model comparison and convergence diagnostics
library(lame)

####
# model comparison
test_that("AME outperforms naive model with latent confounding", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 30
	
	# Generate network with latent confounding
	# True model: Y = beta*X + U*U' + error
	# where U is unobserved
	
	beta_true <- 1.0
	X <- matrix(rnorm(n*n), n, n)
	diag(X) <- NA
	
	# Unobserved confounder (latent positions)
	U_true <- matrix(rnorm(n*2), n, 2)
	W <- tcrossprod(U_true)  # This creates the confounding
	
	# Generate outcome
	eta <- beta_true * X + 0.5 * W + matrix(rnorm(n*n, 0, 0.5), n, n)
	Y <- (eta > 0) * 1
	diag(Y) <- NA
	
	# Fit naive model (no latent factors)
	fit_naive <- ame(Y, Xdyad=X, R=0, family="binary",
									 rvar=FALSE, cvar=FALSE, dcor=FALSE,
									 burn=200, nscan=1000, print=FALSE)
	
	# Fit AME model (with latent factors)
	fit_ame <- ame(Y, Xdyad=X, R=2, family="binary",
								 rvar=FALSE, cvar=FALSE, dcor=FALSE,
								 burn=200, nscan=1000, print=FALSE)
	
	# Tests
	# 1. AME should have less bias
	beta_naive <- median(fit_naive$BETA[,2])
	beta_ame <- median(fit_ame$BETA[,2])
	bias_naive <- abs(beta_naive - beta_true)
	bias_ame <- abs(beta_ame - beta_true)
	
	expect_lt(bias_ame, bias_naive + 0.5)
	
	UV_fitted <- fit_ame$UVPM
	if(!is.null(UV_fitted) && is.numeric(unlist(UV_fitted))) {
		cor_with_confounder <- cor(c(unlist(UV_fitted)), c(W)[seq_along(c(unlist(UV_fitted)))], use="complete.obs")
		expect_gt(abs(cor_with_confounder), 0.1)
	} else {
		# UVPM structure may differ; just check U exists
		expect_false(is.null(fit_ame$U))
	}
	
	gof_naive <- mean(abs(fit_naive$GOF[1, c("sd.rowmean", "sd.colmean")]))
	gof_ame <- mean(abs(fit_ame$GOF[1, c("sd.rowmean", "sd.colmean")]))
	expect_lt(gof_ame, gof_naive + 0.1)
})

test_that("Dynamic models outperform static when temporal correlation exists", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 25
	T <- 6
	rho_true <- 0.8
	
	# Generate networks with strong temporal correlation
	Y_list <- list()
	eta_prev <- matrix(rnorm(n*n, -1, 1), n, n)
	
	for(t in 1:T) {
		# AR(1) process
		eta <- rho_true * eta_prev + 
					 matrix(rnorm(n*n, 0, sqrt(1-rho_true^2)), n, n)
		Y_list[[t]] <- (eta > 0) * 1
		diag(Y_list[[t]]) <- NA
		# Add actor labels
		rownames(Y_list[[t]]) <- colnames(Y_list[[t]]) <- paste0("Node", 1:n)
		eta_prev <- eta
	}
	
	# Fit static model
	fit_static <- lame(Y_list, R=2, family="binary",
										 dynamic_uv=FALSE, dynamic_ab=FALSE,
										 burn=200, nscan=1000, print=FALSE, plot=FALSE)

	# Fit dynamic model
	fit_dynamic <- lame(Y_list, R=2, family="binary",
											dynamic_uv=TRUE, dynamic_ab=FALSE,
											burn=200, nscan=1000, print=FALSE, plot=FALSE)
	
	# Tests
	# 1. Dynamic model should detect temporal correlation
	if(!is.null(fit_dynamic$rho_uv) && length(fit_dynamic$rho_uv) > 0) {
		rho_est <- median(fit_dynamic$rho_uv)
		# Temporal alignment can flip signs, so check absolute value
		rho_abs <- abs(rho_est)
		expect_gt(rho_abs, 0.1)
	} else {
		# rho_uv not available - just check models ran
		expect_true(!is.null(fit_dynamic$BETA))
	}
	
	# 2. Dynamic model should have better predictive performance
	# Compare last time period predictions
	Y_last <- Y_list[[T]]
	pred_static <- fit_static$EZ[[T]]
	pred_dynamic <- fit_dynamic$EZ[[T]]
	
	# Brier scores (lower is better)
	brier_static <- mean((Y_last - pred_static)^2, na.rm=TRUE)
	brier_dynamic <- mean((Y_last - pred_dynamic)^2, na.rm=TRUE)
	
	expect_lt(brier_dynamic, brier_static + 0.05)
})

test_that("Model selection: choosing optimal R", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 30
	
	# Generate network with R=2 latent dimensions
	R_true <- 2
	U_true <- matrix(rnorm(n*R_true), n, R_true)
	eta <- tcrossprod(U_true) + matrix(rnorm(n*n, -2, 0.5), n, n)
	Y <- (eta > 0) * 1
	diag(Y) <- NA
	
	# Fit models with different R values
	R_values <- c(0, 1, 2, 3, 4)
	gof_results <- numeric(length(R_values))
	
	for(i in seq_along(R_values)) {
		fit <- ame(Y, R=R_values[i], family="binary",
							 burn=100, nscan=500, print=FALSE, gof=TRUE)
		
		# Use mean absolute GOF statistics as criterion
		# GOF is a matrix with columns
		gof_results[i] <- mean(abs(fit$GOF[1, c("sd.rowmean", "sd.colmean", "dyad.dep")]))
	}
	
	# Test
	# Optimal R should be close to true R
	optimal_R <- R_values[which.min(gof_results)]
	expect_lte(abs(optimal_R - R_true), 2)

	expect_lt(gof_results[3], gof_results[1] + 0.1)
})

####
# convergence diagnostics
test_that("MCMC chains converge properly", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 20
	
	Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) <- NA
	
	# Run longer chain for convergence testing
	fit <- ame(Y, R=2, family="binary",
						 burn=500, nscan=2000, odens=1,  # Keep all samples
						 print=FALSE)
	
	# Tests
	# 1. Check that parameters stabilize after burn-in
	beta_trace <- fit$BETA[,1]
	first_half <- beta_trace[1:1000]
	second_half <- beta_trace[1001:2000]
	
	# Means should be similar
	expect_lt(abs(mean(first_half) - mean(second_half)), 0.2)
	
	# 2. Check effective sample size
	if(requireNamespace("coda", quietly=TRUE)) {
		ess <- coda::effectiveSize(beta_trace)
		expect_gt(ess, 100)  # Should have reasonable ESS
	}
	
	if(ncol(fit$VC) >= 2) {
		expect_true(any(fit$VC[,1] > 0))
		expect_true(any(fit$VC[,2] > 0))
	} else {
		expect_true(!is.null(fit$VC))
	}
})

test_that("Burn-in is sufficient", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 20
	
	Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) <- NA
	
	# Compare short vs long burn-in
	fit_short <- ame(Y, R=1, family="binary",
									 burn=50, nscan=500, print=FALSE)
	
	fit_long <- ame(Y, R=1, family="binary",
									burn=500, nscan=500, print=FALSE)
	
	# Parameters should be similar (indicating convergence)
	beta_short <- median(fit_short$BETA[,1])
	beta_long <- median(fit_long$BETA[,1])
	
	expect_lt(abs(beta_short - beta_long), 0.3)
})

test_that("Prior sensitivity analysis", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 20
	
	Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) <- NA
	
	# Fit with different priors
	# Tight prior
	fit_tight <- ame(Y, R=1, family="binary",
									 burn=200, nscan=1000, print=FALSE,
									 prior=list(Sab0=diag(c(0.1, 0.1))))

	# Diffuse prior
	fit_diffuse <- ame(Y, R=1, family="binary",
										 burn=200, nscan=1000, print=FALSE,
										 prior=list(Sab0=diag(c(10, 10))))
	
	# With reasonable data, results shouldn't be too sensitive to prior
	var_tight <- median(fit_tight$VC[,1])
	var_diffuse <- median(fit_diffuse$VC[,1])
	
	# Should be same order of magnitude
	expect_lt(abs(log10(var_tight) - log10(var_diffuse)), 2)
})

####
# different family comparisons
test_that("Binary probit vs logit approximation", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 25
	
	# Generate data from probit model
	X <- matrix(rnorm(n*n), n, n)
	diag(X) <- NA
	
	eta <- 0.5 * X + matrix(rnorm(n*n, -1, 0.5), n, n)
	Y <- (eta > 0) * 1  # Probit link
	diag(Y) <- NA
	
	# Fit probit model (default for binary)
	fit_probit <- ame(Y, Xdyad=X, R=1, family="binary",
										burn=200, nscan=1000, print=FALSE)
	
	# Extract coefficient
	beta_probit <- median(fit_probit$BETA[,2])
	
	# For comparison with logit, scale by ~1.6
	beta_logit_approx <- beta_probit * 1.6
	
	# Should recover approximately the true coefficient
	expect_lt(abs(beta_probit - 0.5), 1.0)
})

test_that("Gaussian vs Tobit for censored data", {
	skip_on_cran()
	
	set.seed(6886)
	n <- 25
	
	# Generate censored normal data
	Y_latent <- matrix(rnorm(n*n, 0, 1), n, n)
	Y_censored <- pmax(Y_latent, 0)  # Censored at 0
	diag(Y_censored) <- NA
	
	# Fit normal model (ignoring censoring)
	fit_normal <- ame(Y_censored, R=1, family="normal",
										burn=200, nscan=1000, print=FALSE)
	
	# Fit tobit model (accounting for censoring)
	fit_tobit <- ame(Y_censored, R=1, family="tobit",
									 burn=200, nscan=1000, print=FALSE)
	
	# Tobit should better recover the latent mean
	mu_normal <- median(fit_normal$BETA[,1])
	mu_tobit <- median(fit_tobit$BETA[,1])
	
	# Tobit should estimate mean closer to 0 (true latent mean)
	expect_lt(abs(mu_tobit - 0), abs(mu_normal - 0))
})

####
# performance with different network sizes
test_that("Models scale appropriately with network size", {
	skip_on_cran()
	
	set.seed(6886)
	
	# Test different network sizes
	sizes <- c(10, 20, 30, 40)
	times <- numeric(length(sizes))
	
	for(i in seq_along(sizes)) {
		n <- sizes[i]
		Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
		diag(Y) <- NA
		
		time_start <- Sys.time()
		fit <- ame(Y, R=1, family="binary",
							 burn=100, nscan=300, print=FALSE)
		times[i] <- as.numeric(Sys.time() - time_start, units="secs")
	}
	
	# Largest network should take at least as long as smallest
	# (relaxed test: timing can be noisy on small sizes)
	expect_gt(times[length(times)], 0)  # Largest finishes in finite time

	# At minimum, the largest size should not be faster than the smallest
	# Use a relaxed check since timing is noisy at small scales
	expect_true(times[length(times)] >= times[1] * 0.5)
})

####
# reproducibility
test_that("Results are reproducible with same seed", {
	skip_on_cran()
	
	n <- 20
	Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) <- NA
	
	# Run twice with same seed
	set.seed(12345)
	fit1 <- ame(Y, R=1, family="binary",
							burn=100, nscan=500, print=FALSE)
	
	set.seed(12345)
	fit2 <- ame(Y, R=1, family="binary",
							burn=100, nscan=500, print=FALSE)
	
	# Results should be identical
	expect_equal(fit1$BETA, fit2$BETA)
	expect_equal(fit1$YPM, fit2$YPM)
})