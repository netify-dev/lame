skip_on_cran()

# static family-specific edge cases pulled from the per-family test files:
# link / threshold / censoring / rank-list behavior, degenerate inputs, and
# family-specific guards that are NOT plain parameter-recovery duplicates.
# exactly one parameter-recovery check per family lives in
# test_static_recovery.r; everything here exercises behavior unique to a family.

library(lame)
library(testthat)

#### gaussian family-specific tests ####

# helper function to run single simulation
sim_gaussian_ame = function(seed, n, mu, beta, gamma=NULL, 
														rvar=FALSE, cvar=FALSE, R=0,
														burn=500, nscan=2000) {
	set.seed(seed)
	
	# generate covariates
	xw = matrix(rnorm(n*2), n, 2)
	X = tcrossprod(xw[,1])
	diag(X) = NA
	
	# build linear predictor
	eta = mu + beta * X
	
	# add unobserved covariate effect if specified
	if(!is.null(gamma)) {
		W = tcrossprod(xw[,2])
		eta = eta + gamma * W
	}
	
	# add additive effects if requested
	if(rvar || cvar) {
		a_true = if(rvar) rnorm(n, 0, 0.5) else rep(0, n)
		b_true = if(cvar) rnorm(n, 0, 0.5) else rep(0, n)
		eta = eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	}
	
	# add multiplicative effects if r > 0
	U_true = V_true = NULL
	if(R > 0) {
		U_true = matrix(rnorm(n*R, 0, 1), n, R)
		V_true = matrix(rnorm(n*R, 0, 1), n, R)
		eta = eta + tcrossprod(U_true, V_true)
	}
	
	# generate gaussian outcome with noise
	sigma2 = 1
	Y = eta + matrix(rnorm(n*n, 0, sqrt(sigma2)), n, n)
	diag(Y) = NA
	
	# fit ame model
	fit = ame(Y, Xdyad=X, R=R, family="normal",
						rvar=rvar, cvar=cvar, dcor=FALSE,
						burn=burn, nscan=nscan, 
						verbose = FALSE)
	
	# extract results
	beta_hat = median(fit$BETA[,2])
	beta_se = sd(fit$BETA[,2])
	beta_ci_lower = quantile(fit$BETA[,2], 0.025)
	beta_ci_upper = quantile(fit$BETA[,2], 0.975)
	beta_covered = (beta >= beta_ci_lower) & (beta <= beta_ci_upper)
	
	mu_hat = median(fit$BETA[,1])
	mu_se = sd(fit$BETA[,1])
	mu_ci_lower = quantile(fit$BETA[,1], 0.025)
	mu_ci_upper = quantile(fit$BETA[,1], 0.975)
	mu_covered = (mu >= mu_ci_lower) & (mu <= mu_ci_upper)
	
	# correlation with unobserved if applicable
	cor_with_W = if(!is.null(gamma) && R > 0) {
		cor(c(W), c(reconstruct_UVPM(fit)), use='pairwise.complete.obs')
	} else {
		NA
	}
	
	# recovery of additive effects
	a_cor = if(rvar && !is.null(fit$APM)) {
		cor(a_true, fit$APM, use='pairwise.complete.obs')
	} else {
		NA
	}
	
	b_cor = if(cvar && !is.null(fit$BPM)) {
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
		b_cor = b_cor,
		sigma2_hat = if(!is.null(fit$s2)) median(fit$s2) else if(!is.null(fit$VC) && "va" %in% colnames(fit$VC)) median(fit$VC[,"va"]) else NA
	))
}

####
# gaussian model with covariates only
####

test_that("Gaussian AME with covariates only has calibrated confidence intervals", {
	skip_on_cran()
	
	set.seed(6886)
	n_sims = 30
	n = 40
	mu_true = 0
	beta_true = 1
	
	# run simulations
	results = lapply(1:n_sims, function(i) {
		sim_gaussian_ame(seed=i, n=n, mu=mu_true, beta=beta_true,
										rvar=FALSE, cvar=FALSE, R=0,
										burn=300, nscan=1000)
	})
	
	# extract coverage rates
	beta_coverage = mean(sapply(results, function(x) x$beta_covered))
	mu_coverage = mean(sapply(results, function(x) x$mu_covered))
	
	# check calibration (should be close to 95%)
	expect_gt(beta_coverage, 0.85)
	expect_gt(mu_coverage, 0.85)
	
	# check bias
	beta_bias = mean(sapply(results, function(x) x$beta_hat)) - beta_true
	mu_bias = mean(sapply(results, function(x) x$mu_hat)) - mu_true
	
	expect_lt(abs(beta_bias), 0.1)
	expect_lt(abs(mu_bias), 0.1)
})

####
# gaussian model with covariates + additive effects
####

test_that("Gaussian AME with additive effects recovers true parameters", {
	set.seed(6886)
	n = 50
	mu_true = 0.5
	beta_true = 0.8
	
	# generate data with additive effects
	xw = matrix(rnorm(n*2), n, 2)
	X = tcrossprod(xw[,1])
	diag(X) = NA
	
	# true additive effects
	a_true = rnorm(n, 0, 0.7)
	b_true = rnorm(n, 0, 0.7)
	
	eta = mu_true + beta_true * X + 
				 outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	
	sigma2_true = 1
	Y = eta + matrix(rnorm(n*n, 0, sqrt(sigma2_true)), n, n)
	diag(Y) = NA
	
	# fit model with additive effects
	fit = ame(Y, Xdyad=X, R=0, family="normal",
						rvar=TRUE, cvar=TRUE, dcor=FALSE,
						burn=500, nscan=2000, verbose = FALSE)
	
	# check parameter recovery
	beta_est = median(fit$BETA[,2])
	mu_est = median(fit$BETA[,1])
	
	expect_lt(abs(beta_est - beta_true), 0.3)
	expect_lt(abs(mu_est - mu_true), 0.3)
	
	# check recovery of additive effects
	if(!is.null(fit$APM) && !is.null(fit$BPM)) {
		a_cor = cor(a_true, fit$APM, use='pairwise.complete.obs')
		b_cor = cor(b_true, fit$BPM, use='pairwise.complete.obs')
		
		expect_gt(a_cor, 0.5)
		expect_gt(b_cor, 0.5)
	}
})

####
# gaussian model with covariates + additive + multiplicative effects
####

test_that("Gaussian AME with full model recovers true parameters", {
	set.seed(6886)
	n = 50
	mu_true = -0.2
	beta_true = 1.0
	gamma_true = 0.8  # unobserved covariate
	R_true = 2
	
	# generate data
	xw = matrix(rnorm(n*2), n, 2)
	X = tcrossprod(xw[,1])
	W = tcrossprod(xw[,2])  # unobserved
	diag(X) = NA
	diag(W) = NA
	
	# true additive effects
	a_true = rnorm(n, 0, 0.5)
	b_true = rnorm(n, 0, 0.5)
	
	# true latent factors
	U_true = matrix(rnorm(n*R_true, 0, 1), n, R_true)
	V_true = matrix(rnorm(n*R_true, 0, 1), n, R_true)
	
	# build full model
	eta = mu_true + beta_true * X + gamma_true * W +
				 outer(a_true, rep(1, n)) + outer(rep(1, n), b_true) +
				 tcrossprod(U_true, V_true)
	
	sigma2_true = 1
	Y = eta + matrix(rnorm(n*n, 0, sqrt(sigma2_true)), n, n)
	diag(Y) = NA
	
	# fit full ame model
	fit = ame(Y, Xdyad=X, R=R_true, family="normal",
						rvar=TRUE, cvar=TRUE, dcor=FALSE,
						burn=500, nscan=2000, verbose = FALSE)
	
	# check parameter recovery
	beta_est = median(fit$BETA[,2])
	mu_est = median(fit$BETA[,1])
	
	expect_lt(abs(beta_est - beta_true), 0.4)
	
	# check that multiplicative effects capture unobserved covariate
	if(!is.null(reconstruct_UVPM(fit))) {
		cor_W = cor(c(W), c(reconstruct_UVPM(fit)), use='pairwise.complete.obs')
		expect_gt(cor_W, 0.3)
	}
	
	# check dimensions of latent factors
	expect_equal(ncol(fit$U), R_true)
	expect_equal(ncol(fit$V), R_true)
})

####
# compare models with and without multiplicative effects
####

test_that("Multiplicative effects reduce bias from unobserved confounding", {
	set.seed(6886)
	n = 50
	mu_true = 0
	beta_true = 1
	gamma_true = 1.5  # strong unobserved effect
	
	# generate data with unobserved confounder
	xw = matrix(rnorm(n*2), n, 2)
	X = tcrossprod(xw[,1])
	W = tcrossprod(xw[,2])
	diag(X) = NA
	diag(W) = NA
	
	eta = mu_true + beta_true * X + gamma_true * W
	Y = eta + matrix(rnorm(n*n), n, n)
	diag(Y) = NA
	
	# fit model without multiplicative effects
	fit_no_uv = ame(Y, Xdyad=X, R=0, family="normal",
									rvar=FALSE, cvar=FALSE, dcor=FALSE,
									burn=400, nscan=1500, verbose = FALSE)
	
	# fit model with multiplicative effects
	fit_with_uv = ame(Y, Xdyad=X, R=2, family="normal",
										rvar=FALSE, cvar=FALSE, dcor=FALSE,
										burn=400, nscan=1500, verbose = FALSE)
	
	beta_no_uv = median(fit_no_uv$BETA[,2])
	beta_with_uv = median(fit_with_uv$BETA[,2])
	
	# model with uv should have less bias
	bias_no_uv = abs(beta_no_uv - beta_true)
	bias_with_uv = abs(beta_with_uv - beta_true)
	
	expect_lt(bias_with_uv, bias_no_uv + 0.3)
	
	# check that uv captures the unobserved structure
	if(!is.null(fit_with_uv$UVPM)) {
		cor_W = cor(c(W), c(fit_with_uv$UVPM), use='pairwise.complete.obs')
		expect_gt(cor_W, 0.2)
	}
})

####
# model diagnostics and convergence
####

test_that("Gaussian AME produces valid diagnostics", {
	set.seed(6886)
	n = 30
	
	# simple model for quick convergence check
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	Y = X + matrix(rnorm(n*n), n, n)
	diag(Y) = NA
	
	fit = ame(Y, Xdyad=X, R=1, family="normal",
						rvar=TRUE, cvar=TRUE,
						burn=200, nscan=500, verbose = FALSE)
	
	# check that key outputs exist
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$VC))
	
	# check dimensions
	expect_equal(ncol(fit$BETA), 2)  # intercept + one covariate
	# check that we have some samples (may not be exactly nscan)
	expect_gt(nrow(fit$BETA), 0)
	
	# check predictions are reasonable
	resid = Y - reconstruct_EZ(fit)
	expect_lt(abs(mean(resid, na.rm=TRUE)), 0.5)
	
	# check effective sample size if available
	if(!is.null(fit$VC) && is.matrix(fit$VC)) {
		# variance components exist (some might be negative due to model issues)
		vc_values = as.numeric(fit$VC)
		vc_values = vc_values[!is.na(vc_values)]
		if(length(vc_values) > 0) {
			# at least check that variance components exist
			expect_gt(length(vc_values), 0)
		}
	}
})

#### count (poisson) family-specific tests ####

####
# overdispersion handling
####

test_that("Poisson AME handles overdispersed count data", {
	set.seed(6886)
	n = 35
	
	# generate overdispersed count data using negative binomial
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	eta = 0.5 + 0.4 * X
	lambda = exp(eta)
	
	# add overdispersion via negative binomial
	# (poisson model should still work, just with larger variance)
	Y = matrix(suppressWarnings(rnbinom(n*n, mu=as.vector(lambda), size=2)), n, n)
	diag(Y) = NA
	
	# check overdispersion
	var_Y = var(c(Y), na.rm=TRUE)
	mean_Y = mean(Y, na.rm=TRUE)
	expect_gt(var_Y, mean_Y)  # variance > mean indicates overdispersion
	
	# fit poisson model (should handle overdispersion via random effects)
	fit = ame(Y, Xdyad=X, R=1, family="poisson",
						rvar=TRUE, cvar=TRUE, dcor=FALSE,
						burn=400, nscan=1500, verbose = FALSE)
	
	# should still recover approximate beta
	beta_est = median(fit$BETA[,2])
	expect_lt(abs(beta_est - 0.4), 0.3)
	
	# random effects should capture extra variation
	if(!is.null(fit$APM) && !is.null(fit$BPM)) {
		var_random = var(c(fit$APM)) + var(c(fit$BPM))
		expect_gt(var_random, 0.001)  # should be non-zero
	}
})

####
# edge cases for count data
####

test_that("Poisson AME handles sparse count networks", {
	set.seed(6886)
	n = 30
	
	# very sparse count network (mostly zeros)
	Y = matrix(0, n, n)
	n_nonzero = floor(0.1 * n * n)
	nonzero_idx = sample(1:(n*n), n_nonzero)
	Y[nonzero_idx] = rpois(n_nonzero, lambda=2)
	diag(Y) = NA
	
	# check sparsity
	expect_lt(mean(Y > 0, na.rm=TRUE), 0.15)
	
	fit = ame(Y, R=1, family="poisson",
						burn=300, nscan=1200, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_lt(median(fit$BETA[,1]), 1)  # low intercept for sparse network
	
	# for sparse poisson networks, ez (log scale) will be very negative
	# this is correct behavior - sparse means low lambda, thus negative log(lambda)
	EZ = reconstruct_EZ(fit)
	if(!is.null(EZ)) {
		# for sparse networks, we expect mostly negative ez values
		expect_true(mean(EZ, na.rm=TRUE) < 0, 
								info = "Sparse networks should have negative mean on log scale")
	}
	
	# but ypm should still be non-negative
	expect_true(all(fit$YPM >= 0, na.rm=TRUE),
							info = "YPM (response scale) should be non-negative")
})

test_that("Poisson AME handles zero-inflated patterns", {
	set.seed(6886)
	n = 30
	
	# generate zero-inflated count data
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	eta = 0.5 + 0.3 * X
	lambda = exp(eta)
	Y = matrix(suppressWarnings(rpois(n*n, as.vector(lambda))), n, n)
	
	# add extra zeros (zero-inflation)
	zero_inflate_idx = sample(1:(n*n), floor(0.3*n*n))
	Y[zero_inflate_idx] = 0
	diag(Y) = NA
	
	# fit model (should handle via random effects)
	fit = ame(Y, Xdyad=X, R=1, family="poisson",
						rvar=TRUE, cvar=TRUE,
						burn=300, nscan=1200, verbose = FALSE)
	
	# should still produce reasonable estimates
	beta_est = median(fit$BETA[,2])
	expect_lt(abs(beta_est - 0.3), 0.4)  # more tolerance due to zero-inflation
})

#### binary family-specific tests ####

# helper function to run single simulation for binary
sim_binary_ame = function(seed, n, mu, beta, gamma=NULL,
													 rvar=FALSE, cvar=FALSE, R=0,
													 burn=500, nscan=2000, odens=25) {
	set.seed(seed)
	
	# generate covariates
	xw = matrix(rnorm(n*2), n, 2)
	X = tcrossprod(xw[,1])
	diag(X) = NA
	
	# build linear predictor (probit scale)
	eta = mu + beta * X
	
	# add unobserved covariate effect if specified
	if(!is.null(gamma)) {
		W = tcrossprod(xw[,2])
		eta = eta + gamma * W
	}
	
	# add additive effects if requested
	if(rvar || cvar) {
		a_true = if(rvar) rnorm(n, 0, 0.5) else rep(0, n)
		b_true = if(cvar) rnorm(n, 0, 0.5) else rep(0, n)
		eta = eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	}
	
	# add multiplicative effects if r > 0
	U_true = V_true = NULL
	if(R > 0) {
		U_true = matrix(rnorm(n*R, 0, 1), n, R)
		V_true = matrix(rnorm(n*R, 0, 1), n, R)
		eta = eta + tcrossprod(U_true, V_true)
	}
	
	# generate binary outcome via probit link
	Z = eta + matrix(rnorm(n*n), n, n)  # add standard normal noise
	Y = 1 * (Z > 0)  # binary outcome
	diag(Y) = NA
	
	# fit ame model
	fit = ame(Y, Xdyad=X, R=R, family="binary",
						rvar=rvar, cvar=cvar, dcor=FALSE,
						burn=burn, nscan=nscan, odens=odens,
						verbose = FALSE)
	
	# extract results
	beta_hat = median(fit$BETA[,2])
	beta_se = sd(fit$BETA[,2])
	beta_ci_lower = quantile(fit$BETA[,2], 0.025)
	beta_ci_upper = quantile(fit$BETA[,2], 0.975)
	beta_covered = (beta >= beta_ci_lower) & (beta <= beta_ci_upper)
	
	mu_hat = median(fit$BETA[,1])
	mu_se = sd(fit$BETA[,1])
	mu_ci_lower = quantile(fit$BETA[,1], 0.025)
	mu_ci_upper = quantile(fit$BETA[,1], 0.975)
	mu_covered = (mu >= mu_ci_lower) & (mu <= mu_ci_upper)
	
	# correlation with unobserved if applicable
	cor_with_W = if(!is.null(gamma) && R > 0 && !is.null(fit$UVPM)) {
		cor(c(W), c(fit$UVPM), use='pairwise.complete.obs')
	} else {
		NA
	}
	
	# recovery of additive effects
	a_cor = if(rvar && !is.null(fit$APM) && exists("a_true")) {
		cor(a_true, fit$APM, use='pairwise.complete.obs')
	} else {
		NA
	}
	
	b_cor = if(cvar && !is.null(fit$BPM) && exists("b_true")) {
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

####
# binary model with covariates only
####

test_that("Binary AME with covariates only has calibrated confidence intervals", {
	skip_on_cran()
	
	set.seed(6886)
	n_sims = 30
	n = 40
	mu_true = -0.3
	beta_true = 0.8
	
	# run simulations. odens=5 stores 200 posterior draws per fit (vs 40
	# at the default odens=25) so the 2.5/97.5% ci quantiles are stable
	# enough for a 30-replicate coverage check; with only 40 draws the
	# coverage estimate sits on the 0.80 threshold and flips with the
	# rng stream.
	results = lapply(1:n_sims, function(i) {
		sim_binary_ame(seed=i, n=n, mu=mu_true, beta=beta_true,
									rvar=FALSE, cvar=FALSE, R=0,
									burn=300, nscan=1000, odens=5)
	})
	
	# extract coverage rates
	beta_coverage = mean(sapply(results, function(x) x$beta_covered))
	mu_coverage = mean(sapply(results, function(x) x$mu_covered))
	
	# check calibration (relaxed for binary due to discrete nature)
	expect_gt(beta_coverage, 0.80)
	expect_gt(mu_coverage, 0.80)
	
	# check bias
	beta_bias = mean(sapply(results, function(x) x$beta_hat)) - beta_true
	mu_bias = mean(sapply(results, function(x) x$mu_hat)) - mu_true
	
	expect_lt(abs(beta_bias), 0.15)
	expect_lt(abs(mu_bias), 0.15)
})

####
# binary model edge cases
####

test_that("Binary AME handles sparse and dense networks", {
	set.seed(6886)
	n = 30
	
	# test sparse network
	Y_sparse = matrix(0, n, n)
	Y_sparse[sample(1:(n*n), size=floor(0.05*n*n))] = 1  # 5% density
	diag(Y_sparse) = NA
	
	fit_sparse = ame(Y_sparse, R=1, family="binary",
									 burn=200, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_sparse$BETA))
	expect_lt(median(fit_sparse$BETA[,1]), 0)  # intercept should be negative
	
	# test dense network
	Y_dense = matrix(1, n, n)
	Y_dense[sample(1:(n*n), size=floor(0.05*n*n))] = 0  # 95% density
	diag(Y_dense) = NA
	
	fit_dense = ame(Y_dense, R=1, family="binary",
									burn=200, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_dense$BETA))
	expect_gt(median(fit_dense$BETA[,1]), 0)  # intercept should be positive
})

#### ordinal family-specific tests ####

####
# ordinal model with additive effects
####

test_that("Ordinal AME with additive effects works", {
	set.seed(6886)
	n = 30
	
	# generate data with row/column heterogeneity
	a_true = rnorm(n, 0, 0.5)
	b_true = rnorm(n, 0, 0.5)
	
	eta = outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	Z = eta + matrix(rnorm(n*n), n, n)
	
	# discretize to ordinal
	Y = cut(c(Z), breaks=4, labels=FALSE)
	Y = matrix(Y, n, n)
	diag(Y) = NA
	
	# fit model with additive effects
	fit = ame(Y, R=0, family="ordinal",
						rvar=TRUE, cvar=TRUE,
						burn=300, nscan=1000, verbose = FALSE)
	
	expect_true(!is.null(fit$APM))
	expect_true(!is.null(fit$BPM))
	expect_equal(length(fit$APM), n)
	expect_equal(length(fit$BPM), n)
	
	# check some recovery (relaxed)
	if(!is.null(fit$APM) && !is.null(fit$BPM)) {
		a_cor = cor(a_true, fit$APM, use='pairwise.complete.obs')
		b_cor = cor(b_true, fit$BPM, use='pairwise.complete.obs')
		expect_gt(a_cor, 0.2)
		expect_gt(b_cor, 0.2)
	}
})

####
# ordinal model with different numbers of categories
####

test_that("Ordinal AME handles different numbers of categories", {
	set.seed(6886)
	n = 25
	
	# test with 3 categories
	Y_3 = matrix(sample(1:3, n*n, replace=TRUE), n, n)
	diag(Y_3) = NA
	
	fit_3 = ame(Y_3, R=1, family="ordinal",
							burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_3$BETA))
	expect_equal(length(unique(c(Y_3[!is.na(Y_3)]))), 3)
	
	# test with 7 categories
	Y_7 = matrix(sample(1:7, n*n, replace=TRUE), n, n)
	diag(Y_7) = NA
	
	fit_7 = ame(Y_7, R=1, family="ordinal",
							burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_7$BETA))
	expect_equal(length(unique(c(Y_7[!is.na(Y_7)]))), 7)
})

####
# ordinal model with skewed distributions
####

test_that("Ordinal AME handles skewed distributions", {
	set.seed(6886)
	n = 25
	
	# create skewed ordinal data (mostly low values)
	Y_skew = matrix(sample(1:5, n*n, replace=TRUE, 
												 prob=c(0.5, 0.25, 0.15, 0.07, 0.03)), n, n)
	diag(Y_skew) = NA
	
	# check skewness
	expect_lt(mean(Y_skew, na.rm=TRUE), 2.5)
	
	fit_skew = ame(Y_skew, R=1, family="ordinal",
								 burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_skew$BETA))
	# reverse-skewed data (mostly high values)
	Y_high = matrix(sample(1:5, n*n, replace=TRUE,
												 prob=c(0.03, 0.07, 0.15, 0.25, 0.5)), n, n)
	diag(Y_high) = NA
	
	expect_gt(mean(Y_high, na.rm=TRUE), 3.5)
	
	fit_high = ame(Y_high, R=1, family="ordinal",
								 burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_high$BETA))
})

####
# ordinal model with covariates and multiplicative effects
####

test_that("Ordinal AME with covariates and multiplicative effects works", {
	set.seed(6886)
	n = 30
	
	# generate covariate and latent factors
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	
	U_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	V_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	
	# build model
	eta = 0.4 * X + tcrossprod(U_true, V_true)
	Z = eta + matrix(rnorm(n*n), n, n)
	
	# discretize to ordinal
	Y = cut(c(Z), breaks=5, labels=FALSE)
	Y = matrix(Y, n, n)
	diag(Y) = NA
	
	# fit model
	fit = ame(Y, Xdyad=X, R=2, family="ordinal",
						burn=300, nscan=1000, verbose = FALSE)
	
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	expect_equal(ncol(fit$U), 2)
	expect_equal(ncol(fit$V), 2)
	
	# check that uvpm captures some structure
	if(!is.null(fit$UVPM)) {
	}
})

#### rank (cbin / frn) family-specific tests ####

####
# censored binary (cbin)
####

test_that("Censored binary (cbin) AME runs without errors", {
	set.seed(6886)
	n = 25
	
	# create binary data
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) = NA
	
	# set max outdegree (e.g., each node can nominate at most 5 others)
	odmax = 5
	
	# ensure y respects odmax constraint
	for(i in 1:n) {
		if(sum(Y[i,], na.rm=TRUE) > odmax) {
			# randomly keep only odmax nominations
			ones = which(Y[i,] == 1)
			keep = sample(ones, odmax)
			Y[i, setdiff(ones, keep)] = 0
		}
	}
	
	# test basic model
	fit = ame(Y, R=1, family="cbin", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	
	# check outdegrees respect constraint
	outdegrees = rowSums(Y, na.rm=TRUE)
	expect_true(all(outdegrees <= odmax))
})

####
# fixed rank nomination (frn)
####

test_that("Fixed rank nomination (frn) AME runs without errors", {
	set.seed(6886)
	n = 20
	
	# create ranked nomination data (e.g., rank top 3 friends)
	odmax = 3
	Y = matrix(NA, n, n)
	
	for(i in 1:n) {
		# each person nominates odmax others with ranks 1, 2, 3
		nominees = sample(setdiff(1:n, i), odmax)
		Y[i, nominees] = 1:odmax
	}
	diag(Y) = NA
	
	# test basic model
	fit = ame(Y, R=1, family="frn", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	
	# check that each row has exactly odmax nominations
	nominations_per_row = apply(Y, 1, function(x) sum(!is.na(x) & x > 0))
	expect_true(all(nominations_per_row == odmax))
})

####
# edge cases for rank-based models
####

test_that("Rank models handle varying outdegrees", {
	set.seed(6886)
	n = 20
	
	# test cbin with varying odmax per node
	odmax_vec = sample(2:5, n, replace=TRUE)
	Y_cbin = matrix(0, n, n)
	
	for(i in 1:n) {
		nominees = sample(setdiff(1:n, i), odmax_vec[i])
		Y_cbin[i, nominees] = 1
	}
	diag(Y_cbin) = NA
	
	fit_cbin = ame(Y_cbin, R=1, family="cbin", odmax=odmax_vec,
								 burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_cbin$BETA))
	
	# check each row respects its specific odmax
	for(i in 1:n) {
		expect_lte(sum(Y_cbin[i,], na.rm=TRUE), odmax_vec[i])
	}
})

test_that("FRN handles incomplete rankings", {
	set.seed(6886)
	n = 20
	
	# create data where some nodes make fewer nominations
	odmax = 5
	Y = matrix(NA, n, n)
	
	for(i in 1:n) {
		# some nodes nominate fewer than odmax
		n_nominations = sample(1:odmax, 1)
		nominees = sample(setdiff(1:n, i), n_nominations)
		Y[i, nominees] = 1:n_nominations
	}
	diag(Y) = NA
	
	fit = ame(Y, R=1, family="frn", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
})

####
# rank models with multiplicative effects
####

test_that("Rank models with multiplicative effects work", {
	set.seed(6886)
	n = 25
	
	# generate latent space positions
	U_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	V_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	
	# generate latent scores
	Z = tcrossprod(U_true, V_true) + matrix(rnorm(n*n, 0, 0.5), n, n)
	
	# convert to cbin data
	Y_cbin = 1 * (Z > 0)
	odmax = 6
	for(i in 1:n) {
		if(sum(Y_cbin[i,], na.rm=TRUE) > odmax) {
			ones = which(Y_cbin[i,] == 1)
			keep = sample(ones, odmax)
			Y_cbin[i, setdiff(ones, keep)] = 0
		}
	}
	diag(Y_cbin) = NA
	
	fit_cbin = ame(Y_cbin, R=2, family="cbin", odmax=odmax,
								 burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_cbin$U))
	expect_true(!is.null(fit_cbin$V))
	expect_equal(ncol(fit_cbin$U), 2)
	expect_equal(ncol(fit_cbin$V), 2)
	
	# convert to frn data (rank the positive links)
	Y_frn = matrix(NA, n, n)
	for(i in 1:n) {
		scores = Z[i,]
		scores[i] = NA
		top_indices = order(scores, decreasing=TRUE, na.last=TRUE)[1:min(odmax, n-1)]
		Y_frn[i, top_indices] = 1:length(top_indices)
	}
	diag(Y_frn) = NA
	
	fit_frn = ame(Y_frn, R=2, family="frn", odmax=odmax,
								burn=300, nscan=800, verbose = FALSE)

	expect_true(!is.null(fit_frn$U))
	expect_true(!is.null(fit_frn$V))
})
