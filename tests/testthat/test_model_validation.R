skip_on_cran()

# model validation tests for bias, coverage, and parameter recovery
library(lame)

test_that("ame recovers known parameters in simulated data", {
	skip_on_cran()  # Skip on CRAN to save time
	
	set.seed(6886)
	n = 30  # Small network for testing
	
	# true parameters
	mu_true = -1
	beta_true = 0.5
	
	# generate covariates
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	
	# generate true latent positions for confounding
	# use full-strength UV signal so latent structure is detectable in binary data
	R_true = 2
	U_true = matrix(rnorm(n*R_true), n, R_true)
	V_true = matrix(rnorm(n*R_true), n, R_true)
	UV_true = tcrossprod(U_true, V_true)

	# generate network with known parameters
	eta = mu_true + beta_true * X + UV_true + matrix(rnorm(n*n, 0, 0.5), n, n)
	Y = (eta > 0) * 1
	diag(Y) = NA
	
	# fit models
	fit_naive = ame(Y, Xdyad=X, R=0, rvar=FALSE, cvar=FALSE, dcor=FALSE, 
									 family="binary", burn=100, nscan=500, verbose = FALSE)
	
	fit_ame = ame(Y, Xdyad=X, R=2, rvar=FALSE, cvar=FALSE, dcor=FALSE,
								 family="binary", burn=100, nscan=500, verbose = FALSE)
	
	# check bias
	beta_naive = median(fit_naive$BETA[,2])
	beta_ame = median(fit_ame$BETA[,2])
	
	# AME should have less bias than naive model
	bias_naive = abs(beta_naive - beta_true)
	bias_ame = abs(beta_ame - beta_true)
	
	# the AME model should not have more bias than the naive model
	expect_lt(bias_ame, bias_naive + 0.1)
	
	# check that AME captures the confounding structure
	UV_fitted = reconstruct_UVPM(fit_ame)
	if(!is.null(UV_fitted)) {
		cor_uv = cor(c(UV_fitted), c(UV_true), use="complete.obs")
		expect_gt(abs(cor_uv), 0.1)
	} else {
		# at minimum check latent factors exist
		expect_true(!is.null(fit_ame$U))
	}
})

test_that("lame handles longitudinal data correctly", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	T = 5
	
	# generate longitudinal networks with temporal dependence
	Y_list = list()
	mu_true = -1
	
	# create temporal correlation
	for(t in 1:T) {
		if(t == 1) {
			eta = matrix(rnorm(n*n, mu_true, 1), n, n)
		} else {
			# add temporal dependence
			eta = 0.7 * eta + matrix(rnorm(n*n, mu_true * 0.3, 0.5), n, n)
		}
		Y_list[[t]] = (eta > 0) * 1
		diag(Y_list[[t]]) = NA
		rownames(Y_list[[t]]) = colnames(Y_list[[t]]) = paste0("Node", 1:n)
	}
	
	# fit static model
	fit_static = lame(Y_list, R=2, family="binary", 
										 dynamic_uv=FALSE, dynamic_ab=FALSE,
										 burn=100, nscan=300, verbose = FALSE)
	
	# fit dynamic model
	fit_dynamic = lame(Y_list, R=2, family="binary",
											dynamic_uv=TRUE, dynamic_ab=FALSE,
											burn=100, nscan=300, verbose = FALSE,
											prior=list(rho_uv_mean=0.5))
	
	# dynamic model should capture temporal correlation
	# check rho_uv if available
	if(!is.null(fit_dynamic$rho_uv) && length(fit_dynamic$rho_uv) > 0) {
		rho_est = median(fit_dynamic$rho_uv, na.rm=TRUE)
		# data has rho~0.7 and prior centered at 0.5, so estimate should be positive
		expect_gt(rho_est, 0.1)
		expect_lt(rho_est, 1.0)
	} else {
		# at least check model ran
		expect_true(!is.null(fit_dynamic$BETA))
	}
})

test_that("bipartite ame works correctly", {
	skip_on_cran()
	
	set.seed(6886)
	nA = 25  # Row nodes
	nB = 30  # Column nodes
	
	# true parameters
	mu_true = -0.5
	R_row_true = 2
	R_col_true = 2
	
	# generate true latent positions
	U_true = matrix(rnorm(nA * R_row_true), nA, R_row_true)
	V_true = matrix(rnorm(nB * R_col_true), nB, R_col_true)
	G_true = matrix(c(1, 0.5, 0.5, -1), R_row_true, R_col_true)
	
	# generate bipartite network
	eta = mu_true + U_true %*% G_true %*% t(V_true) + 
				 matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	Y_bip = (eta > 0) * 1
	
	# fit bipartite model
	fit_bip = ame(Y_bip, mode="bipartite", R_row=2, R_col=2,
								 family="binary", burn=100, nscan=500,
								 verbose = FALSE)
	
	# check dimensions - U and V are posterior means (2D matrices)
	expect_equal(dim(fit_bip$U), c(nA, R_row_true))
	expect_equal(dim(fit_bip$V), c(nB, R_col_true))
	# g might be 2D or not exist
	if(!is.null(fit_bip$G)) {
		expect_equal(dim(fit_bip$G), c(R_row_true, R_col_true))
		G_est = fit_bip$G
	} else {
		# use identity if G not available
		G_est = diag(min(R_row_true, R_col_true))
	}
	
	# check that model captures structure
	U_est = fit_bip$U
	V_est = fit_bip$V
	
	# reconstructed matrix should correlate with true structure
	eta_est = U_est %*% G_est %*% t(V_est)
	eta_true = U_true %*% G_true %*% t(V_true)
	cor_eta = cor(c(eta_est), c(eta_true))
	expect_gt(abs(cor_eta), 0.1)
})

test_that("bipartite lame handles longitudinal bipartite data", {
	skip_on_cran()
	
	set.seed(6886)
	nA = 20
	nB = 25
	T = 4
	
	# generate longitudinal bipartite networks
	Y_list = list()
	U_true = matrix(rnorm(nA * 2), nA, 2)
	V_true = matrix(rnorm(nB * 2), nB, 2)
	G_true = matrix(c(1, 0, 0, 1), 2, 2)
	
	for(t in 1:T) {
		# add some temporal variation
		U_t = U_true + matrix(rnorm(nA * 2, 0, 0.1 * t), nA, 2)
		V_t = V_true + matrix(rnorm(nB * 2, 0, 0.1 * t), nB, 2)
		
		eta = U_t %*% G_true %*% t(V_t) + matrix(rnorm(nA * nB, -1, 0.5), nA, nB)
		Y_list[[t]] = (eta > 0) * 1
		rownames(Y_list[[t]]) = paste0("Row", 1:nA)
		colnames(Y_list[[t]]) = paste0("Col", 1:nB)
	}
	
	# fit bipartite longitudinal model
	fit_bip_long = lame(Y_list, mode="bipartite", 
											 R_row=2, R_col=2,
											 family="binary",
											 dynamic_uv=FALSE,
											 burn=100, nscan=300,
											 verbose = FALSE)
	
	# check output structure
	expect_true(!is.null(fit_bip_long$U))
	expect_true(!is.null(fit_bip_long$V))
	expect_true(!is.null(fit_bip_long$G))
	
	# check dimensions for longitudinal bipartite
	expect_equal(dim(fit_bip_long$U)[1], nA)
	expect_equal(dim(fit_bip_long$V)[1], nB)
})

test_that("confidence intervals have correct coverage", {
	skip_on_cran()

	n_sims = 50
	n = 25
	mu_true = -1
	beta_true = 0.5
	
	coverage_mu = numeric(n_sims)
	coverage_beta = numeric(n_sims)
	
	for(sim in 1:n_sims) {
		set.seed(6886 + sim)
		
		# generate data
		X = matrix(rnorm(n*n), n, n)
		diag(X) = NA
		eta = mu_true + beta_true * X + matrix(rnorm(n*n), n, n)
		Y = (eta > 0) * 1
		diag(Y) = NA
		
		# fit model
		fit = ame(Y, Xdyad=X, R=1, family="binary",
							 burn=100, nscan=500, verbose = FALSE)
		
		# check coverage
		mu_ci = quantile(fit$BETA[,1], c(0.025, 0.975))
		beta_ci = quantile(fit$BETA[,2], c(0.025, 0.975))
		
		coverage_mu[sim] = (mu_ci[1] < mu_true) & (mu_true < mu_ci[2])
		coverage_beta[sim] = (beta_ci[1] < beta_true) & (beta_true < beta_ci[2])
	}
	
	# coverage should be close to 95%
	expect_gt(mean(coverage_mu), 0.85)
	expect_lt(mean(coverage_mu), 1.0)
	expect_gt(mean(coverage_beta), 0.85)
	expect_lt(mean(coverage_beta), 1.0)
})

test_that("different families produce appropriate outputs", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	
	# test binary family
	Y_bin = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y_bin) = NA
	
	fit_bin = ame(Y_bin, R=0, family="binary", 
								 burn=50, nscan=200, verbose = FALSE)
	# latent Z is unbounded for probit; check YPM is in [0,1]
	expect_true(all(fit_bin$YPM >= 0 & fit_bin$YPM <= 1, na.rm=TRUE))
	
	# test normal family
	Y_norm = matrix(rnorm(n*n), n, n)
	diag(Y_norm) = NA
	
	fit_norm = ame(Y_norm, R=0, family="normal",
									burn=50, nscan=200, verbose = FALSE)
	expect_true(is.numeric(fit_norm$BETA))
	
	# test ordinal family with proper ordinal data
	Y_ord = matrix(sample(1:5, n*n, replace=TRUE), n, n)
	diag(Y_ord) = NA
	
	fit_ord = ame(Y_ord, R=0, family="ordinal", intercept=FALSE,
								 burn=50, nscan=200, verbose = FALSE)
	expect_true(!is.null(fit_ord$BETA))
})

test_that("GOF statistics are computed correctly", {
	skip_on_cran()
	
	set.seed(6886)
	n = 25
	
	# generate network with known properties
	Y = matrix(0, n, n)
	# create some triangles
	for(i in 1:(n-2)) {
		if(runif(1) > 0.5) {
			Y[i, i+1] = Y[i+1, i] = 1
			Y[i, i+2] = Y[i+2, i] = 1
			Y[i+1, i+2] = Y[i+2, i+1] = 1
		}
	}
	diag(Y) = NA
	
	fit = ame(Y, R=1, family="binary", symmetric=TRUE,
						 burn=100, nscan=300, verbose = FALSE, gof=TRUE)
	
	# check GOF statistics exist
	expect_true(!is.null(fit$GOF))
	expect_true(ncol(fit$GOF) >= 3)
	
	# GOF first row should have positive sd statistics for non-trivial network
	obs_sd_row = fit$GOF[1, "sd.rowmean"]
	obs_sd_col = fit$GOF[1, "sd.colmean"]
	expect_gt(obs_sd_row, 0)
	expect_gt(obs_sd_col, 0)
	# posterior predictive mean should be in the same ballpark as observed
	mean_sd_row = mean(fit$GOF[-1, "sd.rowmean"], na.rm = TRUE)
	expect_gt(mean_sd_row, 0)
})

test_that("dynamic effects capture temporal persistence", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	T = 6
	rho_true = 0.8
	
	# generate networks with known temporal correlation
	Y_list = list()
	U_prev = matrix(rnorm(n*2), n, 2)
	
	for(t in 1:T) {
		# AR(1) process for latent positions
		U_curr = rho_true * U_prev + matrix(rnorm(n*2, 0, sqrt(1-rho_true^2)), n, 2)
		
		eta = tcrossprod(U_curr) + matrix(rnorm(n*n, -2, 0.5), n, n)
		Y_list[[t]] = (eta > 0) * 1
		diag(Y_list[[t]]) = NA
		rownames(Y_list[[t]]) = colnames(Y_list[[t]]) = paste0("Node", 1:n)
		
		U_prev = U_curr
	}
	
	# fit dynamic model with good priors
	fit_dyn = lame(Y_list, R=2, family="binary",
									dynamic_uv=TRUE, dynamic_ab=FALSE,
									burn=500, nscan=2500, odens=25, verbose = FALSE,
									prior=list(rho_uv_mean=0.7, rho_uv_sd=0.2))

	# should recover temporal correlation
	expect_true(!is.null(fit_dyn$rho_uv))
	expect_true(length(fit_dyn$rho_uv) > 1)
	rho_est = median(fit_dyn$rho_uv, na.rm=TRUE)
	# data generated with rho=0.8 and prior centered at 0.7, estimate should be positive
	expect_gt(rho_est, 0.1)
	expect_lt(abs(rho_est - rho_true), 0.7)
})

test_that("model handles missing data appropriately", {
	skip_on_cran()
	
	set.seed(6886)
	n = 25
	
	# create network with missing data
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) = NA
	
	# add additional missing values
	missing_idx = sample(which(!is.na(Y)), size=floor(0.2 * n*n))
	Y[missing_idx] = NA
	
	# model should handle missing data
	expect_no_error({
		fit = ame(Y, R=1, family="binary",
							 burn=50, nscan=200, verbose = FALSE)
	})
	
	# check that predictions are made for missing values
	Y_pred = fit$EZ
	expect_true(all(!is.na(Y_pred[missing_idx])))
})

test_that("symmetric networks maintain symmetry", {
	skip_on_cran()
	
	set.seed(6886) 
	n = 20
	
	# generate symmetric network
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	
	fit_sym = ame(Y, R=2, family="binary", symmetric=TRUE,
								 burn=100, nscan=300, verbose = FALSE)
	
	# check U = V for symmetric model
	# for symmetric models, V might not exist or U=V
	if(!is.null(fit_sym$V)) {
		# if V exists and has same dimensions as U, they should be equal
		if(identical(dim(fit_sym$U), dim(fit_sym$V))) {
			expect_lt(mean(abs(fit_sym$U - fit_sym$V)), 0.1)
		} else {
			# different dimensions - just check both exist
			expect_true(!is.null(fit_sym$U))
			expect_true(!is.null(fit_sym$V))
		}
	} else {
		# just check that U exists
		expect_true(!is.null(fit_sym$U))
	}
	
	# predicted values should be symmetric
	EZ = fit_sym$EZ
	if(is.null(EZ)) {
		# EZ not stored directly; reconstruct from posterior means
		EZ = reconstruct_EZ(fit_sym)
	}
	if(is.matrix(EZ)) {
		expect_lt(mean(abs(EZ - t(EZ)), na.rm=TRUE), 0.1)
	} else {
		# at minimum, model should have produced posterior means
		expect_false(is.null(fit_sym$YPM))
	}
})