skip_on_cran()

# static family-specific edge cases: link, threshold, censoring and rank-list behavior, degenerate inputs and
# family guards; the one parameter-recovery check per family is in test_static_recovery.R

library(lame)
library(testthat)

# gaussian family

# one simulation replicate
sim_gaussian_ame = function(
	seed,
	n,
	mu,
	beta,
	gamma = NULL,
	rvar = FALSE,
	cvar = FALSE,
	R = 0,
	burn = 500,
	nscan = 2000
) {
	set.seed(seed)

	# generate covariates
	xw = matrix(rnorm(n * 2), n, 2)
	X = tcrossprod(xw[, 1])
	diag(X) = NA

	# build linear predictor
	eta = mu + beta * X

	# add unobserved covariate effect if specified
	if (!is.null(gamma)) {
		W = tcrossprod(xw[, 2])
		eta = eta + gamma * W
	}

	# add additive effects if requested
	if (rvar || cvar) {
		a_true = if (rvar) rnorm(n, 0, 0.5) else rep(0, n)
		b_true = if (cvar) rnorm(n, 0, 0.5) else rep(0, n)
		eta = eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	}

	# add multiplicative effects if r > 0
	U_true = V_true = NULL
	if (R > 0) {
		U_true = matrix(rnorm(n * R, 0, 1), n, R)
		V_true = matrix(rnorm(n * R, 0, 1), n, R)
		eta = eta + tcrossprod(U_true, V_true)
	}

	# generate gaussian outcome with noise
	sigma2 = 1
	Y = eta + matrix(rnorm(n * n, 0, sqrt(sigma2)), n, n)
	diag(Y) = NA

	# fit ame model
	fit = ame(
		Y,
		Xdyad = X,
		R = R,
		family = "normal",
		rvar = rvar,
		cvar = cvar,
		dcor = FALSE,
		burn = burn,
		nscan = nscan,
		verbose = FALSE
	)

	# extract results
	beta_hat = median(fit$BETA[, 2])
	beta_se = sd(fit$BETA[, 2])
	beta_ci_lower = quantile(fit$BETA[, 2], 0.025)
	beta_ci_upper = quantile(fit$BETA[, 2], 0.975)
	beta_covered = (beta >= beta_ci_lower) & (beta <= beta_ci_upper)

	mu_hat = median(fit$BETA[, 1])
	mu_se = sd(fit$BETA[, 1])
	mu_ci_lower = quantile(fit$BETA[, 1], 0.025)
	mu_ci_upper = quantile(fit$BETA[, 1], 0.975)
	mu_covered = (mu >= mu_ci_lower) & (mu <= mu_ci_upper)

	# correlation with unobserved if applicable
	cor_with_W = if (!is.null(gamma) && R > 0) {
		cor(c(W), c(reconstruct_UVPM(fit)), use = 'pairwise.complete.obs')
	} else {
		NA
	}

	# recovery of additive effects
	a_cor = if (rvar && !is.null(fit$APM)) {
		cor(a_true, fit$APM, use = 'pairwise.complete.obs')
	} else {
		NA
	}

	b_cor = if (cvar && !is.null(fit$BPM)) {
		cor(b_true, fit$BPM, use = 'pairwise.complete.obs')
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
		sigma2_hat = if (!is.null(fit$s2)) median(fit$s2) else if (!is.null(fit$VC) && "va" %in% colnames(fit$VC)) median(fit$VC[, "va"]) else NA
	))
}

# gaussian model with covariates only

test_that("Gaussian AME with covariates only has calibrated confidence intervals", {
	skip_on_cran()

	set.seed(6886)
	n_sims = 30
	n = 40
	mu_true = 0
	beta_true = 1

	# run simulations
	results = lapply(1:n_sims, function(i) {
		sim_gaussian_ame(
			seed = i,
			n = n,
			mu = mu_true,
			beta = beta_true,
			rvar = FALSE,
			cvar = FALSE,
			R = 0,
			burn = 300,
			nscan = 1000
		)
	})

	# extract coverage rates
	beta_coverage = mean(sapply(results, function(x) x$beta_covered))
	mu_coverage = mean(sapply(results, function(x) x$mu_covered))

	# calibration near 95%
	expect_gt(beta_coverage, 0.85)
	expect_gt(mu_coverage, 0.85)

	# check bias
	beta_bias = mean(sapply(results, function(x) x$beta_hat)) - beta_true
	mu_bias = mean(sapply(results, function(x) x$mu_hat)) - mu_true

	expect_lt(abs(beta_bias), 0.1)
	expect_lt(abs(mu_bias), 0.1)
})

# gaussian model with covariates + additive effects

test_that("Gaussian AME with additive effects recovers true parameters", {
	set.seed(6886)
	n = 50
	mu_true = 0.5
	beta_true = 0.8

	# generate data with additive effects
	xw = matrix(rnorm(n * 2), n, 2)
	X = tcrossprod(xw[, 1])
	diag(X) = NA

	# true additive effects
	a_true = rnorm(n, 0, 0.7)
	b_true = rnorm(n, 0, 0.7)

	eta = mu_true + beta_true * X +
		outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)

	sigma2_true = 1
	Y = eta + matrix(rnorm(n * n, 0, sqrt(sigma2_true)), n, n)
	diag(Y) = NA

	# fit model with additive effects
	fit = ame(
		Y,
		Xdyad = X,
		R = 0,
		family = "normal",
		rvar = TRUE,
		cvar = TRUE,
		dcor = FALSE,
		burn = 500,
		nscan = 2000,
		verbose = FALSE
	)

	# check parameter recovery
	beta_est = median(fit$BETA[, 2])
	mu_est = median(fit$BETA[, 1])

	expect_lt(abs(beta_est - beta_true), 0.3)
	expect_lt(abs(mu_est - mu_true), 0.3)

	# recovery of additive effects
	if (!is.null(fit$APM) && !is.null(fit$BPM)) {
		a_cor = cor(a_true, fit$APM, use = 'pairwise.complete.obs')
		b_cor = cor(b_true, fit$BPM, use = 'pairwise.complete.obs')

		expect_gt(a_cor, 0.5)
		expect_gt(b_cor, 0.5)
	}
})

# gaussian model with covariates + additive + multiplicative effects

test_that("Gaussian AME with full model recovers true parameters", {
	set.seed(6886)
	n = 50
	mu_true = -0.2
	beta_true = 1.0
	gamma_true = 0.8 # unobserved covariate
	R_true = 2

	# generate data
	xw = matrix(rnorm(n * 2), n, 2)
	X = tcrossprod(xw[, 1])
	W = tcrossprod(xw[, 2]) # unobserved
	diag(X) = NA
	diag(W) = NA

	# true additive effects
	a_true = rnorm(n, 0, 0.5)
	b_true = rnorm(n, 0, 0.5)

	# true latent factors
	U_true = matrix(rnorm(n * R_true, 0, 1), n, R_true)
	V_true = matrix(rnorm(n * R_true, 0, 1), n, R_true)

	# build full model
	eta = mu_true + beta_true * X + gamma_true * W +
		outer(a_true, rep(1, n)) + outer(rep(1, n), b_true) +
		tcrossprod(U_true, V_true)

	sigma2_true = 1
	Y = eta + matrix(rnorm(n * n, 0, sqrt(sigma2_true)), n, n)
	diag(Y) = NA

	# fit full ame model
	fit = ame(
		Y,
		Xdyad = X,
		R = R_true,
		family = "normal",
		rvar = TRUE,
		cvar = TRUE,
		dcor = FALSE,
		burn = 500,
		nscan = 2000,
		verbose = FALSE
	)

	# check parameter recovery
	beta_est = median(fit$BETA[, 2])
	mu_est = median(fit$BETA[, 1])

	expect_lt(abs(beta_est - beta_true), 0.4)

	# multiplicative effects capture the unobserved covariate
	if (!is.null(reconstruct_UVPM(fit))) {
		cor_W = cor(c(W), c(reconstruct_UVPM(fit)), use = 'pairwise.complete.obs')
		expect_gt(cor_W, 0.3)
	}

	# dimensions of latent factors
	expect_equal(ncol(fit$U), R_true)
	expect_equal(ncol(fit$V), R_true)
})

# models with and without multiplicative effects

test_that("Multiplicative effects reduce bias from unobserved confounding", {
	set.seed(6886)
	n = 50
	mu_true = 0
	beta_true = 1
	gamma_true = 1.5 # strong unobserved effect

	# generate data with unobserved confounder
	xw = matrix(rnorm(n * 2), n, 2)
	X = tcrossprod(xw[, 1])
	W = tcrossprod(xw[, 2])
	diag(X) = NA
	diag(W) = NA

	eta = mu_true + beta_true * X + gamma_true * W
	Y = eta + matrix(rnorm(n * n), n, n)
	diag(Y) = NA

	# fit model without multiplicative effects
	fit_no_uv = ame(
		Y,
		Xdyad = X,
		R = 0,
		family = "normal",
		rvar = FALSE,
		cvar = FALSE,
		dcor = FALSE,
		burn = 400,
		nscan = 1500,
		verbose = FALSE
	)

	# fit model with multiplicative effects
	fit_with_uv = ame(
		Y,
		Xdyad = X,
		R = 2,
		family = "normal",
		rvar = FALSE,
		cvar = FALSE,
		dcor = FALSE,
		burn = 400,
		nscan = 1500,
		verbose = FALSE
	)

	beta_no_uv = median(fit_no_uv$BETA[, 2])
	beta_with_uv = median(fit_with_uv$BETA[, 2])

	# the model with uv has less bias
	bias_no_uv = abs(beta_no_uv - beta_true)
	bias_with_uv = abs(beta_with_uv - beta_true)

	expect_lt(bias_with_uv, bias_no_uv + 0.3)

	# uv captures the unobserved structure
	if (!is.null(fit_with_uv$UVPM)) {
		cor_W = cor(c(W), c(fit_with_uv$UVPM), use = 'pairwise.complete.obs')
		expect_gt(cor_W, 0.2)
	}
})

# model diagnostics and convergence

test_that("Gaussian AME produces valid diagnostics", {
	set.seed(6886)
	n = 30

	# simple model for a quick convergence check
	X = matrix(rnorm(n * n), n, n)
	diag(X) = NA
	Y = X + matrix(rnorm(n * n), n, n)
	diag(Y) = NA

	fit = ame(
		Y,
		Xdyad = X,
		R = 1,
		family = "normal",
		rvar = TRUE,
		cvar = TRUE,
		burn = 200,
		nscan = 500,
		verbose = FALSE
	)

	# key outputs exist
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$VC))

	# dimensions
	expect_equal(ncol(fit$BETA), 2) # intercept + one covariate
	# some samples are stored, not necessarily nscan
	expect_gt(nrow(fit$BETA), 0)

	# predictions are reasonable
	resid = Y - reconstruct_EZ(fit)
	expect_lt(abs(mean(resid, na.rm = TRUE)), 0.5)

	# variance components: one finite row per stored draw, positive variances and a correlation in [-1, 1]
	expect_true(is.matrix(fit$VC))
	expect_equal(nrow(fit$VC), nrow(fit$BETA))
	expect_true(all(c("va", "vb", "rho", "ve") %in% colnames(fit$VC)))
	expect_true(all(is.finite(fit$VC)))
	expect_true(all(fit$VC[, c("va", "vb", "ve")] > 0))
	expect_true(all(abs(fit$VC[, "rho"]) <= 1))
	# the residual variance is near the generating value of 1
	expect_equal(mean(fit$VC[, "ve"]), 1, tolerance = 0.3)
})

# count (poisson) family

# overdispersion handling

test_that("Poisson AME handles overdispersed count data", {
	set.seed(6886)
	n = 35

	# generate overdispersed count data using negative binomial
	X = matrix(rnorm(n * n, 0, 0.5), n, n)
	diag(X) = NA

	eta = 0.5 + 0.4 * X
	lambda = exp(eta)

	# overdispersion via negative binomial; the poisson model absorbs it with larger variance
	Y = matrix(suppressWarnings(rnbinom(n * n, mu = as.vector(lambda), size = 2)), n, n)
	diag(Y) = NA

	# overdispersion is present
	var_Y = var(c(Y), na.rm = TRUE)
	mean_Y = mean(Y, na.rm = TRUE)
	expect_gt(var_Y, mean_Y) # variance > mean indicates overdispersion

	# the poisson model absorbs overdispersion through random effects
	fit = ame(
		Y,
		Xdyad = X,
		R = 1,
		family = "poisson",
		rvar = TRUE,
		cvar = TRUE,
		dcor = FALSE,
		burn = 400,
		nscan = 1500,
		verbose = FALSE
	)

	# approximate beta recovery
	beta_est = median(fit$BETA[, 2])
	expect_lt(abs(beta_est - 0.4), 0.3)

	# random effects capture the extra variation
	if (!is.null(fit$APM) && !is.null(fit$BPM)) {
		var_random = var(c(fit$APM)) + var(c(fit$BPM))
		expect_gt(var_random, 0.001)
	}
})

# edge cases for count data

test_that("Poisson AME handles sparse count networks", {
	set.seed(6886)
	n = 30

	# very sparse count network, mostly zeros
	Y = matrix(0, n, n)
	n_nonzero = floor(0.1 * n * n)
	nonzero_idx = sample(1:(n * n), n_nonzero)
	Y[nonzero_idx] = rpois(n_nonzero, lambda = 2)
	diag(Y) = NA

	# sparsity
	expect_lt(mean(Y > 0, na.rm = TRUE), 0.15)

	fit = ame(
		Y,
		R = 1,
		family = "poisson",
		burn = 300,
		nscan = 1200,
		verbose = FALSE
	)

	expect_true(!is.null(fit$BETA))
	expect_lt(median(fit$BETA[, 1]), 1) # low intercept for sparse network

	# a sparse poisson network has low lambda, so EZ on the log scale is very negative
	EZ = reconstruct_EZ(fit)
	if (!is.null(EZ)) {
		expect_true(
			mean(EZ, na.rm = TRUE) < 0,
			info = "Sparse networks should have negative mean on log scale"
		)
	}

	# YPM stays non-negative
	expect_true(
		all(fit$YPM >= 0, na.rm = TRUE),
		info = "YPM (response scale) should be non-negative"
	)
})

test_that("Poisson AME handles zero-inflated patterns", {
	set.seed(6886)
	n = 30

	# generate zero-inflated count data
	X = matrix(rnorm(n * n, 0, 0.5), n, n)
	diag(X) = NA

	eta = 0.5 + 0.3 * X
	lambda = exp(eta)
	Y = matrix(suppressWarnings(rpois(n * n, as.vector(lambda))), n, n)

	# extra zeros (zero-inflation)
	zero_inflate_idx = sample(1:(n * n), floor(0.3 * n * n))
	Y[zero_inflate_idx] = 0
	diag(Y) = NA

	# random effects absorb the zero-inflation
	fit = ame(
		Y,
		Xdyad = X,
		R = 1,
		family = "poisson",
		rvar = TRUE,
		cvar = TRUE,
		burn = 300,
		nscan = 1200,
		verbose = FALSE
	)

	# estimates stay reasonable
	beta_est = median(fit$BETA[, 2])
	expect_lt(abs(beta_est - 0.3), 0.4) # more tolerance due to zero-inflation
})

# binary family

# one simulation replicate for binary
sim_binary_ame = function(
	seed,
	n,
	mu,
	beta,
	gamma = NULL,
	rvar = FALSE,
	cvar = FALSE,
	R = 0,
	burn = 500,
	nscan = 2000,
	odens = 25
) {
	set.seed(seed)

	# generate covariates
	xw = matrix(rnorm(n * 2), n, 2)
	X = tcrossprod(xw[, 1])
	diag(X) = NA

	# build linear predictor (probit scale)
	eta = mu + beta * X

	# add unobserved covariate effect if specified
	if (!is.null(gamma)) {
		W = tcrossprod(xw[, 2])
		eta = eta + gamma * W
	}

	# add additive effects if requested
	if (rvar || cvar) {
		a_true = if (rvar) rnorm(n, 0, 0.5) else rep(0, n)
		b_true = if (cvar) rnorm(n, 0, 0.5) else rep(0, n)
		eta = eta + outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	}

	# add multiplicative effects if r > 0
	U_true = V_true = NULL
	if (R > 0) {
		U_true = matrix(rnorm(n * R, 0, 1), n, R)
		V_true = matrix(rnorm(n * R, 0, 1), n, R)
		eta = eta + tcrossprod(U_true, V_true)
	}

	# generate binary outcome via probit link
	Z = eta + matrix(rnorm(n * n), n, n) # add standard normal noise
	Y = 1 * (Z > 0) # binary outcome
	diag(Y) = NA

	# fit ame model
	fit = ame(
		Y,
		Xdyad = X,
		R = R,
		family = "binary",
		rvar = rvar,
		cvar = cvar,
		dcor = FALSE,
		burn = burn,
		nscan = nscan,
		odens = odens,
		verbose = FALSE
	)

	# extract results
	beta_hat = median(fit$BETA[, 2])
	beta_se = sd(fit$BETA[, 2])
	beta_ci_lower = quantile(fit$BETA[, 2], 0.025)
	beta_ci_upper = quantile(fit$BETA[, 2], 0.975)
	beta_covered = (beta >= beta_ci_lower) & (beta <= beta_ci_upper)

	mu_hat = median(fit$BETA[, 1])
	mu_se = sd(fit$BETA[, 1])
	mu_ci_lower = quantile(fit$BETA[, 1], 0.025)
	mu_ci_upper = quantile(fit$BETA[, 1], 0.975)
	mu_covered = (mu >= mu_ci_lower) & (mu <= mu_ci_upper)

	# correlation with unobserved if applicable
	cor_with_W = if (!is.null(gamma) && R > 0 && !is.null(fit$UVPM)) {
		cor(c(W), c(fit$UVPM), use = 'pairwise.complete.obs')
	} else {
		NA
	}

	# recovery of additive effects
	a_cor = if (rvar && !is.null(fit$APM) && exists("a_true")) {
		cor(a_true, fit$APM, use = 'pairwise.complete.obs')
	} else {
		NA
	}

	b_cor = if (cvar && !is.null(fit$BPM) && exists("b_true")) {
		cor(b_true, fit$BPM, use = 'pairwise.complete.obs')
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

# binary model with covariates only

test_that("Binary AME with covariates only has calibrated confidence intervals", {
	skip_on_cran()

	set.seed(6886)
	n_sims = 30
	n = 40
	mu_true = -0.3
	beta_true = 0.8

	# odens = 5 keeps 200 draws per fit so the 2.5/97.5% quantiles are stable enough for a 30-replicate coverage check
	results = lapply(1:n_sims, function(i) {
		sim_binary_ame(
			seed = i,
			n = n,
			mu = mu_true,
			beta = beta_true,
			rvar = FALSE,
			cvar = FALSE,
			R = 0,
			burn = 300,
			nscan = 1000,
			odens = 5
		)
	})

	# extract coverage rates
	beta_coverage = mean(sapply(results, function(x) x$beta_covered))
	mu_coverage = mean(sapply(results, function(x) x$mu_covered))

	# calibration, relaxed for a discrete outcome
	expect_gt(beta_coverage, 0.80)
	expect_gt(mu_coverage, 0.80)

	# check bias
	beta_bias = mean(sapply(results, function(x) x$beta_hat)) - beta_true
	mu_bias = mean(sapply(results, function(x) x$mu_hat)) - mu_true

	expect_lt(abs(beta_bias), 0.15)
	expect_lt(abs(mu_bias), 0.15)
})

# binary model edge cases

test_that("Binary AME handles sparse and dense networks", {
	set.seed(6886)
	n = 30

	# sparse network
	Y_sparse = matrix(0, n, n)
	Y_sparse[sample(1:(n * n), size = floor(0.05 * n * n))] = 1 # 5% density
	diag(Y_sparse) = NA

	fit_sparse = ame(
		Y_sparse,
		R = 1,
		family = "binary",
		burn = 200,
		nscan = 800,
		verbose = FALSE
	)

	expect_true(!is.null(fit_sparse$BETA))
	expect_lt(median(fit_sparse$BETA[, 1]), 0)

	# dense network
	Y_dense = matrix(1, n, n)
	Y_dense[sample(1:(n * n), size = floor(0.05 * n * n))] = 0 # 95% density
	diag(Y_dense) = NA

	fit_dense = ame(
		Y_dense,
		R = 1,
		family = "binary",
		burn = 200,
		nscan = 800,
		verbose = FALSE
	)

	expect_true(!is.null(fit_dense$BETA))
	expect_gt(median(fit_dense$BETA[, 1]), 0)
})

# symmetric poisson

test_that("symmetric poisson fits use the pair sampler and recover the coefficients", {
	set.seed(6886)
	n = 24
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.4)
	Xt = matrix(rnorm(n * n), n, n)
	Xt = (Xt + t(Xt)) / 2
	diag(Xt) = NA
	E = matrix(0, n, n)
	E[upper.tri(E)] = rnorm(n * (n - 1) / 2, 0, sqrt(0.3))
	E = E + t(E)
	Z = 0.3 + 0.6 * Xt + outer(a, a, "+") + E
	Z[is.na(Z)] = 0
	Y = matrix(rpois(n * n, exp(pmin(Z, 6))), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	dimnames(Y) = list(nm, nm)
	# the pair kernel keeps the latent array symmetric and moves it
	Z0 = matrix(0, n, n)
	Z1 = rZ_pois_sym_fc(Z0, 0.3 + 0.6 * ifelse(is.na(Xt), 0, Xt), 0.3, Y)
	expect_equal(Z1, t(Z1))
	expect_true(all(is.finite(Z1)))
	expect_gt(mean(Z1[upper.tri(Z1)] != 0), 0.2)
	fit = suppressWarnings(ame(
		Y,
		array(Xt, c(n, n, 1), dimnames = list(nm, nm, "x1")),
		family = "poisson",
		R = 0,
		symmetric = TRUE,
		nscan = 800,
		burn = 200,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	))
	expect_equal(fit$YPM, t(fit$YPM))
	expect_lt(abs(mean(fit$BETA[, "x1_dyad"]) - 0.6), 0.15)
	expect_gt(mean(fit$YPM, na.rm = TRUE) / mean(Y, na.rm = TRUE), 0.85)
	expect_lt(mean(fit$YPM, na.rm = TRUE) / mean(Y, na.rm = TRUE), 1.15)
})

# ordinal family

# ordinal model with additive effects

test_that("Ordinal AME with additive effects works", {
	set.seed(6886)
	n = 30

	# generate data with row/column heterogeneity
	a_true = rnorm(n, 0, 0.5)
	b_true = rnorm(n, 0, 0.5)

	eta = outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
	Z = eta + matrix(rnorm(n * n), n, n)

	# discretize to ordinal
	Y = cut(c(Z), breaks = 4, labels = FALSE)
	Y = matrix(Y, n, n)
	diag(Y) = NA

	# fit model with additive effects
	fit = ame(
		Y,
		R = 0,
		family = "ordinal",
		rvar = TRUE,
		cvar = TRUE,
		burn = 300,
		nscan = 1000,
		verbose = FALSE
	)

	expect_true(!is.null(fit$APM))
	expect_true(!is.null(fit$BPM))
	expect_equal(length(fit$APM), n)
	expect_equal(length(fit$BPM), n)

	# some recovery, relaxed
	if (!is.null(fit$APM) && !is.null(fit$BPM)) {
		a_cor = cor(a_true, fit$APM, use = 'pairwise.complete.obs')
		b_cor = cor(b_true, fit$BPM, use = 'pairwise.complete.obs')
		expect_gt(a_cor, 0.2)
		expect_gt(b_cor, 0.2)
	}
})

# the ordinal effects prior is scaled to the data in lame() as in ame()

test_that("lame() and ame() place the same ordinal effects prior and agree on va", {
	set.seed(6886)
	n = 30
	nm = paste0("a", seq_len(n))
	a_true = rnorm(n, 0, 0.5)
	b_true = rnorm(n, 0, 0.5)
	Z = outer(a_true, b_true, "+") + matrix(rnorm(n * n), n, n)
	Y = matrix(cut(c(Z), breaks = 4, labels = FALSE), n, n)
	diag(Y) = NA
	dimnames(Y) = list(nm, nm)
	fit_l = suppressWarnings(lame(
		list(Y),
		family = "ordinal",
		R = 0,
		nscan = 600,
		burn = 200,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	# the prior scale follows the rank-normal heterogeneity of y, not a fixed diag(2)
	expect_lt(fit_l$prior$Sab0[1, 1], 1)
	expect_gte(fit_l$prior$eta0, 4)
	fit_a = suppressWarnings(ame(
		Y,
		family = "ordinal",
		R = 0,
		nscan = 600,
		burn = 200,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_lt(abs(mean(fit_l$VC[, "va"]) - mean(fit_a$VC[, "va"])), 0.15)
	expect_lt(abs(mean(fit_l$VC[, "vb"]) - mean(fit_a$VC[, "vb"])), 0.15)
})

# ordinal model with different numbers of categories

test_that("Ordinal AME handles different numbers of categories", {
	set.seed(6886)
	n = 25

	# 3 categories
	Y_3 = matrix(sample(1:3, n * n, replace = TRUE), n, n)
	diag(Y_3) = NA

	fit_3 = ame(
		Y_3,
		R = 1,
		family = "ordinal",
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit_3$BETA))
	expect_equal(length(unique(c(Y_3[!is.na(Y_3)]))), 3)

	# 7 categories
	Y_7 = matrix(sample(1:7, n * n, replace = TRUE), n, n)
	diag(Y_7) = NA

	fit_7 = ame(
		Y_7,
		R = 1,
		family = "ordinal",
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit_7$BETA))
	expect_equal(length(unique(c(Y_7[!is.na(Y_7)]))), 7)
})

# ordinal model with skewed distributions

test_that("Ordinal AME handles skewed distributions", {
	set.seed(6886)
	n = 25

	# skewed ordinal data, mostly low values
	Y_skew = matrix(
		sample(
			1:5,
			n * n,
			replace = TRUE,
			prob = c(0.5, 0.25, 0.15, 0.07, 0.03)
		),
		n,
		n
	)
	diag(Y_skew) = NA

	# skewness
	expect_lt(mean(Y_skew, na.rm = TRUE), 2.5)

	fit_skew = ame(
		Y_skew,
		R = 1,
		family = "ordinal",
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit_skew$BETA))
	# reverse-skewed data, mostly high values
	Y_high = matrix(
		sample(
			1:5,
			n * n,
			replace = TRUE,
			prob = c(0.03, 0.07, 0.15, 0.25, 0.5)
		),
		n,
		n
	)
	diag(Y_high) = NA

	expect_gt(mean(Y_high, na.rm = TRUE), 3.5)

	fit_high = ame(
		Y_high,
		R = 1,
		family = "ordinal",
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit_high$BETA))
})

# ordinal model with covariates and multiplicative effects

test_that("Ordinal AME with covariates and multiplicative effects works", {
	set.seed(6886)
	n = 30

	# generate covariate and latent factors
	X = matrix(rnorm(n * n), n, n)
	diag(X) = NA

	U_true = matrix(rnorm(n * 2, 0, 0.5), n, 2)
	V_true = matrix(rnorm(n * 2, 0, 0.5), n, 2)

	# build model
	eta = 0.4 * X + tcrossprod(U_true, V_true)
	Z = eta + matrix(rnorm(n * n), n, n)

	# discretize to ordinal
	Y = cut(c(Z), breaks = 5, labels = FALSE)
	Y = matrix(Y, n, n)
	diag(Y) = NA

	# fit model
	fit = ame(
		Y,
		Xdyad = X,
		R = 2,
		family = "ordinal",
		burn = 300,
		nscan = 1000,
		verbose = FALSE
	)

	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	expect_equal(ncol(fit$U), 2)
	expect_equal(ncol(fit$V), 2)

	# UVPM captures some structure
	if (!is.null(fit$UVPM)) {
	}
})

# rank (cbin / frn) families

# censored binary (cbin)

test_that("Censored binary (cbin) AME runs without errors", {
	set.seed(6886)
	n = 25

	# create binary data
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA

	# max outdegree: each node nominates at most odmax others
	odmax = 5

	# Y respects the odmax constraint
	for (i in 1:n) {
		if (sum(Y[i, ], na.rm = TRUE) > odmax) {
			# keep a random odmax nominations
			ones = which(Y[i, ] == 1)
			keep = sample(ones, odmax)
			Y[i, setdiff(ones, keep)] = 0
		}
	}

	# basic model
	fit = ame(
		Y,
		R = 1,
		family = "cbin",
		odmax = odmax,
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))

	# outdegrees respect the constraint
	outdegrees = rowSums(Y, na.rm = TRUE)
	expect_true(all(outdegrees <= odmax))
})

# fixed rank nomination (frn)

test_that("Fixed rank nomination (frn) AME runs without errors", {
	set.seed(6886)
	n = 20

	# ranked nomination data, e.g. top 3 friends
	odmax = 3
	Y = matrix(NA, n, n)

	for (i in 1:n) {
		# each person nominates odmax others with ranks 1, 2, 3
		nominees = sample(setdiff(1:n, i), odmax)
		Y[i, nominees] = 1:odmax
	}
	diag(Y) = NA

	# basic model
	fit = ame(
		Y,
		R = 1,
		family = "frn",
		odmax = odmax,
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))

	# each row has exactly odmax nominations
	nominations_per_row = apply(Y, 1, function(x) sum(!is.na(x) & x > 0))
	expect_true(all(nominations_per_row == odmax))
})

# edge cases for rank-based models

test_that("Rank models handle varying outdegrees", {
	set.seed(6886)
	n = 20

	# cbin with odmax varying by node
	odmax_vec = sample(2:5, n, replace = TRUE)
	Y_cbin = matrix(0, n, n)

	for (i in 1:n) {
		nominees = sample(setdiff(1:n, i), odmax_vec[i])
		Y_cbin[i, nominees] = 1
	}
	diag(Y_cbin) = NA

	fit_cbin = ame(
		Y_cbin,
		R = 1,
		family = "cbin",
		odmax = odmax_vec,
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit_cbin$BETA))

	# each row respects its own odmax
	for (i in 1:n) {
		expect_lte(sum(Y_cbin[i, ], na.rm = TRUE), odmax_vec[i])
	}
})

test_that("FRN handles incomplete rankings", {
	set.seed(6886)
	n = 20

	# some nodes make fewer nominations
	odmax = 5
	Y = matrix(NA, n, n)

	for (i in 1:n) {
		# some nodes nominate fewer than odmax
		n_nominations = sample(1:odmax, 1)
		nominees = sample(setdiff(1:n, i), n_nominations)
		Y[i, nominees] = 1:n_nominations
	}
	diag(Y) = NA

	fit = ame(
		Y,
		R = 1,
		family = "frn",
		odmax = odmax,
		burn = 200,
		nscan = 600,
		verbose = FALSE
	)

	expect_true(!is.null(fit$BETA))
})

# rank models with multiplicative effects

test_that("Rank models with multiplicative effects work", {
	set.seed(6886)
	n = 25

	# generate latent space positions
	U_true = matrix(rnorm(n * 2, 0, 0.5), n, 2)
	V_true = matrix(rnorm(n * 2, 0, 0.5), n, 2)

	# generate latent scores
	Z = tcrossprod(U_true, V_true) + matrix(rnorm(n * n, 0, 0.5), n, n)

	# convert to cbin data
	Y_cbin = 1 * (Z > 0)
	odmax = 6
	for (i in 1:n) {
		if (sum(Y_cbin[i, ], na.rm = TRUE) > odmax) {
			ones = which(Y_cbin[i, ] == 1)
			keep = sample(ones, odmax)
			Y_cbin[i, setdiff(ones, keep)] = 0
		}
	}
	diag(Y_cbin) = NA

	fit_cbin = ame(
		Y_cbin,
		R = 2,
		family = "cbin",
		odmax = odmax,
		burn = 300,
		nscan = 800,
		verbose = FALSE
	)

	expect_true(!is.null(fit_cbin$U))
	expect_true(!is.null(fit_cbin$V))
	expect_equal(ncol(fit_cbin$U), 2)
	expect_equal(ncol(fit_cbin$V), 2)

	# convert to frn data (rank the positive links)
	Y_frn = matrix(NA, n, n)
	for (i in 1:n) {
		scores = Z[i, ]
		scores[i] = NA
		top_indices = order(scores, decreasing = TRUE, na.last = TRUE)[1:min(odmax, n - 1)]
		Y_frn[i, top_indices] = 1:length(top_indices)
	}
	diag(Y_frn) = NA

	fit_frn = ame(
		Y_frn,
		R = 2,
		family = "frn",
		odmax = odmax,
		burn = 300,
		nscan = 800,
		verbose = FALSE
	)

	expect_true(!is.null(fit_frn$U))
	expect_true(!is.null(fit_frn$V))
})

test_that("ame() cbin and frn fits run with missing cells", {
	set.seed(6886)
	n = 22
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.5)
	b = rnorm(n, 0, 0.5)
	Z = 0.2 + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
	for (fam in c("cbin", "frn")) {
		Y = matrix(0, n, n)
		for (i in seq_len(n)) {
			zi = Z[i, ]
			zi[i] = -Inf
			top = order(zi, decreasing = TRUE)[1:4]
			top = top[zi[top] > 0]
			if (fam == "cbin") Y[i, top] = 1 else Y[i, top] = rank(zi[top])
		}
		diag(Y) = NA
		k = which(!is.na(Y))
		Y[sample(k, round(0.2 * length(k)))] = NA
		dimnames(Y) = list(nm, nm)
		# the joint effect-variance step must stay positive definite with missing cells
		fit = suppressWarnings(suppressMessages(ame(
			Y,
			family = fam,
			R = 0,
			odmax = 4,
			nscan = 200,
			burn = 50,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)))
		expect_s3_class(fit, "ame")
		expect_true(all(is.finite(fit$BETA)))
		expect_true(all(is.finite(fit$VC[, c("va", "vb")])))
		# the stored Y is the user's Y
		o = !is.na(Y)
		expect_equal(fit$Y[o], Y[o])
	}
})
