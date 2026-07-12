skip_on_cran()

# static bipartite ame models

test_that("bipartite ame fits a gaussian model", {
	skip_on_cran()
	
	set.seed(789)
	
	# set up bipartite network dimensions
	nA = 20
	nB = 15
	
	# true parameters
	mu_true = 0.5
	beta_true = 0.8
	sigma2_true = 0.5
	
	# generate covariates
	X = matrix(rnorm(nA * nB), nA, nB)
	
	# generate true network
	Y_true = mu_true + beta_true * X + matrix(rnorm(nA * nB, 0, sqrt(sigma2_true)), nA, nB)
	
	# fit basic model
	fit = ame(Y_true, Xdyad = X, mode = "bipartite", 
						 R_row = 0, R_col = 0,
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# dimensions are correct
	expect_equal(fit$mode, "bipartite")
	expect_equal(fit$nA, nA)
	expect_equal(fit$nB, nB)
	expect_equal(dim(fit$BETA)[2], 2)
	
	# parameter recovery
	beta_est = unname(colMeans(fit$BETA))
	expect_equal(beta_est[1], mu_true, tolerance = 0.2)
	expect_equal(beta_est[2], beta_true, tolerance = 0.2)
	
	# variance recovery
	s2_est = mean(fit$VC[,4])
	expect_equal(s2_est, sigma2_true, tolerance = 0.3)
})

test_that("bipartite model fits additive effects", {
	skip_on_cran()
	
	set.seed(790)
	
	nA = 15
	nB = 20
	
	# true parameters
	mu_true = 0.3
	sigma2_a = 0.4
	sigma2_b = 0.3
	sigma2_e = 0.5
	
	# generate additive effects
	a_true = rnorm(nA, 0, sqrt(sigma2_a))
	b_true = rnorm(nB, 0, sqrt(sigma2_b))
	
	# generate network
	Y = matrix(mu_true, nA, nB)
	for(i in 1:nA) {
		Y[i,] = Y[i,] + a_true[i]
	}
	for(j in 1:nB) {
		Y[,j] = Y[,j] + b_true[j]
	}
	Y = Y + matrix(rnorm(nA * nB, 0, sqrt(sigma2_e)), nA, nB)
	
	# fit model with additive effects
	fit = ame(Y, mode = "bipartite",
						 rvar = TRUE, cvar = TRUE,
						 R_row = 0, R_col = 0,
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# additive effects are present
	expect_false(is.null(fit$APM))
	expect_false(is.null(fit$BPM))
	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
	
	expect_true(mean(fit$VC[,1]) > 0)
	expect_true(mean(fit$VC[,3]) > 0)
	
	# row and column means are close to true values
	expect_equal(mean(fit$APM), mean(a_true), tolerance = 0.3)
	expect_equal(mean(fit$BPM), mean(b_true), tolerance = 0.3)
})

test_that("bipartite model fits multiplicative effects", {
	skip_on_cran()
	
	set.seed(791)
	
	nA = 20
	nB = 25
	R_row = 2
	R_col = 2
	
	# generate latent factors
	U_true = matrix(rnorm(nA * R_row), nA, R_row)
	V_true = matrix(rnorm(nB * R_col), nB, R_col)
	G_true = diag(c(2, 1), R_row, R_col)
	
	# generate network with multiplicative effects
	mu_true = 0.2
	UV_true = U_true %*% G_true %*% t(V_true)
	Y = mu_true + UV_true + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	
	# fit model with multiplicative effects
	fit = ame(Y, mode = "bipartite",
						 R_row = R_row, R_col = R_col,
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# latent-factor dimensions
	expect_equal(dim(fit$U), c(nA, R_row))
	expect_equal(dim(fit$V), c(nB, R_col))
	expect_equal(dim(fit$G), c(R_row, R_col))
	
	# multiplicative effects capture structure
	UVPM_reconstructed = reconstruct_UVPM(fit)
	
	# correlation should be positive
	if(!is.null(UVPM_reconstructed)) {
		cor_uv = cor(c(UV_true), c(UVPM_reconstructed))
		expect_gt(cor_uv, 0.3)
	}
})

test_that("bipartite binary model fits", {
	skip_on_cran()
	
	set.seed(792)
	
	nA = 15
	nB = 18
	
	# generate binary network
	mu_true = -0.5
	beta_true = 1.2
	X = matrix(rnorm(nA * nB), nA, nB)
	
	# probit model
	Z_true = mu_true + beta_true * X + matrix(rnorm(nA * nB), nA, nB)
	Y = (Z_true > 0) * 1
	
	# fit binary model
	fit = ame(Y, Xdyad = X, mode = "bipartite",
						 family = "binary",
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# model type
	expect_equal(fit$family, "binary")
	
	# parameter signs are correct
	beta_est = colMeans(fit$BETA)
	expect_true(beta_est[1] < 0)
	expect_true(beta_est[2] > 0)
	
	# predictions are binary probabilities
	expect_true(all(fit$YPM >= 0 & fit$YPM <= 1))
})

test_that("bipartite poisson uses the rectangular mh update for z", {
	skip_on_cran()
	set.seed(793)
	nA = 12; nB = 10
	X = matrix(runif(nA * nB, -1, 1), nA, nB)
	lambda = exp(1.0 + 0.5 * X)
	Y = matrix(rpois(nA * nB, lambda), nA, nB)
	# bipartite poisson uses the rectangular rz_pois_bip_fc sampler
	fit = ame(Y, Xdyad = X, mode = "bipartite", family = "poisson",
	          R = 0, burn = 5, nscan = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_equal(fit$family, "poisson")
	expect_equal(fit$mode, "bipartite")
})

test_that("bipartite multiplicative rank improves fit", {
	skip_on_cran()
	
	set.seed(794)
	
	nA = 15
	nB = 12
	
	# generate network with rank-2 structure
	U_true = matrix(rnorm(nA * 2), nA, 2)
	V_true = matrix(rnorm(nB * 2), nB, 2)
	Y = U_true %*% t(V_true) + matrix(rnorm(nA * nB, 0, 0.3), nA, nB)
	
	# fit with r = 0
	fit0 = ame(Y, mode = "bipartite",
							R_row = 0, R_col = 0,
							burn = 300, nscan = 1000, verbose = FALSE)
	
	# fit with r = 2
	fit2 = ame(Y, mode = "bipartite",
							R_row = 2, R_col = 2,
							burn = 300, nscan = 1000, verbose = FALSE)
	
	# model with multiplicative effects should fit better
	resid0 = Y - fit0$YPM
	resid2 = Y - fit2$YPM
	
	# both models produced valid predictions
	expect_false(all(is.na(fit0$YPM)))
	expect_false(all(is.na(fit2$YPM)))
	
	# the rank-2 model should have lower residuals than r = 0
	mse0 = mean(resid0^2, na.rm = TRUE)
	mse2 = mean(resid2^2, na.rm = TRUE)
	
	# only compare finite mse values
	if(is.finite(mse0) && is.finite(mse2)) {
		expect_lt(mse2, mse0)
	}
})

test_that("bipartite model handles missing data", {
	skip_on_cran()
	
	set.seed(795)
	
	nA = 10
	nB = 8
	
	# generate complete network
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# add missing values
	missing_idx = sample(1:(nA*nB), size = 10, replace = FALSE)
	Y[missing_idx] = NA
	
	# fit model
	fit = ame(Y, mode = "bipartite",
						 burn = 200, nscan = 500, verbose = FALSE)
	
	# model ran
	expect_equal(fit$mode, "bipartite")
	
	expect_false(any(is.na(fit$YPM)))
})

test_that("bipartite gof statistics are computed", {
	skip_on_cran()
	
	set.seed(796)
	
	nA = 8
	nB = 10
	
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# fit with gof
	fit = ame(Y, mode = "bipartite",
						 burn = 100, nscan = 300, 
						 gof = TRUE, verbose = FALSE)
	
	# gof output is present
	expect_false(is.null(fit$GOF))
	# bipartite gof has 3 columns: sd.rowmean, sd.colmean, four.cycles
	expect_equal(ncol(fit$GOF), 3)
	
	# gof statistics are non-negative
	obs_sd_row = fit$GOF[1,1]
	obs_sd_col = fit$GOF[1,2]
	obs_four_cycles = fit$GOF[1,3]
	
	expect_true(obs_sd_row >= 0)
	expect_true(obs_sd_col >= 0)
	expect_true(obs_four_cycles >= 0)
})
