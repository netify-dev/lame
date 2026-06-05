skip_on_cran()

# covariate reconstruction checks across families

test_that("Binary AME correctly handles covariates in all functions", {
	set.seed(6886)
	n = 15
	
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	
	# binary probit data with a dyadic covariate
	eta_true = -0.5 + 0.8 * X
	p_true = pnorm(eta_true)
	Y = matrix(NA, n, n)
	non_na_idx = which(!is.na(p_true))
	Y[non_na_idx] = rbinom(length(non_na_idx), 1, p_true[non_na_idx])
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X, R = 0, family = "binary",
						burn = 200, nscan = 500, verbose = FALSE, seed = 456)
	
	expect_false(is.null(fit$X), label = "X should be stored in fit object")
	expect_equal(dim(fit$X), c(n, n, 2))
	
	EZ = reconstruct_EZ(fit)
	expect_false(is.null(EZ))
	
	# reconstructed link tracks the data-generating link
	corr = cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
	expect_gt(corr, 0.5)
	
	pred_link = predict(fit, type = "link")
	expect_equal(pred_link, EZ, tolerance = 1e-10,
							info = "predict(type='link') should match reconstruct_EZ")
	
	pred_resp = predict(fit, type = "response")
	expect_true(all(pred_resp >= 0 & pred_resp <= 1, na.rm = TRUE),
							info = "Binary predictions should be probabilities")
	
	expected_probs = pnorm(EZ)
	mean_diff = mean(abs(fit$YPM - expected_probs), na.rm = TRUE)
	expect_lt(mean_diff, 0.1, label = "YPM should approximate Φ(EZ)")
})

test_that("Poisson AME correctly handles covariates in all functions", {
	set.seed(6886)
	n = 15
	
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	# poisson data with a dyadic covariate
	eta_true = 0.5 + 1.2 * X
	lambda_true = exp(eta_true)
	Y = matrix(NA, n, n)
	non_na_idx = which(!is.na(lambda_true))
	Y[non_na_idx] = rpois(length(non_na_idx), lambda_true[non_na_idx])
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X, R = 0, family = "poisson",
						burn = 200, nscan = 500, odens = 1, verbose = FALSE, seed = 567)
	
	expect_false(is.null(fit$X), label = "X should be stored in fit object")
	
	EZ = reconstruct_EZ(fit)
	expect_false(is.null(EZ))
	
	# reconstructed link tracks the data-generating link
	corr = cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
	expect_gt(corr, 0.5, label = "reconstruct_EZ should correlate with true log(lambda)")
	
	pred_link = predict(fit, type = "link")
	expect_equal(pred_link, EZ, tolerance = 1e-10)
	
	expected_lambda = exp(EZ)
	ratio = mean(fit$YPM / expected_lambda, na.rm = TRUE)
	expect_true(abs(ratio - 1) < 0.15, 
							info = paste("YPM/exp(EZ) ratio should be near 1, got", round(ratio, 3)))
	
	expect_true(all(fit$YPM >= 0, na.rm = TRUE),
							info = "Poisson YPM should be non-negative")
})

test_that("Gaussian AME correctly handles covariates in all functions", {
	set.seed(6886)
	n = 15
	
	X1 = matrix(rnorm(n*n), n, n)
	X2 = matrix(rnorm(n*n), n, n)
	diag(X1) = diag(X2) = NA
	X_array = array(c(X1, X2), dim = c(n, n, 2))
	
	# normal data with two dyadic covariates
	eta_true = 1.0 + 0.5 * X1 - 0.3 * X2
	Y = eta_true + matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X_array, R = 0, family = "normal",
						burn = 200, nscan = 500, verbose = FALSE, seed = 678)
	
	expect_false(is.null(fit$X), label = "X should be stored in fit object")
	expect_equal(dim(fit$X)[3], 3)
	
	EZ = reconstruct_EZ(fit)
	expect_false(is.null(EZ))
	
	# reconstructed link tracks the data-generating link
	corr = cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
	expect_gt(corr, 0.8, label = "reconstruct_EZ should highly correlate with true mean")
	
	mean_diff = mean(abs(fit$YPM - EZ), na.rm = TRUE)
	expect_lt(mean_diff, 0.15, label = "For normal family, YPM should equal EZ")
	
	residuals = Y - EZ
	expect_lt(abs(mean(residuals, na.rm = TRUE)), 0.1)
	
	beta_est = colMeans(fit$BETA)
	expect_equal(length(beta_est), 3)
	expect_true(abs(beta_est[1] - 1.0) < 0.5, label = "Intercept estimate")
	expect_true(abs(beta_est[2] - 0.5) < 0.5, label = "Beta1 estimate")
	expect_true(abs(beta_est[3] - (-0.3)) < 0.5, label = "Beta2 estimate")
})

test_that("Tobit AME correctly handles covariates in all functions", {
	set.seed(6886)
	n = 12
	
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	
	# censored normal data with a dyadic covariate
	eta_true = 0.5 + 0.8 * X
	Y_latent = eta_true + matrix(rnorm(n*n, 0, 0.5), n, n)
	Y = matrix(pmax(0, Y_latent), n, n)
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X, R = 0, family = "tobit",
						burn = 200, nscan = 500, verbose = FALSE, seed = 789)
	
	expect_false(is.null(fit$X), label = "X should be stored in fit object for tobit")
	
	EZ = reconstruct_EZ(fit)
	expect_false(is.null(EZ))
	
	expect_true(all(fit$YPM >= 0, na.rm = TRUE),
							info = "Tobit YPM should be non-negative")
	
	expect_true(any(EZ < 0, na.rm = TRUE) || min(EZ, na.rm = TRUE) < 0.1,
							info = "Tobit EZ (latent) can be negative")
})


test_that("Bipartite models correctly handle covariates across families", {
	set.seed(6886)
	nA = 8
	nB = 10
	
	Xdyad = array(rnorm(nA * nB * 2), c(nA, nB, 2))
	
	# binary bipartite reconstruction
	eta_true = 0.5 + 0.7 * Xdyad[,,1] - 0.4 * Xdyad[,,2]
	p_true = pnorm(eta_true)
	Y_bin = matrix(rbinom(nA * nB, 1, c(p_true)), nA, nB)
	
	fit_bin = ame(Y_bin, Xdyad = Xdyad, mode = "bipartite", family = "binary",
								burn = 200, nscan = 400, verbose = FALSE, seed = 890)
	
	expect_false(is.null(fit_bin$X), label = "Bipartite should store X")
	
	EZ_bin = reconstruct_EZ(fit_bin)
	expect_false(is.null(EZ_bin))
	expect_equal(dim(EZ_bin), c(nA, nB))
	
	corr = cor(c(eta_true), c(EZ_bin))
	expect_gt(corr, 0.4, label = "Bipartite binary EZ should correlate with truth")
	
	# normal bipartite reconstruction
	Y_norm = eta_true + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	
	fit_norm = ame(Y_norm, Xdyad = Xdyad, mode = "bipartite", family = "normal",
								 burn = 200, nscan = 400, verbose = FALSE, seed = 901)
	
	EZ_norm = reconstruct_EZ(fit_norm)
	corr_norm = cor(c(eta_true), c(EZ_norm))
	expect_gt(corr_norm, 0.7, label = "Bipartite normal EZ should highly correlate with truth")
})

test_that("predict() function works correctly with covariates", {
	set.seed(6886)
	n = 12
	
	X = matrix(rnorm(n*n), n, n)
	diag(X) = NA
	Y = 1 + 0.5 * X + matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X, R = 1, family = "normal",
						burn = 200, nscan = 400, verbose = FALSE)
	
	pred_link = predict(fit, type = "link")
	pred_resp = predict(fit, type = "response")

	expect_equal(pred_link, reconstruct_EZ(fit), tolerance = 1e-10)
	# response predictions track stored posterior means
	valid = !is.na(c(pred_resp)) & !is.na(c(fit$YPM))
	expect_gt(cor(c(pred_resp)[valid], c(fit$YPM)[valid]), 0.9)
	
	X_new = array(matrix(rnorm(n*n), n, n), dim = c(n, n, 1))
	diag(X_new[,,1]) = NA
	
	pred_new = predict(fit, newdata = X_new, type = "link")
	expect_false(is.null(pred_new))
	expect_equal(dim(pred_new), c(n, n))
	
	expect_true(mean(abs(pred_new - pred_link), na.rm = TRUE) > 0.1,
							info = "New predictions should differ from original")
})

test_that("simulate() function uses stored X correctly", {
	set.seed(6886)
	n = 10
	
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	Y = 2 + 1.5 * X + matrix(rnorm(n*n, 0, 0.3), n, n)
	diag(Y) = NA
	
	fit = ame(Y, Xdyad = X, R = 0, family = "normal",
						burn = 200, nscan = 400, verbose = FALSE)
	
	# simulate with stored covariates
	warnings_seen = new.env(parent = emptyenv())
	warnings_seen$value = character()
	withCallingHandlers({
		sims = simulate(fit, nsim = 10, seed = 123)
	}, warning = function(w) {
		warnings_seen$value = c(warnings_seen$value, conditionMessage(w))
		invokeRestart("muffleWarning")
	})
	
	expect_false(any(grepl("covariates not stored", warnings_seen$value)),
							 label = "Should not warn about missing X when X is stored")
	
	sim_mean = Reduce(`+`, sims$Y) / length(sims$Y)
	
	expect_true(abs(mean(sim_mean, na.rm = TRUE) - mean(Y, na.rm = TRUE)) < 1,
							info = "Simulated mean should be close to original")
})

test_that("Y is available in fit object for gof()", {
	set.seed(6886)
	n = 10
	Y = matrix(rnorm(n*n), n, n)
	diag(Y) = NA
	
	fit_uni = ame(Y, R = 0, family = "normal",
								burn = 100, nscan = 200, verbose = FALSE)
	
	expect_false(is.null(fit_uni$Y), label = "Unipartite should store Y")
	
	Y_bip = matrix(rnorm(8*10), 8, 10)
	fit_bip = ame(Y_bip, mode = "bipartite", family = "normal",
								burn = 100, nscan = 200, verbose = FALSE)
	
	expect_false(is.null(fit_bip$Y), label = "Bipartite should store Y")
	
	gof_result = gof(fit_uni, nsim = 5, verbose = FALSE)
	expect_false(is.null(gof_result))
	expect_true(nrow(gof_result) == 6)
})
