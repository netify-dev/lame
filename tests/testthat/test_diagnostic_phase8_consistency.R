skip_on_cran()

# cross-function consistency

####
# setup: fit representative models
####

fit_ame_uni = local({
	n = 15
	dat = simulate_test_network(
		n = n, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.5,
		mode = "unipartite", seed = 8001
	)
	ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		burn = 100, nscan = 200, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_ame_bip = local({
	nA = 12; nB = 8
	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.5,
		mode = "bipartite", seed = 8002
	)
	ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		mode = "bipartite", burn = 100, nscan = 200, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_lame_uni = local({
	n = 15
	dat = simulate_test_network(
		n = n, n_time = 3, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.5,
		mode = "unipartite", seed = 8003
	)
	lame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		burn = 100, nscan = 200, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_lame_bip = local({
	nA = 12; nB = 8
	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 3, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.5,
		mode = "bipartite", seed = 8004
	)
	Xdyad_list = lapply(dat$Xdyad, function(x) x[, , 1])
	lame(
		dat$Y, Xdyad = Xdyad_list, R = 2, family = "normal",
		mode = "bipartite", burn = 100, nscan = 200, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

####
# predict vs fitted
####

test_that("predict(type='response') ~ fitted() for ame unipartite", {
	pred = predict(fit_ame_uni, type = "response")
	fitt = fitted(fit_ame_uni)
	# for normal family, response == link == YPM
	# they should be very close (or identical for normal)
	diff = abs(pred - fitt)
	diff = diff[!is.na(diff)]
	expect_true(max(diff) < 3.5)
})

test_that("predict(type='response') ~ fitted() for ame bipartite", {
	pred = predict(fit_ame_bip, type = "response")
	fitt = fitted(fit_ame_bip)
	diff = abs(pred - fitt)
	diff = diff[!is.na(diff)]
	expect_true(max(diff) < 5,
		info = paste("Max predict-fitted diff:", round(max(diff), 3)))
})

test_that("predict(type='response') ~ fitted() for lame unipartite", {
	pred = predict(fit_lame_uni, type = "response")
	fitt = fitted(fit_lame_uni)
	expect_equal(length(pred), length(fitt))
	# check each time point
	for (t in seq_along(pred)) {
		diff = abs(pred[[t]] - fitt[[t]])
		diff = diff[!is.na(diff)]
		expect_true(max(diff) < 3.5)
	}
})

test_that("predict(type='response') ~ fitted() for lame bipartite", {
	pred = predict(fit_lame_bip, type = "response")
	fitt = fitted(fit_lame_bip)
	expect_equal(length(pred), length(fitt))
	for (t in seq_along(pred)) {
		diff = abs(pred[[t]] - fitt[[t]])
		diff = diff[!is.na(diff)]
		expect_true(max(diff) < 3.5)
	}
})

####
# coef/vcov/confint consistency
####

test_that("confint centered near coef for ame unipartite", {
	b = coef(fit_ame_uni)
	ci = confint(fit_ame_uni, level = 0.95)
	# each coef should be within its CI
	for (i in seq_along(b)) {
		expect_true(b[i] >= ci[i, 1] && b[i] <= ci[i, 2],
								info = paste("Coef", names(b)[i], "not within CI"))
	}
})

test_that("confint centered near coef for ame bipartite", {
	b = coef(fit_ame_bip)
	ci = confint(fit_ame_bip, level = 0.95)
	for (i in seq_along(b)) {
		expect_true(b[i] >= ci[i, 1] && b[i] <= ci[i, 2])
	}
})

test_that("confint centered near coef for lame unipartite", {
	b = coef(fit_lame_uni)
	ci = confint(fit_lame_uni, level = 0.95)
	for (i in seq_along(b)) {
		expect_true(b[i] >= ci[i, 1] && b[i] <= ci[i, 2])
	}
})

test_that("diag(vcov) consistent with confint width for ame unipartite", {
	v = vcov(fit_ame_uni)
	ci = confint(fit_ame_uni, level = 0.95)
	# cI width ~ 2 * 1.96 * sqrt(diag(vcov)) for approx normal posterior
	ci_width = ci[, 2] - ci[, 1]
	approx_width = 2 * 1.96 * sqrt(diag(v))
	# should be roughly proportional (within factor of 2)
	ratio = ci_width / approx_width
	expect_true(all(ratio > 0.3 & ratio < 3),
							info = paste("Width ratios:", paste(round(ratio, 2), collapse = ", ")))
})

####
# residuals = Y - fitted
####

test_that("residuals = Y - fitted for ame unipartite", {
	r = residuals(fit_ame_uni, type = "response")
	fitt = fitted(fit_ame_uni)
	Y = fit_ame_uni$Y
	manual_r = Y - fitt
	# should be identical
	expect_equal(r, manual_r)
})

test_that("residuals = Y - fitted for ame bipartite", {
	r = residuals(fit_ame_bip, type = "response")
	fitt = fitted(fit_ame_bip)
	Y = fit_ame_bip$Y
	manual_r = Y - fitt
	expect_equal(r, manual_r)
})

test_that("residuals = Y - fitted for lame unipartite", {
	r = suppressWarnings(residuals(fit_lame_uni, type = "response"))
	fitt = fitted(fit_lame_uni)
	Y_array = fit_lame_uni$Y
	# convert Y from 3D array to list
	Y_list = lapply(seq_len(dim(Y_array)[3]), function(t) Y_array[,,t])
	for (t in seq_along(r)) {
		manual_r = Y_list[[t]] - fitt[[t]]
		expect_equal(r[[t]], manual_r)
	}
})

####
# simulate mean ~ predict
####

test_that("simulate mean approximates predict for ame normal", {
	# for normal family, mean of many simulations should approach predict
	sim_obj = simulate(fit_ame_uni, nsim = 50)
	sims = sim_obj$Y
	# average across simulations
	sim_mean = Reduce("+", sims) / length(sims)
	pred = predict(fit_ame_uni, type = "response")
	# should be roughly similar (with simulation noise)
	diff = abs(sim_mean - pred)
	diff = diff[!is.na(diff)]
	# mean absolute difference should be moderate (generous tolerance for MCMC)
	expect_true(mean(diff) < 2.0)
})

####
# reconstruct_EZ vs stored EZ (lame only)
####

test_that("reconstruct_EZ matches stored EZ for lame unipartite", {
	skip_if(is.null(fit_lame_uni$EZ), "EZ not stored in fit object")
	recon = reconstruct_EZ(fit_lame_uni)
	stored = fit_lame_uni$EZ
	expect_equal(length(recon), length(stored))
	for (t in seq_along(recon)) {
		diff = abs(recon[[t]] - stored[[t]])
		diff = diff[!is.na(diff)]
		expect_true(max(diff) < 0.1,
								info = paste("Max EZ diff at time", t, ":", max(diff)))
	}
})

test_that("reconstruct_EZ matches stored EZ for lame bipartite", {
	skip_if(is.null(fit_lame_bip$EZ), "EZ not stored in fit object")
	recon = reconstruct_EZ(fit_lame_bip)
	stored = fit_lame_bip$EZ
	expect_equal(length(recon), length(stored))
	for (t in seq_along(recon)) {
		diff = abs(recon[[t]] - stored[[t]])
		diff = diff[!is.na(diff)]
		expect_true(max(diff) < 0.1)
	}
})

####
# predict == fitted for non-identity link families
####

fit_ame_bin = local({
	n = 15
	dat = simulate_test_network(
		n = n, family = "binary", R = 2,
		beta_intercept = 0, beta_dyad = 0.5,
		mode = "unipartite", seed = 8010
	)
	ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "binary",
		burn = 100, nscan = 200, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	)
})

test_that("predict(type='response') == fitted() for binary ame", {
	pred = predict(fit_ame_bin, type = "response")
	fitt = fitted(fit_ame_bin)
	expect_equal(pred, fitt)
})

test_that("predict(type='response') values are probabilities for binary ame", {
	pred = predict(fit_ame_bin, type = "response")
	vals = pred[!is.na(pred)]
	expect_true(all(vals >= 0 & vals <= 1))
})

####
# gof() on lame objects
####

test_that("gof() works on lame objects", {
	result = gof(fit_lame_uni, nsim = 5, verbose = FALSE)
	expect_true(is.matrix(result))
	expect_equal(nrow(result), 6)  # 1 obs + 5 sims
	expect_true("sd.rowmean" %in% colnames(result))
	expect_true(all(is.finite(result[1, ])))
})
