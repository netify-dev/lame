skip_on_cran()

# simulation-based parameter recovery tests for unipartite networks.
# these verify that ame() and lame() can recover known parameters
# from simulated unipartite data.

test_that("ame() unipartite normal recovers beta and additive effects", {
	skip_on_cran()

	n = 20
	R = 2
	true_intercept = 1.0
	true_beta = 0.8

	dat = simulate_test_network(
		n = n, n_time = 1,
		family = "normal", R = R,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "unipartite", seed = 123
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "normal",
		burn = 2000, nscan = 8000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat = colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.5,
		info = paste("Intercept:", beta_hat[1], "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.5,
		info = paste("Beta dyad:", beta_hat[2], "vs true:", true_beta))

	# additive effects correlation
	true_a = dat$true_params$a[[1]]
	true_b = dat$true_params$b[[1]]
	# actors may be sorted; align by name
	sorted_names = names(fit$APM)
	orig_names = rownames(dat$Y)
	names(true_a) = orig_names
	names(true_b) = orig_names
	cor_a = cor(fit$APM, true_a[sorted_names])
	cor_b = cor(fit$BPM, true_b[sorted_names])
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# UV product correlation
	true_UV = tcrossprod(dat$true_params$U[[1]], dat$true_params$V[[1]])
	est_UV = fit$U %*% t(fit$V)
	row_order = match(sorted_names, orig_names)
	true_UV_aligned = true_UV[row_order, row_order]
	cor_UV = cor(c(true_UV_aligned), c(est_UV))
	expect_true(cor_UV > 0.2,
		info = paste("UV product correlation:", round(cor_UV, 3)))

	# variance s2 should be near 1
	s2_hat = mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("ame() unipartite binary recovers beta and additive effects", {
	skip_on_cran()

	n = 25
	R = 2

	dat = simulate_test_network(
		n = n, n_time = 1,
		family = "binary", R = R,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "unipartite", seed = 456
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "binary",
		burn = 2000, nscan = 8000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat = colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation
	sorted_names = names(fit$APM)
	orig_names = rownames(dat$Y)
	true_a = dat$true_params$a[[1]]
	true_b = dat$true_params$b[[1]]
	names(true_a) = orig_names
	names(true_b) = orig_names
	cor_a = cor(fit$APM, true_a[sorted_names])
	cor_b = cor(fit$BPM, true_b[sorted_names])
	expect_true(cor_a > 0.2,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.2,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))
})


test_that("lame() unipartite normal recovers beta and additive effects", {
	skip_on_cran()

	n = 15
	n_time = 3
	true_intercept = 0.5
	true_beta = 0.6

	dat = simulate_test_network(
		n = n, n_time = n_time,
		family = "normal", R = 0,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "unipartite", seed = 789
	)

	Xdyad_list = lapply(dat$Xdyad, function(x) x[,,1])

	fit = lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "normal",
		burn = 1000, nscan = 5000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat = colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.6,
		info = paste("Intercept:", round(beta_hat[1], 3), "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.6,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "vs true:", true_beta))

	# additive effects correlation
	sorted_names = names(fit$APM)
	orig_names = rownames(dat$Y[[1]])
	true_a_avg = rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg = rowMeans(do.call(cbind, dat$true_params$b))
	names(true_a_avg) = orig_names
	names(true_b_avg) = orig_names
	cor_a = cor(fit$APM, true_a_avg[sorted_names])
	cor_b = cor(fit$BPM, true_b_avg[sorted_names])
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# s2 should be near 1
	s2_hat = mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1.5,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("lame() unipartite binary recovers beta", {
	skip_on_cran()

	n = 20
	n_time = 3

	dat = simulate_test_network(
		n = n, n_time = n_time,
		family = "binary", R = 0,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "unipartite", seed = 101
	)

	Xdyad_list = lapply(dat$Xdyad, function(x) x[,,1])

	fit = lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "binary",
		burn = 1000, nscan = 5000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat = colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation
	sorted_names = names(fit$APM)
	orig_names = rownames(dat$Y[[1]])
	true_a_avg = rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg = rowMeans(do.call(cbind, dat$true_params$b))
	names(true_a_avg) = orig_names
	names(true_b_avg) = orig_names
	cor_a = cor(fit$APM, true_a_avg[sorted_names])
	cor_b = cor(fit$BPM, true_b_avg[sorted_names])
	expect_true(cor_a > 0.2,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.2,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))
})
