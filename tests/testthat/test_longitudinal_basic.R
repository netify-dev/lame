skip_on_cran()

# longitudinal models via lame()

library(lame)
library(testthat)

test_that("lame() runs with normal family and returns correct structure", {
	skip_on_cran()

	dat = simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 1, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, rho = 0, s2 = 1, seed = 101
	)

	fit = lame(
		Y = dat$Y,
		Xdyad = dat$Xdyad,
		Xrow = dat$Xrow,
		Xcol = dat$Xcol,
		family = "normal",
		R = 0,
		burn = 100,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = TRUE
	)

	# structure checks
	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
	expect_false(is.null(fit$VC))
	expect_false(is.null(fit$APM))
	expect_false(is.null(fit$BPM))
	expect_false(is.null(fit$YPM))

	# BETA has intercept, Xdyad, Xrow and Xcol coefficients: 4 columns
	expect_equal(ncol(fit$BETA), 4)
	# 100 posterior samples (500 / odens=5)
	expect_equal(nrow(fit$BETA), 100)

	# the intercept is in a reasonable range
	beta_intercept_hat = median(fit$BETA[, 1])
	expect_true(
		abs(beta_intercept_hat - 1) < 2,
		info = paste("Intercept estimate:", beta_intercept_hat, "expected ~1")
	)

	# dyad covariate recovery
	beta_dyad_hat = median(fit$BETA[, 2])
	expect_true(
		abs(beta_dyad_hat - 0.5) < 1,
		info = paste("Dyad beta estimate:", beta_dyad_hat, "expected ~0.5")
	)
})

test_that("lame() with covariates recovers known effects", {
	skip_on_cran()

	set.seed(200)
	n = 20
	n_time = 4
	true_beta = 0.8

	# data with a strong known effect
	actors = paste0("a", 1:n)
	Y_list = list()
	Xdyad_list = list()

	for (t in 1:n_time) {
		X = matrix(rnorm(n * n), n, n)
		diag(X) = NA
		rownames(X) = colnames(X) = actors
		Z = true_beta * X + matrix(rnorm(n * n), n, n)
		diag(Z) = NA

		Y_list[[t]] = Z
		Xdyad_list[[t]] = X
	}
	names(Y_list) = names(Xdyad_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		Xdyad = Xdyad_list,
		family = "normal",
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		burn = 200,
		nscan = 1000,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)

	# beta recovery
	beta_hat = median(fit$BETA[, 2])
	beta_ci = quantile(fit$BETA[, 2], c(0.025, 0.975))

	expect_true(
		true_beta >= beta_ci[1] && true_beta <= beta_ci[2],
		info = paste(
			"True beta",
			true_beta,
			"should be in CI [",
			round(beta_ci[1], 3),
			",",
			round(beta_ci[2], 3),
			"]"
		)
	)

	# the point estimate is close
	expect_true(
		abs(beta_hat - true_beta) < 0.3,
		info = paste("Beta hat:", round(beta_hat, 3), "expected:", true_beta)
	)
})

test_that("lame() additive effects correlate with true values", {
	skip_on_cran()

	set.seed(301)
	n = 20
	n_time = 5
	actors = paste0("a", 1:n)

	# strong additive effects with large variance relative to the noise
	true_a = rnorm(n, 0, 2)
	true_b = rnorm(n, 0, 2)

	Y_list = list()
	for (t in 1:n_time) {
		eta = outer(true_a, rep(1, n)) + outer(rep(1, n), true_b)
		Y = eta + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal",
		R = 0,
		rvar = TRUE,
		cvar = TRUE,
		burn = 500,
		nscan = 2000,
		odens = 10,
		verbose = FALSE,
		gof = FALSE
	)

	# additive effects correlate with the truth
	if (!is.null(fit$APM) && length(fit$APM) == n) {
		# match APM to true_a by actor name (lame sorts actors)
		names(true_a) = actors
		apm_matched = fit$APM[actors]
		a_centered = true_a - mean(true_a)
		apm_centered = apm_matched - mean(apm_matched)
		a_cor = cor(a_centered, apm_centered, use = "pairwise.complete.obs")
		expect_true(
			a_cor > 0.3,
			info = paste("Sender effect correlation:", round(a_cor, 3))
		)
	}
	if (!is.null(fit$BPM) && length(fit$BPM) == n) {
		names(true_b) = actors
		bpm_matched = fit$BPM[actors]
		b_centered = true_b - mean(true_b)
		bpm_centered = bpm_matched - mean(bpm_matched)
		b_cor = cor(b_centered, bpm_centered, use = "pairwise.complete.obs")
		expect_true(
			b_cor > 0.3,
			info = paste("Receiver effect correlation:", round(b_cor, 3))
		)
	}
})

test_that("lame() GOF output has correct structure", {
	skip_on_cran()

	dat = simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 0, beta_dyad = 0,
		seed = 105
	)

	fit = lame(
		Y = dat$Y,
		family = "normal",
		R = 0,
		burn = 50,
		nscan = 250,
		odens = 5,
		verbose = FALSE,
		gof = TRUE
	)

	expect_false(
		is.null(fit$GOF),
		info = "GOF should be computed when gof=TRUE"
	)
	if (!is.null(fit$GOF)) {
		# GOF of a lame fit is a list of matrices, one per statistic
		expect_true(is.list(fit$GOF) || is.matrix(fit$GOF))
		if (is.list(fit$GOF)) {
			expect_true(length(fit$GOF) > 0, info = "GOF should have at least one statistic")
		}
	}
})

# further families

test_that("lame() works with ordinal family", {
	skip_on_cran()

	dat = simulate_test_network(
		n = 15, n_time = 3, family = "ordinal",
		R = 0, beta_intercept = 0, beta_dyad = 0.3, seed = 503
	)

	fit = lame(
		Y = dat$Y,
		Xdyad = dat$Xdyad,
		family = "ordinal",
		R = 0,
		burn = 100,
		nscan = 300,
		odens = 3,
		verbose = FALSE,
		gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with all families with R > 0", {
	skip_on_cran()

	families_to_test = c("normal", "binary")

	for (fam in families_to_test) {
		set.seed(510)
		n = 12
		n_time = 3
		actors = paste0("a", 1:n)

		Y_list = list()
		for (t in 1:n_time) {
			if (fam == "normal") {
				Y = matrix(rnorm(n * n), n, n)
			} else {
				Y = matrix(rbinom(n * n, 1, 0.3), n, n)
			}
			diag(Y) = NA
			rownames(Y) = colnames(Y) = actors
			Y_list[[t]] = Y
		}
		names(Y_list) = paste0("t", 1:n_time)

		fit = lame(
			Y = Y_list,
			family = fam,
			R = 2,
			burn = 100,
			nscan = 300,
			odens = 3,
			verbose = FALSE,
			gof = FALSE
		)

		expect_s3_class(fit, "lame")
		expect_false(
			is.null(fit$U),
			info = paste("U should exist for", fam, "with R=2")
		)
		expect_false(
			is.null(fit$V),
			info = paste("V should exist for", fam, "with R=2")
		)
	}
})

test_that("lame() unipartite binary recovers beta and additive effects", {
	skip_on_cran()

	n = 20
	n_time = 3

	dat = simulate_test_network(
		n = n,
		n_time = n_time,
		family = "binary",
		R = 0,
		beta_intercept = 0,
		beta_dyad = 0.8,
		sd_a = 0.4,
		sd_b = 0.4,
		s2 = 1,
		mode = "unipartite",
		seed = 101
	)

	Xdyad_list = lapply(dat$Xdyad, function(x) x[, , 1])

	fit = lame(
		dat$Y,
		Xdyad = Xdyad_list,
		R = 0,
		family = "binary",
		burn = 300,
		nscan = 1200,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 42
	)

	# the 95% interval of beta_dyad covers the generating 0.8 and excludes zero
	ci = quantile(fit$BETA[, 2], c(0.025, 0.975))
	expect_true(ci[1] > 0 && ci[1] < 0.8 && ci[2] > 0.8, info = paste("Beta dyad interval:", paste(round(ci, 3), collapse = " to ")))

	# additive effects correlation
	sorted_names = names(fit$APM)
	orig_names = rownames(dat$Y[[1]])
	true_a_avg = rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg = rowMeans(do.call(cbind, dat$true_params$b))
	names(true_a_avg) = orig_names
	names(true_b_avg) = orig_names
	cor_a = cor(fit$APM, true_a_avg[sorted_names])
	cor_b = cor(fit$BPM, true_b_avg[sorted_names])
	expect_true(
		cor_a > 0.2,
		info = paste("Sender effects correlation:", round(cor_a, 3))
	)
	expect_true(
		cor_b > 0.2,
		info = paste("Receiver effects correlation:", round(cor_b, 3))
	)
})
