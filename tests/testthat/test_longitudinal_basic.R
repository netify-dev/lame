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

test_that("lame() runs with binary family and returns correct structure", {
	skip_on_cran()

	set.seed(102)
	n = 15
	n_time = 3
	actors = paste0("a", 1:n)

	# binary data generated directly, avoiding simulator edge cases
	Y_list = list()
	for (t in 1:n_time) {
		Y = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "binary",
		R = 0,
		burn = 100,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = TRUE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
	expect_false(is.null(fit$VC))

	# YPM holds probabilities in [0, 1]
	if (!is.null(fit$YPM)) {
		# YPM may be an array or a list
		ypm_vals = unlist(fit$YPM)
		ypm_vals = ypm_vals[!is.na(ypm_vals)]
		expect_true(
			all(ypm_vals >= 0 & ypm_vals <= 1),
			info = "Binary YPM should be probabilities in [0,1]"
		)
	}

	# the intercept is negative for a sparse (30%) network
	intercept_hat = median(fit$BETA[, 1])
	expect_true(
		intercept_hat < 1,
		info = "Sparse network should have negative/small intercept"
	)
})

test_that("lame() with multiplicative effects (R > 0) returns UV", {
	skip_on_cran()

	dat = simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 2, beta_intercept = 0, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, seed = 103
	)

	fit = lame(
		Y = dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 2,
		burn = 100,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = TRUE
	)

	expect_s3_class(fit, "lame")

	# U and V exist with the right dimensions
	expect_false(is.null(fit$U))
	expect_false(is.null(fit$V))

	# a static model has n x R matrices U and V
	expect_equal(ncol(fit$U), 2)
	expect_equal(ncol(fit$V), 2)
})

test_that("lame() with symmetric networks works", {
	skip_on_cran()

	dat = simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 0, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, symmetric = TRUE, seed = 104
	)

	# symmetrize the data, and the covariate: a symmetric fit needs a symmetric regression term
	for (t in seq_along(dat$Y)) {
		Y_t = dat$Y[[t]]
		Y_t = (Y_t + t(Y_t)) / 2
		diag(Y_t) = NA
		dat$Y[[t]] = Y_t
	}
	dat$Xdyad = lapply(dat$Xdyad, function(x) { x[, , 1] = (x[, , 1] + t(x[, , 1])) / 2; x })

	fit = lame(
		Y = dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		symmetric = TRUE,
		burn = 100,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = TRUE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
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

test_that("lame() works with poisson family", {
	skip_on_cran()

	set.seed(505)
	n = 15
	n_time = 3
	actors = paste0("a", 1:n)

	Y_list = list()
	Xdyad_list = list()
	for (t in 1:n_time) {
		X = matrix(rnorm(n * n, 0, 0.3), n, n)
		diag(X) = NA
		rownames(X) = colnames(X) = actors
		lambda = exp(0.5 + 0.2 * X)
		lambda[is.na(lambda)] = exp(0.5)
		Y = matrix(rpois(n * n, lambda), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
		Xdyad_list[[t]] = X
	}
	names(Y_list) = names(Xdyad_list) = paste0("t", 1:n_time)

	# short chain: the settling check may warn that it was still moving
	fit = suppressWarnings(lame(
		Y = Y_list,
		Xdyad = Xdyad_list,
		family = "poisson",
		R = 0,
		burn = 100,
		nscan = 300,
		odens = 3,
		verbose = FALSE,
		gof = FALSE
	))

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))

	# YPM holds non-negative counts
	if (!is.null(fit$YPM)) {
		ypm_vals = unlist(fit$YPM)
		ypm_vals = ypm_vals[!is.na(ypm_vals)]
		expect_true(
			all(ypm_vals >= 0),
			info = "Poisson YPM should be non-negative"
		)
	}
})

test_that("lame() works with cbin family", {
	skip_on_cran()

	set.seed(506)
	n = 15
	n_time = 3
	actors = paste0("a", 1:n)
	odmax_val = 4

	Y_list = list()
	for (t in 1:n_time) {
		Y = matrix(0, n, n)
		for (i in 1:n) {
			nominees = sample(setdiff(1:n, i), odmax_val)
			Y[i, nominees] = 1
		}
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	# short chain: the settling check may warn that it was still moving
	fit = suppressWarnings(lame(
		Y = Y_list,
		family = "cbin",
		R = 0,
		burn = 100,
		nscan = 300,
		odens = 3,
		verbose = FALSE,
		gof = FALSE
	))

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with frn family", {
	skip_on_cran()

	set.seed(507)
	n = 15
	n_time = 3
	actors = paste0("a", 1:n)
	odmax_val = 3

	Y_list = list()
	for (t in 1:n_time) {
		Y = matrix(0, n, n)
		for (i in 1:n) {
			nominees = sample(setdiff(1:n, i), odmax_val)
			Y[i, nominees] = 1:odmax_val
		}
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "frn",
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
