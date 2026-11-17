# the running sum of the posterior predictive takes the replicate's mean given each draw

predictive_mean_term = lame:::.predictive_mean_term

test_that("the predictive mean term is the closed-form mean of each family", {
	skip_on_cran()
	set.seed(11)
	n = 5
	EZ = matrix(rnorm(n * n), n, n)
	Ys = 1 * matrix(rbinom(n * n, 1, 0.5), n, n)
	diag(Ys) = NA
	alpha = c(0, 0.7, 1.6)
	levels = c(1, 2, 4, 7)

	off = row(EZ) != col(EZ)
	expect_equal(predictive_mean_term(Ys, EZ, "binary")[off], pnorm(EZ)[off])
	expect_equal(
		predictive_mean_term(Ys, EZ, "binary", s2 = 2.5)[off],
		pnorm(EZ / sqrt(2.5))[off]
	)
	expect_equal(predictive_mean_term(Ys, EZ, "normal")[off], EZ[off])

	cuts = c(-Inf, alpha, Inf)
	ord_mean = EZ * 0
	for (k in seq_along(levels)) {
		ord_mean = ord_mean + levels[k] * (pnorm(cuts[k + 1] - EZ) - pnorm(cuts[k] - EZ))
	}
	expect_equal(
		predictive_mean_term(Ys, EZ, "ordinal", alpha = alpha, levels = levels)[off],
		ord_mean[off]
	)

	# families without a closed form, and rank ordinal, keep the replicate
	for (fam in c("cbin", "frn", "poisson", "ordinal")) {
		expect_identical(predictive_mean_term(Ys, EZ, fam), Ys)
	}

	# structural cells: missing where the replicate is, the replicate where the predictor is not finite
	EZ[1, 2] = NA
	EZ[2, 3] = Inf
	term = predictive_mean_term(Ys, EZ, "binary")
	expect_identical(is.na(term), is.na(Ys))
	expect_identical(term[1, 2], Ys[1, 2])
	expect_identical(term[2, 3], Ys[2, 3])

	# symmetric fits use the symmetrised predictor, slice by slice for a panel
	EZ = array(rnorm(n * n * 2), c(n, n, 2))
	Ya = array(0, c(n, n, 2))
	term = predictive_mean_term(Ya, EZ, "normal", symmetric = TRUE)
	for (k in 1:2) {
		expect_equal(term[, , k], (EZ[, , k] + t(EZ[, , k])) / 2)
	}
})

test_that("the predictive mean term matches the average of the site simulators", {
	skip_on_cran()
	set.seed(12)
	n = 5
	M = 4000
	EZ = matrix(rnorm(n * n, 0, 0.8), n, n)
	alpha = c(0, 0.8)
	levels = c(0, 1, 2)
	mirror = function(Ys) {
		Ys[lower.tri(Ys)] = 0
		Ys + t(Ys)
	}
	sim_mean = function(draw) Reduce(`+`, lapply(seq_len(M), function(m) draw())) / M
	off = row(EZ) != col(EZ)
	EZs = (EZ + t(EZ)) / 2

	# binary directed, binary symmetric, bipartite binary drawn with sd sqrt(s2)
	avg = sim_mean(function() lame:::simY_bin(EZ, 0.4))
	expect_lt(max(abs(avg - predictive_mean_term(avg, EZ, "binary"))[off]), 4 * sqrt(0.25 / M))
	avg = sim_mean(function() mirror(lame:::simY_bin(EZs, 0)))
	term = predictive_mean_term(avg, EZ, "binary", symmetric = TRUE)
	expect_lt(max(abs(avg - term)[off]), 4 * sqrt(0.25 / M))
	EZb = EZ[, 1:3]
	avg = sim_mean(function() 1 * (EZb + matrix(rnorm(n * 3, 0, sqrt(2)), n, 3) > 0))
	expect_lt(max(abs(avg - predictive_mean_term(avg, EZb, "binary", s2 = 2))), 4 * sqrt(0.25 / M))

	# normal with residual variance 1.5
	avg = sim_mean(function() lame:::simY_nrm(EZ, 0.4, 1.5))
	expect_lt(max(abs(avg - predictive_mean_term(avg, EZ, "normal"))[off]), 4 * sqrt(1.5 / M))

	# explicit-cutpoint ordinal, directed and symmetric
	avg = sim_mean(function() lame:::.simY_ord_explicit(EZ, 0.4, alpha, levels))
	term = predictive_mean_term(avg, EZ, "ordinal", alpha = alpha, levels = levels)
	expect_lt(max(abs(avg - term)[off]), 4 * sqrt(1 / M))
	avg = sim_mean(function() mirror(lame:::.simY_ord_explicit(EZs, 0, alpha, levels, TRUE)))
	term = predictive_mean_term(
		avg,
		EZ,
		"ordinal",
		alpha = alpha,
		levels = levels,
		symmetric = TRUE
	)
	expect_lt(max(abs(avg - term)[off]), 4 * sqrt(1 / M))
})

test_that("a short binary fit has fitted probabilities strictly inside (0, 1)", {
	skip_on_cran()
	set.seed(13)
	n = 12
	nm = sprintf("a%02d", seq_len(n))
	Y = 1 * (matrix(rnorm(n * n), n, n) + outer(rnorm(n), rnorm(n), "+") > 0)
	diag(Y) = NA
	dimnames(Y) = list(nm, nm)
	fit = suppressWarnings(ame(
		Y,
		R = 0,
		family = "binary",
		burn = 10,
		nscan = 40,
		odens = 2,
		chains = 1,
		verbose = FALSE
	))
	p = fitted(fit)[row(Y) != col(Y)]
	expect_true(all(p > 0 & p < 1))
	# 20 stored draws: a share of simulated ties would be a multiple of 1/20
	expect_true(any(abs(p * 20 - round(p * 20)) > 1e-6))
})

test_that("YPM keeps the missing-value pattern of the replicates", {
	skip_on_cran()
	set.seed(14)
	n = 10
	nm = sprintf("a%02d", seq_len(n))
	mk = function() {
		Y = 1 * (matrix(rnorm(n * n), n, n) > 0)
		diag(Y) = NA
		dimnames(Y) = list(nm, nm)
		Y
	}
	Y = mk()
	Y[2, 5] = NA
	fit = suppressWarnings(ame(
		Y,
		R = 0,
		family = "binary",
		burn = 10,
		nscan = 20,
		odens = 2,
		chains = 1,
		verbose = FALSE
	))
	expect_identical(unname(is.na(fit$YPM)), row(Y) == col(Y))

	Y_list = list(mk(), mk())
	Y_list[[2]][3, 7] = NA
	fit = suppressWarnings(lame(
		Y_list,
		R = 0,
		family = "binary",
		burn = 10,
		nscan = 20,
		odens = 2,
		chains = 1,
		verbose = FALSE
	))
	for (k in 1:2) {
		expect_identical(unname(is.na(fit$YPM[[k]])), row(Y) == col(Y))
	}
})

test_that("cbin, frn and rank ordinal average ten replicates drawn on their own stream", {
	skip_on_cran()
	set.seed(21)
	n = 6
	EZ = matrix(rnorm(n * n), n, n)
	draw = function() lame:::simY_frn(EZ, 0, 2)
	Ys = draw()
	before = .Random.seed
	term = predictive_mean_term(Ys, EZ, "frn", draw = draw, seed = 5, index = 3)
	# the caller's stream is untouched, also when a replicate fails
	expect_identical(.Random.seed, before)
	expect_error(predictive_mean_term(Ys, EZ, "cbin", draw = function() stop("no draw"), seed = 5, index = 4))
	expect_identical(.Random.seed, before)

	set.seed(lame:::.predictive_stream_seed(5, 3))
	reps = c(list(Ys), lapply(1:9, function(r) draw()))
	expect_equal(term, Reduce(`+`, reps) / 10)
	expect_identical(is.na(term), is.na(Ys))
	expect_identical(predictive_mean_term(Ys, EZ, "frn", draw = draw, seed = 5, index = 3), term)
	expect_false(identical(predictive_mean_term(Ys, EZ, "frn", draw = draw, seed = 5, index = 4), term))
})

test_that("a short cbin fit has YPM on the grid of ten replicates per draw", {
	skip_on_cran()
	set.seed(22)
	n = 12
	nm = sprintf("a%02d", seq_len(n))
	z = matrix(rnorm(n * n), n, n) + outer(rnorm(n), rnorm(n), "+")
	diag(z) = NA
	Y = t(apply(z, 1, function(r) 1 * (rank(-r, na.last = TRUE) <= 3 & r > 0)))
	diag(Y) = NA
	Y[2, 5] = NA
	dimnames(Y) = list(nm, nm)
	fit = suppressMessages(suppressWarnings(ame(
		Y,
		R = 0,
		family = "cbin",
		odmax = 3,
		burn = 10,
		nscan = 40,
		odens = 2,
		chains = 1,
		verbose = FALSE
	)))
	draws = nrow(fit$VC)
	p = fit$YPM[row(Y) != col(Y)]
	# the diagonal stays missing and the hidden cell is predicted, as with one replicate
	expect_identical(unname(is.na(fit$YPM)), row(Y) == col(Y))
	expect_lt(max(abs(p * 10 * draws - round(p * 10 * draws))), 1e-6)
	expect_true(any(abs(p * draws - round(p * draws)) > 1e-6))
})
