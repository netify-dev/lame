test_that("nomination row conditionals retain valid tiny Gaussian variance", {
	skip_on_cran()
	set.seed(47334361)
	n = 3
	Y = matrix(0, n, n)
	diag(Y) = NA
	S = matrix(c(1e-20, 2e-11, 2e-11, 1), 2)
	b = rep(.3, n)
	mu = S[1, 2] / S[2, 2] * b[1]
	v = S[1, 1] - S[1, 2]^2 / S[2, 2]
	args = list(Z = matrix(-100, n, n), Y = Y, a = rep(0, n), b = b, Sab = S, odmax = 3, odobs = rep(0, n), SS = 1)
	for (family in c("cbin", "frn")) {
		input = args
		if (family == "frn") input$YL = matrix(NA_integer_, n, n - 1)
		f = if (family == "cbin") raSab_cbin_fc else raSab_frn_fc
		x = replicate(3000, do.call(f, input)$a[1])
		expect_true(all(is.finite(x)))
		expect_lt(abs((mean(x) - mu) / sqrt(v / length(x))), 6)
		expect_lt(abs((var(x) / v - 1) / sqrt(2 / (length(x) - 1))), 6)
		got = do.call(f, input)
		expect_equal(got$Z - matrix(got$a, n, n), input$Z - matrix(input$a, n, n), tolerance = 2e-14)
	}
})

test_that("checked Gaussian truncation preserves ordinary draws and RNG", {
	checked = getFromNamespace(".lame_checked_truncated_gaussian", "lame")
	old = getFromNamespace("rtnorm_interval_logp", "lame")
	mean = c(-.2, 0, .4)
	lower = c(-2, -Inf, 0)
	upper = c(2, Inf, 3)
	set.seed(47334362)
	before = .Random.seed
	expected = .8 * old(mean / .8, lower / .8, upper / .8)
	RNG = .Random.seed
	assign(".Random.seed", before, envir = .GlobalEnv)
	expect_identical(checked(mean, .8, lower, upper), expected)
	expect_identical(.Random.seed, RNG)
	expect_error(checked(0, 0, -Inf, Inf), class = "lame_truncation_parameter_error")
	expect_identical(.Random.seed, RNG)
	expect_error(checked(0, 1, 2, 1), class = "lame_truncation_support_error")
	expect_identical(.Random.seed, RNG)
})

test_that("latent truncation failures never substitute the previous state", {
	n = 3
	Z = matrix(1, n, n)
	Y = Z
	diag(Y) = NA
	EZ = matrix(-4e307, n, n)
	expect_error(lame:::rZ_bin_sym_fc(Z, EZ, Y), class = "lame_truncation_numerical_error")
	Z = matrix(c(0, 1, 2, 1, 0, 2, 2, 2, 0), n)
	Y = Z
	diag(Y) = NA
	EZ = matrix(4e307, n, n)
	EZ[which(Y == 2)] = -4e307
	diag(EZ) = 0
	expect_error(lame:::rZ_ord_sym_fc(Z, EZ, Y), class = "lame_truncation_numerical_error")
	expect_error(lame:::rZ_ord_fc(Z, EZ, 0, Y), class = "lame_truncation_numerical_error")
})

test_that("rZ_ord_sym_fc reports infeasible windows and samples valid extreme tails", {
	set.seed(6886)
	n = 10
	Y = matrix(sample(1:3, n * n, replace = TRUE), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	# a block of actors far below every cutpoint saturates the normal cdf
	EZ = matrix(0, n, n)
	EZ[1:3, ] = -30
	EZ[, 1:3] = -30
	# a middle-category interval with lower bound 3 and upper bound 1 has no support
	y_bad = matrix(c(NA, 1, 2, 1, NA, 3, 2, 3, NA), 3)
	z_bad = matrix(c(0, 3, 2, 3, 0, 1, 2, 1, 0), 3)
	set.seed(1)
	expect_error(rZ_ord_sym_fc(z_bad, matrix(0, 3, 3), y_bad), class = "lame_truncation_support_error")
	# category-ordered values give valid starting intervals for the tail check
	Z = Y
	diag(Z) = 0
	for (k in 1:30) Z = rZ_ord_sym_fc(Z, EZ, Y)
	expect_true(all(is.finite(Z)))
	expect_equal(Z, t(Z))
})
