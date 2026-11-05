skip_on_cran()

test_that("binary row conditional samples proper far tails and infinite intervals", {
	row_draw = getFromNamespace(".lame_binary_row_truncnorm", "lame")
	set.seed(47334151)
	for (lim in list(c(8, 8.001), c(10, 11), c(-11, -10), c(10, Inf), c(-Inf, -10))) {
		x = row_draw(rep(0, 3000), 1, rep(lim[1], 3000), rep(lim[2], 3000))
		expect_true(all(is.finite(x) & x > lim[1] & x < lim[2]))
		mode = min(max(0, lim[1]), lim[2])
		density = function(z) exp(-mode * z - z^2 / 2)
		lo = lim[1] - mode
		hi = lim[2] - mode
		mass = integrate(density, lo, hi, rel.tol = 1e-11)$value
		cdf = function(z) integrate(density, lo, z, rel.tol = 1e-11)$value / mass
		q = vapply(c(.1, .5, .9), function(p) {
			uniroot(function(z) cdf(z) - p, c(max(lo, -12), min(hi, 12)), tol = 1e-11)$root + mode
		}, 0.0)
		proportion = vapply(q, function(cut) mean(x <= cut), 0.0)
		z = (proportion - c(.1, .5, .9)) / sqrt(c(.1, .5, .9) * c(.9, .5, .1) / length(x))
		expect_lt(max(abs(z)), 6)
	}
	input = list(Z = matrix(.5, 3, 3), Y = matrix(1, 3, 3), a = rep(10.5, 3), b = rep(0, 3), Sab = diag(2), SS = 1)
	diag(input$Y) = NA
	got = do.call(raSab_bin_fc, input)
	expect_true(all(is.finite(got$a) & got$a > 10))
	expect_true(all(is.finite(got$Z)))
	expect_equal(got$Z - matrix(got$a, 3, 3), input$Z - matrix(input$a, 3, 3), tolerance = 2e-14)
	expect_true(all(got$Z[!is.na(input$Y)] > 0))
	expect_true(all(eigen(got$Sab, symmetric = TRUE)$values > 0))
})

test_that("central binary row conditionals preserve complete draws and RNG", {
	row_draw = getFromNamespace(".lame_binary_row_truncnorm", "lame")
	mean = c(-.2, 0, .4)
	lower = c(-2, -Inf, 0)
	upper = c(2, Inf, 3)
	set.seed(47334152)
	before = .Random.seed
	expected = mean + .8 * qnorm(runif(3, pnorm((lower - mean) / .8), pnorm((upper - mean) / .8)))
	RNG = .Random.seed
	assign(".Random.seed", before, envir = .GlobalEnv)
	expect_identical(row_draw(mean, .8, lower, upper), expected)
	expect_identical(.Random.seed, RNG)
	e = new.env(parent = asNamespace("lame"))
	e$.lame_binary_row_truncnorm = function(mean, sd, lower, upper) {
		mean + sd * qnorm(runif(length(mean), pnorm((lower - mean) / sd), pnorm((upper - mean) / sd)))
	}
	old = raSab_bin_fc
	environment(old) = e
	input = list(
		Z = matrix(c(.5, -.5, .7, -.3, .4, -.6, .2, -.8, .9), 3),
		a = c(.1, -.2, .3),
		b = c(-.1, .1, .2),
		Sab = matrix(c(.8, .1, .1, .7), 2),
		SS = 2
	)
	input$Y = (input$Z > 0) * 1
	diag(input$Y) = NA
	before = .Random.seed
	expected = do.call(old, input)
	RNG = .Random.seed
	assign(".Random.seed", before, envir = .GlobalEnv)
	expect_identical(do.call(raSab_bin_fc, input), expected)
	expect_identical(.Random.seed, RNG)
})

test_that("invalid binary row conditionals refuse before consuming RNG", {
	row_draw = getFromNamespace(".lame_binary_row_truncnorm", "lame")
	set.seed(47334153)
	before = .Random.seed
	expect_error(row_draw(0, 0, -Inf, Inf), class = "lame_binary_row_conditional_error")
	expect_identical(.Random.seed, before)
	expect_error(row_draw(0, 1, 1, 1), class = "lame_binary_row_conditional_error")
	expect_identical(.Random.seed, before)
	expect_error(row_draw(1e308, 1, -1e308, 1e308), class = "lame_binary_row_numerical_error")
	expect_identical(.Random.seed, before)
})

test_that("joint binary row and covariance update preserves the original proper prior", {
	set.seed(47334154)
	n = 4
	N = 1800
	eta = 12
	S0 = matrix(c(.8, .13, .13, .6), 2)
	changes = matrix(NA_real_, N, 6)
	carry_error = 0
	all_signs = TRUE
	for (k in seq_len(N)) {
		S = solve(rWishart(1, eta, solve(eta * S0))[, , 1])
		ab = matrix(rnorm(n * 2), n, 2) %*% chol(S)
		a = ab[, 1]
		b = ab[, 2]
		E = matrix(rnorm(n * n), n)
		for (i in seq_len(n - 1)) for (j in (i + 1):n) E[j, i] = .6 * E[i, j] + .8 * E[j, i]
		Z = E + outer(a, b, "+")
		Y = (Z > 0) * 1
		diag(Y) = NA
		got = raSab_bin_fc(Z, Y, a, b, S, S0, eta, SS = 2)
		stat = function(a, S) c(a[1], a[1]^2, a[1] * b[1], S[1, 1], S[1, 2], S[2, 2])
		changes[k, ] = stat(got$a, got$Sab) - stat(a, S)
		carry_error = max(carry_error, abs(got$Z - matrix(got$a, n, n) - (Z - matrix(a, n, n))))
		all_signs = all_signs && all((got$Z[!is.na(Y)] > 0) == (Y[!is.na(Y)] == 1))
	}
	z = colMeans(changes) / (apply(changes, 2, sd) / sqrt(N))
	expect_true(all(is.finite(z)))
	expect_lt(max(abs(z)), 6)
	expect_lt(carry_error, 2e-14)
	expect_true(all_signs)
})
