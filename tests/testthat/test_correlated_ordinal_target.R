test_that("correlated explicit ordinal uses each partner's updated conditional value", {
	z = matrix(c(0, .2, -.8, 0), 2)
	e = matrix(c(0, .3, -.2, 0), 2)
	y = matrix(c(NA, 2, 3, NA), 2)
	for (rho in c(-.7, .6)) {
		set.seed(47328211)
		expected = z
		cuts = c(-Inf, 0, 1, Inf)
		sd = sqrt(1 - rho^2)
		for (j in 1:2) for (i in 1:2) {
			if (i == j) next
			mu = e[i, j] + rho * (expected[j, i] - e[j, i])
			lo = cuts[y[i, j]]; hi = cuts[y[i, j] + 1]
			expected[i, j] = mu + sd * qnorm(pnorm((lo - mu) / sd) + runif(1) *
					(pnorm((hi - mu) / sd) - pnorm((lo - mu) / sd)))
		}
		diag(expected) = rnorm(2, diag(e), sqrt(1 + rho))
		rng = .Random.seed
		set.seed(47328211)
		actual = lame:::rZ_ord_explicit_fc(z, e, y, c(0, 1), rho)
		expect_equal(actual, expected, tolerance = 1e-13)
		expect_identical(.Random.seed, rng)
	}
})

test_that("cutpoints given Z follow the unchanged flat-prior uniform conditional", {
	z = matrix(c(0, -1, .4, -2, 0, 2, .1, 3, 0), 3)
	y = matrix(c(NA, 1, 2, 1, NA, 3, 2, 3, NA), 3)
	set.seed(47328212)
	u = runif(1)
	rng = .Random.seed
	set.seed(47328212)
	a = lame:::.sample_alpha_given_z(c(0, 1), y, z)
	expect_identical(a$alpha, c(0, .4 + (2 - .4) * u))
	expect_identical(.Random.seed, rng)
	set.seed(47328213)
	before = .Random.seed
	z[2, 1] = .2
	expect_error(lame:::.sample_alpha_given_z(c(0, 1), y, z), "category support")
	expect_identical(.Random.seed, before)
})

test_that("correlated explicit ordinal handles partial observations and rejects invalid rho", {
	z = matrix(0, 3, 3)
	e = matrix(0, 3, 3)
	y = matrix(NA_real_, 3, 3); y[1, 2] = 2; y[2, 1] = 3
	for (rho in c(-.99, .99)) {
		set.seed(47328214)
		x = lame:::rZ_ord_explicit_fc(z, e, y, c(0, 1), rho)
		expect_true(all(is.finite(x)))
		expect_true(x[1, 2] > 0 && x[1, 2] < 1 && x[2, 1] > 1)
	}
	before = .Random.seed
	expect_error(lame:::rZ_ord_explicit_fc(z, e, y, c(0, 1), 1), "correlation")
	expect_identical(.Random.seed, before)
})
