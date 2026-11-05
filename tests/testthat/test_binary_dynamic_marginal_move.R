test_that("noncentered AR variance and coordinate Jacobians agree with the original density", {
	set.seed(717634); rho = .7; variance = .2; proposal = .6
	paths = matrix(rnorm(18), 6, 3)
	paths[, 1] = paths[, 1] * sqrt(variance / (1 - rho * rho))
	for (t in 2:3) paths[, t] = rho * paths[, t - 1] + sqrt(variance) * paths[, t]
	innovations = .lame_binary_dynamic_nc_ar_innovations(paths, rho, sqrt(variance))
	vnew = variance * (1 - proposal^2) / (1 - rho^2)
	pnew = .lame_binary_dynamic_nc_ar_paths(innovations, proposal, vnew)
	density = function(p, r, v) sum(dnorm(p[, 1], sd = sqrt(v / (1 - r * r)), log = TRUE)) +
		sum(dnorm(p[, -1] - r * p[, -3], sd = sqrt(v), log = TRUE)) +
		dgamma(1 / v, 2, rate = .2, log = TRUE) - 2 * log(v) + dnorm(r, .7, .2, log = TRUE)
	actual = density(pnew, proposal, vnew) - density(paths, rho, variance) +
		(6 * (3 - 1) / 2 + 1) * log(vnew / variance)
	expected = dnorm(proposal, .7, .2, log = TRUE) - dnorm(rho, .7, .2, log = TRUE) -
		2 * log(vnew / variance) - .2 / vnew + .2 / variance
	expect_equal(actual, expected, tolerance = 1e-10)
	expect_equal(pnew[, 1], paths[, 1], tolerance = 1e-12)
})

test_that("a valid exact AR support endpoint defers the noncentered proposal", {
	set.seed(697184); paths = matrix(rnorm(12, sd = .1), 4, 3)
	rng = .Random.seed
	for (rho in c(-.99, .99)) {
		z = .lame_binary_dynamic_nc_ar_step(paths, rho, .1, function(p) 0, .7, .2, 2, .2)
		expect_identical(z$paths, paths); expect_identical(z$rho, rho)
		expect_identical(z$attempts, 0L); expect_identical(z$accepts, 0L)
	}
	expect_identical(.Random.seed, rng)
	# past the support edge the move is skipped rather than the fit stopped
	z = .lame_binary_dynamic_nc_ar_step(paths, 1, .1, function(p) 0, .7, .2, 2, .2)
	expect_identical(z$paths, paths); expect_identical(z$attempts, 0L)
})
