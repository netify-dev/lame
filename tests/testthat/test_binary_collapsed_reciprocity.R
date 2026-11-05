test_that("binary-pair log probabilities agree with analytic and high-precision values", {
	rho = seq(-.999, .999, length.out = 51)
	expected = log(.25 + asin(rho) / (2 * pi))
	expect_equal(.lame_binary_dynamic_bvn_log(rep(0, 51), rep(0, 51), rho), expected, tolerance = 1e-10)
	expect_equal(
		.lame_binary_dynamic_bvn_log(-40, -30, 0),
		pnorm(-40, log.p = TRUE) + pnorm(-30, log.p = TRUE),
		tolerance = 1e-10
	)
	# independent 256-bit Romberg integration, normalized omitted tail < exp(-120)
	expect_equal(
		.lame_binary_dynamic_bvn_log(-4, -4, -.99),
		-1611.86322643062957071246706,
		tolerance = 1e-10
	)
})

test_that("exact truncated pairs reproduce analytic zero-threshold moments", {
	set.seed(825317)
	for (rho in c(-.95, .95, .99)) {
		P = .25 + asin(rho) / (2 * pi)
		mean = (1 + rho) * dnorm(0) * .5 / P
		second = 1 + rho * sqrt(1 - rho * rho) / (2 * pi * P)
		cross = rho + sqrt(1 - rho * rho) / (2 * pi * P)
		z = .lame_binary_dynamic_biv_positive(rep(0, 12000), rep(0, 12000), rho)
		values = cbind(z, z^2, z[, 1] * z[, 2])
		error = (colMeans(values) - c(mean, mean, second, second, cross)) /
			(apply(values, 2, sd) / sqrt(nrow(values)))
		expect_true(all(z > 0)); expect_lt(max(abs(error)), 6)
		# the marginal-mode envelope keeps this edge case to a practical number of proposals
		expect_lt(attr(z, 'iterations'), 100L)
	}
})

test_that("single observed dyads refresh their partner under the conditional Gaussian law", {
	set.seed(321681)
	eta = array(0, c(2, 2, 1)); eta[1, 2, 1] = .8; eta[2, 1, 1] = -.3
	Y = array(NA_real_, c(2, 2, 1)); Y[1, 2, 1] = 1
	cx = .lame_binary_dynamic_collapsed_rho_context(Y); rho = -.7
	expect_identical(.lame_binary_dynamic_collapsed_rho_loglik(rho, eta, cx), 0)
	z = replicate(8000, {
		zz = .lame_binary_dynamic_collapsed_rho_z(eta, Y, rho, cx)
		c(zz[1, 2, 1] - .8, zz[2, 1, 1] + .3)
	})
	ratio = dnorm(.8) / pnorm(.8); moment1 = 1 - .8 * ratio
	target = c(ratio, rho * ratio, moment1, 1 - rho^2 + rho^2 * moment1, rho * moment1)
	values = rbind(z, z^2, z[1, ] * z[2, ]);mcse = (rowMeans(values) - target) / (apply(values, 1, sd) / sqrt(ncol(values)))
	expect_true(all(z[1, ] > -.8));expect_lt(max(abs(mcse)), 6)
})
