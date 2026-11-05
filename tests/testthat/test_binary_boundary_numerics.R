test_that('observed binary dyad log probabilities remain accurate near reciprocity boundaries', {
	# independent 80/110-digit conditional-integral references, rounded to binary64
	cases = matrix(
		c(
			.36476626058770267, -1.4223376436008821, -.999999941528669, -4782101.588212454,
			-3, -2, -.99999999, -625000028.6325787,
			-3, -2, .99999999, -6.607726221510349
		),
		ncol = 4L,
		byrow = TRUE
	)
	for (k in seq_len(nrow(cases))) {
		a = cases[k, 1L]; b = cases[k, 2L]; rho = cases[k, 3L]; reference = cases[k, 4L]
		tolerance = 3e-11 + 8 * .Machine$double.eps * abs(reference)
		value = .lame_binary_dynamic_bvn_log_quad(a, b, rho)
		expect_lte(abs(value - reference), tolerance)
		expect_lte(abs(.lame_binary_dynamic_bvn_log(a, b, rho) - reference), tolerance)
		expect_lte(abs(.lame_binary_dynamic_bvn_log_quad(b, a, rho) - value), 2 * tolerance)
	}
})

test_that('zero-threshold orthant probabilities use the correct full-support law', {
	for (rho in c(-.99999999999999, -.99999999, -.6, 0, .6, .99999999, .99999999999999)) {
		# atan2 avoids cancellation in pi/2 + asin(rho) at the negative endpoint
		reference = log(atan2(sqrt((1 - rho) * (1 + rho)), -rho) / (2 * pi))
		expect_lte(abs(.lame_binary_dynamic_bvn_log_quad(0, 0, rho) - reference), 3e-11)
	}
})

test_that('normal log-tail displacements preserve sub-ulp physical shifts', {
	# independent 100/140-digit erfc references, reading z+h as a real sum
	cases = matrix(
		c(
			-1e8, 1e-12, .0001,
			-1e8, .1, 9999999.995000001,
			-20, .1, 1.999987663764638,
			-10, .1, 1.0048568501084292
		),
		ncol = 3L,
		byrow = TRUE
	)
	for (k in seq_len(nrow(cases))) {
		value = .lame_normal_log_cdf_shift(cases[k, 1L], cases[k, 2L])
		expect_lte(abs(value - cases[k, 3L]), 3e-13 + 1e-12 * abs(cases[k, 3L]))
	}
})

test_that('normal derivative remains accurate on both Mills-ratio branch boundaries', {
	# independent erfc values at exact binary64 inputs, including nextafter neighbours
	z = c(
		-10.1, -10.000000000000002, -10.0, -9.999999999999998, -9.9, -8.0, -6.0, -5.000000000000001, -5.0, -4.999999999999999, -4.99
	)
	reference = c(
		10.197157545956836, 10.098093233962514, 10.098093233962512, 10.09809323396251, 9.999046787486824, 8.121368112236112, 6.158482604544599, 5.186503967125843, 5.186503967125842, 5.186503967125841, 5.17683147360947
	)
	# libm erfc differs across platforms by a few ulps, so allow 1e-12 relative error
	expect_true(all(abs(.lame_normal_mills_ratio(z) - reference) <= 1e-12 * pmax(1, abs(reference))))
})
