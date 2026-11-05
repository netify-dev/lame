# reference by numerical integration of the defining circular integral after rescaling its peak (no Bessel functions)
.orientation_log_i0_reference = function(k) {
	v = integrate(
		function(u) exp(-2 * k * sin(u / (2 * sqrt(k)))^2),
		0,
		min(pi * sqrt(k), 12),
		rel.tol = 1e-12
	)$value
	log(v / (pi * sqrt(k)))
}

test_that("scaled circular normalization remains accurate in the Bessel tail", {
	x = c(1e4, 1e4 + 1, 1e5, 2e5, 5e5, 1.4e6)
	truth = vapply(x, .orientation_log_i0_reference, numeric(1))
	expect_equal(lame:::.lame_uv_log_i0_scaled(x), truth, tolerance = 2e-14)
	ordinary = c(0, .1, 1, 10, 100, 1e4)
	expect_identical(
		lame:::.lame_uv_log_i0_scaled(ordinary),
		log(besselI(ordinary, 0, expon.scaled = TRUE))
	)
})

test_that("large-concentration orthogonal draws retain both determinant masses", {
	M = diag(c(300000, .35))
	k = c(300000 + .35, 300000 - .35)
	mass = k + vapply(k, .orientation_log_i0_reference, numeric(1))
	prob = plogis(mass[1] - mass[2])
	set.seed(716345)
	q = replicate(6000, as.numeric(lame:::.lame_uv_o2(M)))
	d = q[1, ] * q[4, ] - q[2, ] * q[3, ]
	expect_true(all(is.finite(q)))
	expect_equal(abs(d), rep(1, ncol(q)), tolerance = 1e-12)
	expect_lt(abs(mean(d > 0) - prob), 6 * sqrt(prob * (1 - prob) / ncol(q)))
	# each component has its own sharply concentrated von Mises angle
	scaled_sine = q[2, ] * sqrt(ifelse(d > 0, k[1], k[2]))
	expect_lt(abs(mean(scaled_sine)), 6 * sd(scaled_sine) / sqrt(length(d)))
	reference_second = function(k) {
		weight = function(u) exp(-2 * k * sin(u / (2 * sqrt(k)))^2)
		integrate(
			function(u) k * sin(u / sqrt(k))^2 * weight(u),
			0,
			12,
			rel.tol = 1e-12
		)$value / integrate(weight, 0, 12, rel.tol = 1e-12)$value
	}
	second = prob * reference_second(k[1]) + (1 - prob) * reference_second(k[2])
	expect_lt(
		abs(mean(scaled_sine^2) - second),
		6 * sd(scaled_sine^2) / sqrt(length(d))
	)
	expect_silent(lame:::.lame_uv_o2(diag(c(400000, 100000))))
})
