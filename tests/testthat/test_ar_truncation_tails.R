.ar_tail_pit = function(x, mean, sd, lower, upper) {
	if (upper <= mean) {
		hi = pnorm(upper, mean, sd, log.p = TRUE)
		ratio = exp(pnorm(lower, mean, sd, log.p = TRUE) - hi)
		(exp(pnorm(x, mean, sd, log.p = TRUE) - hi) - ratio) / (1 - ratio)
	} else {
		lo = pnorm(lower, mean, sd, lower.tail = FALSE, log.p = TRUE)
		-expm1(pnorm(x, mean, sd, lower.tail = FALSE, log.p = TRUE) - lo) /
			-expm1(pnorm(upper, mean, sd, lower.tail = FALSE, log.p = TRUE) - lo)
	}
}

.ar_expect_uniform_tail = function(x, mean, sd, lower, upper) {
	expect_true(all(is.finite(x) & x > lower & x < upper))
	expect_equal(length(unique(x)), length(x))
	u = .ar_tail_pit(x, mean, sd, lower, upper)
	expect_true(all(is.finite(u) & u > 0 & u < 1))
	n = length(u)
	expect_lt(abs(mean(u) - .5) / sqrt(1 / (12 * n)), 6)
	expect_lt(abs(mean(u^2) - 1 / 3) / sqrt(4 / (45 * n)), 6)
	expect_lt(abs(mean(u < .1) - .1) / sqrt(.09 / n), 6)
}

test_that("actor AR proposals keep the continuous truncated law beyond probability underflow", {
	set.seed(729531)
	for (mean in c(-3, 3)) {
		# n1 = S1 = 0 drops the stationary correction, isolating the proposal law
		x = replicate(4000, lame:::.lame_rho_uv_mh_stationary(.5, mean, .05^2, 0, 0, 1))
		.ar_expect_uniform_tail(x, mean, .05, -.99, .99)
	}
})

test_that("coefficient AR tails follow the exact conditional without boundary atoms", {
	set.seed(915267)
	for (mean in c(-3, 3)) {
		path = matrix(c(.001, .001), 2, 1)
		variance = 1 / (.001^2 + 1 / .05^2)
		posterior_mean = variance * (.001^2 + mean / .05^2)
		x = replicate(4000, lame:::sample_rho_beta_cpp(path, 1L, 1L, matrix(1), 1, .5, mean, .05, 0, .999)[1])
		.ar_expect_uniform_tail(x, posterior_mean, sqrt(variance), 0, .999)
	}
})

test_that("coefficient AR conditionals preserve small positive variances", {
	set.seed(612783)
	x = replicate(5000, lame:::sample_rho_beta_cpp(matrix(0, 2, 1), 1L, 1L, matrix(1), 1, .5, .4, 1e-8, 0, .999)[1])
	expect_lt(abs(mean(x) - .4) / (1e-8 / sqrt(length(x))), 6)
	expect_lt(abs(var(x) / 1e-16 - 1) / sqrt(2 / (length(x) - 1)), 6)
})
