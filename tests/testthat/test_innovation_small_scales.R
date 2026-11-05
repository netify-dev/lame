test_that("small coefficient innovation scales retain the inverse-Gamma conditional", {
	set.seed(863152);N = 5000L
	sigma = replicate(N, lame:::sample_sigma_beta_cpp(matrix(0, 3, 1), 1L, 1L, matrix(1), .8, 2, 1e-16)[1])
	expect_true(all(is.finite(sigma) & sigma > 0))
	expect_equal(length(unique(sigma)), N)
	# shape 2 + two scalar transitions / 2 = 3; scale stays 1e-16
	u = pgamma(1e-16 / sigma^2, shape = 3)
	expect_lt(abs(mean(u) - .5) / sqrt(1 / (12 * N)), 6)
	expect_lt(abs(mean(u * u) - 1 / 3) / sqrt(4 / (45 * N)), 6)
})

test_that("pooled scale updates retain precision from small innovation variances", {
	set.seed(471382);N = 5000L
	x = replicate(N, lame:::.pool_update_sigma_scale(c(1e-8, 2e-8), alpha = 2, hyper_a = 3, hyper_b = .5))
	u = pgamma(x, shape = 7, rate = .5 + 1e16 + 2.5e15)
	expect_true(all(is.finite(x) & x > 0))
	expect_lt(abs(mean(u) - .5) / sqrt(1 / (12 * N)), 6)
	expect_lt(abs(mean(u * u) - 1 / 3) / sqrt(4 / (45 * N)), 6)
})
