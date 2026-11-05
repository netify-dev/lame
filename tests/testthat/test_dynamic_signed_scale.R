skip_on_cran()

test_that("tiny symmetric coordinates retain their signed factor weights", {
	U = array(c(1, -2, 3) * 1e-160, c(3, 1, 1))
	E = array(c(0, 1.3, -.6, 1.3, 0, .8, -.6, .8, 0), c(3, 3, 1))
	original = serialize(list(U, E), NULL)
	N = 10000
	for (lambda in c(-2, 0, 3)) {
		set.seed(826302)
		draws = replicate(
			N,
			rUV_dynamic_fc(
				U,
				U * lambda,
				E,
				0,
				1,
				1,
				shrink = FALSE,
				symmetric = TRUE
			)$U[, 1, 1]
		)
		# actor 2's conditional after actor 1 is drawn
		precision = 1 + (lambda * draws[1, ])^2 + (lambda * U[3, 1, 1])^2
		linear = E[2, 1, 1] * lambda * draws[1, ] + E[2, 3, 1] * lambda * U[3, 1, 1]
		z = (draws[2, ] - linear / precision) * sqrt(precision)
		expect_lt(abs(mean(z)) * sqrt(N), 4)
		expect_lt(abs(var(z) - 1) * sqrt((N - 1) / 2), 4)
		expect_lt(abs(mean(z * draws[1, ])) * sqrt(N), 4)
	}
	expect_identical(serialize(list(U, E), NULL), original)
})
