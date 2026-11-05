test_that("Gaussian draws preserve small positive and zero variance directions", {
	set.seed(826211)
	x = rmvnorm(40000, c(0, 0, 3), diag(c(1, 1e-18, 0)))
	expect_identical(x[, 3], rep(3, nrow(x)))
	scaled = sweep(x[, 1:2], 2, c(1, 1e-9), "/")
	expect_lt(max(abs(colMeans(scaled))) * sqrt(nrow(x)), 4)
	expect_lt(max(abs(diag(cov(scaled)) - 1)) * sqrt((nrow(x) - 1) / 2), 4)

	direction = c(2, -3, 1)
	set.seed(826212)
	x = rmvnorm(40000, rep(0, 3), tcrossprod(direction))
	expect_lt(max(abs(x - outer(x[, 1] / direction[1], direction))), 1e-6)
	expect_lt(abs(var(x[, 1] / direction[1]) - 1) * sqrt((nrow(x) - 1) / 2), 4)
	expect_identical(
		rmvnorm(2, c(2, 4), matrix(0, 2, 2)),
		matrix(c(2, 2, 4, 4), 2)
	)
})

test_that("ordinary Gaussian draws retain their Cholesky factor and RNG stream", {
	covariance = matrix(c(1, .3, .3, 2), 2)
	factor = chol(covariance)
	set.seed(826213)
	expected = t(t(matrix(rnorm(40), 20, 2) %*% factor) + c(.2, -.4))
	rng = .Random.seed
	set.seed(826213)
	expect_identical(rmvnorm(20, c(.2, -.4), covariance), expected)
	expect_identical(.Random.seed, rng)
	set.seed(826213)
	expect_identical(rmvnorm(20, c(.2, -.4), covariance, factor), expected)
	expect_identical(.Random.seed, rng)
})

test_that("invalid Gaussian covariances fail before consuming randomness", {
	set.seed(826214)
	rng = .Random.seed
	expect_error(rmvnorm(3, c(0, 0), diag(c(1, -1))), "positive semidefinite")
	expect_error(rmvnorm(3, c(0, 0), matrix(c(1, 1, 0, 1), 2)), "symmetric")
	expect_error(rmvnorm(3, c(0, 0), diag(c(1, NA_real_))), "finite covariance")
	expect_error(rmvnorm(3, c(0, 0), diag(3)), "matching mu")
	expect_identical(.Random.seed, rng)
	expect_identical(dim(rmvnorm(0, c(0, 0), diag(2))), c(0L, 2L))
	expect_identical(dim(rmvnorm(4, numeric(), matrix(numeric(), 0, 0))), c(4L, 0L))
})
