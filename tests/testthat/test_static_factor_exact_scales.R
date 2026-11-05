skip_on_cran()

test_that("zero opposing factors retain the exact prior conditional", {
	skip_on_cran()
	set.seed(47315401L)
	n = 4L
	E = array(0, c(n, n, 3L))
	U = V = matrix(0, n, 1L)
	s2 = 1e-8
	# the factor prior of R/uv_prior.R: df = 2R + 2 degrees of freedom and the scale matrix psi I, psi = df
	df = lame:::.lame_uv_prior_df(1L); psi = df
	for (shrink in c(FALSE, TRUE)) for (backend in c("native", "R")) {
		x = replicate(1500L, {
			z = if (backend == "native") {
				lame:::rUV_rep_fc_cpp(
					E,
					U,
					V,
					0,
					s2,
					diag(1 / sqrt(s2), 2),
					1e-6,
					shrink,
					0L,
					df,
					psi
				)
			} else rUV_rep_fc(E, U, V, 0, s2, shrink, df, psi)
			mean(z$U^2)
		})
		# the conditional of a column given a zero partner: psi / (n + df - 2), the mean of 1 / chi-square(n + df)
		target = if (shrink) psi / (n + df - 2) else n
		expect_lt(
			abs(mean(x) - target),
			6 * sd(x) / sqrt(length(x)),
			label = paste(backend, "shrink", shrink)
		)
	}
})

test_that("Wishart draws retain small positive covariance directions", {
	set.seed(47315402L)
	S = diag(c(1, 1e-16))
	x = replicate(2000L, lame:::rwish_cpp(S, 7L)[2, 2] / S[2, 2])
	expect_lt(abs(mean(x) - 7), 6 * sd(x) / sqrt(length(x)))
	expect_error(lame:::rwish_cpp(diag(c(1, -1)), 7L), "positive definite")
	expect_error(lame:::rwish_cpp(matrix(NA_real_, 2, 2), 7L), "finite")
	expect_error(lame:::rwish_cpp(matrix(1, 2, 3), 7L), "square")
	expect_error(lame:::rwish_cpp(diag(2), 0L), "positive")
})

test_that("large factor states retain their inverse-Wishart conditional law", {
	set.seed(47315403L)
	n = 4L
	amplitude = 1e8
	E = array(0, c(n, n, 2L))
	for (backend in c("native", "R")) {
		U = matrix(amplitude, n, 1L)
		V = matrix(0, n, 1L)
		# a zero opposing column leaves a scaled multivariate t prior
		df = lame:::.lame_uv_prior_df(1L); psi = df
		target = (psi + n * amplitude^2) / (n + df - 2)
		x = replicate(1500L, {
			z = if (backend == "native") {
				lame:::rUV_rep_fc_cpp(E, U, V, 0, 1, diag(2), 1e-6, TRUE, 0L, df, psi)
			} else rUV_rep_fc(E, U, V, 0, 1, TRUE, df, psi)
			mean(z$U^2) / target
		})
		expect_lt(abs(mean(x) - 1), 6 * sd(x) / sqrt(length(x)), label = backend)
	}
})

test_that("symmetric roots preserve positive scales and exact zero directions", {
	for (scale in c(1e-24, 1e-16, 1, 1e8)) {
		S = diag(c(scale, 2 * scale))
		H = mhalf(S)
		expect_equal(H %*% H / scale, S / scale, tolerance = 1e-14)
	}
	expect_identical(mhalf(diag(c(1, 1e-16))), diag(c(1, 1e-8)))
	expect_identical(mhalf(matrix(0, 2, 2)), matrix(0, 2, 2))
	H = mhalf(matrix(1, 2, 2))
	expect_equal(H %*% H, matrix(1, 2, 2), tolerance = 1e-14)
	expect_error(mhalf(diag(c(1, -1e-10))), "positive semidefinite")
	expect_error(mhalf(matrix(NA_real_, 2, 2)), "finite")
	expect_error(mhalf(matrix(1, 2, 3)), "square")
})
