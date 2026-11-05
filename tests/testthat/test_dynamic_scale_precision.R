test_that("diagonal innovation precision retains coefficient units", {
	L = diag(c(1, 1e-20))
	expect_identical(lame:::.lame_beta_scale_precision(L), diag(c(1, 1e20)))
	for (L in list(diag(2), diag(c(.001, 10, 500)), matrix(7, 1, 1))) {
		expect_identical(lame:::.lame_beta_scale_precision(L), solve(L))
	}
	set.seed(47327012)
	rng = .Random.seed
	for (bad in list(
		diag(c(0, 1)),
		diag(c(-1, 1)),
		diag(c(Inf, 1)),
		matrix(c(1, .1, .1, 1), 2),
		diag(c(1e-320, 1))
	)) {
		expect_error(lame:::.lame_beta_scale_precision(bad), "innovation")
	}
	expect_identical(.Random.seed, rng)
})

test_that("unequally scaled coefficient paths retain their variance conditional", {
	Lambda = diag(c(1, 1e-20))
	path = cbind(c(0, .5, -.5, 1), 1e-10 * c(0, 1, 2, -1))
	rho = .3
	innov = path[-1, , drop = FALSE] - rho * path[-nrow(path), , drop = FALSE]
	scatter = sum(sweep(innov, 2, sqrt(diag(Lambda)), "/")^2)
	set.seed(47327011)
	expected = sqrt((2 + scatter / 2) / rgamma(4000, 3 + length(innov) / 2, 1))
	rng = .Random.seed
	set.seed(47327011)
	observed = replicate(
		4000,
		lame:::sample_sigma_beta_cpp(
			path,
			c(1L, 1L),
			1L,
			lame:::.lame_beta_scale_precision(Lambda),
			rho,
			3,
			2
		)
	)
	expect_equal(as.numeric(observed), expected, tolerance = 1e-13)
	expect_identical(.Random.seed, rng)
})
