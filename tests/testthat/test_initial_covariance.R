test_that('computed covariance starts remain in their active support', {
	bad = matrix(c(0, .006543598, .006543598, 0), 2)
	x = lame:::.lame_start_sab_covariance(bad)
	expect_equal(diag(x), rep(1e-6, 2));expect_equal(x[1, 2], 0)
	expect_true(all(eigen(x, symmetric = TRUE, only.values = TRUE)$values > 0))
	good = matrix(c(.2, .19999, .19999, .2), 2)
	expect_identical(lame:::.lame_start_sab_covariance(good), good)
	expect_equal(lame:::.lame_start_sab_covariance(bad, c(TRUE, FALSE)), diag(c(1e-6, 0)))
})

test_that('ALS warm starts repair incompatible variance-component cross moments', {
	# the public converter on an ordinary ALS fit
	set.seed(51347);Y = matrix(rnorm(36), 6);diag(Y) = NA
	f = ame_als(Y, R = 0, family = 'normal', verbose = FALSE)
	f$VC[c('va', 'cab', 'vb')] = c(0, .006543598, 0)
	s = als_start_vals(f)$Sab
	expect_true(all(eigen(s, symmetric = TRUE, only.values = TRUE)$values > 0))
	expect_equal(s[1, 2], 0)
})
