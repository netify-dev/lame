# level step of a directed binary dynamic_ab fit with an rw1 intercept path and dyadic correlation
skip_on_cran()

test_that("the level step keeps the exact posterior with an rw1 intercept path and dyadic correlation", {
	cases = list(
		list(variant = 'general', blocks = 'rw1', seed = 21L),
		list(variant = 'general', blocks = 'mixed', seed = 22L),
		list(variant = 'path', blocks = 'mixed', shrink = TRUE, seed = 24L),
		list(variant = 'path', blocks = 'rw1', shrink = FALSE, seed = 25L)
	)
	for (cc in cases) {
		r = do.call(lvl_exact, c(cc, list(rho_d = .4, ndraw = 3e5)))
		expect_gt(r$ess, 500)
		expect_lt(r$maxz, 4.6)
	}
})

test_that("the level-step check detects a mis-stated dimension of the coefficient block", {
	expect_gt(lvl_exact('general', 'rw1', rho_d = .4, ndraw = 3e5, seed = 27L, bad_dim = 2L)$maxz, 5)
})
