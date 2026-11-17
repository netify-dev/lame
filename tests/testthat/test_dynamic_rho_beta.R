# coefficient paths under a dyadic correlation that varies by period
skip_on_cran()

test_that("per-period correlations give the stacked per-period statistics and draws", {
	set.seed(1)
	n = 7L; nt = 4L; pd = 2L
	Xd = lapply(1:nt, function(t) matrix(rnorm(n * n * pd), n * n, pd))
	Xs = lapply(1:nt, function(t) matrix(rnorm(n * n), n * n, 1L))
	Z = lapply(1:nt, function(t) { z = matrix(rnorm(n * n), n); z[2, 5] = NA; z })
	off = lapply(1:nt, function(t) matrix(rnorm(n * n, 0, .2), n))
	rho = c(-.6, .1, .45, .8)
	stats = lame:::beta_dynamic_stats_cpp(Xd, Xs, Z, off, .3, 1.4, rho, FALSE, FALSE, TRUE)
	stacked = list(Omega = array(0, c(pd, pd, nt)), h = matrix(0, nt, pd))
	for (t in 1:nt) {
		one = lame:::beta_dynamic_stats_cpp(Xd[t], Xs[t], Z[t], off[t], .3, 1.4, rho[t], FALSE, FALSE, TRUE)
		stacked$Omega[, , t] = one$Omega[, , 1]
		stacked$h[t, ] = one$h[1, ]
	}
	expect_equal(stats$Omega, stacked$Omega, tolerance = 1e-12)
	expect_equal(stats$h, stacked$h, tolerance = 1e-12)
	args = list(Xd, Xs, Z, off, .3, c(1, .7), c(.3, .2), diag(2), c(0, 0), diag(10, 2), 1.4)
	set.seed(9)
	a = do.call(lame:::sample_beta_dynamic_cpp, c(args, list(rho, FALSE, FALSE, TRUE)))$path
	set.seed(9)
	b = do.call(
		lame:::sample_beta_dynamic_cpp,
		c(args, list(0, FALSE, FALSE, TRUE, stacked$Omega, stacked$h))
	)$path
	expect_identical(a, b)
	# one value per period equal to a scalar is the scalar filter
	set.seed(4)
	a = do.call(lame:::sample_beta_dynamic_cpp, c(args, list(rep(.3, nt), FALSE, FALSE, TRUE)))$path
	set.seed(4)
	b = do.call(lame:::sample_beta_dynamic_cpp, c(args, list(.3, FALSE, FALSE, TRUE)))$path
	expect_identical(a, b)
	expect_error(
		lame:::beta_dynamic_stats_cpp(Xd, Xs, Z, off, .3, 1.4, rho[1:3], FALSE, FALSE, TRUE),
		"one value per period"
	)
})

test_that("dynamic_rho takes the dynamic_ab intercept path and user coefficient blocks", {
	set.seed(21)
	n = 12L; nt = 4L
	nm = sprintf("a%02d", seq_len(n))
	shift = c(-1, -.3, .3, 1)
	X = lapply(seq_len(nt), function(t) array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x")))
	Y = lapply(seq_len(nt), function(t) {
		e = matrix(rnorm(n * n), n); e = e + .5 * t(e)
		m = shift[t] + .5 * X[[t]][, , 1] + outer(rnorm(n, 0, .3), rnorm(n, 0, .3), "+") + e
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	go = function(...) suppressWarnings(suppressMessages(lame(
		Y,
		Xdyad = X,
		R = 0,
		family = "normal",
		dynamic_rho = TRUE,
		burn = 100,
		nscan = 400,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		...
	)))
	f = go(dynamic_ab = TRUE)
	expect_true(isTRUE(f$dynamic_intercept_from_ab))
	expect_length(dim(f$BETA), 3L)
	expect_identical(f$tryErrorChecks$beta, 0)
	expect_gt(stats::cor(apply(f$BETA[, 1, ], 2, mean), shift), .9)
	expect_identical(dim(f$RHO), c(80L, nt))
	g = go(dynamic_beta = "dyad")
	expect_identical(g$dynamic_beta_kind_by_block, c(dyad = "ar1"))
	expect_identical(g$tryErrorChecks$beta, 0)
	expect_identical(g$tryErrorChecks$betaAB, 0)
})
