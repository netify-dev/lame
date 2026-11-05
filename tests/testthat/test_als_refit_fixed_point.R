# a converged empirical-Bayes fit is a fixed point of ame_als_refit(): the count noise variance
# stays at the value the fit found, and the iterations reported are the refit's own

test_that("a refit of a converged static fit on its own data stays at its estimates", {
	skip_on_cran()
	set.seed(12)
	n = 30; nm = paste0("a", 1:n)
	x = matrix(rnorm(n * n), n); dimnames(x) = list(nm, nm)
	a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
	eta = -0.5 + 0.6 * x + outer(a, b, "+")
	Yb = 1 * (eta + matrix(rnorm(n * n), n) > 0); diag(Yb) = NA; dimnames(Yb) = list(nm, nm)
	Yp = matrix(rpois(n * n, exp(eta - 0.5 + rnorm(n * n, 0, .3))), n); diag(Yp) = NA; dimnames(Yp) = list(nm, nm)
	for (fam in c("binary", "poisson")) {
		Y = if (fam == "binary") Yb else Yp
		for (R in c(0, 1)) {
			f = suppressWarnings(suppressMessages(ame_als(Y, Xdyad = x, R = R, family = fam, verbose = FALSE)))
			r = suppressWarnings(suppressMessages(ame_als_refit(f)))
			info = paste(fam, "R =", R)
			expect_lt(max(abs(coef(f) - coef(r))), 1e-4, label = info)
			if (fam == "poisson") expect_identical(r$eb$ve, f$eb$ve, label = info)
			# the refit starts at the variational state, so it runs no more iterations than the fit
			expect_lte(r$eb$variational$iterations, f$eb$variational$iterations, label = info)
			# a second refit does not move it either
			r2 = suppressWarnings(suppressMessages(ame_als_refit(r)))
			expect_lt(max(abs(coef(r) - coef(r2))), 1e-4, label = info)
			# the iterations are the refit's own reweighting steps
			if (fam == "binary") expect_identical(r$iterations, 0L, label = info)
		}
	}
})

test_that("a refit on other data searches the noise variance again", {
	skip_on_cran()
	set.seed(31)
	n = 30; nm = paste0("a", 1:n)
	x = matrix(rnorm(n * n), n); dimnames(x) = list(nm, nm)
	a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
	eta = -0.5 + 0.6 * x + outer(a, b, "+")
	mk = function(sd_e) {
		Y = matrix(rpois(n * n, exp(eta + rnorm(n * n, 0, sd_e))), n); diag(Y) = NA; dimnames(Y) = list(nm, nm); Y
	}
	f = suppressWarnings(suppressMessages(ame_als(mk(.3), Xdyad = x, R = 0, family = "poisson", verbose = FALSE)))
	r = suppressWarnings(suppressMessages(ame_als_refit(f, Y_new = mk(.9))))
	# noisier data pull ve up
	expect_gt(r$eb$ve, f$eb$ve)
})
