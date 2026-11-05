.scale_panel = function(sc) {
	set.seed(5101); n = 15; nm = sprintf("a%02d", 1:n)
	X = lapply(1:3, function(t) array(rnorm(n * n) * sc, c(n, n, 1), dimnames = list(nm, nm, "x")))
	Y = lapply(1:3, function(t) {
		y = 1 * (-0.5 + 0.8 * X[[t]][, , 1] / sc + matrix(rnorm(n * n), n) > 0)
		diag(y) = NA; dimnames(y) = list(nm, nm); y
	})
	list(Y = Y, X = X)
}

test_that("MCMC coefficient draws do not depend on covariate units", {
	skip_on_cran()
	fits = lapply(c(1e3, 1e5), function(sc) {
		d = .scale_panel(sc)
		f = suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			R = 1,
			burn = 50,
			nscan = 200,
			odens = 5,
			seed = 7,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))
		list(slope = f$BETA[, 2] * sc, intercept = f$BETA[, 1])
	})
	expect_equal(fits[[1]]$slope, fits[[2]]$slope, tolerance = 1e-6)
	expect_equal(fits[[1]]$intercept, fits[[2]]$intercept, tolerance = 1e-6)
})

test_that("the ALS route refuses covariates on a very large scale", {
	d = .scale_panel(1e5)
	expect_error(
		suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			verbose = FALSE
		)),
		class = "lame_als_covariate_scale"
	)
	expect_error(
		suppressMessages(ame_als(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = "binary",
			verbose = FALSE
		)),
		class = "lame_als_covariate_scale"
	)
})
