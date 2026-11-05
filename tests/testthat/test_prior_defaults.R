skip_on_cran()

test_that("persistence defaults are scoped to Gaussian bipartite AR1 fits", {
	set.seed(116)
	Y = array(
		rnorm(8 * 7 * 4),
		c(8, 7, 4),
		dimnames = list(paste0("a", 1:8), paste0("b", 1:7), paste0("t", 1:4))
	)
	args = list(
		Y = Y,
		family = "normal",
		mode = "bipartite",
		R = 1,
		dynamic_uv = TRUE,
		burn = 0,
		nscan = 8,
		odens = 1,
		seed = 116,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
	fit = do.call(lame, args)
	expect_equal(fit$prior$rho_uv_mean, .5)
	expect_equal(fit$prior$rho_uv_sd, .2)
	expect_identical(fit$refit_options$prior, list())
	args$prior = list(rho_uv_mean = .9, rho_uv_sd = .1)
	explicit = do.call(lame, args)
	expect_equal(explicit$prior$rho_uv_mean, .9)
	expect_equal(explicit$prior$rho_uv_sd, .1)
	args$prior = NULL
	args$family = "bin"
	args$Y = 1 * (Y > 0)
	binary = do.call(lame, args)
	expect_equal(binary$prior$rho_uv_mean, .9)
	expect_equal(binary$prior$rho_uv_sd, .1)
	args$family = "normal"
	args$mode = "unipartite"
	args$Y = array(
		rnorm(8 * 8 * 4),
		c(8, 8, 4),
		dimnames = list(paste0("a", 1:8), paste0("a", 1:8), paste0("t", 1:4))
	)
	for (tt in 1:4) diag(args$Y[, , tt]) = NA_real_
	directed = do.call(lame, args)
	expect_equal(directed$prior$rho_uv_mean, .9)
	expect_equal(directed$prior$rho_uv_sd, .1)
})
