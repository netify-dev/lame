test_that("checkpoints freeze priors independently of the original expression", {
	set.seed(114)
	Y = array(
		rnorm(8 * 7 * 4),
		c(8, 7, 4),
		dimnames = list(paste0("a", 1:8), paste0("b", 1:7), paste0("t", 1:4))
	)
	requested = list(rho_uv_mean = .2, rho_uv_sd = .15)
	ck = tempfile(fileext = ".rds")
	withr::defer(unlink(ck))
	original = suppressWarnings(suppressMessages(lame(
		Y,
		family = "normal",
		mode = "bipartite",
		R = 1,
		dynamic_uv = TRUE,
		prior = requested,
		burn = 0,
		nscan = 8,
		odens = 1,
		seed = 114,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		checkpoint_path = ck,
		max_seconds = 1e-12
	)))
	rm(requested)
	saved = readRDS(ck)
	expect_identical(saved$resolved_prior, original$prior)
	resumed = suppressMessages(lame_resume(ck, nscan_more = 8))
	reference = suppressMessages(lame_resume(
		ck,
		nscan_more = 8,
		prior = list(rho_uv_mean = .2, rho_uv_sd = .15)
	))
	for (nm in c("BETA", "VC", "U", "V", "G", "rho_uv", "sigma_uv", "start_vals")) {
		expect_identical(resumed[[nm]], reference[[nm]], info = nm)
	}
	expect_identical(resumed$prior, original$prior)
	expect_identical(resumed$prior_supplied, original$prior_supplied)
	expect_identical(resumed$refit_options$prior, original$refit_options$prior)
})

test_that("legacy Gaussian bipartite checkpoints preserve their original defaults", {
	set.seed(115)
	Y = array(
		rnorm(8 * 7 * 4),
		c(8, 7, 4),
		dimnames = list(paste0("a", 1:8), paste0("b", 1:7), paste0("t", 1:4))
	)
	ck = tempfile(fileext = ".rds")
	withr::defer(unlink(ck))
	suppressWarnings(suppressMessages(lame(
		Y,
		family = "normal",
		mode = "bipartite",
		R = 1,
		dynamic_uv = TRUE,
		prior = list(rho_uv_mean = .9, rho_uv_sd = .1),
		burn = 0,
		nscan = 8,
		odens = 1,
		seed = 115,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		checkpoint_path = ck,
		max_seconds = 1e-12
	)))
	legacy = readRDS(ck)
	legacy$resolved_prior = legacy$requested_prior = NULL
	legacy$call_args$prior = NULL
	saveRDS(legacy, ck)
	resumed = suppressMessages(lame_resume(ck, nscan_more = 8))
	reference = suppressMessages(lame_resume(
		ck,
		nscan_more = 8,
		prior = list(rho_uv_mean = .9, rho_uv_sd = .1)
	))
	expect_equal(resumed$prior$rho_uv_mean, .9)
	expect_equal(resumed$prior$rho_uv_sd, .1)
	expect_identical(resumed$refit_options$prior, list())
	expect_false(resumed$prior_supplied)
	for (nm in c("BETA", "VC", "U", "V", "G", "rho_uv", "sigma_uv", "start_vals")) {
		expect_identical(resumed[[nm]], reference[[nm]], info = nm)
	}
})
