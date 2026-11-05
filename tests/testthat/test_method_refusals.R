# generics defined for a fit class only to refuse, and nobs() of a fast fit; run on CRAN
skip_on_cran()

test_that("formula() and logLik() refuse an MCMC lame fit with the reason", {
	fit = structure(list(), class = "lame")
	expect_error(formula(fit), "formula interface")
	expect_error(logLik(fit), "point log-likelihood is not stored")
})

test_that("prior_summary() has no default method and says for which class", {
	expect_error(prior_summary(structure(list(), class = "some_fit")), "some_fit")
	expect_error(prior_summary(1), "numeric")
})

test_that("fitted() and residuals() of a bootstrap point to the underlying fit", {
	b = structure(list(), class = "boot_ame")
	expect_error(fitted(b), "not defined for a")
	expect_error(residuals(b), "underlying")
})

test_that("loo() and waic() of a fast fit refuse and name the MCMC route", {
	fit = structure(list(), class = c("ame_als", "list"))
	expect_error(lame:::loo.ame_als(fit), "point estimator")
	expect_error(lame:::waic.ame_als(fit), "save_log_lik")
	skip_if_not_installed("loo")
	expect_error(loo::loo(fit), "point estimator")
	expect_error(loo::waic(fit), "point estimator")
})

test_that("nobs() of a fast fit counts the finite cells of the stored response", {
	y = matrix(c(NA, 1, 0, 1, NA, 1, 0, Inf, NA), 3)
	expect_identical(nobs(structure(list(Y = y), class = "ame_als")), 5L)
	expect_identical(nobs(structure(list(Y = list(y, y)), class = "ame_als")), 10L)
	expect_error(nobs(structure(list(), class = "ame_als")), "not stored")
	set.seed(3); n = 8
	Y = matrix(rnorm(n * n), n); diag(Y) = NA
	fit = suppressWarnings(suppressMessages(ame_als(Y, family = "normal", R = 0, verbose = FALSE)))
	expect_identical(nobs(fit), as.integer(n * (n - 1)))
})

test_that("lame_resume() points a missing base checkpoint to the per-chain files", {
	d = tempfile("ck"); dir.create(d)
	base = file.path(d, "run.rds")
	expect_error(lame_resume(base), "not found")
	saveRDS(list(), file.path(d, "run_chain1.rds")); saveRDS(list(), file.path(d, "run_chain2.rds"))
	e = expect_error(lame_resume(base), class = "lame_resume_chain_files")
	expect_identical(basename(e$files), c("run_chain1.rds", "run_chain2.rds"))
	unlink(d, recursive = TRUE)
})
