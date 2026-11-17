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
	expect_error(lame:::loo.ame_als(fit), "no per-iteration log-likelihood")
	expect_error(lame:::waic.ame_als(fit), "save_log_lik")
	skip_if_not_installed("loo")
	expect_error(loo::loo(fit), "no per-iteration log-likelihood")
	expect_error(loo::waic(fit), "no per-iteration log-likelihood")
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

test_that("lame_resume() names a saved data object the session no longer holds", {
	skip_on_cran()
	set.seed(4)
	mk = function() lapply(1:3, function(t) {
		y = matrix(rnorm(64), 8, 8, dimnames = list(paste0("a", 1:8), paste0("a", 1:8)))
		diag(y) = NA
		y
	})
	ck = tempfile(fileext = ".rds")
	env = new.env()
	env$panel_gone = mk()
	fit = suppressWarnings(suppressMessages(eval(
		quote(lame(panel_gone, family = "normal", R = 0, nscan = 20, burn = 5, odens = 1, verbose = FALSE, checkpoint_path = ck)),
		env
	)))
	kept = env$panel_gone
	rm("panel_gone", envir = env)
	e = expect_error(lame_resume(fit, nscan_more = 10, .envir = env), class = "lame_resume_missing_object")
	expect_identical(e$objects, "panel_gone")
	expect_match(conditionMessage(e), "panel_gone")
	if (file.exists(ck)) {
		e = expect_error(lame_resume(ck, nscan_more = 10, .envir = env), class = "lame_resume_missing_object")
		expect_identical(e$objects, "panel_gone")
		by_file = suppressWarnings(suppressMessages(lame_resume(ck, nscan_more = 10, Y = kept, .envir = env)))
		expect_s3_class(by_file, "lame")
	}
	# data passed by name replace the saved argument
	again = suppressWarnings(suppressMessages(lame_resume(fit, nscan_more = 10, Y = kept, .envir = env)))
	expect_s3_class(again, "lame")
	expect_identical(nrow(again$BETA), 10L)
	unlink(ck)
})

test_that("print() marks the bipartite innovation sd as set by the prior", {
	skip_on_cran()
	set.seed(9)
	rn = paste0("r", 1:6); cn = paste0("c", 1:7)
	Y = lapply(1:3, function(t) matrix(rnorm(42), 6, 7, dimnames = list(rn, cn)))
	fit = suppressWarnings(suppressMessages(lame(
		Y, family = "normal", mode = "bipartite", R_row = 1, R_col = 1, dynamic_uv = TRUE,
		nscan = 20, burn = 5, odens = 1, verbose = FALSE
	)))
	out = paste(c(capture.output(print(fit)), capture.output(print(fit), type = "message")), collapse = "\n")
	expect_match(out, "set by the prior; G carries the scale", fixed = TRUE)
})
