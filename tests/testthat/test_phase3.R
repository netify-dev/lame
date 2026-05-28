# phase 3 — advanced inference + tidy integration.
#
# Items implemented in phase 3:
#   3.1 hierarchical (rho_b, sigma_b) pooling across blocks (per C2)
#   3.2 gof_temporal() PP temporal trend test
#   3.3 multivariate R-hat for AR(1) paths
#   3.4 time_index + period_exposure (mixed frequency)
#   3.5 tidy.lame + prediction_draws_long
#   3.6 max_seconds + checkpointing + lame_resume
#   3.7 exact rolling-origin LFO at last 3 periods (scoped down)

skip_on_cran()

.fit_phase3 <- function(seed = 7L, ...) {
	data(YX_bin_list, envir = environment())
	set.seed(seed)
	suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad",
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		seed = seed, ...))
}

# -- 3.1 hierarchical pooling ------------------------------------------------

test_that("dynamic_beta_pool='none' is the default and records on the fit", {
	fit <- .fit_phase3()
	expect_identical(fit$dynamic_beta_pool, "none")
})

test_that("dynamic_beta_pool != 'none' fits with hierarchical pooling", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad", dynamic_beta_pool = "both",
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_identical(fit$dynamic_beta_pool, "both")
	expect_true(all(is.finite(fit$BETA)))
})

# -- 3.2 gof_temporal --------------------------------------------------------

test_that("gof_temporal returns a structured result on a dynamic fit", {
	fit <- .fit_phase3()
	res <- suppressWarnings(gof_temporal(fit, stat = "density", n_rep = 30, seed = 42))
	expect_s3_class(res, "gof_temporal")
	expect_true(is.numeric(res$slope_obs))
	expect_length(res$slope_rep, 30L)
	expect_true(is.finite(res$p_pp) && res$p_pp >= 0 && res$p_pp <= 1)
})

test_that("gof_temporal supports reciprocity and transitivity", {
	fit <- .fit_phase3()
	r2 <- suppressWarnings(gof_temporal(fit, stat = "reciprocity",  n_rep = 20, seed = 42))
	r3 <- suppressWarnings(gof_temporal(fit, stat = "transitivity", n_rep = 20, seed = 42))
	expect_s3_class(r2, "gof_temporal")
	expect_s3_class(r3, "gof_temporal")
})

test_that("gof_temporal errors on T < 3", {
	fit <- .fit_phase3()
	bad <- fit
	bad$Y_obs <- fit$Y[, , 1:2]
	expect_error(suppressWarnings(gof_temporal(bad, stat = "density")),
	             "at least 3")
})

# -- 3.3 multivariate R-hat --------------------------------------------------

test_that("rhat_dynamic_beta returns a data frame with one row per coef", {
	fit_list <- lapply(c(1L, 2L), function(s) .fit_phase3(seed = s))
	res <- rhat_dynamic_beta(fit_list)
	expect_s3_class(res, "data.frame")
	expect_true(all(c("coef", "rhat_mvt", "rhat_max_univariate") %in% names(res)))
	expect_equal(nrow(res), dim(fit_list[[1]]$BETA)[2L])
	expect_true(all(is.finite(res$rhat_mvt)))
})

test_that("rhat_dynamic_beta requires at least two fits", {
	fit <- .fit_phase3()
	expect_error(rhat_dynamic_beta(list(fit)), "at least two")
})

test_that("rhat_dynamic_beta requires dynamic_beta", {
	data(YX_bin_list, envir = environment())
	fit_s <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE))
	expect_error(rhat_dynamic_beta(list(fit_s, fit_s)),
	             "3-D")
})

# -- 3.4 mixed frequency -----------------------------------------------------

test_that("time_index = NULL is the default", {
	fit <- .fit_phase3()
	expect_null(fit$time_index)
})

test_that("time_index strictly-increasing validation", {
	data(YX_bin_list, envir = environment())
	expect_error(suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		time_index = c(1, 1, 2, 3),  # non-strictly-increasing
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"strictly increasing")
})

test_that("time_index equal-gap is byte-identical to no time_index", {
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_eq <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		time_index = c(1, 2, 3, 4),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	set.seed(7); fit_no <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_true(identical(fit_eq$BETA, fit_no$BETA))
})

test_that("time_index with unequal gaps fits via gap-aware sampler", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad",
		time_index = c(1, 2, 4, 8),
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_identical(fit$time_index, c(1, 2, 4, 8))
	expect_true(all(is.finite(fit$BETA)))
})

test_that("period_exposure validation", {
	data(YX_bin_list, envir = environment())
	# negative exposure must abort (non-negative requirement)
	expect_error(suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		period_exposure = c(-1, 1, 1, 1),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"non-negative")
})

# -- 3.5 tidy.lame + prediction_draws_long -----------------------------------

test_that("tidy.lame returns data frame with expected columns for dynamic fit", {
	fit <- .fit_phase3()
	td <- tidy(fit)
	expect_s3_class(td, "data.frame")
	expect_true(all(c("term", "period", "estimate", "std.error",
	                  "statistic", "p.value", "conf.low", "conf.high") %in% names(td)))
	expect_equal(nrow(td), dim(fit$BETA)[2L] * dim(fit$BETA)[3L])
})

test_that("tidy.lame on static fit returns one row per coef (no period col)", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_s <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE))
	td_s <- tidy(fit_s)
	expect_equal(nrow(td_s), ncol(fit_s$BETA))
	expect_false("period" %in% names(td_s))
})

test_that("tidy.lame conf.int = FALSE drops CI columns", {
	fit <- .fit_phase3()
	td <- tidy(fit, conf.int = FALSE)
	expect_false("conf.low"  %in% names(td))
	expect_false("conf.high" %in% names(td))
})

test_that("prediction_draws_long returns the right shape", {
	fit <- .fit_phase3()
	pd <- prediction_draws_long(fit, n_draws = 3L, type = "link", seed = 42)
	expect_s3_class(pd, "data.frame")
	# columns follow tidybayes / marginaleffects convention
	# (.draw / .iteration / .value) plus actor names from fit$Y dimnames.
	expect_true(all(c(".chain", ".iteration", ".draw", "period",
	                  "period_label", "i", "j", "actor_i", "actor_j", ".value")
	                %in% names(pd)))
	expect_equal(nrow(pd), 3L * dim(fit$BETA)[3L] * 50L * 50L)
})

# -- 3.6 checkpointing -------------------------------------------------------

test_that("max_seconds aborts the chain cleanly", {
	data(YX_bin_list, envir = environment())
	ck <- tempfile(fileext = ".rds")
	set.seed(7)
	# scale-up nscan and use a tiny max_seconds so the chain definitely
	# hits the budget; on very fast hardware (or hot cache) a short chain
	# can finish before the wall-clock check, which would mask the test.
	fit <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 50000L, burn = 50L, odens = 25L,
		checkpoint_path = ck, checkpoint_every = 5L,
		max_seconds = 0.05,
		verbose = FALSE, plot = FALSE)))
	expect_true(isTRUE(fit$terminated_early),
	            info = sprintf("nscan=50000 with max_seconds=0.05 should terminate early; got terminated_early=%s",
	                           fit$terminated_early))
	expect_true(file.exists(ck))
})

test_that("lame_resume continues from a checkpoint", {
	data(YX_bin_list, envir = environment())
	ck <- tempfile(fileext = ".rds")
	set.seed(7)
	fit1 <- suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 200L, burn = 50L, odens = 5L,
		checkpoint_path = ck, checkpoint_every = 10L,
		max_seconds = 0.3,
		verbose = FALSE, plot = FALSE)))
	skip_if_not(isTRUE(fit1$terminated_early),
	            "first fit did not terminate early; resume not exercised")
	fit2 <- suppressWarnings(suppressMessages(
		lame_resume(ck, nscan_more = 30L)))
	expect_s3_class(fit2, "lame")
	expect_true(nrow(fit2$BETA) >= 1L)
})

test_that("lame_resume errors on missing checkpoint", {
	expect_error(lame_resume("/tmp/nonexistent_lame_checkpoint.rds"),
	             "not found")
})

# -- 3.7 LFO -----------------------------------------------------------------

test_that("lfo returns elpd_lfo and per_period structure", {
	fit <- .fit_phase3()
	res <- suppressWarnings(lfo(fit, periods = c(3L, 4L), refit = TRUE,
	                            nscan = 40L, burn = 10L, odens = 5L))
	expect_s3_class(res, "lfo_lame")
	expect_true(is.finite(res$elpd_lfo))
	expect_equal(nrow(res$per_period), 2L)
	expect_true(all(c("period", "elpd", "n_obs") %in% names(res$per_period)))
})

test_that("lfo errors when periods < 2", {
	fit <- .fit_phase3()
	expect_error(suppressWarnings(lfo(fit, periods = 1L, refit = FALSE)),
	             "leave-out period")
})
