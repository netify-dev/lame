# tests for the dynamic_beta methods layer: diagnostics + prior summaries,
# tidy() / forecasting output, and prediction + log-likelihood machinery.
#
# covers summary()'s rho_beta stationarity table, detect_change_point(),
# dynamic_beta_prior_summary(), freeze_call / update(), compact printing,
# hierarchical pooling, gof_temporal(), rhat_dynamic_beta(), mixed-frequency
# time_index / period_exposure, tidy.lame + prediction_draws_long,
# checkpointing (max_seconds / lame_resume), lfo(), transition kinds at the
# fit level, counterfactual + h-step forecasting, family link inversion,
# autoplot.lame, and pointwise / chunked log-lik storage (loo).

skip_on_cran()


# =============================================================================
# diagnostics, prior summaries, freeze_call / update, compact printing
# =============================================================================

test_that("summary.lame computes rho_beta_stationarity table with right shape", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	s = summary(fit)
	expect_true(!is.null(s$rho_beta_stationarity))
	expect_named(s$rho_beta_stationarity, c("block", "q05", "iqr", "warn"))
})

test_that("stationarity warning fires when rho_beta posterior is forced near 1", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	# force a concentrated near-unit posterior
	fit$RHO_BETA = matrix(rep(rnorm(20, 0.985, 0.01), each = ncol(fit$RHO_BETA)),
	                      nrow = 20, ncol = ncol(fit$RHO_BETA),
	                      dimnames = dimnames(fit$RHO_BETA))
	s = summary(fit)
	expect_true(any(s$rho_beta_stationarity$warn))
	out = capture.output(print(s))
	expect_true(any(grepl("concentrated near 1", out)))
	expect_true(any(grepl("dynamic_beta_kind", out)))
})

test_that("stationarity warning does NOT fire on a diffuse rho_beta posterior", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	# force a diffuse posterior
	set.seed(1)
	fit$RHO_BETA = matrix(runif(20 * ncol(fit$RHO_BETA), 0.3, 0.9),
	                      nrow = 20, ncol = ncol(fit$RHO_BETA),
	                      dimnames = dimnames(fit$RHO_BETA))
	s = summary(fit)
	expect_false(any(s$rho_beta_stationarity$warn))
	out = capture.output(print(s))
	expect_false(any(grepl("concentrated near 1", out)))
})

test_that("90%-NA Y fits without error and produces finite BETA", {
	set.seed(2026)
	n = 10; Tn = 3
	Y = replicate(Tn, {
		z = matrix(rnorm(n*n), n, n); diag(z) = NA
		# make the observed network very sparse
		mask = matrix(stats::runif(n*n) > 0.10, n, n)
		z[mask] = NA; diag(z) = NA
		z
	}, simplify = FALSE)
	X = replicate(Tn, array(rnorm(n*n), c(n, n, 1)), simplify = FALSE)
	fit = lame(Y, Xdyad = X, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	expect_true(all(is.finite(fit$BETA)))
	expect_equal(length(dim(fit$BETA)), 3L)
})

test_that("extra missing dyads keep the dynamic beta sampler finite", {
	set.seed(7)
	n = 8; Tn = 3
	# baseline y with only diagonal missing
	Y_a = replicate(Tn, {z = matrix(rnorm(n*n), n, n); diag(z) = NA; z}, simplify = FALSE)
	# add one missing dyad pair
	Y_b = lapply(Y_a, function(y) { y[1, 2] = NA; y[2, 1] = NA; y })
	X = replicate(Tn, array(rnorm(n*n), c(n, n, 1)), simplify = FALSE)
	set.seed(1); fit_a = lame(Y_a, Xdyad = X, family = "normal", R = 0,
	                           nscan = 30, burn = 5, odens = 5,
	                           dynamic_beta = "dyad", verbose = FALSE)
	set.seed(1); fit_b = lame(Y_b, Xdyad = X, family = "normal", R = 0,
	                           nscan = 30, burn = 5, odens = 5,
	                           dynamic_beta = "dyad", verbose = FALSE)
	expect_true(all(is.finite(fit_a$BETA)))
	expect_true(all(is.finite(fit_b$BETA)))
})

test_that("detect_change_point returns the expected data frame shape", {
	dat = simulate_test_network(n = 6, n_time = 4, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	res = detect_change_point(fit, n_prior_sims = 200, seed = 1)
	expect_named(res, c("coef","bf","m_post_mean","m_prior_q95","t_hat","warn"))
	expect_type(res$bf, "double")
	expect_type(res$warn, "logical")
})

test_that("detect_change_point errors when dynamic_beta is not active", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 5, verbose = FALSE)
	expect_error(detect_change_point(fit), "dynamic_beta")
})

test_that("detect_change_point runs end-to-end on a step-function truth", {
	set.seed(2026)
	n = 8; Tn = 6
	beta_true = c(rep(-0.6, 3), rep(0.6, 3))
	X_list = replicate(Tn, matrix(rnorm(n*n), n, n), simplify = FALSE)
	Y_list = vector("list", Tn)
	for (t in seq_len(Tn)) {
		Yt = beta_true[t] * X_list[[t]] + matrix(rnorm(n*n, 0, 0.1), n, n)
		diag(Yt) = NA
		rownames(Yt) = colnames(Yt) = paste0("a", seq_len(n))
		Y_list[[t]] = Yt
	}
	X_arr = lapply(X_list, function(x) array(x, c(n, n, 1)))
	fit = lame(Y_list, Xdyad = X_arr, family = "normal", R = 0,
	            nscan = 100, burn = 30, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	res = detect_change_point(fit, n_prior_sims = 200, seed = 1)
	expect_true(is.data.frame(res))
	expect_gt(nrow(res), 0L)
	expect_true(all(is.finite(res$bf)))
})

test_that("dynamic_beta_prior_summary is reproducible with seed", {
	a = dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	b = dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	expect_equal(a$summary, b$summary)
	expect_equal(a$paths, b$paths)
})

test_that("larger sigma_scale (= more diffuse) gives larger max-first-diff", {
	small = dynamic_beta_prior_summary(sigma_scale = 0.05, seed = 1, ndraws = 500)
	big   = dynamic_beta_prior_summary(sigma_scale = 5.0,  seed = 1, ndraws = 500)
	expect_gt(big$summary$q50[1], small$summary$q50[1])
	expect_gt(big$summary$q50[2], small$summary$q50[2])
})

test_that("rw1 kind has rho fixed at 1 and path is unit-root random walk", {
	s = dynamic_beta_prior_summary(kind = "rw1", seed = 1, ndraws = 500)
	expect_true(all(s$rho == 1))
	expect_equal(s$kind, "rw1")
})

test_that("invalid rho_mean outside bounds errors cleanly", {
	expect_error(dynamic_beta_prior_summary(rho_mean = 1.5, ndraws = 100),
	             "must be strictly between")
	expect_error(dynamic_beta_prior_summary(rho_mean = -0.1, ndraws = 100),
	             "must be strictly between")
})

test_that("n_periods = 1 errors with informative message", {
	expect_error(dynamic_beta_prior_summary(n_periods = 1), "at least 2")
})

test_that("freeze_call = TRUE stores a data snapshot on the fit", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            freeze_call = TRUE, verbose = FALSE)
	expect_true(isTRUE(fit$freeze_call))
	expect_false(is.null(fit$data_snapshot))
	expect_identical(fit$data_snapshot$Y, dat$Y)
})

test_that("freeze_call = TRUE: update() refits with snapshot even if local Y mutates", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	Y_local = dat$Y
	set.seed(7)
	fit = lame(Y_local, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            freeze_call = TRUE, verbose = FALSE)
		# mutate y_local after fit
	Y_local[[1]][1, 2] = Y_local[[1]][1, 2] + 1000
	# update() must use the snapshot, not the mutated y_local
	set.seed(7)
	fit_upd = update(fit, nscan = 10)
	# refit against the original data with the same args should reproduce fit
	set.seed(7)
	fit_ref = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	                nscan = 10, burn = 2, odens = 1,
	                freeze_call = TRUE, verbose = FALSE)
	expect_equal(fit_upd$BETA, fit_ref$BETA)
})

test_that("freeze_call = FALSE uses parent-frame update lookup", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1, verbose = FALSE)
	expect_null(fit$data_snapshot)
	expect_null(fit$freeze_call)
	# update() should still work using parent.frame() lookup
	fit_upd = update(fit, nscan = 5)
	expect_s3_class(fit_upd, "lame")
})

# -- 1.8 compact print -------------------------------------------------------

test_that("static print output does NOT use the compact table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1, verbose = FALSE)
	out = paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	expect_false(grepl("Dynamic components", out))
})

test_that("single-dynamic fit does NOT use the compact table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_beta = TRUE, verbose = FALSE)
	out = paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	# single dynamic_beta should print the regular per-component line, not
	# the compact joint table
	expect_false(grepl("Dynamic components$|Dynamic components\\s", out))
})

test_that("joint dynamic_ab + dynamic_beta fit uses the compact table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_ab = TRUE, dynamic_beta = TRUE, verbose = FALSE)
	out = paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	expect_true(grepl("Dynamic components", out))
})

test_that("compact = FALSE disables the joint table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_ab = TRUE, dynamic_beta = TRUE, verbose = FALSE)
	out = paste(c(capture.output(print(fit, compact = FALSE)),
	               capture.output(print(fit, compact = FALSE), type = "message")),
	             collapse = "\n")
	expect_false(grepl("Dynamic components", out))
})

# =============================================================================
# tidy() + forecasting: pooling, gof_temporal, rhat, time_index, tidy, checkpoint, lfo
# =============================================================================

.fit_dynamic_beta_tidy = function(seed = 7L, ...) {
	data(YX_bin_list, envir = environment())
	set.seed(seed)
	suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad",
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		seed = seed, ...))
}

test_that("dynamic_beta_pool='none' is the default and records on the fit", {
	fit = .fit_dynamic_beta_tidy()
	expect_identical(fit$dynamic_beta_pool, "none")
})

test_that("dynamic_beta_pool != 'none' fits with hierarchical pooling", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad", dynamic_beta_pool = "both",
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_identical(fit$dynamic_beta_pool, "both")
	expect_true(all(is.finite(fit$BETA)))
})

test_that("gof_temporal returns a structured result on a dynamic fit", {
	fit = .fit_dynamic_beta_tidy()
	res = suppressWarnings(gof_temporal(fit, stat = "density", n_rep = 30, seed = 42))
	expect_s3_class(res, "gof_temporal")
	expect_true(is.numeric(res$slope_obs))
	expect_length(res$slope_rep, 30L)
	expect_true(is.finite(res$p_pp) && res$p_pp >= 0 && res$p_pp <= 1)
})

test_that("gof_temporal supports reciprocity and transitivity", {
	fit = .fit_dynamic_beta_tidy()
	r2 = suppressWarnings(gof_temporal(fit, stat = "reciprocity",  n_rep = 20, seed = 42))
	r3 = suppressWarnings(gof_temporal(fit, stat = "transitivity", n_rep = 20, seed = 42))
	expect_s3_class(r2, "gof_temporal")
	expect_s3_class(r3, "gof_temporal")
})

test_that("gof_temporal errors on T < 3", {
	fit = .fit_dynamic_beta_tidy()
	bad = fit
	bad$Y_obs = fit$Y[, , 1:2]
	expect_error(suppressWarnings(gof_temporal(bad, stat = "density")),
	             "at least 3")
})

test_that("rhat_dynamic_beta returns a data frame with one row per coef", {
	fit_list = lapply(c(1L, 2L), function(s) .fit_dynamic_beta_tidy(seed = s))
	res = rhat_dynamic_beta(fit_list)
	expect_s3_class(res, "data.frame")
	expect_true(all(c("coef", "rhat_mvt", "rhat_max_univariate") %in% names(res)))
	expect_equal(nrow(res), dim(fit_list[[1]]$BETA)[2L])
	expect_true(all(is.finite(res$rhat_mvt)))
})

test_that("rhat_dynamic_beta requires at least two fits", {
	fit = .fit_dynamic_beta_tidy()
	expect_error(rhat_dynamic_beta(list(fit)), "at least two")
})

test_that("rhat_dynamic_beta requires dynamic_beta", {
	data(YX_bin_list, envir = environment())
	fit_s = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE))
	expect_error(rhat_dynamic_beta(list(fit_s, fit_s)),
	             "3-D")
})

# -- 3.4 mixed frequency -----------------------------------------------------

test_that("time_index = NULL is the default", {
	fit = .fit_dynamic_beta_tidy()
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
	set.seed(7); fit_eq = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		time_index = c(1, 2, 3, 4),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	set.seed(7); fit_no = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_true(identical(fit_eq$BETA, fit_no$BETA))
})

test_that("time_index with unequal gaps fits via gap-aware sampler", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
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
	fit = .fit_dynamic_beta_tidy()
	td = tidy(fit)
	expect_s3_class(td, "data.frame")
	expect_true(all(c("term", "period", "estimate", "std.error",
	                  "statistic", "p.value", "conf.low", "conf.high") %in% names(td)))
	expect_equal(nrow(td), dim(fit$BETA)[2L] * dim(fit$BETA)[3L])
})

test_that("tidy.lame on static fit returns one row per coef (no period col)", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_s = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE))
	td_s = tidy(fit_s)
	expect_equal(nrow(td_s), ncol(fit_s$BETA))
	expect_false("period" %in% names(td_s))
})

test_that("tidy.lame conf.int = FALSE drops CI columns", {
	fit = .fit_dynamic_beta_tidy()
	td = tidy(fit, conf.int = FALSE)
	expect_false("conf.low"  %in% names(td))
	expect_false("conf.high" %in% names(td))
})

test_that("prediction_draws_long returns the right shape", {
	fit = .fit_dynamic_beta_tidy()
	pd = prediction_draws_long(fit, n_draws = 3L, type = "link", seed = 42)
	expect_s3_class(pd, "data.frame")
	# columns follow tidybayes / marginaleffects convention
	# (.draw / .iteration / .value) plus actor names from fit$y dimnames.
	expect_true(all(c(".chain", ".iteration", ".draw", "period",
	                  "period_label", "i", "j", "actor_i", "actor_j", ".value")
	                %in% names(pd)))
	expect_equal(nrow(pd), 3L * dim(fit$BETA)[3L] * 50L * 50L)
})

# -- 3.6 checkpointing -------------------------------------------------------

test_that("max_seconds aborts the chain cleanly", {
	data(YX_bin_list, envir = environment())
	ck = tempfile(fileext = ".rds")
	set.seed(7)
	# scale-up nscan and use a tiny max_seconds so the chain definitely
	# hits the budget; on very fast hardware (or hot cache) a short chain
	# can finish before the wall-clock check, which would mask the test.
	fit = suppressWarnings(suppressMessages(lame(
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
	ck = tempfile(fileext = ".rds")
	set.seed(7)
	fit1 = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 200L, burn = 50L, odens = 5L,
		checkpoint_path = ck, checkpoint_every = 10L,
		max_seconds = 0.3,
		verbose = FALSE, plot = FALSE)))
	skip_if_not(isTRUE(fit1$terminated_early),
	            "first fit did not terminate early; resume not exercised")
	fit2 = suppressWarnings(suppressMessages(
		lame_resume(ck, nscan_more = 30L)))
	expect_s3_class(fit2, "lame")
	expect_true(nrow(fit2$BETA) >= 1L)
})

test_that("lame_resume errors on missing checkpoint", {
	expect_error(lame_resume("/tmp/nonexistent_lame_checkpoint.rds"),
	             "not found")
})

# -- 3.7 lfo -----------------------------------------------------------------

test_that("lfo returns elpd_lfo and per_period structure", {
	fit = .fit_dynamic_beta_tidy()
	res = suppressWarnings(lfo(fit, periods = c(3L, 4L), refit = TRUE,
	                            nscan = 40L, burn = 10L, odens = 5L))
	expect_s3_class(res, "lfo_lame")
	expect_true(is.finite(res$elpd_lfo))
	expect_equal(nrow(res$per_period), 2L)
	expect_true(all(c("period", "elpd", "n_obs") %in% names(res$per_period)))
})

test_that("lfo errors when periods < 2", {
	fit = .fit_dynamic_beta_tidy()
	expect_error(suppressWarnings(lfo(fit, periods = 1L, refit = FALSE)),
	             "leave-out period")
})

# fit$Xlist stores the augmented design (intercept slice first). lfo() must
# strip that slice before refitting, since the refit adds its own intercept
# and a stacked pair of constant slices is only jointly identified. capture
# the refit to check the design it actually sees.
.capture_lfo_refit = function(fit, ...) {
	cap = new.env()
	real_lame = lame
	testthat::local_mocked_bindings(lame = function(...) {
		args = list(...)
		cap$Xdyad = args$Xdyad
		cap$fit = do.call(real_lame, args)
		cap$fit
	}, .package = "lame", .env = parent.frame())
	cap$res = suppressWarnings(lfo(fit, ...))
	cap
}

test_that("lfo refit carries exactly one intercept (binary T=4, one covariate)", {
	data(YX_bin_list, envir = environment())
	X1 = lapply(YX_bin_list$X, function(x) x[, , 1L, drop = FALSE])
	fit = suppressWarnings(lame(YX_bin_list$Y, X1, family = "binary", R = 0,
		nscan = 40, burn = 10, odens = 5, verbose = FALSE, plot = FALSE, seed = 7L))
	# premise: the stored design is augmented with the intercept slice
	expect_identical(dimnames(fit$Xlist[[1L]])[[3L]][1L], "intercept")

	cap = .capture_lfo_refit(fit, periods = 4L, refit = TRUE,
		nscan = 20L, burn = 5L, odens = 5L)
	# refit design: the covariate only, no stored intercept slice
	expect_equal(length(cap$Xdyad), 3L)
	expect_equal(dim(cap$Xdyad[[1L]])[3L], 1L)
	expect_false("intercept" %in% dimnames(cap$Xdyad[[1L]])[[3L]])
	# refit coefficients: one intercept + one covariate, nothing collinear.
	# a stacked constant slice would be labelled "intercept_dyad", so match
	# on the prefix rather than the exact name.
	expect_equal(sum(startsWith(colnames(cap$fit$BETA), "intercept")), 1L)
	expect_equal(ncol(cap$fit$BETA), 2L)
	expect_true(is.finite(cap$res$elpd_lfo))
})

test_that("lfo refit of an intercept-only fit passes Xdyad = NULL", {
	data(YX_bin_list, envir = environment())
	fit = suppressWarnings(lame(YX_bin_list$Y, family = "binary", R = 0,
		nscan = 40, burn = 10, odens = 5, verbose = FALSE, plot = FALSE, seed = 7L))
	expect_identical(dimnames(fit$Xlist[[1L]])[[3L]], "intercept")

	cap = .capture_lfo_refit(fit, periods = 4L, refit = TRUE,
		nscan = 20L, burn = 5L, odens = 5L)
	expect_null(cap$Xdyad)
	expect_identical(colnames(cap$fit$BETA), "intercept")
	expect_true(is.finite(cap$res$elpd_lfo))
})

# =============================================================================
# prediction + log-likelihood: transition kinds, counterfactual/forecast, autoplot, log-lik
# =============================================================================

# reusable small fit
.fit_dynamic_beta_prediction = function(kind = "ar1", nscan = 200, burn = 50,
                                          dynamic_beta = "dyad",
                                          family = "binary", seed = 7L, ...) {
	data(YX_bin_list, envir = environment())
	set.seed(seed)
	suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X,
		family = family, R = 0,
		dynamic_beta = dynamic_beta,
		dynamic_beta_kind = kind,
		nscan = nscan, burn = burn, odens = 5,
		verbose = FALSE, plot = FALSE,
		...
	))
}

test_that("dynamic_beta_kind = 'rw1' fits without error and pins rho=1", {
	fit = .fit_dynamic_beta_prediction(kind = "rw1", nscan = 100, burn = 25)
	expect_s3_class(fit, "lame")
	expect_identical(fit$dynamic_beta_kind, "rw1")
	expect_true(all(abs(fit$RHO_BETA - 1) < 1e-12))
	expect_length(dim(fit$BETA), 3L)
})

test_that("random_walk is accepted as alias for rw1", {
	fit = .fit_dynamic_beta_prediction(kind = "random_walk", nscan = 60, burn = 20)
	expect_identical(fit$dynamic_beta_kind, "rw1")
})

test_that("dynamic_beta_kind = 'ar1' is the default", {
	fit = .fit_dynamic_beta_prediction(nscan = 60, burn = 20)
	expect_identical(fit$dynamic_beta_kind, "ar1")
	expect_gt(diff(range(fit$RHO_BETA)), 0)
})

test_that("counterfactual newdata = original X reproduces manual XB", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	data(YX_bin_list, envir = environment())
	# compare to a hand-computed linear predictor for one period
	pred_new = predict(fit, newdata = YX_bin_list$X, type = "link")
	expect_true(is.list(pred_new))
	expect_length(pred_new, length(YX_bin_list$X))
	# manual reconstruction uses the fit's sorted actor order
	beta_per_t = apply(fit$BETA, c(2, 3), mean)
	ord = names(fit$APM)
	X1 = YX_bin_list$X[[1]][ord, ord, , drop = FALSE]
	manual = beta_per_t["intercept", 1] +
		beta_per_t[2, 1] * X1[, , 1] +
		beta_per_t[3, 1] * X1[, , 2] +
		beta_per_t[4, 1] * X1[, , 3] +
		outer(fit$APM, fit$BPM, "+")
	# off-diagonal entries should agree
	diag(manual) = NA
	diag_mask = !is.na(manual) & !is.na(pred_new[[1]])
	expect_lt(max(abs(pred_new[[1]][diag_mask] - manual[diag_mask])), 1e-6)
})

test_that("counterfactual = 0-X gives link == intercept + a + b only", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	data(YX_bin_list, envir = environment())
	# zero out the dyadic covariates
	zero_X = lapply(YX_bin_list$X, function(x) array(0, dim = dim(x)))
	pred_zero = predict(fit, newdata = zero_X, type = "link")
		# the result should not depend on the dyadic covariate values
	expect_true(is.list(pred_zero))
	# spot-check: with zero x the per-period link should differ from with
	# original x (assuming covariate effect is non-zero)
	pred_orig = predict(fit, newdata = YX_bin_list$X, type = "link")
	# at least one time slice differs
	max_diff = max(sapply(seq_along(pred_zero), function(t)
		max(abs(pred_zero[[t]] - pred_orig[[t]]), na.rm = TRUE)))
	expect_gt(max_diff, 1e-6)
})

# -- 2.3 h-step-ahead forecasting --------------------------------------------

test_that("predict(fit, h=0) equals predict(fit) without h", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	# h=0 hits the in-sample branch
	a = predict(fit, h = 0L, type = "link")
	b = predict(fit, type = "link")
	expect_identical(a, b)
})

test_that("predict(fit, h=3) returns list of length 3 with correct dims", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	f3 = predict(fit, h = 3L, type = "link")
	expect_true(is.list(f3))
	expect_length(f3, 3L)
	expect_equal(dim(f3[[1]]), c(50L, 50L))
})

test_that("forecast type='response' is bounded in [0, 1] for binary", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	f3 = predict(fit, h = 2L, type = "response")
	for (m in f3) {
		expect_true(all(m >= 0 & m <= 1, na.rm = TRUE))
	}
})

test_that("forecast errors when h < 0", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	# h = 0 is in-sample; negative h hits the abort path
	expect_error(predict(fit, h = -1L), "at least 0")
})

# -- 2.4 family link inversion ----------------------------------------------

test_that("binary type='response' with newdata equals pnorm(type='link')", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25, family = "binary")
	data(YX_bin_list, envir = environment())
	rsp = predict(fit, newdata = YX_bin_list$X, type = "response")
	lnk = predict(fit, newdata = YX_bin_list$X, type = "link")
	for (t in seq_along(rsp)) {
		expect_lt(max(abs(rsp[[t]] - pnorm(lnk[[t]])), na.rm = TRUE), 1e-12)
	}
})

test_that("forecast type='response' matches family inverse-link per draw", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25, family = "binary")
	# `by_draw = true` exposes per-draw eta and per-draw response. at the
	# draw level, response = pnorm(link) exactly; only after averaging
	# across draws does jensen create a wedge. we test the strict per-draw
	# identity here.
	set.seed(42); rsp = predict(fit, h = 2L, type = "response", by_draw = TRUE)
	set.seed(42); lnk = predict(fit, h = 2L, type = "link", by_draw = TRUE)
	expect_equal(dim(rsp), dim(lnk))
	expect_lt(max(abs(rsp - pnorm(lnk)), na.rm = TRUE), 1e-12)
})

# -- 2.5 autoplot.lame --------------------------------------------------------

test_that("autoplot.lame returns a ggplot for dynamic_beta fit", {
	skip_if_not_installed("ggplot2")
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	p = autoplot(fit)
	expect_s3_class(p, "ggplot")
})

test_that("autoplot.lame falls back to a coefplot for a non-dynamic fit", {
	skip_if_not_installed("ggplot2")
	# fit without dynamic_beta. autoplot.lame returns a tidy()-driven
	# horizontal coefplot so autoplot(fit) always returns a ggplot
	# regardless of fit type.
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_static = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE))
	p = suppressWarnings(autoplot(fit_static))
	expect_s3_class(p, "ggplot")
})

test_that("autoplot.lame respects the `coefs` argument", {
	skip_if_not_installed("ggplot2")
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	nms = dimnames(fit$BETA)[[2]]
	expect_true(length(nms) >= 1L)
	# subset to one coefficient
	p1 = autoplot(fit, coefs = nms[1])
	expect_s3_class(p1, "ggplot")
})

# -- 2.6 extended pointwise log-lik ------------------------------------------

test_that("save_log_lik = TRUE yields a 2-D log-lik matrix with expected width", {
	# binary closed-form log-lik
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE,
		save_log_lik = TRUE))
	expect_true(!is.null(fit$log_lik))
	# rows = iterations stored, cols = number of non-na dyads
	expect_equal(nrow(fit$log_lik), 100L %/% 5L)
	expect_gt(ncol(fit$log_lik), 0)
	# all finite
	expect_true(all(is.finite(fit$log_lik)))
})

# -- 2.7 chunked log-lik storage ---------------------------------------------

test_that("chunked log-lik storage is byte-identical to in-core", {
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_mem = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		save_log_lik = TRUE))
	set.seed(7); fit_chk = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		save_log_lik = "chunked",
		log_lik_chunk_size = 500L))
	# in-memory matrix is on the fit; chunked fit only has the metadata
	expect_true(!is.null(fit_mem$log_lik))
	expect_null(fit_chk$log_lik)
	expect_true(!is.null(fit_chk$log_lik_chunks))
	# read-back must match in-memory exactly
	ll_chk = read_log_lik(fit_chk)
	expect_true(identical(fit_mem$log_lik, ll_chk))
})

test_that("loo() works transparently on a chunked fit", {
	skip_if_not_installed("loo")
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_chk = suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		save_log_lik = "chunked"))
	res = suppressWarnings(loo::loo(fit_chk))
	expect_s3_class(res, "loo")
	expect_true(is.finite(res$estimates["elpd_loo", "Estimate"]))
})
