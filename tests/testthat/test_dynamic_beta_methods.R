# dynamic_beta methods: summary()'s rho_beta stationarity table, detect_change_point(), prior summaries, freeze_call / update(),
# compact printing, pooling, gof_temporal(), rhat_dynamic_beta(), time_index / period_exposure, tidy, checkpointing, lfo(), forecasting, link inversion, autoplot and log-lik storage
skip_on_cran()
# diagnostics, prior summaries, freeze_call / update, compact printing
test_that("stationarity warning fires when rho_beta posterior is forced near 1", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	# a concentrated near-unit posterior
	fit$RHO_BETA = matrix(
		rep(rnorm(20, 0.985, 0.01), each = ncol(fit$RHO_BETA)),
		nrow = 20,
		ncol = ncol(fit$RHO_BETA),
		dimnames = dimnames(fit$RHO_BETA)
	)
	s = summary(fit)
	expect_true(any(s$rho_beta_stationarity$warn))
	out = capture.output(print(s))
	expect_true(any(grepl("concentrated near 1", out)))
	expect_true(any(grepl("dynamic_beta_kind", out)))
})
test_that("stationarity warning does NOT fire on a diffuse rho_beta posterior", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	# a diffuse posterior
	set.seed(1)
	fit$RHO_BETA = matrix(
		runif(20 * ncol(fit$RHO_BETA), 0.3, 0.9),
		nrow = 20,
		ncol = ncol(fit$RHO_BETA),
		dimnames = dimnames(fit$RHO_BETA)
	)
	s = summary(fit)
	expect_false(any(s$rho_beta_stationarity$warn))
	out = capture.output(print(s))
	expect_false(any(grepl("concentrated near 1", out)))
})
test_that("90%-NA Y fits without error and produces finite BETA", {
	set.seed(2026)
	n = 10; Tn = 3
	Y = replicate(Tn, {
		z = matrix(rnorm(n * n), n, n); diag(z) = NA
		# a very sparse observed network
		mask = matrix(stats::runif(n * n) > 0.10, n, n)
		z[mask] = NA; diag(z) = NA
		z
	}, simplify = FALSE)
	X = replicate(Tn, array(rnorm(n * n), c(n, n, 1)), simplify = FALSE)
	fit = lame(
		Y,
		Xdyad = X,
		family = "normal",
		R = 0,
		nscan = 30,
		burn = 5,
		odens = 5,
		dynamic_beta = "dyad",
		verbose = FALSE
	)
	expect_true(all(is.finite(fit$BETA)))
	expect_equal(length(dim(fit$BETA)), 3L)
})
test_that("detect_change_point returns the expected data frame shape", {
	dat = simulate_test_network(n = 6, n_time = 4, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 30,
		burn = 5,
		odens = 5,
		dynamic_beta = "dyad",
		verbose = FALSE
	)
	res = detect_change_point(fit, n_prior_sims = 200, seed = 1)
	expect_named(res, c("coef", "bf", "m_post_mean", "m_prior_q95", "t_hat", "warn"))
	expect_type(res$bf, "double")
	expect_type(res$warn, "logical")
})
test_that("detect_change_point errors when dynamic_beta is not active", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 5,
		verbose = FALSE
	)
	expect_error(detect_change_point(fit), "dynamic_beta")
})
test_that("dynamic_beta_prior_summary is reproducible with seed", {
	a = dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	b = dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	expect_equal(a$summary, b$summary)
	expect_equal(a$paths, b$paths)
})
test_that("larger sigma_scale (= more diffuse) gives larger max-first-diff", {
	small = dynamic_beta_prior_summary(sigma_scale = 0.05, seed = 1, ndraws = 500)
	big = dynamic_beta_prior_summary(sigma_scale = 5.0, seed = 1, ndraws = 500)
	expect_gt(big$summary$q50[1], small$summary$q50[1])
	expect_gt(big$summary$q50[2], small$summary$q50[2])
})
test_that("rw1 kind has rho fixed at 1 and path is unit-root random walk", {
	s = dynamic_beta_prior_summary(kind = "rw1", seed = 1, ndraws = 500)
	expect_true(all(s$rho == 1))
	expect_equal(s$kind, "rw1")
})
test_that("invalid rho_mean outside bounds errors cleanly", {
	expect_error(
		dynamic_beta_prior_summary(rho_mean = 1.5, ndraws = 100),
		"must be strictly between"
	)
	expect_error(
		dynamic_beta_prior_summary(rho_mean = -0.1, ndraws = 100),
		"must be strictly between"
	)
})
test_that("n_periods = 1 errors with informative message", {
	expect_error(dynamic_beta_prior_summary(n_periods = 1), "at least 2")
})
test_that("freeze_call = TRUE: update() refits with snapshot even if local Y mutates", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	Y_local = dat$Y
	set.seed(7)
	fit = lame(
		Y_local,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		freeze_call = TRUE,
		verbose = FALSE
	)
	# change Y_local after the fit
	Y_local[[1]][1, 2] = Y_local[[1]][1, 2] + 1000
	# update() uses the snapshot, not the changed Y_local
	set.seed(7)
	fit_upd = update(fit, nscan = 10)
	# a refit on the original data with the same arguments reproduces the fit
	set.seed(7)
	fit_ref = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		freeze_call = TRUE,
		verbose = FALSE
	)
	expect_equal(fit_upd$BETA, fit_ref$BETA)
})
test_that("freeze_call = FALSE uses parent-frame update lookup", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE
	)
	expect_null(fit$data_snapshot)
	expect_null(fit$freeze_call)
	# update() also works through the parent.frame() lookup
	fit_upd = update(fit, nscan = 5)
	expect_s3_class(fit_upd, "lame")
})
# compact print
test_that("single-dynamic fit does NOT use the compact table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	out = paste(
		c(
			capture.output(print(fit)),
			capture.output(print(fit), type = "message")
		),
		collapse = "\n"
	)
	# a single dynamic_beta prints the regular per-component line, not the compact joint table
	expect_false(grepl("Dynamic components$|Dynamic components\\s", out))
})
test_that("joint dynamic_ab + dynamic_beta fit uses the compact table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		dynamic_ab = TRUE,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	out = paste(
		c(
			capture.output(print(fit)),
			capture.output(print(fit), type = "message")
		),
		collapse = "\n"
	)
	expect_true(grepl("Dynamic components", out))
})
test_that("compact = FALSE disables the joint table", {
	dat = simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		dynamic_ab = TRUE,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	out = paste(
		c(
			capture.output(print(fit, compact = FALSE)),
			capture.output(print(fit, compact = FALSE), type = "message")
		),
		collapse = "\n"
	)
	expect_false(grepl("Dynamic components", out))
})
# tidy() and forecasting: pooling, gof_temporal, rhat, time_index, tidy, checkpoint, lfo
.fit_dynamic_beta_tidy = function(seed = 7L, ...) {
	data(YX_bin_list, envir = environment())
	set.seed(seed)
	suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		dynamic_beta = "dyad",
		nscan = 100,
		burn = 25,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = seed,
		...
	))
}
test_that("dynamic_beta_pool != 'none' fits with hierarchical pooling", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		dynamic_beta = "dyad",
		dynamic_beta_pool = "both",
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_identical(fit$dynamic_beta_pool, "both")
	expect_true(all(is.finite(fit$BETA)))
})
test_that("gof_temporal returns a structured result on a dynamic fit", {
	fit = .fit_dynamic_beta_tidy(posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1))
	res = suppressWarnings(gof_temporal(fit, stat = "density", n_rep = 30, seed = 42))
	expect_s3_class(res, "gof_temporal")
	expect_true(is.numeric(res$slope_obs))
	expect_length(res$slope_rep, 30L)
	expect_true(is.finite(res$p_pp) && res$p_pp >= 0 && res$p_pp <= 1)
})
test_that("gof_temporal supports reciprocity and transitivity", {
	fit = .fit_dynamic_beta_tidy(posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1))
	r2 = suppressWarnings(gof_temporal(fit, stat = "reciprocity", n_rep = 20, seed = 42))
	r3 = suppressWarnings(gof_temporal(fit, stat = "transitivity", n_rep = 20, seed = 42))
	expect_s3_class(r2, "gof_temporal")
	expect_s3_class(r3, "gof_temporal")
})
test_that("gof_temporal errors on T < 3", {
	fit = .fit_dynamic_beta_tidy()
	bad = fit
	bad$Y_obs = fit$Y[, , 1:2]
	expect_error(
		suppressWarnings(gof_temporal(bad, stat = "density")),
		"at least 3"
	)
})
test_that("rhat_dynamic_beta requires at least two fits", {
	fit = .fit_dynamic_beta_tidy()
	expect_error(rhat_dynamic_beta(list(fit)), "at least two")
})
test_that("rhat_dynamic_beta requires dynamic_beta", {
	data(YX_bin_list, envir = environment())
	fit_s = suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 100,
		burn = 25,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	))
	expect_error(
		rhat_dynamic_beta(list(fit_s, fit_s)),
		"3-D"
	)
})
# mixed frequency
test_that("time_index strictly-increasing validation", {
	data(YX_bin_list, envir = environment())
	expect_error(
		suppressWarnings(suppressMessages(lame(
			YX_bin_list$Y,
			YX_bin_list$X,
			family = "binary",
			R = 0,
			time_index = c(1, 1, 2, 3), # not strictly increasing
			nscan = 30,
			burn = 10,
			odens = 5,
			verbose = FALSE,
			plot = FALSE
		))),
		"strictly increasing"
	)
})
test_that("time_index equal-gap is byte-identical to no time_index", {
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_eq = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		time_index = c(1, 2, 3, 4),
		nscan = 30,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	set.seed(7); fit_no = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 30,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_true(identical(fit_eq$BETA, fit_no$BETA))
})
# tidy.lame and prediction_draws_long
test_that("tidy.lame returns data frame with expected columns for dynamic fit", {
	fit = .fit_dynamic_beta_tidy()
	td = tidy(fit)
	expect_s3_class(td, "data.frame")
	expect_true(all(c(
		"term", "period", "estimate", "std.error",
		"statistic", "p.value", "conf.low", "conf.high"
	) %in% names(td)))
	expect_equal(nrow(td), dim(fit$BETA)[2L] * dim(fit$BETA)[3L])
})
test_that("tidy.lame on static fit returns one row per coef (no period col)", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_s = suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 100,
		burn = 25,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	))
	td_s = tidy(fit_s)
	expect_equal(nrow(td_s), ncol(fit_s$BETA))
	expect_false("period" %in% names(td_s))
})
test_that("tidy.lame conf.int = FALSE drops CI columns", {
	fit = .fit_dynamic_beta_tidy()
	td = tidy(fit, conf.int = FALSE)
	expect_false("conf.low" %in% names(td))
	expect_false("conf.high" %in% names(td))
})
test_that("prediction_draws_long returns the right shape", {
	fit = .fit_dynamic_beta_tidy()
	pd = prediction_draws_long(fit, n_draws = 3L, type = "link", seed = 42)
	expect_s3_class(pd, "data.frame")
	# columns follow the tidybayes / marginaleffects convention (.draw / .iteration / .value) plus actor names from the dimnames of fit$Y
	expect_true(all(c(
		".chain", ".iteration", ".draw", "period",
		"period_label", "i", "j", "actor_i", "actor_j", ".value"
	)
		%in% names(pd)))
	expect_equal(nrow(pd), 3L * dim(fit$BETA)[3L] * 50L * 50L)
})
test_that("lame_resume continues from a checkpoint", {
	data(YX_bin_list, envir = environment())
	ck = tempfile(fileext = ".rds")
	set.seed(7)
	fit1 = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 200L,
		burn = 50L,
		odens = 5L,
		checkpoint_path = ck,
		checkpoint_every = 10L,
		max_seconds = 0.3,
		verbose = FALSE,
		plot = FALSE
	)))
	skip_if_not(
		isTRUE(fit1$terminated_early),
		"first fit did not terminate early; resume not exercised"
	)
	# the checkpoint carries the sampler state and the seed, so the
	# continuation restarts where the chain stopped
	ck_obj = readRDS(ck)
	expect_true(all(c("Z", "beta", "a", "b", "U", "V", "rho", "s2", "Sab") %in%
			names(ck_obj$start_vals)))
	expect_equal(
		dim(ck_obj$start_vals$Z),
		c(
			nrow(YX_bin_list$Y[[1]]),
			ncol(YX_bin_list$Y[[1]]),
			length(YX_bin_list$Y)
		)
	)
	expect_equal(ck_obj$seed, 6886)
	fit2 = suppressWarnings(suppressMessages(
		lame_resume(ck, nscan_more = 30L)
	))
	expect_s3_class(fit2, "lame")
	expect_true(nrow(fit2$BETA) >= 1L)
	# the continuation ran from the saved state, without burn-in and with a
	# seed advanced by the completed iterations
	expect_identical(fit2$call$start_vals, as.name(".lame_resume_state"))
	expect_identical(fit2$call$burn, as.integer(ck_obj$burn_remaining))
	expect_identical(fit2$call$seed, as.integer(6886 + ck_obj$iter))
})
test_that("lame_resume errors on missing checkpoint", {
	expect_error(
		lame_resume("/tmp/nonexistent_lame_checkpoint.rds"),
		"not found"
	)
})
# prediction and log-likelihood: transition kinds, counterfactuals, forecasts, autoplot, log-lik
# small fit shared by the tests below
.fit_dynamic_beta_prediction = function(
	kind = "ar1",
	nscan = 200,
	burn = 50,
	dynamic_beta = "dyad",
	family = "binary",
	seed = 7L,
	...
) {
	data(YX_bin_list, envir = environment())
	set.seed(seed)
	suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = family,
		R = 0,
		dynamic_beta = dynamic_beta,
		dynamic_beta_kind = kind,
		nscan = nscan,
		burn = burn,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
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
	# a hand-computed linear predictor for one period
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
	# off-diagonal entries agree
	diag(manual) = NA
	diag_mask = !is.na(manual) & !is.na(pred_new[[1]])
	expect_lt(max(abs(pred_new[[1]][diag_mask] - manual[diag_mask])), 1e-6)
})
# h-step-ahead forecasting
test_that("predict(fit, h=0) equals predict(fit) without h", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	# h = 0 takes the in-sample branch
	a = predict(fit, h = 0L, type = "link")
	b = predict(fit, type = "link")
	expect_identical(a, b)
})
test_that("forecast errors when h < 0", {
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	# h = 0 is in-sample; a negative h aborts
	expect_error(predict(fit, h = -1L), "at least 0")
})
# family link inversion
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
	# by_draw = TRUE exposes per-draw eta and response, where response = pnorm(link) exactly; averaging over draws
	# opens a Jensen gap, so the per-draw identity is tested
	set.seed(42); rsp = predict(fit, h = 2L, type = "response", by_draw = TRUE)
	set.seed(42); lnk = predict(fit, h = 2L, type = "link", by_draw = TRUE)
	expect_equal(dim(rsp), dim(lnk))
	expect_lt(max(abs(rsp - pnorm(lnk)), na.rm = TRUE), 1e-12)
})
# autoplot.lame
test_that("autoplot.lame returns a ggplot for dynamic_beta fit", {
	skip_if_not_installed("ggplot2")
	fit = .fit_dynamic_beta_prediction(nscan = 100, burn = 25)
	p = autoplot(fit)
	expect_s3_class(p, "ggplot")
})
test_that("autoplot.lame falls back to a coefplot for a non-dynamic fit", {
	skip_if_not_installed("ggplot2")
	# a fit without dynamic_beta: autoplot.lame draws a tidy()-driven coefficient plot, so autoplot(fit) is always a ggplot
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_static = suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	))
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

test_that("dynamic beta progress and live plots preserve every sampled draw", {
	seen = list()
	local_mocked_bindings(plot.lame = function(x, ...) {
		seen[[length(seen) + 1L]] <<- x
		invisible(NULL)
	}, .package = "lame")
	for (mode in c("unipartite", "bipartite")) for (rank in 0:1) for (dyn_ab in c(FALSE, TRUE)) {
		dat = simulate_test_network(
			n = 8,
			nA = 8,
			nB = 6,
			n_time = 3,
			family = "binary",
			mode = mode,
			seed = 222
		)
		args = list(
			Y = dat$Y,
			Xdyad = dat$Xdyad,
			family = "binary",
			mode = mode,
			R = rank,
			dynamic_beta = "dyad",
			dynamic_ab = dyn_ab,
			nscan = 20,
			burn = 5,
			odens = 1,
			seed = 771,
			gof = FALSE,
			plot = FALSE,
			verbose = FALSE
		)
		quiet = suppressMessages(do.call(lame, args))
		args$verbose = TRUE; args$plot = TRUE
		start = length(seen)
		output = capture.output({noisy = suppressMessages(do.call(lame, args))})
		label = paste(mode, rank, dyn_ab)
		expect_identical(noisy$BETA, quiet$BETA, info = label)
		expect_identical(noisy$VC, quiet$VC, info = label)
		expect_identical(noisy$YPM, quiet$YPM, info = label)
		expect_equal(length(seen) - start, 2L, info = label)
		for (j in seq_len(length(seen) - start)) {
			shown = seen[[start + j]]
			expect_true(isTRUE(shown$dynamic_beta))
			expect_identical(dim(shown$BETA), c(10L * j, dim(quiet$BETA)[2:3]))
			expect_identical(shown$BETA, quiet$BETA[seq_len(10L * j), , , drop = FALSE])
		}
		if (length(seen) > start) {
			shown = seen[[length(seen)]]
			for (nm in c("APM", "BPM", "EZ", "YPM")) expect_identical(shown[[nm]], quiet[[nm]], info = paste(label, nm))
		}
	}
})
