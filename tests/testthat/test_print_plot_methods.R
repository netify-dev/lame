# print, summary, plot, autoplot, glance, update and sampler_describe methods on small fits that run on CRAN
skip_on_cran()

# everything a method writes, standard output and messages (cli) together
.pm_text = function(expr) {
	msg = capture.output(out <- capture.output(value <- expr), type = "message")
	structure(paste(c(out, msg), collapse = "\n"), value = value)
}

.pm_panel = function(n = 8, TT = 4, seed = 2) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	X = lapply(seq_len(TT), function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, "x1")))
	Y = lapply(seq_len(TT), function(t) {
		m = 0.5 * X[[t]][, , 1] + matrix(rnorm(n * n), n); diag(m) = NA; dimnames(m) = list(nm, nm); m })
	list(Y = Y, X = X)
}
.pm_quiet = function(expr) suppressWarnings(suppressMessages(expr))

test_that("print methods of the post-fit estimators state what they hold", {
	d = .pm_panel()
	adb = als_dynamic_beta(d$Y, d$X, lambda = 1)
	txt = .pm_text(print(adb))
	expect_match(txt, "lambda: 1"); expect_match(txt, "NOT a posterior summary")
	expect_identical(attr(txt, "value"), adb)
	fl = .pm_quiet(lame(
		d$Y,
		Xdyad = d$X,
		family = "normal",
		R = 0,
		burn = 10,
		nscan = 40,
		odens = 2,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1)
	))
	pas = per_actor_slopes(fl, kind = "row", covariate_idx = 1L, lambda = 1)
	txt = .pm_text(print(pas))
	expect_match(txt, "actors x periods: 8 x 4"); expect_match(txt, "2 more rows")
	g = gof_temporal(fl, stat = "mean", n_rep = 8, seed = 4)
	txt = .pm_text(print(g))
	expect_match(txt, "Temporal-trend"); expect_match(txt, "Replicates: 8")
	lf = .pm_quiet(lfo(fl, periods = 4, refit = FALSE))
	txt = .pm_text(print(lf))
	expect_match(txt, "Total elpd_lfo"); expect_match(txt, "Refit per leave-out: FALSE")
	expect_match(txt, "n_obs")
	y4 = d$Y[[4]]
	pit = .pm_quiet(forecast_pit(fl, y4, h = 1))
	txt = .pm_text(print(pit))
	expect_match(txt, "Forecast PIT"); expect_match(txt, "Held-out dyads")
	expect_s3_class(plot(pit), "ggplot")
})

test_that("print methods of the multi-panel and robust fits report their structure", {
	d = .pm_panel()
	fm = .pm_quiet(lame_multi(
		Y_list = list(d$Y, d$Y),
		Xdyad_list = list(d$X, d$X),
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 2,
		verbose = FALSE,
		plot = FALSE
	))
	txt = .pm_text(print(fm))
	expect_match(txt, "K = 2"); expect_match(txt, "Shared beta posterior mean"); expect_match(txt, "x1_dyad")
	fr = .pm_quiet(lame_robust(d$Y, d$X, R = 0, burn = 0, nscan = 4, odens = 1, verbose = FALSE))
	txt = .pm_text(print(fr))
	expect_match(txt, "Robust longitudinal network model"); expect_match(txt, "4 retained draws")
	sr = summary(fr)
	txt = .pm_text(print(sr))
	expect_match(txt, "posterior summary"); expect_match(txt, "x1")
	expect_identical(attr(txt, "value"), sr)
})

test_that("simulation print and summary report the densities of the simulated networks", {
	set.seed(1); n = 10; nm = sprintf("a%02d", 1:n)
	Y = matrix(rbinom(n * n, 1, .3), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	fa = .pm_quiet(
		ame(
			Y,
			family = "binary",
			R = 0,
			burn = 10,
			nscan = 40,
			odens = 2,
			seed = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		)
	)
	s = simulate(fa, nsim = 3, seed = 2)
	expect_s3_class(s, "ame.sim")
	txt = .pm_text(print(s))
	expect_match(txt, "Networks simulated: 3"); expect_match(txt, "binary")
	sm = .pm_quiet(summary(s))
	dens = vapply(s$Y, function(y) mean(y[row(y) != col(y)], na.rm = TRUE), 0)
	expect_equal(sm$densities, dens)
	expect_identical(sm$n_sims, 3L)
	# plot.ame draws the trace plots and returns them invisibly
	expect_invisible(p <- plot(fa))
	expect_s3_class(p, "ggplot")
	expect_s3_class(autoplot(fa), "ggplot")
})

test_that("autoplot and glance of fast fits, and the lame aliases of the ame methods", {
	set.seed(4); n = 8; nm = sprintf("a%02d", 1:n)
	Y = matrix(rnorm(n * n), n, dimnames = list(nm, nm)); diag(Y) = NA
	fals = .pm_quiet(ame_als(Y, family = "normal", R = 0, verbose = FALSE))
	expect_s3_class(autoplot(fals), "ggplot")
	Yl = lapply(1:3, function(t) { m = matrix(rnorm(n * n), n, dimnames = list(nm, nm)); diag(m) = NA; m })
	flals = .pm_quiet(lame(Yl, family = "normal", R = 0, method = "als", verbose = FALSE))
	expect_s3_class(flals, "lame_als")
	expect_s3_class(autoplot(flals), "ggplot")
	expect_identical(lame:::glance.lame_als(flals), lame:::glance.ame_als(flals))
	fl = .pm_quiet(
		lame(
			Yl,
			family = "normal",
			R = 0,
			burn = 5,
			nscan = 20,
			odens = 2,
			seed = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		)
	)
	gl = .pm_quiet(lame:::glance.lame(fl))
	expect_identical(gl, .pm_quiet(lame:::glance.ame(fl)))
	expect_identical(gl$n_periods, 3L)
})

test_that("update() of an ame fit rebuilds the call and refits with the changes", {
	set.seed(5); n = 8; nm = sprintf("a%02d", 1:n)
	Y = matrix(rbinom(n * n, 1, .3), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	fa = .pm_quiet(
		ame(
			Y,
			family = "binary",
			R = 0,
			burn = 10,
			nscan = 40,
			odens = 2,
			seed = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		)
	)
	cl = update(fa, nscan = 20, evaluate = FALSE)
	expect_true(is.call(cl)); expect_equal(cl$nscan, 20); expect_equal(cl$burn, 10)
	f2 = .pm_quiet(update(fa, nscan = 20))
	expect_s3_class(f2, "ame"); expect_equal(nrow(f2$BETA), 10L)
})

test_that("sampler_describe() of an MCMC lame fit names the sampler and the model", {
	set.seed(6); n = 8; nm = sprintf("a%02d", 1:n)
	Yl = lapply(1:3, function(t) { m = matrix(rbinom(n * n, 1, .3), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	fl = .pm_quiet(lame(
		Yl,
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		burn = 5,
		nscan = 20,
		odens = 2,
		seed = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	expect_silent(v <- sampler_describe(fl, verbose = FALSE))
	expect_identical(v, fl)
	txt = .pm_text(sampler_describe(fl))
	expect_match(txt, "Bayesian MCMC"); expect_match(txt, "Periods T: 3"); expect_match(txt, "dynamic_ab")
	expect_match(txt, "Stored draws: 10")
})

test_that("als_dynamic_beta names the coefficient rows after the intercept and the covariate slices", {
	set.seed(12)
	n = 6; TT = 3
	X = lapply(
		seq_len(TT),
		function(t) array(rnorm(n * n * 2), c(n, n, 2), dimnames = list(NULL, NULL, c("trade", "dist")))
	)
	Y = lapply(
		seq_len(TT),
		function(t) { y = X[[t]][, , 1] - X[[t]][, , 2] + matrix(rnorm(n * n), n); diag(y) = NA; y }
	)
	fit = als_dynamic_beta(Y, X, lambda = 1)
	expect_identical(dimnames(fit$beta), list(c("intercept", "trade", "dist"), c("t1", "t2", "t3")))
	expect_identical(rownames(als_dynamic_beta(Y, X, intercept = FALSE)$beta), c("trade", "dist"))
	# unnamed slices are numbered; a single covariate may be a matrix
	Xu = lapply(X, function(x) { dimnames(x) = NULL; x })
	expect_identical(rownames(als_dynamic_beta(Y, Xu)$beta), c("intercept", "x1", "x2"))
	X1 = lapply(X, function(x) x[, , 1])
	expect_identical(rownames(als_dynamic_beta(Y, X1)$beta), c("intercept", "x1"))
	expect_identical(coef(fit), fit$beta)
	expect_error(als_dynamic_beta(Y, c(X[1:2], list(X[[3]][, , 1, drop = FALSE]))), "same number of covariates")
})
