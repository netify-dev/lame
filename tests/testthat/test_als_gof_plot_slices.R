library(lame)
skip_on_cran()

# gof_plot() of a longitudinal ame_als fit uses every time slice, as for an MCMC lame fit

gof_panel = function(n = 14, Tt = 4, seed = 41L) {
	set.seed(seed)
	lapply(seq_len(Tt), function(t) {
		z = matrix(rnorm(n * n, mean = 0.3 * t), n)
		diag(z) = NA_real_
		z
	})
}

test_that("a longitudinal fit's gof plot carries the statistics of every period", {
	y = gof_panel()
	fit = lame_als(y, R = 0, family = "normal", verbose = FALSE)
	expect_no_warning(p <- gof_plot(fit, nsim = 15, seed = 1))
	expect_s3_class(p, "ggplot")
	d = p$data
	expect_setequal(unique(d$time), 1:4)
	# every statistic has a row in every period, with the observed value of that period
	n_stat = length(unique(d$statistic))
	expect_equal(nrow(d), 4L * n_stat)
	expect_true(all(c("observed", "median", "lower", "upper") %in% names(d)))
	expect_true(all(d$lower <= d$median & d$median <= d$upper))
	# the observed statistic of period t is computed on slice t, not on the first slice
	obs_t = vapply(1:4, function(t) gof_stats(fit$Y[, , t], mode = "unipartite", family = "normal")[["sd.rowmean"]], 0)
	got = d$observed[d$statistic == lame:::gof_label("sd.rowmean")]
	expect_equal(got, obs_t)
	expect_gt(length(unique(round(got, 6))), 1L)
})

test_that("the dispatcher and the method give the same plot data", {
	fit = lame_als(gof_panel(), R = 0, family = "normal", verbose = FALSE)
	p1 = gof_plot(fit, nsim = 10, seed = 3)
	p2 = gof_plot.ame_als(fit, nsim = 10, seed = 3)
	expect_equal(p1$data, p2$data)
})

test_that("a binary and a bipartite longitudinal fit are plotted over every period too", {
	set.seed(5)
	yb = lapply(1:3, function(t) {
		m = matrix(rbinom(144, 1, 0.3), 12)
		diag(m) = NA
		m
	})
	fb = lame_als(yb, R = 0, family = "binary", bias_correct = "none", verbose = FALSE)
	p = gof_plot(fb, nsim = 10, seed = 1)
	expect_setequal(unique(p$data$time), 1:3)
	yr = lapply(1:3, function(t) matrix(rnorm(15 * 9), 15, 9))
	fr = lame_als(yr, R = 0, family = "normal", mode = "bipartite", verbose = FALSE)
	pr = gof_plot(fr, nsim = 10, seed = 1)
	expect_setequal(unique(pr$data$time), 1:3)
})

test_that("a cross-sectional fit keeps its histogram of replicates", {
	set.seed(6)
	Y = matrix(rnorm(144), 12)
	diag(Y) = NA
	fit = ame_als(Y, R = 0, family = "normal", verbose = FALSE)
	p = gof_plot(fit, nsim = 20, seed = 1)
	expect_s3_class(p, "ggplot")
	expect_true(all(c("statistic", "observed", "value") %in% names(p$data)))
	expect_equal(nrow(p$data), 20L * length(unique(p$data$statistic)))
})
