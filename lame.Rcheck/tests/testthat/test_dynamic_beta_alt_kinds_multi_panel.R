# dynamic beta, log-likelihood, and multi-panel checks

skip_on_cran()

.fit_dynamic_beta_alt = function(kind = "ar1", ll_method = "observed_exact",
                                  suppress_messages = TRUE, ...) {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	expr = quote(suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad",
		dynamic_beta_kind = kind,
		log_lik_method = ll_method,
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE, ...)))
	if (suppress_messages) suppressMessages(eval(expr)) else eval(expr)
}

test_that("dynamic_beta_kind = 'rw2' fits and is recorded as rw2", {
	fit = .fit_dynamic_beta_alt(kind = "rw2")
	expect_s3_class(fit, "lame")
	expect_identical(fit$dynamic_beta_kind, "rw2")
	expect_true(all(is.finite(fit$BETA)))
	expect_length(dim(fit$BETA), 3L)
})

test_that("rw2 path has modest second differences", {
	set.seed(11)
	fit = .fit_dynamic_beta_alt(kind = "rw2")
	post_mean = apply(fit$BETA, c(2, 3), mean)
	d2 = t(apply(post_mean, 1, function(p) diff(diff(p))))
	expect_lt(max(abs(d2)), 5)
})

test_that("dynamic_beta_kind = 'matern32' fits and is recorded as matern32", {
	fit = .fit_dynamic_beta_alt(kind = "matern32")
	expect_s3_class(fit, "lame")
	expect_identical(fit$dynamic_beta_kind, "matern32")
	expect_true(all(is.finite(fit$BETA)))
	expect_length(dim(fit$BETA), 3L)
})

test_that("invalid dynamic_beta_kind aborts", {
	data(YX_bin_list, envir = environment())
	expect_error(suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_beta = "dyad", dynamic_beta_kind = "magic",
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"dynamic_beta_kind")
})

test_that("log_lik_method default is observed_exact", {
	fit = .fit_dynamic_beta_alt()
	expect_identical(fit$log_lik_method, "observed_exact")
})

test_that("log_lik_method = 'observed_ghk' emits the bias note", {
	expect_message(.fit_dynamic_beta_alt(ll_method = "observed_ghk", suppress_messages = FALSE),
	               "biased downward")
})

test_that("observed_ghk produces a finite log-lik matrix on binary", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		log_lik_method = "observed_ghk",
		save_log_lik = TRUE,
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_true(!is.null(fit[["log_lik"]]))
	expect_true(all(is.finite(fit[["log_lik"]])))
})

test_that("log_lik_method = 'observed_ghk' warns about downward bias", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	expect_message(suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		log_lik_method = "observed_ghk",
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)),
		"biased downward")
})

test_that("log_lik_method = 'augmented' is accepted silently", {
	expect_silent(.fit_dynamic_beta_alt(ll_method = "augmented"))
})

test_that("dynamic_G wired for bipartite + RA/RB > 0", {
	set.seed(7)
	nA = 8; nB = 6; T_per = 3
	Y_list = lapply(seq_len(T_per), function(t) matrix(rnorm(nA*nB), nA, nB))
	X_list = lapply(seq_len(T_per), function(t) array(rnorm(nA*nB*1), c(nA, nB, 1)))
	fit = suppressWarnings(suppressMessages(lame(
		Y_list, X_list, family = "normal", R = 2,
		mode = "bipartite",
		dynamic_G = TRUE,
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_true(isTRUE(fit$dynamic_G))
	expect_equal(dim(fit$G_cube), c(2L, 2L, T_per))
	expect_true(all(is.finite(fit$G_cube)))
})

test_that("dynamic_G on unipartite errors (G = I for unipartite)", {
	data(YX_bin_list, envir = environment())
	expect_error(suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		dynamic_G = TRUE,
		nscan = 30, burn = 10, odens = 5,
		verbose = FALSE, plot = FALSE))),
		"bipartite")
})

test_that("per_actor_slopes returns the expected shape and label", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	pas = per_actor_slopes(fit, kind = "row", covariate_idx = 1L, lambda = 1)
	expect_s3_class(pas, "per_actor_slopes")
	expect_equal(dim(pas$slopes), c(50L, 4L))
	expect_true(all(is.finite(pas$slopes)))
})

test_that("per_actor_slopes errors on bad lambda and bad kind", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 30, burn = 10, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_error(per_actor_slopes(fit, lambda = -1), "non-negative")
	expect_error(per_actor_slopes(fit, kind = "bogus"))
})

test_that("lame_multi with K=1 falls through to standard lame()", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X),
		family = "binary", R = 0,
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_s3_class(fm, "lame_multi")
	expect_equal(fm$K, 1L)
})

test_that("lame_multi with K=2 identical panels has zero deviations", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y, YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X, YX_bin_list$X),
		family = "binary", R = 0,
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_equal(fm$K, 2L)
	expect_length(fm$beta_deviations, 2L)
	# identical panels should have near-zero deviations
	for (k in seq_len(fm$K)) {
		expect_lt(max(abs(fm$beta_deviations[[k]])), 1e-10)
	}
})

test_that("lame_multi gives proper 1/K variance shrinkage with K identical panels", {
	skip_on_cran()
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y, YX_bin_list$Y, YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X, YX_bin_list$X, YX_bin_list$X),
		family = "binary", R = 0,
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_equal(fm$K, 3L)
	expect_lt(abs(fm$var_shrinkage - 1/3), 0.05)
	expect_true(!is.null(fm$BETA_joint))
	expect_true(!is.null(fm$beta_shared_se))
})
