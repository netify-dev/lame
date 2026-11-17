# the fast route has one estimator, the empirical-Bayes fit: the estimator settings of earlier
# releases are not arguments, and a fit saved with one stops at a refit, a bootstrap or a
# covariance with a classed error; a bootstrap reports the analytic error of every quantity it
# gives none for, and a bootstrap that cannot run, or one asked of a dynamic fit, leaves the fit
# with its errors

retired_net = function(fam, seed = 3, n = 25) {
	set.seed(seed)
	x = matrix(rnorm(n * n), n, n)
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	z = -0.3 + 0.5 * x + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
	y = switch(fam, normal = z, binary = 1 * (z > 0), poisson = matrix(rpois(n * n, exp(z / 2)), n, n))
	diag(y) = NA
	dimnames(y) = dimnames(x) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n))
	list(y = y, x = x, w = matrix(rnorm(n), n, 1, dimnames = list(rownames(y), "w")))
}

finite_pos = function(x) length(x) > 0 && all(is.finite(x)) && all(x > 0)

test_that("the removed estimator settings are not arguments", {
	skip_on_cran()
	d = retired_net("binary", n = 12)
	for (nm in c("bias_correct", "non_normal_method")) {
		for (fn in list(ame_als, lame_als, ame, lame, ame_als_refit, ame_als_bootstrap)) expect_false(nm %in% names(formals(fn)))
	}
	expect_error(ame_als(d$y, family = "binary", bias_correct = "jackknife", verbose = FALSE), "unused argument")
	expect_error(ame_als(d$y, family = "binary", non_normal_method = "transform", verbose = FALSE), "unused argument")
	f = ame_als(d$y, family = "binary", verbose = FALSE)
	expect_error(ame_als_refit(f, bias_correct = "jackknife"), "unused argument")
})

test_that("a fit saved with a removed setting stops with a classed error where its estimator is read", {
	skip_on_cran()
	d = retired_net("normal", n = 12)
	f = ame_als(d$y, Xdyad = d$x, verbose = FALSE)
	f$bias_correct = list(method = "jackknife", applied = TRUE)
	expect_error(ame_als_refit(f), class = "lame_als_removed_setting")
	expect_error(ame_als_bootstrap(f, R = 30, verbose = FALSE), class = "lame_als_removed_setting")
	f$eb$vcov = NULL
	expect_error(vcov(f), class = "lame_als_removed_setting")
	expect_true(all(is.finite(coef(f))))
})

test_that("a bootstrap reports the analytic error where it gives none, with the source", {
	skip_on_cran()
	d = retired_net("binary", n = 25)
	f = ame_als(d$y, Xdyad = d$x, Xrow = d$w, R = 2, family = "binary", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 20, seed = 2, verbose = FALSE)
	expect_true(finite_pos(b$se))
	expect_true(finite_pos(b$se_a))
	expect_true(finite_pos(b$se_U[is.finite(as.vector(f$U))]))
	if (!is.null(b$se_source)) {
		expect_true(all(b$se_source$source == "analytic"))
		expect_true(all(nzchar(b$se_source$reason)))
	}
	# the block scheme withholds the node coefficient and va, which keep their analytic errors
	set.seed(5)
	Yl = lapply(1:12, function(s) retired_net("normal", seed = s, n = 12)$y)
	names(Yl) = sprintf("t%02d", 1:12)
	Wl = setNames(rep(list(retired_net("normal", seed = 1, n = 12)$w), 12), names(Yl))
	g = lame_als(Yl, Xrow = Wl, verbose = FALSE)
	expect_no_warning(gb <- lame_als(Yl, Xrow = Wl, verbose = FALSE, bootstrap = 30, bootstrap_type = "block", bootstrap_seed = 1))
	expect_identical(gb$bootstrap$type, "block")
	expect_true(finite_pos(sqrt(diag(vcov(gb)))))
	expect_true(all(c("w_row", "va") %in% gb$bootstrap$se_source$parameter))
	expect_equal(unname(gb$bootstrap$se[["w_row"]]), unname(g$eb$se[["w_row"]]), tolerance = 1e-12)
	expect_equal(unname(gb$bootstrap$vc_se[["va"]]), unname(g$VC_se[["va"]]), tolerance = 1e-12)
	ci = confint(gb)
	expect_true(all(is.finite(ci[c("w_row", "va"), ])))
	out = paste(capture.output(print(summary(gb))), collapse = " ")
	expect_match(out, "Analytic standard errors")
})

test_that("a block bootstrap on a short panel runs the parametric one with a message, not a warning", {
	skip_on_cran()
	Yl = lapply(1:4, function(s) retired_net("normal", seed = s, n = 12)$y)
	names(Yl) = paste0("t", 1:4)
	expect_no_warning(expect_message(
		f <- lame_als(Yl, verbose = FALSE, bootstrap = 30, bootstrap_type = "block", bootstrap_seed = 1),
		class = "lame_als_bootstrap_block_short"
	))
	expect_identical(f$bootstrap$type, "parametric")
	expect_identical(f$bootstrap$requested_type, "block")
	expect_true(finite_pos(sqrt(diag(vcov(f)))))
})

test_that("a bootstrap asked of a dynamic fit, or one that cannot run, leaves the fit with its errors", {
	skip_on_cran()
	Yl = lapply(1:4, function(s) retired_net("binary", seed = s, n = 15)$y)
	names(Yl) = paste0("t", 1:4)
	f0 = lame(Yl, R = 0, family = "binary", method = "als", dynamic_ab = TRUE, verbose = FALSE)
	expect_no_warning(expect_message(
		f1 <- lame(Yl, R = 0, family = "binary", method = "als", dynamic_ab = TRUE, bootstrap = 30, verbose = FALSE),
		class = "lame_als_bootstrap_not_run"
	))
	expect_null(f1$bootstrap)
	expect_equal(coef(f1), coef(f0), tolerance = 1e-12)
	expect_equal(vcov(f1), vcov(f0), tolerance = 1e-12)
	# a count below the minimum stays an error, before anything is fitted
	expect_error(
		lame(Yl, R = 0, family = "binary", method = "als", dynamic_ab = TRUE, bootstrap = 5, verbose = FALSE),
		class = "lame_als_bootstrap_replicates"
	)
	# a static bootstrap that stops with an error after the fit is a message on the fit
	d = retired_net("normal")
	g0 = ame_als(d$y, Xdyad = d$x, verbose = FALSE)
	local_mocked_bindings(ame_als_bootstrap = function(...) cli::cli_abort("too few valid replicates"), .package = "lame")
	expect_message(g1 <- ame_als(d$y, Xdyad = d$x, verbose = FALSE, bootstrap = 30), class = "lame_als_bootstrap_not_run")
	expect_null(g1$bootstrap)
	expect_match(g1$meta$bootstrap_not_run, "too few valid replicates")
	expect_equal(vcov(g1), vcov(g0), tolerance = 1e-12)
})
