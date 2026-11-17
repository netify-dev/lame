library(lame)
skip_on_cran()

# what a bootstrap does not return: parameters a block resample cannot vary carry their analytic errors, labelled;
# block intervals use levels set by the number of blocks

scope_panel = function(n = 12L, Tt = 10L, seed = 5L, node = TRUE) {
	set.seed(seed)
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	xr = rnorm(n); xc = rnorm(n)
	nm = sprintf("a%02d", seq_len(n))
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, n)
		y = 0.3 + 0.6 * x + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.6), n, n)
		if (node) y = y + 0.5 * xr + matrix(-0.4 * xc, n, n, byrow = TRUE)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(nm, nm, "xd"))
	}
	Xr = lapply(seq_len(Tt), function(t) matrix(xr, n, 1, dimnames = list(nm, "xr")))
	Xc = lapply(seq_len(Tt), function(t) matrix(xc, n, 1, dimnames = list(nm, "xc")))
	list(Y = Y, X = X, Xr = if (node) Xr, Xc = if (node) Xc)
}

scope_fit_eb = function(d, ...) {
	suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		Xrow = d$Xr,
		Xcol = d$Xc,
		R = 0,
		family = "normal",
		verbose = FALSE,
		...
	))
}

# the text a print method writes, on one line
scope_text = function(expr) {
	out = utils::capture.output(expr)
	gsub("\\s+", " ", paste(out, collapse = " "))
}

# the block bootstrap

test_that("a block bootstrap of the empirical-Bayes normal fit keeps the analytic errors of its intercept and its shrunk effects", {
	d = scope_panel()
	fit = scope_fit_eb(d)
	expect_null(fit$bias_correct)
	bt = ame_als_bootstrap(fit, R = 30, type = "block", seed = 1, verbose = FALSE)
	expect_setequal(bt$withheld$parameter, c("xr_row", "xc_col", "intercept", "va", "cab", "vb", "ve", "a", "b"))
	# the scheme gives them no replicate error, so they carry the analytic one, labelled
	w = c("intercept", "xr_row", "xc_col")
	expect_equal(unname(bt$se[w]), unname(fit$eb$se[w]), tolerance = 1e-12)
	expect_setequal(bt$se_source$parameter, c(w, "va", "cab", "vb", "ve", "a", "b"))
	expect_true(all(bt$se_source$source == "analytic"))
	expect_equal(unname(bt$vc_se[c("va", "cab", "vb", "ve")]), unname(fit$VC_se[c("va", "cab", "vb", "ve")]), tolerance = 1e-12)
	expect_true(all(is.finite(bt$se["xd_dyad"])) && all(is.finite(c(bt$ci_lo["xd_dyad"], bt$ci_hi["xd_dyad"]))))
	expect_true(all(is.finite(bt$vc_se["rho"])))
	expect_equal(unname(bt$se_a), unname(fit$a_sd), tolerance = 1e-12)
	expect_equal(unname(bt$se_b), unname(fit$b_sd), tolerance = 1e-12)
	expect_true(all(is.finite(confint(bt, which = "a"))))
	expect_true(all(is.finite(confint(bt))))
	expect_true(all(is.finite(confint(bt, which = "beta")["xd_dyad", ])))
	expect_match(scope_text(print(summary(bt))), "Analytic standard errors and intervals for", fixed = TRUE)
	# the parametric scheme of the same fit gives every interval
	bp = ame_als_bootstrap(fit, R = 20, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(all(is.finite(bp$se)) && all(is.finite(bp$se_a)))
})

test_that("the parametric bootstrap of the empirical-Bayes fit returns every interval", {
	d = scope_panel()
	fit = scope_fit_eb(d)
	bp = ame_als_bootstrap(fit, R = 30, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(all(is.finite(bp$se)) && all(is.finite(bp$ci_lo)) && all(is.finite(bp$ci_hi)))
	expect_true(all(is.finite(bp$vc_se[c("va", "cab", "vb", "rho", "ve")])))
	expect_false(anyNA(confint(bp, which = "beta")))
	expect_false(anyNA(vcov(bp)))
	txt = scope_text(print(bp))
	expect_no_match(txt, "No standard error")
})

test_that("the scheme that ran is recorded when a block request is replaced by the parametric one", {
	d = scope_panel(n = 8L, Tt = 6L, seed = 6L)
	fit = suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		R = 0,
		family = "normal",
		verbose = FALSE,
		bootstrap = 30,
		bootstrap_type = "block",
		bootstrap_seed = 1
	))
	expect_identical(fit$bootstrap$type, "parametric")
	expect_identical(fit$meta$bootstrap_type, "parametric")
	expect_null(fit$bootstrap$withheld)
})

test_that("an empirical-Bayes block bootstrap without node covariates has an interval for its intercept and none for its shrunk effects", {
	set.seed(21)
	n = 12L; Tt = 10L
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.7)
	Y = lapply(seq_len(Tt), function(t) {
		y = 1 * (-0.2 + outer(a, a, "+") * 0.5 + matrix(rnorm(n * n), n, n) > 0)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		y
	})
	fit = suppressWarnings(lame_als(Y, R = 0, family = "binary", verbose = FALSE))
	expect_null(fit$bias_correct)
	w = lame:::.ae_boot_withheld(fit, "block")
	# the intercept is that of the centred effects
	expect_false("intercept" %in% w$parameter)
	expect_true(all(c("va", "cab", "vb") %in% w$parameter[w$kind == "vc"]))
	expect_setequal(w$parameter[w$kind == "effects"], c("a", "b"))
	# binary ve is fixed at 1 by identification: it is not a replicate summary to withhold
	expect_false("ve" %in% w$parameter)
	expect_null(lame:::.ae_boot_withheld(fit, "parametric"))

	bt = ame_als_bootstrap(fit, R = 20, type = "block", seed = 1, verbose = FALSE)
	expect_identical(bt$type, "block")
	expect_true(is.finite(bt$se["intercept"]))
	# the shrunk effects keep their analytic posterior sds
	expect_equal(unname(bt$se_a), unname(fit$a_sd), tolerance = 1e-12)
	expect_equal(unname(bt$se_b), unname(fit$b_sd), tolerance = 1e-12)
	expect_true(all(is.finite(confint(bt, which = "a"))) && all(is.finite(confint(bt, which = "b"))))
	expect_true(all(is.finite(confint(bt, which = "vc")["rho", ])))
	expect_no_match(scope_text(print(summary(bt))), "Sender a: no standard error", fixed = TRUE)
	# the parametric scheme gives all of them
	bp = ame_als_bootstrap(fit, R = 12, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(is.finite(bp$se["intercept"]))
	expect_false(anyNA(bp$se_a[is.finite(fit$a)]))
})

test_that("the coefficients of a fit are classed by position", {
	d = scope_panel(n = 8L, Tt = 3L)
	fit = scope_fit_eb(d)
	expect_identical(lame:::.ae_coef_kinds(fit), c("intercept", "dyad", "node", "node"))
	fit0 = suppressWarnings(lame_als(d$Y, R = 0, family = "normal", verbose = FALSE))
	expect_identical(lame:::.ae_coef_kinds(fit0), "intercept")
})

# the levels of a block interval

test_that("block intervals are read at the expanded percentile levels", {
	pr = lame:::.ae_level_probs
	expect_equal(pr(0.95, "parametric", NA), c(0.025, 0.975))
	expect_equal(pr(0.95, "parametric", 10L), c(0.025, 0.975))
	# ten blocks: Phi(sqrt(10 / 9) * qt(0.025, 9))
	expect_equal(pr(0.95, "block", 10L)[1L], stats::pnorm(sqrt(10 / 9) * stats::qt(0.025, 9)))
	expect_lt(pr(0.95, "block", 10L)[1L], 0.01)
	expect_equal(sum(pr(0.95, "block", 10L)), 1)
	# the levels move toward the nominal ones as the number of blocks grows
	expect_gt(pr(0.95, "block", 20L)[1L], pr(0.95, "block", 10L)[1L])
	expect_lt(pr(0.95, "block", 20L)[1L], 0.025)
	expect_equal(pr(0.95, "block", 1e6)[1L], 0.025, tolerance = 1e-4)
	expect_equal(pr(0.9, "block", 10L)[1L], stats::pnorm(sqrt(10 / 9) * stats::qt(0.05, 9)))

	d = scope_panel()
	fit = scope_fit_eb(d)
	bt = ame_als_bootstrap(fit, R = 60, type = "block", seed = 3, verbose = FALSE)
	p10 = lame:::.ae_boot_probs(bt, 0.95)
	expect_equal(p10, pr(0.95, "block", 10L))
	x = bt$coefs[stats::complete.cases(bt$coefs), "xd_dyad"]
	ci = confint(bt, parm = "xd_dyad")
	expect_equal(unname(ci[1, ]), unname(stats::quantile(x, p10)))
	expect_identical(colnames(ci), c("2.5%", "97.5%"))
	# the interval stored on the object is at the same levels, and is wider than the plain percentile one
	expect_equal(unname(c(bt$ci_lo["xd_dyad"], bt$ci_hi["xd_dyad"])), unname(ci[1, ]))
	plain = stats::quantile(x, c(0.025, 0.975))
	expect_lt(ci[1, 1], plain[[1]])
	expect_gt(ci[1, 2], plain[[2]])
	ci9 = confint(bt, parm = "xd_dyad", level = 0.9)
	expect_equal(unname(ci9[1, ]), unname(stats::quantile(x, pr(0.9, "block", 10L))))
	expect_identical(colnames(ci9), c("5%", "95%"))
	# tidy() reads the same levels at another confidence level
	expect_equal(tidy(bt, conf.level = 0.9)$conf.low[2], unname(ci9[1, 1]))
	# the parametric scheme reads the nominal ones
	bp = ame_als_bootstrap(fit, R = 60, type = "parametric", seed = 3, verbose = FALSE)
	xp = bp$coefs[stats::complete.cases(bp$coefs), "xd_dyad"]
	expect_equal(unname(confint(bp, parm = "xd_dyad")[1, ]), unname(stats::quantile(xp, c(0.025, 0.975))))
	expect_true(is.na(bp$n_blocks))
})

# the bootstrap of the default binary and count fits

scope_binary = function(n = 20L, Tt = 3L, seed = 3L) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(Tt), function(t) {
		m = 1 * (matrix(rnorm(n * n), n, n) > 0.3)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	X = lapply(seq_len(Tt), function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, "xd")))
	list(Y = Y, X = X)
}

test_that("a binary fit of the default estimator is bootstrapped", {
	d = scope_binary()
	fit_eb = suppressWarnings(lame_als(d$Y, Xdyad = d$X, family = "binary", verbose = FALSE))
	expect_null(fit_eb$bias_correct)
	b = ame_als_bootstrap(fit_eb, R = 10, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
})

test_that("a Poisson bootstrap of the default estimator has an interval for its intercept", {
	set.seed(8)
	n = 16L; Tt = 3L
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(Tt), function(t) {
		m = matrix(rpois(n * n, 2), n, n)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	X = lapply(seq_len(Tt), function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, "xd")))
	fit_eb = suppressWarnings(lame_als(Y, Xdyad = X, family = "poisson", verbose = FALSE))
	expect_null(fit_eb$bias_correct)
	expect_null(ame_als_bootstrap(fit_eb, R = 10, seed = 1, verbose = FALSE)$withheld)
})

test_that("ame(method = 'als') forwards prior$s20 to the ALS estimator instead of warning that it is ignored", {
	set.seed(31)
	n = 14L
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	expect_no_warning(f1 <- ame(Y, family = "normal", method = "als", prior = list(s20 = 3), verbose = FALSE))
	f2 = ame_als(Y, family = "normal", prior = list(s20 = 3), verbose = FALSE)
	expect_equal(f1$VC, f2$VC)
	f0 = ame(Y, family = "normal", method = "als", verbose = FALSE)
	expect_false(isTRUE(all.equal(f1$VC[["ve"]], f0$VC[["ve"]])))
	# an MCMC-only argument is reported
	expect_warning(ame(Y, family = "normal", method = "als", nscan = 10, verbose = FALSE), "nscan")
})

# the parametric intercept and variance intervals

test_that("a parametric draw of the empirical-Bayes effects carries no intercept projection", {
	d = scope_panel()
	# an empirical-Bayes fit draws its effects uncentred: the intercept is the population one
	eb = scope_fit_eb(d)
	set.seed(1)
	expect_null(attr(lame:::.ae_boot_draw_ab(eb), "intercept_shift"))
})

test_that("the variances va and vb of a parametric bootstrap have a bias-corrected interval", {
	f = lame:::.ae_vc_interval
	expect_equal(f(c(0.2, 0.5), 0.3), c(0.1, 0.4))
	expect_equal(f(c(0.2, 0.9), 0.3), c(0, 0.4))
	d = scope_panel()
	# an empirical-Bayes fit uses the basic interval
	ebf = scope_fit_eb(d)
	be = ame_als_bootstrap(ebf, R = 30, type = "parametric", seed = 4, verbose = FALSE)
	qe = stats::quantile(be$vc_coefs[stats::complete.cases(be$coefs), "va"], c(0.025, 0.975), names = FALSE)
	expect_equal(unname(confint(be, which = "vc")["va", ]), f(qe, ebf$VC[["va"]]))
	# the other components keep their percentile ends
	expect_equal(unname(be$vc_ci_lo["rho"]), unname(stats::quantile(be$vc_coefs[, "rho"], 0.025)))
})
