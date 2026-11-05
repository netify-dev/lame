library(lame)
skip_on_cran()

# what a bootstrap does not return: parameters a block resample cannot vary and fixed-effects fit biases get NA standard errors
# and intervals with the reason; block intervals use levels set by the number of blocks; an uncorrected binary fixed-effects bootstrap is refused

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

# the fixed-effects fit (the normal-family default is the empirical-Bayes one)
scope_fit = function(d, ...) {
	suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		Xrow = d$Xr,
		Xcol = d$Xc,
		R = 0,
		family = "normal",
		bias_correct = "none",
		verbose = FALSE,
		...
	))
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

test_that("a block bootstrap returns no standard error or interval for what slice resampling cannot vary", {
	d = scope_panel()
	fit = scope_fit(d)
	expect_identical(fit$bias_correct$method, "none")
	bt = ame_als_bootstrap(fit, R = 30, type = "block", seed = 1, verbose = FALSE)
	expect_identical(bt$type, "block")
	expect_identical(bt$n_blocks, 10L)

	w = bt$withheld
	expect_setequal(w$parameter, c("xr_row", "xc_col", "va", "cab", "vb", "ve"))
	expect_setequal(w$parameter[w$kind == "coefficient"], c("xr_row", "xc_col"))
	expect_setequal(w$parameter[w$kind == "vc"], c("va", "cab", "vb", "ve"))
	expect_true(all(nzchar(w$reason)))

	# the summaries are NA for those and finite for the others
	expect_true(all(is.na(c(
		bt$se[c("xr_row", "xc_col")],
		bt$ci_lo[c("xr_row", "xc_col")],
		bt$ci_hi[c("xr_row", "xc_col")]
	))))
	expect_true(all(is.finite(c(
		bt$se[c("intercept", "xd_dyad")],
		bt$ci_lo[c("intercept", "xd_dyad")],
		bt$ci_hi[c("intercept", "xd_dyad")]
	))))
	expect_true(all(is.na(c(bt$vc_se[c("va", "cab", "vb", "ve")], bt$vc_ci_lo[c("va", "cab", "vb", "ve")]))))
	expect_true(is.finite(bt$vc_se["rho"]) && is.finite(bt$vc_ci_lo["rho"]) && is.finite(bt$vc_ci_hi["rho"]))
	# the replicates themselves are kept
	expect_true(all(is.finite(bt$coefs[, "xr_row"])))
	expect_true(all(is.finite(bt$vc_coefs[, "va"])))

	# confint, vcov and tidy agree
	ci = confint(bt, which = "beta")
	expect_true(all(is.na(ci[c("xr_row", "xc_col"), ])))
	expect_true(all(is.finite(ci[c("intercept", "xd_dyad"), ])))
	cv = confint(bt, which = "vc")
	expect_true(all(is.na(cv[c("va", "cab", "vb", "ve"), ])))
	expect_true(all(is.finite(cv["rho", ])))
	ci_all = confint(bt)
	expect_true(all(is.finite(ci_all[grepl("^[ab]_", rownames(ci_all)), ])))
	V = vcov(bt)
	expect_true(all(is.na(V[c("xr_row", "xc_col"), ])) && all(is.na(V[, c("xr_row", "xc_col")])))
	expect_true(all(is.finite(V[c("intercept", "xd_dyad"), c("intercept", "xd_dyad")])))
	td = tidy(bt)
	expect_true(all(is.na(td$std.error[td$term %in% c("xr_row", "xc_col")])))
	expect_true(all(is.na(td$conf.low[td$term %in% c("xr_row", "xc_col")])))
	expect_true(all(is.finite(td$std.error[td$term %in% c("intercept", "xd_dyad")])))
	expect_true(all(is.na(tidy(bt, conf.level = 0.9)$conf.high[td$term %in% c("xr_row", "xc_col")])))

	# print and summary say so in plain words and point to the parametric scheme
	txt = scope_text(print(bt))
	expect_match(txt, "Block resample of 10 blocks: the interval ends are the 0.86% and 99.1% quantiles", fixed = TRUE)
	expect_match(txt, "No standard error or interval for xr_row, xc_col", fixed = TRUE)
	expect_match(txt, "No standard error or interval for va, cab, vb", fixed = TRUE)
	expect_match(txt, "No standard error or interval for ve", fixed = TRUE)
	expect_match(txt, "type = \"parametric\"` gives intervals for them", fixed = TRUE)
	expect_no_match(txt, "severely|under-cover")
	txt = scope_text(print(summary(bt)))
	expect_match(txt, "No standard error or interval for xr_row, xc_col", fixed = TRUE)
	expect_no_match(txt, "severely|under-cover")
})

test_that("a block bootstrap of a fixed-effects fit with a multiplicative term keeps the factors and the additive effects", {
	d = scope_panel(n = 12L, Tt = 10L, seed = 9L, node = FALSE)
	fit = suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		R = 1,
		family = "normal",
		bias_correct = "none",
		verbose = FALSE
	))
	bt = ame_als_bootstrap(fit, R = 20, type = "block", seed = 2, verbose = FALSE)
	expect_identical(bt$type, "block")
	expect_setequal(bt$withheld$parameter, c("va", "cab", "vb", "ve"))
	expect_true(all(is.finite(bt$se_a)) && all(is.finite(bt$se_b)))
	expect_true(is.finite(bt$se["intercept"]) && is.finite(bt$se["xd_dyad"]))
	expect_true(any(is.finite(bt$se_U)) || !any(bt$sv_separation$resolved))
	ci = confint(bt, which = "a")
	expect_false(anyNA(ci))
})

test_that("a block bootstrap of the empirical-Bayes normal fit has no interval for its intercept and its shrunk effects", {
	d = scope_panel()
	fit = scope_fit_eb(d)
	expect_identical(fit$bias_correct$method, "eb")
	bt = ame_als_bootstrap(fit, R = 30, type = "block", seed = 1, verbose = FALSE)
	expect_setequal(bt$withheld$parameter, c("xr_row", "xc_col", "intercept", "va", "cab", "vb", "ve", "a", "b"))
	expect_true(all(is.na(bt$se[c("intercept", "xr_row", "xc_col")])))
	expect_true(all(is.finite(bt$se["xd_dyad"])) && all(is.finite(c(bt$ci_lo["xd_dyad"], bt$ci_hi["xd_dyad"]))))
	expect_true(all(is.finite(bt$vc_se["rho"])))
	expect_true(all(is.na(bt$se_a)) && all(is.na(bt$se_b)))
	expect_true(all(is.na(confint(bt, which = "a"))))
	expect_true(all(is.finite(confint(bt, which = "beta")["xd_dyad", ])))
	# the parametric scheme of the same fit gives every interval
	bp = ame_als_bootstrap(fit, R = 20, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(all(is.finite(bp$se)) && all(is.finite(bp$se_a)))
})

test_that("the parametric bootstrap of the same fit returns every interval", {
	d = scope_panel()
	fit = scope_fit(d)
	bp = ame_als_bootstrap(fit, R = 30, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(all(is.finite(bp$se)) && all(is.finite(bp$ci_lo)) && all(is.finite(bp$ci_hi)))
	expect_true(all(is.finite(bp$vc_se[c("va", "cab", "vb", "rho", "ve")])))
	expect_false(anyNA(confint(bp, which = "beta")))
	expect_false(anyNA(vcov(bp)))
	txt = scope_text(print(bp))
	expect_no_match(txt, "No standard error")
})

test_that("a bootstrap of a fit carries its withheld parameters to the fit's summaries", {
	d = scope_panel()
	fit = scope_fit(d, bootstrap = 30, bootstrap_type = "block", bootstrap_seed = 2)
	expect_identical(fit$meta$bootstrap_type, "block")
	expect_setequal(fit$bootstrap$withheld$parameter, c("xr_row", "xc_col", "va", "cab", "vb", "ve"))
	td = tidy(fit)
	expect_identical(td$term, names(coef(fit)))
	expect_true(all(is.na(td$std.error[td$term %in% c("xr_row", "xc_col")])))
	expect_true(all(is.finite(td$std.error[td$term %in% c("intercept", "xd_dyad")])))
	V = vcov(fit)
	expect_true(all(is.na(V[c("xr_row", "xc_col"), ])))
	expect_true(all(is.finite(V[c("intercept", "xd_dyad"), c("intercept", "xd_dyad")])))
	ci = confint(fit)
	expect_true(all(is.na(ci[c("xr_row", "xc_col"), ])))
	sm = summary(fit)
	expect_true(all(is.na(sm$coefficients[c("xr_row", "xc_col"), "Std.Error"])))
	expect_true(all(is.finite(sm$coefficients[c("intercept", "xd_dyad"), "Std.Error"])))
	expect_match(scope_text(print(sm)), "No standard error or interval for xr_row, xc_col", fixed = TRUE)
	# the draws of a parameter without an interval are NA, so a summary of them gives none
	dr = unclass(as_draws(fit))
	expect_true(all(is.na(dr[, "xr_row"])))
	expect_true(all(is.finite(dr[, "xd_dyad"])))
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

test_that("an empirical-Bayes block bootstrap has no interval for its intercept and its shrunk effects", {
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
	expect_identical(fit$bias_correct$method, "eb")
	w = lame:::.ae_boot_withheld(fit, "block")
	expect_true("intercept" %in% w$parameter[w$kind == "coefficient"])
	expect_true(all(c("va", "cab", "vb") %in% w$parameter[w$kind == "vc"]))
	expect_setequal(w$parameter[w$kind == "effects"], c("a", "b"))
	# binary ve is fixed at 1 by identification: it is not a replicate summary to withhold
	expect_false("ve" %in% w$parameter)
	expect_null(lame:::.ae_boot_withheld(fit, "parametric"))

	bt = ame_als_bootstrap(fit, R = 20, type = "block", seed = 1, verbose = FALSE)
	expect_identical(bt$type, "block")
	expect_true(is.na(bt$se["intercept"]))
	expect_true(all(is.na(bt$se_a)) && all(is.na(bt$se_b)))
	expect_true(all(is.na(confint(bt, which = "a"))) && all(is.na(confint(bt, which = "b"))))
	expect_true(all(is.finite(confint(bt, which = "vc")["rho", ])))
	expect_match(scope_text(print(summary(bt))), "Sender a: no standard error", fixed = TRUE)
	# the parametric scheme gives all of them
	bp = ame_als_bootstrap(fit, R = 12, type = "parametric", seed = 1, verbose = FALSE)
	expect_null(bp$withheld)
	expect_true(is.finite(bp$se["intercept"]))
	expect_false(anyNA(bp$se_a[is.finite(fit$a)]))
})

test_that("the coefficients of a fit are classed by position", {
	d = scope_panel(n = 8L, Tt = 3L)
	fit = scope_fit(d)
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
	fit = scope_fit(d)
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

# the bootstrap of an uncorrected binary fixed-effects fit

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

test_that("a bootstrap of an uncorrected binary fixed-effects fit is refused before anything is fitted", {
	d = scope_binary()
	fitted = FALSE
	testthat::local_mocked_bindings(
		.ae_run_budget = function(...) {
			fitted <<- TRUE
			stop("fitted")
		},
		.package = "lame"
	)
	err = tryCatch(
		lame_als(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			bias_correct = "none",
			bootstrap = 30,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_bootstrap_bias")
	expect_false(fitted)
	err = tryCatch(
		ame_als(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = "binary",
			bias_correct = "none",
			bootstrap = 30,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_bootstrap_bias")
	expect_false(fitted)
	expect_match(conditionMessage(err), "empirical-Bayes", fixed = TRUE)
	# the dispatchers inherit it
	err = tryCatch(
		lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			bias_correct = "none",
			bootstrap = 30,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_bootstrap_bias")
	expect_false(fitted)
})

test_that("a bootstrap asked of an uncorrected binary fit afterwards is refused too", {
	d = scope_binary()
	fit = suppressWarnings(lame_als(d$Y, Xdyad = d$X, family = "binary", bias_correct = "none", verbose = FALSE))
	expect_identical(fit$bias_correct$method, "none")
	err = tryCatch(ame_als_bootstrap(fit, R = 20, verbose = FALSE), error = identity)
	expect_s3_class(err, "lame_als_bootstrap_bias")
	expect_error(
		ame_als_bootstrap(fit, R = 20, type = "block", verbose = FALSE),
		class = "lame_als_bootstrap_bias"
	)
	# the point fit stays available
	expect_true(is.finite(coef(fit)["xd_dyad"]))
})

test_that("binary fits with a correction or the default estimator are not refused", {
	d = scope_binary()
	fit_eb = suppressWarnings(lame_als(d$Y, Xdyad = d$X, family = "binary", verbose = FALSE))
	expect_identical(fit_eb$bias_correct$method, "eb")
	b = ame_als_bootstrap(fit_eb, R = 10, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
	# the split-sample correction needs 25 rarer-outcome cells per actor at rank zero
	d30 = scope_binary(n = 30L, Tt = 3L)
	fit_jk = suppressWarnings(lame_als(
		d30$Y,
		Xdyad = d30$X,
		family = "binary",
		bias_correct = "jackknife",
		verbose = FALSE
	))
	expect_identical(fit_jk$bias_correct$method, "jackknife")
	b = ame_als_bootstrap(fit_jk, R = 10, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
	fit_bb = suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		bias_correct = "bootstrap",
		verbose = FALSE
	))
	expect_identical(fit_bb$bias_correct$method, "bootstrap")
	b = ame_als_bootstrap(fit_bb, R = 10, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
})

test_that("a fixed-effects Poisson bootstrap has no interval for its intercept", {
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
	fit = suppressWarnings(lame_als(Y, Xdyad = X, family = "poisson", bias_correct = "none", verbose = FALSE))
	bt = ame_als_bootstrap(fit, R = 12, seed = 1, verbose = FALSE)
	expect_identical(bt$withheld$parameter, "intercept")
	expect_match(bt$withheld$reason, "Poisson", fixed = TRUE)
	expect_true(is.na(bt$se["intercept"]) && is.finite(bt$se["xd_dyad"]))
	expect_true(all(is.na(confint(bt, which = "beta")["intercept", ])))
	expect_true(all(is.finite(bt$vc_se[c("va", "vb", "ve")])))
	# the default estimator has the intercept
	fit_eb = suppressWarnings(lame_als(Y, Xdyad = X, family = "poisson", verbose = FALSE))
	expect_identical(fit_eb$bias_correct$method, "eb")
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

test_that("a parametric draw of the effects records the projection the intercept of a fixed-effects fit carries", {
	d = scope_panel()
	fe = scope_fit(d)
	set.seed(1)
	ab = lame:::.ae_boot_draw_ab(fe)
	expect_true(is.finite(attr(ab, "intercept_shift")) && attr(ab, "intercept_shift") != 0)
	# without node covariates the intercept is the mean of the effects, which are drawn centred
	fe0 = suppressWarnings(lame_als(d$Y, Xdyad = d$X, R = 0, family = "normal", bias_correct = "none", verbose = FALSE))
	expect_equal(attr(lame:::.ae_boot_draw_ab(fe0), "intercept_shift"), 0)
	# an empirical-Bayes fit draws its effects uncentred: the intercept is the population one
	eb = scope_fit_eb(d)
	expect_equal(attr(lame:::.ae_boot_draw_ab(eb), "intercept_shift"), 0)
	# replicate intercepts of the fixed-effects fit are taken relative to it, so their spread is that of the noise,
	# far below the spread added by projecting the drawn effects on the node covariates
	bt = ame_als_bootstrap(fe, R = 60, type = "parametric", seed = 3, verbose = FALSE)
	expect_lt(bt$se[["intercept"]], 0.05)
	expect_gt(bt$se[["xr_row"]], 5 * bt$se[["intercept"]])
})

test_that("the variances va and vb of a parametric bootstrap have a bias-corrected interval", {
	f = lame:::.ae_vc_interval
	expect_equal(f(c(0.2, 0.5), 0.3), c(0.09 / 0.5, 0.09 / 0.2))
	expect_equal(f(c(0.2, 0.5), 0.3, eb = TRUE), c(0.1, 0.4))
	expect_equal(f(c(0.2, 0.9), 0.3, eb = TRUE), c(0, 0.4))
	expect_equal(f(c(-0.1, 0.5), 0.3), c(-0.1, 0.5))
	d = scope_panel()
	fit = scope_fit(d)
	bt = ame_als_bootstrap(fit, R = 60, type = "parametric", seed = 4, verbose = FALSE)
	q = stats::quantile(bt$vc_coefs[stats::complete.cases(bt$coefs), "va"], c(0.025, 0.975), names = FALSE)
	expect_equal(unname(c(bt$vc_ci_lo["va"], bt$vc_ci_hi["va"])), f(q, fit$VC[["va"]]))
	expect_equal(unname(confint(bt, which = "vc")["va", ]), f(q, fit$VC[["va"]]))
	# the other components and the basic type keep their percentile ends
	expect_equal(unname(bt$vc_ci_lo["rho"]), unname(stats::quantile(bt$vc_coefs[, "rho"], 0.025)))
	cb = confint(bt, which = "vc", ci_type = "basic")
	expect_equal(unname(cb["va", ]), pmax(unname(2 * fit$VC[["va"]] - rev(q)), 0))
	# an empirical-Bayes fit uses the basic interval
	ebf = scope_fit_eb(d)
	be = ame_als_bootstrap(ebf, R = 30, type = "parametric", seed = 4, verbose = FALSE)
	qe = stats::quantile(be$vc_coefs[stats::complete.cases(be$coefs), "va"], c(0.025, 0.975), names = FALSE)
	expect_equal(unname(confint(be, which = "vc")["va", ]), f(qe, ebf$VC[["va"]], eb = TRUE))
})
