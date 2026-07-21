# compatibility of lame's s3 methods with legacy fits from the amen package.
# an amen fit (class "ame") stores only BETA/VC/APM/BPM/U/V/UVPM/EZ/YPM/GOF:
# no $Y, $family, $R or $mode. `$Y` partial-matches `$YPM` on such fits, so
# methods that read `object$Y` must use exact lookup, and metadata reporters
# must not fabricate defaults. lame-native fits carry all fields and must be
# unaffected.

# build each fixture once per file; creation stays inside test bodies so
# skips fire before any fitting happens
.amen_compat_cache = new.env(parent = emptyenv())

# skip helper that does not load amen: requireNamespace() (what
# skip_if_not_installed uses) would register amen's summary.ame / plot.ame
# over lame's for the rest of the test session. find.package only checks.
skip_if_no_amen = function() {
	testthat::skip_if(length(find.package("amen", quiet = TRUE)) == 0L,
	                  "amen not installed")
	testthat::skip_if_not_installed("callr")
}

# the amen fit is produced in a child process: loading amen's namespace here
# would register its own summary.ame / plot.ame for class "ame" and shadow
# lame's methods for the rest of the test session. the child returns the bare
# fit object, which is exactly what a porting user loads from an .rds file.
get_amen_fit = function() {
	if (is.null(.amen_compat_cache$afit)) {
		.amen_compat_cache$afit = callr::r(function() {
			set.seed(6886)
			n = 20
			a = rnorm(n, 0, 0.5)
			b = rnorm(n, 0, 0.5)
			Y = 0.4 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 1), n, n)
			diag(Y) = NA
			amen::ame(Y, R = 1, family = "nrm", nscan = 200, burn = 50,
			          odens = 10, plot = FALSE, print = FALSE, gof = TRUE)
		})
	}
	.amen_compat_cache$afit
}

get_lame_ame_fit = function() {
	if (is.null(.amen_compat_cache$lfit)) {
		set.seed(1)
		n = 15
		Y = matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		.amen_compat_cache$lfit = ame(
			Y, R = 1, family = "normal", nscan = 100, burn = 20, odens = 10,
			verbose = FALSE)
	}
	.amen_compat_cache$lfit
}

test_that("residuals() on an amen fit aborts instead of returning zeros", {
	skip_on_cran()
	skip_if_no_amen()
	afit = get_amen_fit()
	# `$Y` partial-matches `$YPM` on an amen fit, which previously produced a
	# silently all-zero residual matrix
	expect_error(residuals(afit), "does not store the observed")
})

test_that("print() on an amen fit completes instead of hard-erroring", {
	skip_on_cran()
	skip_if_no_amen()
	afit = get_amen_fit()
	# x$mode is NULL on an amen fit; `NULL == "bipartite"` under `&&`
	# previously errored with 'missing value where TRUE/FALSE needed'
	expect_output(print(afit), "Additive and Multiplicative Effects")
})

test_that("glance() on an amen fit recovers R from U and keeps family NA", {
	skip_on_cran()
	skip_if_no_amen()
	afit = get_amen_fit()
	# a tidier stays quiet; the NA family is the visible signal
	expect_silent(g <- glance(afit))
	expect_identical(as.integer(g$R), 1L)
	expect_true(is.na(g$family))
	# nobs counts the finite cells of YPM (the NA diagonal is excluded)
	expect_identical(as.integer(g$nobs), sum(is.finite(afit$YPM)))
})

test_that("sampler_describe() on an amen fit reports unknown family, not normal/0", {
	skip_on_cran()
	skip_if_no_amen()
	afit = get_amen_fit()
	out = paste(cli::ansi_strip(
		capture.output(sampler_describe(afit), type = "message")),
		collapse = "\n")
	expect_match(out, "unknown", fixed = TRUE)
	expect_match(out, "Latent rank R: 1", fixed = TRUE)
	expect_no_match(out, "\"normal\"", fixed = TRUE)
})

test_that("lame-native ame fit is unaffected: residuals/print/glance as before", {
	skip_on_cran()
	lfit = get_lame_ame_fit()
	r = residuals(lfit)
	expect_true(is.matrix(r))
	expect_equal(r, lfit$Y - lfit$YPM)
	expect_false(all(r == 0, na.rm = TRUE))
	expect_output(print(lfit), "Additive and Multiplicative Effects")
	g = glance(lfit)
	expect_identical(g$family, "normal")
	expect_identical(as.integer(g$R), 1L)
})

test_that("lame-native lame fit is unaffected: residuals per period as before", {
	skip_on_cran()
	set.seed(2)
	n = 12
	Tt = 3
	Ylist = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n)
		diag(m) = NA
		m
	})
	names(Ylist) = paste0("t", seq_len(Tt))
	fit = lame(Ylist, R = 0, family = "normal", nscan = 60, burn = 20,
	           odens = 10, verbose = FALSE)
	r = residuals(fit)
	expect_true(is.list(r))
	expect_length(r, Tt)
	expect_false(all(vapply(r, function(m) all(m == 0, na.rm = TRUE),
	                        logical(1))))
})

test_that("family alias message fires once per session, not per call", {
	skip_on_cran()
	# reset the show-once flag so this test is order-independent
	.lame_state = lame:::.lame_state
	.lame_state$family_alias_msg_shown = NULL
	set.seed(3)
	n = 12
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	n_alias_msgs = 0L
	for (i in 1:3) {
		withCallingHandlers(
			fit <- ame(Y, R = 0, family = "nrm", nscan = 40, burn = 10,
			           odens = 10, verbose = FALSE),
			message = function(m) {
				if (grepl("amen-style", conditionMessage(m))) {
					n_alias_msgs <<- n_alias_msgs + 1L
				}
				invokeRestart("muffleMessage")
			})
	}
	expect_identical(n_alias_msgs, 1L)
	# the alias still resolves on the silent calls
	expect_identical(fit$family, "normal")
	# lame() shares the flag: no repeat after ame() has shown the note
	Ylist = list(t1 = Y, t2 = Y)
	n_lame_msgs = 0L
	withCallingHandlers(
		lame(Ylist, R = 0, family = "nrm", nscan = 40, burn = 10, odens = 10,
		     verbose = FALSE),
		message = function(m) {
			if (grepl("amen-style", conditionMessage(m))) {
				n_lame_msgs <<- n_lame_msgs + 1L
			}
			invokeRestart("muffleMessage")
		})
	expect_identical(n_lame_msgs, 0L)
})

test_that("vcov.ame_als returns the bootstrap covariance when one is attached", {
	skip_on_cran()
	set.seed(4)
	n = 24
	a = rnorm(n, 0, 0.4)
	b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	Y = 0.3 + 0.6 * Xd + outer(a, b, "+") +
		matrix(rnorm(n * n, 0, 0.7), n, n)
	diag(Y) = NA
	fb = ame_als(Y, Xdyad = Xd, R = 0, family = "normal",
	             bootstrap = 40, bootstrap_seed = 7, verbose = FALSE)
	V = vcov(fb)
	bc = fb$bootstrap$coefs
	# vcov(), summary() and confint() now report the same uncertainty
	expect_equal(V, stats::cov(bc), ignore_attr = TRUE)
	expect_identical(dim(V), c(ncol(bc), ncol(bc)))
	se_v = sqrt(diag(V))
	expect_equal(unname(se_v[names(fb$bootstrap$se)]),
	             unname(fb$bootstrap$se), tolerance = 1e-10)
	sm_se = summary(fb)$coefficients[, "Std.Error"]
	shared = intersect(names(se_v), names(sm_se))
	expect_gt(length(shared), 0L)
	expect_equal(unname(se_v[shared]), unname(sm_se[shared]),
	             tolerance = 1e-8)
	# confint() takes the bootstrap path; every vcov coefficient is covered
	# and its percentile interval is consistent with the bootstrap se
	ci = suppressMessages(confint(fb))
	expect_true(all(names(se_v) %in% rownames(ci)))
	implied_se = (ci[names(se_v), 2] - ci[names(se_v), 1]) /
		(2 * stats::qnorm(0.975))
	rat = implied_se / se_v
	expect_true(all(rat > 0.6 & rat < 1.6))
	# an explicit cluster= is announced, not silently discarded
	expect_message(V2 <- vcov(fb, cluster = "none"),
	               "applies to the sandwich only")
	expect_equal(V2, V)
})

test_that("vcov.ame_als without a bootstrap keeps the sandwich and its warning", {
	skip_on_cran()
	set.seed(5)
	n = 24
	Xd = matrix(rnorm(n * n), n, n)
	Y = 0.3 + 0.6 * Xd + matrix(rnorm(n * n, 0, 0.7), n, n)
	diag(Y) = NA
	f0 = ame_als(Y, Xdyad = Xd, R = 0, family = "normal", verbose = FALSE)
	V = vcov(f0)
	expect_identical(dim(V), c(2L, 2L))
	expect_identical(rownames(V), c("intercept", "dyad1_dyad"))
	# the anti-conservative warning still fires for a discrete-family fit
	Yb = 1 * (Y > median(Y, na.rm = TRUE))
	diag(Yb) = NA
	fp = suppressWarnings(ame_als(Yb, Xdyad = Xd, R = 0, family = "binary",
	                              verbose = FALSE))
	expect_warning(vcov(fp), "conditional sandwich")
})

test_that("reconstruct_EZ and predict(type = 'link') work on a legacy amen fit", {
	skip_on_cran()
	skip_if_no_amen()
	afit = get_amen_fit()
	ez = reconstruct_EZ(afit)
	expect_true(is.matrix(ez))
	expect_true(all(is.finite(ez[!is.na(afit$YPM)])))
	# amen stores EZ directly; the accessor must return it, not crash on
	# the missing $mode metadata
	expect_equal(ez, afit$EZ)
	pr = predict(afit, type = "link")
	expect_true(is.matrix(pr))
	expect_equal(dim(pr), dim(afit$YPM))
})
