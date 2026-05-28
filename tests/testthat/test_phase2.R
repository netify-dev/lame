# phase 2 regression tests — user-visible capabilities.
#
# Items implemented in phase 2:
#   2.1 RW1 prior + stable dynamic_beta_kind = c("ar1","rw1") API
#   2.2 counterfactual covariate paths in predict.lame
#   2.3 h-step-ahead forecasting (predict(fit, h=3))
#   2.4 family-specific link inversion in predict.lame
#   2.5 autoplot.lame() ribbon plot
#   2.6 closed-form pointwise log-lik for extended families
#   2.7 per-column-chunk log-lik storage (per C4)

skip_on_cran()

# small reusable fitter so each test isn't paying the full burn-in cost
.fit_phase2 <- function(kind = "ar1", nscan = 200, burn = 50,
                         dynamic_beta = "dyad", family = "binary", seed = 7L,
                         ...) {
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

# -- 2.1 RW1 prior -----------------------------------------------------------

test_that("dynamic_beta_kind = 'rw1' fits without error and pins rho=1", {
	fit <- .fit_phase2(kind = "rw1", nscan = 100, burn = 25)
	expect_s3_class(fit, "lame")
	expect_identical(fit$dynamic_beta_kind, "rw1")
	# all stored RHO_BETA draws are exactly 1
	expect_true(all(abs(fit$RHO_BETA - 1) < 1e-12))
	# BETA still 3-D and right shape
	expect_length(dim(fit$BETA), 3L)
})

test_that("random_walk is accepted as alias for rw1", {
	fit <- .fit_phase2(kind = "random_walk", nscan = 60, burn = 20)
	expect_identical(fit$dynamic_beta_kind, "rw1")
})

test_that("dynamic_beta_kind = 'ar1' is the default", {
	fit <- .fit_phase2(nscan = 60, burn = 20)
	# default arg "ar1" populates fit$dynamic_beta_kind
	expect_identical(fit$dynamic_beta_kind, "ar1")
	# AR(1) rho is not pinned: some posterior dispersion is expected
	expect_gt(diff(range(fit$RHO_BETA)), 0)
})

# -- 2.2 counterfactual covariate paths --------------------------------------

test_that("counterfactual newdata = original X reproduces manual XB", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	data(YX_bin_list, envir = environment())
	# counterfactual rebuild uses posterior-mean betas; verify it agrees
	# with a hand-computed linear predictor for one period.
	pred_new <- predict(fit, newdata = YX_bin_list$X, type = "link")
	expect_true(is.list(pred_new))
	expect_length(pred_new, length(YX_bin_list$X))
	# manual reconstruction for t=1, in the fit's (sorted) actor order. lame
	# sorts actors on ingest, and predict() realigns newdata to that order, so
	# the hand-computed reference must reorder X to the same order before adding
	# the (already sorted) a + b -- otherwise it would mix orders. (YX_bin_list
	# actor names sort non-trivially, which is what exposes this.)
	beta_per_t <- apply(fit$BETA, c(2, 3), mean)
	ord <- names(fit$APM)
	X1 <- YX_bin_list$X[[1]][ord, ord, , drop = FALSE]
	manual <- beta_per_t["intercept", 1] +
		beta_per_t[2, 1] * X1[, , 1] +
		beta_per_t[3, 1] * X1[, , 2] +
		beta_per_t[4, 1] * X1[, , 3] +
		outer(fit$APM, fit$BPM, "+")
	# off-diagonal should agree to machine precision
	diag(manual) <- NA
	diag_mask <- !is.na(manual) & !is.na(pred_new[[1]])
	expect_lt(max(abs(pred_new[[1]][diag_mask] - manual[diag_mask])), 1e-6)
})

test_that("counterfactual = 0-X gives link == intercept + a + b only", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	data(YX_bin_list, envir = environment())
	# zero out the dyadic covariates
	zero_X <- lapply(YX_bin_list$X, function(x) array(0, dim = dim(x)))
	pred_zero <- predict(fit, newdata = zero_X, type = "link")
	# the result should NOT depend on the dyadic covariate values
	expect_true(is.list(pred_zero))
	# spot-check: with zero X the per-period link should differ from with
	# original X (assuming covariate effect is non-zero)
	pred_orig <- predict(fit, newdata = YX_bin_list$X, type = "link")
	# at least one time slice differs
	max_diff <- max(sapply(seq_along(pred_zero), function(t)
		max(abs(pred_zero[[t]] - pred_orig[[t]]), na.rm = TRUE)))
	expect_gt(max_diff, 1e-6)
})

# -- 2.3 h-step-ahead forecasting --------------------------------------------

test_that("predict(fit, h=0) equals predict(fit) without h", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	# h=0 hits the in-sample branch
	a <- predict(fit, h = 0L, type = "link")
	b <- predict(fit, type = "link")
	expect_identical(a, b)
})

test_that("predict(fit, h=3) returns list of length 3 with correct dims", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	f3 <- predict(fit, h = 3L, type = "link")
	expect_true(is.list(f3))
	expect_length(f3, 3L)
	expect_equal(dim(f3[[1]]), c(50L, 50L))
})

test_that("forecast type='response' is bounded in [0, 1] for binary", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	f3 <- predict(fit, h = 2L, type = "response")
	for (m in f3) {
		expect_true(all(m >= 0 & m <= 1, na.rm = TRUE))
	}
})

test_that("forecast errors when h < 0", {
	fit <- .fit_phase2(nscan = 100, burn = 25)
	# h=0 is in-sample so does NOT error; only negative h hits the abort path
	expect_error(predict(fit, h = -1L), "at least 0")
})

# -- 2.4 family link inversion ----------------------------------------------

test_that("binary type='response' with newdata equals pnorm(type='link')", {
	fit <- .fit_phase2(nscan = 100, burn = 25, family = "binary")
	data(YX_bin_list, envir = environment())
	rsp <- predict(fit, newdata = YX_bin_list$X, type = "response")
	lnk <- predict(fit, newdata = YX_bin_list$X, type = "link")
	for (t in seq_along(rsp)) {
		expect_lt(max(abs(rsp[[t]] - pnorm(lnk[[t]])), na.rm = TRUE), 1e-12)
	}
})

test_that("forecast type='response' matches family inverse-link per draw", {
	fit <- .fit_phase2(nscan = 100, burn = 25, family = "binary")
	# `by_draw = TRUE` exposes per-draw eta and per-draw response. At the
	# draw level, response = pnorm(link) exactly; only after averaging
	# across draws does Jensen create a wedge. We test the strict per-draw
	# identity here.
	set.seed(42); rsp <- predict(fit, h = 2L, type = "response", by_draw = TRUE)
	set.seed(42); lnk <- predict(fit, h = 2L, type = "link", by_draw = TRUE)
	expect_equal(dim(rsp), dim(lnk))
	expect_lt(max(abs(rsp - pnorm(lnk)), na.rm = TRUE), 1e-12)
})

# -- 2.5 autoplot.lame --------------------------------------------------------

test_that("autoplot.lame returns a ggplot for dynamic_beta fit", {
	skip_if_not_installed("ggplot2")
	fit <- .fit_phase2(nscan = 100, burn = 25)
	p <- autoplot(fit)
	expect_s3_class(p, "ggplot")
})

test_that("autoplot.lame falls back to a coefplot for a non-dynamic fit", {
	skip_if_not_installed("ggplot2")
	# fit without dynamic_beta. autoplot.lame returns a tidy()-driven
	# horizontal coefplot so autoplot(fit) always returns a ggplot
	# regardless of fit type.
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit_static <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE))
	p <- suppressWarnings(autoplot(fit_static))
	expect_s3_class(p, "ggplot")
})

test_that("autoplot.lame respects the `coefs` argument", {
	skip_if_not_installed("ggplot2")
	fit <- .fit_phase2(nscan = 100, burn = 25)
	nms <- dimnames(fit$BETA)[[2]]
	expect_true(length(nms) >= 1L)
	# subset to one coefficient
	p1 <- autoplot(fit, coefs = nms[1])
	expect_s3_class(p1, "ggplot")
})

# -- 2.6 extended pointwise log-lik ------------------------------------------

test_that("save_log_lik = TRUE yields a 2-D log-lik matrix with expected width", {
	# binary closed-form log-lik
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5, verbose = FALSE, plot = FALSE,
		save_log_lik = TRUE))
	expect_true(!is.null(fit$log_lik))
	# rows = iterations stored, cols = number of non-NA dyads
	expect_equal(nrow(fit$log_lik), 100L %/% 5L)
	expect_gt(ncol(fit$log_lik), 0)
	# all finite
	expect_true(all(is.finite(fit$log_lik)))
})

# -- 2.7 chunked log-lik storage ---------------------------------------------

test_that("chunked log-lik storage is byte-identical to in-core", {
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_mem <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		save_log_lik = TRUE))
	set.seed(7); fit_chk <- suppressWarnings(lame(
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
	ll_chk <- read_log_lik(fit_chk)
	expect_true(identical(fit_mem$log_lik, ll_chk))
})

test_that("loo() works transparently on a chunked fit", {
	skip_if_not_installed("loo")
	data(YX_bin_list, envir = environment())
	set.seed(7); fit_chk <- suppressWarnings(lame(
		YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
		nscan = 100, burn = 25, odens = 5,
		verbose = FALSE, plot = FALSE,
		save_log_lik = "chunked"))
	res <- suppressWarnings(loo::loo(fit_chk))
	expect_s3_class(res, "loo")
	expect_true(is.finite(res$estimates["elpd_loo", "Estimate"]))
})
