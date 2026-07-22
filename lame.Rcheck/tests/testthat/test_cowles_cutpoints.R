# cowles (1996) z-marginalised mh for explicit ordinal cutpoints.
# tests cover: (a) helper math, (b) byte-identical-default contract,
# (c) end-to-end recovery from a probit-cutpoint dgp, (d) symmetric and
# asymmetric paths, (e) lame() panel path, (f) bipartite fallback.

test_that("log_phi_diff is numerically stable in both tails", {
	# upper tail: both hi, lo > ~7 -- naive log(pnorm(hi) - pnorm(lo)) suffers
	# catastrophic cancellation. reference via survival function:
	# log(s(lo) - s(hi)) is stable.
	hi = 8
	lo = 7.5
	val = lame:::log_phi_diff(hi, lo)
	expect_true(is.finite(val))
	# stable reference: s(lo) - s(hi), both representable as 1e-14 ish
	ref = log(pnorm(lo, lower.tail = FALSE) - pnorm(hi, lower.tail = FALSE))
	expect_lt(abs(val - ref), 1e-8)
	# lower tail: both hi, lo < ~-7 -- mirror image
	hi2 = -7.5
	lo2 = -8
	val2 = lame:::log_phi_diff(hi2, lo2)
	expect_true(is.finite(val2))
	ref2 = log(pnorm(hi2) - pnorm(lo2))
	expect_lt(abs(val2 - ref2), 1e-8)
	# middle: should agree with naive log(pnorm(hi) - pnorm(lo)) when finite
	val3 = lame:::log_phi_diff(0.5, -0.5)
	ref3 = log(pnorm(0.5) - pnorm(-0.5))
	expect_lt(abs(val3 - ref3), 1e-10)
	# infinite bounds: collapses to single cdf
	expect_equal(lame:::log_phi_diff(Inf, 0.5),
	             pnorm(0.5, lower.tail = FALSE, log.p = TRUE),
	             tolerance = 1e-12)
	expect_equal(lame:::log_phi_diff(0.5, -Inf),
	             pnorm(0.5, log.p = TRUE),
	             tolerance = 1e-12)
})

test_that("alpha <-> delta reparameterisation roundtrip is exact", {
	alpha = c(0, 0.4, 1.1, 2.3)
	delta = lame:::.alpha_to_delta(alpha)
	alpha2 = lame:::.delta_to_alpha(delta)
	expect_equal(alpha, alpha2, tolerance = 1e-12)
	# binary edge case
	expect_equal(lame:::.alpha_to_delta(c(0)), numeric(0))
	expect_equal(lame:::.delta_to_alpha(numeric(0)), 0)
})

test_that(".init_alpha_from_data anchors alpha_1 = 0 and is monotone", {
	set.seed(42)
	Y = matrix(sample(1:4, 100, replace = TRUE), 10, 10)
	diag(Y) = NA
	alpha = lame:::.init_alpha_from_data(Y)
	expect_equal(length(alpha), 3L)
	expect_equal(alpha[1L], 0)
	expect_true(all(diff(alpha) > 0))
})

test_that(".init_alpha_from_data aborts on < 2 distinct categories", {
	Y = matrix(1, 6, 6); diag(Y) = NA
	expect_error(lame:::.init_alpha_from_data(Y),
	             "at least 2 distinct")
})

test_that("sample_alpha_cowles accepts/rejects sanely and preserves anchor", {
	set.seed(7)
	n = 20
	K = 4
	true_alpha = c(0, 0.6, 1.4)
	EZ = matrix(rnorm(n*n) * 0.3, n, n); diag(EZ) = 0
	Z = EZ + matrix(rnorm(n*n), n, n)
	Y_int = matrix(findInterval(Z, c(-Inf, true_alpha, Inf)), n, n)
	diag(Y_int) = NA
	alpha = lame:::.init_alpha_from_data(Y_int)
	# 200 sweeps from a reasonable rw proposal sd
	out_alpha = matrix(NA, 200, length(alpha))
	tau = 0.3
	for (s in seq_len(200)) {
		step = lame:::sample_alpha_cowles(alpha, Y_int, EZ, tau, symmetric = FALSE)
		alpha = step$alpha
		out_alpha[s, ] = alpha
	}
	# anchor is preserved across all sweeps
	expect_true(all(out_alpha[, 1L] == 0))
	# monotonicity preserved across all sweeps
	expect_true(all(apply(out_alpha[, -1L, drop = FALSE], 1, function(r) all(diff(r) > 0)) |
	                ncol(out_alpha) == 2L))
	# posterior mean is close-ish to the truth (loose tolerance, short chain)
	pm = colMeans(out_alpha[100:200, , drop = FALSE])
	expect_lt(abs(pm[2L] - true_alpha[2L]), 0.5)
	expect_lt(abs(pm[3L] - true_alpha[3L]), 0.5)
})

test_that("rZ_ord_explicit_fc respects truncation bounds", {
	set.seed(11)
	n = 12
	alpha = c(0, 0.8, 1.5)
	af = c(-Inf, alpha, Inf)
	EZ = matrix(rnorm(n*n) * 0.5, n, n); diag(EZ) = 0
	Y_int = matrix(sample(1:4, n*n, replace = TRUE), n, n); diag(Y_int) = NA
	Z = matrix(0, n, n)
	Z = rZ_ord_explicit_fc(Z, EZ, Y_int, alpha)
	# every cell with y_int = k must have z in [af[k], af[k+1])
	obs = which(!is.na(Y_int))
	yi = Y_int[obs]
	zi = Z[obs]
	lo = af[yi]
	hi = af[yi + 1L]
	expect_true(all(zi > lo - 1e-8))
	expect_true(all(zi < hi + 1e-8))
})

test_that("rZ_ord_sym_explicit_fc returns a symmetric Z", {
	set.seed(13)
	n = 10
	alpha = c(0, 1.0)
	EZ = matrix(rnorm(n*n) * 0.3, n, n)
	EZ = (EZ + t(EZ))/2
	Y_int = matrix(sample(1:3, n*n, replace = TRUE), n, n)
	Y_int = (Y_int + t(Y_int))/2  # rough symmetry
	Y_int[lower.tri(Y_int)] = t(Y_int)[lower.tri(Y_int)]
	diag(Y_int) = NA
	Z = matrix(0, n, n)
	Z = rZ_ord_sym_explicit_fc(Z, EZ, Y_int, alpha)
	# off-diagonal symmetry to machine precision
	off = upper.tri(Z)
	expect_equal(Z[off], t(Z)[off], tolerance = 1e-10)
})

test_that("ordinal_cutpoints = 'data_induced' is byte-identical to default", {
	skip_on_cran()
	set.seed(2026)
	n = 12
	Y = matrix(sample(1:3, n*n, replace = TRUE), n, n); diag(Y) = NA
	f1 = ame(Y, family = "ordinal", R = 0, nscan = 50, burn = 25, odens = 5,
	         rvar = FALSE, cvar = FALSE, intercept = FALSE, seed = 42,
	         verbose = FALSE)
	f2 = ame(Y, family = "ordinal", R = 0, nscan = 50, burn = 25, odens = 5,
	         rvar = FALSE, cvar = FALSE, intercept = FALSE, seed = 42,
	         verbose = FALSE, ordinal_cutpoints = "data_induced")
	expect_identical(f1$BETA, f2$BETA)
	expect_identical(f1$VC, f2$VC)
	expect_identical(f1$YPM, f2$YPM)
	# alpha is null for data_induced
	expect_null(f1$ALPHA)
	expect_null(f2$ALPHA)
})

test_that("ame() explicit-cutpoint path stores ALPHA + ordinal_levels", {
	skip_on_cran()
	set.seed(23)
	n = 18
	EZ = matrix(rnorm(n*n) * 0.4, n, n); diag(EZ) = 0
	true_alpha = c(0, 0.7, 1.5)
	Z = EZ + matrix(rnorm(n*n), n, n)
	Y = matrix(findInterval(Z, c(-Inf, true_alpha, Inf)), n, n); diag(Y) = NA
	fit = ame(Y, family = "ordinal", R = 0, nscan = 80, burn = 40, odens = 5,
	          rvar = FALSE, cvar = FALSE, intercept = FALSE,
	          verbose = FALSE, ordinal_cutpoints = "explicit")
	expect_equal(fit$ordinal_cutpoints, "explicit")
	expect_equal(length(fit$ordinal_levels), 4L)
	expect_false(is.null(fit$ALPHA))
	expect_equal(ncol(fit$ALPHA), 2L)
	expect_equal(colnames(fit$ALPHA), c("alpha_2", "alpha_3"))
	# rough recovery: posterior mean within 0.4 of truth on a short chain
	expect_lt(abs(fit$alpha_post_mean[2L] - true_alpha[2L]), 0.4)
	expect_lt(abs(fit$alpha_post_mean[3L] - true_alpha[3L]), 0.4)
})

test_that("lame() panel explicit-cutpoint path runs and stores ALPHA", {
	skip_on_cran()
	set.seed(31)
	n = 15; Tn = 3
	Y_list = lapply(seq_len(Tn), function(t) {
		EZ = matrix(rnorm(n*n) * 0.3, n, n); diag(EZ) = 0
		Z = EZ + matrix(rnorm(n*n), n, n)
		y = matrix(findInterval(Z, c(-Inf, 0, 0.8, 1.5, Inf)), n, n)
		diag(y) = NA
		y
	})
	fit = lame(Y_list, family = "ordinal", R = 0, nscan = 60, burn = 30, odens = 5,
	           rvar = FALSE, cvar = FALSE, verbose = FALSE,
	           ordinal_cutpoints = "explicit")
	expect_equal(fit$ordinal_cutpoints, "explicit")
	expect_false(is.null(fit$ALPHA))
	expect_equal(ncol(fit$ALPHA), 2L)
})

test_that("ordinal_cutpoints = 'explicit' warns and falls back for bipartite", {
	skip_on_cran()
	set.seed(41)
	nA = 8; nB = 10; Tn = 2
	Y_list = lapply(seq_len(Tn), function(t) {
		mat = matrix(sample(1:3, nA*nB, replace = TRUE), nA, nB)
		rownames(mat) = paste0("a", seq_len(nA))
		colnames(mat) = paste0("b", seq_len(nB))
		mat
	})
	expect_warning({
		fit = lame(Y_list, family = "ordinal", mode = "bipartite",
		           R_row = 0, R_col = 0, nscan = 30, burn = 15, odens = 5,
		           rvar = FALSE, cvar = FALSE, verbose = FALSE,
		           ordinal_cutpoints = "explicit")
	},
		"not available for bipartite ordinal fits")
	# the fit still completes via data-induced fallback
	expect_null(fit$ALPHA)
})

test_that("ordinal_cutpoints = 'explicit' on non-ordinal family warns + ignores", {
	skip_on_cran()
	set.seed(51)
	n = 10
	Y = matrix(rnorm(n*n), n, n); diag(Y) = NA
	expect_warning({
		fit = ame(Y, family = "normal", R = 0, nscan = 20, burn = 10, odens = 5,
		          rvar = FALSE, cvar = FALSE, verbose = FALSE,
		          ordinal_cutpoints = "explicit")
	},
		"only takes effect")
	expect_equal(fit$ordinal_cutpoints, "data_induced")
})
