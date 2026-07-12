# ghk / halton / gauss-hermite pointwise log-lik tests.
# cover the dispatcher, ghk on frn, gauss-hermite on cbin, and the
# reproducibility contract (deterministic given fit seed).

test_that(".halton_seq_randomised respects [0, 1) and decorrelates from base", {
	x = lame:::.halton_seq_randomised(50L, base = 2L, shift = 0.3)
	expect_true(all(x >= 0 & x < 1))
	expect_equal(length(x), 50L)
	# different shifts give different sequences
	y = lame:::.halton_seq_randomised(50L, base = 2L, shift = 0.7)
	expect_false(identical(x, y))
})

test_that(".ghk_row_loglik handles trivial cases cleanly", {
	# d = 0
	expect_equal(lame:::.ghk_row_loglik(numeric(0), integer(0),
	                                     n_mc = 32L), 0)
	# d = 1
	expect_equal(lame:::.ghk_row_loglik(c(0.5), c(1L), n_mc = 32L), 0)
})

test_that(".ghk_row_loglik returns finite log-prob for moderate D", {
	set.seed(91)
	ez = c(2.0, 1.5, 1.0, 0.5, 0.0, -0.5)
	perm = c(1L, 2L, 3L, 4L, 5L, 6L)   # already in descending order
	lp = lame:::.ghk_row_loglik(ez, perm, n_mc = 64L, halton_shift = 0.1)
	expect_true(is.finite(lp))
	expect_lt(lp, 0)  # log of a probability
})

test_that(".ghk_row_loglik falls back gracefully for D > 15", {
	set.seed(13)
	ez = rnorm(20)
	perm = seq_along(ez)
	lp = lame:::.ghk_row_loglik(ez, perm, n_mc = 32L)
	expect_true(is.finite(lp))
	expect_lt(lp, 0)
})

test_that(".cbin_gauss_hermite_loglik returns finite log-prob", {
	set.seed(17)
	# 8-dyad row with 3 "top" dyads
	ez = c(2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5)
	y = c(1, 1, 1, 0, 0, 0, 0, 0)
	lp = lame:::.cbin_gauss_hermite_loglik(y, ez, n_nodes = 10L)
	expect_true(is.finite(lp))
	expect_lt(lp, 0)
})

test_that(".pointwise_loglik_observed_ghk_rank dispatches by family", {
	set.seed(23)
	n_obs = 30
	y_obs = rbinom(n_obs, 1, 0.4)
	ez_obs = rnorm(n_obs, 0, 0.5)
	z_obs = rnorm(n_obs)
	obs_idx = cbind(i = rep(1:5, each = 6),
	                j = rep(1:6, 5),
	                t = rep(1L, n_obs))
	# binary -> closed-form (same as observed_exact)
	lp_bin = lame:::.pointwise_loglik_observed_ghk_rank(y_obs, ez_obs, z_obs,
	                                                   obs_idx, family = "binary",
	                                                   s2 = 1)
	expect_equal(length(lp_bin), n_obs)
	expect_true(all(is.finite(lp_bin)))
	# normal -> closed-form
	y_norm = rnorm(n_obs)
	lp_norm = lame:::.pointwise_loglik_observed_ghk_rank(y_norm, ez_obs, z_obs,
	                                                    obs_idx, family = "normal",
	                                                    s2 = 1)
	expect_equal(length(lp_norm), n_obs)
	expect_true(all(is.finite(lp_norm)))
	# cbin -> gauss-hermite
	lp_cbin = lame:::.pointwise_loglik_observed_ghk_rank(y_obs, ez_obs, z_obs,
	                                                    obs_idx, family = "cbin",
	                                                    s2 = 1)
	expect_true(all(is.finite(lp_cbin)))
})

test_that("v8 GHK on frn produces finite, family-correct log_lik", {
	skip_on_cran()
	set.seed(41)
	n = 10; T = 2
	Xdyad_list = lapply(seq_len(T), function(t) {
		x = array(rnorm(n*n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(T), function(t) {
		Y = matrix(0, n, n)
		for (i in 1:n) {
			nominees = sample(setdiff(1:n, i), 3)
			Y[i, nominees] = sample(1:3)
		}
		diag(Y) = NA
		Y
	})
	fit = suppressMessages(suppressWarnings(
		lame(Y_list, Xdyad = Xdyad_list, family = "frn", R = 1,
		     log_lik_method = "observed_ghk",
		     save_log_lik = TRUE,
		     nscan = 15, burn = 8, odens = 5,
		     verbose = FALSE, seed = 7)))
	expect_false(is.null(fit$log_lik))
	expect_true(all(is.finite(fit$log_lik)))
	expect_true(all(fit$log_lik < 0))
	# log_lik_meta records the deterministic halton shift
	expect_false(is.null(fit$log_lik_meta))
	expect_equal(fit$log_lik_meta$seed, 7)
	expect_true(fit$log_lik_meta$halton_shift >= 0 &&
	            fit$log_lik_meta$halton_shift < 1)
	expect_equal(fit$log_lik_meta$dim_cap, 15L)
})

test_that("v8 GHK is reproducible given seed (halton_shift deterministic)", {
	skip_on_cran()
	set.seed(43)
	n = 8; T = 2
	Xdyad_list = lapply(seq_len(T), function(t) {
		x = array(rnorm(n*n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(T), function(t) {
		Y = matrix(0, n, n)
		for (i in 1:n) {
			nominees = sample(setdiff(1:n, i), 2)
			Y[i, nominees] = sample(1:2)
		}
		diag(Y) = NA
		Y
	})
	f1 = suppressMessages(suppressWarnings(
		lame(Y_list, Xdyad = Xdyad_list, family = "frn", R = 0,
		     log_lik_method = "observed_ghk",
		     save_log_lik = TRUE,
		     nscan = 10, burn = 5, odens = 5,
		     rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 99)))
	f2 = suppressMessages(suppressWarnings(
		lame(Y_list, Xdyad = Xdyad_list, family = "frn", R = 0,
		     log_lik_method = "observed_ghk",
		     save_log_lik = TRUE,
		     nscan = 10, burn = 5, odens = 5,
		     rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 99)))
	expect_equal(f1$log_lik_meta$halton_shift, f2$log_lik_meta$halton_shift)
	# byte-identical log_lik matrices
	expect_identical(f1$log_lik, f2$log_lik)
})

test_that("default log_lik_method = 'observed_exact' is byte-identical", {
	skip_on_cran()
	set.seed(45)
	n = 8; T = 2
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rbinom(n*n, 1, 0.3), n, n)
		diag(mat) = NA
		mat
	})
	f1 = lame(Y_list, family = "binary", R = 0, save_log_lik = TRUE,
	          nscan = 10, burn = 5, odens = 5,
	          rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 13)
	f2 = lame(Y_list, family = "binary", R = 0,
	          log_lik_method = "observed_exact", save_log_lik = TRUE,
	          nscan = 10, burn = 5, odens = 5,
	          rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 13)
	expect_identical(f1$log_lik, f2$log_lik)
	expect_null(f1$log_lik_meta)  # only set when method = observed_ghk
})

test_that("GHK multi-dim Halton dimensions are decorrelated (no lattice)", {
	# each integration dimension must use its own coprime prime base;
	# a single shared base gives cross-dimension correlation ~ -0.71
	h = lame:::.halton_seq_randomised
	bases = lame:::.ghk_prime_bases
	d1 = h(256L, bases[1])
	d2 = h(256L, bases[2])
	expect_lt(abs(cor(d1, d2)), 0.1)
})
