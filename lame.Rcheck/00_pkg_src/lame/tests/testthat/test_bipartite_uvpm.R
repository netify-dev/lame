skip_on_cran()

# bipartite uvpm posterior mean and multi-chain pooling of the
# multiplicative-effects posterior means

# fixed bipartite data with rank-2 structure shared by the tests below
.uvpm_bip_data = function() {
	set.seed(99)
	nA = 15
	nB = 8
	U0 = matrix(rnorm(nA * 2), nA, 2)
	V0 = matrix(rnorm(nB * 2), nB, 2)
	Y = U0 %*% t(V0) + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	rownames(Y) = paste0("r", sprintf("%02d", 1:nA))
	colnames(Y) = paste0("c", sprintf("%02d", 1:nB))
	Y
}

test_that("bipartite ame returns UVPM as a named posterior mean", {
	skip_on_cran()

	Y = .uvpm_bip_data()

	fit_bip = function(seed, nscan) {
		ame(Y, mode = "bipartite", family = "normal", R_row = 2, R_col = 2,
				seed = seed, nscan = nscan, burn = 200, odens = 10,
				verbose = FALSE, gof = FALSE, plot = FALSE)
	}

	seeds = c(1, 2, 3)
	fits_short = lapply(seeds, function(s) fit_bip(s, 500))
	fits_long = lapply(seeds, function(s) fit_bip(s, 5000))

	f1 = fits_long[[1]]
	expect_true(is.matrix(f1$UVPM))
	expect_equal(dim(f1$UVPM), c(nrow(Y), ncol(Y)))
	expect_identical(rownames(f1$UVPM), rownames(Y))
	expect_identical(colnames(f1$UVPM), colnames(Y))

	# a posterior mean tightens across seeds as the chain lengthens; a
	# single retained draw would keep the across-seed spread flat
	sd_across = function(mats) {
		arr = simplify2array(mats)
		mean(apply(arr, c(1, 2), sd))
	}
	sd_short = sd_across(lapply(fits_short, `[[`, "UVPM"))
	sd_long = sd_across(lapply(fits_long, `[[`, "UVPM"))
	expect_lt(sd_long, sd_short)

	# the accessor returns the stored posterior mean
	expect_identical(reconstruct_UVPM(f1), f1$UVPM)
})

test_that("bipartite ame at R_row = R_col = 0 returns an all-zero UVPM", {
	skip_on_cran()

	Y = .uvpm_bip_data()

	fit0 = ame(Y, mode = "bipartite", family = "normal", R_row = 0, R_col = 0,
						 seed = 1, nscan = 300, burn = 100, odens = 10,
						 verbose = FALSE, gof = FALSE, plot = FALSE)

	# amen-compatible: the slot is present and zero, not absent
	expect_false(is.null(fit0$UVPM))
	expect_true(all(fit0$UVPM == 0))
	expect_equal(dim(fit0$UVPM), c(nrow(Y), ncol(Y)))
})

test_that("ame_parallel bipartite pools the stored per-chain UVPMs", {
	skip_on_cran()

	Y = .uvpm_bip_data()

	chains = suppressMessages(
		ame_parallel(Y, n_chains = 2, cores = 1, combine_method = "list",
								 mode = "bipartite", family = "normal",
								 R_row = 2, R_col = 2,
								 nscan = 500, burn = 100, odens = 10,
								 verbose = FALSE, gof = FALSE)
	)
	expect_false(is.null(chains[[1]]$UVPM))
	expect_false(is.null(chains[[2]]$UVPM))

	pooled = suppressMessages(combine_ame_chains(chains))
	avg_uvpm = (chains[[1]]$UVPM + chains[[2]]$UVPM) / 2
	expect_equal(pooled$UVPM, avg_uvpm)

	# the pooled value is the average of stored posterior means, not the
	# average of the chains' final u g v' draws
	last_draw_avg = (chains[[1]]$U %*% chains[[1]]$G %*% t(chains[[1]]$V) +
									 chains[[2]]$U %*% chains[[2]]$G %*% t(chains[[2]]$V)) / 2
	expect_gt(norm(pooled$UVPM - last_draw_avg, "F"), 1e-8)
})

test_that("ame_parallel symmetric unipartite pools ULUPM across chains", {
	skip_on_cran()

	set.seed(7)
	n = 12
	Us = matrix(rnorm(n * 2), n, 2)
	Ys = Us %*% t(Us) + matrix(rnorm(n * n, 0, 0.5), n, n)
	Ys = (Ys + t(Ys)) / 2
	diag(Ys) = NA
	rownames(Ys) = colnames(Ys) = paste0("a", sprintf("%02d", 1:n))

	chains = suppressMessages(
		ame_parallel(Ys, n_chains = 2, cores = 1, combine_method = "list",
								 family = "normal", symmetric = TRUE, R = 2,
								 nscan = 500, burn = 100, odens = 10,
								 verbose = FALSE, gof = FALSE)
	)
	expect_false(is.null(chains[[1]]$ULUPM))
	expect_false(is.null(chains[[2]]$ULUPM))

	pooled = suppressMessages(combine_ame_chains(chains))
	avg_ulupm = (chains[[1]]$ULUPM + chains[[2]]$ULUPM) / 2
	expect_equal(pooled$ULUPM, avg_ulupm)

	# pooling combines both chains rather than passing chain 1 through
	expect_gt(norm(pooled$ULUPM - chains[[1]]$ULUPM, "F"), 1e-8)
})
