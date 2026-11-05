skip_on_cran()

effect_panel = function(bip = FALSE) {
	set.seed(590)
	nr = 8L; nc = if (bip) 7L else nr; nt = 3L
	Y = array(
		rnorm(nr * nc * nt),
		c(nr, nc, nt),
		dimnames = list(paste0("a", seq_len(nr)), paste0(if (bip) "b" else "a", seq_len(nc)), c("P3", "P1", "P2"))
	)
	X = array(
		rnorm(nr * nc * nt),
		c(nr, nc, 1L, nt),
		dimnames = list(dimnames(Y)[[1]], dimnames(Y)[[2]], "x", dimnames(Y)[[3]])
	)
	if (!bip) for (tt in seq_len(nt)) diag(Y[, , tt]) = NA_real_
	list(
		Y = Y,
		Xdyad = X,
		family = "normal",
		mode = if (bip) "bipartite" else "unipartite",
		burn = 20L,
		nscan = 40L,
		odens = 2L,
		seed = 590L,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
}

test_that("saved panel effects reproduce fitted means and joint predictors", {
	skip_on_cran()
	for (bip in c(FALSE, TRUE)) for (dynamic_ab in c(FALSE, TRUE)) for (dynamic_uv in c(FALSE, TRUE)) {
		a = effect_panel(bip); a$R = 1L; a$dynamic_ab = dynamic_ab; a$dynamic_uv = dynamic_uv
		a$dynamic_beta = "x_dyad"
		a$posterior_opts = list(save_UV = TRUE, save_ab = TRUE, thin_ab = 1L)
		f = suppressWarnings(suppressMessages(do.call(lame, a)))
		a$posterior_opts$save_ab = FALSE
		quiet = suppressWarnings(suppressMessages(do.call(lame, a)))
		expect_identical(f$BETA, quiet$BETA)
		expect_identical(f$VC, quiet$VC)
		ab = suppressMessages(simulate_posterior(f, "ab"))
		uv = suppressMessages(simulate_posterior(f, "UV"))
		S = dim(f$BETA)[1]; nr = dim(a$Y)[1]; nc = dim(a$Y)[2]; nt = dim(a$Y)[3]
		if (dynamic_ab) {
			periods = dimnames(a$Y)[[3L]]
			expect_identical(dimnames(f$a_dynamic), list(names(f$APM), periods))
			expect_identical(dimnames(f$b_dynamic), list(names(f$BPM), periods))
			expect_identical(dimnames(f$a_dynamic_sd), dimnames(f$a_dynamic))
			expect_identical(dimnames(f$b_dynamic_sd), dimnames(f$b_dynamic))
			expect_identical(dimnames(f$a_draws)[1:2], dimnames(f$a_dynamic))
			expect_identical(dimnames(f$b_draws)[1:2], dimnames(f$b_dynamic))
			for (tt in seq_along(periods)) {
				expect_identical(unname(ab$a[, periods[tt], , drop = FALSE]), unname(ab$a[, tt, , drop = FALSE]))
				expect_identical(unname(ab$b[, periods[tt], , drop = FALSE]), unname(ab$b[, tt, , drop = FALSE]))
			}
			expect_equal(apply(ab$a, c(1, 2), mean), f$a_dynamic, tolerance = 1e-12)
			expect_equal(apply(ab$b, c(1, 2), mean), f$b_dynamic, tolerance = 1e-12)
			expect_equal(unname(rowMeans(f$a_level_samples)), unname(f$a_level), tolerance = 1e-12)
		} else {
			expect_equal(rowMeans(ab$a), f$APM, tolerance = 1e-12)
			expect_equal(rowMeans(ab$b), f$BPM, tolerance = 1e-12)
		}
		for (tt in seq_len(nt)) {
			eta = array(0, c(nr, nc, S))
			for (d in seq_len(S)) {
				av = if (dynamic_ab) ab$a[, tt, d] else ab$a[, d]
				bv = if (dynamic_ab) ab$b[, tt, d] else ab$b[, d]
				product = if (dynamic_uv) uv[, , tt, d] else uv[, , d]
				eta[, , d] = outer(av, bv, "+") + product
				for (k in seq_len(dim(f$BETA)[2])) eta[, , d] = eta[, , d] + f$Xlist[[tt]][, , k] * f$BETA[d, k, tt]
			}
			target = f$EZ[[tt]]; average = apply(eta, c(1, 2), mean)
			expect_equal(unname(average[is.finite(target)]), unname(target[is.finite(target)]), tolerance = 1e-8)
		}
		qs = suppressMessages(posterior_quantiles(f, "ab", probs = c(.1, .9)))
		expected = apply(ab$a, seq_len(length(dim(ab$a)) - 1L), quantile, probs = c(.1, .9))
		expect_equal(unname(qs$a), unname(expected), tolerance = 0)
	}
})

test_that("additive thinning retains matching draws and handles empty storage", {
	a = effect_panel(); a$R = 0L; a$dynamic_ab = TRUE
	a$posterior_opts = list(save_ab = TRUE, thin_ab = 1L)
	all = suppressMessages(do.call(lame, a))
	a$posterior_opts$thin_ab = 3L
	thin = suppressMessages(do.call(lame, a))
	idx = seq(3L, 20L, 3L)
	expect_equal(thin$ab_draw_index, idx)
	expect_equal(unname(thin$a_draws), unname(all$a_draws[, , idx, drop = FALSE]), tolerance = 0)
	expect_identical(thin$BETA, all$BETA)
	a$max_seconds = 1e-12
	empty = suppressWarnings(suppressMessages(do.call(lame, a)))
	expect_equal(dim(empty$a_draws)[3], 0L)
	expect_length(empty$ab_draw_index, 0L)
	a$max_seconds = Inf; a$posterior_opts$thin_ab = 0
	expect_error(do.call(lame, a), "positive integer")
})


test_that("longitudinal additive extraction requires actual saved draws", {
	fit = structure(list(BETA = matrix(0, 10, 1), APM = c(0, 1), BPM = c(1, 0)), class = "lame")
	expect_error(simulate_posterior(fit, "ab"), "Additive posterior draws were not saved")
})
