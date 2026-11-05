skip_on_cran()

test_that('symmetric saved effects retain actor names and their shared draws', {
	set.seed(517);nm = paste0('n', 1:8)
	Y = lapply(1:3, function(t) {m = matrix(rnorm(64), 8);m = (m + t(m)) / 2;diag(m) = NA;dimnames(m) = list(nm, nm);m})
	names(Y) = c('P3', 'P1', 'P2')
	for (dyn in c(FALSE, TRUE)) for (rank in c(0L, 1L)) {
		f = suppressWarnings(
			suppressMessages(
				lame(
					Y,
					family = 'normal',
					symmetric = TRUE,
					R = rank,
					dynamic_ab = dyn,
					burn = 20,
					nscan = 60,
					odens = 2,
					seed = 517,
					gof = FALSE,
					plot = FALSE,
					verbose = FALSE,
					posterior_opts = list(save_ab = TRUE, thin_ab = 1, save_UV = rank > 0)
				)
			)
		)
		ab = suppressMessages(simulate_posterior(f, 'ab'))
		expect_identical(ab$a, ab$b)
		expect_identical(dimnames(ab$a)[[1]], nm)
		if (dyn) {
			expect_identical(dimnames(f$a_dynamic), list(nm, names(Y)))
			expect_identical(dimnames(f$a_dynamic_sd), dimnames(f$a_dynamic))
			expect_identical(dimnames(f$a_draws)[1:2], dimnames(f$a_dynamic))
			expect_identical(dimnames(ab$a)[1:2], dimnames(f$a_dynamic))
			expect_identical(unname(ab$a[, 'P1', , drop = FALSE]), unname(ab$a[, 2L, , drop = FALSE]))
		}
		if (dyn) { expect_equal(
			apply(ab$a, c(1, 2), mean),
			f$a_dynamic,
			tolerance = 1e-12
		);expect_identical(f$a_level_samples, f$b_level_samples) } else expect_equal(
			rowMeans(ab$a),
			f$APM,
			tolerance = 1e-12
		)
		for (tt in 1:3) {
			aa = if (dyn) ab$a[, tt, ] else ab$a
			expected = mean(f$BETA[, 1]) + outer(rowMeans(aa), rowMeans(aa), '+')
			if (rank > 0) expected = expected + f$ULUPM
			target = f$EZ[[tt]];keep = is.finite(target)
			expect_equal(unname(expected[keep]), unname(target[keep]), tolerance = 1e-8)
		}
	}
})


test_that('dynamic symmetric products can be selected by actor and period names', {
	set.seed(1834)
	actors = c('zeta', 'beta', 'theta', 'alpha', 'gamma', 'delta', 'eta')
	Y = setNames(lapply(seq_len(3), function(t) {
		y = matrix(rnorm(49), 7); y = (y + t(y)) / 2
		diag(y) = NA_real_; dimnames(y) = list(actors, actors); y
	}), c('2001', '2002', '2003'))
	for (family in c('normal', 'binary')) {
		response = if (family == 'binary') lapply(Y, function(y) 1 * (y > 0)) else Y
		fit = suppressMessages(lame(
			response,
			R = 2,
			family = family,
			symmetric = TRUE,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			burn = 20,
			nscan = 40,
			odens = 2,
			seed = 1834,
			gof = FALSE,
			plot = FALSE,
			verbose = FALSE,
			posterior_opts = list(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)
		))
		expect_identical(dimnames(fit$V_draws)[1:3], dimnames(fit$U_draws)[1:3])
		expect_setequal(dimnames(fit$V_draws)[[1]], actors)
		expect_identical(dimnames(fit$V_draws)[[3]], names(Y))
		order = rev(actors); t = '2002'; s = 3L
		named = tcrossprod(fit$U_draws[order, , t, s], fit$V_draws[order, , t, s])
		whole = tcrossprod(fit$U_draws[, , t, s], fit$V_draws[, , t, s])
		expect_equal(named, whole[order, order], tolerance = 1e-12)
		expect_equal(named, t(named), tolerance = 1e-12)
		expect_true(all(is.finite(named)))
	}
})
