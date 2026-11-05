skip_on_cran()

poisson_summary_fixture = function(long = FALSE, bip = FALSE) {
	n = 3L; m = if (bip) 2L else n; S = 8L
	rn = letters[seq_len(n)]; cn = if (bip) LETTERS[seq_len(m)] else rn
	x = array(1, c(n, m, 2), dimnames = list(rn, cn, c('intercept', 'x_dyad')))
	x[, , 2] = matrix(seq_len(n * m) / 10, n, m)
	beta = cbind(intercept = seq(-.2, .2, length.out = S), x_dyad = rep(.4, S))
	y = matrix(2, n, m, dimnames = list(rn, cn)); if (!bip) diag(y) = NA_real_
	f = list(
		BETA = beta,
		VC = cbind(va = rep(0, S), cab = 0, vb = 0, rho = 0, ve = .6),
		family = 'poisson',
		mode = if (bip) 'bipartite' else 'unipartite',
		R = 0L,
		rvar = FALSE,
		cvar = FALSE,
		symmetric = FALSE,
		intercept = TRUE,
		nA = n,
		nB = m,
		APM = setNames(numeric(n), rn),
		BPM = setNames(numeric(m), cn),
		X = x,
		YPM = y
	)
	if (long) {
		f$Xlist = list(x, x); f$YPM = list(y, y)
		f$EZ = lapply(f$Xlist, function(z) .4 * z[, , 2])
		f$n_time = 2L; f$period_exposure = c(1, 2)
		f$dynamic_ab = TRUE; f$rho_ab = rep(.5, S); f$sigma_ab = rep(.1, S)
		f$forecast_state = list(
			a = matrix(0, n, S),
			b = matrix(0, m, S),
			a_level = matrix(0, n, S),
			b_level = matrix(0, m, S),
			draw_index = seq_len(S)
		)
	}
	class(f) = if (long) c('lame', 'ame') else 'ame'
	f
}

test_that('Poisson stored fitted values are returned silently and leave the caller RNG alone', {
	for (long in c(FALSE, TRUE)) for (bip in c(FALSE, TRUE)) {
		f = poisson_summary_fixture(long, bip)
		set.seed(73326); before = .Random.seed
		# no message on a first call or a later one, and no warning
		expect_silent({ a = fitted(f) })
		expect_silent({ b = predict(f) })
		expect_identical(a, f$YPM); expect_identical(b, f$YPM)
		expect_true(identical(.Random.seed, before))
		f$family = 'normal'
		expect_silent(expect_identical(fitted(f), f$YPM))
		expect_silent(expect_identical(predict(f), f$YPM))
	}
	# the package defines no one-time message helper
	expect_false(exists('.warn_poisson_prediction_summary', envir = asNamespace('lame')))
})

test_that('Poisson new-design response predictions are mean counts at the new predictor', {
	for (long in c(FALSE, TRUE)) for (bip in c(FALSE, TRUE)) {
		f = poisson_summary_fixture(long, bip)
		x = f$X[, , 2, drop = FALSE]
		nd = if (long) list(x, x) else x
		expect_silent({ a = predict(f, newdata = nd) })
		eta = predict(f, newdata = nd, type = 'link')
		# exp(eta + ve / 2) with the posterior mean ve = .6, times the exposure of the period
		expect_equal(
			a,
			if (long) Map(function(e, k) exp(e + .3) * c(1, 2)[k], eta, 1:2) else exp(eta + .3),
			tolerance = 1e-12
		)
	}
	# the inverse link alone carries no variance or exposure
	expect_equal(lame:::transform_to_response(matrix(.5, 2, 2), 'poisson'), exp(matrix(.5, 2, 2)))
	# a zero-exposure period is a structural zero
	f = poisson_summary_fixture(TRUE, FALSE); f$period_exposure = c(0, 2)
	x = f$X[, , 2, drop = FALSE]
	a = predict(f, newdata = list(x, x))
	expect_true(all(a[[1]] == 0)); expect_true(all(a[[2]] > 0))
	# the stored predictor path returns the same values when no YPM is stored
	f$YPM = NULL
	b = predict(f)
	expect_equal(b, Map(function(e, k) exp(e + .3) * c(0, 2)[k], f$EZ, 1:2), tolerance = 1e-12)
})

test_that('a Poisson fit without residual-variance draws cannot form a mean count', {
	f = poisson_summary_fixture(); f$VC = f$VC[, setdiff(colnames(f$VC), 've'), drop = FALSE]
	x = f$X[, , 2, drop = FALSE]
	expect_error(predict(f, newdata = x), class = 'lame_poisson_fitted')
	expect_error(lame:::.poisson_fitted_mean(matrix(0, 2, 2), NA_real_), class = 'lame_poisson_fitted')
	expect_error(lame:::.poisson_fitted_mean(matrix(800, 2, 2), 0), 'overflowed')
})

test_that('Poisson forecast points are mean counts at the draw-averaged predictor; draws keep their own means', {
	for (bip in c(FALSE, TRUE)) {
		f = poisson_summary_fixture(TRUE, bip)
		set.seed(63326); before = .Random.seed
		expect_silent({ draw = predict(f, h = 2, by_draw = TRUE, newexposure = c(0, 3), seed = 13) })
		eta = predict(f, h = 2, type = 'link', by_draw = TRUE, seed = 13)
		# each draw's conditional mean: e exp(eta_s + ve_s / 2)
		expect_equal(draw[, , 2, ], exp(eta[, , 2, ] + .3) * 3, tolerance = 1e-12)
		expect_true(all(draw[, , 1, ] == 0))
		# the point forecast: the mean count at the mean predictor and mean variance
		expect_silent({ avg = predict(f, h = 2, newexposure = c(0, 3), seed = 13) })
		expect_equal(avg[[2]], exp(apply(eta[, , 2, ], 1:2, mean) + .3) * 3, tolerance = 1e-12)
		expect_true(all(avg[[1]] == 0))
		expect_identical(dimnames(avg[[2]]), dimnames(draw)[1:2])
		expect_silent({ zero = predict(f, h = 2, newexposure = c(0, 0), seed = 13) })
		expect_true(all(unlist(zero)[is.finite(unlist(zero))] == 0))
		# credible intervals summarize the per-draw conditional means
		expect_silent({ ci = predict(f, h = 2, interval = 'credible', newexposure = c(0, 3), seed = 13) })
		expect_equal(ci[[2]]$median, apply(draw[, , 2, , drop = FALSE], 1:2, median), tolerance = 0)
		# prediction intervals summarize simulated counts
		expect_silent({ pi = predict(f, h = 2, interval = 'prediction', newexposure = c(0, 3), seed = 13) })
		expect_true(all(unlist(pi)[is.finite(unlist(pi))] >= 0))
		after = .Random.seed
		expect_identical(predict(f, h = 2, interval = 'prediction', newexposure = c(0, 3), seed = 13), pi)
		expect_true(identical(.Random.seed, after))
		expect_error(predict(f, h = 2, newexposure = c(-1, 2)), 'non-negative')
	}
})

test_that('joint Poisson count simulation remains available', {
	f = poisson_summary_fixture()
	set.seed(53326); before = .Random.seed
	expect_silent({ z = predict(f, type = 'distribution', n_samples = 50) })
	expect_identical(dim(z), c(3L, 3L, 50L))
	expect_true(all(z[is.finite(z)] >= 0 & z[is.finite(z)] == floor(z[is.finite(z)])))
	expect_false(identical(.Random.seed, before))
})
