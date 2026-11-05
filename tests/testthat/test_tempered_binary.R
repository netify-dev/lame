skip_on_cran()

.tempered_fixture = function() {
	n = 6L
	actors = c('f', 'a', 'd', 'b', 'e', 'c')
	Y = lapply(1:3, function(t) {
		y = outer(1:n, 1:n, function(i, j) as.numeric((i + 2 * j + t) %% 4 == 0))
		y[lower.tri(y)] = t(y)[lower.tri(y)]
		diag(y) = NA_real_
		dimnames(y) = list(actors, actors)
		y
	})
	names(Y) = paste0('t', 1:3)
	Y[[2]][1, 3] = Y[[2]][3, 1] = NA_real_
	X = lapply(
		1:3,
		function(t) array(
			abs(outer(1:n, 1:n, '-')) / n + t / 10,
			c(n, n, 1L),
			dimnames = list(actors, actors, 'distance')
		)
	)
	names(X) = names(Y)
	list(Y = Y, X = X)
}

.tempered_fit = function(...) {
	d = .tempered_fixture()
	args = list(
		Y = d$Y,
		Xdyad = d$X,
		R = 1L,
		family = 'binary',
		symmetric = TRUE,
		dynamic_ab = TRUE,
		sampler = 'tempered',
		seed = 60136L,
		burn = 20L,
		nscan = 40L,
		odens = 4L,
		gof = FALSE,
		verbose = FALSE,
		freeze_call = TRUE,
		prior = list(
			uv_prior = 'unit',
			Sab0 = diag(.3, 2),
			etaab = 6,
			rho_ab_mean = .5,
			rho_ab_sd = .2,
			sigma_ab_shape = 3,
			sigma_ab_scale = .2
		),
		posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE, thin_ab = 1)
	)
	dots = list(...)
	args[names(dots)] = dots
	do.call(lame, args)
}

test_that('tempered saved states reconstruct unique unordered likelihoods', {
	withr::local_seed(4136)
	before = .Random.seed
	f = .tempered_fit(save_log_lik = TRUE)
	expect_identical(.Random.seed, before)
	expect_s3_class(f, 'lame')
	expect_true(all(f$sampler_health$numerical_failures == 0))
	expect_equal(f$g, 36)
	expect_identical(dim(f$BETA), c(10L, 2L))
	expect_identical(dimnames(f$a_draws)[[1]], rownames(f$Y[[1]]))
	index = f$log_lik_index
	expect_true(all(index[, 1] < index[, 2]))
	expect_equal(nobs(f), sum(vapply(f$Y, function(y) sum(upper.tri(y) & is.finite(y)), 0L)))
	expect_equal(nobs(f), ncol(f$log_lik))
	for (s in 1:10) for (t in 1:3) {
		eta = matrix(matrix(f$Xlist[[t]], 36, 2) %*% f$BETA[s, ], 6, 6) +
			outer(f$a_draws[, t, s], f$a_draws[, t, s], '+') +
			tcrossprod(matrix(f$U_samples[, , s], 6, 1), matrix(f$V_samples[, , s], 6, 1))
		j = which(index[, 3] == t)
		cells = index[j, 1:2, drop = FALSE]
		expect_equal(f$log_lik[s, j], pnorm((2 * f$Y[[t]][cells] - 1) * eta[cells], log.p = TRUE), tolerance = 1e-10)
	}
	expect_equal(crossprod(f$U), diag(1), tolerance = 1e-10)
	expect_equal(tcrossprod(f$U, f$V), f$U %*% diag(f$L, nrow = 1L) %*% t(f$U), tolerance = 1e-12)
})

test_that('both tempered priors continue exactly, including coordinate restoration', {
	skip_on_cran()
	for (kind in c('unit', 'hierarchical')) {
		ck = tempfile(fileext = '.rds'); ck2 = tempfile(fileext = '.rds'); ck3 = tempfile(fileext = '.rds')
		withr::defer(unlink(c(ck, ck2, ck3)))
		pr = list(
			uv_prior = kind,
			Sab0 = diag(.3, 2),
			etaab = 6,
			rho_ab_mean = .5,
			rho_ab_sd = .2,
			sigma_ab_shape = 3,
			sigma_ab_scale = .2
		)
		first = .tempered_fit(prior = pr, checkpoint_path = ck)
		full = .tempered_fit(prior = pr, nscan = 80L, checkpoint_path = ck2)
		rest = lame_resume(ck, nscan_more = 40L, checkpoint_path = ck3)
		for (key in c('BETA', 'VC')) expect_identical(rbind(first[[key]], rest[[key]]), full[[key]])
		for (key in c('a_draws', 'b_draws', 'U_samples', 'V_samples'))
			expect_identical(abind::abind(first[[key]], rest[[key]], along = 3), full[[key]])
		a = readRDS(ck2)$state; b = readRDS(ck3)$state
		for (key in c('q', 'stepsizes', 'rng', 'scale_rng')) expect_identical(a[[key]], b[[key]])
		expect_error(lame_resume(ck, seed = 33), 'Start a new fit')
		expect_error(lame_resume(ck, prior = pr), 'Start a new fit')
	}
})

test_that('optional storage and GOF leave forecasts and inference unchanged', {
	full = .tempered_fit()
	lean = .tempered_fit(posterior_opts = NULL)
	thin = .tempered_fit(posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE, thin_ab = 3))
	with_gof = .tempered_fit(gof = TRUE)
	expect_null(lean$a_draws)
	expect_null(lean$U_samples)
	expect_identical(lean$forecast_state, full$forecast_state)
	expect_identical(thin$ab_draw_index, c(1L, 4L, 7L, 10L))
	expect_identical(thin$forecast_state, full$forecast_state)
	expect_identical(with_gof$BETA, full$BETA)
	expect_identical(with_gof$U_samples, full$U_samples)
	expect_length(with_gof$GOF, 5)
	withr::local_seed(995136)
	a = predict(full, h = 2, type = 'response', by_draw = TRUE)
	set.seed(995136)
	b = predict(lean, h = 2, type = 'response', by_draw = TRUE)
	expect_identical(a, b)
})

test_that('tempered inputs align actors and refits preserve the chosen target', {
	d = .tempered_fixture()
	f = .tempered_fit()
	d$Y[[2]] = d$Y[[2]][6:1, 6:1]
	aligned = .tempered_fit(Y = d$Y)
	expect_identical(aligned$BETA, f$BETA)
	expect_identical(aligned$U_samples, f$U_samples)
	cv = lfo(f, periods = 3L, nscan = 40L, burn = 20L, odens = 4L)
	expect_true(is.finite(cv$elpd_lfo))
	expect_equal(cv$per_period$n_obs, 15L)
	pr = capture.output(prior_summary(f))
	expect_true(any(grepl('fixed at one', pr, fixed = TRUE)))
	expect_false(any(grepl('inverse-Wishart', pr, fixed = TRUE)))
})

test_that('tempered sampler refuses unsupported or malformed model settings', {
	expect_error(.tempered_fit(dynamic_ab = FALSE), 'requires binary symmetric')
	expect_error(.tempered_fit(family = 'normal'), 'requires binary symmetric')
	expect_error(.tempered_fit(dynamic_uv = TRUE), 'requires binary symmetric')
	expect_error(.tempered_fit(R = 0), 'positive rank')
	expect_error(.tempered_fit(g = 0), 'finite positive scalar')
	expect_error(.tempered_fit(nscan = 3, odens = 2), 'multiple')
	expect_error(.tempered_fit(seed = Inf), 'integer')
	expect_error(.tempered_fit(max_seconds = 0), 'positive number')
	expect_error(.tempered_fit(save_log_lik = 'chunked'), 'TRUE or FALSE')
	expect_error(.tempered_fit(prior = list(uv_prior = 'typo')), 'uv_prior')
	expect_error(.tempered_fit(prior = list(uv_prior = 'unit', ab_var_cap = 1)), 'ab_var_cap')
	expect_error(.tempered_fit(prior = list(uv_prior = 'unit', Suv0 = 1)), 'documented tempered')
	expect_error(.tempered_fit(posterior_opts = list(save_UV = 'yes')), 'TRUE or FALSE')
	expect_error(.tempered_fit(posterior_opts = list(save_ab = TRUE, thin_ab = 0)), 'integer')
	expect_error(.tempered_fit(posterior_opts = list(misspelled = TRUE)), 'documented storage')
	d = .tempered_fixture(); d$Y[[1]][1, 2] = 1 - d$Y[[1]][2, 1]
	expect_error(.tempered_fit(Y = d$Y), 'symmetric')
	d = .tempered_fixture(); d$X[[1]][1, 2, 1] = Inf
	expect_error(.tempered_fit(Xdyad = d$X), 'finite off the diagonal')
	expect_error(lame(.tempered_fixture()$Y, family = 'binary', prior = list(uv_prior = 'unit')), 'only with')
})


test_that('tempered predictions reject asymmetric covariates on every entry point', {
	f = .tempered_fit()
	d = .tempered_fixture()
	bad = d$X; bad[[1]][1, 2, 1] = bad[[1]][1, 2, 1] + 1
	expect_error(predict(f, newdata = bad, type = 'link'), 'symmetric finite')
	expect_error(predict(f, h = 1, newdata = bad[1], type = 'link'), 'symmetric finite')
	expect_error(lame:::.joint_new_design(f, list(Xdyad = bad), t = 1), 'symmetric finite')
	x = d$X; x[[1]][1, 2, 1] = x[[1]][2, 1, 1] = NA_real_
	expect_error(predict(f, newdata = x, type = 'link'), 'symmetric finite')
	a = predict(f, newdata = d$X, type = 'link')
	expect_true(all(vapply(a, function(x) isTRUE(all.equal(x, t(x))), TRUE)))
	set.seed(9139); a = predict(f, h = 1, newdata = d$X[1], type = 'link', by_draw = TRUE)
	reversed = lapply(d$X[1], function(x) x[6:1, 6:1, , drop = FALSE])
	set.seed(9139); b = predict(f, h = 1, newdata = reversed, type = 'link', by_draw = TRUE)
	expect_identical(a, b)
})

test_that('tempered resume validates the saved RNG and accepts its own sampler selector', {
	skip_on_cran()
	ck = tempfile(fileext = '.rds'); bad = tempfile(fileext = '.rds')
	withr::defer(unlink(c(ck, bad)))
	first = .tempered_fit(checkpoint_path = ck)
	a = lame_resume(ck, nscan_more = 40L)
	b = lame(resume_from = ck, sampler = 'tempered', nscan = 40L)
	expect_identical(a$BETA, b$BETA)
	expect_identical(a$U_samples, b$U_samples)
	expect_error(lame(resume_from = ck, sampler = 'gibbs', nscan = 40L), 'Start a new fit')
	original = readRDS(ck)
	for (kind in c('missing', 'wrong_length', 'wrong_type', 'wrong_generator', 'missing_kind')) {
		z = original
		if (kind == 'missing') z$state$rng = NULL
		if (kind == 'wrong_length') z$state$rng = z$state$rng[1:3]
		if (kind == 'wrong_type') z$state$rng = as.numeric(z$state$rng)
		if (kind == 'wrong_generator') z$state$rng[1] = 10407L
		if (kind == 'missing_kind') z$state$rng_kind = NULL
		saveRDS(z, bad)
		withr::local_seed(92139)
		before = .Random.seed
		expect_error(lame_resume(bad, nscan_more = 40L), 'complete compatible RNG')
		expect_identical(.Random.seed, before)
	}
})


test_that('tempered RNG validation follows every serializable built-in generator', {
	skip_on_cran()
	withr::local_seed(14801)
	original = .Random.seed
	original_kind = RNGkind()
	withr::defer({
		suppressWarnings(do.call(RNGkind, as.list(original_kind)))
		assign('.Random.seed', original, envir = .GlobalEnv)
	})
	generators = c(
		'Wichmann-Hill', 'Marsaglia-Multicarry', 'Super-Duper',
		'Mersenne-Twister', 'Knuth-TAOCP', 'Knuth-TAOCP-2002', "L'Ecuyer-CMRG"
	)
	normals = c('Buggy Kinderman-Ramage', 'Ahrens-Dieter', 'Inversion', 'Kinderman-Ramage')
	binomials = if ('binom.kind' %in% names(formals(RNGkind))) c('Buggy BTPE', 'BTPE') else NA_character_
	for (generator in generators) for (normal in normals)
		for (sample in c('Rounding', 'Rejection')) for (binomial in binomials) {
			kind = c(generator, normal, sample)
			if (!is.na(binomial)) kind = c(kind, binomial)
			# R warns when selecting documented legacy generators; the validator must neither warn nor touch the random stream
			suppressWarnings(do.call(RNGkind, as.list(kind)))
			seed = .Random.seed
			expect_true(lame:::.lame_tempered_validate_rng(seed, kind))
			expect_identical(.Random.seed, seed)
			expect_identical(RNGkind(), kind)
			if (length(kind) == 4L && binomial == 'Buggy BTPE')
				expect_true(lame:::.lame_tempered_validate_rng(seed, kind[1:3]))
			if (length(kind) == 4L && binomial == 'BTPE')
				expect_error(lame:::.lame_tempered_validate_rng(seed, kind[1:3]), 'complete compatible RNG')
		}
})

test_that('tempered checkpoints retain both streams across a caller RNG change', {
	skip_on_cran()
	withr::local_seed(14802)
	original_kind = RNGkind()
	withr::defer(do.call(RNGkind, as.list(original_kind)))
	for (kind in c('unit', 'hierarchical')) {
		RNGkind("L'Ecuyer-CMRG", 'Inversion', 'Rejection')
		ck = tempfile(fileext = '.rds'); ck_full = tempfile(fileext = '.rds'); ck_rest = tempfile(fileext = '.rds')
		withr::defer(unlink(c(ck, ck_full, ck_rest)))
		pr = list(
			uv_prior = kind,
			Sab0 = diag(.3, 2),
			etaab = 6,
			rho_ab_mean = .5,
			rho_ab_sd = .2,
			sigma_ab_shape = 3,
			sigma_ab_scale = .2
		)
		first = .tempered_fit(prior = pr, checkpoint_path = ck)
		full = .tempered_fit(prior = pr, nscan = 80L, checkpoint_path = ck_full)
		RNGkind('Mersenne-Twister', 'Kinderman-Ramage', 'Rejection')
		set.seed(14803)
		caller = .Random.seed; caller_kind = RNGkind()
		rest = lame_resume(ck, nscan_more = 40L, checkpoint_path = ck_rest)
		expect_identical(.Random.seed, caller)
		expect_identical(RNGkind(), caller_kind)
		for (key in c('BETA', 'VC')) expect_identical(rbind(first[[key]], rest[[key]]), full[[key]])
		for (key in c('a_draws', 'b_draws', 'U_samples', 'V_samples'))
			expect_identical(abind::abind(first[[key]], rest[[key]], along = 3), full[[key]])
		a = readRDS(ck_full)$state; b = readRDS(ck_rest)$state
		for (key in c('q', 'stepsizes', 'rng', 'rng_kind', 'scale_rng')) expect_identical(a[[key]], b[[key]])
		saved = readRDS(ck)
		bad = saved
		bad$state$rng_kind = c(bad$state$rng_kind[1:3], 'unknown binomial generator')
		saveRDS(bad, ck_rest)
		expect_error(lame_resume(ck_rest, nscan_more = 40L), 'complete compatible RNG')
		expect_identical(.Random.seed, caller)
		expect_identical(RNGkind(), caller_kind)
		if (kind == 'hierarchical') {
			bad = saved; bad$state$scale_rng = bad$state$scale_rng[1:3]
			saveRDS(bad, ck_rest)
			expect_error(lame_resume(ck_rest, nscan_more = 40L), 'complete compatible RNG')
			expect_identical(.Random.seed, caller)
			expect_identical(RNGkind(), caller_kind)
		}
	}
})
