skip_on_cran()

# deliberately correlated blocks whose contribution cancels draw by draw
joint_fixture = function(long = FALSE, bip = FALSE) {
	n = 4L; m = if (bip) 3L else n; S = 12L; tt = if (long) 2L else 1L
	rn = paste0('r', seq_len(n)); cn = if (bip) paste0('c', seq_len(m)) else rn
	x = array(1, c(n, m, 1), dimnames = list(rn, cn, 'intercept'))
	q = seq_len(S)
	f = list(
		BETA = matrix(q, ncol = 1, dimnames = list(NULL, 'intercept')),
		VC = cbind(va = rep(1, S), cab = 0, vb = 1, rho = 0, ve = 0),
		a_samples = matrix(-q / 2, n, S, byrow = TRUE),
		b_samples = matrix(-q / 4, m, S, byrow = TRUE),
		U_samples = array(rep(q, n), c(S, n, 1)),
		V_samples = array(1, c(m, 1, S)),
		U = matrix(0, n, 1, dimnames = list(rn, NULL)),
		V = matrix(0, m, 1, dimnames = list(cn, NULL)),
		APM = setNames(rep(0, n), rn),
		BPM = setNames(rep(0, m), cn),
		rvar = TRUE,
		cvar = TRUE,
		R = 1L,
		mode = if (bip) 'bipartite' else 'unipartite',
		symmetric = FALSE,
		family = 'normal',
		X = x,
		YPM = matrix(0, n, m, dimnames = list(rn, cn)),
		nA = n,
		nB = m
	)
	f$U_samples = array(rep(-q / 4, each = n), c(n, 1, S))
	if (bip) { f$G_samples = array(rep(2, S), c(1, 1, S)); f$U_samples = f$U_samples / 2 }
	if (long) {
		f$a_draws = array(rep(f$a_samples, each = tt), c(tt, n, S));f$a_draws = aperm(f$a_draws, c(2, 1, 3))
		f$b_draws = array(rep(f$b_samples, each = tt), c(tt, m, S));f$b_draws = aperm(f$b_draws, c(2, 1, 3))
		f$U_draws = array(rep(f$U_samples, each = tt), c(tt, n, 1, S));f$U_draws = aperm(f$U_draws, c(2, 3, 1, 4))
		f$V_draws = array(rep(f$V_samples, each = tt), c(tt, m, 1, S));f$V_draws = aperm(f$V_draws, c(2, 3, 1, 4))
		f$a_samples = f$b_samples = f$U_samples = f$V_samples = NULL
		f$Xlist = rep(list(x), tt);f$YPM = rep(list(f$YPM), tt);f$n_time = tt
		f$dynamic_ab = f$dynamic_uv = TRUE
	}
	class(f) = if (long) c('lame', 'ame') else 'ame'; f
}

test_that('joint predictive outcomes preserve all parameter-block dependence', {
	for (long in c(FALSE, TRUE)) for (bip in c(FALSE, TRUE)) {
		f = joint_fixture(long, bip)
		a = simulate(f, nsim = 31, seed = 71, type = 'posterior', return_latent = TRUE)
		expect_identical(a$method, 'joint posterior')
		expect_length(a$draw_index, 31)
		expect_true(all(unlist(a$Y)[is.finite(unlist(a$Y))] == 0))
		expect_identical(a$Y, a$Z)
		# the default route is the posterior one, and "auto" is its alias
		expect_identical(simulate(f, nsim = 31, seed = 71)$Y, a$Y)
		expect_identical(simulate(f, nsim = 31, seed = 71, type = 'auto')$Y, a$Y)
	}
})

test_that('burning and thinning select shared original iteration indices', {
	f = joint_fixture();f$ab_draw_index = c(2L, 4L, 6L, 8L, 10L, 12L)
	f$a_samples = f$a_samples[, f$ab_draw_index];f$b_samples = f$b_samples[, f$ab_draw_index]
	f$uv_draw_index = c(3L, 6L, 9L, 12L)
	f$U_samples = f$U_samples[, , f$uv_draw_index, drop = FALSE];f$V_samples = f$V_samples[, , f$uv_draw_index, drop = FALSE]
	a = simulate(f, nsim = 50, type = 'posterior', burn_in = 5, thin = 2, seed = 41)
	expect_setequal(a$draw_index, c(6L, 12L))
	expect_true(all(unlist(a$Y)[is.finite(unlist(a$Y))] == 0))
	f$ab_draw_index = NULL
	expect_error(simulate(f, type = 'posterior'), class = 'lame_predictive_draws')
})

test_that('predictive draws have the specified outcomes and shared latent noise', {
	skip_on_cran()
	f = joint_fixture();f$family = 'binary';f$VC[, 'rho'] = .65
	a = simulate(f, nsim = 1200, type = 'posterior', return_latent = TRUE, seed = 171)
	expect_identical(a$Y, lapply(a$Z, function(z) 1 * (z > 0)))
	expect_equal(mean(vapply(a$Y, function(y) y[1, 2], numeric(1))), .5, tolerance = .04)
	z12 = vapply(a$Z, function(z) z[1, 2], numeric(1));z21 = vapply(a$Z, function(z) z[2, 1], numeric(1))
	expect_equal(cor(z12, z21), .65, tolerance = .055)
	p = predict(f, type = 'distribution', n_samples = 20)
	expect_true(all(p[is.finite(p)] %in% c(0, 1)))
	expect_length(attr(p, 'draw_index'), 20)
	f$family = 'poisson';f$VC[, 've'] = .6
	y = simulate(f, nsim = 2500, type = 'posterior', seed = 19)
	v = unlist(y$Y);v = v[is.finite(v)]
	expect_equal(mean(v), exp(.3), tolerance = .04)
	expect_equal(var(v), exp(.3) + (exp(.6) - 1) * exp(.6), tolerance = .15)
})

test_that('a squared bipartite grid retains its diagonal with independent noise', {
	f = joint_fixture();f$mode = 'bipartite';f$G_samples = array(1, c(1, 1, 12));f$nA = f$nB = 4
	f$family = 'normal';f$VC[, 've'] = 2;f$VC[, 'rho'] = .9
	a = simulate(f, nsim = 1500, type = 'posterior', return_latent = TRUE, seed = 32)
	expect_true(all(vapply(a$Y, function(y) all(is.finite(diag(y))), logical(1))))
	expect_equal(var(vapply(a$Y, function(y) y[1, 1], numeric(1))), 2, tolerance = .2)
	expect_lt(abs(cor(vapply(a$Y, function(y) y[1, 2], numeric(1)), vapply(a$Y, function(y) y[2, 1], numeric(1)))), .09)
})

test_that('missing storage is explicit and probability utilities refuse it', {
	f = joint_fixture();f$a_samples = NULL
	expect_error(simulate(f, type = 'posterior'), class = 'lame_predictive_draws')
	expect_error(simulate_posterior(f, 'Y'), class = 'lame_predictive_draws')
	expect_error(predict(f, type = 'distribution'), class = 'lame_predictive_draws')
	# a fit without joint draws is never simulated from a substitute unasked: the
	# default and the "auto" alias abort and name the approximate route
	for (ty in c('posterior', 'auto')) {
		e = expect_error(simulate(f, nsim = 1, type = ty), class = 'lame_simulate_draws')
		expect_s3_class(e, 'lame_predictive_draws')
		expect_match(conditionMessage(e), 'type = "approximate"', fixed = TRUE)
		expect_match(conditionMessage(e), 'posterior_options')
	}
	expect_error(simulate(f), class = 'lame_simulate_draws')
	expect_no_message({x = simulate(f, nsim = 1, type = 'approximate')})
	expect_identical(x$method, 'conditional approximation')
})

test_that('posterior replication refuses legacy period options and names the approximate route', {
	f = joint_fixture(TRUE, FALSE)
	for (ty in c('posterior', 'auto')) {
		e = expect_error(simulate(f, nsim = 1, type = ty, n_time = 3), class = 'lame_simulate_periods')
		expect_match(conditionMessage(e), 'type = "approximate"', fixed = TRUE)
		expect_error(simulate(f, nsim = 1, type = ty, start_from = 'random'), class = 'lame_simulate_periods')
	}
	expect_error(simulate(f, nsim = 1, n_time = 3), class = 'lame_simulate_periods')
})

test_that('the conditional approximation needs one covariate array per period', {
	f = joint_fixture(TRUE);f$dynamic_ab = f$dynamic_uv = FALSE;f$BETA[] = 0;f$R = 0L
	f$VC[, 've'] = 1;f$VC[, 'rho'] = 0
	x = array(1, c(4, 4, 1), dimnames = list(paste0('r', 1:4), paste0('r', 1:4), 'intercept'))
	# a bare array is not a list of periods: its length counts cells
	e = expect_error(
		simulate(f, nsim = 1, type = 'approximate', newdata = list(Xdyad = x)),
		class = 'lame_simulate_xlist'
	)
	expect_match(conditionMessage(e), 'one dyadic covariate array per period')
	expect_error(
		simulate(f, nsim = 1, type = 'approximate', newdata = list(Xdyad = list())),
		class = 'lame_simulate_xlist'
	)
	# a list is cycled or cut to the requested periods
	expect_length(simulate(f, nsim = 1, type = 'approximate', n_time = 3, newdata = list(Xdyad = list(x)))$Y[[1]], 3)
	expect_length(simulate(f, nsim = 1, type = 'approximate', n_time = 1, newdata = list(Xdyad = list(x, x)))$Y[[1]], 1)
})

test_that('actual thinned fits retain indices used by joint prediction', {
	set.seed(36);y = matrix(rnorm(64), 8, 8);diag(y) = NA
	rownames(y) = colnames(y) = paste0('a', 1:8)
	f = suppressMessages(ame(
		y,
		family = 'normal',
		R = 1,
		burn = 10,
		nscan = 48,
		odens = 2,
		seed = 12,
		gof = FALSE,
		verbose = FALSE,
		posterior_opts = list(save_ab = TRUE, save_UV = TRUE, thin_ab = 2, thin_UV = 3)
	))
	expect_identical(as.integer(f$ab_draw_index), seq(2L, 24L, 2L))
	expect_identical(as.integer(f$uv_draw_index), seq(3L, 24L, 3L))
	out = simulate(f, type = 'posterior', nsim = 30, seed = 16)
	expect_true(all(out$draw_index %in% c(6, 12, 18, 24)))
	s = suppressMessages(ame(
		(y + t(y)) / 2,
		family = 'normal',
		R = 1,
		symmetric = TRUE,
		burn = 10,
		nscan = 30,
		odens = 2,
		gof = FALSE,
		verbose = FALSE,
		seed = 11,
		posterior_opts = list(save_ab = TRUE, save_UV = TRUE)
	))
	yy = simulate(s, type = 'posterior', nsim = 3)$Y
	expect_true(all(vapply(yy, function(z) isTRUE(all.equal(z, t(z))), logical(1))))
})

test_that('joint newdata maps named covariates and actors without replacing other draws', {
	f = joint_fixture(TRUE, TRUE)
	for (t in 1:2) f$Xlist[[t]] = array(
		c(rep(1, 12), rep(2, 12)),
		c(4, 3, 2),
		dimnames = list(paste0('r', 1:4), paste0('c', 1:3), c('intercept', 'x_dyad'))
	)
	f$BETA = cbind(f$BETA, x_dyad = rep(3, 12))
	x = array(4, c(4, 3, 1), dimnames = list(paste0('r', 4:1), paste0('c', 3:1), 'x'))
	a = simulate(f, nsim = 3, type = 'posterior', newdata = list(Xdyad = list(x, x)))
	expect_true(all(unlist(a$Y) == 12))
	expect_error(simulate(f, nsim = 1, type = 'posterior', n_time = 3), 'fitted periods')
	expect_error(simulate(f, thin = 0), 'positive integer')
})

test_that('ordinal cutpoints and nomination budgets use the same latent field', {
	f = joint_fixture();f$family = 'ordinal';f$ordinal_cutpoints = 'explicit';f$ordinal_levels = 0:2
	f$ALPHA = matrix(rep(.8, 12), ncol = 1)
	a = simulate(f, nsim = 30, type = 'posterior', return_latent = TRUE, seed = 33)
	for (i in seq_along(a$Y)) expect_equal(
		a$Y[[i]],
		matrix(findInterval(as.vector(a$Z[[i]]), c(0, .8), left.open = TRUE), 4, 4, dimnames = dimnames(a$Z[[i]]))
	)
	f$family = 'frn';f$odmax = 1
	a = simulate(f, nsim = 30, type = 'posterior', return_latent = TRUE, seed = 34)
	for (i in seq_along(a$Y)) {
		expect_true(all(rowSums(a$Y[[i]] > 0, na.rm = TRUE) <= 1))
		expect_true(all(a$Z[[i]][which(a$Y[[i]] > 0)] > 0))
	}
})

test_that('zero exposure retains structural diagonal missingness', {
	f = joint_fixture(TRUE);f$family = 'poisson';f$period_exposure = c(0, 0)
	a = simulate(f, type = 'posterior', nsim = 3)
	expect_true(all(vapply(unlist(a$Y, recursive = FALSE), function(y) all(is.na(diag(y))), logical(1))))
	expect_true(all(unlist(a$Y)[is.finite(unlist(a$Y))] == 0))
})

test_that('fitted-period replications use their matched reciprocity path', {
	f = joint_fixture(TRUE);f$VC[, 've'] = 1;f$VC[, 'rho'] = 0
	f$dynamic_rho = TRUE;f$RHO = cbind(rep(1, 12), rep(-1, 12))
	z = simulate(f, nsim = 8, seed = 417, type = 'posterior')$Y
	for (path in z) {
		expect_equal(path[[1]], t(path[[1]]), tolerance = 1e-12)
		expect_equal(path[[2]], -t(path[[2]]), tolerance = 1e-12)
	}
	# the RHO iterations are selected with BETA, after thinning as well
	f$RHO = cbind(rep(c(-1, 1), 6), rep(c(1, -1), 6))
	z = simulate(f, nsim = 18, seed = 421, type = 'posterior', burn_in = 2, thin = 3)
	for (j in seq_along(z$Y)) for (t in 1:2)
		expect_equal(z$Y[[j]][[t]], f$RHO[z$draw_index[j], t] * t(z$Y[[j]][[t]]), tolerance = 1e-12)
	f$RHO = f$RHO[-1, , drop = FALSE]
	expect_error(simulate(f, type = 'posterior'), class = 'lame_predictive_draws')
})

test_that('future outcomes preserve reciprocity, symmetry and bipartite independence', {
	f = joint_fixture(TRUE);f$dynamic_ab = f$dynamic_uv = FALSE;f$VC[, 've'] = 2
	for (rho in c(-1, 1)) {
		f$VC[, 'rho'] = rho
		z = suppressWarnings(lame:::.forecast_lame(f, h = 2, predictive = TRUE, by_draw = TRUE, seed = 25))
		for (k in 1:2) for (j in seq_len(dim(z)[4]))
			expect_equal(z[, , k, j], rho * t(z[, , k, j]), tolerance = 1e-12)
	}
	f$symmetric = TRUE;f$VC[, 'rho'] = 0
	z = suppressWarnings(lame:::.forecast_lame(f, h = 2, predictive = TRUE, by_draw = TRUE, seed = 31))
	expect_equal(z[, , 1, 1], t(z[, , 1, 1]), tolerance = 1e-12)
	expect_true(all(is.na(diag(z[, , 1, 1]))))
	f$symmetric = FALSE;f$mode = 'bipartite';f$G_samples = array(1, c(1, 1, 12));f$nA = f$nB = 4
	z = suppressWarnings(lame:::.forecast_lame(f, h = 2, predictive = TRUE, by_draw = TRUE, seed = 31))
	expect_true(all(is.finite(z)))
	expect_false(isTRUE(all.equal(z[, , 1, 1], t(z[, , 1, 1]))))
})

test_that('future noise refuses undefined reciprocity and invalid or clipped rates', {
	f = joint_fixture(TRUE);f$dynamic_ab = f$dynamic_uv = FALSE;f$VC[, 've'] = 1
	f$dynamic_rho = TRUE;f$RHO = matrix(.7, 12, 2)
	expect_error(suppressWarnings(lame:::.forecast_lame(f, h = 1, predictive = TRUE)), class = 'lame_forecast_model')
	expect_no_error(suppressWarnings(lame:::.forecast_lame(f, h = 1, type = 'link')))
	f$dynamic_rho = FALSE;f$VC[, 've'] = -1
	expect_error(suppressWarnings(lame:::.forecast_lame(f, h = 1, predictive = TRUE)), 'variance')
	f$family = 'poisson';f$VC[, 've'] = 0;f$BETA = f$BETA + 1000
	expect_error(suppressWarnings(lame:::.forecast_lame(f, h = 1, predictive = TRUE, type = 'response')), 'overflowed')
	expect_error(suppressWarnings(lame:::.forecast_lame(f, h = 1, type = 'response')), 'overflowed')
	for (predictive in c(FALSE, TRUE)) {
		z = suppressWarnings(
			lame:::.forecast_lame(
				f,
				h = 2,
				type = 'response',
				newexposure = c(0, 0),
				predictive = predictive,
				by_draw = TRUE
			)
		)
		expect_true(all(z[is.finite(z)] == 0))
		if (predictive) expect_true(all(is.na(diag(z[, , 1, 1]))))
	}
})

test_that('explicit legacy approximation keeps the fitted reciprocity periods', {
	f = joint_fixture(TRUE);f$dynamic_ab = f$dynamic_uv = FALSE;f$BETA[] = 0;f$R = 0L
	f$VC[, 've'] = 1;f$VC[, 'rho'] = 0;f$dynamic_rho = TRUE;f$RHO = cbind(rep(1, 12), rep(-1, 12))
	z = simulate(f, nsim = 3, seed = 51, type = 'approximate', burn_in = 2, thin = 3)
	expect_identical(z$method, 'conditional approximation')
	for (y in z$Y) {
		expect_equal(y[[1]], t(y[[1]]), tolerance = 1e-12)
		expect_equal(y[[2]], -t(y[[2]]), tolerance = 1e-12)
	}
	expect_error(simulate(f, type = 'approximate', n_time = 3), 'future reciprocity model')
})


test_that('seeded predictive simulations preserve caller RNG on return and error', {
	for (long in c(FALSE, TRUE)) for (bip in c(FALSE, TRUE)) {
		f = joint_fixture(long, bip)
		set.seed(71183);before = .Random.seed
		out = simulate(f, nsim = 5, seed = 3471, type = 'posterior', return_latent = TRUE)
		expect_identical(.Random.seed, before)
		set.seed(3471)
		reference = .simulate_joint(f, .joint_predictive_state(f), 5, if (long) 2L else 1L, TRUE)
		expect_identical(out, reference)
		before = .Random.seed
		expect_error(simulate(f, nsim = 0, seed = 831, type = 'posterior'), 'positive integer')
		expect_identical(.Random.seed, before)
		invisible(simulate(f, nsim = 2, type = 'posterior'))
		expect_false(identical(.Random.seed, before))
		saved = .Random.seed
		tryCatch({
			rm('.Random.seed', envir = globalenv())
			invisible(simulate(f, nsim = 2, seed = 832, type = 'posterior'))
			expect_false(exists('.Random.seed', envir = globalenv(), inherits = FALSE))
		}, finally = assign('.Random.seed', saved, envir = globalenv()))
	}
})
