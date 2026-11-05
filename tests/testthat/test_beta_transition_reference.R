transition_reference_fixture = function(kind = 'rw1', periods = 5L, draws = 30L) {
	nms = c('intercept', 'x_dyad')
	X = array(c(rep(1, 9), rep(.5, 9)), c(3, 3, 2), dimnames = list(letters[1:3], letters[1:3], nms))
	B = array(0, c(draws, 2, periods), dimnames = list(NULL, nms, as.character(seq_len(periods))))
	B[, 2, periods] = 16
	list(
		dynamic_beta = TRUE,
		dynamic_beta_kind = kind,
		BETA = B,
		beta_dynamic_mask = c(FALSE, TRUE),
		beta_dynamic_groups = c('', 'dyad'),
		RHO_BETA = matrix(if (kind == 'rw1') 1 else 0, draws, 1, dimnames = list(NULL, 'dyad')),
		SIGMA_BETA = matrix(2, draws, 1, dimnames = list(NULL, 'dyad')),
		Xlist = rep(list(X), periods),
		prior = list(beta0_mean = 0, beta0_var = 10)
	)
}

test_that('change-point reference has the fitted number of RW1 transitions and scale', {
	skip_on_cran()
	f = transition_reference_fixture(periods = 2L)
	z = detect_change_point(f, n_prior_sims = 40000, seed = 41)
	expect_identical(z$coef, 'x_dyad')
	expect_equal(z$m_post_mean, 4)
	expect_equal(z$t_hat, 2L)
	expect_equal(z$bf, 20)
	expect_equal(z$m_prior_q95, qnorm(.975), tolerance = .025)
	f = transition_reference_fixture(periods = 5L)
	z = detect_change_point(f, n_prior_sims = 40000, seed = 42)
	expect_equal(z$m_prior_q95, qnorm((1 + .95^(1 / 4)) / 2), tolerance = .025)
	expect_equal(z$t_hat, 5L)
})

test_that('change-point reference for AR1 starts from the first state of the matched draw', {
	skip_on_cran()
	# a reference path begins where the fitted path begins: with rho = 0 the sole
	# difference of two periods is the fresh innovation less the first state
	f = transition_reference_fixture(kind = 'ar1', periods = 2L)
	z = detect_change_point(f, n_prior_sims = 40000, seed = 43)
	expect_equal(z$m_prior_q95, qnorm(.975), tolerance = .025)
	# a first state away from zero moves the reference: |N(-b, 1)| with b the
	# first state in units of the innovation scale (4 in this fixture)
	f$BETA[, 2, 1] = 8
	z = detect_change_point(f, n_prior_sims = 40000, seed = 44)
	q = uniroot(function(q) pnorm(q + 2) - pnorm(-q + 2) - .95, c(0, 10))$root
	expect_equal(z$m_prior_q95, q, tolerance = .025)
	expect_gt(z$m_prior_q95, qnorm(.975) + .5)
	expect_equal(z$m_post_mean, 2)
	# the reference does not depend on the prior of the initial state
	f$prior$beta0_var = 1000
	g = detect_change_point(f, n_prior_sims = 40000, seed = 44)
	expect_equal(g$m_prior_q95, z$m_prior_q95)
	expect_false('initial_mean' %in% names(formals(detect_change_point)))
})

test_that('rescaling a covariate preserves standardized RW1 change-point evidence', {
	f = transition_reference_fixture();g = f
	g$Xlist = lapply(g$Xlist, function(x) {x[, , 2] = x[, , 2] * 10; x})
	g$BETA[, 2, ] = g$BETA[, 2, ] / 10
	a = detect_change_point(f, n_prior_sims = 1000, seed = 7)
	b = detect_change_point(g, n_prior_sims = 1000, seed = 7)
	expect_equal(a, b, tolerance = 1e-10)
})

test_that('coefficient reference refuses unsupported or mismatched states', {
	f = transition_reference_fixture()
	for (kind in c('rw2', 'matern32')) {
		f$dynamic_beta_kind = kind
		expect_error(detect_change_point(f), class = 'lame_beta_transition')
	}
	f = transition_reference_fixture();f$time_index = c(1, 2, 4, 5, 7)
	expect_error(detect_change_point(f), 'equally spaced')
	f = transition_reference_fixture();f$SIGMA_BETA = f$SIGMA_BETA[-1, , drop = FALSE]
	expect_error(detect_change_point(f), 'aligned')
	f = transition_reference_fixture();f$SIGMA_BETA[1, 1] = 0
	expect_error(detect_change_point(f), 'positive-scaled')
	f = transition_reference_fixture();f$RHO_BETA[1, 1] = .8
	expect_error(detect_change_point(f), 'Persistence')
	f = transition_reference_fixture()
	expect_error(detect_change_point(f, coefs = 'intercept'), 'dynamic coefficients')
	expect_error(detect_change_point(f, n_prior_sims = 0), 'positive integer')
	expect_error(detect_change_point(f, n_prior_sims = 1.5), 'positive integer')
	f$BETA[1, 2, 1] = NA_real_
	expect_error(detect_change_point(f), 'finite paths')
})

test_that('symmetric collapsed node innovations retain the averaging variance', {
	f = transition_reference_fixture();f$symmetric = TRUE
	dimnames(f$BETA)[[2]][2] = 'x_node'
	f$Xlist = lapply(f$Xlist, function(x) {
		z = array(
			c(x, x[, , 2]),
			c(3, 3, 3),
			dimnames = list(letters[1:3], letters[1:3], c('intercept', 'x_row', 'x_col'))
		)
		z
	})
	s = lame:::.lame_beta_transition_scale(f)
	expect_equal(unname(s$sd), c(0, sqrt(2)))
})

test_that('saved-design innovation scales match the actual sampler geometry', {
	original = lame:::build_beta_state_scales;captured = NULL
	testthat::local_mocked_bindings(build_beta_state_scales = function(...) {
		ans = original(...);captured <<- ans[[1L]];ans
	}, .package = 'lame')
	for (mode in c('directed', 'bipartite', 'symmetric_node')) {
		set.seed(144);n = 5L;TT = 3L;rn = letters[1:n];cn = if (mode == 'bipartite') LETTERS[1:4] else rn
		Y = lapply(seq_len(TT), function(t) {
			y = matrix(rnorm(length(rn) * length(cn)), length(rn), dimnames = list(rn, cn))
			if (mode == 'symmetric_node') y = (y + t(y)) / 2
			if (mode != 'bipartite') diag(y) = NA_real_
			y
		})
		args = list(
			Y = Y,
			family = 'normal',
			R = 0,
			burn = 2,
			nscan = 8,
			odens = 2,
			dynamic_beta_kind = 'rw1',
			verbose = FALSE,
			gof = FALSE,
			seed = 49
		)
		if (mode == 'symmetric_node') {
			args$symmetric = TRUE;args$Xrow = rep(
				list(matrix(seq_len(n) / 3, n, dimnames = list(rn, 'x'))),
				TT
			);args$dynamic_beta = 'x_node'
		} else {
			args$mode = if (mode == 'bipartite') 'bipartite' else 'unipartite'
			args$Xdyad = lapply(
				Y,
				function(y) array(seq_along(y) / 7, c(dim(y), 1), dimnames = c(dimnames(y), list('x')))
			)
			args$dynamic_beta = 'dyad'
		}
		f = suppressWarnings(do.call(lame, args));actual = diag(captured)
		recovered = lame:::.lame_beta_transition_scale(f)$sd[f$beta_dynamic_mask]^2
		expected = if (mode == 'symmetric_node') sum(actual) / 4 else actual
		expect_equal(as.numeric(recovered), as.numeric(expected), tolerance = 1e-12, info = mode)
	}
})
