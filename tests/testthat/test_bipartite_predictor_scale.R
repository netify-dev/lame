# joint scale move of the predictor of a binary bipartite fit

bpls_panel = function(nA = 12L, nB = 10L, N = 3L, seed = 7L, fz = .3) {
	set.seed(seed)
	rn = sprintf('r%02d', seq_len(nA)); cn = sprintf('c%02d', seq_len(nB))
	nz = round(fz * nA); nzb = round(fz * nB)
	a0 = rnorm(nA, 0, .8); b0 = rnorm(nB, 0, .8); a0[seq_len(nz)] = -3.5; b0[seq_len(nzb)] = -3.5
	u = rnorm(nA); v = rnorm(nB)
	Y = X = vector('list', N)
	for (t in seq_len(N)) {
		x = matrix(rnorm(nA * nB), nA, nB)
		z = -1 + .5 * x + outer(a0, b0, '+') + .6 * outer(u, v) + matrix(rnorm(nA * nB), nA, nB)
		y = matrix(as.numeric(z > 0), nA, nB)
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(nA, nB, 1L), dimnames = list(rn, cn, 'x'))
	}
	list(Y = Y, X = X)
}
bpls_fit = function(p, ..., prior = list(), burn = 30, nscan = 120, odens = 2) {
	suppressWarnings(suppressMessages(lame(
		p$Y,
		Xdyad = p$X,
		family = 'binary',
		mode = 'bipartite',
		R_row = 1L,
		R_col = 1L,
		dynamic_uv = TRUE,
		burn = burn,
		nscan = nscan,
		odens = odens,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

test_that('the context lists the observed cells with the sign of each', {
	Y = array(c(1, 0, NA, 1, 0, 0, 1, NA), c(2, 2, 2))
	cx = .lame_bpl_context(Y)
	expect_equal(cx$n1, 3L); expect_equal(cx$n0, 3L)
	expect_equal(cx$sY, c(1, 1, 1, -1, -1, -1))
	expect_equal(Y[cx$idx], c(1, 1, 1, 0, 0, 0))
})

test_that('the acceptance ratio is the posterior ratio of the scaled blocks times the Jacobian', {
	set.seed(2)
	p = 3L; g = 2.5; nA = 5L; nB = 6L
	beta = rnorm(p); a = rnorm(nA); a = a - mean(a); b = rnorm(nB); b = b - mean(b); G = matrix(rnorm(2), 1)
	va = 1.3; vb = .7
	ez = rnorm(40); sY = sample(c(-1, 1), 40, TRUE)
	iV0 = diag(p) / g
	# one proposal: the factor is the first draw of the log-normal step
	set.seed(11); cc = exp(rnorm(1, 0, .3))
	set.seed(11)
	o = .lame_bpl_scale(
		ez,
		sY,
		beta,
		iV0,
		sum(a^2) / va + sum(b^2) / vb,
		nA + nB - 2L,
		0,
		0L,
		sum(G^2),
		length(G),
		tries = 1L,
		sd = .3
	)
	# the log ratio by hand: probit likelihood, Gaussian block priors on their subspaces, Jacobian of the scaling
	lr = sum(pnorm(sY * cc * ez, log.p = TRUE)) - sum(pnorm(sY * ez, log.p = TRUE)) +
		sum(dnorm(cc * beta, 0, sqrt(g), log = TRUE)) - sum(dnorm(beta, 0, sqrt(g), log = TRUE)) + p * log(cc) +
		(-(cc^2 - 1) * sum(a^2) / (2 * va)) + (-(cc^2 - 1) * sum(b^2) / (2 * vb)) + (nA - 1 + nB - 1) * log(cc) +
		sum(dnorm(cc * G, log = TRUE)) - sum(dnorm(G, log = TRUE)) + length(G) * log(cc)
	set.seed(11); rnorm(1, 0, .3); u = log(runif(1))
	expect_equal(o$c != 1, u < lr)
	expect_equal(o$att, 1L)
	expect_equal(o$acc, as.integer(o$c != 1))
})

test_that('the visit respects the bound on the interaction and keeps the counters', {
	set.seed(3)
	ez = rnorm(30); sY = sample(c(-1, 1), 30, TRUE)
	# a factor beyond the room has prior zero: with a room of one the interaction can only shrink
	for (i in 1:20) {
		o = .lame_bpl_scale(ez, sY, c(.1, .2), diag(2), 1, 8L, 0, 0L, 4, 1L, room = 1, tries = 6L, sd = .5)
		expect_lte(o$c, 1); expect_equal(o$att, 6L)
	}
	o = .lame_bpl_scale(ez, sY, c(.1, .2), diag(2), 1, 8L, 0, 0L, 4, 1L, room = 0, tries = 6L, sd = .5)
	expect_equal(o$c, 1); expect_equal(o$acc, 0L)
	o = .lame_bpl_scale(ez, sY, c(.1, .2), diag(2), 1, 8L, 0, 0L, 4, 1L, room = Inf, tries = 40L, sd = .3)
	expect_equal(o$att, 40L); expect_gt(o$acc, 0L); expect_true(is.finite(o$c) && o$c > 0)
	# step adaptation moves the step towards the target acceptance and only when asked
	o1 = .lame_bpl_scale(ez, sY, c(.1, .2), diag(2), 1, 8L, 0, 0L, 4, 1L, tries = 30L, sd = .3, adapt = FALSE)
	o2 = .lame_bpl_scale(ez, sY, c(.1, .2), diag(2), 1, 8L, 0, 0L, 4, 1L, tries = 30L, sd = .3, adapt = TRUE)
	expect_equal(o1$sd, .3); expect_false(o2$sd == .3)
})

test_that('a binary bipartite fit makes the scale move and records its acceptance', {
	skip_on_cran()
	p = bpls_panel()
	f = bpls_fit(p, prior = list(level_move_every = 2), odens = 3)
	expect_true(!is.null(f$level_move))
	expect_equal(unname(f$level_move['attempts']), 8 * 75)
	expect_gt(unname(f$level_move['rate']), 0); expect_lt(unname(f$level_move['rate']), 1)
	expect_true(all(is.finite(f$BETA)) && all(is.finite(f$G_samples)))
	# the adapted step is kept for a continued run, which does not adapt it again
	expect_true(is.finite(f$start_vals$level_move_sd))
	g = bpls_fit(p, prior = list(level_move_every = 2), burn = 0, nscan = 20, start_vals = f$start_vals)
	expect_equal(unname(g$level_move['step_sd']), unname(f$level_move['step_sd']))
	# other families leave it out
	dyn = bpls_fit(p, dynamic_ab = TRUE, prior = list(level_move_every = 2))
	expect_gt(unname(dyn$level_move['attempts']), 0)
	nog = suppressWarnings(suppressMessages(lame(
		p$Y,
		Xdyad = p$X,
		family = 'binary',
		mode = 'bipartite',
		R_row = 0L,
		R_col = 0L,
		burn = 20,
		nscan = 60,
		odens = 2,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = list(level_move_every = 2)
	)))
	expect_gt(unname(nog$level_move['attempts']), 0)
	pois = suppressWarnings(suppressMessages(lame(
		lapply(p$Y, function(y) y * 2),
		Xdyad = p$X,
		family = 'poisson',
		mode = 'bipartite',
		R_row = 0L,
		R_col = 0L,
		burn = 20,
		nscan = 60,
		odens = 2,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	expect_null(pois$level_move)
})

test_that('the visits happen on the schedule of level_move_every', {
	p = bpls_panel()
	f = bpls_fit(p, burn = 0, nscan = 60, prior = list(level_move_every = 10, level_move_tries = 3))
	expect_equal(unname(f$level_move['attempts']), 3 * 6)
})

test_that('the move leaves the exact posterior of a tiny panel invariant', {
	skip_on_cran()
	cases = list(
		list(variant = 'static', ndraw = 3e5),
		list(variant = 'dyn', ndraw = 3e5),
		list(variant = 'nog', ndraw = 2e5)
	)
	for (k in seq_along(cases)) {
		r = do.call(bplt_exact, c(cases[[k]], list(seed = 300L + k)))
		expect_gt(r$ess, 300)
		expect_lt(r$maxz, 4.6)
	}
	# with a bound on the interaction (a proposal beyond it has prior zero)
	r = bplt_exact('static', ndraw = 3e5, seed = 310L, gmax = 1.2)
	expect_lt(r$maxz, 4.6)
})

test_that('the exact-posterior check detects a move that mis-states the prior of a block', {
	skip_on_cran()
	prior = .lame_scale_prior_logratio
	local_mocked_bindings(
		.lame_scale_prior_logratio = function(cc, beta, iV0, q_lev, n_lev, q_dev = 0, n_dev = 0L, q_uv = 0, n_uv = 0L)
			prior(cc, beta, iV0, q_lev, n_lev + 2L, q_dev, n_dev, q_uv, n_uv),
		.package = 'lame'
	)
	expect_gt(bplt_exact('static', ndraw = 3e5, seed = 321L)$maxz, 6)
})
