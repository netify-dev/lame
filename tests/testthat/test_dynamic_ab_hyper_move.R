# collapsed update of the dynamic additive-effect hyperparameters in every mode and family
# (R/dynamic_ab_hyper_move.R; the one-margin target is in R/binary_bipartite_scale_moves.R)
skip_on_cran()

# a stacked target and a symmetric/one-margin context for given totals
hyper_totals = function(n, N, seed, centred = TRUE) {
	set.seed(seed)
	a = matrix(rnorm(n * N), n); b = matrix(rnorm(n * N), n)
	if (centred) { a = sweep(a, 2, colMeans(a)); b = sweep(b, 2, colMeans(b)) }
	list(a = a, b = b)
}

test_that('the stacked density of the totals is the dense Gaussian density', {
	set.seed(3501)
	for (k in 1:40) {
		n = sample(2:7, 1); N = sample(2:6, 1)
		va = exp(rnorm(1)); vb = exp(rnorm(1)); r = runif(1, -.95, .95)
		Sab = matrix(c(va, r * sqrt(va * vb), r * sqrt(va * vb), vb), 2)
		Sig = .lame_bip_ab_cov(exp(rnorm(1)), runif(1, -.9, .9), N)
		tot = hyper_totals(n, N, 3500 + k, centred = k %% 2 == 0)
		df = if (k %% 2 == 0) n - 1L else n
		V2 = kronecker(Sab, matrix(1, N, N)) + kronecker(diag(2), Sig)
		X = cbind(tot$a, tot$b)
		dense = -0.5 * df * as.numeric(determinant(V2)$modulus) - 0.5 * sum(solve(V2) * crossprod(X))
		expect_lt(abs(.lame_uni_ab_loglik(Sab, Sig, .lame_uni_ab_cross(tot$a, tot$b), df) - dense), 1e-9)
	}
	# uncorrelated levels: V2 is block diagonal and the stacked target is the sum of the two one-margin targets
	tot = hyper_totals(6, 4, 3550)
	Sig = .lame_bip_ab_cov(.4, .6, 4)
	stacked = .lame_uni_ab_loglik(diag(c(1.3, .8)), Sig, .lame_uni_ab_cross(tot$a, tot$b), 5L)
	expect_lt(
		abs(
			stacked - (.lame_bip_ab_side_loglik(
				1.3,
				Sig,
				crossprod(tot$a),
				5L
			) + .lame_bip_ab_side_loglik(.8, Sig, crossprod(tot$b), 5L))
		),
		1e-10
	)
})

test_that('integrating the zero-sum levels out of the stacked prior gives the stated density of the totals', {
	# n = 3 actors, T = 2 periods, both margins centred: in an orthonormal basis of the zero-sum space, levels
	# l ~ N(0, Sab (x) I_{n-1}) and deviation rows N(0, Sigma); the density of the totals integrates the 2 (n - 1)
	# level coordinates by product Gauss-Hermite quadrature
	n = 3L; N = 2L; k = n - 1L
	Q = qr.Q(qr(cbind(1, diag(n)[, 1:k])))[, 2:(k + 1)]
	Sab = matrix(c(.7, .4 * sqrt(.7 * 1.3), .4 * sqrt(.7 * 1.3), 1.3), 2)
	Sig = .lame_bip_ab_cov(.6, .55, N); iS = solve(Sig); ldS = as.numeric(determinant(Sig)$modulus)
	tot = hyper_totals(n, N, 3551)
	xi = crossprod(Q, cbind(tot$a, tot$b))
	m = 13L; i = seq_len(m - 1); Jm = matrix(0, m, m); Jm[cbind(i, i + 1)] = sqrt(i); Jm[cbind(i + 1, i)] = sqrt(i)
	e = eigen(Jm, symmetric = TRUE); nodes = e$values; wts = e$vectors[1, ]^2
	grid = as.matrix(expand.grid(rep(list(seq_len(m)), 2L * k)))
	lv = matrix(nodes[grid], nrow(grid)) %*% t(t(chol(kronecker(Sab, diag(k)))))
	w = apply(matrix(wts[grid], nrow(grid)), 1, prod)
	lp = numeric(nrow(grid))
	for (mg in 1:2) for (r in 1:k) {
		dl = matrix(xi[r, (mg - 1) * N + 1:N], nrow(grid), N, byrow = TRUE) - lv[, (mg - 1) * k + r]
		lp = lp - 0.5 * (N * log(2 * pi) + ldS) - 0.5 * rowSums((dl %*% iS) * dl)
	}
	quad = log(sum(w * exp(lp)))
	V2 = kronecker(Sab, matrix(1, N, N)) + kronecker(diag(2), Sig)
	closed = -0.5 * k * (2 * N * log(2 * pi) + as.numeric(determinant(V2)$modulus)) - 0.5 * sum(solve(V2) * crossprod(xi))
	expect_lt(abs(quad - closed), 2e-4)
	# the sampler's target is that density up to a dimension constant
	fast = .lame_uni_ab_loglik(Sab, Sig, .lame_uni_ab_cross(tot$a, tot$b), k)
	expect_lt(abs((fast - closed) - k * N * log(2 * pi)), 1e-9)
})

test_that('the update returns finite, bounded draws and a state it can resume from in every mode', {
	set.seed(3502)
	prior = list(
		etaab = 5,
		Sab0 = matrix(c(.8, .2, .2, 1.1), 2),
		sigma_ab_shape = 2,
		sigma_ab_scale = .3,
		rho_ab_mean = .8,
		rho_ab_sd = .15
	)
	n = 7L; N = 4L
	tot = hyper_totals(n, N, 3503)
	for (mode in c('directed', 'symmetric', 'row', 'col', 'bip')) {
		rv = mode != 'col'; cv = mode %in% c('directed', 'col', 'bip')
		cx = .lame_ab_hyper_context(
			tot$a,
			if (mode == 'bip') tot$b[1:5, ] else tot$b,
			rv,
			cv,
			mode == 'symmetric',
			mode == 'bip',
			TRUE,
			prior,
			1.5
		)
		expect_equal(cx$stacked, mode == 'directed')
		Sab = matrix(c(1, if (mode == 'directed') .2 else 0, if (mode == 'directed') .2 else 0, 1.2), 2)
		state = NULL
		for (v in 1:6) {
			o = .lame_ab_hyper_step(
				tot$a,
				if (mode == 'bip') tot$b[1:5, ] else tot$b,
				Sab,
				.4,
				.7,
				cx,
				state,
				tries = 5L,
				adapt = v <= 4L,
				stat_var_cap = 1.5
			)
			expect_true(all(is.finite(o$Sab)) && o$sigma > 0 && abs(o$rho) < .99)
			expect_lte(o$sigma^2 / (1 - o$rho^2), 1.5 + 1e-12)
			expect_true(isSymmetric(o$Sab))
			if (mode == 'symmetric') expect_equal(o$Sab[1, 1], o$Sab[2, 2])
			if (mode != 'directed') expect_equal(o$Sab[1, 2], 0)
			if (mode == 'row') expect_equal(o$Sab[2, 2], Sab[2, 2])
			if (mode == 'col') expect_equal(o$Sab[1, 1], Sab[1, 1])
			Sab = o$Sab; state = o$state
		}
		expect_equal(unname(state[3]), 6 * 5)
	}
	# a state of another length, such as a checkpoint of another configuration, starts a fresh proposal
	cx = .lame_ab_hyper_context(tot$a, tot$b, TRUE, TRUE, FALSE, FALSE, TRUE, prior)
	o = .lame_ab_hyper_step(tot$a, tot$b, diag(2), .4, .7, cx, state = c(1, 2, 3), tries = 2L)
	expect_equal(unname(o$state[3]), 2)
})

test_that('the stacked update leaves the prior invariant for directed, symmetric and one-margin unipartite networks', {
	skip_on_cran()
	# draw (hyperparameters, levels, deviation paths) from the prior, update given the totals and redraw the levels:
	# the law of the state must not change (Geweke 2004)
	prior = list(
		etaab = 7,
		Sab0 = matrix(c(.6, .2, .2, .9), 2),
		sigma_ab_shape = 2,
		sigma_ab_scale = .3,
		rho_ab_mean = .8,
		rho_ab_sd = .15
	)
	run = function(nrep, mode, n = 5L, N = 3L, centred = TRUE, tries = 3L, seed = 1L, cx_edit = NULL) {
		set.seed(seed)
		rv = mode != 'col'; cv = mode %in% c('directed', 'col'); sy = mode == 'symmetric'
		ar1 = function(s2, rho) { m = matrix(0, n, N); m[, 1] = rnorm(n, 0, sqrt(s2 / (1 - rho^2)))
			for (t in seq_len(N)[-1]) m[, t] = rho * m[, t - 1] + rnorm(n, 0, sqrt(s2)); if (centred) sweep(m, 2, colMeans(m)) else m }
		cx = .lame_ab_hyper_context(matrix(0, n, N), matrix(0, n, N), rv, cv, sy, FALSE, centred, prior)
		if (!is.null(cx_edit)) cx = cx_edit(cx)
		draw = function() {
			repeat { s2 = 1 / rgamma(
				1,
				prior$sigma_ab_shape,
				prior$sigma_ab_scale
			); rho = rnorm(1, prior$rho_ab_mean, prior$rho_ab_sd); if (abs(rho) < .99) break }
			if (mode == 'directed') Sab = solve(rWishart(1, prior$etaab, solve(prior$etaab * prior$Sab0))[, , 1])
			else { v = 1 / rgamma(
				1,
				prior$etaab / 2,
				prior$etaab * prior$Sab0[if (mode == 'col') 2 else 1, if (mode == 'col') 2 else 1] / 2
			)
				Sab = switch(mode, row = diag(c(v, 1e-8)), col = diag(c(1e-8, v)), symmetric = diag(c(v, v))) }
			if (mode == 'directed') { L = matrix(rnorm(2 * n), n) %*% chol(Sab); la = L[, 1]; lb = L[, 2] }
			else { la = if (rv) rnorm(n, 0, sqrt(Sab[1, 1])) else numeric(n); lb = if (sy) la else if (cv) rnorm(n, 0, sqrt(Sab[2, 2])) else numeric(n) }
			if (centred) { la = la - mean(la); lb = if (sy) la else lb - mean(lb) }
			da = if (rv) ar1(s2, rho) else matrix(0, n, N); db = if (sy) da else if (cv) ar1(s2, rho) else matrix(0, n, N)
			list(Sab = Sab, sigma = sqrt(s2), rho = rho, la = la, lb = lb, da = da, db = db)
		}
		stat = function(s) {
			a = s$la + s$da; b = s$lb + s$db
			va = s$Sab[1, 1]; vb = s$Sab[2, 2]; r = s$Sab[1, 2] / sqrt(va * vb)
			c(
				log(va),
				log(vb),
				r,
				r^2,
				log(s$sigma),
				s$rho,
				s$rho^2,
				log(va) * s$rho,
				r * s$rho,
				r * log(va / vb),
				s$la[1]^2,
				s$lb[1]^2,
				s$la[1] * s$lb[1],
				s$la[1] * s$da[1, 1],
				s$da[1, 1]^2,
				s$db[2, N]^2,
				s$da[2, 1] * s$da[2, N],
				sum(s$la^2),
				sum(s$lb^2),
				sum(s$la * s$lb),
				sum(s$da^2),
				sum(s$db^2),
				sum(a[, 1] * a[, N]),
				sum(b[, 1] * b[, N]),
				sum(a[, 1] * b[, N]),
				a[1, 1] * b[1, 1],
				sum(a * b),
				1 / va,
				1 / vb,
				1 / s$sigma^2
			)
		}
		old = new = vector('list', nrep)
		for (r in seq_len(nrep)) {
			s = draw(); old[[r]] = stat(s)
			a = s$la + s$da; b = s$lb + s$db
			o = .lame_ab_hyper_step(a, b, s$Sab, s$sigma, s$rho, cx, NULL, tries)
			s$Sab = o$Sab; s$sigma = o$sigma; s$rho = o$rho
			lev = .lame_ab_level_draw(
				a,
				b,
				s$rho,
				s$sigma,
				s$Sab,
				symmetric = sy,
				center = centred,
				rvar = rv,
				cvar = if (sy) rv else cv
			)
			s$la = lev$a; s$lb = lev$b; s$da = a - s$la; s$db = b - s$lb
			new[[r]] = stat(s)
		}
		d = do.call(rbind, new) - do.call(rbind, old); keep = apply(d, 2, sd) > 0
		max(abs(colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrep))))
	}
	expect_lt(run(3000L, 'directed', seed = 31L), 4.5)
	expect_lt(run(3000L, 'directed', centred = FALSE, seed = 32L), 4.5)
	expect_lt(run(3000L, 'symmetric', seed = 33L), 4.5)
	expect_lt(run(3000L, 'row', seed = 34L), 4.5)
	expect_lt(run(3000L, 'col', seed = 35L), 4.5)
	expect_lt(run(3000L, 'directed', n = 3L, N = 2L, seed = 36L), 4.5)
	# power: a wrong prior scale or a wrong number of rows is detected
	expect_gt(run(6000L, 'directed', seed = 37L, cx_edit = function(cx) { cx$S0 = cx$S0 * 1.5; cx }), 4.5)
	expect_gt(run(6000L, 'directed', seed = 38L, cx_edit = function(cx) { cx$centred = FALSE; cx }), 4.5)
	expect_gt(run(6000L, 'symmetric', seed = 39L, cx_edit = function(cx) { cx$rho_mean = cx$rho_mean - .15; cx }), 4.5)
})

uni_panel = function(family, seed, n = 8L, nt = 3L, sym = FALSE) {
	set.seed(seed)
	nm = sprintf('a%02d', 1:n)
	cell = function() {
		m = switch(
			family,
			binary = matrix(rbinom(n * n, 1, .2), n),
			normal = matrix(rnorm(n * n), n),
			poisson = matrix(rpois(n * n, .5), n),
			ordinal = matrix(sample(1:4, n * n, TRUE, c(.5, .25, .15, .1)), n),
			cbin = { y = matrix(0, n, n); for (i in 1:n) y[i, sample(n, sample(0:3, 1))] = 1; y },
			frn = { y = matrix(0, n, n); for (i in 1:n) { k = sample(0:3, 1); if (k > 0) y[i, sample(n, k)] = seq_len(k) }; y }
		)
		if (sym) m[lower.tri(m)] = t(m)[lower.tri(m)]
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	}
	Y = lapply(1:nt, function(t) cell())
	names(Y) = paste0('p', 1:nt)
	Y
}
uni_run = function(Y, family, prior = list(), sym = FALSE, ...) suppressMessages(lame(
	Y,
	mode = 'unipartite',
	family = family,
	symmetric = sym,
	R = 0L,
	odmax = if (family %in% c('cbin', 'frn')) 3L else NULL,
	dynamic_ab = TRUE,
	burn = 15,
	nscan = 40,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = prior,
	...
))

test_that('lame() runs the update for unipartite directed and symmetric dynamic_ab fits of every family and records it', {
	skip_on_cran()
	for (fam in c('binary', 'normal', 'poisson', 'ordinal', 'cbin', 'frn')) {
		Y = uni_panel(fam, 3600L + nchar(fam))
		on = uni_run(Y, fam)
		expect_equal(unname(on$ab_hyper_move['attempts']), 55 * 8)
		expect_gt(unname(on$ab_hyper_move['accepts']), 0)
		expect_true(
			all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$sigma_ab)) && all(abs(on$rho_ab) < 1)
		)
		off = uni_run(Y, fam, list(ab_hyper_move = FALSE))
		expect_null(off$ab_hyper_move)
		expect_false(isTRUE(all.equal(on$sigma_ab, off$sigma_ab)))
	}
	for (fam in c('binary', 'normal', 'poisson')) {
		Ys = uni_panel(fam, 3700L + nchar(fam), sym = TRUE)
		on = uni_run(Ys, fam, sym = TRUE)
		expect_equal(unname(on$ab_hyper_move['attempts']), 55 * 8)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$sigma_ab)))
		expect_null(uni_run(Ys, fam, list(ab_hyper_move = FALSE), sym = TRUE)$ab_hyper_move)
	}
	# one margin and no margin
	Y = uni_panel('binary', 3800L)
	expect_equal(unname(uni_run(Y, 'binary', rvar = TRUE, cvar = FALSE)$ab_hyper_move['attempts']), 55 * 8)
	expect_equal(unname(uni_run(Y, 'binary', rvar = FALSE, cvar = TRUE)$ab_hyper_move['attempts']), 55 * 8)
	expect_null(uni_run(Y, 'binary', rvar = FALSE, cvar = FALSE)$ab_hyper_move)
	# static effects have no such hyperparameters
	expect_null(
		suppressMessages(lame(
			Y,
			mode = 'unipartite',
			family = 'binary',
			R = 0L,
			dynamic_ab = FALSE,
			burn = 10,
			nscan = 20,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))$ab_hyper_move
	)
})

test_that('lame() runs the update for bipartite dynamic_ab fits of every family', {
	skip_on_cran()
	nA = 9L; nB = 11L; nt = 3L
	for (fam in c('normal', 'poisson', 'ordinal', 'cbin', 'frn')) {
		set.seed(3900L + nchar(fam))
		rn = sprintf('r%02d', 1:nA); cn = sprintf('c%02d', 1:nB)
		Y = lapply(1:nt, function(t) {
			y = switch(
				fam,
				normal = matrix(rnorm(nA * nB), nA),
				poisson = matrix(rpois(nA * nB, .5), nA),
				ordinal = matrix(sample(1:4, nA * nB, TRUE, c(.5, .25, .15, .1)), nA),
				cbin = { y = matrix(0, nA, nB); for (i in 1:nA) y[i, sample(nB, sample(0:3, 1))] = 1; y },
				frn = { y = matrix(0, nA, nB); for (i in 1:nA) { k = sample(0:3, 1); if (k > 0) y[i, sample(nB, k)] = seq_len(k) }; y }
			)
			dimnames(y) = list(rn, cn); y })
		names(Y) = paste0('p', 1:nt)
		run = function(prior) suppressMessages(lame(
			Y,
			mode = 'bipartite',
			family = fam,
			R_row = 0L,
			R_col = 0L,
			odmax = if (fam %in% c('cbin', 'frn')) 3L else NULL,
			dynamic_ab = TRUE,
			burn = 15,
			nscan = 40,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior
		))
		on = run(list())
		expect_equal(unname(on$ab_hyper_move['attempts']), 55 * 8)
		expect_gt(unname(on$ab_hyper_move['accepts']), 0)
		expect_true(
			all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$sigma_ab)) && all(abs(on$rho_ab) < 1)
		)
		expect_null(run(list(ab_hyper_move = FALSE))$ab_hyper_move)
	}
})

test_that('the proposal of the unipartite update survives a checkpoint', {
	skip_on_cran()
	Y = uni_panel('binary', 3950L)
	a = list(
		Y = Y,
		mode = 'unipartite',
		family = 'binary',
		R = 0L,
		dynamic_ab = TRUE,
		burn = 15L,
		nscan = 25L,
		odens = 5L,
		seed = 3L,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	state = f$start_vals$checkpoint_state$values
	expect_true(is.numeric(state$ab_hyper_par))
	expect_equal(length(state$ab_hyper_par), 4L + 5L + 25L)
	expect_equal(unname(state$ab_hyper_par[3]), 40 * 8)
	a$start_vals = f$start_vals; a$burn = 0L; a$nscan = 10L
	g = suppressWarnings(suppressMessages(do.call(lame, a)))
	expect_equal(unname(g$ab_hyper_move['attempts']), 50 * 8)
	expect_equal(g$start_vals$checkpoint_state$values$ab_hyper_par[-(1:4)], state$ab_hyper_par[-(1:4)], tolerance = 0)
})

test_that('the update runs with changing coefficients and with the dynamic correlation of the dyads', {
	skip_on_cran()
	set.seed(3960)
	n = 8L; nt = 3L; nm = sprintf('a%02d', 1:n)
	Y = lapply(
		1:nt,
		function(t) { y = matrix(rbinom(n * n, 1, .25), n); diag(y) = NA; dimnames(y) = list(nm, nm); y }
	); names(Y) = paste0('p', 1:nt)
	X = lapply(1:nt, function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, 'x'))); names(X) = names(Y)
	on = suppressMessages(lame(
		Y,
		Xdyad = X,
		mode = 'unipartite',
		family = 'binary',
		R = 0L,
		dynamic_ab = TRUE,
		dynamic_beta = 'x_dyad',
		burn = 15,
		nscan = 40,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(unname(on$ab_hyper_move['attempts']), 55 * 8)
	expect_true(all(is.finite(on$sigma_ab)) && all(abs(on$rho_ab) < 1))
	on = suppressMessages(lame(
		Y,
		Xdyad = X,
		mode = 'unipartite',
		family = 'binary',
		R = 1L,
		dcor = TRUE,
		dynamic_ab = TRUE,
		burn = 15,
		nscan = 40,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(unname(on$ab_hyper_move['attempts']), 55 * 8)
	expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
})
