# the moves along the ridges of the bipartite sampler (R/binary_bipartite_scale_moves.R)
skip_on_cran()

test_that('the generalised inverse Gaussian draw matches its density', {
	skip_on_cran()
	set.seed(3301)
	for (p in list(c(3.5, 2, 6), c(0.4, 1.5, 0.8))) {
		lambda = p[1]; chi = p[2]; psi = p[3]
		x = vapply(seq_len(3500), function(i) .lame_rgig(lambda, chi, psi), numeric(1))
		u = seq(-14, 14, length.out = 40001)
		lf = lambda * u - 0.5 * (chi * exp(-u) + psi * exp(u)); w = exp(lf - max(lf)); w = w / sum(w)
		m1 = sum(w * exp(u)); m2 = sum(w * exp(2 * u))
		expect_equal(mean(x), m1, tolerance = 0.07)
		expect_equal(sd(x), sqrt(m2 - m1^2), tolerance = 0.1)
		expect_lt(abs(mean(log(x)) - sum(w * u)), 0.05)
	}
	# a lower bound truncates by rejection, a degenerate density has no draw
	x = vapply(seq_len(300), function(i) .lame_rgig(2, 1, 1, lower = 1.5), numeric(1))
	expect_true(all(x >= 1.5))
	expect_true(is.na(.lame_rgig(2, 0, 1)))
	expect_true(is.na(.lame_rgig(2, 1, -1)))
})

test_that('the scale ridge step keeps the product and the bound and draws every coordinate', {
	set.seed(3302)
	for (cf in list(c(8, 12, 1, 1), c(6, 9, 2, 2), c(5, 7, 2, 1), c(10, 6, 1, 3))) {
		U = matrix(rnorm(cf[1] * cf[3]), cf[1]); V = matrix(rnorm(cf[2] * cf[4]), cf[2]); G = matrix(rnorm(cf[3] * cf[4]), cf[3])
		o = .lame_bip_uv_g_scale_step(U, V, G)
		expect_equal(o$U %*% o$G %*% t(o$V), U %*% G %*% t(V), tolerance = 1e-10)
		expect_equal(o$drawn, cf[3] + cf[4])
		expect_false(isTRUE(all.equal(o$G, G)))
	}
	# the prior bound on G holds for every draw
	U = matrix(rnorm(14), 7); V = matrix(rnorm(18), 9); G = matrix(c(.9, -.6, .4, .8), 2)
	for (k in 1:200) { o = .lame_bip_uv_g_scale_step(U, V, G, 1.2); expect_lte(max(abs(o$G)), 1.2 + 1e-12); U = o$U; V = o$V; G = o$G }
})

test_that('the scale ridge step leaves the prior of the positions and the interaction invariant', {
	skip_on_cran()
	# the likelihood is constant along the orbit, so the step is valid for every posterior iff it leaves the prior invariant
	run = function(nrep, nA, nB, RA, RB, g_max = Inf, seed = 1L) {
		set.seed(seed)
		stat = function(U, V, G) {
			M = U %*% G %*% t(V)
			c(
				colSums(U^2),
				colSums(V^2),
				rowSums(G^2),
				colSums(G^2),
				sum(G^2),
				G[1, 1],
				G[1, 1]^2,
				U[1, 1]^2,
				V[nB, RB]^2,
				U[1, 1] * G[1, 1],
				log(sum(U^2)),
				log(sum(G^2)),
				log(sum(V^2)),
				M[1, 1],
				M[2, 3]
			)
		}
		old = new = NULL
		for (r in seq_len(nrep)) {
			repeat { G = matrix(rnorm(RA * RB), RA, RB); if (max(abs(G)) <= g_max) break }
			U = matrix(rnorm(nA * RA), nA); V = matrix(rnorm(nB * RB), nB)
			o = .lame_bip_uv_g_scale_step(U, V, G, g_max)
			old = rbind(old, stat(U, V, G)); new = rbind(new, stat(o$U, o$V, o$G))
		}
		d = new - old; keep = apply(d, 2, sd) > 0
		max(abs(colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrep))))
	}
	expect_lt(run(1000L, 8L, 12L, 1L, 1L, seed = 11L), 4.5)
	expect_lt(run(1000L, 6L, 9L, 2L, 2L, g_max = 1.8, seed = 12L), 4.5)
})

test_that('the density of the totals with the levels integrated out is the stated one', {
	# n = 3 actors, N = 2 periods, centred: the integral over the zero-sum levels of the level prior times the
	# AR(1) prior of the deviations equals |V|^(-(n - 1) / 2) exp(-tr(V^-1 W) / 2) up to a constant in the hyperparameters
	set.seed(3303)
	n = 3L; N = 2L
	a = matrix(rnorm(n * N, 0, 1.2), n); a = sweep(a, 2, colMeans(a))
	Hb = qr.Q(qr(contr.helmert(n)))
	grid = seq(-9, 9, length.out = 361); dz = diff(grid[1:2])
	z = as.matrix(expand.grid(grid, grid))
	log_integral = function(va, s2, rho) {
		Sig = .lame_bip_ab_cov(s2, rho, N); Si = solve(Sig)
		al = z %*% t(Hb) # levels on the zero-sum plane, one row per grid point
		q = rowSums(al^2) / va
		for (i in 1:n) { r = sweep(
			matrix(a[i, ], nrow(al), N, byrow = TRUE),
			1,
			al[, i]
		); q = q + rowSums((r %*% Si) * r) }
		-(n - 1) / 2 * log(va) - (n - 1) / 2 * as.numeric(determinant(Sig)$modulus) + log(sum(exp(-0.5 * q)) * dz^2)
	}
	closed = function(va, s2, rho) { Sig = .lame_bip_ab_cov(s2, rho, N); .lame_bip_ab_side_loglik(va, Sig, crossprod(a), n - 1L) }
	p1 = c(1.3, .09, .7); p2 = c(.6, .25, .2); p3 = c(2.4, .04, .9)
	expect_lt(
		abs(
			(log_integral(p1[1], p1[2], p1[3]) - log_integral(p2[1], p2[2], p2[3])) -
				(closed(p1[1], p1[2], p1[3]) - closed(p2[1], p2[2], p2[3]))
		),
		1e-6
	)
	expect_lt(
		abs(
			(log_integral(p3[1], p3[2], p3[3]) - log_integral(p2[1], p2[2], p2[3])) -
				(closed(p3[1], p3[2], p3[3]) - closed(p2[1], p2[2], p2[3]))
		),
		1e-6
	)
})

test_that('the collapsed update of the additive-effect hyperparameters leaves the prior invariant', {
	skip_on_cran()
	# draw (theta, levels, deviation paths) from the prior, apply the update given the totals and redraw the levels
	# from their conditional: the law of the state must not change (Geweke 2004)
	prior = list(
		etaab = 7,
		Sab0 = diag(c(0.6, 0.9)),
		sigma_ab_shape = 2,
		sigma_ab_scale = 0.3,
		rho_ab_mean = 0.8,
		rho_ab_sd = 0.15
	)
	run = function(nrep, nA = 4L, nB = 5L, N = 3L, centred = TRUE, rvar = TRUE, cvar = TRUE, tries = 3L, seed = 1L) {
		set.seed(seed)
		ar1 = function(n, s2, rho) { m = matrix(0, n, N); m[, 1] = rnorm(n, 0, sqrt(s2 / (1 - rho^2)))
			for (t in seq_len(N)[-1]) m[, t] = rho * m[, t - 1] + rnorm(n, 0, sqrt(s2)); if (centred) sweep(m, 2, colMeans(m)) else m }
		cx = .lame_bip_ab_context(matrix(0, nA, N), matrix(0, nB, N), rvar, cvar, centred, prior)
		draw = function() {
			repeat { va = 1 / rgamma(
				1,
				prior$etaab / 2,
				prior$etaab * prior$Sab0[1, 1] / 2
			); vb = 1 / rgamma(1, prior$etaab / 2, prior$etaab * prior$Sab0[2, 2] / 2)
				s2 = 1 / rgamma(
					1,
					prior$sigma_ab_shape,
					prior$sigma_ab_scale
				); rho = rnorm(1, prior$rho_ab_mean, prior$rho_ab_sd); if (abs(rho) < .99) break }
			la = if (rvar) rnorm(nA, 0, sqrt(va)) else numeric(nA); lb = if (cvar) rnorm(nB, 0, sqrt(vb)) else numeric(nB)
			if (centred) { la = la - mean(la); lb = lb - mean(lb) }
			list(
				va = va,
				vb = vb,
				sigma = sqrt(s2),
				rho = rho,
				la = la,
				lb = lb,
				da = if (rvar) ar1(nA, s2, rho) else matrix(0, nA, N),
				db = if (cvar) ar1(nB, s2, rho) else matrix(0, nB, N)
			)
		}
		stat = function(s) {
			a = s$la + s$da; b = s$lb + s$db
			c(
				log(s$va),
				log(s$vb),
				log(s$sigma),
				s$rho,
				s$rho^2,
				log(s$va) * s$rho,
				s$la[1]^2,
				s$lb[1]^2,
				s$la[1] * s$da[1, 1],
				s$da[1, 1]^2,
				s$db[2, N]^2,
				s$da[2, 1] * s$da[2, N],
				sum(s$la^2),
				sum(s$lb^2),
				sum(s$da^2),
				sum(s$db^2),
				sum(a[, 1] * a[, N]),
				a[1, 1] * a[1, N],
				1 / s$va,
				1 / s$sigma^2
			)
		}
		old = new = NULL
		for (r in seq_len(nrep)) {
			s = draw(); old = rbind(old, stat(s))
			a = s$la + s$da; b = s$lb + s$db
			o = .lame_bip_ab_collapsed_step(a, b, s$va, s$vb, s$sigma, s$rho, cx, NULL, tries = tries)
			s$va = o$va; s$vb = o$vb; s$sigma = o$sigma; s$rho = o$rho
			lev = .lame_ab_level_draw(
				a,
				b,
				s$rho,
				s$sigma,
				diag(c(s$va, s$vb)),
				symmetric = FALSE,
				center = centred,
				rvar = rvar,
				cvar = cvar
			)
			s$la = lev$a; s$lb = lev$b; s$da = a - s$la; s$db = b - s$lb
			new = rbind(new, stat(s))
		}
		d = new - old; keep = apply(d, 2, sd) > 0
		max(abs(colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrep))))
	}
	expect_lt(run(2500L, seed = 21L), 4.5)
	expect_lt(run(2500L, centred = FALSE, seed = 22L), 4.5)
	expect_lt(run(2500L, rvar = TRUE, cvar = FALSE, seed = 23L), 4.5)
	expect_lt(run(2500L, nA = 6L, nB = 3L, N = 2L, seed = 24L), 4.5)
})

test_that('ame() draws along the scale ridge of the positions as a prior switch and records it', {
	set.seed(7401); nA = 9L; nB = 11L
	y = matrix(
		rbinom(nA * nB, 1, .15),
		nA,
		nB,
		dimnames = list(sprintf('r%02d', 1:nA), sprintf('c%02d', 1:nB))
	); y[1:2, ] = 0; y[, 1:3] = 0
	x = array(rnorm(nA * nB), c(nA, nB, 1), dimnames = list(rownames(y), colnames(y), 'x'))
	run = function(prior, R, ...) suppressMessages(ame(
		y,
		Xdyad = x,
		mode = 'bipartite',
		family = 'binary',
		R_row = R,
		R_col = R,
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = prior,
		...
	))
	expect_no_warning({ on = run(list(), 1L) })
	expect_equal(unname(on$uv_g_scale_move['visits']), 80); expect_equal(unname(on$uv_g_scale_move['drawn']), 160)
	expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
	off = run(list(uv_g_scale_move = FALSE), 1L)
	expect_null(off$uv_g_scale_move)
	expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	# without positions there is no ridge
	expect_null(run(list(), 0L)$uv_g_scale_move)
	expect_error(run(list(uv_g_scale_move = NA), 1L), 'uv_g_scale_move must be TRUE or FALSE', fixed = TRUE)
})

test_that('ame() runs the position steps of the actor update without an additive margin', {
	set.seed(7402); nA = 9L; nB = 11L
	y = matrix(
		rbinom(nA * nB, 1, .15),
		nA,
		nB,
		dimnames = list(sprintf('r%02d', 1:nA), sprintf('c%02d', 1:nB))
	); y[1:2, ] = 0; y[, 1:3] = 0
	run = function(R, prior = list(actor_move_every = 2)) suppressMessages(ame(
		y,
		mode = 'bipartite',
		family = 'binary',
		R_row = R,
		R_col = R,
		rvar = FALSE,
		cvar = FALSE,
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = prior
	))
	on = run(1L)
	expect_equal(unname(on$actor_move['visits']), 40)
	expect_equal(unname(on$actor_move['level_calls']), 0)
	expect_equal(unname(on$actor_move['uv_calls']), 40 * (nA + nB))
	expect_equal(unname(on$actor_move['stuck']), 0)
	off = run(1L, list(actor_move = FALSE))
	expect_null(off$actor_move)
	expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	# no margin and no positions: nothing to move
	expect_null(run(0L)$actor_move)
})

lame_panel = function(seed, nA = 8L, nB = 10L, nt = 3L, dens = .12) {
	set.seed(seed)
	rn = sprintf('r%02d', 1:nA); cn = sprintf('c%02d', 1:nB)
	Y = lapply(
		1:nt,
		function(t) { y = matrix(
			rbinom(nA * nB, 1, dens),
			nA,
			nB,
			dimnames = list(rn, cn)
		); y[1:2, ] = 0; y[, 1:3] = 0; y }
	)
	names(Y) = paste0('p', 1:nt)
	X = lapply(1:nt, function(t) array(rnorm(nA * nB), c(nA, nB, 1), dimnames = list(rn, cn, 'x')))
	names(X) = names(Y)
	list(Y = Y, X = X, nA = nA, nB = nB)
}
lame_run = function(p, prior, R, dyn, ...) suppressMessages(lame(
	p$Y,
	Xdyad = p$X,
	mode = 'bipartite',
	family = 'binary',
	R_row = R,
	R_col = R,
	dynamic_ab = dyn,
	burn = 20,
	nscan = 60,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = prior,
	...
))

test_that('lame() draws the additive-effect hyperparameters given the totals as a prior switch and records it', {
	skip_on_cran()
	p = lame_panel(7403L)
	for (R in c(0L, 1L)) {
		expect_no_warning({ on = lame_run(p, list(), R, TRUE) })
		expect_equal(unname(on$ab_hyper_move['attempts']), 80 * 8)
		expect_gt(unname(on$ab_hyper_move['accepts']), 0)
		expect_true(
			all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$sigma_ab)) && all(abs(on$rho_ab) < 1)
		)
		off = lame_run(p, list(ab_hyper_move = FALSE), R, TRUE)
		expect_null(off$ab_hyper_move)
		expect_false(isTRUE(all.equal(on$sigma_ab, off$sigma_ab)))
	}
	# static effects have no hyperparameters of that kind
	expect_null(lame_run(p, list(), 0L, FALSE)$ab_hyper_move)
	expect_error(lame_run(p, list(ab_hyper_move = NA), 0L, TRUE), 'ab_hyper_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(
		lame_run(p, list(uv_g_scale_move = 1), 0L, TRUE),
		'uv_g_scale_move must be TRUE or FALSE',
		fixed = TRUE
	)
})

test_that('lame() draws along the scale ridge of static positions and not of moving ones', {
	skip_on_cran()
	p = lame_panel(7404L, dens = .15)
	for (dyn in c(FALSE, TRUE)) {
		expect_no_warning({ on = lame_run(p, list(), 1L, dyn) })
		expect_equal(unname(on$uv_g_scale_move['visits']), 80); expect_equal(unname(on$uv_g_scale_move['drawn']), 160)
		off = lame_run(p, list(uv_g_scale_move = FALSE), 1L, dyn)
		expect_null(off$uv_g_scale_move)
		expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	}
	expect_null(lame_run(p, list(), 1L, TRUE, dynamic_uv = TRUE)$uv_g_scale_move)
	expect_null(lame_run(p, list(), 0L, TRUE)$uv_g_scale_move)
})

test_that('lame() runs the position steps of the actor update without an additive margin', {
	skip_on_cran()
	p = lame_panel(7405L, dens = .15)
	for (dyn in c(FALSE, TRUE)) {
		on = lame_run(p, list(actor_move_every = 2), 1L, dyn, rvar = FALSE, cvar = FALSE)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['level_calls']), 0)
		expect_equal(unname(on$actor_move['uv_calls']), 40 * (p$nA + p$nB))
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
	}
})

test_that('the proposal of the hyperparameter update survives a checkpoint and a continued run uses it', {
	skip_on_cran()
	p = lame_panel(7406L)
	a = list(
		Y = p$Y,
		Xdyad = p$X,
		mode = 'bipartite',
		family = 'binary',
		dynamic_ab = TRUE,
		burn = 20L,
		nscan = 30L,
		odens = 5L,
		seed = 3L,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	state = f$start_vals$checkpoint_state$values
	expect_true(is.numeric(state$ab_hyper_par))
	expect_equal(unname(state$ab_hyper_par[3]), 50 * 8)
	a$start_vals = f$start_vals; a$burn = 0L; a$nscan = 10L
	a$checkpoint_path = tempfile(fileext = '.rds'); a$max_seconds = 1e-12
	suppressWarnings(suppressMessages(do.call(lame, a)))
	ck = readRDS(a$checkpoint_path); unlink(a$checkpoint_path)
	expect_equal(ck$start_vals$checkpoint_state$values$ab_hyper_par, state$ab_hyper_par, tolerance = 0)
	# a checkpoint without the field restores, and the update starts its proposal afresh
	old = f$start_vals; old$checkpoint_state$values$ab_hyper_par = NULL
	b = a; b$start_vals = old; b$checkpoint_path = NULL; b$max_seconds = NULL
	h = suppressWarnings(suppressMessages(do.call(lame, b)))
	expect_equal(unname(h$ab_hyper_move['attempts']), 10 * 8)
	# the continued run keeps counting and does not adapt again
	a$max_seconds = NULL; a$checkpoint_path = NULL
	g = suppressWarnings(suppressMessages(do.call(lame, a)))
	expect_equal(unname(g$ab_hyper_move['attempts']), 60 * 8)
	expect_equal(g$start_vals$checkpoint_state$values$ab_hyper_par[-(1:4)], state$ab_hyper_par[-(1:4)], tolerance = 0)
})

test_that('lame() draws the additive-effect hyperparameters with changing coefficients as well', {
	skip_on_cran()
	p = lame_panel(7407L, dens = .15)
	on = suppressMessages(lame(
		p$Y,
		Xdyad = p$X,
		mode = 'bipartite',
		family = 'binary',
		dynamic_ab = TRUE,
		dynamic_beta = 'x_dyad',
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(unname(on$ab_hyper_move['attempts']), 80 * 8)
	expect_true(all(is.finite(on$sigma_ab)) && all(abs(on$rho_ab) < 1))
})

test_that('the records of the new updates are kept per chain when chains are pooled', {
	f = list(ab_hyper_move = c(attempts = 8, accepts = 2, rate = .25), uv_g_scale_move = c(visits = 5, drawn = 10))
	pooled = .lame_pool_sampler_diagnostics(f, list(f, f))
	expect_null(pooled$ab_hyper_move); expect_null(pooled$uv_g_scale_move)
	expect_equal(pooled$sampler_diagnostics_by_chain[[2]]$ab_hyper_move, f$ab_hyper_move)
	expect_equal(pooled$sampler_diagnostics_by_chain[[1]]$uv_g_scale_move, f$uv_g_scale_move)
})
