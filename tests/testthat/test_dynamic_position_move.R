skip_on_cran()

test_that('the covariance of a path is the AR(1) law, with local innovation scales when they are given', {
	nt = 5; rho = .7; sg = .6
	cv = .lame_dpm_cov(nt, rho, sg)
	expect_equal(cv, sg^2 / (1 - rho^2) * rho^abs(outer(1:nt, 1:nt, '-')))
	lam = c(1, .5, 2, 1.5, .8)
	cl = .lame_dpm_cov(nt, rho, sg, lam)
	# the precision of the chain: stationary start, transitions N(rho u_{t-1}, sg^2 / lambda_t)
	Pm = matrix(0, nt, nt)
	Pm[1, 1] = (1 - rho^2) / sg^2
	for (t in 2:nt) {
		a = lam[t] / sg^2
		Pm[t, t] = Pm[t, t] + a; Pm[t - 1, t - 1] = Pm[t - 1, t - 1] + a * rho^2
		Pm[t, t - 1] = Pm[t, t - 1] - a * rho; Pm[t - 1, t] = Pm[t - 1, t] - a * rho
	}
	expect_equal(solve(cl), Pm, tolerance = 1e-8)
	# the root reproduces it
	rt = .lame_dpm_root(nt, 2L, rho, sg, lam)
	expect_equal(rt %*% t(rt), kronecker(cl, diag(2)))
})

test_that('the Poisson factor of a path move equals the change of the log likelihood of the cells of its row', {
	set.seed(4)
	nl = 7; nt = 3; rk = 2
	w = array(rnorm(nl * rk * nt), c(nl, rk, nt))
	z = matrix(rnorm(nl * nt, -.5, .5), nl, nt); off = c(0, .3, -.1)
	y = matrix(rpois(nl * nt, exp(z + rep(off, each = nl))), nl, nt)
	ok = matrix(TRUE, nl, nt); ok[2, 1] = FALSE; ok[5, 3] = FALSE
	yl = ifelse(ok, y, 0); el = exp(z + rep(off, each = nl)); el[!ok] = 0
	d = matrix(rnorm(rk * nt, 0, .5), rk, nt)
	dl = vapply(seq_len(nt), function(t) drop(w[, , t] %*% d[, t]), numeric(nl))
	ll = function(zz) sum(dpois(y[ok], exp(zz[ok] + rep(off, each = nl)[ok]), log = TRUE))
	expect_equal(.lame_dpm_poisson_ll(yl, el, w, d), ll(z + dl) - ll(z), tolerance = 1e-10)
	expect_equal(.lame_dpm_poisson_ll(yl, el, w, matrix(0, rk, nt)), 0)
})

test_that('a visit keeps every residual and rebuilds the predictor, in every mode', {
	for (mode in c('dir', 'sym', 'bip')) for (rk in 1:2) {
		s = dpmt_state(mode, rk)
		for (rep in 1:3) {
			o = .lame_dpm_visit(s$cx, s$Z, s$EZ, s$U, s$V, s$G, s$Ld, s$rho_uv, s$sigma_uv)
			expect_lt(max(abs((o$Z - o$EZ) - (s$Z - s$EZ))), 1e-12)
			expect_lt(max(abs(dpmt_ez(mode, s$base, o$U, if (mode == 'sym') NULL else o$V, s$G, s$Ld) - o$EZ)), 1e-12)
			expect_equal(unname(o$stats['dpm_stuck']), 0)
			if (mode == 'sym') expect_lt(max(abs(o$Z - aperm(o$Z, c(2, 1, 3)))), 1e-12)
			expect_gt(max(abs(o$U - s$U)), 0)
			s$U = o$U; s$V = o$V; s$Z = o$Z; s$EZ = o$EZ
		}
		# every actor of every side is visited once (a symmetric network's actors once)
		n_paths = if (mode == 'sym') dim(s$U)[1L] else dim(s$U)[1L] + dim(s$cx$Y)[2L]
		expect_equal(unname(o$stats['dpm_calls']), n_paths)
	}
})

test_that('the t kind reads the local scales, and a visit leaves the coordinates inside the rail', {
	s = dpmt_state('dir', 2L)
	nr = dim(s$U)[1]; nt = dim(s$U)[3]
	lu = matrix(rgamma(nr * nt, 2, 2), nr, nt); lv = matrix(rgamma(nr * nt, 2, 2), nr, nt)
	o = .lame_dpm_visit(
		s$cx,
		s$Z,
		s$EZ,
		s$U,
		s$V,
		rho_uv = s$rho_uv,
		sigma_uv = s$sigma_uv,
		lambda_u = lu,
		lambda_v = lv
	)
	expect_lt(max(abs((o$Z - o$EZ) - (s$Z - s$EZ))), 1e-12)
	expect_gt(max(abs(o$U - s$U)), 0)
	# a rail the start state fits in: no coordinate leaves it, whatever the draws
	rail = max(abs(c(s$U, s$V))) + .2
	for (rep in 1:30) {
		o = .lame_dpm_visit(s$cx, s$Z, s$EZ, s$U, s$V, rho_uv = s$rho_uv, sigma_uv = s$sigma_uv, max_abs = rail)
		expect_lte(max(abs(c(o$U, o$V))), rail)
		s$U = o$U; s$V = o$V; s$Z = o$Z; s$EZ = o$EZ
	}
	# a persistence outside (-1, 1) leaves the state alone
	o = .lame_dpm_visit(s$cx, s$Z, s$EZ, s$U, s$V, rho_uv = 1.2, sigma_uv = s$sigma_uv)
	expect_identical(o$U, s$U); expect_equal(unname(o$stats['dpm_calls']), 0)
})

test_that('lame() moves the position paths of a Poisson fit with the latent array and honours the switches', {
	for (mode in c('dir', 'bip')) {
		p = dpmt_panel(mode)
		on = dpmt_run(p, mode); off = dpmt_run(p, mode, list(actor_move = FALSE))
		# the bipartite scale step is made every sweep (it is not a visit), so it goes with its own switch
		never = dpmt_run(
			p,
			mode,
			c(list(actor_move_every = 1000), if (mode == 'bip') list(bipartite_scale_move = FALSE))
		)
		am = on$actor_move
		n_paths = if (mode == 'bip') 7 + 9 else 16
		expect_equal(unname(am['dpm_calls']), 20 * n_paths)
		expect_equal(unname(am['dpm_stuck']), 0)
		expect_gt(unname(am['dpm_evals']), unname(am['dpm_calls']))
		expect_null(off$actor_move); expect_null(never$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC); expect_identical(never$rho_uv, off$rho_uv)
		expect_false(isTRUE(all.equal(on$rho_uv, off$rho_uv)))
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$U)))
	}
})

test_that('lame() runs the path move for symmetric fits, the t kind, dynamic_ab and a rail, and other families run without it', {
	skip_on_cran()
	p = dpmt_panel('sym')
	sy = dpmt_run(p, 'sym')
	expect_equal(unname(sy$actor_move['dpm_calls']), 20 * 8)
	p = dpmt_panel('dir')
	tt = dpmt_run(p, 'dir', dynamic_uv_kind = 't')
	expect_equal(unname(tt$actor_move['dpm_calls']), 20 * 16); expect_equal(unname(tt$actor_move['dpm_stuck']), 0)
	ab = dpmt_run(p, 'dir', dynamic_ab = TRUE)
	expect_equal(unname(ab$actor_move['dpm_calls']), 20 * 16)
	rl = dpmt_run(p, 'dir', list(uv_max_abs = 1.5), posterior_opts = list(save_UV = TRUE, thin_UV = 1L))
	expect_lte(max(abs(c(rl$U_draws, rl$V_draws))), 1.5)
	# the snap kind and the other families keep the sampler's own updates
	sn = dpmt_run(p, 'dir', dynamic_uv_kind = 'snap')
	expect_false('dpm_calls' %in% names(sn$actor_move))
	pb = dpmt_panel('dir', binary = TRUE)
	bn = dpmt_run(pb, 'dir', family = 'binary', dynamic_uv_kind = 'snap')
	expect_false('dpm_calls' %in% names(bn$actor_move))
	po = lapply(p$Y, function(y) { y[!is.na(y)] = pmin(y[!is.na(y)], 3); y })
	od = suppressMessages(
		suppressWarnings(
			lame(
				po,
				Xdyad = p$X,
				family = 'ordinal',
				dynamic_uv = TRUE,
				R = 1L,
				burn = 10,
				nscan = 30,
				odens = 5,
				seed = 3,
				verbose = FALSE,
				plot = FALSE,
				gof = FALSE
			)
		)
	)
	expect_false('dpm_calls' %in% names(od$actor_move))
})

test_that('the path move leaves the exact posterior of a tiny panel invariant', {
	skip_on_cran()
	cases = list(
		list(mode = 'dir', kind = 'ar1', rk = 1L, ndraw = 2e5),
		list(mode = 'dir', kind = 't', rk = 2L, ndraw = 2e5),
		list(mode = 'sym', kind = 'ar1', rk = 1L, ndraw = 2e5),
		list(mode = 'bip', kind = 'ar1', rk = 2L, ndraw = 2e5),
		list(mode = 'dir', kind = 'ar1', rk = 1L, ndraw = 2e5, rail = 1.2)
	)
	for (k in seq_along(cases)) {
		r = do.call(dpmt_exact, c(cases[[k]], list(seed = 700L + k)))
		expect_gt(r$ess, 400)
		expect_lt(r$maxz, 4.6)
	}
})

test_that('the exact-posterior check detects a move that drops its data factor', {
	skip_on_cran()
	ess = .lame_ess_step
	local_mocked_bindings(
		.lame_ess_step = function(x, mu, root, loglik, ll0) ess(x, mu, root, function(x) 0, 0),
		.package = 'lame'
	)
	expect_gt(dpmt_exact('dir', 'ar1', 1L, ndraw = 2e5, seed = 9L)$maxz, 6)
})
test_that('the exact-posterior check detects a mis-scaled shift', {
	skip_on_cran()
	delta = .lame_dpm_delta
	local_mocked_bindings(.lame_dpm_delta = function(w, d) 1.3 * delta(w, d), .package = 'lame')
	expect_gt(dpmt_exact('dir', 'ar1', 1L, ndraw = 2e5, seed = 10L)$maxz, 6)
})
test_that('the exact-posterior check detects a mis-stated persistence of the prior', {
	skip_on_cran()
	root = .lame_dpm_root
	local_mocked_bindings(
		.lame_dpm_root = function(nt, rk, rho, sigma, lambda = NULL) root(nt, rk, .95 * rho, sigma, lambda),
		.package = 'lame'
	)
	expect_gt(dpmt_exact('dir', 'ar1', 1L, ndraw = 2e5, seed = 11L)$maxz, 4.6)
})

test_that('the probit factor of a path move equals the log likelihood of the observed cells of its row', {
	s = dpmt_state_bin('dir', 2L)
	nr = dim(s$EZ)[1]; nt = dim(s$EZ)[3]
	sd_cell = sqrt(1 - .4)
	off = array(rnorm(length(s$EZ), 0, .3), dim(s$EZ))
	for (i in c(1L, 3L, 5L)) for (row in c(TRUE, FALSE)) {
		sd_ = if (row) s$cx$rows[[i]] else s$cx$cols[[i]]
		line = function(a) matrix(if (row) a[i, , ] else a[, i, ], nr, nt)
		b = (line(s$EZ) + line(off))[sd_$idx]
		d = matrix(rnorm(2 * nt, 0, .5), 2, nt)
		w = if (row) s$V else s$U
		Un = s$U; Vn = s$V
		if (row) Un[i, , ] = Un[i, , ] + d else Vn[i, , ] = Vn[i, , ] + d
		EZn = dpmt_ez('dir', s$base, Un, Vn, NULL, NULL)
		obs = is.finite(if (row) s$Y[i, , ] else s$Y[, i, ]); obs[i, ] = FALSE
		yl = matrix(if (row) s$Y[i, , ] else s$Y[, i, ], nr, nt)
		brute = sum(pnorm((2 * yl[obs] - 1) * (line(EZn) + line(off))[obs] / sd_cell, log.p = TRUE))
		expect_equal(.lame_dpm_probit_ll(b, sd_$sg, sd_$idx, w, d, sd_cell), brute, tolerance = 1e-10)
	}
})

test_that('a binary visit keeps the signs of the data, rebuilds the predictor and moves every path, in every mode', {
	for (mode in c('dir', 'sym', 'bip')) for (rk in 1:2) {
		s = dpmt_state_bin(mode, rk)
		cx = s$cx
		nt = dim(s$U)[3L]
		obs = is.finite(s$Y)
		if (mode != 'bip') for (t in seq_len(nt)) diag(obs[, , t]) = FALSE
		for (rep in 1:3) {
			o = .lame_dpm_visit_binary(cx, s$Z, s$EZ, s$Y, s$rho, s$U, s$V, s$G, s$Ld, s$rho_uv, s$sigma_uv)
			expect_lt(max(abs(dpmt_ez(mode, s$base, o$U, if (mode == 'sym') NULL else o$V, s$G, s$Ld) - o$EZ)), 1e-12)
			expect_true(all(o$Z[obs & s$Y == 1] > 0)); expect_true(all(o$Z[obs & s$Y == 0] <= 0))
			expect_equal(unname(o$stats['dpm_stuck']), 0)
			if (mode == 'sym') expect_lt(max(abs(o$Z - aperm(o$Z, c(2, 1, 3)))), 1e-12)
			expect_gt(max(abs(o$U - s$U)), 0)
			s$U = o$U; s$V = o$V; s$Z = o$Z; s$EZ = o$EZ
		}
		n_paths = if (mode == 'sym') dim(s$U)[1L] else dim(s$U)[1L] + dim(s$Y)[2L]
		expect_equal(unname(o$stats['dpm_calls']), n_paths)
	}
})

test_that('the binary context carries the fields of the driver and reuses an existing context', {
	s = dpmt_state_bin('bip')
	expect_equal(c(s$cx$mode, s$cx$nr, s$cx$nc, s$cx$N), c('bip', '5', '6', '3'))
	base = .lame_actor_move_context(dpmt_state_bin('dir')$Y)
	cx = .lame_dpm_binary_context(NULL, 'dir', base)
	expect_identical(cx$rows, base$rows); expect_equal(c(cx$nr, cx$nc), c(6L, 6L)); expect_identical(cx$mode, 'dir')
})

test_that('lame() moves the position paths of a binary fit with the latent array integrated out and honours the switches', {
	skip_on_cran()
	for (mode in c('dir', 'sym', 'bip')) {
		p = dpmt_panel(mode, binary = TRUE)
		on = dpmt_run(p, mode, family = 'binary'); off = dpmt_run(p, mode, list(actor_move = FALSE), family = 'binary')
		# the draws of the symmetric gauge are made every sweep (they are not visits), so they go with the paths
		never = dpmt_run(
			p,
			mode,
			c(list(actor_move_every = 1000), if (mode == 'sym') list(actor_move_paths = FALSE)),
			family = 'binary'
		)
		am = on$actor_move
		n_paths = if (mode == 'bip') 7 + 9 else if (mode == 'sym') 8 else 16
		expect_equal(unname(am['dpm_calls']), 20 * n_paths)
		expect_equal(unname(am['dpm_stuck']), 0)
		expect_gt(unname(am['dpm_evals']), unname(am['dpm_calls']))
		expect_null(off$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC); expect_identical(never$rho_uv, off$rho_uv)
		expect_false(isTRUE(all.equal(on$rho_uv, off$rho_uv)))
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)) && all(is.finite(on$U)))
	}
})

test_that('lame() runs the binary path move for the t kind and a rail, and leaves a directed fit with the paired-actor step to it', {
	skip_on_cran()
	p = dpmt_panel('dir', binary = TRUE)
	tt = dpmt_run(p, 'dir', family = 'binary', dynamic_uv_kind = 't')
	expect_equal(unname(tt$actor_move['dpm_calls']), 20 * 16); expect_equal(unname(tt$actor_move['dpm_stuck']), 0)
	rl = dpmt_run(
		p,
		'dir',
		list(uv_max_abs = 1.5),
		family = 'binary',
		posterior_opts = list(save_UV = TRUE, thin_UV = 1L)
	)
	expect_lte(max(abs(c(rl$U_draws, rl$V_draws))), 1.5)
	# rank two with dynamic additive effects: the paired-actor slice step of the sampler moves the paths
	sl = dpmt_run(p, 'dir', family = 'binary', rk = 2L, dynamic_ab = TRUE)
	expect_false('dpm_calls' %in% names(sl$actor_move))
})

test_that('the binary path move leaves the exact posterior of a tiny panel invariant', {
	skip_on_cran()
	cases = list(
		list(mode = 'dir', kind = 'ar1', rk = 1L, ndraw = 2e5),
		list(mode = 'dir', kind = 't', rk = 2L, ndraw = 2e5),
		list(mode = 'sym', kind = 'ar1', rk = 1L, ndraw = 2e5),
		list(mode = 'bip', kind = 'ar1', rk = 2L, ndraw = 2e5),
		list(mode = 'dir', kind = 'ar1', rk = 1L, ndraw = 2e5, rail = 1.2)
	)
	for (k in seq_along(cases)) {
		r = do.call(dpmt_exact_binary, c(cases[[k]], list(seed = 800L + k)))
		expect_gt(r$ess, 400)
		expect_lt(r$maxz, 4.6)
	}
})

test_that('the binary exact-posterior check detects a move that drops its data factor, flips the sign of a cell or ignores the dyad factor', {
	skip_on_cran()
	ess = .lame_ess_step
	local_mocked_bindings(
		.lame_ess_step = function(x, mu, root, loglik, ll0) ess(x, mu, root, function(x) 0, 0),
		.package = 'lame'
	)
	expect_gt(dpmt_exact_binary('dir', 'ar1', 1L, ndraw = 2e5, seed = 19L)$maxz, 5)
})
test_that('the binary exact-posterior check detects a flipped sign of the cell factors', {
	skip_on_cran()
	pl = .lame_dpm_probit_ll
	local_mocked_bindings(
		.lame_dpm_probit_ll = function(b, sg, idx, w, d, sd_cell = 1) pl(b, -sg, idx, w, d, sd_cell),
		.package = 'lame'
	)
	expect_gt(dpmt_exact_binary('sym', 'ar1', 1L, ndraw = 2e5, seed = 20L)$maxz, 5)
})
test_that('the binary exact-posterior check detects a move that ignores the shared dyad factor', {
	skip_on_cran()
	local_mocked_bindings(.lame_dyad_factor_draw = function(Z, EZ, rho, upper_idx) array(0, dim(Z)), .package = 'lame')
	expect_gt(dpmt_exact_binary('dir', 'ar1', 1L, ndraw = 2e5, seed = 21L)$maxz, 4.6)
})

test_that('the scale ridge of the symmetric gauge keeps the product and leaves the prior of the paths and weights invariant', {
	skip_on_cran()
	set.seed(11)
	n = 5; nt = 3; rho = .9
	draw = function() {
		u = array(0, c(n, 1, nt)); u[, 1, 1] = rnorm(n)
		for (t in 2:nt) u[, 1, t] = rho * u[, 1, t - 1] + sqrt(1 - rho^2) * rnorm(n)
		list(U = u, Ld = rnorm(1))
	}
	stat = function(s) { u = matrix(s$U[, 1, ], n, nt); ss = sum(u^2); c(
		ss,
		s$Ld,
		s$Ld^2,
		s$Ld * ss,
		sum(u[, 1] * u[, nt]),
		ss^2,
		abs(s$Ld) * ss^2,
		s$Ld * u[1, 1] * u[2, nt]
	) }
	right = function(s) { o = .lame_dpm_sym_scale(s$U, s$Ld, rho); list(U = o$U, Ld = o$Ld) }
	# the last entry of a replicate is the largest change of the product of the weight and the paths
	D = t(replicate(20000, {
		s0 = draw(); s1 = right(s0)
		c(stat(s1) - stat(s0), max(abs(s1$Ld * tcrossprod(s1$U[, 1, 2]) - s0$Ld * tcrossprod(s0$U[, 1, 2]))))
	}))
	expect_lt(max(D[, ncol(D)]), 1e-10)
	D = D[, -ncol(D)]
	z = colMeans(D) / (apply(D, 2, sd) / sqrt(nrow(D)))
	expect_lt(max(abs(z)), 4.2)
	# a kernel that forgets the Jacobian of the map does not leave the prior invariant, and the check sees it
	wrong = function(s) {
		cc = exp(rnorm(1, 0, .4))
		q = sum(s$U[, 1, 1]^2) + sum((s$U[, 1, -1] - rho * s$U[, 1, -nt])^2) / (1 - rho^2)
		if (log(runif(1)) < -(cc^2 - 1) / 2 * q - 0.5 * s$Ld^2 * (cc^-4 - 1)) list(
			U = cc * s$U,
			Ld = s$Ld / cc^2
		) else s
	}
	Dw = t(replicate(20000, { s0 = draw(); stat(wrong(s0)) - stat(s0) }))
	zw = colMeans(Dw) / (apply(Dw, 2, sd) / sqrt(nrow(Dw)))
	expect_gt(max(abs(zw)), 8)
	# the rail rejects a scaled path that leaves the box
	s0 = draw(); big = max(abs(s0$U))
	for (rep in 1:50) { o = .lame_dpm_sym_scale(
		s0$U,
		s0$Ld,
		rho,
		max_abs = big * 1.01
	); expect_lte(max(abs(o$U)), big * 1.01) }
})

test_that('the persistence draw of the symmetric gauge has the conditional of the transitions as its stationary law', {
	skip_on_cran()
	set.seed(12)
	n = 12; nt = 4; rho0 = .94
	u = array(0, c(n, 1, nt)); u[, 1, 1] = rnorm(n)
	for (t in 2:nt) u[, 1, t] = rho0 * u[, 1, t - 1] + sqrt(1 - rho0^2) * rnorm(n)
	m = n * (nt - 1)
	s11 = sum(sapply(2:nt, function(t) sum(u[, 1, t - 1]^2)))
	s12 = sum(sapply(2:nt, function(t) sum(u[, 1, t - 1] * u[, 1, t])))
	s22 = sum(sapply(2:nt, function(t) sum(u[, 1, t]^2)))
	# the exact conditional on a fine grid (normal prior N(.9, .1) truncated to (-.99, .99))
	grid = seq(-.9899, .9899, length.out = 40001)
	lp = -0.5 * m * log(1 - grid^2) - (s22 - 2 * grid * s12 + grid^2 * s11) / (2 * (1 - grid^2)) + dnorm(grid, .9, .1, log = TRUE)
	w = exp(lp - max(lp)); w = w / sum(w)
	gm = sum(w * grid)
	ps = c(.1, .5, .9)
	gq = sapply(ps, function(p) grid[which(cumsum(w) >= p)[1]])
	r = numeric(41000); rc = .9
	for (k in seq_along(r)) { rc = .lame_dpm_rho_unit_mh(u, rc, .9, .1); r[k] = rc }
	r = r[-(1:1000)]
	nb = 50; bl = floor(length(r) / nb)
	bse = function(x) sd(colMeans(matrix(x[seq_len(nb * bl)], bl))) / sqrt(nb)
	expect_lt(abs(mean(r) - gm), 4.5 * bse(r))
	for (k in seq_along(ps)) { ind = as.numeric(r <= gq[k]); expect_lt(abs(mean(ind) - ps[k]), 4.5 * bse(ind)) }
	# the walk moves where the fixed-scale draw of the host hardly does: its acceptance rate is far above the host's
	moves = function(kernel) { x = numeric(4000); rc = .9; for (k in seq_along(x)) { rc = kernel(u, rc, .9, .1); x[k] = rc }; mean(diff(x) != 0) }
	expect_gt(moves(.lame_dpm_rho_unit_mh), 2 * moves(.lame_rho_uv_unit_mh))
})

test_that('lame() draws the scale ridge and the persistence of a symmetric Poisson or binary fit and honours the switches', {
	for (family in c('poisson', 'binary')) {
		p = dpmt_panel('sym', binary = family == 'binary')
		on = dpmt_run(p, 'sym', family = family)
		paths_off = dpmt_run(p, 'sym', list(actor_move_paths = FALSE), family = family)
		off = dpmt_run(p, 'sym', list(actor_move = FALSE), family = family)
		# 40 sweeps, three proposals of rank one each
		expect_gt(
			unname(on$actor_move['dpm_scale_attempts']),
			100
		); expect_lte(unname(on$actor_move['dpm_scale_attempts']), 120)
		expect_equal(unname(on$actor_move['dpm_scale_attempts']) %% 3, 0)
		expect_gt(unname(on$actor_move['dpm_scale_accepts']), 0)
		expect_false(any(grepl('dpm', names(paths_off$actor_move))))
		expect_null(off$actor_move)
		expect_true(all(is.finite(on$rho_uv)) && all(on$rho_uv < .99))
	}
	expect_error(dpmt_run(dpmt_panel('sym'), 'sym', list(actor_move_paths = NA)), 'actor_move_paths')
})

test_that('the bipartite scale step runs for Poisson and ordinal fits with the actor moves and honours its switches', {
	skip_on_cran()
	p = dpmt_panel('bip')
	for (family in c('poisson', 'ordinal')) {
		if (family == 'ordinal') p$Y = lapply(p$Y, function(y) pmin(y, 3))
		on = dpmt_run(p, 'bip', family = family)
		expect_gt(unname(on$bipartite_scale_move['attempts']), 0); expect_lte(
			unname(on$bipartite_scale_move['attempts']),
			40
		)
		expect_equal(unname(on$bipartite_scale_move['numerical_rejections']), 0)
		# no scale step without the actor moves (their draws are then the sampler's own), nor with its own switch off
		expect_null(dpmt_run(p, 'bip', list(actor_move = FALSE), family = family)$bipartite_scale_move)
		expect_null(dpmt_run(p, 'bip', list(bipartite_scale_move = FALSE), family = family)$bipartite_scale_move)
	}
})

test_that('a bipartite visit handles different ranks of the two sides, for the Poisson and the binary family', {
	set.seed(21)
	nr = 5L; nc = 6L; nt = 3L; ra = 1L; rb = 2L
	U = array(
		rnorm(nr * ra * nt, 0, .6),
		c(nr, ra, nt)
	); V = array(rnorm(nc * rb * nt, 0, .6), c(nc, rb, nt)); G = matrix(c(.9, -.5), ra, rb)
	base = array(rnorm(nr * nc * nt, -.4, .4), c(nr, nc, nt))
	EZ = dpmt_ez('bip', base, U, V, G, NULL)
	Yp = array(rpois(nr * nc * nt, .6), c(nr, nc, nt)); Yp[1:2, , ] = 0; Yp[3, 4, 2] = NA
	cp = .lame_ltm_context(Yp, 'poisson', 'bip', offset = c(0, .2, -.3))
	Z = EZ + array(rnorm(length(EZ), 0, .3), dim(EZ))
	for (rep in 1:3) {
		o = .lame_dpm_visit(cp, Z, EZ, U, V, G, NULL, .8, .5)
		expect_lt(max(abs((o$Z - o$EZ) - (Z - EZ))), 1e-12)
		expect_lt(max(abs(dpmt_ez('bip', base, o$U, o$V, G, NULL) - o$EZ)), 1e-12)
		expect_equal(unname(o$stats['dpm_calls']), nr + nc)
		U = o$U; V = o$V; Z = o$Z; EZ = o$EZ
	}
	Yb = array(rbinom(nr * nc * nt, 1, .4), c(nr, nc, nt)); Yb[1:2, , ] = 0; Yb[3, 4, 2] = NA
	cb = .lame_dpm_binary_context(Yb, 'bip')
	Zb = rZ_bin_bip_batch_cpp(array(0, dim(EZ)), EZ, Yb)
	for (rep in 1:3) {
		o = .lame_dpm_visit_binary(cb, Zb, EZ, Yb, 0, U, V, G, NULL, .8, .5)
		expect_lt(max(abs(dpmt_ez('bip', base, o$U, o$V, G, NULL) - o$EZ)), 1e-12)
		obs = is.finite(Yb)
		expect_true(all(o$Z[obs & Yb == 1] > 0)); expect_true(all(o$Z[obs & Yb == 0] <= 0))
		U = o$U; V = o$V; Zb = o$Z; EZ = o$EZ
	}
})
