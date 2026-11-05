# translations of the levels, coefficients, scale and location with the latent array (poisson, ordinal)
skip_on_cran()

# feasible state of a small panel: predictor, latent array, data and levels, centred or free
ltmt_state = function(mode, family, centred, dyn, seed = 7L, n = 6L, N = 3L, nb = 8L, K = 4L, nmiss = 4L) {
	set.seed(seed)
	bip = mode == 'bip'; sym = mode == 'sym'
	nr = n; nc = if (bip) nb else n
	x = array(rnorm(nr * nc * N), c(nr, nc, N))
	if (sym) for (t in 1:N) { m = x[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; x[, , t] = m }
	cen = function(v) if (centred) v - mean(v) else v
	la = cen(rnorm(nr)); lb = if (sym) la else cen(rnorm(nc))
	beta = if (centred) c(-1, .5) else .5
	dev_a = matrix(rnorm(nr * N, 0, .2), nr, N); dev_b = if (sym) dev_a else matrix(rnorm(nc * N, 0, .2), nc, N)
	a_mat = if (dyn) la + dev_a else NULL
	b_mat = if (dyn) { if (sym) a_mat else lb + dev_b } else NULL
	build = function(la, lb, beta, a_mat, b_mat) {
		EZ = array(0, c(nr, nc, N))
		for (t in 1:N) {
			at = if (dyn) a_mat[, t] else la; bt = if (dyn) b_mat[, t] else lb
			EZ[, , t] = (if (centred) beta[1] else 0) + beta[length(beta)] * x[, , t] + outer(at, bt, '+')
		}
		EZ
	}
	EZ = build(la, lb, beta, a_mat, b_mat)
	E = array(rnorm(nr * nc * N, 0, .3), c(nr, nc, N))
	if (sym) for (t in 1:N) { m = E[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; E[, , t] = m }
	Z = EZ + E
	Y = if (family == 'poisson') array(rpois(nr * nc * N, exp(Z)), c(nr, nc, N)) else
		array(
			as.numeric(cut(Z, quantile(Z, seq(0, 1, length.out = K + 1)), include.lowest = TRUE, labels = FALSE)),
			c(nr, nc, N)
		)
	if (!bip) for (t in 1:N) diag(Y[, , t]) = NA
	if (nmiss > 0) Y[sample(length(Y), nmiss)] = NA
	if (sym) for (t in 1:N) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
	list(
		Z = Z,
		EZ = EZ,
		Y = Y,
		la = la,
		lb = lb,
		beta = beta,
		a_mat = a_mat,
		b_mat = b_mat,
		build = build,
		x = x,
		n = nr,
		nc = nc,
		N = N
	)
}
ltmt_sab = function(mode) if (mode == 'bip') diag(c(1.1, .8)) else if (mode == 'sym') matrix(1.1, 2, 2) else matrix(c(1, .3, .3, 1.2), 2)
# every category's cells stay ordered in every period
ltmt_ranks_hold = function(Y, Z) {
	for (t in seq_len(dim(Y)[3L])) {
		y = Y[, , t]; z = Z[, , t]; ok = is.finite(y)
		for (k in sort(unique(y[ok]))[-1L]) if (!(max(z[ok & y < k]) < min(z[ok & y >= k]))) return(FALSE)
	}
	TRUE
}

test_that('the context indexes the cells that carry a data factor', {
	set.seed(3)
	Y = array(rpois(5 * 5 * 3, 1), c(5, 5, 3)); Y[2, 3, 1] = NA; Y[, 4, 2] = NA
	cx = .lame_ltm_context(Y, 'poisson', 'dir', offset = c(0, .5, -.5))
	for (t in 1:3) diag(Y[, , t]) = NA
	expect_equal(length(cx$ok_idx), sum(is.finite(Y)))
	expect_equal(cx$y_ok, Y[cx$ok_idx])
	expect_equal(cx$off_ok, c(0, .5, -.5)[arrayInd(cx$ok_idx, dim(Y))[, 3]])
	# a symmetric pair counts once, with the upper cell's value or else the lower one
	Ys = Y; for (t in 1:3) { m = Ys[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Ys[, , t] = m }
	Ys[1, 2, 1] = NA; Ys[2, 1, 1] = 7
	cs = .lame_ltm_context(Ys, 'poisson', 'sym')
	expect_equal(cs$Y[1, 2, 1], 7); expect_equal(cs$Y[2, 1, 1], 7)
	idx = arrayInd(cs$ok_idx, dim(Ys)); expect_true(all(idx[, 1] < idx[, 2]))
	# ordinal codes start at 1, a missing cell has none, and category cells are listed per period
	Yo = array(sample(c(2, 5, 9), 4 * 4 * 2, TRUE), c(4, 4, 2)); Yo[1, 2, 1] = NA
	co = .lame_ltm_context(Yo, 'ordinal', 'dir')
	expect_equal(co$K, 3L); expect_true(is.na(co$W[1, 2, 1])); expect_true(all(is.na(co$W[cbind(1:4, 1:4, 1)])))
	for (t in 1:2) for (k in 1:3) expect_true(all(co$W[, , t][co$cells[[t]][[k]]] == k))
	# a symmetric ordinal panel with differing triangles is not served; an explicit-cutpoint fit passes its codes
	Yz = Yo[, , 1, drop = FALSE]
	expect_false(.lame_ltm_context(Yz, 'ordinal', 'sym')$ok)
	m = Yz[, , 1]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Yz[, , 1] = m
	expect_true(.lame_ltm_context(Yz, 'ordinal', 'sym')$ok)
	expect_equal(.lame_ltm_context(Yo, 'ordinal', 'dir', Y_int = match(Yo, c(2, 5, 9)))$W, co$W)
})

test_that('the Poisson factor of a translation equals the change of the log likelihood of its cells', {
	for (mode in c('dir', 'sym', 'bip')) {
		s = ltmt_state(mode, 'poisson', TRUE, FALSE, nmiss = 5L)
		off = c(0, .3, -.2)
		cx = .lame_ltm_context(s$Y, 'poisson', mode, offset = off)
		for (row in c(TRUE, FALSE)) for (i in c(1L, 4L)) {
			if (mode == 'sym' && !row) next
			ph = .lame_ltm_pois_phase(cx, s$Z, row)
			ps = c(ys = ph$ys[i], s = ph$S[i])
			# the cells of the actor with data: its row or column, a symmetric pair once
			cells = which(cx$Y == cx$Y, arr.ind = TRUE)
			keep = if (row || mode == 'sym') cells[, 1] == i else cells[, 2] == i
			if (mode == 'sym') keep = cells[, 1] == i | cells[, 2] == i
			cells = cells[keep, , drop = FALSE]
			if (mode == 'sym') cells = cells[cells[, 1] == i, , drop = FALSE]
			ll = function(d) sum(dpois(cx$Y[cells], exp(s$Z[cells] + d + off[cells[, 3]]), log = TRUE))
			for (d in c(-1.2, .4, 1.5)) expect_equal(
				ll(d) - ll(0),
				unname(ps['ys'] * d - ps['s'] * expm1(d)),
				tolerance = 1e-10
			)
		}
	}
})

test_that('a level translation keeps every residual, rebuilds the predictor and the sums, and keeps the ranks', {
	skip_on_cran()
	for (family in c('poisson', 'ordinal')) for (mode in c('dir', 'sym', 'bip')) for (centred in c(TRUE, FALSE)) for (dyn in c(FALSE, TRUE)) {
		s = ltmt_state(mode, family, centred, dyn)
		cx = .lame_ltm_context(s$Y, family, mode, offset = if (family == 'poisson') c(0, .2, -.1) else NULL)
		iV0 = diag(length(s$beta)) / 4
		lev = c(TRUE, mode != 'sym')
		for (rep in 1:5) {
			o = .lame_ltm_visit(
				cx,
				s$Z,
				s$EZ,
				s$la,
				s$lb,
				s$a_mat,
				s$b_mat,
				s$beta,
				ltmt_sab(mode),
				iV0,
				centred = centred,
				lev = lev
			)
			if (mode == 'sym') { o$lb = o$la; o$b_mat = o$a_mat }
			expect_lt(max(abs((o$Z - o$EZ) - (s$Z - s$EZ))), 1e-12)
			expect_lt(max(abs(s$build(o$la, o$lb, o$beta, o$a_mat, o$b_mat) - o$EZ)), 1e-12)
			if (centred) expect_lt(max(abs(sum(o$la)), abs(sum(o$lb))), 1e-12)
			if (family == 'ordinal') expect_true(ltmt_ranks_hold(s$Y, o$Z))
			expect_equal(unname(o$stats['stuck']), 0)
			s$Z = o$Z; s$EZ = o$EZ; s$la = o$la; s$lb = o$lb; s$beta = o$beta; s$a_mat = o$a_mat; s$b_mat = o$b_mat
		}
		# every actor of every margin is visited once (a symmetric network's actors once)
		n_visit = if (mode == 'sym') s$n else s$n + s$nc
		expect_equal(unname(o$stats['level_calls'] + o$stats['skipped']), n_visit)
	}
})

test_that('the cached rank bounds equal the scan through a sequence of translations', {
	skip_on_cran()
	set.seed(11)
	chk = function(mode, nr, nc, N, K = 4L, nmiss = 5L) {
		sym = mode == 'sym'
		Z = array(rnorm(nr * nc * N), c(nr, nc, N))
		if (sym) for (t in 1:N) { m = Z[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Z[, , t] = m }
		Y = array(NA_real_, c(nr, nc, N))
		for (t in 1:N) { z = Z[, , t]; Y[, , t] = as.numeric(
			cut(z, quantile(z, seq(0, 1, length.out = K + 1)), include.lowest = TRUE, labels = FALSE)
		) }
		if (mode != 'bip') for (t in 1:N) diag(Y[, , t]) = NA
		Y[sample(length(Y), nmiss)] = NA
		if (sym) for (t in 1:N) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
		cx = .lame_ltm_context(Y, 'ordinal', mode)
		for (phase in if (sym) 'row' else c('row', 'col')) {
			row = phase == 'row'
			cache = .lame_ltm_ord_cache(cx, Z)
			for (i in sample.int(if (row) nr else nc)) {
				ref = .lame_ltm_ord_bounds(cx, Z, i, row)
				cb = .lame_ltm_ord_bounds_cached(cx, Z, cache, i, row)
				expect_equal(cb$bounds, ref)
				# the bounds are the interval in which the ranks hold with the other cells fixed
				lo = max(ref[1], -2); hi = min(ref[2], 2)
				if (!(lo < hi)) next
				d = runif(1, lo, hi)
				cache = .lame_ltm_ord_update(cache, cb$upd, d)
				if (sym) { Z[i, , ] = Z[i, , ] + d; Z[, i, ] = Z[, i, ] + d } else if (row) Z[i, , ] = Z[i, , ] + d else Z[, i, ] = Z[, i, ] + d
				expect_true(ltmt_ranks_hold(Y, Z))
			}
			cr = .lame_ltm_ord_cache(cx, Z)
			expect_equal(cache$gmax, cr$gmax); expect_equal(cache$gmin, cr$gmin)
		}
	}
	chk('dir', 8, 8, 3); chk(
		'dir', 6, 6, 1, K = 2L, nmiss = 0L
	); chk('sym', 9, 9, 3); chk('bip', 6, 9, 3); chk('bip', 10, 5, 2, K = 5L)
	# just beyond the bounds the ranks break
	s = ltmt_state('dir', 'ordinal', FALSE, FALSE)
	cx = .lame_ltm_context(s$Y, 'ordinal', 'dir')
	for (i in 1:3) {
		bd = .lame_ltm_ord_bounds(cx, s$Z, i, TRUE)
		for (d in c(bd[1] - 1e-6, bd[2] + 1e-6)) if (is.finite(d)) {
			Zs = s$Z; Zs[i, , ] = Zs[i, , ] + d
			expect_false(ltmt_ranks_hold(s$Y, Zs))
		}
		for (d in c(bd[1] + 1e-6, bd[2] - 1e-6)) if (is.finite(d)) {
			Zs = s$Z; Zs[i, , ] = Zs[i, , ] + d
			expect_true(ltmt_ranks_hold(s$Y, Zs))
		}
	}
})

test_that('the explicit-cutpoint bounds keep every value in its interval', {
	s = ltmt_state('dir', 'ordinal', TRUE, FALSE, K = 3L)
	alpha = c(0, .9)
	cx = .lame_ltm_context(s$Y, 'ordinal', 'dir', Y_int = .lame_ltm_codes(s$Y))
	af = c(-Inf, alpha, Inf)
	# a latent array inside the intervals
	set.seed(2)
	Z = s$Z
	for (i in which(is.finite(s$Y))) { k = cx$W[i]; Z[i] = if (k == 1) -abs(rnorm(1)) else if (k == 2) runif(1, 0, .9) else .9 + abs(rnorm(1)) }
	for (i in 1:4) {
		bd = .lame_ltm_ord_bounds(cx, Z, i, TRUE, alpha)
		expect_true(bd[1] < 0 && bd[2] > 0)
		ok = is.finite(s$Y[i, , ])
		for (d in c(bd[1] + 1e-6, bd[2] - 1e-6)) {
			Zs = Z; Zs[i, , ] = Zs[i, , ] + d
			w = cx$W[i, , ]; zi = Zs[i, , ]
			expect_true(all(zi[ok] > af[w[ok]] & zi[ok] <= af[w[ok] + 1L]))
		}
	}
})

test_that('the location move keeps every residual and the ranks, and the strip keeps the cutpoint intervals', {
	for (mode in c('dir', 'sym', 'bip')) for (dyn in c(FALSE, TRUE)) {
		s = ltmt_state(mode, 'ordinal', FALSE, dyn)
		cx = .lame_ltm_context(s$Y, 'ordinal', mode)
		lev = c(TRUE, mode != 'sym')
		for (rep in 1:5) {
			o = .lame_ltm_location(
				cx,
				s$Z,
				s$EZ,
				s$la,
				s$lb,
				s$a_mat,
				s$b_mat,
				s$beta,
				ltmt_sab(mode),
				NULL,
				centred = FALSE,
				lev = lev
			)
			if (mode == 'sym') { o$lb = o$la; o$b_mat = o$a_mat }
			expect_lt(max(abs((o$Z - o$EZ) - (s$Z - s$EZ))), 1e-12)
			expect_lt(max(abs(s$build(o$la, o$lb, o$beta, o$a_mat, o$b_mat) - o$EZ)), 1e-12)
			expect_equal(o$Z, s$Z + o$shift)
			expect_true(ltmt_ranks_hold(s$Y, o$Z))
			s$Z = o$Z; s$EZ = o$EZ; s$la = o$la; s$lb = o$lb; s$a_mat = o$a_mat; s$b_mat = o$b_mat
		}
	}
	# with an intercept the intercept carries the shift; with cutpoints the total shift stays in the strip
	s = ltmt_state('dir', 'ordinal', TRUE, FALSE, K = 3L)
	alpha = c(0, .9)
	cx = .lame_ltm_context(s$Y, 'ordinal', 'dir', Y_int = .lame_ltm_codes(s$Y))
	af = c(-Inf, alpha, Inf)
	Z = s$Z; set.seed(5)
	for (i in which(is.finite(s$Y))) { k = cx$W[i]; Z[i] = if (k == 1) -abs(rnorm(1)) else if (k == 2) runif(1, 0, .9) else .9 + abs(rnorm(1)) }
	EZ = Z - (s$Z - s$EZ)
	for (rep in 1:20) {
		o = .lame_ltm_location(
			cx,
			Z,
			EZ,
			s$la,
			s$lb,
			NULL,
			NULL,
			s$beta,
			ltmt_sab('dir'),
			diag(2) / 4,
			centred = TRUE,
			alpha = alpha
		)
		ok = is.finite(s$Y)
		expect_true(all(o$Z[ok] > af[cx$W[ok]] & o$Z[ok] <= af[cx$W[ok] + 1L]))
		expect_equal(o$beta[1], s$beta[1] + o$shift); expect_equal(o$Z, Z + o$shift)
		Z = o$Z; EZ = o$EZ; s$beta = o$beta
	}
})

test_that('the coefficient move keeps every residual and rebuilds the predictor', {
	for (mode in c('dir', 'sym', 'bip')) for (centred in c(TRUE, FALSE)) {
		s = ltmt_state(mode, 'poisson', centred, FALSE)
		cx = .lame_ltm_context(s$Y, 'poisson', mode, offset = c(0, .2, -.1))
		Xk = if (centred) list(array(1, dim(s$Z)), s$x) else list(s$x)
		for (rep in 1:5) {
			o = .lame_ltm_coef(cx, s$Z, s$EZ, s$beta, Xk, diag(length(s$beta)) / 4)
			expect_lt(max(abs((o$Z - o$EZ) - (s$Z - s$EZ))), 1e-12)
			expect_lt(max(abs(s$build(s$la, s$lb, o$beta, NULL, NULL) - o$EZ)), 1e-12)
			expect_equal(unname(o$stats['coef_calls']), length(s$beta)); expect_equal(unname(o$stats['coef_stuck']), 0)
			s$Z = o$Z; s$EZ = o$EZ; s$beta = o$beta
		}
	}
})

test_that('a position move keeps every residual and rebuilds the predictor', {
	set.seed(31)
	for (mode in c('dir', 'sym', 'bip')) for (R in c(1L, 2L)) {
		s = ltmt_state(mode, 'poisson', FALSE, FALSE)
		nr = s$n; nc = s$nc; N = s$N
		U = matrix(rnorm(nr * R), nr, R); G = matrix(rnorm(R * R, 0, .5), R, R)
		V = if (mode == 'bip') matrix(rnorm(nc * R), nc, R) else matrix(rnorm(nr * R), nr, R)
		Ld = if (mode == 'sym') c(.8, -.5)[seq_len(R)] else NULL
		if (mode == 'sym') V = U %*% diag(Ld, R)
		mult = function(U, V) if (mode == 'bip') U %*% G %*% t(V) else U %*% t(V)
		base = s$EZ
		EZ = base; Z = s$Z
		for (t in 1:N) { EZ[, , t] = EZ[, , t] + mult(U, V); Z[, , t] = Z[, , t] + mult(U, V) }
		cx = .lame_ltm_context(s$Y, 'poisson', mode, offset = c(0, .2, -.1))
		for (rep in 1:4) {
			o = .lame_ltm_uv_visit(
				cx,
				Z,
				EZ,
				U,
				if (mode == 'sym') NULL else V,
				if (mode == 'bip') G else NULL,
				Ld,
				Suv = if (mode == 'dir') diag(2 * R) + .3 else NULL,
				tau = if (mode == 'sym') rep(1.2, R) else NULL
			)
			V1 = if (mode == 'sym') o$U %*% diag(Ld, R) else o$V
			expect_lt(max(abs((o$Z - o$EZ) - (Z - EZ))), 1e-10)
			rebuilt = base; for (t in 1:N) rebuilt[, , t] = rebuilt[, , t] + mult(o$U, V1)
			expect_lt(max(abs(rebuilt - o$EZ)), 1e-10)
			expect_equal(
				unname(o$stats['uv_calls']),
				if (mode == 'sym') nr else nr + nc
			); expect_equal(unname(o$stats['uv_stuck']), 0)
			expect_false(isTRUE(all.equal(o$U, U)))
			Z = o$Z; EZ = o$EZ; U = o$U; V = V1
		}
	}
	# before the second half of burn-in the host runs a diffuse prior and owns directed and symmetric positions
	s = ltmt_state('dir', 'poisson', FALSE, FALSE)
	cx = .lame_ltm_context(s$Y, 'poisson', 'dir')
	o = .lame_ltm_uv_visit(cx, s$Z, s$EZ, matrix(rnorm(6), 6), matrix(rnorm(6), 6), shrink = FALSE)
	expect_equal(unname(o$stats['uv_calls']), 0)
})

test_that('the scale move carries the standardized contrasts of every dyad and keeps (s2, rho) in range', {
	skip_on_cran()
	set.seed(21)
	for (mode in c('dir', 'sym', 'bip')) for (two in unique(c(FALSE, mode == 'dir'))) {
		s = ltmt_state(mode, 'poisson', TRUE, FALSE)
		cx = .lame_ltm_context(s$Y, 'poisson', mode, offset = c(0, .2, -.1))
		s2 = .4; rho = if (two) .3 else 0
		prior = c(shape0 = .5, rate0 = .5)
		lp = .lame_rho_prior_log
		for (rep in 1:8) {
			o = .lame_ltm_scale(
				cx,
				s$Z,
				s$EZ,
				s2,
				rho,
				prior,
				two = two,
				rho_logprior = lp,
				tries = 3L,
				diag_in_target = FALSE
			)
			expect_gt(o$s2, 0); expect_lt(abs(o$rho), 1)
			e0 = s$Z - s$EZ; e1 = o$Z - s$EZ
			if (mode == 'bip') {
				expect_equal(e1 / sqrt(o$s2), e0 / sqrt(s2), tolerance = 1e-10)
			} else if (two) {
				p0 = (e0[cx$pu] + e0[cx$pl]) / 2; q0 = (e0[cx$pu] - e0[cx$pl]) / 2
				p1 = (e1[cx$pu] + e1[cx$pl]) / 2; q1 = (e1[cx$pu] - e1[cx$pl]) / 2
				expect_equal(p1 / sqrt(o$s2 * (1 + o$rho)), p0 / sqrt(s2 * (1 + rho)), tolerance = 1e-10)
				expect_equal(q1 / sqrt(o$s2 * (1 - o$rho)), q0 / sqrt(s2 * (1 - rho)), tolerance = 1e-10)
			} else {
				expect_equal(e1[cx$pu] / sqrt(o$s2), e0[cx$pu] / sqrt(s2), tolerance = 1e-10)
				expect_equal(o$Z[cx$pl][mode == 'sym' | TRUE], if (mode == 'sym') o$Z[cx$pu] else o$Z[cx$pl])
			}
			# the diagonal is left to the host
			expect_equal(o$Z[cx$pd], s$Z[cx$pd])
			s$Z = o$Z; s2 = o$s2; rho = o$rho
		}
	}
	# with every cell missing the data factor is flat, so the move must leave the prior of (s2, rho) invariant
	set.seed(22)
	Y = array(NA_real_, c(4, 4, 2)); cx = .lame_ltm_context(Y, 'poisson', 'dir')
	EZ = array(rnorm(32), c(4, 4, 2))
	draw = function() {
		s2 = 1 / rgamma(1, .5, .5); pw = .lame_rho_prior_power(); rho = 2 * rbeta(1, pw + 1, pw + 1) - 1
		e = array(0, c(4, 4, 2))
		z = matrix(rnorm(2 * length(cx$pu)), ncol = 2) %*% chol(s2 * matrix(c(1, rho, rho, 1), 2))
		e[cx$pu] = z[, 1]; e[cx$pl] = z[, 2]
		list(s2 = s2, rho = rho, Z = EZ + e)
	}
	stat = function(s2, rho, Z, E) c(log(s2), rho, rho^2, log(s2) * rho, ((Z - E)[cx$pu][1])^2 / s2)
	d = replicate(
		3000,
		{ s = draw(); o = .lame_ltm_scale(
			cx,
			s$Z,
			EZ,
			s$s2,
			s$rho,
			c(shape0 = .5, rate0 = .5),
			two = TRUE,
			rho_logprior = .lame_rho_prior_log,
			tries = 4L
		); c(stat(o$s2, o$rho, o$Z, EZ) - stat(s$s2, s$rho, s$Z, EZ)) }
	)
	z = rowMeans(d) / (apply(d, 1, sd) / sqrt(ncol(d)))
	expect_lt(max(abs(z)), 4.5)
})

ltmt_panel = function(mode, family, nt = 3L) {
	set.seed(5)
	nr = if (mode == 'bip') 7L else 8L; nc = if (mode == 'bip') 9L else 8L
	rn = sprintf('r%02d', 1:nr); cn = if (mode == 'bip') sprintf('c%02d', 1:nc) else rn
	Y = lapply(1:nt, function(t) {
		y = matrix(
			if (family == 'poisson') rpois(nr * nc, .4) else sample(1:3, nr * nc, TRUE),
			nr,
			nc,
			dimnames = list(rn, cn)
		)
		if (family == 'poisson') { y[1:2, ] = 0; y[, 1:2] = 0 } else { y[1:2, ] = 1; y[, 1:2] = 1 }
		if (mode == 'sym') y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (mode != 'bip') diag(y) = NA
		y
	})
	X = lapply(
		1:nt,
		function(t) { x = matrix(
			rnorm(nr * nc),
			nr,
			nc,
			dimnames = list(rn, cn)
		); if (mode == 'sym') x[lower.tri(x)] = t(x)[lower.tri(x)]; array(
			x,
			c(nr, nc, 1),
			dimnames = list(rn, cn, 'x')
		) }
	)
	list(Y = Y, X = X)
}
ltmt_run = function(p, mode, family, prior = list(), ...) suppressMessages(suppressWarnings(lame(
	p$Y,
	Xdyad = p$X,
	mode = if (mode == 'bip') 'bipartite' else 'unipartite',
	symmetric = mode == 'sym',
	family = family,
	R = 0L,
	R_row = if (mode == 'bip') 0L else NULL,
	R_col = if (mode == 'bip') 0L else NULL,
	burn = 10,
	nscan = 30,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = modifyList(list(actor_move_every = 2), prior),
	...
)))
# the update records its work, a visit interval beyond the run matches the switch off, and a run with visits differs
ltmt_check = function(mode, family, dyn) {
	p = ltmt_panel(mode, family)
	on = ltmt_run(p, mode, family, dynamic_ab = dyn)
	am = on$actor_move
	expect_equal(unname(am['visits']), 20)
	expect_equal(unname(am['level_calls'] + am['skipped']), 20 * (if (mode == 'sym') 8 else 16))
	expect_equal(unname(am['stuck']), 0)
	if (family == 'poisson') {
		expect_equal(unname(am['coef_calls']), 20 * 2)
		expect_gt(unname(am['scale_attempts']), 0)
		expect_equal(unname(am['rho_attempts'] > 0), mode == 'dir')
	} else expect_equal(unname(am['location_calls']), 20)
	expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
	off = ltmt_run(p, mode, family, prior = list(actor_move = FALSE), dynamic_ab = dyn)
	expect_null(off$actor_move)
	never = ltmt_run(p, mode, family, prior = list(actor_move_every = 1000), dynamic_ab = dyn)
	expect_null(never$actor_move)
	expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
	expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
}

test_that('lame() records the work of the translations and honours the switches', {
	skip_on_cran()
	ltmt_check('dir', 'poisson', FALSE)
	ltmt_check('bip', 'ordinal', TRUE)
	p = ltmt_panel('dir', 'poisson')
	# a fit without additive effects moves the Poisson coefficients and scale
	noab = ltmt_run(p, 'dir', 'poisson', rvar = FALSE, cvar = FALSE)
	expect_equal(unname(noab$actor_move['level_calls']), 0); expect_gt(unname(noab$actor_move['coef_calls']), 0)
	# the option is validated
	expect_error(
		ltmt_run(p, 'dir', 'poisson', prior = list(actor_move = NA)),
		'actor_move must be TRUE or FALSE',
		fixed = TRUE
	)
})

test_that('lame() runs the translations in every mode, family and effect type', {
	skip_on_cran()
	for (family in c('poisson', 'ordinal')) for (mode in c('dir', 'sym', 'bip')) for (dyn in c(FALSE, TRUE)) ltmt_check(mode, family, dyn)
	p = ltmt_panel('dir', 'ordinal')
	ex = ltmt_run(p, 'dir', 'ordinal', ordinal_cutpoints = 'explicit')
	expect_equal(unname(ex$actor_move['visits']), 20); expect_equal(unname(ex$actor_move['location_calls']), 20)
	expect_true(all(is.finite(ex$BETA)))
	# static positions move with the latent array too (unipartite factor prior from mid burn-in)
	for (mode in c('dir', 'sym', 'bip')) {
		p = ltmt_panel(mode, 'poisson')
		pos = suppressMessages(suppressWarnings(lame(
			p$Y,
			Xdyad = p$X,
			mode = if (mode == 'bip') 'bipartite' else 'unipartite',
			symmetric = mode == 'sym',
			family = 'poisson',
			R = if (mode == 'bip') 0L else 1L,
			R_row = if (mode == 'bip') 1L else NULL,
			R_col = if (mode == 'bip') 1L else NULL,
			burn = 20,
			nscan = 30,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = list(actor_move_every = 2)
		)))
		expect_gt(unname(pos$actor_move['uv_calls']), 0); expect_equal(unname(pos$actor_move['uv_stuck']), 0)
		expect_true(all(is.finite(pos$BETA)) && all(is.finite(pos$VC)))
	}
	db = ltmt_run(ltmt_panel('dir', 'poisson'), 'dir', 'poisson', dynamic_beta = TRUE)
	expect_equal(unname(db$actor_move['visits']), 20); expect_equal(unname(db$actor_move['coef_calls']), 0)
	expect_true(all(is.finite(db$BETA)))
})

test_that('ame() runs the bipartite translations and records their work', {
	skip_on_cran()
	for (family in c('poisson', 'ordinal')) {
		p = ltmt_panel('bip', family, nt = 1L)
		run = function(prior) suppressMessages(suppressWarnings(ame(
			p$Y[[1]],
			Xdyad = p$X[[1]],
			mode = 'bipartite',
			family = family,
			R_row = 0L,
			R_col = 0L,
			burn = 10,
			nscan = 40,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior
		)))
		on = run(list(actor_move_every = 2)); off = run(list(actor_move = FALSE)); never = run(list(actor_move_every = 1000))
		expect_equal(unname(on$actor_move['visits']), 25); expect_equal(unname(on$actor_move['stuck']), 0)
		expect_equal(unname(on$actor_move['level_calls'] + on$actor_move['skipped']), 25 * 16)
		expect_null(off$actor_move); expect_null(never$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	}
	# positions move with the latent array in a Poisson fit
	p = ltmt_panel('bip', 'poisson', nt = 1L)
	pos = suppressMessages(suppressWarnings(ame(
		p$Y[[1]],
		Xdyad = p$X[[1]],
		mode = 'bipartite',
		family = 'poisson',
		R_row = 1L,
		R_col = 1L,
		burn = 10,
		nscan = 40,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = list(actor_move_every = 2)
	)))
	expect_gt(unname(pos$actor_move['uv_calls']), 0); expect_equal(unname(pos$actor_move['uv_stuck']), 0)
})

test_that('the kernels leave the exact posterior of a tiny panel invariant', {
	skip_on_cran()
	cases = list(
		list(family = 'poisson', mode = 'dir', comp = TRUE, free = FALSE, kernel = 'visit', ndraw = 2e5),
		list(family = 'poisson', mode = 'sym', comp = FALSE, free = TRUE, kernel = 'visit', ndraw = 2e5),
		list(family = 'poisson', mode = 'bip', comp = TRUE, free = TRUE, kernel = 'coef', ndraw = 2e5),
		list(family = 'poisson', mode = 'dir', comp = TRUE, free = FALSE, kernel = 'scale', ndraw = 3e5),
		list(family = 'poisson', mode = 'bip', comp = TRUE, free = FALSE, kernel = 'scale', ndraw = 3e5),
		list(family = 'poisson', mode = 'dir', comp = TRUE, free = FALSE, kernel = 'uv', ndraw = 3e5),
		list(family = 'poisson', mode = 'sym', comp = TRUE, free = FALSE, kernel = 'uv', ndraw = 3e5),
		list(family = 'ordinal', mode = 'dir', comp = FALSE, free = TRUE, kernel = 'visit', ndraw = 1e6),
		list(family = 'ordinal', mode = 'bip', comp = TRUE, free = FALSE, kernel = 'visit', ndraw = 1e6),
		list(family = 'ordinal', mode = 'sym', comp = FALSE, free = TRUE, kernel = 'location', ndraw = 1e6)
	)
	for (k in seq_along(cases)) {
		r = do.call(ltmt_exact, c(cases[[k]], list(seed = 100L + k)))
		expect_gt(r$ess, 400)
		expect_lt(r$maxz, 4.6)
	}
})

test_that('the exact-posterior check detects a kernel that drops its data factor or mis-states its prior', {
	skip_on_cran()
	phase = .lame_ltm_pois_phase
	local_mocked_bindings(
		.lame_ltm_pois_phase = function(cx, Z, row) { r = phase(cx, Z, row); r$S = r$S * 0; r$ys = r$ys * 0; r },
		.package = 'lame'
	)
	expect_gt(ltmt_exact('poisson', 'dir', TRUE, FALSE, 'visit', ndraw = 2e5, seed = 9L)$maxz, 8)
})
test_that('the exact-posterior check detects a mis-stated level prior and a mis-stated scale prior', {
	skip_on_cran()
	prior = .lame_actor_move_level_prior
	local_mocked_bindings(
		.lame_actor_move_level_prior = function(...) { r = prior(...); c(a1 = r[['a1']], a2 = r[['a2']] * 1.2) },
		.package = 'lame'
	)
	expect_gt(ltmt_exact('poisson', 'dir', TRUE, FALSE, 'visit', ndraw = 2e5, seed = 10L)$maxz, 5)
})
