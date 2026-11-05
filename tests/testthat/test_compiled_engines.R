# compiled kernels of the Poisson translation moves, the symmetric actor-period move and the nomination feasibility check
# match their R versions to rounding (options lame.ltm_engine = "r", lame.pab_engine, lame.nomination_engine) on the same inputs and seeds

ce_with = function(opt, val, code) {
	old = getOption(opt); options(structure(list(val), names = opt)); on.exit(
		options(structure(list(old), names = opt))
	)
	force(code)
}

ce_panel = function(mode, seed = 5L, n = 7L, N = 3L, nb = 9L) {
	set.seed(seed)
	bip = mode == 'bip'; sym = mode == 'sym'; nr = n; nc = if (bip) nb else n
	mk = function(x) { if (sym) for (t in seq_len(dim(x)[3L])) { m = x[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; x[, , t] = m }; x }
	la = rnorm(nr, 0, .5); lb = if (sym) la else rnorm(nc, 0, .5)
	EZ = mk(array(-0.5 + rnorm(nr * nc * N, 0, .2), c(nr, nc, N)))
	for (t in 1:N) EZ[, , t] = EZ[, , t] + outer(la, lb, '+')
	Z = EZ + mk(array(rnorm(nr * nc * N, 0, .3), c(nr, nc, N)))
	Y = mk(array(rpois(nr * nc * N, exp(Z)), c(nr, nc, N)))
	if (!bip) for (t in 1:N) diag(Y[, , t]) = NA
	Y[sample(length(Y), 5)] = NA
	Y = mk(Y)
	list(Z = Z, EZ = EZ, Y = Y, la = la, lb = lb, nr = nr, nc = nc, N = N)
}

test_that('the compiled level translations of a Poisson fit give the draws of the R loop', {
	for (mode in c('dir', 'sym', 'bip')) {
		s = ce_panel(mode)
		cx = lame:::.lame_ltm_context(s$Y, 'poisson', mode, offset = c(0, .2, -.1))
		Sab = if (mode == 'bip') diag(c(1.1, .8)) else if (mode == 'sym') matrix(1.1, 2, 2) else matrix(c(1, .3, .3, 1.2), 2)
		out = lapply(c('r', 'cpp'), function(eng) ce_with('lame.ltm_engine', eng, {
			set.seed(21); st = s
			for (rep in 1:5) {
				o = lame:::.lame_ltm_visit(cx, st$Z, st$EZ, st$la, st$lb, NULL, NULL, 0.5, Sab, NULL, centred = FALSE)
				st$Z = o$Z; st$EZ = o$EZ; st$la = o$la; st$lb = o$lb
			}
			c(o, list(u = runif(1)))
		}))
		expect_equal(out[[1]], out[[2]], tolerance = 1e-10)
	}
})

test_that('the compiled position translations of a Poisson fit give the draws of the R loop', {
	for (mode in c('dir', 'sym', 'bip')) for (R in c(1L, 2L)) {
		s = ce_panel(mode, seed = 8L)
		set.seed(2)
		U = matrix(rnorm(s$nr * R), s$nr, R); G = matrix(rnorm(R * R, 0, .5), R, R)
		V = if (mode == 'bip') matrix(rnorm(s$nc * R), s$nc, R) else matrix(rnorm(s$nr * R), s$nr, R)
		Ld = if (mode == 'sym') c(.8, -.5)[seq_len(R)] else NULL
		if (mode == 'sym') V = U %*% diag(Ld, R)
		cx = lame:::.lame_ltm_context(s$Y, 'poisson', mode, offset = c(0, .2, -.1))
		out = lapply(c('r', 'cpp'), function(eng) ce_with('lame.ltm_engine', eng, {
			set.seed(22); Z = s$Z; EZ = s$EZ; Ui = U; Vi = V
			for (rep in 1:5) {
				o = lame:::.lame_ltm_uv_visit(
					cx,
					Z,
					EZ,
					Ui,
					if (mode == 'sym') NULL else Vi,
					if (mode == 'bip') G else NULL,
					Ld,
					Suv = if (mode == 'dir') diag(2 * R) + .3 else NULL,
					tau = if (mode == 'sym') rep(1.2, R) else NULL
				)
				Z = o$Z; EZ = o$EZ; Ui = o$U; if (mode != 'sym') Vi = o$V
			}
			c(o, list(u = runif(1)))
		}))
		expect_equal(out[[1]], out[[2]], tolerance = 1e-10)
	}
})

test_that('the compiled coefficient translation gives the draws of the R step', {
	s = ce_panel('dir', seed = 9L)
	cx = lame:::.lame_ltm_context(s$Y, 'poisson', 'dir', offset = c(0, .2, -.1))
	cxp = list(ltm = cx, mode = 'dir')
	Xk = list(array(1, dim(s$Z)), array(rnorm(length(s$Z)), dim(s$Z)))
	beta = c(-0.4, 0.3); iV0 = diag(2) / 4
	out = lapply(c('r', 'cpp'), function(eng) ce_with('lame.ltm_engine', eng, {
		set.seed(23); Z = s$Z; EZ = s$EZ; b = beta
		for (rep in 1:6) { o = lame:::.lame_ppm_static(cxp, Z, EZ, b, Xk, iV0); Z = o$Z; EZ = o$EZ; b = o$beta }
		c(o, list(u = runif(1)))
	}))
	expect_equal(out[[1]], out[[2]], tolerance = 1e-10)
})

test_that('the compiled symmetric actor-period loop gives the draws of the R loop', {
	set.seed(1); n = 12
	for (centred in c(TRUE, FALSE)) {
		Zt = matrix(rnorm(n * n, -1, .7), n, n); Zt = (Zt + t(Zt)) / 2; Et = Zt + 0.1
		y = matrix(rpois(n * n, 1), n, n); y = pmax(y, t(y)); y[sample(n * n, 16)] = NA; y = ifelse(is.na(t(y)), NA, y); diag(y) = NA
		eff = rnorm(n); eff[3] = NA; m = rnorm(n, 0, .2); v = rep(0.3, n)
		out = lapply(c('r', 'cpp'), function(eng) ce_with('lame.pab_engine', eng, {
			set.seed(11); Z = Zt; E = Et; e = eff; ev = 0L
			for (s in 1:60) { r = lame:::.lame_pab_sym_period(
				Z,
				E,
				y,
				0.2,
				e,
				m,
				v,
				centred
			); Z = r$Z; E = r$E; e = r$eff; ev = ev + r$evals }
			list(Z = Z, E = E, e = e, ev = ev, u = runif(1))
		}))
		expect_equal(out[[1]], out[[2]], tolerance = 1e-10)
	}
})

test_that('the compiled feasibility check finds the rows of the R check', {
	set.seed(3)
	for (r in 1:400) {
		bip = r %% 2 == 0; frn = r %% 3 == 0; n = sample(4:8, 1); nc = if (bip) sample(4:8, 1) else n; cap = sample(1:3, 1)
		Z = matrix(rnorm(n * nc, -0.3, 1), n, nc); Y = matrix(0, n, nc)
		for (i in 1:n) { k = sample(0:cap, 1); if (k) { cs = sample(
			setdiff(seq_len(nc), if (!bip) i),
			min(k, nc - 1)
		); Y[i, cs] = if (frn) seq_along(cs) else 1 } }
		if (!bip) diag(Y) = NA
		if (r %% 5 == 0) Y[sample(length(Y), 5)] = NA
		if (r %% 7 == 0) Y[sample(n, 1), ] = NA
		YL = NULL
		if (frn) { YL = matrix(NA_integer_, n, cap); for (i in 1:n) for (q in 1:cap) { j = which(Y[i, ] == q); if (length(j)) YL[i, q] = j[1] } }
		odobs = rowSums(!is.na(Y) & Y > 0)
		if (r %% 4 == 0) odobs = odobs + rowSums(lame:::.lame_nom_hidden_cells(Y, bip) & Z > 0)
		a = ce_with('lame.nomination_engine', 'r', lame:::.lame_nom_feasible(Z, Y, YL, rep(cap, n), odobs, bip))
		b = ce_with('lame.nomination_engine', 'cpp', lame:::.lame_nom_feasible(Z, Y, YL, rep(cap, n), odobs, bip))
		expect_identical(a, b)
	}
})
