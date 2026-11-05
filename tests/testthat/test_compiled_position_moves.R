# the compiled visits of the position moves give the draws of their R references from the same seed, to rounding

cpm_same = function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10))

test_that('the compiled visit of the static positions of a cbin / frn panel gives the draws of the R visit', {
	skip_on_cran()
	for (fam in c('cbin', 'frn')) {
		set.seed(1)
		n = 11L; N = 3L; R = 2L; cap = 3L
		Y = array(0, c(n, n, N))
		for (t in 1:N) for (i in 1:n) { k = sample(0:3, 1L); if (k > 0L) Y[i, sample(setdiff(1:n, i), k), t] = if (fam == 'frn') seq_len(k) else 1 }
		for (t in 1:N) diag(Y[, , t]) = NA
		YL = if (fam == 'frn') lapply(
			1:N,
			function(t) t(sapply(1:n, function(i) match(seq_len(cap), Y[i, , t])))
		) else NULL
		U0 = matrix(rnorm(n * R, 0, .7), n, R); V0 = matrix(rnorm(n * R, 0, .7), n, R)
		Z0 = array(0, c(n, n, N))
		for (t in 1:N) {
			EZ = U0 %*% t(V0) - 1
			od = rowSums(Y[, , t] > 0, na.rm = TRUE)
			Z0[, , t] = if (fam == 'cbin') rZ_cbin_fc(matrix(0, n, n), EZ, 0, Y[, , t], cap, od) else rZ_frn_fc(
				matrix(0, n, n),
				EZ,
				0,
				Y[, , t],
				YL[[t]],
				cap,
				od
			)
		}
		ctx = lame:::.lame_spm_context(Y, fam, YL, cap)
		run = function(engine) {
			withr::local_options(lame.spm_engine = engine)
			set.seed(99)
			Z = Z0; U = U0; V = V0
			out = vector('list', 40L)
			for (s in seq_along(out)) {
				o = lame:::.lame_spm_visit(ctx, Z, U, V, shrink = TRUE, uv_prior = list(df = 4, scale = 4 * diag(4)))
				Z = o$Z; U = o$U; V = o$V
				out[[s]] = list(Z = Z, U = U, V = V, st = o$stats)
			}
			out
		}
		a = run('r'); b = run('cpp')
		expect_true(all(mapply(cpm_same, a, b)))
		expect_gt(a[[40L]]$st[['spm_calls']], 0)
	}
})

test_that('the compiled visit of the coordinate paths of a cbin / frn panel gives the draws of the R visit, bipartite and directed', {
	skip_on_cran()
	for (mode in c('bip', 'dir')) for (fam in c('cbin', 'frn')) {
		set.seed(2)
		nr = 9L; nc = if (mode == 'bip') 7L else 9L; N = 3L; cap = 2L; RA = 2L; RB = if (mode == 'bip') 1L else 2L
		Y = array(0, c(nr, nc, N))
		for (t in 1:N) for (i in 1:nr) { k = sample(0:2, 1L); if (k > 0L) Y[i, sample(
			setdiff(1:nc, if (mode == 'dir') i else 0L),
			k
		), t] = if (fam == 'frn') seq_len(k) else 1 }
		if (mode == 'dir') for (t in 1:N) diag(Y[, , t]) = NA
		YL = if (fam == 'frn') lapply(
			1:N,
			function(t) t(sapply(1:nr, function(i) match(seq_len(cap), Y[i, , t])))
		) else NULL
		G = if (mode == 'bip') matrix(rnorm(RA * RB, 1, .3), RA, RB) else NULL
		Uc = array(rnorm(nr * RA * N, 0, .6), c(nr, RA, N)); Vc = array(rnorm(nc * RB * N, 0, .6), c(nc, RB, N))
		Z0 = array(0, c(nr, nc, N))
		for (t in 1:N) {
			Ut = matrix(Uc[, , t], nr, RA); Vt = matrix(Vc[, , t], nc, RB)
			EZ = (if (mode == 'bip') Ut %*% G %*% t(Vt) else Ut %*% t(Vt)) - 1
			od = rowSums(Y[, , t] > 0, na.rm = TRUE)
			Z0[, , t] = if (mode == 'bip') {
				if (fam == 'cbin') rZ_cbin_bip_fc(matrix(0, nr, nc), EZ, Y[, , t], cap, od) else rZ_frn_bip_fc(
					matrix(0, nr, nc),
					EZ,
					Y[, , t],
					YL[[t]],
					cap,
					od
				)
			} else {
				if (fam == 'cbin') rZ_cbin_fc(matrix(0, nr, nc), EZ, 0, Y[, , t], cap, od) else rZ_frn_fc(
					matrix(0, nr, nc),
					EZ,
					0,
					Y[, , t],
					YL[[t]],
					cap,
					od
				)
			}
		}
		ctx = lame:::.lame_npm_context(Y, fam, YL, cap, mode == 'bip')
		run = function(engine) {
			withr::local_options(lame.npm_engine = engine)
			set.seed(98)
			Z = Z0; U = Uc; V = Vc
			out = vector('list', 30L)
			for (s in seq_along(out)) {
				o = lame:::.lame_npm_visit(Z, Y, fam, YL, cap, U, V, G, .8, .6, ctx = ctx)
				Z = o$Z; U = o$U; V = o$V
				out[[s]] = list(Z = Z, U = U, V = V, st = o$stats)
			}
			out
		}
		a = run('r'); b = run('cpp')
		expect_true(all(mapply(cpm_same, a, b)))
	}
})

test_that('the compiled tempered transition gives the draws of the R reference', {
	skip_on_cran()
	set.seed(3)
	n = 8L; N = 3L; R = 2L; rho = .85; off = c(0, .3, -.2); L = c(.8, -.5)
	U = lame:::.lame_uv_unit_paths(array(rnorm(n * R * N), c(n, R, N)), rho)
	Z = Yc = array(0, c(n, n, N))
	for (t in 1:N) {
		Ut = matrix(U[, , t], n, R); z = Ut %*% diag(L) %*% t(Ut) - 1.5 + matrix(rnorm(n * n, 0, .4), n)
		z[lower.tri(z)] = t(z)[lower.tri(z)]; Z[, , t] = z
		y = matrix(rpois(n * n, exp(z + off[t])), n); y[lower.tri(y)] = t(y)[lower.tri(y)]; diag(y) = NA; Yc[, , t] = y
	}
	Yc[1, 2, 1] = Yc[2, 1, 1] = NA; Yc[3, , 2] = Yc[, 3, 2] = NA
	yy = Yc; for (t in 1:N) yy[, , t][lower.tri(yy[, , t], TRUE)] = NA
	dt = lame:::.lame_temper_data(lame:::.lame_ltm_context(yy, 'poisson', 'sym', offset = off), n, N)
	bs = seq(1, .7, length.out = 6)
	run = function(engine, cycles) {
		withr::local_options(lame.temper_engine = engine)
		set.seed(97)
		Zs = Z; Us = U; Ls = L
		out = vector('list', 4L)
		for (s in seq_along(out)) {
			o = lame:::.lame_uv_sym_temper(Us, Ls, Zs, dt, sample.int(2L, 1L), rho, 1, bs, cycles = cycles)
			if (o$accepted) { Zs = o$Z; Us = o$U; Ls = o$L }
			out[[s]] = list(Z = Zs, U = Us, L = Ls, acc = o$accepted, la = o$logacc)
		}
		out
	}
	for (cycles in 1:2) {
		a = run('r', cycles); b = run('cpp', cycles)
		expect_true(all(mapply(cpm_same, a, b)))
	}
})

test_that('the engine option is validated', {
	withr::local_options(lame.spm_engine = 'fortran')
	expect_error(lame:::.lame_move_engine('lame.spm_engine'), 'cpp')
})
