# cbin and frn coordinate paths under dynamic_uv, drawn with the latent rows and columns they bound (R/nomination_path_move.R)

# stationarity on a 3 x 3 two-period panel with one coordinate: AR(1) prior paths and latent values kept when an
# oracle of the family's data rule maps them to the fixed data; one visit must leave every mean unchanged and rows feasible
npm_exact = function(fam, mode, mutant = 'none', target = 2500L, seed = 71L) {
	set.seed(seed)
	n = 3L; N = 2L; rho = .8; sg = .9; Gm = matrix(1.3, 1, 1); bip = mode == 'bip'
	cap = if (fam == 'cbin') 1L else 2L
	Y = array(NA_real_, c(n, n, N))
	setrow = function(t, i, v) { cols = if (bip) 1:n else setdiff(1:n, i); Y[i, cols, t] <<- v }
	if (fam == 'cbin') {
		if (bip) { setrow(1, 1, c(0, 1, 0)); setrow(1, 2, c(0, 0, 0)); setrow(2, 1, c(0, 1, 0)); setrow(2, 2, c(0, 0, 1)) }
		else { setrow(1, 1, c(0, 1)); setrow(1, 2, c(0, 0)); setrow(2, 1, c(0, 1)); setrow(2, 2, c(0, 1)) }
	} else {
		if (bip) { setrow(1, 1, c(2, 1, 0)); setrow(1, 2, c(0, 0, 0)); setrow(2, 1, c(0, 1, 0)); setrow(2, 2, c(1, 0, 2)) }
		else { setrow(1, 1, c(2, 1)); setrow(1, 2, c(0, 0)); setrow(2, 1, c(0, 1)); setrow(2, 2, c(1, 2)) }
	}
	YL = if (fam == 'frn') lapply(1:N, function(t) t(sapply(1:n, function(i) match(seq_len(cap), Y[i, , t])))) else NULL
	oracle_row = function(z, y) {
		k = ncol(z); ok = rep(TRUE, nrow(z))
		if (fam == 'cbin') {
			for (j in 1:k) { mx = z[, j] > 0; for (l in setdiff(1:k, j)) mx = mx & z[, j] > z[, l]; ok = ok & (if (y[j] == 1) mx else !mx) }
		} else {
			nom = lapply(
				1:k,
				function(j) { above = Reduce(
					`+`,
					lapply(setdiff(1:k, j), function(l) (z[, l] > z[, j] & z[, l] > 0) * 1)
				); z[, j] > 0 & above < cap }
			)
			for (j in 1:k) { lower = Reduce(
				`+`,
				lapply(setdiff(1:k, j), function(l) (nom[[l]] & z[, l] < z[, j]) * 1)
			); ok = ok & ifelse(nom[[j]], 1 + lower, 0) == y[j] }
		}
		ok
	}
	feas = function(Zc) {
		ok = rep(TRUE, dim(Zc)[1L])
		for (t in 1:N) for (i in 1:n) { cols = if (bip) 1:n else setdiff(1:n, i); if (any(is.na(Y[i, cols, t]))) next; ok = ok & oracle_row(
			Zc[, i, cols, t],
			Y[i, cols, t]
		) }
		ok
	}
	acc = list(); tot = 0L
	while (tot < target) {
		m = 1e6L; w = sg^2 / (1 - rho^2)
		U = array(rnorm(m * n * N), c(m, n, N)); V = array(rnorm(m * n * N), c(m, n, N))
		U[, , 1] = sqrt(w) * U[, , 1]; V[, , 1] = sqrt(w) * V[, , 1]
		for (t in 2:N) { U[, , t] = rho * U[, , t - 1] + sg * U[, , t]; V[, , t] = rho * V[, , t - 1] + sg * V[, , t] }
		Z = array(0, c(m, n, n, N)); for (t in 1:N) for (i in 1:n) for (j in 1:n) Z[, i, j, t] = U[, i, t] * (if (bip) Gm[1, 1] else 1) * V[, j, t] + rnorm(m)
		k = which(feas(Z))
		acc[[length(acc) + 1L]] = list(
			U = U[k, , , drop = FALSE],
			V = V[k, , , drop = FALSE],
			Z = Z[k, , , , drop = FALSE]
		)
		tot = tot + length(k)
	}
	bind = function(nm) do.call(abind::abind, c(lapply(acc, `[[`, nm), list(along = 1)))
	U = bind('U'); V = bind('V'); Z = bind('Z'); m = dim(U)[1L]
	U1 = U; V1 = V; Z1 = Z
	for (r in seq_len(m)) {
		Uc = array(0, c(n, 1, N)); Vc = array(0, c(n, 1, N)); for (t in 1:N) { Uc[, 1, t] = U[r, , t]; Vc[, 1, t] = V[r, , t] }
		z3 = array(0, c(n, n, N)); for (t in 1:N) z3[, , t] = Z[r, , , t]
		o = lame:::.lame_npm_visit(z3, Y, fam, YL, cap, Uc, Vc, if (bip) Gm else NULL, rho, sg)
		for (t in 1:N) { U1[r, , t] = o$U[, 1, t]; V1[r, , t] = o$V[, 1, t]; Z1[r, , , t] = o$Z[, , t] }
	}
	st = function(U, V, Z) cbind(
		U[, 1, 1],
		U[, 2, 2],
		V[, 1, 1],
		V[, 3, 2],
		U[, 1, 1]^2,
		V[, 2, 2]^2,
		U[, 2, 1] * V[, 3, 1],
		U[, 3, 2] * V[, 1, 2],
		Z[, 1, 2, 1],
		Z[, 2, 3, 2],
		Z[, 3, 1, 1],
		Z[, 1, 1, 1],
		Z[, 2, 2, 2],
		Z[, 1, 2, 1]^2,
		U[, 1, 1] * U[, 1, 2],
		V[, 1, 1] * V[, 1, 2]
	)
	d = st(U1, V1, Z1) - st(U, V, Z)
	list(maxz = max(abs(colMeans(d) / (apply(d, 2L, sd) / sqrt(m)))), feasible = mean(feas(Z1)))
}

test_that('a visit of the paths with their latent rows and columns leaves the posterior of the augmented state unchanged (cbin bipartite, frn bipartite, cbin directed, frn directed)', {
	skip_on_cran()
	skip_if_not_installed('abind')
	for (cfg in list(c('cbin', 'bip'), c('frn', 'bip'), c('cbin', 'dir'), c('frn', 'dir'))) {
		o = npm_exact(cfg[1L], cfg[2L])
		expect_lt(o$maxz, 4.2)
		expect_equal(o$feasible, 1)
	}
})

test_that('a mutant that drops the data condition is caught', {
	skip_on_cran()
	skip_if_not_installed('abind')
	testthat::local_mocked_bindings(.lame_spm_row_ok = function(z, cls, kind, room) TRUE, .package = 'lame')
	expect_gt(npm_exact('cbin', 'bip', target = 1200L)$maxz, 6)
})

test_that('lame() runs the path visit for bipartite and directed cbin / frn fits with dynamic_uv, and the switch omits it', {
	skip_on_cran()
	set.seed(5)
	N = 3L
	rn = sprintf('r%02d', 1:8); cn = sprintf('c%02d', 1:7)
	bip = setNames(
		lapply(
			1:N,
			function(t) { Y = matrix(0, 8, 7); for (i in 1:8) { k = sample(0:2, 1L); if (k > 0L) Y[i, sample(7, k)] = seq_len(k) }; dimnames(Y) = list(rn, cn); Y }
		),
		paste0('p', 1:N)
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				bip,
				family = 'frn',
				mode = 'bipartite',
				R_row = 1L,
				R_col = 1L,
				odmax = 2L,
				dynamic_uv = TRUE,
				burn = 10,
				nscan = 40,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_equal(unname(f$sparse_position_move['visits']), 25)
	f0 = suppressWarnings(
		suppressMessages(
			lame(
				bip,
				family = 'frn',
				mode = 'bipartite',
				R_row = 1L,
				R_col = 1L,
				odmax = 2L,
				dynamic_uv = TRUE,
				burn = 10,
				nscan = 40,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(sparse_position_move = FALSE)
			)
		)
	)
	expect_null(f0$sparse_position_move)
	nm = sprintf('a%02d', 1:8)
	dir = setNames(
		lapply(
			1:N,
			function(t) { Y = matrix(0, 8, 8); for (i in 1:8) { k = sample(0:2, 1L); if (k > 0L) Y[i, sample(setdiff(1:8, i), k)] = 1 }; diag(Y) = NA; dimnames(Y) = list(nm, nm); Y }
		),
		paste0('p', 1:N)
	)
	fd = suppressWarnings(
		suppressMessages(
			lame(
				dir,
				family = 'cbin',
				R = 1L,
				odmax = 2L,
				dynamic_uv = TRUE,
				burn = 10,
				nscan = 40,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_equal(unname(fd$sparse_position_move['visits']), 25)
})
