# the positions of directed cbin, frn and ordinal fits drawn with the latent rows and columns of the actors they bound
# (R/sparse_position_move.R)

# exact-posterior check on a 3-actor panel: prior states (positions from N(0, Suv), latent values around their product)
# are kept when the family's data rule, written from its definition, maps them to the fixed data, so the kept states
# are posterior draws; one visit leaves every mean unchanged and the latent rows feasible
# returns the largest z score of the paired differences of 23 statistics and the share of feasible states after the visit
spm_exact = function(kind, ndraw = 2000L, seed = 3L, Suv = matrix(c(1, .4, .4, 1.3), 2)) {
	set.seed(seed)
	n = 3L
	off = which(row(diag(n)) != col(diag(n)))
	Y = matrix(NA_real_, n, n)
	cap = 2L
	if (kind == 'cbin') { cap = 1L; Y[off] = c(1, 0, 0, 0, 1, 0)[c(1, 2, 3, 4, 5, 6)] }
	if (kind == 'frn') Y[off] = c(2, 0, 0, 0, 1, 1)
	if (kind == 'ordinal') Y[off] = c(1, 3, 2, 3, 3, 2)
	if (kind == 'alpha') Y[off] = c(1, 3, 2, 3, 3, 2)
	alpha = if (kind == 'alpha') c(0, .8) else NULL
	fam = if (kind %in% c('ordinal', 'alpha')) 'ordinal' else kind
	Y3 = array(Y, c(n, n, 1L))
	YL = if (fam == 'frn') list(t(sapply(seq_len(n), function(i) match(seq_len(cap), Y[i, ])))) else NULL
	ctx = lame:::.lame_spm_context(
		Y3,
		fam,
		YL,
		if (fam != 'ordinal') cap else NULL,
		Y_int = if (kind == 'alpha') lame:::.lame_ltm_codes(Y3) else NULL
	)
	cells = lapply(seq_len(n), function(i) which(row(diag(n)) == i & col(diag(n)) != i))
	# the data rule, per row of the latent matrix (rows of `z`: draws; columns: the two off-diagonal cells)
	rule_ok = function(Z) {
		ok = rep(TRUE, dim(Z)[1L])
		if (fam == 'ordinal') {
			for (a in off) for (b in off) {
				if (kind == 'ordinal' && Y[a] < Y[b]) ok = ok & (Z[, a] < Z[, b])
			}
			if (kind == 'alpha') {
				af = c(-Inf, alpha, Inf)
				for (a in off) ok = ok & Z[, a] > af[Y[a]] & Z[, a] < af[Y[a] + 1L]
			}
			return(ok)
		}
		for (i in seq_len(n)) {
			c1 = cells[[i]][1L]; c2 = cells[[i]][2L]
			za = Z[, c1]; zb = Z[, c2]; p = Y[c(c1, c2)]
			if (fam == 'cbin') {
				na = za > 0 & za > zb; nb = zb > 0 & zb > za
				ok = ok & (if (p[1L] == 1) na else !na) & (if (p[2L] == 1) nb else !nb)
			} else {
				pa = za > 0; pb = zb > 0
				ra = ifelse(pa, 1 + (pb & zb < za), 0); rb = ifelse(pb, 1 + (pa & za < zb), 0)
				ok = ok & ra == p[1L] & rb == p[2L]
			}
		}
		ok
	}
	kept = list()
	tot = 0L
	while (tot < ndraw) {
		m = 1e6L
		W = matrix(rnorm(m * n * 2), m * n, 2) %*% chol(Suv)
		U = matrix(W[, 1L], m, n); V = matrix(W[, 2L], m, n)
		Z = matrix(0, m, n * n)
		for (i in seq_len(n)) for (j in seq_len(n)) Z[, i + (j - 1L) * n] = U[, i] * V[, j] + rnorm(m)
		w = which(rule_ok(Z))
		kept[[length(kept) + 1L]] = list(U = U[w, , drop = FALSE], V = V[w, , drop = FALSE], Z = Z[w, , drop = FALSE])
		tot = tot + length(w)
	}
	U = do.call(
		rbind,
		lapply(kept, `[[`, 'U')
	); V = do.call(rbind, lapply(kept, `[[`, 'V')); Z = do.call(rbind, lapply(kept, `[[`, 'Z'))
	m = nrow(U)
	U1 = U; V1 = V; Z1 = Z
	for (r in seq_len(m)) {
		o = lame:::.lame_spm_visit(
			ctx,
			array(Z[r, ], c(n, n, 1L)),
			matrix(U[r, ], n, 1L),
			matrix(V[r, ], n, 1L),
			shrink = TRUE,
			Suv = Suv,
			alpha = alpha
		)
		U1[r, ] = o$U; V1[r, ] = o$V; Z1[r, ] = o$Z
	}
	st = function(U, V, Z) cbind(
		U,
		V,
		U^2,
		V^2,
		U[, 1L] * V[, 2L],
		U[, 2L] * V[, 3L],
		U[, 3L] * V[, 1L],
		Z[, 4L],
		Z[, 8L],
		Z[, 3L],
		Z[, 1L],
		Z[, 5L],
		Z[, 4L]^2,
		Z[, 9L]^2,
		U[, 1L] * Z[, 4L]
	)
	d = st(U1, V1, Z1) - st(U, V, Z)
	z = colMeans(d) / (apply(d, 2L, sd) / sqrt(m))
	list(maxz = max(abs(z)), feasible = mean(rule_ok(Z1)), moved = mean(rowSums(abs(U1 - U)) > 0))
}

test_that('a visit of the positions with their latent rows and columns leaves the posterior of the augmented state unchanged (cbin, frn, ordinal ranks, ordinal cutpoints)', {
	skip_on_cran()
	for (k in c('cbin', 'frn', 'ordinal', 'alpha')) {
		o = spm_exact(k)
		expect_lt(o$maxz, 4.2)
		expect_equal(o$feasible, 1)
		expect_gt(o$moved, .9)
	}
})

test_that('a mutant that drops the data condition, and one that leaves the row aggregates stale, are caught', {
	skip_on_cran()
	# the mutants break the R visit: the compiled visit has its own copy of the conditions (test_compiled_position_moves.R)
	withr::local_options(lame.spm_engine = 'r')
	testthat::local_mocked_bindings(.lame_spm_row_ok = function(z, cls, kind, room) TRUE, .package = 'lame')
	expect_gt(spm_exact('cbin', 800L)$maxz, 6)
	testthat::local_mocked_bindings(.lame_spm_ord_ok = function(...) TRUE, .package = 'lame')
	expect_gt(spm_exact('ordinal', 800L)$maxz, 6)
})

test_that('stale aggregates after a sender move would leave the latent array infeasible', {
	skip_on_cran()
	withr::local_options(lame.spm_engine = 'r')
	testthat::local_mocked_bindings(
		.lame_spm_nom_set_row = function(st, i, znew) { st$Z[i, ] = znew; invisible(NULL) },
		.package = 'lame'
	)
	o = spm_exact('cbin', 2000L)
	expect_lt(o$feasible, 1)
})

test_that('the row condition of a cbin / frn slice is the condition the sweep sets cell by cell', {
	skip_on_cran()
	set.seed(5)
	n = 8L
	for (fam in c('cbin', 'frn')) {
		Y = matrix(0, n, n)
		for (i in seq_len(n)) { k = sample(0:3, 1L); if (k > 0L) Y[i, sample(setdiff(seq_len(n), i), k)] = if (fam == 'frn') seq_len(k) else 1 }
		diag(Y) = NA
		cap = 3L
		YL = if (fam == 'frn') t(sapply(seq_len(n), function(i) match(1:3, Y[i, ]))) else NULL
		Z = matrix(rnorm(n * n), n)
		Z = if (fam == 'cbin') rZ_cbin_fc(Z, matrix(0, n, n), 0, Y, cap, rowSums(Y > 0, na.rm = TRUE)) else rZ_frn_fc(
			Z,
			matrix(0, n, n),
			0,
			Y,
			YL,
			cap,
			rowSums(Y > 0, na.rm = TRUE)
		)
		st = lame:::.lame_spm_nom_slices(
			array(Z, c(n, n, 1L)),
			array(Y, c(n, n, 1L)),
			fam,
			if (fam == 'frn') list(YL) else NULL,
			cap
		)[[1L]]
		for (i in seq_len(n)) expect_true(lame:::.lame_spm_row_ok(Z[i, ], st$cls[i, ], st$kind, st$room_row[i]))
		# one cell of a row moved: inside its limits iff the row is feasible
		for (rep in 1:60) {
			j = sample(n, 1L)
			lim = lame:::.lame_nom_col_limits(st, j)
			i = sample(setdiff(seq_len(n), j), 1L)
			z = Z[i, ] + 0; z[j] = z[j] + rnorm(1L, 0, 2)
			inside = z[j] > lim$lo[i] && z[j] < lim$hi[i]
			expect_identical(lame:::.lame_spm_row_ok(z, st$cls[i, ], st$kind, st$room_row[i]), inside)
		}
	}
})

spm_panel = function(fam, n = 9L, N = 1L, seed = 2L) {
	set.seed(seed)
	one = function() {
		Y = matrix(0, n, n)
		if (fam == 'ordinal') Y = matrix(sample(1:3, n * n, TRUE, c(.6, .3, .1)), n, n)
		else for (i in seq_len(n)) { k = sample(0:2, 1L); if (k > 0L) Y[i, sample(setdiff(seq_len(n), i), k)] = if (fam == 'frn') seq_len(k) else 1 }
		diag(Y) = NA
		dimnames(Y) = list(letters[seq_len(n)], letters[seq_len(n)])
		Y
	}
	if (N == 1L) one() else setNames(lapply(seq_len(N), function(t) one()), paste0('p', seq_len(N)))
}

test_that('ame() and lame() run the visit for directed cbin, frn and ordinal fits and the switch omits it', {
	skip_on_cran()
	for (fam in c('cbin', 'frn', 'ordinal')) {
		od = if (fam == 'ordinal') NULL else 2L
		f = suppressWarnings(
			suppressMessages(
				ame(
					spm_panel(fam),
					family = fam,
					R = 1L,
					odmax = od,
					burn = 10,
					nscan = 30,
					odens = 1,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					seed = 1
				)
			)
		)
		expect_equal(unname(f$sparse_position_move['visits']), 40)
		expect_gt(f$sparse_position_move[['spm_calls']], 0)
		f0 = suppressWarnings(
			suppressMessages(
				ame(
					spm_panel(fam),
					family = fam,
					R = 1L,
					odmax = od,
					burn = 10,
					nscan = 30,
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
		fl = suppressWarnings(
			suppressMessages(
				lame(
					spm_panel(fam, N = 2L),
					family = fam,
					R = 1L,
					odmax = od,
					burn = 10,
					nscan = 30,
					odens = 1,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					seed = 1
				)
			)
		)
		expect_equal(unname(fl$sparse_position_move['visits']), 40)
	}
	expect_error(
		ame(
			spm_panel('cbin'),
			family = 'cbin',
			R = 1L,
			odmax = 2L,
			burn = 5,
			nscan = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(sparse_position_move = 'x')
		),
		'TRUE or FALSE'
	)
})

test_that('the visit is not run for panels it does not serve: hidden cbin cells, symmetric and positionless fits', {
	skip_on_cran()
	Y = spm_panel('cbin')
	Y[2, 3] = NA
	f = suppressWarnings(
		suppressMessages(
			ame(
				Y,
				family = 'cbin',
				R = 1L,
				odmax = 2L,
				burn = 20,
				nscan = 20,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_null(f$sparse_position_move)
	f = suppressWarnings(
		suppressMessages(
			ame(
				spm_panel('cbin'),
				family = 'cbin',
				R = 0L,
				odmax = 2L,
				burn = 20,
				nscan = 20,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_null(f$sparse_position_move)
})

test_that('sparse_position_every sets the sweeps between visits and is validated', {
	skip_on_cran()
	for (fun in c('ame', 'lame')) {
		Y = spm_panel('cbin', N = if (fun == 'lame') 2L else 1L)
		f = suppressWarnings(suppressMessages(do.call(
			fun,
			list(
				Y,
				family = 'cbin',
				R = 1L,
				odmax = 2L,
				burn = 10,
				nscan = 30,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(sparse_position_every = 4)
			)
		)))
		expect_equal(unname(f$sparse_position_move['visits']), 10)
		expect_error(
			do.call(
				fun,
				list(
					Y,
					family = 'cbin',
					R = 1L,
					odmax = 2L,
					burn = 5,
					nscan = 5,
					odens = 1,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = list(sparse_position_every = 0)
				)
			),
			'positive integer'
		)
	}
})
