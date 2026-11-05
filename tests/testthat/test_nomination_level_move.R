# level move of the cbin / frn samplers (R/nomination_actor_move.R, .lame_nom_level_shift): rows at their cap in every
# period translate with the intercept against rows with room, which pin the latent array from both sides. checks what
# the move changes and keeps, the law of the translation, exactness of the kernel under the joint law, the
# intercept-means ridge of bipartite ame(), and the fits that use them

# observed slice from its latent array: row i nominates its min(cap, positive cells) strongest cells, never an absent
# actor; cbin records 1, frn the rank (1 for the weakest nominee)
lv_observe = function(Z, cap, family, bip, absent = integer(0)) {
	Y = matrix(0, nrow(Z), ncol(Z))
	for (i in seq_len(nrow(Z))) {
		z = Z[i, ]
		if (!bip) z[i] = -Inf
		z[absent] = -Inf
		pos = which(z > 0)
		k = min(cap, length(pos))
		if (k > 0) {
			top = pos[order(z[pos], decreasing = TRUE)[seq_len(k)]]
			Y[i, top] = if (family == 'cbin') 1 else rank(z[top])
		}
	}
	if (!bip) diag(Y) = NA
	Y
}

# hide cells and drop actors, then renumber ranks over what is left, as ame() and lame() do
lv_mask = function(Y, family, hidden = integer(0), absent = integer(0)) {
	if (length(absent)) { Y[absent, ] = NA; Y[, absent] = NA }
	if (length(hidden)) Y[hidden] = NA
	if (family == 'frn') for (i in seq_len(nrow(Y))) {
		r = which(!is.na(Y[i, ]) & Y[i, ] > 0)
		if (length(r)) Y[i, r] = rank(Y[i, r])
	}
	Y
}

# whether a latent slice reproduces its data: every observed cell must match the rule when hidden cells are free
# to hold nominations (cbin: the set of nominees; frn: their set and order)
lv_consistent = function(Z, Y, cap, family, bip) {
	hid = lame:::.lame_nom_hidden_cells(Y, bip)
	ok = TRUE
	for (i in seq_len(nrow(Z))) {
		obs = !is.na(Y[i, ])
		if (!bip) obs[i] = FALSE
		el = obs | hid[i, ]
		z = Z[i, ]; z[!el] = -Inf
		pos = which(z > 0)
		k = min(cap, length(pos))
		top = if (k > 0) pos[order(z[pos], decreasing = TRUE)[seq_len(k)]] else integer(0)
		nominated = seq_along(z) %in% top
		ok = ok && all(nominated[obs] == (Y[i, obs] > 0))
		if (family == 'frn') {
			nom = which(obs & Y[i, ] > 0)
			if (length(nom) > 1) ok = ok && identical(order(z[nom]), order(Y[i, nom]))
		}
	}
	ok
}

# one state: intercept and a covariate with a coefficient prior, centred or free levels, latent array and data over
# nt periods, optional hidden cells and absent actors; bip: rectangular, no diagonal, no level covariance
lv_state = function(family, center, bip, hidden = FALSE, absent = FALSE, nt = 2L, nr = 8L, cap = 2, mu = -1.1) {
	nc = if (bip) nr - 2L else nr
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.2, 0.2, 0.5), 2)
	iV0 = matrix(c(0.8, 0.2, 0.2, 0.9), 2)
	beta = c(mu, 0.4)
	ab = matrix(rnorm(2 * nr), nr, 2) %*% t(chol(Sab))
	a = ab[, 1]
	b = if (bip) rnorm(nc, 0, sqrt(Sab[2, 2])) else ab[, 2]
	if (center) { a = a - mean(a); b = b - mean(b) }
	X = array(rnorm(nr * nc * nt), c(nr, nc, nt))
	Z = Y = array(NA_real_, c(nr, nc, nt))
	YL = vector('list', nt)
	for (t in seq_len(nt)) {
		Z[, , t] = beta[1] + beta[2] * X[, , t] + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
		gone = if (absent && t == 2L && !bip) c(2L, 5L) else integer(0)
		y = lv_observe(Z[, , t], cap, family, bip, gone)
		cells = setdiff(
			which(!is.na(y)),
			c(
				outer(gone, seq_len(nr), function(i, j) i + (j - 1) * nr),
				outer(seq_len(nr), gone, function(i, j) i + (j - 1) * nr)
			)
		)
		Y[, , t] = lv_mask(y, family, hidden = if (hidden) sample(cells, 6) else integer(0), absent = gone)
	}
	list(
		Z = Z,
		Y = Y,
		X = X,
		a = a,
		b = b,
		beta = beta,
		Sab = Sab,
		iV0 = iV0,
		cap = cap,
		bip = bip,
		center = center,
		nr = nr,
		nc = nc,
		nt = nt
	)
}
lv_eta = function(s, beta, a, b, t) beta[1] + beta[2] * s$X[, , t] + outer(a, b, '+')

test_that('the level move changes the free rows and the intercept and nothing else, and keeps the data', {
	skip_on_cran()
	set.seed(6601)
	made = 0
	for (family in c('cbin', 'frn')) for (center in c(TRUE, FALSE)) for (bip in c(FALSE, TRUE)) for (case in 1:3) {
		s = lv_state(family, center, bip, hidden = case >= 2, absent = case == 3)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, bip)
		bd = lame:::.lame_nom_shift_bounds(s$Z, cap, masks)
		free = bd$lo < 0 & is.infinite(bd$hi)
		out = lame:::.lame_nom_level_shift(
			s$Z,
			s$a,
			s$b,
			s$Sab,
			cap,
			masks,
			bip,
			rep(TRUE, s$nr),
			s$beta,
			s$iV0,
			center = center
		)
		if (!any(free)) { expect_identical(out$calls, 0); next }
		made = made + 1
		expect_identical(out$calls, 1)
		expect_identical(out$free, sum(free))
		# the free rows move by one common translation and every other row stays put
		d = out$Z - s$Z
		delta = d[which(free)[1], 1, 1]
		expect_true(all(abs(d[free, , ] - delta) < 1e-12))
		expect_identical(out$Z[!free, , ], s$Z[!free, , ])
		expect_gt(abs(delta), 0)
		# the translation respects the lower bound of the free rows
		expect_gt(delta, max(bd$lo[free]) - 1e-12)
		# every residual is unchanged since the predictor moves with the latent value
		for (t in seq_len(s$nt)) expect_equal(
			out$Z[, , t] - lv_eta(s, out$beta, out$a, s$b, t),
			s$Z[, , t] - lv_eta(s, s$beta, s$a, s$b, t)
		)
		# centred levels keep their sum; receivers and other coefficients are untouched
		if (center) expect_equal(sum(out$a), sum(s$a))
		expect_identical(out$beta[-1], s$beta[-1])
		# the pinned rows pay for the translation; free rows rise by delta n_p / n in the centred case
		if (center) expect_equal(out$beta[1] - s$beta[1], delta * sum(free) / s$nr)
		else expect_equal(out$beta[1] - s$beta[1], delta)
		# the latent array reproduces the data, with every row constraint
		for (t in seq_len(s$nt)) expect_true(lv_consistent(out$Z[, , t], s$Y[, , t], s$cap, family, bip))
	}
	expect_gt(made, 20)
	# rows the sparse-actor shrink leaves alone neither move nor count as free
	s = lv_state('cbin', TRUE, FALSE)
	cap = rep(s$cap, s$nr)
	masks = lame:::.lame_nom_shift_masks(s$Y, FALSE)
	out = lame:::.lame_nom_level_shift(
		s$Z,
		s$a,
		s$b,
		s$Sab,
		cap,
		masks,
		FALSE,
		c(FALSE, rep(TRUE, s$nr - 1L)),
		s$beta,
		s$iV0
	)
	expect_identical(out$Z[1, , ], s$Z[1, , ])
	# a panel with no row at its cap has nothing to move
	s0 = lv_state('cbin', TRUE, FALSE, mu = -4)
	m0 = lame:::.lame_nom_shift_masks(s0$Y, FALSE)
	o0 = lame:::.lame_nom_level_shift(
		s0$Z,
		s0$a,
		s0$b,
		s0$Sab,
		rep(s0$cap, s0$nr),
		m0,
		FALSE,
		rep(TRUE, s0$nr),
		s0$beta,
		s0$iV0
	)
	expect_identical(o0$calls, 0)
	expect_identical(o0$Z, s0$Z)
})

test_that('the translation is drawn from the conditional prior cut to what the free rows allow', {
	skip_on_cran()
	set.seed(6602)
	reps = 3000
	for (family in c('cbin', 'frn')) for (center in c(TRUE, FALSE)) for (bip in c(FALSE, TRUE)) {
		s = lv_state(family, center, bip, hidden = TRUE, nt = 1L, mu = -0.6)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, bip)
		bd = lame:::.lame_nom_shift_bounds(s$Z, cap, masks)
		free = bd$lo < 0 & is.infinite(bd$hi)
		expect_true(any(free) && !all(free))
		lower = max(bd$lo[free])
		# law of the translation along its line from the model's definitions: the log prior of the state moved by t,
		# under the centred or free parametrisation, on a grid
		n = s$nr; nf = sum(free)
		w = if (center) ifelse(free, (n - nf) / n, -nf / n) else ifelse(free, 0, -1)
		wb = if (center) nf / n else 1
		r_ab = if (bip) 0 else s$Sab[1, 2] / s$Sab[2, 2]
		m = if (bip) rep(0, n) else s$b * r_ab
		s2 = if (bip) s$Sab[1, 1] else s$Sab[1, 1] - r_ab * s$Sab[1, 2]
		logprior = function(t) {
			bt = s$beta; bt[1] = bt[1] + wb * t
			-0.5 * drop(t(bt) %*% s$iV0 %*% bt) - 0.5 * sum((s$a + w * t - m)^2) / s2
		}
		grid = seq(lower, lower + 12, length.out = 40001)
		lp = vapply(grid, logprior, 0)
		cdf = cumsum(exp(lp - max(lp))); cdf = cdf / cdf[length(cdf)]
		d = replicate(reps, {
			o = lame:::.lame_nom_level_shift(
				s$Z,
				s$a,
				s$b,
				s$Sab,
				cap,
				masks,
				bip,
				rep(TRUE, n),
				s$beta,
				s$iV0,
				center = center
			)
			(o$Z - s$Z)[which(free)[1], 1, 1]
		})
		expect_true(all(d > lower))
		u = approx(grid, cdf, xout = d, rule = 2)$y
		expect_gt(suppressWarnings(ks.test(u, 'punif')$p.value), 1e-3)
	}
})

# one move from the joint law of (coefficients, levels, latent array, data): a kernel invariant for the posterior
# keeps this law, so paired before/after differences of any state function have mean zero (the z values returned)
lv_joint_law = function(
	family,
	center,
	bip,
	reps = 3000,
	hidden = FALSE,
	absent = FALSE,
	seed = 1,
	nr = 6L,
	cap = 2,
	partial = FALSE
) {
	set.seed(seed)
	nc = if (bip) nr - 1L else nr
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.25, 0.25, 0.5), 2)
	iV0 = matrix(c(2, 0.15, 0.15, 0.8), 2)
	Lb = t(chol(solve(iV0))); root = t(chol(Sab))
	nt = 2L
	nab = min(nr, nc)
	funs = function(Z, beta, a, b) c(
		beta,
		beta^2,
		beta[1] * beta[2],
		a,
		a^2,
		beta[1] * a,
		sum(a[seq_len(nab)] * b[seq_len(nab)]),
		as.vector(Z[, , 1]),
		as.vector(Z[, , 1])^2,
		as.numeric(as.vector(Z) > 0),
		rowMeans(Z[, , nt])
	)
	d = NULL
	for (r in seq_len(reps)) {
		beta = drop(Lb %*% rnorm(2))
		ab = matrix(rnorm(2 * nr), nr, 2) %*% t(root)
		a = ab[, 1]; b = if (bip) rnorm(nc, 0, sqrt(Sab[2, 2])) else ab[, 2]
		if (center) { a = a - mean(a); b = b - mean(b) }
		X = array(rnorm(nr * nc * nt), c(nr, nc, nt))
		Z = Y = array(NA_real_, c(nr, nc, nt))
		for (t in seq_len(nt)) {
			Z[, , t] = beta[1] + beta[2] * X[, , t] + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
			gone = if (absent && t == 1L && !bip) nr else integer(0)
			y = lv_observe(Z[, , t], cap, family, bip, gone)
			cells = setdiff(
				which(!is.na(y)),
				c(
					outer(gone, seq_len(nr), function(i, j) i + (j - 1) * nr),
					outer(seq_len(nr), gone, function(i, j) i + (j - 1) * nr)
				)
			)
			Y[, , t] = lv_mask(y, family, hidden = if (hidden) sample(cells, 3) else integer(0), absent = gone)
		}
		masks = lame:::.lame_nom_shift_masks(Y, bip)
		ok = if (partial) c(rep(TRUE, nr - 2L), FALSE, FALSE) else rep(TRUE, nr)
		out = lame:::.lame_nom_level_shift(Z, a, b, Sab, rep(cap, nr), masks, bip, ok, beta, iV0, center = center)
		d = rbind(d, funs(out$Z, out$beta, out$a, b) - funs(Z, beta, a, b))
	}
	keep = apply(d, 2, sd) > 1e-12
	colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrow(d)))
}

test_that('the level move leaves the joint law of coefficients, levels, latent values and data unchanged', {
	skip_on_cran()
	cfg = expand.grid(
		family = c('cbin', 'frn'),
		center = c(TRUE, FALSE),
		bip = c(FALSE, TRUE),
		stringsAsFactors = FALSE
	)
	for (k in seq_len(nrow(cfg))) {
		z = lv_joint_law(
			cfg$family[k],
			cfg$center[k],
			cfg$bip[k],
			hidden = k %% 2 == 0,
			absent = k %% 3 == 0,
			partial = k %% 3 == 1,
			seed = 800 + k
		)
		expect_lt(max(abs(z)), 4.6)
	}
	# power: a move that treats rows with room as free (no upper bound) fails
	bounds = get('.lame_nom_shift_bounds', asNamespace('lame'))
	testthat::local_mocked_bindings(
		.package = 'lame',
		.lame_nom_shift_bounds = function(Z, odmax, masks) {
			bd = bounds(Z, odmax, masks)
			bd$hi[] = Inf
			bd
		}
	)
	expect_gt(max(abs(lv_joint_law('cbin', TRUE, FALSE, reps = 1500, seed = 91))), 6)
})

test_that('the ridge move of bipartite ame() leaves every predictor and the prior of the coefficients and levels unchanged', {
	skip_on_cran()
	set.seed(6603)
	nA = 7; nB = 5; va = 0.7; vb = 0.4
	iV0 = matrix(c(0.6, 0.2, 0.2, 0.9), 2)
	Lb = t(chol(solve(iV0)))
	reps = 15000
	f = function(beta, a, b) c(
		beta,
		beta^2,
		beta[1] * beta[2],
		mean(a),
		mean(a)^2,
		mean(b),
		mean(b)^2,
		a,
		a^2,
		b[1:2],
		beta[1] * mean(a),
		beta[1] * mean(b),
		mean(a) * mean(b),
		sum(a^2),
		sum(b^2)
	)
	d = NULL; worst = 0
	for (r in seq_len(reps)) {
		beta = drop(Lb %*% rnorm(2)); a = rnorm(nA, 0, sqrt(va)); b = rnorm(nB, 0, sqrt(vb))
		o = lame:::.lame_nom_ridge_bip(beta, a, b, iV0, va, vb)
		worst = max(worst, max(abs(outer(o$a, o$b, '+') + o$beta[1] - (outer(a, b, '+') + beta[1]))))
		expect_identical(o$beta[2], beta[2])
		d = rbind(d, f(o$beta, o$a, o$b) - f(beta, a, b))
	}
	expect_lt(worst, 1e-12)
	sdv = apply(d, 2, sd)
	z = colMeans(d) / (sdv / sqrt(reps))
	expect_lt(max(abs(z[sdv > 1e-12])), 4.6)
	# a margin outside the model is left alone
	a0 = rnorm(nA); b0 = rnorm(nB)
	o = lame:::.lame_nom_ridge_bip(c(0.3, -0.2), a0, b0, iV0, va, vb, rvar = FALSE)
	expect_identical(o$a, a0); expect_false(identical(o$b, b0)); expect_identical(o$calls, 1)
	o0 = lame:::.lame_nom_ridge_bip(c(0.3, -0.2), a0, b0, iV0, va, vb, rvar = FALSE, cvar = FALSE)
	expect_identical(o0$calls, 0)
})

# short panels for the fits: levels around an intercept, covariates, hidden cells and, unipartite, absent actors
lv_net = function(family, nt, n, cap, hidden = 0.1, windows = FALSE, bip = FALSE, nb = 8L, seed = 1, shift = -0.7) {
	set.seed(seed)
	nr = n; nc = if (bip) nb else n
	rn = sprintf('a%02d', seq_len(nr)); cn = if (bip) sprintf('b%02d', seq_len(nc)) else rn
	a = rnorm(nr, 0, .6); b = rnorm(nc, 0, .6)
	Y = X = vector('list', nt)
	for (t in seq_len(nt)) {
		keep = seq_len(nr)
		if (windows && !bip && t == 2L) keep = setdiff(keep, c(2L, 5L))
		kc = if (bip) seq_len(nc) else keep
		x = matrix(rnorm(length(keep) * length(kc)), length(keep), length(kc))
		z = shift + 0.5 * x + outer(a[keep], b[kc], '+') + matrix(
			rnorm(length(keep) * length(kc)),
			length(keep),
			length(kc)
		)
		y = lv_observe(z, cap, family, bip)
		cells = which(!is.na(y))
		y = lv_mask(y, family, hidden = if (hidden > 0) sample(cells, round(hidden * length(cells))) else integer(0))
		dimnames(y) = list(rn[keep], cn[kc])
		Y[[t]] = y
		X[[t]] = array(x, c(dim(x), 1L), dimnames = list(rn[keep], cn[kc], 'x1'))
	}
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	list(Y = Y, X = X, cap = cap)
}

test_that('lame() runs the level move for cbin and frn, unipartite and bipartite, and records it', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		d = lv_net(family, 3L, 14L, 2, hidden = 0.1, windows = !bip, bip = bip, nb = 9L, seed = 61)
		run = function(prior, ...) suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = family,
			odmax = d$cap,
			R = 0,
			mode = if (bip) 'bipartite' else 'unipartite',
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		)))
		on = run(list(actor_move_every = 2))
		expect_identical(
			names(on$actor_move),
			c(
				'visits',
				'col_calls',
				'col_skips',
				'level_calls',
				'level_rows',
				if (!bip) c('rho_attempts', 'rho_accepts')
			)
		)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_gt(unname(on$actor_move['level_calls']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_equal(unname(colSums(on$a_samples)), rep(0, ncol(on$a_samples)), tolerance = 1e-8)
		expect_equal(unname(colSums(on$b_samples)), rep(0, ncol(on$b_samples)), tolerance = 1e-8)
		# the default visits every sweep; a schedule that never reaches a visit equals the switch off
		expect_equal(unname(run(list())$actor_move['visits']), 80)
		never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
		expect_null(off$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_identical(never$a_samples, off$a_samples); expect_identical(never$b_samples, off$b_samples)
		# without an intercept the levels are not centred and there is nothing to move against
		ni = run(list(actor_move_every = 2), intercept = FALSE)
		expect_identical(
			names(ni$actor_move),
			c('visits', 'col_calls', 'col_skips', if (!bip) c('rho_attempts', 'rho_accepts'))
		)
	}
})

test_that('lame() translates the levels of a dynamic_ab fit of cbin and frn with their paths', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		d = lv_net(family, 4L, 14L, 2, hidden = 0.1, windows = FALSE, bip = bip, nb = 9L, seed = 62)
		run = function(prior) suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = family,
			odmax = d$cap,
			R = 0,
			mode = if (bip) 'bipartite' else 'unipartite',
			dynamic_ab = TRUE,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE)
		)))
		on = run(list(actor_move_every = 2))
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['col_calls'] + on$actor_move['col_skips']), 40 * (if (bip) 9 else 14))
		expect_gt(unname(on$actor_move['col_calls']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_true(all(is.finite(on$a_dynamic)) && all(is.finite(on$b_dynamic)))
		never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
		expect_null(off$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_identical(never$a_dynamic, off$a_dynamic)
	}
})

test_that('ame() runs the level move and, bipartite, the ridge move for cbin and frn', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		d = lv_net(family, 1L, 14L, 2, hidden = 0.1, bip = bip, nb = 9L, seed = 63)
		run = function(prior, ...) suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = family,
			odmax = d$cap,
			R = 0,
			mode = if (bip) 'bipartite' else 'unipartite',
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		)))
		on = run(list(actor_move_every = 2))
		expect_identical(
			names(on$actor_move),
			c(
				'visits',
				'col_calls',
				'col_skips',
				'level_calls',
				'level_rows',
				if (!bip) c('rho_attempts', 'rho_accepts'),
				if (bip) 'ridge_calls'
			)
		)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_gt(unname(on$actor_move['level_calls']), 0)
		if (bip) expect_equal(unname(on$actor_move['ridge_calls']), 40)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
		expect_null(off$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_identical(never$a_samples, off$a_samples)
		# the ridge does not need the row effects; with neither margin there is nothing to move
		nr_fit = run(list(actor_move_every = 2), rvar = FALSE)
		expect_true(all(is.finite(nr_fit$BETA)))
	}
})

test_that('bipartite ame() draws along the intercept-means ridge for the other families and records it', {
	skip_on_cran()
	set.seed(7720)
	nA = 14; nB = 9
	nms = list(sprintf('r%02d', seq_len(nA)), sprintf('c%02d', seq_len(nB)))
	x = array(rnorm(nA * nB), c(nA, nB, 1), dimnames = c(nms, list('x1')))
	eta = 0.4 + 0.4 * x[, , 1] + outer(rnorm(nA, 0, .6), rnorm(nB, 0, .6), '+')
	for (family in c('normal', 'poisson', 'ordinal')) {
		Y = switch(
			family,
			normal = eta + matrix(rnorm(nA * nB), nA),
			poisson = matrix(rpois(nA * nB, exp(eta - 0.3)), nA),
			ordinal = { z = eta + matrix(rnorm(nA * nB), nA); 0 + (z > 0.4) + (z > 1.4) }
		)
		dimnames(Y) = nms
		run = function(prior, ...) suppressMessages(suppressWarnings(ame(
			Y,
			Xdyad = x,
			family = family,
			mode = 'bipartite',
			R = 0,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		)))
		# the ordinal family fits without an intercept unless asked
		args = if (family == 'ordinal') list(intercept = TRUE) else list()
		on = do.call(run, c(list(list(actor_move_every = 2)), args))
		# Poisson and ordinal fits also run the level translations with the latent array (R/latent_translation_move.R),
		# recorded first; the ridge draws are the last field
		if (family == 'normal') {
			expect_identical(names(on$actor_move), c('visits', 'ridge_calls'))
			expect_equal(unname(on$actor_move), c(40, 40))
		} else {
			expect_identical(names(on$actor_move)[c(1L, length(on$actor_move))], c('visits', 'ridge_calls'))
			expect_equal(unname(on$actor_move[c('visits', 'ridge_calls')]), c(40, 40))
		}
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_equal(unname(do.call(run, c(list(list()), args))$actor_move['ridge_calls']), 80)
		off = do.call(run, c(list(list(actor_move = FALSE)), args))
		expect_null(off$actor_move)
		never = do.call(run, c(list(list(actor_move_every = 10^6)), args))
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC); expect_identical(
			never$a_samples,
			off$a_samples
		)
	}
	# without an intercept there is no ridge; the ordinal fit runs the location draw of the translations
	expect_false('ridge_calls' %in% names(run(list(), intercept = FALSE)$actor_move))
})

# an intercept that is a dynamic_beta path or a static coefficient outside the design prior: the move reads its
# prior as the pair (q, l) of .lame_actor_move_icpt, so a shift by c changes the log prior by -(q c^2 + 2 l c) / 2

lvp_logpath = function(path, rho, sd2, m0, v0) {
	k = length(path)
	lp = -.5 * (path[1] - m0)^2 / v0
	if (k > 1) lp = lp - .5 * sum((path[-1] - rho * path[-k])^2) / sd2
	lp
}
lvp_ic = function(kind, nt = 3L, rho = .6, sd2 = .15, m0 = -.5, v0 = .5, prec = .8) {
	if (kind == 'static') return(list(kind = 'static', value = rnorm(1, 0, 1 / sqrt(prec)), prec = prec))
	p = numeric(nt); p[1] = rnorm(1, m0, sqrt(v0))
	for (t in seq_len(nt)[-1]) p[t] = (if (kind == 'rw1') 1 else rho) * p[t - 1] + rnorm(1, 0, sqrt(sd2))
	list(kind = 'path', path = p, rho = if (kind == 'rw1') 1 else rho, sd2 = sd2, m0 = m0, v0 = v0)
}
lvp_icpt = function(ic, t) if (ic$kind == 'path') ic$path[t] else ic$value

test_that('the level move with the prior of a moving intercept changes the free rows and the path and nothing else', {
	skip_on_cran()
	set.seed(6611)
	made = 0
	for (family in c('cbin', 'frn')) for (kind in c('ar1', 'rw1', 'static')) for (bip in c(FALSE, TRUE)) for (case in 1:3) {
		ic = lvp_ic(kind, 3L)
		s = lv_state(family, TRUE, bip, hidden = case >= 2, absent = case == 3, nt = 3L)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, bip)
		bd = lame:::.lame_nom_shift_bounds(s$Z, cap, masks)
		free = bd$lo < 0 & is.infinite(bd$hi)
		out = lame:::.lame_nom_level_shift(
			s$Z,
			s$a,
			s$b,
			s$Sab,
			cap,
			masks,
			bip,
			rep(TRUE, s$nr),
			s$beta,
			s$iV0,
			ic = ic
		)
		if (!any(free)) { expect_identical(out$calls, 0); expect_identical(out$ic, ic); next }
		made = made + 1
		expect_identical(out$calls, 1)
		d = out$Z - s$Z
		delta = d[which(free)[1], 1, 1]
		expect_true(all(abs(d[free, , ] - delta) < 1e-12))
		expect_identical(out$Z[!free, , ], s$Z[!free, , ])
		expect_gt(delta, max(bd$lo[free]) - 1e-12)
		# the caller's coefficient vector comes back unchanged; the intercept lives in `ic`
		expect_identical(out$beta, s$beta)
		shift = if (kind == 'static') out$ic$value - ic$value else out$ic$path - ic$path
		expect_equal(shift, rep(delta * sum(free) / s$nr, length(shift)), tolerance = 1e-10)
		# every residual is unchanged and the levels sum to zero
		for (t in seq_len(s$nt)) expect_equal(
			out$Z[, , t] - (lvp_icpt(out$ic, t) + outer(out$a, s$b, '+')),
			s$Z[, , t] - (lvp_icpt(ic, t) + outer(s$a, s$b, '+')),
			tolerance = 1e-9
		)
		expect_equal(sum(out$a), sum(s$a), tolerance = 1e-10)
		for (t in seq_len(s$nt)) expect_true(lv_consistent(out$Z[, , t], s$Y[, , t], s$cap, family, bip))
	}
	expect_gt(made, 20)
})

test_that('the translation is drawn from the conditional prior with the path of the intercept in it', {
	skip_on_cran()
	set.seed(6612)
	reps = 3000
	for (family in c('cbin', 'frn')) for (kind in c('ar1', 'rw1', 'static')) for (bip in c(FALSE, TRUE)) {
		ic = lvp_ic(kind, 3L)
		s = lv_state(family, TRUE, bip, hidden = TRUE, nt = 3L, mu = -0.6)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, bip)
		bd = lame:::.lame_nom_shift_bounds(s$Z, cap, masks)
		free = bd$lo < 0 & is.infinite(bd$hi)
		expect_true(any(free) && !all(free))
		lower = max(bd$lo[free])
		n = s$nr; nf = sum(free)
		w = ifelse(free, (n - nf) / n, -nf / n)
		wb = nf / n
		r_ab = if (bip) 0 else s$Sab[1, 2] / s$Sab[2, 2]
		m = if (bip) rep(0, n) else s$b * r_ab
		s2 = if (bip) s$Sab[1, 1] else s$Sab[1, 1] - r_ab * s$Sab[1, 2]
		# log prior along the translation (free rows' levels by w t, the intercept by wb t) from the model's definitions
		logp = function(t) {
			-.5 * sum((s$a + w * t - m)^2) / s2 +
				(if (kind == 'static') -.5 * ic$prec * (ic$value + wb * t)^2 else lvp_logpath(
					ic$path + wb * t,
					ic$rho,
					ic$sd2,
					ic$m0,
					ic$v0
				))
		}
		a2 = -(logp(1) - 2 * logp(0) + logp(-1)); a1 = -(logp(1) - logp(-1)) / 2
		sd_t = 1 / sqrt(a2)
		grid = seq(lower, max(lower, -a1 / a2) + 9 * sd_t, length.out = 40001)
		lp = vapply(grid, logp, 0)
		cdf = cumsum(exp(lp - max(lp))); cdf = cdf / cdf[length(cdf)]
		d = replicate(reps, {
			o = lame:::.lame_nom_level_shift(
				s$Z,
				s$a,
				s$b,
				s$Sab,
				cap,
				masks,
				bip,
				rep(TRUE, n),
				s$beta,
				s$iV0,
				ic = ic
			)
			(o$Z - s$Z)[which(free)[1], 1, 1]
		})
		expect_true(all(d > lower))
		u = approx(grid, cdf, xout = d, rule = 2)$y
		expect_gt(suppressWarnings(ks.test(u, 'punif')$p.value), 1e-3)
	}
})

# one move from the joint law of (intercept path, covariate coefficient, levels, latent array, data): paired
# before/after differences of any state function have mean zero (the z values returned)
lvp_joint_law = function(family, bip, kind, reps = 3000, hidden = FALSE, seed = 1, nr = 6L, cap = 2, ic_used = NULL) {
	set.seed(seed)
	nc = if (bip) nr - 1L else nr
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.25, 0.25, 0.5), 2)
	root = t(chol(Sab))
	nt = 3L; nab = min(nr, nc)
	funs = function(Z, ic, b2, a, b) {
		p = if (ic$kind == 'path') ic$path else rep(ic$value, 2)
		c(
			p,
			p^2,
			mean(p),
			p[1] * p[length(p)],
			b2,
			a,
			a^2,
			mean(p) * a,
			p[1] * a,
			sum(a[seq_len(nab)] * b[seq_len(nab)]),
			as.vector(Z[, , 1]),
			as.numeric(as.vector(Z) > 0),
			rowMeans(Z[, , nt])
		)
	}
	d = NULL
	for (r in seq_len(reps)) {
		ic = lvp_ic(kind, nt)
		b2 = rnorm(1, 0, 1)
		ab = matrix(rnorm(2 * nr), nr, 2) %*% t(root)
		a = ab[, 1]; b = if (bip) rnorm(nc, 0, sqrt(Sab[2, 2])) else ab[, 2]
		a = a - mean(a); b = b - mean(b)
		X = array(rnorm(nr * nc * nt), c(nr, nc, nt))
		Z = Y = array(NA_real_, c(nr, nc, nt))
		for (t in seq_len(nt)) {
			Z[, , t] = lvp_icpt(ic, t) + b2 * X[, , t] + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
			y = lv_observe(Z[, , t], cap, family, bip)
			Y[, , t] = lv_mask(y, family, hidden = if (hidden) sample(which(!is.na(y)), 3) else integer(0))
		}
		masks = lame:::.lame_nom_shift_masks(Y, bip)
		out = lame:::.lame_nom_level_shift(
			Z,
			a,
			b,
			Sab,
			rep(cap, nr),
			masks,
			bip,
			rep(TRUE, nr),
			c(0, b2),
			diag(2),
			ic = if (is.null(ic_used)) ic else ic_used(ic)
		)
		d = rbind(d, funs(out$Z, out$ic, b2, out$a, b) - funs(Z, ic, b2, a, b))
	}
	keep = apply(d, 2, sd) > 1e-12
	colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrow(d)))
}

test_that('the level move with a moving intercept leaves the joint law of path, levels, latent values and data unchanged', {
	skip_on_cran()
	cfg = expand.grid(
		family = c('cbin', 'frn'),
		bip = c(FALSE, TRUE),
		kind = c('ar1', 'rw1', 'static'),
		stringsAsFactors = FALSE
	)
	for (k in seq_len(nrow(cfg))) {
		z = lvp_joint_law(cfg$family[k], cfg$bip[k], cfg$kind[k], hidden = k %% 2 == 0, seed = 2800 + k)
		expect_lt(max(abs(z)), 4.6)
	}
	# power: a move given a ten times too precise prior shifts the intercept too little and fails, for every kind
	tight = function(ic) { if (ic$kind == 'path') { ic$v0 = ic$v0 / 10; ic$sd2 = ic$sd2 / 10 } else ic$prec = ic$prec * 10; ic }
	for (kind in c('ar1', 'rw1', 'static'))
		expect_gt(max(abs(lvp_joint_law('cbin', FALSE, kind, reps = 2000, seed = 2890, ic_used = tight))), 8)
})

test_that('lame() runs the level and receiver moves with dynamic_beta, with and without per-actor slopes, and records them', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		d = lv_net(family, 4L, 14L, 2, hidden = 0.1, windows = !bip, bip = bip, nb = 9L, seed = 64)
		run = function(prior, ...) suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = family,
			odmax = d$cap,
			R = 0,
			mode = if (bip) 'bipartite' else 'unipartite',
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		)))
		specs = list(
			list(dynamic_beta = 'intercept', dynamic_beta_kind = 'ar1'),
			list(dynamic_beta = 'intercept', dynamic_beta_kind = 'rw1'),
			list(dynamic_beta = 'dyad', dynamic_beta_kind = 'ar1'), # the intercept stays a static coefficient
			list(dynamic_beta = TRUE, dynamic_beta_kind = 'rw1'),
			list(dynamic_beta_per_actor = 'row'), # per-actor slopes, static coefficients
			list(dynamic_beta_per_actor = 'col', dynamic_beta = 'intercept', dynamic_beta_kind = 'ar1'),
			# the level and its path move with the intercept path
			list(dynamic_ab = TRUE, dynamic_beta = 'intercept', dynamic_beta_kind = 'rw1')
		)
		for (sp in specs) {
			on = do.call(run, c(list(list(actor_move_every = 2)), sp))
			expect_equal(unname(on$actor_move['visits']), 40)
			expect_gt(unname(on$actor_move['level_calls']), 0)
			expect_gt(unname(on$actor_move['col_calls']), 0)
			# the state handed to the moves stays inside the bounds of its labels, which the collapsed reciprocity update
			# must respect with actor slopes in the latent array
			if (!bip && !is.null(sp$dynamic_beta_per_actor)) expect_equal(unname(on$actor_move['col_skips']), 0)
			expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
			if (is.null(sp$dynamic_ab)) {
				expect_equal(unname(colSums(on$a_samples)), rep(0, ncol(on$a_samples)), tolerance = 1e-8)
				expect_equal(unname(colSums(on$b_samples)), rep(0, ncol(on$b_samples)), tolerance = 1e-8)
			} else expect_true(all(is.finite(on$a_dynamic)) && all(is.finite(on$b_dynamic)))
			never = do.call(
				run,
				c(list(list(actor_move_every = 10^6)), sp)
			); off = do.call(run, c(list(list(actor_move = FALSE)), sp))
			expect_null(off$actor_move)
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
			expect_identical(never$a_samples, off$a_samples); expect_identical(never$b_samples, off$b_samples)
		}
		# without an intercept the dynamic coefficients carry no level to move against, so only the receivers move
		ni = run(list(actor_move_every = 2), intercept = FALSE, dynamic_beta = 'dyad', dynamic_beta_kind = 'ar1')
		expect_identical(
			names(ni$actor_move),
			c('visits', 'col_calls', 'col_skips', if (!bip) c('rho_attempts', 'rho_accepts'))
		)
	}
	# kinds and grids the intercept prior does not cover keep the sampler's own updates
	d = lv_net('cbin', 5L, 12L, 2, hidden = 0.1, seed = 65)
	run2 = function(...) suppressMessages(suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = 'cbin',
		odmax = d$cap,
		R = 0,
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = list(actor_move_every = 2),
		...
	)))
	for (extra in list(
		list(dynamic_beta = 'intercept', dynamic_beta_kind = 'rw2'),
		list(dynamic_beta = 'intercept', time_index = c(1, 2, 4, 5, 6))
	)) {
		fit = do.call(run2, extra)
		expect_false(any(c('level_calls', 'col_calls') %in% names(fit$actor_move)))
		expect_true(all(is.finite(fit$BETA)))
	}
})

test_that('the sender translation of a dynamic_ab fit runs under prior$actor_move and on its schedule', {
	skip_on_cran()
	calls = 0L
	orig = get('.lame_nom_row_shift', asNamespace('lame'))
	testthat::local_mocked_bindings(
		.package = 'lame',
		.lame_nom_row_shift = function(...) { calls <<- calls + 1L; orig(...) }
	)
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		d = lv_net(family, 4L, 14L, 2, hidden = 0.1, windows = FALSE, bip = bip, nb = 9L, seed = 66)
		run = function(prior, ...) { calls <<- 0L; suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = family,
			odmax = d$cap,
			R = 0,
			mode = if (bip) 'bipartite' else 'unipartite',
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE),
			...
		))) }
		run(list(), dynamic_ab = TRUE); expect_equal(calls, 80L)
		run(list(actor_move_every = 4), dynamic_ab = TRUE); expect_equal(calls, 20L)
		run(list(actor_move = FALSE), dynamic_ab = TRUE); expect_equal(calls, 0L)
		run(list(actor_move_every = 10^6), dynamic_ab = TRUE); expect_equal(calls, 0L)
		# a static fit makes the translation every sweep, whatever the switch
		run(list(actor_move = FALSE)); expect_equal(calls, 80L)
	}
})
