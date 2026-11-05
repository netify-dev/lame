# receiver-level move of the cbin / frn samplers (R/nomination_actor_move.R): a receiver's level and latent column translate together
# checked: bounds of a translated column, invariance of residuals, labels and centring, the conditional law of the level, kernel exactness and the fits

# observed slice from its latent array: row i nominates its min(cap, number of positive cells)
# strongest cells (never an absent actor); cbin records 1, frn the rank (1 the weakest nominee)
nm_observe = function(Z, cap, family, bip, absent = integer(0)) {
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

# hide cells and drop actors, then number the ranks consecutively over what is left (as ame() and lame() do)
nm_mask = function(Y, family, hidden = integer(0), absent = integer(0)) {
	if (length(absent)) { Y[absent, ] = NA; Y[, absent] = NA }
	if (length(hidden)) Y[hidden] = NA
	if (family == 'frn') for (i in seq_len(nrow(Y))) {
		r = which(!is.na(Y[i, ]) & Y[i, ] > 0)
		if (length(r)) Y[i, r] = rank(Y[i, r])
	}
	Y
}

# the column holding each rank, row by row
nm_ranks = function(Y) {
	mx = max(c(0, Y), na.rm = TRUE)
	YL = matrix(NA_integer_, nrow(Y), max(mx, 1))
	for (i in seq_len(nrow(Y))) for (r in seq_len(mx)) {
		j = which(Y[i, ] == r)
		if (length(j)) YL[i, r] = j
	}
	YL
}

# is a latent row what the generative rule gives for its labelled cells? y1 the nominees, y0 the
# non-nominees, rank the observed ranks of the nominees (frn)
nm_row_ok = function(z, y1, y0, cap, rank = NULL) {
	cells = which(y1 | y0)
	pos = cells[z[cells] > 0]
	k = min(cap, length(pos))
	top = if (k > 0) pos[order(z[pos], decreasing = TRUE)[seq_len(k)]] else integer(0)
	ok = setequal(top, which(y1))
	if (ok && !is.null(rank)) {
		nom = which(y1 & !is.na(rank))
		if (length(nom) > 1) ok = identical(order(z[nom]), order(rank[nom]))
	}
	ok
}
nm_slice_ok = function(st, Z, Y, family, cap) {
	all(vapply(seq_len(nrow(Z)), function(i) {
		rk = NULL
		if (family == 'frn') { rk = Y[i, ]; rk[!st$isnom[i, ] | is.na(rk) | rk <= 0] = NA }
		nm_row_ok(Z[i, ], st$isnom[i, ], st$cls[i, ] == 0L, cap, rk)
	}, logical(1)))
}

test_that('the bounds of a translated column are exactly the shifts the generative rule allows', {
	skip_on_cran()
	set.seed(20261012)
	cap = 3; checked = 0
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) for (case in 1:4) {
		hidden = case %in% c(2, 4); absent = case >= 3 && !bip
		for (rep in 1:6) {
			nr = 8L; nc = if (bip) 6L else 8L
			a = rnorm(nr, 0, .6); b = rnorm(nc, 0, .6)
			Z = -0.8 + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
			gone = if (absent) sample(nr, 1) else integer(0)
			Y = nm_observe(Z, cap, family, bip, gone)
			cells = which(!is.na(Y))
			Y = nm_mask(
				Y,
				family,
				hidden = if (hidden) sample(
					setdiff(
						cells,
						c(
							outer(gone, seq_len(nr), function(i, j) i + (j - 1) * nr),
							outer(seq_len(nr), gone, function(i, j) i + (j - 1) * nr)
						)
					),
					5
				) else integer(0),
				absent = gone
			)
			YL = if (family == 'frn') nm_ranks(Y) else NULL
			st = lame:::.lame_nom_col_slice(Z, Y, family, YL, rep(cap, nr), bip)
			# the state is consistent with its labels
			expect_true(nm_slice_ok(st, Z, Y, family, cap))
			for (j in sample(nc, 3)) {
				lim = lame:::.lame_nom_col_limits(st, j)
				for (i in seq_len(nr)) {
					if (st$cls[i, j] == -1L) next
					lo = lim$lo[i]; hi = lim$hi[i]
					expect_true(Z[i, j] >= lo - 1e-12 && Z[i, j] <= hi + 1e-12)
					y1 = st$isnom[i, ]; y0 = st$cls[i, ] == 0L
					rk = NULL
					if (family == 'frn') { rk = Y[i, ]; rk[!y1 | is.na(rk) | rk <= 0] = NA }
					inside = if (is.finite(lo) && is.finite(hi)) runif(1, lo, hi) else if (is.finite(hi)) hi - rexp(1) else if (is.finite(lo)) lo + rexp(1) else rnorm(1)
					z2 = Z[i, ]; z2[j] = inside
					expect_true(nm_row_ok(z2, y1, y0, cap, rk))
					for (edge in c(lo, hi)) if (is.finite(edge)) {
						z3 = Z[i, ]; z3[j] = edge + if (edge == lo) -1e-6 else 1e-6
						expect_false(nm_row_ok(z3, y1, y0, cap, rk))
					}
					checked = checked + 1
				}
			}
		}
	}
	expect_gt(checked, 1000)
})

# one state of the move: a latent array consistent with its data, the levels and the variant settings
# (plain: free levels around a fixed intercept; center: levels sum to zero and the intercept moves with them; bip: independent sender and receiver levels)
nm_state = function(family, variant, hidden = FALSE, absent = FALSE, nt = 2L, nr = 8L, cap = 3) {
	bip = grepl('bip', variant); center = grepl('center', variant)
	nc = if (bip) 6L else nr
	a = rnorm(nr, 0, .5); if (center) a = a - mean(a)
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.2, 0.2, 0.5), 2)
	r_ba = if (bip) 0 else Sab[1, 2] / Sab[1, 1]
	s2 = if (bip) Sab[2, 2] else Sab[2, 2] - r_ba * Sab[1, 2]
	b = rnorm(nc, if (bip) 0 else a * r_ba, sqrt(s2)); if (center) b = b - mean(b)
	mu = if (center) 0.3 else -0.6
	Z = Y = array(NA_real_, c(nr, nc, nt))
	YL = vector('list', nt)
	for (t in seq_len(nt)) {
		Z[, , t] = mu + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
		gone = if (absent && t == 2L && !bip) c(2L, 5L) else integer(0)
		y = nm_observe(Z[, , t], cap, family, bip, gone)
		cells = setdiff(
			which(!is.na(y)),
			c(
				outer(gone, seq_len(nr), function(i, j) i + (j - 1) * nr),
				outer(seq_len(nr), gone, function(i, j) i + (j - 1) * nr)
			)
		)
		y = nm_mask(y, family, hidden = if (hidden) sample(cells, 6) else integer(0), absent = gone)
		Y[, , t] = y
		if (family == 'frn') YL[[t]] = nm_ranks(y)
	}
	list(
		Z = Z,
		Y = Y,
		YL = if (family == 'frn') YL else NULL,
		a = a,
		b = b,
		beta = if (center) mu else NULL,
		Sab = Sab,
		iV0 = if (center) matrix(0.8) else NULL,
		cap = cap,
		bip = bip,
		center = center,
		mu = mu,
		nr = nr,
		nc = nc,
		nt = nt
	)
}

test_that('the receiver move keeps the residuals, the labels, the centring and the constraints', {
	skip_on_cran()
	set.seed(5501)
	for (family in c('cbin', 'frn')) for (variant in c('plain', 'center', 'bip', 'bip_center')) for (case in 1:3) {
		s = nm_state(family, variant, hidden = case >= 2, absent = case == 3)
		cap = rep(s$cap, s$nr)
		sh = lame:::.lame_nom_col_shift(
			s$Z,
			s$Y,
			family,
			s$YL,
			cap,
			s$a,
			s$b,
			s$Sab,
			s$bip,
			rep(TRUE, s$nc),
			center = s$center,
			beta = s$beta,
			iV0 = s$iV0
		)
		expect_gt(sh$calls, 0)
		expect_equal(sh$calls + sh$skips, s$nc)
		# every residual is unchanged: the predictor moves with the latent value
		ez_old = (if (s$center) s$beta[1] else s$mu) + outer(s$a, s$b, '+')
		ez_new = (if (s$center) sh$beta[1] else s$mu) + outer(s$a, sh$b, '+')
		for (t in seq_len(s$nt)) expect_equal(sh$Z[, , t] - ez_new, s$Z[, , t] - ez_old)
		if (s$center) expect_equal(sum(sh$b), sum(s$b))
		expect_true(any(sh$b != s$b))
		# the data are still what the latent array says, with the same labels
		for (t in seq_len(s$nt)) {
			st0 = lame:::.lame_nom_col_slice(s$Z[, , t], s$Y[, , t], family, s$YL[[t]], cap, s$bip)
			st1 = lame:::.lame_nom_col_slice(sh$Z[, , t], s$Y[, , t], family, s$YL[[t]], cap, s$bip)
			expect_true(nm_slice_ok(st1, sh$Z[, , t], s$Y[, , t], family, s$cap))
			expect_identical(st1$cls, st0$cls)
		}
		# the columns that may not move stay where they are
		keep = c(TRUE, rep(FALSE, s$nc - 1L))
		sh2 = lame:::.lame_nom_col_shift(s$Z, s$Y, family, s$YL, cap, s$a, s$b, s$Sab, s$bip, keep, center = FALSE)
		expect_identical(sh2$b[-1], s$b[-1])
		expect_identical(sh2$Z[, -1, ], s$Z[, -1, ])
	}
})

test_that('the compiled sweep reproduces the R sweep from the same random numbers', {
	skip_on_cran()
	set.seed(5503)
	for (family in c('cbin', 'frn')) for (variant in c('plain', 'center', 'bip', 'bip_center')) for (case in 1:3) {
		s = nm_state(family, variant, hidden = case >= 2, absent = case == 3, nt = 1L + case %% 3L)
		cap = rep(s$cap, s$nr)
		run = function(engine) {
			set.seed(77 + case)
			lame:::.lame_nom_col_shift(
				s$Z,
				s$Y,
				family,
				s$YL,
				cap,
				s$a,
				s$b,
				s$Sab,
				s$bip,
				rep(TRUE, s$nc),
				center = s$center,
				beta = s$beta,
				iV0 = s$iV0,
				engine = engine
			)
		}
		r = run('r'); cpp = run('cpp')
		expect_gt(r$calls, 0)
		expect_equal(cpp$Z, r$Z, tolerance = 1e-12)
		expect_equal(cpp$b, r$b, tolerance = 1e-12)
		expect_equal(cpp$beta, r$beta, tolerance = 1e-12)
		expect_identical(cpp$calls, r$calls)
		expect_identical(cpp$skips, r$skips)
	}
})

test_that('one column draws its shift from the conditional prior cut to what the column allows', {
	skip_on_cran()
	set.seed(5502)
	reps = 2500
	for (family in c('cbin', 'frn')) for (variant in c('plain', 'center', 'bip_center')) {
		s = nm_state(family, variant, hidden = TRUE, nt = 1L)
		cap = rep(s$cap, s$nr)
		st = lame:::.lame_nom_col_slice(s$Z[, , 1], s$Y[, , 1], family, s$YL[[1]], cap, s$bip)
		r_ba = if (s$bip) 0 else s$Sab[1, 2] / s$Sab[1, 1]
		sd_b = sqrt(if (s$bip) s$Sab[2, 2] else s$Sab[2, 2] - r_ba * s$Sab[1, 2])
		m = if (s$bip) rep(0, s$nc) else s$a * r_ba
		for (j in c(1L, s$nc)) {
			lim = lame:::.lame_nom_col_limits(st, j)
			z = s$Z[, j, 1]
			dlo = max(-Inf, (lim$lo - z)[is.finite(lim$lo)])
			dhi = min(Inf, (lim$hi - z)[is.finite(lim$hi)])
			# the translation has the normal law N(mean, sd) cut to (dlo, dhi): the prior of the shifted level
			# (plain) or the prior of the levels and of the intercept along the compensated direction (centred)
			if (s$center) {
				a2 = (s$nc - 1) / s$nc / sd_b^2 + s$iV0[1, 1] / s$nc^2
				a1 = ((s$b[j] - m[j]) - mean(s$b - m)) / sd_b^2 + s$iV0[1, 1] * s$beta[1] / s$nc
				mean_d = -a1 / a2; sd_d = 1 / sqrt(a2)
			} else { mean_d = m[j] - s$b[j]; sd_d = sd_b }
			lo = (dlo - mean_d) / sd_d; hi = (dhi - mean_d) / sd_d
			mass = pnorm(hi) - pnorm(lo)
			exact_mean = mean_d + sd_d * (dnorm(lo) - dnorm(hi)) / mass
			exact_var = sd_d^2 * (1 + (ifelse(
				is.finite(lo),
				lo * dnorm(lo),
				0
			) - ifelse(is.finite(hi), hi * dnorm(hi), 0)) / mass -
					((dnorm(lo) - dnorm(hi)) / mass)^2)
			ok = seq_len(s$nc) == j
			d = replicate(reps, {
				sh = lame:::.lame_nom_col_shift(
					s$Z,
					s$Y,
					family,
					s$YL,
					cap,
					s$a,
					s$b,
					s$Sab,
					s$bip,
					ok,
					center = s$center,
					beta = s$beta,
					iV0 = s$iV0
				)
				if (s$center) (sh$b[j] - s$b[j]) / (1 - 1 / s$nc) else sh$b[j] - s$b[j]
			})
			expect_true(all(d >= dlo - 1e-9 & d <= dhi + 1e-9))
			expect_lt(abs(mean(d) - exact_mean), 5 * sqrt(exact_var / reps))
			expect_lt(abs(var(d) / exact_var - 1), 0.2)
		}
	}
})

# one sweep from the joint law of (levels, latent array, data) must keep that law, so paired differences of any
# function of (levels, data) before and after have mean zero (the z values returned)
nm_joint_law = function(family, variant, reps = 2500, hidden = FALSE, seed = 1, nr = 5L, cap = 2, engine = 'cpp') {
	set.seed(seed)
	bip = grepl('bip', variant); center = grepl('center', variant)
	nc = if (bip) 4L else nr
	a = c(0.3, -0.2, 0.4, -0.5, 0.1)[seq_len(nr)]; if (center) a = a - mean(a)
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.2, 0.2, 0.5), 2)
	r_ba = if (bip) 0 else Sab[1, 2] / Sab[1, 1]
	s2 = if (bip) Sab[2, 2] else Sab[2, 2] - r_ba * Sab[1, 2]
	m = if (bip) rep(0, nc) else a * r_ba
	iV0 = matrix(0.8)
	H = qr.Q(qr(contr.helmert(nc)))
	funs = function(b, beta0, e) c(b, b^2, b * e, b * (e == 0), if (center) c(beta0, beta0^2, beta0 * b))
	d = NULL
	for (r in seq_len(reps)) {
		b = if (center) m + as.numeric(H %*% rnorm(nc - 1, sd = sqrt(s2))) else rnorm(nc, m, sqrt(s2))
		beta0 = if (center) rnorm(1, 0, 1 / sqrt(iV0[1, 1])) else -0.5
		Z = beta0 + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
		Y = nm_observe(Z, cap, family, bip)
		if (hidden) Y = nm_mask(Y, family, hidden = sample(which(!is.na(Y)), 3)) else Y = nm_mask(Y, family)
		YL = if (family == 'frn') list(nm_ranks(Y)) else NULL
		e = colSums(!is.na(Y) & Y > 0)
		sh = lame:::.lame_nom_col_shift(
			array(Z, c(nr, nc, 1L)),
			array(Y, c(nr, nc, 1L)),
			family,
			YL,
			rep(cap, nr),
			a,
			b,
			Sab,
			bip,
			rep(TRUE, nc),
			center = center,
			beta = if (center) beta0 else NULL,
			iV0 = if (center) iV0 else NULL,
			engine = engine
		)
		after = funs(sh$b, if (center) sh$beta[1] else beta0, e)
		d = rbind(d, after - funs(b, beta0, e))
	}
	colMeans(d) / (apply(d, 2, sd) / sqrt(nrow(d)))
}

test_that('a sweep of receiver moves leaves the joint law of levels, latent values and data unchanged', {
	skip_on_cran()
	cfg = expand.grid(
		family = c('cbin', 'frn'),
		variant = c('plain', 'center', 'bip', 'bip_center'),
		stringsAsFactors = FALSE
	)
	for (k in seq_len(nrow(cfg))) {
		z = nm_joint_law(cfg$family[k], cfg$variant[k], hidden = k %% 2 == 0, seed = 700 + k)
		expect_lt(max(abs(z), na.rm = TRUE), 4.5)
	}
	# the check has power: a move that ignores the room of a row (non-nominees allowed above zero) fails it
	limits = get('.lame_nom_col_limits', asNamespace('lame'))
	testthat::local_mocked_bindings(
		.package = 'lame',
		.lame_nom_col_limits = function(st, j) {
			lim = limits(st, j)
			lim$hi[st$cls[, j] == 0L] = Inf
			lim
		}
	)
	expect_gt(max(abs(nm_joint_law('cbin', 'plain', reps = 1500, seed = 71, engine = 'r')), na.rm = TRUE), 6)
})

# panels of a few periods for the fits: levels around an intercept, covariates, hidden cells and
# (unipartite) actors absent from a period
nm_net = function(family, nt, n, cap, hidden = 0.1, windows = FALSE, bip = FALSE, nb = 8L, seed = 1) {
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
		z = -0.7 + 0.5 * x + outer(a[keep], b[kc], '+') + matrix(
			rnorm(length(keep) * length(kc)),
			length(keep),
			length(kc)
		)
		y = nm_observe(z, cap, family, bip)
		cells = which(!is.na(y))
		y = nm_mask(y, family, hidden = if (hidden > 0) sample(cells, round(hidden * length(cells))) else integer(0))
		dimnames(y) = list(rn[keep], cn[kc])
		Y[[t]] = y
		X[[t]] = array(x, c(dim(x), 1L), dimnames = list(rn[keep], cn[kc], 'x1'))
	}
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	list(Y = Y, X = X, cap = cap)
}

test_that('ame() runs the receiver move for cbin and frn and records it', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) {
		d = nm_net(family, 1L, 12L, 3, hidden = 0.1, seed = 31)
		run = function(prior) suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = family,
			odmax = d$cap,
			R = 0,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = prior,
			posterior_opts = list(save_ab = TRUE)
		)))
		on = run(list(actor_move_every = 2))
		expect_identical(
			names(on$actor_move),
			c('visits', 'col_calls', 'col_skips', 'level_calls', 'level_rows', 'rho_attempts', 'rho_accepts')
		)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['col_calls'] + on$actor_move['col_skips']), 40 * 12)
		expect_gt(unname(on$actor_move['col_calls']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		expect_null(run(list(actor_move = FALSE))$actor_move)
		# the default visits every sweep (20 burn-in and 60 production sweeps)
		expect_equal(unname(run(list())$actor_move['visits']), 80)
		# with the row effects switched off the receivers still move
		nr_fit = suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = family,
			odmax = d$cap,
			R = 0,
			rvar = FALSE,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = list(actor_move_every = 2)
		)))
		expect_equal(unname(nr_fit$actor_move['visits']), 40)
		expect_gt(unname(nr_fit$actor_move['col_calls']), 0)
		# a schedule that never reaches a visit consumes no random numbers: the draws are those of the switch-off
		never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		expect_identical(never$a_samples, off$a_samples); expect_identical(never$b_samples, off$b_samples)
	}
})

test_that('lame() runs the receiver move on changing panels with hidden cells, with and without an intercept', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) {
		d = nm_net(family, 3L, 12L, 3, hidden = 0.1, windows = TRUE, seed = 41)
		for (intercept in c(TRUE, FALSE)) {
			run = function(prior) suppressMessages(suppressWarnings(lame(
				d$Y,
				Xdyad = d$X,
				family = family,
				odmax = d$cap,
				R = 0,
				intercept = intercept,
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
			expect_identical(
				names(on$actor_move),
				c(
					'visits',
					'col_calls',
					'col_skips',
					if (intercept) c('level_calls', 'level_rows'),
					'rho_attempts',
					'rho_accepts'
				)
			)
			expect_equal(unname(on$actor_move['visits']), 40)
			expect_equal(unname(on$actor_move['col_calls'] + on$actor_move['col_skips']), 40 * 12)
			expect_gt(unname(on$actor_move['col_calls']), 0)
			expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
			if (intercept) expect_equal(unname(colSums(on$b_samples)), rep(0, ncol(on$b_samples)), tolerance = 1e-8)
			expect_null(run(list(actor_move = FALSE))$actor_move)
			expect_equal(unname(run(list())$actor_move['visits']), 80)
			never = run(list(actor_move_every = 10^6)); off = run(list(actor_move = FALSE))
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
			expect_identical(never$a_samples, off$a_samples); expect_identical(never$b_samples, off$b_samples)
		}
	}
})

test_that('lame() moves the receivers when the row effects are switched off', {
	skip_on_cran()
	d = nm_net('cbin', 3L, 12L, 3, hidden = 0.1, seed = 43)
	fit = suppressMessages(suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = 'cbin',
		odmax = d$cap,
		R = 0,
		rvar = FALSE,
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = list(actor_move_every = 2)
	)))
	expect_equal(unname(fit$actor_move['visits']), 40)
	expect_gt(unname(fit$actor_move['col_calls']), 0)
	expect_true(all(is.finite(fit$BETA)) && all(is.finite(fit$VC)))
})

test_that('lame() runs the receiver move on bipartite cbin and frn panels', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) {
		d = nm_net(family, 3L, 12L, 3, hidden = 0.1, bip = TRUE, nb = 8L, seed = 51)
		run = function(prior, intercept = TRUE) suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = family,
			odmax = d$cap,
			R = 0,
			mode = 'bipartite',
			intercept = intercept,
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
		for (intercept in c(TRUE, FALSE)) {
			on = run(list(actor_move_every = 2), intercept)
			expect_equal(unname(on$actor_move['visits']), 40)
			expect_equal(unname(on$actor_move['col_calls'] + on$actor_move['col_skips']), 40 * 8)
			expect_gt(unname(on$actor_move['col_calls']), 0)
			expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
			if (intercept) expect_equal(unname(colSums(on$b_samples)), rep(0, ncol(on$b_samples)), tolerance = 1e-8)
			never = run(list(actor_move_every = 10^6), intercept); off = run(list(actor_move = FALSE), intercept)
			expect_null(off$actor_move)
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
			expect_identical(never$b_samples, off$b_samples)
		}
	}
})

# an intercept that is the path of a dynamic_beta fit, or a static coefficient outside the design prior: the sweep reads
# its prior as the pair (q, l) of .lame_actor_move_icpt, so shifting the path by c changes the log prior by -(q c^2 + 2 l c) / 2

db_logpath = function(path, rho, sd2, m0, v0) {
	k = length(path)
	lp = -.5 * (path[1] - m0)^2 / v0
	if (k > 1) lp = lp - .5 * sum((path[-1] - rho * path[-k])^2) / sd2
	lp
}
# an intercept path from its prior
db_ic = function(kind, nt = 3L, rho = .6, sd2 = .15, m0 = -.5, v0 = .5, prec = .8) {
	if (kind == 'static') return(list(kind = 'static', value = rnorm(1, 0, 1 / sqrt(prec)), prec = prec))
	p = numeric(nt); p[1] = rnorm(1, m0, sqrt(v0))
	for (t in seq_len(nt)[-1]) p[t] = (if (kind == 'rw1') 1 else rho) * p[t - 1] + rnorm(1, 0, sqrt(sd2))
	list(kind = 'path', path = p, rho = if (kind == 'rw1') 1 else rho, sd2 = sd2, m0 = m0, v0 = v0)
}
db_icpt = function(ic, t) if (ic$kind == 'path') ic$path[t] else ic$value
# a centred state whose intercept is `ic`: levels summing to zero, latent array and data of nt periods
db_state = function(family, bip, ic, hidden = FALSE, absent = FALSE, nr = 8L, cap = 3, nt = 3L) {
	nc = if (bip) 6L else nr
	a = rnorm(nr, 0, .5); a = a - mean(a)
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.2, 0.2, 0.5), 2)
	r_ba = if (bip) 0 else Sab[1, 2] / Sab[1, 1]
	s2 = if (bip) Sab[2, 2] else Sab[2, 2] - r_ba * Sab[1, 2]
	b = rnorm(nc, if (bip) 0 else a * r_ba, sqrt(s2)); b = b - mean(b)
	Z = Y = array(NA_real_, c(nr, nc, nt))
	YL = vector('list', nt)
	for (t in seq_len(nt)) {
		Z[, , t] = db_icpt(ic, t) + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
		gone = if (absent && t == 2L && !bip) c(2L, 5L) else integer(0)
		y = nm_observe(Z[, , t], cap, family, bip, gone)
		cells = setdiff(
			which(!is.na(y)),
			c(
				outer(gone, seq_len(nr), function(i, j) i + (j - 1) * nr),
				outer(seq_len(nr), gone, function(i, j) i + (j - 1) * nr)
			)
		)
		y = nm_mask(y, family, hidden = if (hidden) sample(cells, 6) else integer(0), absent = gone)
		Y[, , t] = y
		if (family == 'frn') YL[[t]] = nm_ranks(y)
	}
	list(
		Z = Z,
		Y = Y,
		YL = if (family == 'frn') YL else NULL,
		a = a,
		b = b,
		Sab = Sab,
		cap = cap,
		bip = bip,
		nr = nr,
		nc = nc,
		nt = nt,
		ic = ic
	)
}
db_sweep = function(s, ic, ok = rep(TRUE, s$nc), engine = 'cpp', family = if (is.null(s$YL)) 'cbin' else 'frn')
	lame:::.lame_nom_col_shift(
		s$Z,
		s$Y,
		family,
		s$YL,
		rep(s$cap, s$nr),
		s$a,
		s$b,
		s$Sab,
		s$bip,
		ok,
		center = TRUE,
		beta = c(0, 0.4),
		iV0 = diag(c(0.3, 0.9)),
		engine = engine,
		ic = ic
	)

test_that('the receiver sweep with the prior of a moving intercept keeps the residuals and the centring, in both engines', {
	skip_on_cran()
	set.seed(9101)
	for (family in c('cbin', 'frn')) for (kind in c('ar1', 'rw1', 'static')) for (bip in c(FALSE, TRUE)) for (case in 1:3) {
		ic = db_ic(kind)
		s = db_state(family, bip, ic, hidden = case >= 2, absent = case == 3)
		set.seed(500 + case); r = db_sweep(s, ic, engine = 'r', family = family)
		set.seed(500 + case); cpp = db_sweep(s, ic, engine = 'cpp', family = family)
		expect_gt(r$calls, 0)
		expect_equal(cpp$Z, r$Z, tolerance = 1e-12)
		expect_equal(cpp$b, r$b, tolerance = 1e-12)
		expect_identical(cpp$calls, r$calls)
		expect_identical(cpp$skips, r$skips)
		if (kind == 'static') {
			expect_equal(cpp$ic$value, r$ic$value, tolerance = 1e-12)
			shift = cpp$ic$value - ic$value
		} else {
			expect_equal(cpp$ic$path, r$ic$path, tolerance = 1e-12)
			# the whole path moves by one constant
			shift = cpp$ic$path[1] - ic$path[1]
			expect_equal(cpp$ic$path - ic$path, rep(shift, s$nt), tolerance = 1e-10)
		}
		expect_gt(abs(shift), 0)
		# the coefficient vector of the caller is returned as it came: the caller writes the intercept back from `ic`
		expect_identical(cpp$beta, c(0, 0.4))
		# every residual is unchanged and the receiver levels still sum to zero
		for (t in seq_len(s$nt)) expect_equal(
			cpp$Z[, , t] - (db_icpt(cpp$ic, t) + outer(s$a, cpp$b, '+')),
			s$Z[, , t] - (db_icpt(ic, t) + outer(s$a, s$b, '+')),
			tolerance = 1e-9
		)
		expect_equal(sum(cpp$b), sum(s$b), tolerance = 1e-10)
		# the data are still what the latent array says
		cap = rep(s$cap, s$nr)
		for (t in seq_len(s$nt)) {
			st1 = lame:::.lame_nom_col_slice(cpp$Z[, , t], s$Y[, , t], family, s$YL[[t]], cap, bip)
			expect_true(nm_slice_ok(st1, cpp$Z[, , t], s$Y[, , t], family, s$cap))
		}
	}
})

test_that('a static intercept prior object reproduces the sweep with the coefficient prior it stands for', {
	skip_on_cran()
	set.seed(9102)
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) for (engine in c('r', 'cpp')) {
		s = nm_state(family, if (bip) 'bip_center' else 'center', hidden = TRUE)
		cap = rep(s$cap, s$nr)
		set.seed(81)
		plain = lame:::.lame_nom_col_shift(
			s$Z,
			s$Y,
			family,
			s$YL,
			cap,
			s$a,
			s$b,
			s$Sab,
			s$bip,
			rep(TRUE, s$nc),
			center = TRUE,
			beta = s$beta,
			iV0 = s$iV0,
			engine = engine
		)
		set.seed(81)
		ic = list(kind = 'static', value = s$beta, prec = s$iV0[1, 1])
		with_ic = lame:::.lame_nom_col_shift(
			s$Z,
			s$Y,
			family,
			s$YL,
			cap,
			s$a,
			s$b,
			s$Sab,
			s$bip,
			rep(TRUE, s$nc),
			center = TRUE,
			beta = s$beta,
			iV0 = s$iV0,
			engine = engine,
			ic = ic
		)
		expect_equal(with_ic$Z, plain$Z, tolerance = 1e-10)
		expect_equal(with_ic$b, plain$b, tolerance = 1e-10)
		expect_equal(with_ic$ic$value, plain$beta[1], tolerance = 1e-10)
		expect_identical(with_ic$calls, plain$calls)
	}
	# a state outside the bounds of every column leaves the intercept as it came
	s = nm_state('cbin', 'center')
	ic = db_ic('ar1')
	none = lame:::.lame_nom_col_shift(
		s$Z,
		s$Y,
		'cbin',
		NULL,
		rep(s$cap, s$nr),
		s$a,
		s$b,
		s$Sab,
		FALSE,
		rep(FALSE, s$nc),
		center = TRUE,
		beta = 0,
		iV0 = matrix(1),
		ic = ic
	)
	expect_identical(none$ic, ic); expect_identical(none$calls, 0)
})

test_that('one column draws its shift from the conditional prior with the path of the intercept in it', {
	skip_on_cran()
	set.seed(9103)
	reps = 3000
	for (family in c('cbin', 'frn')) for (kind in c('ar1', 'rw1', 'static')) for (bip in c(FALSE, TRUE)) {
		ic = db_ic(kind)
		s = db_state(family, bip, ic, hidden = TRUE)
		cap = rep(s$cap, s$nr)
		nc = s$nc
		r_ba = if (bip) 0 else s$Sab[1, 2] / s$Sab[1, 1]
		sd_b = sqrt(if (bip) s$Sab[2, 2] else s$Sab[2, 2] - r_ba * s$Sab[1, 2])
		m = if (bip) rep(0, nc) else s$a * r_ba
		for (j in c(1L, nc)) {
			# the bounds of the column from the cell classes of every period
			dlo = -Inf; dhi = Inf
			for (t in seq_len(s$nt)) {
				st = lame:::.lame_nom_col_slice(s$Z[, , t], s$Y[, , t], family, s$YL[[t]], cap, bip)
				lim = lame:::.lame_nom_col_limits(st, j)
				z = s$Z[, j, t]
				dlo = max(dlo, (lim$lo - z)[is.finite(lim$lo)]); dhi = min(dhi, (lim$hi - z)[is.finite(lim$hi)])
			}
			w = -rep(1 / nc, nc); w[j] = 1 - 1 / nc
			# the log prior along the translation (b + d w, intercept + d / nc) from the definitions of the model
			logp = function(d) {
				icd = if (kind == 'static') ic$value + d / nc else ic$path + d / nc
				-.5 * sum((s$b + d * w - m)^2) / sd_b^2 + (if (kind == 'static') -.5 * ic$prec * icd^2 else db_logpath(
					icd,
					ic$rho,
					ic$sd2,
					ic$m0,
					ic$v0
				))
			}
			a2 = -(logp(1) - 2 * logp(0) + logp(-1)); a1 = -(logp(1) - logp(-1)) / 2
			mu_d = -a1 / a2; sd_d = 1 / sqrt(a2)
			grid = seq(max(dlo, mu_d - 9 * sd_d), min(dhi, mu_d + 9 * sd_d), length.out = 40001)
			lp = vapply(grid, logp, 0)
			cdf = cumsum(exp(lp - max(lp))); cdf = cdf / cdf[length(cdf)]
			ok = seq_len(nc) == j
			d = replicate(reps, {
				sh = db_sweep(s, ic, ok, family = family)
				(sh$b[j] - s$b[j]) / (1 - 1 / nc)
			})
			expect_true(all(d >= dlo - 1e-9 & d <= dhi + 1e-9))
			u = approx(grid, cdf, xout = d, rule = 2)$y
			expect_gt(suppressWarnings(ks.test(u, 'punif')$p.value), 1e-3)
		}
	}
})

# one sweep from the joint law of (intercept path, levels, latent array, data): paired differences have mean zero
# (the z values returned); `ic_used` is the prior given to the sweep, and a wrong one fails the check
db_joint_law = function(
	family,
	bip,
	kind,
	reps = 2500,
	hidden = FALSE,
	seed = 1,
	nr = 5L,
	cap = 2,
	engine = 'cpp',
	ic_used = NULL
) {
	set.seed(seed)
	nc = if (bip) 4L else nr
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.2, 0.2, 0.5), 2)
	r_ba = if (bip) 0 else Sab[1, 2] / Sab[1, 1]
	s2 = if (bip) Sab[2, 2] else Sab[2, 2] - r_ba * Sab[1, 2]
	a = c(0.3, -0.2, 0.4, -0.5, 0.1)[seq_len(nr)]; a = a - mean(a)
	m = if (bip) rep(0, nc) else a * r_ba
	H = qr.Q(qr(contr.helmert(nc)))
	nt = 3L
	funs = function(b, ic, e) {
		p = if (ic$kind == 'path') ic$path else rep(ic$value, 2)
		c(b, b^2, b * e, b * (e == 0), p, p^2, mean(p), mean(p) * b, p[1] * p[length(p)])
	}
	d = NULL
	for (r in seq_len(reps)) {
		b = m + as.numeric(H %*% rnorm(nc - 1, sd = sqrt(s2)))
		ic = db_ic(kind, nt)
		Z = array(NA_real_, c(nr, nc, nt)); Y = array(NA_real_, c(nr, nc, nt)); YL = vector('list', nt)
		for (t in seq_len(nt)) {
			Z[, , t] = db_icpt(ic, t) + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
			y = nm_observe(Z[, , t], cap, family, bip)
			Y[, , t] = if (hidden) nm_mask(y, family, hidden = sample(which(!is.na(y)), 3)) else nm_mask(y, family)
			if (family == 'frn') YL[[t]] = nm_ranks(Y[, , t])
		}
		e = rowSums(apply(Y, 3, function(y) colSums(!is.na(y) & y > 0)))
		sh = lame:::.lame_nom_col_shift(
			Z,
			Y,
			family,
			YL,
			rep(cap, nr),
			a,
			b,
			Sab,
			bip,
			rep(TRUE, nc),
			center = TRUE,
			beta = c(0, 0),
			iV0 = diag(2),
			engine = engine,
			ic = if (is.null(ic_used)) ic else ic_used(ic)
		)
		d = rbind(d, funs(sh$b, sh$ic, e) - funs(b, ic, e))
	}
	keep = apply(d, 2, sd) > 1e-12
	colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrow(d)))
}

test_that('a sweep of receiver moves leaves the joint law of the intercept path, levels, latent values and data unchanged', {
	skip_on_cran()
	cfg = expand.grid(
		family = c('cbin', 'frn'),
		bip = c(FALSE, TRUE),
		kind = c('ar1', 'rw1', 'static'),
		stringsAsFactors = FALSE
	)
	for (k in seq_len(nrow(cfg))) {
		z = db_joint_law(cfg$family[k], cfg$bip[k], cfg$kind[k], hidden = k %% 2 == 0, seed = 1700 + k)
		expect_lt(max(abs(z), na.rm = TRUE), 4.6)
	}
	# the check has power: a sweep given a prior ten times too precise shifts the intercept too little and fails it, for every kind
	tight = function(ic) { if (ic$kind == 'path') { ic$v0 = ic$v0 / 10; ic$sd2 = ic$sd2 / 10 } else ic$prec = ic$prec * 10; ic }
	for (kind in c('ar1', 'rw1', 'static'))
		expect_gt(
			max(abs(db_joint_law('cbin', FALSE, kind, reps = 2000, seed = 1790, ic_used = tight)), na.rm = TRUE),
			8
		)
})
