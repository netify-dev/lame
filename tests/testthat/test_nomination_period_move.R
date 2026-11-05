# period-level translations of the intercept path of cbin and frn dynamic_beta fits
# (R/nomination_period_move.R, R/nomination_room_move.R): .lame_nom_period_shift is an exact Gibbs
# translation of the periods that rows with room do not tie together, and .lame_nom_room_period is a
# collapsed Metropolis translation of each such period's intercept with the non-nominations integrated out

pm_observe = function(Z, cap, family, bip) {
	Y = matrix(0, nrow(Z), ncol(Z))
	for (i in seq_len(nrow(Z))) {
		z = Z[i, ]
		if (!bip) z[i] = -Inf
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

# does a latent slice reproduce its data under the generative rule? (no hidden cells here)
pm_consistent = function(Z, Y, cap, family, bip) {
	ok = TRUE
	for (i in seq_len(nrow(Z))) {
		obs = !is.na(Y[i, ])
		z = Z[i, ]
		z[!obs] = -Inf
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

# one state: an AR(1) intercept path (first state N(m0, v0), innovations N(0, sd2)), centred levels, latent array and data of nt periods
pm_state = function(family, bip, nt = 4L, nr = 9L, cap = 2, rho = 0, m0 = -0.8, v0 = 0.5, sd2 = 0.16, rho_p = 0.8) {
	nc = if (bip) nr - 2L else nr
	Sab = if (bip) diag(c(0.6, 0.5)) else matrix(c(0.6, 0.25, 0.25, 0.5), 2)
	path = numeric(nt)
	path[1] = rnorm(1, m0, sqrt(v0))
	for (t in seq_len(nt)[-1]) path[t] = rho_p * path[t - 1] + rnorm(1, 0, sqrt(sd2))
	ab = matrix(rnorm(2 * nr), nr, 2) %*% chol(Sab) # rows with covariance Sab
	a = if (bip) rnorm(nr, 0, sqrt(Sab[1, 1])) else ab[, 1]
	b = if (bip) rnorm(nc, 0, sqrt(Sab[2, 2])) else ab[, 2]
	a = a - mean(a)
	b = b - mean(b)
	Z = EZ = Y = array(NA_real_, c(nr, nc, nt))
	for (t in seq_len(nt)) {
		eta = path[t] + outer(a, b, '+')
		e = matrix(rnorm(nr * nc), nr, nc)
		if (!bip && rho != 0) {
			e2 = matrix(rnorm(nr * nr), nr, nr)
			up = upper.tri(e)
			off = !up & row(e) != col(e)
			e[off] = (rho * t(e) + sqrt(1 - rho^2) * e2)[off]
			diag(e) = rnorm(nr, 0, sqrt(1 + rho))
		}
		Z[, , t] = eta + e
		EZ[, , t] = eta
		Y[, , t] = pm_observe(Z[, , t], cap, family, bip)
	}
	list(
		Z = Z,
		EZ = EZ,
		Y = Y,
		a = a,
		b = b,
		path = path,
		Sab = Sab,
		cap = cap,
		bip = bip,
		family = family,
		nr = nr,
		nc = nc,
		nt = nt,
		rho = rho,
		ic = list(kind = 'path', k = 1L, path = path, rho = rho_p, sd2 = sd2, m0 = m0, v0 = v0)
	)
}
pm_eta = function(s, path, a, t) path[t] + outer(a, s$b, '+')

# the components

test_that('the periods and the rows with room fall into components', {
	free = matrix(TRUE, 7, 5)
	free[1, c(1, 2)] = FALSE
	free[2, 2] = FALSE
	free[3, 4] = FALSE
	free[4, c(4, 5)] = FALSE
	cm = lame:::.lame_nom_period_components(free)
	periods = lapply(unname(cm), function(x) x$periods)
	expect_length(cm, 3)
	expect_true(any(vapply(periods, identical, TRUE, c(1L, 2L))))
	expect_true(any(vapply(periods, identical, TRUE, 3L)))
	expect_true(any(vapply(periods, identical, TRUE, c(4L, 5L))))
	for (x in cm) {
		if (identical(x$periods, c(1L, 2L))) expect_identical(x$rows, c(1L, 2L))
		if (identical(x$periods, 3L)) expect_identical(x$rows, integer(0))
		if (identical(x$periods, c(4L, 5L))) expect_identical(x$rows, c(3L, 4L))
	}
	# a row pinned in two periods ties them, a chain ties them all
	free2 = matrix(TRUE, 4, 4)
	free2[1, 1:2] = FALSE
	free2[2, 2:3] = FALSE
	free2[3, 3:4] = FALSE
	expect_length(lame:::.lame_nom_period_components(free2), 1L)
	expect_length(lame:::.lame_nom_period_components(matrix(TRUE, 4, 4)), 4L)
})

# the exact translation

test_that('the period translation keeps the data, every residual and the centring', {
	skip_on_cran()
	set.seed(7101)
	moved = 0
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) for (rep_ in 1:12) {
		s = pm_state(family, bip)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, bip)
		out = lame:::.lame_nom_period_shift(s$Z, s$a, s$b, s$Sab, cap, masks, bip, rep(TRUE, s$nr), s$ic)
		if (out$calls == 0) next
		moved = moved + 1
		expect_equal(sum(out$a), sum(s$a))
		for (t in seq_len(s$nt)) {
			expect_equal(
				out$Z[, , t] - pm_eta(s, out$ic$path, out$a, t),
				s$Z[, , t] - pm_eta(s, s$path, s$a, t),
				tolerance = 1e-10
			)
			expect_true(pm_consistent(out$Z[, , t], s$Y[, , t], s$cap, family, bip))
		}
		# the rows with room (pinned) keep their latent cells wherever the move left their predictor
		bd = lame:::.lame_nom_shift_bounds_t(s$Z, cap, masks)
		pinned = is.finite(bd$hi)
		for (t in seq_len(s$nt)) for (i in which(pinned[, t])) {
			expect_equal(out$Z[i, , t] - s$Z[i, , t], rep(0, s$nc), tolerance = 1e-12)
		}
	}
	expect_gt(moved, 20)
})

test_that('rows that are not free to move keep their level and tie every period together', {
	skip_on_cran()
	set.seed(7104)
	moved = 0
	for (family in c('cbin', 'frn')) for (rep_ in 1:15) {
		s = pm_state(family, FALSE)
		cap = rep(s$cap, s$nr)
		masks = lame:::.lame_nom_shift_masks(s$Y, FALSE)
		ok = rep(TRUE, s$nr)
		ok[c(2, 5)] = FALSE
		out = lame:::.lame_nom_period_shift(s$Z, s$a, s$b, s$Sab, cap, masks, FALSE, ok, s$ic)
		if (out$calls == 0) next
		moved = moved + 1
		# the rows left out are pinned in every period: one component, and their cells stay where they were
		expect_identical(out$comps, 1L)
		expect_identical(out$Z[c(2, 5), , ], s$Z[c(2, 5), , ])
		expect_equal(sum(out$a), sum(s$a))
		for (t in seq_len(s$nt)) {
			expect_equal(
				out$Z[, , t] - pm_eta(s, out$ic$path, out$a, t),
				s$Z[, , t] - pm_eta(s, s$path, s$a, t),
				tolerance = 1e-10
			)
			expect_true(pm_consistent(out$Z[, , t], s$Y[, , t], s$cap, family, FALSE))
		}
	}
	expect_gt(moved, 5)
})

# one move from the joint law, the paired differences of moments before and after (Hotelling T2)
pm_geweke = function(move, family, bip, ndraw, seed, rho = 0, nt = 3L, nr = 9L) {
	set.seed(seed)
	TT = nt
	stat_of = function(s, path, a, Z) {
		mm = function(x) if (length(x)) mean(x) else 0
		c(
			path,
			path^2,
			path[-1] * path[-TT],
			a,
			a^2,
			path[1] * a,
			path[TT] * a,
			vapply(seq_len(TT), function(t) mean(Z[, , t]), 0),
			vapply(seq_len(TT), function(t) mean(Z[, , t]^2), 0),
			vapply(seq_len(TT), function(t) mm(Z[, , t][!is.na(s$Y[, , t]) & s$Y[, , t] == 0]), 0),
			vapply(seq_len(TT), function(t) mm(Z[, , t][!is.na(s$Y[, , t]) & s$Y[, , t] > 0]), 0),
			vapply(
				seq_len(TT),
				function(t) { zz = Z[, , t]; if (bip) 0 else mean(zz[upper.tri(zz)] * t(zz)[upper.tri(zz)]) },
				0
			),
			rowMeans(Z[, , 1]),
			rowMeans(Z[, , TT]),
			as.vector(Z[1:3, 1:2, 1])
		)
	}
	D = vector('list', ndraw)
	for (k in seq_len(ndraw)) {
		s = pm_state(family, bip, nt = nt, nr = nr, rho = rho, m0 = -0.8, v0 = 0.5, sd2 = 0.16)
		out = move(s)
		D[[k]] = stat_of(s, out$path, out$a, out$Z) - stat_of(s, s$path, s$a, s$Z)
	}
	D = do.call(rbind, D)
	N = nrow(D)
	mu = colMeans(D)
	sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), p = pchisq(T2, sum(kp), lower.tail = FALSE), change = mean(abs(D[, seq_len(TT)])))
}
pm_move_period = function(s) {
	cap = rep(s$cap, s$nr)
	out = lame:::.lame_nom_period_shift(
		s$Z,
		s$a,
		s$b,
		s$Sab,
		cap,
		lame:::.lame_nom_shift_masks(s$Y, s$bip),
		s$bip,
		rep(TRUE, s$nr),
		s$ic
	)
	list(Z = out$Z, a = out$a, path = out$ic$path)
}

test_that('one period translation leaves the joint law of path, levels, latent array and data unchanged', {
	skip_on_cran()
	for (cfg in list(list('cbin', FALSE, 3101), list('frn', TRUE, 3102))) {
		r = pm_geweke(pm_move_period, cfg[[1]], cfg[[2]], 4000, cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
		expect_gt(r$change, 0.01)
	}
})

test_that('the exactness check fails for a translation that leaves the path out of the density', {
	skip_on_cran()
	src = deparse(lame:::.lame_nom_period_shift)
	src = sub('sum(du * (path[-1L] - ic$rho * path[-N]))/ic$sd2', '0', src, fixed = TRUE)
	broken = eval(parse(text = paste(src, collapse = '\n')), envir = asNamespace('lame'))
	mv = function(s) {
		cap = rep(s$cap, s$nr)
		out = broken(s$Z, s$a, s$b, s$Sab, cap, lame:::.lame_nom_shift_masks(s$Y, s$bip), s$bip, rep(TRUE, s$nr), s$ic)
		list(Z = out$Z, a = out$a, path = out$ic$path)
	}
	r = pm_geweke(mv, 'cbin', FALSE, 4000, 3103)
	expect_gt(r$maxz, 5)
})

# the collapsed translation

pm_move_room = function(s, fun = lame:::.lame_nom_room_period) {
	cx = lame:::.lame_nom_room_context(s$Y, rep(s$cap, s$nr), s$bip)
	out = fun(s$Z, s$EZ, cx, s$bip, rep(s$rho, s$nt), s$ic)
	list(Z = out$Z, a = s$a, path = out$ic$path)
}

test_that('the collapsed translation keeps the data, the levels and the predictor', {
	skip_on_cran()
	set.seed(7102)
	moved = 0
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) for (rho in c(0, 0.5)) for (rep_ in 1:10) {
		if (bip && rho != 0) next
		s = pm_state(family, bip, rho = rho)
		cx = lame:::.lame_nom_room_context(s$Y, rep(s$cap, s$nr), bip)
		out = lame:::.lame_nom_room_period(s$Z, s$EZ, cx, bip, rep(rho, s$nt), s$ic)
		if (out$stats[['accepts']] == 0) next
		moved = moved + 1
		for (t in seq_len(s$nt)) {
			# the predictor the move carries is that of the shifted state, the levels did not move
			expect_equal(out$EZ[, , t], pm_eta(s, out$ic$path, s$a, t), tolerance = 1e-10)
			expect_true(pm_consistent(out$Z[, , t], s$Y[, , t], s$cap, family, bip))
			# every cell that is not integrated out keeps its residual: the slice translates with the predictor
			w = matrix(FALSE, s$nr, s$nc)
			for (i in which(cx[[t]]$pinned)) w[i, cx[[t]]$Wcols[[i]]] = TRUE
			expect_equal(((out$Z - out$EZ)[, , t] - (s$Z - s$EZ)[, , t])[!w], rep(0, sum(!w)), tolerance = 1e-9)
		}
		expect_identical(out$stats[['stuck']], 0)
	}
	expect_gt(moved, 20)
	# a panel with no row with room has nothing to integrate: the state comes back unchanged
	s = pm_state('cbin', FALSE, cap = 1)
	s$Y[, , ] = ifelse(is.na(s$Y), NA, 1)
	cx = lame:::.lame_nom_room_context(s$Y, rep(1, s$nr), FALSE)
	expect_false(any(vapply(cx, function(ct) any(ct$pinned), NA)))
	out = lame:::.lame_nom_room_period(s$Z, s$EZ, cx, FALSE, rep(0, s$nt), s$ic)
	expect_identical(out$Z, s$Z)
	expect_identical(out$ic$path, s$ic$path)
})

test_that('a period with a row whose hidden cells leave its room open is left alone', {
	set.seed(7103)
	s = pm_state('cbin', FALSE, nt = 3L, nr = 9L, cap = 3)
	# hide one cell of a row with room in period 2 (a row below its cap with an observed non-nomination)
	y = s$Y[, , 2]
	i = which(rowSums(!is.na(y) & y > 0) < 3 & rowSums(!is.na(y) & y == 0) > 1)[1]
	skip_if(is.na(i))
	j = which(!is.na(y[i, ]) & y[i, ] == 0)[1]
	s$Y[i, j, 2] = NA
	cx = lame:::.lame_nom_room_context(s$Y, rep(3, s$nr), FALSE)
	expect_true(cx[[2]]$unsure[i])
	expect_false(cx[[2]]$pinned[i])
	out = lame:::.lame_nom_room_period(s$Z, s$EZ, cx, FALSE, rep(0, s$nt), s$ic)
	expect_identical(out$Z[, , 2], s$Z[, , 2])
	expect_identical(out$ic$path[2], s$ic$path[2])
})

test_that('one collapsed translation leaves the joint law of path, levels, latent array and data unchanged', {
	skip_on_cran()
	for (cfg in list(list('cbin', FALSE, 0, 3201), list('frn', FALSE, 0.4, 3202), list('cbin', TRUE, 0, 3203))) {
		r = pm_geweke(function(s) pm_move_room(s), cfg[[1]], cfg[[2]], 4000, cfg[[4]], rho = cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
		expect_gt(r$change, 0.01)
	}
})

test_that('the exactness check fails for a collapsed translation that drops the integrated cells from its density', {
	skip_on_cran()
	src = deparse(lame:::.lame_nom_room_logf)
	src = sub('if (any(keep)) ', 'if (FALSE) ', src, fixed = TRUE)
	broken = eval(parse(text = paste(src, collapse = '\n')), envir = asNamespace('lame'))
	move = lame:::.lame_nom_room_period
	environment(move) = list2env(list(.lame_nom_room_logf = broken), parent = asNamespace('lame'))
	r = pm_geweke(function(s) pm_move_room(s, move), 'cbin', FALSE, 4000, 3204)
	expect_gt(r$maxz, 6)
})

# the fits

pm_fit = function(family = 'cbin', prior = list(), bip = FALSE, R = 0, dyn = 'intercept', ...) {
	set.seed(5)
	nr = 14L
	nc = if (bip) 9L else nr
	nt = 4L
	w = cumsum(rnorm(nt, 0, 0.3))
	Y = X = vector('list', nt)
	nm_r = sprintf('a%02d', seq_len(nr))
	nm_c = if (bip) sprintf('b%02d', seq_len(nc)) else nm_r
	a = rnorm(nr, 0, 0.7)
	b = rnorm(nc, 0, 0.7)
	for (t in seq_len(nt)) {
		x1 = matrix(rnorm(nr * nc), nr, nc)
		Zg = -1.2 + w[t] + 0.6 * x1 + outer(a, b, '+') + matrix(rnorm(nr * nc), nr, nc)
		if (!bip) diag(Zg) = -Inf
		y = pm_observe(Zg, 3, family, bip)
		dimnames(y) = list(nm_r, nm_c)
		Y[[t]] = y
		X[[t]] = array(x1, c(nr, nc, 1), dimnames = list(nm_r, nm_c, 'x1'))
	}
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	suppressWarnings(suppressMessages(lame(
		Y,
		Xdyad = X,
		family = family,
		odmax = 3,
		mode = if (bip) 'bipartite' else 'unipartite',
		R = R,
		dynamic_beta = dyn,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 3,
		prior = prior,
		...
	)))
}

test_that('the period and room translations run in the fits and are switched off by their own switches', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) for (bip in c(FALSE, TRUE)) {
		fit = pm_fit(family, bip = bip)
		expect_true(all(is.finite(fit$BETA)))
		am = fit$actor_move
		expect_true(
			all(c('period_calls', 'period_comps', 'room_attempts', 'room_accepts', 'room_stuck') %in% names(am))
		)
		expect_gt(am[['period_calls']], 0)
		expect_gt(am[['room_attempts']], 0)
		expect_lte(am[['room_accepts']], am[['room_attempts']])
		expect_identical(am[['room_stuck']], 0)
		expect_false('period_calls' %in% names(pm_fit(family, list(actor_move_period = FALSE), bip = bip)$actor_move))
		expect_false(any(grepl('^room_', names(pm_fit(family, list(actor_move_room = FALSE), bip = bip)$actor_move))))
		expect_null(pm_fit(family, list(actor_move = FALSE), bip = bip)$actor_move)
	}
	expect_error(pm_fit('cbin', list(actor_move_period = 1)), 'actor_move_period must be TRUE or FALSE', fixed = TRUE)
	expect_error(pm_fit('cbin', list(actor_move_room = NA)), 'actor_move_room must be TRUE or FALSE', fixed = TRUE)
})

test_that('the translations run with dynamic_ab and a rank, and not with a static intercept', {
	skip_on_cran()
	fit = pm_fit('frn', dynamic_ab = TRUE)
	expect_true(all(is.finite(fit$BETA)))
	expect_gt(fit$actor_move[['period_calls']], 0)
	fit = pm_fit('cbin', R = 1)
	expect_gt(fit$actor_move[['period_calls']], 0)
	# a static intercept has no path to translate period by period
	fit = pm_fit('cbin', dyn = 'dyad')
	expect_false(any(c('period_calls', 'room_attempts') %in% names(fit$actor_move)))
})
