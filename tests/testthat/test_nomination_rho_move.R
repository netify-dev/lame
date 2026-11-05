# collapsed reciprocity update of the cbin / frn samplers (.lame_nom_rho_collapsed, R/nomination_actor_move.R): rho is drawn
# with the non-nominations and the cells no nomination can sit in integrated out, and those cells are then redrawn given the nominations
# checked: density terms against numerical integration, calibration of the density, what the update keeps, kernel exactness and the fits

# the latent array and data of a draw from the model: dyads with correlation rho (diagonal variance
# 1 + rho), row i nominating its min(cap, number of positive cells) strongest cells, a share of the cells
# hidden and optionally an actor absent from the first period
rm_draw = function(family, n = 6L, nt = 2L, cap = 2, hidden = 0, absent = FALSE, pw = 1, rho = NULL, mu = -1.1) {
	if (is.null(rho)) rho = if (pw == 1) 2 * rbeta(1, 2, 2) - 1 else -cos(pi * runif(1))
	a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
	Z = EZ = Y = array(NA_real_, c(n, n, nt))
	YL = vector('list', nt)
	for (t in seq_len(nt)) {
		eta = mu + outer(a, b, '+') + 0.3 * t / nt
		e1 = matrix(rnorm(n * n), n, n); e2 = matrix(rnorm(n * n), n, n)
		E = e1; low = row(E) > col(E)
		E[low] = (rho * t(e1) + sqrt(1 - rho^2) * e2)[low]
		diag(E) = rnorm(n, 0, sqrt(1 + rho))
		gone = if (absent && t == 1L) n else integer(0)
		z = eta + E
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			zi = z[i, ]; zi[i] = -Inf; zi[gone] = -Inf
			pos = which(zi > 0); k = min(cap, length(pos))
			if (k > 0) { top = pos[order(zi[pos], decreasing = TRUE)[seq_len(k)]]; y[i, top] = if (family == 'cbin') 1 else rank(zi[top]) }
		}
		diag(y) = NA
		if (length(gone)) { y[gone, ] = NA; y[, gone] = NA }
		if (hidden > 0) { cand = which(!is.na(y)); y[cand[runif(length(cand)) < hidden]] = NA }
		if (family == 'frn') for (i in seq_len(n)) { r = which(!is.na(y[i, ]) & y[i, ] > 0); if (length(r)) y[i, r] = rank(y[i, r]) }
		if (family == 'frn') {
			mx = max(c(0, y), na.rm = TRUE); YLt = matrix(NA_integer_, n, max(mx, 1))
			for (i in seq_len(n)) for (r in seq_len(mx)) { j = which(y[i, ] == r); if (length(j)) YLt[i, r] = j }
			YL[[t]] = YLt
		}
		Z[, , t] = z; EZ[, , t] = eta; Y[, , t] = y
	}
	list(rho = rho, Z = Z, EZ = EZ, Y = Y, YL = if (family == 'frn') YL else NULL, n = n, nt = nt, cap = cap)
}

test_that('the terms of the collapsed density are the integrals over the non-nominations', {
	skip_on_cran()
	set.seed(7701)
	for (rep in 1:20) {
		r = runif(1, -0.9, 0.9); c1 = rnorm(1, 1, 1.5); c2 = rnorm(1, 1, 1.5); e = rnorm(1)
		# two non-nominations: P(e1 < c1, e2 < c2) by integration over the first
		ff = integrate(function(x) dnorm(x) * pnorm((c2 - r * x) / sqrt(1 - r^2)), -Inf, c1, rel.tol = 1e-10)$value
		expect_equal(pbivnorm::pbivnorm(c1, c2, r), ff, tolerance = 1e-8)
		# a non-nomination of a nomination's dyad: P(e2 < c2 | e1 = e) is the normal probability of the standardised slack
		bf = integrate(function(x) dnorm((x - r * e) / sqrt(1 - r^2)) / sqrt(1 - r^2), -Inf, c2, rel.tol = 1e-10)$value
		expect_equal(pnorm((c2 - r * e) / sqrt(1 - r^2)), bf, tolerance = 1e-8)
	}
})

test_that('the density of the correlation is calibrated: its distribution function at the true value is uniform', {
	skip_on_cran()
	set.seed(7702)
	grid = seq(-0.9, 0.9, length.out = 721)
	for (cfg in list(
		list('cbin', 1, 0, FALSE),
		list('frn', 1, 0.15, TRUE),
		list('cbin', -0.5, 0.1, TRUE),
		list('frn', -0.5, 0, FALSE)
	)) {
		us = numeric(0)
		while (length(us) < 350) {
			s = rm_draw(cfg[[1]], hidden = cfg[[3]], absent = cfg[[4]], pw = cfg[[2]])
			if (abs(s$rho) >= 0.9) next
			d = lame:::.lame_nom_rho_dyads(s$Z, s$EZ, s$Y, cfg[[1]], s$YL, rep(s$cap, s$n))
			lf = vapply(grid, function(r) lame:::.lame_nom_rho_logf(r, d, cfg[[2]]), 0)
			cdf = cumsum(exp(lf - max(lf))); cdf = cdf / cdf[length(cdf)]
			us = c(us, approx(grid, cdf, xout = s$rho, rule = 2)$y)
		}
		expect_gt(suppressWarnings(ks.test(us, 'punif')$p.value), 1e-3)
	}
})

test_that('the update keeps the nominations and the data and redraws only what it integrated out', {
	skip_on_cran()
	set.seed(7703)
	moved = 0
	# mu: a sparse predictor and a dense one (most rows at their cap, non-nominations far from their ceiling and some above their mean)
	for (family in c('cbin', 'frn')) for (case in 1:3) for (rep in 1:30) {
		s = rm_draw(
			family,
			hidden = if (case >= 2) 0.15 else 0,
			absent = case == 3,
			mu = if (rep %% 3 == 0) 0.9 else -1.1
		)
		out = lame:::.lame_nom_rho_collapsed(s$Z, s$EZ, s$rho, s$Y, family, s$YL, rep(s$cap, s$n), pw = 1)
		if (out$accepts == 0L) { expect_identical(out$Z, s$Z); expect_identical(out$rho, s$rho); next }
		moved = moved + 1
		expect_false(out$rho == s$rho)
		for (t in seq_len(s$nt)) {
			st = lame:::.lame_nom_col_slice(s$Z[, , t], s$Y[, , t], family, s$YL[[t]], rep(s$cap, s$n), FALSE)
			# the nominations (observed or hidden, with their labels) are unchanged
			expect_identical(out$Z[, , t][st$isnom], s$Z[, , t][st$isnom])
			# the data match the latent array, with the same labels, and the ceilings hold
			st1 = lame:::.lame_nom_col_slice(out$Z[, , t], s$Y[, , t], family, s$YL[[t]], rep(s$cap, s$n), FALSE)
			expect_identical(st1$cls, st$cls)
			hid = lame:::.lame_nom_hidden_cells(s$Y[, , t], FALSE)
			for (i in seq_len(s$n)) {
				obs = !is.na(s$Y[i, , t]); obs[i] = FALSE
				el = obs | hid[i, ]
				z = out$Z[i, , t]; z[!el] = -Inf
				pos = which(z > 0); k = min(s$cap, length(pos))
				top = if (k > 0) pos[order(z[pos], decreasing = TRUE)[seq_len(k)]] else integer(0)
				expect_true(all((seq_along(z) %in% top)[obs] == (s$Y[i, obs, t] > 0)))
			}
		}
	}
	expect_gt(moved, 60)
})

# one update from the joint law of (correlation, latent array, data) must keep that law, so paired differences of
# any function of the state before and after have mean zero (the z values returned)
rm_joint_law = function(family, pw, hidden, absent, reps = 2500, seed = 1, n = 6L) {
	set.seed(seed)
	up = upper.tri(matrix(0, n, n))
	funs = function(rho, Z, EZ) {
		E = Z - EZ
		c(
			rho,
			rho^2,
			rho^3,
			mean(vapply(1:2, function(t) mean(E[, , t][up] * t(E[, , t])[up]), 1)),
			mean(vapply(1:2, function(t) mean(E[, , t][up]^2 + t(E[, , t])[up]^2), 1)),
			mean(vapply(1:2, function(t) mean(diag(E[, , t])^2), 1)),
			mean(diag(E[, , 1])),
			as.vector(Z[1:3, 1:3, 1]),
			as.vector(Z[1:3, 1:3, 1])^2,
			as.numeric(as.vector(Z[, , 1]) > 0)[1:30],
			rho * Z[1, 2, 1],
			rho * E[2, 1, 1] * E[1, 2, 1],
			rowMeans(Z[, , 2]),
			# the dyads of the actor that may be absent from the first period (cells no nomination can sit in)
			E[n, 2, 1] * E[2, n, 1],
			E[n, 3, 1] * E[3, n, 1],
			E[n, 2, 1]^2,
			E[2, n, 1]^2,
			rho * E[n, 2, 1] * E[2, n, 1]
		)
	}
	d = NULL
	for (r in seq_len(reps)) {
		s = rm_draw(family, n = n, hidden = hidden, absent = absent, pw = pw)
		o = lame:::.lame_nom_rho_collapsed(s$Z, s$EZ, s$rho, s$Y, family, s$YL, rep(s$cap, n), pw = pw)
		d = rbind(d, funs(o$rho, o$Z, s$EZ) - funs(s$rho, s$Z, s$EZ))
	}
	keep = apply(d, 2, sd) > 1e-12
	colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrow(d)))
}

test_that('the collapsed update leaves the joint law of correlation, latent values and data unchanged', {
	skip_on_cran()
	cfg = list(
		list('cbin', 1, 0, FALSE),
		list('frn', 1, 0.15, TRUE),
		list('cbin', -0.5, 0.15, FALSE),
		list('frn', -0.5, 0, TRUE)
	)
	for (k in seq_along(cfg)) {
		z = rm_joint_law(cfg[[k]][[1]], cfg[[k]][[2]], cfg[[k]][[3]], cfg[[k]][[4]], seed = 900 + k)
		expect_lt(max(abs(z)), 4.6)
	}
	# the check has power: an update that draws the correlation under a flat prior fails it
	logf = get('.lame_nom_rho_logf', asNamespace('lame'))
	testthat::local_mocked_bindings(.package = 'lame', .lame_nom_rho_logf = function(r, d, pw) logf(r, d, 0))
	expect_gt(max(abs(rm_joint_law('cbin', 1, 0, FALSE, reps = 2500, seed = 95))), 5.5)
})

test_that('lame() and ame() run the reciprocity update for cbin and frn and record it', {
	skip_on_cran()
	for (family in c('cbin', 'frn')) {
		set.seed(7710)
		n = 14L; nt = 3L
		nm = sprintf('a%02d', 1:n)
		Y = X = vector('list', nt)
		a = rnorm(n, 0, .6); b = rnorm(n, 0, .6)
		for (t in 1:nt) {
			x = matrix(rnorm(n * n), n, n)
			s = rm_draw(family, n = n, nt = 1L, cap = 2, hidden = 0.1, rho = 0.2)
			y = s$Y[, , 1]; dimnames(y) = list(nm, nm); Y[[t]] = y
			X[[t]] = array(x, c(n, n, 1L), dimnames = list(nm, nm, 'x1'))
		}
		names(Y) = names(X) = as.character(2000 + 1:nt)
		run = function(prior, kind = 'lame', ...) suppressMessages(suppressWarnings(
			if (kind == 'lame') lame(
				Y,
				Xdyad = X,
				family = family,
				odmax = 2,
				R = 0,
				burn = 20,
				nscan = 80,
				odens = 5,
				seed = 3,
				verbose = FALSE,
				plot = FALSE,
				gof = FALSE,
				prior = prior,
				...
			)
			else ame(
				Y[[1]],
				Xdyad = X[[1]],
				family = family,
				odmax = 2,
				R = 0,
				burn = 20,
				nscan = 80,
				odens = 5,
				seed = 3,
				verbose = FALSE,
				gof = FALSE,
				prior = prior,
				...
			)
		))
		for (kind in c('lame', 'ame')) {
			on = run(list(), kind)
			# a visit every tenth sweep: 100 sweeps, 10 visits of 4 proposals
			expect_equal(unname(on$actor_move['rho_attempts']), 40)
			expect_gt(unname(on$actor_move['rho_accepts']), 0)
			expect_true(all(is.finite(on$VC)) && all(abs(on$VC[, 'rho']) < 1))
			# the switches: all moves, the reciprocity update alone, a schedule that is never reached
			expect_null(run(list(actor_move = FALSE), kind)$actor_move)
			expect_false('rho_attempts' %in% names(run(list(actor_move_rho = FALSE), kind)$actor_move))
			never = run(list(actor_move_every = 10^6), kind); off = run(list(actor_move = FALSE), kind)
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		}
		# without a dyadic correlation there is nothing to update
		expect_false('rho_attempts' %in% names(run(list(), 'ame', dcor = FALSE)$actor_move))
	}
})
