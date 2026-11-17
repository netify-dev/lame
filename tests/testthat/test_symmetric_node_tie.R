# the row and column coefficients of a nodal covariate in a symmetric binary or ordinal fit are tied after every draw

stn_sym = function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m }

# small symmetric panel with a time-varying nodal covariate, a dyadic covariate and dynamic additive effects
stn_panel = function(fam = 'binary', n = 30L, Tn = 4L, seed = 1L) {
	set.seed(500 + seed)
	rn = sprintf('a%02d', seq_len(n))
	lev = rnorm(n, 0, .7)
	dev = matrix(0, n, Tn)
	dev[, 1] = rnorm(n, 0, .4)
	for (t in 2:Tn) dev[, t] = .7 * dev[, t - 1] + rnorm(n, 0, .4 * sqrt(1 - .49))
	xs = rnorm(n)
	zz = lapply(seq_len(Tn), function(t) {
		xr = sqrt(.5) * xs + sqrt(.5) * rnorm(n)
		xd = stn_sym(matrix(rnorm(n * n), n, n))
		at = lev + dev[, t]
		list(
			xr = xr,
			xd = xd,
			z = .6 * xd + .35 * outer(xr, xr, '+') + outer(at, at, '+') + stn_sym(matrix(rnorm(n * n), n, n))
		)
	})
	b0 = -unname(quantile(unlist(lapply(zz, function(o) o$z[upper.tri(o$z)])), if (fam == 'poisson') .9 else .92))
	Y = X = Xr = vector('list', Tn)
	for (t in seq_len(Tn)) {
		zt = zz[[t]]$z + b0
		y = if (fam == 'binary') 1 * (zt > 0) else matrix(rpois(n * n, exp(zt - 1)), n, n)
		y = stn_sym(y)
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		Y[[t]] = y
		X[[t]] = array(zz[[t]]$xd, c(n, n, 1), dimnames = list(rn, rn, 'xd'))
		Xr[[t]] = matrix(zz[[t]]$xr, n, 1, dimnames = list(rn, 'xn'))
	}
	names(Y) = names(X) = names(Xr) = paste0('p', seq_len(Tn))
	list(Y = Y, Xdyad = X, Xrow = Xr)
}
stn_fit = function(d, fam, prior = list(), burn = 10L, nscan = 40L, R = 0L, ...) {
	suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$Xdyad,
		Xrow = d$Xrow,
		symmetric = TRUE,
		family = fam,
		R = R,
		burn = burn,
		nscan = nscan,
		odens = 1L,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

test_that('the tie sets the two coefficients of a nodal covariate to their mean and leaves the others', {
	b = c(-1, .2, .8, .5, 1.5)
	expect_equal(.lame_sym_tie(b, 2L, 3L), c(-1, .5, .5, .5, 1.5))
	expect_equal(.lame_sym_tie(b, integer(0), integer(0)), b)
	expect_equal(.lame_sym_tie(c(.2, .8, 1), 1L, 2L), c(.5, .5, 1))
	# two nodal covariates: row slices first, then column slices, as in the design
	expect_equal(.lame_sym_tie(c(0, .1, .3, .9, .5, 2), 2:3, 4:5), c(0, .5, .4, .5, .4, 2))
	# both slices dynamic: paths are tied period by period
	path = cbind(c(1, 2, 3), c(3, 2, 1))
	st = .lame_sym_tie_state(path, c(.5, 1.5), c(2L, 3L), c(1L, 4L), 2L, 3L)
	expect_equal(st$path, cbind(c(2, 2, 2), c(2, 2, 2)))
	expect_equal(st$stat, c(.5, 1.5))
	# both slices static: static entries are tied and a dynamic path is left alone
	st = .lame_sym_tie_state(matrix(1:3, 3, 1), c(.2, .8, 1), 1L, 2:4, 2L, 3L)
	expect_equal(st$stat, c(.5, .5, 1))
	expect_equal(st$path, matrix(1:3, 3, 1))
})


test_that('the exact symmetric draw cannot separate the two slices of a nodal covariate, and the tied draw has the law of their mean', {
	skip_on_cran()
	set.seed(3)
	n = 8L; Tn = 2L; p = 4L
	xr = rnorm(n)
	X = lapply(seq_len(Tn), function(t) {
		a = array(0, c(n, n, p))
		a[, , 1] = 1
		a[, , 2] = matrix(xr, n, n)
		a[, , 3] = matrix(xr, n, n, byrow = TRUE)
		a[, , 4] = stn_sym(matrix(rnorm(n * n), n, n))
		a
	})
	Q = crossprod(matrix(rnorm(16), 4)) / 20 + diag(.05, 4)
	cx = .lame_symmetric_gaussian_context(X, Q, actor = FALSE, center = FALSE)
	# the symmetrised design has the same column for the row and column slices
	expect_equal(cx$D[, 2], cx$D[, 3])
	Z = array(unlist(lapply(seq_len(Tn), function(t) stn_sym(matrix(rnorm(n * n), n, n)))), c(n, n, Tn))
	# posterior of the coefficients given Z, from its definition
	y = unlist(lapply(seq_len(Tn), function(t) (Z[, , t][cx$index] + Z[, , t][cx$mirror]) / 2))
	Lam = crossprod(cx$D) + Q
	Sig = solve(Lam)
	m = drop(Sig %*% crossprod(cx$D, y))
	Tm = diag(p); Tm[2:3, 2:3] = .5
	mu = drop(Tm %*% m)
	V = Tm %*% Sig %*% t(Tm)
	draw = function(tie) replicate(20000, { b = .lame_symmetric_gaussian_draw(Z, 1, cx)$beta; tie(b) })
	one = draw(function(b) .lame_sym_tie(b, 2L, 3L))
	expect_equal(one[2, ], one[3, ])
	emp = rowMeans(one)
	se = sqrt(diag(V) / 20000 * 3)
	expect_lt(max(abs(emp - mu) / se), 4.5)
	expect_lt(max(abs(cov(t(one)) - V) / sqrt(outer(diag(V), diag(V)))), .03)
	# a tie that keeps one of the two rather than their mean has a different law
	bad = draw(function(b) { b[3] = b[2]; b })
	expect_gt(max(abs(rowMeans(bad) - mu) / se), 6)
})


test_that('the tie applies to the binary and ordinal families with a nodal covariate, and to nothing else', {
	skip_on_cran()
	d = stn_panel('binary', n = 12L, Tn = 3L, seed = 2L)
	set.seed(2)
	o = d
	o$Y = lapply(
		d$Y,
		function(y) { v = y[upper.tri(y)] + rnorm(sum(upper.tri(y)), 0, .3); cu = quantile(v, c(.8, .9, .96)); m = matrix(NA_real_, nrow(y), ncol(y)); m[upper.tri(m)] = findInterval(v, cu); m[lower.tri(m)] = t(m)[lower.tri(m)]; dimnames(m) = dimnames(y); m }
	)
	pd = stn_panel('poisson', n = 12L, Tn = 3L, seed = 2L)
	calls = 0L
	local_mocked_bindings(
		.lame_sym_tie = function(beta, rb, cb) { calls <<- calls + 1L; beta },
		.lame_sym_tie_state = function(
			path,
			stat,
			dyn_idx,
			stat_idx,
			rb,
			cb
		) { calls <<- calls + 1L; list(path = path, stat = stat) },
		.package = 'lame'
	)
	count = function(...) { calls <<- 0L; f = stn_fit(...); expect_true(all(is.finite(f$BETA))); calls }
	# binary and ordinal fits with dynamic additive effects, positions or coefficients
	expect_gt(count(d, 'binary', dynamic_ab = TRUE), 0)
	expect_gt(count(o, 'ordinal', dynamic_ab = TRUE), 0)
	expect_gt(count(d, 'binary', dynamic_uv = TRUE, R = 1L), 0)
	expect_gt(count(d, 'binary', dynamic_beta = c('intercept', 'dyad')), 0)
	expect_gt(count(d, 'binary', dynamic_beta = TRUE), 0)
	# the switch, the normal and Poisson families, and a fit without a nodal covariate
	expect_equal(count(d, 'binary', dynamic_ab = TRUE, prior = list(symmetric_node_tie = FALSE)), 0)
	expect_equal(count(pd, 'poisson', dynamic_ab = TRUE), 0)
	dn = d
	dn$Y = lapply(
		d$Y,
		function(y) { m = y + matrix(
			rnorm(length(y), 0, .5),
			nrow(y)
		); m = (m + t(m)) / 2; diag(m) = NA; dimnames(m) = dimnames(y); m }
	)
	expect_equal(count(dn, 'normal', dynamic_ab = TRUE), 0)
	d$Xrow = NULL
	expect_equal(count(d, 'binary', dynamic_ab = TRUE), 0)
})

test_that('the switch is validated in lame()', {
	d = stn_panel('binary', n = 10L, Tn = 2L, seed = 3L)
	expect_error(stn_fit(d, 'binary', prior = list(symmetric_node_tie = 'no')), 'symmetric_node_tie')
	f = stn_fit(d, 'binary', prior = list(symmetric_node_tie = FALSE))
	expect_true(all(is.finite(f$BETA)))
})

test_that('dynamic additive effects of a symmetric binary or ordinal fit with a nodal covariate stay at their scale', {
	skip_on_cran()
	d = stn_panel('binary', n = 30L, Tn = 4L, seed = 1L)
	# innovation sd of the paths in the data is 0.28
	on = stn_fit(d, 'binary', burn = 300L, nscan = 1700L, dynamic_ab = TRUE)
	expect_lt(mean(tail(on$sigma_ab, 100)), 1)
	# every period of the draws x coefficients x periods array
	expect_lt(max(abs(on$BETA[, 'xn_node', ])), 2)
	# without the tie the static coefficient walks along its ridge and the scale runs away; the chain
	# policy refuses that arm, so the fit comes from the condition. both arms drop the intercept
	on0 = stn_fit(d, 'binary', burn = 300L, nscan = 1700L, dynamic_ab = TRUE, intercept = FALSE)
	expect_lt(mean(tail(on0$sigma_ab, 100)), 1)
	expect_lt(max(abs(on0$BETA[, 'xn_node'])), 2)
	off = tryCatch(
		stn_fit(
			d,
			'binary',
			prior = list(symmetric_node_tie = FALSE, chain_extensions = 0),
			burn = 300L,
			nscan = 1700L,
			dynamic_ab = TRUE,
			intercept = FALSE
		),
		lame_chain_unsettled = function(e) e$fit
	)
	expect_gt(mean(tail(off$sigma_ab, 100)), 3)
	# the same for the ordinal family
	set.seed(2)
	o = d
	o$Y = lapply(
		d$Y,
		function(y) { v = y[upper.tri(y)] + rnorm(sum(upper.tri(y)), 0, .3); cu = quantile(v, c(.8, .9, .96)); m = matrix(NA_real_, nrow(y), ncol(y)); m[upper.tri(m)] = findInterval(v, cu); m[lower.tri(m)] = t(m)[lower.tri(m)]; dimnames(m) = dimnames(y); m }
	)
	on = stn_fit(o, 'ordinal', burn = 300L, nscan = 1200L, dynamic_ab = TRUE)
	expect_lt(mean(tail(on$sigma_ab, 100)), 1.5)
})

# the symmetric node tie composed with the scale moves of the binary and ordinal families: one step of the move on a tied state
# leaves the joint law of the tied coefficients, the effects, the latent array and the ties unchanged

# tie_geweke: states are drawn from the prior of the tied model (the dyadic coefficient, the pair of a nodal covariate as the mean of a draw
# of the pair, static additive effects, an unscaled part m standing for dynamic effects), the latent array around the predictor, and the ties
# (binary) from its signs; `visit` is one step of the move; the scaled parameters are multiplied by the factor it returns.
# dcount: the dimension of the scaled parameters the move is told (the tied pair counts once);
# untie: a broken composition that scales only the row slice of the pair
tie_geweke = function(
	kind,
	ndraw,
	seed,
	withm = TRUE,
	dcount = NULL,
	untie = FALSE,
	visit = NULL,
	n = 6L,
	nt = 2L,
	miss = 0.1
) {
	ns = asNamespace('lame')
	set.seed(seed)
	nr = n; nc = n; p = 3L; vb = 0.8; va = 0.5
	X = array(rnorm(nr * nc * p * nt), c(nr, nc, p, nt))
	for (t in 1:nt) {
		d = X[, , 1, t]; d[lower.tri(d)] = t(d)[lower.tri(d)]; X[, , 1, t] = d
		xr = rnorm(nr)
		X[, , 2, t] = matrix(xr, nr, nc)
		X[, , 3, t] = matrix(xr, nr, nc, byrow = TRUE)
	}
	mfix = if (withm) array(rnorm(nr * nc * nt, 0.5, 0.7), c(nr, nc, nt)) else NULL
	if (!is.null(mfix)) for (t in 1:nt) { m_ = mfix[, , t]; m_[lower.tri(m_)] = t(m_)[lower.tri(m_)]; mfix[, , t] = m_ }
	upm = upper.tri(matrix(0, nr, nc))
	E_of = function(s) {
		E = array(0, c(nr, nc, nt))
		for (t in 1:nt) E[, , t] = outer(s$a, s$a, '+') + matrix(
			matrix(X[, , , t], nr * nc, p) %*% s$beta,
			nr,
			nc
		) + if (is.null(mfix)) 0 else mfix[, , t]
		E
	}
	draw_state = function() {
		pr = rnorm(2, 0, sqrt(vb)); mm = mean(pr)
		beta = c(rnorm(1, 0, sqrt(vb)), mm, mm)
		a = rnorm(nr, 0, sqrt(va))
		E = E_of(list(a = a, beta = beta))
		Z = E
		for (t in 1:nt) { z = matrix(rnorm(nr * nc), nr, nc); z[lower.tri(z)] = t(z)[lower.tri(z)]; Z[, , t] = Z[, , t] + z }
		Y = NULL
		if (kind == 'binary') {
			Y = 1 * (Z > 0)
			ms = array(runif(nr * nc * nt) < miss, c(nr, nc, nt))
			for (t in 1:nt) ms[, , t][lower.tri(ms[, , t])] = t(ms[, , t])[lower.tri(ms[, , t])]
			Y[ms] = NA
			for (t in 1:nt) diag(Y[, , t]) = NA
		}
		list(beta = beta, a = a, E = E, Z = Z, Y = Y)
	}
	stat_of = function(s) {
		z = s$Z; eb = s$E
		msk = array(rep(upm, nt), dim(z)); zz = z * msk; ee = eb * msk
		y1 = if (is.null(s$Y)) 0 * z else ifelse(is.na(s$Y), 0, s$Y) * msk
		c(
			s$beta,
			s$beta^2,
			s$beta[1] * s$beta[2],
			sum(s$a^2),
			sum(s$a),
			(s$beta[2] - s$beta[3])^2,
			abs(s$beta[2] - s$beta[3]),
			sapply(
				1:nt,
				function(t) c(
					mean(zz[, , t]),
					mean(zz[, , t]^2),
					mean((zz[, , t] - ee[, , t])^2),
					mean(zz[, , t] * ee[, , t]),
					mean(y1[, , t]),
					mean(y1[, , t] * ee[, , t])
				)
			),
			sum(y1 * ee) * s$beta[1],
			mean(zz^2) * s$beta[2]
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s = draw_state()
		# the priors of the scaled parameters at the tied state: the pair enters through the exchangeable covariance, 2 m^2 / vb
		Apri = sum(s$beta^2) / vb + sum(s$a^2) / va
		dth = (if (is.null(dcount)) 2L else dcount) + nr
		s1 = s
		if (kind == 'binary') {
			cx = ns$.lame_bso_context(s$Y, 'sym')
			o = (if (is.null(visit)) ns$.lame_bso_visit else visit)(s$Z, s$E, mfix, cx, 0, Apri, 0, dth)
		} else {
			o = (if (is.null(visit)) ns$.lame_osm_visit else visit)(s$Z, s$E, mfix, 'sym', rep(0, nt), 1, Apri, 0, dth, affine = TRUE, orbit = TRUE)
		}
		s1$Z = o$Z
		s1$beta = o$c * s$beta
		if (untie) s1$beta[3] = s$beta[3]
		s1$a = o$c * s$a
		s1$E = E_of(s1)
		D[[it]] = stat_of(s1) - stat_of(s)
	}
	D = do.call(rbind, D)
	N = nrow(D); mu = colMeans(D); sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), T2 = T2, df = sum(kp), p = pchisq(T2, sum(kp), lower.tail = FALSE))
}

test_that('the scale moves of the binary and ordinal families keep a tied pair tied and leave the joint law unchanged when the pair counts as one dimension', {
	skip_on_cran()
	for (cfg in list(
		list('binary', 81, TRUE),
		list('binary', 82, FALSE),
		list('ordinal', 83, TRUE),
		list('ordinal', 84, FALSE)
	)) {
		r = tie_geweke(cfg[[1]], 7000, cfg[[2]], withm = cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('counting the pair as two dimensions, or leaving it untied after the move, fails the check', {
	skip_on_cran()
	for (kd in c('binary', 'ordinal')) {
		expect_gt(tie_geweke(kd, 12000, 85, dcount = 3L)$maxz, 4.5)
		expect_gt(tie_geweke(kd, 7000, 86, untie = TRUE)$maxz, 5)
	}
})

test_that('a fit tells the scale moves that a tied pair is one dimension', {
	skip_on_cran()
	stn_p = function(fam, n = 12L, Tn = 3L) {
		set.seed(7)
		rn = sprintf('a%02d', seq_len(n))
		sym = stn_sym
		Y = lapply(
			seq_len(Tn),
			function(t) { z = sym(matrix(rnorm(n * n), n)) + 0.3; y = if (fam == 'binary') 1 * (z > 0.5) else matrix(
				findInterval(z, c(0, .7, 1.4)),
				n
			); y = sym(y); diag(y) = NA; dimnames(y) = list(rn, rn); y }
		)
		X = lapply(
			seq_len(Tn),
			function(t) array(sym(matrix(rnorm(n * n), n)), c(n, n, 1), dimnames = list(rn, rn, 'xd'))
		)
		Xr = lapply(seq_len(Tn), function(t) matrix(rnorm(n), n, 1, dimnames = list(rn, 'xn')))
		names(Y) = names(X) = names(Xr) = paste0('p', seq_len(Tn))
		list(Y = Y, X = X, Xr = Xr)
	}
	got = new.env()
	got$d = numeric(0)
	local_mocked_bindings(
		.lame_bso_visit = function(
			Z,
			E,
			m,
			cx,
			rho,
			A_pri,
			B_pri,
			d_th,
			...
		) { got$d = c(got$d, d_th); list(Z = Z, c = 1, evals = 0L) },
		.lame_osm_visit = function(
			Z,
			E,
			m,
			mode,
			rho,
			s2,
			A_pri,
			B_pri,
			d_th,
			...
		) { got$d = c(got$d, d_th); list(Z = Z, c = 1) },
		.package = 'lame'
	)
	first_d = function(d, fam, tie, ...) {
		got$d = numeric(0)
		suppressWarnings(suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			Xrow = d$Xr,
			symmetric = TRUE,
			family = fam,
			R = 0L,
			burn = 3,
			nscan = 6,
			odens = 1,
			seed = 2,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(symmetric_node_tie = tie, binary_scale_every = 1),
			...
		)))
		got$d[1L]
	}
	for (fam in c('binary', 'ordinal')) {
		d = stn_p(fam)
		expect_equal(first_d(d, fam, FALSE) - first_d(d, fam, TRUE), 1)
		# a tied pair of ar1 paths is one path (three periods) and one long-run mean
		expect_equal(first_d(d, fam, FALSE, dynamic_beta = TRUE) - first_d(d, fam, TRUE, dynamic_beta = TRUE), 4)
		expect_equal(
			first_d(d, fam, FALSE, dynamic_beta = TRUE, dynamic_beta_kind = "rw1") -
				first_d(d, fam, TRUE, dynamic_beta = TRUE, dynamic_beta_kind = "rw1"),
			3
		)
	}
})
