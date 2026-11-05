# translations of the dynamic_beta coefficients of a Poisson fit with the latent array, and the dyadic variance and
# correlation moved with the residuals carried along their conditional law (R/poisson_period_move.R)
skip_on_cran()

# a copy of the function `name` of the package with the text `from` replaced by `to` (one pair or several, applied in turn), on the lines of
# its deparsed source that contain `only` and not `skip`
mut_fn = function(name, from, to, only = NULL, skip = NULL) {
	out = deparse(get(name, asNamespace('lame')), width.cutoff = 500L)
	for (j in seq_along(from)) {
		hit = grepl(from[j], out, fixed = TRUE)
		if (!is.null(only)) hit = hit & grepl(only, out, fixed = TRUE)
		if (!is.null(skip)) hit = hit & !grepl(skip, out, fixed = TRUE)
		if (!any(hit)) stop('mutation target not found: ', from[j])
		out[hit] = sub(from[j], to[j], out[hit], fixed = TRUE)
	}
	fn = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(fn) = asNamespace('lame')
	fn
}

# the path prior

test_that('the path prior of a coefficient is the law the samplers draw from, for every kind and gap pattern', {
	set.seed(3)
	ns = asNamespace('lame')
	nt = 7L
	for (kind in c('ar1', 'rw1', 'rw2', 'matern32')) for (gaps_on in c(FALSE, TRUE)) {
		pos = if (gaps_on) c(1, 2, 4, 5, 8, 9, 10) else seq_len(nt)
		rho = 0.7; sig = 0.4; lam = 1.3; m0 = -0.5; v0 = 2
		pri = ns$.lame_ppm_path_prior(kind, nt, rho, sig, lam, pos, NULL, m0, v0)
		quad = function(x) -0.5 * sum(x * (pri$P %*% x)) + sum(pri$b * x)
		ld = function(x) {
			if (kind %in% c('ar1', 'rw1')) {
				z = ns$.lame_beta_gap_transition(if (kind == 'rw1') 1 else rho, diff(pos))
				-0.5 * (x[1] - m0)^2 / v0 - 0.5 * sum((x[-1] - z$f * x[-nt])^2 / (z$q * sig^2 * lam))
			} else {
				K = ns$.dynamic_beta_kind_to_prec(kind, nt, NULL, pos)
				-0.5 * sum(x * (K %*% x)) / (sig^2 * lam)
			}
		}
		x = rnorm(nt); y = rnorm(nt)
		expect_equal(quad(x) - quad(y), ld(x) - ld(y), tolerance = 1e-9)
	}
})

# the model of the checks

# one state of the model: dynamic paths (the intercept and a dyadic coefficient), a static coefficient with its g-prior, levels, dyadic
# residuals (variance s2, correlation rho), counts from the Poisson law of the latent log rate
pp_state = function(kind, mode, gaps_on, seed, n = 8L, nt = 4L, nc = NULL, p = 2L, asym = FALSE) {
	set.seed(seed)
	ns = asNamespace('lame')
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	nr = n; nc = if (bip) max(3L, n - 2L) else n
	tpos = if (gaps_on) cumsum(c(1, sample(c(1, 2, 3), nt - 1L, replace = TRUE))) else seq_len(nt)
	sigma = c(0.25, 0.15); lam = c(1, 0.8); m0 = c(-1.4, 0.3); v0 = c(0.5, 0.4); rho_b = 0.8
	s2 = 0.3; g = 2.5; rho = if (identical(mode, 'dir')) 0.3 else 0
	# (asym: the design of a symmetric network is not symmetric, as with a row or column covariate; the latent value of a pair has the
	# mean of the two directed predictors and the predictor of the sweep keeps the directed cells)
	Xd = lapply(
		seq_len(nt),
		function(t) { x1 = matrix(rnorm(nr * nc), nr, nc); if (sym && !asym) x1 = (x1 + t(x1)) / 2; cbind(1, as.numeric(x1)) }
	)
	Xs = lapply(
		seq_len(nt),
		function(t) { x2 = matrix(rnorm(nr * nc), nr, nc); if (sym && !asym) x2 = (x2 + t(x2)) / 2; matrix(as.numeric(x2), ncol = 1) }
	)
	Xdz = if (sym) lapply(Xd, function(m) ns$.lame_ppm_sym_cols(m, nr)) else Xd
	Xsz = if (sym) lapply(Xs, function(m) ns$.lame_ppm_sym_cols(m, nr)) else Xs
	path = matrix(0, nt, p)
	for (k in seq_len(p)) {
		if (kind %in% c('ar1', 'rw1')) {
			z = ns$.lame_beta_gap_transition(if (kind == 'rw1') 1 else rho_b, diff(tpos))
			path[1, k] = rnorm(1, m0[k], sqrt(v0[k]))
			for (t in 2:nt) path[t, k] = z$f[t - 1L] * path[t - 1L, k] + rnorm(
				1,
				0,
				sigma[k] * sqrt(lam[k] * z$q[t - 1L])
			)
		} else {
			pri = ns$.lame_ppm_path_prior(kind, nt, rho_b, sigma[k], lam[k], tpos, NULL, m0[k], v0[k])
			path[, k] = backsolve(chol(pri$P), rnorm(nt))
		}
	}
	bst = rnorm(1, 0, sqrt(g * s2)); a = rnorm(nr, 0, 0.3); b = rnorm(nc, 0, 0.3)
	off = rep(0, nt)
	Z = EZ = Y = array(NA_real_, c(nr, nc, nt))
	for (t in seq_len(nt)) {
		eta = matrix(Xd[[t]] %*% path[t, ], nr, nc) + bst * matrix(Xs[[t]][, 1], nr, nc) + outer(a, b, '+')
		mu = if (sym) (eta + t(eta)) / 2 else eta
		e = matrix(rnorm(nr * nc, 0, sqrt(s2)), nr, nc)
		if (identical(mode, 'dir')) { e2 = matrix(
			rnorm(nr * nr, 0, sqrt(s2)),
			nr,
			nr
		); lo = lower.tri(e); e[lo] = (rho * t(e) + sqrt(1 - rho^2) * e2)[lo] }
		if (sym) e[lower.tri(e)] = t(e)[lower.tri(e)]
		zt = mu + e
		y = matrix(rpois(nr * nc, exp(zt + off[t])), nr, nc)
		if (!bip) diag(y) = NA
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		Z[, , t] = zt; EZ[, , t] = eta; Y[, , t] = y
	}
	pr = list(
		kind = kind,
		rho = rep(rho_b, p),
		sigma = sigma,
		lam = lam,
		m0 = m0,
		v0 = v0,
		time_positions = tpos,
		length_scale = NULL
	)
	list(
		path = path,
		bst = bst,
		a = a,
		b = b,
		Z = Z,
		EZ = EZ,
		Y = Y,
		Xd = Xd,
		Xs = Xs,
		Xdz = Xdz,
		Xsz = Xsz,
		pr = pr,
		nr = nr,
		nc = nc,
		nt = nt,
		mode = mode,
		s2 = s2,
		rho = rho,
		g = g
	)
}

test_that('a translation keeps every residual and the predictor consistent with the coefficient state', {
	ns = asNamespace('lame')
	for (cfg in list(
		list('ar1', 'dir', FALSE, FALSE),
		list('rw2', 'bip', FALSE, FALSE),
		list('matern32', 'sym', TRUE, FALSE),
		list('rw1', 'dir', TRUE, FALSE),
		list('ar1', 'sym', FALSE, TRUE)
	)) {
		s = pp_state(cfg[[1]], cfg[[2]], cfg[[3]], 5, asym = cfg[[4]])
		cx = ns$.lame_ppm_context(s$Y, s$mode)
		out = ns$.lame_ppm_visit(cx, s$Z, s$EZ, s$path, s$Xd, s$pr, Xz = s$Xdz)
		expect_gt(out$stats[['period_calls']], 0)
		expect_gt(out$stats[['path_calls']], 0)
		# the residual of a pair is its latent value less the mean of the two directed predictors
		bar = function(a) if (identical(cfg[[2]], 'sym')) { for (t in seq_len(dim(a)[3])) a[, , t] = (a[, , t] + t(a[, , t])) / 2; a } else a
		expect_equal(out$Z - bar(out$EZ), s$Z - bar(s$EZ), tolerance = 1e-10)
		for (t in seq_len(s$nt)) {
			shift = matrix(s$Xd[[t]] %*% (out$path[t, ] - s$path[t, ]), s$nr, s$nc)
			expect_equal(out$EZ[, , t], s$EZ[, , t] + shift, tolerance = 1e-10)
			if (identical(cfg[[2]], 'sym')) expect_equal(out$Z[, , t], t(out$Z[, , t]), tolerance = 1e-12)
		}
		expect_gt(max(abs(out$path - s$path)), 0)
	}
})

# exactness

# paired differences of moments before and after one move from the joint law
pp_test = function(D) {
	N = nrow(D); mu = colMeans(D); sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), T2 = T2, df = sum(kp), p = pchisq(T2, sum(kp), lower.tail = FALSE))
}

pp_geweke_visit = function(
	kind,
	mode,
	gaps_on,
	ndraw,
	seed,
	visit = lame:::.lame_ppm_visit,
	ctx = lame:::.lame_ppm_context,
	static = TRUE,
	asym = FALSE
) {
	ns = asNamespace('lame')
	stat_of = function(s) {
		z = s$Z; y = s$Y; nt = s$nt
		c(
			as.vector(s$path),
			as.vector(s$path)^2,
			s$path[-1L, 1] * s$path[-nt, 1],
			s$path[, 1] * s$path[, 2],
			s$bst,
			s$bst^2,
			s$bst * s$path[, 1],
			sapply(seq_len(nt), function(t) mean(z[, , t])),
			sapply(seq_len(nt), function(t) mean(z[, , t]^2)),
			sapply(seq_len(nt), function(t) mean(y[, , t] * z[, , t], na.rm = TRUE)),
			sapply(seq_len(nt), function(t) mean(z[, , t] * s$Xdz[[t]][, 2]))
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s = pp_state(kind, mode, gaps_on, seed * 100000L + it, asym = asym)
		cx = ctx(s$Y, s$mode)
		out = visit(cx, s$Z, s$EZ, s$path, s$Xd, s$pr, Xz = s$Xdz)
		s1 = s; s1$Z = out$Z; s1$path = out$path
		if (it <= 3) {
			# the translation changes no residual (a pair's: its latent value less the mean of the two directed predictors); the predictor follows the coefficients
			bar = function(a) if (identical(mode, 'sym')) { for (t in seq_len(dim(a)[3])) a[, , t] = (a[, , t] + t(a[, , t])) / 2; a } else a
			if (max(abs((out$Z - bar(out$EZ)) - (s$Z - bar(s$EZ)))) > 1e-8) stop('a residual changed')
			for (t in seq_len(s$nt)) if (max(
				abs(out$EZ[, , t] - s$EZ[, , t] - matrix(s$Xd[[t]] %*% (out$path[t, ] - s$path[t, ]), s$nr, s$nc))
			) > 1e-8) stop('the predictor does not follow the coefficients')
		}
		if (static) {
			Xs_arr = array(vapply(seq_len(s$nt), function(t) s$Xs[[t]][, 1], numeric(s$nr * s$nc)), c(s$nr, s$nc, s$nt))
			Xsz_arr = array(
				vapply(seq_len(s$nt), function(t) s$Xsz[[t]][, 1], numeric(s$nr * s$nc)),
				c(s$nr, s$nc, s$nt)
			)
			sc = ns$.lame_ppm_static(
				cx,
				out$Z,
				out$EZ,
				s$bst,
				list(Xs_arr),
				matrix(1 / (s$g * s$s2)),
				Xz = list(Xsz_arr)
			)
			s1$Z = sc$Z; s1$bst = sc$beta
		}
		D[[it]] = stat_of(s1) - stat_of(s)
	}
	pp_test(do.call(rbind, D))
}

test_that('one translation sweep leaves the joint law of paths, latent array and counts unchanged', {
	skip_on_cran()
	# (rw2 has a nearly flat prior in two directions whose draws are beyond any Poisson rate: its terms are checked against the density above)
	for (cfg in list(
		list('ar1', 'dir', FALSE, 11),
		list('rw1', 'sym', TRUE, 12),
		list('ar1', 'bip', FALSE, 13),
		list('matern32', 'dir', TRUE, 14)
	)) {
		r = pp_geweke_visit(cfg[[1]], cfg[[2]], cfg[[3]], 2000, cfg[[4]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	# a symmetric network whose design is not symmetric (a row or column covariate): the pair moves by the design averaged over it
	r = pp_geweke_visit('ar1', 'sym', FALSE, 2000, 15, asym = TRUE)
	expect_lt(r$maxz, 4.5)
	expect_gt(r$p, 1e-3)
})

test_that('the exactness check fails for translations that do not carry the latent slice with the coefficient', {
	skip_on_cran()
	broken = mut_fn('.lame_ppm_visit', 'Z <- Z + sh', 'Z <- Z')
	expect_error(pp_geweke_visit('ar1', 'dir', FALSE, 50, 21, visit = broken, static = FALSE), 'a residual changed')
	broken = mut_fn(
		'.lame_ppm_visit', 'drop(X[[t]][, moved, drop = FALSE] %*% shifted[t, moved])', 'rep(sum(shifted[t, moved]), nrc)'
	)
	expect_error(
		pp_geweke_visit('ar1', 'dir', FALSE, 50, 22, visit = broken, static = FALSE),
		'the predictor does not follow the coefficients'
	)
	broken = mut_fn('.lame_ppm_visit', 'shifted[t, k] <<- shifted[t, k] + d', 'shifted[t, k] <<- d')
	expect_error(
		pp_geweke_visit('ar1', 'dir', FALSE, 50, 23, visit = broken, static = FALSE),
		'the predictor does not follow the coefficients'
	)
	broken = mut_fn('.lame_ppm_path_prior', 'P[1L, 1L] <- P[1L, 1L] + 1/v0', 'P[1L, 1L] <- P[1L, 1L] + 0')
	env = new.env(parent = asNamespace('lame'))
	broken_visit = lame:::.lame_ppm_visit
	environment(broken_visit) = env
	assign('.lame_ppm_path_prior', broken, envir = env)
	r = pp_geweke_visit('rw1', 'dir', FALSE, 2000, 22, visit = broken_visit, static = FALSE)
	expect_gt(r$maxz, 4.5)
})

# the scale move: (s2, rho) from the prior, the residuals from the model, counts from the Poisson law of the latent log rate
pp_geweke_scale = function(
	mode,
	ndraw,
	seed,
	move = lame:::.lame_ppm_scale,
	miss = 0,
	n = 8L,
	nt = 3L,
	diag_t = FALSE,
	asym = FALSE
) {
	ns = asNamespace('lame')
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym'); two = identical(mode, 'dir')
	nr = n; nc = if (bip) max(3L, n - 2L) else n
	shape0 = 3; rate0 = 0.4; rho_pow = ns$.lame_rho_prior_power()
	set.seed(seed)
	EZ = array(NA_real_, c(nr, nc, nt))
	for (t in seq_len(nt)) { m = matrix(
		rnorm(nr * nc, -0.6, 0.8),
		nr,
		nc
	); if (sym && !asym) m = (m + t(m)) / 2; EZ[, , t] = m }
	# the mean of the latent value of a pair is the mean of the two directed predictors, which differ with a row or column covariate
	EZbar = EZ
	if (sym) for (t in seq_len(nt)) EZbar[, , t] = (EZ[, , t] + t(EZ[, , t])) / 2
	up = upper.tri(matrix(0, nr, nc))
	stat_of = function(s) {
		e = s$Z - EZbar
		per = sapply(
			seq_len(nt),
			function(t) { et = e[, , t]; c(
				mean(et^2),
				if (bip) 0 else mean(et[up] * t(et)[up]),
				mean(s$Y[, , t] * et, na.rm = TRUE),
				mean(et)
			) }
		)
		c(
			log(s$s2),
			s$s2,
			s$rho,
			s$rho^2,
			s$rho * log(s$s2),
			as.vector(per),
			if (diag_t && !bip) sapply(
				seq_len(nt),
				function(t) mean(diag(e[, , t])^2) / (s$s2 * (1 + if (two) s$rho else 0))
			)
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s2 = 1 / rgamma(1, shape0, rate0)
		rho = if (two) 2 * rbeta(1, rho_pow + 1, rho_pow + 1) - 1 else 0
		Z = Y = array(NA_real_, c(nr, nc, nt))
		for (t in seq_len(nt)) {
			e = matrix(rnorm(nr * nc, 0, sqrt(s2)), nr, nc)
			if (two) { e2 = matrix(
				rnorm(nr * nr, 0, sqrt(s2)),
				nr,
				nr
			); lo = lower.tri(e); e[lo] = (rho * t(e) + sqrt(1 - rho^2) * e2)[lo] }
			if (sym) e[lower.tri(e)] = t(e)[lower.tri(e)]
			if (!bip) diag(e) = rnorm(nr, 0, sqrt(if (diag_t && two) s2 * (1 + rho) else s2))
			zt = EZbar[, , t] + e
			y = matrix(rpois(nr * nc, exp(zt)), nr, nc)
			if (miss > 0) { mm = matrix(runif(nr * nc) < miss, nr, nc); if (sym) mm[lower.tri(mm)] = t(mm)[lower.tri(mm)]; y[mm] = NA }
			if (!bip) diag(y) = NA
			if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
			Z[, , t] = zt; Y[, , t] = y
		}
		s = list(s2 = s2, rho = rho, Z = Z, Y = Y)
		cx = ns$.lame_ppm_context(Y, if (bip) 'bip' else if (sym) 'sym' else 'dir')
		out = if (diag_t) move(
			cx,
			Z,
			EZ,
			s2,
			rho,
			c(shape0 = shape0, rate0 = rate0),
			two,
			ns$.lame_rho_prior_log,
			tries = 4L,
			diag_in_target = TRUE
		) else move(cx, Z, EZ, s2, rho, c(shape0 = shape0, rate0 = rate0), two, ns$.lame_rho_prior_log, tries = 4L)
		s1 = s; s1$Z = out$Z; s1$s2 = out$s2; s1$rho = out$rho
		D[[it]] = stat_of(s1) - stat_of(s)
	}
	pp_test(do.call(rbind, D))
}

test_that('the carried residuals leave the joint law of the dyadic variance and correlation, residuals and counts unchanged', {
	skip_on_cran()
	for (cfg in list(list('dir', 31, 0), list('dir', 32, 0.15), list('sym', 33, 0.15), list('bip', 34, 0))) {
		r = pp_geweke_scale(cfg[[1]], 4000, cfg[[2]], miss = cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	# a symmetric network whose directed predictors differ (a row or column covariate in the design): the pair's mean is their average
	r = pp_geweke_scale('sym', 4000, 38, miss = 0.15, asym = TRUE)
	expect_lt(r$maxz, 4.5)
	expect_gt(r$p, 1e-3)
	r = pp_geweke_scale(
		'dir',
		5000,
		35,
		move = mut_fn('.lame_ppm_scale', '(nw$logdet - st$logdet) +', '0 +', only = 'rho_logprior')
	)
	expect_gt(r$maxz, 6)
	r = pp_geweke_scale('dir', 5000, 36, move = mut_fn('.lame_ppm_scale', '(nw$g - st$g) +', '', only = 'rho_logprior'))
	expect_gt(r$maxz, 6)
	# the conditional means only set how large a step is accepted, not what the step leaves invariant
	r = pp_geweke_scale(
		'dir',
		4000,
		37,
		move = mut_fn('.lame_ppm_scale', c('bu <- yu - wu', 'bl <- yl - wl'), c('bu <- 0 * wu', 'bl <- 0 * wl'))
	)
	expect_lt(r$maxz, 4.5)
})

test_that('with the diagonal cells in the variance conditional (ame) the carried residuals leave the joint law unchanged, and the same check fails if they are left alone', {
	skip_on_cran()
	for (cfg in list(list('dir', 51, 0), list('dirind', 52, 0.15))) {
		r = pp_geweke_scale(cfg[[1]], 4000, cfg[[2]], miss = cfg[[3]], diag_t = TRUE)
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	blind = function(cx, Z, EZ, s2, rho, s2_prior, two, rho_logprior, tries = 4L, diag_in_target = TRUE) {
		lame:::.lame_ppm_scale(cx, Z, EZ, s2, rho, s2_prior, two, rho_logprior, tries = tries, diag_in_target = FALSE)
	}
	r = pp_geweke_scale('dir', 6000, 53, move = blind, diag_t = TRUE)
	expect_lt(r$p, 1e-6)
})

test_that('the carried residuals take larger steps than the plain scaling of the residuals', {
	skip_on_cran()
	# a panel with events: the plain carry's log density changes by about sum (y - lambda)^2 s2 (sqrt(c) - 1), the tilted one does not
	ns = asNamespace('lame')
	set.seed(8)
	n = 40L; nt = 2L; s2 = 0.05
	EZ = array(-2.5, c(n, n, nt))
	e = array(rnorm(n * n * nt, 0, sqrt(s2)), c(n, n, nt))
	Y = array(rpois(n * n * nt, exp(EZ + e)), c(n, n, nt))
	for (t in 1:nt) diag(Y[, , t]) = NA
	cx = ns$.lame_ppm_context(Y, 'dir')
	acc_new = acc_old = 0
	for (r in 1:30) {
		out = ns$.lame_ppm_scale(
			cx,
			EZ + e,
			EZ,
			s2,
			0.2,
			c(shape0 = 2, rate0 = 0.1),
			TRUE,
			ns$.lame_rho_prior_log,
			tries = 4L,
			steps = 0.3
		)
		acc_new = acc_new + out$stats[['scale_accepts']] / out$stats[['scale_attempts']]
		old = ns$.lame_ltm_scale(
			cx$ltm,
			EZ + e,
			EZ,
			s2,
			0.2,
			c(shape0 = 2, rate0 = 0.1),
			two = TRUE,
			rho_logprior = ns$.lame_rho_prior_log,
			steps = 0.3,
			tries = 4L
		)
		acc_old = acc_old + old$stats[['scale_accepts']] / old$stats[['scale_attempts']]
	}
	expect_gt(acc_new / 30, acc_old / 30)
})

# the switches and the fits

pp_fit = function(nscan = 30, mode = 'unipartite', kind = 'ar1', ...) {
	dat = simulate_test_network(n = 8, n_time = 4, family = 'poisson', mode = mode, nA = 8, nB = 6, seed = 11)
	suppressWarnings(suppressMessages(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = 'poisson',
		mode = mode,
		R = 0,
		dynamic_beta = c('intercept', 'dyad'),
		dynamic_beta_kind = kind,
		nscan = nscan,
		burn = 10,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5,
		...
	)))
}

test_that('the Poisson switches are validated', {
	for (nm in c('poisson_period_move', 'poisson_static_move')) {
		prior = list(1); names(prior) = nm
		expect_error(pp_fit(prior = prior), paste0(nm, ' must be TRUE or FALSE'), fixed = TRUE)
	}
})

test_that('the translations and the carried residuals run by default and are recorded; the switch omits the translations', {
	skip_on_cran()
	fit = pp_fit()
	rec = fit$poisson_period_move
	expect_true(!is.null(rec))
	expect_equal(unname(rec[['visits']]), 40)
	expect_gt(rec[['period_calls']], 0)
	expect_gt(rec[['path_calls']], 0)
	expect_gt(rec[['scale_attempts']], 0)
	expect_true(all(is.finite(fit$BETA)))
	only_scale = pp_fit(prior = list(poisson_period_move = FALSE))
	expect_equal(only_scale$poisson_period_move[['period_calls']], 0)
	expect_gt(only_scale$poisson_period_move[['scale_attempts']], 0)
})

test_that('every kind, mode and a static panel run with the moves and stay finite', {
	skip_on_cran()
	for (kind in c('ar1', 'rw1', 'rw2', 'matern32')) {
		fit = pp_fit(nscan = 20, kind = kind)
		expect_true(all(is.finite(fit$BETA)))
		expect_gt(fit$poisson_period_move[['period_calls']], 0)
	}
	fit = pp_fit(nscan = 20, mode = 'bipartite', kind = 'rw1')
	expect_true(all(is.finite(fit$BETA)))
	expect_gt(fit$poisson_period_move[['scale_attempts']], 0)
	# no dynamic_beta: the scale move alone
	dat = simulate_test_network(n = 8, n_time = 3, family = 'poisson', seed = 12)
	fit = suppressWarnings(suppressMessages(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = 'poisson',
		R = 0,
		nscan = 20,
		burn = 10,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5
	)))
	expect_equal(fit$poisson_period_move[['period_calls']], 0)
	expect_gt(fit$poisson_period_move[['scale_attempts']], 0)
})

test_that('the translations of the actor levels visit a poisson fit every second sweep by default, every tenth on request', {
	skip_on_cran()
	fit = pp_fit(nscan = 30)
	expect_equal(unname(fit$actor_move[['visits']]), 20)
	fit10 = pp_fit(nscan = 30, prior = list(actor_move_every = 10))
	expect_equal(unname(fit10$actor_move[['visits']]), 4)
})

test_that('the static coefficients of a dynamic_beta poisson fit are translated at the visits of the actor translations', {
	skip_on_cran()
	dat = simulate_test_network(n = 8, n_time = 4, family = 'poisson', mode = 'unipartite', seed = 11)
	fit_static = function(...) suppressWarnings(suppressMessages(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = 'poisson',
		R = 0,
		dynamic_beta = 'intercept',
		nscan = 30,
		burn = 10,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5,
		...
	)))
	p_static = dim(dat$Xdyad[[1]])[3]
	fit = fit_static()
	expect_equal(unname(fit$poisson_period_move[['static_calls']]), 20 * p_static)
	fit10 = fit_static(prior = list(actor_move_every = 10))
	expect_equal(unname(fit10$poisson_period_move[['static_calls']]), 4 * p_static)
	fit_off = fit_static(prior = list(poisson_static_move = FALSE))
	expect_equal(unname(fit_off$poisson_period_move[['static_calls']]), 0)
})

# the compiled kernels

test_that('the compiled carry of the dyads is the map of its definition, and its sums are the ones of the density', {
	set.seed(6)
	nd = 12L
	s2c = 0.21; rhoc = 0.35; s2n = 0.33; rhon = -0.2
	lam_u = rexp(nd, 3); lam_l = rexp(nd, 3)
	ou = rep(TRUE, nd); ol = rep(TRUE, nd); ol[2] = FALSE; ou[5] = FALSE
	yu = ifelse(ou, rpois(nd, lam_u), 0); yl = ifelse(ol, rpois(nd, lam_l), 0)
	wu = ifelse(ou, lam_u, 0); wl = ifelse(ol, lam_l, 0)
	bu = yu - wu; bl = yl - wl
	eu = rnorm(nd, 0, sqrt(s2c)); el = rnorm(nd, 0, sqrt(s2c))
	out = lame:::ppm_dyad_carry_cpp(eu, el, wu, wl, bu, bl, yu, yl, s2c, rhoc, s2n, rhon)
	law = function(i, s2, rho) {
		S = s2 * matrix(c(1, rho, rho, 1), 2)
		A = solve(S) + diag(c(wu[i], wl[i]))
		V = solve(A)
		list(m = drop(V %*% c(bu[i], bl[i])), L = t(chol(V)), V = V)
	}
	logdet = c(0, 0)
	for (i in seq_len(nd)) {
		c0 = law(i, s2c, rhoc); c1 = law(i, s2n, rhon)
		z = solve(c0$L, c(eu[i], el[i]) - c0$m)
		e1 = c1$m + drop(c1$L %*% z)
		expect_equal(c(out$eu[i], out$el[i]), e1, tolerance = 1e-10)
		logdet = logdet + c(log(prod(diag(c1$L))), log(prod(diag(c0$L))))
	}
	expect_equal(out$logdet, logdet[1], tolerance = 1e-10)
	cur = lame:::ppm_dyad_state_cpp(eu, el, wu, wl, yu, yl, s2c, rhoc)
	expect_equal(cur$logdet, logdet[2], tolerance = 1e-10)
	dens = function(e1, e2, s2, rho) -log(2 * pi) - log(s2) - 0.5 * log(1 - rho^2) - (e1^2 - 2 * rho * e1 * e2 + e2^2) / (2 * s2 * (1 - rho^2))
	ld = function(e1, e2, s2, rho) sum(dens(e1, e2, s2, rho))
	# the Gaussian term is the sum of the bivariate log densities up to the constant log(2 pi) per dyad
	expect_equal(out$g - nd * log(2 * pi), ld(out$eu, out$el, s2n, rhon), tolerance = 1e-10)
	expect_equal(cur$g - nd * log(2 * pi), ld(eu, el, s2c, rhoc), tolerance = 1e-10)
	pois = function(a, b) sum(yu * a - wu * exp(a) + yl * b - wl * exp(b))
	expect_equal(out$ll, pois(out$eu, out$el), tolerance = 1e-10)
	expect_equal(cur$ll, pois(eu, el), tolerance = 1e-10)
})

test_that('the compiled carry of single cells and the Poisson sum of a shift match their definitions', {
	set.seed(7)
	n = 15L
	w = rexp(n, 2); y = rpois(n, w); b = y - w; e = rnorm(n, 0, 0.5)
	s2c = 0.25; s2n = 0.4
	out = lame:::ppm_cell_carry_cpp(e, w, b, y, s2c, s2n)
	mc = b / (1 / s2c + w); mn = b / (1 / s2n + w)
	expect_equal(out$e, mn + sqrt((1 / s2c + w) / (1 / s2n + w)) * (e - mc), tolerance = 1e-12)
	expect_equal(out$logdet, -0.5 * sum(log(1 / s2n + w)), tolerance = 1e-12)
	expect_equal(lame:::ppm_cell_state_cpp(e, w, y, s2c)$logdet, -0.5 * sum(log(1 / s2c + w)), tolerance = 1e-12)
	ez = rexp(40); x = rnorm(40)
	expect_equal(lame:::ppm_shift_sum_cpp(ez, x, 0.37), sum(ez * expm1(0.37 * x)), tolerance = 1e-12)
	expect_error(lame:::ppm_shift_sum_cpp(ez, x[-1], 0.1))
})

test_that('the static translation of several coefficients equals its definition with the Poisson factor recomputed from the latent array', {
	ns = asNamespace('lame')
	reference = function(cx, Z, EZ, beta, Xk, iV0) {
		lx = cx$ltm
		for (k in sample.int(length(beta))) {
			x = Xk[[k]]
			xo = x[lx$ok_idx]
			ez = exp(Z[lx$ok_idx] + lx$off_ok)
			ysx = sum(lx$y_ok * xo)
			ll = function(dl) dl * ysx - sum(ez * expm1(dl * xo))
			a1 = drop(iV0 %*% beta)[k]
			a2 = iV0[k, k]
			d = ns$.lame_ess_step(0, -a1 / a2, matrix(1 / sqrt(a2)), ll, 0)$x
			Z = Z + d * x
			EZ = EZ + d * x
			beta[k] = beta[k] + d
		}
		list(Z = Z, EZ = EZ, beta = beta)
	}
	for (mode in c('dir', 'sym', 'bip')) {
		s = pp_state('ar1', mode, FALSE, 17)
		cx = ns$.lame_ppm_context(s$Y, s$mode)
		x1 = array(vapply(seq_len(s$nt), function(t) s$Xs[[t]][, 1], numeric(s$nr * s$nc)), c(s$nr, s$nc, s$nt))
		x2 = array(rnorm(length(x1)), dim(x1))
		x3 = array(rnorm(length(x1)), dim(x1))
		if (identical(mode, 'sym')) for (t in seq_len(s$nt)) { x2[, , t] = (x2[, , t] + t(x2[, , t])) / 2; x3[, , t] = (x3[, , t] + t(x3[, , t])) / 2 }
		beta = c(0.2, -0.3, 0.1)
		iV0 = diag(c(0.5, 0.8, 1.2))
		set.seed(5)
		new = ns$.lame_ppm_static(cx, s$Z, s$EZ, beta, list(x1, x2, x3), iV0)
		set.seed(5)
		ref = reference(cx, s$Z, s$EZ, beta, list(x1, x2, x3), iV0)
		expect_gt(new$stats[['coef_calls']], 2)
		expect_gt(max(abs(new$beta - beta)), 0)
		expect_equal(new$beta, ref$beta, tolerance = 1e-9)
		expect_equal(new$Z, ref$Z, tolerance = 1e-9)
		expect_equal(new$EZ, ref$EZ, tolerance = 1e-9)
		# every residual is unchanged
		expect_equal(new$Z - new$EZ, s$Z - s$EZ, tolerance = 1e-10)
	}
})

# ame()

test_that('ame() Poisson fits move the dyadic variance and correlation with the carried residuals, record it, and the switch omits it', {
	skip_on_cran()
	for (mode in c('dir', 'sym', 'bip')) {
		dat = simulate_test_network(
			n = 9,
			n_time = 1,
			family = 'poisson',
			mode = if (identical(mode, 'bip')) 'bipartite' else 'unipartite',
			nA = 9,
			nB = 6,
			seed = 12
		)
		y = dat$Y; x = dat$Xdyad
		if (identical(mode, 'sym')) { y[lower.tri(y)] = t(y)[lower.tri(y)]; for (k in seq_len(dim(x)[3])) { m = x[, , k]; m[lower.tri(m)] = t(m)[lower.tri(m)]; x[, , k] = m } }
		fit_one = function(...) suppressWarnings(suppressMessages(ame(
			y,
			Xdyad = x,
			family = 'poisson',
			symmetric = identical(mode, 'sym'),
			mode = if (identical(mode, 'bip')) 'bipartite' else 'unipartite',
			dcor = identical(mode, 'dir'),
			burn = 10,
			nscan = 30,
			odens = 1,
			seed = 5,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			...
		)))
		fit = fit_one()
		expect_equal(unname(fit$poisson_period_move[['visits']]), 40)
		expect_gt(fit$poisson_period_move[['scale_attempts']], 0)
		expect_gt(fit$actor_move[['scale_attempts']], 0)
		expect_true(all(is.finite(fit$BETA)) && all(is.finite(fit$VC)))
		expect_equal(unname(fit$actor_move[['visits']]), 20)
		expect_equal(unname(fit_one(prior = list(actor_move_every = 10))$actor_move[['visits']]), 4)
	}
})
