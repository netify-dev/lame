# scale of the linear predictor of a binary probit fit (R/binary_scale_orbit.R): the observed-cell index sets,
# the slice step, the compiled log likelihood, exactness of the move, and the switch and records of the fits
# on a panel that a covariate separates

bso_mut = function(name, from, to) {
	src = deparse(get(name, asNamespace('lame')), width.cutoff = 500L)
	out = sub(from, to, src, fixed = TRUE)
	if (identical(out, src)) stop('mutation target not found: ', from)
	fn = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(fn) = asNamespace('lame')
	fn
}

# the building blocks

test_that('the context lists the observed off-diagonal cells, and a symmetric pair once with the value of its upper cell, else its lower cell', {
	ns = asNamespace('lame')
	Y = array(NA_real_, c(3, 3, 2)); Y[, , 1] = matrix(c(NA, 1, 0, 0, NA, 1, 1, 0, NA), 3); Y[, , 2] = matrix(
		c(NA, NA, 1, 1, NA, 0, 0, 0, NA),
		3
	)
	cx = ns$.lame_bso_context(Y, 'dir')
	expect_equal(length(cx$idx), 6 + 5)
	expect_equal(cx$sg, 2 * Y[cx$idx] - 1)
	expect_equal(length(cx$idg), 6)
	expect_equal(length(cx$ina), 1)
	expect_equal(sort(c(cx$i1, cx$i0)), cx$idx)
	cb = ns$.lame_bso_context(array(c(1, 0, 1, NA, 0, 1), c(2, 3, 1)), 'bip')
	expect_equal(length(cb$idx), 5)
	expect_equal(length(cb$idg), 0)
	# the pair {1, 2} is observed in its lower cell only, {1, 3} in its upper cell, {2, 3} in both
	Ys = array(matrix(c(NA, 1, NA, NA, NA, 0, 1, 0, NA), 3), c(3, 3, 1))
	cs = ns$.lame_bso_context(Ys, 'sym')
	expect_equal(cs$up, c(4, 7, 8))
	expect_equal(cs$lo, c(2, 3, 6))
	expect_equal(cs$sg, c(1, 1, -1))
})

test_that('the slice step leaves a density invariant', {
	skip_on_cran()
	ns = asNamespace('lame')
	set.seed(3)
	m0 = 0.3; s0 = 0.2
	y = replicate(6000, {
		x0 = rnorm(1, m0, s0)
		x0 + ns$.lame_bso_slice(function(u) -(u + x0 - m0)^2 / (2 * s0^2))$u
	})
	expect_equal(mean(y), m0, tolerance = 0.02)
	expect_equal(sd(y), s0, tolerance = 0.05)
	# a skewed density (log of a gamma)
	y = replicate(6000, {
		x0 = log(rgamma(1, 4, 2))
		x0 + ns$.lame_bso_slice(function(u) 4 * (u + x0) - 2 * exp(u + x0), w = 0.3)$u
	})
	expect_equal(mean(y), digamma(4) - log(2), tolerance = 0.03)
	expect_equal(sd(y), sqrt(trigamma(4)), tolerance = 0.05)
	expect_equal(ns$.lame_bso_slice(function(u) NA_real_)$u, 0)
})

test_that('the compiled log likelihood is the sum of the probit log probabilities, and cells far on the correct side add nothing', {
	ns = asNamespace('lame')
	set.seed(4)
	es = rnorm(500, 0, 2); mm = rnorm(500, 0.3, 1); sg = sample(c(-1, 1), 500, replace = TRUE)
	for (k in c(0.4, 1, 2.5)) expect_equal(
		ns$bso_loglik_cpp(es, mm, sg, k, 1 / 0.8),
		sum(pnorm(sg * (k * es + mm) / 0.8, log.p = TRUE)),
		tolerance = 1e-9
	)
	big = c(es, rep(40, 50)); mb = c(mm, rep(0, 50)); sb = c(sg, rep(1, 50))
	expect_equal(
		ns$bso_loglik_cpp(big, mb, sb, 1.3, 1),
		sum(pnorm(sb * (1.3 * big + mb), log.p = TRUE)),
		tolerance = 1e-12
	)
	expect_equal(
		ns$bso_loglik_cpp(c(-5, 5), c(0, 0), c(1, 1), 1, 1),
		pnorm(-5, log.p = TRUE) + pnorm(5, log.p = TRUE),
		tolerance = 1e-12
	)
	expect_error(ns$bso_loglik_cpp(1:3 + 0, 1:2 + 0, c(1, 1, 1), 1, 1))
})

# exactness

# paired differences of moments before and after one move from the joint law of the model (the ties are the signs of the latent array)
bso_geweke = function(
	mode,
	ndraw,
	seed,
	n = 6L,
	nt = 2L,
	pmean = FALSE,
	withm = TRUE,
	visit = lame:::.lame_bso_visit,
	miss = 0.1,
	rho_in = 0.5,
	centred = FALSE,
	count_centred = centred
) {
	ns = asNamespace('lame')
	set.seed(seed)
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym'); dirn = identical(mode, 'dir')
	nr = n; nc = if (bip) max(3L, n - 2L) else n
	p = 2L; vb = 0.8; va = 0.5; vbb = 0.6
	X = array(rnorm(nr * nc * p * nt), c(nr, nc, p, nt))
	rho_t = if (dirn) rho_in else 0
	mfix = if (withm) array(rnorm(nr * nc * nt, 1, 0.7), c(nr, nc, nt)) else NULL
	mu_b = if (pmean) c(0.5, -0.4) else c(0, 0)
	E_of = function(s) {
		E = array(0, c(nr, nc, nt))
		for (t in 1:nt) E[, , t] = outer(s$a, s$b, '+') + matrix(
			matrix(X[, , , t], nr * nc, p) %*% s$beta,
			nr,
			nc
		) + if (is.null(mfix)) 0 else mfix[, , t]
		E
	}
	Ebar_of = function(E) { if (sym) for (t in seq_len(dim(E)[3L])) E[, , t] = (E[, , t] + t(E[, , t])) / 2; E }
	draw_state = function() {
		beta = mu_b + rnorm(p, 0, sqrt(vb)); a = rnorm(nr, 0, sqrt(va)); b = if (sym) a else rnorm(nc, 0, sqrt(vbb))
		# centred: the additive effects sum to zero over the actors, as with an intercept in lame()
		if (centred) { a = a - mean(a); b = if (sym) a else b - mean(b) }
		E = E_of(list(a = a, b = b, beta = beta))
		Z = Ebar_of(E)
		for (t in 1:nt) {
			z = matrix(rnorm(nr * nc), nr, nc)
			if (dirn) { z2 = matrix(rnorm(nr * nr), nr, nr); lo = lower.tri(z); z[lo] = (rho_t * t(z) + sqrt(1 - rho_t^2) * z2)[lo]; diag(z) = rnorm(nr, 0, sqrt(1 + rho_t)) }
			if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
			Z[, , t] = Z[, , t] + z
		}
		Y = 1 * (Z > 0)
		ms = array(runif(nr * nc * nt) < miss, c(nr, nc, nt))
		if (sym) for (t in 1:nt) ms[, , t][lower.tri(ms[, , t])] = t(ms[, , t])[lower.tri(ms[, , t])]
		Y[ms] = NA
		if (!bip) for (t in 1:nt) diag(Y[, , t]) = NA
		list(beta = beta, a = a, b = b, E = E, Z = Z, Y = Y)
	}
	upm = upper.tri(matrix(0, nr, nc))
	stat_of = function(s) {
		z = s$Z; eb = Ebar_of(s$E); y1 = ifelse(is.na(s$Y), 0, s$Y); ob = ifelse(is.na(s$Y), 0, 1)
		if (sym) { msk = array(rep(upm, nt), dim(z)); z = z * msk; eb = eb * msk; y1 = y1 * msk; ob = ob * msk }
		c(
			s$beta,
			s$beta^2,
			s$beta[1] * s$beta[2],
			sum(s$a^2),
			sum(s$b^2),
			sum(s$a),
			sum(s$b),
			sapply(
				1:nt,
				function(t) c(
					mean(z[, , t]),
					mean(z[, , t]^2),
					mean((z[, , t] - eb[, , t])^2),
					mean(z[, , t] * eb[, , t]),
					mean(y1[, , t]),
					mean(y1[, , t] * eb[, , t]),
					mean(y1[, , t] * z[, , t])
				)
			),
			sum(y1 * eb) * s$beta[1],
			sum(ob * eb) * s$beta[2]
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s = draw_state()
		cx = ns$.lame_bso_context(s$Y, mode)
		Apri = sum(s$beta^2) / vb + sum(s$a^2) / va + (if (sym) 0 else sum(s$b^2) / vbb)
		Bpri = if (pmean) sum(mu_b * s$beta) / vb else 0
		dth = p + nr + (if (sym) 0 else nc) - (if (count_centred) (if (sym) 1 else 2) else 0)
		s1 = s
		o = visit(s$Z, s$E, mfix, cx, rho_t, Apri, Bpri, dth)
		s1$Z = o$Z
		s1$beta = o$c * s$beta; s1$a = o$c * s$a; s1$b = if (sym) o$c * s$a else o$c * s$b
		s1$E = E_of(s1)
		# the ties are held: the refreshed latent array agrees with them
		ok = is.na(s$Y) | (1 * (o$Z > 0) == s$Y)
		if (!all(ok)) stop('the refreshed latent array contradicts a tie')
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

test_that('one step of the scale of the linear predictor leaves the joint law of coefficients, effects, latent array and ties unchanged, in every mode', {
	skip_on_cran()
	for (cfg in list(
		list('dir', 51, FALSE, TRUE),
		list('sym', 52, FALSE, TRUE),
		list('bip', 53, FALSE, FALSE),
		list('dir', 54, TRUE, FALSE),
		list('bip', 55, TRUE, TRUE),
		list('sym', 56, TRUE, FALSE)
	)) {
		r = bso_geweke(cfg[[1]], 7000, cfg[[2]], pmean = cfg[[3]], withm = cfg[[4]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('with centred additive effects (an intercept) the scaled dimension is one less per margin, and counting them uncentred fails the check', {
	skip_on_cran()
	for (cfg in list(list('dir', 71), list('bip', 72), list('sym', 73))) {
		r = bso_geweke(cfg[[1]], 7000, cfg[[2]], centred = TRUE)
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	expect_gt(bso_geweke('dir', 7000, 74, centred = TRUE, count_centred = FALSE)$maxz, 6)
})

test_that('the exactness check fails for a move that leaves out the prior, the Jacobian, the dyad factor, the cell scale, the unscaled part or the linear term', {
	skip_on_cran()
	bad = bso_mut('.lame_bso_visit', 'A_pri * k^2/2', '0 * k^2/2')
	expect_gt(bso_geweke('dir', 6000, 61, visit = bad)$maxz, 6)
	bad = bso_mut('.lame_bso_visit', 'd_th * u', '(d_th + 1) * u')
	expect_gt(bso_geweke('dir', 6000, 62, visit = bad)$maxz, 6)
	bad = bso_mut('.lame_bso_visit', '.lame_dyad_factor_draw(Z, E, rho, cx$upper)', 'NULL')
	expect_gt(bso_geweke('dir', 6000, 63, visit = bad)$maxz, 6)
	bad = bso_mut('.lame_bso_visit', 'inv_sd <- 1/sdc', 'inv_sd <- 1')
	expect_gt(bso_geweke('dir', 6000, 64, visit = bad)$maxz, 6)
	bad = bso_mut('.lame_bso_visit', 'es <- e_all - m_all', 'es <- e_all')
	expect_gt(bso_geweke('bip', 6000, 65, visit = bad)$maxz, 6)
	bad = bso_mut('.lame_bso_visit', '+ B_pri * k', '+ 0 * k')
	expect_gt(bso_geweke('dir', 9000, 66, pmean = TRUE, visit = bad)$maxz, 4.5)
})

# the switch and the fits

bso_panel = function(mode = 'dir', n = 8L, nt = 3L, seed = 3L, sep = FALSE) {
	set.seed(seed)
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	nr = n; nc = if (bip) 6L else n
	nmr = sprintf('a%02d', seq_len(nr)); nmc = if (bip) sprintf('b%02d', seq_len(nc)) else nmr
	X = Y = vector('list', nt)
	for (t in seq_len(nt)) {
		x1 = matrix(rnorm(nr * nc), nr, nc); x2 = matrix(rnorm(nr * nc), nr, nc)
		if (sym) { x1 = (x1 + t(x1)) / 2 * sqrt(2); x2 = (x2 + t(x2)) / 2 * sqrt(2) }
		y = if (sep) 1 * (x1 > 0.5) else 1 * (0.7 * x1 - 0.4 * x2 + outer(
			rnorm(nr, 0, .5),
			rnorm(nc, 0, .5),
			'+'
		) + matrix(rnorm(nr * nc), nr, nc) > 0)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		dimnames(y) = list(nmr, nmc)
		if (!bip) diag(y) = NA
		X[[t]] = array(c(x1, x2), c(nr, nc, 2L), dimnames = list(nmr, nmc, c('x1', 'x2'))); Y[[t]] = y
	}
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	list(Y = Y, X = X)
}

bso_fit = function(
	mode = 'dir',
	nt = 3L,
	...,
	prior = list(),
	nscan = 40,
	burn = 10,
	sep = FALSE,
	n = 8L,
	dcor = identical(mode, 'dir')
) {
	dat = bso_panel(mode, nt = nt, sep = sep, n = n)
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	if (nt == 1L) {
		suppressWarnings(suppressMessages(do.call(
			ame,
			c(
				list(
					dat$Y[[1]],
					Xdyad = dat$X[[1]],
					family = 'binary',
					symmetric = sym,
					mode = if (bip) 'bipartite' else 'unipartite',
					dcor = dcor,
					burn = burn,
					nscan = nscan,
					odens = 1,
					seed = 5,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = prior
				),
				list(...)
			)
		)))
	} else {
		suppressWarnings(suppressMessages(do.call(
			lame,
			c(
				list(
					dat$Y,
					Xdyad = dat$X,
					family = 'binary',
					symmetric = sym,
					mode = if (bip) 'bipartite' else 'unipartite',
					dcor = dcor,
					burn = burn,
					nscan = nscan,
					odens = 1,
					seed = 5,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = prior
				),
				list(...)
			)
		)))
	}
}

test_that('binary_scale_move and binary_scale_every are validated', {
	expect_error(
		bso_fit('dir', 1L, prior = list(binary_scale_move = 1)),
		'binary_scale_move must be TRUE or FALSE',
		fixed = TRUE
	)
	expect_error(
		bso_fit('dir', 3L, prior = list(binary_scale_move = NA)),
		'binary_scale_move must be TRUE or FALSE',
		fixed = TRUE
	)
	expect_error(
		bso_fit('dir', 1L, prior = list(binary_scale_every = 0)),
		'binary_scale_every must be a positive whole number',
		fixed = TRUE
	)
	expect_error(
		bso_fit('bip', 3L, prior = list(binary_scale_every = 2.5)),
		'binary_scale_every must be a positive whole number',
		fixed = TRUE
	)
})

test_that('the move runs by default in every mode of ame() and lame(), is recorded, and the switch omits it', {
	skip_on_cran()
	for (cfg in list(
		list('dir', 1L),
		list('sym', 1L),
		list('bip', 1L),
		list('dir', 3L),
		list('sym', 3L),
		list('bip', 3L)
	)) {
		fit = bso_fit(cfg[[1]], cfg[[2]])
		expect_true(!is.null(fit$binary_scale_move))
		expect_equal(unname(fit$binary_scale_move[['visits']]), 10)
		expect_gt(fit$binary_scale_move[['moved']], 0)
		expect_true(all(is.finite(fit$BETA)) && all(is.finite(fit$VC)))
		off = bso_fit(cfg[[1]], cfg[[2]], prior = list(binary_scale_move = FALSE))
		expect_null(off$binary_scale_move)
		every = bso_fit(cfg[[1]], cfg[[2]], prior = list(binary_scale_every = 20))
		expect_equal(unname(every$binary_scale_move[['visits']]), 2)
	}
	# other families are untouched
	dat = simulate_test_network(n = 8, n_time = 3, family = 'normal', seed = 3)
	fn = suppressWarnings(
		suppressMessages(
			lame(
				dat$Y,
				Xdyad = dat$Xdyad,
				family = 'normal',
				nscan = 10,
				burn = 5,
				odens = 1,
				verbose = FALSE,
				plot = FALSE,
				gof = FALSE,
				seed = 5
			)
		)
	)
	expect_null(fn$binary_scale_move)
})

test_that('the move runs with changing coefficients, changing effects, positions, no intercept, actor-specific slopes, absent actors and the joint Hamiltonian configuration', {
	skip_on_cran()
	for (cfg in list(
		list(mode = 'dir', nt = 4L, dynamic_beta = 'x1_dyad', dynamic_beta_kind = 'ar1'),
		list(mode = 'bip', nt = 4L, dynamic_beta = 'x1_dyad', dynamic_beta_kind = 'rw1'),
		list(mode = 'sym', nt = 3L, dynamic_beta = 'x1_dyad'),
		list(mode = 'dir', nt = 3L, dynamic_ab = TRUE),
		list(mode = 'dir', nt = 3L, R = 1),
		list(mode = 'bip', nt = 3L, R = 1),
		list(mode = 'dir', nt = 3L, intercept = FALSE),
		list(mode = 'sym', nt = 3L, dynamic_ab = TRUE),
		list(mode = 'dir', nt = 1L, R = 1),
		list(mode = 'dir', nt = 1L, intercept = FALSE),
		list(mode = 'dir', nt = 4L, dynamic_beta = 'x1_dyad', dynamic_beta_per_actor = 'row'),
		list(mode = 'dir', nt = 4L, dynamic_beta = 'x1_dyad', dynamic_ab = TRUE, R = 1, dcor = FALSE)
	)) {
		args = cfg; args$mode = NULL; args$nt = NULL
		fit = do.call(bso_fit, c(list(cfg$mode, cfg$nt), args))
		expect_gt(fit$binary_scale_move[['moved']], 0)
		expect_true(all(is.finite(fit$BETA)))
	}
})

test_that('on a panel that one covariate separates, the move gives the coefficient far more effective draws than the sweep alone, and a bipartite coefficient is no longer stuck below its posterior', {
	skip_on_cran()
	run = function(mode, arm, n) {
		# the arms are compared at one run length: an arm still moving at the end (the stuck one is) is
		# refused at the first flag rather than continued, and its fit is taken from the condition
		prior = c(if (identical(arm, 'off')) list(binary_scale_move = FALSE) else list(), list(chain_extensions = 0))
		fit = tryCatch(
			bso_fit(mode, 1L, sep = TRUE, n = n, nscan = 1200, burn = 200, prior = prior),
			lame_chain_unsettled = function(e) e$fit
		)
		fit$BETA[, 'x1_dyad']
	}
	on = run('dir', 'on', 20L); off = run('dir', 'off', 20L)
	expect_gt(coda::effectiveSize(on), 4 * coda::effectiveSize(off))
	expect_gt(coda::effectiveSize(on), 30)
	on = run('bip', 'on', 20L); off = run('bip', 'off', 20L)
	expect_gt(mean(on), 1.5 * mean(off))
})

test_that('a covariate that separates the ties is reported as a property, with the sd of its prior, and a covariate that does not is not', {
	dat = bso_panel('dir', n = 10L, nt = 1L, sep = TRUE)
	run = function(Y, X, ...) ame(
		Y,
		Xdyad = X,
		family = 'binary',
		dcor = FALSE,
		burn = 2,
		nscan = 6,
		odens = 1,
		seed = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		...
	)
	w = testthat::capture_warnings(suppressMessages(run(dat$Y[[1]], dat$X[[1]])))
	msg = paste(w, collapse = '\n')
	expect_match(msg, 'x1', fixed = TRUE)
	expect_match(msg, 'separates the ties from the non-ties', fixed = TRUE)
	expect_match(msg, 'every tie lies above every non-tie', fixed = TRUE)
	expect_match(msg, 'a normal with sd', fixed = TRUE)
	expect_false(grepl('drifts|depends on chain length', msg))
	# the sd it states is the sd of the prior the sampler draws the coefficient under
	Y = dat$Y[[1]]; X = dat$X[[1]]
	n = nrow(Y)
	Xd = lame:::precomputeX(lame:::design_array(NULL, NULL, X, TRUE, n, warn = FALSE))
	V0 = solve(lame:::.beta_prior_precision(attributes(Xd)$XX, attributes(Xd)$mX, sum(!is.na(Y))))
	sd_x1 = sqrt(V0[2, 2])
	expect_equal(lame:::.ame_separation_prior_sd(Y, X, NULL, NULL, TRUE, 'unipartite', NA, list(), 1L), sd_x1)
	expect_match(msg, format(signif(sd_x1, 2)), fixed = TRUE)
	# a bipartite coefficient has the independent prior N(0, g s20) with g the number of cells
	expect_equal(
		lame:::.ame_separation_prior_sd(
			matrix(0, 6, 4),
			array(0, c(6, 4, 1)),
			NULL,
			NULL,
			TRUE,
			'bipartite',
			NA,
			list(s20 = 2),
			1L
		),
		sqrt(24 * 2)
	)
	# a covariate that leaves the ties unseparated is not reported
	ns = bso_panel('dir', n = 10L, nt = 1L, sep = FALSE)
	expect_false(any(grepl('separates', testthat::capture_warnings(suppressMessages(run(ns$Y[[1]], ns$X[[1]]))))))
})
