# moves along the scale of the latent array of a rank-likelihood ordinal fit (R/ordinal_scale_move.R): the metric of
# the dyadic covariance, the density of the scale draw, exactness of the two moves under the joint law without data
# (with a broken version failing), the one-pass visit against its definition, and the switch and fit records
skip_on_cran()

osm_mut = function(name, from, to) {
	src = deparse(get(name, asNamespace('lame')), width.cutoff = 500L)
	out = sub(from, to, src, fixed = TRUE)
	if (identical(out, src)) stop('mutation target not found: ', from)
	fn = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(fn) = asNamespace('lame')
	fn
}

# the metric

test_that('the inner product is the metric of the dyadic covariance: dyads with correlation, a diagonal at 1 + rho, pairs, bipartite cells', {
	ns = asNamespace('lame')
	set.seed(1)
	n = 5L
	x = matrix(rnorm(n * n), n); y = matrix(rnorm(n * n), n)
	rho = 0.45; s2 = 1.3
	# the same quadratic form from the covariance of each dyad
	S = s2 * matrix(c(1, rho, rho, 1), 2)
	dot_ref = function(x, y) {
		s = 0
		for (i in 1:(n - 1)) for (j in (i + 1):n) s = s + drop(c(x[i, j], x[j, i]) %*% solve(S) %*% c(y[i, j], y[j, i]))
		s + sum(diag(x) * diag(y)) / ((1 + rho) * s2)
	}
	expect_equal(ns$.lame_osm_dot(x, y, rho, s2, 'dir'), dot_ref(x, y), tolerance = 1e-12)
	expect_equal(ns$.lame_osm_dot(x, x, rho, s2, 'dir'), dot_ref(x, x), tolerance = 1e-12)
	up = upper.tri(x)
	# a symmetric network has one latent value per pair, around the mean of the two directed predictors
	xs = (x + t(x)) / 2; ys = (y + t(y)) / 2
	expect_equal(ns$.lame_osm_dot(x, y, 0, s2, 'sym'), sum(xs[up] * ys[up]) / s2, tolerance = 1e-12)
	expect_equal(ns$.lame_osm_dot(xs, ys, 0, s2, 'sym'), sum(xs[up] * ys[up]) / s2, tolerance = 1e-12)
	xb = matrix(rnorm(12), 3); yb = matrix(rnorm(12), 3)
	expect_equal(ns$.lame_osm_dot(xb, yb, 0, s2, 'bip'), sum(xb * yb) / s2, tolerance = 1e-12)
	expect_equal(ns$.lame_osm_ncell(5, 5, 'dir'), 25)
	expect_equal(ns$.lame_osm_ncell(5, 5, 'sym'), 10)
	expect_equal(ns$.lame_osm_ncell(3, 4, 'bip'), 12)
})

# the scale draw

test_that('the scale draw has the density c^(D - 1) exp(-A c^2 / 2 + B c): gamma when B is zero, the Laplace independence step otherwise', {
	skip_on_cran()
	ns = asNamespace('lame')
	set.seed(2)
	D = 30; A = 25
	for (B in c(0, 12)) {
		cc = replicate(40000, ns$.lame_osm_draw_c(D, A, B))
		grid = seq(0.01, 4, length.out = 20000)
		dens = exp(
			(D - 1) * log(grid) - A * grid^2 / 2 + B * grid - max((D - 1) * log(grid) - A * grid^2 / 2 + B * grid)
		)
		cdf = cumsum(dens) / sum(dens)
		# a Metropolis step from one fixed start is not an independent draw, so only the exact draw is compared with the
		# density by quantiles; the Laplace step is covered by the exactness test below
		if (B == 0) {
			q = c(0.1, 0.25, 0.5, 0.75, 0.9)
			expect_equal(unname(quantile(cc, q)), approx(cdf, grid, xout = q, ties = min)$y, tolerance = 0.01)
		} else {
			expect_true(all(is.finite(cc)) && all(cc > 0))
		}
	}
	expect_equal(ns$.lame_osm_draw_c(1, 5, 0), 1)
	expect_equal(ns$.lame_osm_draw_c(10, 0, 0), 1)
	expect_equal(ns$.lame_osm_draw_c(10, NA, 0), 1)
})

# exactness

# paired differences of moments before and after one move from the joint law of the model without data
osm_geweke = function(
	mode,
	part,
	ndraw,
	seed,
	n = 6L,
	nt = 3L,
	pmean = FALSE,
	withm = TRUE,
	visit = lame:::.lame_osm_visit,
	centred = FALSE,
	count_centred = centred
) {
	set.seed(seed)
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	nr = n; nc = if (bip) max(3L, n - 2L) else n
	p = 2L; vb = 0.8; va = 0.5; vbb = 0.6
	X = array(rnorm(nr * nc * p * nt), c(nr, nc, p, nt))
	# a symmetric network keeps differing directed designs, as a row or column covariate does; a pair has their mean
	rho_t = if (identical(mode, 'dir')) runif(nt, -0.5, 0.7) else rep(0, nt)
	mfix = if (withm) array(rnorm(nr * nc * nt, 0, 0.7), c(nr, nc, nt)) else NULL
	mu_b = if (pmean) c(0.5, -0.4) else c(0, 0)
	upm = upper.tri(matrix(0, nr, nc))
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
		# centred: the additive effects sum to zero over actors, as with an intercept in lame()
		if (centred) { a = a - mean(a); b = if (sym) a else b - mean(b) }
		E = E_of(list(a = a, b = b, beta = beta))
		Z = Ebar_of(E)
		for (t in 1:nt) {
			r = rho_t[t]; z = matrix(rnorm(nr * nc), nr, nc)
			if (identical(mode, 'dir')) { z2 = matrix(rnorm(nr * nr), nr, nr); lo = lower.tri(z); z[lo] = (r * t(z) + sqrt(1 - r^2) * z2)[lo]; diag(z) = rnorm(nr, 0, sqrt(1 + r)) }
			if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
			Z[, , t] = Z[, , t] + z
		}
		list(beta = beta, a = a, b = b, E = E, Z = Z)
	}
	stat_of = function(s) {
		z = s$Z; eb = Ebar_of(s$E); e = z - eb
		if (sym) { z = array(
			vapply(1:nt, function(t) { m_ = z[, , t]; m_[!upm] = 0; m_ }, matrix(0, nr, nc)),
			dim(z)
		); e = array(vapply(1:nt, function(t) { m_ = e[, , t]; m_[!upm] = 0; m_ }, matrix(0, nr, nc)), dim(e)) }
		c(
			s$beta,
			s$beta^2,
			s$beta[1] * s$beta[2],
			sum(s$a^2),
			sum(s$b^2),
			sum(s$a),
			sum(s$b),
			sapply(1:nt, function(t) c(mean(z[, , t]), mean(z[, , t]^2), mean(e[, , t]^2), mean(z[, , t] * eb[, , t]))),
			mean(z^2) * s$beta[1]
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s = draw_state()
		Apri = sum(s$beta^2) / vb + sum(s$a^2) / va + (if (sym) 0 else sum(s$b^2) / vbb)
		Bpri = if (pmean) sum(mu_b * s$beta) / vb else 0
		dth = p + nr + (if (sym) 0 else nc) - (if (count_centred) (if (sym) 1 else 2) else 0)
		s1 = s
		# orbit: the scale of the state; affine: the shift and scale of each period; both: one visit with both
		o = visit(
			s$Z,
			s$E,
			mfix,
			mode,
			rho_t,
			1,
			Apri,
			Bpri,
			dth,
			affine = part %in% c('affine', 'both'),
			orbit = part %in% c('orbit', 'both')
		)
		s1$Z = o$Z
		if (part %in% c('orbit', 'both')) { s1$beta = o$c * s$beta; s1$a = o$c * s$a; s1$b = if (sym) o$c * s$a else o$c * s$b }
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

test_that('one step of the scale of the whole state and of each period leaves the joint law of coefficients, effects and latent array unchanged, in every mode', {
	skip_on_cran()
	for (cfg in list(
		list('dir', 'both', 31, FALSE, TRUE),
		list('sym', 'both', 32, FALSE, TRUE),
		list('bip', 'both', 33, FALSE, FALSE),
		list('dir', 'orbit', 34, TRUE, TRUE),
		list('bip', 'affine', 35, FALSE, TRUE),
		list('sym', 'orbit', 36, TRUE, FALSE)
	)) {
		r = osm_geweke(cfg[[1]], cfg[[2]], 8000, cfg[[3]], pmean = cfg[[4]], withm = cfg[[5]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('the visit written out in terms of the latent array is exact too, and the one-pass visit reproduces it with the same random numbers', {
	ns = asNamespace('lame')
	skip_on_cran()
	r = osm_geweke('dir', 'both', 6000, 37, visit = ns$.lame_osm_visit_ref)
	expect_lt(r$maxz, 4.5)
	expect_gt(r$p, 1e-3)
})

test_that('the one-pass visit equals the definition (the same random numbers give the same latent array and scale) in every mode, with and without an unscaled part', {
	ns = asNamespace('lame')
	set.seed(7)
	for (mode in c('dir', 'sym', 'bip')) for (withm in c(FALSE, TRUE)) {
		nr = 7L; nc = if (mode == 'bip') 5L else 7L; nt = 3L
		# a symmetric network has a symmetric latent array and differing directed predictors; a pair has their mean
		sq = function(a) { if (mode == 'sym') for (t in seq_len(dim(a)[3L])) a[, , t] = (a[, , t] + t(a[, , t])) / 2; a }
		Z = sq(
			array(rnorm(nr * nc * nt, 0, 1.3), c(nr, nc, nt))
		); m = if (withm) array(rnorm(nr * nc * nt, 0, 0.6), c(nr, nc, nt)) else NULL
		E = array(rnorm(nr * nc * nt, 0, 0.8), c(nr, nc, nt)) + if (withm) m else 0
		rho = if (mode == 'dir') c(0.5, -0.3, 0.1) else rep(0, nt)
		for (orbit in c(TRUE, FALSE)) for (affine in c(TRUE, FALSE)) for (B_pri in c(0, 1.5)) {
			set.seed(11); ref = ns$.lame_osm_visit_ref(
				Z,
				E,
				m,
				mode,
				rho,
				1.4,
				9,
				B_pri,
				12,
				affine = affine,
				orbit = orbit
			)
			set.seed(11); one = ns$.lame_osm_visit(
				Z,
				E,
				m,
				mode,
				rho,
				1.4,
				9,
				B_pri,
				12,
				affine = affine,
				orbit = orbit
			)
			expect_equal(one$c, ref$c, tolerance = 1e-10)
			expect_equal(one$Z, ref$Z, tolerance = 1e-10)
		}
		# the sums the kernel returns are the inner products of the metric
		G = ns$osm_gram_cpp(
			Z,
			E,
			if (withm) m else numeric(0),
			c(nr, nc, nt),
			switch(mode, dir = 1L, sym = 2L, bip = 3L)
		)
		kf = if (mode == 'dir') cbind(
			1 / (2 * (1 + rho) * 1.4),
			1 / (2 * (1 - rho) * 1.4),
			1 / ((1 + rho) * 1.4)
		) else cbind(rep(1 / 1.4, nt), 0, 0)
		for (t in seq_len(nt)) {
			expect_equal(
				sum(G[t, 1:3] * kf[t, ]),
				ns$.lame_osm_dot(Z[, , t], Z[, , t], rho[t], 1.4, mode),
				tolerance = 1e-10
			)
			expect_equal(
				sum(G[t, 4:6] * kf[t, ]),
				ns$.lame_osm_dot(Z[, , t], E[, , t], rho[t], 1.4, mode),
				tolerance = 1e-10
			)
			expect_equal(
				sum(G[t, 10:12] * kf[t, ]),
				ns$.lame_osm_dot(E[, , t], E[, , t], rho[t], 1.4, mode),
				tolerance = 1e-10
			)
		}
	}
	expect_error(lame:::osm_gram_cpp(array(0, c(3, 3, 2)), array(0, c(3, 3, 2)), numeric(0), c(3L, 4L, 2L), 1L))
})

test_that('with centred additive effects (an intercept) the scaled dimension is one less per margin, and counting them uncentred fails the check', {
	skip_on_cran()
	for (cfg in list(list('dir', 'both', 81), list('bip', 'both', 82), list('sym', 'orbit', 83))) {
		r = osm_geweke(cfg[[1]], cfg[[2]], 7000, cfg[[3]], centred = TRUE)
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	expect_lt(osm_geweke('dir', 'orbit', 9000, 84, centred = TRUE, count_centred = FALSE)$p, 5e-3)
})

test_that('the exactness check fails for a scale move that leaves out the priors, the Jacobian, the unscaled part or the right metric', {
	skip_on_cran()
	bad_prior = osm_mut('.lame_osm_visit', 'A <- A_pri', 'A <- 0')
	expect_gt(osm_geweke('dir', 'orbit', 6000, 41, visit = bad_prior)$maxz, 6)
	bad_d = osm_mut('.lame_osm_visit', 'nt * ncell + d_th', 'nt * ncell')
	expect_gt(osm_geweke('dir', 'orbit', 6000, 42, visit = bad_d)$maxz, 6)
	bad_m = osm_mut('.lame_osm_visit', 'B <- B_pri + sum(ZM - EM + MM)', 'B <- B_pri')
	expect_gt(osm_geweke('dir', 'orbit', 6000, 43, visit = bad_m)$maxz, 6)
	# mutant: a shift draw that ignores the mean residual
	bad_shift = osm_mut('.lame_osm_visit', 'mu <- -(cc * (sZ[tt] - sE[tt]) - (1 - cc) * sM[tt])/ncell', 'mu <- 0')
	expect_gt(osm_geweke('dir', 'affine', 6000, 44, visit = bad_shift)$maxz, 6)
	# mutant: a period scale that drops the shift from its quadratic term
	bad_a2 = osm_mut('.lame_osm_visit', '2 * cc * d * sZ[tt] * ko[tt]', '0')
	expect_gt(osm_geweke('dir', 'affine', 6000, 46, visit = bad_a2)$maxz, 6)
	# mutant: a period scale whose linear term uses the state scale once rather than twice
	bad_b2 = osm_mut('.lame_osm_visit', 'B2 <- cc^2 * ZE[tt]', 'B2 <- cc * ZE[tt]')
	expect_gt(osm_geweke('dir', 'both', 6000, 47, visit = bad_b2)$maxz, 6)
	# mutant: a metric that confuses the sum and difference of a dyad
	bad_kf = osm_mut(
		'.lame_osm_visit', '1/(2 * (1 + rho) * s2), 1/(2 * (1 - rho) * s2)', '1/(2 * (1 - rho) * s2), 1/(2 * (1 + rho) * s2)'
	)
	expect_gt(osm_geweke('dir', 'orbit', 6000, 45, visit = bad_kf)$maxz, 6)
})

# the switch and the fits

osm_panel = function(mode = 'dir', n = 8L, nt = 3L, seed = 3L, na_actors = 0) {
	set.seed(seed)
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	nr = n; nc = if (bip) 6L else n
	nmr = sprintf('a%02d', seq_len(nr)); nmc = if (bip) sprintf('b%02d', seq_len(nc)) else nmr
	X = Y = vector('list', nt)
	for (t in seq_len(nt)) {
		x1 = matrix(rnorm(nr * nc), nr, nc)
		if (sym) x1 = (x1 + t(x1)) / 2
		z = 0.7 * x1 + outer(rnorm(nr, 0, .5), rnorm(nc, 0, .5), '+') + matrix(rnorm(nr * nc), nr, nc)
		if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
		y = matrix(findInterval(z, quantile(z, c(0.5, 0.8))) + 1L, nr, nc, dimnames = list(nmr, nmc))
		if (!bip) diag(y) = NA
		X[[t]] = array(x1, c(nr, nc, 1L), dimnames = list(nmr, nmc, 'x1')); Y[[t]] = y
	}
	if (na_actors > 0) { gone = seq_len(na_actors); Y[[nt]][gone, ] = NA; if (!bip) Y[[nt]][, gone] = NA }
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	list(Y = Y, X = X)
}

osm_fit = function(mode = 'dir', nt = 3L, ..., prior = list(), nscan = 30, burn = 10) {
	dat = osm_panel(mode, nt = nt, na_actors = list(...)$na_actors %||% 0)
	extra = list(...); extra$na_actors = NULL
	bip = identical(mode, 'bip'); sym = identical(mode, 'sym')
	if (nt == 1L) {
		suppressWarnings(suppressMessages(do.call(
			ame,
			c(
				list(
					dat$Y[[1]],
					Xdyad = dat$X[[1]],
					family = 'ordinal',
					symmetric = sym,
					mode = if (bip) 'bipartite' else 'unipartite',
					dcor = identical(mode, 'dir'),
					burn = burn,
					nscan = nscan,
					odens = 1,
					seed = 5,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = prior
				),
				extra
			)
		)))
	} else {
		suppressWarnings(suppressMessages(do.call(
			lame,
			c(
				list(
					dat$Y,
					Xdyad = dat$X,
					family = 'ordinal',
					symmetric = sym,
					mode = if (bip) 'bipartite' else 'unipartite',
					dcor = identical(mode, 'dir'),
					burn = burn,
					nscan = nscan,
					odens = 1,
					seed = 5,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = prior
				),
				extra
			)
		)))
	}
}

test_that('ordinal_scale_move is validated', {
	expect_error(
		osm_fit('dir', 1L, prior = list(ordinal_scale_move = 1)),
		'ordinal_scale_move must be TRUE or FALSE',
		fixed = TRUE
	)
	expect_error(
		osm_fit('dir', 3L, prior = list(ordinal_scale_move = NA)),
		'ordinal_scale_move must be TRUE or FALSE',
		fixed = TRUE
	)
})

test_that('the moves run by default in every mode of ame() and lame(), are recorded, and the switch omits them', {
	skip_on_cran()
	for (cfg in list(
		list('dir', 1L),
		list('sym', 1L),
		list('bip', 1L),
		list('dir', 3L),
		list('sym', 3L),
		list('bip', 3L)
	)) {
		fit = osm_fit(cfg[[1]], cfg[[2]])
		expect_true(!is.null(fit$ordinal_scale_move))
		expect_equal(unname(fit$ordinal_scale_move[['visits']]), 40)
		expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
		expect_true(all(is.finite(fit$BETA)) && all(is.finite(fit$VC)))
		off = osm_fit(cfg[[1]], cfg[[2]], prior = list(ordinal_scale_move = FALSE))
		expect_null(off$ordinal_scale_move)
		expect_false(isTRUE(all.equal(fit$BETA, off$BETA)))
	}
	# other families are untouched
	dat = simulate_test_network(n = 8, n_time = 3, family = 'binary', seed = 3)
	fb = suppressWarnings(
		suppressMessages(
			lame(
				dat$Y,
				Xdyad = dat$Xdyad,
				family = 'binary',
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
	expect_null(fb$ordinal_scale_move)
})

test_that('the moves run with a dynamic coefficient, dynamic additive effects, positions, explicit cutpoints off, and an actor absent from a period', {
	skip_on_cran()
	fit = osm_fit('dir', 4L, dynamic_beta = 'x1_dyad', dynamic_beta_kind = 'ar1')
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
	expect_true(all(is.finite(fit$BETA)))
	fit = osm_fit('bip', 4L, dynamic_beta = 'x1_dyad', dynamic_beta_kind = 'rw1')
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
	fit = osm_fit('dir', 3L, dynamic_ab = TRUE)
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
	fit = osm_fit('dir', 3L, R = 1)
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
	fit = osm_fit('dir', 3L, na_actors = 2)
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
	expect_true(all(is.finite(fit$BETA)) && all(is.finite(fit$VC)))
	fit = osm_fit('dir', 1L, na_actors = 2)
	expect_gt(fit$ordinal_scale_move[['orbit_moved']], 0)
})
