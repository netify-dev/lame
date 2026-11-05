# ALS start values passed to the MCMC samplers: an ALS fit labels its effects and positions by actor (and period),
# while the sampler state is plain numeric, so the labels stay out of the translation moves that index it by actor
skip_on_cran()

asm_panel = function(mode, family, n = 8L, nb = 6L, N = 3L, seed = 11L) {
	set.seed(seed)
	bip = identical(mode, 'bipartite'); sym = identical(mode, 'symmetric')
	nc = if (bip) nb else n
	rn = sprintf('a%02d', seq_len(n)); cn = if (bip) sprintf('b%02d', seq_len(nc)) else rn
	a0 = rnorm(n, 0, .5); b0 = if (sym) a0 else rnorm(nc, 0, .5)
	Y = X = vector('list', N)
	for (t in seq_len(N)) {
		x = matrix(rnorm(n * nc), n, nc)
		if (sym) x[lower.tri(x)] = t(x)[lower.tri(x)]
		z = 0.2 + 0.4 * x + outer(a0, b0, '+') + matrix(rnorm(n * nc, 0, .5), n, nc)
		if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
		y = switch(family, normal = z, binary = 1 * (z > 0.3), poisson = matrix(rpois(n * nc, exp(z)), n, nc))
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, nc, 1L), dimnames = list(rn, cn, 'x'))
	}
	list(Y = Y, X = X, mode = mode, family = family)
}

# the ALS fit of the panel and the MCMC fit that starts from it (every move visits at every sweep, so that the
# translations see the start values as they are)
asm_args = function(p, dab, duv, R = 1L) {
	a = list(p$Y, Xdyad = p$X, family = p$family, R = R, dynamic_ab = dab, dynamic_uv = duv)
	if (identical(p$mode, 'bipartite')) a$mode = 'bipartite'
	if (identical(p$mode, 'symmetric')) a$symmetric = TRUE
	a
}
asm_als = function(p, dab, duv, bias_correct = NULL) {
	a = c(asm_args(p, dab, duv), list(method = 'als', verbose = FALSE, gof = FALSE, plot = FALSE))
	if (!is.null(bias_correct)) a$bias_correct = bias_correct
	suppressWarnings(suppressMessages(do.call(lame, a)))
}
asm_mcmc = function(p, dab, duv, sv, every = 1L, seed = 1L) {
	a = c(
		asm_args(p, dab, duv),
		list(
			burn = 3,
			nscan = 6,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			seed = seed,
			prior = list(actor_move_every = every),
			start_vals = sv
		)
	)
	suppressWarnings(suppressMessages(do.call(lame, a)))
}
asm_plain = function(sv) {
	for (nm in c('beta', 'a', 'b', 'U', 'V', 'G')) if (!is.null(sv[[nm]])) sv[[nm]] = unname(sv[[nm]])
	sv
}

test_that('a bipartite dynamic_ab fit starts from its ALS fit', {
	p = asm_panel('bipartite', 'poisson')
	fit = asm_als(p, TRUE, FALSE, bias_correct = 'none')
	sv = als_start_vals(fit)
	# the start carries actor labels, as the fit it comes from
	expect_equal(dim(sv$a), c(8L, 3L))
	expect_false(is.null(dimnames(sv$a)))
	m = asm_mcmc(p, TRUE, FALSE, sv, every = 2L)
	expect_equal(dim(m$BETA)[1L], 6L)
	expect_true(all(is.finite(m$BETA)))
	expect_true(all(is.finite(m$VC)))
	# the labels of the actors come from the panel
	expect_identical(names(m$APM), rownames(p$Y[[1L]]))
	expect_identical(names(m$BPM), colnames(p$Y[[1L]]))
})

test_that('a start that carries the paths of one margin only runs', {
	p = asm_panel('bipartite', 'poisson')
	sv = als_start_vals(asm_als(p, TRUE, FALSE, bias_correct = 'none'))
	for (drop in c('a', 'b')) {
		s = sv; s[[drop]] = NULL
		m = asm_mcmc(p, TRUE, FALSE, s, every = 2L)
		expect_equal(dim(m$BETA)[1L], 6L, info = drop)
		expect_true(all(is.finite(m$BETA)) && all(is.finite(m$VC)), info = drop)
	}
})

test_that('a start with one margin missing runs in every mode', {
	skip_on_cran()
	for (mode in c('directed', 'symmetric', 'bipartite')) {
		p = asm_panel(mode, 'normal')
		sv = als_start_vals(asm_als(p, TRUE, FALSE, bias_correct = 'none'))
		for (drop in c('a', 'b', 'ab')) {
			s = sv
			for (d in strsplit(drop, '')[[1L]]) s[[d]] = NULL
			m = asm_mcmc(p, TRUE, FALSE, s, every = 2L)
			expect_true(all(is.finite(m$BETA)), info = paste(mode, drop))
		}
	}
})

test_that('the level prior of a translation does not carry the names of the levels', {
	nm = c(a = .3, b = -.1, c = .5, d = -.7)
	iV = matrix(.5, 1, 1)
	pr = .lame_actor_move_level_prior_bip(2L, nm, .8, iV, 0, TRUE)
	expect_identical(names(pr), c('a1', 'a2'))
	expect_equal(unname(pr), unname(.lame_actor_move_level_prior_bip(2L, unname(nm), .8, iV, 0, TRUE)))
	pr = .lame_actor_move_level_prior_bip(2L, nm, .8, iV, 0, FALSE)
	expect_identical(names(pr), c('a1', 'a2'))
	iS = solve(matrix(c(1, .3, .3, 1.2), 2))
	for (row in c(TRUE, FALSE)) {
		pr = .lame_actor_move_level_prior(3L, row, nm, rev(nm), iS, iV, 0, 4L, TRUE)
		expect_identical(names(pr), c('a1', 'a2'))
		expect_equal(
			unname(pr),
			unname(.lame_actor_move_level_prior(3L, row, unname(nm), unname(rev(nm)), iS, iV, 0, 4L, TRUE))
		)
	}
	pr = .lame_actor_move_level_prior_sym(3L, nm, .9, iV, 0, 4L, TRUE)
	expect_identical(names(pr), c('a1', 'a2'))
})

test_that('start values enter the samplers with plain coefficients and effects', {
	sv = list(
		beta = c(intercept = 1, x = .5),
		a = c(r1 = 1, r2 = 2),
		b = matrix(1:4, 2, dimnames = list(c('c1', 'c2'), c('t1', 't2'))),
		U = matrix(0, 2, 1, dimnames = list(c('r1', 'r2'), NULL)),
		G = matrix(1, 1, 1, dimnames = list('u', 'v')),
		s2 = c(s = 1),
		rho_uv = .5
	)
	pl = .start_vals_plain(sv)
	for (nm in c('beta', 'a', 'b')) {
		expect_null(names(pl[[nm]]))
		expect_null(dimnames(pl[[nm]]))
		expect_equal(as.numeric(pl[[nm]]), as.numeric(sv[[nm]]))
	}
	expect_identical(dim(pl$b), c(2L, 2L))
	# the positions and the interaction matrix, the scalars and everything that is not a numeric slot are left as they are
	expect_identical(pl$U, sv$U)
	expect_identical(pl$G, sv$G)
	expect_identical(pl$s2, sv$s2)
	expect_identical(pl$rho_uv, sv$rho_uv)
	expect_null(.start_vals_plain(NULL))
})

test_that('the chain from named start values is the chain from plain ones', {
	skip_on_cran()
	for (mode in c('bipartite', 'directed')) {
		p = asm_panel(mode, 'poisson')
		sv = als_start_vals(asm_als(p, TRUE, FALSE, bias_correct = 'none'))
		m1 = asm_mcmc(p, TRUE, FALSE, sv, every = 1L)
		m2 = asm_mcmc(p, TRUE, FALSE, asm_plain(sv), every = 1L)
		for (nm in c('BETA', 'VC', 'APM', 'BPM')) expect_identical(m1[[nm]], m2[[nm]])
	}
})

test_that('ALS start values run through every mode, family and dynamic effect', {
	skip_on_cran()
	for (mode in c('directed', 'symmetric', 'bipartite')) for (family in c('normal', 'binary', 'poisson')) {
		p = asm_panel(mode, family)
		for (dab in c(FALSE, TRUE)) for (duv in c(FALSE, TRUE)) {
			fit = asm_als(p, dab, duv, bias_correct = 'none')
			sv = als_start_vals(fit)
			label = paste(mode, family, 'ab', dab, 'uv', duv)
			m = tryCatch(asm_mcmc(p, dab, duv, sv), error = function(e) e)
			expect_false(
				inherits(m, 'error'),
				info = if (inherits(m, 'error')) paste(label, conditionMessage(m)) else label
			)
			if (inherits(m, 'error')) next
			expect_true(all(is.finite(m$BETA)), info = label)
			expect_true(all(is.finite(m$VC)), info = label)
			expect_equal(dim(m$BETA)[1L], 6L, info = label)
		}
	}
})

test_that('ALS start values run through the default estimators of the count and binary families', {
	skip_on_cran()
	for (mode in c('directed', 'bipartite')) for (family in c('binary', 'poisson')) {
		p = asm_panel(mode, family)
		for (dab in c(FALSE, TRUE)) {
			sv = als_start_vals(asm_als(p, dab, FALSE))
			label = paste(mode, family, 'ab', dab)
			m = tryCatch(asm_mcmc(p, dab, FALSE, sv), error = function(e) e)
			expect_false(
				inherits(m, 'error'),
				info = if (inherits(m, 'error')) paste(label, conditionMessage(m)) else label
			)
		}
	}
})

test_that('an ALS start also enters a fit with changing coefficients', {
	skip_on_cran()
	for (mode in c('directed', 'bipartite')) {
		p = asm_panel(mode, 'poisson')
		sv = als_start_vals(asm_als(p, FALSE, FALSE, bias_correct = 'none'))
		a = c(
			asm_args(p, FALSE, FALSE),
			list(
				dynamic_beta = TRUE,
				burn = 3,
				nscan = 6,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(actor_move_every = 1),
				start_vals = sv
			)
		)
		m = tryCatch(suppressWarnings(suppressMessages(do.call(lame, a))), error = function(e) e)
		expect_false(inherits(m, 'error'), info = if (inherits(m, 'error')) paste(mode, conditionMessage(m)) else mode)
		if (!inherits(m, 'error')) expect_true(all(is.finite(m$BETA)))
	}
})

test_that('ame() starts a bipartite fit from an ALS fit in every family', {
	skip_on_cran()
	for (family in c('normal', 'binary', 'poisson')) {
		# a panel large enough that the binary ALS fit does not separate
		pf = asm_panel('bipartite', family, n = 16L, nb = 12L)
		y = pf$Y[[1L]]; x = pf$X[[1L]]
		fit = suppressWarnings(suppressMessages(ame(
			y,
			Xdyad = x,
			family = family,
			R = 1,
			mode = 'bipartite',
			method = 'als',
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			bias_correct = 'none'
		)))
		sv = als_start_vals(fit)
		# an ALS start carries no latent array: the sampler draws its own
		expect_null(sv$Z)
		m = suppressWarnings(suppressMessages(ame(
			y,
			Xdyad = x,
			family = family,
			R = 1,
			mode = 'bipartite',
			burn = 3,
			nscan = 6,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			seed = 1,
			prior = list(actor_move_every = 1),
			start_vals = sv
		)))
		expect_true(all(is.finite(m$BETA)), info = family)
	}
})

test_that('a bipartite ame() start whose slots do not fit the panel is reported and the defaults are used', {
	skip_on_cran()
	pf = asm_panel('bipartite', 'poisson')
	y = pf$Y[[1L]]; x = pf$X[[1L]]
	sv = list(beta = c(0, 0, 0, 0), a = rep(0, 3L), U = matrix(0, 8L, 2L))
	expect_warning(
		m <- suppressMessages(ame(
			y,
			Xdyad = x,
			family = 'poisson',
			R = 1,
			mode = 'bipartite',
			burn = 3,
			nscan = 6,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			seed = 1,
			start_vals = sv
		)),
		'do not match this bipartite panel'
	)
	expect_true(all(is.finite(m$BETA)))
	# a start that names only some slots keeps the defaults of the others
	m2 = suppressMessages(ame(
		y,
		Xdyad = x,
		family = 'poisson',
		R = 1,
		mode = 'bipartite',
		burn = 3,
		nscan = 6,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 1,
		start_vals = list(beta = c(0, .3))
	))
	expect_true(all(is.finite(m2$BETA)))
})
