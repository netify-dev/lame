# exact draw of the sender-receiver shift a + c, b - c of fits without an intercept

abs_sim = function(n = 12, TT = 3, seed = 11, nB = NULL) {
	set.seed(seed)
	nb = if (is.null(nB)) n else nB
	a = rnorm(n, 0.5, 1)
	b = rnorm(nb, -0.3, 1)
	lapply(seq_len(TT), function(t) {
		y = outer(a, b, '+') + matrix(rnorm(n * nb), n, nb)
		if (is.null(nB)) diag(y) = NA
		col_names = if (is.null(nB)) paste0('r', seq_len(n)) else paste0('c', seq_len(nb))
		dimnames(y) = list(paste0('r', seq_len(n)), col_names)
		y
	})
}

abs_fit = function(Y, burn = 20, nscan = 40, seed = 3, prior = list(), ...) {
	suppressWarnings(suppressMessages(lame(
		Y,
		family = 'normal',
		R = 0,
		burn = burn,
		nscan = nscan,
		odens = 1,
		seed = seed,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

# z scores of prior-started chains alternating a draw of a given b with the shift
abs_law_z = function(shift, paired, nA, nB, Sab, reps = 4000, steps = 10, seed = 1) {
	set.seed(seed)
	L = t(chol(Sab))
	out = matrix(NA_real_, reps, 4)
	for (r in seq_len(reps)) {
		if (paired) {
			x = L %*% matrix(rnorm(2 * nA), 2)
			a = x[1, ]
			b = x[2, ]
		} else {
			a = rnorm(nA, 0, sqrt(Sab[1, 1]))
			b = rnorm(nB, 0, sqrt(Sab[2, 2]))
		}
		for (k in seq_len(steps)) {
			a = if (paired) {
				Sab[1, 2] / Sab[2, 2] * b + rnorm(nA, 0, sqrt(Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2]))
			} else {
				rnorm(nA, 0, sqrt(Sab[1, 1]))
			}
			cc = shift(a, b, Sab, paired)
			a = a + cc
			b = b - cc
		}
		out[r, ] = c(a[1], b[1], sum(a), sum(b))
	}
	tot = rowSums(out[, 3:4])
	out = cbind(out, tot)
	v = c(
		Sab[1, 1],
		Sab[2, 2],
		nA * Sab[1, 1],
		nB * Sab[2, 2],
		if (paired) nA * sum(Sab) else nA * Sab[1, 1] + nB * Sab[2, 2]
	)
	c(colMeans(out) / sqrt(v / reps), (apply(out, 2, var) / v - 1) / sqrt(2 / (reps - 1)))
}

test_that('the shift draw keeps the prior of the effects, and a wrong-signed mean does not', {
	skip_on_cran()
	draw = lame:::.lame_ab_shift_draw
	broken = function(a, b, Sab, paired) -draw(a, b, Sab, paired)
	S_pair = matrix(c(1, .5, .5, 2), 2)
	S_bip = diag(c(1, 2.5))
	# a shorter law on CRAN
	law = if (identical(Sys.getenv('NOT_CRAN'), 'true')) c(4000, 10) else c(800, 3)
	z = function(shift, paired, nA, nB, Sab) abs_law_z(shift, paired, nA, nB, Sab, reps = law[1], steps = law[2])
	expect_true(all(abs(z(draw, TRUE, 6, 6, S_pair)) < 4))
	expect_true(all(abs(z(draw, FALSE, 5, 8, S_bip)) < 4))
	expect_gt(max(abs(z(broken, TRUE, 6, 6, S_pair))), 10)
	expect_gt(max(abs(z(broken, FALSE, 5, 8, S_bip))), 10)
	# with a diagonal covariance the paired and the independent forms agree
	set.seed(2)
	a = rnorm(7)
	b = rnorm(7)
	set.seed(5)
	c1 = draw(a, b, diag(c(.7, 1.9)), TRUE)
	set.seed(5)
	c2 = draw(a, b, diag(c(.7, 1.9)), FALSE)
	expect_equal(c1, c2)
	# an undefined conditional leaves the effects
	expect_identical(draw(c(a, NA), b, diag(2), FALSE), 0)
})

test_that('lame() draws the shift without an intercept, records it, and the switch and an intercept omit it', {
	skip_on_cran()
	Y = abs_sim()
	for (dyn in c(FALSE, TRUE)) {
		on = abs_fit(Y, intercept = FALSE, dynamic_ab = dyn)
		expect_equal(unname(on$ab_shift_move['visits']), 60)
		expect_equal(unname(on$ab_shift_move['moved']), 60)
		off = abs_fit(Y, intercept = FALSE, dynamic_ab = dyn, prior = list(ab_shift_move = FALSE))
		expect_null(off$ab_shift_move)
		expect_null(abs_fit(Y, intercept = TRUE, dynamic_ab = dyn)$ab_shift_move)
	}
	Yb = abs_sim(8, 3, nB = 6)
	bip = abs_fit(Yb, intercept = FALSE, mode = 'bipartite')
	expect_equal(unname(bip$ab_shift_move['visits']), 60)
	# a symmetric fit has a = b and no shift
	Ys = lapply(Y, function(y) { y = (y + t(y)) / 2; y })
	expect_null(abs_fit(Ys, intercept = FALSE, symmetric = TRUE)$ab_shift_move)
	expect_error(abs_fit(Y, intercept = FALSE, prior = list(ab_shift_move = NA)), 'ab_shift_move')
})

# cached posterior summaries of four seeds with the move and four without
abs_cache = new.env()
abs_runs = function(cfg) {
	if (!is.null(abs_cache[[cfg]])) return(abs_cache[[cfg]])
	dyn = cfg == 'dynamic'
	bip = cfg == 'bipartite'
	Y = if (bip) abs_sim(15, 4, nB = 10) else abs_sim(15, 4)
	one = function(seed, on) {
		f = abs_fit(
			Y,
			burn = 200,
			nscan = if (dyn) 1500 else 1000,
			seed = seed,
			intercept = FALSE,
			dynamic_ab = dyn,
			mode = if (bip) 'bipartite' else 'unipartite',
			posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
			prior = if (on) list() else list(ab_shift_move = FALSE)
		)
		A = if (dyn) f$a_level_samples else f$a_samples
		B = if (dyn) f$b_level_samples else f$b_samples
		rms = function(M) apply(M, 2, function(x) sqrt(mean(x^2)))
		# bulk ess of the rms of the effects, the worst period under dynamic_ab
		ess_rms = function(D, S) {
			if (!dyn) return(posterior::ess_bulk(rms(S)))
			min(apply(D, 2, function(Mt) posterior::ess_bulk(rms(Mt))))
		}
		# a_i + b_j in the first period
		ab = if (dyn) {
			outer(rowMeans(f$a_draws[, 1, ]), rowMeans(f$b_draws[, 1, ]), '+')
		} else {
			outer(rowMeans(A), rowMeans(B), '+')
		}
		list(
			theta = c(
				rowMeans(A),
				rowMeans(B),
				sum_a = mean(colSums(A)),
				sum_b = mean(colSums(B)),
				va = mean(f$VC[, 'va']),
				vb = mean(f$VC[, 'vb'])
			),
			ab = ab,
			ess_a = ess_rms(f$a_draws, A),
			ess_b = ess_rms(f$b_draws, B)
		)
	}
	out = list(on = lapply(1:4, one, on = TRUE), off = lapply(1:4, one, on = FALSE))
	abs_cache[[cfg]] = out
	out
}

test_that('the move leaves the posterior of the effects and their variances unchanged', {
	skip_on_cran()
	skip_if_not_installed('posterior')
	for (cfg in c('static', 'dynamic')) {
		r = abs_runs(cfg)
		on = sapply(r$on, `[[`, 'theta')
		off = sapply(r$off, `[[`, 'theta')
		welch = (rowMeans(on) - rowMeans(off)) / sqrt(apply(on, 1, var) / 4 + apply(off, 1, var) / 4)
		# tight bound on the quantities the move acts on, wider on the actor means
		moved = c('sum_a', 'sum_b', 'va', 'vb')
		expect_true(all(abs(welch[moved]) < 3), info = cfg)
		expect_true(all(abs(welch) < 4), info = cfg)
		ab_on = Reduce(`+`, lapply(r$on, `[[`, 'ab')) / 4
		ab_off = Reduce(`+`, lapply(r$off, `[[`, 'ab')) / 4
		expect_lt(max(abs(ab_on - ab_off)), 0.05)
	}
})

test_that('the move mixes the size of the effects where the conditional draws cross the shift slowly', {
	skip_on_cran()
	skip_if_not_installed('posterior')
	for (cfg in c('dynamic', 'bipartite')) {
		r = abs_runs(cfg)
		for (nm in c('ess_a', 'ess_b')) {
			e_on = sapply(r$on, `[[`, nm)
			e_off = sapply(r$off, `[[`, nm)
			expect_gt(median(e_on), 3 * median(e_off), label = paste(cfg, nm))
		}
	}
	# static directed effects mix about as well with the move
	r = abs_runs('static')
	expect_gt(
		median(sapply(r$on, `[[`, 'ess_a')),
		0.7 * median(sapply(r$off, `[[`, 'ess_a'))
	)
})

# ame(): bipartite fits take the move, directed unipartite fits draw the effects jointly

ame_abs_sim = function(nA = 12, nB = 9, seed = 11) {
	set.seed(seed)
	a = rnorm(nA, 0.5, 1)
	b = rnorm(nB, -0.3, 1)
	y = outer(a, b, '+') + matrix(rnorm(nA * nB), nA, nB)
	dimnames(y) = list(paste0('r', seq_len(nA)), paste0('c', seq_len(nB)))
	y
}

ame_abs_fit = function(
	Y,
	burn = 20,
	nscan = 40,
	seed = 3,
	prior = list(),
	mode = 'bipartite',
	family = 'normal',
	R = 0,
	...
) {
	suppressWarnings(suppressMessages(ame(
		Y,
		family = family,
		mode = mode,
		R = R,
		burn = burn,
		nscan = nscan,
		odens = 1,
		seed = seed,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		prior = prior,
		...
	)))
}

test_that('ame() draws the shift for bipartite fits without an intercept and records it', {
	skip_on_cran()
	Y = ame_abs_sim()
	on = ame_abs_fit(Y, intercept = FALSE)
	expect_equal(unname(on$ab_shift_move['visits']), 60)
	expect_equal(unname(on$ab_shift_move['moved']), 60)
	off = ame_abs_fit(Y, intercept = FALSE, prior = list(ab_shift_move = FALSE))
	expect_null(off$ab_shift_move)
	expect_null(ame_abs_fit(Y, intercept = TRUE)$ab_shift_move)
	expect_null(ame_abs_fit(Y, intercept = FALSE, rvar = FALSE)$ab_shift_move)
	# the directed unipartite draw is joint
	Yu = abs_sim(10, 1)[[1]]
	expect_null(ame_abs_fit(Yu, intercept = FALSE, mode = 'unipartite')$ab_shift_move)
	expect_error(ame_abs_fit(Y, intercept = FALSE, prior = list(ab_shift_move = NA)), 'ab_shift_move')
	# a known prior name
	expect_no_warning(
		suppressMessages(ame(
			Y,
			mode = 'bipartite',
			R = 0,
			burn = 5,
			nscan = 10,
			odens = 1,
			chains = 1,
			verbose = FALSE,
			gof = FALSE,
			intercept = FALSE,
			prior = list(ab_shift_move = TRUE)
		)),
		message = 'ab_shift_move'
	)
})

test_that('every family of the bipartite ame() path without an intercept runs with the move', {
	skip_on_cran()
	set.seed(21)
	nA = 12; nB = 9; cap = 3
	z = outer(rnorm(nA, 0, .6), rnorm(nB, 0, .6), '+') + matrix(rnorm(nA * nB), nA, nB)
	nom = function(fam) {
		y = matrix(0, nA, nB)
		for (i in seq_len(nA)) {
			pos = which(z[i, ] > 0)
			k = min(cap, length(pos))
			if (k > 0) {
				top = pos[order(z[i, pos], decreasing = TRUE)[seq_len(k)]]
				y[i, top] = if (fam == 'cbin') 1 else rank(z[i, top])
			}
		}
		y
	}
	data = list(
		normal = z,
		binary = 1 * (z > 0),
		poisson = matrix(rpois(nA * nB, exp(0.5 + z / 2)), nA, nB),
		ordinal = matrix(as.integer(cut(z, c(-Inf, -.5, .5, 1.2, Inf))), nA, nB),
		cbin = nom('cbin'),
		frn = nom('frn')
	)
	for (fam in names(data)) {
		y = data[[fam]]
		dimnames(y) = list(paste0('r', seq_len(nA)), paste0('c', seq_len(nB)))
		args = list(y, family = fam, intercept = FALSE, R = 1)
		if (fam %in% c('cbin', 'frn')) args$odmax = cap
		f = do.call(ame_abs_fit, args)
		expect_equal(unname(f$ab_shift_move['visits']), 60, label = fam)
		expect_true(all(is.finite(f$APM)) && all(is.finite(f$BPM)), label = fam)
	}
})

# four seeds with the move and four without on one 12 x 9 network
ame_abs_cache = new.env()
ame_abs_runs = function() {
	if (!is.null(ame_abs_cache$r)) return(ame_abs_cache$r)
	Y = ame_abs_sim()
	one = function(seed, on) {
		f = ame_abs_fit(
			Y,
			burn = 200,
			nscan = 2000,
			seed = seed,
			intercept = FALSE,
			posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
			prior = if (on) list() else list(ab_shift_move = FALSE)
		)
		A = f$a_samples
		B = f$b_samples
		rms = function(M) apply(M, 2, function(x) sqrt(mean(x^2)))
		list(
			theta = c(
				sum_a = mean(colSums(A)),
				sum_b = mean(colSums(B)),
				va = mean(f$VC[, 'va']),
				vb = mean(f$VC[, 'vb'])
			),
			ab = outer(rowMeans(A), rowMeans(B), '+'),
			ess_a = posterior::ess_bulk(rms(A)),
			ess_b = posterior::ess_bulk(rms(B))
		)
	}
	ame_abs_cache$r = list(on = lapply(1:4, one, on = TRUE), off = lapply(1:4, one, on = FALSE))
	ame_abs_cache$r
}

test_that('the ame() move leaves the posterior of the bipartite effects unchanged', {
	skip_on_cran()
	skip_if_not_installed('posterior')
	r = ame_abs_runs()
	on = sapply(r$on, `[[`, 'theta')
	off = sapply(r$off, `[[`, 'theta')
	welch = (rowMeans(on) - rowMeans(off)) / sqrt(apply(on, 1, var) / 4 + apply(off, 1, var) / 4)
	expect_true(all(abs(welch) < 3))
	ab_on = Reduce(`+`, lapply(r$on, `[[`, 'ab')) / 4
	ab_off = Reduce(`+`, lapply(r$off, `[[`, 'ab')) / 4
	expect_lt(max(abs(ab_on - ab_off)), 0.05)
})

test_that('the ame() move mixes the size of the bipartite effects', {
	skip_on_cran()
	skip_if_not_installed('posterior')
	r = ame_abs_runs()
	for (nm in c('ess_a', 'ess_b')) {
		expect_gt(
			median(sapply(r$on, `[[`, nm)),
			2 * median(sapply(r$off, `[[`, nm)),
			label = nm
		)
	}
})
