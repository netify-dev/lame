# the scale split of the static factors of a directed network (R/uv_scale_split.R)
skip_on_cran()

usp_mut = function(from, to) {
	f = lame:::.lame_uv_split_draw
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# rows from the matrix-t prior (covariance from IW(df, Psi), rows N(0, Suv)); one step of the move keeps their law
# statistics: squared factor norms, their inner product, the log determinant of the cross product, and products of these
usp_law = function(
	ndraw,
	seed,
	n = 7L,
	R = 2L,
	hier = TRUE,
	df = 2L * R + 2L,
	Psi = NULL,
	draw = lame:::.lame_uv_split_draw
) {
	set.seed(seed)
	Pg = Psi
	if (length(Psi) == 1L) Pg = Psi * diag(2L * R)
	if (is.null(Psi)) { M = crossprod(matrix(rnorm((2 * R)^2), 2 * R)); Psi = (M / (2 * R) + diag(.6, 2 * R)) * df }
	out = matrix(0, ndraw, 10)
	stat = function(U, V) { S = crossprod(cbind(U, V)); a = sum(U^2); b = sum(V^2); c(
		a,
		b,
		sum(U * V),
		as.numeric(determinant(S)$modulus),
		a * b,
		log(a / b),
		U[1, 1] * V[2, 1],
		log(a / b)^2,
		a^2,
		b^2
	) }
	for (r in seq_len(ndraw)) {
		if (hier) { Pg = if (is.null(Pg)) Psi else Pg; Suv = solve(
			rWishart(1, df, solve(Pg))[, , 1]
		) ; W = matrix(rnorm(n * 2 * R), n) %*% chol(Suv) }
		# the covariance held at n s I, s = Psi / df the unit of the prior
		else W = matrix(rnorm(n * 2 * R, 0, sqrt(n * (if (length(Psi) == 1L) Psi else mean(diag(Psi))) / df)), n)
		U = W[, 1:R, drop = FALSE]; V = W[, R + 1:R, drop = FALSE]
		o = draw(U, V, df, Psi, hier)
		out[r, ] = stat(o$U, o$V) - stat(U, V)
		# the product is untouched
		if (max(abs(o$U %*% t(o$V) - U %*% t(V))) > 1e-9) stop('the product moved')
	}
	z = colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))
	ev = eigen(cov(out), symmetric = TRUE)
	m = drop(crossprod(ev$vectors, colMeans(out)))
	T2 = ndraw * sum(m^2 / ev$values)
	list(maxz = max(abs(z)), p = pchisq(T2, 10, lower.tail = FALSE))
}

test_that('one step of the scale split leaves the law of the factors under the matrix-t prior unchanged, for a matrix scale, a scalar scale and the fixed covariance', {
	skip_on_cran()
	for (cfg in list(list(8001, 2L, TRUE), list(8002, 1L, TRUE), list(8003, 3L, TRUE), list(8004, 2L, FALSE))) {
		r = usp_law(7000, cfg[[1]], R = cfg[[2]], hier = cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	r = usp_law(7000, 8005, Psi = 9)
	expect_lt(r$maxz, 4.5)
	# the fixed covariance of a response in large units (s = 400): the scale the sampler holds, not n I
	r = usp_law(7000, 8006, hier = FALSE, Psi = 400 * 6)
	expect_lt(r$maxz, 4.5)
	expect_gt(r$p, 1e-3)
})

test_that('the check fails for a split drawn with the wrong degrees of freedom, scale, Jacobian or covariance', {
	skip_on_cran()
	expect_gt(
		usp_law(
			7000,
			8011,
			draw = usp_mut('-(n + df) * sum(log(diag(ch)))', '-(n + df + 4) * sum(log(diag(ch)))')
		)$maxz,
		5
	)
	expect_gt(usp_law(7000, 8012, draw = usp_mut('Psi + S * outer(d, d)', '2 * Psi + S * outer(d, d)'))$maxz, 5)
	expect_gt(
		usp_law(
			7000,
			8013,
			draw = usp_mut('-(n + df) * sum(log(diag(ch)))', '-(n + df) * sum(log(diag(ch))) + 2 * u')
		)$maxz,
		5
	)
	expect_gt(
		usp_law(
			7000,
			8014,
			draw = usp_mut('-(n + df) * sum(log(diag(ch)))', '-(n + df) * sum(log(diag(ch))) - 0.5 * n * u^2')
		)$maxz,
		5
	)
	expect_gt(usp_law(7000, 8015, hier = FALSE, draw = usp_mut('(2 * n * s)', '(8 * n * s)'))$maxz, 5)
	expect_gt(usp_law(7000, 8016, hier = FALSE, Psi = 400 * 6, draw = usp_mut('(2 * n * s)', '(2 * n)'))$maxz, 5)
})

test_that('the move leaves a degenerate state alone', {
	ns = asNamespace('lame')
	U = matrix(c(1, NA, 2, 3), 2)
	expect_identical(ns$.lame_uv_split_draw(U, U, 6, 1)$c, 1)
	expect_identical(ns$.lame_uv_split_draw(matrix(0, 1, 2), matrix(0, 1, 2), 6, 1)$c, 1)
})

usp_panel = function(n = 24L, Tn = 2L, seed = 3L) {
	set.seed(seed)
	rn = sprintf('a%02d', seq_len(n))
	U = matrix(rnorm(n * 2, 0, .8), n); V = matrix(rnorm(n * 2, 0, .8), n)
	Y = lapply(
		seq_len(Tn),
		function(t) { y = U %*% t(V) + matrix(rnorm(n * n), n); diag(y) = NA; dimnames(y) = list(rn, rn); y }
	)
	names(Y) = paste0('p', seq_len(Tn))
	Y
}
usp_fit = function(Y, ..., prior = list(), burn = 20, nscan = 60) suppressWarnings(suppressMessages(lame(
	Y,
	family = 'normal',
	R = 2L,
	burn = burn,
	nscan = nscan,
	odens = 1,
	seed = 5,
	verbose = FALSE,
	gof = FALSE,
	plot = FALSE,
	prior = prior,
	posterior_opts = list(save_UV = TRUE),
	...
)))

test_that('lame() and ame() draw the split by default, record it, and the switch omits it', {
	Y = usp_panel()
	on = usp_fit(Y)
	off = usp_fit(Y, prior = list(uv_split_move = FALSE))
	expect_equal(unname(on$uv_split_move['visits']), 80)
	expect_gt(unname(on$uv_split_move['moved']), 70)
	expect_null(off$uv_split_move)
	expect_false(identical(on$BETA, off$BETA))
	expect_error(usp_fit(Y, prior = list(uv_split_move = 'no')), 'uv_split_move')
	a = suppressWarnings(
		suppressMessages(
			ame(
				Y[[1]],
				family = 'normal',
				R = 2L,
				burn = 20,
				nscan = 60,
				odens = 1,
				seed = 5,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE
			)
		)
	)
	expect_equal(unname(a$uv_split_move['visits']), 80)
	a0 = suppressWarnings(
		suppressMessages(
			ame(
				Y[[1]],
				family = 'normal',
				R = 2L,
				burn = 20,
				nscan = 60,
				odens = 1,
				seed = 5,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				prior = list(uv_split_move = FALSE)
			)
		)
	)
	expect_null(a0$uv_split_move)
	expect_error(
		suppressWarnings(
			ame(
				Y[[1]],
				family = 'normal',
				R = 2L,
				burn = 5,
				nscan = 10,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				prior = list(uv_split_move = 'no')
			)
		),
		'uv_split_move'
	)
	# moving positions draw the split of their paths; hold off the drift check so the record is not extended
	local_mocked_bindings(
		.lame_chain_moving = function(fit, ...) data.frame(parameter = character(0), text = character(0)),
		.package = 'lame'
	)
	d = usp_fit(Y, dynamic_uv = TRUE)
	expect_equal(unname(d$uv_split_move['visits']), 80)
})

test_that('the split of the scale mixes at nominal speed with the move', {
	skip_on_cran()
	Y = usp_panel(30L, 2L, 4L)
	ess = function(x) { x = x - mean(x); a = acf(
		x,
		lag.max = 300,
		plot = FALSE
	)$acf[-1]; k = which(a < .05)[1]; if (is.na(k)) k = 300; length(x) / (1 + 2 * sum(a[seq_len(k - 1L)])) }
	split = function(f) vapply(
		seq_len(dim(f$U_samples)[3]),
		function(s) log(sqrt(sum(f$U_samples[, , s]^2)) / sqrt(sum(f$V_samples[, , s]^2))),
		0
	)
	on = split(usp_fit(Y, burn = 200, nscan = 2000))
	off = split(usp_fit(Y, prior = list(uv_split_move = FALSE), burn = 200, nscan = 2000))
	expect_gt(ess(on), 3 * ess(off))
	expect_gt(ess(on), 500)
	se = sqrt(var(on) / ess(on) + var(off) / ess(off))
	expect_lt(abs(mean(on) - mean(off)) / se, 4.5)
})

# the column scales and shears of the static gauge U A, V A^-T, and the split of a symmetric fit
ugd_mut = function(fn, from, to) {
	f = get(fn, envir = asNamespace('lame'))
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# three steps of the move from a prior draw keep its law: log norms and correlations of the columns, the split
# and the log determinant, with their squares; returns the largest z of the mean change
ugd_law = function(ndraw, seed, n = 7L, R = 2L, hier = TRUE, df = 2L * R + 2L, Psi = df * 1.3, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_gauge_draw
	set.seed(seed)
	Pg = if (length(Psi) == 1L) Psi * diag(2L * R) else Psi
	stat = function(U, V) {
		S = crossprod(cbind(U, V))
		C = cov2cor(S)
		x = c(log(diag(S)), C[upper.tri(C)], log(sum(U^2) / sum(V^2)), as.numeric(determinant(S[1:R, 1:R])$modulus))
		c(x, x^2)
	}
	out = NULL
	for (r in seq_len(ndraw)) {
		if (hier) {
			Suv = solve(rWishart(1, df, solve(Pg))[, , 1])
			W = matrix(rnorm(n * 2 * R), n) %*% chol(Suv)
		} else W = matrix(rnorm(n * 2 * R, 0, sqrt(n * Psi / df)), n)
		U = W[, 1:R, drop = FALSE]; V = W[, R + 1:R, drop = FALSE]
		o = list(U = U, V = V)
		for (k in 1:3) o = draw(o$U, o$V, df, Psi, hier)
		if (max(abs(o$U %*% t(o$V) - U %*% t(V))) > 1e-8) stop('the product moved')
		out = rbind(out, stat(o$U, o$V) - stat(U, V))
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))))
}

test_that('the column scales and shears leave the law of the factors under the matrix-t prior unchanged', {
	skip_on_cran()
	expect_lt(ugd_law(3000, 8101), 4.5)
	expect_lt(ugd_law(3000, 8102, R = 3L), 4.5)
	expect_lt(ugd_law(3000, 8103, hier = FALSE), 4.5)
	# the checks fail for the wrong degrees of freedom, a wrong measure on the shear and the wrong held scale
	expect_gt(ugd_law(3000, 8104, draw = ugd_mut(
		'.lame_uv_gauge_draw',
		'-0.5 * (n + df)',
		'-0.5 * (n + df + 4)'
	)), 5)
	expect_gt(ugd_law(3000, 8105, draw = ugd_mut(
		'.lame_uv_gauge_draw',
		'lp(crossprod(B, S %*% B))',
		'lp(crossprod(B, S %*% B)) - abs(h)'
	)), 5)
	expect_gt(ugd_law(3000, 8106, hier = FALSE, draw = ugd_mut('.lame_uv_gauge_draw', '(2 * n * s)', '(8 * n * s)')), 5)
})

test_that('the gauge draw keeps the product, leaves rank one and a degenerate state alone', {
	set.seed(8107)
	U = matrix(rnorm(20), 10); V = matrix(rnorm(20), 10)
	o = lame:::.lame_uv_gauge_draw(U, V, 6, 6)
	expect_equal(o$U %*% t(o$V), U %*% t(V), tolerance = 1e-10)
	expect_equal(o$moved, 4L)
	# rank one is the scale split alone and takes no random numbers
	seed_before = .Random.seed
	expect_identical(lame:::.lame_uv_gauge_draw(U[, 1, drop = FALSE], V[, 1, drop = FALSE], 4, 4)$moved, 0L)
	expect_identical(.Random.seed, seed_before)
	U[2, 1] = NA
	expect_identical(lame:::.lame_uv_gauge_draw(U, V, 6, 6)$moved, 0L)
})

# symmetric factors: column precision Gamma(1, 1/2) or held at 1 / n, eigen-weights N(0, l_var)
usd_law = function(ndraw, seed, n = 8L, R = 2L, hier = TRUE, l_var = 1.7, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_sym_split_draw
	set.seed(seed)
	stat = function(U, V) {
		x = c(log(colSums(U^2)), log(colSums(V^2)), sum(U[, 1] * V[, 1]))
		c(x, x^2)
	}
	out = NULL
	for (r in seq_len(ndraw)) {
		sdu = if (hier) 1 / sqrt(rexp(R, 0.5)) else rep(sqrt(n), R)
		U = sweep(matrix(rnorm(n * R), n), 2, sdu, '*')
		V = sweep(U, 2, rnorm(R, 0, sqrt(l_var)), '*')
		o = list(U = U, V = V)
		for (k in 1:3) o = draw(o$U, o$V, l_var, hier)
		if (max(abs(o$U %*% t(o$V) - U %*% t(V))) > 1e-8) stop('the product moved')
		out = rbind(out, stat(o$U, o$V) - stat(U, V))
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))))
}

test_that('the split of a symmetric fit leaves the law of the factors and eigen-weights unchanged', {
	skip_on_cran()
	expect_lt(usd_law(3000, 8111), 4.5)
	expect_lt(usd_law(3000, 8112, hier = FALSE), 4.5)
	expect_lt(usd_law(3000, 8113, R = 1L, l_var = 25), 4.5)
	# the checks fail for a wrong Jacobian and a wrong exponent of the integrated precision
	expect_gt(usd_law(3000, 8114, draw = ugd_mut('.lame_uv_sym_split_draw', '(n - 2) * u', '(n - 1) * u')), 5)
	expect_gt(usd_law(3000, 8115, draw = ugd_mut('.lame_uv_sym_split_draw', '-0.5 * (n + 2)', '-0.5 * (n + 4)')), 5)
})

test_that('a dynamic_rho fit with moving actor effects and static positions draws the whole gauge and settles', {
	skip_on_cran()
	# the sweep cell uni_normal_ab_drho_none_int at a small size: directed normal panel, a dyadic correlation path
	set.seed(8121)
	n = 16L; N = 4L; rn = sprintf('a%02d', seq_len(n))
	U = matrix(rnorm(n * 2, 0, .8), n); V = matrix(rnorm(n * 2, 0, .8), n)
	Y = lapply(seq_len(N), function(t) {
		rt = c(.2, .4, .5, .3)[t]
		z1 = matrix(rnorm(n * n), n); z2 = matrix(rnorm(n * n), n)
		e = z1; lo = lower.tri(e); e[lo] = rt * t(z1)[lo] + sqrt(1 - rt^2) * z2[lo]
		y = -.5 + rnorm(n, 0, .3) + U %*% t(V) + .5 * e
		diag(y) = NA; dimnames(y) = list(rn, rn); y
	})
	names(Y) = paste0('t', seq_len(N))
	calls = 0L
	gauge_draw = lame:::.lame_uv_gauge_draw
	trace_gauge = function(...) { calls <<- calls + 1L; gauge_draw(...) }
	local_mocked_bindings(.lame_uv_gauge_draw = trace_gauge, .package = 'lame')
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = 'normal',
		R = 2L,
		dynamic_ab = TRUE,
		dynamic_rho = TRUE,
		burn = 300,
		nscan = 1200,
		odens = 10,
		seed = 7,
		chains = 2,
		cores = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	expect_s3_class(fit, 'lame')
	expect_gt(calls, 0L)
	rh = fit$rhat
	expect_true('uv_split' %in% rh$parameter)
	expect_lt(rh$rhat[rh$parameter == 'uv_split'], 1.05)
})

test_that('with the switch off the gauge and the symmetric split are never drawn', {
	stop_draw = function(...) stop('drawn with the switch off')
	local_mocked_bindings(
		.lame_uv_gauge_draw = stop_draw,
		.lame_uv_sym_split_draw = stop_draw,
		.lame_uv_split_draw = stop_draw,
		.package = 'lame'
	)
	Y = usp_panel(12L, 2L, 5L)
	expect_null(usp_fit(Y, prior = list(uv_split_move = FALSE), burn = 5, nscan = 10)$uv_split_move)
	Ys = lapply(Y, function(y) { y = (y + t(y)) / 2; y })
	expect_null(usp_fit(Ys, symmetric = TRUE, prior = list(uv_split_move = FALSE), burn = 5, nscan = 10)$uv_split_move)
	expect_s3_class(
		suppressWarnings(suppressMessages(ame(
			Ys[[1]],
			family = 'normal',
			R = 2L,
			symmetric = TRUE,
			burn = 5,
			nscan = 10,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(uv_split_move = FALSE)
		))),
		'ame'
	)
})

test_that('symmetric lame() and ame() draw the split of their columns by default', {
	Ys = lapply(usp_panel(12L, 2L, 6L), function(y) { y = (y + t(y)) / 2; y })
	f = usp_fit(Ys, symmetric = TRUE, burn = 10, nscan = 20)
	expect_equal(unname(f$uv_split_move['visits']), 30)
	expect_gt(unname(f$uv_split_move['moved']), 20)
	a = suppressWarnings(suppressMessages(ame(
		Ys[[1]],
		family = 'normal',
		R = 2L,
		symmetric = TRUE,
		burn = 10,
		nscan = 20,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	expect_equal(unname(a$uv_split_move['visits']), 30)
})

test_that('the shears of the bipartite gauge keep the law of the positions and the interaction, and the check has power', {
	skip_on_cran()
	bip_law = function(ndraw, seed, nA = 6L, nB = 8L, RA = 2L, RB = 3L, g_var = 2, draw = NULL) {
		if (is.null(draw)) draw = lame:::.lame_bip_uv_g_scale_step
		set.seed(seed)
		stat = function(U, V, G) {
			x = c(log(colSums(U^2)), log(colSums(V^2)), sum(U[, 1] * U[, 2]), sum(V[, 1] * V[, 3]), sum(G[1, ] * G[2, ]), G[1, 1])
			c(x, x^2)
		}
		out = NULL
		for (r in seq_len(ndraw)) {
			U = matrix(rnorm(nA * RA), nA); V = matrix(rnorm(nB * RB), nB); G = matrix(rnorm(RA * RB, 0, sqrt(g_var)), RA)
			o = list(U = U, V = V, G = G)
			for (k in 1:3) o = draw(o$U, o$V, o$G, g_var = g_var)
			out = rbind(out, stat(o$U, o$V, o$G) - stat(U, V, G))
		}
		max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))))
	}
	expect_lt(bip_law(3000, 8131), 4.5)
	expect_gt(bip_law(3000, 8132, draw = ugd_mut(
		'.lame_bip_uv_g_scale_step',
		'- sum(Gm[k, ] * Gm[l, ])/g_var)',
		'+ sum(Gm[k, ] * Gm[l, ])/g_var)'
	)), 5)
})
