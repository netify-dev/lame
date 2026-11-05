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
	# symmetric, bipartite and moving positions have no split to draw
	d = usp_fit(Y, dynamic_uv = TRUE)
	expect_null(d$uv_split_move)
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
