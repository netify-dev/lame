# the public routing expressions, evaluated on controlled states
skip_on_cran()

.find_conditional_route = function(x, predicate) {
	if (!is.call(x)) return(NULL)
	if (predicate(x)) return(x)
	for (y in as.list(x)) {
		z = tryCatch(.find_conditional_route(y, predicate), error = function(e) NULL)
		if (!is.null(z)) return(z)
	}
	NULL
}

test_that("disabled AME margins are absent from the joint Gaussian law", {
	skip_on_cran()
	ns = asNamespace("lame")
	expr = .find_conditional_route(
		body(get("ame_unipartite", ns)),
		function(x)
			identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("tmp")) &&
			grepl("rbeta_ab_fc", paste(deparse(x), collapse = " "), fixed = TRUE)
	)
	expect_false(is.null(expr))
	draw = function(rvar, cvar, x, Sab, s2) {
		n = nrow(x); Z = matrix(0, n, n); UV_cache = Z
		X_precomp = precomputeX(array(x, c(n, n, 1L)))
		g = n * n; rho = 0; symmetric = FALSE; sym_gaussian_context = NULL
		eval(expr); c(tmp$beta, if (rvar) tmp$a else if (cvar) tmp$b)
	}
	environment(draw) = list2env(list(expr = expr), parent = ns)
	n = 3L; s2 = 1e-12; samples = 2000L
	for (flags in list(c(FALSE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) {
		rvar = flags[1]; cvar = flags[2]
		x = if (rvar) matrix(rep(c(1, 2, -1), each = n), n) else
			matrix(rep(c(1, 2, -1), n), n)
		Sab = diag(c(if (rvar) .8 * s2 else 1e-8, if (cvar) .6 * s2 else 1e-8))
		q = 1L + n * as.integer(rvar || cvar)
		D = matrix(0, n * n, q); D[, 1] = as.numeric(x)
		ij = expand.grid(i = seq_len(n), j = seq_len(n))
		if (q > 1L) for (cell in seq_len(n * n))
			D[cell, 1L + if (rvar) ij$i[cell] else ij$j[cell]] = 1
		P = c(
			sum(x * x) * (1 + 1e-6) / (n * n),
			if (q > 1L) rep(1 / if (rvar) Sab[1, 1] else Sab[2, 2], n)
		)
		S = solve(crossprod(D) / s2 + diag(P, q))
		set.seed(751 + q + as.integer(rvar))
		z = t(matrix(replicate(samples, draw(rvar, cvar, x, Sab, s2)), nrow = q))
		se = sqrt((outer(diag(S), diag(S)) + S * S) / (samples - 1L))
		expect_lt(max(abs(cov(z) - S) / se), 7)
		expect_lt(max(abs(colMeans(z)) / sqrt(diag(S) / samples)), 7)
	}
})

test_that("a valid beta-only AME fit emits no actor-rank warning", {
	Y = matrix(c(0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0), 5)
	expect_warning(
		ame(
			Y,
			family = "binary",
			R = 0,
			rvar = FALSE,
			cvar = FALSE,
			dcor = FALSE,
			burn = 1,
			nscan = 4,
			odens = 1,
			gof = FALSE,
			verbose = FALSE,
			seed = 752
		),
		NA
	)
})

test_that("dynamic actors use the prior-preserving bipartite beta conditional", {
	ns = asNamespace("lame")
	expr = .find_conditional_route(
		body(get("lame", ns)),
		function(x)
			identical(x[[1]], as.name("if")) &&
			identical(paste(deparse(x[[2]]), collapse = ""), "p_bip > 0") &&
			grepl("resid_b", paste(deparse(x), collapse = " "), fixed = TRUE)
	)
	expect_false(is.null(expr))
	draw = function(x) {
		nA = 2L; nB = 3L; N = 1L; p_bip = 2L; beta = numeric(p_bip)
		Z = array(0, c(nA, nB, N)); EZ_uv_only = Z
		a_mat = matrix(0, nA, N); b_mat = matrix(0, nB, N)
		Xlist = list(x); s2 = .7; g = 12; eval(expr); beta
	}
	environment(draw) = list2env(list(expr = expr), parent = ns)
	set.seed(753)
	z = replicate(2000, draw(array(1e4, c(2, 3, 2))))
	contrast = (z[1, ] - z[2, ]) / sqrt(2)
	# identical likelihood columns leave this direction at its N(0, g * s2) prior
	expect_lt(abs(var(contrast) / (12 * .7) - 1) / sqrt(2 / 1999), 7)
	set.seed(754); before = .Random.seed
	expect_error(draw(array(1e16, c(2, 3, 2))), "unresolved|rescal", ignore.case = TRUE)
	expect_identical(.Random.seed, before)
})
