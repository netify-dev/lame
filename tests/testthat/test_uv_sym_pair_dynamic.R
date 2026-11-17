# the pair gauge of symmetric dynamic position paths against their weights (R/uv_gauge_dynamic.R)
skip_on_cran()

usp_mut = function(fn, from, to) {
	f = get(fn, envir = asNamespace('lame'))
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# gauge-sensitive statistics of the paths: column norms and correlation over all periods and in each period
usp_stat = function(X) {
	S = crossprod(matrix(aperm(X, c(1, 3, 2)), ncol = dim(X)[2]))
	per = unlist(lapply(seq_len(dim(X)[3]), function(t) {
		St = crossprod(matrix(X[, , t], dim(X)[1]))
		c(log(diag(St)), St[1, 2] / sqrt(St[1, 1] * St[2, 2]))
	}))
	x = c(log(diag(S)), cov2cor(S)[upper.tri(S)], per)
	c(x, x^2)
}

# unit-variance AR(1) paths from their prior, moved three times; the statistics keep their law
usp_law = function(ndraw, seed, L = c(1.3, -0.7), n = 6L, N = 3L, rho = 0.8, max_abs = Inf, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_sym_pair_dynamic_draw
	set.seed(seed)
	R = length(L)
	out = NULL
	for (r in seq_len(ndraw)) {
		repeat {
			U = array(0, c(n, R, N))
			U[, , 1] = rnorm(n * R)
			for (t in 2:N) U[, , t] = rho * U[, , t - 1] + sqrt(1 - rho^2) * rnorm(n * R)
			if (max(abs(U)) <= max_abs) break
		}
		o = list(U = U)
		for (k in 1:3) o = draw(o$U, L, rho, max_abs = max_abs)
		for (t in seq_len(N)) {
			p0 = U[, , t] %*% diag(L) %*% t(U[, , t])
			p1 = o$U[, , t] %*% diag(L) %*% t(o$U[, , t])
			if (max(abs(p1 - p0)) > 1e-8) stop('the product moved')
		}
		if (max(abs(o$U)) > max_abs) stop('the bound was left')
		s1 = usp_stat(o$U) - usp_stat(U)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(nrow(out)))))
}

test_that('the hyperbolic and circular pair maps keep U L U\' and the law of the paths', {
	ns = asNamespace('lame')
	for (l in list(c(1.4, -0.6), c(0.8, 0.3), c(-2, -0.5))) {
		A = ns$.lame_uv_sym_pair_map(0.37, l)
		expect_equal(A %*% diag(l) %*% t(A), diag(l), tolerance = 1e-12)
		expect_equal(det(A), 1, tolerance = 1e-12)
	}
	expect_lt(usp_law(2000, 9301), 4.5)
	expect_lt(usp_law(2000, 9302, L = c(1, 0.5)), 4.5)
	expect_lt(usp_law(2000, 9303, L = c(1, -0.6, 0.4)), 4.5)
	expect_lt(usp_law(2000, 9304, max_abs = 2), 4.5)
	# the checks fail without the transition into the period and with the innovations of the tail held
	expect_gt(usp_law(2000, 9305, draw = usp_mut('.lame_uv_sym_pair_dynamic_draw', '- 2 * sum(A * C)', '')), 5)
	expect_gt(usp_law(2000, 9306, draw = usp_mut('.lame_uv_sym_pair_dynamic_draw', 'M <- M + crossprod(e)', 'M <- M')), 5)
})

test_that('lame() runs the pair gauge with the switch on and a symmetric rho .9 panel passes the two-chain check', {
	set.seed(9311)
	n = 20L; N = 4L; rho = 0.9
	rn = sprintf('a%02d', 1:n)
	U = array(0, c(n, 2, N))
	U[, , 1] = rnorm(2 * n, 0, 0.8)
	for (t in 2:N) U[, , t] = rho * U[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * n, 0, 0.8)
	Y = lapply(1:N, function(t) {
		e = matrix(rnorm(n * n), n)
		e[lower.tri(e)] = t(e)[lower.tri(e)]
		y = U[, , t] %*% diag(c(1, -1)) %*% t(U[, , t]) + 0.5 * e
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		y
	})
	names(Y) = paste0('t', 1:N)
	fit = function(prior, chains = 2, nscan = 1500) suppressMessages(lame(
		Y, family = 'normal', R = 2L, symmetric = TRUE, dynamic_uv = TRUE, burn = 300, nscan = nscan, odens = 15,
		seed = 6886, verbose = FALSE, gof = FALSE, plot = FALSE, prior = prior, chains = chains, cores = 1
	))
	f = fit(list())
	expect_s3_class(f, 'lame')
	expect_gt(f$uv_sym_pair_move[['moved']], 0)
	expect_null(fit(list(uv_sym_pair_move = FALSE), chains = 1, nscan = 150)$uv_sym_pair_move)
	expect_error(fit(list(uv_sym_pair_move = NA), chains = 1, nscan = 150), 'uv_sym_pair_move')
})
