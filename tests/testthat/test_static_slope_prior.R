# the prior of a static slope does not depend on whether another coefficient has a path

test_that('the static block keeps the marginal of the panel prior, free of s2', {
	set.seed(21)
	X = matrix(rnorm(60 * 4), 60, 4)
	X[, 1] = 1
	P = crossprod(X) / 50
	S = c(2, 4)
	P0 = lame:::.lame_static_panel_precision(P, S)
	expect_equal(P0, solve(solve(P)[S, S]), tolerance = 1e-12)
	expect_identical(lame:::.lame_static_panel_precision(P, 1:4), (P + t(P)) / 2)
	# every static coefficient of a unipartite fit takes P0 whatever s2 is
	expect_identical(lame:::.lame_static_prior_prec(P0, c(9, 9), S, 0.3), P0)
	expect_identical(lame:::.lame_static_s2_terms(P0, c(1, 2), c(9, 9), S), c(df = 0, ss = 0))
	# bipartite static coefficients stay N(0, g s2) each
	expect_identical(lame:::.lame_static_prior_prec(NULL, c(2, 8), 1:2, 0.5), diag(1 / (c(2, 8) * 0.5), 2))
	expect_identical(
		lame:::.lame_static_s2_terms(NULL, c(1, 2), c(2, 8), 1:2),
		lame:::.lame_s2_beta_terms(c(1, 2), c(2, 8))
	)
	# a static intercept moves under its conditional given the other static coefficients
	P3 = lame:::.lame_static_panel_precision(P, 1:3)
	b = c(0.4, -1, 2)
	ic = lame:::.lame_static_icpt_prior(b, 1L, P3, rep(9, 3), 1:3, 1)
	V = solve(P3)
	expect_equal(ic$prec, 1 / (V[1, 1] - V[1, -1] %*% solve(V[-1, -1], V[-1, 1]))[1, 1], tolerance = 1e-12)
	expect_equal(ic$mean, drop(V[1, -1] %*% solve(V[-1, -1], b[-1])), tolerance = 1e-12)
	lq = lame:::.lame_actor_move_icpt(ic, NULL, NULL)
	expect_equal(unname(lq['l']), drop(P3 %*% b)[1], tolerance = 1e-12)
	expect_equal(unname(lq['q']), P3[1, 1], tolerance = 1e-12)
	# the scale terms and the ridge state read the same precision
	path = matrix(rnorm(6), 3, 2)
	terms = lame:::.lame_path_scale_terms(b[2:3], path, c(1, 0.5), c(0.3, 0.4), diag(2), diag(2), c(9, 9), 2, P0)
	base = lame:::.lame_path_scale_terms(numeric(0), path, c(1, 0.5), c(0.3, 0.4), diag(2), diag(2), numeric(0), 2)
	expect_equal(terms$q - base$q, drop(crossprod(b[2:3], P0 %*% b[2:3])), tolerance = 1e-12)
	st = lame:::.lame_cov_ridge_dyn_state(4, S, c(1, 3), b[2:3], path, c(9, 9), 2, diag(c(5, 6)), TRUE, P0)
	expect_identical(st$P0[S, S], P0)
	expect_identical(st$P0[1, 1], 1 / 5)
})

lsp_panel = function(n = 10L, N = 3L, family = 'normal', seed = 5, scale_y = 1) {
	set.seed(seed)
	ids = sprintf('a%02d', seq_len(n))
	Y = X = vector('list', N)
	for (t in seq_len(N)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(ids, ids, 'x'))
		y = 0.2 * t + 0.5 * x[, , 1] + matrix(rnorm(n * n), n, n)
		if (family == 'binary') y = 1 * (y > 0)
		diag(y) = NA
		dimnames(y) = list(ids, ids)
		Y[[t]] = scale_y * y
		X[[t]] = x
	}
	names(Y) = names(X) = paste0('t', seq_len(N))
	list(Y = Y, X = X)
}
lsp_fit = function(d, ...) suppressWarnings(suppressMessages(lame(
	d$Y,
	Xdyad = d$X,
	R = 0,
	burn = 10,
	nscan = 20,
	odens = 2,
	chains = 1,
	verbose = FALSE,
	gof = FALSE,
	plot = FALSE,
	prior = list(chain_extensions = 0),
	...
)))

test_that('a long-run mean or path level takes the prior of its coefficient when static', {
	set.seed(22)
	X = matrix(rnorm(60 * 4), 60, 4)
	X[, 1] = 1
	P = crossprod(X) / 50
	v = lame:::.lame_static_panel_var(P, c(2, 4))
	expect_equal(v, diag(solve(P))[c(2, 4)], tolerance = 1e-12)
	# the variance coefficient 2 has under the marginal of a static block holding it
	P0 = lame:::.lame_static_panel_precision(P, 1:3)
	expect_equal(solve(P0)[2, 2], v[1], tolerance = 1e-12)
})

test_that('a unipartite fit resolves the same static slope prior with and without paths', {
	skip_on_cran()
	for (family in c('normal', 'binary')) {
		d = lsp_panel(family = family)
		static = lsp_fit(d, family = family)
		dyn_ab = lsp_fit(d, family = family, dynamic_ab = TRUE)
		dyn_icpt = lsp_fit(d, family = family, dynamic_beta = 'intercept')
		expect_identical(dyn_ab$g, static$g)
		expect_identical(dyn_icpt$g, static$g)
		expect_identical(static$coefficient_prior$kind, 'panel_design')
		expect_match(dyn_ab$coefficient_prior$description, 'panel design prior')
	}
	# the default g of the normal family carries the unit of the response
	d = lsp_panel()
	d10 = lsp_panel(scale_y = 10)
	expect_equal(lsp_fit(d10, dynamic_ab = TRUE)$g, 100 * lsp_fit(d, dynamic_ab = TRUE)$g, tolerance = 1e-10)
})

test_that('a bipartite fit resolves the same static slope prior with and without paths', {
	skip_on_cran()
	set.seed(8)
	rows = sprintf('r%02d', 1:8)
	cols = sprintf('c%02d', 1:6)
	X = lapply(1:3, function(t) array(rnorm(48, 0, 3), c(8, 6, 1), dimnames = list(rows, cols, 'x')))
	Y = lapply(1:3, function(t) matrix(0.4 * X[[t]][, , 1] + rnorm(48), 8, 6, dimnames = list(rows, cols)))
	names(X) = names(Y) = paste0('t', 1:3)
	fit = function(...) suppressWarnings(suppressMessages(lame(
		Y,
		Xdyad = X,
		mode = 'bipartite',
		R_row = 0,
		R_col = 0,
		family = 'normal',
		burn = 10,
		nscan = 20,
		odens = 2,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = list(chain_extensions = 0),
		...
	)))
	static = fit()
	dyn = fit(dynamic_ab = TRUE)
	expect_equal(dyn$g, static$g, tolerance = 1e-12)
	# the slope's g is per unit of its covariate
	expect_lt(static$g[2], static$g[1] / 5)
})
