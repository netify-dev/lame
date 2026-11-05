test_that("dynamic bipartite warm starts preserve the fitted G interaction", {
	skip_on_cran()
	set.seed(1135008L);nr = 24L;nc = 20L;nt = 6L
	a = matrix(rnorm(nr * nt, sd = .4), nr);b = matrix(rnorm(nc * nt, sd = .4), nc)
	u = matrix(rnorm(nr, sd = .6), nr);v = matrix(rnorm(nc, sd = .6), nc)
	Y = X = vector('list', nt)
	for (t in 1:nt) {
		x = matrix(rnorm(nr * nc), nr);eta = -.5 + .5 * x + outer(a[, t], b[, t], '+') + tcrossprod(u, v)
		y = matrix(rbinom(nr * nc, 1, pnorm(eta)), nr)
		dimnames(y) = list(paste0('r', 1:nr), paste0('c', 1:nc))
		Y[[t]] = y;X[[t]] = array(x, c(nr, nc, 1), dimnames = c(dimnames(y), list('x')))
	}
	f = lame:::lame_dynamic_als(
		Y,
		Xdyad = X,
		R = 1,
		family = 'binary',
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		mode = 'bipartite',
		verbose = FALSE,
		seed = 441L,
		bias_correct = 'none'
	)
	seed = lame:::.ae_static_seed(f)
	expect_identical(seed$G, f$G)
	expect_gt(abs(f$G[1, 1] - 1), .1)
	g = ame_als_refit(f, max_iter = 1000, tol = 1e-10)
	expect_equal(g$irls_trace[1], f$objective, tolerance = 1e-9)
	expect_lte(g$objective, f$objective + 1e-8 * max(1, abs(f$objective)))
})


test_that("unequal bipartite ranks remain valid warm starts", {
	seed = list(
		a = numeric(5),
		b = numeric(6),
		coefficients = c(intercept = 0),
		U = matrix(0, 5, 1),
		V = matrix(0, 6, 2),
		G = matrix(c(2, -3), 1, 2)
	)
	prep = list(nr = 5L, nc = 6L)
	usable = lame:::.lda_usable_start(seed, prep, 2L, FALSE, R_row = 1L, R_col = 2L)
	expect_identical(usable, seed)
})

test_that("Poisson warm seeds reconstruct the fitted count means", {
	set.seed(813121L);nr = 18L;nc = 16L
	Y = lapply(1:4, function(t) matrix(rpois(nr * nc, exp(.2 + rnorm(nr * nc, sd = .7))), nr))
	f = suppressWarnings(lame:::lame_dynamic_als(
		Y,
		R = 0,
		family = 'poisson',
		dynamic_ab = TRUE,
		mode = 'bipartite',
		verbose = FALSE,
		bias_correct = 'none'
	))
	seed = lame:::.ae_static_seed(f)
	for (t in 1:4) {
		eta = seed$coef_path['intercept', t] + outer(seed$a_path[, t], seed$b_path[, t], '+')
		expect_equal(exp(eta), unname(f$fitted[[t]]), tolerance = 1e-10, ignore_attr = TRUE)
	}
	g = suppressWarnings(ame_als_refit(f, max_iter = 200, tol = 1e-8))
	expect_equal(g$irls_trace[1], f$objective, tolerance = 1e-8)
	expect_lte(g$objective, f$objective + 1e-8 * max(1, abs(f$objective)))
})
