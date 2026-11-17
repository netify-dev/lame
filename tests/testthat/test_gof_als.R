# gof() of an ALS fit: replicates at the point estimate, the ones gof_plot() draws

skip_on_cran()

gal_panel = function(n = 10L, N = 3L, seed = 12, changing = FALSE) {
	set.seed(seed)
	ids = sprintf('a%02d', seq_len(n))
	a = rnorm(n, 0, 0.5)
	Y = X = vector('list', N)
	for (t in seq_len(N)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(ids, ids, 'x1'))
		y = outer(a, a, '+') + 0.5 * x[, , 1] + matrix(rnorm(n * n), n, n)
		diag(y) = NA
		dimnames(y) = list(ids, ids)
		keep = if (changing && t == 2) ids[-1] else ids
		Y[[t]] = y[keep, keep]
		X[[t]] = x[keep, keep, , drop = FALSE]
	}
	names(Y) = names(X) = paste0('t', seq_len(N))
	list(Y = Y, X = X)
}

test_that('gof() of a cross-sectional ALS fit uses its point-estimate replicates', {
	d = gal_panel(N = 1L)
	fit = suppressWarnings(suppressMessages(ame_als(d$Y[[1]], Xdyad = d$X[[1]], R = 0, verbose = FALSE)))
	G = gof(fit, nsim = 8, verbose = FALSE)
	expect_identical(dim(G), c(9L, 5L))
	expect_identical(attr(G, 'replication_method'), 'point estimate')
	expect_equal(G[1, ], gof_stats(d$Y[[1]], mode = 'unipartite'), tolerance = 1e-12)
	sims = simulate(fit, nsim = 8, seed = 6886)
	manual = t(vapply(sims$Y, function(y) {
		y = if (is.list(y)) y[[1]] else y
		y[!is.finite(d$Y[[1]])] = NA
		gof_stats(y, mode = 'unipartite')
	}, numeric(5)))
	expect_equal(unname(G[-1, ]), unname(manual), tolerance = 1e-12)
	# the method dispatched through ame()
	fit2 = suppressWarnings(suppressMessages(ame(d$Y[[1]], Xdyad = d$X[[1]], R = 0, method = 'als', verbose = FALSE)))
	expect_true(all(is.finite(gof(fit2, nsim = 4, verbose = FALSE))))
})

test_that('gof() of a longitudinal ALS fit with changing actors runs per period', {
	d = gal_panel(changing = TRUE)
	fit = suppressWarnings(suppressMessages(lame(d$Y, Xdyad = d$X, R = 0, method = 'als', verbose = FALSE)))
	G = gof(fit, nsim = 5, verbose = FALSE)
	expect_identical(dim(G), c(6L, 5L))
	obs = Reduce('+', lapply(d$Y, gof_stats, mode = 'unipartite')) / 3
	expect_equal(G[1, ], obs, tolerance = 1e-12)
	expect_true(all(is.finite(G)))
})
