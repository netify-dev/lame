# gof() on panels whose actor set changes, period by period on the actors present
skip_on_cran()

gca_panel = function(n = 8L, N = 4L, family = 'normal', symmetric = FALSE, seed = 3) {
	set.seed(seed)
	ids = sprintf('a%02d', seq_len(n))
	Y = vector('list', N)
	for (t in seq_len(N)) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(ids, ids))
		if (family == 'binary') y = 1 * (y > 0.3)
		if (symmetric) y[lower.tri(y)] = t(y)[lower.tri(y)]
		diag(y) = NA
		# one actor absent from the second period and another from the last
		keep = setdiff(ids, c(if (t == 2) ids[1], if (t == N) ids[2]))
		Y[[t]] = y[keep, keep]
	}
	names(Y) = paste0('t', seq_len(N))
	Y
}
gca_bip = function(nA = 7L, nB = 6L, N = 3L, seed = 4) {
	set.seed(seed)
	rows = sprintf('r%02d', seq_len(nA))
	cols = sprintf('c%02d', seq_len(nB))
	Y = lapply(seq_len(N), function(t) matrix(rnorm(nA * nB), nA, nB, dimnames = list(rows, cols)))
	Y[[2]] = Y[[2]][-1, , drop = FALSE]
	Y[[3]] = Y[[3]][, -2, drop = FALSE]
	names(Y) = paste0('t', seq_len(N))
	Y
}
gca_fit = function(Y, ...) suppressWarnings(suppressMessages(lame(
	Y,
	burn = 20,
	nscan = 60,
	odens = 2,
	seed = 7,
	chains = 1,
	verbose = FALSE,
	plot = FALSE,
	prior = list(chain_extensions = 0),
	...
)))
# observed statistics the sampler stored, averaged over the periods
gca_sampler_obs = function(fit) vapply(fit$GOF, function(g) mean(g[, 1]), numeric(1))

test_that('gof() of a directed panel with changing actors matches the sampler per period', {
	Y = gca_panel()
	for (dyn in list(list(), list(dynamic_ab = TRUE), list(R = 1, dynamic_uv = TRUE))) {
		fit = do.call(gca_fit, c(list(Y), modifyList(list(R = 0, family = 'normal'), dyn)))
		G = gof(fit, nsim = 6, verbose = FALSE)
		expect_identical(dim(G), c(7L, 5L))
		expect_equal(G[1, names(fit$GOF)], gca_sampler_obs(fit), tolerance = 1e-12)
		expect_true(all(is.finite(G)))
		# the observed periods on their own actors give the same check
		expect_equal(gof(fit, Y = Y, nsim = 6, verbose = FALSE), G, tolerance = 1e-12)
	}
	# each replicate is evaluated on the actors present in its period
	fit = gca_fit(Y, R = 0, family = 'normal')
	sims = simulate(fit, nsim = 6, seed = 6886, type = 'posterior')
	manual = vapply(sims$Y, function(yy) mean(vapply(seq_along(Y), function(t) {
		y = yy[[t]][rownames(Y[[t]]), colnames(Y[[t]])]
		y[!is.finite(Y[[t]])] = NA
		gof_stats(y, mode = 'unipartite')[['sd.rowmean']]
	}, numeric(1))), numeric(1))
	expect_equal(unname(gof(fit, nsim = 6, verbose = FALSE)[-1, 'sd.rowmean']), manual, tolerance = 1e-12)
})

test_that('gof() of a symmetric panel with changing actors matches the sampler per period', {
	Y = gca_panel(family = 'binary', symmetric = TRUE)
	fit = gca_fit(Y, R = 0, family = 'binary', symmetric = TRUE)
	G = gof(fit, nsim = 6, verbose = FALSE)
	expect_equal(G[1, names(fit$GOF)], gca_sampler_obs(fit), tolerance = 1e-12)
	expect_true(all(is.finite(G)))
})

test_that('gof() of a bipartite panel with changing rows and columns matches the sampler per period', {
	Y = gca_bip()
	fit = gca_fit(Y, mode = 'bipartite', R_row = 1, R_col = 1, family = 'normal')
	G = gof(fit, nsim = 6, verbose = FALSE)
	expect_identical(colnames(G), names(fit$GOF))
	expect_equal(G[1, ], gca_sampler_obs(fit), tolerance = 1e-12)
	expect_true(all(is.finite(G)))
	# an actor outside the predictive grid with observed cells is still refused
	bad = Y
	rownames(bad[[1]])[1] = 'unknown'
	expect_error(gof(fit, Y = bad, nsim = 2, verbose = FALSE), class = 'lame_gof_alignment')
})
