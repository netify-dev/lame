skip_on_cran()

.weighted_cycle_reference = function(Y) {
	observed = is.finite(Y); v = Y[observed]; E = Y
	if (length(v) < 4L || stats::sd(v) == 0) return(NA_real_)
	E[observed] = (v - mean(v)) / stats::sd(v)
	products = numeric()
	for (i in seq_len(nrow(Y) - 1L)) for (k in (i + 1L):nrow(Y))
		for (j in seq_len(ncol(Y) - 1L)) for (l in (j + 1L):ncol(Y)) {
			z = E[cbind(c(i, i, k, k), c(j, l, j, l))]
			if (all(is.finite(z))) products = c(products, prod(z))
		}
	if (length(products)) mean(products) else NA_real_
}

test_that('Gaussian rectangle dependence agrees with explicit enumeration', {
	set.seed(735821)
	for (k in seq_len(20L)) {
		Y = matrix(rnorm(30), 5L, 6L);Y[sample.int(30, 5)] = NA
		x = gof_stats_bipartite(Y, family = 'normal')['four.cycle.dep']
		expect_equal(unname(x), .weighted_cycle_reference(Y), tolerance = 1e-12)
		expect_equal(unname(x), .lame_bipartite_four_cycle_dependence(t(Y)), tolerance = 1e-12)
		expect_equal(unname(x), .lame_bipartite_four_cycle_dependence(7 - 3 * Y), tolerance = 1e-12)
		expect_equal(unname(x), .lame_bipartite_four_cycle_dependence(Y * 1e200), tolerance = 1e-12)
		expect_equal(unname(x), .lame_bipartite_four_cycle_dependence(Y * 1e-200), tolerance = 1e-12)
	}
})

test_that('fixed nonzero masks no longer guarantee a Gaussian predictive match', {
	Y = matrix(c(-2, -1, 1, 2, 3, -3, 4, -4, 2, -2, 1, -1), 3L, 4L)
	Z = matrix(c(1, 4, -2, -3, 1, 2, 3, 4, -4, -2, -1, -3), 3L, 4L)
	Y[1, 1] = Z[1, 1] = NA
	expect_equal(gof_stats_bipartite(Y)['four.cycles'], gof_stats_bipartite(Z)['four.cycles'])
	a = gof_stats(Y, mode = 'bipartite', family = 'normal')
	b = gof_stats(Z, mode = 'bipartite', family = 'normal')
	expect_true('four.cycle.dep' %in% names(a));expect_false('four.cycles' %in% names(a))
	expect_gt(abs(a['four.cycle.dep'] - b['four.cycle.dep']), .01)
	empty = matrix(NA_real_, 3, 4);empty[cbind(c(1, 1, 2, 3), c(1, 2, 3, 4))] = 1:4
	expect_true(is.na(gof_stats_bipartite(empty, family = 'normal')['four.cycle.dep']))
	expect_true(is.na(gof_stats_bipartite(matrix(1, 3, 4), family = 'normal')['four.cycle.dep']))
	B = matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0), 3, 4)
	expect_identical(gof_stats_bipartite(B, family = 'binary'), gof_stats_bipartite(B))
})

test_that('Gaussian bipartite fits and post-hoc checks use observed weighted rectangles', {
	set.seed(912354);Y = matrix(rnorm(42), 6, 7, dimnames = list(paste0('r', 1:6), paste0('c', 1:7)));Y[1, 1] = NA
	f = ame(
		Y,
		mode = 'bipartite',
		R = 1,
		family = 'normal',
		burn = 20,
		nscan = 40,
		odens = 5,
		posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE, thin_UV = 1, thin_ab = 1),
		gof = TRUE,
		verbose = FALSE,
		seed = 79136
	)
	expect_true('four.cycle.dep' %in% colnames(f$GOF))
	expect_equal(unname(f$GOF[1, 'four.cycle.dep']), .weighted_cycle_reference(Y), tolerance = 1e-12)
	expect_gt(sd(f$GOF[-1, 'four.cycle.dep']), 0)
	G = gof(f, nsim = 6, verbose = FALSE, custom_gof = function(y) c(n_observed = sum(is.finite(y))))
	expect_true('four.cycle.dep' %in% colnames(G));expect_true(all(G[, 'n_observed'] == 41))
	expect_identical(
		attr(G, 'replication_method'),
		'joint posterior'
	);expect_length(attr(G, 'draw_index'), nrow(G) - 1L)
	missing = f;missing$U_samples = missing$V_samples = missing$G_samples = NULL
	expect_error(gof(missing, nsim = 2, verbose = FALSE), class = 'lame_predictive_draws')
	expect_equal(unname(G[1, 'four.cycle.dep']), .weighted_cycle_reference(Y), tolerance = 1e-12)
	expect_error(gof_plot(f, statistics = 'four.cycles'), class = 'lame_gof_uninformative')
	expect_s3_class(gof_plot(f, statistics = 'four.cycle.dep'), 'ggplot')
})

test_that('longitudinal weighted GOF preserves periods and uses complete rectangles', {
	set.seed(216783);Y = lapply(
		1:2,
		function(t) matrix(rnorm(42), 6, 7, dimnames = list(paste0('r', 1:6), paste0('c', 1:7)))
	)
	Y[[1]][1, 1] = NA
	f = lame(
		Y,
		mode = 'bipartite',
		R = 0,
		family = 'normal',
		burn = 20,
		nscan = 40,
		odens = 5,
		posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
		gof = TRUE,
		verbose = FALSE,
		plot = FALSE,
		seed = 71931
	)
	expect_true('four.cycle.dep' %in% names(f$GOF))
	expect_equal(unname(f$GOF$four.cycle.dep[, 1]), vapply(Y, .weighted_cycle_reference, numeric(1)), tolerance = 1e-12)
	expect_true(all(apply(f$GOF$four.cycle.dep[, -1, drop = FALSE], 1, sd) > 0))
	p = gof_plot(f, statistics = 'four.cycle.dep');expect_s3_class(p, 'ggplot')
	expect_equal(sort(unique(p$data$time)), 1:2)
	expect_s3_class(gof_plot(f, type = 'static', statistics = 'four.cycle.dep'), 'ggplot')
})

test_that('ALS weighted GOF compares the same observed cells', {
	set.seed(116273);Y = matrix(rnorm(42), 6, 7, dimnames = list(paste0('r', 1:6), paste0('c', 1:7)));Y[1, 1] = NA
	f = ame_als(Y, mode = 'bipartite', R = 0, family = 'normal', verbose = FALSE)
	p = gof_plot(f, nsim = 5, seed = 92736)
	sims = simulate(f, nsim = 5, seed = 92736)$Y
	expected = vapply(sims, function(y) {y[!is.finite(Y)] = NA_real_;.weighted_cycle_reference(y)}, numeric(1))
	z = p$data[p$data$statistic == 'four.cycle.dep', ]
	expect_equal(z$value, unname(expected), tolerance = 1e-12)
	expect_equal(z$observed, rep(.weighted_cycle_reference(Y), 5L), tolerance = 1e-12)
	expect_true(all(is.finite(z$value)));expect_gt(sd(z$value), 0)
})
