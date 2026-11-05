test_that('temporal density excludes missing cells and retains bipartite diagonals', {
	y = matrix(c(NA, 1, NA, NA), 2)
	expect_equal(.pick_temporal_stat('density')(y), 1)
	square = diag(1, 3)
	expect_equal(.pick_temporal_stat('density', 'bipartite')(square), 1 / 3)
	expect_equal(.pick_temporal_stat('mean', 'bipartite')(square), 1 / 3)
	expect_equal(.pick_temporal_stat('density')(square), 0)
	expect_equal(.pick_temporal_stat('mean')(square), 0)
})

test_that('temporal transitivity agrees with explicit observed-triple enumeration', {
	set.seed(491706)
	for (i in 1:20) {
		y = matrix(sample(c(NA, 0, 1), 25, replace = TRUE), 5);diag(y) = NA
		complete_paths = closures = 0
		for (a in 1:5) for (b in setdiff(1:5, a)) for (c in setdiff(1:5, c(a, b))) {
			v = c(y[a, b], y[b, c], y[a, c])
			if (all(is.finite(v)) && all(v[1:2] != 0)) {complete_paths = complete_paths + 1;closures = closures + (v[3] != 0)}
		}
		expected = if (complete_paths) closures / complete_paths else NA_real_
		expect_equal(.pick_temporal_stat('transitivity')(y), expected)
	}
})

test_that('temporal prediction uses joint rows and identical observed-cell masks', {
	set.seed(617225)
	Y = lapply(
		1:3,
		function(t) matrix(
			rnorm(42, mean = t / 2),
			6,
			7,
			dimnames = list(paste0('r', 1:6), paste0('c', 1:7))
		)
	)
	Y[[1]][1:2, 2] = NA;Y[[2]] = Y[[2]][c(4, 2, 6, 1), c(7, 3, 4, 1, 6), drop = FALSE]
	f = lame(
		Y,
		mode = 'bipartite',
		family = 'normal',
		R = 0,
		burn = 20,
		nscan = 40,
		odens = 5,
		posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
		gof = FALSE,
		seed = 87325,
		verbose = FALSE,
		plot = FALSE
	)
	set.seed(73952);rng = .Random.seed
	g = gof_temporal(f, stat = 'mean', n_rep = 8, seed = 4826)
	expect_identical(.Random.seed, rng)
	sim = simulate(f, nsim = 8, seed = 4826, type = 'posterior')
	expected = t(vapply(sim$Y, function(yy) vapply(1:3, function(t) {
		observed = f$Y[, , t];x = yy[[t]][rownames(observed), colnames(observed), drop = FALSE]
		x[!is.finite(observed)] = NA_real_;mean(x, na.rm = TRUE)
	}, numeric(1)), numeric(3)))
	expect_equal(g$stat_rep_by_t, expected, tolerance = 1e-12)
	expect_equal(unname(g$stat_obs_by_t), unname(apply(f$Y, 3, mean, na.rm = TRUE)))
	expect_equal(g$slope_rep, apply(expected, 1, function(x) unname(coef(lm(x ~ I(1:3)))[2])), tolerance = 1e-12)
	expect_identical(g$draw_index, sim$draw_index)
	expect_identical(g$replication_method, 'joint posterior')
	f$a_samples = f$a_draws = NULL
	expect_error(gof_temporal(f, stat = 'mean', n_rep = 2, seed = 4826), class = 'lame_predictive_draws')
	expect_error(gof_temporal(f, stat = 'reciprocity', n_rep = 2), class = 'lame_gof_uninformative')
	expect_error(gof_temporal(f, stat = 'density', n_rep = 2), class = 'lame_gof_uninformative')
	expect_error(gof_temporal(f, n_rep = 1), 'at least 2')
})

test_that('constant observed trends and tied predictive slopes are compared correctly', {
	y = matrix(0, 3, 3);diag(y) = NA
	f = structure(list(Y = rep(list(y), 3), family = 'normal', mode = 'unipartite'), class = 'lame')
	make_sim = function(slopes) list(
		Y = lapply(slopes, function(s) lapply(c(-1, 0, 1), function(t) y + s * t)),
		method = 'joint posterior',
		draw_index = seq_along(slopes)
	)
	for (scale in c(1, 1e-100)) {
		testthat::with_mocked_bindings({
			g = gof_temporal(f, stat = 'mean', n_rep = 4)
			expect_equal(g$slope_obs, 0);expect_false(g$degenerate);expect_equal(g$p_pp, 1)
			expect_equal(g$slope_rep, c(-1, 0, 0, 1) * scale)
		}, simulate = function(object, nsim, seed, type) {
			expect_identical(type, 'posterior');make_sim(c(-1, 0, 0, 1) * scale)
		}, .package = 'stats')
	}
	testthat::with_mocked_bindings({
		expect_warning(g <- gof_temporal(f, stat = 'mean', n_rep = 4), class = 'lame_gof_uninformative')
		expect_true(g$degenerate);expect_true(is.na(g$p_pp));expect_match(g$reason, 'identical')
	}, simulate = function(...) make_sim(rep(0, 4)), .package = 'stats')
})
