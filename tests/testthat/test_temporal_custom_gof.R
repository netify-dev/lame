.temporal_gof_fixture = function() {
	set.seed(817254)
	Y = lapply(
		1:3,
		function(t) matrix(
			rnorm(56, mean = t / 3),
			7,
			8,
			dimnames = list(paste0('r', 1:7), paste0('c', 1:8))
		)
	)
	Y[[1]][1, 2] = NA;Y[[2]] = Y[[2]][c(3, 7, 4, 2, 6), c(8, 3, 5, 1, 6, 2), drop = FALSE]
	Y[[2]][2, 4] = NA;Y[[3]][3, 1] = NA;names(Y) = c('early', 'middle', 'late')
	f = lame(
		Y,
		mode = 'bipartite',
		R = 0,
		family = 'normal',
		burn = 20,
		nscan = 40,
		odens = 5,
		posterior_opts = posterior_options(save_ab = TRUE, thin_ab = 1),
		gof = FALSE,
		verbose = FALSE,
		plot = FALSE,
		seed = 32617
	)
	list(fit = f, Y = Y)
}

test_that('custom GOF uses identical per-period matrix summaries and exact nsim', {
	z = .temporal_gof_fixture();f = z$fit;Y = z$Y
	stat = function(y) { stopifnot(is.matrix(y));c(
		cells = length(y),
		observed = sum(is.finite(y)),
		average = mean(y, na.rm = TRUE)
	) }
	G = gof(f, Y = Y, nsim = 30, custom_gof = stat, verbose = FALSE)
	expected_obs = Reduce('+', lapply(Y, stat)) / 3
	expect_equal(unname(G[1, names(expected_obs)]), unname(expected_obs), tolerance = 1e-12)
	expect_equal(nrow(G), 31L);expect_true(all(G[, 'cells'] == expected_obs['cells']))
	expect_true(all(G[, 'observed'] == expected_obs['observed']))
	sims = simulate(f, nsim = 30, seed = 6886, type = 'posterior')
	manual = vapply(sims$Y, function(yy) mean(vapply(1:3, function(t) {
		y = yy[[t]][rownames(Y[[t]]), colnames(Y[[t]]), drop = FALSE]
		y[!is.finite(Y[[t]])] = NA_real_;mean(y, na.rm = TRUE)
	}, numeric(1))), numeric(1))
	expect_equal(unname(G[-1, 'average']), manual, tolerance = 1e-12)
	expect_identical(attr(G, 'draw_index'), sims$draw_index)
	expect_true(anyDuplicated(attr(G, 'draw_index')) > 0)
	L = gof(f, Y = Y, nsim = 4, custom_gof = list(cells = function(y) length(y)), verbose = FALSE)
	expect_true(all(L[, 'cells'] == expected_obs['cells']))
})

test_that('named actor and period permutations preserve observed-cell GOF', {
	z = .temporal_gof_fixture();Y = z$Y
	permuted = lapply(Y, function(y) y[rev(seq_len(nrow(y))), rev(seq_len(ncol(y))), drop = FALSE])
	permuted = permuted[c(3, 1, 2)]
	a = gof(z$fit, Y = Y, nsim = 5, verbose = FALSE)
	b = gof(z$fit, Y = permuted, nsim = 5, verbose = FALSE)
	expect_equal(a, b, tolerance = 1e-12)
	bad = Y;rownames(bad[[1]])[1] = 'unknown'
	expect_error(gof(z$fit, Y = bad, nsim = 2, verbose = FALSE), class = 'lame_gof_alignment')
	expect_error(gof(z$fit, Y = Y[1:2], nsim = 2, verbose = FALSE), class = 'lame_gof_alignment')
	expect_error(gof(z$fit, nsim = 0, verbose = FALSE), 'positive integer')
	expect_error(gof(z$fit, nsim = 1.5, verbose = FALSE), 'positive integer')
})

test_that('GOF alignment keeps singleton axes and refuses positional recycling', {
	x = array(1:6, c(1, 3, 2), dimnames = list('row', c('a', 'b', 'c'), c('first', 'second')))
	y = .lame_gof_periods(x)
	expect_identical(dim(y[[1]]), c(1L, 3L));expect_identical(names(y), c('first', 'second'))
	expect_equal(.lame_gof_align_replicate(matrix(1:6, 2, 3), matrix(0, 2, 3)), matrix(1:6, 2, 3))
	expect_error(.lame_gof_align_replicate(matrix(1:6, 2, 3), matrix(0, 1, 3)), class = 'lame_gof_alignment')
})


test_that('explicit rank zero overrides compatibility coordinate placeholders', {
	z = .temporal_gof_fixture();f = z$fit
	expect_equal(f$R, 0);expect_gt(dim(f$U)[2], 0)
	state = .joint_predictive_state(f)
	expect_null(state$U);expect_null(state$V)
	row = state$index[1];ab = match(row, state$ab_index)
	independent = outer(f$a_samples[, ab], f$b_samples[, ab], '+') + f$BETA[row, 1]
	expect_equal(unname(.joint_predictor(f, state, row, 1)), unname(independent), tolerance = 1e-12)
	sim = simulate(f, nsim = 3, seed = 219731, type = 'posterior')
	expect_identical(sim$method, 'joint posterior');expect_length(sim$Y, 3)
})
