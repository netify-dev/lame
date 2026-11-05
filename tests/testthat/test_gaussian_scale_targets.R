test_that('bipartite factor draws retain small and large prior variances', {
	set.seed(492161)
	U = array(0, c(2, 2, 1));V = array(0, c(3, 2, 1));E = array(0, c(2, 3, 1));G = matrix(0, 2, 2)
	for (sig in c(1e-8, 1e4)) {
		z = replicate(4000, lame:::rUV_dynamic_bip_fc_cpp(U, V, E, G, 0, sig, 1)$U[1, , 1]) / sig
		expect_lt(max(abs(rowMeans(z))), 6 / sqrt(ncol(z)))
		expect_lt(max(abs(rowMeans(z^2) - 1)), 6 * sqrt(2 / ncol(z)))
	}
	bad = U;bad[1] = NaN;copy = bad
	expect_error(lame:::rUV_dynamic_bip_fc_cpp(bad, V, E, G, .5, .2, 1), 'Invalid bipartite')
	expect_identical(bad, copy)
})

test_that('reciprocal coefficient conditionals retain small observation variance', {
	set.seed(563175)
	Z = matrix(0, 2, 2);diag(Z) = NA;ve = 1e-8;rho = .5
	args = list(
		Xdyn_list = list(matrix(1, 4, 1)),
		Xstat_list = list(matrix(0, 4, 0)),
		Z_list = list(Z),
		offset_list = list(matrix(0, 2, 2)),
		beta_static = numeric(0),
		rho_by_coef = .5,
		sigma_by_coef = .2,
		Lambda = matrix(1),
		beta0_mean = 0,
		beta0_cov = matrix(1),
		s2 = ve,
		dyad_rho = rho,
		bipartite = FALSE,
		symmetric = FALSE,
		use_dyad_rho = TRUE
	)
	variance = 1 / (1 + 2 / (ve * (1 + rho)))
	d = replicate(4000, do.call(lame:::sample_beta_dynamic_cpp, args)$path[1, 1]) / sqrt(variance)
	expect_lt(abs(mean(d)), 6 / sqrt(length(d)))
	expect_lt(abs(mean(d^2) - 1), 6 * sqrt(2 / length(d)))
	args$beta0_cov = matrix(-1)
	expect_error(do.call(lame:::sample_beta_dynamic_cpp, args), 'coefficient predictive covariance')
})

test_that('coefficient smoothing retains tiny positive innovation noise', {
	set.seed(748912)
	args = list(
		Xdyn_list = rep(list(matrix(0, 4, 1)), 2),
		Xstat_list = rep(list(matrix(0, 4, 0)), 2),
		Z_list = rep(list(matrix(NA_real_, 2, 2)), 2),
		offset_list = rep(list(matrix(0, 2, 2)), 2),
		beta_static = numeric(0),
		rho_by_coef = 1,
		sigma_by_coef = 1e-8,
		Lambda = matrix(1),
		beta0_mean = 0,
		beta0_cov = matrix(1),
		s2 = 1,
		dyad_rho = 0,
		bipartite = FALSE,
		symmetric = FALSE,
		use_dyad_rho = FALSE
	)
	d = replicate(4000, diff(do.call(lame:::sample_beta_dynamic_cpp, args)$path[, 1])) / 1e-8
	expect_lt(abs(mean(d)), 6 / sqrt(length(d)))
	expect_lt(abs(mean(d^2) - 1), 6 * sqrt(2 / length(d)))
})

test_that('additive scales and symmetric factor loadings do not acquire floors', {
	set.seed(917425)
	zero = matrix(0, 3, 3);rho = .4;scale = 1e-16;shape = 2 + nrow(zero) * ncol(zero)
	sd = replicate(
		4000,
		lame:::.lame_sample_sigma_ab(
			zero,
			zero,
			rho,
			1,
			FALSE,
			prior_shape = 2,
			prior_scale = scale
		)
	)
	precision = scale / sd^2
	expect_lt(abs(mean(precision) - shape), 6 * sqrt(shape / length(sd)))
	expect_true(all(sd < 1e-6))
	U = array(1, c(3, 1, 1));E = array(0, c(3, 3, 1));ve = 1e-16;variance = 1 / (1 + 3 / ve)
	z = replicate(4000, lame:::.lame_sym_dyn_L_draw(E, U, 0, ve)) / sqrt(variance)
	expect_lt(abs(mean(z)), 6 / sqrt(length(z)))
	expect_lt(abs(mean(z^2) - 1), 6 * sqrt(2 / length(z)))
})
