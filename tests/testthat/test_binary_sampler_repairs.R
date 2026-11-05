skip_on_cran()

.binary_repair_panel = function() {
	withr::local_seed(71440)
	actors = paste0('n', c(7, 2, 5, 1, 4, 6, 3))
	lapply(seq_len(3), function(t) {
		z = matrix(rnorm(49), 7, dimnames = list(actors, actors))
		diag(z) = NA_real_
		z
	})
}

test_that('disabled dynamic actor margins remain exactly zero', {
	y = .binary_repair_panel()
	for (family in c('binary', 'normal')) for (margin in c('row', 'col')) {
		f = lame(
			if (family == 'binary') lapply(y, function(z) 1 * (z > 0)) else y,
			family = family,
			R = 0,
			dynamic_ab = TRUE,
			rvar = margin == 'row',
			cvar = margin == 'col',
			seed = 7886,
			burn = 20,
			nscan = 40,
			odens = 10,
			gof = FALSE,
			verbose = FALSE
		)
		sv = f$start_vals$checkpoint_state$values
		fields = if (margin == 'row') c('b', 'b_lev', 'b_mat') else c('a', 'a_lev', 'a_mat')
		expect_true(all(vapply(fields, function(k) all(sv[[k]] == 0), TRUE)))
		expect_true(all(unlist(f$tryErrorChecks) == 0))
	}
	expect_identical(
		lame:::.lame_ab_cond_prior(diag(c(.7, 0)), 1L, numeric(5)),
		list(mean = NULL, var = .7)
	)
	expect_error(lame:::.lame_ab_cond_prior(matrix(c(1, .1, .1, 0), 2), 1L, numeric(5)), 'disabled')
})

test_that('the exact symmetric block agrees with independent Gaussian conditioning', {
	withr::local_seed(71441)
	n = 5L; p = 2L
	x = array(rnorm(n * n * p), c(n, n, p))
	Q = matrix(c(1, .2, .2, .8), 2)
	cx = lame:::.lame_symmetric_gaussian_context(list(x), Q, actor = TRUE, center = TRUE)
	ij = which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	H = qr.Q(qr(contr.helmert(n)))
	D = t(vapply(seq_len(nrow(ij)), function(k) {
		i = ij[k, 1]; j = ij[k, 2]
		c((x[i, j, ] + x[j, i, ]) / 2, H[i, ] + H[j, ])
	}, numeric(p + n - 1)))
	V = matrix(0, p + n - 1, p + n - 1)
	V[1:p, 1:p] = solve(Q); V[p + seq_len(n - 1), p + seq_len(n - 1)] = diag(.7, n - 1)
	z = matrix(rnorm(n * n), n); z = (z + t(z)) / 2
	A = V %*% t(D); marginal = diag(nrow(D)) + D %*% A
	reference_mean = drop(A %*% solve(marginal, z[upper.tri(z)]))
	reference_cov = V - A %*% solve(marginal, t(A))
	observed = lame:::.lame_symmetric_gaussian_parameters(array(z, c(n, n, 1)), .7, cx)
	expect_equal(observed$mean, reference_mean, tolerance = 1e-11)
	expect_equal(chol2inv(observed$root), reference_cov, tolerance = 1e-11)
	diag(z) = 1e9
	changed = lame:::.lame_symmetric_gaussian_parameters(array(z, c(n, n, 1)), .7, cx)
	expect_identical(observed, changed)
})

test_that('new directed binary moves run and expose numerical failures', {
	y = lapply(.binary_repair_panel(), function(z) 1 * (z > 0))
	withr::local_seed(71442)
	before = .Random.seed
	f = lame(
		y,
		family = 'binary',
		R = 0,
		dynamic_ab = TRUE,
		seed = 7886,
		burn = 20,
		nscan = 100,
		odens = 10,
		gof = FALSE,
		verbose = FALSE
	)
	expect_identical(.Random.seed, before)
	expect_equal(unname(f$joint_levels_move['attempts']), 12)
	expect_equal(unname(f$joint_levels_move['numerical_failures']), 0)
	expect_true(all(unlist(f$tryErrorChecks) == 0))
	a = ame(
		y[[1]],
		family = 'binary',
		R = 0,
		dcor = TRUE,
		seed = 7886,
		burn = 25,
		nscan = 100,
		odens = 10,
		gof = FALSE,
		verbose = FALSE
	)
	expect_equal(unname(a$rho_collapsed_move['attempts']), 20)
	expect_true(all(unlist(a$tryErrorChecks) == 0))
	expect_identical(.Random.seed, before)
})

test_that('symmetric binary fits use the exact block with static and dynamic actors', {
	y = lapply(.binary_repair_panel(), function(z) 1 * (z + t(z) > 0))
	for (dynamic in c(FALSE, TRUE)) {
		f = lame(
			y,
			family = 'binary',
			symmetric = TRUE,
			dynamic_ab = dynamic,
			R = 1,
			seed = 7886,
			burn = 30,
			nscan = 100,
			odens = 10,
			gof = FALSE,
			verbose = FALSE
		)
		expect_true(all(is.finite(f$BETA)))
		expect_true(all(unlist(f$tryErrorChecks) == 0))
		expect_null(f$joint_levels_move)
	}
	a = ame(
		y[[1]],
		family = 'binary',
		symmetric = TRUE,
		R = 1,
		seed = 7886,
		burn = 30,
		nscan = 100,
		odens = 10,
		gof = FALSE,
		verbose = FALSE
	)
	expect_true(all(is.finite(a$BETA)))
	expect_true(all(unlist(a$tryErrorChecks) == 0))
	expect_null(a$rho_collapsed_move)
})
