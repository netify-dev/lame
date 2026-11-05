test_that('unresolved factor precision fails explicitly without changing caller state', {
	U = array(0, c(3, 2, 2));V = array(1e10, dim(U));E = array(0, c(3, 3, 2));original = list(U = U, V = V, E = E)
	expect_error(
		lame:::rUV_dynamic_fc_cpp(U, V, E, .8, .5, 1, TRUE, FALSE),
		'Gaussian precision|chol.*decomposition failed'
	)
	expect_error(
		lame:::rUV_dynamic_snap_fc_cpp(U, V, E, .8, .5, 1, 2, .1, matrix(0, 3, 2), matrix(0, 3, 2), TRUE, FALSE),
		'Gaussian precision|chol.*decomposition failed'
	)
	expect_error(
		lame:::rUV_dynamic_t_fc_cpp(U, V, E, .8, .5, 1, 5, matrix(1, 3, 2), matrix(1, 3, 2), TRUE, FALSE),
		'Gaussian precision|chol.*decomposition failed'
	)
	expect_identical(list(U = U, V = V, E = E), original)
})

test_that('snap selection respects zero and unit prior probabilities', {
	set.seed(731876);U = array(rnorm(3 * 2 * 3) * 100, c(3, 2, 3));V = U;E = array(0, c(3, 3, 3))
	for (p in c(0, 1)) {
		z = lame:::rUV_dynamic_snap_fc_cpp(U, V, E, .8, .5, 1, 2, p, matrix(0, 3, 3), matrix(0, 3, 3), TRUE, FALSE)
		expect_true(all(z$delta_u[, -1] == p));expect_true(all(z$delta_v[, -1] == p))
		expect_true(all(is.finite(c(z$U, z$V))))
	}
})
