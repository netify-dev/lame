test_that('factor balancing remains finite at subnormal squared norms', {
	skip_on_cran()
	set.seed(811800L);nr = 10L;nc = 9L;rank = 2L
	W = matrix(runif(nr * nc, .3, .7), nr, nc);W[sample(length(W), 20)] = 0;W[1:2, ] = 0;W[, nc] = 0
	Y = matrix(rnorm(nr * nc), nr, nc);obs = W > 0
	U = matrix(rnorm(nr * rank), nr, rank);V = matrix(rnorm(nc * rank), nc, rank)
	for (exponent in c(0, 1, 10, seq(140, 163, by = .5))) for (order in 1:2) {
		u = U * 10^-exponent;v = V * 10^-exponent;if (order == 2) v = v[, 2:1]
		new = .ae_lowrank_em_als(Y, W, obs, rank, u, v, c(.167, .168), c(.168, .167), 1, 6)
		native = als_factor_update_cpp(Y, W, u, v, c(.167, .168), c(.168, .167), 1, 6, 15L, 1e-6)
		expect_true(all(is.finite(c(new$U, new$V, new$Omat, new$sig2_u, new$sig2_v))))
		expect_equal(new$Omat, tcrossprod(native$U, native$V), tolerance = 1e-12)
		penalty = function(z) sum(sweep(z$U^2, 2, z$sig2_u, '/')) + sum(sweep(z$V^2, 2, z$sig2_v, '/'))
		expect_lte(penalty(new), penalty(native) + 1e-12 * max(1, penalty(native)))
	}
})
