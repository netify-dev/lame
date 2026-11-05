test_that("innovation variance counts coordinates on each side of a bipartite model", {
	set.seed(782361)
	for (ranks in list(c(1L, 2L), c(2L, 1L), c(2L, 3L), c(2L, 2L))) {
		U = array(rnorm(3 * ranks[1] * 4, sd = .7), c(3L, ranks[1], 4L))
		V = array(rnorm(4 * ranks[2] * 4, sd = .7), c(4L, ranks[2], 4L))
		rho = .4
		for (weighted in FALSE) {
			wu = matrix(seq(.5, 1.5, length.out = 12), 3, 4)
			wv = matrix(seq(.4, 1.8, length.out = 16), 4, 4)
			wu[1, 2] = wv[2, 3] = 0
			if (!weighted) { wu[wu > 0] = 1; wv[wv > 0] = 1 }
			# an independent vector of whitened Gaussian coordinates
			whiten = function(z, weight) {
				ans = sqrt(1 - rho^2) * as.vector(z[, , 1])
				for (tt in 2:4) {
					keep = which(weight[, tt] > 0)
					current = matrix(z[, , tt], dim(z)[1], dim(z)[2])[keep, , drop = FALSE]
					previous = matrix(z[, , tt - 1L], dim(z)[1], dim(z)[2])[keep, , drop = FALSE]
					ans = c(ans, as.vector(sweep(current - rho * previous, 1, sqrt(weight[keep, tt]), '*')))
				}
				ans
			}
			coordinates = c(whiten(U, wu), whiten(V, wv))
			shape = 2 + length(coordinates) / 2
			rate = 1 + sum(coordinates^2) / 2
			set.seed(53917)
			reference = sqrt(1 / rgamma(512, shape = shape, rate = rate))
			set.seed(53917)
			actual = replicate(
				512,
				.lame_sample_sigma_uv_drift_only(
					U,
					V,
					rho,
					.7,
					FALSE,
					delta_u = 1 - wu,
					delta_v = 1 - wv,
					prior_shape = 2,
					prior_scale = 1
				)
			)
			expect_equal(actual, reference, tolerance = 1e-12, info = paste(ranks, weighted, collapse = '/'))
		}
	}
})
