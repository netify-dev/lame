test_that("Poisson marginal log scores resolve broad latent variances", {
	for (v in c(4, 100, 1e4)) {
		for (y in c(0, 1, 3)) {
			eta = 4
			f = function(z) dpois(y, exp(eta + sqrt(v) * z)) * dnorm(z)
			oracle = log(sum(vapply(seq_len(4), function(j) {
				cuts = c(-Inf, -1, 0, 1, Inf)
				integrate(
					f,
					cuts[j],
					cuts[j + 1],
					rel.tol = 1e-11,
					abs.tol = 1e-12,
					subdivisions = 1000
				)$value
			}, 0.)))
			expect_lt(abs(lame:::.lame_pois_marginal_loglik(y, eta, v) - oracle), 1e-8)
		}
	}
})

test_that("large Poisson counts agree with an independent Gamma mixture identity", {
	for (y in c(1e4, 1e6, 1e12)) {
		eta = log(y) - .1; v = 4
		# p(y) = E[N(log(G); eta, v)]/y for G ~ Gamma(shape=y, rate=1).
		oracle = log(
			integrate(
				function(u) dnorm(log(qgamma(u, y)), eta, sqrt(v)),
				0,
				1,
				rel.tol = 1e-11,
				abs.tol = 1e-12,
				subdivisions = 1000
			)$value
		) - log(y)
		expect_lt(abs(lame:::.lame_pois_marginal_loglik(y, eta, v) - oracle), 1e-8)
	}
})

test_that("Poisson integration is deterministic and keeps positive tiny variance", {
	set.seed(473285)
	rng = .Random.seed
	for (v in c(1e-20, 1e-12)) {
		ans = lame:::.lame_pois_marginal_loglik(c(0, 1, 2), rep(log(2), 3), v)
		oracle = vapply(0:2, function(y) log(integrate(function(z) {
			dpois(y, 2 * exp(sqrt(v) * z)) * dnorm(z)
		}, -10, 10, rel.tol = 1e-12, abs.tol = 1e-13)$value), 0.)
		expect_equal(ans, oracle, tolerance = 1e-12)
	}
	expect_identical(.Random.seed, rng)
	expect_error(lame:::.lame_pois_marginal_loglik(1, 0, 0), "positive")
	expect_error(lame:::.lame_pois_marginal_loglik(.5, 0, 1), "integer")
})
