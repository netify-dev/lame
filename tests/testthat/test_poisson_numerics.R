skip_on_cran()

test_that("Poisson log masses preserve extreme rates and cancellation near the mode", {
	pm = lame:::.poisson_log_mass
	expect_equal(pm(c(0, 1, 2), -1000), c(0, -1000, -2000 - log(2)), tolerance = 1e-12)
	expect_equal(pm(c(0, 1), Inf), c(-Inf, -Inf))
	expect_equal(pm(c(0, 1), -Inf), c(0, -Inf))
	for (y in c(1, 7, 1e8, 1e12)) for (delta in c(-.1, -1e-8, 0, 1e-8, .1)) {
		eta = log(y) + delta
		expect_equal(pm(y, eta), dpois(y, exp(eta), log = TRUE), tolerance = 1e-8)
	}
	z = matrix(c(0, -1000, log(1e12), 0), 2);y = matrix(c(NA, 1, 1e12, NA), 2)
	expected = matrix(c(0, sum(pm(y[c(2, 3)], z[c(2, 3)])), sum(pm(y[c(2, 3)], z[c(2, 3)])), 0), 2)
	expect_equal(lame:::ldZgbme_cpp(z, y, z, 0, 1), expected, tolerance = 1e-10)
	expect_identical(lame:::ldZgbme_cpp(z, y, z, 0, 1), lame:::ldZgbme_pois_exposure_cpp(z, y, z, 0, 1, 0))
})

test_that("Poisson simulation preserves large rates and explicit two-mode diagonals", {
	set.seed(66431L);ez = matrix(log(2e6), 20, 20)
	y = simY_pois(ez, bipartite = TRUE)
	expect_true(all(is.finite(y)));expect_gt(min(y), 1.9e6)
	ordinary = simY_pois(ez);expect_true(all(is.na(diag(ordinary))))
	expect_gt(min(ordinary[is.finite(ordinary)]), 1.9e6)
	expect_equal(lame:::simY_pois_exposure(matrix(1000, 3, 3), 0, TRUE), matrix(0, 3, 3))
	expect_error(lame:::simY_pois_exposure(ez, -1, TRUE), 'finite and nonnegative')
	expect_error(simY_pois(matrix(1000, 3, 3), bipartite = TRUE), 'overflowed')
	tiny = lame:::simY_pois_exposure(matrix(720, 3, 3), 1e-10, TRUE)
	expect_true(all(is.finite(tiny)));expect_gt(min(tiny), 1e300)
})

test_that("count dispersion uses the uncapped fitted mean", {
	for (mu in c(20, 1e12)) {
		y = array(rep(c(0, 2 * mu), 50), c(10, 10, 1));eta = list(matrix(log(mu), 10, 10))
		variance = lame:::.ae_pois_latent_var(y, eta, array(TRUE, dim(y)), 0)
		expect_equal(variance, log(2 - 1 / mu), tolerance = 1e-9)
	}
	expect_equal(lame:::.poisson_inverse_log1p(log1p(2e6)), 2e6, tolerance = 1e-10)
})

test_that("latent count updates move and preserve the extreme-tail Gaussian reference", {
	set.seed(66432L);n = 4096L;mu = -1000;ez = matrix(mu, 2, 2);y = matrix(c(NA, 1, 1, NA), 2)
	before = matrix(rnorm(2 * n, mu + 1), n, 2);after = matrix(0, n, 2)
	for (k in 1:n) {
		z = ez;z[c(2, 3)] = before[k, ];z = lame:::rZ_pois_fc_cpp(z, ez, 0, 1, y);after[k, ] = z[c(2, 3)]
	}
	delta = after - before
	expect_gt(mean(delta != 0), .1)
	expect_lt(abs(mean(delta)), 5 * sd(as.vector(delta)) / sqrt(length(delta)))
	for (mode in c('symmetric', 'bipartite')) {
		initial = rnorm(n, mu + 1);updated = numeric(n)
		for (k in 1:n) {
			if (mode == 'symmetric') {z = matrix(initial[k], 2, 2);updated[k] = lame:::rZ_pois_sym_fc(z, ez, 1, y)[1, 2]} else updated[k] = lame:::rZ_pois_bip_fc(
				matrix(initial[k], 1, 1),
				matrix(mu, 1, 1),
				1,
				matrix(1, 1, 1)
			)[1, 1]
		}
		delta = updated - initial
		expect_gt(mean(delta != 0), .1)
		expect_lt(abs(mean(delta)), 5 * sd(delta) / sqrt(n))
	}
})

test_that("square bipartite observation helpers retain every eligible cell", {
	ez = matrix(1000, 3, 3)
	expect_true(all(simY_bin(ez, 0, bipartite = TRUE) == 1))
	expect_equal(simY_nrm(ez, 0, 0, bipartite = TRUE), ez)
	expect_true(all(is.finite(simY_ord(ez, 0, matrix(rep(1:3, 3), 3), bipartite = TRUE))))
	ranks = simY_frn(ez, 0, 3, bipartite = TRUE)
	expect_true(all(is.finite(ranks)));expect_true(all(rowSums(ranks > 0) == 3))
	expect_true(all(is.na(diag(simY_frn(ez, 0, 3)))))
})


test_that("public conditional simulation retains square bipartite observations", {
	n = 3L;S = 6L;nm = list(paste0('r', 1:n), paste0('c', 1:n))
	x = array(1, c(n, n, 1), dimnames = c(nm, list('intercept')))
	f = list(
		BETA = matrix(rep(0, S), ncol = 1, dimnames = list(NULL, 'intercept')),
		VC = cbind(va = rep(.1, S), vb = .1, rho = .8, ve = .01),
		R = 0L,
		APM = setNames(rep(0, n), nm[[1]]),
		BPM = setNames(rep(0, n), nm[[2]]),
		nA = n,
		nB = n,
		mode = 'bipartite',
		symmetric = FALSE,
		X = x,
		Y = matrix(rep(0:2, n), n, n, dimnames = nm),
		YPM = matrix(0, n, n, dimnames = nm)
	)
	for (long in c(FALSE, TRUE)) for (start in c('posterior', 'random')) for (family in c(
		'normal', 'binary', 'poisson', 'ordinal', 'frn', 'cbin'
	)) {
		fit = f;fit$family = family
		fit$BETA[] = if (family == 'poisson') log(2e6) else 0
		if (long) {fit$Xlist = rep(list(x), 2);fit$Y = rep(list(f$Y), 2);fit$YPM = rep(list(f$YPM), 2);fit$n_time = 2L;fit$dynamic_ab = fit$dynamic_uv = FALSE}
		class(fit) = if (long) c('lame', 'ame') else 'ame'
		set.seed(66433);saved = .Random.seed
		draws = if (long) simulate(
			fit,
			nsim = 4,
			seed = 101,
			type = 'approximate',
			start_from = start
		) else simulate(fit, nsim = 4, seed = 101, type = 'approximate')
		expect_identical(.Random.seed, saved)
		ys = if (long) unlist(draws$Y, recursive = FALSE) else draws$Y
		for (y in ys) {expect_true(all(is.finite(y)));if (family == 'poisson') expect_gt(max(y), 1.2e6)}
	}
})
