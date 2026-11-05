.als_prediction_fixture = function(link = 'probit', mode = 'unipartite', rho = .6, n = 20L) {
	f = list(
		EZ = list(matrix(-1, n, n)),
		dims = list(n_row = n, n_col = n),
		n_time = 1L,
		family = 'binary',
		VC = c(rho = rho, ve = 1),
		symmetric = FALSE,
		mode = mode,
		link = link,
		row_names = as.character(seq_len(n)),
		col_names = as.character(seq_len(n))
	)
	class(f) = 'ame_als';f
}

test_that("binary simulations preserve marginal links and reciprocal probabilities", {
	skip_on_cran()
	n = 20L;nsim = 1200L
	for (mode in c('unipartite', 'bipartite')) for (link in c('probit', 'logit')) {
		rho = if (mode == 'bipartite') 0 else .6;f = .als_prediction_fixture(link, mode, rho, n)
		p = if (link == 'probit') pnorm(-1) else plogis(-1);q = qnorm(p)
		joint = if (rho == 0) p * p else integrate(
			function(z) dnorm(z) * pnorm((q - rho * z) / sqrt(1 - rho^2)),
			-Inf,
			q,
			rel.tol = 1e-12
		)$value
		set.seed(998L);before = .Random.seed;sims = simulate(f, nsim = nsim, seed = 918372L)$Y
		expect_identical(.Random.seed, before)
		average = vapply(sims, function(y) mean(y, na.rm = TRUE), numeric(1))
		pair = vapply(sims, function(y) mean((y * t(y))[upper.tri(y)], na.rm = TRUE), numeric(1))
		expect_lt(abs(mean(average) - p) / (sd(average) / sqrt(nsim)), 5)
		expect_lt(abs(mean(pair) - joint) / (sd(pair) / sqrt(nsim)), 5)
		if (mode == 'bipartite') expect_true(all(vapply(sims, function(y) all(is.finite(diag(y))), logical(1))))
	}
})

test_that("binary pair likelihood agrees with the analytic zero-threshold estimate", {
	estimate = get('.ae_binary_pair_rho', asNamespace('lame'))
	n = 10L;pairs = which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)[1:20, , drop = FALSE]
	for (same in c(0L, 2L, 10L, 18L, 20L)) {
		first = c(rep(0, same / 2), rep(1, same / 2), rep(0, (20 - same) / 2), rep(1, (20 - same) / 2))
		second = c(rep(0, same / 2), rep(1, same / 2), rep(1, (20 - same) / 2), rep(0, (20 - same) / 2))
		y = matrix(NA_real_, n, n);y[pairs] = first;y[pairs[, 2:1]] = second
		yy = array(y, c(n, n, 1));eta = array(0, dim(yy))
		reference = sin(pi * (same / 20 - .5))
		for (link in c('probit', 'logit')) {
			z = estimate(yy, eta, link)
			expect_equal(z$rho, reference, tolerance = 1e-6)
			expect_identical(z$identified, TRUE)
			expect_equal(z$n_pairs, 20)
			expect_identical(z$boundary, same %in% c(0L, 20L))
		}
		yy[9, 10, 1] = 1
		expect_equal(estimate(yy, eta, 'probit')$rho, reference, tolerance = 1e-6)
	}
})

test_that("endpoint interval probabilities agree with independent quadrature", {
	lo = c(-12, -2, -.2, 1, 11);hi = c(-11, -1, .3, 2, 12)
	reference = vapply(
		seq_along(lo),
		function(k) log(integrate(dnorm, lo[k], hi[k], rel.tol = 1e-12, abs.tol = 1e-100)$value),
		numeric(1)
	)
	expect_equal(get('.ae_normal_interval_log', asNamespace('lame'))(lo, hi), reference, tolerance = 1e-10)
	expect_identical(get('.ae_normal_interval_log', asNamespace('lame'))(c(1, 2), c(1, 1)), c(-Inf, -Inf))
})

test_that("unidentified binary dependence is explicit and perfect dependence is retained", {
	y = array(NA_real_, c(3, 3, 1));y[1, 2, 1] = 1
	z = get('.ae_binary_pair_rho', asNamespace('lame'))(y, array(0, dim(y)), 'probit')
	expect_false(z$identified);expect_true(is.na(z$rho));expect_equal(z$n_pairs, 0)
	f = .als_prediction_fixture(n = 2L);f$rho_estimation = z;f$VC['rho'] = NA_real_
	expect_error(simulate(f, seed = 91), class = 'lame_als_prediction')
	f$rho_estimation = list(identified = TRUE);f$EZ = list(matrix(0, 2, 2))
	for (rho in c(-1, 1)) {
		f$VC['rho'] = rho;sims = simulate(f, nsim = 200, seed = 5121)$Y
		one = vapply(sims, function(y) y[1, 2], numeric(1));two = vapply(sims, function(y) y[2, 1], numeric(1))
		expect_true(if (rho == 1) all(one == two) else all(one + two == 1))
	}
})

test_that("the corrected fit estimates reciprocity at its reported predictor", {
	set.seed(781148L);n = 30L
	y = lapply(1:4, function(t) {z = matrix(rbinom(n * n, 1, .5), n);diag(z) = NA_real_;z})
	f = lame_als(y, R = 0, family = 'binary', dcor = FALSE, verbose = FALSE)
	expect_true(f$converged)
	reference = get('.ae_binary_pair_rho', asNamespace('lame'))(f$Y, array(unlist(f$EZ), dim(f$Y)), f$link)
	expect_equal(unname(f$VC['rho']), reference$rho, tolerance = 1e-12)
	expect_equal(f$rho_estimation, reference, tolerance = 1e-12)
	expect_true(is.finite(f$rho_working))
	bad = f$Y;bad[1, 2, 1] = 2
	expect_error(ame_als_refit(f, Y_new = bad), class = 'lame_als_data')
	f$VC['rho'] = 1;expect_equal(als_start_vals(f)$rho, .99)
	f$VC['rho'] = -1;expect_equal(als_start_vals(f)$rho, -.99)
})

test_that("prior summaries describe IRLS factor regularization", {
	f = .als_prediction_fixture();f$R = 2L;f$non_normal_method = 'irls'
	output = paste(capture.output(prior_summary(f)), collapse = ' ')
	expect_match(output, 'regularized latent factors')
	expect_false(grepl('There are no priors', output, fixed = TRUE))
})
