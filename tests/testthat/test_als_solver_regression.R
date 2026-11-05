test_that("balanced cross-section partitions repair an avoidable boundary half", {
	skip_on_cran()
	set.seed(62977L);n = 120L;y = matrix(rbinom(n * n, 1, .5), n);diag(y) = NA_real_
	y[, 2] = 0;y[2, 2] = NA_real_;y[c(1, 3), 2] = 1
	yy = array(y, c(n, n, 1));rng = .Random.seed
	sp = lame:::.ae_jk_split(yy, 'binary', FALSE, FALSE, R = 0L)
	expect_null(sp$problem);expect_gt(sp$checkerboard_attempt, 1L)
	expect_lte(sp$checkerboard_attempt, 32L);expect_identical(.Random.seed, rng)
	for (mask in sp$hide) {
		half = yy;half[mask] = NA_real_;state = lame:::.ae_jk_actor_state(half, 'binary', FALSE)
		expect_false(any(state$row$bound));expect_false(any(state$col$bound))
		expect_equal(mask[, , 1], t(mask[, , 1]));expect_lte(diff(range(rowSums(mask[, , 1]))), 1)
	}
	expect_identical(sp, lame:::.ae_jk_split(yy, 'binary', FALSE, FALSE, R = 0L))
	y[3, 2] = 0
	refused = lame:::.ae_jk_split(array(y, c(n, n, 1)), 'binary', FALSE, FALSE, R = 0L)
	expect_match(refused$problem, 'all 32 balanced checkerboard partitions fail')
	expect_true(all(refused$tried$no_tie > 0))
})

test_that("factor variance ordering cannot falsely stop the outer likelihood", {
	skip_on_cran()
	for (case in list(
		list(group = 8L, rep = 13L, n = 90L, r = 8L, sd = .4, link = 'probit'),
		list(group = 3L, rep = 24L, n = 60L, r = 2L, sd = .6, link = 'logit')
	)) {
		set.seed(1120000L + 1000L * case$group + case$rep);n = case$n;r = case$r
		a = rnorm(n, sd = .5);b = rnorm(n, sd = .5)
		u = matrix(rnorm(n * r, sd = case$sd), n, r);v = matrix(rnorm(n * r, sd = case$sd), n, r)
		Y = X = vector('list', 6L)
		for (t in 1:6) {
			xx = matrix(rnorm(n * n), n);eta = -.5 + .5 * xx + outer(a, b, '+') + tcrossprod(u, v)
			threshold = if (case$link == 'probit') eta else qnorm(plogis(eta, log.p = TRUE), log.p = TRUE)
			g = matrix(rnorm(n * n), n);noise = (sqrt(.7) * (g + t(g)) + sqrt(.3) * (g - t(g))) / sqrt(2)
			yy = 1 * (threshold + noise > 0);diag(yy) = NA_real_;dimnames(yy) = list(paste0('r', 1:n), paste0('r', 1:n))
			Y[[t]] = yy;X[[t]] = array(xx, c(n, n, 1), dimnames = c(dimnames(yy), list('x')))
		}
		fit = lame_als(
			Y,
			Xdyad = X,
			R = r,
			family = 'binary',
			link = case$link,
			bias_correct = 'jackknife',
			verbose = FALSE
		)
		expect_true(fit$converged)
		certificates = c(list(fit$irls_diagnostics), fit$bias_correct$irls_diagnostics)
		expect_length(certificates, 3L)
		expect_true(all(vapply(certificates, function(z) isTRUE(z$pass), logical(1))))
	}
})
