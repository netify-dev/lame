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
			verbose = FALSE
		)
		expect_true(fit$converged)
		expect_true(isTRUE(fit$irls_diagnostics$pass))
	}
})
