test_that("stalled directed count paths reach their unchanged score criterion", {
	skip_on_cran()
	set.seed(1984001);n = 60L;nt = 8L;rank = 1L
	paths = function(n, r, sd, rho) {x = array(0, c(n, r, nt));x[, , 1] = rnorm(
		n * r,
		0,
		sd / sqrt(1 - rho^2)
	);for (t in 2:nt) x[, , t] = rho * x[, , t - 1] + rnorm(n * r, 0, sd);x}
	a = matrix(rep(rnorm(n, 0, .3), nt), n) + matrix(paths(n, 1, .15, .7), n, nt)
	b = matrix(rep(rnorm(n, 0, .3), nt), n) + matrix(paths(n, 1, .15, .7), n, nt)
	u = paths(n, rank, .3, .7);v = paths(n, rank, .3, .7);Y = X = vector('list', nt)
	for (t in seq_len(nt)) {
		x = matrix(rnorm(n * n), n);eta = .5 + .5 * x + outer(a[, t], b[, t], '+') + tcrossprod(
			matrix(u[, , t], n, rank),
			matrix(v[, , t], n, rank)
		)
		g = matrix(rnorm(n * n), n);z = (sqrt(.7) * (g + t(g)) + sqrt(.3) * (g - t(g))) / sqrt(2)
		y = matrix(rpois(n * n, exp(eta + sqrt(.4) * z)), n);diag(y) = NA_real_
		dimnames(y) = list(
			paste0('a', seq_len(n)),
			paste0('a', seq_len(n))
		);Y[[t]] = y;X[[t]] = array(x, c(n, n, 1), dimnames = c(dimnames(y), list('x')))
	}
	fit = lame(
		Y,
		Xdyad = X,
		R = rank,
		family = 'poisson',
		method = 'als',
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		als_max_iter = 200L,
		als_tol = 1e-6,
		verbose = FALSE
	)
	expect_true(fit$converged)
	expect_true(fit$diagnostics$stationarity$pass)
	if (!is.null(fit$diagnostics$joint_refinement)) {
		expect_true(fit$diagnostics$joint_refinement$accepted)
		expect_lte(
			sum(vapply(fit$diagnostics$joint_refinement$cycles, `[[`, integer(1), 'budget')),
			fit$diagnostics$joint_refinement$max_iter
		)
		expect_lte(fit$diagnostics$joint_refinement$objective_after, fit$diagnostics$joint_refinement$objective_before)
		expect_equal(fit$objective, fit$diagnostics$joint_refinement$objective_after, tolerance = 1e-10)
	}
	expect_equal(tail(fit$dev_history, 1), fit$objective, tolerance = 1e-10)
	expect_true(is.finite(fit$VC['ve']) && fit$VC['ve'] > 0 && fit$VC['ve'] < 1)
})
