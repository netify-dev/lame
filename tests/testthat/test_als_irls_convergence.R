skip_on_cran()

test_that("factor balancing preserves predictions and minimizes the gauge penalty", {
	set.seed(983116L);records = list()
	for (rank in c(1L, 2L, 4L)) for (case in 1:25) {
		nr = 9L;nc = 7L;u = matrix(rnorm(nr * rank), nr);v = matrix(rnorm(nc * rank), nc);su = exp(rnorm(rank));sv = exp(rnorm(rank))
		if (case %% 5 == 0) {u[, 1] = 0;v[, 1] = 0};if (rank > 1 && case %% 7 == 0) {u[, 2] = u[, 1];v[, 2] = v[, 1]}
		z = get(
			'.ae_factor_balance',
			asNamespace('lame')
		)(u, v, su, sv);before = sum(colSums(u^2) / su) + sum(colSums(v^2) / sv);after = sum(colSums(z$U^2) / su) + sum(colSums(z$V^2) / sv)
		product_error = max(
			abs(tcrossprod(u, v) - tcrossprod(z$U, z$V))
		);expect_lt(product_error, 1e-11);expect_lte(after, before + 1e-9)
		bound = 2 * sum(
			svd(tcrossprod(u, v), nu = 0, nv = 0)$d[1:rank] / sort(sqrt(su * sv), decreasing = TRUE)
		);expect_equal(after, bound, tolerance = 1e-9)
		gauge_obj = function(h) {h = matrix(h, rank); uu = z$U %*% h;vv = z$V %*% t(solve(h));sum(colSums(uu^2) / su) + sum(colSums(vv^2) / sv)}
		ident = c(diag(rank));eps = 1e-5;grad = vapply(
			seq_along(ident),
			function(k) {a = b = ident;a[k] = a[k] + eps;b[k] = b[k] - eps;(gauge_obj(a) - gauge_obj(b)) / (2 * eps)},
			numeric(1)
		)
		expect_lt(max(abs(grad)), 1e-5)
		for (j in 1:20) {h = diag(rank) + matrix(rnorm(rank^2, sd = .1), rank);expect_gte(gauge_obj(h), after - 1e-8)}
		records[[length(records) + 1L]] = data.frame(
			rank = rank,
			case = case,
			product_error = product_error,
			gradient_error = max(abs(grad)),
			penalty_gain = before - after
		)
	}
})

test_that("a Poisson fit stopped after one iteration is refused", {
	set.seed(995217L)
	nr = 24L;nc = 20L;nt = 2L;a = rnorm(nr, sd = .4);b = rnorm(nc, sd = .4)
	x = lapply(seq_len(nt), function(t) matrix(rnorm(nr * nc), nr))
	y = lapply(x, function(xx) matrix(rpois(nr * nc, exp(.2 + 1.1 * xx + outer(a, b, '+'))), nr))
	for (t in seq_len(nt)) y[[t]][1, 2] = NA_real_
	xx = lapply(x, function(z) array(z, c(nr, nc, 1)))
	err = tryCatch(
		lame_als(
			y,
			Xdyad = xx,
			mode = 'bipartite',
			family = 'poisson',
			R = 0,
			max_iter = 1,
			tol = 1e-14,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, 'lame_als_convergence')
})

test_that("a symmetric zero product retains its updated variance scale", {
	ns = asNamespace('lame');n = 8L;rank = 2L;zero = matrix(0, n, n)
	weights = matrix(1, n, n);diag(weights) = 0
	z = get('.ae_lowrank_update_pen', ns)(zero, matrix(1, n, n), weights, weights > 0,
		rank, TRUE, matrix(0, n, rank), matrix(0, n, rank), rep(.2, rank), rep(.2, rank), .1, 6)
	expect_equal(z$sig2_u, rep(.6 / (6 + 2 * n * rank), rank), tolerance = 1e-12)
	expect_equal(z$Omat, zero, tolerance = 1e-12)
})

test_that("capped IRLS increments retain the family likelihood score", {
	ns = asNamespace('lame')
	eta = c(-30, -4, 0, 6, 25);y = c(0, 100, 1, 1e5, 1e12)
	z = get('.ae_irls_poisson', ns)(y, eta);score = y - exp(eta)
	allowance = 32 * .Machine$double.eps * pmax(1, abs(eta)) * z$w + 1e-12 * pmax(1, abs(score))
	expect_true(all(abs(z$w * (z$z - eta) - score) <= allowance))
	eta = rep(c(-40, -10, -3, 0, 3, 10, 40), 2);y = rep(0:1, each = 7)
	for (link in c('probit', 'logit')) {
		work = get(paste0('.ae_irls_binary_', link), ns)(y, eta)
		loglik = function(e) if (link == 'probit') pnorm((2 * y - 1) * e, log.p = TRUE) else -log1p(exp(-abs(e))) - pmax(-(2 * y - 1) * e, 0)
		h = 1e-5;reference = (loglik(eta + h) - loglik(eta - h)) / (2 * h)
		expect_equal(work$w * (work$z - eta), reference, tolerance = 1e-7)
		expect_equal(get('.ae_irls_deviance', ns)(y, eta, 'binary', link), -2 * sum(loglik(eta)), tolerance = 1e-12)
	}
})

test_that("Poisson fitted means and predictions can exceed one million", {
	y = matrix(2e6, 10, 10);diag(y) = NA_real_
	f = ame_als(y, R = 0, family = 'poisson', verbose = FALSE)
	expect_true(f$converged)
	sim = unlist(simulate(f, nsim = 20, seed = 927114L)$Y);sim = sim[is.finite(sim)]
	expect_lt(abs(mean(sim) - 2e6) / sqrt(2e6 / length(sim)), 6)
	expect_error(get('.als_pois_rate', asNamespace('lame'))(1000), class = 'lame_als_prediction')
})
