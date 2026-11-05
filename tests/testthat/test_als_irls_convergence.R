skip_on_cran()

test_that("IRLS score certificates agree with independent likelihood gradients", {
	ns = asNamespace('lame');set.seed(952819L)
	for (mode in c(
		'directed', 'bipartite', 'symmetric'
	)) for (kind in c('probit', 'logit', 'poisson')) for (rank in 0:2) {
		nr = 5L;nc = if (mode == 'bipartite') 4L else nr;nt = 2L;sym = mode == 'symmetric';family = if (kind == 'poisson') 'poisson' else 'binary';link = if (kind == 'poisson') 'log' else kind
		Y = X = vector('list', nt)
		for (tt in 1:nt) {
			yy = matrix(if (family == 'binary') rbinom(nr * nc, 1, .5) else rpois(nr * nc, 2), nr)
			xx = matrix(rnorm(nr * nc), nr)
			if (sym) {yy[lower.tri(yy)] = t(yy)[lower.tri(yy)];xx = (xx + t(xx)) / 2}
			if (mode != 'bipartite') diag(yy) = NA_real_
			yy[1, 2] = NA_real_;if (sym) yy[2, 1] = NA_real_
			dimnames(yy) = list(
				paste0('r', 1:nr),
				if (mode == 'bipartite') paste0('c', 1:nc) else paste0('r', 1:nc)
			);Y[[tt]] = yy;X[[tt]] = array(xx, c(nr, nc, 1), dimnames = c(dimnames(yy), list('x')))
		}
		prep = get('.ae_prepare', ns)(Y, X, NULL, NULL, if (mode == 'bipartite') 'bipartite' else 'unipartite', TRUE, 'quiet')
		f = list(
			mu = -.2,
			beta = .4,
			a = rnorm(nr, sd = .2),
			b = rnorm(nc, sd = .2),
			U = matrix(rnorm(nr * rank, sd = .4), nr, rank),
			V = matrix(rnorm(nc * rank, sd = .4), nc, rank),
			sig2_u = seq(.6, 1.2, length.out = rank),
			sig2_v = seq(.8, 1.4, length.out = rank)
		)
		if (sym) {f$b = f$a;f$U = if (rank) qr.Q(
			qr(matrix(rnorm(nr * rank), nr))
		) else matrix(0, nr, 0);f$L = seq(.5, -.3, length.out = rank);f$sig2_u = f$sig2_v = rep(.7, rank);f$V = f$U;f$Omat = if (rank) f$U %*% (f$L * t(f$U)) else matrix(0, nr, nc)} else {f$Omat = tcrossprod(f$U, f$V)}
		certificate = get('.ae_irls_stationarity', ns)(prep, f, family, link, rank, sym)
		if (sym) {sg = sign(f$L);B = if (rank) sweep(f$U, 2, sqrt(abs(f$L)), '*') else matrix(0, nr, 0);par = c(f$mu, f$beta, f$a, c(B))} else {par = c(
			f$mu,
			f$beta,
			f$a,
			f$b,
			c(f$U),
			c(f$V)
		)}
		objective = function(z) {
			mu = z[1];beta = z[2];a = z[seq_len(nr) + 2L];off = 2L + nr
			if (sym) {b = a;B = matrix(
				z[seq_len(nr * rank) + off],
				nr,
				rank
			);O = if (rank) B %*% (sg * t(B)) else matrix(0, nr, nc);pen = if (rank) -4 * sum(
				dnorm(B, sd = sqrt(f$sig2_u[1]), log = TRUE)
			) else 0} else {
				b = z[seq_len(nc) + off];off = off + nc;U = matrix(
					z[seq_len(nr * rank) + off],
					nr,
					rank
				);off = off + nr * rank;V = matrix(z[seq_len(nc * rank) + off], nc, rank);O = tcrossprod(U, V)
				pen = if (rank) -2 * sum(
					vapply(
						seq_len(rank),
						function(k) sum(
							dnorm(U[, k], sd = sqrt(f$sig2_u[k]), log = TRUE)
						) + sum(dnorm(V[, k], sd = sqrt(f$sig2_v[k]), log = TRUE)),
						numeric(1)
					)
				) else 0
			}
			value = pen
			for (tt in 1:nt) {eta = mu + beta * prep$X[, , 1, tt] + outer(a, b, '+') + O;keep = is.finite(prep$Yarr[, , tt]);yy = prep$Yarr[, , tt][keep]
				ll = if (family == 'poisson') dpois(yy, exp(eta[keep]), log = TRUE) else dbinom(
					yy,
					1,
					if (kind == 'probit') pnorm(eta[keep]) else plogis(eta[keep]),
					log = TRUE
				)
				value = value - 2 * sum(ll)
			}
			value
		}
		step = 1e-5;numeric_grad = vapply(
			seq_along(par),
			function(k) {hi = lo = par;hi[k] = hi[k] + step;lo[k] = lo[k] - step;(objective(hi) - objective(lo)) / (2 * step)},
			numeric(1)
		)
		error = max(abs(numeric_grad - certificate$gradient));expect_lt(error, 1e-5)

	}
})

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

test_that("Poisson IRLS reaches the independent fixed-effects GLM solution", {
	set.seed(995217L)
	nr = 24L;nc = 20L;nt = 2L;a = rnorm(nr, sd = .4);b = rnorm(nc, sd = .4)
	x = lapply(seq_len(nt), function(t) matrix(rnorm(nr * nc), nr))
	y = lapply(x, function(xx) matrix(rpois(nr * nc, exp(.2 + 1.1 * xx + outer(a, b, '+'))), nr))
	for (t in seq_len(nt)) y[[t]][1, 2] = NA_real_
	xx = lapply(x, function(z) array(z, c(nr, nc, 1)))
	f = lame_als(y, Xdyad = xx, mode = 'bipartite', family = 'poisson', R = 0, bias_correct = 'none', verbose = FALSE)
	df = do.call(
		rbind,
		lapply(
			seq_len(nt),
			function(t) data.frame(
				y = c(y[[t]]),
				x = c(x[[t]]),
				sender = factor(c(row(y[[t]]))),
				receiver = factor(c(col(y[[t]])))
			)
		)
	)
	oracle = glm(
		y ~ x + sender + receiver,
		data = df,
		family = poisson(),
		control = glm.control(epsilon = 1e-12, maxit = 100)
	)
	expect_true(f$converged)
	expect_true(f$irls_diagnostics$pass)
	expect_equal(unname(f$beta[1]), unname(coef(oracle)['x']), tolerance = 1e-5)
	expect_equal(f$deviance, deviance(oracle), tolerance = 1e-7)
	expected_weights = array(unlist(f$fitted), dim(f$obs_weights));keep = is.finite(f$obs_weights)
	expect_equal(unname(f$obs_weights[keep]), unname(expected_weights[keep]), tolerance = 1e-8)
	expect_true(all(is.na(f$obs_weights[!is.finite(f$Y)])))
	expect_true(all(f$irls_diagnostics$conditional_gain >= -1e-8))
	err = tryCatch(
		lame_als(
			y,
			Xdyad = xx,
			mode = 'bipartite',
			family = 'poisson',
			R = 0,
			bias_correct = 'none',
			max_iter = 1,
			tol = 1e-14,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, 'lame_als_convergence')
	expect_false(err$state$pass)
	expect_identical(err$kind, 'irls')
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

test_that("large count residuals retain the exact Poisson mean target", {
	n = 60L;y = matrix(1, n, n);diag(y) = NA_real_
	y[cbind(seq_len(n), c(2:n, 1L))] = 10000
	f = ame_als(y, R = 0, family = 'poisson', bias_correct = 'none', tol = 1e-10, verbose = FALSE)
	expect_true(f$converged)
	expect_equal(mean(unlist(f$fitted), na.rm = TRUE), mean(y, na.rm = TRUE), tolerance = 1e-9)
	expect_true(f$irls_diagnostics$pass)
})

test_that("Poisson fitted means and predictions can exceed one million", {
	y = matrix(2e6, 10, 10);diag(y) = NA_real_
	f = ame_als(y, R = 0, family = 'poisson', bias_correct = 'none', verbose = FALSE)
	expect_true(f$converged)
	expect_equal(mean(unlist(f$fitted), na.rm = TRUE), 2e6, tolerance = 1e-9)
	expect_equal(mean(unlist(predict(f, type = 'response')), na.rm = TRUE), 2e6, tolerance = 1e-9)
	expect_lt(abs(f$deviance), 1e-6)
	sim = unlist(simulate(f, nsim = 20, seed = 927114L)$Y);sim = sim[is.finite(sim)]
	expect_lt(abs(mean(sim) - 2e6) / sqrt(2e6 / length(sim)), 6)
	expect_error(get('.als_pois_rate', asNamespace('lame'))(1000), class = 'lame_als_prediction')
})
