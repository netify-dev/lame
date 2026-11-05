test_that("proper high-rank factor priors preserve independent joint draws", {
	skip_on_cran()
	n = 5L;nt = 2L;s2 = .8;rho = .45
	# Geweke prior-predictive check of the joint law of factors, data and covariance for both implementations,
	# under the default prior of R/uv_prior.R and two other (df, psi); bounded/log statistics avoid undefined IW moments
	checks = list();N = 1200L
	for (rank in 3:4) for (pr in c('default', 'previous', 'tight')) for (backend in c('compiled', 'R')) {
		df = switch(
			pr,
			default = 2 * rank + 2,
			previous = max(rank + 2, 2 * rank),
			tight = 2 * rank + 6
		);psi = switch(pr, default = df, previous = 1, tight = 0.5)
		set.seed(69351 + rank + 10 * match(pr, c('default', 'previous', 'tight')));delta = matrix(0, N, 6)
		for (k in seq_len(N)) {
			S = solve(
				rWishart(1, df, diag(2 * rank) / psi)[, , 1]
			);W = matrix(
				rnorm(n * 2 * rank),
				n,
				2 * rank
			) %*% chol(S);U = W[, seq_len(rank), drop = FALSE];V = W[, rank + seq_len(rank), drop = FALSE];eta = tcrossprod(U, V);E = array(0, c(n, n, nt))
			for (t in seq_len(nt)) {z = matrix(0, n, n);diag(z) = rnorm(n, sd = sqrt(s2 * (1 + rho)));for (i in 1:(n - 1)) for (j in (i + 1):n) {a = rnorm(2);z[i, j] = sqrt(s2) * a[1];z[j, i] = sqrt(s2) * (rho * a[1] + sqrt(1 - rho^2) * a[2])};E[, , t] = eta + z}
			stat = function(U, V) {p = tcrossprod(U, V);c(
				log1p(sum(U^2) + sum(V^2)),
				atan(p[1, 2]),
				atan(p[1, 2]) * atan(E[1, 2, 1]),
				mean(atan(p)),
				mean(atan(p) * atan(E[, , 1])),
				as.numeric(abs(U[1, 1]) < 1)
			)}
			z = if (backend == 'compiled') lame:::rUV_rep_fc_cpp(
				E,
				U,
				V,
				rho,
				s2,
				diag(2),
				1e-6,
				TRUE,
				0:(rank - 1),
				df,
				psi
			) else rUV_rep_fc(E, U, V, rho, s2, TRUE, df, psi)
			delta[k, ] = stat(z$U, z$V) - stat(U, V)
		}
		z = colMeans(delta) / (apply(delta, 2, sd) / sqrt(N));expect_true(
			all(is.finite(z)),
			info = paste(rank, pr, backend)
		);expect_true(max(abs(z)) < 6, info = paste(rank, pr, backend));checks[[paste(rank, pr, backend)]] = z
	}

})

test_that("large prior factor draws retain a finite precision update", {
	x = readRDS(test_path("fixtures", "static_factor_precision_extreme.rds"))
	# captured from a proper rank-four prior draw on which the I + W'W inversion breaks down
	set.seed(535315)
	original = x
	ih = lame:::mhalf(solve(x$s2 * matrix(c(1, x$rho, x$rho, 1), 2)))
	for (k in 1:10) {
		a = rUV_rep_fc(x$E, x$U, x$V, x$rho, x$s2, TRUE)
		b = lame:::rUV_rep_fc_cpp(
			x$E,
			x$U,
			x$V,
			x$rho,
			x$s2,
			ih,
			1e-6,
			TRUE,
			0:3,
			lame:::.lame_uv_prior_df(4L),
			lame:::.lame_uv_prior_df(4L)
		)
		expect_true(all(is.finite(c(a$U, a$V, b$U, b$V))))
	}
	expect_identical(x, original)
})

test_that("the cross-sectional factor update (rSuv_fc then rUV_fc) preserves the joint law under the shared prior", {
	skip_on_cran()
	# Geweke check of the ame() kernel: covariance from its inverse-Wishart conditional, then rows given the covariance;
	# one sweep must leave the prior-predictive law of (Suv, U, V) unchanged
	n = 6L;s2 = .8;rho = .45;N = 1500L
	for (rank in 1:3) for (pr in c('default', 'other')) {
		df = if (pr == 'default') lame:::.lame_uv_prior_df(rank) else 2 * rank + 5;s = if (pr == 'default') 1 else 0.3
		set.seed(7115 + rank + (pr == 'other') * 10);delta = matrix(0, N, 5)
		for (k in seq_len(N)) {
			S = solve(
				rWishart(1, df, diag(2 * rank) / (df * s))[, , 1]
			);W = matrix(
				rnorm(n * 2 * rank),
				n,
				2 * rank
			) %*% chol(S);U = W[, seq_len(rank), drop = FALSE];V = W[, rank + seq_len(rank), drop = FALSE]
			z = matrix(0, n, n);diag(z) = rnorm(n, sd = sqrt(s2 * (1 + rho)));for (i in 1:(n - 1)) for (j in (i + 1):n) {a = rnorm(2);z[i, j] = sqrt(s2) * a[1];z[j, i] = sqrt(s2) * (rho * a[1] + sqrt(1 - rho^2) * a[2])}
			E = tcrossprod(U, V) + z
			stat = function(U, V, S) {p = tcrossprod(U, V);c(
				log(sum(diag(S))),
				log1p(sum(U^2) + sum(V^2)),
				atan(p[1, 2]),
				atan(p[1, 2]) * atan(E[1, 2]),
				mean(atan(p) * atan(E))
			)}
			S2 = rSuv_fc(U, V, Suv0 = diag(2 * rank) * s, kappa0 = df);z2 = rUV_fc(E, U, V, S2, rho, s2, offset = 0)
			delta[k, ] = stat(z2$U, z2[['V']], S2) - stat(U, V, S)
		}
		zz = colMeans(delta) / (apply(delta, 2, sd) / sqrt(N));expect_true(
			all(is.finite(zz)),
			info = paste(rank, pr)
		);expect_true(max(abs(zz)) < 6, info = paste(rank, pr))
	}
})
