skip_on_cran()

test_that('correlated additive priors retain their near-boundary conditional', {
	S = matrix(c(1, .9999, .9999, 1), 2);lev = c(-1, 0, 1)
	for (k in 1:2) {p = lame:::.lame_ab_cond_prior(S, k, lev);expect_equal(p$mean, .9999 * lev);expect_equal(p$var, 1 - .9999^2)}
	p = lame:::.lame_ab_cond_prior(S * 1e-12, 1, lev)
	expect_equal(p$mean, .9999 * lev);expect_lt(abs(p$var / (1e-12 * (1 - .9999^2)) - 1), 1e-10)
})

test_that('additive AR conditionals preserve stationary one-period and small scales', {
	for (nt in c(1L, 3L)) for (v in c(.2, 1e-14)) {
		rho = .8;M = matrix(seq_len(3 * nt) / 20, 3, nt)
		K = v / (1 - rho^2) * rho^abs(outer(seq_len(nt), seq_len(nt), '-'));Q = solve(K)
		for (t in seq_len(nt)) {
			z = lame:::.lame_ab_ar1_cond(M, t, nt, rho, sqrt(v))
			mu = if (nt == 1L) rep(0, 3) else -drop(M[, -t, drop = FALSE] %*% Q[-t, t]) / Q[t, t]
			expect_equal(z$m, mu, tolerance = 1e-10);expect_lt(abs(z$v * Q[t, t] - 1), 1e-10)
		}
		st = lame:::.lame_ab_level_path_stats(M, rho, sqrt(v))
		expect_lt(max(abs(st$prec / sum(Q) - 1)), 1e-10)
		expect_lt(max(abs(st$num / drop(M %*% rowSums(Q)) - 1)), 1e-10)
	}
})

test_that('tiny joint additive levels retain their coupled constrained Gaussian law', {
	set.seed(145867);n = 3L;N = 6000L
	S = matrix(c(1, .8, .8, 2), 2) * 1e-10;pa = c(1, 3, 8) * 1e10;pb = c(4, 2, 1) * 1e10;ha = c(1, -2, .5) * 1e5;hb = c(-1, .5, 1) * 1e5
	Q = kronecker(diag(n), solve(S)) + diag(c(rbind(pa, pb)))
	basis = kronecker(qr.Q(qr(contr.helmert(n))), diag(2));V = solve(crossprod(basis, Q %*% basis));mu = drop(
		basis %*% (V %*% crossprod(basis, c(rbind(ha, hb))))
	);C = basis %*% V %*% t(basis)
	x = replicate(N, {z = lame:::.lame_ab_level_joint(ha, pa, hb, pb, S);c(rbind(z$a, z$b))})
	expect_lt(max(abs((rowMeans(x) - mu) / sqrt(diag(C) / N))), 6)
	se = sqrt((outer(diag(C), diag(C)) + C^2) / N)
	expect_lt(max(abs((cov(t(x)) - C) / se)), 6)
	expect_lt(max(abs(colSums(x[c(1, 3, 5), ]))), 1e-15)
	expect_lt(max(abs(colSums(x[c(2, 4, 6), ]))), 1e-15)
})

test_that('native additive draws retain small stationary prior variances', {
	set.seed(763127);n = 3L;N = 5000L;rho = .8;variance = 1e-14
	Y = array(NA_real_, c(n, n, 1));X = array(0, dim(Y));target = (1 - 1 / n) * variance / (1 - rho^2)
	z = replicate(
		N,
		lame:::sample_dynamic_ab_cpp(
			matrix(0, n, 1),
			matrix(0, n, 1),
			Y,
			X,
			rho,
			sqrt(variance),
			1,
			FALSE,
			0,
			TRUE,
			TRUE,
			FALSE
		)$a[, 1]
	)
	expect_lt(max(abs(rowMeans(z) / sqrt(target / N))), 6)
	expect_lt(max(abs((rowMeans(z * z) - target) / sqrt(2 * target^2 / N))), 6)
	expect_lt(max(abs(colSums(z))), 1e-15)
})
