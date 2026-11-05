# independent prior draws and observation densities expose absolute numerical floors
skip_on_cran()

.gibbs_scale_z = function(draws, mean, covariance) {
	N = nrow(draws)
	se = sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) / (N - 1))
	c(
		mean = max(abs(colMeans(draws) - mean) / sqrt(diag(covariance) / N)),
		covariance = max(abs(cov(draws) - covariance) / se)
	)
}

test_that("small inverse-gamma conditionals retain their exact draws and RNG", {
	set.seed(25911)
	reference = 1 / rgamma(4000, 5, 1e-12)
	rng = .Random.seed
	set.seed(25911)
	draws = replicate(4000, lame:::.lame_safe_inv_gamma(5, 1e-12))
	expect_identical(draws, reference)
	expect_identical(.Random.seed, rng)
	expect_false(any(draws == 1e-8))
	a = c(1, -1, 2) * 1e-7
	set.seed(25912)
	reference = 1 / rgamma(4000, 5, (7e-14 + sum(a^2)) / 2)
	set.seed(25912)
	draws = replicate(4000, rSab_fc(a, a, diag(1e-14, 2), 7, symmetric = TRUE)[1, 1])
	expect_identical(draws, reference)
})

test_that("small centered inverse-Wishart updates preserve their joint prior", {
	skip_on_cran()
	set.seed(25913)
	n = 3L; eta = 7; scale = 1e-14
	feature = function(S) c(
		S[1, 1] / scale,
		S[2, 2] / scale,
		S[1, 2] / scale,
		as.numeric(determinant(S / scale, logarithm = TRUE)$modulus)
	)
	delta = t(replicate(4000, {
		S = solve(rWishart(1, eta, diag(2) / (eta * scale))[, , 1])
		ab = matrix(rnorm(2 * n), n, 2) %*% chol(S)
		ab = sweep(ab, 2, colMeans(ab))
		feature(rSab_fc(ab[, 1], ab[, 2], diag(scale, 2), eta, center = TRUE)) - feature(S)
	}))
	expect_lt(max(abs(colMeans(delta) / (apply(delta, 2, sd) / sqrt(nrow(delta))))), 6)
})

test_that("mixed positive actor scales remain supported in both actor kernels", {
	skip_on_cran()
	n = 3L; Sab = diag(c(1, 1e-12))
	ij = expand.grid(i = seq_len(n), j = seq_len(n))
	D = matrix(0, n * n, 2 * n)
	for (k in seq_len(nrow(ij))) D[k, c(ij$i[k], n + ij$j[k])] = 1
	S = solve(crossprod(D) + diag(rep(1 / diag(Sab), each = n)))
	set.seed(25914)
	z = t(replicate(4000, {
		x = suppressWarnings(rbeta_ab_fc(matrix(0, n, n), Sab, 0))
		c(x$a, x$b)
	}))
	expect_lt(max(.gibbs_scale_z(z, numeric(2 * n), S)), 6)
	set.seed(25916)
	z = t(replicate(4000, {
		x = lame:::.lame_static_actor_draw(array(0, c(n, n, 1)), Sab, 0, 1)
		c(as.numeric(x$a), as.numeric(x$b))
	}))
	expect_lt(max(.gibbs_scale_z(z, numeric(2 * n), S)), 6)
})

test_that("replicated R actors retain small SPD and disabled PSD directions", {
	n = 3L; X = array(numeric(0), c(n, n, 0, 1)); Y = array(0, c(n, n, 1))
	set.seed(25915)
	z = t(replicate(4000, {
		x = rbeta_ab_rep_fc(Y, diag(1e-12, 2), 0, X)
		c(x$a, x$b)
	}))
	ij = expand.grid(i = seq_len(n), j = seq_len(n))
	D = matrix(0, n * n, 2 * n)
	for (k in seq_len(nrow(ij))) D[k, c(ij$i[k], n + ij$j[k])] = 1
	S = solve(crossprod(D) + diag(1e12, 2 * n))
	expect_lt(max(.gibbs_scale_z(z, numeric(2 * n), S)), 6)
	set.seed(25941)
	b = t(replicate(4000, {
		x = rbeta_ab_rep_fc(Y, diag(c(0, 1)), 0, X)
		stopifnot(all(x$a == 0))
		x$b
	}))
	expect_lt(max(.gibbs_scale_z(b, numeric(n), diag(1 / (n + 1), n))), 6)
})

test_that("bipartite coefficient prior survives collinear scaled covariates", {
	set.seed(25930)
	X = list(array(1e8, c(2, 3, 2))); Y = array(0, c(2, 3, 1))
	z = replicate(
		4000,
		as.numeric(
			lame:::rbeta_ab_bip_gibbs_cpp(
				Y,
				X,
				matrix(0, 2, 3),
				numeric(2),
				numeric(3),
				.7,
				12,
				1,
				1,
				FALSE,
				FALSE
			)$beta
		)
	)
	contrast = (z[1, ] - z[2, ]) / sqrt(2)
	expect_lt(abs(mean(contrast)) / sqrt(8.4 / length(contrast)), 6)
	expect_lt(abs(var(contrast) / 8.4 - 1) / sqrt(2 / (length(contrast) - 1)), 6)
	set.seed(25999); rng = .Random.seed
	expect_error(
		lame:::rbeta_ab_bip_gibbs_cpp(
			Y,
			list(array(1e16, c(2, 3, 2))),
			matrix(0, 2, 3),
			numeric(2),
			numeric(3),
			1,
			12,
			1,
			1,
			FALSE,
			FALSE
		),
		"numerically unresolved"
	)
	expect_identical(.Random.seed, rng)
})

test_that("weak interaction priors receive no numerical precision ridge", {
	set.seed(25935)
	R = array(0, c(1, 1, 1)); U = V = R
	z = replicate(4000, as.numeric(lame:::sample_G_bip_cpp(R, U, V, 1e-10, 1)))
	expect_lt(abs(var(z) * 1e-10 - 1) / sqrt(2 / (length(z) - 1)), 6)
})

test_that("singular replicated covariance obeys its active Gaussian law", {
	n = 3L; X = array(numeric(0), c(n, n, 0, 1)); Y = array(0, c(n, n, 1))
	for (center in c(FALSE, TRUE)) {
		Q = diag(1 + 2 * n, n) + matrix(2, n, n)
		S = solve(Q)
		if (center) S = (diag(n) - matrix(1 / n, n, n)) / (1 + 2 * n)
		set.seed(25942 + as.integer(center))
		z = t(replicate(4000, {
			x = rbeta_ab_rep_fc(Y, matrix(1, 2, 2), 0, X, center = center)
			stopifnot(max(abs(x$a - x$b)) < 1e-12)
			if (center) stopifnot(abs(sum(x$a)) < 1e-12)
			x$a
		}))
		expect_lt(max(.gibbs_scale_z(z, numeric(n), S)), 6)
	}
})

test_that("nonfinite additive covariance draws fail explicitly", {
	fn = rSab_fc
	mock = new.env(parent = environment(fn))
	mock$rgamma = function(...) 0
	environment(fn) = mock
	expect_error(fn(c(1, -1), c(1, -1), symmetric = TRUE), "conditional draw")
})

test_that("factor recovery retains positive scales below a direct inverse cutoff", {
	n = 3L; Sab = diag(c(1, 1e-20)); Y = array(0, c(n, n, 1))
	X = array(numeric(0), c(n, n, 0, 1))
	ij = expand.grid(i = seq_len(n), j = seq_len(n))
	D = matrix(0, n * n, 2 * n)
	for (k in seq_len(nrow(ij))) D[k, c(ij$i[k], n + ij$j[k])] = c(1, 1e-10)
	S = solve(diag(2 * n) + crossprod(D))
	set.seed(25947)
	z = t(replicate(4000, {
		x = rbeta_ab_rep_fc(Y, Sab, 0, X)
		c(x$a, x$b / 1e-10)
	}))
	expect_lt(max(.gibbs_scale_z(z, numeric(2 * n), S)), 6)
	x = rbeta_ab_rep_fc(Y, matrix(0, 2, 2), 0, X)
	expect_identical(x$a, numeric(n))
	expect_identical(x$b, numeric(n))
})
