.static_marginal_fixture = function(rank = 1L, uv_prior = NULL) {
	n = 5L; nt = 3L
	X = lapply(
		seq_len(nt),
		function(t) array(c(rep(1, n * n), sin(seq_len(n * n) / 7 + t), cos(seq_len(n * n) / 9 + t)), c(n, n, 3L))
	)
	Y = array(0, c(n, n, nt))
	for (t in seq_len(nt)) {
		Y[, , t] = outer(seq_len(n), seq_len(n), function(i, j) (2 * i + j + t) %% 3 == 0)
		diag(Y[, , t]) = NA_real_
	}
	Y[2, 3, 1] = NA_real_
	cx = lame:::.lame_binary_static_hmc_context(Y, X, c(1L, 2L), 3L, rank, uv_prior)
	h = list(
		g = c(1.3, 2.1),
		rho_beta = .6,
		var_beta = .4,
		iC0 = matrix(1 / 2.3),
		iSab = solve(matrix(c(.8, .2, .2, .5), 2))
	)
	size = 2L + nt + 2L * (n - 1L) + 2L * n * rank
	q = seq_len(size) / size / 3 - .15
	list(Y = Y, X = X, cx = cx, h = h, q = q, n = n, nt = nt, rank = rank)
}

test_that('static actor observed density agrees with the original independent target', {
	# the unit prior used when the host passes none, and the data-scaled prior of a sparse network
	# (R/uv_prior.R: round((2R + 2) m) degrees of freedom, scale matrix df s I)
	for (rank in 1:2) for (k in c(0, 9)) {
		df0 = 2 * rank + 2 + k; psi0 = if (k > 0) df0 * .3 else 2 * rank + 2
		z = .static_marginal_fixture(rank, if (k > 0) list(df = df0, psi = psi0))
		n = z$n; nt = z$nt; H = qr.Q(qr(contr.helmert(n)))
		reference = function(q) {
			beta = q[1:2]; path = q[2 + seq_len(nt)]
			level = matrix(q[2 + nt + seq_len(2 * (n - 1))], n - 1, 2)
			start = 2 + nt + 2 * (n - 1)
			U = matrix(q[start + seq_len(n * rank)], n, rank)
			V = matrix(q[start + n * rank + seq_len(n * rank)], n, rank)
			actor = H %*% level
			answer = 0
			for (t in seq_len(nt)) for (i in seq_len(n)) for (j in seq_len(n)) {
				y = z$Y[i, j, t]
				if (i == j || is.na(y)) next
				eta = sum(z$X[[t]][i, j, 1:2] * beta) + z$X[[t]][i, j, 3] * path[t] +
					actor[i, 1] + actor[j, 2] + sum(U[i, ] * V[j, ])
				answer = answer + pnorm((2 * y - 1) * eta, log.p = TRUE)
			}
			answer = answer - sum(beta^2 / z$h$g) / 2 - path[1]^2 / 2.3 / 2
			for (t in 2:nt) answer = answer - (path[t] - .6 * path[t - 1])^2 / .4 / 2
			for (i in seq_len(n - 1)) answer = answer - drop(t(level[i, ]) %*% z$h$iSab %*% level[i, ]) / 2
			W = cbind(U, V)
			answer - (n + df0) / 2 * as.numeric(
				determinant(psi0 * diag(2 * rank) + crossprod(W), logarithm = TRUE)$modulus
			)
		}
		got = lame:::.lame_binary_static_hmc_evaluate(z$q, z$cx, z$h)
		expect_equal(got$logp, reference(z$q), tolerance = 1e-12)
		finite_difference = vapply(seq_along(z$q), function(j) {
			plus = minus = z$q; plus[j] = plus[j] + 1e-5; minus[j] = minus[j] - 1e-5
			(reference(plus) - reference(minus)) / 2e-5
		}, numeric(1))
		expect_equal(got$gradient, finite_difference, tolerance = 1e-7)
	}
})

test_that('bounded static actor proposals reverse with the original phase-space state', {
	for (rank in 1:2) {
		z = .static_marginal_fixture(rank)
		mass = lame:::.lame_binary_static_hmc_mass(z$cx, z$h)
		for (scale in c(1, 30)) {
			q = z$q * scale; p = cos(seq_along(q)) * sqrt(mass)
			a = lame:::.lame_binary_static_hmc_integrate(q, p, z$cx, z$h, mass, .08, 7L)
			expect_false(is.null(a))
			b = lame:::.lame_binary_static_hmc_integrate(a$q, -a$p, z$cx, z$h, mass, .08, 7L)
			expect_equal(b$q, q, tolerance = 1e-10)
			expect_equal(b$p, -p, tolerance = 1e-10)
		}
	}
})

test_that('static actor joint moves cover changing intercepts and pooled coefficient paths', {
	skip_on_cran()
	withr::local_seed(14901)
	n = 7L; actors = paste0('a', c(7, 2, 5, 1, 4, 6, 3))
	Y = lapply(seq_len(3), function(t) {
		y = matrix(rbinom(n * n, 1, .3), n, n, dimnames = list(actors, actors))
		diag(y) = NA_real_; y
	})
	X = lapply(
		seq_len(3),
		function(t) array(
			rnorm(n * n * 2),
			c(n, n, 2L),
			dimnames = list(actors, actors, c('x1', 'x2'))
		)
	)
	caller = .Random.seed
	for (rank in 1:2) for (kind in c('ar1', 'rw1')) {
		f = lame(
			Y,
			Xdyad = X,
			family = 'binary',
			R = rank,
			dynamic_ab = FALSE,
			dynamic_beta = c('intercept', 'x1_dyad'),
			dynamic_beta_kind = kind,
			dynamic_beta_pool = 'both',
			dcor = FALSE,
			seed = 7886,
			burn = 100,
			nscan = 100,
			odens = 10,
			gof = FALSE,
			verbose = FALSE
		)
		expect_identical(.Random.seed, caller)
		expect_equal(unname(f$marginal_move['attempts']), 8)
		expect_equal(unname(f$marginal_move['numerical_rejections']), 0)
		expect_true(all(unlist(f$tryErrorChecks) == 0))
		expect_true(all(is.finite(f$BETA)))
	}
})
