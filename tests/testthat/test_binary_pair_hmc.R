.pair_hmc_fixture = function() {
	set.seed(918365)
	n = 3L; nt = 2L
	Y = array(NA_real_, c(n, n, nt))
	Y[cbind(c(1, 2, 2, 3, 1, 3), c(2, 1, 3, 2, 3, 1), c(1, 1, 1, 1, 2, 2))] = c(1, 0, 0, 0, 0, 1)
	Y[2, 1, 2] = 0
	X = lapply(1:nt, function(t) {
		x = array(rnorm(n * n * 2), c(n, n, 2)); x[, , 1] = 1; x
	})
	P = matrix(c(3, .2, .2, 2), 2); Sab = matrix(c(.2, .03, .03, .15), 2)
	cx = .lame_binary_dynamic_pair_hmc_context(Y, X, 2L)
	h = .lame_binary_dynamic_pair_hmc_hyper(cx, P, Sab, .6, .25, .7, .2, -.4)
	K = function(r, s) s * s / (1 - r * r) * r^abs(outer(1:nt, 1:nt, '-'))
	ia = solve(K(.6, .25)); iu = solve(K(.7, .2)); il = solve(Sab)
	density = function(q) {
		z = .lame_binary_dynamic_pair_hmc_decode(q, cx, h); lp = 0
		for (t in 1:nt) {
			eta = matrix(matrix(X[[t]], n * n, 2) %*% z$beta, n, n) +
				outer(z$level[, 1] + z$adev[, t], z$level[, 2] + z$bdev[, t], '+') +
				z$U[, , t] %*% t(z$V[, , t])
			for (i in 1:(n - 1)) for (j in (i + 1):n) {
				yy = c(Y[i, j, t], Y[j, i, t]); v = c(eta[i, j], eta[j, i]); ok = is.finite(yy)
				if (all(ok)) {
					sign = 2 * yy - 1; r = h$rho * prod(sign)
					lp = lp + log(
						integrate(
							function(z) dnorm(z) *
								pnorm((sign[2] * v[2] - r * z) / sqrt(1 - r * r)),
							lower = -Inf,
							upper = sign[1] * v[1],
							rel.tol = 1e-12
						)$value
					)
				} else if (any(ok)) lp = lp + pnorm((2 * yy[ok] - 1) * v[ok], log.p = TRUE)
			}
		}
		u = matrix(z$U, ncol = nt); v = matrix(z$V, ncol = nt)
		lp - .5 * (drop(z$beta %*% P %*% z$beta) + sum(z$level * (z$level %*% il)) +
				sum(z$adev * (z$adev %*% ia)) + sum(z$bdev * (z$bdev %*% ia)) +
				sum(u * (u %*% iu)) + sum(v * (v %*% iu)))
	}
	list(cx = cx, h = h, density = density, q = rnorm(cx$size, sd = .3))
}

test_that('observed-pair Hamiltonian density uses reciprocal and missing-pair likelihoods', {
	f = .pair_hmc_fixture(); q2 = f$q + rnorm(length(f$q), sd = .1)
	value = function(q) .lame_binary_dynamic_pair_hmc_evaluate(q, f$cx, f$h, NULL)$logp
	expect_equal(value(f$q) - value(q2), f$density(f$q) - f$density(q2), tolerance = 1e-10)
	eps = 1e-5
	grad = vapply(seq_along(f$q), function(j) {
		plus = minus = f$q; plus[j] = plus[j] + eps; minus[j] = minus[j] - eps
		(f$density(plus) - f$density(minus)) / (2 * eps)
	}, numeric(1))
	actual = .lame_binary_dynamic_pair_hmc_evaluate(f$q, f$cx, f$h, NULL)$gradient
	expect_lt(max(abs(actual - grad)), 1e-6)
})

test_that('observed-pair likelihood gradients remain accurate below log probability -1600', {
	Y = array(NA_real_, c(2, 2, 1)); Y[1, 2, 1] = Y[2, 1, 1] = 1
	cx = .lame_binary_dynamic_collapsed_rho_context(Y); eta = array(-4, c(2, 2, 1))
	a = .lame_binary_dynamic_pair_likelihood_gradient(eta, -.99, cx); eps = 1e-5
	numerical = vapply(c(2, 3), function(j) {
		plus = minus = eta; plus[j] = plus[j] + eps; minus[j] = minus[j] - eps
		(.lame_binary_dynamic_pair_likelihood_gradient(plus, -.99, cx)$logp -
				.lame_binary_dynamic_pair_likelihood_gradient(minus, -.99, cx)$logp) / (2 * eps)
	}, numeric(1))
	expect_lt(a$logp, -1600); expect_true(all(is.finite(a$gradient)))
	expect_lt(max(abs(a$gradient[c(2, 3)] - numerical)), 1e-5)
})
