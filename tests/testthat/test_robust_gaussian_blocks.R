# independent density Hessians, not copies of the compiled conditionals
skip_on_cran()

robust_quadratic_reference = function(negative_log_density, dimension) {
	unit = diag(dimension)
	f = negative_log_density
	f0 = f(rep(0, dimension))
	h = vapply(seq_len(dimension), function(i) (f(-unit[i, ]) - f(unit[i, ])) / 2, numeric(1))
	Q = matrix(0, dimension, dimension)
	for (i in seq_len(dimension)) for (j in seq_len(dimension))
		Q[i, j] = f(unit[i, ] + unit[j, ]) - f(unit[i, ]) - f(unit[j, ]) + f0
	S = solve(Q)
	list(mean = drop(S %*% h), covariance = S)
}

robust_check_gaussian = function(draws, target) {
	N = nrow(draws)
	S = target$covariance
	expect_lt(max(abs(colMeans(draws) - target$mean) / sqrt(diag(S) / N)), 6)
	se = sqrt((outer(diag(S), diag(S)) + S^2) / (N - 1))
	expect_lt(max(abs(cov(draws) - S) / se), 6)
}

test_that("weighted coefficient and intercept draws match independent Gaussian densities", {
	set.seed(74191)
	n = 3L
	nt = 3L
	Y = array(rnorm(n * n * nt), c(n, n, nt))
	W = array(runif(length(Y), .1, 2), dim(Y))
	W[c(1, 5, 9, 18)] = 0
	v = .7
	rho = -.35
	iv = .4
	X = matrix(rnorm(length(Y) * 2), length(Y), 2)
	precision = c(.3, .8)
	f = function(b) sum(W * (Y - array(X %*% b, dim(Y)))^2) / (2 * v) + sum(precision * b^2) / 2
	target = robust_quadratic_reference(f, 2)
	draws = t(replicate(4000, as.numeric(lame:::robust_weighted_beta(Y, W, X, v, precision))))
	robust_check_gaussian(draws, target)
	f = function(mu) {
		ar = (1 - rho^2) * mu[1]^2 + sum((mu[-1] - rho * mu[-nt])^2)
		ar / (2 * iv) + sum(W * (Y - array(rep(mu, each = n * n), dim(Y)))^2) / (2 * v)
	}
	target = robust_quadratic_reference(f, nt)
	draws = t(replicate(4000, as.numeric(lame:::robust_weighted_mu(Y, W, v, rho, iv))))
	robust_check_gaussian(draws, target)
})

test_that("weighted actor paths obey the exact constrained Gaussian posterior", {
	set.seed(52291)
	n = 3L
	nt = 2L
	Y = array(rnorm(n * n * nt), c(n, n, nt))
	W = array(runif(length(Y), .3, 2), dim(Y))
	v = .7
	rho = .4
	iv = .3
	Sab = matrix(c(.8, .2, .2, .5), 2)
	other = c(-.2, .1, .1)
	q = nt + 1L
	for (tau in list(NULL, c(.2, 3))) {
		f = function(x) {
			X = matrix(x, n, q, byrow = TRUE)
			a = X[, 1:nt, drop = FALSE]
			level = X[, q]
			value = 0
			for (i in seq_len(n)) {
				delta = a[i, ] - level[i]
				value = value + (if (is.null(tau)) ((1 - rho^2) * delta[1]^2 + (delta[2] - rho * delta[1])^2) else (tau[1] * (1 - rho^2) * delta[1]^2 + tau[2] * (delta[2] - rho * delta[1])^2)) / (2 * iv)
				lev = c(level[i], other[i])
				value = value + drop(t(lev) %*% solve(Sab, lev)) / 2
				for (t in seq_len(nt)) value = value + sum(W[i, , t] * (Y[i, , t] - a[i, t])^2) / (2 * v)
			}
			value
		}
		target = robust_quadratic_reference(f, n * q)
		S = target$covariance
		L = kronecker(matrix(1, 1, n), diag(q))
		SL = S %*% t(L)
		target$covariance = S - SL %*% solve(L %*% SL, t(SL))
		target$mean = drop(target$mean - SL %*% solve(L %*% SL, L %*% target$mean))
		draws = t(replicate(4000, {
			z = lame:::robust_weighted_actor(Y, W, v, rho, iv, Sab, other, TRUE, tau)
			as.numeric(t(cbind(z$path, z$level)))
		}))
		expect_lt(max(abs(draws %*% t(L))), 1e-10)
		robust_check_gaussian(draws, target)
		# the receiver conditional is the sender conditional of the transposed data
		set.seed(327)
		row = lame:::robust_weighted_actor(Y, W, v, rho, iv, Sab, other, TRUE, tau)
		set.seed(327)
		column = lame:::robust_weighted_actor(
			aperm(Y, c(2, 1, 3)),
			aperm(W, c(2, 1, 3)),
			v,
			rho,
			iv,
			Sab[2:1, 2:1],
			other,
			FALSE,
			tau
		)
		expect_equal(row, column)
	}
})

test_that("weighted latent trajectories match an independent AR Gaussian density", {
	set.seed(91367)
	n = 3L
	nt = 2L
	rank = 2L
	Y = array(rnorm(n * n * nt), c(n, n, nt))
	W = array(runif(length(Y), .1, 2), dim(Y))
	W[c(1, 5, 9, 18)] = 0
	U = V = array(rnorm(n * rank * nt), c(n, rank, nt))
	v = .7
	rho = -.4
	iv = .3
	f = function(z) {
		u = matrix(z, rank, nt)
		prior = ((1 - rho^2) * sum(u[, 1]^2) + sum((u[, 2] - rho * u[, 1])^2)) / (2 * iv)
		prior + sum(vapply(
			seq_len(nt),
			function(t)
				sum(W[1, , t] * (Y[1, , t] - V[, , t] %*% u[, t])^2) / (2 * v),
			numeric(1)
		))
	}
	target = robust_quadratic_reference(f, rank * nt)
	draws = t(replicate(4000, as.numeric(lame:::robust_weighted_uv(Y, W, U, V, v, rho, iv)$U[1, , ])))
	robust_check_gaussian(draws, target)
})
