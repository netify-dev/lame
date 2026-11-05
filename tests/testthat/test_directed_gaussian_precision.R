test_that("directed effects use the exact joint Gaussian under near reciprocal data", {
	n = 5L
	rho = .9999
	Sab = diag(c(.3, .2))
	X = precomputeX(design_array(intercept = TRUE, n = n))
	Z = matrix(seq_len(n * n) / (n * n), n)
	# cell design in physical beta, sender and receiver coordinates
	D = cbind(1, kronecker(matrix(1, n, 1), diag(n)), kronecker(diag(n), matrix(1, n, 1)))
	W = mhalf(solve(matrix(c(1, rho, rho, 1), 2)))
	perm = as.vector(t(matrix(seq_len(n * n), n)))
	Dw = W[1, 1] * D + W[1, 2] * D[perm, ]
	zw = W[1, 1] * as.vector(Z) + W[1, 2] * as.vector(t(Z))
	Q = crossprod(Dw) + diag(c(2, rep(1 / Sab[1, 1], n), rep(1 / Sab[2, 2], n)))
	C = solve(Q)
	mu = drop(C %*% crossprod(Dw, zw))
	set.seed(47327902)
	N = 3000L
	x = replicate(N, {
		a = rbeta_ab_fc(Z, Sab, rho, X, iV0 = matrix(2, 1, 1))
		c(a$beta, a$a, a$b)
	})
	mean_z = (rowMeans(x) - mu) / sqrt(diag(C) / N)
	cov_z = (cov(t(x)) - C) / sqrt((C^2 + outer(diag(C), diag(C))) / (N - 1))
	expect_lt(max(abs(mean_z)), 4.5)
	expect_lt(max(abs(cov_z)), 4.5)
})
