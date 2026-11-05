test_that("symmetric factor draws preserve very small conditional variances", {
	withr::local_seed(47318003L)
	n = 4L; U = V = matrix(1, n, 1L); E = matrix(0, n, n); s2 = 1e-14
	variance = 1 / (1 / n + (n - 1L) / s2)
	draws = replicate(1600L, lame:::rUV_sym_fc_cpp(E, U, V, s2, FALSE, rep(1L, 4L))$U[2L, 1L]) / sqrt(variance)
	expect_lt(abs(mean(draws)), 6 / sqrt(length(draws)))
	expect_lt(abs(mean(draws^2) - 1), 6 * sqrt(2 / length(draws)))
})

test_that("a zero first actor does not erase a nonzero symmetric factor weight", {
	withr::local_seed(47318004L)
	n = 4L; U = matrix(c(0, 1, 2, 3), n, 1L); V = 2 * U
	E = matrix(.2, n, n); diag(E) = 0; s2 = .7
	precision = 1 / n + 4 * sum(U[-1L, 1L]^2) / s2
	mean = (2 * sum(U[-1L, 1L] * E[1L, -1L]) / s2) / precision
	draws = replicate(2400L, lame:::rUV_sym_fc_cpp(E, U, V, s2, FALSE, 0L)$U[1L, 1L])
	standardized = (draws - mean) * sqrt(precision)
	expect_lt(abs(base::mean(standardized)), 6 / sqrt(length(draws)))
	expect_lt(abs(base::mean(standardized^2) - 1), 6 * sqrt(2 / length(draws)))
})

test_that("symmetric row updates retain correlated Gaussian conditional draws", {
	withr::local_seed(47318005L)
	U = matrix(c(.2, 1, -.4, .6, 1.1, .7, -.6, .9, .2, -.8), 5L, 2L)
	L = c(.8, -1.2); V = sweep(U, 2L, L, "*")
	E = outer(1:5, 1:5, function(i, j) sin(i + j)); s2 = .6
	D = sweep(U[-3L, , drop = FALSE], 2L, L, "*")
	Q = diag(1 / 5, 2L) + crossprod(D) / s2
	mean = solve(Q, crossprod(D, E[3L, -3L]) / s2)
	draws = replicate(2400L, lame:::rUV_sym_fc_cpp(E, U, V, s2, FALSE, 2L)$U[3L, ])
	standardized = chol(Q) %*% (draws - drop(mean))
	expect_lt(max(abs(rowMeans(standardized))), 6 / sqrt(ncol(draws)))
	expect_lt(max(abs(tcrossprod(standardized) / ncol(draws) - diag(2L))), 6 * sqrt(2 / ncol(draws)))
})
