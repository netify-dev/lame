test_that('tempered decoding builds the intercept path of each kind', {
	for (rank in 1:2) for (kind in c('unit', 'hierarchical')) for (path_kind in c('rw1', 'rw2', 'matern32')) {
		n = 5L; nt = 4L
		Y = lapply(1:nt, function(t) {
			y = outer(1:n, 1:n, function(i, j) as.numeric((i + j + t) %% 3 == 0))
			diag(y) = NA_real_; y
		})
		input = lame:::.lame_tempered_input(Y, NULL, rank)
		prior = lame:::.lame_tempered_prior(
			list(
				uv_prior = kind,
				Sab0 = diag(.4, 2),
				etaab = 6,
				rho_ab_mean = .5,
				rho_ab_sd = .2,
				sigma_ab_shape = 3,
				sigma_ab_scale = .3
			),
			input$Y
		)
		design = lame:::.lame_tempered_design(input, NA)
		cx = lame:::.lame_tempered_context(design$Y, design$Xlist, rank, design$Q, prior, path_kind)
		# the gaussian block ends with the nt - 1 intercept increments; log sigma_beta is last
		d = n * rank + rank + cx$p + (n - 1L) * (nt + 1L) + (nt - 1L)
		q = sin(seq_len(d + 4L)) / 5
		q[d + 1:4] = c(-1, -1, atanh(.3 / .99), log(.4))
		# intercept path: zero in the first period, then the scaled change of the kind
		z = lame:::.lame_tempered_decode(q, cx, design$Q)
		w = q[d - (nt - 2L):0]
		change = switch(path_kind, rw1 = cumsum(w), rw2 = cumsum(cumsum(w)), drop(cx$intercept_basis %*% w))
		expect_equal(z$intercept_path, c(0, change) * .4 * sqrt(cx$lam0))
	}
})

test_that('the tempered intercept bases give the covariance of each kind', {
	nt = 6L
	ls = 2.5
	A = lame:::.lame_tempered_intercept_basis('matern32', nt, ls)
	a = sqrt(3) * abs(outer(1:nt, 1:nt, '-')) / ls
	S = (1 + a) * exp(-a)
	D = cbind(-1, diag(nt - 1L))
	expect_equal(tcrossprod(A), D %*% S %*% t(D))
	expect_identical(A[upper.tri(A)], rep(0, sum(upper.tri(A))))
	expect_equal(lame:::.lame_tempered_intercept_basis('matern32', nt), lame:::.lame_tempered_intercept_basis('matern32', nt, 2))
	# rw2: the second differences of the change are the standard coordinates
	w = rnorm(nt - 1L)
	change = c(0, drop(lame:::.lame_tempered_intercept_basis('rw2', nt) %*% w))
	expect_equal(diff(change, differences = 2), w[-1])
	expect_equal(change[2], w[1])
})
