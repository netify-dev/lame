test_that('signed symmetric start recovers the sign of each leading eigenvalue', {
	set.seed(4401); n = 15L
	u = qr.Q(qr(matrix(rnorm(n * 2L), n, 2L)))
	E = u %*% diag(c(-6, 3)) %*% t(u)
	Z = array(E, c(n, n, 3L)); for (t in 1:3) diag(Z[, , t]) = NA
	s = .lame_signed_sym_start(Z, 2L)
	L0 = colSums(s$U * s$V) / colSums(s$U^2)
	expect_equal(sort(L0), c(-1, 1))
	# the leading factor is the one with the largest absolute eigenvalue
	expect_equal(unname(L0[1]), -1)
	expect_equal(abs(cor(s$U[, 1], u[, 1])), 1, tolerance = 0.05)
	expect_null(.lame_signed_sym_start(array(0, c(n, n, 2L)), 2L))
})

test_that('a static symmetric binary panel recovers a heterophilous dimension', {
	skip_on_cran()
	# started at U = V = 0 this chain stays in the all-positive sector and fits the truth poorly
	set.seed(1005); n = 18L; nt = 9L; nm = sprintf('c%02d', 1:n)
	a = rnorm(n, 0, 0.5); u = matrix(rnorm(n * 2L), n, 2L)
	EZ = -1 + outer(a, a, '+') + u %*% diag(c(1.2, -0.8)) %*% t(u)
	Y = lapply(seq_len(nt), function(t) {
		e = matrix(rnorm(n * n), n, n); e[lower.tri(e)] = t(e)[lower.tri(e)]
		y = 1 * (EZ + e > 0); diag(y) = NA; dimnames(y) = list(nm, nm); y
	})
	# the start values are under test, so the collapsed actor update is off (its recovery is checked in
	# test_binary_actor_collapsed_move.R); the intercept still drifts here, so the chain policy continues the chain
	fit = suppressMessages(lame(
		Y,
		family = 'binary',
		symmetric = TRUE,
		R = 2L,
		burn = 200L,
		nscan = 1000L,
		odens = 10L,
		seed = 12L,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = list(actor_move = FALSE)
	))
	expect_gte(fit$auto_extend$extensions, 1L)
	expect_equal(sort(sign(fit$L)), c(-1, 1))
	truth = EZ[upper.tri(EZ)]
	fitted = fit$EZ[[1]][upper.tri(EZ)]
	expect_gt(cor(truth, fitted), 0.85)
})
