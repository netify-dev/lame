# the noise-variance correction inverts the additive-block coefficient matrix, singular in the a_i + d / b_j - d gauge
# and once more per extra component of a disconnected graph; the result is the minimum-norm inverse whether or not Cholesky succeeds

# weights of a bipartite panel whose rows and columns fall into disjoint blocks
block_weights = function(sizes_r, sizes_c, seed = 1) {
	set.seed(seed)
	W = matrix(0, sum(sizes_r), sum(sizes_c))
	r0 = 0; c0 = 0
	for (k in seq_along(sizes_r)) {
		ri = r0 + seq_len(sizes_r[k]); ci = c0 + seq_len(sizes_c[k])
		W[ri, ci] = matrix(runif(length(ri) * length(ci), 0.5, 3), length(ri), length(ci))
		r0 = r0 + sizes_r[k]; c0 = c0 + sizes_c[k]
	}
	W
}
additive_matrix = function(W) {
	rbind(cbind(diag(rowSums(W), nrow(W)), W), cbind(t(W), diag(colSums(W), ncol(W))))
}
# centered traces of the sender and receiver blocks of the inverse
block_traces = function(S, na) {
	nb = nrow(S) - na
	c(
		a = lame:::.ae_centered_trace(S[seq_len(na), seq_len(na), drop = FALSE]),
		b = lame:::.ae_centered_trace(S[na + seq_len(nb), na + seq_len(nb), drop = FALSE])
	)
}
permute_design = function(K, na, perm_a, perm_b) {
	idx = c(perm_a, na + perm_b)
	K[idx, idx]
}

test_that("a disconnected additive block gives the minimum-norm inverse, whatever the rounding", {
	W = block_weights(c(5, 4), c(4, 6), seed = 11)
	na = nrow(W)
	K = additive_matrix(W)
	null_vec = c(rep(1, na), rep(-1, ncol(W)))
	reference = block_traces(MASS::ginv(K), na)
	S = lame:::.ae_gauge_inverse(K, null_vec)
	expect_true(all(is.finite(S)))
	expect_lt(max(abs(S)), 10)
	expect_equal(block_traces(S, na), reference, tolerance = 1e-8)
	# the same system reordered and with rounding-level noise: the factorization may change, the traces do not
	set.seed(12)
	for (rep in 1:8) {
		perm_a = sample(na); perm_b = sample(ncol(W))
		Kp = permute_design(K, na, perm_a, perm_b)
		noise = matrix(rnorm(length(Kp)), nrow(Kp)); noise = (noise + t(noise)) / 2
		Kp = Kp + 1e-15 * mean(diag(Kp)) * noise
		nvp = c(rep(1, na), rep(-1, ncol(W)))
		Sp = lame:::.ae_gauge_inverse(Kp, nvp)
		expect_lt(max(abs(Sp)), 10)
		expect_equal(block_traces(Sp, na), reference, tolerance = 1e-6)
	}
})

test_that("a symmetric additive block that is singular gives the minimum-norm inverse", {
	# complete bipartite adjacency between two groups of three: the alternating vector is a null direction of diag(rowsums) + w
	A = matrix(0, 6, 6); A[1:3, 4:6] = 1; A[4:6, 1:3] = 1
	K = diag(rowSums(A)) + A
	S = lame:::.ae_gauge_inverse(K)
	expect_lt(max(abs(S)), 10)
	expect_equal(
		lame:::.ae_centered_trace(S),
		lame:::.ae_centered_trace(MASS::ginv(K)),
		tolerance = 1e-8
	)
})

test_that("a connected additive block keeps the plain inverse", {
	set.seed(21)
	W = matrix(runif(7 * 5, 0.5, 3), 7, 5)
	K = additive_matrix(W)
	null_vec = c(rep(1, 7), rep(-1, 5))
	Kg = K + (mean(diag(K)) / nrow(K)) * tcrossprod(null_vec)
	expect_identical(lame:::.ae_gauge_inverse(K, null_vec), chol2inv(chol(Kg)))
})

test_that("the noise variances of a disconnected design are finite and match the pseudo-inverse", {
	W = block_weights(c(6, 5), c(5, 6), seed = 31)
	nr = nrow(W); nc = ncol(W); Tt = 2L
	mask = array(rep(W > 0, Tt), c(nr, nc, Tt))
	Wobs = array(rep(W, Tt), c(nr, nc, Tt))
	ve = 0.8
	# unit weights (no working weights) and the working weights of an IRLS fit
	totals = list(unit = (W > 0) * Tt, working = W * Tt)
	for (kind in names(totals)) {
		reference = ve * block_traces(MASS::ginv(additive_matrix(totals[[kind]])), nr)
		weights = if (kind == "unit") NULL else Wobs
		got = lame:::.ae_effect_noise_var(mask, weights, ve, FALSE, rep(1, nr), rep(1, nc))
		expect_true(all(is.finite(unlist(got))))
		expect_lt(max(unlist(got)), 10)
		expect_equal(unname(unlist(got)), unname(reference), tolerance = 1e-8)
	}
})

test_that("a fit to a disconnected panel reports finite noise variances", {
	skip_on_cran()
	set.seed(41)
	n = 16L
	nm = sprintf("a%02d", seq_len(n))
	Yl = lapply(1:3, function(t) {
		y = outer(rnorm(n, 0, 0.5), rnorm(n, 0, 0.5), "+") + matrix(rnorm(n * n), n, n)
		diag(y) = NA
		y[1:8, 9:16] = NA; y[9:16, 1:8] = NA
		dimnames(y) = list(nm, nm)
		y
	})
	fit = suppressWarnings(suppressMessages(lame_als(
		Yl,
		R = 0,
		family = "normal",
		verbose = FALSE,
		bias_correct = "none"
	)))
	expect_true(all(is.finite(fit$vc_noise)))
	expect_lt(max(fit$vc_noise), 1)
	expect_true(all(is.finite(fit$VC[c("va", "vb")])))
	expect_lt(max(fit$VC[c("va", "vb")]), 10)
})
