skip_on_cran()

library(lame)

# the bipartite cross-section sampler stores its draws in the identified form of the one-mode sampler and
# of lame(): positions centred over the actors, the means of the product in the effects and the intercept

bi_net = function(fam = "normal", nA = 16, nB = 11, seed = 1) {
	set.seed(seed)
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	U = matrix(rnorm(nA * 2, 0.5, 0.8), nA); V = matrix(rnorm(nB * 2, -0.5, 0.8), nB)
	x = matrix(rnorm(nA * nB), nA, nB)
	z = 0.3 + 0.5 * x + outer(a, b, "+") + tcrossprod(U, V) + matrix(rnorm(nA * nB), nA, nB)
	y = switch(fam, normal = z, binary = 1 * (z > 0), ordinal = matrix(findInterval(z, c(-1, 0.5, 2)), nA, nB))
	dimnames(y) = list(sprintf("r%02d", seq_len(nA)), sprintf("c%02d", seq_len(nB)))
	list(Y = y, X = array(x, c(nA, nB, 1), dimnames = c(dimnames(y), list("x"))))
}
bi_fit = function(d, fam, ...) suppressWarnings(ame(
	d$Y, Xdyad = d$X, family = fam, mode = "bipartite", burn = 50, nscan = 200, odens = 5, chains = 1,
	verbose = FALSE, gof = FALSE,
	posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1L, thin_UV = 1L), ...
))

test_that("stored draws are centred and keep the predictor of the sampled state", {
	for (fam in c("normal", "binary", "ordinal")) for (rc in list(c(2, 2), c(2, 1))) {
		lab = paste(fam, rc[1], rc[2])
		f = bi_fit(bi_net(fam, seed = 3), fam, R_row = rc[1], R_col = rc[2])
		expect_lt(max(abs(colMeans(f$U))), 1e-10, label = lab)
		expect_lt(max(abs(colMeans(f$V))), 1e-10, label = lab)
		expect_lt(max(abs(rowMeans(f$UVPM))), 1e-10, label = lab)
		expect_lt(max(abs(colMeans(f$UVPM))), 1e-10, label = lab)
		# the last stored draw is the last sweep, whose sampled state start_vals keeps
		ns = ncol(f$a_samples); sv = f$start_vals
		a_s = f$a_samples[, ns]; b_s = f$b_samples[, ns]
		U_s = matrix(f$U_samples[, , ns], nrow(f$U)); V_s = matrix(f$V_samples[, , ns], nrow(f$V))
		G_s = matrix(f$G_samples[, , ns], rc[1], rc[2])
		expect_lt(max(abs(colMeans(U_s))), 1e-10, label = lab)
		has_int = "intercept" %in% colnames(f$BETA)
		expect_identical(has_int, fam != "ordinal", label = lab)
		b_last = f$BETA[nrow(f$BETA), ]
		shown = (if (has_int) b_last[[1]] else 0) + outer(a_s, b_s, "+") + U_s %*% G_s %*% t(V_s)
		drawn = (if (has_int) sv$beta[1] else 0) + outer(sv$a, sv$b, "+") + sv$U %*% sv$G %*% t(sv$V)
		expect_equal(unname(shown), unname(drawn), tolerance = 1e-10, label = lab)
		slopes = if (has_int) -1 else seq_along(b_last)
		expect_equal(unname(b_last[slopes]), sv$beta[slopes], label = lab)
		expect_gt(max(abs(colMeans(sv$U))), 1e-6, label = lab)
		expect_equal(unname(f$APM), rowMeans(f$a_samples), tolerance = 1e-10, label = lab)
		expect_equal(unname(f$BPM), rowMeans(f$b_samples), tolerance = 1e-10, label = lab)
	}
})

test_that("a fit without one margin keeps the sampled coordinates", {
	f = bi_fit(bi_net("normal", seed = 4), "normal", R = 2, cvar = FALSE)
	expect_gt(max(abs(colMeans(f$U))), 1e-6)
	expect_identical(f$U, f$start_vals$U)
})
