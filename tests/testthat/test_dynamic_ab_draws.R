library(lame)
skip_on_cran()

# brute-force constrained conditional of a 2n-vector of levels: build the
# joint normal, condition analytically on both margins summing to zero
constrained_moments = function(num_a, prec_a, num_b, prec_b, Sab) {
	n = length(num_a); iS = solve(Sab)
	Q = matrix(0, 2 * n, 2 * n); ia = seq_len(n); ib = n + ia
	diag(Q)[ia] = prec_a + iS[1, 1]; diag(Q)[ib] = prec_b + iS[2, 2]
	Q[cbind(ia, ib)] = iS[1, 2]; Q[cbind(ib, ia)] = iS[1, 2]
	V = solve(Q); m = V %*% c(num_a, num_b)
	A = matrix(0, 2 * n, 2); A[ia, 1] = 1; A[ib, 2] = 1
	VA = V %*% A; AVA = t(A) %*% VA
	list(
		mean = as.numeric(m - VA %*% solve(AVA, t(A) %*% m)),
		cov = V - VA %*% solve(AVA) %*% t(VA)
	)
}

test_that("the joint level draw is the constrained conditional with the full Sab", {
	set.seed(4)
	n = 6; Sab = matrix(c(0.9, -0.35, -0.35, 0.6), 2, 2)
	pa = runif(n, 0.5, 4); pb = runif(n, 0.5, 4)
	na_ = rnorm(n, 0, 2); nb_ = rnorm(n, 0, 2)
	ref = constrained_moments(na_, pa, nb_, pb, Sab)
	B = 40000
	draws = matrix(NA_real_, B, 2 * n)
	for (b in seq_len(B)) {
		d = lame:::.lame_ab_level_joint(na_, pa, nb_, pb, Sab, symmetric = FALSE, center = TRUE)
		draws[b, ] = c(d$a, d$b)
	}
	expect_lt(max(abs(rowSums(draws[, 1:n]))), 1e-10)
	expect_lt(max(abs(rowSums(draws[, n + 1:n]))), 1e-10)
	z_mean = (colMeans(draws) - ref$mean) / sqrt(pmax(diag(ref$cov), 1e-12) / B)
	expect_lt(max(abs(z_mean)), 4.5)
	# the sender-receiver cross moment carries Sab[1, 2]; the diagonal-prior
	# construction would put it near zero
	emp = cov(draws)
	expect_lt(max(abs(emp - ref$cov)) / max(abs(ref$cov)), 0.08)
	# without the constraint the margins do not sum to zero
	d0 = lame:::.lame_ab_level_joint(na_, pa, nb_, pb, Sab, symmetric = FALSE, center = FALSE)
	expect_gt(abs(sum(d0$a)), 1e-6)
})

test_that("the conditional prior of one margin given the other has the right moments", {
	Sab = matrix(c(0.9, -0.35, -0.35, 0.6), 2, 2)
	b_lev = c(0.5, -0.2, 0.1)
	cp = lame:::.lame_ab_cond_prior(Sab, 1L, b_lev)
	expect_equal(cp$mean, Sab[1, 2] / Sab[2, 2] * b_lev)
	expect_equal(cp$var, Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2])
	cp2 = lame:::.lame_ab_cond_prior(Sab, 2L, b_lev)
	expect_equal(cp2$var, Sab[2, 2] - Sab[1, 2]^2 / Sab[1, 1])
	# a diagonal Sab (bipartite) falls back to the marginal prior
	cpd = lame:::.lame_ab_cond_prior(diag(c(0.9, 0.6)), 1L, b_lev)
	expect_null(cpd$mean); expect_equal(cpd$var, 0.9)
})

test_that("the deviation kernel draws each period on the sum-to-zero subspace when asked", {
	set.seed(9)
	n = 7; TT = 3
	Z = array(rnorm(n * n * TT), c(n, n, TT)); EZ = array(0, c(n, n, TT))
	for (t in seq_len(TT)) diag(Z[, , t]) = NA
	a0 = matrix(0, n, TT); b0 = matrix(0, n, TT)
	for (dyad_rho in c(0, 0.4)) {
		out = lame:::sample_dynamic_ab_cpp(a0, b0, Z, EZ, 0.7, 0.4, 1, FALSE, dyad_rho, center = TRUE)
		expect_true(all(is.finite(out$a)) && all(is.finite(out$b)))
		expect_lt(max(abs(colSums(out$a))), 1e-8)
		expect_lt(max(abs(colSums(out$b))), 1e-8)
		out0 = lame:::sample_dynamic_ab_cpp(a0, b0, Z, EZ, 0.7, 0.4, 1, FALSE, dyad_rho, center = FALSE)
		expect_gt(max(abs(colSums(out0$a))), 1e-6)
	}
	# symmetric: one path, drawn jointly
	Zs = Z; for (t in seq_len(TT)) { Zs[, , t] = (Zs[, , t] + t(Zs[, , t])) / 2; diag(Zs[, , t]) = NA }
	outs = lame:::sample_dynamic_ab_cpp(a0, a0, Zs, EZ, 0.7, 0.4, 1, TRUE, 0, center = TRUE)
	expect_equal(outs$a, outs$b)
	expect_lt(max(abs(colSums(outs$a))), 1e-8)
	# a missing actor-period stays out of the constraint and is not drawn
	Zm = Z; Zm[3, , 2] = NA; Zm[, 3, 2] = NA
	outm = lame:::sample_dynamic_ab_cpp(a0, b0, Zm, EZ, 0.7, 0.4, 1, FALSE, 0, center = TRUE)
	expect_true(all(is.finite(outm$a)))
})

test_that("the likelihood-route level draw centres by conditioning", {
	set.seed(3)
	n = 6; TT = 2
	resid = array(rnorm(n * n * TT), c(n, n, TT)); for (t in 1:TT) diag(resid[, , t]) = NA
	dev = matrix(0, n, TT); other = matrix(0, n, TT)
	for (rho in c(0, 0.3)) {
		lev = lame:::.lame_ab_level_draw_lik(
			resid,
			dev,
			other,
			0.8,
			1,
			1L,
			rho = rho,
			skip_diag = TRUE,
			mean_level = rep(0.1, n),
			center = TRUE
		)
		expect_length(lev, n)
		expect_lt(abs(sum(lev)), 1e-10)
	}
})

test_that("bipartite levels are drawn per margin with no pair coupling", {
	set.seed(6)
	na = 7; nb = 5
	d = expect_silent(lame:::.lame_ab_level_joint(
		rnorm(na),
		runif(na, 1, 3),
		rnorm(nb),
		runif(nb, 1, 3),
		diag(c(0.9, 0.6)),
		symmetric = FALSE,
		center = TRUE
	))
	expect_length(d$a, na); expect_length(d$b, nb)
	expect_lt(abs(sum(d$a)), 1e-10); expect_lt(abs(sum(d$b)), 1e-10)
})

test_that("summary of a dynamic_ab fit carries one logical dynamic_ab entry", {
	set.seed(44)
	n = 10
	nm = paste0("a", seq_len(n))
	Yl = lapply(1:3, function(t) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm))
		diag(y) = NA
		y
	})
	fit = suppressWarnings(lame(
		Yl,
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		burn = 10,
		nscan = 30,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		chains = 1
	))
	s = summary(fit)
	expect_equal(sum(names(s) == "dynamic_ab"), 1L)
	expect_identical(s$dynamic_ab, TRUE)
	expect_match(
		paste(capture.output(print(s)), collapse = "\n"),
		"Dynamic additive effects: enabled"
	)
	# an object without the flag reads as FALSE rather than NULL
	fit$dynamic_ab = NULL
	expect_identical(summary(fit)$dynamic_ab, FALSE)
})
