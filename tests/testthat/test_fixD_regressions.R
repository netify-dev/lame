# MCMC correctness checks: unit variance in the symmetric ordinal Z kernel, single counting of dyads in
# symmetric lame(), single storage of posterior_opts draws in ame(), and Sab sampling for poisson; direct
# conditional laws plus short recovery fits against simulated parameters

test_that("#1 symmetric ordinal latent draws have unit variance", {
	set.seed(33501)
	n = 60L
	eta = matrix(rnorm(n * n), n, n)
	y = matrix(1, n, n)
	diag(y) = NA
	z = matrix(0, n, n)
	draw = rZ_ord_sym_fc(z, eta, y)
	upper = upper.tri(draw)
	residual = draw[upper] - ((eta + t(eta)) / 2)[upper]
	m = length(residual)
	# one rank class has no truncation bounds, so the unordered dyads are independent unit normals
	expect_lt(abs(mean(residual)), 5 / sqrt(m))
	expect_lt(abs(var(residual) - 1), 5 * sqrt(2 / (m - 1)))
	expect_equal(draw, t(draw))
})

test_that("symmetric ordinal recovers a known unit-scale dyadic coefficient", {
	skip_on_cran()
	set.seed(3); n = 30
	Xd = matrix(rnorm(n * n), n, n); Xd = (Xd + t(Xd)) / sqrt(2)
	Em = matrix(rnorm(n * n), n, n); Em = (Em + t(Em)) / sqrt(2)
	Zl = 1.0 * Xd + Em
	Y = matrix(
		cut(
			Zl,
			quantile(Zl, seq(0, 1, length = 7), na.rm = TRUE),
			labels = FALSE,
			include.lowest = TRUE
		),
		n,
		n
	)
	diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	Xa = array(Xd, c(n, n, 1), dimnames = list(rownames(Y), colnames(Y), "x1"))
	fs = ame(
		Y,
		Xdyad = Xa,
		R = 0,
		symmetric = TRUE,
		family = "ordinal",
		burn = 300,
		nscan = 1500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 1
	)
	# the simulated residual has variance one, so the coefficient is one
	estimate = mean(fs$BETA[, "x1_dyad"])
	expect_gt(estimate, 0.85)
	expect_lt(estimate, 1.2)
})

test_that("#11 lame() symmetric dyadic-coef sd matches ame() (not sqrt(2) smaller)", {
	skip_on_cran()
	set.seed(7); n = 35
	Xd = matrix(rnorm(n * n), n, n); Xd = (Xd + t(Xd)) / 2
	E = matrix(rnorm(n * n), n, n); E = (E + t(E)) / sqrt(2)
	av = rnorm(n, 0, 0.4)
	Y = 0.3 + 0.8 * Xd + outer(av, av, "+") + E
	Y = (Y + t(Y)) / 2; diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	Xa = array(Xd, c(n, n, 1), dimnames = list(rownames(Y), colnames(Y), "x1"))
	fa = ame(
		Y,
		Xdyad = Xa,
		R = 0,
		symmetric = TRUE,
		family = "normal",
		burn = 400,
		nscan = 3000,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 5
	)
	fl = lame(
		list(t1 = Y),
		Xdyad = list(t1 = Xa),
		R = 0,
		symmetric = TRUE,
		family = "normal",
		burn = 400,
		nscan = 3000,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 5
	)
	sd_ratio = sd(fa$BETA[, "x1_dyad"]) / sd(fl$BETA[, "x1_dyad"])
	# double counting gives about 1.41; single counting about 1.0
	expect_gt(sd_ratio, 0.8)
	expect_lt(sd_ratio, 1.25)
	# additive variance agrees between the two entry points
	expect_equal(mean(fl$VC[, "va"]), mean(fa$VC[, "va"]), tolerance = 0.4)
})

test_that("#12 ame() stores each posterior draw once (no duplication, no crash)", {
	skip_on_cran()
	set.seed(1); n = 18
	Y = matrix(rnorm(n), n) %*% t(matrix(rnorm(n), n)) + matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	f = ame(
		Y,
		R = 1,
		family = "normal",
		burn = 60,
		nscan = 240,
		odens = 4,
		verbose = FALSE,
		gof = FALSE,
		seed = 2,
		posterior_opts = list(save_UV = TRUE, save_ab = TRUE, thin_UV = 1, thin_ab = 1)
	)
	n_eff = nrow(f$BETA)
	k = dim(f$U_samples)[3]
	expect_equal(k, n_eff)
	# consecutive pairs are distinct draws
	frac_dup = mean(sapply(
		seq(1, k - 1, 2),
		function(i) isTRUE(all.equal(f$U_samples[, , i], f$U_samples[, , i + 1]))
	))
	expect_lt(frac_dup, 0.1)

	# symmetric fits with partial posterior_opts store V_samples (U L per draw) so U L U' can be rebuilt per draw
	Ys = Y; Ys[lower.tri(Ys)] = t(Ys)[lower.tri(Ys)]
	fsym = expect_no_error(
		ame(
			Ys,
			R = 1,
			symmetric = TRUE,
			family = "normal",
			burn = 40,
			nscan = 80,
			odens = 2,
			verbose = FALSE,
			gof = FALSE,
			seed = 3,
			posterior_opts = list(save_UV = TRUE, thin_UV = 1)
		)
	)
	expect_false(is.null(fsym$V_samples))
	expect_equal(dim(fsym$V_samples), dim(fsym$U_samples))
	s1 = fsym$U_samples[, , 1] %*% t(fsym$V_samples[, , 1])
	expect_equal(s1, t(s1), tolerance = 1e-8)
})

test_that("#19 ame() samples Sab for poisson (additive variance not frozen)", {
	skip_on_cran()
	set.seed(3); n = 28
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	X = array(rnorm(n * n), c(n, n, 1))
	eta = 0.2 + 0.5 * X[, , 1] + outer(a, b, "+")
	Y = matrix(rpois(n * n, exp(pmin(eta, 5))), n, n); diag(Y) = NA
	nm = sprintf("a%02d", 1:n); rownames(Y) = colnames(Y) = nm
	dimnames(X) = list(nm, nm, "x")
	f = ame(
		Y,
		Xdyad = X,
		R = 0,
		family = "poisson",
		nscan = 1200,
		burn = 400,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 7
	)
	# a sampled Sab gives a non-constant va column
	expect_gt(sd(f$VC[, "va"]), 1e-4)
	# and va is well away from its start value of about 0.05
	expect_gt(mean(f$VC[, "va"]), 0.1)
})

test_that("ame() returns the posterior-mean multiplicative matrix (amen parity)", {
	skip_on_cran()
	set.seed(4); n = 14
	nm = paste0("a", 1:n)

	Ys = matrix(rbinom(n * n, 1, 0.3), n, n)
	Ys[lower.tri(Ys)] = t(Ys)[lower.tri(Ys)]
	diag(Ys) = NA
	dimnames(Ys) = list(nm, nm)
	fs = ame(
		Ys,
		R = 2,
		symmetric = TRUE,
		family = "binary",
		burn = 20,
		nscan = 100,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)
	expect_false(is.null(fs$ULUPM))
	expect_equal(dim(fs$ULUPM), c(n, n))
	expect_equal(rownames(fs$ULUPM), nm)
	expect_equal(fs$ULUPM, t(fs$ULUPM), tolerance = 1e-10)

	Ya = matrix(rnorm(n * n), n, n)
	diag(Ya) = NA
	dimnames(Ya) = list(nm, nm)
	fa = ame(
		Ya,
		R = 2,
		family = "normal",
		burn = 20,
		nscan = 60,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)
	expect_false(is.null(fa$UVPM))
	expect_equal(dim(fa$UVPM), c(n, n))

	# R = 0 returns a zero matrix rather than dropping the slot
	f0 = ame(
		Ys,
		R = 0,
		symmetric = TRUE,
		family = "binary",
		burn = 20,
		nscan = 60,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)
	expect_false(is.null(f0$ULUPM))
	expect_true(all(f0$ULUPM == 0))
})
