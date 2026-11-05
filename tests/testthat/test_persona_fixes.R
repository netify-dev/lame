skip_on_cran()

# input alignment, sampler correctness, reporting and support checks across families and modes

test_that("ame() realigns named-but-permuted covariates to Y", {
	set.seed(1)
	n = 12
	actors = paste0("a", sprintf("%02d", seq_len(n)))
	xv = rnorm(n)
	xr = matrix(xv, n, 1, dimnames = list(actors, "xv"))
	Y = outer(2 * xv, rep(0, n), "+") + matrix(rnorm(n * n, sd = 0.1), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = actors
	perm = sample(n)

	fit_aligned = ame(
		Y,
		Xrow = xr,
		R = 0,
		family = "normal",
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	fit_perm = suppressMessages(ame(
		Y,
		Xrow = xr[perm, , drop = FALSE],
		R = 0,
		family = "normal",
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	# same seed and realigned design give identical draws
	expect_equal(coef(fit_perm), coef(fit_aligned))

	# a named covariate whose actors differ from Y's aborts
	xr_bad = xr
	rownames(xr_bad) = paste0("z", seq_len(n))
	expect_error(
		ame(
			Y,
			Xrow = xr_bad,
			R = 0,
			family = "normal",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		),
		"actor names"
	)
})

test_that("ame() realigns permuted Xdyad by dimnames", {
	set.seed(2)
	n = 10
	actors = paste0("a", sprintf("%02d", seq_len(n)))
	Xd = array(
		rnorm(n * n),
		c(n, n, 1),
		dimnames = list(actors, actors, "x")
	)
	Y = 3 * Xd[, , 1] + matrix(rnorm(n * n, sd = 0.1), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = actors
	perm = sample(n)
	Xd_perm = Xd[perm, perm, , drop = FALSE]

	fit = suppressMessages(ame(
		Y,
		Xdyad = Xd_perm,
		R = 0,
		family = "normal",
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	expect_gt(unname(coef(fit)["x_dyad"]), 2)
})

test_that("saved log_lik maps every coefficient to its own design slice", {
	set.seed(3)
	n = 10
	Xd = array(rnorm(n * n), c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	Y = 1.5 + 2 * Xd[, , 1] + matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	fit = ame(
		Y,
		Xdyad = Xd,
		R = 0,
		family = "normal",
		rvar = FALSE,
		cvar = FALSE,
		dcor = FALSE,
		nscan = 20,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		save_log_lik = TRUE
	)

	obs_idx = which(!is.na(Y))
	s = nrow(fit$BETA)
	eta = matrix(0, n, n)
	for (k in seq_along(fit$BETA[s, ])) {
		eta = eta + fit$BETA[s, k] * fit$X[, , k]
	}
	ve = fit$VC[s, ncol(fit$VC)]
	ll_manual = dnorm(Y[obs_idx], mean = eta[obs_idx], sd = sqrt(ve), log = TRUE)
	expect_equal(unname(fit$log_lik[s, ]), unname(ll_manual), tolerance = 1e-8)
})

test_that("exactly symmetric directed normal data reject a proved improper model", {
	set.seed(4)
	n = 8
	A = matrix(rnorm(n * n), n, n)
	Ys = (A + t(A)) / 2
	diag(Ys) = NA
	dimnames(Ys) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))
	before = .Random.seed
	expect_message(
		expect_error(
			ame(
				Ys,
				R = 0,
				nscan = 60,
				burn = 15,
				odens = 5,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE
			),
			class = "lame_improper_normal_reciprocity"
		),
		"exactly symmetric"
	)
	expect_identical(.Random.seed, before)
	fit_sym = ame(
		Ys,
		symmetric = TRUE,
		R = 0,
		nscan = 6,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	expect_s3_class(fit_sym, "ame")
	expect_message({fit_fixed = ame(
		Ys,
		dcor = FALSE,
		R = 0,
		nscan = 6,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)}, "exactly symmetric")
	expect_s3_class(fit_fixed, "ame")
})

test_that("lame() rejects duplicated actor names and rectangular unipartite Y", {
	Y1 = matrix(rnorm(36), 6, 6)
	rownames(Y1) = colnames(Y1) = c("a1", "a1", "a2", "a3", "a4", "a5")
	expect_error(
		lame(
			list(Y1, Y1),
			family = "normal",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		"duplicated actor name"
	)

	Yr = matrix(
		rnorm(54),
		6,
		9,
		dimnames = list(paste0("r", 1:6), paste0("c", 1:9))
	)
	expect_error(
		lame(
			list(Yr, Yr),
			family = "normal",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		"bipartite"
	)
})

test_that("lame() validates R like ame() does", {
	Y1 = matrix(rnorm(36), 6, 6)
	rownames(Y1) = colnames(Y1) = paste0("a", 1:6)
	expect_error(
		lame(
			list(Y1, Y1),
			R = -2,
			family = "normal",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		"non-negative integer"
	)
	expect_error(
		lame(
			list(Y1, Y1),
			R = 1.7,
			family = "normal",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		"non-negative integer"
	)
})

test_that("abbreviated model.name argument cannot silently change the fit", {
	Y = matrix(rbinom(100, 1, 0.3), 10, 10)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:10)
	# lame: an abbreviated spelling aborts with guidance
	expect_error(
		lame(
			list(Y, Y),
			model = "bin",
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		"partially matched"
	)
})

test_that("constant discrete Y aborts instead of crashing in chol", {
	Y1 = matrix(1, 8, 8)
	diag(Y1) = NA
	expect_error(
		ame(
			Y1,
			family = "binary",
			R = 0,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		),
		"constant"
	)
})

test_that("compute_mcmc_diagnostics flags identical chains", {
	Y = matrix(rnorm(64), 8, 8)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:8)
	fit1 = ame(
		Y,
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	fit2 = ame(
		Y,
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	expect_warning(
		compute_mcmc_diagnostics(list(fit1, fit2)),
		"identical"
	)
})

test_that("procrustes_align does not partial-match fit$GOF into G", {
	skip_on_cran()
	skip_if_not_installed("abind")
	data(YX_bin_list)
	Y_bin = lapply(YX_bin_list$Y[1:3], function(y) {
		y = 1 * (y > 0)
		y
	})
	fit = suppressWarnings(lame(
		Y_bin,
		R = 1,
		family = "binary",
		dynamic_uv = TRUE,
		nscan = 10,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = TRUE
	))
	pa = procrustes_align(fit)
	expect_false(identical(pa$G, fit$GOF))
})

test_that("netify symmetric attribute propagates through the bridge", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	df = expand.grid(
		i = paste0("a", 1:8),
		j = paste0("a", 1:8),
		stringsAsFactors = FALSE
	)
	df = df[df$i < df$j, ]
	set.seed(5)
	df$y = rbinom(nrow(df), 1, 0.4)
	net = netify::netify(
		df,
		actor1 = "i",
		actor2 = "j",
		weight = "y",
		missing_to_zero = FALSE
	)

	fit = suppressMessages(ame(
		net,
		R = 0,
		nscan = 40,
		burn = 10,
		odens = 2,
		gof = FALSE,
		verbose = FALSE,
		plot = FALSE
	))
	expect_true(fit$symmetric)
	expect_error(
		suppressMessages(ame(
			net,
			symmetric = FALSE,
			R = 0,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		)),
		"conflicts"
	)
})

test_that("lame() accepts a list of data.frame Xrow covariates", {
	set.seed(6)
	n = 6
	actors = paste0("a", seq_len(n))
	Y1 = matrix(
		rnorm(n * n),
		n,
		n,
		dimnames = list(actors, actors)
	)
	diag(Y1) = NA
	xr = data.frame(xv = rnorm(n))
	fit = lame(
		list(Y1, Y1),
		Xrow = list(xr, xr),
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(fit, "lame")
})

test_that("4-D amen-layout Xdyad is accepted as the per-period list it stands for", {
	set.seed(6886)
	n = 5
	actors = paste0("a", seq_len(n))
	Yl = lapply(1:3, function(t) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(actors, actors))
		diag(y) = NA
		y
	})
	Xd4 = array(rnorm(n * n * 2 * 3), c(n, n, 2, 3), dimnames = list(actors, actors, c("x1", "x2"), NULL))
	Xl = lapply(1:3, function(t) Xd4[, , , t])
	f4 = suppressMessages(lame(
		Yl,
		Xdyad = Xd4,
		family = "normal",
		nscan = 10,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	))
	fl = suppressMessages(lame(
		Yl,
		Xdyad = Xl,
		family = "normal",
		nscan = 10,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	))
	expect_s3_class(f4, "lame")
	expect_equal(f4$BETA, fl$BETA)
})

# statistical and computational checks

test_that("lame() coefficients are scale-equivariant (no ridge attenuation)", {
	set.seed(11)
	n = 15; Tt = 3
	x1 = array(rnorm(n * n * Tt), c(n, n, Tt))
	mk = function(s) lapply(1:Tt, function(t) {
		y = (1 + 1 * x1[, , t] + rnorm(n * n)) * s; diag(y) = NA
		rownames(y) = colnames(y) = sprintf("a%02d", 1:n); y })
	Xl = lapply(
		1:Tt,
		function(t) array(
			x1[, , t],
			c(n, n, 1),
			dimnames = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n), "x1")
		)
	)
	f1 = lame(
		mk(1),
		Xdyad = Xl,
		family = "normal",
		R = 0,
		nscan = 150,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 3
	)
	f2 = lame(
		mk(1e3),
		Xdyad = Xl,
		family = "normal",
		R = 0,
		nscan = 150,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 3
	)
	b1 = coef(f1)["x1_dyad"]; b2 = coef(f2)["x1_dyad"] / 1e3
	expect_gt(b1, 0.7)
	expect_lt(abs(b2 - b1) / abs(b1), 0.15)
})

test_that("ame() inference is scale-invariant (g-prior + adaptive Sab0)", {
	skip_on_cran()
	set.seed(12)
	n = 15
	x = matrix(rnorm(n * n), n, n)
	Y = 0.5 + x + rnorm(n * n); diag(Y) = NA
	Xd = array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	f1 = ame(
		Y,
		Xdyad = Xd,
		family = "normal",
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 3
	)
	f2 = ame(
		Y * 1e4,
		Xdyad = Xd,
		family = "normal",
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 3
	)
	expect_lt(max(abs(apply(f2$BETA, 2, mean) / 1e4 - apply(f1$BETA, 2, mean))), 0.05)
	# rho is not manufactured at odd scales
	f3 = ame(
		Y * 1e-2,
		Xdyad = Xd,
		family = "normal",
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 3
	)
	expect_lt(abs(mean(f3$VC[, "rho"]) - mean(f1$VC[, "rho"])), 0.15)
})

test_that("dynamic_uv latent trajectories are not frozen (sigma_uv IG rate)", {
	skip_on_cran()
	set.seed(13)
	Tt = 4; nn = 15
	U0 = matrix(rnorm(nn), nn, 1)
	Yl = list()
	for (t in 1:Tt) {
		U0 = 0.9 * U0 + matrix(rnorm(nn, 0, 0.5), nn, 1)
		z = -0.5 + U0 %*% t(U0) + matrix(rnorm(nn * nn), nn, nn)
		y = 1 * (z > 0); diag(y) = NA
		rownames(y) = colnames(y) = sprintf("a%02d", 1:nn)
		Yl[[t]] = y
	}
	fd = lame(
		Yl,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 3
	)
	expect_gt(mean(abs(fd$U[, , Tt] - fd$U[, , 1])), 0.03)
})

test_that("poisson posterior predictive includes the overdispersion layer", {
	skip_on_cran()
	set.seed(14)
	n = 20
	x = matrix(rnorm(n * n), n, n)
	z = 1 + 0.5 * x + rnorm(n * n)
	Yp = matrix(rpois(n * n, exp(z)), n, n); diag(Yp) = NA
	fp = ame(
		Yp,
		Xdyad = array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x")),
		family = "poisson",
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		gof = TRUE,
		plot = FALSE,
		seed = 3,
		save_log_lik = TRUE
	)
	# YPM is the marginal mean, not biased low by exp(-s2/2)
	expect_gt(mean(fp$YPM, na.rm = TRUE) / mean(Yp, na.rm = TRUE), 0.85)
	expect_true(all(is.finite(fp$log_lik)))
})

test_that("posterior_opts save_UV_draws stores dynamic latent trajectory draws", {
	skip_on_cran()
	set.seed(15)
	n = 12; Tt = 3
	U0 = matrix(rnorm(n), n, 1)
	Yl = list()
	for (t in 1:Tt) {
		U0 = 0.9 * U0 + matrix(rnorm(n, 0, 0.5), n, 1)
		z = -0.5 + U0 %*% t(U0) + matrix(rnorm(n * n), n, n)
		y = 1 * (z > 0); diag(y) = NA
		rownames(y) = colnames(y) = sprintf("a%02d", 1:n)
		Yl[[t]] = y
	}
	fd = lame(
		Yl,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 3,
		posterior_opts = list(save_UV_draws = TRUE)
	)
	# [actor, dim, period, draw] with nscan/odens = 8 stored draws
	expect_equal(dim(fd$U_draws), c(12, 1, 3, 8))
	expect_equal(dim(fd$V_draws), c(12, 1, 3, 8))
	expect_true(all(is.finite(fd$U_draws[, , , 8])))
	expect_true(all(is.finite(fd$V_draws[, , , 8])))
	expect_equal(rownames(fd$U_draws), sprintf("a%02d", 1:n))
	# with storage off the draw slots are absent
	f0 = lame(
		Yl,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 3,
		posterior_opts = posterior_options()
	)
	expect_null(f0$U_draws)
	expect_null(f0$V_draws)
})

# exact truncated-normal z augmentation

test_that("bipartite binary ame() recovers the probit slope (no clip attenuation)", {
	# the z update draws the truncated-normal full conditional, so the probit slope is not attenuated
	set.seed(21)
	nA = 30; nB = 25
	X = array(rnorm(nA * nB), c(nA, nB, 1))
	Z0 = -0.5 + 1.0 * X[, , 1] + matrix(rnorm(nA * nB), nA, nB)
	Y = 1 * (Z0 > 0)
	rownames(Y) = paste0("r", sprintf("%02d", 1:nA))
	colnames(Y) = paste0("c", sprintf("%02d", 1:nB))
	dimnames(X) = list(rownames(Y), colnames(Y), "x1")
	glm_slope = coef(glm(c(Y) ~ c(X[, , 1]), binomial("probit")))[2]
	fit = ame(
		Y,
		Xdyad = X,
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		family = "binary",
		mode = "bipartite",
		burn = 200,
		nscan = 600,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)
	slope = colMeans(fit$BETA)[2]
	expect_lt(abs(slope / glm_slope - 1), 0.15)
})

test_that("rZ_bin_bip_batch_cpp never draws on the wrong side of the truncation", {
	skip_on_cran()
	# cells with extreme ez get z on the correct side of zero
	set.seed(5)
	ez_vals = c(-40, -10, -8, -6, 0, 6, 8, 10, 40)
	nE = length(ez_vals)
	EZ = array(rep(ez_vals, 2), c(nE, 2, 1))
	Y = array(c(rep(0, nE), rep(1, nE)), c(nE, 2, 1))
	Z = array(0, c(nE, 2, 1))
	for (r in 1:200) {
		z = rZ_bin_bip_batch_cpp(Z, EZ, Y)
		expect_true(all(is.finite(z)))
		expect_true(all(z[, 1, 1] <= 0)) # y = 0 cells
		expect_true(all(z[, 2, 1] >= 0)) # y = 1 cells
	}
})

test_that("rZ_bin_fused_cpp has no frozen sign-violating cells at extreme ez", {
	set.seed(3)
	n = 10
	EZ = matrix(rnorm(n * n, 0, 1.5), n, n)
	EZ[1, 2] = 8; EZ[2, 1] = -8
	Y = 1 * (EZ + matrix(rnorm(n * n), n, n) > 0)
	Y[1, 2] = 0; Y[2, 1] = 1
	Z = EZ + matrix(rnorm(n * n), n, n)
	k12 = k21 = numeric(300)
	for (s in 1:300) {
		Z = rZ_bin_fused_cpp(Z, EZ, 0.55, Y)
		k12[s] = Z[1, 2]; k21[s] = Z[2, 1]
	}
	# both cells move while keeping their signs
	expect_true(all(k12 <= 0))
	expect_true(all(k21 >= 0))
	expect_gt(sd(k12), 0)
	expect_gt(sd(k21), 0)
})

test_that("rtnorm_interval_logp stays inside deep-tail intervals", {
	set.seed(7)
	ez = rep(c(-12, 0, 12), each = 2000)
	z1 = lame:::rtnorm_interval_logp(ez, 0, Inf)
	expect_true(all(is.finite(z1)) && all(z1 >= 0))
	z0 = lame:::rtnorm_interval_logp(ez, -Inf, 0)
	expect_true(all(is.finite(z0)) && all(z0 <= 0))
	z2 = lame:::rtnorm_interval_logp(ez, 1, 2)
	expect_true(all(is.finite(z2)) && all(z2 >= 1) && all(z2 <= 2))
	# moment check at ez = 0 on (0, Inf): mean 0.7979, var 1 - 2/pi
	z = lame:::rtnorm_interval_logp(rep(0, 20000), 0, Inf)
	expect_lt(abs(mean(z) - 0.7979), 0.02)
	expect_lt(abs(var(z) - (1 - 2 / pi)), 0.02)
})

test_that("rZ_nrm_fc imputes mixed pairs from the reciprocal-conditional law", {
	set.seed(2)
	n = 4
	Y = matrix(0, n, n); diag(Y) = NA; Y[1, 2] = NA; Y[2, 1] = 3
	Z = Y; Z[is.na(Z)] = 0; Z[2, 1] = 3
	d = replicate(3000, rZ_nrm_fc(Z, matrix(0, n, n), 0.9, 1, Y)[1, 2])
	# full conditional: N(0.9 * 3, 1 - 0.81) = N(2.7, 0.436^2)
	expect_lt(abs(mean(d) - 2.7), 0.05)
	expect_lt(abs(sd(d) - sqrt(0.19)), 0.04)
	# rectangular (bipartite) input runs on the t(A) path
	Yb = matrix(rnorm(12), 3, 4); Yb[1, 2] = NA
	expect_no_error(rZ_nrm_fc(Yb, matrix(0, 3, 4), 0.5, 1, Yb))
})

test_that("rZ_nrm_batch_cpp both-missing pairs keep within-pair correlation rho", {
	set.seed(4)
	n = 6; rho = 0.8
	Y = array(rnorm(n * n), c(n, n, 1))
	Y[1, 2, 1] = NA; Y[2, 1, 1] = NA
	EZ = array(0, c(n, n, 1))
	Z = Y; Z[is.na(Z)] = 0
	d = replicate(4000, {
		z = rZ_nrm_batch_cpp(Z, EZ, rho, 1, Y)$Z
		c(z[1, 2, 1], z[2, 1, 1])
	})
	# the stationary within-pair correlation is rho
	expect_gt(cor(d[1, ], d[2, ]), rho - 0.05)
	expect_lt(abs(mean(d[1, ])), 0.06)
	expect_lt(abs(sd(d[1, ]) - 1), 0.05)
})

test_that("rs2_rep_fc gamma shape and rate match the C++ kernel target", {
	skip_on_cran()
	# shape (m + 1) / 2 and an explicit fixed positive prior rate in both kernels
	set.seed(1)
	n = 3; E = array(rnorm(n * n), c(n, n, 1)); rho = 0.3
	rhoMat = solve(matrix(c(1, rho, rho, 1), 2, 2)); H = mhalf(rhoMat)
	Et = E[, , 1]
	EM = cbind(Et[upper.tri(Et)], t(Et)[upper.tri(Et)]) %*% H
	ss = sum(EM^2); m = length(EM); s20 = 0.7; rate = (ss + s20) / 2
	set.seed(99)
	d = replicate(4e4, rs2_rep_fc_cpp(E, rhoMat, s20))
	expect_lt(abs(mean(1 / d) * rate - (m + 1) / 2), 0.05)
	# the exported R reference targets the same fixed-prior posterior
	set.seed(100)
	dR = replicate(4e4, rs2_rep_fc(E, rho, s20))
	expect_lt(abs(mean(1 / dR) * rate - (m + 1) / 2), 0.05)
	# all-zero residuals do not yield NaN
	E0 = array(0, c(n, n, 1))
	expect_true(is.finite(rs2_rep_fc(E0, rho)))
})

test_that("poisson marginal log-lik uses adaptive GH accurate in the tails (defect 23)", {
	skip_on_cran()
	# the adaptive Gauss-Hermite rule centred at the per-observation Laplace mode matches integrate() to < 1e-4,
	# including tail cells with |log y - eta| >> sqrt(s2)
	exact = function(y, eta, s2) {
		z0 = log(y + 0.5); sdp = 1 / sqrt(y + 1 / s2)
		f = function(z) exp(
			dpois(y, exp(z), log = TRUE) +
				dnorm(z, eta, sqrt(s2), log = TRUE) + 200
		)
		log(
			integrate(
				f,
				z0 - 40 * sdp,
				z0 + 40 * sdp,
				rel.tol = 1e-13
			)$value
		) - 200
	}
	# the internal poisson branch through a synthetic fit whose cell eta is the design slice and ve the single VC column
	grid = expand.grid(
		y = c(0, 1, 5, 25, 100),
		eta = c(-2, 0, 2, 4),
		s2 = c(0.25, 1, 3)
	)
	max_err = 0
	for (s2 in unique(grid$s2)) {
		sub = grid[grid$s2 == s2, ]
		ny = nrow(sub)
		Yg = matrix(sub$y, ny, 1)
		fit = list(
			BETA = matrix(1, 1, 1, dimnames = list(NULL, "x")),
			VC = matrix(s2, 1, 1, dimnames = list(NULL, "ve")),
			X = array(
				sub$eta,
				c(ny, 1, 1),
				dimnames = list(NULL, NULL, "x")
			)
		)
		ll = lame:::.ame_compute_log_lik_post_hoc(
			fit,
			Y = Yg,
			family = "poisson"
		)
		ex = mapply(exact, sub$y, sub$eta, s2)
		max_err = max(max_err, max(abs(as.numeric(ll[1, ]) - ex)))
	}
	expect_lt(max_err, 1e-4)

	# divergent-draw regime (eta >> log y); a non-square Y keeps the observation from being dropped as a diagonal cell
	fit_div = list(
		BETA = matrix(1, 1, 1, dimnames = list(NULL, "x")),
		VC = matrix(1, 1, 1, dimnames = list(NULL, "ve")),
		X = array(800, c(1, 2, 1), dimnames = list(NULL, NULL, "x"))
	)
	ll_div = lame:::.ame_compute_log_lik_post_hoc(
		fit_div,
		Y = matrix(3, 1, 2),
		family = "poisson"
	)
	expect_equal(as.numeric(ll_div[1, 1]), -315459.72, tolerance = 1e-2)

	# WAIC on a poisson panel stays finite and sane
	set.seed(14)
	n = 20
	x = matrix(rnorm(n * n), n, n)
	z = 1 + 0.5 * x + rnorm(n * n)
	Yp = matrix(rpois(n * n, exp(z)), n, n); diag(Yp) = NA
	fp = ame(
		Yp,
		Xdyad = array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x")),
		family = "poisson",
		R = 0,
		dcor = FALSE,
		nscan = 100,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 3,
		save_log_lik = TRUE
	)
	expect_true(all(is.finite(fp$log_lik)))
	w = suppressWarnings(waic(fp))
	expect_true(is.finite(w$estimates["waic", "Estimate"]))
})

test_that("#8 deep-tail truncated draws stay finite for the ordinal kernel", {
	# an (lb, ub) interval more than 8 sd from ez stays finite through the log.p-scale inverse cdf
	set.seed(11)
	n = 6
	# an independent random start with an inverted middle window is refused
	Yo = matrix(sample(0:2, n * n, replace = TRUE), n, n); diag(Yo) = NA
	invalid = matrix(rnorm(n * n), n, n)
	before = .Random.seed
	expect_error(
		rZ_ord_fc(invalid, matrix(15, n, n), rho = 0, Y = Yo),
		class = "lame_truncation_support_error"
	)
	# a supported category-ordered state gives a genuine deep-tail conditional
	assign(".Random.seed", before, envir = .GlobalEnv)
	Zo = Yo
	diag(Zo) = 0
	Zo = rZ_ord_fc(Zo, matrix(15, n, n), rho = 0, Y = Yo)
	expect_true(all(is.finite(Zo[row(Zo) != col(Zo)])))
	levels = sort(unique(Yo[is.finite(Yo)]))
	for (k in seq_len(length(levels) - 1L)) {
		expect_lt(
			max(Zo[!is.na(Yo) & Yo == levels[k]]),
			min(Zo[!is.na(Yo) & Yo == levels[k + 1L]])
		)
	}
})

test_that("#4 observed_ghk poisson pointwise log-lik == overdispersed marginal", {
	# the GHK poisson branch uses the adaptive Gauss-Hermite marginal of the default log_lik path, not the
	# equidispersed plug-in dpois(round(y), exp(eta))
	set.seed(3)
	n_obs = 40
	y_obs = rpois(n_obs, 3)
	ez_obs = rnorm(n_obs, 1, 0.7)
	obs_idx = cbind(i = rep(1:8, each = 5), j = rep(1:5, 8), t = rep(1L, n_obs))
	s2 = 0.8
	ghk = lame:::.pointwise_loglik_observed_ghk(
		y_obs,
		ez_obs,
		rnorm(n_obs),
		obs_idx,
		family = "poisson",
		s2 = s2
	)
	marg = lame:::.lame_pois_marginal_loglik(round(y_obs), ez_obs, s2)
	expect_equal(ghk, marg)
	# the rank dispatcher sends poisson through the same path
	ghk_rank = lame:::.pointwise_loglik_observed_ghk_rank(
		y_obs,
		ez_obs,
		rnorm(n_obs),
		obs_idx,
		family = "poisson",
		s2 = s2
	)
	expect_equal(ghk_rank, marg)
})

test_that("#9/#13 summary(lame)$coefficients aliases $beta", {
	set.seed(913)
	n = 8; Tt = 3
	Yl = lapply(seq_len(Tt), function(t) { m = matrix(rnorm(n * n), n, n); diag(m) = NA; m })
	Xl = lapply(
		seq_len(Tt),
		function(t)
			array(rnorm(n * n), c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	)
	lf = lame(
		Yl,
		Xdyad = Xl,
		family = "normal",
		R = 0,
		nscan = 30,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	sl = summary(lf)
	expect_false(is.null(sl$coefficients))
	expect_identical(sl$coefficients, sl$beta)
	# `$coef` resolves by partial matching (no exact `coef` element)
	expect_identical(sl$coef, sl$beta)
})

test_that("#5 dynamic_beta drift percentage does not explode on a zero-crossing", {
	# per-period means -0.3, -0.1, 0.1, 0.3: time average near 0, drift range near 0.6, so a percentage of the mean is undefined
	niter = 40; Tt = 4
	BETA = array(
		0,
		c(niter, 2, Tt),
		dimnames = list(NULL, c("intercept", "x_dyad"), NULL)
	)
	traj = c(-0.3, -0.1, 0.1, 0.3)
	set.seed(5)
	for (t in seq_len(Tt)) {
		BETA[, 1, t] = rnorm(niter, 1.0, 0.02) # bounded away from 0
		BETA[, 2, t] = rnorm(niter, traj[t], 0.01) # crosses zero
	}
	fake = structure(
		list(
			BETA = BETA,
			VC = matrix(1, niter, 1, dimnames = list(NULL, "ve")),
			n_time = Tt,
			family = "normal",
			mode = "unipartite",
			dynamic_beta = TRUE,
			beta_dynamic_mask = c(FALSE, TRUE),
			RHO_BETA = matrix(0.8, niter, 1, dimnames = list(NULL, "x_dyad")),
			call = quote(lame())
		),
		class = "lame"
	)
	sf = summary(fake)
	# a zero-crossing coefficient has an NA percentage
	expect_true(is.na(sf$beta_dynamic_drift_pct[["x_dyad"]]))
	# nothing explodes anywhere
	expect_true(all(
		is.na(sf$beta_dynamic_drift_pct) |
			sf$beta_dynamic_drift_pct < 1000
	))
	# the absolute drift range is reported (about 0.6)
	expect_equal(unname(sf$beta_dynamic_drift[["x_dyad"]]), 0.6, tolerance = 0.05)
	# print runs and never emits a runaway percentage
	out = capture.output(print(sf))
	expect_false(any(grepl("[0-9]{4,}", out)))
})

test_that("#7 mismatched Xdyad actor dimension gives a clear error, not base R crash", {
	set.seed(7)
	n = 20
	Yb = matrix(rbinom(n * n, 1, 0.3), n, n); diag(Yb) = NA
	Xbad = array(rnorm(19 * 19), c(19, 19, 1), dimnames = list(NULL, NULL, "x"))
	# binary fits report the dimension error up front, before the separation scan
	expect_error(
		ame(
			Yb,
			Xdyad = Xbad,
			family = "binary",
			R = 0,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		),
		"do not match"
	)
	Yn = matrix(rnorm(n * n), n, n); diag(Yn) = NA
	expect_error(
		ame(
			Yn,
			Xdyad = Xbad,
			family = "normal",
			R = 0,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE
		),
		"do not match"
	)
})

test_that("#6 predict.ame(newdata) carries dimnames and realigns by actor name", {
	set.seed(6)
	n = 10
	actors = paste0("a", sprintf("%02d", seq_len(n)))
	Xd = array(rnorm(n * n), c(n, n, 1), dimnames = list(actors, actors, "x"))
	Y = 2 * Xd[, , 1] + matrix(rnorm(n * n, sd = 0.2), n, n); diag(Y) = NA
	rownames(Y) = colnames(Y) = actors
	fit = ame(
		Y,
		Xdyad = Xd,
		R = 1,
		family = "normal",
		nscan = 50,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	nd = array(rnorm(n * n), c(n, n, 1), dimnames = list(actors, actors, "x"))
	p_ord = predict(fit, newdata = nd, type = "link")
	# output carries the fit's actor dimnames
	expect_false(is.null(rownames(p_ord)))
	expect_identical(rownames(p_ord), rownames(fit$YPM))
	expect_identical(colnames(p_ord), colnames(fit$YPM))
	# reordered but named newdata gives the same prediction
	perm = sample(n)
	nd_perm = nd[perm, perm, , drop = FALSE]
	p_perm = suppressMessages(predict(fit, newdata = nd_perm, type = "link"))
	expect_equal(p_perm, p_ord, tolerance = 1e-10)
})

test_that("#12 identical chains suppress the green 'converged' success tick", {
	set.seed(12)
	n = 8
	Yn = matrix(rnorm(n * n), n, n); diag(Yn) = NA
	Xn = array(rnorm(n * n), c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	# same seed gives byte-identical BETA across chains
	mk = function() ame(
		Yn,
		Xdyad = Xn,
		family = "normal",
		R = 0,
		seed = 42,
		nscan = 30,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	c1 = mk(); c2 = mk()
	msgs = character(0)
	out = character(0)
	withCallingHandlers(
		{ out = capture.output(compute_mcmc_diagnostics(list(c1, c2))) },
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
		warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
	)
	allmsg = paste(c(msgs, out), collapse = "\n")
	expect_true(grepl("identical", allmsg, ignore.case = TRUE))
	expect_false(grepl("All parameters converged", allmsg))
})

# adversarial inputs

test_that("dynamic_beta + latent factors + additive effects works (a/b via joint sampler)", {
	skip_on_cran()
	set.seed(404); n = 18; Tp = 5
	bt = seq(-0.8, 0.8, length.out = Tp)
	Yl = Xl = vector("list", Tp)
	au = rnorm(n, 0, 0.5); bu = rnorm(n, 0, 0.5); uu = rnorm(n, 0, 0.6); vu = rnorm(n, 0, 0.6)
	for (t in 1:Tp) {
		xt = matrix(rnorm(n * n), n, n)
		et = 0.2 + bt[t] * xt + outer(au, bu, "+") + outer(uu, vu) +
			matrix(rnorm(n * n, 0, 0.6), n, n)
		diag(et) = NA; rownames(et) = colnames(et) = paste0("a", 1:n)
		Yl[[t]] = et
		Xl[[t]] = array(xt, c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	}
	# directed normal with R > 0 and additive effects: the joint rbeta_ab sampler (whitened unit prior floor)
	# keeps the additive variance bounded
	fa = lame(
		Yl,
		Xdyad = Xl,
		family = "normal",
		R = 1,
		dynamic_beta = "dyad",
		nscan = 400,
		burn = 150,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(fa, "lame")
	expect_true(all(is.finite(fa$VC[, "va"])) && max(fa$VC[, "va"]) < 40)
	# a symmetric latent space is supported too, with a symmetric covariate
	Ys = lapply(Yl, function(m) { s = (m + t(m)) / 2; diag(s) = NA; s })
	Xs = lapply(Xl, function(x) { x[, , 1] = (x[, , 1] + t(x[, , 1])) / 2; x })
	fs = lame(
		Ys,
		Xdyad = Xs,
		family = "normal",
		R = 1,
		dynamic_beta = "dyad",
		symmetric = TRUE,
		nscan = 400,
		burn = 150,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(fs, "lame")
	expect_true(max(fs$VC[, "va"]) < 40)
	# ordinal is supported in the same combination; a noisy ordinal panel recovers cleanly
	Yo = lapply(Yl, function(m) {
		z = matrix(
			cut(as.numeric(m), c(-Inf, -0.5, 0.5, Inf), labels = FALSE),
			nrow(m),
			ncol(m)
		)
		diag(z) = NA; dimnames(z) = dimnames(m); z
	})
	fo = lame(
		Yo,
		Xdyad = Xl,
		family = "ordinal",
		R = 1,
		dynamic_beta = "dyad",
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(fo, "lame")
	expect_true(all(is.finite(fo$VC[, "va"])))
	# R = 0 dynamic_beta and additive-off fits both run
	f0 = lame(
		Yl,
		Xdyad = Xl,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		nscan = 50,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(f0, "lame")
	f1 = lame(
		Yl,
		Xdyad = Xl,
		family = "normal",
		R = 1,
		dynamic_beta = "dyad",
		rvar = FALSE,
		cvar = FALSE,
		nscan = 50,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE
	)
	expect_s3_class(f1, "lame")
})

test_that("dynamic_uv R=1 post-processing (simulate/gof/uv_plot) does not crash", {
	skip_on_cran()
	set.seed(2); nn = 15; U0 = matrix(rnorm(nn), nn, 1); Yl = list()
	for (t in 1:4) {
		U0 = 0.9 * U0 + matrix(rnorm(nn, 0, 0.5), nn, 1)
		z = -0.4 + U0 %*% t(U0) + matrix(rnorm(nn * nn), nn, nn)
		y = 1 * (z > 0); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", 1:nn); Yl[[t]] = y
	}
	f = suppressWarnings(lame(
		Yl,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		nscan = 100,
		burn = 40,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 7,
		posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)
	))
	expect_no_error(suppressWarnings(suppressMessages(simulate(f, nsim = 1))))
	expect_no_error(suppressWarnings(suppressMessages(gof(f))))
	expect_no_error(suppressWarnings(suppressMessages(uv_plot(f))))
})

test_that("plot() works on a dynamic_beta fit (3-D BETA)", {
	skip_on_cran()
	set.seed(1); n = 12
	Xl = lapply(
		1:4,
		function(t) array(
			rnorm(n * n),
			c(n, n, 1),
			dimnames = list(NULL, NULL, "x")
		)
	)
	bt = seq(-0.8, 0.8, length.out = 4)
	Yl = lapply(1:4, function(t) {
		z = bt[t] * Xl[[t]][, , 1] + matrix(rnorm(n * n), n, n)
		y = 1 * (z > 0); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", 1:n); y
	})
	f = suppressWarnings(lame(
		Yl,
		Xdyad = Xl,
		family = "binary",
		R = 0,
		dynamic_beta = "dyad",
		nscan = 100,
		burn = 40,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 7
	))
	pf = tempfile(fileext = ".pdf"); pdf(pf)
	expect_no_error(suppressWarnings(suppressMessages(plot(f))))
	dev.off(); unlink(pf)
})

test_that("bipartite save_UV also stores G_samples for reconstruction", {
	set.seed(6); nA = 10; nB = 8
	Yb = matrix(rnorm(nA * nB), nA, nB)
	rownames(Yb) = paste0("r", 1:nA); colnames(Yb) = paste0("c", 1:nB)
	f = suppressWarnings(ame(
		Yb,
		R = 2,
		mode = "bipartite",
		family = "normal",
		nscan = 100,
		burn = 40,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		seed = 7,
		posterior_opts = list(save_UV = TRUE, thin_UV = 1)
	))
	expect_false(is.null(f$G_samples))
	expect_equal(dim(f$G_samples)[1:2], c(2L, 2L))
	expect_equal(dim(f$G_samples)[3], dim(f$U_samples)[3])
})

test_that("symmetric intercept-less fits keep coefficient names (dynamic_beta selectable)", {
	skip_on_cran()
	set.seed(11); n = 12; Tp = 3
	X1 = matrix(rnorm(n * n), n, n); Xs = (X1 + t(X1)) / 2
	Xd = array(Xs, c(n, n, 1), dimnames = list(NULL, NULL, "x"))
	Yl = lapply(seq_len(Tp), function(t) {
		z = 0.5 * Xs + matrix(rnorm(n * n), n, n); z = (z + t(z)) / 2
		y = matrix(cut(as.numeric(z), c(-Inf, -0.5, 0.5, Inf), labels = FALSE), n, n)
		diag(y) = NA; rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	names(Yl) = paste0("t", seq_len(Tp))
	# symmetric ordinal defaults to intercept = FALSE and keeps the coefficient name (bnames[-numeric(0)] edge case)
	f = suppressWarnings(lame(
		Yl,
		Xdyad = lapply(seq_len(Tp), function(t) Xd),
		family = "ordinal",
		symmetric = TRUE,
		R = 0,
		nscan = 40,
		burn = 15,
		odens = 2,
		verbose = FALSE,
		plot = FALSE
	))
	expect_identical(colnames(f$BETA), "x_dyad")
	# block-based dynamic_beta selection works on that fit spec
	fd = suppressWarnings(lame(
		Yl,
		Xdyad = lapply(seq_len(Tp), function(t) Xd),
		family = "ordinal",
		symmetric = TRUE,
		R = 0,
		dynamic_beta = "dyad",
		nscan = 40,
		burn = 15,
		odens = 2,
		verbose = FALSE,
		plot = FALSE
	))
	expect_true(is.matrix(coef(fd)) && "x_dyad" %in% rownames(coef(fd)))
})

test_that("removed families give the curated unknown-family error", {
	n = 10
	Y = matrix(rbinom(n * n, 1, 0.3), n, n); diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	for (fam in c("tobit", "rrl")) {
		expect_error(
			ame(Y, family = fam, nscan = 10, burn = 5),
			"not[[:space:]]+a[[:space:]]+recognised[[:space:]]+family"
		)
		expect_error(
			lame(list(t1 = Y), family = fam, nscan = 10, burn = 5),
			"not[[:space:]]+a[[:space:]]+recognised[[:space:]]+family"
		)
	}
})

test_that("bipartite R_row/R_col fills UVPM like R does", {
	skip_on_cran()
	set.seed(7); nA = 12; nB = 9
	U0 = matrix(rnorm(nA * 2, 0, 0.8), nA, 2); V0 = matrix(rnorm(nB * 2, 0, 0.8), nB, 2)
	Yl = list(); Tt = 3
	for (t in seq_len(Tt)) {
		eta = -0.4 + U0 %*% t(V0)
		Y = matrix(rbinom(nA * nB, 1, pnorm(eta)), nA, nB)
		rownames(Y) = paste0("r", sprintf("%02d", seq_len(nA)))
		colnames(Y) = paste0("c", sprintf("%02d", seq_len(nB)))
		Yl[[paste0("t", t)]] = Y
	}
	f1 = lame(
		Y = Yl,
		family = "binary",
		mode = "bipartite",
		R_row = 2,
		R_col = 2,
		nscan = 60,
		burn = 20,
		odens = 5,
		seed = 11,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)
	f2 = lame(
		Y = Yl,
		family = "binary",
		mode = "bipartite",
		R = 2,
		nscan = 60,
		burn = 20,
		odens = 5,
		seed = 11,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)
	expect_gt(max(abs(unlist(f1$UVPM))), 0)
	expect_equal(f1$UVPM, f2$UVPM)
	expect_identical(f1$R, f2$R)
	# a one-sided rank drops the bilinear term and says so
	expect_warning(
		lame(
			Y = Yl,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			nscan = 20,
			burn = 5,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		"both[[:space:]]+ranks"
	)
})

test_that("bipartite unequal-rank fits survive the full post-processing surface", {
	skip_on_cran()
	set.seed(9); nA = 12L; nB = 9L; Tt = 3
	U0 = matrix(rnorm(nA * 2, 0, 0.8), nA, 2); V0 = matrix(rnorm(nB * 2, 0, 0.8), nB, 2)
	Yl = list()
	for (t in seq_len(Tt)) {
		Y = matrix(rbinom(nA * nB, 1, pnorm(-0.3 + U0 %*% t(V0))), nA, nB)
		rownames(Y) = paste0("r", sprintf("%02d", seq_len(nA)))
		colnames(Y) = paste0("c", sprintf("%02d", seq_len(nB)))
		Yl[[paste0("t", t)]] = Y
	}
	f21 = lame(
		Y = Yl,
		family = "binary",
		mode = "bipartite",
		R_row = 2,
		R_col = 1,
		nscan = 40,
		burn = 10,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)
	# rank-1 side: predict(h > 0) keeps the latent slice (static fit, so the constant-path advisory warning is expected)
	pr = suppressWarnings(predict(f21, h = 2, n_draws = 5, seed = 1))
	expect_length(pr, 2)
	expect_identical(dim(pr[[1]]), c(nA, nB))
	# forecasts carry actor names
	expect_identical(rownames(pr[[1]]), rownames(Yl[[1]]))
	expect_identical(colnames(pr[[1]]), colnames(Yl[[1]]))
	# plot() does not index a second latent column that V lacks
	plot_file = tempfile(fileext = ".pdf")
	grDevices::pdf(plot_file)
	tryCatch(expect_no_error(plot(f21)), finally = {
		grDevices::dev.off()
		unlink(plot_file)
	})
	# reconstruct_UVPM returns the stored posterior mean for lame fits
	expect_identical(reconstruct_UVPM(f21), f21$UVPM)
	expect_true(is.matrix(f21$UVPM) && max(abs(f21$UVPM)) > 0)
	# unequal ranks raise no vector-recycling warnings
	f23 = expect_no_warning(
		lame(
			Y = Yl,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 3,
			nscan = 40,
			burn = 10,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_true(max(abs(unlist(f23$UVPM))) > 0)
})

test_that("uv_plot trajectory handles unequal bipartite ranks", {
	skip_on_cran()
	set.seed(10); nA = 10; nB = 8; Tt = 3
	Yl = list()
	for (t in seq_len(Tt)) {
		Y = matrix(rbinom(nA * nB, 1, 0.35), nA, nB)
		rownames(Y) = paste0("r", sprintf("%02d", seq_len(nA)))
		colnames(Y) = paste0("c", sprintf("%02d", seq_len(nB)))
		Yl[[paste0("t", t)]] = Y
	}
	f31 = suppressWarnings(lame(
		Y = Yl,
		family = "binary",
		mode = "bipartite",
		R_row = 3,
		R_col = 1,
		dynamic_uv = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		seed = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_s3_class(suppressMessages(uv_plot(f31, plot_type = "trajectory")), "ggplot")
	expect_s3_class(suppressMessages(uv_plot(f31, plot_type = "faceted")), "ggplot")
	# snapshot slices keep actor labels
	p = suppressMessages(uv_plot(f31, plot_type = "snapshot"))
	expect_s3_class(p, "ggplot")
})

test_that("forecasting works from a changing-composition fit", {
	skip_on_cran()
	set.seed(11); actors = paste0("a", sprintf("%02d", 1:12))
	Yl = list()
	for (t in 1:4) {
		# period 1 holds only 9 of the 12 actors
		act_t = if (t == 1) actors[1:9] else actors
		n_t = length(act_t)
		Y = matrix(
			rbinom(n_t * n_t, 1, 0.3),
			n_t,
			n_t,
			dimnames = list(act_t, act_t)
		)
		diag(Y) = NA
		Yl[[paste0("t", t)]] = Y
	}
	fit = suppressWarnings(lame(
		Y = Yl,
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		seed = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	pr = predict(fit, h = 1, n_draws = 5, seed = 2)
	expect_identical(dim(pr[[1]]), c(12L, 12L))
	expect_identical(rownames(pr[[1]]), sort(actors))
})

test_that("bipartite dynamic_ab recovers full-strength coefficients", {
	skip_on_cran()
	set.seed(501); nA = 12; nB = 9; Tt = 3
	ra = sprintf("r%02d", 1:nA); ca = sprintf("c%02d", 1:nB)
	X = lapply(
		1:Tt,
		function(t) array(
			rnorm(nA * nB),
			c(nA, nB, 1),
			dimnames = list(ra, ca, "x1")
		)
	)
	av = rnorm(nA, 0, 0.3); bv = rnorm(nB, 0, 0.3)
	Y = lapply(1:Tt, function(t) {
		m = -0.3 + 0.5 * X[[t]][, , 1] + outer(av, rep(1, nB)) +
			outer(rep(1, nA), bv) + matrix(rnorm(nA * nB), nA, nB)
		dimnames(m) = list(ra, ca); m
	})
	fd = lame(
		Y,
		Xdyad = X,
		family = "normal",
		mode = "bipartite",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 600,
		burn = 150,
		odens = 5,
		seed = 501,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	f0 = lame(
		Y,
		Xdyad = X,
		family = "normal",
		mode = "bipartite",
		R = 0,
		nscan = 600,
		burn = 150,
		odens = 5,
		seed = 501,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	# the dynamic_ab chain centres near the static estimate, not at half of it
	expect_gt(mean(fd$BETA[, "x1_dyad"]) / mean(f0$BETA[, "x1_dyad"]), 0.8)
	expect_gt(acf(fd$BETA[, "x1_dyad"], plot = FALSE)$acf[2], -0.3)
})

test_that("simulate() supports cbin fits", {
	skip_on_cran()
	set.seed(1); n = 14
	X = matrix(rnorm(n * n), n, n)
	Z = -0.3 + 0.5 * X + matrix(rnorm(n * n), n, n)
	Y = matrix(0, n, n)
	for (i in 1:n) {
		zi = Z[i, ]; zi[i] = -Inf; pos = which(zi > 0)
		m = min(length(pos), 4)
		if (m > 0) Y[i, pos[order(zi[pos], decreasing = TRUE)][1:m]] = 1
	}
	diag(Y) = NA; rownames(Y) = colnames(Y) = paste0("a", 1:n)
	Xa = array(X, c(n, n, 1), dimnames = list(rownames(Y), colnames(Y), "x"))
	fit = ame(
		Y,
		Xdyad = Xa,
		family = "cbin",
		odmax = rep(4, n),
		nscan = 100,
		burn = 30,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	s = simulate(fit, nsim = 2)
	expect_length(s$Y, 2)
	expect_true(all(unlist(s$Y) %in% c(0, 1, NA)))
	fl = lame(
		list(t1 = Y, t2 = Y),
		Xdyad = list(Xa, Xa),
		family = "cbin",
		odmax = rep(4, n),
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)
	sl = simulate(fl, nsim = 2)
	expect_length(sl$Y, 2)
})
