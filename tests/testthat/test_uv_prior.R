# the prior of the multiplicative effects of a unipartite network: one definition (R/uv_prior.R) that ame(), lame(),
# the exact moves and the empirical-Bayes estimators read
skip_on_cran()

test_that("the degrees of freedom are the smallest for which the prior has a finite mean", {
	for (R in 1:6) {
		df = lame:::.lame_uv_prior_df(R)
		expect_identical(df, 2L * R + 2L)
		# an inverse-Wishart prior on a 2R x 2R covariance has a finite mean from df = 2R + 2
		expect_gt(df, 2 * R + 1)
		expect_lte(df - 1L, 2 * R + 1)
	}
})

test_that("the scale is the data scale for normal and one where the link fixes it", {
	set.seed(5)
	n = 30
	Y = matrix(rnorm(n * n, 3, 2), n) + outer(rnorm(n, 0, 1.5), rnorm(n, 0, 1.5), "+"); diag(Y) = NA
	for (fam in c(
		"binary", "ordinal", "cbin", "frn", "poisson"
	)) expect_identical(lame:::.lame_uv_prior_scale(Y, fam), 1)
	# typical size of the additive effects: root mean variance of the row and column means of the centred response
	# (u_i'v_j has the units of the response, so a scale in y^2 would give y and 1000 y different priors)
	Ec = Y - mean(Y, na.rm = TRUE)
	a = rowMeans(Ec, na.rm = TRUE); b = colMeans(Ec, na.rm = TRUE)
	expect_equal(lame:::.lame_uv_prior_scale(Y, "normal"), sqrt(mean(c(var(a), var(b)))))
	expect_equal(
		lame:::.lame_uv_prior_scale(1e3 * Y, "normal"),
		1e3 * lame:::.lame_uv_prior_scale(Y, "normal"),
		tolerance = 1e-12
	)
	# a three-dimensional array (lame()) gives the same scale to a single period
	expect_equal(
		lame:::.lame_uv_prior_scale(array(Y, c(n, n, 1L)), "normal"),
		lame:::.lame_uv_prior_scale(Y, "normal"),
		tolerance = 1e-12
	)
	# a response without variation between actors falls back on the variance of the response
	Y0 = matrix(rnorm(n * n), n); Y0 = Y0 - rowMeans(Y0); Y0 = t(t(Y0) - colMeans(Y0)); diag(Y0) = NA
	expect_gt(lame:::.lame_uv_prior_scale(Y0, "normal"), 0)
	pr = lame:::.lame_uv_prior(Y, "normal", 2L)
	expect_equal(pr$kappa0, 6L); expect_equal(pr$Suv0, diag(lame:::.lame_uv_prior_scale(Y, "normal"), 4))
	expect_equal(pr$psi, 6 * pr$scale)
	expect_null(lame:::.lame_uv_prior(Y, "binary", 0L))
})

test_that("ame() and lame() fit under the same prior of the factors, for every family", {
	skip_on_cran()
	set.seed(11)
	n = 20
	mk = function(fam) {
		Z = outer(rnorm(n, 0, .5), rnorm(n, 0, .5), "+") + matrix(rnorm(n * n), n)
		Y = switch(
			fam,
			binary = 1 * (Z > 0),
			normal = 2 * Z + 1,
			poisson = matrix(rpois(n * n, exp(Z / 2)), n),
			ordinal = matrix(as.integer(cut(Z, c(-Inf, -.7, 0, .7, Inf))), n)
		)
		diag(Y) = NA; dimnames(Y) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n))); Y
	}
	for (fam in c("binary", "normal", "poisson", "ordinal")) {
		Y = mk(fam)
		fa = suppressWarnings(
			ame(Y, R = 2, family = fam, nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE, dcor = FALSE)
		)
		fl = suppressWarnings(
			lame(
				list(Y),
				R = 2,
				family = fam,
				nscan = 60,
				burn = 20,
				odens = 5,
				verbose = FALSE,
				plot = FALSE,
				dcor = FALSE
			)
		)
		want = lame:::.lame_uv_prior(Y, fam, 2L)
		expect_equal(fa$prior$kappa0, want$kappa0, info = fam)
		expect_equal(fa$prior$Suv0, want$Suv0, info = fam)
		expect_equal(fl$static_uv_prior$df, want$kappa0, info = fam)
		# the inverse-Wishart scale matrix of lame() is the one ame() builds as kappa0 * Suv0
		expect_equal(fl$static_uv_prior$scale, fa$prior$kappa0 * fa$prior$Suv0, info = fam)
	}
	# a supplied Suv0 or kappa0 is honoured by ame()
	Y = mk("binary")
	f = suppressWarnings(ame(
		Y,
		R = 1,
		family = "binary",
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		dcor = FALSE,
		prior = list(Suv0 = diag(2) * .3, kappa0 = 9)
	))
	expect_equal(f$prior$kappa0, 9); expect_equal(f$prior$Suv0, diag(2) * .3)
})

test_that("the moves and the exact densities that integrate the factor covariance out read the same prior", {
	for (R in 1:3) {
		df = lame:::.lame_uv_prior_df(R)
		u = lame:::.lame_uv_prior_unit(R)
		expect_equal(u$df, df); expect_equal(u$scale, diag(df, 2L * R))
		# the draw of the collapsed actor move under the default (lame()) is the inverse-Wishart of the shared prior
		set.seed(40 + R); n = 9L
		U = matrix(rnorm(n * R), n); V = matrix(rnorm(n * R), n); W = cbind(U, V)
		draws = replicate(3000, lame:::.lame_actor_move_suv(U, V, n, R, TRUE))
		mean_theory = (u$scale + crossprod(W)) / (n + u$df - 2 * R - 1)
		expect_lt(max(abs(apply(draws, c(1, 2), mean) - mean_theory)) / max(abs(mean_theory)), .1)
	}
	# the scale move of the static positions: the density ratio of the shared prior
	set.seed(3); n = 6L; R = 2L; W = matrix(rnorm(n * 2 * R), n); f = .7
	ev = eigen(crossprod(W), symmetric = TRUE, only.values = TRUE)$values
	df = lame:::.lame_uv_prior_df(R)
	lp = function(W) -(n + df) / 2 * as.numeric(determinant(df * diag(2 * R) + crossprod(W), logarithm = TRUE)$modulus)
	oracle = lp(sqrt(f) * W) - lp(W) + n * R * log(f)
	expect_equal(lame:::.lame_static_scale_logratio(f, 0, 0, ev, n, R, TRUE), oracle, tolerance = 1e-10)
})

test_that("the variational prior of a directed network is the shared one", {
	set.seed(21)
	n = 30
	u = rnorm(n); v = rnorm(n)
	Y = 1 * (matrix(rnorm(n * n), n) + outer(
		rnorm(n, 0, .5),
		rnorm(n, 0, .5),
		"+"
	) + 1.2 * outer(u, v) > 0.5); diag(Y) = NA
	dimnames(Y) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))
	f = suppressWarnings(ame_als(Y, R = 1, family = "binary", dcor = FALSE, verbose = FALSE))
	st = lame:::.ae_mf_object_state(f)
	cx = lame:::.ae_mf_context(st$prep, st$ab, 1L, "probit")
	# the data-scaled prior of ame() and lame() (R/uv_prior.R), read off the response the estimator fits
	want = lame:::.lame_uv_prior(st$prep$Yarr, "binary", 1L)
	expect_equal(cx$kappa0_joint, want$df)
	expect_equal(cx$psi_joint, want$df * want$scale)
	expect_equal(want$scale, lame:::.lame_binary_probit_scale(Y)$vscale, tolerance = 1e-12)
})

test_that("the link-fixed families take the data scale of their additive effects and the information multiplier", {
	set.seed(8)
	n = 30
	Z = outer(rnorm(n, 0, .6), rnorm(n, 0, .6), "+") + matrix(rnorm(n * n), n)
	Yb = 1 * (Z > 1.2); diag(Yb) = NA
	# binary: the variance of the actor means of the standardised residuals of the minority class (diagonal at zero)
	# and the multiplier 0.25 / mean(p (1 - p)) of the probit fit
	ybar = mean(Yb, na.rm = TRUE); E = (Yb - ybar) / dnorm(qnorm(ybar)); diag(E) = 0
	a = rowMeans(E, na.rm = TRUE); b = colMeans(E, na.rm = TRUE)
	vs = mean(diag(cov(cbind(a, b)))); P = pnorm(qnorm(ybar) + outer(a, b, "+")); m = 0.25 / mean(P * (1 - P))
	pr = lame:::.lame_uv_prior(Yb, "binary", 2L)
	expect_equal(pr$scale, vs, tolerance = 1e-12); expect_equal(pr$df, as.integer(round(6 * m)))
	expect_equal(pr$Suv0, diag(vs, 4), tolerance = 1e-12); expect_equal(pr$psi, pr$df * pr$scale)
	expect_gt(pr$df, 6L)
	# the class coded the other way round is the same network
	expect_equal(lame:::.lame_uv_prior(1 - Yb, "binary", 2L)[c("df", "scale")], pr[c("df", "scale")], tolerance = 1e-12)
	# a panel pools its periods; one period is the matrix
	expect_equal(
		lame:::.lame_uv_prior(array(Yb, c(n, n, 1L)), "binary", 2L)[c("df", "scale")],
		pr[c("df", "scale")],
		tolerance = 1e-12
	)
	Y2 = array(c(Yb, Yb), c(n, n, 2L))
	expect_equal(lame:::.lame_uv_prior(Y2, "binary", 2L)[c("df", "scale")], pr[c("df", "scale")], tolerance = 1e-12)
	# ordinal, cbin and frn: the probit moments of the additive prior; poisson: the log-rate moments
	Yo = matrix(as.integer(cut(Z, c(-Inf, -.7, 0, .7, Inf))), n); diag(Yo) = NA
	sc = lame:::.lame_probit_ab_scale(Yo, "ordinal")
	po = lame:::.lame_uv_prior(Yo, "ordinal", 2L)
	expect_equal(po$scale, sc$vscale); expect_equal(po$df, 6L)
	Yc = Yb; Yc[is.na(Yc)] = NA
	sc = lame:::.lame_probit_ab_scale(Yc, "cbin", rep(n, n))
	pc = lame:::.lame_uv_prior(Yc, "cbin", 1L, rep(n, n))
	expect_equal(pc$scale, sc$vscale); expect_equal(pc$df, as.integer(round(4 * sc$vdfmlt)))
	Yp = matrix(rpois(n * n, exp(Z / 2)), n); diag(Yp) = NA
	pp = lame:::.lame_uv_prior(Yp, "poisson", 2L)
	expect_equal(pp$scale, lame:::.ame_pois_effect_scale(Yp)); expect_equal(pp$df, 6L)
	# the moves of a symmetric fit keep the unit prior; a network without variation falls back on it
	expect_equal(lame:::.lame_uv_prior(Yb, "binary", 2L, directed = FALSE)[c("df", "scale")], list(df = 6L, scale = 1))
	Y0 = matrix(0, n, n); diag(Y0) = NA
	expect_equal(lame:::.lame_uv_prior(Y0, "binary", 2L)[c("df", "scale")], list(df = 6L, scale = 1))
})
