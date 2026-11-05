skip_on_cran()

library(lame)

# dynamic empirical-Bayes estimator of count panels (lame_dynamic_als(), family = "poisson"; R/lame_dynamic_als_eb_pln.R,
# R/lame_dynamic_als_eb_plnpair.R, src/dynamic_eb_pln.cpp): Poisson-lognormal cells, ve as a hyperparameter, first-order posterior path means, directed dyadic correlation

# small count panel from the lame dynamic_ab model: y | z ~ Poisson(exp(z)), z = eta + e, e ~ N(0, ve), eta = mu + x'beta + a_it + b_jt
# with a level plus a zero-mean AR(1) deviation per actor; the noises of a directed dyad correlate at rho_d; `zero` actors are silent
pln_sim = function(
	n = 14,
	nb = n,
	Tt = 4,
	mode = "directed",
	seed = 1,
	mu = -0.3,
	beta = c(0.4, -0.3),
	va = 0.5,
	vb = 0.5,
	rho = 0.7,
	sigma = 0.3,
	ve = 0.3,
	rho_d = 0,
	R = 0,
	suv = 0.8,
	zero = 0L,
	drop = 0
) {
	set.seed(seed)
	sym = identical(mode, "symmetric")
	bip = identical(mode, "bipartite")
	if (!bip) nb = n
	ar1 = function(m) {
		out = matrix(0, m, Tt)
		out[, 1] = rnorm(m, 0, sigma / sqrt(1 - rho^2))
		for (t in 2:Tt) out[, t] = rho * out[, t - 1] + rnorm(m, 0, sigma)
		out
	}
	alpha = rnorm(n, 0, sqrt(va))
	gam = if (sym) alpha else rnorm(nb, 0, sqrt(vb))
	d = ar1(n)
	e = if (sym) d else ar1(nb)
	if (R > 0) {
		U0 = matrix(rnorm(n * R, 0, sqrt(suv)), n, R)
		V0 = if (sym) U0 else matrix(rnorm(nb * R, 0, sqrt(suv)), nb, R)
	}
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	cn = if (bip) paste0("b", sprintf("%02d", seq_len(nb))) else rn
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		xt = array(
			rnorm(n * nb * length(beta)),
			c(n, nb, length(beta)),
			dimnames = list(rn, cn, paste0("x", seq_along(beta)))
		)
		if (sym) for (k in seq_along(beta)) {
			m = xt[, , k]
			m[lower.tri(m)] = t(m)[lower.tri(m)]
			xt[, , k] = m
		}
		eta = mu + outer(alpha + d[, t], gam + e[, t], "+")
		for (k in seq_along(beta)) eta = eta + beta[k] * xt[, , k]
		if (R > 0) eta = eta + tcrossprod(U0, V0)
		err = matrix(rnorm(n * nb, 0, sqrt(ve)), n, nb)
		if (rho_d != 0 && !sym && !bip) {
			lo = lower.tri(err)
			err[lo] = rho_d * t(err)[lo] + sqrt(1 - rho_d^2) * matrix(rnorm(n * nb, 0, sqrt(ve)), n, nb)[lo]
		}
		if (sym) err[lower.tri(err)] = t(err)[lower.tri(err)]
		y = matrix(rpois(n * nb, exp(eta + err)), n, nb)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (zero > 0) { y[seq_len(zero), ] = 0L; if (!bip) y[, seq_len(zero)] = 0L }
		dimnames(y) = list(rn, cn)
		if (!bip) diag(y) = NA
		Y[[t]] = y
		X[[t]] = xt
	}
	if (drop > 0) for (t in seq_len(Tt)) {
		gone = sample(n, round(drop * n))
		Y[[t]][gone, ] = NA
		if (!bip) Y[[t]][, gone] = NA
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, truth = list(mu = mu, beta = beta, ve = ve, rho = rho, sigma = sigma, rho_d = rho_d))
}

pln_fit = function(d, ...) {
	suppressWarnings(lame:::lame_dynamic_als(d$Y, d$X, family = "poisson", dynamic_ab = TRUE, verbose = FALSE, ...))
}

# the cell

test_that("the count cell's derivatives are the cumulants of the posterior of the noise", {
	ns = asNamespace("lame")
	ve = 0.4
	eta = c(-1.5, -0.3, 0.4, 1.2, 2.5, -0.3, 0.8)
	y = c(0, 1, 3, 2, 12, 0, 7)
	K = ns$.lde_pln_kin(ve)
	cell = function(e) ns$.lde_pln_terms(y, e, ve, lgamma(y + 1), K)
	h = 1e-3
	up = cell(eta + h)
	dn = cell(eta - h)
	m = cell(eta)
	# the score is the derivative of the log likelihood, the information minus the derivative of the score, and the first and
	# second derivatives of the information are the next two derivatives of the log likelihood with the sign of h' = -l'''
	expect_equal(m$s, (up$l - dn$l) / (2 * h), tolerance = 1e-6)
	expect_equal(m$w, -(up$s - dn$s) / (2 * h), tolerance = 1e-6)
	expect_equal(m$hp, (up$w - dn$w) / (2 * h), tolerance = 1e-5)
	expect_equal(m$hpp, (up$hp - dn$hp) / (2 * h), tolerance = 1e-4)
	# the information is positive: the cell is log concave in the predictor
	expect_true(all(m$w > 0))
	# the log likelihood is the integral of the lognormal-Poisson probability
	pr = vapply(
		seq_along(y),
		function(i) stats::integrate(
			function(z) stats::dpois(y[i], exp(eta[i] + z)) * stats::dnorm(z, 0, sqrt(ve)),
			-12,
			12,
			rel.tol = 1e-11
		)$value,
		0
	)
	expect_equal(as.vector(exp(m$l)), pr, tolerance = 1e-8)
})

test_that("the compiled count kernels agree with their R references", {
	ns = asNamespace("lame")
	ve = 0.3
	set.seed(3)
	y = array(c(rpois(30, 2), NA, 0), c(4, 8))
	m = array(rnorm(32, 0, 1), c(4, 8))
	v = array(runif(32, 0, 0.5), c(4, 8))
	lgy = lgamma(pmax(y, 0, na.rm = TRUE) + 1)
	on.exit(options(lame.als_engine = NULL))
	a = ns$.lde_pln_terms(y, m, ve, lgy)
	options(lame.als_engine = "r")
	b = ns$.lde_pln_terms(y, m, ve, lgy)
	for (f in names(a)) expect_equal(a[[f]], b[[f]], tolerance = 1e-9, info = f)
	options(lame.als_engine = NULL)
	e1 = ns$.lde_pln_expect(y, m, v, ve, lgy)
	options(lame.als_engine = "r")
	e2 = ns$.lde_pln_expect(y, m, v, ve, lgy)
	for (f in names(e1)) expect_equal(e1[[f]], e2[[f]], tolerance = 1e-9, info = f)
})

test_that("the expected count cell is the plain cell without noise and the integral over the predictor with it", {
	ns = asNamespace("lame")
	ve = 0.35
	y = array(c(0, 1, 4, 9), c(2, 2))
	m = array(c(-0.7, 0.2, 0.9, 1.6), c(2, 2))
	lgy = lgamma(y + 1)
	z0 = ns$.lde_pln_expect(y, m, array(0, c(2, 2)), ve, lgy)
	p0 = ns$.lde_pln_terms(y, m, ve, lgy)
	expect_equal(z0$es, p0$s, tolerance = 1e-10)
	expect_equal(z0$ew, p0$w, tolerance = 1e-10)
	expect_equal(z0$eh, p0$hp, tolerance = 1e-10)
	expect_equal(z0$ehh, p0$hpp, tolerance = 1e-10)
	expect_equal(z0$el, p0$l, tolerance = 1e-10)
	expect_equal(z0$lm, p0$l, tolerance = 1e-10)
	v = array(c(0.2, 0.5, 0.05, 0.9), c(2, 2))
	zv = ns$.lde_pln_expect(y, m, v, ve, lgy)
	# the marginal log likelihood: the predictor integrated out of the likelihood is the lognormal-Poisson probability with the
	# variance of the noise plus that of the predictor, and the same by quadrature over the predictor
	pr = vapply(
		seq_along(y),
		function(i) stats::integrate(
			function(z) stats::dpois(y[i], exp(m[i] + z)) * stats::dnorm(z, 0, sqrt(ve + v[i])),
			-14,
			14,
			rel.tol = 1e-11
		)$value,
		0
	)
	expect_equal(exp(as.vector(zv$lm)), pr, tolerance = 1e-8)
	ex = vapply(
		seq_along(y),
		function(i) stats::integrate(
			function(z) as.vector(
				ns$.lde_pln_terms(rep(y[i], length(z)), m[i] + z, ve, rep(lgy[i], length(z)))$s
			) * stats::dnorm(z, 0, sqrt(v[i])),
			-8 * sqrt(v[i]),
			8 * sqrt(v[i]),
			rel.tol = 1e-9
		)$value,
		0
	)
	expect_equal(as.vector(zv$es), ex, tolerance = 1e-7)
})

test_that("the expected count cell is accurate where the predictor is very uncertain", {
	ns = asNamespace("lame")
	ve = 0.4
	# the integral over the predictor by a dense rule against the rule the expected cell chooses by the largest variance (9, 21 or 41
	# nodes): the score, the information and the second moment of the noise
	for (v in c(1, 6, 25)) for (y in c(0, 3)) {
		m = 0.5
		lgy = lgamma(y + 1)
		z = seq(-9, 9, length.out = 4001)
		wz = stats::dnorm(z) * (z[2] - z[1])
		nd = ns$.lde_pln_terms(rep(y, length(z)), m + sqrt(v) * z, ve, rep(lgy, length(z)), 28L)
		ap = ns$.lde_pln_expect(array(y, c(1, 1)), array(m, c(1, 1)), array(v, c(1, 1)), ve, array(lgy, c(1, 1)))
		expect_equal(as.vector(ap$es), sum(wz * nd$s), tolerance = 1e-3)
		expect_equal(as.vector(ap$ew), sum(wz * nd$w), tolerance = 1e-3)
		expect_equal(as.vector(ap$ee), sum(wz * nd$ee), tolerance = 1e-3)
	}
	expect_identical(length(ns$.lde_pln_outer_rule(2)$z[[1]]), 3L)
	expect_identical(length(ns$.lde_pln_outer_rule(6)$z[[1]]), 21L)
	expect_identical(length(ns$.lde_pln_outer_rule(25)$z[[1]]), 41L)
})

test_that("the second moment of the noise is the cumulant form of the cell", {
	ns = asNamespace("lame")
	ve = 0.4
	y = c(0, 2, 5, 1)
	eta = c(-1, 0.3, 1.1, 0.5)
	nd = ns$.lde_pln_terms(y, eta, ve, lgamma(y + 1))
	# E[e^2] = Var e + (E e)^2 with E e = ve * score and Var e = ve - ve^2 * w
	expect_equal(nd$ee, ve - ve^2 * nd$w + ve^2 * nd$s^2, tolerance = 1e-9)
})

# the fit

test_that("a dynamic count fit is an empirical-Bayes fit of the Poisson-lognormal model by default", {
	skip_on_cran()
	d = pln_sim(n = 16, Tt = 4, seed = 11)
	f = pln_fit(d)
	expect_s3_class(f, "lame_dynamic_als")
	expect_identical(f$family, "poisson")
	expect_true(f$converged)
	expect_identical(f$bias_correct$method, "eb")
	expect_true(f$bias_correct$default)
	# the noise variance is a fitted hyperparameter, the path persistence is estimated, and the predictor is on the latent scale
	ve = unname(f$VC["ve"])
	expect_gt(ve, 0.1)
	expect_lt(ve, 0.7)
	expect_identical(f$s2, ve)
	expect_identical(f$pois_latent_var, ve)
	expect_identical(f$eb$ve, ve)
	expect_false(isTRUE(all.equal(f$rho_ab, 0.8)))
	expect_gt(f$rho_ab, 0.2)
	expect_lt(f$rho_ab, 0.98)
	expect_equal(f$fitted[[1]], exp(f$EZ[[1]] + ve / 2))
	# coefficients near the truth, Laplace standard errors finite and the right size
	expect_lt(max(abs(f$coefficients[-1] - d$truth$beta)), 0.15)
	expect_true(all(is.finite(f$eb$se)))
	# the intercept is reported with the effects at mean zero, so its standard error is that of the identified combination
	expect_lt(f$eb$se[["intercept"]], 3 * f$eb$se[["x1_dyad"]])
	expect_equal(f$eb$se, sqrt(diag(vcov(f))), tolerance = 1e-8)
	expect_equal(f$deviance, -2 * sum(vapply(seq_along(d$Y), function(t) {
		ok = is.finite(d$Y[[t]])
		sum(stats::dpois(d$Y[[t]][ok], f$fitted[[t]][ok], log = TRUE))
	}, 0)), tolerance = 0.01)
	# the public route is the same fit
	fp = suppressWarnings(
		lame(d$Y, Xdyad = d$X, family = "poisson", dynamic_ab = TRUE, method = "als", verbose = FALSE)
	)
	expect_s3_class(fp, "lame_dynamic_als")
	expect_identical(fp$bias_correct$method, "eb")
	expect_identical(fp$coefficients, f$coefficients)
	expect_identical(fp$eb$se, f$eb$se)
})

test_that("the explicit estimators of the earlier releases are still available and the request for empirical Bayes is honoured", {
	skip_on_cran()
	d = pln_sim(n = 14, Tt = 4, seed = 12)
	fn = pln_fit(d, bias_correct = "none")
	expect_identical(fn$bias_correct$method, "none")
	# the penalised fit fixes the persistence of the paths at the prior mean
	expect_equal(fn$rho_ab, 0.8)
	fe = pln_fit(d, bias_correct = "eb")
	expect_identical(fe$bias_correct$method, "eb")
	expect_false(isTRUE(all.equal(fe$rho_ab, 0.8)))
	fj = pln_fit(d, bias_correct = "jackknife")
	expect_true(isTRUE(fj$bias_correct$applied))
	# the empirical-Bayes estimator also serves the normal family
	fnn = lame:::lame_dynamic_als(
		lapply(d$Y, function(y) y + 0.5),
		d$X,
		family = "normal",
		dynamic_ab = TRUE,
		bias_correct = "eb",
		verbose = FALSE
	)
	expect_identical(fnn$bias_correct$method, "eb")
	expect_error(
		lame:::lame_dynamic_als(d$Y, d$X, family = "poisson", dynamic_ab = TRUE, bootstrap = 5L, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
})

test_that("directed, bipartite and undirected count panels, with and without missing cells, fit and report", {
	skip_on_cran()
	for (mode in c("directed", "bipartite", "symmetric")) {
		d = pln_sim(n = 14, nb = 10, Tt = 4, mode = mode, seed = 21, drop = 0.1)
		f = suppressWarnings(lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "poisson",
			dynamic_ab = TRUE,
			verbose = FALSE,
			mode = if (mode == "bipartite") "bipartite" else "unipartite",
			symmetric = mode == "symmetric"
		))
		expect_true(f$converged, info = mode)
		expect_gt(f$VC[["ve"]], 0.05)
		expect_true(all(is.finite(f$eb$se)), info = mode)
		expect_lt(max(abs(f$coefficients[-1] - d$truth$beta)), 0.25)
	}
})

test_that("moving coefficients of a count panel are fitted with their paths", {
	d = pln_sim(n = 14, Tt = 5, seed = 13)
	f = pln_fit(d, dynamic_beta = c("intercept", "x1_dyad"), dynamic_beta_kind = "rw1")
	expect_true(f$converged)
	expect_true(f$dynamic_beta)
	expect_equal(dim(f$coef_path), c(3, 5))
	expect_true(all(is.finite(f$eb$se_path)))
	expect_equal(unname(f$sigma_beta[1]) > 0, TRUE)
})

test_that("an actor without a count in any period has finite, shrunk effects", {
	d = pln_sim(n = 16, Tt = 4, seed = 14, zero = 2L, mu = -0.2)
	f = pln_fit(d)
	expect_true(f$converged)
	expect_true(all(is.finite(f$a_dynamic)) && all(is.finite(f$b_dynamic)))
	# the silent actors sit well below the others but are bounded by the prior (the fixed-effects fit runs to minus infinity)
	zi = 1:2
	expect_lt(max(f$a_level[zi]), min(f$a_level[-zi]) + 1)
	expect_gt(min(c(f$a_level[zi], f$b_level[zi])), -9)
	expect_true(all(f$a_level[zi] < 0))
})

test_that("a panel with an all-NA period or an all-zero period fits", {
	d = pln_sim(n = 12, Tt = 4, seed = 15)
	Y0 = d$Y
	Y0[[2]][!is.na(Y0[[2]])] = 0
	f0 = pln_fit(list(Y = Y0, X = d$X))
	expect_true(f0$converged)
	expect_true(all(is.finite(f0$eb$se)))
	Yn = d$Y
	Yn[[3]][] = NA_real_
	fn = pln_fit(list(Y = Yn, X = d$X))
	expect_true(fn$converged)
	expect_true(all(is.finite(fn$eb$se)))
	expect_true(all(is.finite(fn$a_dynamic)))
})

# the correlation of a directed panel

test_that("the dyadic correlation of a directed count panel is estimated at the fitted state and enters the standard errors", {
	skip_on_cran()
	d = pln_sim(n = 24, Tt = 4, seed = 31, rho_d = 0.6, ve = 0.5, mu = 0)
	f = pln_fit(d)
	expect_true(f$converged)
	pr = f$eb$pair
	expect_false(is.null(pr))
	expect_gt(pr$n_pairs, 500)
	expect_gt(pr$rho, 0.3)
	expect_lt(pr$rho, 0.85)
	expect_equal(unname(f$VC["rho"]), pr$rho)
	# dcor = FALSE leaves the cells independent: no pair record, and the effects and the noise are those of independent cells
	fi = pln_fit(d, dcor = FALSE)
	expect_null(fi$eb$pair)
	expect_equal(fi$coefficients, f$coefficients)
	expect_equal(fi$VC[["ve"]], f$VC[["ve"]])
	# the pair information enters the standard errors, which differ from those of independent cells, and a replay of the fit keeps the switch
	expect_true(all(is.finite(fi$eb$se)))
	expect_gt(max(abs(f$eb$se / fi$eb$se - 1)), 1e-3)
	gi = suppressWarnings(ame_als_refit(fi, verbose = FALSE))
	expect_null(gi$eb$pair)
	# an undirected or a bipartite panel has no pairs of cells
	db = pln_sim(n = 12, nb = 9, Tt = 3, mode = "bipartite", seed = 32)
	fb = suppressWarnings(
		lame:::lame_dynamic_als(db$Y, db$X, family = "poisson", dynamic_ab = TRUE, mode = "bipartite", verbose = FALSE)
	)
	expect_null(fb$eb$pair)
})

test_that("an independent panel has a correlation near zero", {
	skip_on_cran()
	d = pln_sim(n = 24, Tt = 4, seed = 33, rho_d = 0, ve = 0.5, mu = 0)
	f = pln_fit(d)
	expect_lt(abs(f$eb$pair$rho), 0.3)
})

# latent structure

test_that("static factors and moving positions of a count panel are fitted", {
	skip_on_cran()
	d = pln_sim(n = 16, Tt = 4, seed = 41, R = 1, suv = 1.2, mu = -0.2)
	f1 = pln_fit(d, R = 1)
	expect_true(f1$converged)
	expect_identical(f1$R, 1L)
	expect_equal(dim(f1$U), c(16, 1))
	expect_true(all(is.finite(f1$eb$se)))
	expect_gt(cor(as.vector(f1$UVPM), as.vector(tcrossprod(f1$U, f1$V))), 0.99)
	f2 = pln_fit(d, R = 1, dynamic_uv = TRUE)
	expect_true(f2$converged)
	expect_true(f2$dynamic_uv)
	expect_equal(dim(f2$U), c(16, 1, 4))
	# the heavy-tailed transition of the positions and the switches of the additive margins
	f3 = pln_fit(d, R = 1, dynamic_uv = TRUE, dynamic_uv_kind = "t")
	expect_true(f3$converged)
	expect_identical(f3$dynamic_uv_kind, "t")
	d0 = pln_sim(n = 14, Tt = 3, seed = 42)
	fr = pln_fit(d0, rvar = FALSE)
	expect_true(fr$converged)
	expect_true(all(fr$a_dynamic == 0))
	expect_false(any(fr$b_dynamic == 0))
	expect_false(unname(fr$eb$effects[["sender"]]))
})

# reporting

test_that("the methods of a count fit work and read the latent variance", {
	d = pln_sim(n = 14, Tt = 3, seed = 51)
	f = pln_fit(d)
	expect_true(any(grepl("Poisson-lognormal", capture.output(print(f)))))
	s = capture.output(summary(f))
	expect_true(any(grepl("ve", s)))
	expect_equal(dim(confint(f)), c(3, 2))
	pr = predict(f)
	expect_equal(pr[[1]], f$fitted[[1]])
	expect_equal(predict(f, type = "link")[[1]], f$EZ[[1]])
	sm = simulate(f, nsim = 1, seed = 1)
	expect_equal(length(sm$Y), 1L)
	expect_error(logLik(f))
})

test_that("a refit from the fitted state is at its fixed point", {
	skip_on_cran()
	d = pln_sim(n = 14, Tt = 4, seed = 52)
	f = pln_fit(d)
	g = suppressWarnings(ame_als_refit(f, verbose = FALSE))
	expect_equal(g$iterations, 0L)
	expect_equal(g$coefficients, f$coefficients, tolerance = 1e-5)
	expect_equal(g$VC[["ve"]], f$VC[["ve"]], tolerance = 1e-4)
})

test_that("the warm seed of a count fit is the latent predictor it fitted", {
	set.seed(813122L)
	nr = 14L; nc = 11L
	Y = lapply(1:3, function(t) matrix(rpois(nr * nc, exp(0.2 + rnorm(nr * nc, sd = 0.7))), nr))
	f = suppressWarnings(
		lame:::lame_dynamic_als(Y, R = 0, family = "poisson", dynamic_ab = TRUE, mode = "bipartite", verbose = FALSE)
	)
	expect_identical(f$bias_correct$method, "eb")
	seed = lame:::.ae_static_seed(f)
	for (t in 1:3) {
		eta = seed$coef_path["intercept", t] + outer(seed$a_path[, t], seed$b_path[, t], "+")
		expect_equal(eta, unname(f$EZ[[t]]), tolerance = 1e-10, ignore_attr = TRUE)
		expect_equal(exp(eta + f$VC[["ve"]] / 2), unname(f$fitted[[t]]), tolerance = 1e-10, ignore_attr = TRUE)
	}
	# the penalised fit's seed adds back the half variance its predictor subtracts
	fn = suppressWarnings(lame:::lame_dynamic_als(
		Y,
		R = 0,
		family = "poisson",
		dynamic_ab = TRUE,
		mode = "bipartite",
		verbose = FALSE,
		bias_correct = "none"
	))
	sn = lame:::.ae_static_seed(fn)
	eta = sn$coef_path["intercept", 1] + outer(sn$a_path[, 1], sn$b_path[, 1], "+")
	expect_equal(exp(eta), unname(fn$fitted[[1]]), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the paths are reported at the mode shifted by minus half of the vector of the Laplace correction", {
	d = pln_sim(n = 14, Tt = 4, seed = 53, zero = 1L, mu = -1)
	ns = asNamespace("lame")
	f = pln_fit(d)
	# the reported linear predictor of a quiet actor sits below its mode: the shift is mode - u / 2 with u >= 0 where the
	# information of a cell rises with the predictor
	prep = ns$.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ctx = ns$.lde_context(prep, FALSE, "log", list(), R = 0L, dynamic_ab = TRUE, dcor = FALSE)
	fit = ns$.lde_fit(ctx, max_iter = 300L, tol = 5e-5)
	ms = ns$.lde_pln_mean_shift(ctx, fit$ev)
	expect_equal(dim(ms$x), dim(fit$st$x))
	expect_true(all(is.finite(ms$x)) && all(is.finite(ms$lev)))
	expect_gt(ms$lev[1], 0)
	# the reported levels are centred over the actors, so the check is on differences
	expect_equal(
		unname(f$a_level[1] - f$a_level[2]),
		unname((fit$st$lev[1] - 0.5 * ms$lev[1]) - (fit$st$lev[2] - 0.5 * ms$lev[2])),
		tolerance = 1e-6
	)
})

test_that("the second moment of the noise is taken about the first-order posterior mean of the predictors", {
	skip_on_cran()
	ns = asNamespace("lame")
	d = pln_sim(n = 14, Tt = 4, seed = 53, mu = -2.5)
	prep = ns$.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	for (rv in c(TRUE, FALSE)) {
		ctx = ns$.lde_context(prep, FALSE, "log", list(), R = 0L, dynamic_ab = TRUE, dcor = FALSE, rvar = rv)
		fit = ns$.lde_fit(ctx, max_iter = 300L, tol = 5e-5)
		eta = ns$.lde_eta(ctx, fit$st)
		ds = ns$.lde_pln_dshift(ctx, fit$st, eta, fit$ev$lap, fit$ev$fac)
		expect_identical(dim(ds), dim(ctx$y))
		expect_true(all(is.finite(ds[ctx$obs])))
		# it is the change of the predictors under the shift of the paths and levels by minus half of the vector of the Laplace correction
		ms = ns$.lde_pln_mean_shift(ctx, fit$ev)
		stm = fit$st
		stm$x = fit$st$x - 0.5 * ms$x
		stm$lev = fit$st$lev - 0.5 * ms$lev
		if (any(ctx$fz)) {
			stm$x[ctx$fz, ] = 0
			stm$lev[ctx$fz] = 0
		}
		expect_equal(ds, ns$.lde_eta(ctx, stm) - eta, tolerance = 1e-12)
		# the posterior mean of a sparse panel's predictors sits below its mode
		expect_lt(mean(ds[ctx$obs]), 0)
		if (!rv) {
			# a disabled sender margin has no shift: the shift of a cell is the receiver's (its effects sit at 1e-7 or less, not at zero)
			spread = apply(ds, c(2, 3), function(v) diff(range(v[is.finite(v)])))
			expect_lt(max(spread[is.finite(spread)]), 1e-5)
		} else {
			# the statistic the variance's update reads rises with the centre: the second moment of the noise about the lower predictors is larger
			s0 = ns$.lde_pln_ee_sum(ctx, eta, fit$ev$lap$vcell, fit$st$hyp$ve)
			s1 = ns$.lde_pln_ee_sum(ctx, eta + ds, fit$ev$lap$vcell, fit$st$hyp$ve)
			expect_gt(s1, s0)
			# the fitted state carries the shift of its last evaluation
			expect_equal(fit$st$dshift, ds, tolerance = 1e-3)
		}
	}
	# a binary fit has no count cell and no shift
	Yb = lapply(d$Y, function(y) { y[!is.na(y)] = as.numeric(y[!is.na(y)] > 0); y })
	prb = ns$.ae_prepare(Yb, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	cb = ns$.lde_context(prb, FALSE, "probit", list(), R = 0L, dynamic_ab = TRUE, dcor = FALSE)
	fb = ns$.lde_fit(cb, max_iter = 100L, tol = 1e-4)
	expect_null(fb$st$dshift)
})

test_that("the noise variance of a sparse panel with many periods is not underestimated", {
	skip_on_cran()
	# 90 percent zero cells and 16 periods: the data inform ve weakly and its EM image contracts at a rate near one,
	# so the fixed point is sensitive to the statistic the update reads
	ve = vapply(
		1:2,
		function(s) pln_fit(pln_sim(n = 24, Tt = 16, seed = s, mu = -3.2, ve = 0.3), max_iter = 1000L)$VC[["ve"]],
		0
	)
	expect_gt(min(ve), 0.24)
	expect_gt(mean(ve), 0.27)
	expect_lt(mean(ve), 0.45)
})

test_that("a refit of a sparse count fit is at its fixed point", {
	skip_on_cran()
	d = pln_sim(n = 20, Tt = 5, seed = 54, mu = -3)
	f = pln_fit(d, max_iter = 1000L)
	expect_true(f$converged)
	g = suppressWarnings(ame_als_refit(f, verbose = FALSE))
	expect_identical(g$iterations, 0L)
	expect_equal(g$VC[["ve"]], f$VC[["ve"]], tolerance = 1e-4)
	expect_equal(g$coefficients, f$coefficients, tolerance = 1e-5)
})

test_that("a count fit agrees with the long run of the sampler on a small panel", {
	skip_on_cran()
	skip_if(identical(Sys.getenv("LAME_SLOW_TESTS"), ""), "the sampler comparison is run with LAME_SLOW_TESTS set")
	d = pln_sim(n = 14, Tt = 4, seed = 61)
	f = pln_fit(d)
	m = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		dynamic_ab = TRUE,
		burn = 1500,
		nscan = 12000,
		odens = 10,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	B = m$BETA
	bm = colMeans(B)
	bs = apply(B, 2, sd)
	expect_lt(max(abs((f$coefficients[-1] - bm[-1]) / bs[-1])), 0.5)
	expect_equal(f$VC[["ve"]], mean(m$VC[, "ve"]), tolerance = 0.25)
	# the reported linear predictor is the posterior mean of the sampler's to a fraction of the cells' posterior sd
	dz = unlist(f$EZ) - unlist(m$EZ)
	expect_lt(sqrt(mean(dz^2, na.rm = TRUE)), 0.2)
})

# the stopping rule, the overrides and the priors

test_that("the tail of an iteration is read for its distance to the limit", {
	ns = asNamespace("lame")
	# a geometric tail: the moves sum to |move| r / (1 - r)
	for (r in c(0.3, 0.9, 0.995)) {
		h = 2 - 0.1 * cumsum(r^(0:8))
		tl = ns$.lde_tail(h)
		expect_equal(tl$rate, r, tolerance = 1e-8)
		expect_equal(tl$dist, abs(2 - 0.1 / (1 - r) - tail(h, 1)), tolerance = 1e-6)
		expect_true(tl$mono)
		expect_true(tl$stable)
	}
	# too few moves to read: unknown
	expect_identical(ns$.lde_tail(c(1, 0.9, 0.85))$dist, Inf)
	# moves that grow are no tail
	expect_gt(ns$.lde_tail(c(1, 0.99, 0.97, 0.94, 0.9))$dist, 1)
	# a sequence at its limit has distance zero
	expect_equal(ns$.lde_tail(rep(1, 6))$dist, 0)
	# the rate of a sum of two tails rises as the faster one dies out: the distance read from the largest recent ratio is a lower bound,
	# and the tail is not called stable until the ratios agree
	k = 0:12
	h = -cumsum(0.05 * 0.5^k + 0.01 * 0.97^k)
	early = ns$.lde_tail(h[1:6])
	late = ns$.lde_tail(h[1:13])
	expect_false(early$stable)
	expect_lt(early$dist, abs(-(0.05 / 0.5 + 0.01 / 0.03) - tail(h[1:6], 1)))
	expect_true(late$rate > early$rate)
	expect_false(ns$.lde_tail_done(early, 0.0071))
	expect_true(ns$.lde_tail_done(ns$.lde_tail(c(1, 1 - 1e-6, 1 - 1.5e-6, 1 - 1.75e-6, 1 - 1.875e-6)), 0.0071))
})

test_that("a count fit stops at the limit of the noise variance's iteration and not where its moves are small", {
	skip_on_cran()
	# a sparse panel whose 90 percent zero cells inform ve weakly, with an EM image that contracts at a rate near one;
	# the fit at the default tolerance is within a percent of the limit of the iteration
	d = pln_sim(n = 24, Tt = 5, seed = 71, mu = -3, ve = 0.3)
	expect_gt(mean(unlist(d$Y) == 0, na.rm = TRUE), 0.85)
	f = pln_fit(d, max_iter = 1000L)
	g = pln_fit(d, tol = 1e-9, max_iter = 3000L)
	expect_true(f$converged)
	expect_true(g$converged)
	expect_equal(f$VC[["ve"]], g$VC[["ve"]], tolerance = 0.01)
	expect_lt(abs(f$coefficients[["intercept"]] - g$coefficients[["intercept"]]), 0.02 * f$eb$se[["intercept"]])
})

test_that("an empirical-Bayes count fit can be refit as one of the penalised estimators and replays its switches", {
	skip_on_cran()
	d = pln_sim(n = 14, Tt = 4, seed = 72)
	f = pln_fit(d)
	for (bc in c("none", "jackknife")) {
		g = suppressWarnings(ame_als_refit(f, bias_correct = bc, verbose = FALSE))
		h = pln_fit(d, bias_correct = bc)
		expect_identical(g$bias_correct$method, bc)
		expect_true(g$converged)
		expect_equal(g$coefficients, h$coefficients, tolerance = 1e-4)
		expect_equal(g$VC[["ve"]], h$VC[["ve"]], tolerance = 1e-3)
	}
})

test_that("the prior of the noise variance is the sampler's and the prior of a bipartite panel reads the row actors", {
	skip_on_cran()
	d = pln_sim(n = 14, Tt = 4, seed = 73)
	f1 = pln_fit(d)
	f25 = pln_fit(d, prior = list(s20 = 25))
	expect_identical(f1$eb$ve_prior, c(nu0 = 1, s20 = 1))
	expect_identical(f25$eb$ve_prior, c(nu0 = 1, s20 = 25))
	expect_gt(f25$VC[["ve"]], f1$VC[["ve"]])
	expect_error(pln_fit(d, prior = list(s20 = -1)))
	# lame() takes the degrees of freedom of the additive effects from the number of row actors
	db = pln_sim(n = 10, nb = 60, Tt = 3, mode = "bipartite", seed = 74)
	fb = suppressWarnings(
		lame:::lame_dynamic_als(db$Y, db$X, family = "poisson", dynamic_ab = TRUE, mode = "bipartite", verbose = FALSE)
	)
	expect_equal(fb$eb$eta0, round(4 + 3 * 10 / 100))
})

test_that("a binary fit's empirical-Bayes record has no count entries", {
	d = pln_sim(n = 12, Tt = 3, seed = 75)
	Yb = lapply(d$Y, function(y) { y[!is.na(y)] = as.numeric(y[!is.na(y)] > 0); y })
	fb = suppressWarnings(lame:::lame_dynamic_als(Yb, d$X, family = "binary", dynamic_ab = TRUE, verbose = FALSE))
	expect_false(any(c("ve", "ve_prior") %in% names(fb$eb)))
	expect_true(all(c("ve", "ve_prior") %in% names(pln_fit(d)$eb)))
})
