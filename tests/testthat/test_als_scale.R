skip_on_cran()

library(lame)

# scale of the factor prior in the variational binary ALS fit: the cavity weight of a cell, the cavity covariances
# of the factor rows, the variational scale update, the undirected and bipartite priors, and the coefficient uncertainty the scale carries

sc_panel = function(
	kind = "directed",
	n = 14,
	Tt = 1,
	R = 1,
	seed = 1,
	mu = -0.5,
	beta = 0.7,
	sd_uv = 0.8
) {
	set.seed(seed)
	sym = kind == "symmetric"; bip = kind == "bipartite"
	nc = if (bip) max(6, round(0.7 * n)) else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("c%02d", seq_len(nc)) else rn
	a = rnorm(n, 0, 0.6); b = if (sym) a else rnorm(nc, 0, 0.6)
	U = matrix(rnorm(n * max(R, 1), 0, sd_uv), n)[, seq_len(R), drop = FALSE]
	V = if (sym) U else matrix(rnorm(nc * max(R, 1), 0, sd_uv), nc)[, seq_len(R), drop = FALSE]
	O = if (R > 0) tcrossprod(U, V) else matrix(0, n, nc)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * nc), n, nc); e = matrix(rnorm(n * nc), n, nc)
		if (sym) { x[lower.tri(x)] = t(x)[lower.tri(x)]; e[lower.tri(e)] = t(e)[lower.tri(e)] }
		y = 1 * (mu + beta * x + outer(a, b, "+") + O + e > 0)
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn); Y[[t]] = y
		X[[t]] = array(x, c(n, nc, 1), dimnames = list(rn, cn, "x_dyad"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, n = n, nc = nc, Tt = Tt, kind = kind)
}

sc_fit = function(d, R = 1, ...) {
	md = if (d$kind == "bipartite") "bipartite" else "unipartite"
	args = list(
		R = R,
		family = "binary",
		mode = md,
		symmetric = d$kind == "symmetric",
		verbose = FALSE,
		...
	)
	suppressWarnings(suppressMessages(
		if (d$Tt == 1L) do.call(ame_als, c(list(d$Y[[1]], Xdyad = d$X[[1]]), args))
		else do.call(lame_als, c(list(d$Y, Xdyad = d$X), args))
	))
}

# the mode of a panel and the engine state started from it
sc_engine = function(d, R = 1, ...) {
	md = if (d$kind == "bipartite") "bipartite" else "unipartite"
	Y = if (d$Tt == 1L) d$Y[[1]] else d$Y
	X = if (d$Tt == 1L) d$X[[1]] else d$X
	f = suppressWarnings(suppressMessages(lame:::.ae_run(
		Y,
		X,
		NULL,
		NULL,
		R,
		"binary",
		md,
		d$kind == "symmetric",
		200,
		1e-6,
		FALSE,
		6886,
		longitudinal = d$Tt > 1L,
		call = quote(sc_engine()),
		lowrank_method = "mm",
		non_normal_method = "irls",
		link = "probit",
		linear_solver = "eigen",
		multistart = "none",
		bias_correct = NULL,
		eb_posterior = "mode"
	)))
	st = lame:::.ae_eb_state(f)
	res = lame:::.ae_mf_fit(st$prep, st$fit, st$ab, st$R, d$kind == "symmetric", "probit", st$w, ...)
	res$prep = st$prep; res$ab = st$ab
	res
}

test_that("the cavity weight is the curvature of the log likelihood integrated over the predictor", {
	for (link in c("probit", "logit")) for (v in c(0, 0.4, 2)) for (m in c(-3, -0.5, 1.2)) for (yy in 0:1) {
		y = array(yy, c(1, 1, 1))
		got = as.numeric(lame:::.ae_mf_info_w(y, array(m, c(1, 1, 1)), array(v, c(1, 1, 1)), link))
		sg = 2 * yy - 1
		pr = function(e) if (link == "probit") pnorm(sg * e) else plogis(sg * e)
		lz = function(mm) if (v == 0) log(pr(mm)) else log(
			integrate(function(e) pr(e) * dnorm(e, mm, sqrt(v)), -Inf, Inf, rel.tol = 1e-12)$value
		)
		h = 1e-3
		want = -(lz(m + h) - 2 * lz(m) + lz(m - h)) / h^2
		expect_equal(got, want, tolerance = if (link == "probit") 2e-3 else 1e-2, info = paste(link, v, m, yy))
		expect_gte(got, 0)
		# the probit weight of a cell cannot exceed 1 / (1 + v), the information of a known unit noise
		if (link == "probit") expect_lte(got, 1 / (1 + v) + 1e-12)
	}
	# missing cells carry no weight
	y = array(c(1, NA, 0, 1), c(2, 2, 1))
	w = lame:::.ae_mf_info_w(y, array(0, c(2, 2, 1)), array(0.3, c(2, 2, 1)), "probit")
	expect_identical(w[is.na(y)], 0)
	expect_true(all(w[!is.na(y)] > 0))
})

test_that("the cavity covariance of a directed factor row is the inverse of its information form and exceeds the mean-field one", {
	d = sc_panel("directed", n = 20, R = 2, seed = 4)
	res = sc_engine(d, 2, max_iter = 600L, tol = 1e-7)
	cx = res$cx; S = res$S
	C = lame:::.ae_mf_info_sum(cx, res$ev)
	cv = lame:::.ae_mf_cavity_cov(cx, S, C)
	pb = lame:::.ae_mf_prior_blocks(S$prior, 2L)
	for (i in c(1, 7, 20)) {
		A = pb$uu + Reduce(`+`, lapply(seq_len(cx$nc), function(j) C[i, j] * tcrossprod(S[["V"]][j, ])))
		expect_equal(cv$Sig[, , i], solve(A), tolerance = 1e-10)
		# the mean-field covariance also carries the variance of the other side, so the cavity covariance is not below it
		expect_true(all(eigen(cv$Sig[, , i] - S$Sig[, , i], symmetric = TRUE, only.values = TRUE)$values > -1e-10))
	}
	for (j in c(2, 11)) {
		A = pb$vv + Reduce(`+`, lapply(seq_len(cx$nr), function(i) C[i, j] * tcrossprod(S$U[i, ])))
		expect_equal(cv$Om[, , j], solve(A), tolerance = 1e-10)
	}
	# the undirected row leaves its own cell out (the diagonal of the weights)
	ds = sc_panel("symmetric", n = 16, R = 1, seed = 5)
	rs = sc_engine(ds, 1, max_iter = 600L, tol = 1e-7)
	Cs = lame:::.ae_mf_info_sum(rs$cx, rs$ev)
	cs = lame:::.ae_mf_cavity_cov(rs$cx, rs$S, Cs)
	Fs = rs$S$F * rs$S$sg
	Cn = Cs; diag(Cn) = 0
	for (i in c(1, 9)) {
		A = 1 / rs$S$prior$tau + sum(Cn[i, ] * Fs[, 1]^2)
		expect_equal(cs$Sig[1, 1, i], 1 / A, tolerance = 1e-10)
	}
})

test_that("the directed scale is the variational update from the cavity second moments and the cavity fit is its fixed point", {
	d = sc_panel("directed", n = 24, R = 2, seed = 6)
	res = sc_engine(d, 2, max_iter = 2000L, tol = 1e-8, scale = "cavity")
	expect_identical(res$termination, "converged")
	S = res$S; cx = res$cx
	C = lame:::.ae_mf_info_sum(cx, res$ev)
	tg = lame:::.ae_mf_scale_target(cx, S, C)
	pr = S$prior
	cv = lame:::.ae_mf_cavity_cov(cx, S, C)
	S2 = S; S2$Sig = cv$Sig; S2$Om = cv$Om
	sm = lame:::.ae_mf_second_moments(S2, 2L)
	# (Psi + SE) / (df + n), reduced to the isotropic family
	M = (pr$psi + sm$SE) / (sm$nu + pr$kappa0)
	expect_equal(tg$s, mean(diag(M)), tolerance = 1e-12)
	expect_equal(tg$rho, mean(diag(M[1:2, 3:4])) / mean(diag(M)), tolerance = 1e-12)
	# options(lame.eb_scale_stat = "mean") sets the posterior mean of the covariance, (Psi + SE) / (df + n - 2 R - 1)
	op = options(lame.eb_scale_stat = "mean")
	tm = lame:::.ae_mf_scale_target(cx, S, C)
	expect_equal(tm$s, mean(diag((pr$psi + sm$SE) / (sm$nu + pr$kappa0 - 5))), tolerance = 1e-12)
	expect_gt(tm$s, tg$s)
	# a supplied prior of df <= 2 R + 1 would make that update expand; the divisor keeps one pseudo-observation beyond the rows
	S3 = S; S3$prior$kappa0 = 3
	t3 = lame:::.ae_mf_scale_target(cx, S3, C)
	expect_equal(t3$s, mean(diag((pr$psi + sm$SE) / (sm$nu + 1))), tolerance = 1e-12)
	options(op)
	expect_equal(lame:::.ae_mf_scale_target(cx, S, C)$s, tg$s, tolerance = 1e-12)
	expect_error({ options(lame.eb_scale_stat = "median"); lame:::.ae_mf_scale_target(cx, S, C) }, "mean")
	options(op)
	# the iteration stopped where the prior is the update
	expect_equal(res$scale$s, tg$s, tolerance = 2e-3)
	expect_equal(diag(pr$Suv), rep(res$scale$s, 4), tolerance = 1e-12)
	expect_identical(res$scale$rule, "cavity")
	# a fit at a fixed prior keeps it
	S1 = S; S1$prior = lame:::.ae_mf_prior_shift(S$prior, 1.7)
	r1 = lame:::.ae_mf_fit(
		res$prep,
		NULL,
		res$ab,
		2L,
		FALSE,
		"probit",
		NULL,
		max_iter = 400L,
		tol = 1e-7,
		start = S1,
		scale = "fixed"
	)
	expect_equal(r1$S$prior$Suv, S1$prior$Suv, tolerance = 1e-14)
	expect_identical(r1$scale$rule, "fixed")
})

test_that("the cavity rule lifts the scale and the product above the mean-field update where the data are weak", {
	skip_on_cran()
	d = sc_panel("directed", n = 40, R = 2, seed = 7, mu = -2.2, sd_uv = 0.8)
	cav = sc_fit(d, R = 2)
	old = options(lame.eb_scale = "plugin")
	on.exit(options(old), add = TRUE)
	plug = sc_fit(d, R = 2)
	options(old)
	expect_identical(cav$eb$variational$scale$rule, "cavity")
	expect_identical(plug$eb$variational$scale$rule, "plugin")
	expect_gt(cav$eb$variational$scale$s, mean(diag(plug$eb$suv)))
	# the variational means (the product of this weak network is withheld from the report)
	expect_gt(sum(lame:::.ae_product_restore(cav)$UVPM^2), sum(lame:::.ae_product_restore(plug)$UVPM^2))
	expect_true(cav$eb$variational$converged)
	# the reported covariances are the cavity ones, the state of the updates is kept beside them
	fp = cav$eb$factor_posterior
	expect_true(all(fp$Sig[1, 1, ] >= fp$Sig_mf[1, 1, ] - 1e-12))
	expect_identical(fp$prior$kind, "joint")
})

test_that("an unknown scale rule in the option is refused", {
	d = sc_panel("directed", n = 16, R = 1, seed = 3)
	old = options(lame.eb_scale = "mean")
	on.exit(options(old), add = TRUE)
	expect_error(sc_fit(d, R = 1), "cavity")
	expect_error(sc_fit(d, R = 0), "cavity")
	options(old)
	old = options(lame.eb_scale_uncertainty = "no")
	expect_error(sc_fit(d, R = 1), "TRUE")
})

test_that("a pure-noise panel fitted with a rank it does not have does not run away", {
	skip_on_cran()
	d = sc_panel("directed", n = 40, Tt = 2, R = 0, seed = 8, mu = -3)
	for (r in 1:3) {
		fit = sc_fit(d, R = r)
		v = fit$eb$variational
		expect_true(v$converged, info = r)
		expect_lt(v$scale$s, 1.5)
		expect_gt(v$scale$s, 0.05)
	}
})

test_that("the posterior-mean update of the scale settles on a pure-noise panel at every rank and falls back on the plug-in where its steps stall", {
	skip_on_cran()
	d = sc_panel("directed", n = 40, Tt = 2, R = 0, seed = 8, mu = -3)
	pl = sc_fit(d, R = 3)
	op = options(lame.eb_scale_stat = "mean")
	on.exit(options(op), add = TRUE)
	for (r in 1:3) {
		fit = sc_fit(d, R = r)
		v = fit$eb$variational
		expect_true(v$converged, info = paste("rank", r))
		expect_gt(v$scale$s, 0.05)
		expect_lt(v$scale$s, 1.5)
	}
	# at rank 3 the fit settles at the fixed point of the fit without the pair stage (stalled rules are tested in
	# test_als_retirement.R and test_als_multistart_selection.R)
	expect_equal(fit$eb$variational$scale$s, sc_fit(d, R = 3, dcor = FALSE)$eb$variational$scale$s, tolerance = 0.01)
	expect_gt(fit$eb$variational$scale$s, 1.2 * pl$eb$variational$scale$s)
	# where the steps do stall (forced here) the fit is repeated with the plug-in update, from the same start
	real = lame:::.ae_mf_fit
	calls = character(0)
	testthat::local_mocked_bindings(.ae_mf_fit = function(..., scale_stat = NULL) {
		calls <<- c(calls, if (is.null(scale_stat)) "default" else scale_stat)
		res = real(..., scale_stat = scale_stat)
		# only the first iteration stalls; the later calls are the two fits of the scale's uncertainty
		if (length(calls) == 1L) { res$converged = FALSE; res$termination = "no improving step" }
		res
	}, .package = "lame")
	fb = sc_fit(d, R = 3)
	expect_identical(calls[1:2], c("default", "plugin"))
	expect_equal(fb$eb$variational$scale$s, pl$eb$variational$scale$s, tolerance = 1e-8)
})

test_that("the undirected prior is the sampler's: the density of the scale is normal and its plug-in is the inverse mean precision", {
	# p(tau) = tau^-2 I(1 / (2 tau)) is a density
	pt = function(tau) exp(lame:::.ae_mf_sym_logprior(tau))
	tot = integrate(function(x) pt(exp(x)) * exp(x), -25, 25, rel.tol = 1e-10)$value
	expect_equal(tot, 1, tolerance = 1e-6)
	# I(c) is continuous across the switch to the asymptotic series and positive
	cc = c(11.99, 12.01)
	Ic = lame:::.ae_mf_sym_I(cc)
	expect_true(all(Ic > 0))
	expect_lt(abs(diff(log(Ic)) - (log(cc[2]^-2) - log(cc[1]^-2))), 1e-3)
	# the plug-in scale of a column is 1 / E[1 / tau | m]
	for (mm in c(2, 9)) for (nn in c(10, 30)) {
		w = function(tau) exp(lame:::.ae_mf_sym_logprior(tau) - 0.5 * nn * log(tau) - mm / (2 * tau))
		z = integrate(function(x) w(exp(x)) * exp(x), -25, 25, rel.tol = 1e-10)$value
		zi = integrate(function(x) w(exp(x)) * exp(x) / exp(x), -25, 25, rel.tol = 1e-10)$value
		expect_equal(lame:::.ae_mf_sym_tau_plug(mm, nn), z / zi, tolerance = 1e-3)
	}
})

test_that("the bipartite prior is the sampler's: the sender scales are the eigenvalues of a Wishart(R, I) matrix and the plug-in is the inverse mean precision of the conditional", {
	# one column: the eigenvalue of a Wishart(1, 1), a chi-square(1)
	for (mm in c(3, 12)) for (nn in c(10, 25)) {
		w = function(s) s^(-(nn + 1) / 2) * exp(-s / 2 - mm / (2 * s))
		z = integrate(function(x) w(exp(x)) * exp(x), -25, 25, rel.tol = 1e-10)$value
		zi = integrate(function(x) w(exp(x)) * exp(x) / exp(x), -25, 25, rel.tol = 1e-10)$value
		expect_equal(lame:::.ae_mf_bip_scale_plug(mm, nn, 1), z / zi, tolerance = 1e-3)
	}
	# two columns: the conditional of one is repelled from the other
	m = c(9, 2); n = 20; s = c(1.4, 0.3)
	got = lame:::.ae_mf_bip_scale_plug(m, n, s)
	for (k in 1:2) {
		l = 3 - k
		w = function(x) x^(-(n + 1) / 2) * exp(-x / 2 - m[k] / (2 * x)) * abs(x - s[l])
		cut = log(s[l])
		pieces = function(f) integrate(f, -25, cut, rel.tol = 1e-10)$value + integrate(f, cut, 25, rel.tol = 1e-10)$value
		z = pieces(function(x) w(exp(x)) * exp(x))
		zi = pieces(function(x) w(exp(x)))
		expect_equal(got[k], z / zi, tolerance = 2e-3)
	}
	# the density the conditional comes from is that of the eigenvalues of a Wishart(2, I): the moments
	# E tr W = df p = 4 and E det W = df (df - 1) = 2 of its joint density with the product of differences
	lt = seq(-14, 6, length.out = 700)
	dl = lt[2] - lt[1]
	lam = exp(lt)
	f = outer(lam, lam, function(a, b) (a * b)^(-0.5) * exp(-(a + b) / 2) * abs(a - b) * a * b)
	tot = sum(f) * dl^2
	expect_equal(sum(f * outer(lam, lam, "+")) * dl^2 / tot, 4, tolerance = 5e-3)
	expect_equal(sum(f * outer(lam, lam, "*")) * dl^2 / tot, 2, tolerance = 5e-3)
	# a column the rows say nothing about (a sum of squares far below the number of rows) is not held at the prior's scale
	expect_lt(lame:::.ae_mf_bip_scale_plug(c(40, 0.2), 60, c(2, 0.5))[2], 0.1)
})

test_that("a variational fit of every kind records its scale rule and a finite uncertainty of the scale", {
	skip_on_cran()
	for (kind in c("directed", "bipartite", "symmetric")) {
		d = sc_panel(kind, n = 22, R = 1, seed = 9)
		fit = sc_fit(d)
		sc = fit$eb$variational$scale
		expect_identical(sc$rule, "cavity", info = kind)
		expect_true(is.finite(sc$s) && sc$s > 0, info = kind)
		expect_true(is.finite(sc$sd_log) && sc$sd_log > 0, info = kind)
		expect_true(fit$eb$variational$converged, info = kind)
	}
})

test_that("the scale adds its uncertainty to the covariance of the coefficients", {
	skip_on_cran()
	# at mu = -1.5 the data-scaled factor prior (R/uv_prior.R) leaves the variational factors of this panel away from zero
	d = sc_panel("directed", n = 36, R = 2, seed = 10, mu = -1.5)
	with = sc_fit(d, R = 2)
	old = options(lame.eb_scale_uncertainty = FALSE)
	on.exit(options(old), add = TRUE)
	without = sc_fit(d, R = 2)
	options(old)
	expect_true(is.na(without$eb$variational$scale$sd_log))
	expect_gt(with$eb$variational$scale$sd_log, 0)
	expect_equal(coef(with), coef(without), tolerance = 1e-8)
	expect_true(all(with$eb$se >= without$eb$se - 1e-12))
	expect_gt(sum(with$eb$se^2), sum(without$eb$se^2))
	# the conditional covariance is the same
	expect_equal(with$eb$vcov_cond, without$eb$vcov_cond, tolerance = 1e-8)
})

test_that("the curvature of the hyperprior is the second derivative of its log density", {
	skip_on_cran()
	d = sc_panel("directed", n = 16, R = 2, seed = 11)
	res = sc_engine(d, 2, max_iter = 600L, tol = 1e-7, scale = "cavity")
	S = res$S; p = 4L
	s = mean(diag(S$prior$Suv))
	expect_equal(lame:::.ae_mf_scale_prior_curvature(res$cx, S), p / (2 * s), tolerance = 1e-12)
	db = sc_panel("bipartite", n = 16, R = 2, seed = 11)
	rb = sc_engine(db, 2, max_iter = 600L, tol = 1e-7, scale = "cavity")
	expect_equal(lame:::.ae_mf_scale_prior_curvature(rb$cx, rb$S), sum(rb$S$prior$su) / 2, tolerance = 1e-12)
	# undirected: differences of the log density in the log of a common multiplier
	ds = sc_panel("symmetric", n = 16, R = 2, seed = 11)
	rs = sc_engine(ds, 2, max_iter = 600L, tol = 1e-7, scale = "cavity")
	tau = rs$S$prior$tau
	lp = function(x) sum(lame:::.ae_mf_sym_logprior(tau * exp(x)) + log(tau * exp(x)))
	h = 1e-2
	expect_equal(
		lame:::.ae_mf_scale_prior_curvature(rs$cx, rs$S),
		-(lp(h) - 2 * lp(0) + lp(-h)) / h^2,
		tolerance = 2e-2
	)
})

test_that("a refit on the same data stays at the scale of the fit, a refit on new data re-estimates it", {
	skip_on_cran()
	d = sc_panel("directed", n = 24, Tt = 2, R = 1, seed = 12)
	fit = sc_fit(d, R = 1)
	re = ame_als_refit(fit, verbose = FALSE)
	expect_identical(re$eb$variational$scale$rule, "cavity")
	expect_equal(re$eb$variational$scale$s, fit$eb$variational$scale$s, tolerance = 1e-3)
	# the same fit on fewer cells: the refit and a cold fit agree on the scale and the coefficients
	set.seed(5)
	Y2 = simplify2array(d$Y)
	ix = which(!is.na(Y2)); Y2[sample(ix, round(0.3 * length(ix)))] = NA
	Y2l = lapply(seq_len(dim(Y2)[3]), function(t) Y2[, , t]); names(Y2l) = names(d$Y)
	re2 = suppressWarnings(suppressMessages(ame_als_refit(fit, Y_new = Y2, verbose = FALSE)))
	cold = suppressWarnings(suppressMessages(lame_als(Y2l, Xdyad = d$X, R = 1, family = "binary", verbose = FALSE)))
	expect_gt(abs(re2$eb$variational$scale$s - fit$eb$variational$scale$s), 1e-4)
	expect_equal(re2$eb$variational$scale$s, cold$eb$variational$scale$s, tolerance = 2e-2)
	expect_equal(coef(re2), coef(cold), tolerance = 2e-2)
})
