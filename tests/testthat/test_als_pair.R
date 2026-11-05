skip_on_cran()

library(lame)

# dyadic correlation of the binary probit empirical-Bayes fit of a directed network: the bivariate probit probability
# of a dyad, its expectation under the variational posterior, the additive-block information with the pairs' cross
# weights, and the correlation estimated with the other hyperparameters

pr_panel = function(n = 24, Tt = 1, R = 1, rho = 0.6, seed = 1, mu = -0.8, beta = 0.7, sd_uv = 0.8) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	U = matrix(rnorm(n * max(R, 1), 0, sd_uv), n)[, seq_len(R), drop = FALSE]
	V = matrix(rnorm(n * max(R, 1), 0, sd_uv), n)[, seq_len(R), drop = FALSE]
	O = if (R > 0) tcrossprod(U, V) else matrix(0, n, n)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, n)
		e = matrix(rnorm(n * n), n, n); z2 = matrix(rnorm(n * n), n, n)
		e[lower.tri(e)] = (rho * t(e) + sqrt(1 - rho^2) * z2)[lower.tri(e)]
		y = 1 * (mu + beta * x + outer(a, b, "+") + O + e > 0)
		diag(y) = NA
		dimnames(y) = list(rn, rn); Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(rn, rn, "x_dyad"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, n = n, Tt = Tt)
}

pr_fit = function(d, R = 1, dcor = TRUE, ...) {
	suppressWarnings(suppressMessages(
		if (d$Tt == 1L) ame_als(d$Y[[1]], d$X[[1]], R = R, family = "binary", verbose = FALSE, dcor = dcor, ...)
		else lame_als(d$Y, d$X, R = R, family = "binary", verbose = FALSE, dcor = dcor, ...)
	))
}

test_that("the derivatives of the log bivariate normal probability match finite differences", {
	a = c(-1.2, 0.3, -3, 2); b = c(0.5, -0.8, -2.5, 1.5); r = c(0.6, -0.4, 0.8, 0.2); h = 1e-5
	d = lame:::.ae_pair_logp(a, b, r)
	fd = function(f, x) (f(x + h) - f(x - h)) / (2 * h)
	expect_equal(d$La, fd(function(x) lame:::.ae_pair_logp(x, b, r)$L, a), tolerance = 1e-6)
	expect_equal(d$Lb, fd(function(x) lame:::.ae_pair_logp(a, x, r)$L, b), tolerance = 1e-6)
	expect_equal(d$Laa, fd(function(x) lame:::.ae_pair_logp(x, b, r)$La, a), tolerance = 1e-6)
	expect_equal(d$Lbb, fd(function(x) lame:::.ae_pair_logp(a, x, r)$Lb, b), tolerance = 1e-6)
	expect_equal(d$Lab, fd(function(x) lame:::.ae_pair_logp(a, x, r)$La, b), tolerance = 1e-6)
	expect_equal(d$Lr, fd(function(x) lame:::.ae_pair_logp(a, b, x)$L, r), tolerance = 1e-6)
	expect_equal(d$Lrr, fd(function(x) lame:::.ae_pair_logp(a, b, x)$Lr, r), tolerance = 1e-5)
})

test_that("the log bivariate normal cdf is accurate where pbivnorm is not", {
	# the lower tail with a negative correlation: pbivnorm returns values off by orders of magnitude there
	ab = rbind(c(-8, -8, -0.5), c(-10, -5, -0.9), c(-6, -6, -0.3), c(-8, -8, 0.8), c(-3, -8, 0.8))
	lg = lame:::.ae_pair_logphi2(ab[, 1], ab[, 2], ab[, 3])
	ref = vapply(seq_len(nrow(ab)), function(k) {
		kap = sqrt(1 - ab[k, 3]^2); a = ab[k, 1]
		lf = function(x) dnorm(x, log = TRUE) + pnorm((ab[k, 2] - ab[k, 3] * x) / kap, log.p = TRUE)
		m = max(vapply(seq(a - 14, a, length.out = 2001), lf, 0))
		m + log(integrate(function(x) exp(lf(x) - m), -Inf, a, rel.tol = 1e-12, subdivisions = 500L)$value)
	}, 0)
	expect_equal(lg, ref, tolerance = 1e-7)
	# and agrees with pbivnorm where that is accurate
	a = c(-1, 0.5, 1.5); b = c(0.2, -0.7, 1); r = c(0.5, -0.3, 0.9)
	expect_equal(lame:::.ae_pair_logphi2(a, b, r), log(pbivnorm::pbivnorm(a, b, r)), tolerance = 1e-8)
})

test_that("the expected pair problem is that of two independent cells at rho = 0 and is symmetric in the cells", {
	m1 = c(-2, -0.5, 1); m2 = c(-1, 0.4, -3); v1 = c(0.3, 0.5, 0.8); v2 = c(0.6, 0.2, 0.4)
	for (y in list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))) {
		s1 = rep(2 * y[1] - 1, 3); s2 = rep(2 * y[2] - 1, 3)
		q = lame:::.ae_pair_quad(m1, m2, v1, v2, rep(0, 3), s1, s2, 0, K = 9L)
		c1 = lame:::.ae_mf_cells(array(y[1], c(3, 1, 1)), array(m1, c(3, 1, 1)), array(v1, c(3, 1, 1)), "probit", 9L)
		c2 = lame:::.ae_mf_cells(array(y[2], c(3, 1, 1)), array(m2, c(3, 1, 1)), array(v2, c(3, 1, 1)), "probit", 9L)
		expect_equal(q$ell, as.vector(c1$ell + c2$ell), tolerance = 1e-9)
		expect_equal(q$w11, as.vector(c1$w), tolerance = 1e-9)
		expect_equal(q$w22, as.vector(c2$w), tolerance = 1e-9)
		expect_equal(q$psi1, as.vector(c1$w * (c1$z - m1)), tolerance = 1e-9)
		expect_equal(q$w12, rep(0, 3), tolerance = 1e-9)
	}
	# the two cells exchanged (the rule builds the second predictor given the first, so the nodes differ
	# and the two agree to the accuracy of the rule, not to rounding)
	q1 = lame:::.ae_pair_quad(m1, m2, v1, v2, c(0.1, 0.05, 0.1), c(1, -1, 1), c(-1, -1, 1), 0.5, K = 7L)
	q2 = lame:::.ae_pair_quad(m2, m1, v2, v1, c(0.1, 0.05, 0.1), c(-1, -1, 1), c(1, -1, 1), 0.5, K = 7L)
	expect_equal(q1$ell, q2$ell, tolerance = 1e-5)
	expect_equal(q1$psi1, q2$psi2, tolerance = 1e-5)
	expect_equal(q1$w11, q2$w22, tolerance = 2e-3)
	expect_equal(q1$w12, q2$w12, tolerance = 2e-3)
})

test_that("five nodes a side reproduce a thirteen-node rule in the quantities the fit uses", {
	set.seed(3)
	n = 60
	m1 = runif(n, -5, 0.5); m2 = runif(n, -5, 0.5); v1 = runif(n, 0.05, 1); v2 = runif(n, 0.05, 1)
	cv = runif(n, -0.4, 0.4) * sqrt(v1 * v2)
	for (y in list(c(0, 0), c(1, 0), c(1, 1))) for (rho in c(0.3, 0.8, -0.5)) {
		s1 = rep(2 * y[1] - 1, n); s2 = rep(2 * y[2] - 1, n)
		ref = lame:::.ae_pair_quad(m1, m2, v1, v2, cv, s1, s2, rho, K = 13L)
		q = lame:::.ae_pair_quad(m1, m2, v1, v2, cv, s1, s2, rho, K = 5L)
		for (nm in c("ell", "psi1", "psi2", "w11", "w22", "w12")) {
			expect_lt(sum(abs(q[[nm]] - ref[[nm]])) / sum(abs(ref[[nm]])), 1e-2)
		}
	}
})

test_that("the information of the additive block with the cross weights is the hessian of the expected log likelihood", {
	d = pr_panel(n = 12, R = 1, rho = 0.7, seed = 4)
	fit = pr_fit(d, R = 1)
	ms = lame:::.ae_mf_object_state(fit)
	cx = lame:::.ae_mf_context(ms$prep, ms$ab, ms$R, "probit", pair = list(rho = fit$eb$pair$rho, update = FALSE))
	S = ms$S; nr = cx$nr; p = cx$p
	gradfun = function(S) {
		ev = lame:::.ae_mf_evaluate(cx, S)
		psi = ev$wk$w * (ev$wk$z - ev$eta); psi[!is.finite(psi)] = 0
		c(
			sum(psi),
			vapply(seq_len(p), function(k) sum(psi * array(cx$X[, , k, ], dim(psi))), 0),
			rowSums(psi, dims = 1L),
			colSums(rowSums(psi, dims = 2L))
		)
	}
	setx = function(S, x) { S$mu = x[1]; S$beta = x[1 + seq_len(p)]; S$a = x[1 + p + seq_len(nr)]; S$b = x[1 + p + nr + seq_len(nr)]; S }
	x0 = c(S$mu, S$beta, S$a, S$b)
	h = 1e-4
	Hfd = vapply(seq_along(x0), function(k) {
		xp = x0; xm = x0; xp[k] = xp[k] + h; xm[k] = xm[k] - h
		(gradfun(setx(S, xp)) - gradfun(setx(S, xm))) / (2 * h)
	}, numeric(length(x0)))
	Hfd = (Hfd + t(Hfd)) / 2
	ev = lame:::.ae_mf_evaluate(cx, S)
	Warr = ev$wk$w
	Zst = ifelse(Warr > 0, ev$pair$G / Warr, 0)
	st = lame:::ae_bcd_stats_cpp(Warr, Zst, cx$X, cx$nr, cx$nc, cx$Tt, cx$p)
	ctx = lame:::.ae_eb_context(cx$ab, st, rowSums(Warr, dims = 2L), cx$p)
	ctx2 = lame:::.ae_pair_add_cross(ctx, lame:::.ae_pair_cross(cx, ev$pair$cw))
	idx = c(ctx$i_th, ctx$i_s)
	expect_lt(max(abs(Hfd + ctx2$H[idx, idx])), 1e-4 * max(abs(ctx2$H[idx, idx])))
	# the diagonal weights alone are not the Hessian: the cross weights matter at this correlation
	expect_gt(max(abs(Hfd + ctx$H[idx, idx])), 100 * max(abs(Hfd + ctx2$H[idx, idx])))
})

test_that("the posterior of rho is the prior without information and the peak with a lot", {
	expect_equal(lame:::.ae_pair_rho_post(0.9, 1e-6)$mean, 0, tolerance = 1e-3)
	expect_equal(lame:::.ae_pair_rho_post(0.4, 1e6)$mean, 0.4, tolerance = 1e-3)
	expect_lt(lame:::.ae_pair_rho_post(0.9, 30)$mean, 0.9)
	expect_gt(lame:::.ae_pair_rho_post(0.9, 30)$mean, 0.7)
	expect_equal(lame:::.ae_pair_rho_post(-0.3, 200)$sd, 1 / sqrt(200), tolerance = 0.05)
	# the prior is the one ame() and lame() share (R/rho_prior.R): 2 Beta(2, 2) - 1 has the variance 1/5
	expect_equal(lame:::.lame_rho_prior_power(), 1)
	expect_equal(lame:::.ae_pair_rho_post(0, 1e-9)$sd, sqrt(1 / 5), tolerance = 1e-3)
	# and its exponent is the one the helper gives: the arc-sine prior (exponent -1/2, variance 1/2) pulls a
	# correlation near the boundary less far in
	expect_equal(lame:::.ae_pair_rho_post(0, 1e-9, power = -0.5)$sd, sqrt(1 / 2), tolerance = 1e-3)
	expect_lt(lame:::.ae_pair_rho_post(0.9, 30)$mean, lame:::.ae_pair_rho_post(0.9, 30, power = -0.5)$mean)
})

test_that("the correlation of the fast fit reads its prior from the helper the samplers read", {
	d = pr_panel(n = 12, R = 0, rho = 0.5, seed = 12)
	fit_rho = function() pr_fit(d, R = 0)$eb$pair$rho
	base = fit_rho()
	# a different exponent of the helper moves the estimate (the pairs of a small network say little about rho)
	testthat::local_mocked_bindings(.lame_rho_prior_power = function() -0.5, .package = "lame")
	expect_gt(abs(fit_rho() - base), 1e-3)
})

test_that("ame_als() and lame_als() have the one prior of the samplers", {
	d = pr_panel(n = 20, R = 0, rho = 0.5, seed = 12)
	a = suppressWarnings(suppressMessages(ame_als(d$Y[[1]], d$X[[1]], R = 0, family = "binary", verbose = FALSE)))
	l = suppressWarnings(suppressMessages(lame_als(d$Y, d$X, R = 0, family = "binary", verbose = FALSE)))
	expect_identical(a$eb$pair$rho, l$eb$pair$rho)
	expect_identical(coef(a), coef(l))
	# a refit goes back to the same value
	expect_equal(ame_als_refit(l, verbose = FALSE)$eb$pair$rho, l$eb$pair$rho, tolerance = 1e-3)
})

test_that("a fit with a dyadic correlation estimates it, reports it as the rho of the model and refits with it", {
	skip_on_cran()
	d = pr_panel(n = 40, R = 1, rho = 0.7, seed = 5)
	fit = pr_fit(d, R = 1)
	expect_true(fit$converged)
	expect_false(is.null(fit$eb$pair))
	expect_equal(fit$VC[["rho"]], fit$eb$pair$rho)
	expect_equal(fit$rho_estimation$method, "variational pair likelihood")
	expect_gt(fit$eb$pair$rho, 0.4)
	expect_lt(abs(fit$eb$pair$rho - 0.7), 3 * fit$eb$pair$sd + 0.1)
	expect_true(is.finite(fit$eb$pair$sd) && fit$eb$pair$sd > 0)
	expect_true(all(is.finite(fit$eb$se)))
	# the independent-cell fit (dcor = FALSE) has no pair model
	ind = pr_fit(d, R = 1, dcor = FALSE)
	expect_null(ind$eb$pair)
	expect_null(ind$rho_estimation$method)
	# the refit starts from the pair model and returns to it
	re = ame_als_refit(fit, verbose = FALSE)
	expect_false(is.null(re$eb$pair))
	expect_equal(re$eb$pair$rho, fit$eb$pair$rho, tolerance = 1e-3)
	expect_equal(coef(re), coef(fit), tolerance = 1e-3)
})

test_that("without reciprocity the fit finds a small correlation and nearly the coefficients of the independent-cell fit", {
	skip_on_cran()
	d = pr_panel(n = 40, R = 1, rho = 0, seed = 6)
	fit = pr_fit(d, R = 1)
	ind = pr_fit(d, R = 1, dcor = FALSE)
	expect_lt(abs(fit$eb$pair$rho), 3 * fit$eb$pair$sd)
	expect_equal(coef(fit), coef(ind), tolerance = 0.1)
})

test_that("dcor is a probit, directed, empirical-Bayes feature and asks for nothing else", {
	d = pr_panel(n = 16, R = 0, rho = 0.5, seed = 7)
	# the default is on for a directed network and is not an error for the fits that have no pair model
	expect_silent(
		suppressMessages(ame_als(d$Y[[1]], d$X[[1]], R = 0, family = "binary", bias_correct = "none", verbose = FALSE))
	)
	# the logit link has no dyadic correlation of its own
	expect_error(
		ame_als(d$Y[[1]], d$X[[1]], R = 0, family = "binary", link = "logit", dcor = TRUE, verbose = FALSE),
		"probit"
	)
	ind = suppressWarnings(
		suppressMessages(ame_als(d$Y[[1]], d$X[[1]], R = 0, family = "binary", link = "logit", verbose = FALSE))
	)
	expect_null(ind$eb$pair)
	expect_error(ame_als(d$Y[[1]], d$X[[1]], R = 0, family = "binary", dcor = NA, verbose = FALSE), "TRUE")
})

test_that("a panel of several periods with missing cells is fitted and the reported correlation is the model's", {
	d = pr_panel(n = 20, Tt = 3, R = 1, rho = 0.5, seed = 8)
	set.seed(2)
	for (t in seq_along(d$Y)) { ix = which(!is.na(d$Y[[t]])); d$Y[[t]][sample(ix, round(0.2 * length(ix)))] = NA }
	fit = pr_fit(d, R = 1)
	expect_true(fit$converged)
	expect_lt(fit$eb$pair$n_pairs, 3L * 20L * 19L / 2L)
	expect_gt(fit$eb$pair$n_pairs, 0L)
	expect_true(all(is.finite(fit$eb$se)))
	expect_s3_class(fit, "lame_als")
})

test_that("the bound of a fit with a dyadic correlation is recomputed from the object", {
	d = pr_panel(n = 20, R = 1, rho = 0.5, seed = 9)
	for (R in 0:1) {
		fit = pr_fit(d, R = R)
		expect_false(is.null(fit$eb$pair))
		expect_equal(
			lame:::.ae_mf_elbo_object(fit)$elbo,
			fit$eb$variational$elbo,
			tolerance = 1e-6,
			info = paste("R =", R)
		)
	}
})

test_that("the coarse stage of the quadrature reaches the fixed point of the rule of the fit", {
	d = pr_panel(n = 30, R = 1, rho = 0.6, seed = 10)
	mode = suppressWarnings(suppressMessages(lame:::.ae_run(
		d$Y[[1]],
		d$X[[1]],
		NULL,
		NULL,
		1L,
		"binary",
		"unipartite",
		FALSE,
		200,
		1e-6,
		FALSE,
		6886,
		longitudinal = FALSE,
		call = quote(pr_mode()),
		lowrank_method = "als",
		non_normal_method = "irls",
		link = "probit",
		linear_solver = "eigen",
		multistart = "none",
		bias_correct = NULL,
		eb_posterior = "mode"
	)))
	st = lame:::.ae_eb_state(mode)
	go = function(pair) lame:::.ae_mf_fit(
		st$prep,
		st$fit,
		st$ab,
		st$R,
		FALSE,
		"probit",
		st$w,
		max_iter = 1000L,
		tol = 1e-6,
		scale = "cavity",
		pair = pair
	)
	a = go(list())
	b = go(list(coarse = FALSE))
	expect_true(a$converged && b$converged)
	expect_equal(a$S$beta, b$S$beta, tolerance = 2e-3)
	expect_equal(a$S$rho, b$S$rho, tolerance = 5e-3)
	expect_equal(a$elbo$elbo, b$elbo$elbo, tolerance = 1e-5)
	expect_identical(a$cx$pair_quad, 5L)
})

test_that("ame() and lame() with method = 'als' pass dcor to the fast estimator", {
	d = pr_panel(n = 20, R = 0, rho = 0.5, seed = 14, Tt = 2)
	a = suppressWarnings(
		suppressMessages(ame(d$Y[[1]], Xdyad = d$X[[1]], R = 0, family = "binary", method = "als", verbose = FALSE))
	)
	b = suppressWarnings(
		suppressMessages(
			ame(d$Y[[1]], Xdyad = d$X[[1]], R = 0, family = "binary", method = "als", dcor = FALSE, verbose = FALSE)
		)
	)
	expect_false(is.null(a$eb$pair))
	expect_null(b$eb$pair)
	l = suppressWarnings(
		suppressMessages(lame(d$Y, Xdyad = d$X, R = 0, family = "binary", method = "als", verbose = FALSE))
	)
	m = suppressWarnings(
		suppressMessages(
			lame(d$Y, Xdyad = d$X, R = 0, family = "binary", method = "als", dcor = FALSE, verbose = FALSE)
		)
	)
	expect_false(is.null(l$eb$pair))
	expect_null(m$eb$pair)
})

test_that("the pair quadrature in blocks of pairs is the quadrature of the whole vector", {
	set.seed(7)
	n = 12000
	m1 = rnorm(n, -2, 0.8); m2 = rnorm(n, -2, 0.8); v1 = runif(n, 0.05, 0.5); v2 = runif(n, 0.05, 0.5)
	cv = runif(n, -0.1, 0.2) * sqrt(v1 * v2)
	s1 = sample(c(-1, 1), n, TRUE, prob = c(0.8, 0.2)); s2 = sample(c(-1, 1), n, TRUE, prob = c(0.8, 0.2))
	a = lame:::.ae_pair_quad(m1, m2, v1, v2, cv, s1, s2, 0.4, K = 3L, chunk = 5000L)
	b = lame:::.ae_pair_quad(m1, m2, v1, v2, cv, s1, s2, 0.4, K = 3L, chunk = 1e9)
	expect_identical(a, b)
	expect_length(a$ell, n)
})
