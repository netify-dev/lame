skip_on_cran()

library(lame)

# empirical-Bayes estimator of the binary ALS fit (bias_correct = "eb", the
# default of a binary IRLS fit)

# a small directed binary panel with known effects
eb_panel = function(
	n = 12,
	Tt = 1,
	seed = 1,
	mu = -0.5,
	beta = 0.7,
	sd_ab = 0.6,
	R = 0,
	sd_uv = 0.6,
	sym = FALSE,
	quiet_actor = NULL
) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, sd_ab); b = if (sym) a else rnorm(n, 0, sd_ab)
	U = matrix(rnorm(n * max(R, 1), 0, sd_uv), n, max(R, 1))[, seq_len(R), drop = FALSE]
	V = if (sym) U else matrix(rnorm(n * max(R, 1), 0, sd_uv), n, max(R, 1))[, seq_len(R), drop = FALSE]
	O = if (R > 0) tcrossprod(U, V) else matrix(0, n, n)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, n)
		e = matrix(rnorm(n * n), n, n)
		if (sym) { x[lower.tri(x)] = t(x)[lower.tri(x)]; e[lower.tri(e)] = t(e)[lower.tri(e)] }
		y = 1 * (mu + beta * x + outer(a, b, "+") + O + e > 0)
		if (!is.null(quiet_actor)) y[quiet_actor, ] = 0
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(nm, nm, "x_dyad"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, a = a, b = b, beta = beta, mu = mu, n = n, Tt = Tt, O = O)
}

# the posterior mode with the Laplace correction, the stage a variational fit starts from
mode_fit = function(d, ..., R = 0, symmetric = FALSE, link = "probit") {
	Y = if (d$Tt == 1L) d$Y[[1]] else d$Y
	X = if (d$Tt == 1L) d$X[[1]] else d$X
	suppressWarnings(suppressMessages(lame:::.ae_run(
		Y,
		X,
		NULL,
		NULL,
		R,
		"binary",
		"unipartite",
		symmetric,
		200,
		1e-6,
		FALSE,
		6886,
		longitudinal = d$Tt > 1L,
		call = quote(mode_fit()),
		lowrank_method = "mm",
		non_normal_method = "irls",
		link = link,
		linear_solver = "eigen",
		multistart = "none",
		bias_correct = NULL,
		eb_posterior = "mode",
		...
	)))
}

eb_fit = function(d, ..., R = 0, symmetric = FALSE) {
	args = list(R = R, family = "binary", verbose = FALSE, symmetric = symmetric, ...)
	suppressWarnings(suppressMessages(
		if (d$Tt == 1L) do.call(ame_als, c(list(d$Y[[1]], Xdyad = d$X[[1]]), args))
		else do.call(lame_als, c(list(d$Y, Xdyad = d$X), args))
	))
}

test_that("a binary irls fit is empirical Bayes by default and on request", {
	d = eb_panel()
	fit = eb_fit(d)
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(fit$bias_correct$default)
	expect_true(isTRUE(fit$converged))
	expect_true(all(is.finite(fit$coefficients)))
	expect_identical(
		names(fit$eb),
		c(
			"sab", "Sab0", "eta0", "gamma", "sab_change",
			"sig2_u", "sig2_v", "mu_int", "a_raw", "b_raw",
			"log_marginal", "vcov", "vcov_cond", "se",
			"se_method", "vcov_error", "estimator", "variational", "laplace", "pair"
		)
	)
	expect_true(is.finite(fit$eb$log_marginal))
	# the estimator that produced the values is on the record
	expect_identical(fit$eb$estimator, "variational")
	expect_true(fit$eb$variational$converged)
	expect_identical(mode_fit(d)$eb$estimator, "laplace")

	asked = eb_fit(d, bias_correct = "eb")
	expect_false(asked$bias_correct$default)
	expect_identical(asked$bias_correct$requested, "eb")
	expect_identical(coef(asked), coef(fit))

	# the options that are not the default keep the fixed-effects fit
	none = eb_fit(d, bias_correct = "none")
	expect_identical(none$bias_correct$method, "none")
	expect_null(none$eb)
	expect_false(isTRUE(all.equal(coef(none), coef(fit), tolerance = 1e-3)))
	expect_error(eb_fit(d, bias_correct = "jackknife"), class = "lame_jk_split")
})

test_that("longitudinal and bipartite binary fits are empirical Bayes by default", {
	d = eb_panel(n = 10, Tt = 3, seed = 2)
	fit = eb_fit(d)
	expect_s3_class(fit, "lame_als")
	expect_identical(fit$bias_correct$method, "eb")

	set.seed(3)
	nr = 10; nc = 7
	rn = sprintf("r%02d", seq_len(nr)); cn = sprintf("c%02d", seq_len(nc))
	a = rnorm(nr, 0, .6); b = rnorm(nc, 0, .6)
	x = matrix(rnorm(nr * nc), nr, nc)
	y = 1 * (-0.3 + 0.7 * x + outer(a, b, "+") + matrix(rnorm(nr * nc), nr, nc) > 0)
	dimnames(y) = list(rn, cn)
	xa = array(x, c(nr, nc, 1), dimnames = list(rn, cn, "x_dyad"))
	bf = suppressWarnings(suppressMessages(
		ame_als(y, Xdyad = xa, family = "binary", mode = "bipartite", verbose = FALSE)
	))
	expect_identical(bf$bias_correct$method, "eb")
	expect_length(bf$a, nr)
	expect_length(bf$b, nc)
	expect_equal(dim(bf$eb$vcov), c(2L, 2L))
	# independent sender and receiver variances
	expect_equal(bf$eb$sab[1, 2], 0)
})

test_that("the empirical-Bayes estimator is refused where it is not defined", {
	set.seed(4)
	# (a normal response has the estimator, test_als_normal_eb.R)
	Yp = matrix(rpois(144, 2), 12, 12); diag(Yp) = NA
	expect_error(
		ame_als(Yp, family = "poisson", non_normal_method = "transform", bias_correct = "eb", verbose = FALSE),
		class = "lame_als_eb"
	)
	d = eb_panel()
	expect_error(
		eb_fit(d, bias_correct = "eb", non_normal_method = "transform"),
		class = "lame_als_eb"
	)
	# the moving positions of a logit link have no empirical-Bayes estimator on the dynamic route (the
	# sampler's model is probit)
	dd = eb_panel(n = 10, Tt = 3, seed = 5)
	expect_error(
		suppressWarnings(suppressMessages(
			lame:::lame_dynamic_als(
				dd$Y,
				dd$X,
				R = 1,
				family = "binary",
				dynamic_uv = TRUE,
				dynamic_ab = TRUE,
				link = "logit",
				bias_correct = "eb",
				verbose = FALSE
			)
		)),
		class = "lame_als_eb"
	)
	fit = eb_fit(d)
	# a transformed count fit has no empirical-Bayes refit either
	expect_error(
		ame_als_refit(
			suppressWarnings(
				ame_als(Yp, family = "poisson", non_normal_method = "transform", verbose = FALSE)
			),
			bias_correct = "eb"
		),
		class = "lame_als_eb"
	)
	expect_identical(ame_als_refit(fit)$bias_correct$method, "eb")
})

test_that("an actor without a tie has a finite effect", {
	d = eb_panel(n = 14, seed = 6, quiet_actor = 3)
	fit = eb_fit(d)
	expect_null(fit$not_estimated)
	expect_true(is.finite(fit$a[3]))
	expect_lt(fit$a[3], 0)
	expect_true(all(is.finite(fit$a)) && all(is.finite(fit$b)))
	expect_true(all(is.finite(sqrt(diag(vcov(fit))))))

	# the fixed-effects fit has no finite value for that sender
	none = eb_fit(d, bias_correct = "none")
	expect_true(is.na(none$a[3]))
	expect_true("a03" %in% none$not_estimated$actor)
})

test_that("the covariance of a variational fit is the inverse penalised information at its expected weights", {
	d = eb_panel(n = 8, seed = 7)
	fit = eb_fit(d, link = "logit")
	n = d$n
	y = d$Y[[1]]; x = d$X[[1]][, , 1]
	cells = which(!is.na(y), arr.ind = TRUE)
	Z = cbind(
		1,
		x[cells],
		sapply(seq_len(n), function(k) as.numeric(cells[, 1] == k)),
		sapply(seq_len(n), function(k) as.numeric(cells[, 2] == k))
	)
	P = solve(fit$eb$sab)
	pen = matrix(0, 2 + 2 * n, 2 + 2 * n)
	ia = 2 + seq_len(n); ib = 2 + n + seq_len(n)
	pen[ia, ia] = P[1, 1] * diag(n); pen[ib, ib] = P[2, 2] * diag(n)
	pen[ia, ib] = P[1, 2] * diag(n); pen[ib, ia] = P[2, 1] * diag(n)
	w = fit$obs_weights[, , 1][cells]
	V = solve(crossprod(Z, w * Z) + pen)[1:2, 1:2]
	expect_equal(unname(vcov(fit)), unname(V), tolerance = 1e-8)
	# the weights are the expectation over the predictor's distribution, not the
	# weights at the fitted point
	plug = plogis(fit$EZ[[1]][cells])
	expect_true(all(is.finite(w) & w > 0))
	expect_gt(max(abs(w - plug * (1 - plug))), 1e-4)
})

test_that("the Laplace covariance of the mode is the inverse penalised information", {
	# a logit link makes the expected and observed information equal, so the
	# covariance is the inverse numerical Hessian of the penalised negative
	# log-likelihood at the fitted point
	d = eb_panel(n = 8, seed = 7)
	fit = mode_fit(d, link = "logit")
	n = d$n
	y = d$Y[[1]]; x = d$X[[1]][, , 1]
	obs = !is.na(y)
	P = solve(fit$eb$sab)
	neg_log_post = function(par) {
		eta = par[1] + par[2] * x + outer(par[2 + seq_len(n)], par[2 + n + seq_len(n)], "+")
		al = cbind(par[2 + seq_len(n)], par[2 + n + seq_len(n)])
		0.5 * sum((al %*% P) * al) - sum(dbinom(y[obs], 1, plogis(eta[obs]), log = TRUE))
	}
	par = c(fit$mu, fit$beta, fit$a, fit$b)
	H = stats::optimHess(par, neg_log_post, control = list(ndeps = rep(1e-5, length(par))))
	V = solve(H)[1:2, 1:2]
	expect_equal(unname(vcov(fit)), unname(V), tolerance = 2e-3)
	expect_equal(dimnames(vcov(fit)), list(c("intercept", "x_dyad"), c("intercept", "x_dyad")))
	# conditional on the factors there are none, so the two covariances agree
	expect_equal(unname(fit$eb$vcov), unname(fit$eb$vcov_cond), tolerance = 1e-10)
})

test_that("print, summary, vcov, confint and tidy report the Laplace standard errors", {
	d = eb_panel(seed = 8)
	fit = eb_fit(d)
	msg = gsub("\\s+", " ", paste(capture.output(print(fit)), collapse = " "))
	expect_match(msg, "empirical Bayes", ignore.case = TRUE)
	expect_match(msg, "population intercept")
	expect_false(grepl("Descriptive summaries", msg))
	smsg = gsub("\\s+", " ", paste(capture.output(print(summary(fit))), collapse = " "))
	expect_false(grepl("Descriptive summaries", smsg))
	expect_match(smsg, "Laplace")

	sm = summary(fit)
	expect_true(all(c("Estimate", "Std.Error") %in% colnames(sm$coefficients)))
	expect_true(all(is.finite(sm$coefficients[, "Std.Error"])))
	V = vcov(fit)
	expect_equal(V, t(V))
	expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
	expect_equal(unname(sqrt(diag(V))), unname(fit$eb$se))
	ci = confint(fit)
	expect_equal(rownames(ci), names(coef(fit)))
	expect_true(all(ci[, 1] < coef(fit) & coef(fit) < ci[, 2]))
	expect_equal(unname(ci[, 2] - coef(fit)), unname(qnorm(0.975) * sqrt(diag(V))))
	td = tidy(fit)
	expect_identical(unique(td$se_source), "laplace")
	expect_identical(unique(glance(fit)$se_source), "laplace")
})

test_that("the effects are not re-centred and the intercept is the population intercept", {
	d = eb_panel(seed = 9)
	fit = eb_fit(d)
	x = d$X[[1]][, , 1]
	eta = unname(fit$coefficients["intercept"]) + unname(fit$coefficients["x_dyad"]) * x +
		outer(fit$a, fit$b, "+")
	diag(eta) = NA
	expect_equal(unname(fit$EZ[[1]]), unname(eta), tolerance = 1e-10)
	# at the fixed point of the variational equations the intercept is stationary,
	# so the posterior-mean effects of a prior with mean zero sum to zero (the
	# sum of the effect equations is the intercept equation); the posterior mode
	# with the Laplace correction is not forced to
	expect_lt(abs(sum(fit$a)) + abs(sum(fit$b)), 1e-8)
	md = mode_fit(d)
	expect_gt(abs(sum(md$a)) + abs(sum(md$b)), 1e-6)
	# the fixed-effects fit centres them
	none = eb_fit(d, bias_correct = "none")
	expect_lt(abs(mean(none$a)), 1e-8)
})

test_that("a refit and the parametric bootstrap use the empirical-Bayes estimator", {
	skip_on_cran()
	d = eb_panel(n = 14, seed = 10)
	fit = eb_fit(d)
	re = ame_als_refit(fit, verbose = FALSE)
	expect_identical(re$bias_correct$method, "eb")
	expect_equal(coef(re), coef(fit), tolerance = 1e-3)
	expect_equal(re$eb$sab, fit$eb$sab, tolerance = 1e-3)
	bt = suppressWarnings(suppressMessages(
		ame_als_bootstrap(fit, R = 10, seed = 1, verbose = FALSE)
	))
	expect_s3_class(bt, "boot_ame")
	expect_gte(bt$n_valid, 10L)
	expect_true(all(is.finite(bt$se)))
})

test_that("multistart keeps the empirical-Bayes estimator", {
	skip_on_cran()
	d = eb_panel(n = 14, seed = 11, R = 1)
	fit = eb_fit(d, R = 1, multistart = "cheap")
	expect_identical(fit$bias_correct$method, "eb")
	expect_length(fit$multistart_sse, 4L)
	expect_true(any(fit$multistart_converged))
})

test_that("a symmetric binary fit is empirical Bayes with one effect per actor", {
	d = eb_panel(n = 14, seed = 12, sym = TRUE)
	fit = eb_fit(d, symmetric = TRUE)
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(isTRUE(fit$converged))
	expect_equal(fit$a, fit$b)
	# a scalar variance: the covariance is diagonal with equal entries
	expect_equal(fit$eb$sab[1, 1], fit$eb$sab[2, 2])
	expect_equal(fit$eb$sab[1, 2], 0)
	expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

test_that("an undirected fit with a latent rank is empirical Bayes by default and on request", {
	skip_on_cran()
	d = eb_panel(n = 14, seed = 13, sym = TRUE, R = 1)
	fit = eb_fit(d, R = 1, symmetric = TRUE)
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(fit$bias_correct$default)
	expect_true(isTRUE(fit$converged))
	asked = eb_fit(d, R = 1, symmetric = TRUE, bias_correct = "eb")
	expect_identical(asked$bias_correct$requested, "eb")
	expect_identical(coef(asked), coef(fit))
	expect_equal(dim(fit$U), c(14L, 1L))
	# one scale for the whole factor, stored once per column
	expect_length(fit$eb$sig2_u, 1L)
	expect_true(is.finite(fit$eb$sig2_u) && fit$eb$sig2_u > 0)
	expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
	# a refit keeps the estimator
	re = ame_als_refit(fit, verbose = FALSE)
	expect_identical(re$bias_correct$method, "eb")
	expect_equal(coef(re), coef(fit), tolerance = 1e-2)
	# the fixed-effects rule is available on request
	none = eb_fit(d, R = 1, symmetric = TRUE, bias_correct = "none")
	expect_identical(none$bias_correct$method, "none")
	expect_null(none$eb)
	# the starts of a multistart fit are empirical-Bayes fits too
	ms = eb_fit(d, R = 1, symmetric = TRUE, multistart = "cheap")
	expect_identical(ms$bias_correct$method, "eb")
	expect_length(ms$multistart_sse, 4L)
	expect_true(any(ms$multistart_converged))
})

test_that("the factor scale of an undirected fit carries the posterior second moment", {
	# no weight and no factor: the update is the prior expectation of the squared factor plus the
	# sampler's hyperprior (scale matrix suv0, kappa0 degrees of freedom), so suv0 / kappa0 is its fixed point;
	# without the second moment it falls to the hyperprior alone (the absorbing zero of the factor)
	nr = 20; R = 2
	zero = matrix(0, nr, nr)
	obs = matrix(FALSE, nr, nr)
	kappa0 = 2 * R + 2; suv0 = 1
	with = lame:::.ae_lowrank_em_mm(
		zero,
		zero,
		zero,
		obs,
		R,
		suv0,
		kappa0,
		cold_start = FALSE,
		second_moment = TRUE,
		sig2_in = 0.4
	)
	expect_equal(with$sig2_u, rep((suv0 + nr * R * 0.4) / (nr * R + kappa0), R))
	fixed = lame:::.ae_lowrank_em_mm(
		zero,
		zero,
		zero,
		obs,
		R,
		suv0,
		kappa0,
		cold_start = FALSE,
		second_moment = TRUE,
		sig2_in = suv0 / kappa0
	)
	expect_equal(fixed$sig2_u, rep(suv0 / kappa0, R))
	without = lame:::.ae_lowrank_em_mm(
		zero,
		zero,
		zero,
		obs,
		R,
		suv0,
		kappa0,
		cold_start = FALSE
	)
	expect_equal(without$sig2_u, rep(kappa0 * suv0 / (kappa0 + 2 * nr * R), R))
	expect_lt(without$sig2_u[1], 0.2 * with$sig2_u[1])
	# the trace of one row: sum over rows of tr((sum_j c_ij f_j f_j' + I / s)^-1)
	set.seed(1)
	fm = matrix(rnorm(nr * R), nr, R)
	cm = matrix(runif(nr * nr, 0, 2), nr, nr); cm = (cm + t(cm)) / 2; diag(cm) = 0
	direct = sum(vapply(seq_len(nr), function(i) {
		g = crossprod(fm, cm[i, ] * fm) + diag(1 / 0.7, R)
		sum(diag(solve(g)))
	}, 0))
	expect_equal(lame:::.ae_sym_trace(fm, cm, 0.7), direct)
})

test_that("the conditional trace of the factor rows is the trace of the diagonal blocks of the factor information", {
	skip_on_cran()
	# the trace the scale update adds is the one of the information the Laplace
	# covariance is built from, on the same scale, for both factor kinds
	for (sym in c(FALSE, TRUE)) {
		d = eb_panel(n = 10, seed = 72, mu = -0.6, R = 2, sym = sym, sd_uv = 0.8)
		# the fit's state: the variational means (a weak product is withheld from the report)
		fit = lame:::.ae_product_restore(eb_fit(d, R = 2, symmetric = sym))
		lp = lame:::.ae_eb_laplace_state(fit, inverse = FALSE)
		Hf = lp$jp$Hf; lay = lp$jp$lay; wsum = lp$jp$wsum
		rows = function(i, side) vapply(1:2, function(r) if (side == "u") lay$iu(r)[i] else lay$iv(r)[i], 0L)
		if (sym) {
			fm = sweep(fit$U, 2, sqrt(abs(fit$L)), "*")
			tr_h = sum(vapply(1:10, function(i) { ix = rows(i, "u"); sum(diag(solve(Hf[ix, ix]))) }, 0))
			expect_equal(lame:::.ae_sym_trace(fm, wsum, fit$eb$sig2_u[1]), tr_h, tolerance = 1e-10)
		} else {
			tr_u = sum(vapply(1:10, function(i) { ix = rows(i, "u"); sum(diag(solve(Hf[ix, ix]))) }, 0))
			tr_v = sum(vapply(1:10, function(i) { ix = rows(i, "v"); sum(diag(solve(Hf[ix, ix]))) }, 0))
			one = function(w, f, s) sum(vapply(seq_len(nrow(w)), function(i) {
				sum(diag(solve(crossprod(f, w[i, ] * f) + diag(1 / s, 2))))
			}, 0))
			expect_equal(one(wsum, fit$V, fit$eb$sig2_u), tr_u, tolerance = 1e-10)
			expect_equal(one(t(wsum), fit$U, fit$eb$sig2_v), tr_v, tolerance = 1e-10)
		}
	}
})

test_that("an undirected rank-one term is not shrunk out of the fit", {
	skip_on_cran()
	# in small networks the product of the fit follows the truth
	found = vapply(1:4, function(s) {
		d = eb_panel(n = 30, seed = 200 + s, mu = -0.5, R = 1, sym = TRUE, sd_uv = 0.9)
		# the variational means: the engine does not shrink the term out (whether it is reported is the information's)
		fit = lame:::.ae_product_restore(eb_fit(d, R = 1, symmetric = TRUE))
		truth = d$O
		off = row(truth) != col(truth)
		c(norm = sqrt(sum(fit$UVPM[off]^2)), cor = cor(fit$UVPM[off], truth[off]))
	}, numeric(2))
	expect_true(all(found["norm", ] > 0.5))
	expect_gt(mean(found["cor", ]), 0.5)
})

# predictive probabilities

# the variance of the predictor of a cell is g' S g under the joint Laplace
# covariance S of everything the predictor depends on: with no latent rank the
# predictor is linear in the parameters and the formula is exact
test_that("the predictive variance is the variance of the predictor under the joint Laplace covariance", {
	d = eb_panel(n = 12, seed = 32, mu = -0.7)
	fit = eb_fit(d)
	pr = lame:::.ae_eb_predictive(fit)
	lp = lame:::.ae_eb_laplace_state(fit)
	S = lp$jp$S; ctx = lp$jp$ctx; lay = lp$jp$lay
	x = d$X[[1]][, , 1]
	for (cell in list(c(2, 5), c(7, 3), c(11, 1))) {
		i = cell[1]; j = cell[2]
		g = numeric(nrow(S))
		g[ctx$i_th] = c(1, x[i, j])
		g[lay$ia[i]] = 1; g[lay$ib[j]] = 1
		expect_equal(pr$var[i, j, 1], drop(crossprod(g, S %*% g)), tolerance = 1e-10)
	}
	# with a latent rank the arrays are positive on the observed cells, and an
	# undirected network has the same predictor in both orders
	for (sym in c(FALSE, TRUE)) {
		dr = eb_panel(n = 12, seed = 31, mu = -0.7, R = 1, sym = sym)
		fr = eb_fit(dr, R = 1, symmetric = sym)
		pr = lame:::.ae_eb_predictive(fr)
		expect_equal(dim(pr$var), c(12L, 12L, 1L))
		off = row(pr$var[, , 1]) != col(pr$var[, , 1])
		expect_true(all(pr$var[, , 1][off] > 0))
		expect_true(all(is.finite(pr$mean[, , 1][off])))
		if (sym) {
			expect_equal(pr$var[, , 1], t(pr$var[, , 1]), tolerance = 1e-10)
			expect_equal(pr$mean[, , 1], t(pr$mean[, , 1]), tolerance = 1e-10)
		}
	}
})

# the posterior mean of the effects and factors is the mode moved by minus the
# covariance times q, the gradient of half the log determinant of their
# information; q is checked against central differences of that log determinant
test_that("the displacement of the integrated parameters is the gradient of half the log determinant of their information", {
	skip_on_cran()
	fisher_w = function(eta) {
		tm = lame:::.ae_probit_terms(eta, NULL)
		exp(2 * tm$ld - tm$lp1 - tm$lp0)
	}
	as_array = function(obj) {
		eta = simplify2array(obj$EZ)
		if (length(dim(eta)) == 2L) array(eta, c(dim(eta), 1L)) else eta
	}
	prepare = function(obj) {
		w = fisher_w(as_array(obj))
		obj$obs_weights = w
		obj
	}
	half_ld = function(obj) {
		lp = lame:::.ae_eb_laplace_state(obj, inverse = FALSE)
		ctx = lp$jp$ctx
		iz = setdiff(seq_len(nrow(lp$jp$Hf)), c(ctx$i_th, ctx$i_g))
		sum(log(diag(chol(lp$jp$Hf[iz, iz, drop = FALSE]))))
	}
	# the identity is one of the fitted state: the directed product of this ten-actor panel is below the withholding threshold
	# under the data-scaled factor prior (R/uv_prior.R), so it is read without the rule
	old = options(lame.eb_withhold_product = FALSE)
	on.exit(options(old), add = TRUE)
	for (sym in c(FALSE, TRUE)) {
		d = eb_panel(n = 10, seed = 51, mu = -0.6, R = 1, sym = sym, sd_uv = 1.3)
		# the gradient is that of the information of independent cells: a fit with a dyadic correlation adds the cross weights of its pairs
		fit = prepare(eb_fit(d, R = 1, symmetric = sym, dcor = FALSE))
		expect_true(!sym || abs(fit$L[1]) > 0.5)
		lp = lame:::.ae_eb_laplace_state(fit)
		jp = lp$jp; ctx = jp$ctx; lay = jp$lay
		iz = setdiff(seq_len(nrow(jp$Hf)), c(ctx$i_th, ctx$i_g))
		posz = integer(nrow(jp$Hf)); posz[iz] = seq_along(iz)
		Sz = solve(jp$Hf[iz, iz])
		q = lame:::.ae_eb_logdet_grad(lp, lame:::.ae_eb_slots(lp), Sz, posz)
		Lr = if (sym) abs(fit$L[1]) else 1
		sg = if (sym && fit$L[1] < 0) -1 else 1
		Fm = if (sym) fit$U[, 1] * sqrt(Lr) else NULL
		bump = function(obj, what, i, eps) {
			for (t in seq_along(obj$EZ)) {
				ez = obj$EZ[[t]]
				if (what == "a") {
					ez[i, ] = ez[i, ] + eps
					if (sym) ez[, i] = ez[, i] + eps
				} else if (what == "b") {
					ez[, i] = ez[, i] + eps
				} else if (what == "u") {
					if (sym) {
						dz = sg * eps * Fm
						ez[i, ] = ez[i, ] + dz; ez[, i] = ez[, i] + dz
					} else ez[i, ] = ez[i, ] + eps * obj$V[, 1]
				} else {
					ez[, i] = ez[, i] + eps * obj$U[, 1]
				}
				obj$EZ[[t]] = ez
			}
			if (what == "a") { obj$eb$a_raw[i] = obj$eb$a_raw[i] + eps; if (sym) obj$eb$b_raw[i] = obj$eb$a_raw[i] }
			if (what == "b") obj$eb$b_raw[i] = obj$eb$b_raw[i] + eps
			if (what == "u") obj$U[i, 1] = obj$U[i, 1] + eps / (if (sym) sqrt(Lr) else 1)
			if (what == "v") obj$V[i, 1] = obj$V[i, 1] + eps
			prepare(obj)
		}
		pos_of = function(what, i) switch(what, a = lay$ia[i], b = lay$ib[i], u = lay$iu(1)[i], v = lay$iv(1)[i])
		eps = 1e-5
		for (what in if (sym) c("a", "u") else c("a", "b", "u", "v")) {
			for (i in c(2, 5)) {
				fd = (half_ld(bump(fit, what, i, eps)) - half_ld(bump(fit, what, i, -eps))) / (2 * eps)
				expect_equal(
					unname(q[posz[pos_of(what, i)]]),
					fd,
					tolerance = 1e-5,
					info = paste(if (sym) "symmetric" else "directed", what, i)
				)
			}
		}
	}
})

test_that("the displacement of the effects agrees with the importance-sampled posterior mean", {
	skip_on_cran()
	# a logit link makes the observed and the Fisher information equal: the
	# first-order Laplace displacement of the mean can be compared with the
	# mean of the posterior of the effects, drawn by importance sampling from
	# the Laplace Gaussian
	d = eb_panel(n = 8, seed = 61, mu = -1.6, sd_ab = 0.8)
	fit = mode_fit(d, link = "logit")
	n = 8; y = d$Y[[1]]; x = d$X[[1]][, , 1]; obs = !is.na(y)
	lp = lame:::.ae_eb_laplace_state(fit); jp = lp$jp; ctx = jp$ctx; lay = jp$lay
	iz = setdiff(seq_len(nrow(jp$Hf)), c(ctx$i_th, ctx$i_g))
	posz = integer(nrow(jp$Hf)); posz[iz] = seq_along(iz)
	it = c(ctx$i_th, ctx$i_g)
	Sf = jp$S
	Sz = Sf[iz, iz] - Sf[iz, it] %*% solve(Sf[it, it], Sf[it, iz])
	q = lame:::.ae_eb_logdet_grad(lp, lame:::.ae_eb_slots(lp), Sz, posz)
	shift = -as.vector(Sz %*% q)
	mu = fit$eb$mu_int; beta = unname(coef(fit)["x_dyad"])
	P = solve(fit$eb$sab)
	logpost = function(zz) {
		a = zz[1:n]; b = zz[n + 1:n]
		eta = mu + beta * x + outer(a, b, "+")
		al = cbind(a, b)
		sum(dbinom(y[obs], 1, plogis(eta[obs]), log = TRUE)) - 0.5 * sum((al %*% P) * al)
	}
	zhat = c(fit$eb$a_raw, fit$eb$b_raw)
	ord = c(posz[lay$ia], posz[lay$ib])
	ch = chol(Sz[ord, ord])
	set.seed(1)
	ns = 20000
	draws = matrix(rnorm(ns * 2 * n), ns) %*% ch + rep(zhat, each = ns)
	lw = apply(draws, 1, logpost) - (-0.5 * rowSums(((draws - rep(zhat, each = ns)) %*% solve(ch))^2))
	w = exp(lw - max(lw)); w = w / sum(w)
	is_shift = colSums(draws * w) - zhat
	expect_gt(cor(shift[ord], is_shift), 0.95)
	expect_gt(sum(shift[ord] * is_shift) / sum(is_shift^2), 0.8)
	expect_lt(sum(shift[ord] * is_shift) / sum(is_shift^2), 1.25)
})

test_that("an empirical-Bayes fit predicts with the posterior of its predictor", {
	d = eb_panel(n = 16, seed = 33, mu = -1.2)
	fit = eb_fit(d)
	plug = pnorm(fit$EZ[[1]])
	expect_equal(unname(fitted(fit, plug_in = TRUE)), unname(plug), tolerance = 1e-12)
	pred = fitted(fit)
	expect_equal(dim(pred), dim(plug))
	expect_true(all(is.na(diag(pred))))
	off = row(pred) != col(pred)
	expect_true(all(pred[off] > 0 & pred[off] < 1))
	# the predictive probabilities are not the plug-in ones, and they keep the
	# rare events away from zero: the plug-in value of a rare cell is lower
	expect_gt(max(abs(pred[off] - plug[off])), 1e-4)
	lowest = off & plug < quantile(plug[off], 0.05)
	expect_gt(mean(pred[lowest]), mean(plug[lowest]))
	expect_equal(unname(predict(fit)), unname(pred), tolerance = 1e-12)
	expect_equal(unname(predict(fit, plug_in = TRUE)), unname(plug), tolerance = 1e-12)
	# the residuals follow the fitted values
	expect_equal(unname(residuals(fit)), unname(d$Y[[1]] - pred), tolerance = 1e-12)
	expect_equal(unname(residuals(fit, plug_in = TRUE)), unname(d$Y[[1]] - plug), tolerance = 1e-12)
	expect_equal(
		unname(residuals(eb_fit(d, bias_correct = "none"))),
		unname(d$Y[[1]] - fitted(eb_fit(d, bias_correct = "none"))),
		tolerance = 1e-12
	)
	# the link scale is the posterior mean of the predictor: for a variational fit
	# the fitted predictor itself, for the posterior mode with the Laplace
	# correction it sits below the mode on average (the effects of rare cells
	# are skewed)
	lk = predict(fit, type = "link")
	expect_equal(unname(lk[off]), unname(fit$EZ[[1]][off]), tolerance = 1e-10)
	md = mode_fit(d)
	expect_lt(mean(predict(md, type = "link")[off]), mean(md$EZ[[1]][off]))
	expect_equal(unname(predict(fit, type = "link", plug_in = TRUE)), unname(fit$EZ[[1]]), tolerance = 1e-12)
	# the mean is linear in the covariates: shifting the design shifts it by
	# the coefficient
	xa = d$X[[1]]
	lk2 = predict(fit, newdata = xa + 0.5, type = "link")
	expect_equal(unname(lk2[off] - lk[off]), rep(0.5 * unname(coef(fit)["x_dyad"]), sum(off)), tolerance = 1e-10)
	expect_equal(unname(predict(fit, newdata = xa)), unname(pred), tolerance = 1e-12)
})

test_that("longitudinal, bipartite and latent-rank fits predict, and a fit without the posterior falls back with a warning", {
	d = eb_panel(n = 10, Tt = 3, seed = 34, R = 1)
	fit = eb_fit(d, R = 1)
	pr = fitted(fit)
	expect_type(pr, "list")
	expect_length(pr, 3L)
	expect_equal(names(pr), names(fit$EZ))
	expect_true(all(vapply(pr, function(m) all(m[row(m) != col(m)] > 0 & m[row(m) != col(m)] < 1), TRUE)))
	expect_equal(unname(fitted(fit, plug_in = TRUE)[[2]]), unname(pnorm(fit$EZ[[2]])), tolerance = 1e-12)

	set.seed(35)
	nr = 12; nc = 8
	rn = sprintf("r%02d", seq_len(nr)); cn = sprintf("c%02d", seq_len(nc))
	a = rnorm(nr, 0, .6); b = rnorm(nc, 0, .6)
	x = matrix(rnorm(nr * nc), nr, nc)
	y = 1 * (-0.5 + 0.7 * x + outer(a, b, "+") + matrix(rnorm(nr * nc), nr, nc) > 0)
	dimnames(y) = list(rn, cn)
	xa = array(x, c(nr, nc, 1), dimnames = list(rn, cn, "x_dyad"))
	for (rank in 0:1) {
		bf = suppressWarnings(suppressMessages(
			ame_als(y, Xdyad = xa, R = rank, family = "binary", mode = "bipartite", verbose = FALSE)
		))
		pb = fitted(bf)
		expect_equal(dim(pb), c(nr, nc))
		expect_false(anyNA(pb))
		expect_true(all(pb > 0 & pb < 1))
	}

	# a fit that does not carry the state the posterior is built from
	old = fit
	old$eb$mu_int = NULL
	expect_warning(p_old <- fitted(old, plug_in = FALSE), NA)
	expect_equal(unname(p_old[[1]]), unname(pnorm(fit$EZ[[1]])), tolerance = 1e-12)
})

test_that("the held-out cells of an empirical-Bayes fit are scored by the predictive probabilities", {
	skip_on_cran()
	# a sparse panel: the plug-in value of an event the fit has not seen is
	# often far below the rate of such cells; the predictive probability is
	# not, and it scores better
	scores = vapply(1:6, function(s) {
		d = eb_panel(n = 40, seed = 400 + s, mu = -2.4)
		y = d$Y[[1]]
		set.seed(s)
		held = sample(which(!is.na(y)), round(0.2 * sum(!is.na(y))))
		yt = y; yt[held] = NA
		dt = d; dt$Y[[1]] = yt
		fit = eb_fit(dt)
		ll = function(p) -mean(ifelse(y[held] == 1, log(p[held]), log1p(-p[held])))
		c(plug = ll(fitted(fit, plug_in = TRUE)), pred = ll(fitted(fit)))
	}, numeric(2))
	expect_lt(mean(scores["pred", ]), mean(scores["plug", ]))
})

# the evidence

test_that("the Laplace evidence prefers the rank the data carry", {
	skip_on_cran()
	rank_by_evidence = function(d, sym = FALSE) {
		ev = vapply(0:2, function(r) eb_fit(d, R = r, symmetric = sym)$eb$log_marginal, 0)
		which.max(ev) - 1L
	}
	with_rank = vapply(
		1:3,
		function(s) rank_by_evidence(eb_panel(n = 40, seed = 500 + s, mu = -0.8, R = 1, sd_uv = 1)),
		0L
	)
	without = vapply(1:3, function(s) rank_by_evidence(eb_panel(n = 40, seed = 600 + s, mu = -0.8)), 0L)
	expect_true(all(with_rank >= 1L))
	expect_true(all(without == 0L))
	# the stored value is the one rebuilt from the fitted object
	for (sym in c(FALSE, TRUE)) {
		fit = eb_fit(eb_panel(n = 14, seed = 36, R = 1, sym = sym), R = 1, symmetric = sym)
		expect_equal(fit$eb$log_marginal, lame:::.ae_eb_evidence(fit)$log_marginal)
		expect_true(is.finite(fit$eb$log_marginal))
	}
})

test_that("an undirected term with a negative eigenvalue is recovered with its sign", {
	skip_on_cran()
	n = 30; nm = sprintf("a%02d", seq_len(n))
	res = vapply(1:3, function(s) {
		set.seed(300 + s)
		a = rnorm(n, 0, .6); u = rnorm(n, 0, 1)
		x = matrix(rnorm(n * n), n); x[lower.tri(x)] = t(x)[lower.tri(x)]
		e = matrix(rnorm(n * n), n); e[lower.tri(e)] = t(e)[lower.tri(e)]
		O = -tcrossprod(u)
		y = 1 * (-1 + 0.7 * x + outer(a, a, "+") + O + e > 0); diag(y) = NA
		dimnames(y) = list(nm, nm)
		xa = array(x, c(n, n, 1), dimnames = list(nm, nm, "x_dyad"))
		fit = suppressWarnings(suppressMessages(
			ame_als(y, Xdyad = xa, R = 1, family = "binary", symmetric = TRUE, verbose = FALSE)
		))
		off = row(O) != col(O)
		c(L = fit$L[1], cor = cor(fit$UVPM[off], O[off]))
	}, numeric(2))
	expect_true(all(res["L", ] < 0))
	expect_gt(mean(res["cor", ]), 0.5)
})

# slower checks

test_that("the empirical-Bayes fit removes the fixed-effects bias of the slope", {
	skip_on_cran()
	slopes = vapply(1:16, function(s) {
		d = eb_panel(n = 30, seed = 100 + s, mu = -1)
		c(
			none = unname(coef(eb_fit(d, bias_correct = "none"))["x_dyad"]),
			eb = unname(coef(eb_fit(d))["x_dyad"])
		)
	}, numeric(2))
	expect_gt(mean(slopes["none", ]) - 0.7, 0.03)
	expect_lt(abs(mean(slopes["eb", ]) - 0.7), 0.03)
})

test_that("the Laplace standard errors agree with the parametric bootstrap", {
	skip_on_cran()
	d = eb_panel(n = 30, seed = 21, mu = -1, R = 1)
	fit = eb_fit(d, R = 1)
	bt = suppressWarnings(suppressMessages(
		ame_als_bootstrap(fit, R = 40, seed = 2, verbose = FALSE)
	))
	ratio = bt$se[c("intercept", "x_dyad")] / fit$eb$se[c("intercept", "x_dyad")]
	expect_true(all(ratio > 0.6 & ratio < 1.6))
})

test_that("node covariates are fixed effects of the prior mean", {
	skip_on_cran()
	d = eb_panel(n = 30, seed = 22, mu = -1)
	set.seed(23)
	nm = rownames(d$Y[[1]])
	xr = matrix(rnorm(30), 30, 1, dimnames = list(nm, "xr"))
	xc = matrix(rnorm(30), 30, 1, dimnames = list(nm, "xc"))
	fit = suppressWarnings(suppressMessages(
		ame_als(d$Y[[1]], Xdyad = d$X[[1]], Xrow = xr, Xcol = xc, family = "binary", verbose = FALSE)
	))
	expect_identical(fit$bias_correct$method, "eb")
	expect_identical(names(coef(fit)), c("intercept", "x_dyad", "xr_row", "xc_col"))
	expect_equal(dim(vcov(fit)), c(4L, 4L))
	eta = fit$coefficients["intercept"] + fit$coefficients["x_dyad"] * d$X[[1]][, , 1] +
		fit$coefficients["xr_row"] * outer(xr[, 1], rep(1, 30)) +
		fit$coefficients["xc_col"] * outer(rep(1, 30), xc[, 1]) + outer(fit$a, fit$b, "+")
	diag(eta) = NA
	expect_equal(unname(fit$EZ[[1]]), unname(eta), tolerance = 1e-10)
})

test_that("the explicit bootstrap and jackknife corrections are unchanged", {
	skip_on_cran()
	d = eb_panel(n = 60, seed = 31, mu = 0)
	bt = eb_fit(d, bias_correct = "bootstrap")
	expect_identical(bt$bias_correct$method, "bootstrap")
	expect_false(is.null(bt$raw))
	expect_null(bt$eb)
	jk = eb_fit(d, bias_correct = "jackknife")
	expect_identical(jk$bias_correct$method, "jackknife")
	expect_false(is.null(jk$raw))
})

test_that("the longitudinal empirical-Bayes fit with a latent factor converges", {
	skip_on_cran()
	d = eb_panel(n = 30, Tt = 4, seed = 41, mu = -1, R = 1, sd_uv = 0.8)
	fit = eb_fit(d, R = 1)
	expect_true(isTRUE(fit$converged))
	expect_identical(fit$irls_diagnostics$termination, "converged")
	expect_lt(fit$irls_diagnostics$max_score, 1e-3)
	# the Laplace covariance integrates the factors out, so it is the wider
	expect_true(all(sqrt(diag(fit$eb$vcov)) >= sqrt(diag(fit$eb$vcov_cond)) - 1e-10))
})

test_that("sparse panels with a latent factor reach the full certificate", {
	skip_on_cran()
	# sparse directed probit panels from the estimator's simulation design (sender and
	# receiver sd 0.6, correlation 0.3, one factor of variance 0.6)
	sim = function(n, mu, seed) {
		set.seed(seed)
		ab = matrix(rnorm(2 * n), n, 2) %*% chol(0.36 * matrix(c(1, 0.3, 0.3, 1), 2))
		U = matrix(rnorm(n, 0, sqrt(0.6)), n, 1); V = matrix(rnorm(n, 0, sqrt(0.6)), n, 1)
		nm = sprintf("a%03d", seq_len(n))
		x = matrix(rnorm(n * n), n, n)
		y = 1 * (mu + 0.7 * x + outer(ab[, 1], ab[, 2], "+") + tcrossprod(U, V) +
				matrix(rnorm(n * n), n, n) > 0)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		list(Y = y, X = array(x, c(n, n, 1), dimnames = list(nm, nm, "x_dyad")))
	}
	for (case in list(
		c(30, -1, 7003058),
		c(50, -2.2, 7019029),
		c(50, -2.2, 7019044),
		c(50, -2.2, 7019062),
		c(50, -2.2, 7019064)
	)) {
		d = sim(case[1], case[2], case[3])
		fit = suppressMessages(ame_als(
			d$Y,
			Xdyad = d$X,
			R = 1,
			family = "binary",
			verbose = FALSE
		))
		expect_true(fit$converged)
		expect_identical(fit$irls_diagnostics$termination, "converged")
		expect_lte(fit$irls_diagnostics$max_score, fit$irls_diagnostics$threshold)
	}
})

test_that("a panel whose mode does not converge is refused with the convergence class", {
	skip_on_cran()
	d = eb_panel(n = 30, Tt = 2, seed = 51, mu = -1, R = 1)
	err = tryCatch(
		suppressMessages(lame_als(d$Y, Xdyad = d$X, R = 1, family = "binary", max_iter = 3, verbose = FALSE)),
		error = function(e) e
	)
	expect_s3_class(err, "lame_als_convergence")
	msg = conditionMessage(err)
	# the state and the remedy
	expect_match(msg, "max_iter")
	expect_match(msg, "smaller rank")
	expect_match(msg, "rank 1")
})

test_that("the warm-start seed of an empirical-Bayes fit carries its variance components", {
	d = eb_panel(n = 14, Tt = 3, seed = 21)
	fit = eb_fit(d)
	seed = lame:::.ae_static_seed(fit)
	expect_true(all(c("a_level", "b_level", "eb") %in% names(seed)))
	expect_identical(seed$eb$sab, fit$eb$sab)
})

test_that("a dynamic empirical-Bayes fit is refit from its own seed without iterating", {
	skip_on_cran()
	d = eb_panel(n = 14, Tt = 4, seed = 22)
	fit = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		method = "als",
		dynamic_ab = TRUE,
		verbose = FALSE
	)))
	# the dynamic empirical-Bayes estimator is a separate piece of the package; skip where it is absent
	skip_if(!inherits(fit, "lame_dynamic_als") || !identical(fit$bias_correct$method, "eb"))
	seed = lame:::.ae_static_seed(fit)
	expect_true(is.numeric(seed$a_level) && is.numeric(seed$b_level) && is.list(seed$eb))
	for (bc in list(NULL, "eb")) {
		g = suppressWarnings(suppressMessages(ame_als_refit(fit, bias_correct = bc)))
		expect_identical(g$bias_correct$method, "eb")
		expect_equal(g$iterations, 0)
		expect_equal(unname(coef(g)), unname(coef(fit)))
	}
})
