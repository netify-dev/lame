skip_on_cran()

library(lame)

# dynamic empirical-Bayes estimator for normal panels (R/lame_dynamic_als_eb_nrm.R): a gaussian cell in the binary and
# count engine, the residual variance as a hyperparameter, the dyadic correlation as an exact bivariate-normal pair,
# and the expected cell under position noise

# small normal panel from the dynamic_ab model: y = mu + x'beta + a_it + b_jt + e, e ~ N(0, ve), each actor a level
# plus a zero-mean AR(1) deviation, and dyadic residual correlation rho_d in directed panels
nrm_sim = function(
	n = 14,
	nb = n,
	Tt = 4,
	mode = "directed",
	seed = 1,
	mu = 0.3,
	beta = c(0.4, -0.3),
	va = 0.5,
	vb = 0.5,
	rho = 0.7,
	sigma = 0.3,
	ve = 1,
	rho_d = 0,
	R = 0,
	suv = 1,
	dyn_uv = FALSE,
	drop = 0,
	scale = 1
) {
	set.seed(seed)
	sym = identical(mode, "symmetric")
	bip = identical(mode, "bipartite")
	if (!bip) nb = n
	ar1 = function(m, Tt) {
		out = matrix(0, m, Tt)
		out[, 1] = rnorm(m, 0, sigma / sqrt(1 - rho^2))
		for (t in 2:Tt) out[, t] = rho * out[, t - 1] + rnorm(m, 0, sigma)
		out
	}
	alpha = rnorm(n, 0, sqrt(va))
	gam = if (sym) alpha else rnorm(nb, 0, sqrt(vb))
	d = ar1(n, Tt)
	e = if (sym) d else ar1(nb, Tt)
	if (R > 0) {
		if (dyn_uv) {
			U = array(0, c(n, R, Tt)); V = array(0, c(nb, R, Tt))
			for (r in seq_len(R)) { U[, r, ] = ar1(n, Tt) * 0.8; V[, r, ] = if (sym) U[, r, ] else ar1(nb, Tt) * 0.8 }
		} else {
			U0 = matrix(rnorm(n * R, 0, sqrt(suv)), n, R)
			V0 = if (sym) U0 else matrix(rnorm(nb * R, 0, sqrt(suv)), nb, R)
		}
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
		if (R > 0) eta = eta + if (dyn_uv) tcrossprod(
			matrix(U[, , t], n, R),
			matrix(V[, , t], nb, R)
		) else tcrossprod(U0, V0)
		err = matrix(rnorm(n * nb, 0, sqrt(ve)), n, nb)
		if (rho_d != 0 && !sym && !bip) {
			lo = lower.tri(err)
			err[lo] = rho_d * t(err)[lo] + sqrt(1 - rho_d^2) * matrix(rnorm(n * nb, 0, sqrt(ve)), n, nb)[lo]
		}
		if (sym) err[lower.tri(err)] = t(err)[lower.tri(err)]
		y = scale * (eta + err)
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

nrm_fit = function(d, ...) {
	suppressWarnings(lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, verbose = FALSE, ...))
}

# panel context for the internal functions
nrm_ctx = function(d, dcor = FALSE, R = 0L, ...) {
	ns = asNamespace("lame")
	prep = ns$.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ns$.lde_context(prep, FALSE, "identity", list(), R = R, dynamic_ab = TRUE, dcor = dcor, ...)
}

# the cell

test_that("the normal cell's working problem is the gaussian problem, for a single cell and for a pair", {
	ns = asNamespace("lame")
	d = nrm_sim(n = 8, Tt = 3, seed = 2, rho_d = 0.5)
	ctx = nrm_ctx(d, dcor = TRUE)
	set.seed(3)
	eta = array(rnorm(length(ctx$y)), dim(ctx$y))
	hyp = list(ve = 0.7, rho_d = 0.4)
	wk = ns$.lde_work_nrm(ctx, eta, hyp)
	obs = ctx$obs
	r = ctx$y - eta
	pc = ctx$pc
	# an unpaired cell is the independent gaussian cell: weight 1 / ve, score r / ve, response y
	single = obs
	single[c(pc$one, pc$two)] = FALSE
	expect_equal(wk$w[single], rep(1 / 0.7, sum(single)))
	expect_equal(wk$score[single], r[single] / 0.7)
	expect_equal(wk$z[single], ctx$y[single])
	# a pair: bivariate normal precision w = 1 / (ve (1 - rho^2)), cross weight -rho w, score w r1 - rho w r2
	w = 1 / (0.7 * (1 - 0.4^2))
	expect_equal(wk$w[pc$one], rep(w, length(pc$one)))
	expect_equal(wk$cw[pc$one], rep(-0.4 * w, length(pc$one)))
	expect_equal(wk$score[pc$one], w * r[pc$one] - 0.4 * w * r[pc$two])
	expect_equal(wk$z[pc$one], eta[pc$one] + wk$score[pc$one] / wk$w[pc$one])
	# the weight's derivatives in the predictor vanish: a gaussian cell has constant information
	expect_true(all(wk$hprime == 0) && all(wk$hpp == 0))
	# the deviance is -2 times the log density of independent cells and bivariate normal pairs
	S = 0.7 * matrix(c(1, 0.4, 0.4, 1), 2)
	ll_pair = vapply(seq_along(pc$one), function(k) {
		rr = c(r[pc$one[k]], r[pc$two[k]])
		-0.5 * (sum(rr * solve(S, rr)) + log(det(2 * pi * S)))
	}, 0)
	ll_single = stats::dnorm(r[single], 0, sqrt(0.7), log = TRUE)
	expect_equal(wk$dev, -2 * (sum(ll_pair) + sum(ll_single)))
	# under noise v the expected deviance adds v / ve (v1 + v2 for a pair); the marginal deviance uses variance ve + v
	v = array(0.2, dim(ctx$y))
	wn = ns$.lde_work_nrm(ctx, eta, hyp, v)
	expect_equal(wn$w, wk$w)
	expect_equal(wn$score, wk$score)
	expect_equal(wn$dev - wk$dev, 0.2 * sum(single) / 0.7 + 0.2 * 2 * length(pc$one) * (1 / (0.7 * (1 - 0.4^2))))
	Sn = S + diag(0.2, 2)
	llm_pair = vapply(seq_along(pc$one), function(k) {
		rr = c(r[pc$one[k]], r[pc$two[k]])
		-0.5 * (sum(rr * solve(Sn, rr)) + log(det(2 * pi * Sn)))
	}, 0)
	expect_equal(wn$dev_marg, -2 * (sum(llm_pair) + sum(stats::dnorm(r[single], 0, sqrt(0.9), log = TRUE))))
})

# the posterior at fixed hyperparameters is exact

# dense posterior mean and covariance of (coefficients, levels, paths) at a fit's hyperparameters
nrm_dense = function(ctx, hy, Re_fun = NULL) {
	ns = asNamespace("lame")
	na = ctx$na; Tt = ctx$Tt; q = ctx$q; p = ctx$mth
	lam = ns$.lde_level_precision(hy$sab, ctx$kind, na, ctx$nb)
	s = 1 / hy$s2
	rho = hy$rho
	Pth = ns$.lde_theta_prior(ctx, hy)
	# stationary AR(1) deviation prior: precision diagonal (1, 1 + rho^2, ..., 1 + rho^2, 1), -rho off it, over s2
	Qd = matrix(0, Tt, Tt)
	diag(Qd) = c(1, rep(1 + rho^2, max(Tt - 2, 0)), 1)[seq_len(Tt)]
	if (Tt >= 2) for (t in 1:(Tt - 1)) Qd[t, t + 1] = Qd[t + 1, t] = -rho
	Qd = Qd * s
	nu = p + q + q * Tt
	iu_th = seq_len(p); iu_lev = p + seq_len(q); iu_x = function(t) p + q + (t - 1) * q + seq_len(q)
	P = matrix(0, nu, nu)
	P[iu_th, iu_th] = Pth
	P[iu_lev, iu_lev] = lam
	Dm = matrix(0, q * Tt, q + q * Tt)
	for (t in 1:Tt) for (e in 1:q) { r = (e - 1) * Tt + t; Dm[r, e] = -1; Dm[r, q + (t - 1) * q + e] = 1 }
	Pd = crossprod(Dm, kronecker(diag(q), Qd) %*% Dm)
	ix = c(iu_lev, p + q + seq_len(q * Tt))
	P[ix, ix] = P[ix, ix] + Pd
	A = NULL; yv = numeric(0); cells = NULL
	for (t in 1:Tt) {
		ob = which(as.vector(ctx$obs[, , t]))
		ii = ((ob - 1) %% na) + 1; jj = ((ob - 1) %/% na) + 1
		Ar = matrix(0, length(ob), nu)
		Ar[, iu_th[ctx$thidx[[t]]]] = ctx$Xc[[t]][ob, , drop = FALSE]
		for (k in seq_along(ob)) { Ar[k, iu_x(t)[ii[k]]] = 1; Ar[k, iu_x(t)[if (ctx$sym) ii[k] else na + jj[k]]] = 1 }
		A = rbind(A, Ar)
		yv = c(yv, as.vector(ctx$y[, , t])[ob])
		cells = rbind(cells, cbind(t = t, i = ii, j = jj))
	}
	Re = diag(1 / hy$ve, nrow(A))
	if (isTRUE(ctx$dcor) && ctx$pc$n_pairs > 0L) {
		key = paste(cells[, "t"], cells[, "i"], cells[, "j"])
		inf = ns$.lde_nrm_pair_info(hy$ve, hy$rho_d)
		for (k in seq_len(nrow(cells))) if (cells[k, "i"] < cells[k, "j"]) {
			m = match(paste(cells[k, "t"], cells[k, "j"], cells[k, "i"]), key)
			if (!is.na(m)) { Re[k, k] = Re[m, m] = inf$w; Re[k, m] = Re[m, k] = inf$c }
		}
	}
	H = crossprod(A, Re %*% A) + P
	list(
		u = solve(H, crossprod(A, Re %*% yv)),
		H = H,
		A = A,
		iu_th = iu_th,
		iu_lev = iu_lev,
		iu_x = iu_x,
		p = p,
		nu = nu,
		cells = cells
	)
}

test_that("the fit is the exact posterior mean at its hyperparameters, with and without the dyadic correlation", {
	skip_on_cran()
	ns = asNamespace("lame")
	for (dc in c(FALSE, TRUE)) {
		d = nrm_sim(n = 8, Tt = 3, seed = 4 + dc, rho_d = 0.5, drop = 0.1)
		ctx = nrm_ctx(d, dcor = dc)
		fit = ns$.lde_fit(ctx, max_iter = 600L, tol = 1e-10)
		expect_true(fit$converged)
		st = fit$st
		dn = nrm_dense(ctx, st$hyp)
		expect_lt(max(abs(dn$u[dn$iu_th] - st$theta)), 1e-6)
		expect_lt(max(abs(dn$u[dn$iu_lev] - st$lev)), 1e-5)
		expect_lt(max(abs(unlist(lapply(1:ctx$Tt, function(t) dn$u[dn$iu_x(t)])) - as.vector(st$x))), 1e-6)
		# the posterior variance of a cell's additive predictor given the coefficients is the engine's vcell
		ix = unlist(lapply(1:ctx$Tt, function(t) dn$iu_x(t)))
		Vc = solve(dn$H[-dn$iu_th, -dn$iu_th])
		tt = dn$cells[, "t"]; ii = dn$cells[, "i"]; jj = dn$cells[, "j"]
		vd = vapply(seq_len(nrow(dn$cells)), function(k) {
			a = rep(0, dn$nu - dn$p)
			a[(dn$iu_x(tt[k]) - dn$p)[ii[k]]] = 1
			a[(dn$iu_x(tt[k]) - dn$p)[ctx$na + jj[k]]] = 1
			drop(a %*% Vc %*% a)
		}, 0)
		vc_eng = vapply(seq_len(nrow(dn$cells)), function(k) fit$ev$lap$vcell[ii[k], jj[k], tt[k]], 0)
		expect_equal(vc_eng, vd, tolerance = 1e-8)
		if (dc) {
			# and the covariance of the two predictors of a dyad
			pc = ctx$pc
			cvd = vapply(seq_len(pc$n_pairs), function(k) {
				a1 = rep(0, dn$nu - dn$p); a2 = a1
				a1[(dn$iu_x(pc$t[k]) - dn$p)[pc$i[k]]] = 1; a1[(dn$iu_x(pc$t[k]) - dn$p)[ctx$na + pc$j[k]]] = 1
				a2[(dn$iu_x(pc$t[k]) - dn$p)[pc$j[k]]] = 1; a2[(dn$iu_x(pc$t[k]) - dn$p)[ctx$na + pc$i[k]]] = 1
				drop(a1 %*% Vc %*% a2)
			}, 0)
			expect_equal(fit$ev$lap$cvpair, cvd, tolerance = 1e-8)
		}
	}
})

test_that("the residual variance and the dyadic correlation are at the fixed points of their EM images", {
	ns = asNamespace("lame")
	d = nrm_sim(n = 10, Tt = 4, seed = 6, rho_d = 0.5)
	ctx = nrm_ctx(d, dcor = TRUE)
	fit = ns$.lde_fit(ctx, max_iter = 800L, tol = 1e-10)
	expect_true(fit$converged)
	hy = fit$st$hyp
	ev = fit$ev
	eta = ns$.lde_eta(ctx, fit$st)
	# ve is the posterior mean of the sampler's inverse-gamma conditional given the expected sum of squares
	q = ns$.lde_nrm_q_sum(ctx, eta, ev$lap$vcell, ev$lap$cvpair, hy$rho_d)
	pr = ctx$prior$ve_prior
	expect_equal(unname(pr[["nu0"]]), 1)
	expect_equal(unname(pr[["s20"]]), stats::var(ctx$y[ctx$obs]))
	expect_equal(hy$ve, (pr[["nu0"]] * pr[["s20"]] + q) / (sum(ctx$obs) + pr[["nu0"]] - 2), tolerance = 1e-7)
	# rho_d is the root of the cubic of the pairs' expected complete-data log posterior
	rs = ns$.lde_nrm_rho_step(ctx, hy, eta, ev$lap)
	expect_equal(rs$rho, hy$rho_d, tolerance = 1e-7)
	pc = ctx$pc
	r = ctx$y - eta
	S11 = sum(r[pc$one]^2 + ev$lap$vcell[pc$one] + r[pc$two]^2 + ev$lap$vcell[pc$two])
	S12 = sum(r[pc$one] * r[pc$two] + ev$lap$cvpair)
	M = pc$n_pairs
	lp = function(rho) -((M - 2) / 2) * log(1 - rho^2) - (S11 - 2 * rho * S12) / (2 * hy$ve * (1 - rho^2))
	expect_gte(lp(hy$rho_d), max(lp(hy$rho_d + c(-0.01, 0.01))))
})

# the route

test_that("a dynamic normal fit is an empirical-Bayes fit of the Gaussian model by default and the penalised fit is still available", {
	d = nrm_sim(n = 14, Tt = 4, seed = 11, rho_d = 0.4)
	f = nrm_fit(d)
	expect_identical(f$family, "normal")
	expect_identical(f$bias_correct$method, "eb")
	expect_true(f$converged)
	expect_true(all(c("ve", "ve_prior") %in% names(f$eb)))
	expect_equal(unname(f$VC[["ve"]]), f$eb$ve)
	expect_equal(f$s2, f$eb$ve)
	expect_gt(f$rho_ab, 0)
	expect_lt(f$rho_ab, 0.99)
	# a directed panel estimates the dyadic correlation; bipartite and undirected panels have none
	expect_gt(unname(f$VC[["rho"]]), 0.1)
	expect_equal(unname(f$VC[["rho"]]), f$eb$pair$rho)
	expect_true(is.finite(f$eb$pair$sd))
	expect_identical(f$rho_estimation$n_pairs, f$eb$pair$n_pairs)
	f0 = nrm_fit(d, dcor = FALSE)
	expect_null(f0$eb$pair)
	# an explicit request for the estimator gives the default
	g = nrm_fit(d, bias_correct = "eb")
	expect_equal(g$coefficients, f$coefficients)
	# the penalised fit keeps its fixed persistence and moment estimates
	p = nrm_fit(d, bias_correct = "none")
	expect_identical(p$bias_correct$method, "none")
	expect_null(p$eb)
	expect_equal(p$rho_ab, 0.8)
	# bootstrapping the empirical-Bayes estimator is refused; the penalised one allows it
	expect_error(nrm_fit(d, bootstrap = 5L), class = "lame_dynamic_eb_unsupported")
	expect_error(nrm_fit(d, lambda_ab = 1), class = "lame_dynamic_eb_unsupported")
	# the estimator covers normal, binary and count panels
	expect_error(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "normal",
			dynamic_ab = TRUE,
			bias_correct = "jackknife",
			verbose = FALSE
		),
		class = "lame_jk_split"
	)
})

test_that("directed, bipartite and undirected normal panels, with and without missing cells and moving coefficients, fit and report", {
	skip_on_cran()
	cases = list(
		nrm_sim(n = 12, Tt = 4, seed = 21, rho_d = 0.4),
		nrm_sim(n = 12, Tt = 4, seed = 22, drop = 0.15),
		nrm_sim(n = 12, nb = 8, Tt = 4, mode = "bipartite", seed = 23),
		nrm_sim(n = 12, Tt = 4, mode = "symmetric", seed = 24)
	)
	modes = c("unipartite", "unipartite", "bipartite", "unipartite")
	for (k in seq_along(cases)) {
		d = cases[[k]]
		f = suppressWarnings(
			lame:::lame_dynamic_als(
				d$Y,
				d$X,
				family = "normal",
				dynamic_ab = TRUE,
				mode = modes[k],
				symmetric = k == 4L,
				verbose = FALSE
			)
		)
		expect_true(f$converged, info = k)
		expect_true(all(is.finite(f$coefficients)) && all(is.finite(f$eb$se)), info = k)
		expect_gt(f$VC[["ve"]], 0.5)
		expect_lt(f$VC[["ve"]], 1.6)
	}
	d = nrm_sim(n = 12, Tt = 5, seed = 25)
	f = nrm_fit(d, dynamic_beta = c("intercept", "x1_dyad"), dynamic_beta_kind = "rw1")
	expect_true(f$converged)
	expect_true(all(f$beta_dynamic_mask[c("intercept", "x1_dyad")]))
	expect_true(all(is.finite(f$coef_path)))
})

# recovery, scale, methods

test_that("a normal fit recovers the model of a simulated panel", {
	skip_on_cran()
	d = nrm_sim(n = 30, Tt = 5, seed = 31, rho_d = 0.5, ve = 1, va = 0.5, vb = 0.5)
	f = nrm_fit(d)
	expect_true(f$converged)
	expect_lt(abs(f$VC[["ve"]] - 1), 0.1)
	expect_lt(abs(f$VC[["rho"]] - 0.5), 0.1)
	expect_lt(abs(f$sigma_ab - 0.3), 0.06)
	expect_lt(abs(f$rho_ab - 0.7), 0.2)
	expect_lt(max(abs(f$coefficients[-1] - c(0.4, -0.3))), 0.05)
	# standard errors reflect the slopes' sampling, so the truth is within three
	expect_lt(max(abs(f$coefficients[-1] - c(0.4, -0.3)) / f$eb$se[-1]), 3)
})

test_that("the fit is equivariant to the unit of the response", {
	skip_on_cran()
	d = nrm_sim(n = 14, Tt = 4, seed = 41, rho_d = 0.3)
	d10 = nrm_sim(n = 14, Tt = 4, seed = 41, rho_d = 0.3, scale = 10)
	f = nrm_fit(d)
	g = nrm_fit(d10)
	expect_equal(g$coefficients, 10 * f$coefficients, tolerance = 2e-3)
	expect_equal(g$VC[["ve"]], 100 * f$VC[["ve"]], tolerance = 2e-3)
	expect_equal(g$VC[["va"]], 100 * f$VC[["va"]], tolerance = 2e-3)
	expect_equal(g$sigma_ab, 10 * f$sigma_ab, tolerance = 2e-3)
	expect_equal(g$rho_ab, f$rho_ab, tolerance = 2e-3)
	expect_equal(g$VC[["rho"]], f$VC[["rho"]], tolerance = 2e-3)
})

test_that("the floors of the variances are relative to the response: tiny and huge units fit as the unit response does", {
	skip_on_cran()
	d = nrm_sim(n = 14, Tt = 4, seed = 41, rho_d = 0.3)
	f = nrm_fit(d)
	expect_true(f$converged)
	for (s in c(1e-6, 1e6)) {
		g = nrm_fit(nrm_sim(n = 14, Tt = 4, seed = 41, rho_d = 0.3, scale = s))
		expect_true(g$converged)
		# slopes and variances only: the coefficient prior g = n^2 max(1, var(y)) follows the unit only above one,
		# so the intercept of a tiny-unit response differs slightly
		expect_equal(unname(g$coefficients[-1]), s * unname(f$coefficients[-1]), tolerance = 2e-3)
		expect_equal(g$VC[["ve"]], s^2 * f$VC[["ve"]], tolerance = 2e-3)
		expect_equal(g$VC[["va"]], s^2 * f$VC[["va"]], tolerance = 5e-3)
		expect_equal(g$sigma_ab, s * f$sigma_ab, tolerance = 5e-3)
		expect_equal(g$rho_ab, f$rho_ab, tolerance = 5e-3)
		expect_equal(g$VC[["rho"]], f$VC[["rho"]], tolerance = 2e-3)
	}
})

test_that("the prior of the factor covariance reads the shared definition of R/uv_prior.R", {
	d = nrm_sim(n = 10, Tt = 3, seed = 61)
	ctx = nrm_ctx(d, R = 1L)
	prep = asNamespace("lame")$.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet"
	)
	expect_equal(ctx$suv0, asNamespace("lame")$.lame_uv_prior_scale(prep$Yarr, "normal"))
	# and the sampler's typical additive effect size: root mean variance of the row and column means of the centred response
	Ec = prep$Yarr - mean(prep$Yarr, na.rm = TRUE)
	expect_equal(
		ctx$suv0,
		sqrt(mean(c(stats::var(apply(Ec, 1L, mean, na.rm = TRUE)), stats::var(apply(Ec, 2L, mean, na.rm = TRUE)))))
	)
})

test_that("the methods of a normal empirical-Bayes fit work, a refit stops at once and the prior of the residual variance is the sampler's", {
	d = nrm_sim(n = 12, Tt = 4, seed = 51, rho_d = 0.3)
	f = nrm_fit(d)
	expect_s3_class(f, "lame_dynamic_als")
	expect_output(print(f))
	expect_output(print(summary(f)))
	expect_equal(dim(vcov(f)), c(3L, 3L))
	ci = confint(f)
	expect_true(all(ci[, 1] < f$coefficients & f$coefficients < ci[, 2]))
	expect_equal(predict(f)[[1]], f$EZ[[1]])
	expect_equal(fitted(f)[[1]], f$EZ[[1]])
	expect_length(simulate(f, nsim = 1, seed = 1)$Y, 1L)
	g = suppressWarnings(ame_als_refit(f, verbose = FALSE))
	expect_identical(g$iterations, 0L)
	expect_equal(g$VC[["ve"]], f$VC[["ve"]], tolerance = 1e-6)
	expect_equal(g$coefficients, f$coefficients, tolerance = 1e-6)
	expect_equal(unname(f$eb$ve_prior[["s20"]]), stats::var(unlist(f$Y), na.rm = TRUE))
	f5 = nrm_fit(d, prior = list(s20 = 25))
	expect_equal(unname(f5$eb$ve_prior[["s20"]]), 25)
	expect_gt(f5$VC[["ve"]], f$VC[["ve"]])
})

test_that("static factors and moving positions of a normal panel are fitted", {
	skip_on_cran()
	d = nrm_sim(n = 16, Tt = 4, seed = 61, R = 1, suv = 1.5)
	f = nrm_fit(d, R = 1)
	expect_true(f$converged)
	expect_equal(dim(f$U), c(16L, 1L))
	expect_gt(f$VC[["ve"]], 0.6)
	expect_lt(f$VC[["ve"]], 1.4)
	# the product explains what the additive effects do not
	f0 = nrm_fit(d)
	expect_lt(f$VC[["ve"]], f0$VC[["ve"]])
	dm = nrm_sim(n = 14, Tt = 4, seed = 62, R = 1, dyn_uv = TRUE)
	g = nrm_fit(dm, R = 1, dynamic_uv = TRUE)
	expect_true(g$converged)
	expect_equal(dim(g$U), c(14L, 1L, 4L))
	expect_true(all(is.finite(g$coefficients)))
})

test_that("the ridge of a failed cholesky factor is relative to the information of a cell of unit variance", {
	# a normal cell has information 1 / var(y), so the ridge on a singular matrix scales with the unit of the response
	S = matrix(1, 2L, 2L) * 1e-6
	ch_unit = lame:::.lde_chol(S)
	ch_rel = lame:::.lde_chol(S, unit = 1e-6)
	expect_equal(max(abs(crossprod(ch_unit) - S)), 1e-10, tolerance = 1e-6)
	expect_equal(max(abs(crossprod(ch_rel) - S)), 1e-16, tolerance = 1e-3)
	# a factorisable matrix is untouched at any unit
	P = diag(c(2, 3))
	expect_identical(lame:::.lde_chol(P, unit = 1e-6), chol(P))
	# the unit is the response's for a normal panel and one for other families
	expect_equal(lame:::.lde_cunit(list(prior = list(vunit = 1e6))), 1e-6)
	expect_identical(lame:::.lde_cunit(list(prior = list(vunit = 1))), 1)
	expect_identical(lame:::.lde_cunit(list(prior = list())), 1)
	d = nrm_sim(n = 8, Tt = 3, seed = 3)
	f = nrm_fit(d)
	expect_true(f$converged)
})
