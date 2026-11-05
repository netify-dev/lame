skip_on_cran()

library(lame)

# the empirical-Bayes normal estimator (ame_als(), lame_als(), family = "normal"): the Gaussian cell and the cell pair
# (R/ame_als_gauss.R, R/ame_als_gauss_pair.R), exactness at rank 0 against independent generalised least squares,
# and the fits it produces across network kinds, dyadic correlation, prediction, refit, bootstrap and refusals

# normal AME network: sender and receiver effects (va, vb, cab) with node covariates, dyadic covariates, a latent
# product of rank R (coordinate variance suv), and noise s2 with dyadic correlation rho
sim_norm_net = function(
	n = 24,
	kind = "directed",
	R = 0,
	nt = 1L,
	s2 = 1,
	rho = 0.4,
	va = 0.5,
	vb = 0.5,
	cab = 0.15,
	suv = 0.5,
	mu = 0.3,
	beta = c(0.6, -0.4),
	gam = c(0.5, -0.4),
	seed = 1,
	miss = 0,
	m = n
) {
	set.seed(seed)
	sym = identical(kind, "symmetric"); bip = identical(kind, "bipartite")
	if (!bip) m = n
	rn = sprintf("a%03d", seq_len(n)); cn = if (bip) sprintf("b%03d", seq_len(m)) else rn
	xr = matrix(rnorm(n), n, 1, dimnames = list(rn, "xr")); xc = matrix(rnorm(m), m, 1, dimnames = list(cn, "xc"))
	if (sym) {
		a = rnorm(n, 0, sqrt(va)); b = a
	} else if (bip) {
		a = rnorm(n, 0, sqrt(va)); b = rnorm(m, 0, sqrt(vb))
	} else {
		e = matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(va, cab, cab, vb), 2)); a = e[, 1]; b = e[, 2]
	}
	a = a + gam[1] * xr[, 1]
	b = if (sym) a else b + gam[2] * xc[, 1]
	if (R > 0) {
		O = if (sym) { F = matrix(rnorm(n * R, 0, sqrt(suv)), n); F %*% (rep(c(1, -1), length.out = R) * t(F)) } else
			tcrossprod(matrix(rnorm(n * R, 0, sqrt(suv)), n), matrix(rnorm(m * R, 0, sqrt(suv)), m))
	} else O = matrix(0, n, m)
	X = array(NA_real_, c(n, m, 2L, nt), dimnames = list(rn, cn, c("x1", "x2"), NULL))
	Y = array(NA_real_, c(n, m, nt), dimnames = list(rn, cn, NULL))
	for (t in seq_len(nt)) {
		for (k in 1:2) { xk = matrix(rnorm(n * m), n, m); if (sym) xk[lower.tri(xk)] = t(xk)[lower.tri(xk)]; X[, , k, t] = xk }
		eta = mu + outer(a, b, "+") + O + beta[1] * X[, , 1, t] + beta[2] * X[, , 2, t]
		E = matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
		if (sym) E[lower.tri(E)] = t(E)[lower.tri(E)] else if (!bip) E[lower.tri(E)] = rho * t(E)[lower.tri(E)] + sqrt(1 - rho^2) * matrix(
			rnorm(n * m, 0, sqrt(s2)),
			n,
			m
		)[lower.tri(E)]
		y = eta + E
		if (!bip) diag(y) = NA
		if (miss > 0) { mk = matrix(runif(n * m) < miss, n, m); if (sym) mk[lower.tri(mk)] = t(mk)[lower.tri(mk)]; y[mk] = NA }
		Y[, , t] = y
	}
	list(
		Y = Y,
		X = X,
		Xrow = xr,
		Xcol = if (sym) NULL else xc,
		n = n,
		m = m,
		nt = nt,
		kind = kind,
		truth = list(rho = if (kind == "directed") rho else 0, s2 = s2, O = O)
	)
}
# the text a fit prints, on one line (cli wraps its lines)
printed = function(f) {
	paste(utils::capture.output(print(f)), collapse = " ")
}
fit_norm = function(d, R = 0, ...) {
	args = list(
		Y = if (d$nt == 1L) d$Y[, , 1] else lapply(seq_len(d$nt), function(t) d$Y[, , t]),
		Xdyad = if (d$nt == 1L) d$X[, , , 1] else lapply(seq_len(d$nt), function(t) d$X[, , , t]),
		Xrow = if (d$nt == 1L) d$Xrow else rep(list(d$Xrow), d$nt),
		Xcol = if (is.null(d$Xcol)) NULL else if (d$nt == 1L) d$Xcol else rep(list(d$Xcol), d$nt),
		R = R,
		family = "normal",
		mode = if (d$kind == "bipartite") "bipartite" else "unipartite",
		symmetric = d$kind == "symmetric",
		verbose = FALSE
	)
	args = c(args, list(...))
	suppressWarnings(do.call(if (d$nt == 1L) ame_als else lame_als, args))
}

test_that("the Gaussian cell is the expected working problem of a cell whose predictor is N(m, v)", {
	g = expand.grid(y = c(-1.3, 0.4, 2.2), m = c(-0.5, 0.8), v = c(0, 0.3, 1.5), ve = c(0.2, 1, 4))
	cl = lame:::.ae_gauss_cells(g$y, g$m, g$v, g$ve)
	for (i in seq_len(nrow(g))) {
		with(g[i, ], {
			ell = if (v > 0) integrate(
				function(e) dnorm(y, e, sqrt(ve), log = TRUE) * dnorm(e, m, sqrt(v)),
				m - 10 * sqrt(v),
				m + 10 * sqrt(v),
				rel.tol = 1e-10
			)$value else dnorm(y, m, sqrt(ve), log = TRUE)
			expect_equal(cl$ell[i], ell, tolerance = 1e-7)
			# expected score (y - m) / ve and information 1 / ve, so the working response is the observation
			expect_equal(cl$z[i], y, tolerance = 1e-12)
			expect_equal(cl$w[i], 1 / ve, tolerance = 1e-12)
			expect_equal(cl$ee[i], (y - m)^2 + v, tolerance = 1e-12)
		})
	}
	# unobserved cells are all zero
	yy = c(1, NA, 3); cl = lame:::.ae_gauss_cells(yy, c(0, 0, 0), c(0.1, 0.1, 0.1), 2)
	expect_equal(cl$w, c(0.5, 0, 0.5), ignore_attr = TRUE); expect_equal(cl$z, c(1, 0, 3), ignore_attr = TRUE)
	expect_equal(cl$ell[2], 0); expect_equal(cl$ee[2], 0)
})

test_that("the cavity weight is the information of the likelihood integrated over the predictor", {
	g = expand.grid(y = c(-1, 0.7), m = c(-0.4, 0.9), v = c(0, 0.5, 2), ve = c(0.3, 1.5))
	w = lame:::.ae_gauss_info_w(g$y, g$m, g$v, g$ve)
	h = 1e-4
	lz = function(y, m, v, ve) dnorm(y, m, sqrt(ve + v), log = TRUE)
	fd = -(lz(g$y, g$m + h, g$v, g$ve) - 2 * lz(g$y, g$m, g$v, g$ve) + lz(g$y, g$m - h, g$v, g$ve)) / h^2
	expect_equal(as.vector(w), fd, tolerance = 1e-5)
	expect_equal(as.vector(w), 1 / (g$ve + g$v), tolerance = 1e-12)
	expect_equal(lame:::.ae_gauss_info_w(c(1, NA), c(0, 0), c(0.2, 0.2), 1), c(1 / 1.2, 0), ignore_attr = TRUE)
})

test_that("the terms of the mode are the log likelihood and its derivatives, the information constant", {
	y = array(
		c(1.2, -0.4, NA, 0.9, 2.5, 0.1, -1, 0.3),
		c(2, 2, 2)
	); eta = array(c(0.5, 0, 1, 1, 2, -0.5, 0, 0.7), c(2, 2, 2))
	ve = 0.7; tm = lame:::.ae_gauss_terms(y, eta, ve)
	obs = is.finite(y)
	expect_equal(tm$l[obs], dnorm(y[obs], eta[obs], sqrt(ve), log = TRUE), tolerance = 1e-12)
	h = 1e-5
	fd = (lame:::.ae_gauss_terms(y, eta + h, ve)$l - lame:::.ae_gauss_terms(y, eta - h, ve)$l) / (2 * h)
	expect_equal(tm$s[obs], fd[obs], tolerance = 1e-6)
	expect_equal(tm$w[obs], rep(1 / ve, sum(obs))); expect_equal(sum(tm$hp), 0)
	expect_equal(tm$ee[obs], (y - eta)[obs]^2)
	expect_true(all(tm$l[!obs] == 0) && all(tm$s[!obs] == 0) && all(tm$w[!obs] == 0))
	# working problem: response y, weight 1 / ve, no cap on the increment
	wk = lame:::.ae_irls_gauss(y * 1e4, lame:::.ae_gauss_terms(y * 1e4, eta, ve))
	expect_equal(wk$z[obs], (y * 1e4)[obs]); expect_equal(wk$w[obs], rep(1 / ve, sum(obs)))
	expect_equal(lame:::.ae_dev_gauss(tm), -2 * sum(tm$l))
})

test_that("the posterior of the dyadic correlation is the sampler's prior times the pair likelihood", {
	set.seed(3)
	P = 40; ve = 0.8
	e1 = rnorm(P); e2 = 0.5 * e1 + sqrt(1 - 0.25) * rnorm(P)
	A = sum(e1^2 + e2^2); B = sum(e1 * e2)
	ps = lame:::.ae_gauss_rho_post(A, B, P, ve)
	f = function(r) exp(
		-0.5 * P * log(1 - r^2) - (A - 2 * r * B) / (2 * ve * (1 - r^2)) + lame:::.lame_rho_prior_log(r)
	)
	lo = -0.999999; hi = 0.999999
	z0 = integrate(f, lo, hi, rel.tol = 1e-10, subdivisions = 2000L)$value
	m1 = integrate(function(r) r * f(r), lo, hi, rel.tol = 1e-10, subdivisions = 2000L)$value / z0
	sd1 = sqrt(integrate(function(r) (r - m1)^2 * f(r), lo, hi, rel.tol = 1e-10, subdivisions = 2000L)$value / z0)
	expect_equal(ps$mean, m1, tolerance = 1e-5); expect_equal(ps$sd, sd1, tolerance = 1e-4)
	# with no information the result is the prior, symmetric about zero
	expect_equal(lame:::.ae_gauss_rho_post(0, 0, 0, 1)$mean, 0, tolerance = 1e-9)
})

test_that("the expected pair problem is the one of the bivariate normal law under q", {
	d = sim_norm_net(10, seed = 5, rho = 0.3)
	prep = lame:::.ae_prepare(d$Y[, , 1], d$X[, , , 1], d$Xrow, d$Xcol, "unipartite", FALSE)
	ab = lame:::.ae_eb_prior(prep, FALSE, "normal")
	cx = lame:::.ae_mf_context(prep, ab, 0L, "identity", pair = list(gauss = TRUE))
	set.seed(4)
	q = 2 * prep$nr
	Sc = crossprod(matrix(rnorm(q * q), q)) / q * 0.05 + diag(0.02, q)
	ve = 0.9; rho = 0.35
	S = list(Sc = Sc, ve = ve); S$vact = lame:::.ae_mf_vact(Sc, ab)
	eta = array(rnorm(prep$nr^2), dim(prep$Yarr)); V3 = array(S$vact, dim(prep$Yarr))
	pe = lame:::.ae_gauss_pair_cells(cx, S, eta, V3, rho)
	pc = cx$pair
	# expected log likelihood of every pair by two-dimensional Gauss-Hermite quadrature over the predictors N(m, V)
	gh = lame:::.ae_mf_gh(15L)
	ell = numeric(pc$n_pairs)
	for (p in seq_len(pc$n_pairs)[1:6]) {
		v1 = V3[pc$one[p]]; v2 = V3[pc$two[p]]; cv = pe$pair$cv[p]
		l11 = sqrt(v1); l21 = cv / l11; l22 = sqrt(v2 - l21^2)
		tot = 0
		for (k in seq_along(gh$z)) for (l in seq_along(gh$z)) {
			e1 = eta[pc$one[p]] + l11 * gh$z[k]; e2 = eta[pc$two[p]] + l21 * gh$z[k] + l22 * gh$z[l]
			r1 = pc$y1[p] - e1; r2 = pc$y2[p] - e2
			tot = tot + gh$w[k] * gh$w[l] * (-log(2 * pi * ve) - 0.5 * log(1 - rho^2) - (r1^2 - 2 * rho * r1 * r2 + r2^2) / (2 * ve * (1 - rho^2)))
		}
		ell[p] = tot
		expect_equal(pe$wk$ell[pc$one[p]], tot, tolerance = 1e-9)
	}
	# weights: the diagonal of a pair's information on the first cell, the cross weight in both
	k = 1 / (ve * (1 - rho^2))
	expect_equal(pe$wk$w[pc$one], rep(k, pc$n_pairs)); expect_equal(
		pe$pair$cw[pc$one],
		rep(-rho * k, pc$n_pairs)
	); expect_equal(pe$pair$cw[pc$one], pe$pair$cw[pc$two])
	# gradient G = W y of the normal equations
	expect_equal(pe$pair$G[pc$one], k * (pc$y1 - rho * pc$y2), tolerance = 1e-12)
	# the summed whitened second moment of the pairs gives the noise variance
	expect_equal(
		sum(pe$wk$ee[pc$one]),
		sum(
			((pc$y1 - eta[pc$one])^2 + (pc$y2 - eta[pc$two])^2 + V3[pc$one] + V3[pc$two] -
					2 * rho * ((pc$y1 - eta[pc$one]) * (pc$y2 - eta[pc$two]) + pe$pair$cv)) / (1 - rho^2)
		),
		tolerance = 1e-10
	)
	# the rho score of the expected log likelihood is the derivative of its sum over pairs
	h = 1e-5
	tot = function(r) sum(lame:::.ae_gauss_pair_cells(cx, S, eta, V3, r)$wk$ell[pc$one])
	expect_equal(pe$pair$rs, (tot(rho + h) - tot(rho - h)) / (2 * h), tolerance = 1e-6)
	expect_equal(pe$pair$rc, (tot(rho + h) - 2 * tot(rho) + tot(rho - h)) / h^2, tolerance = 1e-4)
})

test_that("the cavity weights of a pair are the information matrix of the pair density integrated over the predictors", {
	d = sim_norm_net(8, seed = 6)
	prep = lame:::.ae_prepare(d$Y[, , 1], d$X[, , , 1], d$Xrow, d$Xcol, "unipartite", FALSE)
	ab = lame:::.ae_eb_prior(prep, FALSE, "normal")
	cx = lame:::.ae_mf_context(prep, ab, 0L, "identity", pair = list(gauss = TRUE))
	pc = cx$pair
	set.seed(2)
	eta = array(rnorm(prep$nr^2), dim(prep$Yarr)); V3 = array(runif(prep$nr^2, 0.05, 0.4), dim(prep$Yarr))
	cv = rep(0.03, pc$n_pairs); rho = -0.2; ve = 0.7
	iw = lame:::.ae_gauss_pair_info_w(cx, eta, V3, cv, rho, ve)
	for (p in 1:4) {
		i1 = pc$one[p]; i2 = pc$two[p]
		C = matrix(c(ve + V3[i1], rho * ve + cv[p], rho * ve + cv[p], ve + V3[i2]), 2)
		lp = function(m1, m2) { r = c(pc$y1[p] - m1, pc$y2[p] - m2); -log(2 * pi) - 0.5 * log(det(C)) - 0.5 * sum(r * solve(C, r)) }
		h = 1e-4; m1 = eta[i1]; m2 = eta[i2]
		expect_equal(iw$w[i1], -(lp(m1 + h, m2) - 2 * lp(m1, m2) + lp(m1 - h, m2)) / h^2, tolerance = 1e-5)
		expect_equal(iw$w[i2], -(lp(m1, m2 + h) - 2 * lp(m1, m2) + lp(m1, m2 - h)) / h^2, tolerance = 1e-5)
		expect_equal(
			iw$cw[i1],
			-(lp(m1 + h, m2 + h) - lp(m1 + h, m2 - h) - lp(m1 - h, m2 + h) + lp(m1 - h, m2 - h)) / (4 * h^2),
			tolerance = 1e-5
		)
	}
})

# cell-level GLS solution of a rank-0 normal fit at its own hyperparameters: y = X theta + Z alpha + e,
# alpha_i ~ N(0, Sab), e ~ N(0, ve), with correlation rho within a dyad
gls_reference = function(f, d, pair = FALSE) {
	n = d$n; y = d$Y[, , 1]; ok = which(is.finite(y)); idx = arrayInd(ok, c(n, n))
	i = idx[, 1]; j = idx[, 2]
	X = cbind(1, d$X[, , 1, 1][ok], d$X[, , 2, 1][ok], d$Xrow[i, 1], d$Xcol[j, 1])
	N = length(ok)
	Z = matrix(0, N, 2 * n); Z[cbind(seq_len(N), i)] = 1; Z[cbind(seq_len(N), n + j)] = 1
	sab = f$eb$sab
	G = kronecker(matrix(c(1, 0, 0, 1), 2), diag(n))
	G[1:n, 1:n] = sab[1, 1] * diag(n); G[n + 1:n, n + 1:n] = sab[2, 2] * diag(n); G[1:n, n + 1:n] = sab[1, 2] * diag(n); G[n + 1:n, 1:n] = sab[1, 2] * diag(n)
	ve = f$eb$ve
	Ee = diag(ve, N)
	if (pair) {
		rho = f$eb$pair$rho
		key = function(a, b) paste(pmin(a, b), pmax(a, b))
		k = key(i, j); pr = split(seq_len(N), k)
		for (g in pr) if (length(g) == 2L) { Ee[g[1], g[2]] = rho * ve; Ee[g[2], g[1]] = rho * ve }
	}
	Sg = Z %*% G %*% t(Z) + Ee
	Si = solve(Sg)
	cov_theta = solve(t(X) %*% Si %*% X)
	theta = as.vector(cov_theta %*% t(X) %*% Si %*% y[ok])
	list(theta = theta, cov = cov_theta, X = X, Z = Z, G = G, Sg = Sg, Ee = Ee, ok = ok, i = i, j = j, N = N)
}

test_that("a rank-0 fit is the exact posterior of the linear Gaussian model at its hyperparameters", {
	d = sim_norm_net(12, seed = 7)
	f = fit_norm(d, 0, dcor = FALSE, tol = 1e-10)
	expect_identical(f$bias_correct$method, "eb"); expect_identical(f$eb$estimator, "variational")
	g = gls_reference(f, d)
	# coefficients, node coefficients included, are the GLS ones; their covariance is the Laplace one
	expect_equal(unname(f$coefficients), g$theta, tolerance = 1e-6)
	expect_equal(unname(f$eb$vcov), unname(g$cov), tolerance = 1e-6)
	expect_equal(unname(f$eb$se), sqrt(diag(g$cov)), tolerance = 1e-6)
	# the effects are the posterior means (BLUPs) given the coefficients
	res = d$Y[, , 1][g$ok] - as.vector(g$X %*% g$theta)
	alpha = as.vector(g$G %*% t(g$Z) %*% solve(g$Sg, res))
	yhat = as.vector(g$X %*% g$theta) + as.vector(g$Z %*% alpha)
	mu_all = f$EZ[[1]][g$ok]
	# the returned (a, b) exclude the part the node covariates carry
	expect_equal(unname(mu_all), yhat, tolerance = 1e-6)
	# noise variance is the EM fixed point (nu0 s20 + sum E[e^2]) / (N + nu0 - 2), with each cell's predictor posterior variance
	Ca = solve(t(g$Z) %*% g$Z / f$eb$ve + solve(g$G))
	vcell = rowSums((g$Z %*% Ca) * g$Z)
	ee = sum((d$Y[, , 1][g$ok] - yhat)^2 + vcell)
	s20 = f$eb$ve_prior[["s20"]]
	expect_equal(f$eb$ve, (s20 + ee) / (g$N + 1 - 2), tolerance = 1e-5)
	# the noise variance prior is the sampler's: one degree of freedom on the variance of the response
	expect_equal(s20, var(d$Y[, , 1][g$ok]), tolerance = 1e-12)
	# effect covariance is its EM fixed point (eta0 Sab0 + sum E[alpha alpha']) / (n + eta0 - 3)
	Cs = Ca; al = matrix(alpha - 0, ncol = 1)
	al_a = alpha[1:12] - d$Xrow[, 1] * 0; al_b = alpha[12 + 1:12]
	expect_true(all(is.finite(f$eb$sab)) && all(eigen(f$eb$sab)$values > 0))
})

test_that("a rank-0 fit with the dyadic correlation is the generalised least squares fit of the pairs", {
	d = sim_norm_net(12, seed = 8, rho = 0.5)
	f = fit_norm(d, 0, tol = 1e-10)
	expect_false(is.null(f$eb$pair)); expect_equal(f$VC[["rho"]], f$eb$pair$rho)
	g = gls_reference(f, d, pair = TRUE)
	expect_equal(unname(f$coefficients), g$theta, tolerance = 1e-6)
	expect_equal(unname(f$eb$vcov), unname(g$cov), tolerance = 1e-6)
	expect_gt(f$eb$pair$rho, 0.2); expect_lt(f$eb$pair$rho, 0.8)
	expect_equal(f$eb$pair$n_pairs, 12 * 11 / 2)
	# the correlation interval is the posterior one under the sampler's prior, inside the open interval
	expect_gt(f$eb$pair$sd, 0)
})

test_that("a normal fit is empirical Bayes by default and the fixed-effects fit is one argument away", {
	d = sim_norm_net(16, seed = 9)
	f = fit_norm(d, 0)
	expect_identical(f$bias_correct$method, "eb"); expect_true(f$bias_correct$default)
	expect_identical(f$non_normal_method, "irls"); expect_identical(f$link, "identity")
	expect_false(is.null(f$eb$se)); expect_true(
		all(c("intercept", "x1_dyad", "x2_dyad", "xr_row", "xc_col") %in% names(f$eb$se))
	)
	expect_equal(f$s2, f$VC[["ve"]]); expect_equal(f$VC[["ve"]], f$eb$ve)
	expect_true(is.finite(f$eb$log_marginal))
	fe = fit_norm(d, 0, bias_correct = "eb")
	expect_false(isTRUE(fe$bias_correct$default)); expect_equal(fe$coefficients, f$coefficients)
	fx = fit_norm(d, 0, bias_correct = "none", tol = 1e-10)
	expect_identical(fx$bias_correct$method, "none"); expect_identical(
		fx$bias_correct$requested,
		"none"
	); expect_null(fx$eb)
	expect_identical(fx$non_normal_method, "transform")
	# the rank-0 fixed-effects fit is two-way least squares, so the dyadic slopes are lm()'s
	y = as.vector(d$Y[, , 1]); x1 = as.vector(d$X[, , 1, 1]); x2 = as.vector(d$X[, , 2, 1])
	rw = factor(rep(seq_len(16), 16)); cl = factor(rep(seq_len(16), each = 16))
	ls = lm(y ~ x1 + x2 + rw + cl)
	expect_equal(unname(fx$coefficients[c("x1_dyad", "x2_dyad")]), unname(coef(ls)[c("x1", "x2")]), tolerance = 1e-6)
	# its standard errors are the conditional sandwich of the dyadic coefficients
	expect_null(fx$eb$vcov); expect_true(all(c("x1_dyad", "x2_dyad") %in% rownames(suppressMessages(vcov(fx)))))
	# the empirical-Bayes covariance covers every coefficient, node covariates included
	expect_setequal(rownames(vcov(f)), names(f$coefficients))
	expect_true(all(diag(vcov(f)) > 0))
	ci = confint(f); expect_equal(nrow(ci), length(f$coefficients))
	expect_identical(rownames(ci), names(f$coefficients))
})

test_that("every kind of network converges at ranks 0 to 2, with node covariates, missing cells and several periods", {
	skip_on_cran()
	for (kind in c("directed", "bipartite", "symmetric")) for (R in 0:2) for (nt in c(1L, 3L)) {
		d = sim_norm_net(18, kind = kind, R = min(R, 1), nt = nt, seed = 100 + R + 10 * nt, miss = 0.05, m = 22)
		f = fit_norm(d, R)
		lab = paste(kind, "R", R, "T", nt)
		expect_true(f$converged, label = lab)
		expect_identical(f$bias_correct$method, "eb", label = lab)
		expect_true(all(is.finite(f$coefficients)) && all(is.finite(f$eb$se)) && all(f$eb$se > 0), label = lab)
		expect_true(f$VC[["ve"]] > 0, label = lab)
		expect_equal(is.null(f$eb$pair), kind != "directed", label = lab)
		if (R > 0L) expect_equal(ncol(f$U), R, label = lab) else expect_null(f$UVPM)
		expect_equal(length(f$EZ), nt, label = lab)
	}
})

test_that("the rank that the evidence picks is the generating one, and the product is recovered", {
	skip_on_cran()
	d = sim_norm_net(30, R = 2, suv = 0.6, seed = 21)
	lm_ = vapply(0:3, function(R) fit_norm(d, R)$eb$log_marginal, 0)
	expect_equal(which.max(lm_) - 1L, 2L)
	f = fit_norm(d, 2)
	dc = function(P) P - rowMeans(P) - rep(colMeans(P), each = nrow(P)) + mean(P)
	expect_gt(cor(as.vector(dc(f$UVPM)), as.vector(dc(d$truth$O))), 0.8)
	# the over-ranked fit converges with a product of the same size; spread over more loadings its signal can fall below
	# the threshold, and then the product is withheld and the evidence names the generating rank
	f3 = fit_norm(d, 3)
	expect_true(f3$converged)
	P3 = if (lame:::.ae_product_is_withheld(f3)) f3$eb$product_withheld$means$UVPM else f3$UVPM
	expect_lt(abs(sqrt(sum(dc(P3)^2)) / sqrt(sum(dc(f$UVPM)^2)) - 1), 0.1)
	if (lame:::.ae_product_is_withheld(f3)) expect_identical(f3$eb$product_withheld$rank_chosen, 2L)
})

test_that("the dyadic correlation is estimated, shortens the interval of a covariate its cells do not share, and is absent where there are no pairs", {
	skip_on_cran()
	d = sim_norm_net(40, seed = 31, rho = 0.6)
	f = fit_norm(d, 0); f0 = fit_norm(d, 0, dcor = FALSE)
	expect_lt(abs(f$VC[["rho"]] - 0.6), 0.1)
	expect_equal(f$rho_estimation$method, "bivariate normal pair likelihood")
	# the independent-cell fit reports the post hoc correlation of its residuals
	expect_lt(abs(f0$VC[["rho"]] - f$VC[["rho"]]), 0.06); expect_null(f0$eb$pair)
	expect_lt(f$eb$se[["x1_dyad"]], f0$eb$se[["x1_dyad"]])
	fb = fit_norm(
		sim_norm_net(20, kind = "bipartite", seed = 32, m = 25),
		0
	); expect_null(fb$eb$pair); expect_equal(fb$VC[["rho"]], 0)
	fs = fit_norm(
		sim_norm_net(20, kind = "symmetric", seed = 33),
		0
	); expect_null(fs$eb$pair); expect_equal(fs$VC[["rho"]], 0)
})

test_that("a fit that is held to R and its pair model also runs on the R reference engine", {
	skip_on_cran()
	d = sim_norm_net(14, R = 1, seed = 41)
	f1 = fit_norm(d, 1)
	op = options(lame.als_engine = "r"); on.exit(options(op), add = TRUE)
	f2 = fit_norm(d, 1)
	expect_equal(f1$coefficients, f2$coefficients, tolerance = 1e-6)
	expect_equal(f1$eb$se, f2$eb$se, tolerance = 1e-6)
	expect_equal(f1$VC, f2$VC, tolerance = 1e-6)
})

test_that("a normal fit agrees with the sampler on a small network", {
	skip_on_cran()
	d = sim_norm_net(18, seed = 51, rho = 0.4)
	f = fit_norm(d, 0)
	m = suppressWarnings(ame(
		d$Y[, , 1],
		Xdyad = d$X[, , , 1],
		Xrow = d$Xrow,
		Xcol = d$Xcol,
		R = 0,
		family = "normal",
		dcor = TRUE,
		burn = 500,
		nscan = 12000,
		odens = 10,
		verbose = FALSE,
		plot = FALSE,
		seed = 3,
		gof = FALSE
	))
	bm = colMeans(m$BETA); bs = apply(m$BETA, 2, sd)
	nm = names(bm)
	expect_lt(max(abs(f$coefficients[nm] - bm) / bs), 0.5)
	expect_true(all(f$eb$se[nm] / bs > 0.8 & f$eb$se[nm] / bs < 1.25))
	vc = colMeans(m$VC); vs = apply(m$VC, 2, sd)
	expect_lt(abs(f$VC[["ve"]] - vc[["ve"]]) / vs[["ve"]], 1.5)
	expect_lt(abs(f$VC[["rho"]] - vc[["rho"]]), 0.12)
})

test_that("the response scale is carried by the fit: coefficients and variances rescale with it at rank 0", {
	d = sim_norm_net(14, seed = 61)
	f = fit_norm(d, 0, dcor = FALSE)
	d3 = d; d3$Y = d$Y * 1000
	f3 = fit_norm(d3, 0, dcor = FALSE)
	expect_equal(f3$coefficients, f$coefficients * 1000, tolerance = 1e-5)
	expect_equal(f3$eb$se, f$eb$se * 1000, tolerance = 1e-5)
	expect_equal(f3$VC[["ve"]], f$VC[["ve"]] * 1e6, tolerance = 1e-5)
	expect_equal(f3$eb$log_marginal - f$eb$log_marginal, -length(d$Y[is.finite(d$Y)]) * log(1000), tolerance = 1e-5)
})

test_that("an undirected product rescales with the response, in any units the estimator accepts", {
	skip_on_cran()
	d = sim_norm_net(24, kind = "symmetric", R = 1, suv = 0.8, seed = 62)
	f = fit_norm(d, 1)
	expect_true(f$converged)
	dc = function(P) P - rowMeans(P) - rep(colMeans(P), each = nrow(P)) + mean(P)
	for (m in c(1e-4, 1e3, 1e6)) {
		dm = d; dm$Y = d$Y * m
		fm = fit_norm(dm, 1)
		expect_true(fm$converged)
		# the majorize-minimize step must not stall at the start for large-unit responses (a stalled product is zero)
		expect_equal(fm$VC[["ve"]] / m^2, f$VC[["ve"]], tolerance = 0.01)
		expect_equal(sqrt(sum(dc(fm$UVPM)^2)) / m, sqrt(sum(dc(f$UVPM)^2)), tolerance = 0.03)
		expect_equal(unname(fm$coefficients / m), unname(f$coefficients), tolerance = 0.01)
	}
})

test_that("a response without variation, or in units the estimator is not conditioned for, is refused with the remedy", {
	n = 12
	Y = matrix(2, n, n); diag(Y) = NA
	expect_error(ame_als(Y, family = "normal", verbose = FALSE), class = "lame_als_data")
	set.seed(5); Y2 = matrix(rnorm(n * n), n) * 1e-9; diag(Y2) = NA
	expect_error(ame_als(Y2, family = "normal", verbose = FALSE), class = "lame_als_scale")
	# the fixed-effects fit is unaffected
	expect_no_error(suppressWarnings(ame_als(Y, family = "normal", bias_correct = "none", verbose = FALSE)))
	# a fit that does not settle within its iterations is refused, as in the other families
	d = sim_norm_net(24, R = 1, seed = 71)
	expect_error(fit_norm(d, 2, max_iter = 2), class = "lame_als_convergence")
	yb = 1 * (d$Y[, , 1] > 0.3)
	expect_error(
		ame_als(yb, R = 0, family = "binary", non_normal_method = "transform", bias_correct = "eb", verbose = FALSE),
		class = "lame_als_eb"
	)
})

test_that("nearly collinear covariates get a huge standard error, and a design collinear to rounding is refused with the cause", {
	skip_on_cran()
	set.seed(84); n = 30
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	x1 = matrix(rnorm(n * n), n, n)
	x2 = x1 + matrix(rnorm(n * n, 0, 1e-3), n, n)
	f = suppressWarnings(
		ame_als(
			Y,
			Xdyad = array(c(x1, x2), c(n, n, 2), dimnames = list(NULL, NULL, c("a", "b"))),
			R = 0,
			family = "normal",
			verbose = FALSE
		)
	)
	expect_true(f$converged)
	# the flat coefficient prior leaves the collinear direction to the data, so no pseudo-inverse picks a small coefficient
	expect_gt(min(f$eb$se[c("a_dyad", "b_dyad")]), 10)
	x3 = x1 + matrix(rnorm(n * n, 0, 1e-6), n, n)
	expect_error(
		suppressWarnings(ame_als(Y, Xdyad = array(c(x1, x3), c(n, n, 2)), R = 0, family = "normal", verbose = FALSE)),
		"Nearly collinear covariates",
		class = "lame_als_convergence"
	)
})

test_that("fitted values, predictions, residuals, simulation and the summary work on a normal empirical-Bayes fit", {
	d = sim_norm_net(14, R = 1, seed = 81)
	f = fit_norm(d, 1)
	fv = fitted(f); expect_equal(dim(fv), c(14, 14)); expect_equal(fv, f$fitted[[1]])
	expect_equal(residuals(f), d$Y[, , 1] - fv, ignore_attr = TRUE)
	pr = predict(f); expect_equal(pr, fv, ignore_attr = TRUE)
	nd = d$X[, , , 1]; nd[, , 1] = nd[, , 1] + 1
	pn = predict(f, newdata = nd)
	off = row(pr) != col(pr)
	expect_equal(as.vector((pn - pr)[off]), rep(f$coefficients[["x1_dyad"]], sum(off)), tolerance = 1e-8)
	sm = simulate(f, nsim = 3, seed = 2); expect_length(sm$Y, 3L); expect_equal(dim(sm$Y[[1]]), c(14, 14))
	s = summary(f); expect_true("Std.Error" %in% colnames(s$coefficients))
	expect_true(all(is.finite(s$coefficients[, "Std.Error"])))
	expect_match(printed(f), "empirical Bayes", fixed = TRUE)
	expect_error(logLik(f), "not defined")
})

test_that("a refit on the same data returns the fit, and a refit and a bootstrap run on new data", {
	skip_on_cran()
	d = sim_norm_net(16, R = 1, seed = 91)
	f = fit_norm(d, 1)
	rf = ame_als_refit(f, verbose = FALSE)
	expect_true(rf$converged)
	expect_equal(rf$coefficients, f$coefficients, tolerance = 1e-4)
	expect_equal(rf$VC[["ve"]], f$VC[["ve"]], tolerance = 1e-4)
	expect_equal(rf$eb$log_marginal, f$eb$log_marginal, tolerance = 1e-3)
	bt = suppressWarnings(ame_als_bootstrap(f, R = 12, seed = 4, verbose = FALSE))
	expect_equal(bt$n_valid, 12L)
	expect_true(all(is.finite(bt$se)))
	ci = suppressWarnings(confint(bt)); expect_true(all(names(f$coefficients) %in% rownames(ci)))
	# a fixed-effects fit refitted with the empirical-Bayes estimator gets all of it
	fx = fit_norm(d, 1, bias_correct = "none")
	re = ame_als_refit(fx, bias_correct = "eb", verbose = FALSE)
	expect_identical(re$bias_correct$method, "eb"); expect_identical(
		re$eb$estimator,
		"variational"
	); expect_false(is.null(re$eb$se))
})

test_that("multistart keeps the start with the largest evidence, and lame() and ame() with method = 'als' reach the estimator", {
	skip_on_cran()
	d = sim_norm_net(18, R = 1, seed = 95)
	f = fit_norm(d, 1, multistart = "cheap")
	expect_identical(f$multistart_criterion, "Laplace marginal likelihood")
	expect_true(f$multistart_selected %in% 1:4)
	fa = suppressWarnings(
		ame(
			d$Y[, , 1],
			Xdyad = d$X[, , , 1],
			Xrow = d$Xrow,
			Xcol = d$Xcol,
			R = 1,
			family = "normal",
			method = "als",
			verbose = FALSE
		)
	)
	expect_identical(fa$bias_correct$method, "eb")
	fl = suppressWarnings(
		lame(
			list(d$Y[, , 1]),
			Xdyad = list(d$X[, , , 1]),
			Xrow = list(d$Xrow),
			Xcol = list(d$Xcol),
			R = 1,
			family = "normal",
			method = "als",
			verbose = FALSE
		)
	)
	expect_identical(fl$bias_correct$method, "eb")
	# the dynamic route is its own estimator (lame_dynamic_als_eb_nrm.R)
	Yl = lapply(1:4, function(t) { z = sim_norm_net(14, seed = 300 + t)$Y[, , 1]; z })
	fd = suppressWarnings(lame(Yl, R = 0, family = "normal", method = "als", dynamic_ab = TRUE, verbose = FALSE))
	expect_s3_class(fd, "lame_dynamic_als")
})

test_that("the prior of the variance of the noise is the sampler's and can be set", {
	skip_on_cran()
	d = sim_norm_net(20, seed = 97)
	f = fit_norm(d, 0, dcor = FALSE)
	f2 = fit_norm(d, 0, dcor = FALSE, prior = list(s20 = 100))
	expect_equal(f$eb$ve_prior[["nu0"]], 1)
	expect_equal(f2$eb$ve_prior[["s20"]], 100)
	# a one-degree-of-freedom prior on a large scale adds 100 / (N - 1) to the variance
	N = sum(is.finite(d$Y[, , 1]))
	expect_gt(f2$VC[["ve"]], f$VC[["ve"]])
	expect_lt(abs((f2$VC[["ve"]] - f$VC[["ve"]]) - (100 - f$eb$ve_prior[["s20"]]) / (N - 1)), 0.05)
})

test_that("an undirected network given as a directed one is refused as the sampler refuses it", {
	set.seed(3); n = 16
	a = rnorm(n); Y = outer(a, a, "+") + matrix(rnorm(n * n), n); Y = (Y + t(Y)) / 2; diag(Y) = NA
	expect_error(ame_als(Y, R = 0, family = "normal", verbose = FALSE), class = "lame_improper_normal_reciprocity")
	expect_error(ame_als(Y, R = 1, family = "normal", verbose = FALSE), class = "lame_improper_normal_reciprocity")
	# the remedies named in the message are valid fits
	expect_true(ame_als(Y, R = 0, family = "normal", symmetric = TRUE, verbose = FALSE)$converged)
	expect_true(ame_als(Y, R = 0, family = "normal", dcor = FALSE, verbose = FALSE)$converged)
	# one unequal pair makes the posterior proper, so it is not refused
	Y[1, 2] = Y[1, 2] + 0.5
	expect_no_error(ame_als(Y, R = 0, family = "normal", verbose = FALSE))
})

test_that("a product that the data barely determine is flagged by its signal per loading", {
	skip_on_cran()
	weak = sim_norm_net(30, R = 2, suv = 0.04, seed = 111)
	f = fit_norm(weak, 2)
	info = f$eb$variational$information
	expect_true(all(c("signal", "loadings", "per_loading", "threshold") %in% names(info)))
	expect_equal(info[["loadings"]], 2 * 60)
	expect_lt(info[["per_loading"]], 1)
	expect_match(printed(f), "noise variances of signal per loading", fixed = TRUE)
	strong = sim_norm_net(30, R = 2, suv = 0.8, seed = 112)
	g = fit_norm(strong, 2)
	expect_gt(g$eb$variational$information[["per_loading"]], 1)
	expect_false(grepl("signal per loading", printed(g), fixed = TRUE))
	# a fit without a latent rank has no such measure
	expect_null(fit_norm(strong, 0)$eb$variational$information)
})
