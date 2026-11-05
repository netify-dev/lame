skip_on_cran()

library(lame)

# variational estimator of the binary ALS fit: a Gaussian variational posterior of actor effects and latent factors,
# with the variance of the product carried through the likelihood

# a small binary panel of the three kinds, with known effects and a rank-R product
vb_panel = function(
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

vb_fit = function(d, R = 1, ...) {
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

# the posterior mode with the Laplace correction, the stage the variational fit starts from
vb_mode = function(d, R = 1) {
	md = if (d$kind == "bipartite") "bipartite" else "unipartite"
	Y = if (d$Tt == 1L) d$Y[[1]] else d$Y
	X = if (d$Tt == 1L) d$X[[1]] else d$X
	suppressWarnings(suppressMessages(lame:::.ae_run(
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
		call = quote(vb_mode()),
		lowrank_method = "mm",
		non_normal_method = "irls",
		link = "probit",
		linear_solver = "eigen",
		multistart = "none",
		bias_correct = NULL,
		eb_posterior = "mode"
	)))
}

# the engine on a mode fit: the state, the context and the iteration record
vb_engine = function(d, R = 1, ...) {
	f = vb_mode(d, R)
	st = lame:::.ae_eb_state(f)
	lame:::.ae_mf_fit(st$prep, st$fit, st$ab, st$R, d$kind == "symmetric", "probit", st$w, ...)
}

test_that("a variational fit records its estimator, its iteration and its lower bound", {
	skip_on_cran()
	for (kind in c("directed", "bipartite", "symmetric")) {
		d = vb_panel(kind, seed = 3)
		fit = vb_fit(d)
		expect_identical(fit$eb$estimator, "variational", info = kind)
		v = fit$eb$variational
		expect_true(v$converged, info = kind)
		expect_identical(v$termination, "converged", info = kind)
		expect_true(
			all(is.finite(unlist(v[c("residual", "step", "elbo", "ell", "kl_actors", "kl_factors")]))),
			info = kind
		)
		expect_identical(v$start, "posterior mode")
		# rank evidence is the Laplace marginal likelihood of the mode; the variational lower bound is kept beside it
		md = vb_mode(d)
		expect_equal(fit$eb$log_marginal, md$eb$log_marginal, tolerance = 1e-10, info = kind)
		expect_equal(fit$eb$log_marginal, v$mode_log_marginal, info = kind)
		expect_equal(lame:::.ae_eb_evidence(fit)$log_marginal, fit$eb$log_marginal, info = kind)
		expect_equal(lame:::.ae_mf_elbo_object(fit)$elbo, v$elbo, tolerance = 1e-6, info = kind)
		expect_lt(v$elbo, v$ell)
		expect_true(fit$converged, info = kind)
		fp = fit$eb$factor_posterior
		expect_identical(dim(fp$Sig), c(1L, 1L, as.integer(d$n)), info = kind)
		expect_true(all(fp$Sig > 0), info = kind)
		# the mode it starts from is the Laplace estimator
		expect_identical(vb_mode(d)$eb$estimator, "laplace", info = kind)
		expect_null(vb_mode(d)$eb$variational, info = kind)
	}
	# a fit without a latent rank is variational too, without factors
	fit0 = vb_fit(vb_panel("directed", seed = 4), R = 0)
	expect_identical(fit0$eb$estimator, "variational")
	expect_null(fit0$eb$factor_posterior)
	expect_true(fit0$eb$variational$converged)
})

test_that("the factor covariances of a directed fit are positive definite and the prior couples a sender and a receiver", {
	d = vb_panel("directed", n = 16, R = 2, seed = 5)
	fit = vb_fit(d, R = 2)
	fp = fit$eb$factor_posterior
	expect_identical(dim(fp$Sig), c(2L, 2L, 16L))
	for (i in seq_len(16)) {
		expect_true(all(eigen(fp$Sig[, , i], symmetric = TRUE, only.values = TRUE)$values > 0))
		expect_true(all(eigen(fp$Om[, , i], symmetric = TRUE, only.values = TRUE)$values > 0))
	}
	expect_identical(fp$prior$kind, "joint")
	expect_identical(dim(fit$eb$suv), c(4L, 4L))
	expect_true(all(eigen(fit$eb$suv, symmetric = TRUE, only.values = TRUE)$values > 0))
	# a bipartite fit has two independent sets of actors, so no coupling
	fb = vb_fit(vb_panel("bipartite", n = 16, R = 2, seed = 5), R = 2)
	expect_identical(fb$eb$factor_posterior$prior$kind, "diag")
	expect_null(fb$eb$suv)
})

test_that("the quadrature gives the expected score, weight and log likelihood of a gaussian predictor", {
	skip_on_cran()
	for (link in c("probit", "logit")) for (v in c(0, 0.3, 1)) for (m in c(-2, 0.4, 1.5)) for (yy in 0:1) {
		cl = lame:::.ae_mf_cells(array(yy, c(1, 1, 1)), array(m, c(1, 1, 1)), array(v, c(1, 1, 1)), link)
		sg = 2 * yy - 1
		psi = function(e) if (link == "probit") sg * exp(
			dnorm(e, log = TRUE) - pnorm(sg * e, log.p = TRUE)
		) else sg * (1 - plogis(sg * e))
		wt = function(e) if (link == "probit") { lam = exp(
			dnorm(e, log = TRUE) - pnorm(sg * e, log.p = TRUE)
		); lam * (lam + sg * e) } else plogis(e) * plogis(-e)
		ll = function(e) if (link == "probit") pnorm(sg * e, log.p = TRUE) else plogis(sg * e, log.p = TRUE)
		ex = function(f) if (v == 0) f(m) else integrate(
			function(e) f(e) * dnorm(e, m, sqrt(v)),
			-Inf,
			Inf,
			rel.tol = 1e-12
		)$value
		Es = ex(psi); Ew = ex(wt)
		info = paste(link, v, m, yy)
		expect_equal(as.numeric(cl$w), Ew, tolerance = 1e-4, info = info)
		expect_equal(as.numeric(cl$ell), ex(ll), tolerance = 1e-4, info = info)
		expect_equal(as.numeric(cl$z), m + Es / Ew, tolerance = 1e-4, info = info)
	}
	# the working increment is capped at 30, and the weight keeps the score
	cl = lame:::.ae_mf_cells(array(1, c(1, 1, 1)), array(-40, c(1, 1, 1)), array(0.1, c(1, 1, 1)), "logit")
	expect_equal(as.numeric(cl$z), -40 + 30)
	# unobserved cells carry nothing
	cl = lame:::.ae_mf_cells(array(NA_real_, c(1, 1, 1)), array(0.2, c(1, 1, 1)), array(0.1, c(1, 1, 1)), "probit")
	expect_identical(c(as.numeric(cl$z), as.numeric(cl$w), as.numeric(cl$ell)), c(0, 0, 0))
})

test_that("the variance of the product is the variance under the independent gaussian posteriors", {
	set.seed(1)
	R = 2; n = 3; N = 4e5
	U = matrix(rnorm(n * R), n); V = matrix(rnorm(n * R), n)
	psd = function() { A = matrix(rnorm(4), 2); crossprod(A) / 4 + diag(0.05, 2) }
	Sig = array(unlist(lapply(seq_len(n), function(i) psd())), c(2, 2, n))
	Om = array(unlist(lapply(seq_len(n), function(i) psd())), c(2, 2, n))
	pv = lame:::.ae_mf_prod_var(U, V, Sig, Om)
	draw = function(m, S) m + t(chol(S)) %*% matrix(rnorm(2 * N), 2)
	i = 1; j = 2
	pr = colSums(draw(U[i, ], Sig[, , i]) * draw(V[j, ], Om[, , j]))
	expect_equal(pv[i, j], var(pr), tolerance = 0.02)
	# an undirected network: the product of two rows through the sign pattern
	sg = c(1, -1)
	pvs = lame:::.ae_mf_prod_var_sym(U, sg, Sig)
	pr = colSums(draw(U[i, ], Sig[, , i]) * (sg * draw(U[j, ], Sig[, , j])))
	expect_equal(pvs[i, j], var(pr), tolerance = 0.02)
	expect_equal(pvs, t(pvs))
})

test_that("a variational fit is a stationary point of its objective and the mode is not", {
	d = vb_panel("directed", n = 16, R = 1, seed = 5)
	res = vb_engine(d, 1, tol = 1e-9, max_iter = 3000L)
	expect_true(res$converged)
	cx = res$cx
	st = lame:::.ae_eb_state(vb_mode(d, 1))
	S0 = lame:::.ae_mf_init(cx, st$fit, st$w)
	val = function(S) {
		ev = lame:::.ae_mf_evaluate(cx, S)
		lame:::.ae_mf_merit(cx, S, ev, S$sab)
	}
	sym2 = function(k) { A = matrix(rnorm(k * k), k); (A + t(A)) / 2 }
	# largest directional derivative over random directions of the means and covariances, by central differences
	slope = function(S, what, seed) {
		set.seed(seed)
		eps = 1e-5
		bump = function(S, e) {
			if (what == "mean") {
				S$mu = S$mu + e * dm[1]; S$beta = S$beta + e * dm[2]
				S$a = S$a + e * da; S$b = S$b + e * db; S$U = S$U + e * dU; S$V = S$V + e * dV
			} else {
				S$Sc = S$Sc + e * dSc; S$vact = lame:::.ae_mf_vact(S$Sc, cx$ab)
				S$Sig = S$Sig + e * dSig; S$Om = S$Om + e * dOm
			}
			S
		}
		dm = rnorm(2); da = rnorm(16); db = rnorm(16); dU = matrix(rnorm(16), 16); dV = matrix(rnorm(16), 16)
		dSc = 0.01 * sym2(32); dSig = array(0.01 * rnorm(16), c(1, 1, 16)); dOm = array(0.01 * rnorm(16), c(1, 1, 16))
		(val(bump(S, eps)) - val(bump(S, -eps))) / (2 * eps)
	}
	for (what in c("mean", "cov")) {
		g_vb = max(abs(vapply(1:4, function(s) slope(res$S, what, s), 0)))
		g_mode = max(abs(vapply(1:4, function(s) slope(S0, what, s), 0)))
		cat(what, g_vb, g_mode, "\n")
		expect_lt(g_vb, 1e-3 * g_mode, label = paste("gradient at the variational fit,", what))
		expect_gt(g_mode, 1e-3)
	}
})

test_that("the divergences in the lower bound are those of the gaussians they name", {
	kl = function(m1, S1, S0) {
		k = length(m1)
		0.5 * (sum(diag(solve(S0) %*% S1)) + as.numeric(t(m1) %*% solve(S0) %*% m1) - k +
				as.numeric(
					determinant(S0, logarithm = TRUE)$modulus
				) - as.numeric(determinant(S1, logarithm = TRUE)$modulus))
	}
	# a directed fit: the rows (u_i, v_i) against N(0, Suv), the sender and the receiver independent under q
	d = vb_panel("directed", n = 14, R = 2, seed = 19)
	res = vb_engine(d, 2, tol = 1e-7, max_iter = 2000L)
	S = res$S; Suv = S$prior$Suv
	by_hand = sum(vapply(seq_len(14), function(i) {
		Q = matrix(0, 4, 4); Q[1:2, 1:2] = S$Sig[, , i]; Q[3:4, 3:4] = S$Om[, , i]
		kl(c(S$U[i, ], S$V[i, ]), Q, Suv)
	}, 0))
	expect_equal(lame:::.ae_mf_kl_factors(S, 2L), by_hand, tolerance = 1e-10)
	# a bipartite fit: per-column scales
	db = vb_panel("bipartite", n = 14, R = 2, seed = 19)
	rb = vb_engine(db, 2, tol = 1e-7, max_iter = 2000L)
	Sb = rb$S
	by_hand = sum(vapply(seq_len(dim(Sb$Sig)[3]), function(i) kl(Sb$U[i, ], Sb$Sig[, , i], diag(Sb$prior$su)), 0)) +
		sum(vapply(seq_len(dim(Sb$Om)[3]), function(j) kl(Sb$V[j, ], Sb$Om[, , j], diag(Sb$prior$sv)), 0))
	expect_equal(lame:::.ae_mf_kl_factors(Sb, 2L), by_hand, tolerance = 1e-10)
	# an undirected fit: rows against N(0, tau I)
	ds = vb_panel("symmetric", n = 16, R = 2, seed = 19)
	rs = vb_engine(ds, 2, tol = 1e-7, max_iter = 2000L)
	Ss = rs$S
	by_hand = sum(vapply(seq_len(16), function(i) kl(Ss$F[i, ], Ss$Sig[, , i], diag(Ss$prior$tau, 2)), 0))
	expect_equal(lame:::.ae_mf_kl_factors(Ss, 2L), by_hand, tolerance = 1e-10)
	# the effects
	ab = res$cx$ab
	Lam = lame:::.ae_eb_precision(S$sab, ab)
	expect_equal(
		lame:::.ae_mf_kl_actors(S$Sc, c(S$a, S$b), 0, Lam),
		kl(c(S$a, S$b), S$Sc, solve(Lam)),
		tolerance = 1e-8
	)
	# the lower bound is the expected log likelihood (nu = 1/2 on the ordered cells of an undirected network) less both
	el = res$elbo
	expect_equal(el$elbo, el$ell - el$kl_actors - el$kl_factors)
	expect_equal(rs$elbo$ell, 0.5 * rs$ev$ell)
})

test_that("the fixed point does not depend on the damping or on the number of sweeps", {
	skip_on_cran()
	d = vb_panel("directed", n = 16, R = 2, seed = 6)
	r1 = vb_engine(d, 2, tol = 1e-8, max_iter = 4000L)
	r2 = vb_engine(d, 2, tol = 1e-8, max_iter = 4000L, step_max = 0.4, inner = 2L)
	expect_true(r1$converged && r2$converged)
	expect_equal(r2$S$beta, r1$S$beta, tolerance = 1e-5)
	expect_equal(r2$elbo$elbo, r1$elbo$elbo, tolerance = 1e-6)
	off = row(r1$ev$Pm) != col(r1$ev$Pm)
	expect_lt(max(abs(r2$ev$Pm - r1$ev$Pm)[off]), 5e-4)
	# damping takes more iterations to reach the same point
	expect_gt(r2$iterations, r1$iterations)
})

test_that("rescaling a factor column leaves the likelihood alone and the closed-form gauge maximises the objective along it", {
	for (kind in c("bipartite", "directed")) {
		d = vb_panel(kind, n = 18, R = 2, seed = 7)
		res = vb_engine(d, 2, tol = 1e-7, max_iter = 2000L)
		cx = res$cx; S = res$S
		# scale the first column by a factor that unbalances the two sides
		unbal = function(S, c1) {
			cc = c(c1, 1)
			S$U = sweep(S$U, 2L, cc, "*"); S$V = sweep(S$V, 2L, 1 / cc, "*")
			S$Sig = S$Sig * array(tcrossprod(cc), dim(S$Sig)); S$Om = S$Om * array(tcrossprod(1 / cc), dim(S$Om))
			S
		}
		Su = unbal(S, 3)
		e0 = lame:::.ae_mf_evaluate(cx, S); e1 = lame:::.ae_mf_evaluate(cx, Su)
		expect_equal(e1$ell, e0$ell, tolerance = 1e-10)
		expect_equal(e1$Pm, e0$Pm, tolerance = 1e-10)
		term = function(S) lame:::.ae_mf_factor_term(S, 2L)
		Sg = lame:::.ae_mf_gauge(cx, Su)
		grid = vapply(exp(seq(-1.5, 1.5, length.out = 41)), function(c1) term(unbal(Su, c1)), 0)
		expect_gte(term(Sg), max(grid) - 1e-6)
		# and the reported fit sits at the gauge optimum already
		expect_equal(term(lame:::.ae_mf_gauge(cx, S)), term(S), tolerance = 1e-6)
	}
})

test_that("a refinement that cannot converge is refused with the convergence class, the state and the remedy", {
	d = vb_panel("directed", n = 20, R = 1, seed = 2)
	# the mode does not settle in two steps
	err = tryCatch(
		suppressMessages(ame_als(d$Y[[1]], Xdyad = d$X[[1]], R = 1, family = "binary", verbose = FALSE, max_iter = 2)),
		error = function(e) e
	)
	expect_s3_class(err, "lame_als_convergence")
	expect_match(conditionMessage(err), "posterior mode did not reach")
	expect_match(conditionMessage(err), "smaller rank")
	# the refinement itself: a record that stopped at its iteration limit is refused with its state
	rec = list(termination = "iteration limit", iterations = 40L, residual = 3.2e-3, tolerance = 1e-6)
	err2 = tryCatch(lame:::.ae_mf_refuse(rec, 3L), error = function(e) e)
	expect_s3_class(err2, "lame_als_convergence")
	msg = conditionMessage(err2)
	expect_match(msg, "variational refinement of the latent factors did not converge")
	expect_match(msg, "40 iterations")
	expect_match(msg, "rank 3")
	expect_match(msg, "max_iter")
	rec$termination = "no improving step"
	expect_match(
		conditionMessage(tryCatch(lame:::.ae_mf_refuse(rec), error = function(e) e)),
		"no damped step improved"
	)
})

test_that("the effects of a variational fit sum to zero without node covariates and the predictor is the fitted one", {
	d = vb_panel("directed", n = 14, R = 1, seed = 8)
	fit = vb_fit(d)
	expect_lt(abs(sum(fit$a)) + abs(sum(fit$b)), 1e-7)
	pr = lame:::.ae_eb_predictive(fit)
	off = row(fit$EZ[[1]]) != col(fit$EZ[[1]])
	expect_equal(unname(pr$mean[, , 1][off]), unname(fit$EZ[[1]][off]), tolerance = 1e-10)
	# the variance of the predictor contains the variance of the product
	post = fit$eb$factor_posterior
	st = lame:::.ae_product_restore(fit)
	vp = lame:::.ae_mf_prod_var(st$U, st$V, post$Sig, post$Om)
	expect_true(all(pr$var[, , 1][off] > vp[off]))
	expect_equal(
		unname(fitted(fit)[off]),
		unname(pnorm(pr$mean[, , 1] / sqrt(1 + pr$var[, , 1]))[off]),
		tolerance = 1e-10
	)
	expect_equal(unname(fitted(fit, plug_in = TRUE)[off]), unname(pnorm(fit$EZ[[1]])[off]), tolerance = 1e-10)
})

test_that("an undirected variational fit keeps the sign pattern of the mode and a symmetric product", {
	skip_on_cran()
	d = vb_panel("symmetric", n = 24, R = 2, seed = 9)
	fit = vb_fit(d, R = 2)
	md = vb_mode(d, 2)
	expect_identical(fit$eb$estimator, "variational")
	expect_identical(sign(fit$L), sign(md$L))
	expect_true(isSymmetric(unname(fit$UVPM), tol = 1e-8))
	expect_equal(unname(fit$UVPM), unname(fit$U %*% (fit$L * t(fit$U))), tolerance = 1e-8)
	expect_equal(crossprod(fit$U), diag(2), tolerance = 1e-8, ignore_attr = TRUE)
	fp = fit$eb$factor_posterior
	expect_identical(fp$prior$kind, "columns")
	expect_equal(unname(fit$UVPM), unname(tcrossprod(sweep(fp$F, 2L, fp$sg, "*"), fp$F)), tolerance = 1e-8)
	# a rank with a negative eigenvalue
	set.seed(10)
	n = 24; nm = sprintf("a%02d", 1:n)
	u1 = rnorm(n, 0, 0.9); u2 = rnorm(n, 0, 0.9)
	O = tcrossprod(u1) - tcrossprod(u2)
	x = matrix(rnorm(n * n), n); x[lower.tri(x)] = t(x)[lower.tri(x)]
	e = matrix(rnorm(n * n), n); e[lower.tri(e)] = t(e)[lower.tri(e)]
	y = 1 * (-0.3 + 0.7 * x + O + e > 0); diag(y) = NA; dimnames(y) = list(nm, nm)
	dn = list(
		Y = list(t1 = y),
		X = list(t1 = array(x, c(n, n, 1), dimnames = list(nm, nm, "x_dyad"))),
		n = n,
		nc = n,
		Tt = 1L,
		kind = "symmetric"
	)
	fn = vb_fit(dn, R = 2)
	expect_identical(sort(sign(fn$L)), c(-1, 1))
	expect_identical(sign(fn$L), sign(vb_mode(dn, 2)$L))
})

test_that("the variational and the Laplace fit agree without a latent rank", {
	for (kind in c("directed", "bipartite", "symmetric")) {
		d = vb_panel(kind, n = 30, R = 0, seed = 11)
		vb = vb_fit(d, R = 0); md = vb_mode(d, 0)
		expect_equal(unname(coef(vb)), unname(coef(md)), tolerance = 0.06, info = kind)
		expect_equal(unname(vb$eb$sab), unname(md$eb$sab), tolerance = 0.15, info = kind)
	}
})

test_that("a refit of a variational fit starts from the fit and agrees with it", {
	skip_on_cran()
	d = vb_panel("directed", n = 18, R = 1, seed = 12)
	fit = vb_fit(d)
	re = ame_als_refit(fit, verbose = FALSE)
	expect_identical(re$eb$estimator, "variational")
	expect_identical(re$eb$variational$start, "previous variational fit")
	expect_lt(re$eb$variational$iterations, 12L)
	expect_equal(coef(re), coef(fit), tolerance = 1e-4)
	expect_equal(re$eb$log_marginal, fit$eb$log_marginal, tolerance = 1e-6)
	expect_equal(re$eb$variational$elbo, fit$eb$variational$elbo, tolerance = 1e-6)
})

test_that("the parametric bootstrap refits with the variational estimator", {
	skip_on_cran()
	d = vb_panel("directed", n = 20, R = 1, seed = 13)
	fit = vb_fit(d)
	bt = suppressWarnings(suppressMessages(ame_als_bootstrap(fit, R = 10, seed = 1, verbose = FALSE)))
	expect_gte(bt$n_valid, 10L)
	expect_true(all(is.finite(bt$se)))
	ratio = bt$se[c("intercept", "x_dyad")] / fit$eb$se[c("intercept", "x_dyad")]
	expect_true(all(ratio > 0.4 & ratio < 2.5))
})

test_that("the lower bound prefers the rank the data carry", {
	skip_on_cran()
	d = vb_panel("directed", n = 40, R = 2, seed = 14, sd_uv = 1)
	lb = vapply(0:3, function(r) vb_fit(d, R = r)$eb$log_marginal, 0)
	expect_gt(lb[3], lb[1] + 5)
	expect_gt(lb[3], lb[2])
	# the extra column of rank 3 is pruned: it does not raise the bound by much
	expect_lt(lb[4] - lb[3], 2)
	# no rank in the data: the bound does not reward a rank
	d0 = vb_panel("directed", n = 40, R = 0, seed = 15)
	lb0 = vapply(0:2, function(r) vb_fit(d0, R = r)$eb$log_marginal, 0)
	expect_lt(max(lb0[2:3]) - lb0[1], 2)
})

test_that("a rank above the generating one does not inflate the slope as the mode does", {
	skip_on_cran()
	d = vb_panel("directed", n = 40, R = 0, seed = 16, mu = -1)
	vb = vb_fit(d, R = 2); md = vb_mode(d, 2); r0 = vb_fit(d, R = 0)
	expect_lt(abs(coef(vb)["x_dyad"] - coef(r0)["x_dyad"]), 0.75 * abs(coef(md)["x_dyad"] - coef(r0)["x_dyad"]))
})

test_that("a longitudinal and a bipartite variational fit predict and summarise", {
	for (kind in c("directed", "bipartite")) {
		d = vb_panel(kind, n = 14, Tt = 3, R = 1, seed = 17)
		fit = vb_fit(d)
		expect_identical(fit$eb$estimator, "variational")
		pr = fitted(fit)
		expect_length(pr, 3L)
		expect_true(
			all(vapply(pr, function(m) all(is.finite(m[!is.na(m)]) & m[!is.na(m)] > 0 & m[!is.na(m)] < 1), TRUE))
		)
		msg = paste(utils::capture.output(print(summary(fit))), collapse = "")
		expect_match(msg, "variational", ignore.case = TRUE)
		expect_true(all(is.finite(summary(fit)$coefficients[, "Std.Error"])))
	}
})

test_that("the mode is available as the estimator and is unchanged by the refinement code", {
	d = vb_panel("directed", n = 14, R = 1, seed = 18)
	md = vb_mode(d, 1)
	expect_identical(md$eb$estimator, "laplace")
	expect_true(isTRUE(md$irls_diagnostics$pass))
	vb = vb_fit(d, R = 1)
	# the certificate of the variational fit is the one of the mode it started from
	expect_true(isTRUE(vb$irls_diagnostics$pass))
	expect_equal(vb$irls_diagnostics$max_score, md$irls_diagnostics$max_score)
	# the mode and the variational fit differ in the product they report
	expect_gt(max(abs(lame:::.ae_product_restore(vb)$UVPM - md$UVPM), na.rm = TRUE), 1e-4)
})
