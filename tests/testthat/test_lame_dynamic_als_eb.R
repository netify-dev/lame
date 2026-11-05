skip_on_cran()

library(lame)

# a small binary panel from the lame dynamic_ab model: a level plus a zero-mean AR(1) deviation per actor, probit links
eb_sim = function(
	n = 14,
	nb = n,
	Tt = 4,
	mode = "directed",
	seed = 1,
	mu = -1,
	beta = c(0.5, -0.4),
	va = 0.5,
	vb = 0.5,
	rho = 0.7,
	sigma = 0.3,
	dcor = 0
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
		err = matrix(rnorm(n * nb), n, nb)
		# the two errors of a dyad are correlated (directed network; the draws exist only when dcor is set)
		if (dcor != 0 && !sym && !bip) {
			lo = lower.tri(err)
			err[lo] = dcor * t(err)[lo] + sqrt(1 - dcor^2) * matrix(rnorm(n * nb), n, nb)[lo]
		}
		z = eta + err
		if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
		y = 1 * (z > 0)
		dimnames(y) = list(rn, cn)
		if (!bip) diag(y) = NA
		Y[[t]] = y
		X[[t]] = xt
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

# the same with static factors: a rank R product of static positions added to the predictor (directed network);
# the draws follow the validation file's generator so that its panels (seed 31) are reproduced
eb_sim_fac = function(
	n = 20,
	Tt = 4,
	seed = 31,
	mu = -1.2,
	R = 1,
	suv = 0.5,
	beta = c(0.6, -0.5),
	va = 0.6,
	vb = 0.6,
	rho = 0.7,
	sigma = 0.35,
	dcor = 0
) {
	set.seed(seed)
	ar1 = function(m) {
		out = matrix(0, m, Tt)
		out[, 1] = rnorm(m, 0, sigma / sqrt(1 - rho^2))
		for (t in 2:Tt) out[, t] = rho * out[, t - 1] + rnorm(m, 0, sigma)
		out
	}
	alpha = rnorm(n, 0, sqrt(va))
	gam = rnorm(n, 0, sqrt(vb))
	d = ar1(n)
	e = ar1(n)
	U = matrix(rnorm(n * R, 0, sqrt(suv)), n, R)
	V = matrix(rnorm(n * R, 0, sqrt(suv)), n, R)
	rn = sprintf("a%03d", seq_len(n))
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		xt = array(
			rnorm(n * n * length(beta)),
			c(n, n, length(beta)),
			dimnames = list(rn, rn, paste0("x", seq_along(beta)))
		)
		eta = mu + outer(alpha + d[, t], gam + e[, t], "+")
		for (k in seq_along(beta)) eta = eta + beta[k] * xt[, , k]
		eta = eta + tcrossprod(U, V)
		err = matrix(rnorm(n * n), n, n)
		# the two errors of a dyad are correlated (the draws exist only when dcor is set)
		if (dcor != 0) {
			lo = lower.tri(err)
			err[lo] = dcor * t(err)[lo] + sqrt(1 - dcor^2) * matrix(rnorm(n * n), n, n)[lo]
		}
		y = 1 * (eta + err > 0)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
		X[[t]] = xt
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

eb_fit = function(d, ...) {
	suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		verbose = FALSE,
		...
	)))
}

# a state of a small directed panel with a missing cell in two periods (its partner is a cell without a pair), at which
# the iteration is not at a solution: random paths, a coefficient vector, a correlation of 0.45
pair_state = function(n = 9, Tt = 3, seed = 9) {
	d = eb_sim(n = n, Tt = Tt, seed = seed, beta = 0.5)
	d$Y[[1]][2, 5] = NA
	d$Y[[2]][7, 3] = NA
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	ctx = lame:::.lde_context(prep, FALSE, "probit", coef_prior = "flat", centre = FALSE, dcor = TRUE)
	st = lame:::.lde_init(ctx)
	set.seed(100 + seed)
	st$x[] = rnorm(length(st$x), 0, 0.3)
	st$lev = rowMeans(st$x)
	st$theta[1:2] = c(-1, 0.4)
	st$hyp$rho_d = 0.45
	list(ctx = ctx, st = st, prep = prep)
}

test_that("the block-tridiagonal solver reproduces the dense information matrix", {
	set.seed(3)
	Tt = 5L; qx = 6L; qw = 4L; qs = qx + qw; nl = qx; m = 3L; r = nl + m
	mk = function(k) crossprod(matrix(rnorm(k * k), k)) * 0.3 + diag(k) * 4
	A = lapply(seq_len(Tt), function(t) mk(qs))
	b = c(rep(0.9, qx), rep(0.5, qw))
	C = lapply(seq_len(Tt), function(t) matrix(rnorm(qs * r), qs, r) * 0.15)
	Hzz = crossprod(matrix(rnorm(r * r), r)) * 0.5 + diag(r) * 12
	fac = lame:::.lde_factor(A, b, C, Hzz, qx, nl)
	Hd = lame:::.lde_dense_matrix(A, b, C, Hzz)
	Hi = solve(Hd)
	blk = function(t) (t - 1L) * qs + seq_len(qs)
	zi = Tt * qs + seq_len(r)
	expect_equal(
		fac$logdet_xx + fac$logdet_z,
		as.numeric(determinant(Hd)$modulus),
		tolerance = 1e-10
	)
	rs = lapply(seq_len(Tt), function(t) rnorm(qs))
	rz = rnorm(r)
	sol = lame:::.lde_solve(fac, rs, rz)
	ref = as.vector(Hi %*% c(unlist(rs), rz))
	expect_equal(unlist(sol$s), ref[seq_len(Tt * qs)], tolerance = 1e-10)
	expect_equal(sol$z, ref[zi], tolerance = 1e-10)
	post = lame:::.lde_posterior(fac)
	ix = seq_len(qx)
	for (t in seq_len(Tt)) {
		expect_equal(post$var_s[, t], diag(Hi[blk(t), blk(t)]), tolerance = 1e-10)
		expect_equal(
			post$cov_xl[, t],
			diag(Hi[blk(t)[ix], zi[seq_len(nl)]]),
			tolerance = 1e-10
		)
		if (t < Tt) {
			expect_equal(
				post$cov_lag[, t],
				diag(Hi[blk(t), blk(t + 1L)]),
				tolerance = 1e-10
			)
		}
	}
	expect_equal(post$cov_zz, Hi[zi, zi], tolerance = 1e-10)
})

test_that("the laplace correction is the derivative of the laplace-marginal objective", {
	# logit link: the observed and the expected information coincide, so the
	# estimating equations are exactly the gradient of the marginal objective
	for (mode in c("directed", "bipartite", "symmetric")) {
		d = eb_sim(n = 10, nb = 8, Tt = 3, mode = mode, seed = 3, beta = 0.5)
		prep = lame:::.ae_prepare(
			d$Y,
			d$X,
			NULL,
			NULL,
			if (mode == "bipartite") "bipartite" else "unipartite",
			longitudinal = TRUE,
			node_time_policy = "quiet",
			symmetric = mode == "symmetric"
		)
		ctx = lame:::.lde_context(
			prep,
			mode == "symmetric",
			"logit",
			coef_prior = "flat",
			centre = FALSE
		)
		hyp = list(sab = diag(c(0.5, 0.5)), s2 = 0.2, rho = 0.6)
		st0 = lame:::.lde_init(ctx)
		st0$hyp = hyp
		inner = function(theta) {
			st = st0
			st$theta = theta
			for (k in 1:200) {
				eta = lame:::.lde_eta(ctx, st)
				wk = lame:::.lde_work(ctx$y, eta, ctx$link)
				sysm = lame:::.lde_system(ctx, st, wk)
				facx = lame:::.lde_factor(sysm$Ax, sysm$bx, sysm$Cx, sysm$Hzx, ctx$q, ctx$q)
				rs = lapply(seq_len(ctx$Tt), function(t) {
					xi = matrix(ctx$Xc[[t]] %*% theta, ctx$na, ctx$nb)
					lame:::.lde_stack_sums(wk$w[, , t] * (wk$z[, , t] - xi), ctx)
				})
				sol = lame:::.lde_solve_cond(facx, rs)
				dx = max(abs(do.call(cbind, sol$s) - st$x), abs(sol$lev - st$lev))
				st$x = do.call(cbind, sol$s)
				st$lev = sol$lev
				if (dx < 1e-12) break
			}
			st
		}
		marg = function(theta) {
			st = inner(theta)
			eta = lame:::.lde_eta(ctx, st)
			wk = lame:::.lde_work(ctx$y, eta, ctx$link)
			sysm = lame:::.lde_system(ctx, st, wk)
			facx = lame:::.lde_factor(sysm$Ax, sysm$bx, sysm$Cx, sysm$Hzx, ctx$q, ctx$q)
			dd = st$x - st$lev
			Qd = lame:::.lde_q_times(t(dd), hyp$rho)
			pen = sum(t(Qd) * dd) / hyp$s2 + sum(st$lev * (sysm$Lam %*% st$lev))
			ctx$nu * wk$dev + pen + facx$logdet_xx + facx$logdet_l
		}
		theta0 = c(-0.7, 0.45)
		st = inner(theta0)
		eta = lame:::.lde_eta(ctx, st)
		wk = lame:::.lde_work(ctx$y, eta, ctx$link)
		sysm = lame:::.lde_system(ctx, st, wk)
		facx = lame:::.lde_factor(sysm$Ax, sysm$bx, sysm$Cx, sysm$Hzx, ctx$q, ctx$q)
		post = lame:::.lde_posterior(facx, cond = TRUE)
		lap = lame:::.lde_laplace(ctx, facx, post$cond, wk, st)
		score = numeric(ctx$pall)
		for (t in seq_len(ctx$Tt)) {
			score = score + ctx$nu * as.vector(crossprod(ctx$Xc[[t]], as.vector(wk$score[, , t])))
		}
		analytic = -2 * (score - lap$corr)
		eps = 1e-5
		numeric_g = vapply(seq_along(theta0), function(k) {
			e = numeric(length(theta0))
			e[k] = eps
			(marg(theta0 + e) - marg(theta0 - e)) / (2 * eps)
		}, numeric(1))
		expect_equal(analytic, numeric_g, tolerance = 1e-5, info = mode)
		# the correction is not negligible: the uncorrected score is a different number
		expect_gt(max(abs(lap$corr)), 0.1)
	}
})

test_that("a dynamic binary ALS fit is an empirical-Bayes fit by default", {
	d = eb_sim(n = 14, Tt = 4, seed = 2)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	))
	expect_s3_class(fit, "lame_dynamic_als")
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(isTRUE(fit$bias_correct$default))
	expect_true(isTRUE(fit$bias_correct$dynamic))
	expect_true(fit$converged)
	expect_true(all(is.finite(fit$a_dynamic)))
	expect_true(all(is.finite(fit$eb$se)))
	expect_gt(fit$rho_ab, -1)
	expect_lt(fit$rho_ab, 1)
	expect_gt(fit$sigma_ab, 0)
	expect_true(all(eigen(fit$eb$sab, only.values = TRUE)$values > 0))
	expect_equal(unname(fit$VC["va"]), fit$eb$sab[1, 1])
	expect_equal(dim(fit$a_dynamic), c(14L, 4L))
	# the covariance of the coefficients is what vcov() and confint() report
	expect_equal(vcov(fit), fit$eb$vcov)
	expect_equal(unname(sqrt(diag(vcov(fit)))), unname(fit$eb$se))
	ci = confint(fit)
	expect_true(all(ci[, 1] < fit$coefficients & fit$coefficients < ci[, 2]))
	out = capture.output(print(fit))
	expect_true(any(grepl("empirical", out)))
	out = capture.output(print(summary(fit)))
	expect_true(any(grepl("rho_ab", out)))
	# the explicit request is the same fit, recorded as requested
	fe = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		bias_correct = "eb",
		verbose = FALSE
	))
	expect_identical(fe$bias_correct$requested, "eb")
	expect_false(isTRUE(fe$bias_correct$default))
	expect_equal(fe$coefficients, fit$coefficients)
	# the uncorrected penalised fit is available by name
	fn = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		bias_correct = "none",
		verbose = FALSE
	))
	expect_identical(fn$bias_correct$method, "none")
	expect_null(fn$eb)
})

test_that("undirected, bipartite, factor and moving-coefficient fits are empirical Bayes too", {
	skip_on_cran()
	spec = list(
		sym = list(mode = "symmetric", args = list(symmetric = TRUE)),
		bip = list(mode = "bipartite", args = list(mode = "bipartite")),
		fac = list(mode = "directed", args = list(R = 1)),
		mov = list(
			mode = "directed",
			args = list(
				dynamic_beta = "x1_dyad",
				dynamic_beta_kind = "rw1"
			)
		)
	)
	for (nm in names(spec)) {
		d = eb_sim(n = 10, nb = 7, Tt = 3, mode = spec[[nm]]$mode, seed = 8)
		fit = suppressWarnings(do.call(
			lame,
			c(
				list(
					d$Y,
					Xdyad = d$X,
					family = "binary",
					dynamic_ab = TRUE,
					method = "als",
					verbose = FALSE
				),
				spec[[nm]]$args
			)
		))
		expect_identical(fit$bias_correct$method, "eb", info = nm)
		expect_true(all(is.finite(fit$eb$se)), info = nm)
		expect_true(all(is.finite(fit$a_dynamic)), info = nm)
		expect_equal(unname(sqrt(diag(vcov(fit)))), unname(fit$eb$se), info = nm)
	}
	# a moving coefficient reports a path and its standard errors per period
	expect_equal(dim(fit$coef_path), c(3L, 3L))
	expect_gt(diff(range(fit$coef_path["x1_dyad", ])), 0)
	expect_equal(diff(range(fit$coef_path["x2_dyad", ])), 0)
	expect_equal(dim(fit$eb$se_path), c(3L, 3L))
	expect_true(all(fit$eb$se_path > 0))
})

# the fitter, called directly: the front door forwards these switches only where lame() supports them

test_that("a disabled margin is a point mass at zero, as in the sampler", {
	skip_on_cran()
	d = eb_sim(n = 10, Tt = 3, seed = 12)
	f0 = eb_fit(d)
	fr = eb_fit(d, rvar = FALSE)
	expect_true(fr$converged)
	expect_true(all(fr$a_dynamic == 0) && all(fr$a_level == 0))
	expect_equal(unname(fr$VC["va"]), 0)
	expect_equal(unname(fr$VC["cab"]), 0)
	expect_false(any(fr$b_dynamic == 0))
	expect_equal(unname(fr$eb$effects), c(FALSE, TRUE))
	fc = eb_fit(d, cvar = FALSE)
	expect_true(all(fc$b_dynamic == 0) && all(fc$a_dynamic != 0))
	# the fit is not the full fit with the effects thrown away
	expect_gt(max(abs(fr$coefficients - f0$coefficients)), 1e-3)
	# an undirected network has one effect per actor, switched by nvar (rvar and
	# cvar do not apply to it)
	ds = eb_sim(n = 10, Tt = 3, seed = 13, mode = "symmetric")
	fs = lame:::lame_dynamic_als(
		ds$Y,
		ds$X,
		family = "binary",
		symmetric = TRUE,
		dynamic_ab = TRUE,
		nvar = FALSE,
		verbose = FALSE
	)
	expect_true(all(fs$a_dynamic == 0))
	fs2 = lame:::lame_dynamic_als(
		ds$Y,
		ds$X,
		family = "binary",
		symmetric = TRUE,
		dynamic_ab = TRUE,
		rvar = FALSE,
		verbose = FALSE
	)
	expect_true(all(fs2$a_dynamic != 0))
	# bipartite margins
	db = eb_sim(n = 10, nb = 7, Tt = 3, mode = "bipartite", seed = 14)
	fb = suppressWarnings(lame:::lame_dynamic_als(
		db$Y,
		db$X,
		family = "binary",
		mode = "bipartite",
		dynamic_ab = TRUE,
		cvar = FALSE,
		verbose = FALSE
	))
	expect_true(all(fb$b_dynamic == 0) && all(fb$a_dynamic != 0))
	# every margin off: a probit regression with no path law to report
	fn = eb_fit(d, rvar = FALSE, cvar = FALSE)
	expect_true(fn$converged)
	expect_null(fn$rho_ab)
	expect_false(fn$dynamic_ab)
	expect_true(all(is.finite(fn$eb$se)))
})

test_that("switches a route cannot honour are named, never dropped", {
	d = eb_sim(n = 10, Tt = 3, seed = 15)
	catch = function(expr) {
		msgs = character()
		val = withCallingHandlers(expr, warning = function(w) {
			msgs <<- c(msgs, conditionMessage(w))
			invokeRestart("muffleWarning")
		})
		list(val = val, msgs = msgs)
	}
	# the empirical-Bayes route models the dyadic correlation of a directed network (lame()'s default): no message,
	# a record of the pair likelihood; dcor = FALSE is the independent-cell fit
	r = catch(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		dcor = FALSE,
		verbose = FALSE
	))
	expect_false(any(grepl("ignored MCMC-only", r$msgs)))
	expect_s3_class(r$val, "lame_dynamic_als")
	expect_null(r$val$eb$pair)
	r = catch(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		dcor = TRUE,
		verbose = FALSE
	))
	expect_false(any(grepl("ignored MCMC-only", r$msgs)))
	expect_false(is.null(r$val$eb$pair))
	# the penalised route honours none of them
	r = catch(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		bias_correct = "none",
		rvar = FALSE,
		dcor = TRUE,
		dynamic_beta_pool = "both",
		verbose = FALSE
	))
	expect_true(any(grepl("ignored MCMC-only arguments.*rvar", r$msgs)))
	expect_true(any(grepl("ignored MCMC-only arguments.*dcor", r$msgs)))
	expect_true(any(grepl("dynamic_beta_pool", r$msgs)))
	# honoured switches are silent
	r = catch(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		rvar = FALSE,
		cvar = TRUE,
		dynamic_beta_pool = "both",
		verbose = FALSE
	))
	expect_false(any(grepl("ignored MCMC-only", r$msgs)))
	# a malformed switch is an error
	expect_error(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			dynamic_ab = TRUE,
			rvar = "no",
			verbose = FALSE
		),
		"rvar"
	)
})

test_that("unequal row and column ranks of a bipartite network fit the interaction of rank min(R_row, R_col)", {
	# a small panel's product is withheld below its information threshold, so the rule is off here
	withr::local_options(list(lame.eb_withhold_product = FALSE))
	skip_on_cran()
	d = eb_sim(n = 12, nb = 8, Tt = 3, mode = "bipartite", seed = 16)
	args = list(d$Y, d$X, family = "binary", mode = "bipartite", dynamic_ab = TRUE, verbose = FALSE)
	f1 = suppressWarnings(do.call(lame:::lame_dynamic_als, c(args, list(R = 1))))
	f21 = suppressWarnings(do.call(lame:::lame_dynamic_als, c(args, list(R_row = 2, R_col = 1))))
	f12 = suppressWarnings(do.call(lame:::lame_dynamic_als, c(args, list(R_row = 1, R_col = 2))))
	# the interaction is the one of rank 1: the same fit, factors padded with zeros
	expect_equal(f21$coefficients, f1$coefficients)
	expect_equal(f12$coefficients, f1$coefficients)
	expect_equal(dim(f21$U), c(12L, 2L)); expect_equal(dim(f21$V), c(8L, 1L))
	expect_equal(dim(f12$U), c(12L, 1L)); expect_equal(dim(f12$V), c(8L, 2L))
	expect_true(all(f21$U[, 2] == 0))
	expect_equal(unname(f21$G), matrix(c(1, 0), 2, 1))
	expect_equal(c(f21$R_row, f21$R_col), c(2L, 1L))
	expect_equal(f21$UVPM, f1$UVPM)
	# as in the sampler, one rank cannot be positive and the other zero
	expect_error(
		suppressWarnings(do.call(lame:::lame_dynamic_als, c(args, list(R_row = 1, R_col = 0)))),
		"Only one of"
	)
})

test_that("a second-order random-walk coefficient path is an empirical-Bayes path", {
	d = eb_sim(n = 10, Tt = 5, seed = 17)
	f = eb_fit(d, dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw2")
	expect_true(f$converged)
	expect_identical(f$bias_correct$method, "eb")
	expect_equal(dim(f$coef_path), c(3L, 5L))
	expect_gt(diff(range(f$coef_path["x1_dyad", ])), 0)
	expect_true(all(is.finite(f$eb$se_path)))
	expect_equal(unname(f$rho_beta), 1)
	expect_gt(unname(f$sigma_beta), 0)
	# the second-difference prior needs three periods, as in the penalised route
	expect_error(
		eb_fit(
			list(Y = d$Y[1:2], X = d$X[1:2]),
			dynamic_beta = "x1_dyad",
			dynamic_beta_kind = "rw2"
		),
		"at least 3 time periods"
	)
	# a rougher path than a random walk would give is shrunk further by the
	# second difference: the path of a straight trend is recovered closely
	set.seed(3)
	trend = 0.2 * (seq_len(5) - 3)
	Y = d$Y; X = d$X
	for (t in seq_len(5)) {
		eta = -0.8 + trend[t] * X[[t]][, , 1] - 0.4 * X[[t]][, , 2]
		y = 1 * (eta + matrix(rnorm(100), 10) > 0); diag(y) = NA
		dimnames(y) = dimnames(d$Y[[t]]); Y[[t]] = y
	}
	g = eb_fit(list(Y = Y, X = X), dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw2")
	expect_true(g$converged)
	expect_gt(stats::cor(g$coef_path["x1_dyad", ], trend), 0.5)
})

test_that("dynamic_beta_pool shares the innovation scales across the blocks", {
	skip_on_cran()
	d = eb_sim(n = 8, Tt = 3, seed = 18)
	args = list(dynamic_beta = c("intercept", "x1_dyad"), dynamic_beta_kind = "rw1")
	f0 = do.call(eb_fit, c(list(d), args))
	fs = do.call(eb_fit, c(list(d), args, list(dynamic_beta_pool = "sigma")))
	expect_null(f0$eb$beta_pool$sigma_share)
	expect_true(is.finite(fs$eb$beta_pool$sigma_share) && fs$eb$beta_pool$sigma_share > 0)
	expect_identical(fs$dynamic_beta_pool, "sigma")
	# the hierarchy moves the scales
	expect_gt(max(abs(fs$sigma_beta - f0$sigma_beta)), 1e-4)
	# one moving block has nothing to pool against
	f1 = eb_fit(d, dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw1")
	f1p = eb_fit(d, dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw1", dynamic_beta_pool = "both")
	expect_equal(f1p$coefficients, f1$coefficients)
})

test_that("an exhausted iteration budget is refused, with the unfinished state in the condition", {
	d = eb_sim(n = 10, Tt = 3, seed = 6)
	err = expect_error(
		lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			dynamic_ab = TRUE,
			method = "als",
			als_max_iter = 2,
			verbose = FALSE
		),
		class = "lame_als_convergence"
	)
	# the message states the state at the stop
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(
		msg,
		"did not reach the empirical-Bayes convergence criterion: it stopped after 2 passes (iteration limit)",
		fixed = TRUE
	)
	expect_match(msg, "largest scaled residual of the empirical-Bayes equations is")
	expect_false(grepl("caution|read the", msg))
	# the last iterate is a field of the condition
	fit = err$fit
	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$converged)
	expect_false(fit$diagnostics$stationarity$pass)
	expect_gt(fit$diagnostics$stationarity$max_score, fit$diagnostics$stationarity$threshold)
	expect_equal(fit$iterations, 2L)
	expect_lte(length(fit$dev_history), 2L)
	# the same panel at the default budget is a fit
	ok = expect_no_warning(
		lame(d$Y, Xdyad = d$X, family = "binary", dynamic_ab = TRUE, method = "als", verbose = FALSE)
	)
	expect_true(ok$converged)
})

test_that("an actor without an event has finite effects under the random-effect prior", {
	d = eb_sim(n = 14, Tt = 4, seed = 5)
	for (t in seq_along(d$Y)) {
		d$Y[[t]][1, ] = ifelse(is.na(d$Y[[t]][1, ]), NA, 0)
		d$Y[[t]][, 2] = ifelse(is.na(d$Y[[t]][, 2]), NA, 0)
	}
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	))
	expect_true(fit$converged)
	expect_true(all(is.finite(fit$a_dynamic)))
	expect_true(all(is.finite(fit$b_dynamic)))
	# a sender without a tie is pulled below the average sender, but not to minus infinity
	expect_lt(mean(fit$a_dynamic[1, ]), 0)
	expect_gt(min(fit$a_dynamic[1, ]), -6)
	expect_lt(mean(fit$b_dynamic[2, ]), 0)
})

# moving latent positions: the building blocks

test_that("the position update reproduces the dense posterior of every actor's path", {
	set.seed(11)
	n = 3L; m = 5L; R = 2L; Tt = 4L
	W = lapply(seq_len(Tt), function(t) matrix(runif(n * m, 0.05, 0.4), n, m))
	Rr = lapply(seq_len(Tt), function(t) matrix(rnorm(n * m), n, m))
	Vh = lapply(seq_len(Tt), function(t) matrix(rnorm(m * R, 0, 0.6), m, R))
	Om = lapply(seq_len(Tt), function(t) {
		a = array(0, c(R, R, m))
		for (j in seq_len(m)) a[, , j] = crossprod(matrix(rnorm(4, 0, 0.2), 2)) + diag(0.05, R)
		a
	})
	s2 = 0.12; rho = 0.8
	out = lame:::.lde_pos_side(W, Rr, Vh, Om, s2, rho, n, R, Tt)
	P0 = kronecker(lame:::.lde_pos_q(Tt, rho) / s2, diag(R))
	kl = 0
	for (i in seq_len(n)) {
		P = P0
		h = numeric(Tt * R)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			for (j in seq_len(m)) {
				P[ix, ix] = P[ix, ix] + W[[t]][i, j] * (tcrossprod(Vh[[t]][j, ]) + Om[[t]][, , j])
				h[ix] = h[ix] + W[[t]][i, j] * Rr[[t]][i, j] * Vh[[t]][j, ]
			}
		}
		S = solve(P)
		mu = as.vector(S %*% h)
		expect_equal(as.vector(out$mean[i, , ]), mu, tolerance = 1e-10)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			expect_equal(out$sig[, , i, t], S[ix, ix], tolerance = 1e-10)
			if (t < Tt) expect_equal(out$cl[, , i, t], S[ix, ix + R], tolerance = 1e-10)
		}
		# the KL divergence of N(mu, S) from the prior N(0, P0^-1), by its definition
		Sp = solve(P0)
		kl = kl + 0.5 * (sum(diag(solve(Sp, S))) + sum(mu * (P0 %*% mu)) - length(mu) +
				as.numeric(determinant(Sp)$modulus) - as.numeric(determinant(S)$modulus))
	}
	pos = list(qu = out, qv = out)
	expect_equal(lame:::.lde_pos_kl2(pos, R, Tt, lame:::.lde_pos_hp(s2 / (1 - rho^2), rho)), 4 * kl, tolerance = 1e-10)
	# the stationary AR(1) precision is the inverse of its covariance
	cv = outer(1:5, 1:5, function(a, b) 0.7^abs(a - b)) / (1 - 0.7^2)
	expect_equal(solve(lame:::.lde_pos_q(5L, 0.7)), cv, tolerance = 1e-10)
	# the EM update of the persistence at a held scale maximises its objective
	St = lame:::.lde_pos_stats(pos, R, Tt)
	pr = lame:::.lde_pos_prior(list())
	S0 = 0.5
	cnt = St$nser * St$Tt
	f = function(rho) {
		s2 = S0 * (1 - rho^2)
		-(cnt / 2) * log(S0) - ((cnt - St$nser) / 2) * log(1 - rho^2) -
			(St$S1 / (2 * S0) + (St$Sn - 2 * rho * St$Sc + rho^2 * St$Sl) / (2 * S0 * (1 - rho^2))) -
			(pr$shape + 1) * log(s2) - pr$scale / s2 - (rho - pr$rho_mean)^2 / (2 * pr$rho_sd^2)
	}
	g = seq(-0.99, 0.99, by = 1e-4)
	rr = lame:::.lde_pos_rho(St, S0, pr, rho_ref = 0.5)
	expect_equal(rr$rho, g[which.max(vapply(g, f, numeric(1)))], tolerance = 1e-3)
	expect_equal(rr$gain, 2 * (f(rr$rho) - f(0.5)), tolerance = 1e-8)
})

test_that("a sweep of both sides with unequal actor sets is the dense posterior of each side given the other", {
	set.seed(31)
	na = 3L; nb = 5L; R = 2L; Tt = 3L
	obs = array(runif(na * nb * Tt) > 0.2, c(na, nb, Tt))
	ctx = list(Tt = Tt, na = na, nb = nb, R = R, sym = FALSE, obs = obs)
	wk = list(w = array(runif(na * nb * Tt, 0.05, 0.4), c(na, nb, Tt)), z = array(rnorm(na * nb * Tt), c(na, nb, Tt)))
	eta_add = array(rnorm(na * nb * Tt, 0, 0.3), c(na, nb, Tt))
	mk = function(n) {
		sig = array(0, c(R, R, n, Tt))
		for (i in seq_len(n)) for (t in seq_len(Tt)) sig[, , i, t] = crossprod(matrix(rnorm(4, 0, 0.2), 2)) + diag(0.05, R)
		list(mean = array(rnorm(n * R * Tt, 0, 0.5), c(n, R, Tt)), sig = sig)
	}
	pos = list(qu = mk(na), qv = mk(nb))
	hp = lame:::.lde_pos_hp(0.5, 0.8)
	out = lame:::.lde_pos_sweep(ctx, pos, hp, wk, eta_add)
	P0 = kronecker(lame:::.lde_pos_q(Tt, hp$rho) / hp$s2, diag(R))
	W = ifelse(obs, wk$w, 0)
	Rr = ifelse(obs, wk$z - eta_add, 0)
	dense = function(i, own_is_row, other) {
		P = P0
		h = numeric(Tt * R)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			m = dim(other$mean)[1L]
			for (j in seq_len(m)) {
				w = if (own_is_row) W[i, j, t] else W[j, i, t]
				r = if (own_is_row) Rr[i, j, t] else Rr[j, i, t]
				v = other$mean[j, , t]
				P[ix, ix] = P[ix, ix] + w * (tcrossprod(v) + other$sig[, , j, t])
				h[ix] = h[ix] + w * r * v
			}
		}
		S = solve(P)
		list(mean = as.vector(S %*% h), S = S)
	}
	# the row side is given the incoming column side, the column side the updated row side
	for (i in seq_len(na)) {
		d = dense(i, TRUE, pos$qv)
		expect_equal(as.vector(out$qu$mean[i, , ]), d$mean, tolerance = 1e-10)
		for (t in seq_len(Tt)) expect_equal(
			out$qu$sig[, , i, t],
			d$S[(t - 1L) * R + 1:R, (t - 1L) * R + 1:R],
			tolerance = 1e-10
		)
	}
	for (j in seq_len(nb)) {
		d = dense(j, FALSE, out$qu)
		expect_equal(as.vector(out$qv$mean[j, , ]), d$mean, tolerance = 1e-10)
		for (t in seq_len(Tt)) expect_equal(
			out$qv$sig[, , j, t],
			d$S[(t - 1L) * R + 1:R, (t - 1L) * R + 1:R],
			tolerance = 1e-10
		)
	}
})

test_that("the reported positions of an undirected network keep the product and the sign structure, and report is idempotent", {
	set.seed(17)
	n = 7L; R = 3L; Tt = 4L
	sg = c(1, -1, 1)
	pos = list(qu = list(mean = array(rnorm(n * R * Tt), c(n, R, Tt)), sig = array(0, c(R, R, n, Tt))))
	for (i in seq_len(n)) for (t in seq_len(Tt)) pos$qu$sig[, , i, t] = crossprod(matrix(rnorm(9, 0, 0.3), 3)) + diag(0.1, R)
	ctx = list(Tt = Tt, R = R, sgn = sg, sym = TRUE)
	rp = lame:::.lde_pos_report_sym(ctx, pos)
	prod_of = function(U, V, t) tcrossprod(U[, , t], V[, , t])
	for (t in seq_len(Tt)) {
		before = tcrossprod(sweep(pos$qu$mean[, , t], 2L, sg, "*"), pos$qu$mean[, , t])
		expect_equal(prod_of(rp$U, rp$V, t), before, tolerance = 1e-10)
	}
	# the rotation is orthogonal and mixes only columns of one sign
	expect_equal(crossprod(rp$rotation), diag(R), tolerance = 1e-10)
	expect_equal(rp$rotation[c(1L, 3L), 2L], c(0, 0))
	expect_equal(rp$rotation[2L, c(1L, 3L)], c(0, 0))
	# reported again, the reported positions stay where they are
	pos2 = list(qu = list(mean = rp$U, sig = array(0, c(R, R, n, Tt))))
	for (i in seq_len(n)) for (t in seq_len(Tt)) pos2$qu$sig[, , i, t] = crossprod(
		rp$rotation,
		pos$qu$sig[, , i, t] %*% rp$rotation
	)
	rp2 = lame:::.lde_pos_report_sym(ctx, pos2)
	expect_equal(rp2$U, rp$U, tolerance = 1e-8)
})

test_that("the position noise enters the cell problem as the expected score, weight and log likelihood", {
	y = matrix(c(1, 0, 1, 0, NA, 1), 3)
	eta = matrix(c(-1.3, 0.4, -2.2, 0.9, 0.2, -45), 3)
	# no noise: the plain problem, the capped increment of the extreme cell included
	w0 = lame:::.lde_work(y, eta, "probit")
	wz = lame:::.lde_work_gh(y, eta, matrix(0, 3, 2))
	for (nm in c("z", "score", "w", "hprime", "hpp", "dev")) expect_equal(wz[[nm]], w0[[nm]], tolerance = 1e-10)
	expect_equal(wz$dev_marg, w0$dev, tolerance = 1e-10)
	# an unobserved cell carries nothing
	expect_equal(c(wz$w[2, 2], wz$score[2, 2], wz$hprime[2, 2], wz$hpp[2, 2]), numeric(4))
	# with noise: expectations over N(0, v) of the cell's score, Fisher weight, its derivatives and log likelihood
	# (the 9-node rule is accurate to 1e-4 at these variances)
	y2 = matrix(c(1, 0, 1, 0), 2)
	eta2 = matrix(c(-1.3, 0.4, -2.2, 0.9), 2)
	v2 = matrix(c(0.3, 0.1, 0.5, 0.2), 2)
	wg = lame:::.lde_work_gh(y2, eta2, v2)
	ex = function(k, f) {
		integrate(function(d) {
			f(lame:::.lde_probit_raw(rep(y2[k], length(d)), eta2[k] + d)) * dnorm(d, 0, sqrt(v2[k]))
		}, -Inf, Inf, rel.tol = 1e-12)$value
	}
	for (k in 1:4) {
		expect_equal(wg$score[k], ex(k, function(r) r$score), tolerance = 1e-4)
		expect_equal(wg$w[k], ex(k, function(r) pmax(r$wraw, 1e-8)), tolerance = 1e-4)
		expect_equal(wg$hprime[k], ex(k, function(r) r$hprime), tolerance = 1e-4)
		expect_equal(wg$hpp[k], ex(k, function(r) r$hpp), tolerance = 1e-3)
	}
	expect_equal(wg$dev, -2 * sum(vapply(1:4, function(k) ex(k, function(r) r$ell), numeric(1))), tolerance = 1e-6)
	# the working response is the incoming predictor plus the expected score over the expected weight
	expect_equal(wg$z, eta2 + wg$score / wg$w, tolerance = 1e-8)
	# the blocks in which the reference processes a large panel change no digit, a partial last block included
	wr = lame:::.lde_work_gh_ref(y2, eta2, v2)
	wb = lame:::.lde_work_gh_ref(y2, eta2, v2, block = 3L)
	for (nm in c(
		"z", "w", "score", "hprime", "hpp", "dev", "dev_marg"
	)) expect_equal(wb[[nm]], wr[[nm]], tolerance = 1e-14)
	# a cell whose outcome is certain keeps a positive weight and no shift (the floor), and everything stays finite
	we = lame:::.lde_work_gh(matrix(1, 1, 1), matrix(40, 1, 1), matrix(0.1, 1, 1))
	expect_true(all(is.finite(c(we$z, we$w, we$score, we$hprime, we$hpp, we$dev))))
	expect_equal(c(we$w), 1e-8)
	expect_equal(c(we$z), 40)
	# the marginal deviance integrates the noise out of the likelihood, not of its logarithm: exact for the probit,
	# and never above the deviance of the bound (Jensen)
	lm = function(k) {
		log(integrate(function(d) {
			p = pnorm(eta2[k] + d)
			(if (y2[k] == 1) p else 1 - p) * dnorm(d, 0, sqrt(v2[k]))
		}, -Inf, Inf, rel.tol = 1e-12)$value)
	}
	expect_equal(wg$dev_marg, -2 * sum(vapply(1:4, lm, numeric(1))), tolerance = 1e-8)
	expect_lt(wg$dev_marg, wg$dev)
})

test_that("the compiled cell problem equals the R reference, extreme predictors, missing cells and zero noise included", {
	set.seed(31)
	n = 400L
	y = array(rbinom(n, 1, 0.3), c(20, 10, 2))
	y[sample(n, 15)] = NA
	cases = list(
		list(eta = rnorm(n, -1.5, 1.5), v = rexp(n, 3)),
		list(eta = runif(n, -12, 12), v = runif(n, 0, 4)),
		list(eta = runif(n, -45, 45), v = runif(n, 0, 20)),
		list(eta = rnorm(n, -1.5, 1.5), v = rep(0, n)),
		list(eta = rnorm(n, -1.5, 1.5), v = rep(1e-12, n))
	)
	for (cs in cases) {
		eta = array(cs$eta, dim(y))
		v = array(cs$v, dim(y))
		a = lame:::.lde_work_gh_ref(y, eta, v)
		b = lame:::.lde_work_gh(y, eta, v)
		for (nm in c("z", "w", "score", "hprime", "hpp")) expect_equal(b[[nm]], a[[nm]], tolerance = 1e-10)
		expect_equal(b$dev, a$dev, tolerance = 1e-10)
		expect_equal(b$dev_marg, a$dev_marg, tolerance = 1e-10)
		expect_identical(dim(b$z), dim(y))
		# the missing cells carry nothing
		expect_true(all(b$w[is.na(y)] == 0 & b$score[is.na(y)] == 0))
	}
	# a scalar noise variance is a variance for every cell
	expect_equal(
		lame:::.lde_work_gh(y, array(-1, dim(y)), 0.3)$w,
		lame:::.lde_work_gh(y, array(-1, dim(y)), array(0.3, dim(y)))$w
	)
	# what is not wanted is not formed, and what is wanted does not change
	eta = array(rnorm(n, -1.5, 1.5), dim(y))
	v = array(rexp(n, 3), dim(y))
	full = lame:::.lde_work_gh(y, eta, v)
	mid = lame:::.lde_work_gh(y, eta, v, what = 1L)
	low = lame:::.lde_work_gh(y, eta, v, what = 0L)
	for (nm in c("z", "w")) {
		expect_identical(low[[nm]], full[[nm]])
		expect_identical(mid[[nm]], full[[nm]])
	}
	for (nm in c("score", "hprime", "hpp", "dev_marg")) expect_identical(mid[[nm]], full[[nm]])
	expect_true(is.na(mid$dev) && is.na(low$dev) && is.na(low$dev_marg))
	expect_true(all(low$score == 0))
})

test_that("the logit link: the compiled cell problem equals the R reference, is the plain logit problem without noise, and integrates the probability", {
	set.seed(32)
	n = 300L
	y = array(rbinom(n, 1, 0.4), c(15, 10, 2))
	y[sample(n, 12)] = NA
	for (cs in list(
		list(eta = rnorm(n, 0, 1.5), v = rexp(n, 2)),
		list(eta = runif(n, -35, 35), v = runif(n, 0, 20)),
		list(eta = rnorm(n, 0, 1.5), v = rep(0, n))
	)) {
		eta = array(cs$eta, dim(y))
		v = array(cs$v, dim(y))
		a = lame:::.lde_work_gh_ref(y, eta, v, link = "logit")
		b = lame:::.lde_work_gh(y, eta, v, link = "logit")
		for (nm in c("z", "w", "score", "hprime", "hpp")) expect_equal(b[[nm]], a[[nm]], tolerance = 1e-10)
		expect_equal(b$dev, a$dev, tolerance = 1e-10)
		expect_equal(b$dev_marg, a$dev_marg, tolerance = 1e-10)
	}
	# without noise it is the plain problem of the logit link
	ob = is.finite(y)
	eta = array(rnorm(n), dim(y))
	g0 = lame:::.lde_work_gh(y, eta, 0, link = "logit")
	pl = lame:::.lde_work(y[ob], eta[ob], "logit")
	expect_equal(g0$z[ob], pl$z, tolerance = 1e-10)
	expect_equal(g0$w[ob], pl$w, tolerance = 1e-10)
	expect_equal(g0$dev, pl$dev, tolerance = 1e-10)
	# the derivatives of the expected weight are those of the weight (cells whose increment is not capped)
	h = 1e-5
	v = array(runif(n, 0.2, 3), dim(y))
	w0 = lame:::.lde_work_gh(y, eta, v, what = 1L, link = "logit")
	wp = lame:::.lde_work_gh(y, eta + h, v, what = 1L, link = "logit")
	wm = lame:::.lde_work_gh(y, eta - h, v, what = 1L, link = "logit")
	ok = ob & abs(w0$z - eta) < 29
	expect_equal(w0$hprime[ok], ((wp$w - wm$w) / (2 * h))[ok], tolerance = 1e-6)
	expect_equal(w0$hpp[ok], ((wp$hprime - wm$hprime) / (2 * h))[ok], tolerance = 1e-6)
	# the marginal likelihood of a cell is the integral of its probability over the predictor
	m = 0.7; vv = 5
	f = function(x) stats::plogis(m + sqrt(vv) * x) * stats::dnorm(x)
	one = lame:::.lde_work_gh(
		array(1, c(1, 1, 1)),
		array(m, c(1, 1, 1)),
		array(vv, c(1, 1, 1)),
		what = 1L,
		link = "logit"
	)
	expect_equal(exp(-0.5 * one$dev_marg), stats::integrate(f, -Inf, Inf)$value, tolerance = 1e-7)
})

test_that("a cell whose predictor has a large variance is integrated by a larger rule, and the rule is the 9-node one where the variance is small", {
	# the weight, the score and the derivative of the weight of one cell by numerical integration of the
	# probit terms over the predictor's law
	ex = function(y, eta, v, f) {
		integrate(
			function(d) f(lame:::.lde_probit_raw(rep(y, length(d)), eta + d)) * dnorm(d, 0, sqrt(v)),
			-Inf,
			Inf,
			rel.tol = 1e-12,
			subdivisions = 2000L
		)$value
	}
	for (v in c(0.5, 10)) {
		truth_w = ex(1, 1.5, v, function(r) pmax(r$wraw, 1e-8))
		truth_s = ex(1, 1.5, v, function(r) r$score)
		ld = lame:::.lde_work_gh(matrix(1), matrix(1.5), matrix(v))
		r9 = lame:::.lde_work_gh(matrix(1), matrix(1.5), matrix(v), per_sd = 0)
		if (v < 1) {
			# the base rule below a variance of 1.27: the ladder is the 9-node rule
			expect_identical(ld$w, r9$w)
			expect_lt(abs(ld$w / truth_w - 1), 1e-3)
		} else {
			# at a variance of 10 the 9-node rule is wrong by a fifth in the weight, the ladder by under 0.5 percent
			expect_gt(abs(r9$w / truth_w - 1), 0.05)
			expect_lt(abs(ld$w / truth_w - 1), 5e-3)
			expect_lt(abs(ld$score / truth_s - 1), 5e-3)
		}
	}
	# the rule of a cell depends on its own variance only
	y = matrix(1, 1, 2)
	both = lame:::.lde_work_gh(y, matrix(1.5, 1, 2), matrix(c(0.5, 10), 1, 2))
	expect_identical(both$w[1], c(lame:::.lde_work_gh(matrix(1), matrix(1.5), matrix(0.5))$w))
	expect_identical(both$w[2], c(lame:::.lde_work_gh(matrix(1), matrix(1.5), matrix(10))$w))
})

test_that("the compiled position sweep equals the R reference for directed and bipartite panels", {
	set.seed(4)
	for (cfg in list(
		c(na = 12, nb = 12, R = 2, Tt = 4),
		c(na = 9, nb = 9, R = 1, Tt = 2),
		c(na = 10, nb = 7, R = 3, Tt = 5),
		c(na = 6, nb = 11, R = 2, Tt = 3)
	)) {
		na = cfg[["na"]]; nb = cfg[["nb"]]; R = cfg[["R"]]; Tt = cfg[["Tt"]]
		dims = c(na, nb, Tt)
		side = function(m) {
			sg = array(0, c(R, R, m, Tt))
			cl = array(0, c(R, R, m, max(Tt - 1L, 1L)))
			for (i in seq_len(m)) for (t in seq_len(Tt)) sg[, , i, t] = crossprod(matrix(rnorm(R * R), R)) * 0.1 + diag(0.05, R)
			for (i in seq_len(m)) for (t in seq_len(Tt - 1L)) cl[, , i, t] = matrix(rnorm(R * R), R) * 0.02
			list(mean = array(rnorm(m * R * Tt), c(m, R, Tt)), sig = sg, cl = cl)
		}
		ctx = list(sym = FALSE, na = na, nb = nb, R = R, Tt = Tt, obs = array(TRUE, dims))
		ctx$obs[sample(prod(dims), 10)] = FALSE
		wk = list(w = array(runif(prod(dims), 0.05, 0.4), dims), z = array(rnorm(prod(dims)), dims))
		wk$w[sample(prod(dims), 5)] = NA
		eta = array(rnorm(prod(dims), -1, 0.5), dims)
		pos = list(qu = side(na), qv = side(nb))
		hp = lame:::.lde_pos_hp(0.7, 0.85)
		a = lame:::.lde_pos_sweep_ref(ctx, pos, hp, wk, eta)
		b = lame:::.lde_pos_sweep(ctx, pos, hp, wk, eta)
		for (sd in c("qu", "qv")) {
			for (f in c(
				"mean", "sig", "cl", "P0", "ldp", "trp", "ldp0"
			)) expect_equal(b[[sd]][[f]], a[[sd]][[f]], tolerance = 1e-10)
			expect_identical(dim(b[[sd]]$sig), dim(a[[sd]]$sig))
			expect_identical(dim(b[[sd]]$cl), dim(a[[sd]]$cl))
		}
	}
})

test_that("the product variance by one matrix product per period equals the sum of outer products", {
	set.seed(6)
	R = 2L; Tt = 3L
	mk = function(m) {
		sg = array(0, c(R, R, m, Tt))
		for (i in seq_len(m)) for (t in seq_len(Tt)) sg[, , i, t] = crossprod(matrix(rnorm(4), 2)) * 0.2 + diag(0.1, 2)
		list(mean = array(rnorm(m * R * Tt), c(m, R, Tt)), sig = sg)
	}
	pos = list(qu = mk(7L), qv = mk(5L))
	expect_equal(lame:::.lde_pos_prodvar(pos, R, Tt), lame:::.lde_pos_prodvar_ref(pos, R, Tt), tolerance = 1e-12)
})

test_that("the structured border complement and the triangular quadratic forms equal the dense ones", {
	set.seed(12)
	Tt = 4L; qx = 6L; m = 3L; r = qx + m
	mk = function(k) crossprod(matrix(rnorm(k * k), k)) * 0.3 + diag(k) * 4
	A = lapply(seq_len(Tt), function(t) mk(qx))
	b = rep(0.7, qx)
	cd = lapply(seq_len(Tt), function(t) rnorm(qx, -0.3, 0.1))
	C = lapply(seq_len(Tt), function(t) cbind(diag(cd[[t]]), matrix(rnorm(qx * m), qx, m) * 0.2))
	Hzz = crossprod(matrix(rnorm(r * r), r)) * 0.5 + diag(r) * 12
	fa = lame:::.lde_factor(A, b, C, Hzz, qx, qx)
	fb = lame:::.lde_factor(A, b, C, Hzz, qx, qx, cd)
	expect_equal(fb$Sz, fa$Sz, tolerance = 1e-12)
	expect_equal(fb$Szi, fa$Szi, tolerance = 1e-10)
	Hi = solve(lame:::.lde_dense_matrix(A, b, C, Hzz))
	blk = function(t) (t - 1L) * qx + seq_len(qx)
	zi = Tt * qx + seq_len(r)
	il = seq_len(qx)
	# the conditional covariance of the paths given the border outside the levels, by its definition
	pf = lame:::.lde_posterior(fb, cond = TRUE)
	pn = lame:::.lde_posterior(fb, cond = TRUE, moments = FALSE)
	for (t in seq_len(Tt)) {
		Hl = lame:::.lde_dense_matrix(A, b, lapply(C, function(Ct) Ct[, il, drop = FALSE]), Hzz[il, il])
		ref = solve(Hl)[blk(t), blk(t)]
		expect_equal(pf$cond[[t]], ref, tolerance = 1e-10)
		expect_equal(pn$cond[[t]], ref, tolerance = 1e-10)
		expect_equal(pf$var_s[, t], diag(Hi[blk(t), blk(t)]), tolerance = 1e-10)
	}
	# the moments are skipped when they are not wanted
	expect_null(pn$var_s)
	expect_null(pn$cov_zz)
	# the conditional solve with the kept factor of the level block is the solve with its inverse
	rs = lapply(seq_len(Tt), function(t) rnorm(qx))
	sc = lame:::.lde_solve_cond(fb, rs)
	fc = fb
	fc$chl = NULL
	sd = lame:::.lde_solve_cond(fc, rs)
	expect_equal(sc$lev, sd$lev, tolerance = 1e-12)
	expect_equal(unlist(sc$s), unlist(sd$s), tolerance = 1e-12)
})

test_that("a gate of plain passes between extrapolations changes the speed of the iteration, not its fixed point", {
	skip_on_cran()
	d = eb_sim(n = 12, Tt = 4, seed = 8, mu = -0.8)
	prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ctx = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE)
	fit = lame:::.lde_fit(ctx, lame:::.lde_init(ctx), max_iter = 8L, tol = 5e-5)
	st = fit$st
	st$hyp$uv = lame:::.lde_pos_hp(0.5, st$hyp$uv$rho)
	r0 = lame:::.lde_fit_core_dyn(fit$ctx, st, max_iter = 300L, tol = 1e-4, freeze_ab = TRUE, gate = 0L)
	r3 = lame:::.lde_fit_core_dyn(fit$ctx, st, max_iter = 300L, tol = 1e-4, freeze_ab = TRUE, gate = 3L)
	expect_true(r0$converged && r3$converged)
	# no gate: extrapolated passes follow each other; a gate of 3: at least 3 plain passes between two
	j0 = which(r0$trace$jump)
	j3 = which(r3$trace$jump)
	expect_gte(length(j3), 2L)
	expect_equal(min(diff(j0)), 1L)
	expect_gte(min(diff(j3)), 4L)
	expect_equal(r3$ev$objective, r0$ev$objective, tolerance = 1e-6)
	expect_equal(r3$st$theta, r0$st$theta, tolerance = 1e-4)
})

test_that("the reported positions are the fitted ones in the gauge of the average product", {
	set.seed(5)
	n = 6L; R = 2L; Tt = 4L
	mk = function(n) {
		m = array(rnorm(n * R * Tt), c(n, R, Tt))
		sg = array(0, c(R, R, n, Tt))
		for (i in seq_len(n)) for (t in seq_len(Tt)) sg[, , i, t] = crossprod(matrix(rnorm(4, 0, 0.3), 2)) + diag(0.1, 2)
		list(mean = m, sig = sg)
	}
	pos = list(qu = mk(6L), qv = mk(5L))
	ctx = list(Tt = Tt, R = R)
	rp = lame:::.lde_pos_report(ctx, pos)
	expect_equal(crossprod(rp$rotation), diag(R), tolerance = 1e-10)
	for (t in seq_len(Tt)) {
		expect_equal(
			tcrossprod(rp$U[, , t], rp$V[, , t]),
			tcrossprod(pos$qu$mean[, , t], pos$qv$mean[, , t]),
			tolerance = 1e-10
		)
	}
	expect_true(all(rp$U_sd > 0))
	# the posterior sd of a rotated coordinate follows from the rotated covariance
	Q = rp$rotation
	expect_equal(rp$U_sd[2, , 3], sqrt(diag(crossprod(Q, pos$qu$sig[, , 2, 3] %*% Q))), tolerance = 1e-10)
})

test_that("the variance of the product of two independent gaussian positions is the sum of the three covariance terms", {
	set.seed(23)
	R = 2L; Tt = 2L
	spd = function() { a = matrix(rnorm(4, 0, 0.5), 2); crossprod(a) + diag(0.2, 2) }
	Su = spd(); Sv = spd()
	mu = c(0.7, -0.4); mv = c(-0.3, 0.9)
	pos = list(
		qu = list(mean = array(rep(mu, each = 1), c(1, R, Tt)), sig = array(Su, c(R, R, 1, Tt))),
		qv = list(mean = array(rep(mv, each = 1), c(1, R, Tt)), sig = array(Sv, c(R, R, 1, Tt)))
	)
	pos$qu$mean[1, , ] = mu; pos$qv$mean[1, , ] = mv
	vv = lame:::.lde_pos_prodvar(pos, R, Tt)[1L, 1L, 1L]
	nsim = 4e5
	u = matrix(rnorm(nsim * 2), nsim) %*% chol(Su) + rep(mu, each = nsim)
	v = matrix(rnorm(nsim * 2), nsim) %*% chol(Sv) + rep(mv, each = nsim)
	expect_equal(vv, var(rowSums(u * v)), tolerance = 0.02)
	expect_equal(vv, sum(outer(mu, mu) * Sv) + sum(Su * outer(mv, mv)) + sum(Su * Sv), tolerance = 1e-10)
})

test_that("the undirected position sweep updates every actor against the others' current q", {
	set.seed(21)
	n = 5L; R = 2L; Tt = 3L
	sg = c(1, -1)
	ctx = list(
		Tt = Tt,
		na = n,
		nb = n,
		R = R,
		sym = TRUE,
		sgn = sg,
		obs = array(TRUE, c(n, n, Tt))
	)
	for (t in 1:Tt) diag(ctx$obs[, , t]) = FALSE
	wsym = function(a) { for (t in 1:Tt) { m = a[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; a[, , t] = m }; a }
	wk = list(
		w = wsym(array(runif(n * n * Tt, 0.05, 0.4), c(n, n, Tt))),
		z = wsym(array(rnorm(n * n * Tt), c(n, n, Tt)))
	)
	eta_add = array(0, c(n, n, Tt))
	mk = function() {
		sig = array(0, c(R, R, n, Tt)); cl = array(0, c(R, R, n, Tt - 1L))
		for (i in 1:n) for (t in 1:Tt) sig[, , i, t] = crossprod(matrix(rnorm(4, 0, 0.2), 2)) + diag(0.3, R)
		for (i in 1:n) for (t in 1:(Tt - 1L)) cl[, , i, t] = diag(0.2, R)
		list(mean = array(rnorm(n * R * Tt, 0, 0.5), c(n, R, Tt)), sig = sig, cl = cl)
	}
	qin = mk()
	hp = lame:::.lde_pos_hp(0.6, 0.8)
	pos = list(qu = qin, qv = lame:::.lde_pos_mirror(qin, sg))
	out = lame:::.lde_pos_sweep_sym(ctx, pos, hp, wk, eta_add)
	P0 = kronecker(lame:::.lde_pos_q(Tt, hp$rho) / hp$s2, diag(R))
	for (i in 1:n) {
		# actors before i at their updated q, the others at the incoming q
		P = P0
		h = numeric(Tt * R)
		for (j in setdiff(1:n, i)) {
			src = if (j < i) out$qu else qin
			for (t in 1:Tt) {
				ix = (t - 1L) * R + 1:R
				mj = sg * src$mean[j, , t]
				Ej = outer(sg, sg) * (tcrossprod(src$mean[j, , t]) + src$sig[, , j, t])
				P[ix, ix] = P[ix, ix] + wk$w[i, j, t] * Ej
				h[ix] = h[ix] + wk$w[i, j, t] * wk$z[i, j, t] * mj
			}
		}
		S = solve(P)
		expect_equal(as.vector(out$qu$mean[i, , ]), as.vector(S %*% h), tolerance = 1e-10)
		expect_equal(out$qu$sig[, , i, 2], S[R + 1:R, R + 1:R], tolerance = 1e-10)
	}
	# the receiver side is the mirror of the sender side, so the product is F S F'
	O = lame:::.lde_pos_product(out, Tt)
	for (t in 1:Tt) {
		Ft = matrix(out$qu$mean[, , t], n)
		expect_equal(O[, , t], Ft %*% diag(sg, R) %*% t(Ft), tolerance = 1e-10)
		expect_equal(O[, , t], t(O[, , t]), tolerance = 1e-10)
	}
	# one side's statistics and KL for an undirected network
	expect_equal(lame:::.lde_pos_stats(out, R, Tt, sym = TRUE)$nser, n * R)
	expect_equal(
		lame:::.lde_pos_kl2(out, R, Tt, hp, sym = TRUE),
		lame:::.lde_pos_kl2(list(qu = out$qu, qv = out$qu), R, Tt, hp) / 2,
		tolerance = 1e-10
	)
})

test_that("the scale search finds an interior minimum and extends its bracket to a minimum outside it", {
	f = function(s0) function(S) 3 * (log(S) - log(s0))^2 + 5
	g1 = lame:::.lde_golden_log(f(1), 0.08, 4)
	expect_false(g1$edge)
	expect_equal(g1$S[which.min(g1$value)], 1, tolerance = 0.12)
	# a minimiser above the bracket: found by extending the range upward
	g2 = lame:::.lde_golden_log(f(30), 0.08, 4)
	expect_false(g2$edge)
	expect_equal(g2$S[which.min(g2$value)], 30, tolerance = 0.12)
	expect_gt(g2$range[2], 30)
	# and below it
	g3 = lame:::.lde_golden_log(f(0.01), 0.08, 4)
	expect_false(g3$edge)
	expect_equal(g3$S[which.min(g3$value)], 0.01, tolerance = 0.12)
	# a function that keeps falling is reported as bounded, not as an interior minimum
	g4 = lame:::.lde_golden_log(function(S) -log(S), 0.08, 4)
	expect_true(g4$edge)
})

test_that("the numeric leaves of a state list are read and put back in order, and nothing else changes", {
	x = list(a = c(1.5, 2.5, 3.5), b = list(c = matrix(1:4 + 0, 2), d = NULL, e = "text"), f = 7, g = TRUE)
	v = lame:::.lde_leaves(x)
	expect_equal(v, c(1.5, 2.5, 3.5, 1, 2, 3, 4, 7))
	pos = new.env()
	pos$i = 0L
	y = lame:::.lde_leaves_put(x, v * 2, pos)
	expect_equal(pos$i, length(v))
	expect_equal(y$a, x$a * 2)
	expect_equal(y$b$c, x$b$c * 2)
	expect_equal(y$f, 14)
	expect_identical(y$b$e, "text")
	expect_true(y$g)
	expect_true("d" %in% names(y$b) && is.null(y$b$d))
})

test_that("the squared extrapolation reaches the limit of a slowly contracting tail", {
	A = diag(c(0.98, 0.3, 0.1))
	b = c(1, -2, 0.5)
	fixed = solve(diag(3) - A, b)
	F = function(x) as.vector(A %*% x + b)
	x = c(0, 0, 0)
	for (k in 1:20) x = F(x)
	y1 = F(x)
	y2 = F(y1)
	xs = lame:::.lde_squarem(x, y1, y2)
	expect_lt(max(abs(xs - fixed)), 1e-5)
	expect_gt(max(abs(y2 - fixed)), 10)
	# a halved extrapolation lies between the plain double step and the full one
	xh = lame:::.lde_squarem(x, y1, y2, shrink = 0.5)
	expect_lt(max(abs(xh - fixed)), max(abs(y2 - fixed)))
	expect_gt(max(abs(xh - fixed)), max(abs(xs - fixed)))
	# no curvature, no extrapolation
	expect_null(lame:::.lde_squarem(c(1, 1), c(2, 2), c(3, 3)))
})

test_that("the change of the objective still to come is the sum of a geometric tail", {
	obj = 100 + 5 * 0.9^(0:30)
	k = 12
	expect_equal(lame:::.lde_dyn_remaining(obj[(k - 3):k]), obj[k] - 100, tolerance = 1e-6)
	expect_equal(lame:::.lde_dyn_remaining(obj[1:3]), Inf)
	expect_equal(lame:::.lde_dyn_remaining(rep(1, 6)), 0)
	# a change that has stopped shrinking is not settled, nor one that follows a change of zero
	expect_gt(lame:::.lde_dyn_remaining(c(1, 2, 3, 4)), 100)
	expect_gt(lame:::.lde_dyn_remaining(c(1, 2, 2, 5)), 100)
	expect_gt(lame:::.lde_dyn_remaining(c(1, 1, 1, 6)), 100)
})

test_that("mixing two position posteriors mixes means and covariance blocks, and the mirror follows", {
	set.seed(9)
	n = 5L; R = 2L; Tt = 3L
	mk = function() list(
		mean = array(rnorm(n * R * Tt), c(n, R, Tt)),
		sig = array(rnorm(R * R * n * Tt), c(R, R, n, Tt)),
		cl = array(rnorm(R * R * n * (Tt - 1L)), c(R, R, n, Tt - 1L)),
		ldp = 1.5
	)
	p0 = list(qu = mk(), qv = mk()); p1 = list(qu = mk(), qv = mk())
	m0 = lame:::.lde_pos_mix(p0, p1, 0)
	m1 = lame:::.lde_pos_mix(p0, p1, 1)
	for (f in c("mean", "sig", "cl")) {
		expect_equal(m0$qu[[f]], p0$qu[[f]]); expect_equal(m1$qv[[f]], p1$qv[[f]])
	}
	mh = lame:::.lde_pos_mix(p0, p1, 0.25)
	expect_equal(mh$qu$mean, 0.75 * p0$qu$mean + 0.25 * p1$qu$mean)
	expect_equal(mh$qv$cl, 0.75 * p0$qv$cl + 0.25 * p1$qv$cl)
	# an undirected network has one posterior: the receiver side is the sign mirror of the sender side
	sg = c(1, -1)
	ms = lame:::.lde_pos_mix(p0, p1, 0.4, sym = TRUE, sgn = sg)
	expect_equal(ms$qv, lame:::.lde_pos_mirror(ms$qu, sg))
})

test_that("the balance of the two sides keeps the product, minimises the kl divergence along its valley and is idempotent", {
	set.seed(3)
	hp = lame:::.lde_pos_hp(0.6, 0.8)
	# one side: the precision of every actor is the prior plus a block diagonal, as in the sweeps; the dense
	# inverses are kept as the reference
	mkside = function(n, R, Tt, scale) {
		P0 = kronecker(lame:::.lde_pos_q(Tt, hp$rho) / hp$s2, diag(R))
		mean_ = array(rnorm(n * R * Tt, 0, scale), c(n, R, Tt))
		sig = array(0, c(R, R, n, Tt)); cl = array(0, c(R, R, n, max(Tt - 1L, 1L)))
		sgs = vector("list", n); ldp = 0; trp = 0
		for (i in seq_len(n)) {
			P = P0
			for (t in seq_len(Tt)) {
				ix = (t - 1L) * R + seq_len(R)
				P[ix, ix] = P[ix, ix] + crossprod(matrix(rnorm(R * R), R)) * 2 + diag(0.2, R)
			}
			Sg = solve(P); sgs[[i]] = Sg
			ldp = ldp + as.numeric(determinant(P)$modulus); trp = trp + sum(P0 * Sg)
			for (t in seq_len(Tt)) {
				ix = (t - 1L) * R + seq_len(R)
				sig[, , i, t] = Sg[ix, ix]
				if (t < Tt) cl[, , i, t] = Sg[ix, ix + R]
			}
		}
		list(
			q = list(
				mean = mean_,
				sig = sig,
				cl = cl,
				ldp = ldp,
				trp = trp,
				ldp0 = n * as.numeric(determinant(P0)$modulus),
				P0 = P0
			),
			sgs = sgs
		)
	}
	dense_kl2 = function(side, B) {
		qs = side$q; n = dim(qs$mean)[1L]; R = dim(qs$mean)[2L]; Tt = dim(qs$mean)[3L]
		D = kronecker(diag(Tt), B); tot = 0
		for (i in seq_len(n)) {
			Sg = D %*% side$sgs[[i]] %*% t(D); m = D %*% as.vector(qs$mean[i, , ])
			tot = tot + sum(qs$P0 * Sg) + sum(m * (qs$P0 %*% m)) - R * Tt -
				as.numeric(determinant(qs$P0)$modulus) - as.numeric(determinant(Sg)$modulus)
		}
		tot
	}
	# the transformation recovered from the means (one R x R matrix for every actor and period)
	recover = function(q0, q1) {
		n = dim(q0$mean)[1L]; R = dim(q0$mean)[2L]; Tt = dim(q0$mean)[3L]
		X0 = matrix(aperm(q0$mean, c(1, 3, 2)), n * Tt, R); X1 = matrix(aperm(q1$mean, c(1, 3, 2)), n * Tt, R)
		t(solve(crossprod(X0), crossprod(X0, X1)))
	}
	for (cfg in list(c(8, 8, 2, 4), c(8, 6, 2, 3), c(7, 9, 1, 4))) {
		nu = cfg[1]; nv = cfg[2]; R = cfg[3]; Tt = cfg[4]
		su = mkside(nu, R, Tt, 1.0); sv = mkside(nv, R, Tt, 0.3)
		pos = list(qu = su$q, qv = sv$q)
		g = lame:::.lde_pos_gauge(pos, hp, R, Tt)
		# the product and its variance do not change
		expect_equal(lame:::.lde_pos_product(g, Tt), lame:::.lde_pos_product(pos, Tt), tolerance = 1e-12)
		expect_equal(lame:::.lde_pos_prodvar(g, R, Tt), lame:::.lde_pos_prodvar(pos, R, Tt), tolerance = 1e-12)
		# the stored pieces of the KL divergence are those of the dense computation, before and after
		expect_equal(
			lame:::.lde_pos_kl2(pos, R, Tt, hp),
			dense_kl2(su, diag(R)) + dense_kl2(sv, diag(R)),
			tolerance = 1e-10
		)
		Bu = recover(pos$qu, g$qu); Bv = recover(pos$qv, g$qv)
		expect_equal(Bu %*% t(Bv), diag(R), tolerance = 1e-10)
		d1 = dense_kl2(su, Bu) + dense_kl2(sv, Bv)
		expect_equal(lame:::.lde_pos_kl2(g, R, Tt, hp), d1, tolerance = 1e-10)
		expect_lt(d1, lame:::.lde_pos_kl2(pos, R, Tt, hp))
		# a minimum: every perturbation of the transformation raises it
		for (k in 1:10) {
			Bp = diag(R) + matrix(rnorm(R * R, 0, 0.05), R)
			expect_gt(dense_kl2(su, Bp %*% Bu) + dense_kl2(sv, t(solve(Bp)) %*% Bv), d1)
		}
		# idempotent
		g2 = lame:::.lde_pos_gauge(g, hp, R, Tt)
		expect_equal(g2$qu$mean, g$qu$mean, tolerance = 1e-10)
		expect_equal(g2$qv$sig, g$qv$sig, tolerance = 1e-10)
	}
	# the minimiser for sides with as many actors is the geometric mean: G M_u G = M_v
	Mu = crossprod(matrix(rnorm(9), 3)) + diag(3); Mv = crossprod(matrix(rnorm(9), 3)) + diag(3)
	G = lame:::.lde_pos_gauge_g(Mu, Mv, 0)
	expect_equal(G %*% Mu %*% G, Mv, tolerance = 1e-10)
	# and with unequal numbers G M_u G - c G = M_v
	G = lame:::.lde_pos_gauge_g(Mu, Mv, 5)
	expect_equal(G %*% Mu %*% G - 5 * G, Mv, tolerance = 1e-10)
})

test_that("the precision of an ar(1) path with weighted steps is the inverse of its covariance, and the sweep with weights is the dense posterior", {
	set.seed(17)
	Tt = 5L; rho = 0.7
	w = c(1, runif(Tt - 1L, 0.2, 3))
	# the covariance of the path: the first period is stationary, the step into period t has the variance 1 / w_t
	B = diag(Tt)
	cv = matrix(0, Tt, Tt)
	for (a in seq_len(Tt)) for (b in seq_len(Tt)) {
		# the covariance of x_a and x_b by the recursion x_t = rho x_t-1 + e_t / sqrt(w_t)
		v = numeric(Tt); v[1] = 1 / (1 - rho^2)
		for (t in 2:Tt) v[t] = rho^2 * v[t - 1] + 1 / w[t]
		if (a <= b) cv[a, b] = rho^(b - a) * v[a] else cv[a, b] = rho^(a - b) * v[b]
	}
	Q = lame:::.lde_pos_qw(Tt, rho, w)
	expect_equal(Q, solve(cv), tolerance = 1e-10)
	expect_equal(determinant(Q)$modulus[[1]], log(1 - rho^2) + sum(log(w[-1L])), tolerance = 1e-10)
	expect_equal(lame:::.lde_pos_qw(Tt, rho, rep(1, Tt)), lame:::.lde_pos_q(Tt, rho), tolerance = 1e-14)
	# the diagonal and off-diagonal of every actor at once
	lam = cbind(1, matrix(runif(3 * (Tt - 1L), 0.2, 3), 3))
	qc = lame:::.lde_pos_qcoef(3L, Tt, rho, lam)
	for (i in 1:3) {
		Qi = lame:::.lde_pos_qw(Tt, rho, lam[i, ])
		expect_equal(qc$d[i, ], diag(Qi), tolerance = 1e-14)
		expect_equal(qc$e[i, ], Qi[cbind(1:(Tt - 1L), 2:Tt)], tolerance = 1e-14)
	}
	# the sweep: the dense posterior of every actor under its own prior precision, in R and compiled
	n = 4L; m = 5L; R = 2L; Tt = 4L
	Wl = lapply(seq_len(Tt), function(t) matrix(runif(n * m, 0.05, 0.4), n, m))
	Rl = lapply(seq_len(Tt), function(t) matrix(rnorm(n * m), n, m))
	Vh = lapply(seq_len(Tt), function(t) matrix(rnorm(m * R, 0, 0.6), m, R))
	Om = lapply(seq_len(Tt), function(t) {
		a = array(0, c(R, R, m))
		for (j in seq_len(m)) a[, , j] = crossprod(matrix(rnorm(4, 0, 0.2), 2)) + diag(0.05, R)
		a
	})
	lam = cbind(1, matrix(runif(n * (Tt - 1L), 0.2, 3), n))
	s2 = 0.12
	out = lame:::.lde_pos_side(Wl, Rl, Vh, Om, s2, rho, n, R, Tt, lam = lam)
	for (i in seq_len(n)) {
		P = kronecker(lame:::.lde_pos_qw(Tt, rho, lam[i, ]) / s2, diag(R))
		h = numeric(Tt * R)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			for (j in seq_len(m)) {
				P[ix, ix] = P[ix, ix] + Wl[[t]][i, j] * (tcrossprod(Vh[[t]][j, ]) + Om[[t]][, , j])
				h[ix] = h[ix] + Wl[[t]][i, j] * Rl[[t]][i, j] * Vh[[t]][j, ]
			}
		}
		S = solve(P)
		expect_equal(as.vector(out$mean[i, , ]), as.vector(S %*% h), tolerance = 1e-10)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			expect_equal(out$sig[, , i, t], S[ix, ix], tolerance = 1e-10)
			if (t < Tt) expect_equal(out$cl[, , i, t], S[ix, ix + R], tolerance = 1e-10)
		}
	}
	# the compiled sweep with weights equals the reference
	dims = c(n, m, Tt)
	ctx = list(sym = FALSE, na = n, nb = m, R = R, Tt = Tt, obs = array(TRUE, dims))
	side = function(k) {
		sg = array(0, c(R, R, k, Tt)); cl = array(0, c(R, R, k, Tt - 1L))
		for (i in seq_len(k)) for (t in seq_len(Tt)) sg[, , i, t] = crossprod(matrix(rnorm(4), 2)) * 0.1 + diag(0.05, R)
		list(mean = array(rnorm(k * R * Tt), c(k, R, Tt)), sig = sg, cl = cl)
	}
	pos = list(qu = side(n), qv = side(m))
	pos$qu$lam = lam
	pos$qv$lam = cbind(1, matrix(runif(m * (Tt - 1L), 0.2, 3), m))
	wk = list(w = array(runif(prod(dims), 0.05, 0.4), dims), z = array(rnorm(prod(dims)), dims))
	eta = array(rnorm(prod(dims), -1, 0.5), dims)
	hp = lame:::.lde_pos_hp(0.7, 0.85)
	a = lame:::.lde_pos_sweep_ref(ctx, pos, hp, wk, eta)
	b = lame:::.lde_pos_sweep(ctx, pos, hp, wk, eta)
	for (sd in c("qu", "qv")) {
		for (f in c("mean", "sig", "cl", "ldp", "trp", "lam")) expect_equal(
			b[[sd]][[f]],
			a[[sd]][[f]],
			tolerance = 1e-10
		)
	}
})

test_that("the mixing weights take the variational update, and the kl divergence of the scale mixture is the gaussian part and the gamma part", {
	set.seed(19)
	n = 4L; R = 2L; Tt = 4L; nu = 4
	rho = 0.8; s2 = 0.3
	hp = lame:::.lde_pos_hp(s2 / (1 - rho^2), rho)
	# an actor's joint posterior of the path, dense; the blocks the sweeps keep are its diagonal and lag-one blocks
	mk = function() {
		mean_ = array(rnorm(n * R * Tt, 0, 0.6), c(n, R, Tt))
		sig = array(0, c(R, R, n, Tt)); cl = array(0, c(R, R, n, Tt - 1L)); ful = vector("list", n)
		for (i in seq_len(n)) {
			A = matrix(rnorm((R * Tt)^2, 0, 0.2), R * Tt); P = crossprod(A) + diag(0.3, R * Tt)
			ful[[i]] = P
			for (t in seq_len(Tt)) {
				ix = (t - 1L) * R + seq_len(R)
				sig[, , i, t] = P[ix, ix]
				if (t < Tt) cl[, , i, t] = P[ix, ix + R]
			}
		}
		list(mean = mean_, sig = sig, cl = cl, ful = ful)
	}
	q = mk()
	# E |u_t - rho u_t-1|^2 by the definition: the mean square plus the trace of the covariance of the step
	ee = lame:::.lde_pos_eesq(q, rho)
	for (i in seq_len(n)) for (t in 2:Tt) {
		ix1 = (t - 2L) * R + seq_len(R); ix2 = (t - 1L) * R + seq_len(R)
		Dm = cbind(-rho * diag(R), diag(R)); idx = c(ix1, ix2)
		step_cov = Dm %*% q$ful[[i]][idx, idx] %*% t(Dm)
		step_mean = q$mean[i, , t] - rho * q$mean[i, , t - 1L]
		expect_equal(ee$ee[i, t - 1L], sum(step_mean^2) + sum(diag(step_cov)), tolerance = 1e-10)
	}
	expect_equal(
		ee$e1,
		vapply(seq_len(n), function(i) sum(q$mean[i, , 1L]^2) + sum(diag(q$ful[[i]][1:R, 1:R])), 0),
		tolerance = 1e-10
	)
	# the weights: the sampler's full conditional gamma((nu + R) / 2, (nu + |step|^2 / s2) / 2) with the step replaced by its expectation
	lu = lame:::.lde_pos_lambda(q, hp, nu, R)
	expect_equal(lu$lam[, 1L], rep(1, n))
	expect_equal(lu$lam[, -1L], (nu + R) / (nu + ee$ee / s2), tolerance = 1e-12)
	expect_equal(lu$elog[, -1L], digamma((nu + R) / 2) - log(lu$b[, -1L]), tolerance = 1e-12)
	# a step that is large against the innovation scale gets a small weight
	expect_true(all(diff(lu$lam[order(ee$ee[, 1L]), 2L]) <= 0))
	# the KL divergence: the Gaussian part, E log |P0| at E log lambda and the energy at E lambda, plus the gamma KL
	q$lam = lu$lam; q$elog = lu$elog; q$b = lu$b; q$ldp = 5.3
	pos = list(qu = q, qv = q)
	kl = lame:::.lde_pos_kl2(pos, R, Tt, hp, t_nu = nu)
	one = 0
	for (i in seq_len(n)) {
		P0 = kronecker(lame:::.lde_pos_qw(Tt, rho, q$lam[i, ]) / s2, diag(R))
		Sg = q$ful[[i]]; m = as.vector(q$mean[i, , ])
		# E log |P0(lambda)| = log |P0(E lambda)| + R sum_t (E log lambda_t - log E lambda_t)
		ldp0 = as.numeric(determinant(P0)$modulus) + R * sum(q$elog[i, -1L] - log(q$lam[i, -1L]))
		one = one + sum(P0 * Sg) + sum(m * (P0 %*% m)) - R * Tt - ldp0
	}
	gam = 0
	for (i in seq_len(n)) for (t in 2:Tt) {
		a = (nu + R) / 2; b = q$b[i, t]
		f = function(x) dgamma(x, a, b) * (dgamma(x, a, b, log = TRUE) - dgamma(x, nu / 2, nu / 2, log = TRUE))
		gam = gam + integrate(f, 0, Inf, rel.tol = 1e-10)$value
	}
	# twice the KL divergence of each of the two identical sides: the Gaussian part, the log determinant of q's precision and twice the gamma KL
	expect_equal(kl, 2 * (one + q$ldp + 2 * gam), tolerance = 1e-6)
	# and with every weight one the Gaussian KL: the same function without weights
	q1 = q; q1$lam = NULL
	expect_equal(
		lame:::.lde_pos_kl2(list(qu = q1, qv = q1), R, Tt, hp),
		2 * lame:::.lde_pos_kl2(list(qu = q1), R, Tt, hp, sym = TRUE),
		tolerance = 1e-12
	)
})

test_that("the student-t transition is fitted, reports its mixing weights and tends to the gaussian fit as its degrees of freedom grow", {
	skip_on_cran()
	d = eb_sim(n = 10, Tt = 4, seed = 8, mu = -0.8)
	args = list(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		method = "als",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		verbose = FALSE
	)
	ft = suppressWarnings(suppressMessages(do.call(lame, c(args, list(dynamic_uv_kind = "t")))))
	expect_true(isTRUE(ft$converged))
	expect_identical(ft$dynamic_uv_kind, "t")
	expect_identical(ft$eb$uv_kind, "t")
	expect_equal(ft$eb$t_nu, 4)
	# the weights: one in the first period, positive, and not all one
	expect_equal(dim(ft$lambda_u), c(10, 4))
	expect_equal(dim(ft$lambda_v), c(10, 4))
	expect_true(all(ft$lambda_u[, 1L] == 1) && all(ft$lambda_u > 0))
	expect_gt(diff(range(ft$lambda_u[, -1L])), 0.01)
	# with a very large number of degrees of freedom the scale mixture is the Gaussian transition
	f1 = suppressWarnings(suppressMessages(do.call(lame, c(args, list(dynamic_uv_kind = "ar1")))))
	expect_null(f1$lambda_u)
	fl = suppressWarnings(
		suppressMessages(do.call(lame, c(args, list(dynamic_uv_kind = "t", prior = list(t_nu = 1e6)))))
	)
	expect_equal(unname(fl$coefficients), unname(f1$coefficients), tolerance = 1e-5)
	expect_equal(fl$eb$uv_scale, f1$eb$uv_scale, tolerance = 1e-4)
	expect_equal(fl$rho_uv, f1$rho_uv, tolerance = 1e-4)
	expect_equal(range(fl$lambda_u), c(1, 1), tolerance = 1e-4)
	# a restart record of one kind does not seed the other
	prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	c1 = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE)
	c2 = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE, dynamic_uv_kind = "t")
	expect_false(identical(lame:::.lde_signature(c1), lame:::.lde_signature(c2)))
	# an undirected network shares the weights of its two sides
	ds = eb_sim(n = 10, Tt = 4, mode = "symmetric", seed = 9, mu = -0.8)
	fs = suppressWarnings(suppressMessages(lame(
		ds$Y,
		Xdyad = ds$X,
		family = "binary",
		method = "als",
		symmetric = TRUE,
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = "t",
		verbose = FALSE
	)))
	expect_equal(fs$lambda_u, fs$lambda_v)
})

test_that("the balance of the two sides is the minimum of the kl divergence under weighted steps", {
	set.seed(23)
	R = 2L; Tt = 4L; rho = 0.8; s2 = 0.4
	hp = lame:::.lde_pos_hp(s2 / (1 - rho^2), rho)
	mkside = function(n, scale) {
		lam = cbind(1, matrix(runif(n * (Tt - 1L), 0.3, 2.5), n))
		mean_ = array(rnorm(n * R * Tt, 0, scale), c(n, R, Tt))
		sig = array(0, c(R, R, n, Tt)); cl = array(0, c(R, R, n, Tt - 1L)); ful = vector("list", n)
		for (i in seq_len(n)) {
			A = matrix(rnorm((R * Tt)^2, 0, 0.3), R * Tt); P = crossprod(A) + diag(0.2, R * Tt)
			ful[[i]] = P
			for (t in seq_len(Tt)) {
				ix = (t - 1L) * R + seq_len(R)
				sig[, , i, t] = P[ix, ix]
				if (t < Tt) cl[, , i, t] = P[ix, ix + R]
			}
		}
		list(mean = mean_, sig = sig, cl = cl, lam = lam, elog = log(lam), b = matrix(1, n, Tt), ldp = 0, ful = ful)
	}
	# the energy of a side after the change u -> B u, by the dense matrices, at fixed weights
	energy = function(side, B) {
		D = kronecker(diag(Tt), B); tot = 0
		for (i in seq_len(dim(side$mean)[1L])) {
			P0 = kronecker(lame:::.lde_pos_qw(Tt, rho, side$lam[i, ]) / s2, diag(R))
			Sg = D %*% side$ful[[i]] %*% t(D); m = D %*% as.vector(side$mean[i, , ])
			tot = tot + sum(P0 * Sg) + sum(m * (P0 %*% m))
		}
		tot
	}
	su = mkside(7L, 1); sv = mkside(9L, 0.3)
	pos = list(qu = su, qv = sv)
	g = lame:::.lde_pos_gauge(pos, hp, R, Tt)
	recover = function(q0, q1) {
		n = dim(q0$mean)[1L]
		X0 = matrix(aperm(q0$mean, c(1, 3, 2)), n * Tt, R); X1 = matrix(aperm(q1$mean, c(1, 3, 2)), n * Tt, R)
		t(solve(crossprod(X0), crossprod(X0, X1)))
	}
	Bu = recover(su, g$qu); Bv = recover(sv, g$qv)
	expect_equal(Bu %*% t(Bv), diag(R), tolerance = 1e-10)
	# the change of the energy plus the change of the entropy: -log |B| per period and actor on each side, with the
	# sign that the sender (positive) and receiver (negative) log determinants carry
	ent = function(Bp) {
		Bu2 = Bp %*% Bu; Bv2 = t(solve(Bp)) %*% Bv
		energy(su, Bu2) + energy(sv, Bv2) - 2 * 7 * Tt * log(abs(det(Bu2))) - 2 * 9 * Tt * log(abs(det(Bv2)))
	}
	f0 = ent(diag(R))
	for (k in 1:10) {
		Bp = diag(R) + matrix(rnorm(R * R, 0, 0.05), R)
		expect_gt(ent(Bp), f0 - 1e-9)
	}
})

test_that("a position update that goes back and forth is damped and the iteration converges", {
	skip_on_cran()
	# a strongly identified small panel: the plain map of the sweeps alternates and does not settle (its merit swings
	# by tens of percent from pass to pass, far above the threshold)
	set.seed(9)
	n = 10L; Tt = 4L
	ar1 = function() {
		out = matrix(0, n, Tt)
		out[, 1] = rnorm(n, 0, 1 / sqrt(1 - 0.95^2))
		for (t in 2:Tt) out[, t] = 0.95 * out[, t - 1] + rnorm(n)
		out
	}
	u = ar1(); v = ar1()
	nm = paste0("a", sprintf("%02d", seq_len(n)))
	Y = lapply(seq_len(Tt), function(t) {
		y = 1 * (-0.5 + tcrossprod(u[, t], v[, t]) + matrix(rnorm(n * n), n) > 0)
		dimnames(y) = list(nm, nm); diag(y) = NA
		y
	})
	names(Y) = paste0("t", seq_len(Tt))
	prep = lame:::.ae_prepare(Y, NULL, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ctx = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE)
	fit = lame:::.lde_fit(ctx, lame:::.lde_init(ctx), max_iter = 8L, tol = 5e-5)
	# held at the stationary variance 10 and persistence 0.95 of the strongly identified panel
	fit$st$hyp$uv = lame:::.lde_pos_hp(10, 0.95)
	# the first two plain moves point in opposite directions
	e0 = lame:::.lde_eval_dyn(fit$ctx, fit$st, 5e-5)
	e1 = lame:::.lde_eval_dyn(fit$ctx, e0$cand, 5e-5)
	d1 = lame:::.lde_dyn_vec(e0$cand) - lame:::.lde_dyn_vec(fit$st)
	d2 = lame:::.lde_dyn_vec(e1$cand) - lame:::.lde_dyn_vec(e0$cand)
	expect_lt(sum(d1 * d2) / sqrt(sum(d1^2) * sum(d2^2)), -0.5)
	# the iteration detects it at the first pass (the half step, no extrapolation tried), damps the whole state
	# (the positions included) and converges
	r = lame:::.lde_fit_core_dyn(fit$ctx, fit$st, max_iter = 300L, tol = 5e-5)
	expect_equal(r$trace$step[2], 0.5)
	expect_false(r$trace$jump[2])
	expect_true(r$converged)
	expect_lt(r$ev$cert$merit, sqrt(5e-5))
	expect_true(any(r$trace$step < 1, na.rm = TRUE))
	# the move vector of the relaxation carries the means of the positions, and nothing else changes
	hv = lame:::.lde_hvec(fit$ctx, fit$st)
	st0 = fit$st
	st0$pos = NULL
	hv0 = lame:::.lde_hvec(fit$ctx, st0)
	expect_equal(length(hv) - length(hv0), length(fit$st$pos$qu$mean) + length(fit$st$pos$qv$mean))
	expect_equal(hv[seq_along(hv0)], hv0)
})

test_that("moving positions: a small directed panel is fitted by the default route and reported", {
	# a small panel's product is withheld below its information threshold, so the rule is off here
	withr::local_options(list(lame.eb_withhold_product = FALSE))
	skip_on_cran()
	d = eb_sim(n = 10, Tt = 4, seed = 8, mu = -0.8)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		method = "als",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		verbose = FALSE
	))
	expect_true(isTRUE(fit$dynamic_uv))
	expect_identical(fit$dynamic_uv_kind, "ar1")
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(fit$converged)
	expect_equal(dim(fit$U), c(10L, 1L, 4L))
	expect_equal(dim(fit$V), c(10L, 1L, 4L))
	expect_equal(dim(fit$U_sd), dim(fit$U))
	expect_equal(dim(fit$Oarr), c(10L, 10L, 4L))
	expect_equal(fit$UVPM, fit$Oarr)
	# the reported product is the product of the reported positions in every period
	for (t in 1:4) expect_equal(fit$Oarr[, , t], tcrossprod(fit$U[, , t], fit$V[, , t]), ignore_attr = TRUE)
	# the position hyperparameters and the scale search are reported
	expect_true(is.finite(fit$rho_uv) && abs(fit$rho_uv) < 1)
	expect_true(is.finite(fit$sigma_uv) && fit$sigma_uv > 0)
	expect_identical(fit$eb$rho_uv, fit$rho_uv)
	expect_true(fit$eb$uv_scale > 0.08 && fit$eb$uv_scale < 4)
	expect_true(is.data.frame(fit$eb$uv_search) && nrow(fit$eb$uv_search) >= 4L)
	expect_true(all(is.finite(fit$eb$se)))
	# the stationarity check covers the position step
	expect_true("positions" %in% names(fit$diagnostics$stationarity$by_block))
	expect_true(fit$diagnostics$stationarity$pass)
	# the linear predictor carries the product
	eta = fit$EZ[[3]]
	expect_equal(
		eta[2, 5],
		fit$Oarr[2, 5, 3] + fit$a_dynamic[2, 3] + fit$b_dynamic[5, 3] +
			sum(fit$coef_path[, 3] * c(1, d$X[[3]][2, 5, ])),
		tolerance = 1e-8
	)
})

test_that("verbose = TRUE reports the position scale search as it goes", {
	skip_on_cran()
	d = eb_sim(n = 8, Tt = 3, seed = 6, mu = -0.8)
	msgs = testthat::capture_messages(suppressWarnings(
		lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			R = 1,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			verbose = TRUE
		)
	))
	expect_true(any(grepl("Position scale S", msgs, fixed = TRUE)))
	expect_gte(sum(grepl("Position scale S", msgs, fixed = TRUE)), 6L)
})

test_that("a position scale at the edge of the searched range is refused, with the state in the condition", {
	skip_on_cran()
	d = eb_sim(n = 8, Tt = 3, seed = 6, mu = -0.8)
	real = lame:::.lde_pos_profile
	# the upper edge: the objective keeps falling as the scale grows
	testthat::local_mocked_bindings(
		.lde_pos_profile = function(...) { pr = real(...); pr$at_edge = TRUE; pr$S = pr$range[2L]; pr },
		.package = "lame"
	)
	err = expect_error(
		suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			R = 1,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			verbose = FALSE
		)),
		class = "lame_dynamic_uv_scale_edge"
	)
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(msg, "does not identify the stationary variance of the moving positions")
	expect_match(msg, "keeps falling up to the largest scale searched")
	expect_match(msg, "reproduce the panel exactly")
	expect_false(grepl("with care", msg))
	expect_s3_class(err$fit, "lame_dynamic_als")
	expect_true(is.finite(err$fit$eb$uv_scale))
})

test_that("a position scale at the lower edge is refused with its own statement", {
	skip_on_cran()
	d = eb_sim(n = 8, Tt = 3, seed = 6, mu = -0.8)
	real = lame:::.lde_pos_profile
	testthat::local_mocked_bindings(
		.lde_pos_profile = function(...) { pr = real(...); pr$at_edge = TRUE; pr$S = pr$range[1L]; pr },
		.package = "lame"
	)
	err = expect_error(
		suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			R = 1,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			verbose = FALSE
		)),
		class = "lame_dynamic_uv_scale_edge"
	)
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(msg, "lowest at the smallest scale searched")
	expect_match(msg, "no variation that the positions can resolve")
})

test_that("a scale search that did not settle near its minimum is refused, with the state in the condition", {
	skip_on_cran()
	d = eb_sim(n = 8, Tt = 3, seed = 6, mu = -0.8)
	real = lame:::.lde_pos_profile
	testthat::local_mocked_bindings(
		.lde_pos_profile = function(...) { pr = real(...); pr$unsettled = TRUE; pr$grid$converged[which.min(
			abs(log(pr$grid$S) - log(pr$S))
		)] = FALSE; pr },
		.package = "lame"
	)
	err = expect_error(
		suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			R = 1,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			verbose = FALSE
		)),
		class = "lame_dynamic_uv_scale_unsettled"
	)
	# a failure of the iteration is also a convergence failure
	expect_s3_class(err, "lame_als_convergence")
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(msg, "did not settle at 1 of the 5 position scales nearest the minimum of the search")
	expect_s3_class(err$fit, "lame_dynamic_als")
})

test_that("a restart record belongs to the problem it was made for", {
	d = eb_sim(n = 8, Tt = 3, seed = 6, mu = -0.8)
	prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ctx = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE)
	sg = lame:::.lde_signature(ctx)
	expect_identical(sg, lame:::.lde_signature(ctx))
	for (what in list(
		list(link = "logit"),
		list(beta_pool = "sigma"),
		list(coef_prior = "flat"),
		list(R = 2L),
		list(Tt = 4L)
	)) {
		other = ctx
		other[names(what)] = what
		expect_false(identical(sg, lame:::.lde_signature(other)), info = names(what))
	}
})

test_that("specifications the empirical-Bayes estimator cannot validate are refused by name", {
	d = eb_sim(n = 10, Tt = 3, seed = 4)
	args = list(d$Y, Xdyad = d$X, family = "binary", method = "als", verbose = FALSE)
	# the heavy-tailed kind of moving latent positions is not offered for a bipartite network (the interaction
	# U G V' leaves the size of the positions to G); the error names the way out
	db = eb_sim(n = 8, nb = 6, Tt = 3, mode = "bipartite", seed = 4)
	expect_error(
		lame(
			db$Y,
			Xdyad = db$X,
			family = "binary",
			mode = "bipartite",
			method = "als",
			R = 1,
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			dynamic_ab = TRUE,
			verbose = FALSE
		),
		"bipartite"
	)
	expect_error(
		lame:::.lde_check_supported(
			structure(list(bip = TRUE), class = "list"),
			FALSE,
			1L,
			1L,
			1L,
			TRUE,
			TRUE,
			"t",
			list(any = FALSE),
			"ar1",
			NULL,
			NULL,
			NULL,
			list(),
			"probit"
		),
		class = "lame_dynamic_eb_unsupported"
	)
	# the cell problem under the positions' noise is built for the probit link of the sampler
	expect_error(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			R = 1,
			dynamic_uv = TRUE,
			dynamic_ab = TRUE,
			link = "logit",
			verbose = FALSE
		),
		class = "lame_dynamic_eb_unsupported"
	)
	# a Matern path has no dynamic ALS route at all, for any family
	expect_error(
		do.call(
			lame,
			c(
				args,
				list(
					dynamic_ab = TRUE, dynamic_beta = "x1_dyad",
					dynamic_beta_kind = "matern32"
				)
			)
		),
		"matern32"
	)
	# fixed smoothness penalties belong to the uncorrected model
	expect_error(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			dynamic_ab = TRUE,
			lambda_ab = 2,
			verbose = FALSE
		),
		class = "lame_dynamic_eb_unsupported"
	)
	# a bootstrap is not offered; the estimator carries Laplace standard errors
	expect_error(
		do.call(lame, c(args, list(dynamic_ab = TRUE, bootstrap = 5L))),
		class = "lame_dynamic_eb_unsupported"
	)
	# the other corrections stay static or Poisson only
	expect_error(
		do.call(lame, c(args, list(dynamic_ab = TRUE, bias_correct = "bootstrap"))),
		class = "lame_jk_split"
	)
	expect_error(
		do.call(lame, c(args, list(dynamic_ab = TRUE, bias_correct = "jackknife"))),
		class = "lame_jk_split"
	)
	# the empirical-Bayes estimator is defined for the normal, binary and Poisson likelihoods: a normal panel is a fit, and the
	# bootstrap and the fixed penalties are refused for it as for the others
	dp = d
	dp$Y = lapply(dp$Y, function(y) y * 2 + 0.5)
	fn = lame:::lame_dynamic_als(dp$Y, dp$X, family = "normal", dynamic_ab = TRUE, bias_correct = "eb", verbose = FALSE)
	expect_true(fn$converged)
	expect_identical(fn$family, "normal")
	expect_true(isTRUE(fn$bias_correct$dynamic))
	expect_error(
		lame:::lame_dynamic_als(dp$Y, dp$X, family = "normal", dynamic_ab = TRUE, bootstrap = 5L, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
	expect_error(
		lame:::lame_dynamic_als(dp$Y, dp$X, family = "normal", dynamic_ab = TRUE, lambda_ab = 2, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
})

test_that("the cavity rule's scale update is the root of its residual, and an uninformed prior goes to psi / df", {
	set.seed(5)
	na = 12L; R = 2L
	ctx = list(R = R, sym = FALSE, kind = "directed", na = na, nb = na, suv0 = 1, fscale = "cavity")
	mkq = function(m) list(
		mean = array(m, c(na, R, 1L)),
		sig = array(rep(diag(0.1, R), na), c(R, R, na, 1L)),
		cl = array(0, c(R, R, na, 1L))
	)
	pos = list(qu = mkq(matrix(rnorm(na * R, 0, 0.7), na, R)), qv = mkq(matrix(rnorm(na * R, 0, 0.7), na, R)))
	prior = lame:::.lde_fac_prior_init(ctx, pos$qu$mean[, , 1L], pos$qv$mean[, , 1L])
	expect_identical(prior$kind, "joint")
	expect_equal(prior$Suv[1L, 1L], prior$Suv[2L, 2L])
	C = matrix(runif(na * na, 0.1, 0.5), na, na)
	diag(C) = 0
	# repeated updates reach the root: the plug-in at the prior returned reproduces its scale
	for (i in 1:8) prior = lame:::.lde_fac_scale_update(ctx, pos, prior, C)$prior
	S = lame:::.lde_fac_state(ctx, pos, prior)
	tg = lame:::.ae_mf_scale_target(list(R = R, sym = FALSE), S, C)
	expect_equal(tg$s, prior$Suv[1L, 1L], tolerance = 1e-6)
	expect_equal(tg$rho, prior$Suv[1L, R + 1L] / prior$Suv[1L, 1L], tolerance = 1e-3)
	# no information about any row and zero means: the update has the fixed point psi / df (suv0 = 1: the null
	# branch, where a panel without structure must not run away)
	pos0 = list(qu = mkq(matrix(0, na, R)), qv = mkq(matrix(0, na, R)))
	pr0 = lame:::.lde_fac_prior_init(ctx, pos0$qu$mean[, , 1L] + 1, pos0$qv$mean[, , 1L] + 1)
	for (i in 1:8) pr0 = lame:::.lde_fac_scale_update(ctx, pos0, pr0, C * 0)$prior
	expect_equal(pr0$Suv[1L, 1L], 1, tolerance = 1e-6)
	# the undirected and bipartite kinds: a scale for each column
	ctx_s = list(R = R, sym = TRUE, kind = "symmetric", na = na, nb = na, suv0 = 1, fscale = "cavity", sgn = c(1, -1))
	ps = lame:::.lde_fac_prior_init(ctx_s, pos$qu$mean[, , 1L], pos$qu$mean[, , 1L])
	expect_identical(ps$kind, "columns")
	expect_length(ps$tau, R)
	ctx_b = list(R = R, sym = FALSE, kind = "bipartite", na = na, nb = na, suv0 = 1, fscale = "cavity")
	pb = lame:::.lde_fac_prior_init(ctx_b, pos$qu$mean[, , 1L], pos$qv$mean[, , 1L])
	expect_identical(pb$kind, "diag")
	expect_equal(pb$sv, rep(1, R))
})

test_that("a fit with static factors leaves its growth phase and stays at its fixed point", {
	skip_on_cran()
	# a panel whose first extrapolation comes while the factor scale is still growing out of its start
	d = eb_sim_fac(n = 20, Tt = 4, seed = 31, mu = -1.2, R = 1)
	fit = suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		dynamic_ab = TRUE,
		R = 1L,
		verbose = FALSE,
		max_iter = 1000L
	)))
	expect_true(fit$converged)
	expect_lt(fit$iterations, 150L)
	# the record of a converged fit is its own fixed point
	g = suppressWarnings(suppressMessages(ame_als_refit(fit)))
	expect_equal(g$iterations, 0)
	expect_equal(unname(coef(g)), unname(coef(fit)), tolerance = 1e-8)
	# no extrapolation made the objective worse than the plain pass by more than the range it had covered
	tr = lame:::.lde_fit(
		lame:::.lde_context(
			lame:::.ae_prepare(
				d$Y,
				d$X,
				NULL,
				NULL,
				"unipartite",
				longitudinal = TRUE,
				node_time_policy = "quiet"
			),
			FALSE,
			"probit",
			dynamic_ab = TRUE,
			R = 1L
		),
		max_iter = 300L
	)$trace
	expect_lt(max(tr$objective) - min(tr$objective), 0.25 * abs(tr$objective[1L]))
})

test_that("static factors: the priors are the sampler's, the scale carries its uncertainty, and options(lame.eb_scale) selects the rule", {
	skip_on_cran()
	options_w = options(lame.eb_withhold_intercept = FALSE)
	on.exit(options(options_w), add = TRUE)
	old = options(lame.eb_scale = "cavity", lame.eb_scale_uncertainty = TRUE)
	on.exit(options(old), add = TRUE)
	d = eb_sim_fac(n = 14, Tt = 3, seed = 8, R = 1)
	fc = eb_fit(d, R = 1)
	expect_identical(fc$eb$factor_scale$rule, "cavity")
	Suv = fc$eb$factor_prior$Suv
	expect_equal(Suv[1L, 1L], Suv[2L, 2L])
	expect_null(attr(fc$eb$factor_prior, "acc"))
	expect_true(is.finite(fc$eb$factor_scale$sd_log) && fc$eb$factor_scale$sd_log > 0)
	# the standard errors carry the uncertainty of the scale
	options(lame.eb_scale_uncertainty = FALSE)
	fn = eb_fit(d, R = 1)
	expect_true(is.na(fn$eb$factor_scale$sd_log))
	expect_true(all(fn$eb$se <= fc$eb$se + 1e-12) && any(fn$eb$se < fc$eb$se - 1e-8))
	expect_equal(unname(coef(fn)), unname(coef(fc)), tolerance = 1e-8)
	# the plug-in rule
	options(lame.eb_scale = "plugin", lame.eb_scale_uncertainty = TRUE)
	fp = eb_fit(d, R = 1)
	expect_identical(fp$eb$factor_scale$rule, "plugin")
	expect_identical(fp$eb$factor_prior$kind, "joint")
	expect_true(is.na(fp$eb$factor_scale$sd_log))
	options(lame.eb_scale = "nonsense")
	expect_error(eb_fit(d, R = 1), "cavity")
	# bipartite: the receiver factor stays at the unit scale; undirected: a scale for every column of F
	options(lame.eb_scale = "cavity")
	db = eb_sim(n = 12, nb = 8, Tt = 3, mode = "bipartite", seed = 4)
	fb = suppressMessages(suppressWarnings(lame(
		db$Y,
		Xdyad = db$X,
		family = "binary",
		method = "als",
		mode = "bipartite",
		dynamic_ab = TRUE,
		R = 1,
		verbose = FALSE
	)))
	expect_identical(fb$eb$factor_prior$kind, "diag")
	expect_equal(fb$eb$factor_prior$sv, 1)
	ds = eb_sim(n = 12, Tt = 3, mode = "symmetric", seed = 5)
	fs = suppressMessages(suppressWarnings(lame(
		ds$Y,
		Xdyad = ds$X,
		family = "binary",
		method = "als",
		symmetric = TRUE,
		dynamic_ab = TRUE,
		R = 2,
		verbose = FALSE
	)))
	expect_identical(fs$eb$factor_prior$kind, "columns")
	expect_length(fs$eb$factor_prior$tau, 2L)
})

test_that("a fit with static factors and a logit link is a logit fit", {
	skip_on_cran()
	d = eb_sim_fac(n = 16, Tt = 3, seed = 12, R = 1)
	fp = eb_fit(d, R = 1, link = "probit")
	fl = eb_fit(d, R = 1, link = "logit")
	expect_true(fl$converged)
	# the logit coefficients are about 1.7 times the probit ones, and not the same fit relabelled
	r = coef(fl) / coef(fp)
	expect_true(all(r > 1.3 & r < 2.2))
	expect_gt(max(abs(coef(fl) - coef(fp))), 0.1)
})

# the dyadic correlation of a directed network (lame()'s dcor = TRUE)

test_that("the pair node of the shared quadrature: the Fisher information, its derivatives and its derivatives in the correlation", {
	set.seed(1)
	n = 150
	a = rnorm(n, -1, 1.5); b = rnorm(n, -1, 1.5); r = runif(n, -0.9, 0.9)
	s1 = sample(c(-1, 1), n, TRUE); s2 = sample(c(-1, 1), n, TRUE)
	nd = lame:::.ae_pair_node_fisher_r(a, b, r, s1, s2, TRUE)
	ref = lame:::.ae_pair_logp(s1 * a, s2 * b, s1 * s2 * r)
	expect_equal(nd$lo, ref$L, tolerance = 1e-9)
	expect_equal(nd$psi1, s1 * ref$La, tolerance = 1e-9)
	expect_equal(nd$psi2, s2 * ref$Lb, tolerance = 1e-9)
	# the Fisher information is the expectation over the four outcomes of the outer product of the score
	I11 = I22 = I12 = numeric(n)
	for (g1 in c(1, -1)) for (g2 in c(1, -1)) {
		dd = lame:::.ae_pair_logp(g1 * a, g2 * b, g1 * g2 * r)
		l1 = g1 * dd$La; l2 = g2 * dd$Lb
		I11 = I11 + exp(dd$L) * l1^2; I22 = I22 + exp(dd$L) * l2^2; I12 = I12 + exp(dd$L) * l1 * l2
	}
	expect_equal(nd$I11, I11, tolerance = 1e-8)
	expect_equal(nd$I22, I22, tolerance = 1e-8)
	expect_lt(max(abs(nd$I12 - I12)), 1e-9)
	expect_true(all(nd$I12^2 <= nd$I11 * nd$I22))
	# its derivatives in the two predictors
	h = 1e-5
	for (nm in c("I11", "I22", "I12")) {
		fa = (lame:::.ae_pair_node_fisher_r(a + h, b, r, s1, s2, FALSE)[[nm]] - lame:::.ae_pair_node_fisher_r(a - h, b, r, s1, s2, FALSE)[[nm]]) / (2 * h)
		fb = (lame:::.ae_pair_node_fisher_r(a, b + h, r, s1, s2, FALSE)[[nm]] - lame:::.ae_pair_node_fisher_r(a, b - h, r, s1, s2, FALSE)[[nm]]) / (2 * h)
		expect_lt(max(abs(nd[[paste0(nm, "a")]] - fa)), 1e-5 * max(abs(fa)))
		expect_lt(max(abs(nd[[paste0(nm, "b")]] - fb)), 1e-5 * max(abs(fb)))
	}
	# the derivatives in the correlation: the score of the log probability and its curvature, the derivatives of the two
	# scores and of the Fisher information
	nr = lame:::.ae_pair_node_fisher_r(a, b, r, s1, s2, TRUE, TRUE)
	up = lame:::.ae_pair_node_fisher_r(a, b, r + h, s1, s2, FALSE)
	dn = lame:::.ae_pair_node_fisher_r(a, b, r - h, s1, s2, FALSE)
	for (nm in c("psi1", "psi2", "I11", "I22", "I12")) {
		fr = (up[[nm]] - dn[[nm]]) / (2 * h)
		expect_lt(max(abs(nr[[switch(nm, psi1 = "psi1r", psi2 = "psi2r", paste0(nm, "r"))]] - fr)), 1e-5 * max(abs(fr)))
	}
	expect_lt(max(abs(nr$s0 - (up$lo - dn$lo) / (2 * h))), 1e-5 * max(abs(nr$s0)))
	expect_lt(max(abs(nr$rc0 - (up$lo - 2 * nd$lo + dn$lo) / h^2)), 1e-4 * max(abs(nr$rc0)))
	# without a correlation the pair is two independent cells: the cell's Fisher weight and no cross weight
	n0 = lame:::.ae_pair_node_fisher_r(a, b, rep(0, n), s1, s2, TRUE)
	expect_equal(n0$I11, dnorm(a)^2 / (pnorm(a) * pnorm(-a)), tolerance = 1e-9)
	expect_equal(n0$I22, dnorm(b)^2 / (pnorm(b) * pnorm(-b)), tolerance = 1e-9)
	expect_equal(n0$I12, rep(0, n), tolerance = 1e-12)
	# predictors far in the tails: finite, with a positive weight
	ex = lame:::.ae_pair_node_fisher_r(
		c(-14, -9, 7, 12, 0),
		c(-14, 4, -9, 12, -12),
		rep(-0.8, 5),
		c(1, 1, -1, 1, -1),
		c(1, -1, -1, 1, 1),
		TRUE
	)
	expect_true(all(vapply(ex, function(v) all(is.finite(v)), TRUE)))
	expect_true(all(ex$I11 > 0 & ex$I22 > 0))
})

test_that("the information of the pairs of cells is the explicit sum over the dyads, and the right-hand side is the pair's gradient and information", {
	su = pair_state(n = 8, Tt = 3, seed = 5)
	ctx = su$ctx; st = su$st
	expect_true(ctx$dcor)
	eta = lame:::.lde_eta(ctx, st)
	wk = lame:::.lde_work_plain(ctx, eta, st$hyp)
	sysm = lame:::.lde_system(ctx, st, wk)
	na = ctx$na; q = ctx$q; m = ctx$mth
	Hth = matrix(0, m, m)
	for (t in seq_len(ctx$Tt)) {
		Kx = matrix(0, q, q); Gx = matrix(0, m, q); rx = numeric(q)
		g = function(i, j) { v = numeric(q); v[i] = 1; v[na + j] = 1; v }
		th = function(i, j) { v = numeric(m); v[ctx$thidx[[t]]] = ctx$Xc[[t]][i + (j - 1) * na, ]; v }
		for (i in seq_len(na)) for (j in seq_len(na)) {
			if (!ctx$obs[i, j, t]) next
			W = wk$w[i, j, t]
			Kx = Kx + W * tcrossprod(g(i, j)); Gx = Gx + W * tcrossprod(th(i, j), g(i, j)); Hth = Hth + W * tcrossprod(th(i, j))
			rx = rx + g(i, j) * W * wk$z[i, j, t]
			cc = wk$cw[i, j, t]
			if (cc != 0) {
				Kx = Kx + cc * tcrossprod(g(i, j), g(j, i)); Gx = Gx + cc * tcrossprod(th(i, j), g(j, i))
				Hth = Hth + cc * tcrossprod(th(i, j), th(j, i))
				rx = rx + g(i, j) * cc * eta[j, i, t]
			}
		}
		expect_equal(sysm$K[[t]], Kx, tolerance = 1e-12)
		expect_equal(sysm$G[[t]], Gx, tolerance = 1e-12)
		expect_equal(unname(sysm$rx[, t]), rx, tolerance = 1e-12)
	}
	expect_equal(sysm$Hth - sysm$Pth, Hth, tolerance = 1e-12)
	# cells whose partner is missing are independent cells, with no cross weight
	t_miss = 1L
	expect_equal(wk$cw[2, 5, t_miss], 0)
	expect_equal(wk$cw[5, 2, t_miss], 0)
	expect_gt(wk$w[5, 2, t_miss], 0)
	expect_gte(min(eigen(sysm$K[[1]], symmetric = TRUE, only.values = TRUE)$values), -1e-12)
})

test_that("the information of the pairs of cells carries the node covariates of both cells as the explicit sum over the dyads", {
	n = 7; Tt = 2
	d = eb_sim(n = n, Tt = Tt, seed = 15, beta = 0.5)
	rn = rownames(d$Y[[1]])
	Xr = lapply(seq_len(Tt), function(t) matrix(sin(seq_len(n) + t), n, 1, dimnames = list(rn, "wr")))
	Xc = lapply(seq_len(Tt), function(t) matrix(cos(seq_len(n) * t), n, 1, dimnames = list(rn, "wc")))
	d$Y[[1]][3, 6] = NA
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		Xr,
		Xc,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	ctx = lame:::.lde_context(prep, FALSE, "probit", coef_prior = "flat", centre = FALSE, dcor = TRUE)
	expect_gt(ctx$pall, 2L)
	st = lame:::.lde_init(ctx)
	set.seed(21)
	st$x[] = rnorm(length(st$x), 0, 0.3)
	st$lev = rowMeans(st$x)
	st$theta[] = seq(-0.5, 0.5, length.out = length(st$theta))
	st$hyp$rho_d = -0.35
	eta = lame:::.lde_eta(ctx, st)
	wk = lame:::.lde_work_plain(ctx, eta, st$hyp)
	sysm = lame:::.lde_system(ctx, st, wk)
	na = ctx$na; q = ctx$q; m = ctx$mth
	Hth = matrix(0, m, m)
	for (t in seq_len(ctx$Tt)) {
		Kx = matrix(0, q, q); Gx = matrix(0, m, q); rx = numeric(q)
		g = function(i, j) { v = numeric(q); v[i] = 1; v[na + j] = 1; v }
		th = function(i, j) { v = numeric(m); v[ctx$thidx[[t]]] = ctx$Xc[[t]][i + (j - 1) * na, ]; v }
		for (i in seq_len(na)) for (j in seq_len(na)) {
			if (!ctx$obs[i, j, t]) next
			W = wk$w[i, j, t]
			Kx = Kx + W * tcrossprod(g(i, j)); Gx = Gx + W * tcrossprod(th(i, j), g(i, j)); Hth = Hth + W * tcrossprod(th(i, j))
			rx = rx + g(i, j) * W * wk$z[i, j, t]
			cc = wk$cw[i, j, t]
			if (cc != 0) {
				Kx = Kx + cc * tcrossprod(g(i, j), g(j, i)); Gx = Gx + cc * tcrossprod(th(i, j), g(j, i))
				Hth = Hth + cc * tcrossprod(th(i, j), th(j, i))
				rx = rx + g(i, j) * cc * eta[j, i, t]
			}
		}
		expect_equal(sysm$K[[t]], Kx, tolerance = 1e-12)
		expect_equal(sysm$G[[t]], Gx, tolerance = 1e-12)
		expect_equal(unname(sysm$rx[, t]), rx, tolerance = 1e-12)
	}
	expect_equal(sysm$Hth - sysm$Pth, Hth, tolerance = 1e-12)
})

test_that("the deviance of the additive predictor alone is the sum over the cells and the pairs", {
	su = pair_state(n = 7, Tt = 3, seed = 6)
	ctx = su$ctx; st = su$st
	eta = lame:::.lde_eta(ctx, st)
	wk = lame:::.lde_work_plain(ctx, eta, st$hyp)
	expect_equal(lame:::.lde_dev_plain(ctx, eta, st$hyp), wk$dev, tolerance = 1e-9)
	# a pair is counted once, with the bivariate probit probability
	pc = ctx$pc
	yy = ctx$y
	yy[c(pc$one, pc$two)] = NA
	by_hand = lame:::.lde_work(yy, eta, "probit")$dev -
		2 * sum(
			lame:::.ae_pair_logphi2(
				pc$s1 * eta[pc$one],
				pc$s2 * eta[pc$two],
				pc$s1 * pc$s2 * rep(st$hyp$rho_d, pc$n_pairs)
			)
		)
	expect_equal(wk$dev, by_hand, tolerance = 1e-9)
	# without a dyadic correlation both are the independent cells'
	ctx0 = ctx; ctx0$dcor = FALSE; ctx0$pc = NULL
	expect_equal(lame:::.lde_dev_plain(ctx0, eta, st$hyp), lame:::.lde_work(ctx$y, eta, "probit")$dev)
})

test_that("the laplace correction of the pairs: the direct term differentiates the log determinant, the shift of the mode is the dense solve", {
	su = pair_state(n = 9, Tt = 3, seed = 9)
	ctx = su$ctx; st = su$st
	ev0 = lame:::.lde_eval(ctx, st, candidate = FALSE)
	# the weights are not capped at this state, so the correction is formed for the uncapped information
	expect_lt(max(abs(ev0$wk$score[ctx$obs] / ev0$wk$w[ctx$obs])), 30)
	logdet = function(ctx, st, delta) {
		s = st
		s$theta[1] = s$theta[1] + delta
		ev = lame:::.lde_eval(ctx, s, candidate = FALSE)
		ev$fac$logdet_xx + ev$fac$logdet_l
	}
	h = 1e-5
	# the direct term: the sum over the cells of T_m is the derivative of the log determinant in a common shift of the predictors
	expect_equal(sum(ev0$lap$Tcell[ctx$obs]), (logdet(ctx, st, h) - logdet(ctx, st, -h)) / (2 * h), tolerance = 1e-7)
	# the same without a dyadic correlation (the independent cells' h' v)
	ctx0 = ctx; ctx0$dcor = FALSE; ctx0$pc = NULL
	st0 = st; st0$hyp$rho_d = NULL
	ev00 = lame:::.lde_eval(ctx0, st0, candidate = FALSE)
	expect_equal(
		sum(ev00$lap$Tcell[ctx0$obs]),
		(logdet(ctx0, st0, h) - logdet(ctx0, st0, -h)) / (2 * h),
		tolerance = 1e-7
	)
	# the shift of the mode: u = M^-1 E'T with the dense information M of the paths and the levels, and the partner's cross weight
	q = ctx$q; na = ctx$na; Tt = ctx$Tt
	il = seq_len(q)
	sm = ev0$sysm
	Hd = lame:::.lde_dense_matrix(
		sm$A,
		sm$b,
		lapply(sm$C, function(Ct) Ct[, il, drop = FALSE]),
		sm$Hzz[il, il, drop = FALSE]
	)
	rhs = c(
		unlist(lapply(seq_len(Tt), function(t) c(rowSums(ev0$lap$Tcell[, , t]), colSums(ev0$lap$Tcell[, , t])))),
		numeric(q)
	)
	u = solve(Hd, rhs)
	kappa = array(0, dim(ev0$wk$w))
	for (t in seq_len(Tt)) {
		ut = u[(t - 1L) * q + seq_len(q)]
		uc = outer(ut[seq_len(na)], ut[na + seq_len(na)], "+")
		kappa[, , t] = 0.5 * (ev0$lap$Tcell[, , t] - ev0$wk$w[, , t] * uc - ev0$wk$cw[, , t] * t(uc))
	}
	kappa[!ctx$obs] = 0
	expect_equal(ev0$lap$kappa, kappa, tolerance = 1e-9)
})

test_that("the score of the correlation is the slope of the engine's profile objective, and the log determinant and the shift of the mode matter", {
	skip_on_cran()
	d = eb_sim(n = 18, Tt = 4, seed = 52, dcor = 0.7)
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	ctx = lame:::.lde_context(prep, FALSE, "probit", dcor = TRUE)
	fit = lame:::.lde_fit(ctx, max_iter = 300L, tol = 1e-8)
	ctx = fit$ctx
	rho0 = fit$st$hyp$rho_d
	ctxh = ctx
	ctxh$rho_hold = TRUE
	prior = function(r) -2 * lame:::.lame_rho_prior_log(r)
	# the objective at the converged fit with the correlation held (the profile), without the prior of the correlation
	at = function(r) {
		st = fit$st
		st$hyp$rho_d = r
		rr = lame:::.lde_fit_core(ctxh, st, max_iter = 500L, tol = 1e-9)
		expect_true(rr$converged)
		list(obj = rr$ev$objective - prior(r), st = rr$st)
	}
	h = 0.05
	up = at(rho0 + h); dn = at(rho0 - h); mid = at(rho0)
	fd = -(up$obj - dn$obj) / (2 * h) / 2
	ev = lame:::.lde_eval(ctx, mid$st, candidate = FALSE)
	pr = ev$wk$pair
	expect_lt(abs(ev$lap$rho_score$rs - fd), 0.25 * abs(fd))
	# the plug-in score (the pairs at the mode alone) is several units away: the correction is the larger part
	expect_gt(abs(sum(pr$s0) - fd), 5 * abs(ev$lap$rho_score$rs - fd))
	# the curvature of the expected log probability is the node's, summed
	expect_equal(ev$lap$rho_score$rc, sum(pr$rc0))
	expect_lt(ev$lap$rho_score$rc, 0)
	# a held correlation takes no step, and the record of it is the held value
	ev_h = lame:::.lde_eval(ctxh, mid$st, candidate = FALSE)
	expect_null(ev_h$rho_step)
	expect_equal(ev_h$em_self$rho_d, mid$st$hyp$rho_d)
})

test_that("the cross information of the pairs enters the factor blocks as the explicit sum over the dyads", {
	set.seed(3)
	d = eb_sim(n = 7, Tt = 2, seed = 4)
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	ctx = lame:::.lde_context(prep, FALSE, "probit", R = 2L, coef_prior = "flat", centre = FALSE, dcor = TRUE)
	na = ctx$na; R = ctx$R; q = ctx$q
	W = matrix(runif(na * na, 0.1, 0.5), na, na); diag(W) = 0
	Ct = matrix(runif(na * na, -0.1, 0.1), na, na); Ct = Ct + t(Ct); diag(Ct) = 0
	U = matrix(rnorm(na * R), na, R); V = matrix(rnorm(na * R), na, R)
	Sig = array(0, c(R, R, na)); Om = array(0, c(R, R, na))
	wb = lame:::.lde_fac_wblocks(ctx, W, U, V, Sig, Om, Ct)
	ix = lame:::.lde_widx(ctx)
	nf = 2L * na * R
	# the gradient of the predictor of the cell (i, j) in the paths and in the factors (u_i: v_j, v_j: u_i)
	gp = function(i, j) { v = numeric(q); v[i] = 1; v[na + j] = 1; v }
	gf = function(i, j) {
		v = numeric(nf)
		for (r in seq_len(R)) { v[ix$u[[r]][i]] = V[j, r]; v[ix$v[[r]][j]] = U[i, r] }
		v
	}
	Kxw = matrix(0, q, nf); Kww = matrix(0, nf, nf)
	for (i in seq_len(na)) for (j in seq_len(na)) if (i != j) {
		Kxw = Kxw + W[i, j] * tcrossprod(gp(i, j), gf(i, j)) + Ct[i, j] * tcrossprod(gp(i, j), gf(j, i))
		Kww = Kww + W[i, j] * tcrossprod(gf(i, j)) + Ct[i, j] * tcrossprod(gf(i, j), gf(j, i))
	}
	expect_equal(wb$Kxw, Kxw, tolerance = 1e-12)
	expect_equal(wb$Kww, Kww, tolerance = 1e-12)
	# the block between a coefficient and the factors: sum c_ij x_ij,k gf(j, i), the way the joint covariance forms it
	Xk = matrix(rnorm(na * na), na, na)
	cross = numeric(nf)
	for (i in seq_len(na)) for (j in seq_len(na)) if (i != j) cross = cross + Xk[i, j] * Ct[i, j] * gf(j, i)
	expect_equal(lame:::.lde_wsum(ctx, t(Xk * Ct), U, V), cross, tolerance = 1e-12)
	# the cavity weights of the pairs are the static estimator's, and the independent cells' without a correlation
	ev = lame:::.lde_eval(ctx, lame:::.lde_init(ctx), candidate = FALSE)
	hyp = ev$st$hyp; hyp$rho_d = 0.4
	V3 = array(0.3, dim(ctx$y))
	cw1 = lame:::.lde_cavity_w(ctx, ev$eta, V3, ev$lap, hyp)
	cw2 = lame:::.ae_pair_info_w(list(pair = ctx$pc, y = ctx$y, link = "probit"), ev$eta, V3, ev$lap$cvpair, 0.4)
	expect_equal(cw1$w, cw2$w); expect_equal(cw1$cw, cw2$cw)
	ctx0 = ctx; ctx0$dcor = FALSE
	expect_equal(
		lame:::.lde_cavity_w(ctx0, ev$eta, V3, ev$lap, hyp)$w,
		lame:::.ae_mf_info_w(ctx$y, ev$eta, V3, "probit")
	)
	expect_null(lame:::.lde_cavity_w(ctx0, ev$eta, V3, ev$lap, hyp)$cw)
})

test_that("a directed fit models the dyadic correlation by default and recovers it; dcor = FALSE is the independent-cell fit", {
	skip_on_cran()
	d = eb_sim(n = 24, Tt = 4, seed = 41, dcor = 0.6)
	f1 = eb_fit(d)
	f0 = eb_fit(d, dcor = FALSE)
	expect_true(f1$converged && f0$converged)
	expect_false(is.null(f1$eb$pair))
	expect_null(f0$eb$pair)
	pr = f1$eb$pair
	expect_equal(unname(f1$VC["rho"]), pr$rho)
	expect_true(f1$rho_estimation$identified)
	expect_equal(
		pr$n_pairs,
		sum(
			vapply(
				seq_along(d$Y),
				function(t) { y = d$Y[[t]]; sum(is.finite(y[upper.tri(y)]) & is.finite(t(y)[upper.tri(y)])) },
				1
			)
		)
	)
	expect_true(pr$rho > 0 && abs(pr$rho - 0.6) < 4 * pr$sd)
	expect_gt(pr$sd, 0.03)
	expect_lt(pr$sd, 0.2)
	expect_null(pr$prior)
	# the independent-cell fit reports the correlation of the pairs of its cells
	expect_false(isTRUE(all.equal(unname(f0$VC["rho"]), pr$rho, tolerance = 1e-6)))
	# the score is stationary at the fitted correlation (the posterior mean of the quadratic model is the value itself)
	expect_lt(abs(pr$score / pr$curvature), 0.1)
	# a deviance that counts a pair once is below the independent cells' at the same predictor
	expect_lt(f1$deviance, f0$deviance)
	# coefficients agree: the correlation changes them by a fraction of their standard errors
	expect_lt(max(abs(coef(f1) - coef(f0)) / f0$eb$se), 1)
	# print mentions it
	expect_true(any(grepl("Dyadic correlation", testthat::capture_messages(lame:::.lde_print_eb(f1)))))
})

test_that("without reciprocity the fit with the dyadic correlation is the independent-cell fit to a fraction of a standard error", {
	skip_on_cran()
	d = eb_sim(n = 24, Tt = 4, seed = 42)
	f1 = eb_fit(d)
	f0 = eb_fit(d, dcor = FALSE)
	expect_lt(abs(f1$eb$pair$rho), 4 * f1$eb$pair$sd)
	expect_lt(max(abs(coef(f1) - coef(f0)) / f0$eb$se), 0.3)
	expect_lt(max(abs(f1$eb$se / f0$eb$se - 1)), 0.05)
})

test_that("the pair problem at a correlation held at zero is the independent-cell fit, on a panel with events and on one without", {
	skip_on_cran()
	d = eb_sim(n = 14, Tt = 3, seed = 5, dcor = 0.5)
	for (zero in c(FALSE, TRUE)) {
		if (zero) d$Y = lapply(d$Y, function(y) { y[!is.na(y)] = 0; y })
		prep = lame:::.ae_prepare(
			d$Y,
			d$X,
			NULL,
			NULL,
			"unipartite",
			longitudinal = TRUE,
			node_time_policy = "quiet",
			symmetric = FALSE
		)
		ctx0 = lame:::.lde_context(prep, FALSE, "probit", dcor = FALSE)
		ctx1 = lame:::.lde_context(prep, FALSE, "probit", dcor = TRUE)
		ctx1$rho_hold = TRUE
		f0 = lame:::.lde_fit(ctx0, max_iter = 300L, tol = 5e-5)
		f1 = lame:::.lde_fit(ctx1, max_iter = 300L, tol = 5e-5)
		expect_true(f0$converged && f1$converged)
		expect_equal(f1$iterations, f0$iterations)
		expect_equal(f1$st$theta, f0$st$theta, tolerance = 1e-8)
		expect_equal(f1$st$x, f0$st$x, tolerance = 1e-8)
		expect_equal(f1$ev$objective, f0$ev$objective, tolerance = 1e-6)
		expect_equal(f1$st$hyp$rho_d, 0)
	}
	# an uninformative panel (no events: nothing to say about the correlation) returns it to the prior's mean and converges
	f = suppressWarnings(
		lame:::.lde_fit(lame:::.lde_context(prep, FALSE, "probit", dcor = TRUE), max_iter = 300L, tol = 5e-5)
	)
	expect_true(f$converged)
	expect_lt(abs(f$st$hyp$rho_d), 1e-3)
	# and so does a panel in which every dyad has exactly one tie (the log likelihood of the correlation is not concave)
	d$Y = lapply(d$Y, function(y) { y[upper.tri(y)] = 1; y[lower.tri(y)] = 0; y })
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	f = suppressWarnings(
		lame:::.lde_fit(lame:::.lde_context(prep, FALSE, "probit", dcor = TRUE), max_iter = 300L, tol = 5e-5)
	)
	expect_true(f$converged)
})

test_that("the fixed point of the dyadic correlation does not depend on where the iteration starts it", {
	skip_on_cran()
	d = eb_sim(n = 20, Tt = 4, seed = 47, dcor = 0.5)
	prep = lame:::.ae_prepare(
		d$Y,
		d$X,
		NULL,
		NULL,
		"unipartite",
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = FALSE
	)
	ctx = lame:::.lde_context(prep, FALSE, "probit", dcor = TRUE)
	fits = lapply(c(0, 0.9, -0.8), function(r0) {
		st0 = lame:::.lde_init(ctx)
		st0$hyp$rho_d = r0
		lame:::.lde_fit(ctx, st0, max_iter = 600L, tol = 5e-5)
	})
	expect_true(all(vapply(fits, function(f) f$converged, TRUE)))
	rho = vapply(fits, function(f) f$st$hyp$rho_d, 1)
	expect_lt(diff(range(rho)), 2e-3)
	th = vapply(fits, function(f) f$st$theta, numeric(ctx$mth))
	expect_lt(max(apply(th, 1L, function(v) diff(range(v)))), 5e-3)
})

test_that("a refit from the record of a fit with a dyadic correlation is at its fixed point, and the problem is part of the record", {
	skip_on_cran()
	d = eb_sim(n = 18, Tt = 3, seed = 43, dcor = 0.5)
	f1 = eb_fit(d)
	f2 = suppressWarnings(suppressMessages(ame_als_refit(f1)))
	expect_equal(f2$iterations, 0L)
	expect_equal(unname(f2$VC["rho"]), unname(f1$VC["rho"]))
	expect_equal(coef(f2), coef(f1))
	# a record of the independent-cell problem is not the record of this one
	f0 = eb_fit(d, dcor = FALSE)
	f3 = suppressWarnings(
		suppressMessages(
			lame:::lame_dynamic_als(d$Y, d$X, family = "binary", dynamic_ab = TRUE, verbose = FALSE, start = f0)
		)
	)
	expect_gt(f3$iterations, 0L)
	expect_false(is.null(f3$eb$pair))
})

test_that("the dyadic correlation of a directed fit with missing cells and a changing composition, and the networks and links that have none", {
	skip_on_cran()
	d = eb_sim(n = 14, Tt = 4, seed = 44, dcor = 0.5)
	Y2 = d$Y
	for (t in 1:4) { gone = ((t * 3):(t * 3 + 1)) %% 14 + 1; Y2[[t]][gone, ] = NA; Y2[[t]][, gone] = NA }
	f = eb_fit(list(Y = Y2, X = d$X))
	expect_true(f$converged)
	expect_gt(f$eb$pair$n_pairs, 0L)
	Y3 = lapply(
		1:4,
		function(t) { keep = setdiff(
			rownames(d$Y[[t]]),
			rownames(d$Y[[t]])[((t - 1) * 2 + 1):((t - 1) * 2 + 2)]
		); d$Y[[t]][keep, keep] }
	)
	X3 = lapply(
		1:4,
		function(t) { keep = setdiff(
			rownames(d$Y[[t]]),
			rownames(d$Y[[t]])[((t - 1) * 2 + 1):((t - 1) * 2 + 2)]
		); d$X[[t]][keep, keep, , drop = FALSE] }
	)
	names(Y3) = names(X3) = names(d$Y)
	f3 = eb_fit(list(Y = Y3, X = X3))
	expect_true(f3$converged)
	# the logit link has no such correlation: an explicit request is an error, the default is the independent cells
	expect_error(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			dynamic_ab = TRUE,
			link = "logit",
			dcor = TRUE,
			verbose = FALSE
		),
		"probit"
	)
	fl = eb_fit(d, link = "logit")
	expect_null(fl$eb$pair)
	# undirected and bipartite networks have no pair of cells: the request is accepted and the record is absent
	ds = eb_sim(n = 10, Tt = 3, mode = "symmetric", seed = 45)
	fs = suppressWarnings(
		suppressMessages(
			lame(
				ds$Y,
				Xdyad = ds$X,
				family = "binary",
				method = "als",
				symmetric = TRUE,
				dynamic_ab = TRUE,
				dcor = TRUE,
				verbose = FALSE
			)
		)
	)
	expect_null(fs$eb$pair)
	db = eb_sim(n = 10, nb = 7, Tt = 3, mode = "bipartite", seed = 46)
	fb = suppressWarnings(
		suppressMessages(
			lame(
				db$Y,
				Xdyad = db$X,
				family = "binary",
				method = "als",
				mode = "bipartite",
				dynamic_ab = TRUE,
				dcor = TRUE,
				verbose = FALSE
			)
		)
	)
	expect_null(fb$eb$pair)
	# the front door: the default of lame() for a directed network, and the switch
	fd = suppressWarnings(
		suppressMessages(lame(d$Y, Xdyad = d$X, family = "binary", method = "als", dynamic_ab = TRUE, verbose = FALSE))
	)
	expect_false(is.null(fd$eb$pair))
	fd0 = suppressWarnings(
		suppressMessages(
			lame(d$Y, Xdyad = d$X, family = "binary", method = "als", dynamic_ab = TRUE, dcor = FALSE, verbose = FALSE)
		)
	)
	expect_null(fd0$eb$pair)
})

test_that("the pair quadrature of the engine: the node rule per pair, and the compiled kernel equals the R reference", {
	set.seed(4)
	n = 100
	m1 = rnorm(n, -1.5, 1.4); m2 = rnorm(n, -1.5, 1.4)
	v1 = c(rep(0, 10), runif(n - 10, 0.01, 2)); v2 = pmax(v1 * runif(n, 0.5, 1.5), 0)
	s1 = sample(c(-1, 1), n, TRUE, prob = c(0.8, 0.2)); s2 = sample(c(-1, 1), n, TRUE, prob = c(0.8, 0.2))
	# one node without noise, then the smallest rung with nine nodes a standard deviation of the larger noise
	kk = lame:::.lde_pair_rule(v1, v2)
	expect_true(all(kk[v1 == 0 & v2 == 0] == 1L))
	expect_true(all(kk[pmax(v1, v2) > 1e-12] >= pmin(61, 9 * sqrt(pmax(v1, v2)[pmax(v1, v2) > 1e-12]))))
	expect_true(all(kk %in% c(1L, 5L, 9L, 15L, 21L, 31L, 41L, 61L)))
	for (rho in c(-0.45, 0.3)) {
		k = lame:::.lde_pair_quad(m1, m2, v1, v2, s1, s2, rho, TRUE, rho_deriv = TRUE)
		r = withr::with_options(
			list(lame.als_engine = "r"),
			lame:::.lde_pair_quad(m1, m2, v1, v2, s1, s2, rho, TRUE, rho_deriv = TRUE)
		)
		expect_setequal(names(k), names(r))
		expect_identical(as.integer(k$K), as.integer(r$K))
		for (f in setdiff(names(r), "K")) expect_lt(
			max(abs(k[[f]] - r[[f]])),
			1e-9 * max(1, max(abs(r[[f]]))),
			label = f
		)
	}
	k = lame:::.lde_pair_quad(m1, m2, v1, v2, s1, s2, 0.2, FALSE)
	expect_false("I11a" %in% names(k))
	# the number of threads changes nothing
	q1 = withr::with_options(
		list(lame.als_threads = 1L),
		lame:::.lde_pair_quad(m1, m2, v1, v2, s1, s2, 0.3, TRUE, rho_deriv = TRUE)
	)
	q3 = withr::with_options(
		list(lame.als_threads = 3L),
		lame:::.lde_pair_quad(m1, m2, v1, v2, s1, s2, 0.3, TRUE, rho_deriv = TRUE)
	)
	expect_identical(q3, q1)
})

# a call of the pair quadrature with wide noise, the grid path forced (the planner would take it where the noise is wide)
pair_wide_call = function(seed = 21, n = 70, sd_range = c(1, 5), m_sd = 2.5) {
	set.seed(seed)
	list(
		m1 = rnorm(n, -1, m_sd),
		m2 = rnorm(n, -1, m_sd),
		v1 = runif(n, sd_range[1], sd_range[2])^2,
		v2 = runif(n, sd_range[1], sd_range[2])^2,
		s1 = sample(c(-1, 1), n, TRUE),
		s2 = sample(c(-1, 1), n, TRUE)
	)
}
pair_wide_quad = function(cl, rho, ..., tier = 2L, deriv = TRUE, rho_deriv = TRUE) {
	withr::with_options(
		c(list(lame.pair_grid = "force"), list(...)),
		lame:::.lde_pair_quad(cl$m1, cl$m2, cl$v1, cl$v2, cl$s1, cl$s2, rho, deriv, rho_deriv, tier = tier)
	)
}
pair_rel_diff = function(x, y) {
	fs = setdiff(intersect(names(x), names(y)), "K")
	max(vapply(fs, function(f) max(abs(x[[f]] - y[[f]]) / pmax(1, abs(y[[f]]))), 0))
}

test_that("the table of the pair quadrature on the shared grid equals the R reference, the nodes are right at extreme predictors, and the number of threads changes nothing", {
	# the compiled table and contraction against the R reference on a small grid
	rho = 0.35; cc = 3; du = 0.12; klo = -22L; khi = 22L
	tb = lame:::lde_pair_table_cpp(rho, cc, du, klo, khi, 15L, 1L)
	tr = lame:::.lde_pair_table_r(rho, cc, du, klo, khi, 15L)
	expect_lt(max(abs(tb$table - tr$table) / pmax(1, abs(tr$table))), 1e-7)
	set.seed(5)
	n = 40
	m1 = rnorm(n, -1, 3); m2 = rnorm(n, -1, 3); sd1 = runif(n, 0.6, 5); sd2 = runif(n, 0.6, 5); v = sample(0:3, n, TRUE)
	for (dr in list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE))) {
		a = lame:::lde_pair_apply_cpp(tb, m1, m2, sd1, sd2, v, 7.5, dr[1], dr[2], 1L)
		b = lame:::.lde_pair_apply_r(tr, m1, m2, sd1, sd2, v, 7.5, dr[1], dr[2])
		expect_setequal(names(a), names(b))
		expect_identical(a$ok, b$ok)
		expect_identical(as.integer(a$K), as.integer(b$K))
		expect_lt(pair_rel_diff(a, b), 1e-9)
	}
	expect_identical(
		lame:::lde_pair_apply_cpp(tb, m1, m2, sd1, sd2, v, 7.5, TRUE, TRUE, 3L),
		lame:::lde_pair_apply_cpp(tb, m1, m2, sd1, sd2, v, 7.5, TRUE, TRUE, 1L)
	)
	# a window outside the grid, and a predictor without noise, are not computed
	oo = lame:::lde_pair_apply_cpp(
		tb,
		c(0, 500, 0),
		c(0, 0, 0),
		c(1, 1, 1),
		c(1, 1, 0),
		c(0L, 0L, 0L),
		7.5,
		FALSE,
		FALSE,
		1L
	)
	expect_identical(oo$ok, c(1L, 0L, 0L))
	expect_identical(oo$ell[2:3], c(0, 0))
	# the nodes at extreme predictors (a predictor of -40 has a probability below the smallest double: the logarithms come from the
	# log-domain integral, the subtraction of the margins must not be taken where both terms have underflowed)
	g = c(-70, -45, -39, -38, -20, -5, -0.5, 0, 0.5, 5, 20, 38, 39, 45, 70)
	A = rep(g, each = length(g)); B = rep(g, times = length(g))
	for (rho in c(-0.7, 0.3)) for (v in 0:3) {
		cp = lame:::lde_pair_direct_cpp(
			A,
			B,
			0 * A,
			0 * A,
			rep(v, length(A)),
			0L * seq_along(A),
			0L * seq_along(A),
			1L + 0L * seq_along(A),
			1L + 0L * seq_along(A),
			rho,
			3,
			0.05,
			7.5,
			TRUE,
			TRUE,
			1L
		)
		rr = lame:::.lde_pair_node_r(A, B, rho, bitwShiftL(1L, v))
		obs = rr[, 12 + 7 * v + 1:7]
		cmp = cbind(cp$ell, cp$psi1, cp$psi2, cp$s0, cp$rc0, cp$psi1r, cp$psi2r)
		expect_true(all(is.finite(cmp)))
		expect_lt(max(abs(cmp - obs) / pmax(1, abs(obs))), 1e-7)
		expect_lt(
			max(
				abs(
					cbind(cp$I11, cp$I22, cp$I12, cp$I11a, cp$I22b, cp$I11r) - rr[, c(1, 2, 3, 4, 7, 10)]
				) / pmax(1, abs(rr[, c(1, 2, 3, 4, 7, 10)]))
			),
			1e-9
		)
	}
})

test_that("lde_pair_table_blend_cpp is the linear combination of tables, whatever the number of threads", {
	set.seed(3)
	a = rnorm(10000); b = rnorm(10000); cc = rnorm(10000)
	w = c(0.5, -0.25, 0.75)
	r1 = lame:::lde_pair_table_blend_cpp(list(a, b, cc), w, 1L)
	expect_equal(r1, 0.5 * a - 0.25 * b + 0.75 * cc, tolerance = 1e-14)
	expect_identical(lame:::lde_pair_table_blend_cpp(list(a, b, cc), w, 3L), r1)
	expect_identical(lame:::lde_pair_table_blend_cpp(list(a), 1), a)
	expect_error(lame:::lde_pair_table_blend_cpp(list(a, b[-1L]), c(1, 0)), "sizes")
	expect_error(lame:::lde_pair_table_blend_cpp(list(a, b), 1), "lengths")
})

test_that("the family of tables in the correlation: members are shared inside the window, the interpolation reproduces the exact table, a window is rebuilt outside it, a fit leaves no cache", {
	skip_on_cran()
	cl = pair_wide_call()
	lame:::.lde_pair_reset()
	r0 = pair_wide_quad(cl, 0.3)
	f0 = lame:::.lde_pair_env$fam
	expect_gt(f0$n, 1L)
	expect_identical(f0$tier, 2L)
	expect_equal(f0$rc, 0.3)
	expect_lte(f0$h, 0.015)
	expect_equal(sum(abs(f0$bw)), f0$n - 1)
	# a nearby correlation: the same members, the table is their combination; the answer is the exact table's
	rho1 = 0.3 + 0.6 * f0$h
	r1 = pair_wide_quad(cl, rho1)
	f1 = lame:::.lde_pair_env$fam
	expect_identical(f1$tabs, f0$tabs)
	expect_equal(f1$cur$rho, rho1)
	lame:::.lde_pair_reset()
	e1 = pair_wide_quad(cl, rho1, lame.pair_family = FALSE)
	expect_equal(lame:::.lde_pair_env$fam$n, 1L)
	expect_lt(pair_rel_diff(r1, e1), 2e-9)
	# the second call at the same correlation does not combine again (the table of the evaluation is kept)
	lame:::.lde_pair_reset()
	pair_wide_quad(cl, 0.3)
	tab_a = lame:::.lde_pair_env$fam$cur$tab
	pair_wide_quad(cl, 0.3, deriv = FALSE, rho_deriv = FALSE)
	expect_identical(lame:::.lde_pair_env$fam$cur$tab, tab_a)
	# the correlation outside the window: a new family about it; the coarse tier makes a smaller family that the fine tier does not use
	rho2 = 0.3 + 3 * f0$h
	pair_wide_quad(cl, rho2)
	expect_equal(lame:::.lde_pair_env$fam$rc, rho2)
	lame:::.lde_pair_reset()
	pair_wide_quad(cl, 0.3, tier = 1L)
	f_t1 = lame:::.lde_pair_env$fam
	expect_identical(f_t1$tier, 1L)
	expect_lt(f_t1$n, f0$n)
	pair_wide_quad(cl, 0.3, tier = 2L)
	expect_identical(lame:::.lde_pair_env$fam$tier, 2L)
	expect_gt(lame:::.lde_pair_env$fam$n, f_t1$n)
	lame:::.lde_pair_reset()
	expect_null(lame:::.lde_pair_env$fam)
	expect_null(lame:::.lde_pair_env$hint)
	# a grid that no table fits in memory is refused by name, before anything is allocated
	big = list(N = 3000L, cc = 3, du = 0.0125, klo = -1500L, khi = 1499L)
	expect_error(lame:::.lde_pair_family(big, 0.3, 15L, 1L, 2L), "too wide")
})

test_that("a fit with the family of tables is the fit with one exact table a correlation, leaves no cache, and is the same fit twice in a session", {
	skip_on_cran()
	d = eb_sim(n = 10, Tt = 3, seed = 3, dcor = 0.5)
	f = eb_fit(d, R = 1, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_null(lame:::.lde_pair_env$fam)
	expect_null(lame:::.lde_pair_env$hint)
	fe = withr::with_options(list(lame.pair_family = FALSE), eb_fit(d, R = 1, dynamic_uv = TRUE))
	expect_lt(max(abs(f$coefficients - fe$coefficients)), 1e-6)
	expect_lt(abs(f$VC["rho"] - fe$VC["rho"]), 1e-6)
	expect_lt(abs(log(f$eb$uv_scale / fe$eb$uv_scale)), 1e-4)
	f2 = eb_fit(d, R = 1, dynamic_uv = TRUE)
	expect_identical(f2$coefficients, f$coefficients)
	expect_identical(f2$iterations, f$iterations)
})

test_that("the pair quadrature is accurate across the noise range and the family does not change it", {
	skip_on_cran()
	# against a finer grid (half the spacing) and against 199-point Gauss-Hermite where that rule is accurate (standard deviations up to 2)
	for (sdv in c(0.3, 1, 3, 8, 15)) {
		set.seed(round(10 * sdv))
		n = 30
		cl = list(
			m1 = rnorm(n, -2, 2),
			m2 = rnorm(n, -2, 2),
			v1 = rep(sdv^2, n),
			v2 = rep((0.5 + runif(1)) * sdv, n)^2,
			s1 = sample(c(-1, 1), n, TRUE),
			s2 = sample(c(-1, 1), n, TRUE)
		)
		lame:::.lde_pair_reset()
		q = withr::with_options(
			list(lame.pair_grid = "force"),
			lame:::.lde_pair_quad(cl$m1, cl$m2, cl$v1, cl$v2, cl$s1, cl$s2, 0.45, TRUE, TRUE)
		)
		lame:::.lde_pair_reset()
		qf = withr::with_options(
			list(lame.pair_grid = "force", lame.pair_du_scale = 0.5, lame.pair_family = FALSE),
			lame:::.lde_pair_quad(cl$m1, cl$m2, cl$v1, cl$v2, cl$s1, cl$s2, 0.45, TRUE, TRUE)
		)
		expect_lt(pair_rel_diff(q, qf), 1e-8, label = paste("fine grid, sd", sdv))
		if (sdv <= 2) {
			gh = lame:::.ae_pair_quad_fisher(cl$m1, cl$m2, cl$v1, cl$v2, 0, cl$s1, cl$s2, 0.45, 199L, TRUE, TRUE)
			expect_lt(pair_rel_diff(q, gh), 1e-8, label = paste("gauss-hermite 199, sd", sdv))
		}
	}
	lame:::.lde_pair_reset()
})

test_that("the planner: no wide noise is the shared rule, the grid can be switched off, and the coarse tier is a coarser rule that agrees to 1e-4", {
	skip_on_cran()
	cl = pair_wide_call(sd_range = c(0, 0.05))
	lame:::.lde_pair_reset()
	expect_null(lame:::.lde_pair_plan(cl$m1, cl$m2, sqrt(cl$v1), sqrt(cl$v2), 2L * (cl$s1 < 0) + (cl$s2 < 0), 0.3, 19L))
	k = lame:::.lde_pair_quad(cl$m1, cl$m2, cl$v1, cl$v2, cl$s1, cl$s2, 0.3, TRUE, TRUE)
	sh = lame:::.ae_pair_quad_fisher(
		cl$m1,
		cl$m2,
		cl$v1,
		cl$v2,
		0,
		cl$s1,
		cl$s2,
		0.3,
		lame:::.lde_pair_rule(cl$v1, cl$v2),
		TRUE,
		TRUE
	)
	for (f in setdiff(names(sh), "K")) expect_identical(k[[f]], sh[[f]], label = f)
	cl = pair_wide_call()
	off = withr::with_options(
		list(lame.pair_grid = FALSE),
		lame:::.lde_pair_quad(cl$m1, cl$m2, cl$v1, cl$v2, cl$s1, cl$s2, 0.3, TRUE, TRUE)
	)
	sh = lame:::.ae_pair_quad_fisher(
		cl$m1,
		cl$m2,
		cl$v1,
		cl$v2,
		0,
		cl$s1,
		cl$s2,
		0.3,
		lame:::.lde_pair_rule(cl$v1, cl$v2),
		TRUE,
		TRUE
	)
	for (f in setdiff(names(sh), "K")) expect_identical(off[[f]], sh[[f]], label = f)
	vr = 2L * (cl$s1 < 0) + (cl$s2 < 0)
	lame:::.lde_pair_reset()
	pl = lame:::.lde_pair_plan(cl$m1, cl$m2, sqrt(cl$v1), sqrt(cl$v2), vr, 0.3, 19L, force = TRUE)
	expect_true(any(pl$tab))
	expect_true(all(pl$tab | seq_along(vr) %in% pl$direct))
	expect_identical(pl$vmask, 15L)
	fine = pair_wide_quad(cl, 0.3, tier = 2L)
	lame:::.lde_pair_reset()
	coarse = pair_wide_quad(cl, 0.3, tier = 1L)
	expect_lt(pair_rel_diff(coarse, fine), 2e-4)
	# the coarse tier plans a coarser grid and a narrower window than the fine tier
	lame:::.lde_pair_reset()
	p1 = lame:::.lde_pair_plan(cl$m1, cl$m2, sqrt(cl$v1), sqrt(cl$v2), vr, 0.3, 19L, tier = 1L, force = TRUE)
	lame:::.lde_pair_reset()
	p2 = lame:::.lde_pair_plan(cl$m1, cl$m2, sqrt(cl$v1), sqrt(cl$v2), vr, 0.3, 19L, tier = 2L, force = TRUE)
	expect_gte(p1$du, p2$du)
	expect_lt(p1$cutoff, p2$cutoff)
	lame:::.lde_pair_reset()
})

test_that("the quadratic search of the scale finds a parabola from three evaluations and its vertex beyond them, and reports no vertex as an edge", {
	calls = 0
	mk = function(S0, curv = 5) function(S) { calls <<- calls + 1; 3 + curv * (log(S) - log(S0))^2 }
	# the minimiser between the points: three evaluations
	r = lame:::.lde_quad_search_log(mk(2.2), 2, w = 0.12)
	expect_false(r$edge)
	expect_equal(calls, 3)
	expect_equal(sort(r$S), sort(2 * exp(c(-0.12, 0, 0.12))))
	# the minimiser beyond them: the vertex is evaluated, then the three nearest the best one close the search
	calls = 0
	r = lame:::.lde_quad_search_log(mk(3.5), 2, w = 0.12)
	expect_false(r$edge)
	expect_equal(calls, 4)
	expect_equal(max(r$S), 3.5, tolerance = 1e-8)
	# no convexity: no vertex (the search steps outward and gives up)
	r = lame:::.lde_quad_search_log(function(S) -log(S), 2, w = 0.12, max_extra = 2L)
	expect_true(r$edge)
	# a minimiser beyond the span (here a factor 4 from the centre): the search is held to the span and reports the edge
	r = lame:::.lde_quad_search_log(mk(8), 2, w = 0.12, span = 0.7)
	expect_true(r$edge)
	expect_lte(max(r$S), 2 * exp(0.7) * (1 + 1e-9))
})

test_that("the profile of the pair likelihood by the quadratic search agrees with the golden search and the full search", {
	skip_on_cran()
	# a panel with moving positions of standard deviation 1.2 (a scale of 2.5) and errors correlated within the dyad
	set.seed(8)
	n = 16; Tt = 3
	U = array(rnorm(n * Tt, 0, 1.2), c(n, 1, Tt)); V = array(rnorm(n * Tt, 0, 1.2), c(n, 1, Tt))
	rn = sprintf("a%02d", seq_len(n))
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		xt = array(rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x1"))
		eta = -1 + 0.5 * xt[, , 1] + tcrossprod(U[, , t], V[, , t])
		err = matrix(rnorm(n * n), n, n)
		lo = lower.tri(err)
		err[lo] = 0.5 * t(err)[lo] + sqrt(1 - 0.25) * matrix(rnorm(n * n), n, n)[lo]
		y = 1 * (eta + err > 0)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y; X[[t]] = xt
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	d = list(Y = Y, X = X)
	fq = eb_fit(d, R = 1, dynamic_uv = TRUE)
	fg = withr::with_options(list(lame.dcor_golden_search = TRUE), eb_fit(d, R = 1, dynamic_uv = TRUE))
	ff = withr::with_options(list(lame.dcor_full_search = TRUE), eb_fit(d, R = 1, dynamic_uv = TRUE))
	for (f in list(fq, fg, ff)) expect_true(f$converged)
	expect_lt(abs(log(fq$eb$uv_scale / ff$eb$uv_scale)), 0.03)
	expect_lt(abs(log(fg$eb$uv_scale / ff$eb$uv_scale)), 0.03)
	expect_lt(max(abs(fq$coefficients - ff$coefficients)), 5e-3)
	expect_lt(abs(fq$VC["rho"] - ff$VC["rho"]), 5e-3)
	# the search of the pair likelihood is the quadratic one: a handful of evaluations about the independent cells' scale
	expect_lte(nrow(fq$eb$uv_search), 7L)
	expect_gte(nrow(fq$eb$uv_search), 3L)
	# a quadratic search that finds no vertex hands the panel to the search over [0.08, 4], before any polish is paid for: the same fit
	testthat::local_mocked_bindings(.lde_quad_search_log = function(at, centre, ...) {
		S = centre * exp(c(0, -0.12, 0.12))
		list(S = S, value = vapply(S, at, 0), edge = TRUE, range = range(S))
	}, .package = "lame")
	fm = eb_fit(d, R = 1, dynamic_uv = TRUE)
	expect_true(fm$converged)
	expect_identical(fm$eb$uv_search, ff$eb$uv_search)
	expect_identical(fm$coefficients, ff$coefficients)
})

test_that("static factors with the dyadic correlation: the fit converges, recovers the correlation, and the joint covariance carries the cross terms", {
	skip_on_cran()
	options_w = options(lame.eb_withhold_intercept = FALSE)
	on.exit(options(options_w), add = TRUE)
	d = eb_sim_fac(n = 20, Tt = 3, seed = 33, R = 1, suv = 0.6, dcor = 0.6)
	f1 = eb_fit(d, R = 1)
	f0 = eb_fit(d, R = 1, dcor = FALSE)
	expect_true(f1$converged && f0$converged)
	pr = f1$eb$pair
	expect_false(is.null(pr))
	expect_true(pr$rho > 0 && abs(pr$rho - 0.6) < 4 * pr$sd)
	expect_true(all(is.finite(f1$eb$se)) && all(f1$eb$se > 0))
	# the cross information of the pairs makes the standard errors smaller or equal, never wild
	expect_true(all(f1$eb$se / f0$eb$se > 0.8 & f1$eb$se / f0$eb$se < 1.1))
	# the joint covariance of the coefficients is symmetric positive definite
	V = f1$eb$vcov
	expect_equal(V, t(V), tolerance = 1e-10)
	expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), 0)
	# the fit is its own fixed point
	f2 = suppressWarnings(suppressMessages(ame_als_refit(f1)))
	expect_equal(f2$iterations, 0L)
	expect_equal(unname(f2$VC["rho"]), unname(f1$VC["rho"]))
	# the cavity rule's scale uncertainty refits hold the correlation and the fit reports the scale's standard deviation
	expect_true(is.finite(f1$eb$factor_scale$sd_log))
})

test_that("moving positions with the dyadic correlation: the fit converges and the correlation is the panel's", {
	skip_on_cran()
	source_vsim = function() {
		set.seed(8)
		n = 16; Tt = 3
		U = array(rnorm(n * Tt, 0, 0.6), c(n, 1, Tt)); V = array(rnorm(n * Tt, 0, 0.6), c(n, 1, Tt))
		rn = sprintf("a%02d", seq_len(n))
		Y = X = vector("list", Tt)
		for (t in seq_len(Tt)) {
			xt = array(rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x1"))
			eta = -1 + 0.5 * xt[, , 1] + tcrossprod(U[, , t], V[, , t])
			err = matrix(rnorm(n * n), n, n)
			lo = lower.tri(err)
			err[lo] = 0.5 * t(err)[lo] + sqrt(1 - 0.25) * matrix(rnorm(n * n), n, n)[lo]
			y = 1 * (eta + err > 0)
			dimnames(y) = list(rn, rn)
			diag(y) = NA
			Y[[t]] = y; X[[t]] = xt
		}
		names(Y) = names(X) = paste0("t", seq_len(Tt))
		list(Y = Y, X = X)
	}
	d = source_vsim()
	f = eb_fit(d, R = 1, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_false(is.null(f$eb$pair))
	expect_true(abs(f$eb$pair$rho) < 0.98 && is.finite(f$eb$pair$sd))
	expect_equal(unname(f$VC["rho"]), f$eb$pair$rho)
	expect_true(all(is.finite(f$eb$se)))
	f0 = eb_fit(d, R = 1, dynamic_uv = TRUE, dcor = FALSE)
	expect_null(f0$eb$pair)
})

test_that("the pair problem under wide noise (strongly identified positions) is the independent cells' at a zero correlation, with no blow-up at extreme nodes", {
	set.seed(12)
	n = 8; Tt = 1
	y = array(rbinom(n * n * Tt, 1, 0.4), c(n, n, Tt))
	for (i in seq_len(n)) y[i, i, ] = NA
	pc = lame:::.lde_pair_context(y)
	ctx = list(pc = pc, rho_hold = NULL, dcor = TRUE, y = y)
	eta = array(rnorm(n * n * Tt, 0, 2.5), dim(y))
	V = array(runif(n * n * Tt, 9, 25), dim(y))
	wp = lame:::.lde_work_pair(ctx, y, eta, V, 0, 1L, "probit")
	wi = lame:::.lde_work_gh(y, eta, V, what = 1L, link = "probit")
	ok = is.finite(y)
	# nothing blows up: the extreme nodes (|eta| up to 40) carry a probability near 1e-300 and the log-domain integral must be right there
	for (f in c("z", "w", "score", "hprime")) expect_true(all(is.finite(wp[[f]][ok])), label = f)
	expect_lt(max(abs(wp$score - wi$score)[ok]), 2e-3)
	expect_lt(max(abs(wp$w - wi$w)[ok]), 5e-2 * max(wi$w))
	expect_lt(max(abs(wp$cw)), 5e-2 * max(wi$w))
	# the log probability at extreme predictors, either order of the limits
	expect_equal(
		lame:::.ae_pair_logphi2(31, -5.56, 0),
		pnorm(-5.56, log.p = TRUE) + pnorm(31, log.p = TRUE),
		tolerance = 1e-4
	)
	expect_equal(
		lame:::.ae_pair_logphi2(-5.56, 31, 0),
		pnorm(-5.56, log.p = TRUE) + pnorm(31, log.p = TRUE),
		tolerance = 1e-4
	)
	expect_equal(
		lame:::.ae_pair_logphi2(25, -9, 0),
		pnorm(-9, log.p = TRUE) + pnorm(25, log.p = TRUE),
		tolerance = 1e-4
	)
})
