# collapsed per-actor update of the binary samplers (R/binary_actor_collapsed_move.R): directed and symmetric,
# lame() (centred levels) and ame() (unconstrained levels), complete panels and panels with absent actors
skip_on_cran()

# tiny binary panel: the last actor has no events; with `absent`, actor 3 misses the second period (or the whole cross-section),
# one cell is missing, and a symmetric panel has one pair observed only in its lower cell and one only in its upper cell
am_panel = function(n, nt, sym, absent) {
	Y = array(NA_real_, c(n, n, nt))
	for (t in seq_len(nt)) {
		y = matrix(rbinom(n * n, 1, .3), n)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		y[n, ] = 0; y[, n] = 0; diag(y) = NA
		Y[, , t] = y
	}
	if (absent) {
		if (nt >= 2L) { Y[3, , 2] = NA; Y[, 3, 2] = NA } else if (!sym) { Y[3, , 1] = NA; Y[, 3, 1] = NA }
		Y[2, 1, 1] = NA
		if (sym) Y[1, 4, 1] = NA
	}
	Y
}

test_that('the collapsed actor visit keeps the centring, the predictor and the signs', {
	skip_on_cran()
	set.seed(481516); n = 5L; nt = 3L
	for (sym in c(FALSE, TRUE)) for (absent in c(FALSE, TRUE)) for (dyn in c(FALSE, TRUE)) for (rank in c(0L, 2L)) {
		Y = am_panel(n, nt, sym, absent)
		if (absent) { Y[4, , 3] = NA; Y[, 4, 3] = NA }
		x = array(rnorm(n * n * nt), c(n, n, nt))
		if (sym) for (t in 1:nt) { m = x[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; x[, , t] = m }
		cx = if (sym) .lame_actor_move_context_sym(Y) else .lame_actor_move_context(Y)
		U = if (rank) matrix(rnorm(n * rank, sd = .5), n) else NULL
		Ld = if (rank) c(1.2, -.7)[seq_len(rank)] else NULL
		V = if (rank) (if (sym) U %*% diag(Ld, rank) else matrix(rnorm(n * rank, sd = .5), n)) else NULL
		la = rnorm(n, sd = .4); la = la - mean(la)
		lb = if (sym) la else { v = rnorm(n, sd = .4); v - mean(v) }
		dev = function() { d = matrix(rnorm(n * nt, sd = .2), n); sweep(d, 2, colMeans(d)) }
		am = if (dyn) matrix(la, n, nt) + dev() else NULL
		bm = if (dyn) (if (sym) am else matrix(lb, n, nt) + dev()) else NULL
		beta = c(-1, .4)
		predictor = function(beta, A, B, U, V) {
			z = array(0, c(n, n, nt))
			for (t in 1:nt) {
				z[, , t] = beta[1] + beta[2] * x[, , t] + outer(A[, t], B[, t], '+')
				if (!is.null(U)) z[, , t] = z[, , t] + tcrossprod(U, V)
			}
			z
		}
		A = if (dyn) am else matrix(la, n, nt); B = if (sym) A else if (dyn) bm else matrix(lb, n, nt)
		EZ = predictor(beta, A, B, U, V); Z = EZ + array(rnorm(n * n * nt), c(n, n, nt))
		if (!sym) { Z[cx$i1] = abs(Z[cx$i1]); Z[cx$i0] = -abs(Z[cx$i0]) }
		Sab = matrix(c(1, .3, .3, .8), 2); iV0 = matrix(c(1.2, .2, .2, 1), 2)
		for (k in 1:6) {
			if (sym) {
				o = .lame_actor_move_visit_sym(Z, EZ, Y, .9, beta, la, am, U, Ld, 1.1, iV0, k %% 2 == 0, cx)
				Z = o$Z; EZ = o$EZ; beta = o$beta; la = o$la; lb = la; am = o$a_mat; bm = am
				if (rank) { U = o$U; V = U %*% diag(Ld, rank) }
			} else {
				o = .lame_actor_move_visit(Z, EZ, -.35, beta, la, lb, am, bm, U, V, Sab, iV0, k %% 2 == 0, cx)
				Z = o$Z; EZ = o$EZ; beta = o$beta; la = o$la; lb = o$lb; am = o$a_mat; bm = o$b_mat; U = o$U; V = o$V
			}
			A = if (dyn) am else matrix(la, n, nt); B = if (sym) A else if (dyn) bm else matrix(lb, n, nt)
			expect_lt(max(abs(EZ - predictor(beta, A, B, U, V))), 1e-10)
			expect_lt(abs(sum(la)), 1e-10); expect_lt(abs(sum(lb)), 1e-10)
			if (dyn) { expect_lt(max(abs(colSums(am))), 1e-10); expect_lt(max(abs(colSums(bm))), 1e-10) }
			if (!sym) { expect_true(all(Z[cx$i1] > 0)); expect_true(all(Z[cx$i0] < 0)) }
			expect_true(all(is.finite(Z)))
		}
		expect_gt(o$stats[['level_calls']], 0); expect_equal(o$stats[['stuck']], 0)
		if (rank) expect_gt(o$stats[['uv_calls']], 0) else expect_equal(o$stats[['uv_calls']], 0)
	}
})

test_that('absent actor-periods and one-sided missing pairs only remove cells from the likelihood', {
	set.seed(5); n = 5L; nt = 2L
	Y = am_panel(n, nt, FALSE, TRUE)
	cx = .lame_actor_move_context(Y)
	obs = is.finite(Y); for (t in 1:nt) diag(obs[, , t]) = FALSE
	for (i in seq_len(n)) {
		expect_equal(length(cx$rows[[i]]$idx), sum(obs[i, , ]))
		expect_equal(length(cx$cols[[i]]$idx), sum(obs[, i, ]))
	}
	# the absent actor-period has no cell in its row or column of that period, the others keep theirs
	row3 = matrix(obs[3, , ], n, nt); expect_equal(sum(row3[, 2]), 0); expect_gt(sum(row3[, 1]), 0)
	expect_equal(sum(cx$rows[[3]]$idx > n), sum(row3[, 2:nt]))
	# symmetric: a pair is observed through its upper cell, or through its lower cell when that is missing
	Ys = am_panel(n, 1L, TRUE, TRUE)
	cs = .lame_actor_move_context_sym(Ys)
	expect_true(1 %in% cs$rows[[4]]$partner) # pair {1, 4}: upper cell missing, lower cell observed
	expect_true(1 %in% cs$rows[[2]]$partner) # pair {1, 2}: lower cell missing, upper cell observed
	Ysm = Ys; Ysm[4, 1, 1] = NA # the pair {1, 4} with neither cell
	expect_false(1 %in% .lame_actor_move_context_sym(Ysm)$rows[[4]]$partner)
})

test_that('the elliptical slice step draws the Gaussian posterior it is given', {
	set.seed(31337)
	# prior N(mu, root root'), Gaussian likelihood exp(-(x - c)' P (x - c) / 2): the posterior is Gaussian in closed form
	mu = c(1, -1); root = t(chol(matrix(c(4, 1.2, 1.2, 1), 2))); S0 = tcrossprod(root)
	cen = c(3, 0); P = matrix(c(3, -1, -1, 2), 2)
	loglik = function(x) -.5 * drop(t(x - cen) %*% P %*% (x - cen))
	post_cov = solve(solve(S0) + P); post_mean = drop(post_cov %*% (solve(S0) %*% mu + P %*% cen))
	x = mu; ll = loglik(x); nstep = 10000L; keep = matrix(NA_real_, nstep, 2); stuck = 0L
	for (k in seq_len(nstep)) { o = .lame_ess_step(x, mu, root, loglik, ll); x = o$x; ll = o$ll; keep[k, ] = x; stuck = stuck + o$stuck }
	expect_equal(stuck, 0L)
	ess = apply(keep, 2, coda::effectiveSize)
	expect_lt(max(abs(colMeans(keep) - post_mean) / (sqrt(diag(post_cov)) / sqrt(ess))), 4.5)
	expect_lt(max(abs(cov(keep) / post_cov - 1)), .1)
})

test_that('the level translation prior equals the shifted joint Gaussian prior', {
	set.seed(1103); n = 6L; i = 4L
	iS = solve(matrix(c(1.3, -.4, -.4, .7), 2)); iV0 = crossprod(matrix(rnorm(9), 3)) + diag(.5, 3)
	for (centred in c(TRUE, FALSE)) for (row in c(TRUE, FALSE)) {
		la = rnorm(n); lb = rnorm(n); beta = c(-.7, .3, 1.1)
		if (centred) { la = la - mean(la); lb = lb - mean(lb) }
		joint = function(delta) {
			w = numeric(n); w[i] = 1; if (centred) w = w - 1 / n
			a = la + if (row) delta * w else 0; b = lb + if (row) 0 else delta * w
			be = beta; if (centred) be[1] = be[1] + delta / n
			-.5 * (sum((cbind(a, b) %*% iS) * cbind(a, b)) + if (centred) drop(t(be) %*% iV0 %*% be) else 0)
		}
		pr = .lame_actor_move_level_prior(i, row, la, lb, iS, iV0, beta, n, centred)
		for (d in c(-2, -.5, .3, 1.7))
			expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
	}
	# symmetric: one effect per actor, entering the row and the column, the intercept moves by 2 delta / n
	va = 1.7
	for (centred in c(TRUE, FALSE)) {
		la = rnorm(n); if (centred) la = la - mean(la); beta = c(-.7, .3, 1.1)
		joint = function(delta) {
			w = numeric(n); w[i] = 1; if (centred) w = w - 1 / n
			a = la + delta * w; be = beta; if (centred) be[1] = be[1] + 2 * delta / n
			-.5 * (sum(a^2) / va + if (centred) drop(t(be) %*% iV0 %*% be) else 0)
		}
		pr = .lame_actor_move_level_prior_sym(i, la, va, iV0, beta, n, centred)
		for (d in c(-2, -.5, .3, 1.7))
			expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
	}
})

test_that('the factor priors of the positions are the sampler priors', {
	set.seed(77); n = 9L; R = 2L
	U = matrix(rnorm(n * R), n); V = matrix(rnorm(n * R), n)
	# ame: inverse Wishart with kappa0 + n degrees of freedom and scale kappa0 Suv0 + W'W, at every sweep
	pri = list(df = 6, scale = 6 * diag(2 * R))
	draws = replicate(4000, .lame_actor_move_suv(U, V, n, R, FALSE, pri))
	W = cbind(U, V); mean_theory = (pri$scale + crossprod(W)) / (n + pri$df - 2 * R - 1)
	expect_lt(max(abs(apply(draws, c(1, 2), mean) / mean_theory - 1)), .1)
	# lame: the diffuse n I before the second half of burn-in
	expect_equal(.lame_actor_move_suv(U, V, n, R, FALSE), diag(n, 2 * R))
	# symmetric: Gamma((2 + n) / 2, (1 + sum u^2) / 2) column precisions, 1 / n before the second half of burn-in
	expect_equal(.lame_actor_move_tau(U, n, R, FALSE), rep(1 / n, R))
	tau = replicate(4000, .lame_actor_move_tau(U, n, R, TRUE))
	expect_lt(max(abs(rowMeans(tau) / ((2 + n) / (1 + colSums(U^2))) - 1)), .06)
})

# exact-posterior check of the visit alone on a tiny model with absent actor-periods: the chain of visits (fixed hyperparameters)
# against self-normalised importance sampling from the prior; the chain-mean error is the spread of independent chains
am_kernel_z = function(mode, R, rho = 0, dyn = FALSE, seed = 1L, nch = 12L, nvis = 1500L, nmc = 30000L) {
	sym = mode %in% c('lame_sym', 'ame_sym'); centred = mode %in% c('lame_dir', 'lame_sym')
	n = 4L; N = if (centred) 2L else 1L
	set.seed(seed)
	Y = am_panel(n, N, sym, TRUE)
	xx = array(rnorm(n * n * N), c(n, n, N))
	if (sym) for (t in seq_len(N)) { m = xx[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; xx[, , t] = m }
	b2 = .5; Sab = matrix(c(1.2, .4, .4, .8), 2); va = 1; iV0 = matrix(c(1.5, .3, .3, 1.1), 2); iS = solve(Sab)
	Ld = c(1.3, -.8)[seq_len(R)]; tau_fix = c(1.6, .7)[seq_len(R)]
	Suv_fix = NULL
	if (R > 0 && !sym) { B = matrix(rnorm((2 * R)^2), 2 * R); Suv_fix = (crossprod(B) / (2 * R) + diag(.3, 2 * R)) * .8 }
	if (R > 0) testthat::local_mocked_bindings(
		.lame_actor_move_tau = function(U, n, R, shrink) tau_fix,
		.lame_actor_move_suv = function(U, V, n, R, shrink, prior = NULL) Suv_fix,
		.package = 'lame'
	)
	H = qr.Q(qr(contr.helmert(n)))
	if (dyn) {
		da = matrix(rnorm(n * N, sd = .3), n); da = sweep(da, 2, colMeans(da))
		db = matrix(rnorm(n * N, sd = .3), n); db = sweep(db, 2, colMeans(db))
	} else da = db = matrix(0, n, N)
	cx = if (sym) .lame_actor_move_context_sym(Y) else .lame_actor_move_context(Y)
	mkEZ = function(b0, la, lb, U, V) {
		z = array(0, c(n, n, N))
		for (t in 1:N) {
			z[, , t] = b0 + b2 * xx[, , t] + outer(la + da[, t], (if (sym) la + da[, t] else lb + db[, t]), '+')
			if (R > 0) z[, , t] = z[, , t] + tcrossprod(U, V)
		}
		z
	}
	iu = which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	loglik = function(eta) {
		ll = 0
		for (t in 1:N) {
			e1 = eta[cbind(iu[, 1], iu[, 2], t)]; e2 = eta[cbind(iu[, 2], iu[, 1], t)]
			y1 = Y[cbind(iu[, 1], iu[, 2], t)]; y2 = Y[cbind(iu[, 2], iu[, 1], t)]
			if (sym) {
				y = ifelse(is.finite(y1), y1, y2); ok = is.finite(y)
				ll = ll + sum(pnorm((2 * y[ok] - 1) * (e1[ok] + e2[ok]) / 2, log.p = TRUE)); next
			}
			both = is.finite(y1) & is.finite(y2); o1 = is.finite(y1) & !is.finite(y2); o2 = !is.finite(y1) & is.finite(y2)
			s1 = 2 * y1[both] - 1; s2 = 2 * y2[both] - 1
			ll = ll + sum(log(pmax(pbivnorm::pbivnorm(s1 * e1[both], s2 * e2[both], s1 * s2 * rho), 1e-300))) +
				sum(
					pnorm((2 * y1[o1] - 1) * e1[o1], log.p = TRUE)
				) + sum(pnorm((2 * y2[o2] - 1) * e2[o2], log.p = TRUE))
		}
		ll
	}
	b0_mean = -iV0[1, 2] * b2 / iV0[1, 1]
	unpack = function(zz) {
		k = 0; take = function(m) { out = zz[k + seq_len(m)]; k <<- k + m; out }
		b0 = if (centred) take(1) else .3
		if (centred) { la = drop(H %*% take(n - 1)); lb = if (sym) la else drop(H %*% take(n - 1)) } else { la = take(n); lb = if (sym) la else take(n) }
		U = V = NULL
		if (R > 0) { U = matrix(take(n * R), n, R); V = if (sym) U %*% diag(Ld, R) else matrix(take(n * R), n, R) }
		list(b0 = b0, la = la, lb = lb, U = U, V = V)
	}
	pack = function(b0, la, lb, U, V) c(
		if (centred) b0,
		if (centred) crossprod(H, la) else la,
		if (!sym) { if (centred) crossprod(H, lb) else lb },
		if (R > 0) c(U),
		if (R > 0 && !sym) c(V)
	)
	logprior = function(zz) {
		k = 0; take = function(m) { out = zz[k + seq_len(m)]; k <<- k + m; out }
		lp = 0
		if (centred) lp = lp - .5 * iV0[1, 1] * (take(1) - b0_mean)^2
		if (sym) lp = lp - .5 * sum(take(if (centred) n - 1 else n)^2) / va
		else { za = take(if (centred) n - 1 else n); zb = take(if (centred) n - 1 else n); lp = lp - .5 * sum(
			(cbind(za, zb) %*% iS) * cbind(za, zb)
		) }
		if (R > 0) {
			Uc = matrix(take(n * R), n, R)
			if (sym) lp = lp - .5 * sum(sweep(Uc^2, 2, tau_fix, '*'))
			else { W_ = cbind(Uc, matrix(take(n * R), n, R)); lp = lp - .5 * sum((W_ %*% solve(Suv_fix)) * W_) }
		}
		lp
	}
	run_chain = function() {
		st = list(b0 = if (centred) -1 else .3, la = c(.2, -.3, .4, -.3), lb = c(-.1, .2, .1, -.2), U = NULL, V = NULL)
		if (centred) { st$la = st$la - mean(st$la); st$lb = st$lb - mean(st$lb) }
		if (R > 0) { st$U = matrix(rnorm(n * R, sd = .5), n); st$V = if (sym) st$U %*% diag(Ld, R) else matrix(rnorm(n * R, sd = .5), n) }
		if (sym) st$lb = st$la
		beta = c(st$b0, b2)
		EZ = mkEZ(st$b0, st$la, st$lb, st$U, st$V); Z = EZ + array(rnorm(n * n * N), c(n, n, N))
		if (sym) for (t in 1:N) Z[, , t] = rZ_bin_sym_fc(Z[, , t], EZ[, , t], Y[, , t], .9)
		else { Z[cx$i1] = abs(Z[cx$i1]); Z[cx$i0] = -abs(Z[cx$i0]) }
		keep = matrix(NA_real_, nvis, length(pack(st$b0, st$la, st$lb, st$U, st$V)))
		for (k in seq_len(nvis + 100L)) {
			if (sym) {
				o = .lame_actor_move_visit_sym(
					Z,
					EZ,
					Y,
					.9,
					beta,
					st$la,
					if (dyn) st$la + da else NULL,
					st$U,
					if (R > 0) Ld else NULL,
					va,
					iV0,
					TRUE,
					cx,
					centred = centred
				)
				Z = o$Z; EZ = o$EZ; beta = o$beta; st$la = o$la; st$lb = st$la; if (R > 0) { st$U = o$U; st$V = st$U %*% diag(Ld, R) }
			} else {
				o = .lame_actor_move_visit(
					Z,
					EZ,
					rho,
					beta,
					st$la,
					st$lb,
					if (dyn) st$la + da else NULL,
					if (dyn) st$lb + db else NULL,
					st$U,
					st$V,
					Sab,
					iV0,
					TRUE,
					cx,
					centred = centred
				)
				Z = o$Z; EZ = o$EZ; beta = o$beta; st$la = o$la; st$lb = o$lb; if (R > 0) { st$U = o$U; st$V = o$V }
			}
			st$b0 = beta[1]
			if (k > 100L) keep[k - 100L, ] = pack(st$b0, st$la, st$lb, st$U, st$V)
		}
		keep
	}
	keeps = lapply(seq_len(nch), function(i) run_chain())
	pool = do.call(rbind, keeps); d0 = ncol(pool)
	sub = pool[seq(1, nrow(pool), by = 8), , drop = FALSE]
	mu_p = colMeans(sub); Lp = t(chol(cov(sub) * 1.3 + diag(1e-8, d0))); dfp = 6
	draws = matrix(NA_real_, nmc, d0); logw = numeric(nmc)
	for (m in seq_len(nmc)) {
		g = rchisq(1, dfp); e = rnorm(d0)
		zz = mu_p + drop(Lp %*% e) / sqrt(g / dfp)
		s_ = unpack(zz)
		logw[m] = loglik(mkEZ(s_$b0, s_$la, s_$lb, s_$U, s_$V)) + logprior(zz) - (-(dfp + d0) / 2 * log1p(sum(e^2) / g))
		draws[m, ] = zz
	}
	w = exp(logw - max(logw)); w = w / sum(w)
	zfun = function(Fc, Fd) {
		ref = colSums(Fd * w); cm = t(sapply(Fc, colMeans))
		(colMeans(cm) - ref) / sqrt(apply(cm, 2, var) / length(Fc) + colSums(sweep(Fd, 2, ref)^2 * w^2))
	}
	act_of = function(v) { s_ = unpack(v); c(b0 = s_$b0, a = s_$la, b = if (!sym) s_$lb) }
	acts = lapply(keeps, function(k) t(apply(k, 1, act_of))); ad = t(apply(draws, 1, act_of))
	za = zfun(acts, ad); za = za[apply(ad, 2, sd) > 0]; if (!centred) za = za[-1]
	zi = NULL
	if (R > 0) {
		inv_fun = function(v) { s_ = unpack(v); P_ = tcrossprod(s_$U, s_$V)
			if (sym) c(P_[upper.tri(P_)], rowSums(s_$U^2)) else c(
				P_[upper.tri(P_)],
				P_[lower.tri(P_)],
				rowSums(s_$U^2),
				rowSums(s_$V^2)
			) }
		zi = zfun(lapply(keeps, function(k) t(apply(k, 1, inv_fun))), t(apply(draws, 1, inv_fun)))
	}
	list(levels = max(abs(za)), invariants = if (is.null(zi)) 0 else max(abs(zi)), ess_is = 1 / sum(w^2))
}

test_that('the collapsed actor visit leaves the exact conditional posterior invariant', {
	skip_on_cran()
	skip_if_not_installed('pbivnorm')
	cfg = list(
		list('lame_dir', 0L, -.4, TRUE), # directed, centred levels, dynamic paths, dyadic correlation
		list('lame_sym', 1L, 0, TRUE), # symmetric, one effect per actor, eigenmodel positions
		list('ame_dir', 1L, -.4, FALSE), # cross-section, unconstrained levels, hierarchical factor prior
		list('ame_sym', 0L, 0, FALSE)
	)
	for (k in seq_along(cfg)) {
		z = am_kernel_z(cfg[[k]][[1]], cfg[[k]][[2]], rho = cfg[[k]][[3]], dyn = cfg[[k]][[4]], seed = 100L + k)
		expect_gt(z$ess_is, 2000)
		expect_lt(z$levels, 4.5)
		expect_lt(z$invariants, 4.5)
	}
})

test_that('lame() runs the collapsed actor update on incomplete and symmetric panels and records its work', {
	skip_on_cran()
	set.seed(7305); nm = sprintf('a%02d', 1:9)
	# changing composition: the actors present differ across the three periods
	sets = list(1:8, c(2:9), c(1:7, 9))
	mk = function(s, sym) {
		y = matrix(rbinom(length(s)^2, 1, .15), length(s), dimnames = list(nm[s], nm[s]))
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		diag(y) = NA; y[nrow(y), ] = ifelse(is.na(y[nrow(y), ]), NA, 0); y[, nrow(y)] = ifelse(is.na(y[, nrow(y)]), NA, 0); y
	}
	for (sym in c(FALSE, TRUE)) {
		Y = lapply(sets, mk, sym = sym)
		# a symmetric fit needs a symmetric covariate
		X = lapply(sets, function(s) {
			x = matrix(rnorm(length(s)^2), length(s))
			if (sym) x = (x + t(x)) / 2
			array(x, c(length(s), length(s), 1), dimnames = list(nm[s], nm[s], 'x'))
		})
		run = function(prior, dyn, R) suppressMessages(suppressWarnings(lame(
			Y,
			Xdyad = X,
			family = 'binary',
			R = R,
			symmetric = sym,
			dynamic_ab = dyn,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior
		)))
		for (dyn in c(FALSE, TRUE)) for (R in c(0, 1)) {
			on = run(list(actor_move_every = 2), dyn, R)
			expect_equal(unname(on$actor_move['visits']), 40)
			expect_equal(unname(on$actor_move['level_calls']), 40 * (if (sym) 1 else 2) * 9)
			expect_equal(unname(on$actor_move['uv_calls'] > 0), R > 0)
			expect_equal(unname(on$actor_move['stuck']), 0)
			if (!sym) {
				expect_equal(unname(on$actor_move['rho_attempts']), 40 * 4)
				expect_equal(
					unname(run(list(actor_move_every = 2, actor_move_rho = FALSE), dyn, R)$actor_move['rho_attempts']),
					0
				)
			}
			expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
			expect_null(run(list(actor_move = FALSE), dyn, R)$actor_move)
			# a schedule that never reaches a visit consumes no random numbers: the draws are those of the switch-off
			never = run(list(actor_move_every = 10^6), dyn, R); off = run(list(actor_move = FALSE), dyn, R)
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
		}
	}
	expect_error(run(list(actor_move = NA), FALSE, 0), 'actor_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(
		run(list(actor_move_every = 0), FALSE, 0),
		'actor_move_every must be a positive whole number',
		fixed = TRUE
	)
	expect_error(run(list(actor_move_rho = 1), FALSE, 0), 'actor_move_rho must be TRUE or FALSE', fixed = TRUE)
})

test_that('ame() runs the collapsed actor update, directed and symmetric, and records its work', {
	skip_on_cran()
	set.seed(6012); n = 10; nm = sprintf('a%02d', 1:n)
	for (sym in c(FALSE, TRUE)) {
		Y = matrix(rbinom(n * n, 1, .15), n, dimnames = list(nm, nm))
		if (sym) Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
		diag(Y) = NA; Y[n, ] = ifelse(is.na(Y[n, ]), NA, 0); Y[, n] = ifelse(is.na(Y[, n]), NA, 0)
		Y[3, ] = NA; Y[, 3] = NA; Y[2, 5] = NA; if (sym) Y[5, 2] = NA
		xd = array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, 'x'))
		if (sym) { m = xd[, , 1]; m[lower.tri(m)] = t(m)[lower.tri(m)]; xd[, , 1] = m }
		fit = function(prior, R) suppressMessages(ame(
			Y,
			Xdyad = xd,
			family = 'binary',
			R = R,
			symmetric = sym,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = prior
		))
		run = function(prior, R) suppressWarnings(fit(prior, R))
		for (R in c(0, 2)) {
			on = run(list(actor_move_every = 2), R)
			expect_equal(unname(on$actor_move['visits']), 40)
			expect_equal(unname(on$actor_move['level_calls']), 40 * (if (sym) 1 else 2) * n)
			expect_equal(unname(on$actor_move['uv_calls'] > 0), R > 0)
			expect_equal(unname(on$actor_move['stuck']), 0)
			expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
			expect_null(run(list(actor_move = FALSE), R)$actor_move)
			# the default is on
			expect_gt(unname(run(list(), R)$actor_move['visits']), 0)
			never = run(list(actor_move_every = 10^6), R); off = run(list(actor_move = FALSE), R)
			expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
			expect_identical(never$APM, off$APM); expect_identical(never$UVPM, off$UVPM)
		}
	}
	# only the binary family has the update; the prior names are checked
	yn = Y; yn[is.finite(yn)] = rnorm(sum(is.finite(yn)))
	fn = suppressMessages(
		suppressWarnings(ame(yn, family = 'normal', burn = 5, nscan = 20, odens = 5, verbose = FALSE, gof = FALSE))
	)
	expect_null(fn$actor_move)
	expect_error(run(list(actor_move = NA), 0), 'actor_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(run(list(actor_move = 1), 0), 'actor_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(run(list(actor_move_every = 0), 0), 'actor_move_every must be a positive whole number', fixed = TRUE)
	expect_warning(fit(list(actor_moves = TRUE), 0), 'Unknown')
	expect_no_warning(fit(list(actor_move = FALSE, actor_move_every = 5), 0))
})

test_that('the default symmetric update recovers a heterophilous dimension', {
	skip_on_cran()
	# the panel of test_symmetric_signed_start.R, fitted with the collapsed actor update on (the default)
	set.seed(1005); n = 18L; nt = 9L; nm = sprintf('c%02d', 1:n)
	a = rnorm(n, 0, 0.5); u = matrix(rnorm(n * 2L), n, 2L)
	EZ = -1 + outer(a, a, '+') + u %*% diag(c(1.2, -0.8)) %*% t(u)
	Y = lapply(seq_len(nt), function(t) {
		e = matrix(rnorm(n * n), n, n); e[lower.tri(e)] = t(e)[lower.tri(e)]
		y = 1 * (EZ + e > 0); diag(y) = NA; dimnames(y) = list(nm, nm); y
	})
	fit = lame(
		Y,
		family = 'binary',
		symmetric = TRUE,
		R = 2L,
		burn = 200L,
		nscan = 1000L,
		odens = 10L,
		seed = 12L,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)
	expect_gt(fit$actor_move[['uv_calls']], 0)
	expect_equal(sort(sign(fit$L)), c(-1, 1))
	expect_gt(cor(EZ[upper.tri(EZ)], fit$EZ[[1]][upper.tri(EZ)]), 0.80)
})
