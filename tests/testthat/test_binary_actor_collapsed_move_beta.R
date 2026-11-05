# collapsed per-actor update of the binary samplers with dynamic coefficients and fewer additive margins
# (R/binary_actor_collapsed_move.R): the prior of a moving intercept under a level translation, the exactness of the visit, and the fits

# tiny binary panel: the last actor has no events, actor 3 is missing from the second period, and one cell is missing
amb_tiny_panel = function(n, nt, sym = FALSE) {
	Y = array(NA_real_, c(n, n, nt))
	for (t in seq_len(nt)) {
		y = matrix(rbinom(n * n, 1, .3), n)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		y[n, ] = 0; y[, n] = 0; diag(y) = NA
		Y[, , t] = y
	}
	if (nt >= 2L) { Y[3, , 2] = NA; Y[, 3, 2] = NA }
	Y[2, 1, 1] = NA
	if (sym) Y[1, 4, 1] = NA # a pair observed in one cell only
	Y
}

# log density, up to a constant, of an intercept path under the coefficient sampler's prior:
# beta_1 ~ N(m0, v0), beta_t = rho beta_{t-1} + e_t, e_t ~ N(0, sd2)
amb_logpath = function(path, rho, sd2, m0, v0) {
	k = length(path)
	lp = -.5 * (path[1] - m0)^2 / v0
	if (k > 1) lp = lp - .5 * sum((path[-1] - rho * path[-k])^2) / sd2
	lp
}

test_that('a shift of the intercept path changes its log prior by -(q c^2 + 2 l c) / 2', {
	set.seed(2301)
	for (rho in c(1, .6, -.3)) for (k in c(2L, 5L)) {
		path = rnorm(k); m0 = .4; v0 = 3; sd2 = .25
		ic = list(kind = 'path', path = path, rho = rho, sd2 = sd2, m0 = m0, v0 = v0)
		lq = .lame_actor_move_icpt(ic, NULL, NULL)
		for (cc in c(-1.7, -.2, .5, 2.3))
			expect_equal(
				amb_logpath(path + cc, rho, sd2, m0, v0) - amb_logpath(path, rho, sd2, m0, v0),
				-(lq[["q"]] * cc^2 + 2 * lq[["l"]] * cc) / 2,
				tolerance = 1e-10
			)
		# a random walk changes only the prior of its first state
		if (rho == 1) expect_equal(lq[['q']], 1 / v0)
		# the shift moves the whole path
		expect_equal(.lame_actor_move_icpt_shift(ic, NULL, .7)$ic$path, path + .7)
		# the coefficient sampler's own quadratic form (R/binary_marginal_move.R), for the zero prior mean
		ic0 = ic; ic0$m0 = 0
		lq0 = .lame_actor_move_icpt(ic0, NULL, NULL)
		qform = function(p) .lame_path_scale_terms(
			numeric(0),
			matrix(p, ncol = 1),
			rho,
			sqrt(sd2),
			diag(1),
			diag(v0, 1),
			1,
			1
		)$q
		for (cc in c(-1.7, .5)) expect_equal(
			-.5 * (qform(path + cc) - qform(path)),
			-(lq0[['q']] * cc^2 + 2 * lq0[['l']] * cc) / 2,
			tolerance = 1e-12
		)
	}
	# a static intercept: N(0, 1 / prec)
	ic = list(kind = 'static', value = .8, prec = .35)
	lq = .lame_actor_move_icpt(ic, NULL, NULL)
	for (cc in c(-1, .6)) expect_equal(
		-.5 * .35 * ((.8 + cc)^2 - .8^2),
		-(lq[["q"]] * cc^2 + 2 * lq[["l"]] * cc) / 2,
		tolerance = 1e-10
	)
	expect_equal(.lame_actor_move_icpt_shift(ic, NULL, .2)$ic$value, 1)
	# the coefficient prior of a static fit is unchanged
	iV0 = crossprod(matrix(rnorm(9), 3)) + diag(.5, 3); beta = c(-.7, .3, 1.1)
	lq = .lame_actor_move_icpt(NULL, iV0, beta)
	expect_equal(unname(lq[['l']]), drop(iV0 %*% beta)[1]); expect_equal(unname(lq[['q']]), iV0[1, 1])
	expect_equal(.lame_actor_move_icpt_shift(NULL, beta, .4)$beta[1], beta[1] + .4)
})

test_that('the level translation prior with a moving intercept equals the shifted joint Gaussian prior', {
	set.seed(2302); n = 6L; i = 4L; N = 4L
	iS = solve(matrix(c(1.3, -.4, -.4, .7), 2))
	for (kind in c('path', 'static')) for (row in c(TRUE, FALSE)) {
		la = rnorm(n); lb = rnorm(n); la = la - mean(la); lb = lb - mean(lb)
		ic = if (kind == 'path') list(kind = 'path', path = rnorm(N), rho = .7, sd2 = .4, m0 = .2, v0 = 5)
		else list(kind = 'static', value = .6, prec = .3)
		logic = function(sh) if (kind == 'path') amb_logpath(sh$path, .7, .4, .2, 5) else -.5 * .3 * sh$value^2
		joint = function(delta) {
			w = rep(-1 / n, n); w[i] = 1 - 1 / n
			a = la + if (row) delta * w else 0; b = lb + if (row) 0 else delta * w
			-.5 * sum((cbind(a, b) %*% iS) * cbind(a, b)) + logic(.lame_actor_move_icpt_shift(ic, NULL, delta / n)$ic)
		}
		pr = .lame_actor_move_level_prior(i, row, la, lb, iS, NULL, NULL, n, TRUE, ic)
		for (d in c(-2, -.5, .3, 1.7))
			expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
	}
	va = 1.7
	for (kind in c('path', 'static')) {
		la = rnorm(n); la = la - mean(la)
		ic = if (kind == 'path') list(kind = 'path', path = rnorm(N), rho = 1, sd2 = .4, m0 = 0, v0 = 5)
		else list(kind = 'static', value = .6, prec = .3)
		logic = function(sh) if (kind == 'path') amb_logpath(sh$path, 1, .4, 0, 5) else -.5 * .3 * sh$value^2
		joint = function(delta) {
			w = rep(-1 / n, n); w[i] = 1 - 1 / n
			-.5 * sum((la + delta * w)^2) / va + logic(.lame_actor_move_icpt_shift(ic, NULL, 2 * delta / n)$ic)
		}
		pr = .lame_actor_move_level_prior_sym(i, la, va, NULL, NULL, n, TRUE, ic)
		for (d in c(-2, -.5, .3, 1.7))
			expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
	}
})

# exact-posterior check of the visit alone on a tiny model: the chain of visits (hyperparameters fixed) against
# self-normalised importance sampling from the exact likelihood and prior, with the spread of independent chains as the error
#   icpt  'iV0' the static coefficient prior, 'static' a static intercept of a dynamic_beta fit, 'rw1' / 'ar1' an
#         intercept path (the visit moves the path by a constant, so the harness holds the path's shape fixed and
#         checks the constant, with the levels), 'none' no intercept (uncentred levels, as ame())
#   lev   the additive margins in the model (sender, receiver; a symmetric fit has the sender entry only); a margin
#         out of the model is zero
#   R     rank of the static positions (0 or 1); dyn: dynamic_ab deviations (fixed); sym: a symmetric network
amb_kernel_z = function(
	icpt,
	lev = c(TRUE, TRUE),
	R = 0L,
	dyn = FALSE,
	seed = 1L,
	nch = 10L,
	nvis = 1200L,
	nmc = 25000L,
	rho = -.4,
	sym = FALSE
) {
	n = 4L; N = if (icpt %in% c('rw1', 'ar1')) 3L else 2L
	shp = c(.5, -.2, -.3)
	centred = icpt != 'none'
	set.seed(seed)
	Y = amb_tiny_panel(n, N, sym)
	xx = array(rnorm(n * n * N), c(n, n, N))
	if (sym) for (t in seq_len(N)) { m = xx[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; xx[, , t] = m }
	b2 = .5; Sab = matrix(c(1.2, .4, .4, .8), 2)
	if (!all(lev)) Sab = diag(c(1.2, .8)); iS = solve(Sab)
	va = 1
	iV0 = matrix(c(1.5, .3, .3, 1.1), 2); prec = .7; rho_p = if (icpt == 'ar1') .6 else 1; sd2 = .3; m0 = .2; v0 = 4
	b0_mean = -iV0[1, 2] * b2 / iV0[1, 1]
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
		db = if (sym) da else { m_ = matrix(rnorm(n * N, sd = .3), n); sweep(m_, 2, colMeans(m_)) }
	} else da = db = matrix(0, n, N)
	cx = if (sym) .lame_actor_move_context_sym(Y) else .lame_actor_move_context(Y)
	nb0 = if (centred) 1L else 0L
	mkEZ = function(b0, la, lb, U, V) {
		z = array(0, c(n, n, N))
		for (t in 1:N) {
			z[, , t] = (if (nb0 == 0L) .3 else if (icpt %in% c('rw1', 'ar1')) b0 + shp[t] else b0) + b2 * xx[, , t] +
				outer(la + da[, t], (if (sym) la + da[, t] else lb + db[, t]), '+')
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
	nl = if (centred) n - 1L else n
	unpack = function(zz) {
		k = 0; take = function(m) { out = zz[k + seq_len(m)]; k <<- k + m; out }
		b0 = if (nb0 > 0L) take(nb0) else NULL
		la = if (lev[1]) { z_ = take(nl); if (centred) drop(H %*% z_) else z_ } else numeric(n)
		lb = if (sym) la else if (lev[2]) { z_ = take(nl); if (centred) drop(H %*% z_) else z_ } else numeric(n)
		U = V = NULL
		if (R > 0) { U = matrix(take(n * R), n, R); V = if (sym) U %*% diag(Ld, R) else matrix(take(n * R), n, R) }
		list(b0 = b0, la = la, lb = lb, U = U, V = V)
	}
	pack = function(b0, la, lb, U, V) c(
		b0,
		if (lev[1]) { if (centred) crossprod(H, la) else la },
		if (lev[2] && !sym) { if (centred) crossprod(H, lb) else lb },
		if (R > 0) c(U),
		if (R > 0 && !sym) c(V)
	)
	logprior = function(zz) {
		s_ = unpack(zz); lp = 0
		if (icpt == 'iV0') lp = lp - .5 * iV0[1, 1] * (s_$b0 - b0_mean)^2
		if (icpt == 'static') lp = lp - .5 * prec * s_$b0^2
		if (icpt %in% c('rw1', 'ar1')) lp = lp + amb_logpath(shp + s_$b0, rho_p, sd2, m0, v0)
		k = 0; take = function(m) { out = zz[k + seq_len(m)]; k <<- k + m; out }
		if (nb0 > 0L) take(nb0)
		if (sym) { if (lev[1]) lp = lp - .5 * sum(take(nl)^2) / va }
		else if (all(lev)) { za = take(nl); zb = take(nl); lp = lp - .5 * sum((cbind(za, zb) %*% iS) * cbind(za, zb)) }
		else if (lev[1]) lp = lp - .5 * sum(take(nl)^2) / Sab[1, 1]
		else if (lev[2]) lp = lp - .5 * sum(take(nl)^2) / Sab[2, 2]
		if (R > 0) {
			Uc = matrix(take(n * R), n, R)
			if (sym) lp = lp - .5 * sum(sweep(Uc^2, 2, tau_fix, '*'))
			else { W_ = cbind(Uc, matrix(take(n * R), n, R)); lp = lp - .5 * sum((W_ %*% solve(Suv_fix)) * W_) }
		}
		lp
	}
	run_chain = function() {
		st = list(
			b0 = if (nb0 == 0L) NULL else -1,
			la = if (lev[1]) c(.2, -.3, .4, -.3) else numeric(n),
			lb = if (lev[2]) c(-.1, .2, .1, -.2) else numeric(n),
			U = NULL,
			V = NULL
		)
		if (centred) { if (lev[1]) st$la = st$la - mean(st$la); if (lev[2]) st$lb = st$lb - mean(st$lb) }
		if (sym) st$lb = st$la
		if (R > 0) { st$U = matrix(rnorm(n * R, sd = .5), n); st$V = if (sym) st$U %*% diag(Ld, R) else matrix(rnorm(n * R, sd = .5), n) }
		beta = c(if (icpt == 'iV0') st$b0 else 0, b2)
		EZ = mkEZ(st$b0, st$la, st$lb, st$U, st$V); Z = EZ + array(rnorm(n * n * N), c(n, n, N))
		if (sym) for (t in 1:N) Z[, , t] = rZ_bin_sym_fc(Z[, , t], EZ[, , t], Y[, , t], .9)
		else { Z[cx$i1] = abs(Z[cx$i1]); Z[cx$i0] = -abs(Z[cx$i0]) }
		keep = matrix(NA_real_, nvis, length(pack(st$b0, st$la, st$lb, st$U, st$V)))
		for (k in seq_len(nvis + 100L)) {
			ic = switch(
				icpt,
				static = list(kind = 'static', value = st$b0, prec = prec),
				rw1 = ,
				ar1 = list(kind = 'path', path = shp + st$b0, rho = rho_p, sd2 = sd2, m0 = m0, v0 = v0),
				NULL
			)
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
					centred = centred,
					ic = ic,
					lev = lev[1]
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
					centred = centred,
					ic = ic,
					lev = lev
				)
				Z = o$Z; EZ = o$EZ; beta = o$beta; st$la = o$la; st$lb = o$lb; if (R > 0) { st$U = o$U; st$V = o$V }
			}
			if (icpt == 'iV0') st$b0 = beta[1] else if (icpt == 'static') st$b0 = o$ic$value else if (nb0 > 0L) st$b0 = mean(o$ic$path)
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
	act_of = function(v) { s_ = unpack(v); c(b0 = s_$b0, a = if (lev[1]) s_$la, b = if (lev[2] && !sym) s_$lb) }
	za = 0
	if (any(lev) || nb0 > 0L) {
		acts = lapply(keeps, function(k) t(apply(k, 1, act_of))); ad = t(apply(draws, 1, act_of))
		za = zfun(acts, ad); za = za[apply(ad, 2, sd) > 0]
	}
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

test_that('the collapsed actor visit with a moving or static dynamic_beta intercept leaves the exact posterior invariant', {
	skip_on_cran()
	skip_if_not_installed('pbivnorm')
	cfg = list(
		list('rw1', c(TRUE, TRUE), 0L, FALSE), # random-walk intercept path
		list('ar1', c(TRUE, TRUE), 0L, FALSE), # AR(1) intercept path
		list('rw1', c(TRUE, TRUE), 0L, TRUE), # with dynamic additive effects
		list('static', c(TRUE, TRUE), 0L, FALSE), # static intercept next to moving coefficients
		list('ar1', c(TRUE, FALSE), 1L, FALSE), # one margin, static positions
		list('rw1', c(TRUE, TRUE), 1L, FALSE, TRUE), # symmetric: the intercept moves by 2 delta / n
		list('static', c(TRUE, TRUE), 0L, TRUE, TRUE)
	)
	for (k in seq_along(cfg)) {
		z = amb_kernel_z(
			cfg[[k]][[1]],
			cfg[[k]][[2]],
			cfg[[k]][[3]],
			cfg[[k]][[4]],
			seed = 200L + k,
			sym = length(cfg[[k]]) >= 5
		)
		expect_gt(z$ess_is, 1500)
		expect_lt(z$levels, 4.5)
		expect_lt(z$invariants, 4.5)
	}
})

test_that('the collapsed actor visit with fewer additive margins leaves the exact posterior invariant', {
	skip_on_cran()
	skip_if_not_installed('pbivnorm')
	cfg = list(
		list('none', c(FALSE, FALSE), 1L), # no margin: the position steps alone
		list('none', c(TRUE, FALSE), 1L), # sender margin only
		list('none', c(FALSE, TRUE), 0L), # receiver margin only, no positions
		list('iV0', c(FALSE, TRUE), 1L), # centred fit, receiver margin only
		list('none', c(FALSE, FALSE), 1L, TRUE)
	) # symmetric without the actor effect: positions alone
	for (k in seq_along(cfg)) {
		z = amb_kernel_z(
			cfg[[k]][[1]],
			cfg[[k]][[2]],
			cfg[[k]][[3]],
			FALSE,
			seed = 300L + k,
			sym = length(cfg[[k]]) >= 4
		)
		expect_gt(z$ess_is, 1500)
		expect_lt(z$levels, 4.5)
		expect_lt(z$invariants, 4.5)
	}
})

# a small directed panel of a few periods for the fits
amb_panel = function(n = 12L, nt = 5L, seed = 11L) {
	set.seed(seed); nm = sprintf('a%02d', 1:n)
	Y = lapply(1:nt, function(t) {
		y = matrix(rbinom(n * n, 1, .16), n, dimnames = list(nm, nm)); diag(y) = NA
		y[n, ] = ifelse(is.na(y[n, ]), NA, 0); y[, n] = ifelse(is.na(y[, n]), NA, 0); y
	})
	X = lapply(1:nt, function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, 'x')))
	names(Y) = names(X) = as.character(2000 + 1:nt)
	list(Y = Y, X = X)
}

test_that('lame() runs the collapsed actor update with dynamic coefficients and records its work', {
	skip_on_cran()
	d = amb_panel()
	run = function(prior, ...) suppressMessages(suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = 'binary',
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = prior,
		...
	)))
	for (kind in c('rw1', 'ar1')) for (blocks in list(
		'intercept',
		c('intercept', 'x_dyad'),
		'x_dyad'
	)) for (dab in c(FALSE, TRUE)) {
		on = run(list(actor_move_every = 2), R = 1, dynamic_beta = blocks, dynamic_beta_kind = kind, dynamic_ab = dab)
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['level_calls']), 40 * 2 * 12)
		expect_gt(unname(on$actor_move['uv_calls']), 0)
		expect_equal(unname(on$actor_move['stuck']), 0)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		off = run(list(actor_move = FALSE), R = 1, dynamic_beta = blocks, dynamic_beta_kind = kind, dynamic_ab = dab)
		expect_null(off$actor_move)
		# a schedule that never reaches a visit consumes no random numbers: the draws are those of the switch-off
		never = run(
			list(actor_move_every = 10^6),
			R = 1,
			dynamic_beta = blocks,
			dynamic_beta_kind = kind,
			dynamic_ab = dab
		)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
	}
	# the intercept path moves with the levels, so the stored intercept differs between the arms while every cell's predictor is unchanged
	on = run(list(actor_move_every = 1), R = 0, dynamic_beta = 'intercept', dynamic_beta_kind = 'rw1')
	expect_true(all(is.finite(on$BETA)))
	# kinds and grids the update does not cover keep the sampler's own updates
	expect_null(
		run(list(actor_move_every = 2), R = 0, dynamic_beta = 'intercept', dynamic_beta_kind = 'rw2')$actor_move
	)
	expect_null(
		run(list(actor_move_every = 2), R = 0, dynamic_beta = 'intercept', time_index = c(1, 2, 4, 5, 6))$actor_move
	)
})

test_that('ame() and lame() run the position steps alone when the model has no additive margin', {
	skip_on_cran()
	d = amb_panel(nt = 3L)
	for (rv in list(c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE))) {
		fl = suppressMessages(suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = 'binary',
			R = 1,
			rvar = rv[1],
			cvar = rv[2],
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = list(actor_move_every = 2)
		)))
		fa = suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = 'binary',
			R = 1,
			rvar = rv[1],
			cvar = rv[2],
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = list(actor_move_every = 2)
		)))
		for (f in list(fl, fa)) {
			expect_gt(unname(f$actor_move['uv_calls']), 0)
			expect_equal(unname(f$actor_move['level_calls']), if (any(rv)) 40 * 12 else 0)
			expect_true(all(is.finite(f$BETA)) && all(is.finite(f$VC)))
		}
		off = suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = 'binary',
			R = 1,
			rvar = rv[1],
			cvar = rv[2],
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = list(actor_move = FALSE)
		)))
		never = suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = 'binary',
			R = 1,
			rvar = rv[1],
			cvar = rv[2],
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE,
			prior = list(actor_move_every = 10^6)
		)))
		expect_null(off$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC)
	}
	# without any margin and without positions there is nothing to move
	expect_null(
		suppressMessages(suppressWarnings(ame(
			d$Y[[1]],
			Xdyad = d$X[[1]],
			family = 'binary',
			R = 0,
			rvar = FALSE,
			cvar = FALSE,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			gof = FALSE
		)))$actor_move
	)
	# symmetric, without the actor effect
	set.seed(4); n = 10; nm = sprintf('a%02d', 1:n)
	ys = matrix(
		rbinom(n * n, 1, .2),
		n,
		dimnames = list(nm, nm)
	); ys[lower.tri(ys)] = t(ys)[lower.tri(ys)]; diag(ys) = NA
	fs = suppressMessages(suppressWarnings(ame(
		ys,
		family = 'binary',
		R = 1,
		symmetric = TRUE,
		nvar = FALSE,
		burn = 40,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		prior = list(actor_move_every = 2)
	)))
	expect_equal(unname(fs$actor_move['level_calls']), 0)
	expect_gt(unname(fs$actor_move['uv_calls']), 0)
})
