skip_on_cran()

test_that('the bipartite context indexes the observed cells of each actor', {
	ctx = .lame_actor_move_context_bip
	set.seed(417); nA = 4L; nB = 5L; nt = 3L
	Y = array(rbinom(nA * nB * nt, 1, .4), c(nA, nB, nt))
	Y[2, 3, 1] = NA; Y[3, , 2] = NA; Y[, 5, 3] = NA # a missing cell, an actor absent from a period, another absent from another
	cx = ctx(Y)
	expect_equal(c(cx$nA, cx$nB, cx$N), c(nA, nB, nt))
	for (i in 1:nA) {
		s = cx$rows[[i]]; slab = matrix(Y[i, , ], nB, nt)
		expect_equal(s$idx, which(is.finite(slab)))
		expect_equal(s$sg, 2 * slab[s$idx] - 1)
	}
	for (j in 1:nB) {
		s = cx$cols[[j]]; slab = matrix(Y[, j, ], nA, nt)
		expect_equal(s$idx, which(is.finite(slab)))
		expect_equal(s$sg, 2 * slab[s$idx] - 1)
	}
	# the cells of the absent actors are not in the index sets
	expect_false(any(cx$rows[[3]]$idx %in% (nB + seq_len(nB))))
	expect_equal(length(cx$rows[[3]]$idx), 2 * nB - 1) # period 1 in full, period 2 absent, period 3 without column 5
	# a cross-section is a one-period panel, an actor with no observed cell has an empty factor
	y = Y[, , 1]; y[1, ] = NA
	c1 = ctx(y); expect_equal(c1$N, 1L); expect_length(c1$rows[[1]]$idx, 0L)
	# a value that is neither 0 nor 1 is unconstrained in the latent draw, so it carries no factor
	y2 = Y[, , 1]; y2[3, 2] = 2
	expect_false(2L %in% ctx(y2)$rows[[3]]$idx); expect_equal(length(ctx(y2)$rows[[3]]$idx), nB - 1L)
})

test_that('the bipartite level translation prior equals the shifted joint Gaussian prior', {
	set.seed(1103)
	for (nA in c(4L, 7L)) for (centred in c(TRUE, FALSE)) for (compensate in c(TRUE, FALSE)) {
		x = rnorm(nA, sd = .6); if (centred) x = x - mean(x)
		v = 1.4; beta = c(-.7, .3, 1.1); iV0 = crossprod(matrix(rnorm(9), 3)) + diag(.5, 3)
		for (i in c(1L, nA)) {
			joint = function(delta) {
				w = rep(0, nA); w[i] = 1
				if (compensate) w = w - 1 / nA
				be = beta; if (compensate) be[1] = be[1] + delta / nA
				# the levels' own prior on the whole space or on the zero-sum hyperplane (the same density there)
				-.5 * (sum((x + delta * w)^2) / v + drop(t(be) %*% iV0 %*% be))
			}
			pr = .lame_actor_move_level_prior_bip(i, x, v, iV0, beta, compensate)
			for (d in c(-2, -.5, .3, 1.7))
				expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
		}
	}
})

# the log density, up to a constant, of an intercept path under the coefficient sampler's prior:
# beta_1 ~ N(m0, v0), beta_t = rho beta_{t-1} + e_t, e_t ~ N(0, sd2)
bipmb_logpath = function(path, rho, sd2, m0, v0) {
	k = length(path)
	lp = -.5 * (path[1] - m0)^2 / v0
	if (k > 1) lp = lp - .5 * sum((path[-1] - rho * path[-k])^2) / sd2
	lp
}

test_that('the bipartite level translation prior with a moving intercept equals the shifted joint Gaussian prior', {
	set.seed(1104)
	N = 4L; v = 1.4
	for (nA in c(4L, 7L)) for (kind in c('path', 'static')) for (rho in c(1, .6)) {
		x = rnorm(nA, sd = .6); x = x - mean(x)
		ic = if (kind == 'path') list(
			kind = 'path',
			path = rnorm(N),
			rho = rho,
			sd2 = .35,
			m0 = .3,
			v0 = 4
		) else list(kind = 'static', value = .6, prec = .3)
		logic = function(sh) if (kind == 'path') bipmb_logpath(sh$path, rho, .35, .3, 4) else -.5 * .3 * sh$value^2
		for (i in c(1L, nA)) {
			joint = function(delta) {
				w = rep(-1 / nA, nA); w[i] = 1 - 1 / nA
				-.5 * sum((x + delta * w)^2) / v + logic(.lame_actor_move_icpt_shift(ic, NULL, delta / nA)$ic)
			}
			pr = .lame_actor_move_level_prior_bip(i, x, v, NULL, NULL, TRUE, ic)
			for (d in c(-2, -.5, .3, 1.7))
				expect_equal(joint(d) - joint(0), -(pr[['a2']] * d^2 + 2 * pr[['a1']] * d) / 2, tolerance = 1e-12)
		}
	}
	# a fit without an intercept to compensate and a fit with static coefficients are unchanged by the argument
	iV0 = crossprod(matrix(rnorm(9), 3)) + diag(.5, 3); beta = c(-.7, .3, 1.1); x = rnorm(5, sd = .5)
	expect_identical(
		.lame_actor_move_level_prior_bip(2L, x, 1.2, iV0, beta, TRUE, NULL),
		.lame_actor_move_level_prior_bip(2L, x, 1.2, iV0, beta, TRUE)
	)
	ic = list(kind = 'path', path = rnorm(3), rho = .5, sd2 = .3, m0 = 0, v0 = 3)
	expect_identical(
		.lame_actor_move_level_prior_bip(2L, x, 1.2, iV0, beta, FALSE, ic),
		c(a1 = x[2] / 1.2, a2 = 1 / 1.2)
	)
})

test_that('the bipartite visit keeps the predictor, the centring and the signs', {
	skip_on_cran()
	ctx = .lame_actor_move_context_bip; visit = .lame_actor_move_visit_bip
	set.seed(481516); nA = 5L; nB = 6L; nt = 3L
	Y = array(rbinom(nA * nB * nt, 1, .25), c(nA, nB, nt))
	Y[nA, , ] = 0; Y[, nB, ] = 0; Y[2, 3, 1] = NA; Y[1, , 2] = NA
	x = array(rnorm(nA * nB * nt), c(nA, nB, nt))
	cx = ctx(Y)
	for (dyn in c(FALSE, TRUE)) for (rank in c(0L, 2L)) for (centred in c(TRUE, FALSE)) {
		U = if (rank) matrix(rnorm(nA * rank, sd = .5), nA) else NULL; V = if (rank) matrix(rnorm(nB * rank, sd = .5), nB) else NULL
		G = if (rank) matrix(rnorm(rank * rank, sd = .7), rank) else NULL
		la = rnorm(nA, sd = .4); lb = rnorm(nB, sd = .4)
		if (centred) { la = la - mean(la); lb = lb - mean(lb) }
		dev_a = matrix(rnorm(nA * nt, sd = .2), nA); dev_a = sweep(dev_a, 2, colMeans(dev_a))
		dev_b = matrix(rnorm(nB * nt, sd = .2), nB); dev_b = sweep(dev_b, 2, colMeans(dev_b))
		am = if (dyn) matrix(la, nA, nt) + dev_a else NULL
		bm = if (dyn) matrix(lb, nB, nt) + dev_b else NULL
		beta = c(-1, .4)
		predictor = function(beta, A, B, U, V) { z = array(0, c(nA, nB, nt)); for (t in 1:nt) {
			z[, , t] = beta[1] + beta[2] * x[, , t] + outer(A[, t], B[, t], '+'); if (!is.null(U)) z[, , t] = z[, , t] + U %*% G %*% t(V) }; z }
		A = if (dyn) am else matrix(la, nA, nt); B = if (dyn) bm else matrix(lb, nB, nt)
		EZ = predictor(beta, A, B, U, V)
		iV0 = matrix(c(1.2, .2, .2, 1), 2); sa0 = sum(la); sb0 = sum(lb)
		for (k in 1:6) {
			o = visit(
				EZ,
				beta,
				la,
				lb,
				am,
				bm,
				U,
				V,
				G,
				1.3,
				.8,
				iV0,
				TRUE,
				cx,
				dev = if (dyn) list(rho = .6, sigma = .4),
				free_means = !centred
			)
			EZ = o$EZ; beta = o$beta; la = o$la; lb = o$lb; am = o$a_mat; bm = o$b_mat; U = o$U; V = o$V
			A = if (dyn) am else matrix(la, nA, nt); B = if (dyn) bm else matrix(lb, nB, nt)
			expect_lt(max(abs(EZ - predictor(beta, A, B, U, V))), 1e-10)
			# the translation sums to zero: the centring of the levels (and of every period's path) is kept; the levels that are not
			# centred have the ridge step, which moves their means (with the intercept) and nothing else
			if (centred) { expect_lt(abs(sum(la) - sa0), 1e-10); expect_lt(abs(sum(lb) - sb0), 1e-10)
				if (dyn) { expect_lt(
					max(abs(colSums(am) - sa0)),
					1e-10
				); expect_lt(max(abs(colSums(bm) - sb0)), 1e-10) } }
		}
		expect_gt(o$stats[['level_calls']], 0); expect_equal(o$stats[['stuck']], 0)
		if (rank) expect_gt(o$stats[['uv_calls']], 0) else expect_equal(o$stats[['uv_calls']], 0)
		# the deviation paths of actor pairs are drawn under dynamic effects only: nA %/% 2 + nB %/% 2 pairs per visit
		if (dyn) expect_equal(
			o$stats[['dev_calls']],
			nA %/% 2 + nB %/% 2
		) else expect_false('dev_calls' %in% names(o$stats))
	}
	# the ridge step leaves the predictor unchanged and the identified sum intercept + mean(a) + mean(b) too
	la = rnorm(nA, sd = .4); lb = rnorm(nB, sd = .4); beta = c(-1, .4)
	EZ = array(0, c(nA, nB, nt)); for (t in 1:nt) EZ[, , t] = beta[1] + beta[2] * x[, , t] + outer(la, lb, '+')
	cx0 = ctx(
		array(NA_real_, c(nA, nB, nt))
	) # no observed cell: the other blocks draw from their priors, the ridge from its conditional
	o = visit(
		EZ,
		beta,
		la,
		lb,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1.3,
		.8,
		matrix(c(1.2, .2, .2, 1), 2),
		TRUE,
		cx0,
		free_means = TRUE
	)
	expect_equal(o$stats[['ridge_calls']], 1)
	expect_lt(
		max(
			abs(
				o$EZ - array(
					sapply(1:nt, function(t) o$beta[1] + o$beta[2] * x[, , t] + outer(o$la, o$lb, '+')),
					c(nA, nB, nt)
				)
			)
		),
		1e-10
	)
	# a margin without effects is left alone, and without an intercept a level moves by itself
	o = visit(
		array(0, c(nA, nB, nt)),
		c(.2),
		rep(0, nA),
		rnorm(nB),
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		1,
		1,
		matrix(1),
		FALSE,
		cx,
		FALSE,
		TRUE
	)
	expect_equal(o$la, rep(0, nA)); expect_equal(o$beta, .2)
	expect_equal(o$stats[['level_calls']], nB)
})

test_that('the bipartite visit leaves the prior invariant: joint-distribution test', {
	skip_on_cran()
	# draw (theta, y) from the joint, apply one visit to theta given y; theta' must have the law of theta,
	# so every statistic has the same mean before and after (Geweke 2004), whatever the data drawn
	ctx = .lame_actor_move_context_bip; visit = .lame_actor_move_visit_bip
	run = function(nrep, nA, nB, nt, R, mode, dyn, seed, ar = c(.7, .5), margins = TRUE) {
		set.seed(seed)
		va = 1.3; vb = .8; g = 3
		x = array(rnorm(nA * nB * nt), c(nA, nB, nt))
		mask = array(TRUE, c(nA, nB, nt)); mask[sample(
			length(mask),
			floor(.2 * length(mask))
		)] = FALSE; mask[2, , nt] = FALSE; mask[, nB, 1] = FALSE
		Ha = qr.Q(qr(contr.helmert(nA))); Hb = qr.Q(qr(contr.helmert(nB)))
		intercept = mode != 'noint'; centred = mode == 'centred'; p = if (intercept) 2L else 1L
		pred = function(beta, A, B, U, V, G) { z = array(0, c(nA, nB, nt)); for (t in 1:nt) {
			z[, , t] = (if (intercept) beta[1] else 0) + beta[p] * x[, , t] + outer(A[, t], B[, t], '+')
			if (R > 0) z[, , t] = z[, , t] + U %*% G %*% t(V) }; z }
		sts = function(beta, la, lb, A, B, U, V, G) {
			eta = pred(beta, A, B, U, V, G)
			s = c(
				beta,
				beta^2,
				la,
				la^2,
				lb,
				lb^2,
				la[1] * lb[1],
				beta[1] * la[1],
				beta[p] * lb[2],
				A[1, nt],
				A[nA, 1]^2,
				B[2, nt] * A[1, 1],
				A[2, 2] * A[3, 1],
				B[1, nt]^2,
				sum(A[, nt]^2),
				sum(B[, 1]^2),
				A[2, nt] * B[3, nt],
				eta[1, 1, 1],
				eta[nA, nB, nt],
				eta[2, 2, 1]^2,
				sum(eta),
				sum(eta^2)
			)
			if (R > 0) s = c(
				s,
				U[1, 1],
				U[1, 1]^2,
				V[nB, R],
				V[2, 1]^2,
				U[1, 1] * V[1, 1],
				sum(U^2),
				sum(V^2),
				G[1, 1] * U[2, 1],
				sum((U %*% G)^2)
			)
			s }
		old = new = NULL
		for (r in seq_len(nrep)) {
			beta = rnorm(p, 0, sqrt(g))
			la = if (centred) drop(sqrt(va) * Ha %*% rnorm(nA - 1)) else rnorm(nA, 0, sqrt(va))
			lb = if (centred) drop(sqrt(vb) * Hb %*% rnorm(nB - 1)) else rnorm(nB, 0, sqrt(vb))
			# without an additive margin only the positions are visited
			if (!margins) { la = numeric(nA); lb = numeric(nB) }
			if (dyn) {
				# independent stationary AR(1) deviation paths per actor, conditioned to sum to zero in every period (the actors share
				# one covariance, so the conditional is the per-period centring)
				ar1 = function(n) { m = matrix(0, n, nt); m[, 1] = rnorm(n, 0, ar[2] / sqrt(1 - ar[1]^2))
					for (t in seq_len(nt)[-1]) m[, t] = ar[1] * m[, t - 1] + rnorm(n, 0, ar[2]); sweep(m, 2, colMeans(m)) }
				da = ar1(nA); db = ar1(nB)
				am = matrix(la, nA, nt) + da; bm = matrix(lb, nB, nt) + db } else am = bm = NULL
			A = if (dyn) am else matrix(la, nA, nt); B = if (dyn) bm else matrix(lb, nB, nt)
			if (R > 0) { U = matrix(rnorm(nA * R), nA); V = matrix(rnorm(nB * R), nB); G = matrix(rnorm(R * R), R) } else U = V = G = NULL
			eta = pred(beta, A, B, U, V, G)
			Y = 1 * (eta + array(rnorm(length(eta)), dim(eta)) > 0); Y[!mask] = NA
			cx = ctx(Y)
			st = visit(
				eta,
				beta,
				la,
				lb,
				am,
				bm,
				U,
				V,
				G,
				va,
				vb,
				diag(p) / g,
				intercept,
				cx,
				move_a = margins,
				move_b = margins,
				dev = if (dyn) list(rho = ar[1], sigma = ar[2]),
				free_means = (mode == 'free')
			)
			A2 = if (dyn) st$a_mat else matrix(st$la, nA, nt); B2 = if (dyn) st$b_mat else matrix(st$lb, nB, nt)
			old = rbind(old, sts(beta, la, lb, A, B, U, V, G)); new = rbind(
				new,
				sts(st$beta, st$la, st$lb, A2, B2, st$U, st$V, G)
			)
		}
		d = new - old; keep = apply(d, 2, sd) > 0
		z = colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrep))
		max(abs(z)) }
	expect_lt(run(3000L, 3L, 4L, 2L, 0L, 'centred', FALSE, 1L), 4.5)
	expect_lt(run(3000L, 4L, 5L, 3L, 1L, 'centred', TRUE, 2L), 4.5)
	expect_lt(run(3000L, 3L, 4L, 2L, 1L, 'free', FALSE, 3L), 4.5)
	expect_lt(run(3000L, 3L, 4L, 2L, 0L, 'noint', FALSE, 4L), 4.5)
	expect_lt(run(3000L, 4L, 5L, 2L, 1L, 'centred', FALSE, 5L, margins = FALSE), 4.5)
	expect_lt(run(3000L, 4L, 5L, 2L, 2L, 'noint', FALSE, 6L, margins = FALSE), 4.5)
})

test_that('the bipartite visit with a moving or static dynamic_beta intercept leaves the prior invariant: joint-distribution test', {
	skip_on_cran()
	# as above, with the intercept a path (beta_1 ~ N(m0, v0), beta_t = rho beta_{t-1} + e_t, e_t ~ N(0, sd2)) or a static
	# coefficient next to moving ones (N(0, g)); the statistics include the path
	ctx = .lame_actor_move_context_bip; visit = .lame_actor_move_visit_bip
	run = function(nrep, nA, nB, nt, R, kind, dyn, seed, rho_b = .6, move_b = TRUE, ic_edit = NULL) {
		set.seed(seed)
		va = 1.3; vb = .8; g = 3; sd2 = .35; m0 = .3; v0 = 4; prec = 1 / g
		x = array(rnorm(nA * nB * nt), c(nA, nB, nt))
		mask = array(TRUE, c(nA, nB, nt)); mask[sample(
			length(mask),
			floor(.2 * length(mask))
		)] = FALSE; mask[2, , nt] = FALSE; mask[, nB, 1] = FALSE
		Ha = qr.Q(qr(contr.helmert(nA))); Hb = qr.Q(qr(contr.helmert(nB)))
		pred = function(path, b2, A, B, U, V, G) { z = array(0, c(nA, nB, nt)); for (t in 1:nt) {
			z[, , t] = path[t] + b2 * x[, , t] + outer(A[, t], B[, t], '+')
			if (R > 0) z[, , t] = z[, , t] + U %*% G %*% t(V) }; z }
		sts = function(path, b2, la, lb, A, B, U, V, G) {
			eta = pred(path, b2, A, B, U, V, G)
			s = c(
				path,
				path^2,
				b2,
				b2^2,
				la,
				la^2,
				lb,
				lb^2,
				la[1] * lb[1],
				path[1] * la[1],
				path[nt] * lb[2],
				mean(path) * b2,
				b2 * lb[2],
				A[1, nt],
				A[nA, 1]^2,
				B[2, nt] * A[1, 1],
				A[2, 2] * A[3, 1],
				B[1, nt]^2,
				sum(A[, nt]^2),
				sum(B[, 1]^2),
				A[2, nt] * B[3, nt],
				eta[1, 1, 1],
				eta[nA, nB, nt],
				eta[2, 2, 1]^2,
				sum(eta),
				sum(eta^2)
			)
			if (R > 0) s = c(
				s,
				U[1, 1],
				U[1, 1]^2,
				V[nB, R],
				V[2, 1]^2,
				U[1, 1] * V[1, 1],
				sum(U^2),
				sum(V^2),
				G[1, 1] * U[2, 1],
				sum((U %*% G)^2)
			)
			s }
		old = new = vector('list', nrep)
		for (r in seq_len(nrep)) {
			path = if (kind == 'static') rep(
				rnorm(1, 0, sqrt(1 / prec)),
				nt
			) else { p_ = numeric(nt); p_[1] = rnorm(1, m0, sqrt(v0))
				for (t in seq_len(nt)[-1]) p_[t] = rho_b * p_[t - 1] + rnorm(1, 0, sqrt(sd2)); p_ }
			b2 = rnorm(1, 0, sqrt(g))
			la = drop(sqrt(va) * Ha %*% rnorm(nA - 1)); lb = drop(sqrt(vb) * Hb %*% rnorm(nB - 1))
			if (dyn) {
				ar1 = function(n) { m = matrix(0, n, nt); m[, 1] = rnorm(n, 0, .5 / sqrt(1 - .7^2))
					for (t in seq_len(nt)[-1]) m[, t] = .7 * m[, t - 1] + rnorm(n, 0, .5); sweep(m, 2, colMeans(m)) }
				am = matrix(la, nA, nt) + ar1(nA); bm = matrix(lb, nB, nt) + ar1(nB) } else am = bm = NULL
			A = if (dyn) am else matrix(la, nA, nt); B = if (dyn) bm else matrix(lb, nB, nt)
			if (R > 0) { U = matrix(rnorm(nA * R), nA); V = matrix(rnorm(nB * R), nB); G = matrix(rnorm(R * R), R) } else U = V = G = NULL
			eta = pred(path, b2, A, B, U, V, G)
			Y = 1 * (eta + array(rnorm(length(eta)), dim(eta)) > 0); Y[!mask] = NA
			cx = ctx(Y)
			ic = if (kind == 'static') list(kind = 'static', value = path[1], prec = prec) else list(
				kind = 'path',
				path = path,
				rho = rho_b,
				sd2 = sd2,
				m0 = m0,
				v0 = v0
			)
			if (!is.null(ic_edit)) ic = ic_edit(ic)
			st = visit(
				eta,
				c(mean(path), b2),
				la,
				lb,
				am,
				bm,
				U,
				V,
				G,
				va,
				vb,
				diag(2) / g,
				TRUE,
				cx,
				move_b = move_b,
				dev = if (dyn) list(rho = .7, sigma = .5),
				ic = ic
			)
			path2 = if (kind == 'static') rep(st$ic$value, nt) else st$ic$path
			A2 = if (dyn) st$a_mat else matrix(st$la, nA, nt); B2 = if (dyn) st$b_mat else matrix(st$lb, nB, nt)
			old[[r]] = sts(path, b2, la, lb, A, B, U, V, G); new[[r]] = sts(
				path2,
				b2,
				st$la,
				st$lb,
				A2,
				B2,
				st$U,
				st$V,
				G
			)
		}
		d = do.call(rbind, new) - do.call(rbind, old); keep = apply(d, 2, sd) > 0
		max(abs(colMeans(d)[keep] / (apply(d, 2, sd)[keep] / sqrt(nrep)))) }
	expect_lt(run(3000L, 3L, 4L, 3L, 0L, 'rw1', FALSE, 11L, rho_b = 1), 4.5)
	expect_lt(run(3000L, 4L, 5L, 3L, 1L, 'ar1', TRUE, 12L), 4.5)
	expect_lt(run(3000L, 3L, 4L, 2L, 1L, 'static', FALSE, 13L), 4.5)
	expect_lt(run(3000L, 4L, 5L, 3L, 0L, 'ar1', FALSE, 14L, move_b = FALSE), 4.5)
	expect_lt(run(3000L, 4L, 5L, 4L, 1L, 'rw1', TRUE, 15L, rho_b = 1), 4.5)
	# the test has power: a wrong prior of the first state or of the transitions is detected
	expect_gt(run(6000L, 3L, 4L, 3L, 0L, 'ar1', FALSE, 16L, ic_edit = function(ic) { ic$v0 = ic$v0 / 8; ic }), 4.5)
	expect_gt(run(6000L, 3L, 4L, 3L, 0L, 'ar1', FALSE, 17L, ic_edit = function(ic) { ic$sd2 = ic$sd2 / 8; ic }), 4.5)
})

test_that('the bipartite visit leaves the exact conditional posterior invariant', {
	skip_on_cran()
	ctx = .lame_actor_move_context_bip; visit = .lame_actor_move_visit_bip
	set.seed(2718); nA = 3L; nB = 4L; nt = 2L
	Y = array(rbinom(nA * nB * nt, 1, .3), c(nA, nB, nt)); Y[nA, , ] = 0; Y[, nB, ] = 0; Y[2, 3, 1] = NA
	x = array(rnorm(nA * nB * nt), c(nA, nB, nt)); b2 = .5
	va = 1.2; vb = .8; iV0 = matrix(c(1.5, .3, .3, 1.1), 2)
	cx = ctx(Y); sg = 2 * Y - 1; ob = is.finite(Y)
	Ha = qr.Q(qr(contr.helmert(nA))); Hb = qr.Q(qr(contr.helmert(nB)))
	pred = function(b0, a, b) { z = array(0, c(nA, nB, nt)); for (t in 1:nt) z[, , t] = b0 + b2 * x[, , t] + outer(a, b, '+'); z }
	loglik = function(eta) sum(pnorm(sg[ob] * eta[ob], log.p = TRUE))
	beta = c(-1, b2); a = c(.2, -.3, .1); a = a - mean(a); b = c(-.1, .2, .1, -.2); b = b - mean(b)
	EZ = pred(beta[1], a, b)
	nv = 15000L; keep = matrix(NA_real_, nv, 1 + nA + nB)
	for (k in seq_len(nv)) {
		# the second coefficient is held at b2 by the test: only the intercept and the levels move
		o = visit(EZ, beta, a, b, NULL, NULL, NULL, NULL, NULL, va, vb, iV0, TRUE, cx)
		EZ = o$EZ; beta = o$beta; a = o$la; b = o$lb; keep[k, ] = c(beta[1], a, b)
	}
	# reference: importance sampling from the exact Gaussian prior on the centred subspaces, the intercept
	# conditional on the second coefficient, weighted by the exact likelihood
	nmc = 60000L; m0 = -iV0[1, 2] * b2 / iV0[1, 1]; s0 = 1 / sqrt(iV0[1, 1])
	dr = matrix(NA_real_, nmc, 1 + nA + nB); lw = numeric(nmc)
	for (m in seq_len(nmc)) { aa = drop(
		sqrt(va) * Ha %*% rnorm(nA - 1)
	); bb = drop(sqrt(vb) * Hb %*% rnorm(nB - 1)); b0 = rnorm(1, m0, s0)
		dr[m, ] = c(b0, aa, bb); lw[m] = loglik(pred(b0, aa, bb)) }
	w = exp(lw - max(lw)); w = w / sum(w); ref = colSums(dr * w)
	se_is = sqrt(colSums(sweep(dr, 2, ref)^2 * w^2))
	ess = pmax(apply(keep, 2, function(v) coda::effectiveSize(v)), 50)
	se_chain = apply(keep, 2, sd) / sqrt(ess)
	z = (colMeans(keep) - ref) / sqrt(se_is^2 + se_chain^2)
	expect_lt(max(abs(z)), 4)
	expect_lt(max(abs(apply(keep, 2, sd) / sqrt(colSums(sweep(dr, 2, ref)^2 * w)) - 1)), .15)
})

test_that('the bipartite visit with a moving or static dynamic_beta intercept leaves the exact conditional posterior invariant', {
	skip_on_cran()
	# a tiny model; the visit moves the intercept path by a constant, so the path's shape is held and the constant is the unknown
	ctx = .lame_actor_move_context_bip; visit = .lame_actor_move_visit_bip
	set.seed(2719); nA = 3L; nB = 4L; nt = 3L
	Y = array(rbinom(nA * nB * nt, 1, .3), c(nA, nB, nt)); Y[nA, , ] = 0; Y[, nB, ] = 0; Y[2, 3, 1] = NA
	x = array(rnorm(nA * nB * nt), c(nA, nB, nt)); b2 = .5
	va = 1.2; vb = .8; shp = c(.5, -.2, .3)
	cx = ctx(Y); sg = 2 * Y - 1; ob = is.finite(Y)
	Ha = qr.Q(qr(contr.helmert(nA))); Hb = qr.Q(qr(contr.helmert(nB)))
	pred = function(path, a, b) { z = array(0, c(nA, nB, nt)); for (t in 1:nt) z[, , t] = path[t] + b2 * x[, , t] + outer(a, b, '+'); z }
	loglik = function(eta) sum(pnorm(sg[ob] * eta[ob], log.p = TRUE))
	for (kind in c('rw1', 'ar1', 'static')) {
		rho_b = if (kind == 'ar1') .6 else 1; sd2 = .3; m0 = .2; v0 = 4; prec = .7
		mk_ic = function(b0) if (kind == 'static') list(kind = 'static', value = b0, prec = prec) else list(
			kind = 'path',
			path = shp + b0,
			rho = rho_b,
			sd2 = sd2,
			m0 = m0,
			v0 = v0
		)
		shape = function(b0) if (kind == 'static') rep(b0, nt) else shp + b0
		# the prior of the constant: Gaussian, read off the log prior of the path by differences
		lp_b0 = function(b0) if (kind == 'static') -.5 * prec * b0^2 else bipmb_logpath(shp + b0, rho_b, sd2, m0, v0)
		q = -(lp_b0(1) + lp_b0(-1) - 2 * lp_b0(0)); l = -(lp_b0(1) - lp_b0(-1)) / 2; mu0 = -l / q
		a = c(.2, -.3, .1); a = a - mean(a); b = c(-.1, .2, .1, -.2); b = b - mean(b); b0 = -1
		EZ = pred(shape(b0), a, b)
		nv = 15000L; keep = matrix(NA_real_, nv, 1 + nA + nB)
		for (k in seq_len(nv)) {
			o = visit(
				EZ,
				c(mean(shape(b0)), b2),
				a,
				b,
				NULL,
				NULL,
				NULL,
				NULL,
				NULL,
				va,
				vb,
				diag(2),
				TRUE,
				cx,
				ic = mk_ic(b0)
			)
			EZ = o$EZ; a = o$la; b = o$lb
			b0 = if (kind == 'static') o$ic$value else o$ic$path[1] - shp[1]
			keep[k, ] = c(b0, a, b)
		}
		nmc = 60000L
		dr = matrix(NA_real_, nmc, 1 + nA + nB); lw = numeric(nmc)
		for (m in seq_len(nmc)) { aa = drop(
			sqrt(va) * Ha %*% rnorm(nA - 1)
		); bb = drop(sqrt(vb) * Hb %*% rnorm(nB - 1)); bb0 = rnorm(1, mu0, 1 / sqrt(q))
			dr[m, ] = c(bb0, aa, bb); lw[m] = loglik(pred(shape(bb0), aa, bb)) }
		w = exp(lw - max(lw)); w = w / sum(w); ref = colSums(dr * w)
		se_is = sqrt(colSums(sweep(dr, 2, ref)^2 * w^2))
		ess = pmax(apply(keep, 2, function(v) coda::effectiveSize(v)), 50)
		se_chain = apply(keep, 2, sd) / sqrt(ess)
		z = (colMeans(keep) - ref) / sqrt(se_is^2 + se_chain^2)
		expect_lt(max(abs(z)), 4)
		expect_lt(max(abs(apply(keep, 2, sd) / sqrt(colSums(sweep(dr, 2, ref)^2 * w)) - 1)), .15)
	}
})

test_that('ame() runs the bipartite actor update as a prior switch and records its work', {
	set.seed(7305); nA = 9L; nB = 11L
	y = matrix(
		rbinom(nA * nB, 1, .12),
		nA,
		nB,
		dimnames = list(sprintf('r%02d', 1:nA), sprintf('c%02d', 1:nB))
	); y[1:2, ] = 0; y[, 1:3] = 0
	x = array(rnorm(nA * nB), c(nA, nB, 1), dimnames = list(rownames(y), colnames(y), 'x'))
	# ame() warns about prior names it does not know until it lists the actor-move switches; nothing else may warn
	run = function(prior, R, ...) suppressMessages(withCallingHandlers(
		ame(
			y,
			Xdyad = x,
			mode = 'bipartite',
			family = 'binary',
			R_row = R,
			R_col = R,
			burn = 20,
			nscan = 60,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = prior,
			...
		),
		warning = function(w) if (grepl('Unknown.*prior', conditionMessage(w))) invokeRestart('muffleWarning')
	))
	for (R in c(0L, 1L)) {
		expect_no_warning({ on = run(list(actor_move_every = 2), R) })
		expect_equal(unname(on$actor_move['visits']), 40)
		expect_equal(unname(on$actor_move['level_calls']), 40 * (nA + nB))
		expect_equal(unname(on$actor_move['uv_calls'] > 0), R > 0)
		expect_equal(unname(on$actor_move['uv_calls']), if (R > 0) 40 * (nA + nB) else 0)
		expect_equal(unname(on$actor_move['stuck']), 0)
		# the free (not centred) levels also share a ridge with the intercept, drawn once per visit
		expect_equal(unname(on$actor_move['ridge_calls']), 40)
		expect_true(all(is.finite(on$BETA)) && all(is.finite(on$VC)))
		off = run(list(actor_move = FALSE), R)
		expect_null(off$actor_move)
		# a visit interval beyond the run never visits: the draws are those of the switch off
		never = run(list(actor_move_every = 1000), R)
		expect_null(never$actor_move)
		expect_identical(never$BETA, off$BETA); expect_identical(never$VC, off$VC); expect_identical(never$APM, off$APM)
		expect_false(isTRUE(all.equal(on$BETA, off$BETA)))
	}
	expect_error(run(list(actor_move = NA), 0L), 'actor_move must be TRUE or FALSE', fixed = TRUE)
	expect_error(run(list(actor_move_every = 0), 0L), 'actor_move_every must be a positive whole number', fixed = TRUE)
	# other families run no level or position step; with an intercept they draw along the ridge, which has no likelihood term
	yn = y + matrix(rnorm(nA * nB), nA)
	nm = suppressWarnings(
		ame(
			yn,
			mode = 'bipartite',
			family = 'normal',
			burn = 10,
			nscan = 20,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_identical(names(nm$actor_move), c('visits', 'ridge_calls'))
	expect_equal(unname(nm$actor_move), c(30, 30))
	expect_null(
		suppressWarnings(ame(
			yn,
			mode = 'bipartite',
			family = 'normal',
			burn = 10,
			nscan = 20,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			prior = list(actor_move = FALSE)
		))$actor_move
	)
})
