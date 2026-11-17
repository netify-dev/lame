# R functions only the tests use: reference forms of the ALS evidence, the variational objective and
# the position sweeps of the dynamic empirical-Bayes fit, which the tests compare the estimators'
# compiled or structured paths with

# the evidence of a fitted object (the value stored in fit$eb$log_marginal)
.ae_eb_evidence = function(object) {
	if (identical(object$eb$estimator, "variational")) return(object$eb$laplace)
	lp = .ae_eb_laplace_state(object, inverse = FALSE)
	jp = lp$jp; ctx = jp$ctx
	iz = setdiff(seq_len(nrow(jp$Hf)), c(ctx$i_th, ctx$i_g))
	.ae_eb_evidence_core(
		lp$prep,
		lp$ab,
		lp$fit,
		lp$R,
		lp$eta0,
		lp$link,
		2 * sum(log(diag(chol(jp$Hf[iz, iz, drop = FALSE])))),
		2 * sum(log(diag(chol(jp$H0[ctx$i_s, ctx$i_s, drop = FALSE]))))
	)
}

# the evidence lower bound of a variational fit, recomputed from the object at
# its stored state (the value stored in fit$eb$variational$elbo)
.ae_mf_elbo_object = function(object) {
	ms = .ae_mf_object_state(object)
	link = .ae_default(object$link, "probit")
	cx = .ae_mf_context(
		ms$prep,
		ms$ab,
		ms$R,
		link,
		pair = if (!is.null(object$eb$pair)) list(rho = object$eb$pair$rho, quad = object$eb$pair$quad, update = FALSE)
	)
	ev = .ae_mf_evaluate(cx, ms$S)
	.ae_mf_elbo(cx, ms$S, ev)
}

# q of one side's positions given the other's. W and Rr are per-period lists of the (n x m)
# weights and working residuals (zero off the observed cells), Vh and Om the other side's
# (m x R) means and (R x R x m) covariances. returns the means (n x R x Tt), the diagonal
# covariance blocks of each actor's stacked path (R x R x n x Tt), the lag-one blocks
# (R x R x n x (Tt - 1)) and the KL pieces sum_i log|P_i| (ldp), sum_i tr(P0 Sigma_i) (trp)
# and n log|P0| (ldp0)
.lde_pos_side = function(W, Rr, Vh, Om, s2, rho, n, R, Tt, lam = NULL) {
	P0 = kronecker(.lde_pos_q(Tt, rho) / s2, diag(R))
	Lam = array(0, dim = c(R, R, n, Tt))
	H = array(0, dim = c(n, R, Tt))
	for (t in seq_len(Tt)) {
		for (k in seq_len(R)) for (l in k:R) {
			ev = Vh[[t]][, k] * Vh[[t]][, l] + Om[[t]][k, l, ]
			v = as.vector(W[[t]] %*% ev)
			Lam[k, l, , t] = v
			Lam[l, k, , t] = v
		}
		H[, , t] = (W[[t]] * Rr[[t]]) %*% Vh[[t]]
	}
	mean_ = array(0, dim = c(n, R, Tt))
	sig = array(0, dim = c(R, R, n, Tt))
	cl = array(0, dim = c(R, R, n, max(Tt - 1L, 1L)))
	ldp = 0
	trp = 0
	ldp0 = n * as.numeric(determinant(P0)$modulus)
	for (i in seq_len(n)) {
		# the prior precision of this actor's path (its transition weights, if the transition is a scale mixture)
		P0i = if (is.null(lam)) P0 else kronecker(.lde_pos_qw(Tt, rho, lam[i, ]) / s2, diag(R))
		Pi = P0i
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			Pi[ix, ix] = Pi[ix, ix] + Lam[, , i, t]
		}
		ch = .lde_chol(Pi)
		ldp = ldp + 2 * sum(log(diag(ch)))
		Sg = chol2inv(ch)
		trp = trp + sum(P0i * Sg)
		mean_[i, , ] = matrix(as.vector(Sg %*% as.vector(H[i, , ])), R, Tt)
		for (t in seq_len(Tt)) {
			ix = (t - 1L) * R + seq_len(R)
			sig[, , i, t] = Sg[ix, ix]
			if (t < Tt) cl[, , i, t] = Sg[ix, ix + R]
		}
	}
	list(mean = mean_, sig = sig, cl = cl, ldp = ldp, trp = trp, ldp0 = ldp0, P0 = P0)
}

# R reference of the sweep, used by the tests
.lde_pos_sweep_ref = function(ctx, pos, hp, wk, eta_add) {
	if (isTRUE(ctx$sym)) return(.lde_pos_sweep_sym(ctx, pos, hp, wk, eta_add))
	Tt = ctx$Tt; na = ctx$na; nb = ctx$nb; R = ctx$R
	obs = ctx$obs
	W = lapply(seq_len(Tt), function(t) {
		w = wk$w[, , t]
		w[!obs[, , t] | !is.finite(w)] = 0
		w
	})
	Rr = lapply(seq_len(Tt), function(t) {
		r = wk$z[, , t] - eta_add[, , t]
		r[!obs[, , t] | !is.finite(r)] = 0
		r
	})
	Vh = lapply(seq_len(Tt), function(t) matrix(pos$qv$mean[, , t], nb, R))
	Om = lapply(seq_len(Tt), function(t) array(pos$qv$sig[, , , t], c(R, R, nb)))
	qu = .lde_pos_side(W, Rr, Vh, Om, hp$s2, hp$rho, na, R, Tt, lam = pos$qu$lam)
	Uh = lapply(seq_len(Tt), function(t) matrix(qu$mean[, , t], na, R))
	Omu = lapply(seq_len(Tt), function(t) array(qu$sig[, , , t], c(R, R, na)))
	qv = .lde_pos_side(lapply(W, t), lapply(Rr, t), Uh, Omu, hp$s2, hp$rho, nb, R, Tt, lam = pos$qv$lam)
	qu$lam = pos$qu$lam
	qv$lam = pos$qv$lam
	list(qu = qu, qv = qv)
}

# R reference, one outer product per pair and period (used by the tests)
.lde_pos_prodvar_ref = function(pos, R, Tt) {
	Um = pos$qu$mean; Vm = pos$qv$mean
	na = dim(Um)[1L]; nb = dim(Vm)[1L]
	out = array(0, dim = c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		for (k in seq_len(R)) for (l in seq_len(R)) {
			su = pos$qu$sig[k, l, , t]; sv = pos$qv$sig[k, l, , t]
			out[, , t] = out[, , t] + outer(su, sv) + outer(Um[, k, t] * Um[, l, t], sv) +
				outer(su, Vm[, k, t] * Vm[, l, t])
		}
	}
	out
}

# expected working problem of a binary (probit or logit) cell whose predictor is N(eta, vnoise),
# vnoise the variance of the product under q, as in .ae_mf_cells(): E log p and its derivative in
# eta, E psi, and the expected Fisher weight with its first two derivatives (for the Laplace
# correction), all by Gauss-Hermite quadrature on the ladder above. the working response is
# eta + E psi / E w (the increment capped as in .lde_work()); vnoise = 0 gives the plain problem.
#   dev:      -2 E log p, the deviance of the expected log likelihood
#   dev_marg: -2 log E p = -2 log Phi(s c eta), c = (1 + vnoise)^-1/2, the marginal deviance used by
#             the profile over the position scale, since E log p is optimised by the collapse
#             (by quadrature for the logit)
.lde_work_gh_ref = function(
	y,
	eta,
	vnoise,
	K = 9L,
	delta_cap = formals(.lde_work)$delta_cap,
	wfloor = formals(.lde_work)$wfloor,
	block = 20000L,
	per_sd = .lde_gh_per_sd,
	link = "probit"
) {
	logit = identical(link, "logit")
	ix = which(is.finite(y))
	nobs = length(ix)
	rules = .lde_gh_rules(K)
	es = ew = eh = ehh = el = dmarg = numeric(nobs)
	vv = pmax(vnoise[ix], 0)
	# the rule of every cell: the smallest rung with at least per_sd * sd nodes, the last one beyond it
	lev = pmin(findInterval(per_sd * sqrt(vv), rules$K, left.open = TRUE) + 1L, length(rules$K))
	lev[is.na(lev)] = 1L
	for (r in seq_along(rules$K)) {
		cells = which(lev == r)
		nd = list(z = rules$z[[r]], w = rules$w[[r]])
		# cells in blocks so the cells x nodes matrices stay small; the node terms are those of
		# .lde_probit_raw(), with the outcome picked by index
		for (b0 in seq.int(1L, max(length(cells), 1L), by = block)) {
			if (!length(cells)) break
			bi = cells[b0:min(b0 + block - 1L, length(cells))]
			yy = y[ix[bi]]
			m = eta[ix[bi]]
			v = vv[bi]
			e = outer(sqrt(v), nd$z) + m
			one = yy == 1
			if (logit) {
				# the score of an outcome of 1 is p0 and minus that of an outcome of 0 is p1 (the roles of lam1 and
				# lam0), so the weight is p0 p1 = lam1 lam0; w' = w (p0 - p1), w'' = w' (p0 - p1) - 2 w^2
				lp1 = stats::plogis(e, log.p = TRUE)
				lp0 = stats::plogis(e, lower.tail = FALSE, log.p = TRUE)
				lam1 = stats::plogis(e, lower.tail = FALSE)
				lam0 = stats::plogis(e)
				wr = lam1 * lam0
				mm = lam1 - lam0
				hpr = wr * mm
				hppr = hpr * mm - 2 * wr * wr
			} else {
				lp1 = stats::pnorm(e, log.p = TRUE)
				lp0 = stats::pnorm(e, lower.tail = FALSE, log.p = TRUE)
				ld = -0.5 * e * e - 0.91893853320467267
				lam1 = exp(ld - lp1)
				lam0 = exp(ld - lp0)
				wr = lam1 * lam0
				mm = lam0 - lam1 - 2 * e
				hpr = wr * mm
				hppr = hpr * mm + wr * (lam0 * lam0 + lam1 * lam1 + e * (lam1 - lam0) - 2)
			}
			sc = lam0
			sc[one, ] = lam1[one, ]
			lg = lp0
			lg[one, ] = lp1[one, ]
			es[bi] = as.vector(sc %*% nd$w) * (2 * one - 1)
			ew[bi] = as.vector(pmax(wr, wfloor) %*% nd$w)
			eh[bi] = as.vector(hpr %*% nd$w)
			ehh[bi] = as.vector(hppr %*% nd$w)
			el[bi] = as.vector(lg %*% nd$w)
			# marginal deviance: the noise integrated out of the likelihood, Phi(c eta) (by quadrature for the logit)
			if (logit) {
				pm = lam1
				pm[one, ] = lam0[one, ]
				dmarg[bi] = log(as.vector(pm %*% nd$w))
			} else {
				dmarg[bi] = .lde_probit_raw(yy, m / sqrt(1 + v))$ell
			}
		}
	}
	delta = es / ew
	big = abs(delta) > delta_cap
	if (any(big)) {
		delta[big] = sign(delta[big]) * delta_cap
		ew[big] = abs(es[big]) / delta_cap
	}
	z = w = score = hprime = hpp = array(0, dim(y))
	z[ix] = eta[ix] + delta
	w[ix] = ew
	score[ix] = es
	hprime[ix] = ifelse(is.finite(eh), eh, 0)
	hpp[ix] = ifelse(is.finite(ehh), ehh, 0)
	list(z = z, w = w, score = score, hprime = hprime, hpp = hpp, dev = -2 * sum(el), dev_marg = -2 * sum(dmarg))
}

# the full information matrix built densely, in the order (s_1, ..., s_Tt, z)
.lde_dense_matrix = function(A, b, C, Hzz) {
	Tt = length(A); qs = nrow(A[[1L]]); r = ncol(Hzz)
	n = Tt * qs + r
	H = matrix(0, n, n)
	blk = function(t) (t - 1L) * qs + seq_len(qs)
	zi = Tt * qs + seq_len(r)
	for (t in seq_len(Tt)) {
		H[blk(t), blk(t)] = A[[t]]
		if (t < Tt) {
			H[blk(t), blk(t + 1L)] = diag(-b, qs)
			H[blk(t + 1L), blk(t)] = diag(-b, qs)
		}
		H[blk(t), zi] = C[[t]]
		H[zi, blk(t)] = t(C[[t]])
	}
	H[zi, zi] = Hzz
	H
}

# one side: the precision of every actor is the prior plus a block diagonal, as in the sweeps; the dense
# inverses are kept as the reference
.lde_pos_side_ref = function(n, R, Tt, hp, scale) {
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
