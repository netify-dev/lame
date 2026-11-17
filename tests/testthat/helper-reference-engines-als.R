# R references of the compiled kernels of the ALS empirical-Bayes estimators (src/ae_kernels.cpp,
# src/ae_pln_kernels.cpp, src/dynamic_eb_*.cpp). the tests compare each kernel with its reference, and
# local_als_r_engine() runs whole fits on the references in place of the kernels

.ae_pair_cache = new.env(parent = emptyenv())

# 64-point Gauss-Legendre rule on (-1, 1)
.ae_pair_gl = function() {
	hit = .ae_pair_cache[["gl"]]
	if (!is.null(hit)) return(hit)
	n = 64L
	i = seq_len(n - 1L)
	J = matrix(0, n, n)
	J[cbind(i, i + 1L)] = J[cbind(i + 1L, i)] = i / sqrt(4 * i^2 - 1)
	e = eigen(J, symmetric = TRUE)
	out = list(x = e$values, w = 2 * e$vectors[1L, ]^2)
	.ae_pair_cache[["gl"]] = out
	out
}

# log of the bivariate normal cdf Phi2(a, b; r). pbivnorm is inaccurate below a probability of 1e-6, so small
# probabilities are integrated in the log domain by Gauss-Legendre over (a - 14, a), a the smaller limit
.ae_pair_logphi2_r = function(a, b, r) {
	p = pbivnorm::pbivnorm(a, b, r)
	lp = log(pmax(p, 1e-300))
	bad = which(!is.finite(p) | p < 1e-6)
	if (length(bad)) {
		gl = .ae_pair_gl()
		# the integrand's mass sits below the smaller limit, so the window is placed under it
		ab = pmin(a[bad], b[bad]); bb = pmax(a[bad], b[bad]); rb = r[bad]
		w = 14
		xs = ab - outer(rep(1, length(ab)), w * (1 - (gl$x + 1) / 2))
		lw = matrix(log(gl$w * w / 2), length(ab), length(gl$x), byrow = TRUE)
		lf = stats::dnorm(xs, log = TRUE) + stats::pnorm((bb - rb * xs) / sqrt(1 - rb^2), log.p = TRUE) + lw
		mx = apply(lf, 1L, max)
		lp[bad] = mx + log(rowSums(exp(lf - mx)))
	}
	lp
}

# log P, P = Phi2(a, b; r), and its derivatives: gradient (La, Lb), Hessian (Laa, Lbb, Lab),
# score of r (Lr = phi2 / P) and its derivative (Lrr)
.ae_pair_logp_r = function(a, b, r) {
	kap = sqrt(1 - r^2)
	lP = .ae_pair_logphi2_r(a, b, r)
	ua = (a - r * b) / kap; ub = (b - r * a) / kap
	la = stats::dnorm(a, log = TRUE); lb = stats::dnorm(b, log = TRUE)
	lna = stats::pnorm(ua, log.p = TRUE); lnb = stats::pnorm(ub, log.p = TRUE)
	lda = stats::dnorm(ua, log = TRUE); ldb = stats::dnorm(ub, log = TRUE)
	ga = exp(la + lnb - lP); gb = exp(lb + lna - lP)
	hab = exp(la + ldb - log(kap) - lP)
	paa = -a * ga - (r / kap) * exp(la + ldb - lP)
	pbb = -b * gb - (r / kap) * exp(lb + lda - lP)
	nn = a^2 - 2 * r * a * b + b^2
	dl2 = (r * (1 - r^2) + a * b * (1 - r^2) - r * nn) / (1 - r^2)^2
	list(
		L = lP,
		La = ga,
		Lb = gb,
		Laa = paa - ga^2,
		Lbb = pbb - gb^2,
		Lab = hab - ga * gb,
		Lr = hab,
		Lrr = hab * dl2 - hab^2
	)
}

# expected log probability of a pair, its score psi and information W in the two predictors, and the score
# and curvature of rho, for (eta1, eta2) bivariate normal (m1, m2, v1, v2, cv), by a tensor Gauss-Hermite
# rule of K nodes a side; s1, s2 the signs 2 y - 1
.ae_pair_quad_r = function(m1, m2, v1, v2, cv, s1, s2, rho, K = 7L, chunk = 8192L) {
	n = length(m1)
	if (n <= chunk) return(.ae_pair_quad_block(m1, m2, v1, v2, cv, s1, s2, rho, K))
	# blocks of pairs keep the per-node temporaries in cache; the result is identical
	idx = split(seq_len(n), ceiling(seq_len(n) / chunk))
	res = lapply(idx, function(i) .ae_pair_quad_block(m1[i], m2[i], v1[i], v2[i], cv[i], s1[i], s2[i], rho, K))
	nm = names(res[[1L]])
	stats::setNames(lapply(nm, function(k) unlist(lapply(res, `[[`, k), use.names = FALSE)), nm)
}

.ae_pair_quad_block = function(m1, m2, v1, v2, cv, s1, s2, rho, K) {
	nodes = .ae_mf_gh(K)
	n = length(m1)
	sd1 = sqrt(pmax(v1, 0))
	l21 = ifelse(sd1 > 0, cv / pmax(sd1, 1e-300), 0)
	l22 = sqrt(pmax(v2 - l21^2, 0))
	r = s1 * s2 * rho
	ss = s1 * s2
	ell = psi1 = psi2 = w11 = w22 = w12 = rs = rc = numeric(n)
	for (k in seq_len(K)) {
		e1 = m1 + sd1 * nodes$z[k]
		for (l in seq_len(K)) {
			wt = nodes$w[k] * nodes$w[l]
			e2 = m2 + l21 * nodes$z[k] + l22 * nodes$z[l]
			d = .ae_pair_logp_r(s1 * e1, s2 * e2, r)
			ell = ell + wt * d$L
			psi1 = psi1 + wt * s1 * d$La
			psi2 = psi2 + wt * s2 * d$Lb
			w11 = w11 - wt * d$Laa
			w22 = w22 - wt * d$Lbb
			w12 = w12 - wt * ss * d$Lab
			rs = rs + wt * ss * d$Lr
			rc = rc + wt * d$Lrr
		}
	}
	list(ell = ell, psi1 = psi1, psi2 = psi2, w11 = w11, w22 = w22, w12 = w12, rs = rs, rc = rc)
}

# the four outcome probabilities P(s1, s2) = Phi2(s1 a, s2 b; s1 s2 r) as an n x 4 matrix, from one bivariate
# probability and the margins; any below 1e-5 is recomputed directly, since the subtraction loses it
.ae_pair_probs_r = function(a, b, r) {
	pa = stats::pnorm(a)
	pb = stats::pnorm(b)
	p11 = exp(.ae_pair_logphi2_r(a, b, r))
	P = cbind(p11, pa - p11, pb - p11, 1 - pa - pb + p11)
	sg1 = c(1, 1, -1, -1)
	sg2 = c(1, -1, 1, -1)
	for (j in 2:4) {
		bad = which(!(P[, j] >= 1e-5))
		if (length(bad)) {
			P[bad, j] = exp(.ae_pair_logphi2_r(sg1[j] * a[bad], sg2[j] * b[bad], sg1[j] * sg2[j] * r[bad]))
		}
	}
	P
}

# per-node terms of the pair problem at (a, b) with correlation r and observed signs (s1, s2): log probability,
# score, Fisher information (I11, I22, I12) over the four outcomes, with `deriv` its derivatives in a and b and
# with `rho_deriv` those in r. probabilities are floored at 1e-150 and non-finite derivatives set to zero
.ae_pair_node_fisher_r = function(a, b, r, s1, s2, deriv = TRUE, rho_deriv = FALSE) {
	n = length(a)
	kap = sqrt(1 - r^2)
	ta = (a - r * b) / kap
	tb = (b - r * a) / kap
	fa = stats::dnorm(a)
	fb = stats::dnorm(b)
	P = .ae_pair_probs_r(a, b, r)
	Pc = pmax(P, 1e-150)
	sg1 = c(1, 1, -1, -1)
	sg2 = c(1, -1, 1, -1)
	ss = sg1 * sg2
	# A_s = Phi(s2 tb), B_s = Phi(s1 ta)
	ptb = stats::pnorm(tb)
	ntb = stats::pnorm(-tb)
	pta = stats::pnorm(ta)
	nta = stats::pnorm(-ta)
	A = cbind(ptb, ntb, ptb, ntb)
	B = cbind(pta, pta, nta, nta)
	F11 = rowSums(A^2 / Pc)
	F22 = rowSums(B^2 / Pc)
	F12 = rowSums(rep(ss, each = n) * A * B / Pc)
	I11 = fa^2 * F11
	I22 = fb^2 * F22
	I12 = fa * fb * F12
	# the observed pair: its column of the matrix
	k = 1L + 2L * (s1 < 0) + (s2 < 0)
	ix = cbind(seq_len(n), k)
	Po = Pc[ix]
	psi1 = s1 * fa * A[ix] / Po
	psi2 = s2 * fb * B[ix] / Po
	out = list(lo = log(Po), psi1 = psi1, psi2 = psi2, I11 = I11, I22 = I22, I12 = I12)
	if (!deriv) return(out)
	ft = stats::dnorm(tb)
	fta = stats::dnorm(ta)
	# derivatives of A_s, B_s, P_s in a and b: dt_b/da = -r/kap, dt_b/db = 1/kap, dt_a/da = 1/kap, dt_a/db = -r/kap,
	# dP_s/da = s1 phi(a) A_s, dP_s/db = s2 phi(b) B_s
	sg1m = rep(sg1, each = n)
	sg2m = rep(sg2, each = n)
	Aa = (-r / kap * ft) * sg2m
	Ab = (ft / kap) * sg2m
	Ba = (fta / kap) * sg1m
	Bb = (-r / kap * fta) * sg1m
	Pa = (fa * sg1m) * A
	Pb = (fb * sg2m) * B
	P2 = Pc^2
	F11a = rowSums(2 * A * Aa / Pc - A^2 * Pa / P2)
	F11b = rowSums(2 * A * Ab / Pc - A^2 * Pb / P2)
	F22a = rowSums(2 * B * Ba / Pc - B^2 * Pa / P2)
	F22b = rowSums(2 * B * Bb / Pc - B^2 * Pb / P2)
	ssm = rep(ss, each = n)
	F12a = rowSums(ssm * ((Aa * B + A * Ba) / Pc - A * B * Pa / P2))
	F12b = rowSums(ssm * ((Ab * B + A * Bb) / Pc - A * B * Pb / P2))
	out$I11a = -2 * a * fa^2 * F11 + fa^2 * F11a
	out$I11b = fa^2 * F11b
	out$I22a = fb^2 * F22a
	out$I22b = -2 * b * fb^2 * F22 + fb^2 * F22b
	out$I12a = -a * fa * fb * F12 + fa * fb * F12a
	out$I12b = -b * fa * fb * F12 + fa * fb * F12b
	for (f in c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b")) out[[f]][!is.finite(out[[f]])] = 0
	if (!rho_deriv) return(out)
	# the probability of an outcome has P_r = s phi2(a, b; r) (Plackett, s = s1 s2 of the outcome, phi2 the bivariate
	# normal density, phi(a) phi(tb) / kap) and the mixed derivatives P_ar = s phi(a) phi(tb) (r b - a) / kap^3,
	# P_br = s phi(b) phi(ta) (r a - b) / kap^3 (dt_b/dr = (r b - a) / kap^3, dt_a/dr = (r a - b) / kap^3)
	phi2 = fa * ft / kap
	Pr = matrix(ssm * phi2, n, 4L)
	Par = matrix(ssm * (fa * ft * (r * b - a) / kap^3), n, 4L)
	Pbr = matrix(ssm * (fb * fta * (r * a - b) / kap^3), n, 4L)
	out$I11r = rowSums(2 * Pa * Par / Pc - Pa^2 * Pr / P2)
	out$I22r = rowSums(2 * Pb * Pbr / Pc - Pb^2 * Pr / P2)
	out$I12r = rowSums((Par * Pb + Pa * Pbr) / Pc - Pa * Pb * Pr / P2)
	Po_r = Pr[ix]
	out$s0 = Po_r / Po
	out$psi1r = Par[ix] / Po - Pa[ix] * Po_r / Po^2
	out$psi2r = Pbr[ix] / Po - Pb[ix] * Po_r / Po^2
	Lr = phi2 / Po
	nn = a^2 - 2 * r * a * b + b^2
	dl2 = s1 * s2 * (r * (1 - r^2) + a * b * (1 - r^2) - r * nn) / (1 - r^2)^2
	out$rc0 = Lr * dl2 - Lr^2
	for (f in c("I11r", "I22r", "I12r", "s0", "psi1r", "psi2r", "rc0")) out[[f]][!is.finite(out[[f]])] = 0
	out
}

# expected Fisher information problem of a pair under q, (eta_1, eta_2) bivariate normal, by a tensor
# Gauss-Hermite rule of K[i] nodes a side (K may differ by pair); returns ell, the expected score and information,
# with `deriv` their derivatives in the two means and with `rho_deriv` (needs `deriv`) those in the correlation
.ae_pair_quad_fisher_r = function(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv = FALSE, rho_deriv = FALSE) {
	n = length(m1)
	rho_deriv = deriv && rho_deriv
	v1 = rep_len(v1, n); v2 = rep_len(v2, n); cv = rep_len(cv, n); s1 = rep_len(s1, n); s2 = rep_len(s2, n)
	rho = rep_len(rho, n); K = rep_len(as.integer(K), n)
	nm = c(
		"ell",
		"psi1",
		"psi2",
		"I11",
		"I22",
		"I12",
		if (deriv) c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b"),
		if (rho_deriv) c("s0", "rc0", "psi1r", "psi2r", "I11r", "I22r", "I12r")
	)
	out = stats::setNames(replicate(length(nm), numeric(n), simplify = FALSE), nm)
	for (kk in sort(unique(K))) {
		ix = which(K == kk)
		nd = .ae_mf_gh(kk)
		sd1 = sqrt(pmax(v1[ix], 0))
		l21 = ifelse(sd1 > 0, cv[ix] / pmax(sd1, 1e-300), 0)
		l22 = sqrt(pmax(v2[ix] - l21^2, 0))
		acc = stats::setNames(replicate(length(nm), numeric(length(ix)), simplify = FALSE), nm)
		for (k in seq_len(kk)) {
			a = m1[ix] + sd1 * nd$z[k]
			for (l in seq_len(kk)) {
				wt = nd$w[k] * nd$w[l]
				b = m2[ix] + l21 * nd$z[k] + l22 * nd$z[l]
				nd1 = .ae_pair_node_fisher_r(a, b, rho[ix], s1[ix], s2[ix], deriv, rho_deriv)
				for (f in nm) acc[[f]] = acc[[f]] + wt * nd1[[switch(f, ell = "lo", f)]]
			}
		}
		for (f in nm) out[[f]][ix] = acc[[f]]
	}
	out
}

# the sums of the cross weights that the additive block needs: the dyad-level sums over time
# of cw, of cw [1, x] (gmc: the cross analogue of gm), the (theta, theta) block
# sum cw x_ij x_ji' and the (theta, actor) blocks. cw is symmetric in (i, j) within a period
.ae_pair_cross_r = function(cx, cw) {
	nr = cx$nr; Tt = cx$Tt; p = cx$p
	cw[!is.finite(cw)] = 0
	Csum = rowSums(cw, dims = 2L)
	gmc = .ae_eb_hprime_sums(cw, cx$X, nr, nr, Tt, p)
	pa = p + 1L
	x1 = function(t, k) {
		if (k == 1L) return(matrix(1, nr, nr))
		x = cx$X[, , k - 1L, t]
		x[!is.finite(x)] = 0
		x
	}
	xx = matrix(0, pa, pa)
	Mk = lapply(seq_len(pa), function(k) matrix(0, nr, nr))
	for (t in seq_len(Tt)) {
		ct = cw[, , t]
		for (k in seq_len(pa)) {
			xk = x1(t, k) * ct
			Mk[[k]] = Mk[[k]] + xk
			for (l in k:pa) {
				v = sum(xk * t(x1(t, l)))
				xx[k, l] = xx[k, l] + v
				if (l != k) xx[l, k] = xx[l, k] + v
			}
		}
	}
	list(
		Csum = Csum,
		gmc = gmc,
		xx = xx,
		Mk = Mk,
		d = rowSums(Csum),
		th_a = t(vapply(Mk, colSums, numeric(nr))),
		th_b = t(vapply(Mk, rowSums, numeric(nr)))
	)
}

# reference implementation of the cell (as pln_node() in src/ae_pln.h): log likelihood, score,
# information, its derivative, and the noise posterior's second moment and mean; processed in chunks
.ae_pln_node_r = function(y, eta, ve, K = .ae_pln_kin(max(ve)), chunk = 20000L) {
	n = length(y)
	ve = rep_len(ve, n)
	out = vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix = ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] = .ae_pln_node_block(y[ix], eta[ix], ve[ix], K)
	}
	if (!length(out)) return(list(
		l = numeric(),
		s = numeric(),
		w = numeric(),
		hp = numeric(),
		ee = numeric(),
		me = numeric()
	))
	if (length(out) == 1L) return(out[[1L]])
	nm = names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

.ae_pln_node_block = function(y, eta, ve, K) {
	n = length(y)
	eta = pmin(pmax(eta, -700), 700)
	ve = pmin(pmax(ve, 1e-8), 1e6)
	iv = 1 / ve
	ly = ifelse(y > 0, log(y), -Inf)
	# posterior mode of the noise by a monotone Newton iteration from the right of the root
	zeta = pmax(eta, ly)
	u = y + eta / ve
	act = seq_len(n)
	for (it in seq_len(100L)) {
		if (!length(act)) break
		zt = zeta[act]
		ez = exp(zt)
		step = (ez + zt / ve[act] - u[act]) / (ez + iv[act])
		zt = zt - step
		zeta[act] = zt
		act = act[!(abs(step) <= 1e-13 * pmax(1, abs(zt)))]
	}
	Xs = exp(zeta)
	s = 1 / sqrt(Xs + iv)
	estar = zeta - eta
	rule = .ae_pln_rule(K)
	Z = matrix(rule$z, n, K, byrow = TRUE)
	tt = s * Z
	# expm1(t) - t - t^2 / 2, by its series where the subtraction loses every digit
	r = expm1(tt) - tt - 0.5 * tt * tt
	sm = abs(tt) < 0.1
	if (any(sm)) {
		t1 = tt[sm]
		r[sm] = t1 * t1 * t1 / 6 * (1 + t1 / 4 * (1 + t1 / 5 * (1 + t1 / 6 * (1 + t1 / 7 * (1 + t1 / 8 *
				(1 + t1 / 9 * (1 + t1 / 10 * (1 + t1 / 11))))))))
	}
	lw = matrix(rule$lw, n, K, byrow = TRUE) - Xs * r
	ex = expm1(tt)
	mx = lw[cbind(seq_len(n), max.col(lw, ties.method = "first"))]
	p = exp(lw - mx)
	tot = rowSums(p)
	p = p / tot
	# a node whose weight underflowed carries nothing
	ex[p == 0] = 0
	ez1 = rowSums(p * Z)
	d = Z - ez1
	vz = rowSums(p * d * d)
	k3 = rowSums(p * d * d * d)
	eex = rowSums(p * ex)
	cze = rowSums(p * d * (ex - eex))
	me = estar + s * ez1
	gs = y * estar - Xs - estar * estar / (2 * ve)
	list(
		l = y * eta - lgamma(y + 1) - 0.5 * log(2 * pi * ve) + gs + log(s) + 0.918938533204672741780329736406 +
			mx + log(tot),
		s = me / ve,
		w = Xs * s * cze / ve,
		hp = -(s * s * s) * k3 / (ve * ve * ve),
		ee = s * s * vz + me * me,
		me = me
	)
}

# the same by the compiled kernel (ae_pln_node_cpp); zero off the observed cells and where
# the predictor is not finite
.ae_pln_terms_r = function(y, eta, ve, K = .ae_pln_kin(ve)) {
	ix = which(is.finite(y) & is.finite(eta))
	nd = .ae_pln_node_r(y[ix], eta[ix], ve, K)
	fill = function(v) {
		o = array(0, .ae_dims(y))
		o[ix] = v
		o
	}
	list(
		l = fill(nd$l),
		s = fill(nd$s),
		w = fill(nd$w),
		hp = fill(nd$hp),
		ee = fill(nd$ee),
		me = fill(nd$me),
		ve = ve
	)
}

# expected working problem of cells with predictor N(m, v), as .ae_mf_cells_r() for a binary cell,
# plus `ee`, the expected second moment of the noise; the response increment is capped at 30
.ae_pln_cells_r = function(y, m, v, ve, K = 9L, Kin = .ae_pln_kin(ve)) {
	obs = is.finite(y) & is.finite(m)
	ix = which(obs)
	mm = m[ix]
	vv = pmax(rep_len(v, length(y))[ix], 0)
	yy = y[ix]
	rules = .ae_pln_outer(K)
	rule_of = findInterval(vv, rules$vcut, left.open = TRUE) + 1L
	es = ew = el = ee = numeric(length(ix))
	for (r in seq_along(rules$z)) {
		sel = which(rule_of == r)
		if (!length(sel)) next
		sdv = sqrt(vv[sel])
		for (k in seq_along(rules$z[[r]])) {
			nd = .ae_pln_node_r(yy[sel], mm[sel] + sdv * rules$z[[r]][k], ve, Kin)
			es[sel] = es[sel] + rules$w[[r]][k] * nd$s
			ew[sel] = ew[sel] + rules$w[[r]][k] * nd$w
			el[sel] = el[sel] + rules$w[[r]][k] * nd$l
			ee[sel] = ee[sel] + rules$w[[r]][k] * nd$ee
		}
	}
	delta = es / pmax(ew, 1e-12)
	big = abs(delta) > 30
	if (any(big)) {
		delta[big] = sign(delta[big]) * 30
		ew[big] = abs(es[big]) / 30
	}
	z = w = ell = eev = array(0, .ae_dims(y))
	z[ix] = mm + delta
	w[ix] = ew
	ell[ix] = el
	eev[ix] = ee
	list(z = z, w = w, ell = ell, ee = eev)
}

# cavity weight of the cells for the covariance of the factor rows (see .ae_mf_info_w_r): the
# information of a Poisson-lognormal cell with variance v + ve at mean m
.ae_pln_info_w_r = function(y, m, v, ve, Kin = .ae_pln_kin(max(ve + pmax(v, 0), na.rm = TRUE))) {
	ix = which(is.finite(y) & is.finite(m))
	w = .ae_pln_node_r(y[ix], m[ix], ve + pmax(rep_len(v, length(y))[ix], 0), Kin)$w
	w[!is.finite(w) | w < 0] = 0
	out = array(0, .ae_dims(y))
	out[ix] = w
	out
}

# R reference for vectors of pairs: counts y1, y2, means m1, m2 and latent covariance (s11, s22, s12);
# returns the log probability L, the score (L1, L2) and the information (w11, w22, w12)
.ae_pln2_node_r = function(y1, y2, m1, m2, s11, s22, s12, K = .ae_pln2_kin(max(s11, s22)), chunk = 4000L) {
	n = length(y1)
	s11 = rep_len(s11, n); s22 = rep_len(s22, n); s12 = rep_len(s12, n)
	out = vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix = ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] = .ae_pln2_block(y1[ix], y2[ix], m1[ix], m2[ix], s11[ix], s22[ix], s12[ix], K)
	}
	if (!length(out)) return(
		list(L = numeric(), L1 = numeric(), L2 = numeric(), w11 = numeric(), w22 = numeric(), w12 = numeric())
	)
	if (length(out) == 1L) return(out[[1L]])
	nm = names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

.ae_pln2_block = function(y1, y2, m1, m2, s11, s22, s12, K) {
	n = length(y1)
	m1 = pmin(pmax(m1, -700), 700); m2 = pmin(pmax(m2, -700), 700)
	s11 = pmin(pmax(s11, 1e-8), 1e6); s22 = pmin(pmax(s22, 1e-8), 1e6)
	# hold the correlation away from one so the covariance stays invertible
	s12 = pmax(pmin(s12, 0.999999 * sqrt(s11 * s22)), -0.999999 * sqrt(s11 * s22))
	det = s11 * s22 - s12 * s12
	p11 = s22 / det; p22 = s11 / det; p12 = -s12 / det
	# start from the single-cell modes: e = zeta - m, exp(zeta) + zeta / s = y + m / s
	mode1 = function(y, m, s) {
		zeta = pmax(m, ifelse(y > 0, log(y), -Inf)); u = y + m / s
		act = seq_len(n)
		for (it in seq_len(100L)) {
			if (!length(act)) break
			zt = zeta[act]; ez = exp(zt)
			step = (ez + zt / s[act] - u[act]) / (ez + 1 / s[act])
			zt = zt - step; zeta[act] = zt
			act = act[!(abs(step) <= 1e-13 * pmax(1, abs(zt)))]
		}
		zeta - m
	}
	e1 = mode1(y1, m1, s11); e2 = mode1(y2, m2, s22)
	phi = function(e1, e2) y1 * e1 - exp(m1 + e1) + y2 * e2 - exp(m2 + e2) - 0.5 * (p11 * e1 * e1 + 2 * p12 * e1 * e2 + p22 * e2 * e2)
	act = seq_len(n)
	for (it in seq_len(100L)) {
		if (!length(act)) break
		a = act
		X1 = exp(m1[a] + e1[a]); X2 = exp(m2[a] + e2[a])
		g1 = y1[a] - X1 - (p11[a] * e1[a] + p12[a] * e2[a]); g2 = y2[a] - X2 - (p12[a] * e1[a] + p22[a] * e2[a])
		h11 = X1 + p11[a]; h22 = X2 + p22[a]; h12 = p12[a]
		dt = h11 * h22 - h12 * h12
		d1 = (h22 * g1 - h12 * g2) / dt; d2 = (h11 * g2 - h12 * g1) / dt
		# damp the step until the objective does not fall
		t = rep(1, length(a)); ok = rep(FALSE, length(a))
		n1 = e1[a]; n2 = e2[a]
		for (h in 0:30) {
			todo = which(!ok)
			if (!length(todo)) break
			c1 = e1[a][todo] + t[todo] * d1[todo]; c2 = e2[a][todo] + t[todo] * d2[todo]
			aa = a[todo]
			fc = y1[aa] * c1 - exp(m1[aa] + c1) + y2[aa] * c2 - exp(m2[aa] + c2) - 0.5 * (p11[aa] * c1 * c1 + 2 * p12[aa] * c1 * c2 + p22[aa] * c2 * c2)
			fo = y1[aa] * e1[aa] - exp(m1[aa] + e1[aa]) + y2[aa] * e2[aa] - exp(m2[aa] + e2[aa]) - 0.5 * (p11[aa] * e1[aa] * e1[aa] + 2 * p12[aa] * e1[aa] * e2[aa] + p22[aa] * e2[aa] * e2[aa])
			good = is.finite(fc) & fc >= fo - 1e-12 * (1 + abs(fo))
			n1[todo[good]] = c1[good]; n2[todo[good]] = c2[good]; ok[todo[good]] = TRUE
			t[todo[!good]] = t[todo[!good]] / 2
		}
		e1[a] = n1; e2[a] = n2
		act = a[!(pmax(abs(t * d1), abs(t * d2)) <= 1e-12 * pmax(1, abs(n1), abs(n2)))]
	}
	X1 = exp(m1 + e1); X2 = exp(m2 + e2)
	h11 = X1 + p11; h22 = X2 + p22; h12 = p12
	dh = h11 * h22 - h12 * h12
	a11 = h22 / dh; a22 = h11 / dh; a12 = -h12 / dh # Laplace covariance
	c11 = sqrt(a11); c21 = a12 / c11; c22 = sqrt(pmax(a22 - c21 * c21, 1e-300))
	rule = .ae_pln_rule(K)
	z = rule$z
	zi = rep(z, each = K); zj = rep(z, times = K) # product grid, i slowest
	lwg = rep(rule$lw, each = K) + rep(rule$lw, times = K)
	d1m = outer(c11, zi) # n x K^2
	d2m = outer(c21, zi) + outer(c22, zj)
	r = function(t) {
		out = expm1(t) - t - 0.5 * t * t
		sm = abs(t) < 0.1
		if (any(sm)) { t1 = t[sm]; out[sm] = t1 * t1 * t1 / 6 * (1 + t1 / 4 * (1 + t1 / 5 * (1 + t1 / 6 * (1 + t1 / 7 * (1 + t1 / 8 * (1 + t1 / 9 * (1 + t1 / 10 * (1 + t1 / 11)))))))) }
		out
	}
	lw = matrix(lwg, n, K * K, byrow = TRUE) - X1 * r(d1m) - X2 * r(d2m)
	mx = lw[cbind(seq_len(n), max.col(lw, ties.method = "first"))]
	pw = exp(lw - mx); tot = rowSums(pw); pw = pw / tot
	x1 = expm1(d1m); x2 = expm1(d2m)
	x1[pw == 0] = 0; x2[pw == 0] = 0
	E1 = rowSums(pw * d1m); E2 = rowSums(pw * d2m)
	q1 = d1m - E1; q2 = d2m - E2
	ex1 = rowSums(pw * x1); ex2 = rowSums(pw * x2)
	u1 = x1 - ex1; u2 = x2 - ex2
	# m_kj = Cov(e_k, X_j) = X_j Cov(d_k, expm1(d_j))
	M11 = X1 * rowSums(pw * q1 * u1); M12 = X2 * rowSums(pw * q1 * u2)
	M21 = X1 * rowSums(pw * q2 * u1); M22 = X2 * rowSums(pw * q2 * u2)
	Ee1 = e1 + E1; Ee2 = e2 + E2
	w11 = p11 * M11 + p12 * M21; w12a = p11 * M12 + p12 * M22
	w12b = p12 * M11 + p22 * M21; w22 = p12 * M12 + p22 * M22
	gs = y1 * e1 - X1 + y2 * e2 - X2 - 0.5 * (p11 * e1 * e1 + 2 * p12 * e1 * e2 + p22 * e2 * e2)
	L = y1 * m1 + y2 * m2 - lgamma(y1 + 1) - lgamma(y2 + 1) - 0.5 * log(det / (a11 * a22 - a12 * a12)) + gs + mx + log(tot)
	list(
		L = L,
		L1 = p11 * Ee1 + p12 * Ee2,
		L2 = p12 * Ee1 + p22 * Ee2,
		w11 = w11,
		w22 = w22,
		w12 = 0.5 * (w12a + w12b)
	)
}

# expected working problem of the observed cells of `y` with predictor N(m, v): response z = m + E psi / w,
# weight w = -E psi' and expected log likelihood ell, zero elsewhere; the increment is capped at 30 as in IRLS
.ae_mf_cells_r = function(y, m, v, link, K = 9L) {
	obs = is.finite(y)
	ix = which(obs)
	sgn = 2 * y[ix] - 1
	mm = m[ix]
	nodes = .ae_mf_gh(K)
	sdv = sqrt(pmax(v[ix], 0))
	probit = identical(link, "probit")
	es = ew = el = 0
	for (k in seq_along(nodes$z)) {
		e = mm + sdv * nodes$z[k]
		se = sgn * e
		if (probit) {
			lps = stats::pnorm(se, log.p = TRUE)
			lam = exp(-0.5 * e * e - 0.91893853320467267 - lps)
			psi = sgn * lam
			w = lam * (lam + se)
		} else {
			lps = stats::plogis(se, log.p = TRUE)
			p1 = exp(lps)
			psi = sgn * (1 - p1)
			w = p1 * (1 - p1)
		}
		es = es + nodes$w[k] * psi
		ew = ew + nodes$w[k] * w
		el = el + nodes$w[k] * lps
	}
	delta = es / pmax(ew, 1e-12)
	big = abs(delta) > 30
	if (any(big)) {
		delta[big] = sign(delta[big]) * 30
		ew[big] = abs(es[big]) / 30
	}
	z = w = ell = array(0, dim(y))
	z[ix] = mm + delta
	w[ix] = ew
	ell[ix] = el
	list(z = z, w = w, ell = ell)
}

# cavity information weight of a cell for the factor-row covariances, nu = -d^2 log Z / dm^2 with
# Z(m, v) = E p(y | eta), eta ~ N(m, v); the expected weight E w(eta) would count the other side's variance twice
.ae_mf_info_w_r = function(y, m, v, link, K = 9L) {
	obs = is.finite(y)
	ix = which(obs)
	sgn = 2 * y[ix] - 1
	mm = m[ix]
	vv = pmax(v[ix], 0)
	w = if (identical(link, "probit")) {
		sq = sqrt(1 + vv)
		t = sgn * mm / sq
		lam = exp(stats::dnorm(t, log = TRUE) - stats::pnorm(t, log.p = TRUE))
		out = lam * (lam + t)
		far = t < -30
		if (any(far)) out[far] = 1 - 1 / t[far]^2 + 6 / t[far]^4
		out / (1 + vv)
	} else {
		nodes = .ae_mf_gh(K)
		sdv = sqrt(vv)
		z0 = z1 = z2 = 0
		for (k in seq_along(nodes$z)) {
			p1 = stats::plogis(sgn * (mm + sdv * nodes$z[k]))
			z0 = z0 + nodes$w[k] * p1
			z1 = z1 + nodes$w[k] * sgn * p1 * (1 - p1)
			z2 = z2 + nodes$w[k] * p1 * (1 - p1) * (1 - 2 * p1)
		}
		(z1 / z0)^2 - z2 / z0
	}
	w[!is.finite(w) | w < 0] = 0
	out = array(0, dim(y))
	out[ix] = w
	out
}

# R reference for the general rank, one system at a time
.ae_mf_row_solve_r = function(Gflat, L, rhs) {
	n = nrow(rhs); R = ncol(rhs)
	x = matrix(0, n, R); Sig = array(0, c(R, R, n))
	for (i in seq_len(n)) {
		g = matrix(Gflat[i, ], R, R) + L
		g = (g + t(g)) / 2
		Si = .ae_safe_inverse(g)
		Sig[, , i] = Si
		x[i, ] = Si %*% rhs[i, ]
	}
	list(x = x, Sig = Sig)
}

# normal log density and both log tail probabilities at a probit predictor,
# computed once per iterate; `outcome` is the log probability of each cell's
# observed outcome and `one` the cells with y = 1 (NULL: not needed)
.ae_probit_terms_r = function(eta, one) {
	lp1 = stats::pnorm(eta, log.p = TRUE)
	lp0 = stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
	outcome = NULL
	if (!is.null(one)) {
		outcome = lp0
		outcome[one] = lp1[one]
	}
	list(ld = stats::dnorm(eta, log = TRUE), lp1 = lp1, lp0 = lp0, outcome = outcome)
}

# sum_t h' [1, x] for every dyad (n_row * n_col x (p + 1)); h' can be negative, so
# the R reference sums its positive and negative parts apart
.ae_eb_hprime_sums_r = function(hprime, X, nr, nc, Tt, p) {
	hp = hprime
	hp[!is.finite(hp)] = 0
	pos = ae_bcd_stats_cpp(pmax(hp, 0), hp, X, nr, nc, Tt, p)$gm
	neg = ae_bcd_stats_cpp(pmax(-hp, 0), hp, X, nr, nc, Tt, p)$gm
	pos - neg
}

# d h / d eta of the probit Fisher weight from the probit terms (the reference of ae_probit_hprime_cpp)
ae_probit_hprime_r = function(eta, ld, lp1, lp0) {
	lam1 = exp(ld - lp1)
	lam0 = exp(ld - lp0)
	lam1 * lam0 * (lam0 - lam1 - 2 * eta)
}

# runs the ALS estimators of the calling frame on the R references: every compiled quadrature kernel of the
# static and dynamic empirical-Bayes fits is replaced by its reference, and the table family of the pair grid
# is off (a reference table is one exact table)
local_als_r_engine = function(env = parent.frame()) {
	row_solve = get(".ae_mf_row_solve", asNamespace("lame"))
	testthat::local_mocked_bindings(
		.ae_mf_cells = function(y, m, v, link, K = 9L, ve = NULL, kin = NULL) {
			if (identical(link, "log")) return(.ae_pln_cells(y, m, v, ve, K, if (is.null(kin)) .ae_pln_kin(ve) else kin))
			if (identical(link, "identity")) return(.ae_gauss_cells(y, m, v, ve))
			.ae_mf_cells_r(y, m, v, link, K)
		},
		.ae_mf_info_w = function(y, m, v, link, K = 9L, ve = NULL) {
			if (identical(link, "log")) return(.ae_pln_info_w(y, m, v, ve))
			if (identical(link, "identity")) return(.ae_gauss_info_w(y, m, v, ve))
			.ae_mf_info_w_r(y, m, v, link, K)
		},
		.ae_mf_row_solve = function(Gflat, L, rhs) {
			if (ncol(rhs) <= 2L) row_solve(Gflat, L, rhs) else .ae_mf_row_solve_r(Gflat, L, rhs)
		},
		.ae_pair_logphi2 = function(a, b, r) .ae_pair_logphi2_r(a, b, rep_len(r, length(a))),
		.ae_pair_logp = function(a, b, r) .ae_pair_logp_r(a, b, rep_len(r, length(a))),
		.ae_pair_quad = function(m1, m2, v1, v2, cv, s1, s2, rho, K = 7L, chunk = 8192L) {
			.ae_pair_quad_r(m1, m2, v1, v2, cv, s1, s2, rho, K, chunk)
		},
		.ae_pair_quad_fisher = function(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv = FALSE, rho_deriv = FALSE) {
			.ae_pair_quad_fisher_r(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv, rho_deriv)
		},
		.ae_pair_cross = function(cx, cw) .ae_pair_cross_r(cx, cw),
		.ae_pln2_node = function(y1, y2, m1, m2, s11, s22, s12, K = .ae_pln2_kin(max(s11, s22))) {
			.ae_pln2_node_r(y1, y2, m1, m2, s11, s22, s12, K)
		},
		.ae_pln_terms = function(y, eta, ve, K = .ae_pln_kin(ve)) .ae_pln_terms_r(y, eta, ve, K),
		.ae_pln_cells = function(y, m, v, ve, K = 9L, Kin = .ae_pln_kin(ve)) .ae_pln_cells_r(y, m, v, ve, K, Kin),
		.ae_pln_info_w = function(y, m, v, ve, Kin = .ae_pln_kin(max(ve + pmax(v, 0), na.rm = TRUE))) {
			.ae_pln_info_w_r(y, m, v, ve, Kin)
		},
		.ae_probit_terms = function(eta, one) .ae_probit_terms_r(eta, one),
		ae_probit_hprime_cpp = function(eta, ld, lp1, lp0, threads = 1L) ae_probit_hprime_r(eta, ld, lp1, lp0),
		.ae_eb_hprime_sums = function(hprime, X, nr, nc, Tt, p) .ae_eb_hprime_sums_r(hprime, X, nr, nc, Tt, p),
		.ae_bvn_log = function(a, b, rho) .lame_binary_dynamic_bvn_log(a, b, rho),
		.lde_pln_terms = function(y, eta, ve, lgy = NULL, K = .lde_pln_kin(ve)) .lde_pln_terms_r(y, eta, ve, lgy, K),
		.lde_pln_expect = function(y, m, v, ve, lgy = NULL) .lde_pln_expect_ref(y, m, v, ve, lgy),
		.lde_pair_family_on = function() FALSE,
		lde_pair_table_cpp = function(r, cc, du, klo, khi, vmask, threads = 1L) .lde_pair_table_r(r, cc, du, klo, khi, vmask),
		lde_pair_apply_cpp = function(tab, m1, m2, sd1, sd2, variant, cutoff, deriv, rho_deriv, threads = 1L) {
			.lde_pair_apply_r(tab, m1, m2, sd1, sd2, variant, cutoff, deriv, rho_deriv)
		},
		lde_pair_direct_cpp = function(m1, m2, sd1, sd2, variant, type1, type2, K1, K2, rho, cc, du, cutoff, deriv, rho_deriv, threads = 1L) {
			.lde_pair_direct_r(m1, m2, sd1, sd2, variant, type1, type2, K1, K2, rho, cc, du, cutoff, deriv, rho_deriv)
		},
		.package = "lame",
		.env = env
	)
}
