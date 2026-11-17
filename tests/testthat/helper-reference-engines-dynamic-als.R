# R references of the compiled pair-grid and Poisson-lognormal kernels of the dynamic ALS
# empirical-Bayes estimator (src/dynamic_eb_pair.cpp, src/dynamic_eb_pln.cpp)

# the 12 + 4 x 7 functions at the nodes (a, b), n x 40: the information functions of .ae_pair_node_fisher_r() and, for each
# outcome (+,+) (+,-) (-,+) (-,-), the log probability of the observed pair and its derivatives from the log-domain terms
# of the bivariate probit (no floor of the probabilities: they grow in the tails)
.lde_pair_node_r = function(a, b, r, vmask = 15L) {
	n = length(a)
	nd = .ae_pair_node_fisher_r(a, b, rep_len(r, n), rep(1, n), rep(1, n), TRUE, TRUE)
	out = matrix(0, n, 40L)
	info = c("I11", "I22", "I12", "I11a", "I11b", "I22a", "I22b", "I12a", "I12b", "I11r", "I22r", "I12r")
	for (k in seq_along(info)) out[, k] = nd[[info[k]]]
	sg1 = c(1, 1, -1, -1)
	sg2 = c(1, -1, 1, -1)
	for (j in 1:4) {
		if (!bitwAnd(vmask, bitwShiftL(1L, j - 1L))) next
		s1 = sg1[j]; s2 = sg2[j]
		A = s1 * a; B = s2 * b; R = rep_len(s1 * s2 * r, n)
		d = .ae_pair_logp_r(A, B, R)
		om = 1 - R^2
		LAR = d$Lr * (R * B - A) / om - d$La * d$Lr
		LBR = d$Lr * (R * A - B) / om - d$Lb * d$Lr
		cols = cbind(d$L, s1 * d$La, s2 * d$Lb, s1 * s2 * d$Lr, d$Lrr, s2 * LAR, s1 * LBR)
		cols[!is.finite(cols)] = 0
		out[, 12L + 7L * (j - 1L) + 1:7] = cols
	}
	out
}

.lde_pair_table_r = function(rho, cc, du, klo, khi, vmask = 15L) {
	k = klo:khi
	x = cc * sinh(k * du)
	N = length(k)
	a = rep(x, each = N) # index i of the first cell varies slowest: [i * N + j]
	b = rep(x, times = N)
	nodes = .lde_pair_node_r(a, b, rho, vmask)
	list(
		table = as.vector(nodes),
		N = N,
		cc = cc,
		du = du,
		klo = klo,
		khi = khi,
		rho = rho,
		vmask = vmask,
		cols = nodes
	)
}

# the contraction of the columns (the same output as lde_pair_apply_cpp())
.lde_pair_apply_r = function(tab, m1, m2, sd1, sd2, variant, cutoff = 7.5, deriv = TRUE, rho_deriv = FALSE) {
	N = tab$N
	x = tab$cc * sinh((tab$klo:tab$khi) * tab$du)
	jac = tab$cc * cosh((tab$klo:tab$khi) * tab$du) * tab$du
	n = length(m1)
	rho_deriv = deriv && rho_deriv
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
	# output column -> table column
	info = c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7, 8, 9, NA, NA, NA, NA, 10, 11, 12)
	obs = c(1, 2, 3, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4, 5, 6, 7, NA, NA, NA)
	out = stats::setNames(replicate(length(nm), numeric(n), simplify = FALSE), nm)
	ok = integer(n)
	K = integer(n)
	cols = matrix(tab$table, N * N, 40L)
	for (p in seq_len(n)) {
		w1 = jac * stats::dnorm((x - m1[p]) / sd1[p]) / sd1[p]
		w2 = jac * stats::dnorm((x - m2[p]) / sd2[p]) / sd2[p]
		i1 = which(abs(x - m1[p]) <= cutoff * sd1[p])
		i2 = which(abs(x - m2[p]) <= cutoff * sd2[p])
		if (!length(i1) || !length(i2) || m1[p] - cutoff * sd1[p] < x[1L] || m1[p] + cutoff * sd1[p] > x[N] ||
				m2[p] - cutoff * sd2[p] < x[1L] || m2[p] + cutoff * sd2[p] > x[N]) next
		wt = as.vector(outer(w1[i1], w2[i2]))
		ix = as.vector(outer((i1 - 1L) * N, i2, "+"))
		for (c in seq_along(nm)) {
			tc = if (!is.na(info[c])) info[c] else 12L + 7L * variant[p] + obs[c]
			out[[c]][p] = sum(wt * cols[ix, tc])
		}
		ok[p] = 1L
		K[p] = max(length(i1), length(i2))
	}
	c(out, list(ok = ok, K = K))
}

# the direct kernel's R reference: the product of the rules a predictor, the functions of .lde_pair_node_r() at the nodes
.lde_pair_direct_r = function(
	m1,
	m2,
	sd1,
	sd2,
	variant,
	type1,
	type2,
	K1,
	K2,
	rho,
	cc,
	du,
	cutoff = 7.5,
	deriv = TRUE,
	rho_deriv = FALSE
) {
	n = length(m1)
	rho_deriv = deriv && rho_deriv
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
	info = c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7, 8, 9, NA, NA, NA, NA, 10, 11, 12)
	obs = c(1, 2, 3, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4, 5, 6, 7, NA, NA, NA)
	out = stats::setNames(replicate(length(nm), numeric(n), simplify = FALSE), nm)
	ok = integer(n)
	K = integer(n)
	rule = function(m, sd, type, k) {
		if (!(sd > 0) || (type == 0L && k <= 1L)) return(list(x = m, w = 1))
		if (type == 0L) {
			g = .ae_mf_gh(k)
			return(list(x = m + sd * g$z, w = g$w))
		}
		i0 = ceiling(asinh((m - cutoff * sd) / cc) / du)
		i1 = floor(asinh((m + cutoff * sd) / cc) / du)
		u = (i0:i1) * du
		x = cc * sinh(u)
		list(x = x, w = cc * cosh(u) * du * stats::dnorm((x - m) / sd) / sd)
	}
	for (p in seq_len(n)) {
		r1 = rule(m1[p], sd1[p], type1[p], K1[p])
		r2 = rule(m2[p], sd2[p], type2[p], K2[p])
		a = rep(r1$x, times = length(r2$x))
		b = rep(r2$x, each = length(r1$x))
		w = as.vector(outer(r1$w, r2$w))
		cols = .lde_pair_node_r(a, b, rho, bitwShiftL(1L, variant[p]))
		for (c in seq_along(nm)) {
			tc = if (!is.na(info[c])) info[c] else 12L + 7L * variant[p] + obs[c]
			out[[c]][p] = sum(w * cols[, tc])
		}
		ok[p] = 1L
		K[p] = max(length(r1$x), length(r2$x))
	}
	c(out, list(ok = ok, K = K))
}

# R reference for pln_node() of src/ae_pln.h with the fourth cumulant, in chunks to bound memory
.lde_pln_node_block = function(y, eta, ve, K) {
	n = length(y)
	eta = pmin(pmax(eta, -700), 700)
	ve = pmin(pmax(ve, 1e-8), 1e6)
	iv = 1 / ve
	ly = ifelse(y > 0, log(y), -Inf)
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
	ex[p == 0] = 0
	ez1 = rowSums(p * Z)
	d = Z - ez1
	vz = rowSums(p * d * d)
	k3 = rowSums(p * d * d * d)
	k4 = rowSums(p * d * d * d * d) - 3 * vz * vz
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
		hpp = -(s * s * s * s) * k4 / (ve * ve * ve * ve),
		ee = s * s * vz + me * me,
		me = me
	)
}

.lde_pln_node_r = function(y, eta, ve, K = .lde_pln_kin(max(ve)), chunk = 20000L) {
	n = length(y)
	ve = rep_len(ve, n)
	out = vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix = ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] = .lde_pln_node_block(y[ix], eta[ix], ve[ix], K)
	}
	if (!length(out)) return(list(
		l = numeric(),
		s = numeric(),
		w = numeric(),
		hp = numeric(),
		hpp = numeric(),
		ee = numeric(),
		me = numeric()
	))
	if (length(out) == 1L) return(out[[1L]])
	nm = names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

# R reference of lde_pln_expect_cpp(): the cell terms averaged over the predictor
# m + sqrt(v) z by an outer Gauss-Hermite rule; `lm` is the log expected probability of the count
.lde_pln_expect_r = function(y, m, v, ve, Kin, K2) {
	d = .ae_dims(y)
	ix = which(is.finite(y) & is.finite(m))
	mm = m[ix]
	vv = pmax(rep_len(v, length(y))[ix], 0)
	yy = y[ix]
	rules = .lde_pln_outer_rule(suppressWarnings(max(vv, 0)))
	rule_of = findInterval(vv, rules$vcut, left.open = TRUE) + 1L
	n = length(ix)
	acc = list(es = numeric(n), ew = numeric(n), eh = numeric(n), ehh = numeric(n), el = numeric(n), ee = numeric(n))
	for (r in seq_along(rules$z)) {
		sel = which(rule_of == r)
		if (!length(sel)) next
		sdv = sqrt(vv[sel])
		for (k in seq_along(rules$z[[r]])) {
			nd = .lde_pln_node_r(yy[sel], mm[sel] + sdv * rules$z[[r]][k], ve, Kin)
			wk = rules$w[[r]][k]
			acc$es[sel] = acc$es[sel] + wk * nd$s
			acc$ew[sel] = acc$ew[sel] + wk * nd$w
			acc$eh[sel] = acc$eh[sel] + wk * nd$hp
			acc$ehh[sel] = acc$ehh[sel] + wk * nd$hpp
			acc$el[sel] = acc$el[sel] + wk * nd$l
			acc$ee[sel] = acc$ee[sel] + wk * nd$ee
		}
	}
	acc$lm = .lde_pln_node_r(yy, mm, ve + vv, K2)$l
	lapply(acc, function(a) { o = array(0, d); o[ix] = a; o })
}

# cell terms of .lde_pln_terms() on the observed entries by the R reference, zero elsewhere
.lde_pln_terms_r = function(y, eta, ve, lgy = NULL, K = .lde_pln_kin(ve)) {
	d = .ae_dims(y)
	ix = which(is.finite(y) & is.finite(eta))
	nd = .lde_pln_node_r(y[ix], eta[ix], ve, K)
	lapply(nd, function(v) { o = array(0, d); o[ix] = v; o })
}

# expected cells of .lde_pln_expect() by the R reference, with the rules the compiled path uses
.lde_pln_expect_ref = function(y, m, v, ve, lgy = NULL) {
	d = .ae_dims(y)
	v = array(pmax(rep_len(as.numeric(v), length(y)), 0), d)
	Kin = .lde_pln_kin(ve)
	vmax = suppressWarnings(max(v[is.finite(y) & is.finite(m)], 0))
	K2 = .ae_pln_kin(1.5 * (ve + vmax))
	.lde_pln_expect_r(y, m, v, ve, Kin, K2)
}
