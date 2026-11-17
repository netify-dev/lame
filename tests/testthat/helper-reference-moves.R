# R references only the tests use for sampler moves with compiled kernels: the coefficient translation,
# the ordinal scale move, the replicated Gaussian block and the hidden cells of censored nominations

# Poisson coefficient translations: coefficient k and the latent array move by
# delta x_k together, drawn by an elliptical slice step
.lame_ltm_coef = function(cx, Z, EZ, beta, Xk, iV0) {
	stats = c(coef_calls = 0, coef_evals = 0, coef_stuck = 0)
	for (k in sample.int(length(beta))) {
		x = Xk[[k]]
		xo = x[cx$ok_idx]
		if (!all(is.finite(xo)) || !any(xo != 0)) next
		ez = exp(Z[cx$ok_idx] + cx$off_ok)
		ysx = sum(cx$y_ok * xo)
		ll = function(dl) dl * ysx - sum(ez * expm1(dl * xo))
		a1 = drop(iV0 %*% beta)[k]
		a2 = iV0[k, k]
		if (!is.finite(a1) || !is.finite(a2) || a2 <= 0) next
		ans = .lame_ess_step(0, -a1 / a2, matrix(1 / sqrt(a2)), ll, 0)
		stats["coef_calls"] = stats["coef_calls"] + 1
		stats["coef_evals"] = stats["coef_evals"] + ans$evals
		stats["coef_stuck"] = stats["coef_stuck"] + ans$stuck
		d = ans$x
		if (d == 0) next
		Z = Z + d * x
		EZ = EZ + d * x
		beta[k] = beta[k] + d
	}
	list(Z = Z, EZ = EZ, beta = beta, stats = stats)
}

# inner product of two matrices in the metric of the covariance of the latent array;
# `up` the upper-triangle mask (square networks)
.lame_osm_dot = function(x, y, rho, s2, mode, up = NULL) {
	if (identical(mode, "bip")) return(sum(x * y) / s2)
	if (is.null(up)) up = upper.tri(x)
	if (identical(mode, "sym")) {
		xs = (x + t(x)) / 2
		ys = (y + t(y)) / 2
		return(sum(xs[up] * ys[up]) / s2)
	}
	xt = t(x)
	yt = t(y)
	sum(((x + xt) * (y + yt))[up]) / (2 * (1 + rho) * s2) +
		sum(((x - xt) * (y - yt))[up]) / (2 * (1 - rho) * s2) +
		sum(diag(x) * diag(y)) / ((1 + rho) * s2)
}

# the shift and the scale of the latent array of one period given its predictor
.lame_osm_period = function(Zt, Et, rho, s2, mode) {
	nr = nrow(Zt)
	nc = ncol(Zt)
	up = if (identical(mode, "bip")) NULL else upper.tri(Zt)
	ncell = .lame_osm_ncell(nr, nc, mode)
	if (identical(mode, "sym")) {
		Zt = (Zt + t(Zt)) / 2
		Et = (Et + t(Et)) / 2
	}
	# shift: the density of Z + d around E is Gaussian in d
	res = Zt - Et
	if (identical(mode, "sym")) {
		prec = ncell / s2
		mu = -sum(res[up]) / ncell
	} else if (identical(mode, "bip")) {
		prec = ncell / s2
		mu = -sum(res) / ncell
	} else {
		prec = ncell / ((1 + rho) * s2)
		mu = -sum(res) / ncell
	}
	Zt = Zt + stats::rnorm(1L, mu, 1 / sqrt(prec))
	# scale
	A = .lame_osm_dot(Zt, Zt, rho, s2, mode, up)
	B = .lame_osm_dot(Zt, Et, rho, s2, mode, up)
	cc = .lame_osm_draw_c(ncell, A, B)
	if (cc != 1) Zt = cc * Zt
	Zt
}

# reference visit written directly on the latent array (the definition that
# .lame_osm_visit computes from one-pass sums); arguments as for .lame_osm_visit
.lame_osm_visit_ref = function(Z, E, m, mode, rho, s2, A_pri, B_pri, d_th, affine = TRUE, orbit = TRUE) {
	if (!all(is.finite(Z)) || !all(is.finite(E))) return(list(Z = Z, c = 1, E = E))
	nr = dim(Z)[1L]
	nc = dim(Z)[2L]
	nt = dim(Z)[3L]
	rho = rep_len(rho, nt)
	up = if (identical(mode, "bip")) NULL else upper.tri(matrix(0, nr, nc))
	ncell = .lame_osm_ncell(nr, nc, mode)
	if (identical(mode, "sym")) {
		sym_t = function(a) {
			for (tt in seq_len(dim(a)[3L])) a[, , tt] = (a[, , tt] + t(a[, , tt])) / 2
			a
		}
		E = sym_t(E)
		if (!is.null(m)) m = sym_t(m)
	}
	# the orbit of the whole state: u = Z - E_scaled, the unscaled part m
	cc = 1
	if (orbit) {
		Es = if (is.null(m)) E else E - m
		A = A_pri
		B = B_pri
		for (tt in seq_len(nt)) {
			u = Z[, , tt] - Es[, , tt]
			A = A + .lame_osm_dot(u, u, rho[tt], s2, mode, up)
			if (!is.null(m)) B = B + .lame_osm_dot(u, m[, , tt], rho[tt], s2, mode, up)
		}
		cc = .lame_osm_draw_c(nt * ncell + d_th, A, B)
		if (cc != 1) {
			Z = cc * Z
			E = if (is.null(m)) cc * E else cc * Es + m
		}
	}
	if (affine) {
		for (tt in seq_len(nt)) Z[, , tt] = .lame_osm_period(Z[, , tt], E[, , tt], rho[tt], s2, mode)
	}
	list(Z = Z, c = cc, E = E)
}

# R reference of cbin_implied_labels_cpp, used in the tests
.lame_cbin_implied_labels_r = function(Z, Y, odmax, hid) {
	Ya = Y
	cap = rep_len(odmax, nrow(Y))
	n1 = rowSums(!is.na(Y) & Y > 0)
	for (i in which(rowSums(hid) > 0L)) {
		h = which(hid[i, ])
		zh = Z[i, h]
		lab = numeric(length(h))
		r = max(cap[i] - n1[i], 0)
		pos = which(zh > 0)
		if (r > 0 && length(pos)) {
			top = pos[order(zh[pos], decreasing = TRUE)][seq_len(min(r, length(pos)))]
			lab[top] = 1
		}
		Ya[i, h] = lab
	}
	Ya
}

# log of the standard normal mass on (lo, hi), stable in both tails
.lame_log_norm_mass = function(lo, hi) {
	out = rep(-Inf, length(lo))
	ok = hi > lo
	up = ok & lo >= 0
	dn = ok & !up
	if (any(up)) {
		a = stats::pnorm(lo[up], lower.tail = FALSE, log.p = TRUE)
		b = stats::pnorm(hi[up], lower.tail = FALSE, log.p = TRUE)
		out[up] = a + log1p(-exp(b - a))
	}
	if (any(dn)) {
		a = stats::pnorm(hi[dn], log.p = TRUE)
		b = stats::pnorm(lo[dn], log.p = TRUE)
		out[dn] = a + log1p(-exp(b - a))
	}
	out
}

.lame_cbin_hidden_moves_r = function(Z, EZ, rho, Y, odmax, hid, bip) {
	n_r = nrow(Z); n_c = ncol(Z)
	sq = !bip
	cond_sd = if (sq) sqrt(1 - rho^2) else 1
	cap = rep_len(odmax, n_r)
	o1 = !is.na(Y) & Y > 0
	o0 = !is.na(Y) & Y == 0
	if (sq) { diag(o1) = FALSE; diag(o0) = FALSE }
	n1 = rowSums(o1)
	r = cap - n1
	z1 = Z; z1[!o1] = Inf
	m1 = apply(z1, 1L, min)
	z0 = Z; z0[!o0] = -Inf
	M0 = apply(z0, 1L, max)
	P0 = rowSums(o0 & Z > 0)
	# the hidden cells of each row, in a random order
	hl = lapply(seq_len(n_r), function(i) {
		h = which(hid[i, ])
		if (length(h) > 1L) h[sample.int(length(h))] else h
	})
	n_round = max(lengths(hl))
	if (!n_round) return(Z)
	for (k in seq_len(n_round)) {
		rows = which(lengths(hl) >= k)
		cols = vapply(rows, function(i) hl[[i]][k], integer(1))
		blocks = if (sq) list(rows < cols, rows > cols) else list(rep(TRUE, length(rows)))
		for (bl in blocks) {
			if (!any(bl)) next
			ii = rows[bl]; jj = cols[bl]
			pos = ii + (jj - 1L) * n_r
			x_old = Z[pos]
			# the other hidden cells of each row above 0, M0 and m1
			cnt = function(t) rowSums(hid[ii, , drop = FALSE] & (Z[ii, , drop = FALSE] > t)) - (x_old > t)
			a0 = cnt(0); aM = cnt(M0[ii]); am = cnt(m1[ii])
			M = M0[ii]; m = m1[ii]; rr = r[ii]; cp = cap[ii]; nn = n1[ii]; p0 = P0[ii]
			Mp = pmax(M, 0)
			lo = cbind(-Inf, 0, Mp, m)
			hi = cbind(0, Mp, m, Inf)
			feas = matrix(FALSE, length(ii), 4L)
			for (g in 1:4) {
				gt0 = g >= 2L
				gtM = g >= 3L & M > 0
				gtm = g == 4L
				K = nn + a0 + p0 + gt0
				ok_le = K <= cp & p0 == 0
				ok_gt = K > cp & (am + gtm) <= rr & (M <= 0 | (aM + gtM) >= rr)
				feas[, g] = (ok_le | ok_gt) & hi[, g] > lo[, g]
			}
			mu = EZ[pos]
			if (sq && rho != 0) {
				mp = jj + (ii - 1L) * n_r
				mu = mu + rho * (Z[mp] - EZ[mp])
			}
			lm = .lame_log_norm_mass((lo - mu) / cond_sd, (hi - mu) / cond_sd)
			lm = matrix(lm, length(ii), 4L)
			lm[!feas] = -Inf
			mx = apply(lm, 1L, max)
			good = is.finite(mx)
			if (!any(good)) next
			w = exp(lm[good, , drop = FALSE] - mx[good])
			w = w / rowSums(w)
			u = stats::runif(sum(good))
			g_pick = 1L + rowSums(u > t(apply(w, 1L, cumsum)))
			g_pick = pmin(g_pick, 4L)
			sel = cbind(seq_len(sum(good)), g_pick)
			lb = lo[good, , drop = FALSE][sel]; ub = hi[good, , drop = FALSE][sel]
			val = mu[good] + cond_sd * rtnorm_interval_logp(
				rep(0, sum(good)),
				(lb - mu[good]) / cond_sd,
				(ub - mu[good]) / cond_sd
			)
			bad = !is.finite(val)
			if (any(bad)) val[bad] = x_old[good][bad]
			Z[pos[good]] = val
		}
	}
	Z
}
