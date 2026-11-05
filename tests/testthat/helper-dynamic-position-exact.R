# exact-posterior check of the position-path move on a tiny Poisson panel: prior draws weighted by the Poisson
# mass are posterior draws, so one visit must keep every weighted mean; returns z scores and the effective size
# mode: 'dir', 'sym' or 'bip'; kind: 'ar1' or 't'; rk: rank of the positions; rail: coordinate box that truncates the prior
dpmt_exact = function(mode, kind = 'ar1', rk = 1L, ndraw = 2e5, seed = 1L, rail = Inf, nbatch = 5e4) {
	set.seed(seed)
	nr = if (mode == 'bip') 2L else if (mode == 'sym') 4L else 3L
	nc = if (mode == 'bip') 3L else nr
	N = 3L
	sq = mode != 'bip'
	rho_uv = 0.8
	sigma_uv = if (mode == 'sym') sqrt(1 - rho_uv^2) else 0.4
	nu = 5; s2 = 0.5; rho = if (mode == 'dir') 0.4 else 0
	offs = c(0, .3, -.2)
	Gm = if (mode == 'bip') matrix(c(0.7, -0.4, 0.5, 0.3)[seq_len(rk * rk)], rk, rk) else NULL
	Ld = if (mode == 'sym') c(0.9, -0.6)[seq_len(rk)] else NULL
	sym_fill = function(a) { for (t in seq_len(dim(a)[3L])) { m = a[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; a[, , t] = m }; a }
	X = array(rnorm(nr * nc * N), c(nr, nc, N)); if (mode == 'sym') X = sym_fill(X)
	ncell = nr * nc * N
	ci = arrayInd(seq_len(ncell), c(nr, nc, N))
	base = array(-.3 + .4 * X, c(nr, nc, N))
	Y = array(rpois(ncell, exp(-.5 + .3 * X)), c(nr, nc, N)); Y[1, , ] = 0; Y[, 1, ] = 0
	if (mode == 'sym') Y = sym_fill(Y)
	if (sq) for (t in 1:N) diag(Y[, , t]) = NA
	Y[sample(which(is.finite(Y)), 1)] = NA
	if (mode == 'sym') Y = sym_fill(Y)
	cx = .lame_ltm_context(Y, 'poisson', mode, offset = offs)
	ar_paths = function(m, n, lam = NULL) {
		e = array(rnorm(m * n * rk * N), c(m, n, rk, N))
		out = array(0, c(m, n, rk, N))
		out[, , , 1] = e[, , , 1] * sigma_uv / sqrt(1 - rho_uv^2)
		for (t in 2:N) out[, , , t] = rho_uv * out[, , , t - 1] + sigma_uv * e[, , , t] / (if (is.null(lam)) 1 else sqrt(rep(lam[, , t], rk)))
		out
	}
	draw = function(m) {
		lam_u = lam_v = NULL
		if (kind == 't') {
			lam_u = array(rgamma(m * nr * N, nu / 2, rate = nu / 2), c(m, nr, N))
			lam_v = array(rgamma(m * nc * N, nu / 2, rate = nu / 2), c(m, nc, N))
		}
		Um = ar_paths(m, nr, lam_u)
		Vm = if (mode == 'sym') NULL else ar_paths(m, nc, lam_v)
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			mult = if (mode == 'sym') rowSums(sapply(seq_len(rk), function(r) Ld[r] * Um[, i, r, t] * Um[, j, r, t]))
			else if (mode == 'bip') Reduce(
				'+',
				lapply(
					seq_len(rk),
					function(r) Reduce('+', lapply(seq_len(rk), function(q) Gm[r, q] * Um[, i, r, t] * Vm[, j, q, t]))
				)
			)
			else rowSums(sapply(seq_len(rk), function(r) Um[, i, r, t] * Vm[, j, r, t]))
			ez[, cc] = base[i, j, t] + mult
		}
		sd_z = sqrt(s2)
		E = matrix(rnorm(m * ncell), m, ncell) * sd_z
		if (mode == 'dir') {
			c1 = sqrt((1 + rho) / 2); c2 = sqrt((1 - rho) / 2)
			for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
				cij = which(
					ci[, 1] == i & ci[, 2] == j & ci[, 3] == t
				); cji = which(ci[, 1] == j & ci[, 2] == i & ci[, 3] == t)
				g1 = rnorm(m); g2 = rnorm(m)
				E[, cij] = sd_z * (c1 * g1 + c2 * g2); E[, cji] = sd_z * (c1 * g1 - c2 * g2)
			}
		}
		if (mode == 'sym') for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr)
			E[, which(
				ci[, 1] == j & ci[, 2] == i & ci[, 3] == t
			)] = E[, which(ci[, 1] == i & ci[, 2] == j & ci[, 3] == t)]
		if (sq) { dgc = which(ci[, 1] == ci[, 2]); E[, dgc] = matrix(
			rnorm(m * length(dgc)),
			m,
			length(dgc)
		) * sqrt((1 + rho) * s2) }
		Zm = ez + E
		lw = rep(0, m)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			if ((sq && i == j) || (mode == 'sym' && i > j)) next
			y = if (mode == 'sym') { yu = Y[i, j, t]; if (is.finite(yu)) yu else Y[j, i, t] } else Y[i, j, t]
			if (is.finite(y)) lw = lw + y * (Zm[, cc] + offs[t]) - exp(Zm[, cc] + offs[t])
		}
		inbox = rep(TRUE, m)
		if (is.finite(rail)) { inbox = apply(abs(Um), 1, max) <= rail; if (!is.null(Vm)) inbox = inbox & (apply(abs(Vm), 1, max) <= rail) }
		list(U = Um, V = Vm, ez = ez, Zm = Zm, w = exp(lw - max(lw)) * inbox, lu = lam_u, lv = lam_v)
	}
	state = function(bt, r) list(
		Z = array(bt$Zm[r, ], c(nr, nc, N)),
		EZ = array(bt$ez[r, ], c(nr, nc, N)),
		U = array(bt$U[r, , , ], c(nr, rk, N)),
		V = if (is.null(bt$V)) NULL else array(bt$V[r, , , ], c(nc, rk, N)),
		lu = if (is.null(bt$lu)) NULL else matrix(bt$lu[r, , ], nr, N),
		lv = if (is.null(bt$lv)) NULL else matrix(bt$lv[r, , ], nc, N)
	)
	stat = function(s) {
		Z = s$Z; EZ = s$EZ; U = s$U; V = if (is.null(s$V)) s$U else s$V
		r12 = Z[1, 2, 1] - EZ[1, 2, 1]; r21 = Z[2, 1, 1] - EZ[2, 1, 1]
		c(
			U[1, 1, 1],
			U[2, 1, 2],
			U[1, 1, 1] * U[1, 1, 2],
			U[2, 1, 1] * U[2, 1, N],
			sum(U^2),
			sum(V^2),
			U[1, 1, N],
			V[2, 1, 1],
			V[1, 1, 2] * V[1, 1, N],
			U[1, 1, 1] * V[1, 1, 1],
			Z[1, 2, 1],
			Z[1, 2, 1]^2,
			Z[2, 1, 2] * Z[1, 2, 2],
			Z[1, 1, 1],
			Z[1, 1, 1]^2,
			mean(Z),
			mean(EZ),
			EZ[1, 2, 1] * EZ[2, 1, 1],
			r12 * r21,
			r12^2,
			EZ[1, 2, 3],
			Z[2, 3, 3],
			mean(EZ[1, , ]),
			mean(EZ[, 1, ]),
			sum(EZ^2),
			EZ[1, 1, 1] * sum(U[1, , 1]^2),
			EZ[2, 2, 2] * sum(U[2, , 2]^2)
		)
	}
	kern = function(s) {
		o = .lame_dpm_visit(
			cx,
			s$Z,
			s$EZ,
			s$U,
			s$V,
			Gm,
			Ld,
			rho_uv,
			sigma_uv,
			max_abs = rail,
			lambda_u = s$lu,
			lambda_v = s$lv
		)
		s$Z = o$Z; s$EZ = o$EZ; s$U = o$U; if (mode != 'sym') s$V = o$V
		s
	}
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		keep = which(bt$w > 1e-12 * max(bt$w))
		if (length(keep) > 4000L) {
			w0 = quantile(bt$w[keep], 1 - 2 * 4000 / length(keep))
			pk = pmin(1, bt$w[keep] / w0); sel = runif(length(keep)) < pk
			wk = (bt$w[keep] / pk)[sel]; keep = keep[sel]
		} else wk = bt$w[keep]
		Dl[[length(Dl) + 1L]] = do.call(
			rbind,
			lapply(keep, function(r) { s0 = state(bt, r); stat(kern(s0)) - stat(s0) })
		)
		Wl[[length(Wl) + 1L]] = wk
	}
	D = do.call(rbind, Dl); W = unlist(Wl); wn = W / sum(W)
	dm = colSums(D * wn); se = sqrt(colSums(wn^2 * sweep(D, 2, dm)^2))
	z = ifelse(se > 0, dm / se, 0)
	list(z = z, maxz = max(abs(z)), ess = sum(W)^2 / sum(W^2))
}

# the same check for a binary panel (.lame_dpm_visit_binary): paths and dyad factors weighted by the probit mass, the
# latent array drawn from its conditional; `rho` is the dyadic correlation (the diagonal variance of a symmetric panel)
dpmt_exact_binary = function(mode, kind = 'ar1', rk = 1L, ndraw = 2e5, seed = 1L, rail = Inf, nbatch = 5e4) {
	set.seed(seed)
	nr = if (mode == 'bip') 2L else if (mode == 'sym') 4L else 3L
	nc = if (mode == 'bip') 3L else nr
	N = 3L
	sq = mode != 'bip'
	rho_uv = 0.8
	sigma_uv = if (mode == 'sym') sqrt(1 - rho_uv^2) else 0.4
	nu = 5; rho = if (mode == 'dir') 0.4 else 0.3
	sdc = if (mode == 'dir') sqrt(1 - abs(rho)) else 1
	Gm = if (mode == 'bip') matrix(c(0.7, -0.4, 0.5, 0.3)[seq_len(rk * rk)], rk, rk) else NULL
	Ld = if (mode == 'sym') c(0.9, -0.6)[seq_len(rk)] else NULL
	sym_fill = function(a) { for (t in seq_len(dim(a)[3L])) { m = a[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; a[, , t] = m }; a }
	X = array(rnorm(nr * nc * N), c(nr, nc, N)); if (mode == 'sym') X = sym_fill(X)
	ncell = nr * nc * N
	ci = arrayInd(seq_len(ncell), c(nr, nc, N))
	base = array(-.3 + .4 * X, c(nr, nc, N))
	Y = array(rbinom(ncell, 1, pnorm(base)), c(nr, nc, N)); Y[1, , ] = 0; Y[, 1, ] = 0
	if (mode == 'sym') Y = sym_fill(Y)
	if (sq) for (t in 1:N) diag(Y[, , t]) = NA
	Y[sample(which(is.finite(Y)), 1)] = NA
	if (mode == 'sym') Y = sym_fill(Y)
	cx = .lame_dpm_binary_context(Y, mode)
	ar_paths = function(m, n, lam = NULL) {
		e = array(rnorm(m * n * rk * N), c(m, n, rk, N))
		out = array(0, c(m, n, rk, N))
		out[, , , 1] = e[, , , 1] * sigma_uv / sqrt(1 - rho_uv^2)
		for (t in 2:N) out[, , , t] = rho_uv * out[, , , t - 1] + sigma_uv * e[, , , t] / (if (is.null(lam)) 1 else sqrt(rep(lam[, , t], rk)))
		out
	}
	draw = function(m) {
		lam_u = lam_v = NULL
		if (kind == 't') {
			lam_u = array(rgamma(m * nr * N, nu / 2, rate = nu / 2), c(m, nr, N))
			lam_v = array(rgamma(m * nc * N, nu / 2, rate = nu / 2), c(m, nc, N))
		}
		Um = ar_paths(m, nr, lam_u)
		Vm = if (mode == 'sym') NULL else ar_paths(m, nc, lam_v)
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			mult = if (mode == 'sym') Reduce(
				'+',
				lapply(seq_len(rk), function(r) Ld[r] * Um[, i, r, t] * Um[, j, r, t])
			)
			else if (mode == 'bip') Reduce(
				'+',
				lapply(
					seq_len(rk),
					function(r) Reduce('+', lapply(seq_len(rk), function(q) Gm[r, q] * Um[, i, r, t] * Vm[, j, q, t]))
				)
			)
			else Reduce('+', lapply(seq_len(rk), function(r) Um[, i, r, t] * Vm[, j, r, t]))
			ez[, cc] = base[i, j, t] + mult
		}
		off = matrix(0, m, ncell)
		if (mode == 'dir') {
			k = sqrt(abs(rho)); sg = if (rho < 0) -1 else 1
			for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
				cij = which(
					ci[, 1] == i & ci[, 2] == j & ci[, 3] == t
				); cji = which(ci[, 1] == j & ci[, 2] == i & ci[, 3] == t)
				f = rnorm(m); off[, cij] = k * f; off[, cji] = sg * k * f
			}
		}
		lw = rep(0, m)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			if ((sq && i == j) || (mode == 'sym' && i > j)) next
			y = if (mode == 'sym') { yu = Y[i, j, t]; if (is.finite(yu)) yu else Y[j, i, t] } else Y[i, j, t]
			if (is.finite(y)) lw = lw + pnorm((2 * y - 1) * (ez[, cc] + off[, cc]) / sdc, log.p = TRUE)
		}
		inbox = rep(TRUE, m)
		if (is.finite(rail)) { inbox = apply(abs(Um), 1, max) <= rail; if (!is.null(Vm)) inbox = inbox & (apply(abs(Vm), 1, max) <= rail) }
		list(U = Um, V = Vm, ez = ez, off = off, w = exp(lw - max(lw)) * inbox, lu = lam_u, lv = lam_v)
	}
	state = function(bt, r) {
		EZ = array(bt$ez[r, ], c(nr, nc, N)); Z0 = array(0, c(nr, nc, N))
		Z = if (mode == 'dir') .lame_binary_z_given_factor(
			Z0,
			EZ,
			array(bt$off[r, ], c(nr, nc, N)),
			rho,
			cx$i1,
			cx$i0,
			cx$ina,
			cx$idg
		)
		else if (mode == 'sym') { for (t in 1:N) Z0[, , t] = rZ_bin_sym_fc(Z0[, , t], EZ[, , t], Y[, , t], rho); Z0 }
		else rZ_bin_bip_batch_cpp(Z0, EZ, Y)
		list(
			Z = Z,
			EZ = EZ,
			U = array(bt$U[r, , , ], c(nr, rk, N)),
			V = if (is.null(bt$V)) NULL else array(bt$V[r, , , ], c(nc, rk, N)),
			lu = if (is.null(bt$lu)) NULL else matrix(bt$lu[r, , ], nr, N),
			lv = if (is.null(bt$lv)) NULL else matrix(bt$lv[r, , ], nc, N)
		)
	}
	stat = function(s) {
		Z = s$Z; EZ = s$EZ; U = s$U; V = if (is.null(s$V)) s$U else s$V
		r12 = Z[1, 2, 1] - EZ[1, 2, 1]; r21 = Z[2, 1, 1] - EZ[2, 1, 1]
		c(
			U[1, 1, 1],
			U[2, 1, 2],
			U[1, 1, 1] * U[1, 1, 2],
			U[2, 1, 1] * U[2, 1, N],
			sum(U^2),
			sum(V^2),
			U[1, 1, N],
			V[2, 1, 1],
			V[1, 1, 2] * V[1, 1, N],
			U[1, 1, 1] * V[1, 1, 1],
			Z[1, 2, 1],
			Z[1, 2, 1]^2,
			Z[2, 1, 2] * Z[1, 2, 2],
			Z[1, 1, 1],
			Z[1, 1, 1]^2,
			mean(Z),
			mean(EZ),
			EZ[1, 2, 1] * EZ[2, 1, 1],
			r12 * r21,
			r12^2,
			EZ[1, 2, 3],
			Z[2, 3, 3],
			mean(EZ[1, , ]),
			mean(EZ[, 1, ]),
			sum(EZ^2),
			mean(Z > 0),
			sum(Z[, , 2] * EZ[, , 2]),
			EZ[1, 1, 1] * sum(U[1, , 1]^2),
			EZ[2, 2, 2] * sum(U[2, , 2]^2)
		)
	}
	kern = function(s) {
		o = .lame_dpm_visit_binary(
			cx,
			s$Z,
			s$EZ,
			Y,
			if (mode == 'bip') 0 else rho,
			s$U,
			s$V,
			Gm,
			Ld,
			rho_uv,
			sigma_uv,
			max_abs = rail,
			lambda_u = s$lu,
			lambda_v = s$lv
		)
		s$Z = o$Z; s$EZ = o$EZ; s$U = o$U; if (mode != 'sym') s$V = o$V
		s
	}
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		keep = which(bt$w > 1e-12 * max(bt$w))
		if (length(keep) > 4000L) {
			w0 = quantile(bt$w[keep], 1 - 2 * 4000 / length(keep))
			pk = pmin(1, bt$w[keep] / w0); sel = runif(length(keep)) < pk
			wk = (bt$w[keep] / pk)[sel]; keep = keep[sel]
		} else wk = bt$w[keep]
		Dl[[length(Dl) + 1L]] = do.call(
			rbind,
			lapply(keep, function(r) { s0 = state(bt, r); stat(kern(s0)) - stat(s0) })
		)
		Wl[[length(Wl) + 1L]] = wk
	}
	D = do.call(rbind, Dl); W = unlist(Wl); wn = W / sum(W)
	dm = colSums(D * wn); se = sqrt(colSums(wn^2 * sweep(D, 2, dm)^2))
	z = ifelse(se > 0, dm / se, 0)
	list(z = z, maxz = max(abs(z)), ess = sum(W)^2 / sum(W^2))
}
