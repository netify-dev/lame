# composition of two moves of one Poisson sweep with moving positions and coefficients: the position-path move
# (R/dynamic_position_move.R) and the coefficient-path translations (R/poisson_period_move.R). on a tiny panel prior draws
# are weighted by the Poisson mass of the observed cells, and the two moves in sequence must leave every weighted mean
# unchanged; a second move that reads the predictor from before the first one fails the check

ppm_dpm_exact = function(mode, ord, ndraw, seed, asym = FALSE, nbatch = 5e4) {
	ns = asNamespace('lame')
	ltm_context = ns$.lame_ltm_context; ppm_context = ns$.lame_ppm_context
	gap = ns$.lame_beta_gap_transition; sym_cols = ns$.lame_ppm_sym_cols
	rk = 1L
	set.seed(seed)
	nr = if (mode == "bip") 2L else if (mode == "sym") 4L else 3L
	nc = if (mode == "bip") 3L else nr
	N = 3L
	sq = mode != "bip"
	rho_uv = 0.8
	sigma_uv = if (mode == "sym") sqrt(1 - rho_uv^2) else 0.4
	s2 = 0.5; rho = if (mode == "dir") 0.4 else 0
	offs = c(0, .3, -.2)
	Gm = if (mode == "bip") matrix(0.7, 1, 1) else NULL
	Ld = if (mode == "sym") 0.9 else NULL
	# coefficient paths for an intercept and a covariate, AR(1)
	pc = 2L; sigma_b = c(0.25, 0.15); lam_b = c(1, 0.8); m0 = c(-1.2, 0.3); v0 = c(0.5, 0.4); rho_b = 0.8; tpos = seq_len(N)
	tr = gap(rho_b, diff(tpos))
	sym_fill = function(a) { for (t in seq_len(dim(a)[3L])) { m = a[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; a[, , t] = m }; a }
	X = array(rnorm(nr * nc * N), c(nr, nc, N)); if (mode == "sym" && !asym) X = sym_fill(X)
	Xd = lapply(seq_len(N), function(t) cbind(1, as.numeric(X[, , t])))
	Xdz = if (mode == "sym") lapply(Xd, function(m) sym_cols(m, nr)) else Xd
	ncell = nr * nc * N
	ci = arrayInd(seq_len(ncell), c(nr, nc, N))
	Y = array(rpois(ncell, exp(-1.2 + .3 * X)), c(nr, nc, N)); Y[1, , ] = 0; Y[, 1, ] = 0
	if (mode == "sym") Y = sym_fill(Y)
	if (sq) for (t in 1:N) diag(Y[, , t]) = NA
	Y[sample(which(is.finite(Y)), 1)] = NA
	if (mode == "sym") Y = sym_fill(Y)
	cx = ltm_context(Y, "poisson", mode, offset = offs)
	cxp = ppm_context(Y, mode, offset = offs, ltm = cx)
	pr = list(
		kind = "ar1",
		rho = rep(rho_b, pc),
		sigma = sigma_b,
		lam = lam_b,
		m0 = m0,
		v0 = v0,
		time_positions = tpos,
		length_scale = NULL
	)
	ar_paths = function(m, n) {
		e = array(rnorm(m * n * rk * N), c(m, n, rk, N))
		out = array(0, c(m, n, rk, N))
		out[, , , 1] = e[, , , 1] * sigma_uv / sqrt(1 - rho_uv^2)
		for (t in 2:N) out[, , , t] = rho_uv * out[, , , t - 1] + sigma_uv * e[, , , t]
		out
	}
	cell_of = function(i, j, t) which(ci[, 1] == i & ci[, 2] == j & ci[, 3] == t)
	draw = function(m) {
		Um = ar_paths(m, nr)
		Vm = if (mode == "sym") NULL else ar_paths(m, nc)
		bp = array(0, c(m, N, pc))
		for (k in seq_len(pc)) {
			bp[, 1, k] = rnorm(m, m0[k], sqrt(v0[k]))
			for (t in 2:N) bp[, t, k] = tr$f[t - 1L] * bp[, t - 1L, k] + rnorm(
				m,
				0,
				sigma_b[k] * sqrt(lam_b[k] * tr$q[t - 1L])
			)
		}
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			mult = if (mode == "sym") Ld * Um[, i, 1, t] * Um[, j, 1, t] else if (mode == "bip") Gm[1, 1] * Um[, i, 1, t] * Vm[, j, 1, t] else Um[, i, 1, t] * Vm[, j, 1, t]
			ez[, cc] = bp[, t, 1] * Xd[[t]][(j - 1L) * nr + i, 1] + bp[, t, 2] * Xd[[t]][(j - 1L) * nr + i, 2] + mult
		}
		mu = ez
		if (mode == "sym") for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
			cu = cell_of(i, j, t); cl = cell_of(j, i, t)
			mu[, cu] = mu[, cl] = (ez[, cu] + ez[, cl]) / 2
		}
		sd_z = sqrt(s2)
		E = matrix(rnorm(m * ncell), m, ncell) * sd_z
		if (mode == "dir") {
			c1 = sqrt((1 + rho) / 2); c2 = sqrt((1 - rho) / 2)
			for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr) {
				cij = cell_of(i, j, t); cji = cell_of(j, i, t)
				g1 = rnorm(m); g2 = rnorm(m)
				E[, cij] = sd_z * (c1 * g1 + c2 * g2); E[, cji] = sd_z * (c1 * g1 - c2 * g2)
			}
		}
		if (mode == "sym") for (t in 1:N) for (i in 1:(nr - 1)) for (j in (i + 1):nr) E[, cell_of(j, i, t)] = E[, cell_of(i, j, t)]
		if (sq) { dgc = which(ci[, 1] == ci[, 2]); E[, dgc] = matrix(
			rnorm(m * length(dgc)),
			m,
			length(dgc)
		) * sqrt((1 + rho) * s2) }
		Zm = mu + E
		lw = rep(0, m)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			if ((sq && i == j) || (mode == "sym" && i > j)) next
			y = if (mode == "sym") { yu = Y[i, j, t]; if (is.finite(yu)) yu else Y[j, i, t] } else Y[i, j, t]
			if (is.finite(y)) lw = lw + y * (Zm[, cc] + offs[t]) - exp(Zm[, cc] + offs[t])
		}
		list(U = Um, V = Vm, bp = bp, ez = ez, Zm = Zm, w = exp(lw - max(lw)))
	}
	state = function(bt, r) list(
		Z = array(bt$Zm[r, ], c(nr, nc, N)),
		EZ = array(bt$ez[r, ], c(nr, nc, N)),
		U = array(bt$U[r, , , ], c(nr, rk, N)),
		V = if (is.null(bt$V)) NULL else array(bt$V[r, , , ], c(nc, rk, N)),
		path = matrix(bt$bp[r, , ], N, pc)
	)
	bar = function(a) { if (mode == "sym") for (t in seq_len(dim(a)[3L])) a[, , t] = (a[, , t] + t(a[, , t])) / 2; a }
	stat = function(s) {
		Z = s$Z; EZ = bar(s$EZ); U = s$U; V = if (is.null(s$V)) s$U else s$V; P = s$path
		r12 = Z[1, 2, 1] - EZ[1, 2, 1]; r21 = Z[2, 1, 1] - EZ[2, 1, 1]
		c(
			U[1, 1, 1],
			U[2, 1, 2],
			U[1, 1, 1] * U[1, 1, 2],
			sum(U^2),
			sum(V^2),
			U[1, 1, N],
			V[2, 1, 1],
			U[1, 1, 1] * V[1, 1, 1],
			P[1, 1],
			P[N, 1],
			P[2, 2],
			P[1, 1] * P[N, 1],
			P[2, 1] * P[2, 2],
			P[1, 1]^2,
			sum(P[, 2]^2),
			P[1, 1] * U[1, 1, 1],
			P[2, 2] * sum(U[, 1, 2]^2),
			Z[1, 2, 1],
			Z[1, 2, 1]^2,
			Z[2, 1, 2] * Z[1, 2, 2],
			Z[1, 1, 1],
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
			EZ[1, 1, 1] * P[2, 1],
			Z[1, 2, 2] * P[N, 2]
		)
	}
	k_dpm = function(s) {
		o = ns$.lame_dpm_visit(cx, s$Z, s$EZ, s$U, s$V, Gm, Ld, rho_uv, sigma_uv)
		s$Z = o$Z; s$EZ = o$EZ; s$U = o$U; if (mode != "sym") s$V = o$V
		s
	}
	k_ppm = function(s, EZ_in = s$EZ) {
		o = ns$.lame_ppm_visit(cxp, s$Z, EZ_in, s$path, Xd, pr, Xz = Xdz)
		s$Z = o$Z; s$EZ = o$EZ; s$path = o$path
		s
	}
	kern = switch(
		ord,
		dpm = k_dpm,
		ppm = k_ppm,
		dpm_ppm = function(s) k_ppm(k_dpm(s)),
		ppm_dpm = function(s) k_dpm(k_ppm(s)),
		stale = function(s) { e0 = s$EZ; s1 = k_dpm(s); k_ppm(s1, EZ_in = e0) },
		stop("unknown order")
	)
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

test_that('the position-path move and the translations of the coefficient paths compose exactly, in either order, on bipartite, directed and symmetric panels', {
	skip_on_cran()
	for (cfg in list(
		list('bip', 'dpm_ppm', 1, FALSE),
		list('bip', 'ppm_dpm', 2, FALSE),
		list('dir', 'dpm_ppm', 3, FALSE),
		list('sym', 'ppm_dpm', 4, TRUE)
	)) {
		r = ppm_dpm_exact(cfg[[1]], cfg[[2]], 6e4, cfg[[3]], asym = cfg[[4]])
		expect_gt(r$ess, 2500)
		expect_lt(r$maxz, 4.5)
	}
})

test_that('the composition check fails when the second move reads the predictor from before the first', {
	skip_on_cran()
	r = ppm_dpm_exact('bip', 'stale', 6e4, 5)
	expect_gt(r$maxz, 6)
})
