# exact-posterior check of the level step of a directed binary fit with an rw1 intercept path.
# prior states weighted by the probit mass of the observed pairs; one visit keeps every weighted mean.
# variant 'general' has no positions, 'path' rank-one positions (shrink: matrix-t prior or burn-in gaussian);
# blocks 'rw1' adds a static coefficient, 'mixed' an ar1 one; `bad_dim` mis-states the block dimension
lvl_exact = function(
	variant = 'general',
	blocks = 'rw1',
	rho_d = 0.4,
	shrink = TRUE,
	ndraw = 2e5,
	seed = 1L,
	tries = 3L,
	nbatch = 1e5L,
	bad_dim = 0L
) {
	set.seed(seed)
	n = 3L; N = 2L; R = 1L
	Sab = matrix(c(1, .3, .3, .8), 2); rho_ab = .7; sigma_ab = .5
	v0 = 2; sb = .6
	g = 3; s2 = 1
	v0x = 1.5; rho_x = .6; sx = .4; lam_x = 1.3
	df = 2L * R + 2L; psi = df
	X = array(rnorm(n * n * N), c(n, n, N))
	Y = array(rbinom(n * n * N, 1, .45), c(n, n, N))
	for (t in 1:N) diag(Y[, , t]) = NA
	dg = array(FALSE, c(n, n, N))
	for (t in 1:N) diag(dg[, , t]) = TRUE
	i1 = which(Y == 1 & !dg)
	i0 = which(Y == 0 & !dg)
	ina = which(is.na(Y) & !dg)
	idg = which(dg)
	obs = c(i1, i0); sY = c(rep(1, length(i1)), rep(-1, length(i0)))
	upper = upper.tri(matrix(0, n, n))
	pr = which(upper, arr.ind = TRUE)
	# the pair cells, period by period
	c1 = unlist(lapply(1:N, function(t) pr[, 1] + (pr[, 2] - 1) * n + (t - 1) * n * n))
	c2 = unlist(lapply(1:N, function(t) pr[, 2] + (pr[, 1] - 1) * n + (t - 1) * n * n))
	s1 = 2 * Y[c1] - 1; s2c = 2 * Y[c2] - 1
	cs = chol(Sab)
	draw = function(m) {
		lv = array(rnorm(m * n * 2), c(m, n, 2))
		for (k in 1:m) lv[k, , ] = lv[k, , ] %*% cs
		a = lv[, , 1] - rowMeans(lv[, , 1]); b = lv[, , 2] - rowMeans(lv[, , 2])
		mk = function() {
			d = array(0, c(m, n, N))
			d[, , 1] = matrix(rnorm(m * n), m, n) * sigma_ab / sqrt(1 - rho_ab^2)
			for (t in 2:N) d[, , t] = rho_ab * d[, , t - 1] + matrix(rnorm(m * n), m, n) * sigma_ab
			for (t in 1:N) d[, , t] = d[, , t] - rowMeans(d[, , t])
			d
		}
		da = mk(); db = mk()
		ic = matrix(0, m, N)
		ic[, 1] = rnorm(m, 0, sqrt(v0))
		for (t in 2:N) ic[, t] = ic[, t - 1] + rnorm(m, 0, sb)
		xc = if (blocks == 'rw1') matrix(rnorm(m, 0, sqrt(g * s2)), m, N) else {
			x = matrix(0, m, N); x[, 1] = rnorm(m, 0, sqrt(v0x))
			for (t in 2:N) x[, t] = rho_x * x[, t - 1] + rnorm(m, 0, sx * sqrt(lam_x))
			x
		}
		U = V = NULL
		if (variant == 'path') {
			W = array(0, c(m, n, 2))
			for (k in 1:m) {
				Sig = if (shrink) solve(stats::rWishart(1, df, diag(2) / psi)[, , 1]) else diag(2) * n
				W[k, , ] = matrix(rnorm(n * 2), n, 2) %*% chol(Sig)
			}
			U = W[, , 1]; V = W[, , 2]
		}
		ez = array(0, c(m, n * n * N))
		for (t in 1:N) for (j in 1:n) for (i in 1:n) {
			cc = i + (j - 1) * n + (t - 1) * n * n
			ez[, cc] = ic[, t] + xc[, t] * X[i, j, t] + a[, i] + da[, i, t] + b[, j] + db[, j, t] +
				(if (variant == 'path') U[, i] * V[, j] else 0)
		}
		lw = rowSums(log(pmax(pbivnorm::pbivnorm(
			as.vector(t(t(ez[, c1, drop = FALSE]) * s1)),
			as.vector(t(t(ez[, c2, drop = FALSE]) * s2c)),
			rep(s1 * s2c * rho_d, each = m)
		), 1e-300)) |> matrix(m))
		list(
			a = a,
			b = b,
			da = da,
			db = db,
			ic = ic,
			xc = xc,
			U = U,
			V = V,
			ez = ez,
			w = exp(lw - max(lw))
		)
	}
	# draw the latent pairs given the state and the data
	draw_z = function(ez) {
		Z = array(ez, c(n, n, N))
		for (k in seq_along(c1)) {
			e1 = ez[c1[k]]; e2 = ez[c2[k]]
			hit = FALSE
			for (tr in 1:2000) {
				z = rnorm(2); z2 = rho_d * z[1] + sqrt(1 - rho_d^2) * z[2]
				z1 = e1 + z[1]; z2 = e2 + z2
				if ((z1 > 0) == (s1[k] > 0) && (z2 > 0) == (s2c[k] > 0)) { hit = TRUE; break }
			}
			# rare pair: gibbs sweeps of the truncated bivariate normal from a feasible point
			if (!hit) {
				tn = function(m, sd, pos) {
					u = if (pos) runif(1, pnorm(0, m, sd), 1) else runif(1, 0, pnorm(0, m, sd))
					qnorm(min(max(u, 1e-300), 1 - 1e-16), m, sd)
				}
				z1 = if (s1[k] > 0) abs(e1) + .1 else -abs(e1) - .1
				z2 = if (s2c[k] > 0) abs(e2) + .1 else -abs(e2) - .1
				cs1 = sqrt(1 - rho_d^2)
				for (sw in 1:300) {
					z1 = tn(e1 + rho_d * (z2 - e2), cs1, s1[k] > 0)
					z2 = tn(e2 + rho_d * (z1 - e1), cs1, s2c[k] > 0)
				}
			}
			Z[c1[k]] = z1; Z[c2[k]] = z2
		}
		Z
	}
	state = function(bt, r) list(
		a = bt$a[r, ],
		b = bt$b[r, ],
		da = bt$da[r, , ],
		db = bt$db[r, , ],
		ic = bt$ic[r, ],
		xc = bt$xc[r, ],
		U = if (variant == 'path') bt$U[r, ] else NULL,
		V = if (variant == 'path') bt$V[r, ] else NULL,
		ez = array(bt$ez[r, ], c(n, n, N)),
		Z = draw_z(bt$ez[r, ])
	)
	stat = function(s) {
		c(
			s$ic, s$ic[1]^2, s$ic[2] - s$ic[1], (s$ic[2] - s$ic[1])^2, s$xc, s$xc[1]^2, s$ic[1] * s$xc[1],
			s$a[1], sum(s$a^2), s$b[2], sum(s$b^2), s$a[1] * s$b[1],
			s$da[1, 1], sum(s$da^2), sum(s$db^2), s$da[1, 1] * s$da[1, 2],
			if (variant == 'path') c(s$U[1] * s$V[2], sum(s$U^2), sum(s$V^2)) else NULL,
			s$ez[1, 2, 1], s$ez[1, 2, 1]^2, mean(s$ez), sum(s$ez^2),
			s$Z[1, 2, 1], s$Z[2, 1, 1] * s$Z[1, 2, 1], s$Z[2, 3, 2], s$Z[2, 3, 2]^2
		)
	}
	kern = function(s) {
		EZ = s$ez; Z = s$Z
		off = .lame_dyad_factor_draw(Z, EZ, rho_d, upper)
		sdc = sqrt(1 - abs(rho_d)); ez_o = EZ[obs] / sdc; off_o = off[obs] / sdc
		ll_cur = sum(pnorm(sY * (ez_o + off_o), log.p = TRUE))
		iS = solve(Sab); AB = cbind(s$a, s$b)
		q_lev = sum(AB * (AB %*% iS)); n_lev = 2L * n
		dss = .lame_ab_dev_ss(s$da, s$db, rho_ab); q_dev = dss$ss / sigma_ab^2; n_dev = dss$cnt
		path = if (blocks == 'rw1') matrix(s$ic, N, 1) else cbind(s$ic, s$xc)
		pt = .lame_path_scale_terms(
			if (blocks == 'rw1') s$xc[1] else numeric(0),
			path,
			if (blocks == 'rw1') 1 else c(1, rho_x),
			if (blocks == 'rw1') sb else c(sb, sx),
			if (blocks == 'rw1') diag(1, 1) else diag(c(1, lam_x)),
			if (blocks == 'rw1') diag(v0, 1) else diag(c(v0, v0x)),
			if (blocks == 'rw1') g else numeric(0),
			s2
		)
		pt$dimension = pt$dimension + bad_dim
		if (variant == 'path') {
			uv_eigen = eigen(crossprod(cbind(s$U, s$V)), symmetric = TRUE, only.values = TRUE)$values
			uv_eigen = pmax(uv_eigen, 0)
		}
		ctot = 1
		for (k in seq_len(tries)) {
			cc = ctot * exp(rnorm(1, 0, .3))
			ll_new = sum(pnorm(sY * (cc * ez_o + off_o), log.p = TRUE))
			lr = ll_new - ll_cur + if (variant == 'path') {
				.lame_static_scale_logratio(
					cc / ctot,
					ctot^2 * (pt$q + q_lev + q_dev),
					pt$dimension + n_lev - 2L + max(n_dev - 2L * N, 0L),
					ctot * uv_eigen,
					n,
					R,
					shrink,
					df,
					psi
				)
			} else {
				.lame_level_coef_logratio(
					cc / ctot,
					ctot,
					numeric(0),
					matrix(0, 0, 0),
					pt,
					q_lev,
					n_lev - 2L,
					q_dev,
					max(n_dev - 2L * N, 0L)
				)
			}
			if (is.finite(lr) && log(runif(1)) < lr) { ctot = cc; ll_cur = ll_new }
		}
		if (ctot != 1) {
			s$ic = ctot * s$ic; s$xc = ctot * s$xc
			s$a = ctot * s$a; s$b = ctot * s$b; s$da = ctot * s$da; s$db = ctot * s$db
			if (variant == 'path') { s$U = sqrt(ctot) * s$U; s$V = sqrt(ctot) * s$V }
			EZ = ctot * EZ
			s$Z = .lame_binary_z_given_factor(Z, EZ, off, rho_d, i1, i0, ina, idg)
			s$ez = EZ
		}
		s
	}
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		ok = which(is.finite(bt$w))
		keep = ok[bt$w[ok] > 1e-12 * max(bt$w[ok])]
		if (length(keep) > 6000L) {
			w0 = quantile(bt$w[keep], max(0, 1 - 2 * 6000 / length(keep)))
			pk = pmin(1, bt$w[keep] / w0); sel = runif(length(keep)) < pk
			wk = (bt$w[keep] / pk)[sel]; keep = keep[sel]
		} else wk = bt$w[keep]
		Dl[[length(Dl) + 1L]] = do.call(rbind, lapply(keep, function(r) {
			s0 = state(bt, r)
			stat(kern(s0)) - stat(s0)
		}))
		Wl[[length(Wl) + 1L]] = wk
	}
	D = do.call(rbind, Dl); W = unlist(Wl); wn = W / sum(W)
	dm = colSums(D * wn); se = sqrt(colSums(wn^2 * sweep(D, 2, dm)^2))
	z = ifelse(se > 0, dm / se, 0)
	list(z = z, maxz = max(abs(z)), ess = sum(W)^2 / sum(W^2), n = length(W))
}
