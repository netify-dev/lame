# exact-posterior check of the predictor scale step with an rw1 intercept path, as bplt_exact.
# blocks 'rw1' adds a static coefficient, 'mixed' an ar1 one; `bad_dim` mis-states the path dimension
bplt_path_exact = function(
	blocks = 'rw1',
	ndraw = 3e5,
	seed = 1L,
	tries = 3L,
	nbatch = 1e5L,
	bad_dim = 0L
) {
	set.seed(seed)
	nA = 3L; nB = 4L; N = 2L
	g_prior = 3; s2 = 1; va = 1.2; vb = 0.8; rho_ab = 0.7; sigma_ab = 0.5
	v0 = 2; sb = .6; v0x = 1.5; rho_x = .6; sx = .4; lam_x = 1.3
	X = array(rnorm(nA * nB * N), c(nA, nB, N))
	u = rnorm(nA); v = rnorm(nB)
	Y = array(rbinom(nA * nB * N, 1, pnorm(-.3 + .5 * X)), c(nA, nB, N))
	Y[1, , ] = 0
	Y[, 1, ] = 0
	Y[sample(length(Y), 2)] = NA
	cx = .lame_bpl_context(Y)
	ncell = nA * nB * N
	ci = arrayInd(seq_len(ncell), c(nA, nB, N))
	draw = function(m) {
		ic = matrix(0, m, N)
		ic[, 1] = rnorm(m, 0, sqrt(v0))
		for (t in 2:N) ic[, t] = ic[, t - 1] + rnorm(m, 0, sb)
		xc = if (blocks == 'rw1') matrix(rnorm(m, 0, sqrt(g_prior * s2)), m, N) else {
			x = matrix(0, m, N); x[, 1] = rnorm(m, 0, sqrt(v0x))
			for (t in 2:N) x[, t] = rho_x * x[, t - 1] + rnorm(m, 0, sx * sqrt(lam_x))
			x
		}
		a = matrix(rnorm(m * nA, 0, sqrt(va)), m, nA); a = a - rowMeans(a)
		b = matrix(rnorm(m * nB, 0, sqrt(vb)), m, nB); b = b - rowMeans(b)
		G = rnorm(m)
		mk = function(n) {
			d = array(0, c(m, n, N))
			d[, , 1] = matrix(rnorm(m * n), m, n) * sigma_ab / sqrt(1 - rho_ab^2)
			for (t in 2:N) d[, , t] = rho_ab * d[, , t - 1] + matrix(rnorm(m * n), m, n) * sigma_ab
			for (t in 1:N) d[, , t] = d[, , t] - rowMeans(d[, , t])
			d
		}
		da = mk(nA); db = mk(nB)
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			ez[, cc] = ic[, t] + xc[, t] * X[i, j, t] + a[, i] + b[, j] +
				da[, i, t] + db[, j, t] + G * u[i] * v[j]
		}
		lw = rep(0, m)
		for (cc in seq_len(ncell)) {
			y = Y[ci[cc, 1], ci[cc, 2], ci[cc, 3]]
			if (is.finite(y)) lw = lw + pnorm((2 * y - 1) * ez[, cc], log.p = TRUE)
		}
		list(ic = ic, xc = xc, a = a, b = b, G = G, da = da, db = db, ez = ez, w = exp(lw - max(lw)))
	}
	state = function(bt, r) list(
		ic = bt$ic[r, ],
		xc = bt$xc[r, ],
		a = bt$a[r, ],
		b = bt$b[r, ],
		G = bt$G[r],
		da = bt$da[r, , ],
		db = bt$db[r, , ],
		ez = array(bt$ez[r, ], c(nA, nB, N))
	)
	stat = function(s) c(
		s$ic, s$ic[1]^2, (s$ic[2] - s$ic[1])^2, s$xc, s$xc[1]^2, s$ic[1] * s$xc[1],
		s$a[1], sum(s$a^2), s$b[1], sum(s$b^2), s$G, s$G^2, s$G * s$ic[2],
		s$da[1, 1], sum(s$da^2), sum(s$db^2),
		s$ez[1, 1, 1], s$ez[1, 1, 1]^2, mean(s$ez), sum(s$ez^2), mean(s$ez[2, , ])
	)
	kern = function(s) {
		one = 1 - min(rho_ab^2, 0.9801)
		ss = one * (sum(s$da[, 1]^2) + sum(s$db[, 1]^2)) +
			sum((s$da[, 2:N, drop = FALSE] - rho_ab * s$da[, 1:(N - 1), drop = FALSE])^2) +
			sum((s$db[, 2:N, drop = FALSE] - rho_ab * s$db[, 1:(N - 1), drop = FALSE])^2)
		path = if (blocks == 'rw1') matrix(s$ic, N, 1) else cbind(s$ic, s$xc)
		pt = .lame_path_scale_terms(
			if (blocks == 'rw1') s$xc[1] else numeric(0),
			path,
			if (blocks == 'rw1') 1 else c(1, rho_x),
			if (blocks == 'rw1') sb else c(sb, sx),
			if (blocks == 'rw1') diag(1, 1) else diag(c(1, lam_x)),
			if (blocks == 'rw1') diag(v0, 1) else diag(c(v0, v0x)),
			if (blocks == 'rw1') g_prior else numeric(0),
			s2
		)
		pt$dimension = pt$dimension + bad_dim
		ca = .lame_bpl_coef_args(numeric(0), matrix(0, 0, 0), s$G^2, 1L, pt)
		o = .lame_bpl_scale(
			s$ez[cx$idx],
			cx$sY,
			ca$beta,
			ca$iV0,
			sum(s$a^2) / va + sum(s$b^2) / vb,
			nA + nB - 2L,
			ss / sigma_ab^2,
			(nA + nB) * N - 2L * N,
			ca$q_g,
			ca$n_g,
			tries = tries,
			sd = 0.3,
			adapt = FALSE
		)
		if (o$c != 1) {
			for (nm in c('ic', 'xc', 'a', 'b', 'G', 'da', 'db', 'ez')) s[[nm]] = o$c * s[[nm]]
		}
		s
	}
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		ok = which(is.finite(bt$w))
		keep = ok[bt$w[ok] > 1e-12 * max(bt$w[ok])]
		if (length(keep) > 8000L) {
			w0 = quantile(bt$w[keep], max(0, 1 - 2 * 8000 / length(keep)))
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
