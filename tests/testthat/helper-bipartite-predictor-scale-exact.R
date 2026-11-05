# exact-posterior check of the predictor scale move of a binary bipartite fit on a tiny panel: prior draws weighted
# by the probit mass of the observed cells; one visit must leave every weighted mean unchanged. returns z scores of
# the before/after differences and the effective size. variant: 'static', 'dyn' (dynamic_ab) or 'nog' (no G)
bplt_exact = function(variant = 'static', ndraw = 1e5, seed = 1L, tries = 3L, gmax = Inf, nbatch = 1e5L) {
	set.seed(seed)
	nA = 3L; nB = 4L; N = 2L
	dyn = variant == 'dyn'; hasG = variant != 'nog'
	g_prior = 3; s2 = 1; va = 1.2; vb = 0.8; rho_ab = 0.7; sigma_ab = 0.5; p = 2L
	X = array(rnorm(nA * nB * N), c(nA, nB, N))
	u = rnorm(nA); v = rnorm(nB)
	Y = array(rbinom(nA * nB * N, 1, pnorm(-.3 + .5 * X)), c(nA, nB, N)); Y[1, , ] = 0; Y[, 1, ] = 0
	Y[sample(length(Y), 2)] = NA
	cx = .lame_bpl_context(Y)
	iV0 = diag(p) / (g_prior * s2)
	ncell = nA * nB * N
	ci = arrayInd(seq_len(ncell), c(nA, nB, N))
	draw = function(m) {
		B = matrix(rnorm(m * p), m, p) * sqrt(g_prior * s2)
		a = matrix(rnorm(m * nA, 0, sqrt(va)), m, nA); a = a - rowMeans(a)
		b = matrix(rnorm(m * nB, 0, sqrt(vb)), m, nB); b = b - rowMeans(b)
		G = if (hasG) rnorm(m) else rep(0, m)
		if (is.finite(gmax)) G = ifelse(abs(G) <= gmax, G, NA)
		da = db = NULL
		if (dyn) {
			mk = function(n) {
				d = array(0, c(m, n, N)); d[, , 1] = matrix(rnorm(m * n), m, n) * sigma_ab / sqrt(1 - rho_ab^2)
				for (t in 2:N) d[, , t] = rho_ab * d[, , t - 1] + matrix(rnorm(m * n), m, n) * sigma_ab
				for (t in 1:N) d[, , t] = d[, , t] - rowMeans(d[, , t])
				d
			}
			da = mk(nA); db = mk(nB)
		}
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) {
			i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
			ez[, cc] = B[, 1] + B[, 2] * X[i, j, t] + a[, i] + b[, j] + (if (dyn) da[, i, t] + db[, j, t] else 0) + G * u[i] * v[j]
		}
		lw = rep(0, m)
		for (cc in seq_len(ncell)) {
			y = Y[ci[cc, 1], ci[cc, 2], ci[cc, 3]]
			if (is.finite(y)) lw = lw + pnorm((2 * y - 1) * ez[, cc], log.p = TRUE)
		}
		lw[is.na(lw)] = -Inf
		list(B = B, a = a, b = b, G = G, da = da, db = db, ez = ez, w = exp(lw - max(lw[is.finite(lw)])))
	}
	state = function(bt, r) list(
		B = bt$B[r, ],
		a = bt$a[r, ],
		b = bt$b[r, ],
		G = bt$G[r],
		da = if (dyn) bt$da[r, , ] else NULL,
		db = if (dyn) bt$db[r, , ] else NULL,
		ez = array(bt$ez[r, ], c(nA, nB, N))
	)
	stat = function(s) {
		ez = s$ez
		c(
			s$B[1],
			s$B[2],
			s$B[1]^2,
			s$B[1] * s$B[2],
			s$a[1],
			s$a[1]^2,
			sum(s$a^2),
			s$b[1],
			sum(s$b^2),
			s$a[1] * s$b[1],
			s$G,
			s$G^2,
			s$B[1] * s$G,
			s$B[2] * s$a[1],
			ez[1, 1, 1],
			ez[1, 1, 1]^2,
			ez[2, 3, 2] * ez[1, 1, 1],
			mean(ez),
			sum(ez^2),
			s$G * s$a[2],
			s$G * s$b[3],
			mean(ez[2, , ]),
			mean(ez[, 2, ]),
			if (dyn) c(s$da[1, 1], sum(s$da^2), s$db[2, 1] * s$db[2, 2], sum(s$db^2), s$da[1, 1] * s$a[1]) else NULL
		)
	}
	kern = function(s) {
		q_dev = 0; n_dev = 0L
		if (dyn) {
			one = 1 - min(rho_ab^2, 0.9801)
			ss = one * (sum(s$da[, 1]^2) + sum(s$db[, 1]^2)) + sum(
				(s$da[, 2:N, drop = FALSE] - rho_ab * s$da[, 1:(N - 1), drop = FALSE])^2
			) +
				sum((s$db[, 2:N, drop = FALSE] - rho_ab * s$db[, 1:(N - 1), drop = FALSE])^2)
			q_dev = ss / sigma_ab^2; n_dev = (nA + nB) * N - 2L * N
		}
		o = .lame_bpl_scale(
			s$ez[cx$idx],
			cx$sY,
			s$B,
			iV0,
			sum(s$a^2) / va + sum(s$b^2) / vb,
			nA + nB - 2L,
			q_dev,
			n_dev,
			if (hasG) s$G^2 else 0,
			if (hasG) 1L else 0L,
			room = if (hasG && is.finite(gmax) && abs(s$G) > 0) gmax / abs(s$G) else Inf,
			tries = tries,
			sd = 0.3,
			adapt = FALSE
		)
		if (o$c != 1) {
			s$B = o$c * s$B; s$a = o$c * s$a; s$b = o$c * s$b
			if (hasG) s$G = o$c * s$G
			if (dyn) { s$da = o$c * s$da; s$db = o$c * s$db }
			s$ez = o$c * s$ez
		}
		s
	}
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		ok = which(is.finite(bt$w) & !is.na(bt$G))
		keep = ok[bt$w[ok] > 1e-12 * max(bt$w[ok])]
		if (length(keep) > 8000L) {
			w0 = quantile(bt$w[keep], max(0, 1 - 2 * 8000 / length(keep)))
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
	list(z = z, maxz = max(abs(z)), ess = sum(W)^2 / sum(W^2), n = length(W))
}
