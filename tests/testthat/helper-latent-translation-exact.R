# exact-posterior check of the translation kernels on a tiny panel: prior draws of (coefficients, levels, latent array)
# weighted by the data factor, so one application of a kernel must leave every weighted mean unchanged; returns z scores
# and the effective size. family: 'poisson', 'ordinal' (rank likelihood) or 'ordinal_alpha' (cutpoints 0 and .8);
# mode: 'dir', 'sym' or 'bip'; comp: intercept with levels centred (free = FALSE) or free; kernel: 'visit',
# 'location', 'coef', 'scale' or 'uv' (static positions, Poisson)
ltmt_exact = function(family, mode, comp, free, kernel = 'visit', ndraw = 4e5, seed = 1L, hook = NULL, nbatch = 1e5) {
	set.seed(seed)
	fam = if (family == 'ordinal_alpha') 'ordinal' else family
	nr = if (mode == 'bip') 2L else if (mode == 'sym') 4L else 3L
	nc = if (mode == 'bip') 3L else nr
	N = if (fam == 'poisson') 2L else 1L
	sq = mode != 'bip'
	Sab = if (mode == 'bip') diag(c(1.1, .8)) else if (mode == 'sym') matrix(1.1, 2, 2) else matrix(c(1, .3, .3, 1.4), 2)
	rho = if (mode == 'dir') .4 else 0
	s2 = if (fam == 'poisson') .5 else 1
	s20 = 1; V0 = 4; p = if (comp) 2L else 1L; iV0 = diag(p) / V0
	Suv_t = matrix(c(1, .4, .4, 1.3), 2); Ld_t = 0.8; tau_t = 1.3; G_t = 0.7
	offs = if (fam == 'poisson') c(0, .3)[seq_len(N)] else rep(0, N)
	alpha = if (family == 'ordinal_alpha') c(0, .8) else NULL
	sym_fill = function(a) { for (t in seq_len(dim(a)[3L])) { m = a[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; a[, , t] = m }; a }
	X = array(rnorm(nr * nc * N), c(nr, nc, N)); if (mode == 'sym') X = sym_fill(X)
	ncell = nr * nc * N
	ci = arrayInd(seq_len(ncell), c(nr, nc, N))
	offd = if (sq) ci[, 1] != ci[, 2] else rep(TRUE, ncell)
	if (fam == 'poisson') {
		Y = array(rpois(ncell, exp(-.5 + .3 * X)), c(nr, nc, N)); Y[1, , ] = 0; Y[, 1, ] = 0
		if (mode == 'sym') Y = sym_fill(Y)
	} else {
		repeat {
			Y = array(NA_real_, c(nr, nc, N)); Y[which(offd)] = sample(rep(1:3, length.out = sum(offd)))
			if (mode == 'sym') Y = sym_fill(Y)
			if (length(unique(Y[is.finite(Y)])) == 3L) break
		}
	}
	if (sq) for (t in 1:N) diag(Y[, , t]) = NA
	if (fam == 'poisson') { Y[sample(which(is.finite(Y)), 1)] = NA; if (mode == 'sym') Y = sym_fill(Y) }
	cx = .lame_ltm_context(
		Y,
		fam,
		mode,
		offset = if (fam == 'poisson') offs else NULL,
		Y_int = if (!is.null(alpha)) .lame_ltm_codes(Y) else NULL
	)
	codes = if (!is.null(alpha)) as.vector(.lame_ltm_codes(Y)) else NULL
	draw = function(m) {
		scale_k = kernel == 'scale'
		s2v = if (scale_k) 1 / rgamma(m, 0.5, rate = s20 / 2) else rep(s2, m)
		rhov = if (scale_k && mode == 'dir') 2 * rbeta(
			m,
			.lame_rho_prior_power() + 1,
			.lame_rho_prior_power() + 1
		) - 1 else rep(rho, m)
		B = matrix(rnorm(m * p), m, p) * sqrt(if (scale_k && mode == 'bip') V0 * s2v else V0)
		U1 = matrix(rnorm(m * nr), m, nr); U2 = matrix(rnorm(m * nc), m, nc)
		Um = Vm = NULL
		if (kernel == 'uv') {
			if (mode == 'dir') { Z2 = matrix(rnorm(m * nr * 2), m * nr, 2) %*% chol(Suv_t); Um = matrix(Z2[, 1], m, nr); Vm = matrix(Z2[, 2], m, nr) }
			else if (mode == 'sym') { Um = matrix(rnorm(m * nr, 0, 1 / sqrt(tau_t)), m, nr); Vm = Um * Ld_t }
			else { Um = matrix(rnorm(m * nr), m, nr); Vm = matrix(rnorm(m * nc), m, nc) }
		}
		if (mode == 'dir') { la = sqrt(Sab[1, 1]) * U1; lb = Sab[1, 2] / sqrt(Sab[1, 1]) * U1 + sqrt(
			Sab[2, 2] - Sab[1, 2]^2 / Sab[1, 1]
		) * U2 }
		else if (mode == 'sym') { la = sqrt(Sab[1, 1]) * U1; lb = la }
		else { la = sqrt(Sab[1, 1]) * U1; lb = sqrt(Sab[2, 2]) * U2 }
		if (!free) { la = la - rowMeans(la); lb = if (mode != 'sym') lb - rowMeans(lb) else la }
		ez = matrix(0, m, ncell)
		for (cc in seq_len(ncell)) ez[, cc] = (if (comp) B[, 1] else 0) + B[, p] * X[ci[cc, 1], ci[cc, 2], ci[cc, 3]] + la[, ci[cc, 1]] + lb[, ci[cc, 2]] +
			(if (kernel == 'uv') (if (mode == 'bip') G_t else 1) * Um[, ci[cc, 1]] * Vm[, ci[cc, 2]] else 0)
		sd_z = sqrt(s2v)
		E = matrix(rnorm(m * ncell), m, ncell) * sd_z
		if (mode == 'dir') {
			c1 = sqrt((1 + rhov) / 2); c2 = sqrt((1 - rhov) / 2)
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
		if (sq) { dg = which(ci[, 1] == ci[, 2]); E[, dg] = matrix(
			rnorm(m * length(dg)),
			m,
			length(dg)
		) * sqrt((1 + rhov) * s2v) }
		Zm = ez + E
		if (fam == 'poisson') {
			lw = rep(0, m)
			for (cc in seq_len(ncell)) {
				i = ci[cc, 1]; j = ci[cc, 2]; t = ci[cc, 3]
				if ((sq && i == j) || (mode == 'sym' && i > j)) next
				y = if (mode == 'sym') { yu = Y[i, j, t]; if (is.finite(yu)) yu else Y[j, i, t] } else Y[i, j, t]
				if (is.finite(y)) lw = lw + y * (Zm[, cc] + offs[t]) - exp(Zm[, cc] + offs[t])
			}
			w = exp(lw - max(lw))
		} else {
			ok = rep(TRUE, m)
			if (is.null(alpha)) {
				cells = which(is.finite(as.vector(Y)))
				for (t in 1:N) { ct = cells[ci[cells, 3] == t]; for (a_ in ct) for (b_ in ct) if (as.vector(Y)[a_] < as.vector(Y)[b_]) ok = ok & (Zm[, a_] < Zm[, b_]) }
			} else {
				af = c(-Inf, alpha, Inf)
				for (cc in which(is.finite(as.vector(Y)))) ok = ok & Zm[, cc] > af[codes[cc]] & Zm[, cc] <= af[codes[cc] + 1L]
			}
			w = as.numeric(ok)
		}
		list(B = B, la = la, lb = lb, ez = ez, Zm = Zm, w = w, s2 = s2v, rho = rhov, Um = Um, Vm = Vm)
	}
	state = function(bt, r) list(
		Z = array(bt$Zm[r, ], c(nr, nc, N)),
		EZ = array(bt$ez[r, ], c(nr, nc, N)),
		la = bt$la[r, ],
		lb = bt$lb[r, ],
		beta = bt$B[r, ],
		s2 = bt$s2[r],
		rho = bt$rho[r],
		U = if (!is.null(bt$Um)) matrix(bt$Um[r, ], ncol = 1L),
		V = if (!is.null(bt$Vm)) matrix(bt$Vm[r, ], ncol = 1L)
	)
	stat = function(s) {
		Z = s$Z; EZ = s$EZ; a1 = s$la[1]; b1 = s$lb[1]; a2 = s$la[2]; b2 = s$lb[min(2, length(s$lb))]
		r1 = Z[1, 2, 1] - EZ[1, 2, 1]; r2 = Z[2, 1, 1] - EZ[2, 1, 1]
		c(
			s$beta[1],
			s$beta[p],
			a1,
			b1,
			a2,
			b2,
			a1 * b1,
			a1 * a2,
			a1^2,
			b1^2,
			sum(s$la^2),
			sum(s$lb^2),
			if (mode == 'bip') 0 else sum(s$la * s$lb),
			Z[1, 2, 1],
			Z[1, 2, 1]^2,
			Z[2, 3, 1],
			Z[1, 2, 1] * Z[2, 1, 1],
			Z[1, 1, 1],
			Z[1, 1, 1]^2,
			mean(Z),
			mean(EZ[1, , ]),
			mean(EZ[, 1, ]),
			EZ[1, 2, 1] * EZ[2, 1, 1],
			r1 * r2,
			r1^2,
			s$beta[p] * a1,
			if (comp) c(s$beta[1] * a1, s$beta[1]^2) else c(0, 0),
			if (kernel == 'uv') c(
				sum(drop(s$U)^2),
				sum(drop(s$V)^2),
				drop(s$U)[1] * drop(s$V)[2],
				drop(s$U)[1],
				drop(s$V)[1],
				sum(outer(drop(s$U), drop(s$V))^2),
				drop(s$U)[2] * drop(s$U)[1]
			) else NULL,
			if (kernel == 'scale') c(
				log(s$s2),
				s$s2,
				s$rho,
				s$rho^2,
				log(s$s2) * s$rho,
				1 / s$s2,
				(Z[1, 2, 1] - EZ[1, 2, 1])^2 / s$s2
			) else NULL
		)
	}
	lev = c(TRUE, mode != 'sym')
	kern = function(s) {
		if (kernel == 'location') {
			o = .lame_ltm_location(
				cx,
				s$Z,
				s$EZ,
				s$la,
				s$lb,
				NULL,
				NULL,
				s$beta,
				Sab,
				iV0,
				centred = comp,
				lev = lev,
				alpha = alpha
			)
			s[c('Z', 'EZ', 'la', 'lb', 'beta')] = o[c('Z', 'EZ', 'la', 'lb', 'beta')]
		} else if (kernel == 'uv') {
			o = .lame_ltm_uv_visit(
				cx,
				s$Z,
				s$EZ,
				s$U,
				if (mode == 'sym') NULL else s$V,
				if (mode == 'bip') matrix(G_t, 1, 1) else NULL,
				if (mode == 'sym') Ld_t else NULL,
				Suv = if (mode == 'dir') Suv_t else NULL,
				tau = if (mode == 'sym') tau_t else NULL
			)
			s[c('Z', 'EZ', 'U')] = o[c('Z', 'EZ', 'U')]
			s$V = if (mode == 'sym') o$U * Ld_t else o$V
		} else if (kernel == 'coef') {
			o = .lame_ltm_coef(cx, s$Z, s$EZ, s$beta, if (comp) list(array(1, c(nr, nc, N)), X) else list(X), iV0)
			s[c('Z', 'EZ', 'beta')] = o[c('Z', 'EZ', 'beta')]
		} else if (kernel == 'scale') {
			df_b = if (mode == 'bip') p else 0; ss_b = if (mode == 'bip') sum(s$beta^2) / V0 else 0
			o = .lame_ltm_scale(
				cx,
				s$Z,
				s$EZ,
				s$s2,
				s$rho,
				c(shape0 = (1 + df_b) / 2, rate0 = (s20 + ss_b) / 2),
				two = mode == 'dir',
				rho_logprior = .lame_rho_prior_log,
				tries = 4L
			)
			s[c('Z', 's2', 'rho')] = o[c('Z', 's2', 'rho')]
			if (sq) for (t in 1:N) { zt = s$Z[, , t]; diag(zt) = rnorm(
				nr,
				diag(s$EZ[, , t]),
				sqrt((1 + s$rho) * s$s2)
			); s$Z[, , t] = zt }
		} else {
			o = .lame_ltm_visit(
				cx,
				s$Z,
				s$EZ,
				s$la,
				s$lb,
				NULL,
				NULL,
				s$beta,
				Sab,
				iV0,
				centred = comp,
				lev = lev,
				alpha = alpha
			)
			s[c('Z', 'EZ', 'la', 'lb', 'beta')] = o[c('Z', 'EZ', 'la', 'lb', 'beta')]
		}
		if (mode == 'sym') s$lb = s$la
		s
	}
	if (!is.null(hook)) hook()
	Dl = list(); Wl = list(); done = 0
	while (done < ndraw) {
		bt = draw(nbatch); done = done + nbatch
		keep = which(bt$w > 1e-12 * max(bt$w))
		if (fam == 'poisson' && length(keep) > 6000L) {
			w0 = quantile(bt$w[keep], 1 - 2 * 6000 / length(keep))
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
