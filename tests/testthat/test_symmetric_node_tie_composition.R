# the symmetric node tie composed with the scale moves of the binary and ordinal families: one step of the move on a tied state
# leaves the joint law of the tied coefficients, the effects, the latent array and the ties unchanged
skip_on_cran()

# tie_geweke: states are drawn from the prior of the tied model (the dyadic coefficient, the pair of a nodal covariate as the mean of a draw
# of the pair, static additive effects, an unscaled part m standing for dynamic effects), the latent array around the predictor, and the ties
# (binary) from its signs; `visit` is one step of the move; the scaled parameters are multiplied by the factor it returns.
# dcount: the dimension of the scaled parameters the move is told (the tied pair counts once);
# untie: a broken composition that scales only the row slice of the pair
tie_geweke = function(
	kind,
	ndraw,
	seed,
	withm = TRUE,
	dcount = NULL,
	untie = FALSE,
	visit = NULL,
	n = 6L,
	nt = 2L,
	miss = 0.1
) {
	ns = asNamespace('lame')
	set.seed(seed)
	nr = n; nc = n; p = 3L; vb = 0.8; va = 0.5
	X = array(rnorm(nr * nc * p * nt), c(nr, nc, p, nt))
	for (t in 1:nt) {
		d = X[, , 1, t]; d[lower.tri(d)] = t(d)[lower.tri(d)]; X[, , 1, t] = d
		xr = rnorm(nr)
		X[, , 2, t] = matrix(xr, nr, nc)
		X[, , 3, t] = matrix(xr, nr, nc, byrow = TRUE)
	}
	mfix = if (withm) array(rnorm(nr * nc * nt, 0.5, 0.7), c(nr, nc, nt)) else NULL
	if (!is.null(mfix)) for (t in 1:nt) { m_ = mfix[, , t]; m_[lower.tri(m_)] = t(m_)[lower.tri(m_)]; mfix[, , t] = m_ }
	upm = upper.tri(matrix(0, nr, nc))
	E_of = function(s) {
		E = array(0, c(nr, nc, nt))
		for (t in 1:nt) E[, , t] = outer(s$a, s$a, '+') + matrix(
			matrix(X[, , , t], nr * nc, p) %*% s$beta,
			nr,
			nc
		) + if (is.null(mfix)) 0 else mfix[, , t]
		E
	}
	draw_state = function() {
		pr = rnorm(2, 0, sqrt(vb)); mm = mean(pr)
		beta = c(rnorm(1, 0, sqrt(vb)), mm, mm)
		a = rnorm(nr, 0, sqrt(va))
		E = E_of(list(a = a, beta = beta))
		Z = E
		for (t in 1:nt) { z = matrix(rnorm(nr * nc), nr, nc); z[lower.tri(z)] = t(z)[lower.tri(z)]; Z[, , t] = Z[, , t] + z }
		Y = NULL
		if (kind == 'binary') {
			Y = 1 * (Z > 0)
			ms = array(runif(nr * nc * nt) < miss, c(nr, nc, nt))
			for (t in 1:nt) ms[, , t][lower.tri(ms[, , t])] = t(ms[, , t])[lower.tri(ms[, , t])]
			Y[ms] = NA
			for (t in 1:nt) diag(Y[, , t]) = NA
		}
		list(beta = beta, a = a, E = E, Z = Z, Y = Y)
	}
	stat_of = function(s) {
		z = s$Z; eb = s$E
		msk = array(rep(upm, nt), dim(z)); zz = z * msk; ee = eb * msk
		y1 = if (is.null(s$Y)) 0 * z else ifelse(is.na(s$Y), 0, s$Y) * msk
		c(
			s$beta,
			s$beta^2,
			s$beta[1] * s$beta[2],
			sum(s$a^2),
			sum(s$a),
			(s$beta[2] - s$beta[3])^2,
			abs(s$beta[2] - s$beta[3]),
			sapply(
				1:nt,
				function(t) c(
					mean(zz[, , t]),
					mean(zz[, , t]^2),
					mean((zz[, , t] - ee[, , t])^2),
					mean(zz[, , t] * ee[, , t]),
					mean(y1[, , t]),
					mean(y1[, , t] * ee[, , t])
				)
			),
			sum(y1 * ee) * s$beta[1],
			mean(zz^2) * s$beta[2]
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		s = draw_state()
		# the priors of the scaled parameters at the tied state: the pair enters through the exchangeable covariance, 2 m^2 / vb
		Apri = sum(s$beta^2) / vb + sum(s$a^2) / va
		dth = (if (is.null(dcount)) 2L else dcount) + nr
		s1 = s
		if (kind == 'binary') {
			cx = ns$.lame_bso_context(s$Y, 'sym')
			o = (if (is.null(visit)) ns$.lame_bso_visit else visit)(s$Z, s$E, mfix, cx, 0, Apri, 0, dth)
		} else {
			o = (if (is.null(visit)) ns$.lame_osm_visit else visit)(s$Z, s$E, mfix, 'sym', rep(0, nt), 1, Apri, 0, dth, affine = TRUE, orbit = TRUE)
		}
		s1$Z = o$Z
		s1$beta = o$c * s$beta
		if (untie) s1$beta[3] = s$beta[3]
		s1$a = o$c * s$a
		s1$E = E_of(s1)
		D[[it]] = stat_of(s1) - stat_of(s)
	}
	D = do.call(rbind, D)
	N = nrow(D); mu = colMeans(D); sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), T2 = T2, df = sum(kp), p = pchisq(T2, sum(kp), lower.tail = FALSE))
}

test_that('the scale moves of the binary and ordinal families keep a tied pair tied and leave the joint law unchanged when the pair counts as one dimension', {
	skip_on_cran()
	for (cfg in list(
		list('binary', 81, TRUE),
		list('binary', 82, FALSE),
		list('ordinal', 83, TRUE),
		list('ordinal', 84, FALSE)
	)) {
		r = tie_geweke(cfg[[1]], 7000, cfg[[2]], withm = cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('counting the pair as two dimensions, or leaving it untied after the move, fails the check', {
	skip_on_cran()
	for (kd in c('binary', 'ordinal')) {
		expect_gt(tie_geweke(kd, 12000, 85, dcount = 3L)$maxz, 4.5)
		expect_gt(tie_geweke(kd, 7000, 86, untie = TRUE)$maxz, 5)
	}
})

test_that('a fit tells the scale moves that a tied pair is one dimension', {
	stn_p = function(fam, n = 12L, Tn = 3L) {
		set.seed(7)
		rn = sprintf('a%02d', seq_len(n))
		sym = function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m }
		Y = lapply(
			seq_len(Tn),
			function(t) { z = sym(matrix(rnorm(n * n), n)) + 0.3; y = if (fam == 'binary') 1 * (z > 0.5) else matrix(
				findInterval(z, c(0, .7, 1.4)),
				n
			); y = sym(y); diag(y) = NA; dimnames(y) = list(rn, rn); y }
		)
		X = lapply(
			seq_len(Tn),
			function(t) array(sym(matrix(rnorm(n * n), n)), c(n, n, 1), dimnames = list(rn, rn, 'xd'))
		)
		Xr = lapply(seq_len(Tn), function(t) matrix(rnorm(n), n, 1, dimnames = list(rn, 'xn')))
		names(Y) = names(X) = names(Xr) = paste0('p', seq_len(Tn))
		list(Y = Y, X = X, Xr = Xr)
	}
	got = new.env()
	got$d = numeric(0)
	local_mocked_bindings(
		.lame_bso_visit = function(
			Z,
			E,
			m,
			cx,
			rho,
			A_pri,
			B_pri,
			d_th,
			...
		) { got$d = c(got$d, d_th); list(Z = Z, c = 1, evals = 0L) },
		.lame_osm_visit = function(
			Z,
			E,
			m,
			mode,
			rho,
			s2,
			A_pri,
			B_pri,
			d_th,
			...
		) { got$d = c(got$d, d_th); list(Z = Z, c = 1) },
		.package = 'lame'
	)
	first_d = function(d, fam, tie, ...) {
		got$d = numeric(0)
		suppressWarnings(suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			Xrow = d$Xr,
			symmetric = TRUE,
			family = fam,
			R = 0L,
			burn = 3,
			nscan = 6,
			odens = 1,
			seed = 2,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(symmetric_node_tie = tie, binary_scale_every = 1),
			...
		)))
		got$d[1L]
	}
	for (fam in c('binary', 'ordinal')) {
		d = stn_p(fam)
		expect_equal(first_d(d, fam, FALSE) - first_d(d, fam, TRUE), 1)
		expect_equal(first_d(d, fam, FALSE, dynamic_beta = TRUE) - first_d(d, fam, TRUE, dynamic_beta = TRUE), 3)
	}
})
