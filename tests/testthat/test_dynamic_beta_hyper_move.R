# hyperparameters of the dynamic_beta paths (R/dynamic_beta_hyper_move.R): exact conditionals of the innovation scale and
# the persistence every sweep, and an interweaving step on (scale, persistence, path) with the standardised innovations fixed
# checked here: the compiled sufficient statistics, the maps, the exactness of the step, and the switches
skip_on_cran()

# the sufficient statistics

test_that('the compiled sufficient statistics are those of the R observation model', {
	set.seed(3)
	cmp = function(n_r, n_c, bip, sym, rho, nas) {
		nt = 3L
		Xd = Xs = Z = off = vector('list', nt)
		for (t in seq_len(nt)) {
			Xd[[t]] = matrix(rnorm(n_r * n_c * 2), n_r * n_c, 2)
			Xs[[t]] = matrix(rnorm(n_r * n_c), n_r * n_c, 1)
			Zt = matrix(rnorm(n_r * n_c), n_r, n_c)
			if (nas > 0) Zt[sample(length(Zt), round(nas * length(Zt)))] = NA
			Z[[t]] = Zt
			off[[t]] = matrix(rnorm(n_r * n_c, 0, 0.3), n_r, n_c)
		}
		bs = rnorm(1)
		s2 = 1.7
		use_rho = (!bip) && (!sym) && rho != 0
		cc = lame:::beta_dynamic_stats_cpp(Xd, Xs, Z, off, bs, s2, rho, bip, sym, use_rho)
		for (t in seq_len(nt)) {
			r = Z[[t]] - off[[t]] - matrix(Xs[[t]] %*% bs, n_r, n_c)
			zz = lame:::.lame_beta_observation_stats(Xd[[t]], r, s2, if (use_rho) rho else 0, bip, sym)
			expect_equal(cc$Omega[, , t], zz$Q, tolerance = 1e-10)
			expect_equal(unname(cc$h[t, ]), unname(zz$h), tolerance = 1e-10)
		}
	}
	cmp(10, 10, FALSE, FALSE, 0, 0)
	cmp(10, 10, FALSE, FALSE, 0.4, 0.05)
	cmp(10, 10, FALSE, FALSE, -0.3, 0.2)
	cmp(9, 6, TRUE, FALSE, 0, 0.1)
	cmp(10, 10, FALSE, TRUE, 0.8, 0)
})

test_that('the filter reads the statistics it is handed and draws the same path', {
	set.seed(2)
	n = 12; nt = 4L; p = 3L
	Xd = Xs = Z = off = vector('list', nt)
	for (t in seq_len(nt)) {
		Xd[[t]] = cbind(1, matrix(rnorm(n * n * 2), n * n, 2))
		Xs[[t]] = matrix(rnorm(n * n), n * n, 1)
		Z[[t]] = matrix(rnorm(n * n), n, n)
		off[[t]] = matrix(rnorm(n * n, 0, 0.3), n, n)
	}
	args = list(
		Xdyn_list = Xd,
		Xstat_list = Xs,
		Z_list = Z,
		offset_list = off,
		beta_static = 0.3,
		rho_by_coef = rep(0.8, p),
		sigma_by_coef = rep(0.2, p),
		Lambda = diag(1, p),
		beta0_mean = rep(0, p),
		beta0_cov = diag(10, p),
		s2 = 1.3,
		dyad_rho = 0.3,
		bipartite = FALSE,
		symmetric = FALSE,
		use_dyad_rho = TRUE
	)
	for (use_rho in c(TRUE, FALSE)) {
		args$use_dyad_rho = use_rho
		st = do.call(
			lame:::beta_dynamic_stats_cpp,
			args[c(
				'Xdyn_list', 'Xstat_list', 'Z_list', 'offset_list', 'beta_static', 's2', 'dyad_rho', 'bipartite', 'symmetric', 'use_dyad_rho'
			)]
		)
		set.seed(9)
		a = do.call(lame:::sample_beta_dynamic_cpp, args)
		set.seed(9)
		b = do.call(lame:::sample_beta_dynamic_cpp, c(args, list(Omega_in = st$Omega, h_in = st$h)))
		expect_identical(a$path, b$path)
	}
	# statistics of the wrong shape are refused
	expect_error(do.call(lame:::sample_beta_dynamic_cpp, c(args, list(Omega_in = st$Omega[, , -1], h_in = st$h))))
})

# the maps

test_that('the path map keeps the standardised innovations and inverts itself', {
	set.seed(4)
	lam = c(0.7, 1.9)
	for (gaps in list(rep(1, 7), c(1, 1, 2, 1, 3, 1, 1))) {
		nt = length(gaps) + 1L
		x = matrix(rnorm(nt * 2), nt, 2)
		for (kind in c('ar1', 'rw1', 'rw2', 'matern32')) {
			rho0 = if (kind == 'ar1') 0.7 else 1
			same = lame:::.lame_beta_nc_path(x, kind, 0.3, rho0, 0.3, rho0, lam, gaps)
			expect_equal(same, x, tolerance = 1e-12)
			rho1 = if (kind == 'ar1') 0.55 else 1
			moved = lame:::.lame_beta_nc_path(x, kind, 0.3, rho0, 0.5, rho1, lam, gaps)
			expect_gt(max(abs(moved - x)), 1e-3)
			back = lame:::.lame_beta_nc_path(moved, kind, 0.5, rho1, 0.3, rho0, lam, gaps)
			expect_equal(back, x, tolerance = 1e-10)
			if (kind %in% c('ar1', 'rw1')) {
				# the first state does not move; the innovations in units of the scale do not change
				expect_identical(moved[1, ], x[1, ])
				z0 = lame:::.lame_beta_gap_transition(rho0, gaps)
				z1 = lame:::.lame_beta_gap_transition(rho1, gaps)
				w0 = (x[-1, ] - z0$f * x[-nt, ]) / (0.3 * sqrt(outer(z0$q, lam)))
				w1 = (moved[-1, ] - z1$f * moved[-nt, ]) / (0.5 * sqrt(outer(z1$q, lam)))
				expect_equal(w1, w0, tolerance = 1e-10)
			}
			if (kind == 'matern32') expect_equal(moved, x * (0.5 / 0.3), tolerance = 1e-12)
			if (kind == 'rw2') {
				# the part the second differences leave flat (constants and linear functions of time) stays where it is
				flat = lame:::.lame_beta_rw2_flat(x, c(0, cumsum(gaps)))
				expect_equal(lame:::.lame_beta_rw2_flat(moved, c(0, cumsum(gaps))), flat, tolerance = 1e-10)
				expect_equal(moved - flat, (x - flat) * (0.5 / 0.3), tolerance = 1e-10)
			}
		}
	}
})

# exactness

# one step from the joint law of (hyperparameters, path, data): sigma_g^2 ~ IG(shape, scale), rho_g truncated normal
# (logistic normal and a shared Gamma scale when pooled), the path from the kind's prior, y = path + noise of precision Omega_t
# returns the paired differences of moments before and after the step
hm_geweke = function(kind, pool, gaps_on, ndraw, seed, step = lame:::.lame_beta_nc_step, nt = 8L) {
	set.seed(seed)
	ns = asNamespace('lame')
	p = 3L; group_id = c(1L, 1L, 2L); G = 2L; lam = c(0.8, 1.5, 1)
	shape = 2; scale = 0.15; a_hyp = 2; b_hyp = 0.3; rho_mean = 0.8; rho_sd = 0.15; lo = 0; hi = 0.999; rho_mu = 0.5; rho_tau = 0.8
	m0 = 0; v0 = 4
	pos = if (gaps_on) c(1, 2, 3, 5, 6, 9, 10, 11)[seq_len(nt)] else seq_len(nt)
	gaps = diff(pos)
	Om = array(0, c(p, p, nt))
	for (t in seq_len(nt)) { A = matrix(rnorm(p * p), p); Om[, , t] = 5 * (diag(p) * 0.6 + 0.4 * tcrossprod(A) / p + 0.1) }
	rtn = function(mu, sd) repeat { x = rnorm(1, mu, sd); if (x > lo && x < hi) return(x) }
	stat_of = function(s) c(
		s$sig,
		log(s$sig),
		s$sig^2,
		if (kind == 'ar1') c(s$rho, s$rho^2, s$rho * s$sig) else numeric(0),
		as.vector(s$path[c(1, nt %/% 2, nt), ]),
		as.vector(s$path[c(1, nt %/% 2, nt), ])^2,
		sum(s$path[, 1]^2),
		sum(diff(s$path[, 3])^2),
		s$sig[1] * s$path[nt, 1],
		sum((s$path - s$y)^2),
		colSums((s$path - s$y)^2),
		sum(s$path * s$y),
		log(s$sig) * sum((s$path[, 1] - s$y[, 1])^2)
	)
	D = vector('list', ndraw)
	acc = att = 0
	for (k in seq_len(ndraw)) {
		S = if (pool %in% c('sigma', 'both')) rgamma(1, a_hyp, b_hyp) else scale
		sig = sqrt(1 / rgamma(G, shape, S))
		rho = if (pool %in% c('rho', 'both') && kind == 'ar1') lo + (hi - lo) * plogis(rnorm(G, rho_mu, rho_tau)) else if (kind == 'ar1') replicate(G, rtn(rho_mean, rho_sd)) else rep(1, G)
		path = matrix(0, nt, p)
		for (j in seq_len(p)) {
			g = group_id[j]
			if (kind %in% c('ar1', 'rw1')) {
				z = ns$.lame_beta_gap_transition(rho[g], gaps)
				path[1, j] = rnorm(1, m0, sqrt(v0))
				for (t in 2:nt) path[t, j] = z$f[t - 1] * path[t - 1, j] + sig[g] * sqrt(lam[j] * z$q[t - 1]) * rnorm(1)
			} else {
				K = ns$.dynamic_beta_kind_to_prec(kind, nt, NULL, pos)
				path[, j] = backsolve(chol(K / (sig[g]^2 * lam[j])), rnorm(nt))
			}
		}
		y = path
		for (t in seq_len(nt)) y[t, ] = path[t, ] + drop(t(chol(solve(Om[, , t]))) %*% rnorm(p))
		h = t(vapply(seq_len(nt), function(t) drop(Om[, , t] %*% y[t, ]), numeric(p)))
		s = list(sig = sig, rho = rho, path = path, y = y)
		pscale = if (pool %in% c('sigma', 'both')) ns$.pool_update_sigma_scale(
			sig,
			alpha = shape,
			hyper_a = a_hyp,
			hyper_b = b_hyp
		) else scale
		out = step(
			path,
			list(Omega = Om, h = h),
			group_id,
			G,
			kind,
			lam,
			sig,
			rho,
			shape,
			pscale,
			rho_mean,
			rho_sd,
			lo,
			hi,
			rho_pool = if (pool %in% c('rho', 'both') && kind == 'ar1') list(mu = rho_mu, tau = rho_tau) else NULL,
			movable_rho = identical(kind, 'ar1'),
			gaps = gaps,
			tries = 6L
		)
		acc = acc + out$accepts
		att = att + out$attempts
		s1 = s; s1$path = out$path; s1$sig = out$sigma; s1$rho = out$rho
		D[[k]] = stat_of(s1) - stat_of(s)
	}
	D = do.call(rbind, D)
	N = nrow(D)
	mu = colMeans(D)
	sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), T2 = T2, df = sum(kp), p = pchisq(T2, sum(kp), lower.tail = FALSE), accept = acc / att)
}

test_that('one interweaving step leaves the joint law of hyperparameters, path and data unchanged', {
	skip_on_cran()
	for (cfg in list(
		list('ar1', 'none', FALSE, 11),
		list('ar1', 'both', TRUE, 12),
		list('rw1', 'sigma', FALSE, 13),
		list('matern32', 'none', TRUE, 14)
	)) {
		r = hm_geweke(cfg[[1]], cfg[[2]], cfg[[3]], 4000, cfg[[4]])
		expect_gt(r$accept, 0.05)
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
	# rw2: the part of the path its second differences leave flat stays fixed, its (ridge) prior enters the ratio
	r = hm_geweke('rw2', 'none', FALSE, 4000, 15)
	expect_gt(r$accept, 0.05)
	expect_lt(r$maxz, 4.5)
})

test_that('the exactness check fails for a step that accepts every proposal', {
	skip_on_cran()
	src = deparse(lame:::.lame_beta_nc_step)
	src = sub('log(stats::runif(1L)) < ', 'TRUE || ', src, fixed = TRUE)
	broken = eval(parse(text = paste(src, collapse = '\n')), envir = asNamespace('lame'))
	r = hm_geweke('ar1', 'none', FALSE, 3000, 21, step = broken)
	expect_gt(r$maxz, 6)
})

# the switches and the fits

hm_fit = function(nscan = 80, ...) {
	data(YX_bin_list, envir = environment())
	suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = 'binary',
		R = 0,
		dynamic_beta = 'dyad',
		nscan = nscan,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5,
		...
	)))
}

test_that('beta_hyper_nc is validated', {
	expect_error(hm_fit(prior = list(beta_hyper_nc = NA)), 'beta_hyper_nc must be TRUE or FALSE', fixed = TRUE)
})

test_that('the hyperparameters are drawn every sweep', {
	skip_on_cran()
	on = hm_fit()
	sg_on = as.numeric(as.matrix(on$SIGMA_BETA)[, 1])
	# every sweep: a new value at (nearly) every stored sweep
	expect_gt(length(unique(sg_on)), 0.9 * length(sg_on))
	# the interweaving step is recorded when it ran and only then
	expect_true(!is.null(on$beta_hyper_move))
	expect_equal(unname(on$beta_hyper_move[['sweeps']]), 100)
	expect_gt(on$beta_hyper_move[['accepts']], 0)
	expect_lt(on$beta_hyper_move[['rate']], 1)
	expect_null(hm_fit(prior = list(beta_hyper_nc = FALSE))$beta_hyper_move)
	# without the step the schedule is still every sweep
	sg_h = as.numeric(as.matrix(hm_fit(prior = list(beta_hyper_nc = FALSE))$SIGMA_BETA)[, 1])
	expect_gt(length(unique(sg_h)), 0.9 * length(sg_h))
})

test_that('every kind, pooling and gapped periods fit with the step', {
	skip_on_cran()
	data(YX_bin_list, envir = environment())
	for (kind in c('ar1', 'rw1', 'rw2', 'matern32')) for (pool in c('none', 'both')) {
		fit = suppressWarnings(suppressMessages(lame(
			YX_bin_list$Y,
			YX_bin_list$X,
			family = 'binary',
			R = 0,
			dynamic_beta = c('intercept', 'dyad'),
			dynamic_beta_kind = kind,
			dynamic_beta_pool = pool,
			nscan = 40,
			burn = 15,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			seed = 5
		)))
		expect_true(all(is.finite(fit$BETA)))
		expect_true(all(is.finite(as.matrix(fit$SIGMA_BETA))))
		expect_gt(fit$beta_hyper_move[['attempts']], 0)
	}
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = 'binary',
		R = 0,
		dynamic_beta = 'dyad',
		time_index = c(1, 2, 4, 8),
		nscan = 40,
		burn = 15,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5
	)))
	expect_true(all(is.finite(fit$BETA)))
	expect_gt(fit$beta_hyper_move[['attempts']], 0)
})
