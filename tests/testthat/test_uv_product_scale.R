# the size of the product of dynamic positions with the latent array integrated out or carried (R/uv_product_scale.R)
skip_on_cran()

ups_mut = function(fn, from, to) {
	f = get(fn, envir = asNamespace('lame'))
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# one joint draw of the parameters, the latent array and the data, so the state is a posterior draw
# `stationary`: the inverse gamma (shape, scale) is on the stationary variance w, sigma^2 = w (1 - rho_uv^2)
ups_state = function(family, mode, kind = 'ar1', n = 6L, nB = 5L, R = 2L, N = 3L, rho_uv = 0.7, shape = 3, scale = 2,
	g_var = 1.5, s2 = 0.3, rho = 0.4, nu = 4, kappa = 1.2, pi_snap = 0.3, stationary = FALSE) {
	bip = mode == 'bip'
	nr = n
	nc = if (bip) nB else n
	sigma = sqrt((if (stationary) 1 - rho_uv^2 else 1) / rgamma(1, shape, rate = scale))
	path = function(m) {
		lam = matrix(1, m, N)
		del = matrix(0, m, N)
		# kind t: the first period and every transition with its own weight
		if (kind == 't') lam[] = rgamma(m * N, nu / 2, nu / 2)
		if (kind == 'snap') del[, -1] = rbinom(m * (N - 1L), 1, pi_snap)
		X = array(0, c(m, R, N))
		X[, , 1] = rnorm(m * R, 0, sigma / sqrt(1 - rho_uv^2)) / sqrt(lam[, 1])
		for (t in 2:N) {
			drift = rho_uv * X[, , t - 1] + matrix(rnorm(m * R), m) * sigma / sqrt(lam[, t])
			reset = rnorm(m * R, 0, if (is.null(kappa)) sigma / sqrt(1 - rho_uv^2) else kappa)
			X[, , t] = ifelse(matrix(del[, t], m, R) == 1, reset, drift)
		}
		list(X = X, lam = lam, del = del)
	}
	pu = path(nr)
	pv = path(nc)
	G = if (bip) matrix(rnorm(R * R, 0, sqrt(g_var)), R) else NULL
	base = array(-0.3 + rnorm(nr * nc * N, 0, 0.3), c(nr, nc, N))
	EZ = base + lame:::.lame_uvps_product(pu$X, pv$X, G)
	sd_e = if (family == 'normal') sqrt(s2) else 1
	if (family == 'poisson') {
		Z = EZ + rnorm(length(EZ), 0, sqrt(s2))
		Y = array(rpois(length(Z), exp(Z)), dim(Z))
	} else if (bip) {
		Z = EZ + rnorm(length(EZ), 0, sd_e)
		Y = if (family == 'normal') Z else 1 * (Z > 0)
	} else {
		Z = EZ
		for (t in 1:N) {
			e1 = matrix(rnorm(n * n), n); e2 = matrix(rnorm(n * n), n)
			e = e1
			lo = lower.tri(e)
			e[lo] = rho * t(e1)[lo] + sqrt(1 - rho^2) * e2[lo]
			diag(e) = rnorm(n, 0, sqrt(1 + rho))
			Z[, , t] = EZ[, , t] + sd_e * e
		}
		Y = if (family == 'normal') Z else 1 * (Z > 0)
	}
	if (!bip) for (t in 1:N) diag(Y[, , t]) = NA
	Y[sample(length(Y), 4L)] = NA
	list(Y = Y, Z = Z, EZ = EZ, base = base, U = pu$X, V = pv$X, G = G, sigma = sigma, lam_u = pu$lam, lam_v = pv$lam,
		del_u = pu$del, del_v = pv$del)
}

ups_law = function(ndraw, seed, family, mode, kind = 'ar1', draw = NULL, rho = 0.4, kappa = 1.2, max_abs = Inf) {
	if (is.null(draw)) draw = lame:::.lame_uvps_visit
	set.seed(seed)
	bip = mode == 'bip'
	rho_d = if (bip || family == 'poisson') 0 else rho
	stat = function(st) {
		P = lame:::.lame_uvps_product(st$U, st$V, st$G)
		split = log(apply(st$U^2, 3, sum) / apply(st$V^2, 3, sum))
		x = c(log(apply(P^2, 3, sum)), split, mean(st$Z), if (!bip) log(st$sigma), if (bip) log(sum(st$G^2)), if (kind == 't') c(log(sum(st$lam_u[, -1])), log(sum(st$lam_u[, 1]))))
		c(x, x^2)
	}
	out = NULL
	for (r in seq_len(ndraw)) {
		repeat {
			st = ups_state(family, mode, kind, rho = rho, kappa = kappa)
			if (max(abs(st$U), abs(st$V)) <= max_abs) break
		}
		cx = lame:::.lame_uvps_context(st$Y, family, mode)
		o = st
		for (k in 1:2) {
			v = draw(
				cx, o$Z, o$base + lame:::.lame_uvps_product(o$U, o$V, o$G), o$Y, o$U, o$V, 0.7, o$sigma, 3, 2,
				rho = rho_d, G = o$G, g_var = 1.5, max_abs = max_abs, kind = kind,
				lambda_u = if (kind == 't') o$lam_u, lambda_v = if (kind == 't') o$lam_v,
				kappa = if (kind == 'snap') kappa, delta_u = o$del_u, delta_v = o$del_v
			)
			o$Z = v$Z; o$U = v$U; o$V = v$V; o$G = v$G; o$sigma = v$sigma
			if (kind == 't') { o$lam_u = v$lambda_u; o$lam_v = v$lambda_v }
		}
		if (max(abs(o$U), abs(o$V)) > max_abs) stop('the bound was left')
		if (family == 'binary') {
			ok = is.finite(o$Y)
			if (any((o$Z[ok] > 0) != (o$Y[ok] == 1))) stop('the latent array left the sign region')
		}
		s1 = stat(o) - stat(st)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(nrow(out)))))
}

test_that('the size moves keep the posterior of bipartite Poisson and binary fits', {
	expect_lt(ups_law(1500, 9201, 'poisson', 'bip'), 4.5)
	expect_lt(ups_law(1500, 9202, 'binary', 'bip'), 4.5)
	expect_lt(ups_law(1500, 9203, 'poisson', 'bip', max_abs = 3), 4.5)
	# the check fails with a wrong jacobian of G
	expect_gt(ups_law(1500, 9204, 'poisson', 'bip', draw = ups_mut('.lame_uvps_visit', '+ dG * u', '+ 3 * dG * u')), 5)
})

test_that('the size moves keep the posterior of bipartite snap fits', {
	expect_lt(ups_law(1500, 9205, 'poisson', 'bip', kind = 'snap'), 4.5)
	expect_lt(ups_law(1500, 9206, 'binary', 'bip', kind = 'snap', kappa = NULL), 4.5)
	# the check fails with the resets of the bipartite paths weighted as drifts
	expect_gt(ups_law(1500, 9207, 'poisson', 'bip', kind = 'snap', kappa = 0.4, draw = ups_mut(
		'.lame_uvps_visit',
		's2/.lame_snap_kappa(kappa, rho_uv, sigma_uv)^2',
		'1'
	)), 5)
})

test_that('the size moves keep the posterior of directed fits of every kind', {
	expect_lt(ups_law(1500, 9211, 'binary', 'dir'), 4.5)
	expect_lt(ups_law(1500, 9212, 'poisson', 'dir'), 4.5)
	expect_lt(ups_law(1500, 9213, 'binary', 'dir', kind = 't'), 4.5)
	expect_lt(ups_law(1500, 9214, 'poisson', 'dir', kind = 'snap'), 4.5)
	expect_lt(ups_law(1500, 9218, 'binary', 'dir', kind = 'snap', kappa = NULL), 4.5)
	# the checks fail without the jacobian or the path prior of one period, with a wrong split under kind t and
	# with the resets of kind snap weighted as drifts
	expect_gt(ups_law(1500, 9215, 'binary', 'dir', draw = ups_mut('.lame_uvps_visit', '+ dt * u/2', '')), 5)
	expect_gt(ups_law(1500, 9216, 'binary', 'dir', kind = 't', draw = ups_mut('.lame_uvps_visit', 'lp(exp(u/2))', '0')), 5)
	expect_gt(ups_law(1500, 9217, 'binary', 'dir', kind = 't', draw = ups_mut('.lame_uvps_visit', 'fv(exp(-v))', 'fv(exp(v))')), 5)
	expect_gt(ups_law(1500, 9219, 'poisson', 'dir', kind = 'snap', kappa = 0.4, draw = ups_mut(
		'.lame_uvps_visit',
		's2/.lame_snap_kappa(kappa, rho_uv, sigma_uv)^2',
		'1'
	)), 5)
})

# the persistence drawn with the paths following: rho from its prior, then the paths and the data
upnc_law = function(ndraw, seed, family, mode, kind = 'ar1', draw = NULL, kappa = 1.2, max_abs = Inf, sigma_max = Inf,
	stationary = FALSE, stationary_draw = stationary) {
	if (is.null(draw)) draw = lame:::.lame_uv_nc_rho_draw
	set.seed(seed)
	bip = mode == 'bip'
	rho_d = if (bip || family == 'poisson') 0 else 0.4
	stat = function(st, r) {
		P = lame:::.lame_uvps_product(st$U, st$V, st$G)
		x = c(r, log(st$sigma), log(apply(P^2, 3, sum)))
		c(x, x^2)
	}
	out = NULL
	for (i in seq_len(ndraw)) {
		repeat {
			r = rnorm(1, 0.7, 0.3)
			if (abs(r) >= 0.99) next
			st = ups_state(family, mode, kind, rho_uv = r, rho = rho_d, kappa = kappa, stationary = stationary)
			if (max(abs(st$U), abs(st$V)) <= max_abs && st$sigma <= sigma_max) break
		}
		cx = lame:::.lame_uvps_context(st$Y, family, mode)
		o = st
		ro = r
		for (k in 1:2) {
			v = draw(
				cx, o$Z, o$base + lame:::.lame_uvps_product(o$U, o$V, o$G), o$Y, o$U, o$V, ro, o$sigma, 0.7, 0.3, 3, 2,
				rho = rho_d, s2 = 0.3, G = o$G, sigma_max = sigma_max, max_abs = max_abs, kind = kind,
				lambda_u = if (kind == 't') o$lam_u, lambda_v = if (kind == 't') o$lam_v,
				kappa = if (kind == 'snap') kappa, delta_u = o$del_u, delta_v = o$del_v, stationary = stationary_draw
			)
			o$Z = v$Z; o$U = v$U; o$V = v$V; o$sigma = v$sigma; ro = v$rho
		}
		if (max(abs(o$U), abs(o$V)) > max_abs || o$sigma > sigma_max) stop('the bound was left')
		if (family == 'normal') {
			ok = is.finite(o$Y)
			if (any(o$Z[ok] != o$Y[ok])) stop('an observed cell moved')
		}
		s1 = stat(o, ro) - stat(st, r)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[i, ] = s1
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(nrow(out)))))
}

test_that('the persistence with the paths following keeps the posterior of each family, mode and kind', {
	expect_lt(upnc_law(1500, 9221, 'normal', 'dir'), 4.5)
	expect_lt(upnc_law(1500, 9222, 'binary', 'dir', kind = 't'), 4.5)
	expect_lt(upnc_law(1500, 9223, 'poisson', 'dir', kind = 'snap'), 4.5)
	expect_lt(upnc_law(1500, 9224, 'normal', 'dir', kind = 'snap', kappa = NULL), 4.5)
	expect_lt(upnc_law(1500, 9225, 'poisson', 'bip', max_abs = 3, sigma_max = 1.2), 4.5)
	expect_lt(upnc_law(1500, 9226, 'binary', 'bip'), 4.5)
	# the check fails with a wrong prior of sigma_uv along the ridge
	expect_gt(upnc_law(1500, 9227, 'normal', 'dir', draw = ups_mut('.lame_uv_nc_rho_draw', '- shape * log(om)', '- 4 * shape * log(om)')), 5)
})

test_that('the persistence with the paths following keeps the posterior of bipartite snap fits and under the prior on the stationary variance', {
	expect_lt(upnc_law(1500, 9228, 'poisson', 'bip', kind = 'snap'), 4.5)
	expect_lt(upnc_law(1500, 9229, 'binary', 'bip', kind = 'snap', kappa = NULL), 4.5)
	expect_lt(upnc_law(1500, 9230, 'normal', 'bip', kind = 'snap', stationary = TRUE), 4.5)
	expect_lt(upnc_law(1500, 9232, 'poisson', 'dir', stationary = TRUE), 4.5)
	expect_lt(upnc_law(1500, 9233, 'binary', 'dir', kind = 'snap', kappa = NULL, stationary = TRUE), 4.5)
	# the checks fail with the resets of the bipartite paths drawn as drifts and with the innovation prior read
	# under the prior on the stationary variance
	as_drift = function(...) {
		a = list(...)
		a$kind = 'ar1'
		do.call(lame:::.lame_uv_nc_rho_draw, a)
	}
	expect_gt(upnc_law(6000, 9234, 'binary', 'bip', kind = 'snap', kappa = 0.4, draw = as_drift), 5)
	expect_gt(upnc_law(1500, 9235, 'normal', 'dir', stationary = TRUE, stationary_draw = FALSE), 5)
})

test_that('lame() runs the size and persistence draws with their switches on and a bipartite Poisson panel passes the two-chain check', {
	set.seed(9231)
	nA = 12L; nB = 10L; N = 4L; rho = 0.9
	ar = function(m) {
		X = array(0, c(m, 2, N))
		X[, , 1] = rnorm(2 * m, 0, 0.6)
		for (t in 2:N) X[, , t] = rho * X[, , t - 1] + sqrt(1 - rho^2) * rnorm(2 * m, 0, 0.6)
		X
	}
	U = ar(nA); V = ar(nB)
	Y = lapply(1:N, function(t) {
		y = matrix(rpois(nA * nB, exp(U[, , t] %*% t(V[, , t]) + rnorm(nA * nB, 0, 0.3))), nA)
		dimnames(y) = list(sprintf('r%02d', 1:nA), sprintf('c%02d', 1:nB))
		y
	})
	names(Y) = paste0('t', 1:N)
	fit = function(prior, chains = 2, nscan = 1500) suppressMessages(lame(
		Y, family = 'poisson', mode = 'bipartite', R_row = 2L, R_col = 2L, dynamic_uv = TRUE, burn = 300, nscan = nscan,
		odens = 15, seed = 6886, verbose = FALSE, gof = FALSE, plot = FALSE, prior = prior, chains = chains, cores = 1
	))
	f = fit(list())
	expect_s3_class(f, 'lame')
	expect_gt(f$uv_size_move[['moved']], 0)
	expect_gt(f$uv_rho_nc_move[['moved']], 0)
	g = fit(list(uv_size_move = FALSE, uv_rho_nc_move = FALSE), chains = 1, nscan = 150)
	expect_null(g$uv_size_move)
	expect_null(g$uv_rho_nc_move)
	expect_error(fit(list(uv_size_move = 'yes'), chains = 1, nscan = 150), 'uv_size_move')
})
