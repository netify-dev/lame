# hyperparameters of the actor-specific slopes (R/per_actor_hyper_move.R): exact conditionals of the innovation
# variance and persistence at every sweep, plus an interweaving step on (variance, persistence, paths)
# with the standardised innovations held fixed

# the maps

test_that('the path map keeps the standardised innovations, inverts itself and keeps the zero sum over actors', {
	set.seed(4)
	n = 6L; nt = 5L
	w = matrix(rnorm(n * nt), n, nt)
	w = sweep(w, 2, colMeans(w))
	th = lame:::.lame_pa_rebuild(w, 0.6, 0.35)
	expect_equal(colSums(th), rep(0, nt), tolerance = 1e-10)
	expect_equal(lame:::.lame_pa_innovations(th, 0.6, 0.35), w, tolerance = 1e-10)
	moved = lame:::.lame_pa_rebuild(w, -0.3, 0.8)
	expect_equal(colSums(moved), rep(0, nt), tolerance = 1e-10)
	expect_equal(lame:::.lame_pa_innovations(moved, -0.3, 0.8), w, tolerance = 1e-10)
	# the stationary first state: variance sigma^2 / (1 - rho^2) per unit innovation
	expect_equal(th[, 1], 0.35 / sqrt(1 - 0.6^2) * w[, 1], tolerance = 1e-12)
})

test_that('the likelihood of the paths is the one the path sweep reads', {
	set.seed(5)
	n = 4L; nt = 3L
	th = matrix(rnorm(n * nt), n, nt); H = matrix(rgamma(n * nt, 2), n, nt); h = matrix(rnorm(n * nt), n, nt)
	expect_equal(lame:::.lame_pa_loglik(th, H, h, NULL), sum(h * th) - 0.5 * sum(H * th^2))
	P = lapply(seq_len(nt), function(t) { m = matrix(rnorm(n * n, 0, 0.3), n); diag(H[, t]) + crossprod(m) })
	ref = sum(vapply(seq_len(nt), function(t) sum(h[, t] * th[, t]) - 0.5 * drop(t(th[, t]) %*% P[[t]] %*% th[, t]), 0))
	expect_equal(lame:::.lame_pa_loglik(th, H, h, P), ref)
})

# exactness

# one step from the joint law: sigma^2 ~ IG(shape, scale), (rho + 1) / 2 ~ Beta(a, b) on (-0.99, 0.99),
# centred stationary AR(1) actor paths and Gaussian statistics h = P theta + noise
pa_geweke = function(n, nt, coupled, ndraw, seed, step = lame:::.lame_pa_nc_step) {
	set.seed(seed)
	shape = 3; scale = 0.4; ra = 2; rb = 2
	H = matrix(rgamma(n * nt, 2, 1), n, nt)
	Pl = NULL
	if (coupled) { Pl = lapply(
		seq_len(nt),
		function(t) { m = matrix(rnorm(n * n, 0, 0.25), n); diag(H[, t]) + crossprod(m) }
	); for (t in seq_len(nt)) H[, t] = diag(Pl[[t]]) }
	Lp = if (coupled) lapply(Pl, function(P) t(chol(P))) else NULL
	stat_of = function(s) {
		th = s$th
		c(
			log(s$s2),
			s$s2,
			s$rho,
			s$rho^2,
			s$rho * log(s$s2),
			colSums(th^2),
			colSums(th[, -1, drop = FALSE] * th[, -nt, drop = FALSE]),
			colSums(s$h * th),
			th[1, ],
			th[1, ]^2
		)
	}
	D = vector('list', ndraw)
	worst_sum = 0
	for (it in seq_len(ndraw)) {
		s2 = 1 / rgamma(1, shape, scale)
		repeat { r = 2 * rbeta(1, ra, rb) - 1; if (abs(r) < 0.99) break }
		th = matrix(0, n, nt)
		for (i in seq_len(n)) { th[i, 1] = rnorm(1, 0, sqrt(s2 / (1 - r^2))); for (t in seq_len(nt)[-1]) th[i, t] = r * th[i, t - 1] + rnorm(1, 0, sqrt(s2)) }
		th = sweep(th, 2, colMeans(th))
		h = matrix(0, n, nt)
		for (t in seq_len(nt)) h[, t] = if (coupled) drop(
			Pl[[t]] %*% th[, t] + Lp[[t]] %*% rnorm(n)
		) else H[, t] * th[, t] + sqrt(H[, t]) * rnorm(n)
		s = list(s2 = s2, rho = r, th = th, h = h)
		out = step(th, H, h, Pl, r, s2, shape, scale, ra, rb, tries = 4L)
		s1 = s; s1$th = out$theta; s1$s2 = out$sigma2; s1$rho = out$rho
		worst_sum = max(worst_sum, max(abs(colSums(s1$th))))
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
	list(maxz = max(abs(z)), T2 = T2, df = sum(kp), p = pchisq(T2, sum(kp), lower.tail = FALSE), zero_sum = worst_sum)
}

test_that('one interweaving step leaves the joint law of hyperparameters, paths and statistics unchanged', {
	skip_on_cran()
	for (cfg in list(
		list(6L, 4L, FALSE, 11),
		list(6L, 4L, TRUE, 12),
		list(9L, 5L, FALSE, 13),
		list(5L, 1L, FALSE, 14)
	)) {
		r = pa_geweke(cfg[[1]], cfg[[2]], cfg[[3]], 3500, cfg[[4]])
		expect_lt(r$zero_sum, 1e-8)
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('the exactness check fails for a step that leaves out the stationary first state or accepts every proposal', {
	skip_on_cran()
	mut = function(name, from, to) {
		src = deparse(get(name, asNamespace('lame')))
		out = sub(from, to, src, fixed = TRUE)
		stopifnot(!identical(out, src))
		fn = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
		environment(fn) = asNamespace('lame')
		fn
	}
	env = new.env(parent = asNamespace('lame'))
	assign(
		'.lame_pa_rebuild',
		mut('.lame_pa_rebuild', 'th[, 1L] <- sigma/sqrt(1 - rho^2) * w[, 1L]', 'th[, 1L] <- sigma * w[, 1L]'),
		envir = env
	)
	step = lame:::.lame_pa_nc_step
	environment(step) = env
	r = pa_geweke(6L, 4L, FALSE, 2500, 21, step = step)
	expect_gt(r$maxz, 6)
	step = mut('.lame_pa_nc_step', 'log(stats::runif(1L)) < ', 'TRUE || ')
	r = pa_geweke(6L, 4L, FALSE, 2500, 22, step = step)
	expect_gt(r$maxz, 6)
})

# the switches and the fits

pa_fit = function(nscan = 60, burn = 20, kind = 'row', ident = 'center', mode = 'unipartite', ...) {
	dat = simulate_test_network(n = 8, n_time = 4, family = 'normal', mode = mode, nA = 8, nB = 6, seed = 3)
	suppressWarnings(suppressMessages(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = 'normal',
		mode = mode,
		R = 0,
		dynamic_beta_per_actor = kind,
		per_actor_identifiability = ident,
		keep_per_actor = 'summary',
		nscan = nscan,
		burn = burn,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5,
		...
	)))
}

test_that('actor_hyper_nc is validated', {
	expect_error(pa_fit(prior = list(actor_hyper_nc = NA)), 'actor_hyper_nc must be TRUE or FALSE', fixed = TRUE)
})

test_that('the hyperparameters are drawn at every sweep from the first', {
	skip_on_cran()
	on = pa_fit()
	# a new value at nearly every stored sweep
	expect_gt(length(unique(on$SIGMA_ACTOR)), 0.9 * length(on$SIGMA_ACTOR))
	# the step is recorded when it ran and only then
	expect_true(!is.null(on$pa_hyper_move))
	expect_equal(unname(on$pa_hyper_move[['sweeps']]), 80)
	expect_gt(on$pa_hyper_move[['accepts']], 0)
	expect_lt(on$pa_hyper_move[['rate']], 1)
	expect_null(pa_fit(prior = list(actor_hyper_nc = FALSE))$pa_hyper_move)
	every = pa_fit(prior = list(actor_hyper_nc = FALSE))
	expect_gt(length(unique(every$SIGMA_ACTOR)), 0.9 * length(every$SIGMA_ACTOR))
	# the hyperparameters start moving at the first sweep: a fit with no burn-in has a draw different from the start values
	nb = pa_fit(nscan = 20, burn = 0)
	expect_false(isTRUE(all.equal(nb$SIGMA_ACTOR[1], sqrt(0.1))))
})

test_that('every identifiability, margin and mode runs with the step and keeps the paths centred', {
	skip_on_cran()
	for (cfg in list(
		list('row', 'center', 'unipartite'),
		list('col', 'exact_center', 'unipartite'),
		list('row', 'center', 'bipartite')
	)) {
		fit = pa_fit(nscan = 30, burn = 10, kind = cfg[[1]], ident = cfg[[2]], mode = cfg[[3]])
		expect_true(all(is.finite(fit$SIGMA_ACTOR)))
		expect_gt(fit$pa_hyper_move[['attempts']], 0)
		expect_lt(max(abs(colSums(fit$theta_actor_mean))), 1e-8)
	}
})
