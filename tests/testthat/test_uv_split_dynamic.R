# exact draw of the scale split of dynamic position paths
skip_on_cran()

usd_mut = function(from, to) {
	f = lame:::.lame_uv_split_dynamic_draw
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# precision weights w and persistence indicators k of the transitions under each kind
usd_weights = function(n, N, kind, sigma, nu = 4, kappa = 1.3, pi_snap = 0.3) {
	w = matrix(1, n, N)
	k = matrix(1, n, N)
	if (N > 1L && kind == 't') w[, -1] = rgamma(n * (N - 1L), nu / 2, nu / 2)
	if (N > 1L && kind == 'snap') {
		d = matrix(rbinom(n * (N - 1L), 1, pi_snap), n)
		k[, -1] = 1 - d
		w[, -1] = ifelse(d == 1, sigma^2 / kappa^2, 1)
	}
	list(w = w, k = k)
}

# paths from the prior: x_1 ~ N(0, sigma^2 / (w_1 (1 - rho^2))), x_t ~ N(k_t rho x_{t-1}, sigma^2 / w_t)
usd_paths = function(n, R, N, rho, sigma, wk) {
	X = array(0, c(n, R, N))
	X[, , 1] = rnorm(n * R) * sigma / sqrt(wk$w[, 1] * (1 - rho^2))
	if (N > 1L) for (t in 2:N) {
		X[, , t] = wk$k[, t] * rho * X[, , t - 1] + rnorm(n * R) * sigma / sqrt(wk$w[, t])
	}
	X
}

# quadratic form of the prior
usd_q = function(X, wk, rho) {
	n = dim(X)[1]
	N = dim(X)[3]
	sl = function(t) matrix(X[, , t], n)
	q = sum(wk$w[, 1] * (1 - rho^2) * rowSums(sl(1)^2))
	if (N > 1L) for (t in 2:N) q = q + sum(wk$w[, t] * rowSums((sl(t) - wk$k[, t] * rho * sl(t - 1))^2))
	q
}

# z scores and hotelling p value of one step of the move from prior paths, cut at max_abs by rejection
usd_law = function(
	ndraw,
	seed,
	nA = 7L,
	RA = 2L,
	nB = 5L,
	RB = 1L,
	N = 3L,
	rho = 0.8,
	sigma = 0.6,
	kind = 'ar1',
	max_abs = Inf,
	draw = lame:::.lame_uv_split_dynamic_draw
) {
	set.seed(seed)
	G = matrix(rnorm(RA * RB), RA, RB)
	stat = function(U, V, wu, wv) {
		su = apply(U^2, 3, sum)
		sv = apply(V^2, 3, sum)
		c(
			su,
			sv,
			0.5 * log(su / sv),
			usd_q(U, wu, rho),
			usd_q(V, wv, rho),
			sum(U[, , 1] * U[, , min(2L, N)])
		)
	}
	out = NULL
	moved = 0
	top = 0
	for (r in seq_len(ndraw)) {
		repeat {
			wu = usd_weights(nA, N, kind, sigma)
			wv = usd_weights(nB, N, kind, sigma)
			U = usd_paths(nA, RA, N, rho, sigma, wu)
			V = usd_paths(nB, RB, N, rho, sigma, wv)
			if (max(abs(U), abs(V)) <= max_abs) break
		}
		o = draw(U, V, rho, sigma, wu$w, wv$w, wu$k, wv$k, max_abs = max_abs)
		s1 = stat(o$U, o$V, wu, wv) - stat(U, V, wu, wv)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
		moved = moved + o$moved
		top = max(top, abs(o$U), abs(o$V))
		# the product of every period is untouched
		for (t in seq_len(N)) {
			p0 = matrix(U[, , t], nA) %*% G %*% t(matrix(V[, , t], nB))
			p1 = matrix(o$U[, , t], nA) %*% G %*% t(matrix(o$V[, , t], nB))
			if (max(abs(p1 - p0)) > 1e-9) stop('the product moved')
		}
	}
	z = colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))
	ev = eigen(cov(out), symmetric = TRUE)
	m = drop(crossprod(ev$vectors, colMeans(out)))
	T2 = ndraw * sum(m^2 / ev$values)
	list(maxz = max(abs(z)), p = pchisq(T2, ncol(out), lower.tail = FALSE), moved = moved / ndraw, top = top)
}

test_that('the broken versions of the move change the law of the paths', {
	# jacobian without the periods, or without the second side
	expect_gt(usd_law(7000, 9101, draw = usd_mut('N * lam * u', 'lam * u'))$maxz, 5)
	f = usd_mut('lam <- dU[1L] * dU[2L] - dV[1L] * dV[2L]', 'lam <- dU[1L] * dU[2L]')
	expect_gt(usd_law(7000, 9102, draw = f)$maxz, 5)
	# first period at the transition precision, not the stationary one
	expect_gt(usd_law(7000, 9103, draw = usd_mut('p1 <- 1 - rho * rho', 'p1 <- 1'))$maxz, 5)
	# innovation variance left out
	expect_gt(usd_law(7000, 9104, draw = usd_mut('/(2 * s2)', '/2'))$maxz, 5)
	# pull of the neighbouring periods left out of the per-period step
	expect_gt(usd_law(7000, 9105, draw = usd_mut('- 2 * bu * exp(u)', ''))$maxz, 5)
	# t and snap weights ignored
	f = usd_mut('w_u <- one(w_u, dU[1L])', 'w_u <- one(NULL, dU[1L])')
	expect_gt(usd_law(7000, 9106, N = 4L, kind = 't', draw = f)$maxz, 5)
	f = usd_mut('k_u <- one(k_u, dU[1L])', 'k_u <- one(NULL, dU[1L])')
	expect_gt(usd_law(7000, 9107, N = 4L, kind = 'snap', draw = f)$maxz, 5)
})

test_that('one step of the split leaves the law of the paths under their AR(1) prior unchanged', {
	cfgs = list(
		# bipartite, every kind of weight
		list(seed = 9001),
		list(seed = 9002, nA = 4L, RA = 1L, nB = 6L, RB = 2L, N = 2L, rho = 0.5, sigma = 1.4),
		list(seed = 9003, N = 4L, kind = 't'),
		list(seed = 9004, N = 4L, kind = 'snap'),
		# unipartite (U and V of the same size, no jacobian)
		list(seed = 9005, nA = 6L, RA = 2L, nB = 6L, RB = 2L, N = 4L, rho = 0.9, sigma = 0.5),
		list(seed = 9006, nA = 6L, RA = 2L, nB = 6L, RB = 2L, N = 4L, rho = 0.9, sigma = 0.5, kind = 't'),
		list(seed = 9007, nA = 6L, RA = 2L, nB = 6L, RB = 2L, N = 4L, rho = 0.9, sigma = 0.5, kind = 'snap')
	)
	for (cfg in cfgs) {
		r = do.call(usd_law, c(list(ndraw = 7000), cfg))
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
		expect_gt(r$moved, 1)
	}
})

test_that('the move respects the coordinate bound and keeps the law of the bounded prior', {
	r = usd_law(7000, 9011, max_abs = 2.2)
	expect_lt(r$maxz, 4.5)
	expect_gt(r$p, 1e-3)
	expect_lte(r$top, 2.2)
	# the bound ignored: draws leave the support
	f = usd_mut('lm <- lims(max(abs(U_cube)), max(abs(V_cube)))', 'lm <- c(-Inf, Inf)')
	expect_gt(usd_law(3000, 9012, max_abs = 1.6, draw = f)$top, 1.6)
})

usd_still = function(fit, ...) data.frame(parameter = character(0), text = character(0))
usd_uni = function(n = 10L, N = 3L, seed = 1, family = 'normal') {
	set.seed(seed)
	nm = sprintf('a%02d', seq_len(n))
	U = matrix(rnorm(n * 2), n)
	V = matrix(rnorm(n * 2), n)
	Y = vector('list', N)
	for (t in seq_len(N)) {
		U = 0.9 * U + sqrt(0.19) * matrix(rnorm(n * 2), n)
		V = 0.9 * V + sqrt(0.19) * matrix(rnorm(n * 2), n)
		y = U %*% t(V) + matrix(rnorm(n * n), n)
		y = switch(
			family,
			normal = y,
			binary = 1 * (y > 0),
			poisson = matrix(rpois(n * n, exp(y / 2)), n),
			ordinal = matrix(as.numeric(cut(y, c(-Inf, -1, 0, 1, Inf))), n)
		)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
	}
	names(Y) = paste0('p', seq_len(N))
	Y
}
usd_bip = function(nA = 10L, nB = 8L, N = 3L, seed = 2) {
	set.seed(seed)
	U = matrix(rnorm(nA * 2), nA)
	V = matrix(rnorm(nB * 2), nB)
	Y = vector('list', N)
	for (t in seq_len(N)) {
		U = 0.9 * U + sqrt(0.19) * matrix(rnorm(nA * 2), nA)
		V = 0.9 * V + sqrt(0.19) * matrix(rnorm(nB * 2), nB)
		y = U %*% t(V) + matrix(rnorm(nA * nB), nA)
		dimnames(y) = list(sprintf('r%02d', seq_len(nA)), sprintf('c%02d', seq_len(nB)))
		Y[[t]] = y
	}
	names(Y) = paste0('p', seq_len(N))
	Y
}
usd_fit = function(Y, ..., bip = FALSE, prior = list(), burn = 20, nscan = 60, chains = 1) {
	args = list(
		Y = Y,
		dynamic_uv = TRUE,
		burn = burn,
		nscan = nscan,
		odens = 1,
		seed = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		chains = chains,
		...
	)
	rank_args = if (bip) list(mode = 'bipartite', R_row = 2L, R_col = 2L) else list(R = 2L)
	args = c(args, rank_args)
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

test_that('lame() draws the split of dynamic paths in every mode, family and kind, records it, and the switch omits it', {
	local_mocked_bindings(.lame_chain_moving = usd_still, .package = 'lame')
	cases = list(
		list(Y = usd_bip(), bip = TRUE),
		list(Y = usd_uni()),
		list(Y = usd_uni(family = 'binary'), family = 'binary'),
		list(Y = usd_uni(family = 'poisson'), family = 'poisson'),
		list(Y = usd_uni(family = 'ordinal'), family = 'ordinal'),
		list(Y = usd_uni(), dynamic_uv_kind = 't'),
		list(Y = usd_uni(), dynamic_uv_kind = 'snap')
	)
	for (cs in cases) {
		on = do.call(usd_fit, cs)
		expect_equal(unname(on$uv_split_move['visits']), 80)
		expect_gt(unname(on$uv_split_move['moved']), 0)
		off = do.call(usd_fit, c(cs, list(prior = list(uv_split_move = FALSE))))
		expect_null(off$uv_split_move)
		# an ordinal fit has no intercept, so its variance components carry the difference
		draws = function(f) if (ncol(f$BETA)) f$BETA else f$VC
		expect_false(identical(draws(on), draws(off)))
	}
	# a symmetric fit splits each column of its paths against its weight (R/uv_gauge_dynamic.R)
	Ys = lapply(usd_uni(), function(y) { y = (y + t(y)) / 2; y })
	expect_equal(unname(usd_fit(Ys, symmetric = TRUE)$uv_split_move['visits']), 80)
	expect_null(usd_fit(Ys, symmetric = TRUE, prior = list(uv_split_move = FALSE))$uv_split_move)
})

test_that('a bipartite fit with unequal ranks and absent actors draws the split', {
	local_mocked_bindings(.lame_chain_moving = usd_still, .package = 'lame')
	Y = usd_bip()
	# a row actor absent from the second period, a column actor from the third
	Y[[2]] = Y[[2]][-1, , drop = FALSE]
	Y[[3]] = Y[[3]][, -2, drop = FALSE]
	fit_bip = function(prior = list()) suppressWarnings(suppressMessages(lame(
		Y = Y,
		mode = 'bipartite',
		R_row = 2L,
		R_col = 1L,
		dynamic_uv = TRUE,
		burn = 20,
		nscan = 60,
		odens = 1,
		seed = 5,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior
	)))
	on = fit_bip()
	expect_identical(c(on$R_row, on$R_col), c(2L, 1L))
	expect_equal(unname(on$uv_split_move['visits']), 80)
	expect_gt(unname(on$uv_split_move['moved']), 0)
	expect_true(all(is.finite(on$BETA)))
	off = fit_bip(list(uv_split_move = FALSE))
	expect_null(off$uv_split_move)
	expect_false(identical(on$BETA, off$BETA))
})

test_that('a two-chain fit keeps the record of each chain, and a resumed fit draws the split', {
	local_mocked_bindings(
		.lame_chain_moving = usd_still,
		.lame_gate_table = function(fits) data.frame(parameter = character(0), rhat = numeric(0)),
		.package = 'lame'
	)
	f2 = usd_fit(usd_bip(), bip = TRUE, chains = 2)
	expect_identical(f2$chains, 2L)
	rec = lapply(f2$sampler_diagnostics_by_chain, `[[`, 'uv_split_move')
	expect_length(rec, 2L)
	for (r in rec) {
		expect_equal(unname(r['visits']), 80)
		expect_gt(unname(r['moved']), 0)
	}
	f1 = usd_fit(usd_uni())
	r1 = suppressWarnings(suppressMessages(lame_resume(f1, nscan_more = 30)))
	expect_equal(unname(r1$uv_split_move['visits']), 30)
	expect_gt(unname(r1$uv_split_move['moved']), 0)
})

test_that('the split of dynamic paths mixes at nominal speed with the move', {
	local_mocked_bindings(.lame_chain_moving = usd_still, .package = 'lame')
	Y = usd_bip(20L, 15L, 5L, seed = 4)
	ess = function(x) {
		x = x - mean(x)
		a = acf(x, lag.max = 300, plot = FALSE)$acf[-1]
		k = which(a < .05)[1]
		if (is.na(k)) k = 300
		length(x) / (1 + 2 * sum(a[seq_len(k - 1L)]))
	}
	split = function(f) {
		d = dim(f$U_draws)
		sapply(seq_len(d[3]), function(t) vapply(
			seq_len(d[4]),
			function(s) {
				log(sqrt(sum(f$U_draws[, , t, s]^2)) / sqrt(sum(f$V_draws[, , t, s]^2)))
			},
			0
		))
	}
	on = split(usd_fit(Y, bip = TRUE, burn = 200, nscan = 2000))
	off = split(usd_fit(
		Y,
		bip = TRUE,
		prior = list(uv_split_move = FALSE),
		burn = 200,
		nscan = 2000
	))
	e_on = apply(on, 2, ess)
	e_off = apply(off, 2, ess)
	expect_true(all(e_on > 3 * e_off))
	expect_true(all(e_on > 500))
})

test_that('the bipartite article fit settles after at most one longer block', {
	Y_list = readRDS(test_path('fixtures', 'bipartite_article_Y_list.rds'))
	fit = suppressWarnings(suppressMessages(lame(
		Y = Y_list,
		mode = 'bipartite',
		R_row = 2,
		R_col = 2,
		family = 'normal',
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		burn = 500,
		nscan = 2000,
		odens = 20,
		chains = 2,
		verbose = FALSE,
		plot = FALSE
	)))
	# 100 stored draws per chain are too few to compare some series, so one longer block may follow
	expect_lte(fit$auto_extend$extensions, 1L)
	expect_true(all(fit$rhat$informative, na.rm = TRUE))
	# every monitored series under its gate threshold (1.05, or the noise bound of a short run)
	expect_true(all(fit$rhat$rhat <= fit$rhat$threshold, na.rm = TRUE))
	expect_lt(max(fit$rhat$rhat, na.rm = TRUE), 1.06)
})
