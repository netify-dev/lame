# the column scales and shears of the common gauge of dynamic position paths (R/uv_gauge_dynamic.R)
skip_on_cran()

ugy_mut = function(fn, from, to) {
	f = get(fn, envir = asNamespace('lame'))
	out = deparse(f, width.cutoff = 500L)
	hit = grepl(from, out, fixed = TRUE)
	if (!any(hit)) stop('mutation target not found: ', from)
	out[hit] = sub(from, to, out[hit], fixed = TRUE)
	g = eval(parse(text = paste(out, collapse = '\n')), envir = asNamespace('lame'))
	environment(g) = asNamespace('lame')
	g
}

# precision weights and persistence indicators of the transitions of each kind
ugy_weights = function(n, N, kind, sigma, nu = 4, kappa = 1.3, pi_snap = 0.3) {
	w = matrix(1, n, N)
	k = matrix(1, n, N)
	if (kind == 't') w[, -1] = rgamma(n * (N - 1L), nu / 2, nu / 2)
	if (kind == 'snap') {
		d = matrix(rbinom(n * (N - 1L), 1, pi_snap), n)
		k[, -1] = 1 - d
		w[, -1] = ifelse(d == 1, sigma^2 / kappa^2, 1)
	}
	list(w = w, k = k)
}

# paths from the prior: x_1 ~ N(0, sigma^2 / (w_1 (1 - rho^2))), x_t ~ N(k_t rho x_{t-1}, sigma^2 / w_t)
ugy_paths = function(n, R, N, rho, sigma, wk) {
	X = array(0, c(n, R, N))
	X[, , 1] = rnorm(n * R) * sigma / sqrt(wk$w[, 1] * (1 - rho^2))
	for (t in 2:N) X[, , t] = wk$k[, t] * rho * X[, , t - 1] + rnorm(n * R) * sigma / sqrt(wk$w[, t])
	X
}

# gauge-sensitive statistics of a side and their squares
ugy_side = function(X) {
	S = crossprod(matrix(aperm(X, c(1, 3, 2)), ncol = dim(X)[2]))
	C = cov2cor(S)
	c(log(diag(S)), C[upper.tri(C)])
}
ugy_z = function(out) max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(nrow(out)))))

ugy_law = function(ndraw, seed, n = 6L, R = 2L, N = 3L, rho = 0.8, sigma = 0.6, kind = 'ar1', max_abs = Inf, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_gauge_dynamic_draw
	set.seed(seed)
	stat = function(U, V) { x = c(ugy_side(U), ugy_side(V), log(sum(U^2) / sum(V^2))); c(x, x^2) }
	out = NULL
	for (r in seq_len(ndraw)) {
		repeat {
			wu = ugy_weights(n, N, kind, sigma); wv = ugy_weights(n, N, kind, sigma)
			U = ugy_paths(n, R, N, rho, sigma, wu); V = ugy_paths(n, R, N, rho, sigma, wv)
			if (max(abs(U), abs(V)) <= max_abs) break
		}
		o = list(U = U, V = V)
		for (k in 1:3) o = draw(o$U, o$V, rho, sigma, wu$w, wv$w, wu$k, wv$k, max_abs = max_abs)
		for (t in seq_len(N)) if (max(abs(o$U[, , t] %*% t(o$V[, , t]) - U[, , t] %*% t(V[, , t]))) > 1e-9) stop('the product moved')
		if (max(abs(o$U), abs(o$V)) > max_abs) stop('the bound was left')
		s1 = stat(o$U, o$V) - stat(U, V)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	ugy_z(out)
}

test_that('the column scales and shears of directed paths keep the law of the paths under each kind and the bound', {
	expect_lt(ugy_law(2500, 9101), 4.5)
	expect_lt(ugy_law(2500, 9102, R = 3L, kind = 't'), 4.5)
	expect_lt(ugy_law(2500, 9103, kind = 'snap'), 4.5)
	expect_lt(ugy_law(2500, 9104, max_abs = 2.5), 4.5)
	# the checks fail for a wrong sign in the shear mean and a wrong exponent of the column scale
	expect_gt(ugy_law(2500, 9105, draw = ugy_mut(
		'.lame_uv_gauge_dynamic_draw',
		'-(Mu[k, l] - Mv[k, l])/prec',
		'-(Mu[k, l] + Mv[k, l])/prec'
	)), 5)
	expect_gt(ugy_law(2500, 9106, draw = ugy_mut(
		'.lame_uv_gauge_dynamic_draw',
		'exp(-2 * u) * b)/(2 * s2)',
		'exp(-2 * u) * b)/(2 * s2) + 2 * u'
	)), 5)
})

ugy_bip_law = function(ndraw, seed, nA = 6L, nB = 5L, RA = 2L, RB = 2L, N = 3L, rho = 0.7, sigma = 0.7, g_var = 1.5,
	g_max = Inf, max_abs = Inf, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_gauge_bip_dynamic_draw
	set.seed(seed)
	stat = function(U, V, G) {
		x = c(ugy_side(U), ugy_side(V), log(rowSums(G^2)), log(colSums(G^2)), if (RA > 1L) sum(G[1, ] * G[2, ]))
		c(x, x^2)
	}
	one = function(n) list(w = matrix(1, n, N), k = matrix(1, n, N))
	out = NULL
	for (r in seq_len(ndraw)) {
		repeat {
			U = ugy_paths(nA, RA, N, rho, sigma, one(nA)); V = ugy_paths(nB, RB, N, rho, sigma, one(nB))
			G = matrix(rnorm(RA * RB, 0, sqrt(g_var)), RA)
			if (max(abs(U), abs(V)) <= max_abs && max(abs(G)) <= g_max) break
		}
		o = list(U = U, V = V, G = G)
		for (k in 1:3) o = draw(o$U, o$V, o$G, rho, sigma, g_var = g_var, g_max = g_max, max_abs = max_abs)
		for (t in seq_len(N)) {
			p0 = U[, , t] %*% G %*% t(V[, , t]); p1 = o$U[, , t] %*% o$G %*% t(o$V[, , t])
			if (max(abs(p1 - p0)) > 1e-9) stop('the product moved')
		}
		if (max(abs(o$G)) > g_max || max(abs(o$U), abs(o$V)) > max_abs) stop('the bound was left')
		s1 = stat(o$U, o$V, o$G) - stat(U, V, G)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	ugy_z(out)
}

test_that('the gauge of bipartite paths against G keeps the law of the paths and G, also cut at the bounds', {
	expect_lt(ugy_bip_law(2500, 9111), 4.5)
	expect_lt(ugy_bip_law(2500, 9112, RA = 3L, RB = 1L), 4.5)
	expect_lt(ugy_bip_law(2500, 9113, g_max = 2, max_abs = 2.5), 4.5)
	expect_gt(ugy_bip_law(2500, 9114, draw = ugy_mut(
		'.lame_uv_gauge_bip_dynamic_draw',
		'lam <- N * dim(X)[1L] - ncol(Gm)',
		'lam <- N * dim(X)[1L]'
	)), 5)
	expect_gt(ugy_bip_law(2500, 9115, draw = ugy_mut(
		'.lame_uv_gauge_bip_dynamic_draw',
		'- sum(Gm[k, ] * Gm[l, ])/g_var)',
		'+ sum(Gm[k, ] * Gm[l, ])/g_var)'
	)), 5)
})

ugy_sym_law = function(ndraw, seed, n = 6L, R = 2L, N = 3L, rho = 0.8, l_var = 2, max_abs = Inf, draw = NULL) {
	if (is.null(draw)) draw = lame:::.lame_uv_sym_split_dynamic_draw
	set.seed(seed)
	sigma = sqrt(1 - rho^2)
	stat = function(U, L) { x = c(ugy_side(U), log(L^2)); c(x, x^2) }
	one = list(w = matrix(1, n, N), k = matrix(1, n, N))
	out = NULL
	for (r in seq_len(ndraw)) {
		repeat {
			U = ugy_paths(n, R, N, rho, sigma, one)
			if (max(abs(U)) <= max_abs) break
		}
		L = rnorm(R, 0, sqrt(l_var))
		o = list(U = U, Ld = L)
		for (k in 1:3) o = draw(o$U, o$Ld, rho, sigma, l_var, max_abs = max_abs)
		for (t in seq_len(N)) {
			p0 = U[, , t] %*% diag(L) %*% t(U[, , t]); p1 = o$U[, , t] %*% diag(o$Ld) %*% t(o$U[, , t])
			if (max(abs(p1 - p0)) > 1e-9) stop('the product moved')
		}
		s1 = stat(o$U, o$Ld) - stat(U, L)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	ugy_z(out)
}

test_that('the split of symmetric paths against their weights keeps the law of the paths and the weights', {
	expect_lt(ugy_sym_law(2500, 9121), 4.5)
	expect_lt(ugy_sym_law(2500, 9122, max_abs = 2), 4.5)
	expect_gt(ugy_sym_law(2500, 9123, draw = ugy_mut(
		'.lame_uv_sym_split_dynamic_draw',
		'lam <- dU[1L] * dU[3L] - 2',
		'lam <- dU[1L] * dU[3L]'
	)), 5)
})

test_that('the gauge draws keep a degenerate state and lame() runs them only with the switch on', {
	ns = asNamespace('lame')
	U = array(rnorm(12), c(3, 2, 2)); U[1] = NA
	expect_identical(ns$.lame_uv_gauge_dynamic_draw(U, U, .5, 1)$moved, 0L)
	expect_identical(ns$.lame_uv_gauge_dynamic_draw(U[, 1, , drop = FALSE], U[, 1, , drop = FALSE], .5, 1)$moved, 0L)
	set.seed(9131)
	n = 10L; rn = sprintf('a%02d', 1:n)
	Y = lapply(1:3, function(t) { y = matrix(rnorm(n * n), n) + tcrossprod(matrix(rnorm(2 * n), n), matrix(rnorm(2 * n), n)); diag(y) = NA; dimnames(y) = list(rn, rn); y })
	names(Y) = paste0('t', 1:3)
	calls = 0L
	gd = lame:::.lame_uv_gauge_dynamic_draw
	local_mocked_bindings(
		.lame_uv_gauge_dynamic_draw = function(...) { calls <<- calls + 1L; gd(...) },
		.lame_chain_moving = function(fit, ...) data.frame(parameter = character(0), text = character(0)),
		.package = 'lame'
	)
	fit = function(prior) suppressWarnings(suppressMessages(lame(
		Y, family = 'normal', R = 2L, dynamic_uv = TRUE, burn = 5, nscan = 10, odens = 1, seed = 2,
		verbose = FALSE, gof = FALSE, plot = FALSE, prior = prior
	)))
	fit(list())
	expect_gt(calls, 0L)
	calls = 0L
	fit(list(uv_split_move = FALSE))
	expect_equal(calls, 0L)
})
