# exact draw of the orientation of dynamic position paths, rotations and reflections
skip_on_cran()

# copies of the draw and its helpers with one function mutated
urd_mut = function(from, to, fn = '.lame_uv_rotation_dynamic_draw') {
	ns = asNamespace('lame')
	fe = new.env(parent = ns)
	helpers = grep('^\\.lame_uvrd_', ls(ns, all.names = TRUE), value = TRUE)
	for (nm in c('.lame_uv_rotation_dynamic_draw', helpers)) {
		f = get(nm, envir = ns)
		if (nm == fn) {
			out = deparse(f, width.cutoff = 500L)
			for (i in seq_along(from)) {
				hit = grepl(from[i], out, fixed = TRUE)
				if (!any(hit)) stop('mutation target not found: ', from[i])
				out[hit] = sub(from[i], to[i], out[hit], fixed = TRUE)
			}
			f = eval(parse(text = paste(out, collapse = '\n')))
		}
		environment(f) = fe
		assign(nm, f, envir = fe)
	}
	get('.lame_uv_rotation_dynamic_draw', envir = fe)
}

# precision weights w and persistence indicators k of the transitions under each kind
urd_kappa = 1.3
urd_pi = 0.3
urd_snap_wk = function(d, sigma, kappa = urd_kappa) {
	k = 1 - d
	k[, 1] = 1
	list(w = ifelse(k == 1, 1, sigma^2 / kappa^2), k = k, d = d)
}
urd_weights = function(n, N, kind, sigma, nu = 4) {
	w = matrix(1, n, N)
	k = matrix(1, n, N)
	if (kind == 't') w[, -1] = rgamma(n * (N - 1L), nu / 2, nu / 2)
	if (kind %in% c('snap', 'snapc')) {
		d = cbind(0, matrix(rbinom(n * (N - 1L), 1, urd_pi), n))
		return(urd_snap_wk(d, sigma))
	}
	list(w = w, k = k)
}

# paths from the prior: x_1 ~ N(0, sigma^2 / (w_1 (1 - rho^2))), x_t ~ N(k_t rho x_{t-1}, sigma^2 / w_t)
urd_paths = function(n, R, N, rho, sigma, wk) {
	X = array(0, c(n, R, N))
	X[, , 1] = rnorm(n * R) * sigma / sqrt(wk$w[, 1] * (1 - rho^2))
	for (t in 2:N) {
		X[, , t] = wk$k[, t] * rho * X[, , t - 1] + rnorm(n * R) * sigma / sqrt(wk$w[, t])
	}
	X
}

# quadratic form of the prior
urd_q = function(X, wk, rho) {
	n = dim(X)[1]
	N = dim(X)[3]
	sl = function(t) matrix(X[, , t], n)
	q = sum(wk$w[, 1] * (1 - rho^2) * rowSums(sl(1)^2))
	for (t in 2:N) q = q + sum(wk$w[, t] * rowSums((sl(t) - wk$k[, t] * rho * sl(t - 1))^2))
	q
}

# statistics that see the orientation: per-period spread, per-transition alignment, turn and handedness
urd_stat = function(X, wk, rho) {
	n = dim(X)[1]
	R = dim(X)[2]
	N = dim(X)[3]
	sl = function(t) matrix(X[, , t], n)
	out = urd_q(X, wk, rho)
	for (t in seq_len(N)) {
		out = c(out, sum(sl(t)[, 1]^2 - sl(t)[, 2]^2), sum(sl(t)[, 1] * sl(t)[, 2]))
	}
	for (t in 2:N) {
		C = crossprod(sl(t - 1), sl(t))
		out = c(out, sum(diag(C)), sign(det(C)))
		for (p in 1:(R - 1)) out = c(out, C[p, p + 1] - C[p + 1, p])
	}
	out
}

# z scores of one step of the move from prior paths, cut at max_abs by rejection.
# kind snap holds the resets fixed, snapc lets the move sum them out and draw them again
urd_law = function(
	ndraw,
	seed,
	nA = 7L,
	nB = 6L,
	R = 2L,
	N = 3L,
	rho = 0.8,
	sigma = 0.6,
	kind = 'ar1',
	max_abs = Inf,
	draw = lame:::.lame_uv_rotation_dynamic_draw
) {
	set.seed(seed)
	out = NULL
	moved = 0
	top = 0
	flips = 0
	for (r in seq_len(ndraw)) {
		repeat {
			wu = urd_weights(nA, N, kind, sigma)
			wv = urd_weights(nB, N, kind, sigma)
			U = urd_paths(nA, R, N, rho, sigma, wu)
			V = urd_paths(nB, R, N, rho, sigma, wv)
			if (max(abs(U), abs(V)) <= max_abs) break
		}
		if (kind == 'snapc') {
			snap = list(pi = urd_pi, kappa = urd_kappa, delta_u = wu$d, delta_v = wv$d)
			o = draw(U, V, rho, sigma, max_abs = max_abs, snap = snap)
			wu1 = urd_snap_wk(o$delta_u, sigma)
			wv1 = urd_snap_wk(o$delta_v, sigma)
			s0 = c(urd_stat(U, wu, rho), urd_stat(V, wv, rho), colSums(wu$d)[-1], colSums(wv$d)[-1])
			s1 = c(
				urd_stat(o$U, wu1, rho),
				urd_stat(o$V, wv1, rho),
				colSums(wu1$d)[-1],
				colSums(wv1$d)[-1]
			) - s0
		} else {
			o = draw(U, V, rho, sigma, wu$w, wv$w, wu$k, wv$k, max_abs = max_abs)
			s0 = c(urd_stat(U, wu, rho), urd_stat(V, wv, rho))
			s1 = c(urd_stat(o$U, wu, rho), urd_stat(o$V, wv, rho)) - s0
		}
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
		moved = moved + o$moved
		top = max(top, abs(o$U), abs(o$V))
		# count transitions whose handedness changed, which only a reflection does
		hand = function(X, t) sign(det(crossprod(matrix(X[, , t - 1], nA), matrix(X[, , t], nA))))
		for (t in 2:N) flips = flips + (hand(o$U, t) != hand(U, t))
		# the product of every period is untouched
		for (t in seq_len(N)) {
			p0 = tcrossprod(matrix(U[, , t], nA), matrix(V[, , t], nB))
			p1 = tcrossprod(matrix(o$U[, , t], nA), matrix(o$V[, , t], nB))
			if (max(abs(p1 - p0)) > 1e-12) stop('the product moved')
		}
	}
	# drop statistics the draws never changed
	sdv = apply(out, 2, sd)
	z = colMeans(out[, sdv > 0, drop = FALSE]) / (sdv[sdv > 0] / sqrt(ndraw))
	list(maxz = max(abs(z)), moved = moved / ndraw, top = top, flips = flips / ndraw)
}

test_that('one step of the rotation leaves the law of the paths under their AR(1) prior unchanged', {
	cfgs = list(
		list(seed = 9201),
		list(seed = 9202, N = 4L, kind = 't'),
		list(seed = 9203, N = 4L, kind = 'snap'),
		list(seed = 9204, N = 4L, kind = 'snapc'),
		list(seed = 9205, R = 3L, N = 4L, rho = 0.9, sigma = 0.5),
		list(seed = 9206, R = 3L, N = 4L, kind = 't'),
		list(seed = 9207, R = 3L, N = 4L, kind = 'snap'),
		list(seed = 9208, R = 3L, N = 4L, kind = 'snapc')
	)
	for (cfg in cfgs) {
		r = do.call(urd_law, c(list(ndraw = 4000), cfg))
		expect_lt(r$maxz, 4)
		expect_gt(r$moved, 1)
	}
})

test_that('the broken versions of the move change the law of the paths', {
	# M_t with the wrong sign, transposed
	f = urd_mut('M <- rho/s2', 'M <- -rho/s2', '.lame_uvrd_step')
	expect_gt(urd_law(3000, 9211, draw = f)$maxz, 5)
	f = urd_mut(
		'B <- M[pq, pq, drop = FALSE]',
		'B <- t(M[pq, pq, drop = FALSE])',
		'.lame_uvrd_step'
	)
	expect_gt(urd_law(3000, 9213, draw = f)$maxz, 5)
	# summed-out resets: mirror always taken, mirror not self-inverse, resets not redrawn, wrong reset sd
	f = urd_mut('log(stats::runif(1)) >= ll_p - ll', 'FALSE', '.lame_uvrd_step')
	expect_gt(urd_law(3000, 9217, N = 4L, kind = 'snapc', draw = f)$maxz, 5)
	f = urd_mut('tcrossprod(refl, P)', 'refl', '.lame_uvrd_mirror')
	expect_gt(urd_law(3000, 9214, N = 4L, kind = 'snapc', draw = f)$maxz, 5)
	f = urd_mut(
		'st$du[, t] <- .lame_uvrd_reset',
		'st$du[, t] <- st$du[, t] + 0 * .lame_uvrd_reset',
		'.lame_uvrd_step'
	)
	expect_gt(urd_law(3000, 9218, N = 4L, kind = 'snapc', draw = f)$maxz, 5)
	f = urd_mut('k2 <- as.numeric(snap$kappa)^2', 'k2 <- 4 * as.numeric(snap$kappa)^2')
	expect_gt(urd_law(3000, 9219, N = 4L, kind = 'snapc', draw = f)$maxz, 5)
	# snap resets and t weights ignored
	f = urd_mut('k_u <- one(k_u, dU[1L])', 'k_u <- one(NULL, dU[1L])')
	expect_gt(urd_law(3000, 9215, N = 4L, kind = 'snap', draw = f)$maxz, 5)
	f = urd_mut('w_u <- one(w_u, dU[1L])', 'w_u <- one(NULL, dU[1L])')
	expect_gt(urd_law(3000, 9216, N = 4L, R = 3L, kind = 't', draw = f)$maxz, 5)
})

test_that('the move reaches the reflections, which the conditional draws cannot', {
	# the move flips the handedness of transitions; a mutant confined to rotations never does
	r = urd_law(1500, 9231, nA = 3L, nB = 3L, N = 4L, rho = 0.5)
	expect_gt(r$flips, 0.2)
	expect_gt(urd_law(1500, 9232, nA = 3L, nB = 3L, N = 4L, rho = 0.5, kind = 'snapc')$flips, 0.2)
	f = urd_mut(
		'.lame_uv_o2(B)',
		'{ E2 <- .lame_uv_o2(B); if (det(E2) < 0) diag(2) else E2 }',
		'.lame_uvrd_step'
	)
	expect_identical(urd_law(1500, 9231, nA = 3L, nB = 3L, N = 4L, rho = 0.5, draw = f)$flips, 0)
})

test_that('the move respects the coordinate bound and keeps the law of the bounded prior', {
	r = urd_law(3000, 9221, max_abs = 2.5)
	expect_lt(r$maxz, 4)
	expect_lte(r$top, 2.5)
	r = urd_law(3000, 9223, N = 4L, kind = 'snapc', max_abs = 2.5)
	expect_lt(r$maxz, 4)
	expect_lte(r$top, 2.5)
	f = urd_mut('if (top > max_abs)', 'if (FALSE)', '.lame_uvrd_step')
	expect_gt(urd_law(1000, 9222, max_abs = 2.5, draw = f)$top, 2.5)
})

urd_still = function(fit, ...) data.frame(parameter = character(0), text = character(0))
urd_uni = function(n = 10L, N = 3L, seed = 1, R = 2L, sd = 1, noise = 1) {
	set.seed(seed)
	nm = sprintf('a%02d', seq_len(n))
	U = matrix(rnorm(n * R, 0, sd), n)
	V = matrix(rnorm(n * R, 0, sd), n)
	Y = vector('list', N)
	for (t in seq_len(N)) {
		if (t > 1) {
			U = 0.9 * U + sqrt(0.19) * matrix(rnorm(n * R, 0, sd), n)
			V = 0.9 * V + sqrt(0.19) * matrix(rnorm(n * R, 0, sd), n)
		}
		y = U %*% t(V) + matrix(rnorm(n * n, 0, noise), n)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
	}
	names(Y) = paste0('p', seq_len(N))
	Y
}
urd_fit = function(Y, ..., prior = list(), burn = 20, nscan = 60, R = 2L, odens = 1, seed = 5) {
	suppressWarnings(suppressMessages(lame(
		Y = Y,
		R = R,
		dynamic_uv = TRUE,
		burn = burn,
		nscan = nscan,
		odens = odens,
		seed = seed,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

test_that('lame() draws the orientation of directed paths under every kind, records it, and the switch omits it', {
	local_mocked_bindings(.lame_chain_moving = urd_still, .package = 'lame')
	Y = urd_uni()
	for (kind in c('ar1', 't', 'snap')) {
		on = urd_fit(Y, dynamic_uv_kind = kind)
		expect_equal(unname(on$uv_rotation_move['visits']), 80)
		expect_gt(unname(on$uv_rotation_move['moved']), 0)
		off = urd_fit(Y, dynamic_uv_kind = kind, prior = list(uv_rotation_move = FALSE))
		expect_null(off$uv_rotation_move)
		expect_false(identical(on$VC, off$VC))
	}
	f3 = urd_fit(urd_uni(R = 3L), R = 3L, dynamic_uv_kind = 'snap')
	expect_gt(unname(f3$uv_rotation_move['moved']), 0)
	# rank-one, symmetric and bipartite fits have no rotation move
	expect_null(urd_fit(Y, R = 1L)$uv_rotation_move)
	Ys = lapply(Y, function(y) (y + t(y)) / 2)
	expect_null(urd_fit(Ys, symmetric = TRUE)$uv_rotation_move)
	Yb = lapply(Y, function(y) {
		y = y[1:8, 3:10]
		dimnames(y) = list(sprintf('r%02d', 1:8), sprintf('c%02d', 1:8))
		y
	})
	fb = suppressWarnings(suppressMessages(lame(
		Y = Yb,
		mode = 'bipartite',
		R_row = 2L,
		R_col = 2L,
		dynamic_uv = TRUE,
		burn = 20,
		nscan = 60,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	expect_null(fb$uv_rotation_move)
	expect_error(urd_fit(Y, prior = list(uv_rotation_move = 'yes')), 'uv_rotation_move')
})

test_that('a two-chain fit keeps the record of each chain', {
	local_mocked_bindings(
		.lame_chain_moving = urd_still,
		.lame_gate_table = function(fits) data.frame(parameter = character(0), rhat = numeric(0)),
		.package = 'lame'
	)
	f2 = urd_fit(urd_uni(), dynamic_uv_kind = 'snap', chains = 2)
	rec = lapply(f2$sampler_diagnostics_by_chain, `[[`, 'uv_rotation_move')
	expect_length(rec, 2L)
	for (r in rec) {
		expect_equal(unname(r['visits']), 80)
		expect_gt(unname(r['moved']), 0)
	}
})

test_that('a snap chain does not stay in a reflected period with the move', {
	# rank-two AR(1) data: with the move the chain keeps the handedness of every transition, pi_snap stays at
	# its prior level and the split of every period mixes. (Without the rotation and whole-path moves the chain
	# of the snap build stayed in a reflected tail at this seed; with the gauge moves of the paths it no longer
	# does at seeds 1 to 4 or 10, so that arm is no control and is not run.)
	Y = urd_uni(15L, 5L, seed = 5, sd = 0.8, noise = 0.5)
	ess = function(x) {
		x = x - mean(x)
		a = acf(x, lag.max = 300, plot = FALSE)$acf[-1]
		k = which(a < .05)[1]
		if (is.na(k)) k = 300
		length(x) / (1 + 2 * sum(a[seq_len(k - 1L)]))
	}
	run = function(on) {
		f = urd_fit(
			Y,
			dynamic_uv_kind = 'snap',
			burn = 300,
			nscan = 1500,
			seed = 10,
			prior = list(snap_kappa = 4, uv_rotation_move = on, actor_move_paths = on, chain_extensions = 0)
		)
		d = dim(f$U_draws)
		reflected = function(t, s) det(crossprod(f$U_draws[, , t - 1, s], f$U_draws[, , t, s])) < 0
		log_split = function(t, s) log(sqrt(sum(f$U_draws[, , t, s]^2)) / sqrt(sum(f$V_draws[, , t, s]^2)))
		refl = sapply(2:d[3], function(t) mean(sapply(seq_len(d[4]), function(s) reflected(t, s))))
		split = sapply(seq_len(d[3]), function(t) sapply(seq_len(d[4]), function(s) log_split(t, s)))
		list(
			refl = refl,
			pi = mean(f$pi_snap),
			ess_pi = ess(f$pi_snap),
			ess_split = apply(split, 2, ess)
		)
	}
	on = run(TRUE)
	expect_identical(max(on$refl), 0)
	expect_lt(on$pi, 0.015)
	expect_gt(on$ess_pi, 500)
	expect_true(all(on$ess_split > 500))
})

test_that('a rank-two snap fit reports, monitors and plots rho_uv, sigma_uv and pi_snap', {
	Y = urd_uni(12L, 4L, seed = 5, sd = 0.8, noise = 0.5)
	f = urd_fit(
		Y,
		dynamic_uv_kind = 'snap',
		chains = 2,
		burn = 100,
		nscan = 600,
		odens = 10,
		# a persistent 4-period panel, fitted under a prior that says so
		prior = list(snap_kappa = 4, rho_uv_mean = 0.9, rho_uv_sd = 0.1)
	)
	hyper = c('rho_uv', 'sigma_uv', 'pi_snap')
	s = summary(f)
	expect_equal(s$rho_uv, mean(f$rho_uv))
	expect_equal(s$sigma_uv, mean(f$sigma_uv))
	expect_equal(s$pi_snap, mean(f$pi_snap))
	out = gsub('\\s+', ' ', paste(capture.output(print(s)), collapse = ' '))
	expect_false(grepl('not reported', out))
	# the two-chain gate monitors the three
	expect_true(all(hyper %in% f$rhat$parameter))
	expect_true(all(f$rhat$rhat[f$rhat$parameter %in% hyper] < 1.05))
	expect_true(all(hyper %in% colnames(lame:::.lame_gate_series(f))))
	p = trace_plot(f, params = 'dynamic')
	d = if (inherits(p, 'patchwork')) p[[1]]$data else p$data
	expect_setequal(
		unique(as.character(d$parameter)),
		c('Position persistence (rho_uv)', 'Position innovation (sigma_uv)', 'Position reset rate (pi_snap)')
	)
})
