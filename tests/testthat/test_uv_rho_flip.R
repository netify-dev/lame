# exact mirror move of the persistence of dynamic latent paths (R/uv_rho_flip.R)

# stationary AR(1) paths of one side
urf_paths = function(n, R, N, rho, sigma) {
	X = array(0, c(n, R, N))
	X[, , 1] = rnorm(n * R, 0, sigma / sqrt(1 - rho^2))
	for (t in seq_len(N)[-1]) X[, , t] = rho * X[, , t - 1] + rnorm(n * R, 0, sigma)
	X
}

# log density of the stationary AR(1) prior of the paths
urf_log_prior = function(X, rho, sigma) {
	N = dim(X)[3]
	lp = sum(dnorm(X[, , 1], 0, sigma / sqrt(1 - rho^2), log = TRUE))
	for (t in seq_len(N)[-1]) lp = lp + sum(dnorm(X[, , t], rho * X[, , t - 1], sigma, log = TRUE))
	lp
}

# alternating period signs
urf_mirror = function(X) sweep(X, 3, rep_len(c(1, -1), dim(X)[3]), '*')

test_that('the mirror keeps every product and the path prior, and the move takes a mirror state back', {
	set.seed(11)
	U = urf_paths(6, 2, 5, 0.85, 0.4)
	V = urf_paths(5, 2, 5, 0.85, 0.4)
	G = matrix(rnorm(4), 2, 2)
	Um = urf_mirror(U)
	Vm = urf_mirror(V)
	for (t in 1:5) expect_equal(Um[, , t] %*% G %*% t(Vm[, , t]), U[, , t] %*% G %*% t(V[, , t]))
	expect_equal(urf_log_prior(Um, -0.85, 0.4), urf_log_prior(U, 0.85, 0.4))
	# from the mirror the prior ratio exceeds one, so every uniform takes the move
	back = lame:::.lame_uv_rho_flip(Um, Vm, -0.85, 0.5, 0.2, 0.999)
	expect_true(back$moved)
	expect_identical(back$rho, 0.85)
	expect_identical(back$U, U)
	expect_identical(back[["V"]], V)
	# from the positive state the move is taken with probability exp(-2 rho m / s^2)
	lr = -2 * 0.85 * 0.5 / 0.2^2
	expect_false(lame:::.lame_uv_rho_flip(U, V, 0.85, 0.5, 0.2, exp(lr) * 1.01)$moved)
	expect_true(lame:::.lame_uv_rho_flip(U, V, 0.85, 0.5, 0.2, exp(lr) * 0.99)$moved)
	# a symmetric fit carries one set of paths
	sym = lame:::.lame_uv_rho_flip(Um, NULL, -0.85, 0.9, 0.1, 0.5)
	expect_identical(sym$U, U)
	expect_null(sym[["V"]])
	# without a uniform the move is not tried
	for (u in list(NULL, NA_real_)) {
		skip_move = lame:::.lame_uv_rho_flip(Um, Vm, -0.85, 0.5, 0.2, u)
		expect_false(skip_move$visited)
		expect_identical(skip_move$U, Um)
		expect_identical(skip_move$rho, -0.85)
	}
})

test_that('a recycled metropolis uniform is uniform on either side of the decision', {
	skip_on_cran()
	set.seed(12)
	u = runif(20000)
	for (la in log(c(0.03, 0.3, 0.9))) {
		accepted = log(u) < la
		r = mapply(lame:::.lame_recycle_u, u, la, accepted)
		for (side in list(accepted, !accepted)) {
			expect_gt(suppressWarnings(ks.test(r[side], 'punif'))$p.value, 1e-3)
		}
	}
	expect_identical(lame:::.lame_recycle_u(0.4, 0.5, TRUE), 0.4)
})

# chain without data: exact path draws given the persistence, the package draw of the persistence
# with its recycled uniform, then the mirror move. rho then follows its prior
urf_prior_chain = function(n_iter, kind, prior_mean, prior_sd, flip = lame:::.lame_uv_rho_flip) {
	nA = 3; nB = 2; R = 1; N = 4
	rho = prior_mean
	sigma = if (kind == 'unit') sqrt(1 - rho^2) else 0.5
	out = numeric(n_iter)
	for (it in seq_len(n_iter)) {
		U = urf_paths(nA, R, N, rho, sigma)
		V = if (kind == 'unit') NULL else urf_paths(nB, R, N, rho, sigma)
		d = switch(
			kind,
			slice = lame:::.lame_uv_rho_w_step(U, V, rho, sigma, prior_mean, prior_sd, 2, 1),
			unit = lame:::.lame_rho_uv_unit_mh(U, rho, prior_mean, prior_sd, with_u = TRUE),
			drift = lame:::.lame_sample_rho_uv_drift_only(
				U,
				V,
				sigma,
				rho,
				FALSE,
				matrix(0, nA, N),
				matrix(0, nB, N),
				prior_mean = prior_mean,
				prior_sd = prior_sd,
				with_u = TRUE
			)
		)
		if (kind == 'slice') sigma = d$sigma
		rho = flip(U, V, d$rho, prior_mean, prior_sd, d$u)$rho
		if (kind == 'unit') sigma = sqrt(1 - rho^2)
		out[it] = rho
	}
	out
}

# standardized gap between the chain's share of negative rho and the prior's, batch-means error
urf_gap = function(x, prior_mean, prior_sd) {
	lo = pnorm(-0.99, prior_mean, prior_sd)
	p0 = (pnorm(0, prior_mean, prior_sd) - lo) / (pnorm(0.99, prior_mean, prior_sd) - lo)
	neg = matrix(x < 0, ncol = 50)
	(mean(neg) - p0) / (sd(colMeans(neg)) / sqrt(50))
}

# the move with its log ratio halved
urf_wrong_flip = function() {
	f = lame:::.lame_uv_rho_flip
	src = deparse(f, width.cutoff = 500L)
	hit = grepl('-2 * rho * prior_mean/prior_sd^2', src, fixed = TRUE)
	stopifnot(sum(hit) == 1L)
	src[hit] = sub('-2 * rho * prior_mean/prior_sd^2', '-rho * prior_mean/prior_sd^2', src[hit], fixed = TRUE)
	g = eval(parse(text = paste(src, collapse = '\n')))
	environment(g) = asNamespace('lame')
	g
}

test_that('without data the persistence keeps its prior under each draw and the move, and a wrong ratio does not', {
	skip_on_cran()
	m = 0.15
	s = 0.3
	for (kind in c('slice', 'unit', 'drift')) {
		set.seed(13)
		expect_lt(abs(urf_gap(urf_prior_chain(5000, kind, m, s), m, s)), 4)
		set.seed(13)
		expect_gt(abs(urf_gap(urf_prior_chain(5000, kind, m, s, urf_wrong_flip()), m, s)), 6)
	}
})

# directed bipartite gaussian panel with AR(1) positions
urf_bip = function(nA = 10L, nB = 8L, N = 4L, seed = 3) {
	set.seed(seed)
	U = matrix(rnorm(nA * 2), nA)
	V = matrix(rnorm(nB * 2), nB)
	Y = vector('list', N)
	for (t in seq_len(N)) {
		if (t > 1) {
			U = 0.9 * U + sqrt(0.19) * matrix(rnorm(nA * 2), nA)
			V = 0.9 * V + sqrt(0.19) * matrix(rnorm(nB * 2), nB)
		}
		Y[[t]] = matrix(
			U %*% t(V) + rnorm(nA * nB, 0, 0.5),
			nA,
			nB,
			dimnames = list(sprintf('r%02d', seq_len(nA)), sprintf('c%02d', seq_len(nB)))
		)
	}
	names(Y) = paste0('t', seq_len(N))
	Y
}

urf_fit = function(Y, ..., prior = list(), burn = 30, nscan = 30) {
	suppressWarnings(suppressMessages(lame(
		Y = Y,
		mode = 'bipartite',
		R_row = 2,
		R_col = 2,
		family = 'normal',
		dynamic_uv = TRUE,
		burn = burn,
		nscan = nscan,
		odens = 1,
		seed = 7,
		chains = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

test_that('a chain started at the mirror of the persistence returns to positive rho_uv, and records it', {
	skip_on_cran()
	local_mocked_bindings(
		.lame_chain_moving = function(fit, ...) data.frame(parameter = character(0), text = character(0)),
		.package = 'lame'
	)
	Y = urf_bip()
	f0 = urf_fit(Y)
	expect_gt(min(f0$rho_uv), 0)
	expect_equal(unname(f0$uv_rho_flip_move['visits']), 60)
	# the mirror of the final state: negative persistence, alternating period signs
	sv = f0$start_vals
	sv$checkpoint_state = NULL
	sv$rho_uv = -abs(sv$rho_uv)
	sv$U = urf_mirror(sv$U)
	sv$V = urf_mirror(sv$V)
	on = urf_fit(Y, start_vals = sv, burn = 0)
	expect_gt(min(on$rho_uv), 0)
	expect_gte(unname(on$uv_rho_flip_move['moved']), 1)
	# the conditional draws alone keep the sampler at the mirror; its draws are stored with rho_uv >= 0
	off = urf_fit(Y, start_vals = sv, burn = 0, prior = list(uv_rho_flip_move = FALSE))
	expect_lt(off$start_vals$rho_uv, 0)
	expect_gte(min(off$rho_uv), 0)
	expect_null(off$uv_rho_flip_move)
	expect_error(urf_fit(Y, prior = list(uv_rho_flip_move = 'yes')), 'uv_rho_flip_move')
})
