# dense Gaussian precision and determinant oracles, independent of the sampler's AR gradients and SVD factors
skip_on_cran()

binary_move_fixture = function(rank, uv_prior = NULL) {
	n = 4L; nt = 3L; k = n - 1L
	Y = array(NA_real_, c(n, n, nt))
	Y[cbind(c(1, 2, 3, 4, 1, 2), c(2, 3, 1, 2, 4, 1), c(1, 1, 2, 2, 3, 3))] = c(1, 0, 0, 1, 0, 0)
	X = lapply(seq_len(nt), function(t) {
		x = array(rnorm(n * n * 3), c(n, n, 3)); x[, , 1] = 1; x
	})
	h = list(
		g = 1.2,
		rho_beta = c(1, .4),
		var_beta = c(.2, .3),
		iC0 = solve(matrix(c(.6, .1, .1, .8), 2)),
		iSab = solve(matrix(c(.5, .15, .15, .4), 2)),
		rho_ab = .65,
		var_ab = .4
	)
	dims = list(
		static = 1L,
		path = c(nt, 2L),
		level = c(k, 2L),
		adev = c(nt, k),
		bdev = c(nt, k),
		U = c(n, rank),
		V = c(n, rank)
	)
	cursor = 0L
	indices = lapply(dims, function(d) {
		ix = cursor + seq_len(prod(d)); cursor <<- cursor + prod(d)
		if (length(d) > 1L) dim(ix) = d
		ix
	})
	restore = function(q) lapply(indices, function(ix) {
		z = q[ix]; if (!is.null(dim(ix))) dim(z) = dim(ix); z
	})
	Q = matrix(0, cursor, cursor)
	Q[indices$static, indices$static] = 1 / h$g
	add_ar = function(ix, rho, variance, initial) {
		Q[ix[1], ix[1]] <<- Q[ix[1], ix[1]] + initial
		for (t in 2:length(ix)) {
			v = numeric(cursor); v[ix[t]] = 1; v[ix[t - 1]] = -rho
			Q <<- Q + tcrossprod(v) / variance
		}
	}
	Q[indices$path[1, ], indices$path[1, ]] = h$iC0
	for (j in 1:2) add_ar(indices$path[, j], h$rho_beta[j], h$var_beta[j], 0)
	for (i in seq_len(k)) {
		Q[indices$level[i, ], indices$level[i, ]] = h$iSab
		for (block in c("adev", "bdev")) add_ar(
			indices[[block]][, i],
			h$rho_ab,
			h$var_ab,
			(1 - h$rho_ab^2) / h$var_ab
		)
	}
	H = qr.Q(qr(contr.helmert(n)))
	density = function(q) {
		z = restore(q); levels = H %*% z$level
		a = H %*% t(z$adev); b = H %*% t(z$bdev); uv = tcrossprod(z$U, z$V)
		lp = 0
		for (t in seq_len(nt)) {
			eta = matrix(matrix(X[[t]], n * n, 3) %*% c(z$static, z$path[t, ]), n, n) +
				outer(levels[, 1] + a[, t], levels[, 2] + b[, t], "+") + uv
			observed = is.finite(Y[, , t])
			lp = lp + sum(pnorm((2 * Y[, , t][observed] - 1) * eta[observed], log.p = TRUE))
		}
		W = cbind(z$U, z$V)
		# factor covariance prior of R/uv_prior.R: df degrees of freedom, scale psi I (both 2R + 2 by default)
		df0 = if (is.null(uv_prior)) 2 * rank + 2 else uv_prior$df
		psi0 = if (is.null(uv_prior)) 2 * rank + 2 else uv_prior$psi
		lp - drop(q %*% Q %*% q) / 2 - (n + df0) / 2 *
			as.numeric(determinant(psi0 * diag(2 * rank) + crossprod(W), logarithm = TRUE)$modulus)
	}
	prior_draw = function() {
		path = matrix(0, nt, 2)
		path[1, ] = drop(rnorm(2) %*% chol(solve(h$iC0)))
		for (t in 2:nt) path[t, ] = h$rho_beta * path[t - 1, ] + rnorm(2) * sqrt(h$var_beta)
		ar = function() {
			z = matrix(0, nt, k); z[1, ] = rnorm(k, sd = sqrt(h$var_ab / (1 - h$rho_ab^2)))
			for (t in 2:nt) z[t, ] = h$rho_ab * z[t - 1, ] + rnorm(k, sd = sqrt(h$var_ab))
			z
		}
		S = solve(rWishart(1, 2 * rank + 2, diag(2 * rank) / (2 * rank + 2))[, , 1])
		W = matrix(rnorm(n * 2 * rank), n, 2 * rank) %*% chol(S)
		unlist(
			list(
				rnorm(1, sd = sqrt(h$g)),
				path,
				matrix(rnorm(k * 2), k, 2) %*% chol(solve(h$iSab)),
				ar(),
				ar(),
				W[, seq_len(rank), drop = FALSE],
				W[, rank + seq_len(rank), drop = FALSE]
			),
			use.names = FALSE
		)
	}
	list(
		context = lame:::.lame_binary_hmc_context(Y, X, 1, 2:3, rank, uv_prior),
		hyper = h,
		density = density,
		draw = prior_draw,
		restore = restore,
		size = cursor
	)
}

test_that("binary marginal gradients and trajectories agree with independent density", {
	set.seed(341535)
	# the unit prior and the data-scaled prior of a sparse network (R/uv_prior.R: df = round((2R + 2) m), psi = df s)
	for (rank in 1:2) for (k in c(0, 7)) {
		f = binary_move_fixture(
			rank,
			if (k > 0) list(df = 2 * rank + 2 + k, psi = (2 * rank + 2 + k) * .3)
		); q = rnorm(f$size, sd = .3)
		ans = lame:::.lame_binary_hmc_evaluate(q, f$context, f$hyper)
		q2 = rnorm(f$size, sd = .3)
		expect_equal(
			ans$logp - lame:::.lame_binary_hmc_evaluate(q2, f$context, f$hyper)$logp,
			f$density(q) - f$density(q2),
			tolerance = 1e-10
		)
		eps = 1e-5
		gradient = vapply(seq_along(q), function(i) {
			a = b = q; a[i] = a[i] + eps; b[i] = b[i] - eps
			(f$density(a) - f$density(b)) / (2 * eps)
		}, numeric(1))
		expect_equal(ans$gradient, gradient, tolerance = 1e-6)
		mass = lame:::.lame_binary_hmc_mass(f$context, f$hyper)
		p = rnorm(length(q)) * sqrt(mass)
		forward = lame:::.lame_binary_hmc_integrate(q, p, f$context, f$hyper, mass, .12, 7)
		backward = lame:::.lame_binary_hmc_integrate(
			forward$q,
			-forward$p,
			f$context,
			f$hyper,
			mass,
			.12,
			7
		)
		expect_equal(backward$q, q, tolerance = 1e-10)
		expect_equal(backward$p, -p, tolerance = 1e-10)
	}
})

test_that("the complete binary marginal update preserves independent prior draws", {
	skip_on_cran()
	set.seed(493817)
	for (rank in 1:2) {
		f = binary_move_fixture(rank); f$context$obs[] = FALSE
		statistic = function(q) {
			z = f$restore(q)
			c(
				z$static,
				z$path[3, 1],
				z$level[1, 1],
				z$adev[2, 1],
				atan(sum(z$U[1, ] * z$V[2, ])),
				log1p(sum(z$U^2) + sum(z$V^2))
			)
		}
		delta = replicate(800, {
			q = f$draw()
			z = lame:::.lame_binary_hmc_step(q, f$context, f$hyper, .15, L = 4)
			statistic(z$q) - statistic(q)
		})
		z = rowMeans(delta) / (apply(delta, 1, sd) / sqrt(ncol(delta)))
		expect_true(all(is.finite(z)))
		expect_lt(max(abs(z)), 6)
	}
})

test_that("static factor scaling includes the covariance and coefficient path priors", {
	set.seed(759381)
	nt = 4L; p = 2L; path = matrix(rnorm(nt * p), nt, p)
	static = rnorm(3); rho = c(-.3, 1); sigma = c(.2, .7)
	Lambda = diag(c(.3, 2)); C0 = matrix(c(2, .4, .4, 1), 2); g = c(.5, 2, 4); s2 = 1.3
	Q = matrix(0, nt * p, nt * p); first = 1 + (0:(p - 1)) * nt
	Q[first, first] = solve(C0)
	for (j in 1:p) for (t in 2:nt) {
		v = numeric(nt * p); v[(j - 1) * nt + t] = 1; v[(j - 1) * nt + t - 1] = -rho[j]
		Q = Q + tcrossprod(v) / (sigma[j]^2 * Lambda[j, j])
	}
	terms = lame:::.lame_path_scale_terms(static, path, rho, sigma, Lambda, C0, g, s2)
	expect_equal(terms$q, sum(static^2 / (g * s2)) + drop(c(path) %*% Q %*% c(path)), tolerance = 1e-10)
	expect_equal(terms$dimension, length(static) + length(path))
	for (rank in 1:2) for (shrink in c(FALSE, TRUE)) {
		n = 5L; W = matrix(rnorm(n * 2 * rank), n, 2 * rank); x = rnorm(17); f = .6
		prior = function(x, W) sum(dnorm(x, log = TRUE)) + if (shrink) {
			-(n + 2 * rank + 2) / 2 * as.numeric(
				determinant((2 * rank + 2) * diag(2 * rank) + crossprod(W), logarithm = TRUE)$modulus
			)
		} else sum(dnorm(W, sd = sqrt(n), log = TRUE))
		oracle = prior(f * x, sqrt(f) * W) - prior(x, W) + (length(x) + length(W) / 2) * log(f)
		expect_equal(
			lame:::.lame_static_scale_logratio(
				f,
				sum(x^2),
				length(x),
				eigen(crossprod(W), symmetric = TRUE, only.values = TRUE)$values,
				n,
				rank,
				shrink
			),
			oracle,
			tolerance = 1e-10
		)
	}
})

test_that("binary joint moves retain aligned draws and resumable tuning", {
	skip_on_cran()
	set.seed(361853); n = 8L; nt = 3L
	X = lapply(
		seq_len(nt),
		function(t) array(
			rnorm(n * n),
			c(n, n, 1),
			dimnames = list(as.character(1:n), as.character(1:n), "x")
		)
	)
	Y = lapply(seq_len(nt), function(t) {
		z = matrix(
			rbinom(n * n, 1, .2),
			n,
			n,
			dimnames = list(as.character(1:n), as.character(1:n))
		)
		diag(z) = NA; z[1, 2] = 1; z[2, 1] = 0; z
	})
	args = list(
		Y = Y,
		Xdyad = X,
		R = 2,
		family = "binary",
		dcor = FALSE,
		dynamic_ab = TRUE,
		dynamic_beta = "x_dyad",
		dynamic_beta_kind = "rw1",
		burn = 80,
		nscan = 160,
		odens = 5,
		seed = 71833,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		posterior_opts = list(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1),
		checkpoint_path = tempfile(fileext = ".rds"),
		checkpoint_every = 20
	)
	checkpoint_path = args$checkpoint_path
	on.exit(unlink(checkpoint_path), add = TRUE)
	rng = .Random.seed
	fit = suppressMessages(do.call(lame, args))
	expect_identical(.Random.seed, rng)
	expect_gt(unname(fit$marginal_move["attempts"]), 0)
	state = lame:::.joint_predictive_state(fit)
	means = lapply(
		seq_len(nt),
		function(t) Reduce(
			`+`,
			lapply(
				state$index,
				function(s) lame:::.joint_predictor(fit, state, s, t)
			)
		) / length(state$index)
	)
	reported = unlist(fit$EZ); expected = unlist(means); keep = is.finite(reported)
	expect_equal(unname(expected[keep]), unname(reported[keep]), tolerance = 1e-10)
	checkpoint = readRDS(args$checkpoint_path)
	step = checkpoint$start_vals$checkpoint_state$values$marginal_sd
	args$checkpoint_path = NULL; args$checkpoint_every = NULL
	args$start_vals = checkpoint$start_vals; args$burn = 0; args$nscan = 40
	resumed = suppressMessages(do.call(lame, args))
	expect_identical(unname(resumed$marginal_move["step_sd"]), step)
	args$start_vals$checkpoint_state$values$marginal_sd = NULL
	args$start_vals$marginal_move_sd = NULL
	legacy = suppressMessages(do.call(lame, args))
	expect_identical(unname(legacy$marginal_move["step_sd"]), .12)
	args$start_vals = NULL; args$prior = list(marginal_move = FALSE)
	off = suppressMessages(do.call(lame, args))
	expect_null(off$marginal_move)
	args$start_vals = checkpoint$start_vals
	args$nscan = 10; args$odens = 1; args$prior = list(level_move_every = 1)
	first = suppressMessages(do.call(lame, args))
	expect_gt(unname(first$level_move["attempts"]), 0)
	args$prior = list(marginal_move = NA)
	expect_error(suppressMessages(do.call(lame, args)), "marginal_move must be TRUE or FALSE", fixed = TRUE)
})

test_that("the joint scale move runs for binary fits with static latent positions", {
	set.seed(9101); n = 10; nm = sprintf("a%02d", 1:n)
	Y = lapply(1:3, function(t) { y = matrix(rbinom(n * n, 1, 0.2), n, dimnames = list(nm, nm)); diag(y) = NA; y })
	X = lapply(1:3, function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, "x")))
	fit = lame(
		Y,
		Xdyad = X,
		family = "binary",
		R = 2,
		dynamic_ab = TRUE,
		burn = 20,
		nscan = 60,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		prior = list(level_move_every = 1)
	)
	expect_gt(unname(fit$level_move["attempts"]), 0)
})
