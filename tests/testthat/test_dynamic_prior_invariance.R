# independent draws from joint priors detect incompatible Gibbs conditionals
skip_on_cran()

prior_invariance_z = function(difference) {
	colMeans(difference) / (apply(difference, 2, sd) / sqrt(nrow(difference)))
}

test_that("centered actor covariance updates preserve the constrained prior", {
	skip_on_cran()
	set.seed(39218)
	n = 3L
	eta = 7
	feature = function(S) c(S[1, 1], S[2, 2], S[1, 2], determinant(S, logarithm = TRUE)$modulus)
	difference = t(replicate(4000, {
		S = solve(rWishart(1, eta, diag(2) / eta)[, , 1])
		ab = matrix(rnorm(2 * n), n, 2) %*% chol(S)
		ab = sweep(ab, 2, colMeans(ab))
		next_S = rSab_fc(ab[, 1], ab[, 2], Sab0 = diag(2), eta0 = eta, center = TRUE)
		feature(next_S) - feature(S)
	}))
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})

test_that("centered AR actor hyperparameters preserve their joint prior", {
	skip_on_cran()
	set.seed(719361)
	n = 3L
	nt = 3L
	path = function(rho, variance) {
		x = matrix(0, n, nt)
		x[, 1] = rnorm(n, sd = sqrt(variance / (1 - rho^2)))
		for (t in 2:nt) x[, t] = rho * x[, t - 1L] + rnorm(n, sd = sqrt(variance))
		sweep(x, 2, colMeans(x))
	}
	prior = list(sigma_ab_shape = 3, sigma_ab_scale = .4, rho_ab_mean = .8, rho_ab_sd = .15)
	difference = t(replicate(5000, {
		rho = qnorm(runif(1, pnorm(-.99, .8, .15), pnorm(.99, .8, .15)), .8, .15)
		variance = 1 / rgamma(1, 3, .4)
		z = lame:::.lame_ab_hyper_update(
			path(rho, variance),
			path(rho, variance),
			rho,
			sqrt(variance),
			prior,
			center = TRUE
		)
		c(
			rho = z$rho_ab - rho,
			variance = z$sigma_ab^2 - variance,
			log_variance = log(z$sigma_ab^2) - log(variance)
		)
	}))
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})

test_that("AR coefficient and residual-scale proposals respect joint variance support", {
	expect_identical(lame:::.lame_ar_rho_support(.95, .5, 1, 2), .5)
	expect_identical(lame:::.lame_ar_rho_support(.3, .5, 1, 2), .3)
	expect_identical(lame:::.lame_ar_s2_support(.1, 2, rho_ab = .5, sigma_ab = 1, ab_cap = 1), 2)
	expect_identical(lame:::.lame_ar_s2_support(1.5, 2, rho_ab = .5, sigma_ab = 1, ab_cap = 1), 1.5)
	expect_error(lame:::.lame_ar_start_support(.95, 1, 2, resume = TRUE), "Checkpoint")
})

test_that("truncated innovation draws retain tail resolution without a boundary atom", {
	set.seed(91731)
	cap = .02
	draws = replicate(1000, lame:::.lame_rinvgamma_sd_trunc(2, 1, cap, cap / 2))
	expect_true(all(is.finite(draws) & draws > 0 & draws < cap))
	expect_length(unique(draws), length(draws))
	# conditional survival probabilities are Uniform(0,1), even when the
	# untruncated survival probability is below floating-point range
	log_tail = pgamma(1 / cap^2, 2, 1, lower.tail = FALSE, log.p = TRUE)
	uniform = exp(pgamma(1 / draws^2, 2, 1, lower.tail = FALSE, log.p = TRUE) - log_tail)
	expect_lt(abs(mean(uniform) - .5) / sqrt(1 / (12 * length(draws))), 6)
})

test_that("bounded coordinate conditionals preserve an exact fixed-data posterior", {
	set.seed(51991)
	N = 6000L
	bound = 1
	y = 1.5
	s2 = .5
	log_density = function(u, v) -.5 * (u^2 + v^2) - (y - u * v)^2 / (2 * s2)
	# exact independent starts by accept/reject from a uniform rectangle;
	# the density peaks at u = v = +/-1
	starts = matrix(numeric(), 0, 2)
	while (nrow(starts) < N) {
		u = runif(20000, -bound, bound)
		v = runif(20000, -bound, bound)
		keep = log(runif(length(u))) <= log_density(u, v) - log_density(1, 1)
		starts = rbind(starts, cbind(u[keep], v[keep]))
	}
	feature = function(u, v) c(u, v, u^2, v^2, u * v, u^2 * v^2, u^2 - v^2)
	difference = matrix(NA_real_, N, 7)
	for (i in seq_len(N)) {
		u = array(starts[i, 1], c(1, 1, 1))
		v = array(starts[i, 2], c(1, 1, 1))
		z = lame:::rUV_dynamic_fc_cpp(
			u,
			v,
			array(y, c(1, 1, 1)),
			0,
			1,
			s2,
			FALSE,
			FALSE,
			max_abs = bound
		)
		stopifnot(max(abs(z$U), abs(z$V)) <= bound)
		difference[i, ] = feature(z$U, z$V) - feature(u, v)
	}
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})

test_that("binary scale-move dimensions agree with constrained density and Jacobian", {
	set.seed(21871)
	n = 4L
	nt = 3L
	rank = 2L
	H = qr.Q(qr(contr.helmert(n)))
	level = matrix(rnorm((n - 1L) * 2), n - 1L, 2)
	dev = array(rnorm((n - 1L) * 2 * nt), c(n - 1L, 2, nt))
	uv = array(rnorm(n * rank * nt * 2), c(n, rank, nt, 2))
	beta = c(.3, -.7)
	P = diag(c(.4, .8))
	cscale = 1.2
	negative_log_prior = function(b, a, d, u) {
		sum(b * (P %*% b)) / 2 + sum(a^2) / 2 + sum(d^2) / 2 + sum(u^2) / 2
	}
	exact = negative_log_prior(beta, level, dev, uv) -
		negative_log_prior(cscale * beta, cscale * level, cscale * dev, sqrt(cscale) * uv) +
		(length(beta) + length(level) + length(dev) + length(uv) / 2) * log(cscale)
	actual = lame:::.lame_scale_prior_logratio(
		cscale,
		beta,
		P,
		sum((H %*% level)^2),
		length(level),
		sum(dev^2),
		length(dev),
		sum(uv^2),
		length(uv)
	)
	expect_equal(actual, exact, tolerance = 1e-12)
	# when deviations and their SD move together, their normalized Gaussian
	# density and coordinate Jacobian cancel at the constrained dimension
	variance = .7
	shape = 3
	rate = .4
	log_ig = function(v) dgamma(1 / v, shape, rate, log = TRUE) - 2 * log(v)
	log_gaussian = function(d, v) sum(dnorm(d, sd = sqrt(v), log = TRUE))
	exact = log_gaussian(cscale * dev, cscale^2 * variance) - log_gaussian(dev, variance) +
		log_ig(cscale^2 * variance) - log_ig(variance) + (length(dev) + 2) * log(cscale)
	expect_equal(lame:::.lame_scale_ig_logratio(variance, cscale^2 * variance, shape, rate), exact)
})

test_that("the static joint draw conditions actor means without changing the intercept prior", {
	set.seed(52591)
	n = 3L
	nt = 2L
	N = 4000L
	Y = array(2, c(n, n, nt))
	xflat = array(1, c(n * n, 1, nt))
	draws = t(replicate(N, {
		z = lame:::.lame_replicated_gaussian_draw(
			Y,
			xflat,
			diag(2),
			0,
			1,
			g = 1,
			center = TRUE
		)
		c(z$beta, sum(z$a), sum(z$b))
	}))
	# the centered actor columns are orthogonal to the all-ones intercept
	precision = n * n * nt + n * n
	target_mean = n * n * nt * 2 / precision
	target_variance = 1 / precision
	expect_lt(abs(mean(draws[, 1]) - target_mean) / sqrt(target_variance / N), 6)
	expect_lt(abs(var(draws[, 1]) - target_variance) / (target_variance * sqrt(2 / (N - 1))), 6)
	expect_lt(max(abs(draws[, 2:3])), 1e-10)
})

test_that("rectangular centered actor sweeps preserve the proper joint prior", {
	set.seed(91821)
	nr = 3L
	nc = 4L
	nt = 2L
	ve = .7
	g = 1.5
	va = .6
	vb = .9
	X = lapply(seq_len(nt), function(t) array(1, c(nr, nc, 1)))
	feature = function(beta, a, b) c(beta^2, mean(a^2), mean(b^2), beta * a[1], a[1] * b[1])
	difference = t(replicate(5000, {
		beta = rnorm(1, sd = sqrt(g * ve))
		a = rnorm(nr, sd = sqrt(va)); a = a - mean(a)
		b = rnorm(nc, sd = sqrt(vb)); b = b - mean(b)
		Z = array(rep(beta + outer(a, b, "+"), nt), c(nr, nc, nt)) +
			array(rnorm(nr * nc * nt, sd = sqrt(ve)), c(nr, nc, nt))
		z = lame:::rbeta_ab_bip_gibbs_cpp(
			Z,
			X,
			matrix(0, nr, nc),
			a,
			b,
			ve,
			g,
			va,
			vb,
			TRUE,
			TRUE,
			center = TRUE
		)
		stopifnot(abs(sum(z$a)) < 1e-10, abs(sum(z$b)) < 1e-10)
		feature(z$beta, z$a, z$b) - feature(beta, a, b)
	}))
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})

test_that("disabled dynamic actor margins remain outside the conditional model", {
	set.seed(47812)
	Z = array(rnorm(3 * 3 * 2), c(3, 3, 2))
	zero = matrix(0, 3, 2)
	for (role in c("sender", "receiver")) {
		z = lame:::sample_dynamic_ab_cpp(
			zero,
			zero,
			Z,
			array(0, dim(Z)),
			.4,
			.6,
			.7,
			FALSE,
			.5,
			center = TRUE,
			rvar = role == "sender",
			cvar = role == "receiver"
		)
		expect_equal(if (role == "sender") z$b else z$a, zero)
		expect_lt(max(abs(colSums(z$a)), abs(colSums(z$b))), 1e-10)
		# transposition and swapping active roles describe the same model
		set.seed(185)
		row = lame:::sample_dynamic_ab_cpp(
			zero,
			zero,
			Z,
			array(0, dim(Z)),
			.4,
			.6,
			.7,
			FALSE,
			.5,
			TRUE,
			TRUE,
			FALSE
		)
		set.seed(185)
		column = lame:::sample_dynamic_ab_cpp(
			zero,
			zero,
			aperm(Z, c(2, 1, 3)),
			array(0, dim(Z)),
			.4,
			.6,
			.7,
			FALSE,
			.5,
			TRUE,
			FALSE,
			TRUE
		)
		expect_equal(row$a, column$b)
	}
})

test_that("centered nomination contrasts preserve their prior and all residuals", {
	skip_on_cran()
	set.seed(31331)
	n = 5L
	nt = 2L
	S = matrix(c(.8, .25, .25, .6), 2)
	Y = array(NA_real_, c(n, n, nt))
	masks = lame:::.lame_nom_shift_masks(Y, FALSE)
	Z = array(rnorm(n * n * nt), c(n, n, nt))
	difference = t(replicate(3000, {
		ab = matrix(rnorm(2 * n), n, 2) %*% chol(S)
		ab = sweep(ab, 2, colMeans(ab))
		a = ab[, 1]; b = ab[, 2]
		z = lame:::.lame_nom_row_shift(
			Z,
			a,
			b,
			S,
			rep(n, nt),
			masks,
			FALSE,
			rep(TRUE, n),
			center = TRUE
		)
		stopifnot(
			abs(sum(z$a)) < 1e-10,
			max(abs(sweep(z$Z, 1, z$a) - sweep(Z, 1, a))) < 1e-10
		)
		c(mean(z$a^2) - mean(a^2), mean(z$a * b) - mean(a * b))
	}))
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})

test_that("dynamic bipartite conditionals use full paths in static actor and G updates", {
	set.seed(51629)
	Y = lapply(1:3, function(t) matrix(rnorm(56), 7, 8))
	actual_G = get("sample_G_bip_cpp", asNamespace("lame"))
	actual_ab = get("rbeta_ab_bip_gibbs_cpp", asNamespace("lame"))
	g_calls = ab_calls = 0L
	testthat::local_mocked_bindings(
		sample_G_bip_cpp = function(R, U, V, lambdaG, s2) {
			frame = parent.frame()
			expect_identical(U, get("U_cube", frame))
			expect_identical(V, get("V_cube", frame))
			g_calls <<- g_calls + 1L
			actual_G(R, U, V, lambdaG, s2)
		},
		rbeta_ab_bip_gibbs_cpp = function(
			Z,
			Xlist,
			UV_eff,
			a_current,
			b_current,
			s2,
			g_prior,
			va,
			vb,
			rvar,
			cvar,
			center = FALSE
		) {
			frame = parent.frame()
			want = get("Z", frame)
			U = get("U_cube", frame)
			V = get("V_cube", frame)
			G = get("G", frame)
			for (t in 1:3) want[, , t] = want[, , t] - U[, , t] %*% G %*% t(V[, , t])
			expect_equal(Z, want)
			expect_equal(UV_eff, matrix(0, 7, 8))
			ab_calls <<- ab_calls + 1L
			actual_ab(
				Z,
				Xlist,
				UV_eff,
				a_current,
				b_current,
				s2,
				g_prior,
				va,
				vb,
				rvar,
				cvar,
				center
			)
		},
		.package = "lame"
	)
	f = suppressMessages(lame(
		Y,
		family = "normal",
		mode = "bipartite",
		R = 2,
		dynamic_uv = TRUE,
		dynamic_ab = FALSE,
		burn = 0,
		nscan = 4,
		odens = 1,
		seed = 87731,
		verbose = FALSE,
		gof = FALSE
	))
	expect_equal(g_calls, 4L)
	expect_equal(ab_calls, 4L)
	expect_true(all(is.finite(f$VC)))
})

test_that("actor-to-factor composition refreshes collapsed artificial diagonals", {
	skip_on_cran()
	set.seed(73218)
	n = 3L
	nt = 2L
	ar = .4
	ve = .3
	path = function(sig, center) {
		x = matrix(0, n, nt)
		x[, 1] = rnorm(n, sd = sig / sqrt(1 - ar^2))
		x[, 2] = ar * x[, 1] + rnorm(n, sd = sig)
		if (center) x = sweep(x, 2, colMeans(x))
		x
	}
	uv_mean = function(u, v) array(
		unlist(lapply(
			1:nt,
			function(t)
				outer(u[, , t], v[, , t])
		)),
		c(n, n, nt)
	)
	ab_mean = function(a, b) array(
		unlist(lapply(
			1:nt,
			function(t)
				outer(a[, t], b[, t], "+")
		)),
		c(n, n, nt)
	)
	feature = function(a, b, u, v) c(
		mean(u^2),
		mean(v^2),
		mean(uv_mean(u, v)^2),
		mean(ab_mean(a, b) * uv_mean(u, v))
	)
	difference = t(replicate(4000, {
		a = path(.7, TRUE); b = path(.7, TRUE)
		u = array(path(.6, FALSE), c(n, 1, nt))
		v = array(path(.6, FALSE), c(n, 1, nt))
		uv = uv_mean(u, v)
		mean = ab_mean(a, b) + uv
		Z = mean + array(rnorm(n * n * nt, sd = sqrt(ve)), dim(mean))
		before = feature(a, b, u, v)
		ab = lame:::sample_dynamic_ab_cpp(a, b, Z, uv, ar, .7, ve, FALSE, 0, TRUE)
		offset = ab_mean(ab$a, ab$b)
		Z = lame:::.lame_refresh_diagonal(Z, lame:::.lame_diagonal(offset + uv), 0, ve)
		z = lame:::rUV_dynamic_fc_cpp(u, v, Z - offset, ar, .6, ve, FALSE, FALSE)
		feature(ab$a, ab$b, z$U, z$V) - before
	}))
	expect_lt(max(abs(prior_invariance_z(difference))), 6)
})


test_that("static actors receive fresh diagonals after collapsed coefficient updates", {
	set.seed(7843)
	Y = lapply(1:3, function(t) {
		y = matrix(rnorm(36), 6, 6)
		diag(y) = NA_real_
		y
	})
	X = lapply(
		Y,
		function(y) array(
			rnorm(36),
			c(6, 6, 1),
			dimnames = list(NULL, NULL, "x")
		)
	)
	state = new.env(parent = emptyenv())
	state$coefficient = FALSE
	state$refreshed = FALSE
	state$calls = 0L
	actual_beta = get("sample_beta_dynamic_cpp", asNamespace("lame"))
	actual_refresh = get(".lame_refresh_diagonal", asNamespace("lame"))
	actual_ab = get(".lame_static_actor_draw", asNamespace("lame"))
	testthat::local_mocked_bindings(
		sample_beta_dynamic_cpp = function(...) {
			state$coefficient = TRUE
			state$refreshed = FALSE
			actual_beta(...)
		},
		.lame_refresh_diagonal = function(...) {
			if (state$coefficient) state$refreshed = TRUE
			actual_refresh(...)
		},
		.lame_static_actor_draw = function(...) {
			expect_true(state$coefficient && state$refreshed)
			state$calls = state$calls + 1L
			state$coefficient = FALSE
			actual_ab(...)
		},
		.package = "lame"
	)
	f = suppressMessages(lame(
		Y,
		Xdyad = X,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		dynamic_ab = FALSE,
		dcor = FALSE,
		burn = 0,
		nscan = 4,
		odens = 1,
		seed = 917,
		verbose = FALSE,
		gof = FALSE
	))
	expect_equal(state$calls, 4L)
	expect_true(all(is.finite(f$BETA)))
})


test_that("actor conditionals retain residual correlation near either boundary", {
	set.seed(5431)
	n = 3L
	Y = matrix(rnorm(n * n), n)
	ij = which(row(Y) != col(Y))
	xy = arrayInd(ij, c(n, n))
	X = diag(n)[xy[, 1], ]
	N = 1800L
	for (rho in c(-.997, .997)) {
		noise = diag(.7, length(ij))
		for (i in seq_along(ij)) for (j in seq_along(ij)) {
			if (all(xy[i, ] == rev(xy[j, ]))) noise[i, j] = .7 * rho
		}
		W = solve(noise)
		V = solve(diag(1 / .8, n) + crossprod(X, W %*% X))
		mu = drop(V %*% crossprod(X, W %*% Y[ij]))
		u = rowSums(V)
		mu = mu - u * sum(mu) / sum(u)
		V = V - outer(u, u) / sum(u)
		for (kind in c("path", "level")) {
			x = t(replicate(N, if (kind == "path") {
				lame:::sample_dynamic_ab_cpp(
					matrix(0, n, 1),
					matrix(0, n, 1),
					array(Y, c(n, n, 1)),
					array(0, c(n, n, 1)),
					0,
					sqrt(.8),
					.7,
					FALSE,
					rho,
					TRUE,
					TRUE,
					FALSE
				)$a[, 1]
			} else {
				lame:::.lame_ab_level_draw_lik(
					array(Y, c(n, n, 1)),
					matrix(0, n, 1),
					matrix(0, n, 1),
					.8,
					.7,
					1,
					rho,
					TRUE,
					center = TRUE
				)
			}))
			expect_lt(max(abs(colMeans(x) - mu) / sqrt(diag(V) / N)), 6)
			expect_lt(max(abs(cov(x) - V) / sqrt((V^2 + outer(diag(V), diag(V))) / (N - 1))), 6)
		}
	}
})
