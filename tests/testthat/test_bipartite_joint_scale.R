test_that("bipartite scale acceptance includes the complete density and Jacobian", {
	set.seed(518207)
	rho = .3
	for (rank in 1:2) for (sigma in c(1e-8, .5, 2)) {
		U = array(rnorm(3 * rank * 3, sd = sigma), c(3L, rank, 3L))
		V = array(rnorm(4 * rank * 3, sd = sigma), c(4L, rank, 3L))
		G = matrix(rnorm(rank^2), rank)
		shape = 2
		scale = sigma^2
		density = function(U, V, G, sigma) {
			path = function(z) {
				ans = sum(dnorm(z[, , 1], sd = sigma / sqrt(1 - rho^2), log = TRUE))
				for (t in 2:3) ans = ans + sum(dnorm(z[, , t] - rho * z[, , t - 1], sd = sigma, log = TRUE))
				ans
			}
			path(U) + path(V) + sum(dnorm(G, log = TRUE)) - (shape + 1) * log(sigma^2) - scale / sigma^2
		}
		for (delta in c(-.7, -1e-7, 0, .3, .8)) {
			multiplier = exp(delta)
			direct = density(U * multiplier, V * multiplier, G / multiplier^2, sigma * multiplier) -
				density(U, V, G, sigma) + (length(U) + length(V) - 2 * length(G) + 2) * delta
			got = .lame_bipartite_scale_log_ratio(delta, sigma, G, shape, scale)
			# subtracting complete log densities loses relative precision near zero
			expect_lt(abs(got - direct), 1e-10)
			expect_equal(
				got,
				-.lame_bipartite_scale_log_ratio(
					-delta,
					sigma * multiplier,
					G / multiplier^2,
					shape,
					scale
				),
				tolerance = 1e-10
			)
		}
	}
	expect_equal(.lame_bipartite_scale_log_ratio(0, .5, matrix(0, 1, 1), 2, 1), 0)
	expect_equal(
		.lame_bipartite_scale_log_ratio(.1, 1e-140, matrix(1e-200, 1, 1), 2, 1e-280),
		-6 * .1 - expm1(-.2),
		tolerance = 1e-12
	)
	expect_identical(.lame_bipartite_scale_log_ratio(.1, 1, matrix(1e200, 1, 1), 2, 1), Inf)
	expect_identical(.lame_bipartite_scale_log_ratio(-.1, 1, matrix(1e200, 1, 1), 2, 1), -Inf)
})

test_that("joint scale moves preserve products and every physical prior bound", {
	set.seed(946022)
	for (rank in 1:2) {
		z = list(
			U = array(rnorm(3 * rank * 3, sd = .4), c(3L, rank, 3L)),
			V = array(rnorm(4 * rank * 3, sd = .4), c(4L, rank, 3L)),
			G = matrix(.4, rank, rank),
			sigma = .4
		)
		products = function(z) vapply(
			1:3,
			function(t) as.vector(
				matrix(z$U[, , t], 3, rank) %*% z$G %*% t(matrix(z$V[, , t], 4, rank))
			),
			numeric(12)
		)
		original = products(z)
		accepts = bounds = numerical = 0L
		discrepancy = 0
		for (k in 1:500) {
			z = do.call(
				.lame_bipartite_scale_step,
				c(
					z[c('U', 'V', 'G', 'sigma')],
					list(rho = .3, shape = 2, scale = 1, max_abs = 2, max_G = 1.5, sigma_max = 1.2, stat_var_cap = 1.4)
				)
			)
			accepts = accepts + z$accepts
			bounds = bounds + z$bound_rejections
			numerical = numerical + z$numerical_rejections
			discrepancy = max(discrepancy, max(abs(products(z) - original)))
			stopifnot(
				max(abs(z$U), abs(z$V)) <= 2,
				max(abs(z$G)) <= 1.5,
				z$sigma <= 1.2,
				z$sigma^2 / (1 - .3^2) <= 1.4
			)
		}
		expect_lt(discrepancy, 1e-11)
		expect_gt(accepts, 20L)
		expect_gt(bounds, 0L)
		expect_identical(numerical, 0L)
	}
})

test_that("Gaussian bipartite fits expose scale health and retain a public opt-out", {
	skip_on_cran()
	set.seed(685341)
	rn = paste0('r', 1:5);cn = paste0('c', 1:6)
	Y = lapply(1:3, function(t) matrix(rnorm(30), 5, 6, dimnames = list(rn, cn)))
	args = list(
		Y = Y,
		mode = 'bipartite',
		family = 'normal',
		R = 1,
		dynamic_uv = TRUE,
		rvar = FALSE,
		cvar = FALSE,
		burn = 500,
		nscan = 2000,
		odens = 1,
		seed = 44351,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		posterior_opts = list(save_UV = TRUE, thin_UV = 1)
	)
	f = do.call(lame, args)
	expect_identical(unname(f$bipartite_scale_move['attempts']), 2250)
	expect_gt(unname(f$bipartite_scale_move['accepts']), 100)
	expect_identical(unname(f$bipartite_scale_move['numerical_rejections']), 0)
	expect_true(all(f$sampler_health$numerical_failures == 0))
	expect_equal(nrow(f$VC), 2000)
	args$prior = list(bipartite_scale_move = FALSE)
	args$nscan = 20;args$burn = 20
	off = do.call(lame, args)
	expect_null(off$bipartite_scale_move)
	args$prior$bipartite_scale_move = 1
	expect_error(do.call(lame, args), 'bipartite_scale_move must be TRUE or FALSE', fixed = TRUE)
	pooled = .lame_pool_sampler_diagnostics(list(), list(f, f))
	expect_null(pooled$bipartite_scale_move)
	expect_equal(pooled$sampler_diagnostics_by_chain[[2]]$bipartite_scale_move, f$bipartite_scale_move)
	expect_error(
		.lame_report_bipartite_scale_health(list(
			bipartite_scale_move = c(attempts = 10, numerical_rejections = 1)
		)),
		class = 'lame_sampler_failure'
	)
	expect_no_error(.lame_report_bipartite_scale_health(list(
		bipartite_scale_move = c(attempts = 10, numerical_rejections = 0, bound_rejections = 20)
	)))
})
