test_that('actor-slope prior controls are validated and default draws are unchanged', {
	skip_on_cran()
	set.seed(974612);n = 4L;nt = 3L;X = lapply(1:nt, function(t) array(rnorm(n * n), c(n, n, 1)));Y = lapply(
		1:nt,
		function(t) {x = matrix(rnorm(n * n), n);diag(x) = NA;x}
	)
	args = list(
		Y = Y,
		Xdyad = X,
		family = 'normal',
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		dynamic_beta_per_actor = 'row',
		burn = 40L,
		nscan = 160L,
		odens = 4L,
		seed = 91746,
		verbose = FALSE,
		gof = FALSE
	)
	# default innovation scale is the variance of a normal response (R/normal_scale.R)
	defaults = list(
		sigma_actor_shape = 2,
		sigma_actor_scale = lame:::.lame_unit_var(Y, 'normal'),
		rho_actor_shape1 = 2,
		rho_actor_shape2 = 2
	)
	a = do.call(lame, args);b = do.call(lame, c(args, list(prior = defaults)))
	for (nm in c(
		'BETA', 'VC', 'THETA_ACTOR', 'RHO_ACTOR', 'SIGMA_ACTOR', 'EZ', 'YPM'
	)) expect_identical(a[[nm]], b[[nm]])
	changed = list(sigma_actor_shape = 3, sigma_actor_scale = .1, rho_actor_shape1 = 1.5, rho_actor_shape2 = 2.5)
	c = do.call(lame, c(args, list(prior = changed)))
	expect_equal(c$per_actor_prior$innovation_variance$shape, 3)
	expect_equal(c$per_actor_prior$innovation_variance$scale, .1)
	expect_equal(c$per_actor_prior$rho$shape1, 1.5)
	expect_equal(c$per_actor_prior$rho$shape2, 2.5)
	expect_false(identical(c$SIGMA_ACTOR, a$SIGMA_ACTOR))
	for (nm in names(defaults)) for (value in list(0, -1, Inf, NA_real_, c(1, 2), 'x')) {
		p = setNames(list(value), nm)
		expect_error(do.call(lame, c(args, list(prior = p))), 'finite positive number')
	}
})
