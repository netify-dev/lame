test_that('Hamiltonian acceptance is stable at finite extreme momentum scales', {
	set.seed(882614)
	for (i in 1:100) {
		old = rnorm(10);new = rnorm(10);mass = exp(rnorm(10));a = rnorm(1);b = rnorm(1)
		expect_identical(.lame_hmc_log_acceptance(a, b, new, old, mass), a - b - .5 * sum((new^2 - old^2) / mass))
	}
	expect_equal(.lame_hmc_log_acceptance(-3, -2, 2e154, 3e154, 1e308), 1.5, tolerance = 1e-14)
	expect_equal(.lame_hmc_log_acceptance(-3, -2, 2e154, 3e154, 1e300), 249999999, tolerance = 1e-14)
	expect_false(is.finite(.lame_hmc_log_acceptance(0, 0, 1e200, 0, 1)))
})

test_that('Hamiltonian steps distinguish recoverable and unresolved energy overflow', {
	for (prefix in c('.lame_binary_hmc', '.lame_binary_dynamic_ar_hmc')) {
		ns = asNamespace('lame');e = new.env(parent = ns);f = get(paste0(prefix, '_step'), ns);environment(f) = e
		assign(paste0(prefix, '_evaluate'), function(...) list(logp = 0, gradient = 0), e)
		args = list(q = 0, cx = NULL, h = NULL, step = .1, L = 1)
		if ('offset' %in% names(formals(f))) args$offset = NULL
		assign(paste0(prefix, '_mass'), function(...) 1e308, e)
		assign(
			paste0(prefix, '_integrate'),
			function(...) list(q = 0, p = 1e155, evaluation = list(logp = 0, gradient = 0), failed = FALSE),
			e
		)
		set.seed(84152);old_normal = rnorm(1);set.seed(84152);z = do.call(f, args)
		expect_false(z$nonfinite);expect_equal(log(z$prob), .5 * (old_normal^2 - 100), tolerance = 1e-12)
		assign(paste0(prefix, '_mass'), function(...) 1, e)
		assign(
			paste0(prefix, '_integrate'),
			function(...) list(q = 0, p = 1e200, evaluation = list(logp = 0, gradient = 0), failed = FALSE),
			e
		)
		z = do.call(f, args)
		expect_true(z$nonfinite);expect_false(z$accepted);expect_identical(z$q, 0);expect_equal(z$prob, 0)
	}
})
