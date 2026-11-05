skip_on_cran()

.capture_sampler_health_messages = function(expr) {
	warnings = messages = character()
	value = withCallingHandlers(
		expr,
		warning = function(e) { warnings <<- c(warnings, conditionMessage(e)); invokeRestart('muffleWarning') },
		message = function(e) { messages <<- c(messages, conditionMessage(e)); invokeRestart('muffleMessage') }
	)
	list(value = value, warnings = warnings, messages = messages)
}

test_that('a few numerical errors cannot be hidden by many valid bound rejections', {
	h = .lame_sampler_health(list(UV = 101L, G = 200L, beta = 1L), 100L, 199L)
	expect_identical(h$numerical_failures, c(1, 1, 1))
	expect_identical(h$prior_bound_rejections, c(100, 199, 0))
	e = expect_error(.lame_report_sampler_health(h, 100000L), class = 'lame_sampler_failure')
	expect_true(
		all(vapply(c('UV', 'G', 'beta'), function(nm) grepl(nm, conditionMessage(e), fixed = TRUE), logical(1)))
	)
	expect_identical(e$health, h)
	expect_error(.lame_sampler_health(list(UV = 1L), 2L), 'Invalid sampler')
})

test_that('valid prior-bound rejections and zero counts are not failures and say nothing while fitting', {
	h = .lame_sampler_health(list(UV = 3L, G = 7L, beta = 0L), 3L, 7L)
	a = .capture_sampler_health_messages(.lame_report_sampler_health(h, 10L))
	expect_length(a$warnings, 0L); expect_length(a$messages, 0L); expect_true(all(h$numerical_failures == 0))
	z = .capture_sampler_health_messages(.lame_report_sampler_health(
		.lame_sampler_health(list(UV = 0L, G = 0L, beta = 0L)),
		10L
	))
	expect_length(z$warnings, 0L); expect_length(z$messages, 0L)
	# the summary states the rejections as a property when they are on one percent or more of the sweeps
	line = lame:::.lame_bound_share_line(list(sampler_health = h, sampler_health_iterations = 10L))
	expect_match(line, 'prior\\$uv_max_abs'); expect_match(line, 'prior\\$G_max'); expect_match(line, '30%')
	expect_match(line, 'Metropolis rejection')
	expect_null(lame:::.lame_bound_share_line(list(sampler_health = h, sampler_health_iterations = 10000L)))
	expect_null(lame:::.lame_bound_share_line(list()))
})

.health_fit_data = function() {
	set.seed(316924)
	lapply(
		1:2,
		function(t) matrix(
			rnorm(30),
			6L,
			5L,
			dimnames = list(paste0('r', 1:6), paste0('c', 1:5))
		)
	)
}

test_that('a caught native sampler error aborts the fit and carries its counts', {
	local_mocked_bindings(rUV_dynamic_bip_fc_cpp = function(...) stop('injected native failure'), .package = 'lame')
	e = expect_error(
		lame(
			Y = .health_fit_data(),
			mode = 'bipartite',
			R = 1L,
			family = 'normal',
			dynamic_uv = TRUE,
			burn = 0L,
			nscan = 10L,
			odens = 1L,
			verbose = FALSE,
			gof = FALSE,
			seed = 29164L
		),
		class = 'lame_sampler_failure'
	)
	h = e$health
	expect_gt(h$numerical_failures[h$step == 'UV'], 0)
	expect_equal(h$prior_bound_rejections[h$step == 'UV'], 0)
})

test_that('actual bounded Gaussian proposals have exact separate fit counts', {
	native = get('rUV_dynamic_bip_fc_cpp', asNamespace('lame')); visits = 0L
	local_mocked_bindings(rUV_dynamic_bip_fc_cpp = function(...) {
		z = native(...)
		visits <<- visits + as.integer(z$bound_rejections > 0)
		z
	}, .package = 'lame')
	a = .capture_sampler_health_messages(lame(
		Y = .health_fit_data(),
		mode = 'bipartite',
		R = 1L,
		family = 'normal',
		dynamic_uv = TRUE,
		burn = 0L,
		nscan = 20L,
		odens = 1L,
		prior = list(uv_max_abs = .05),
		verbose = FALSE,
		gof = FALSE,
		seed = 29164L
	))
	f = a$value; h = f$sampler_health
	expect_gt(visits, 0L)
	expect_equal(h$prior_bound_rejections[h$step == 'UV'], visits)
	expect_equal(f$tryErrorChecks$UV, visits)
	expect_equal(h$numerical_failures[h$step == 'UV'], 0)
	expect_false(any(grepl('Numerical sampling failures', a$warnings, fixed = TRUE)))
})

test_that('pooled fits retain a failure from a later chain and sum bound counts', {
	make = function(uv, g, uv_bound, g_bound, iters) list(
		BETA = matrix(1:8, ncol = 1),
		VC = matrix(1:8, ncol = 1),
		tryErrorChecks = list(UV = uv, G = g),
		sampler_health = .lame_sampler_health(list(UV = uv, G = g), uv_bound, g_bound),
		sampler_health_iterations = iters,
		marginal_move = c(attempts = 10, accepts = 5, step_sd = .1)
	)
	a = make(2L, 0L, 2L, 0L, 20L); b = make(3L, 4L, 2L, 4L, 40L)
	f = combine_ame_chains(list(a, b), diagnostics = FALSE)
	expect_equal(f$sampler_health$numerical_failures, c(1, 0))
	expect_equal(f$sampler_health$prior_bound_rejections, c(4, 4))
	expect_equal(f$tryErrorChecks, list(UV = 5, G = 4))
	expect_equal(f$refusal_share, c(UV = 5 / 60, G = 4 / 60))
	expect_equal(f$sampler_health_iterations, 60)
	expect_null(f$marginal_move)
	expect_identical(f$sampler_diagnostics_by_chain[[2L]]$marginal_move, b$marginal_move)
	b$sampler_health = NULL
	legacy = combine_ame_chains(list(a, b), diagnostics = FALSE)
	expect_null(legacy$sampler_health)
	expect_identical(legacy$sampler_diagnostics_by_chain[[1L]]$sampler_health, a$sampler_health)
	expect_equal(legacy$tryErrorChecks, list(UV = 5, G = 4))
})


test_that('joint numerical rejections abort the fit even when conditional health is clean', {
	f = list(
		marginal_move = c(numerical_rejections = 1),
		dynamic_marginal_move = c(numerical_rejections = 2, bound_rejections = 90)
	)
	e = expect_error(.lame_report_joint_sampler_health(f), class = 'lame_sampler_failure')
	expect_match(conditionMessage(e), 'marginal_move')
	expect_match(conditionMessage(e), 'dynamic_marginal_move')
	f$marginal_move['numerical_rejections'] = 0
	f$dynamic_marginal_move['numerical_rejections'] = 0
	expect_no_error(.lame_report_joint_sampler_health(f))
})

test_that('an actual joint trajectory numerical rejection aborts the fit', {
	set.seed(361853); n = 6L; nt = 2L
	X = lapply(
		seq_len(nt),
		function(t) array(
			rnorm(n * n),
			c(n, n, 1),
			dimnames = list(as.character(1:n), as.character(1:n), 'x')
		)
	)
	Y = lapply(seq_len(nt), function(t) {
		z = matrix(rbinom(n * n, 1, .2), n, n, dimnames = list(as.character(1:n), as.character(1:n)))
		diag(z) = NA; z[1, 2] = 1; z[2, 1] = 0; z
	})
	local_mocked_bindings(.lame_binary_hmc_step = function(q, ...) {
		list(q = q, accepted = FALSE, prob = 0, nonfinite = TRUE)
	}, .package = 'lame')
	e = expect_error(
		suppressMessages(lame(
			Y = Y,
			Xdyad = X,
			R = 1L,
			family = 'binary',
			dcor = FALSE,
			dynamic_ab = TRUE,
			dynamic_beta = 'x_dyad',
			dynamic_beta_kind = 'rw1',
			burn = 0L,
			nscan = 20L,
			odens = 5L,
			seed = 71833L,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		class = 'lame_sampler_failure'
	)
	expect_match(conditionMessage(e), 'Joint Hamiltonian trajectories', fixed = TRUE)
	expect_equal(unname(e$failures['marginal_move']), 1)
})

.separation_fit = function(eta, family = 'binary', R = 2L) {
	n = 6L; nt = 2L; nm = as.character(1:n)
	EZ = lapply(seq_len(nt), function(t) {
		m = matrix(eta[t], n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	Y = array(1, c(n, n, nt), dimnames = list(nm, nm, NULL))
	for (t in seq_len(nt)) diag(Y[, , t]) = NA
	list(family = family, R = R, EZ = EZ, YPM = EZ, Y = Y)
}

test_that('a binary fit with a saturated linear predictor warns of separation', {
	a = .capture_sampler_health_messages(.lame_separation_check(.separation_fit(c(0.5, 12))))
	expect_equal(a$value, 0.5)
	expect_length(a$warnings, 1L)
	expect_match(a$warnings, 'separated')
	expect_match(a$warnings, 'smaller', fixed = TRUE)
	r0 = .capture_sampler_health_messages(.lame_separation_check(.separation_fit(c(-12, -12), R = 0L)))
	expect_match(r0$warnings, 'all present or all absent', fixed = TRUE)
})

test_that('ordinary binary and non-binary fits do not warn of separation', {
	a = .capture_sampler_health_messages(.lame_separation_check(.separation_fit(c(-3, 4))))
	expect_equal(a$value, 0); expect_length(a$warnings, 0L)
	b = .capture_sampler_health_messages(.lame_separation_check(.separation_fit(c(20, 20), family = 'normal')))
	expect_length(b$warnings, 0L)
	# unobserved cells do not count toward the share
	f = .separation_fit(c(0.5, 12)); f$Y[, , 2] = NA
	expect_length(.capture_sampler_health_messages(.lame_separation_check(f))$warnings, 0L)
})

test_that('separation share matches actors by name when the actor set changes', {
	nm = as.character(1:6)
	Y = array(1, c(6L, 6L, 2L), dimnames = list(nm, nm, NULL))
	for (t in 1:2) diag(Y[, , t]) = NA
	Y[6, , 2] = NA; Y[, 6, 2] = NA # actor 6 absent in period 2
	Y[1, 2, 2] = NA # one missing cell among present actors
	e1 = matrix(0, 6L, 6L, dimnames = list(nm, nm)); diag(e1) = NA
	e2 = matrix(0, 5L, 5L, dimnames = list(nm[1:5], nm[1:5])); diag(e2) = NA
	e2['1', '2'] = 20 # saturated only at the missing cell
	f = list(family = 'binary', R = 2L, EZ = list(e1, e2), YPM = list(e1, e2), Y = Y)
	a = .capture_sampler_health_messages(.lame_separation_check(f))
	expect_equal(a$value, 0)
	expect_length(a$warnings, 0L)
})
