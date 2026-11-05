.ame_health_fixture = function() {
	matrix(c(NA, 1, 0, 1, 0, NA, 1, 0, 1, 0, NA, 0, 0, 1, 1, NA), 4L)
}
.ame_health_capture = function(expr) {
	warnings = character()
	fit = withCallingHandlers(
		expr,
		warning = function(w) {warnings <<- c(warnings, conditionMessage(w));invokeRestart('muffleWarning')}
	)
	list(fit = fit, warnings = warnings)
}
.ame_health_fit = function() {
	ame(
		Y = .ame_health_fixture(),
		family = 'binary',
		R = 0L,
		dcor = TRUE,
		burn = 10L,
		nscan = 50L,
		odens = 5L,
		verbose = FALSE,
		gof = FALSE,
		seed = 47337341L
	)
}

test_that('AME reports full-sweep clean numerical health', {
	got = .ame_health_capture(.ame_health_fit())
	expect_true(all(got$fit$sampler_health$numerical_failures == 0))
	expect_true(all(got$fit$sampler_health$prior_bound_rejections == 0))
	expect_equal(got$fit$sampler_health_iterations, 60L)
	expect_false(any(grepl('Numerical sampling failures', got$warnings, fixed = TRUE)))
})

test_that('caught AME collapsed reciprocity errors abort the fit with their counts', {
	local_mocked_bindings(.ame_binary_collapsed_rho_mh = function(...) stop('injected MH failure'), .package = 'lame')
	e = expect_error(.ame_health_fit(), class = 'lame_sampler_failure')
	h = e$health
	expect_equal(h$numerical_failures[h$step == 'rho_collapsed'], 2L)
	expect_match(gsub('[[:space:]]+', ' ', conditionMessage(e)), 'of 60 sweeps')
})

test_that('caught AME utility refresh errors abort the fit with their counts', {
	local_mocked_bindings(
		.ame_binary_collapsed_rho_mh = function(rho, eta, cx, tries = 4L) list(
			rho = if (rho < .9) rho + .01 else rho - .01,
			accepts = 1L,
			attempts = tries
		),
		.lame_binary_dynamic_collapsed_rho_z = function(...) stop('injected Z failure'),
		.package = 'lame'
	)
	e = expect_error(.ame_health_fit(), class = 'lame_sampler_failure')
	h = e$health
	expect_equal(h$numerical_failures[h$step == 'rho_collapsed'], 2L)
	expect_true(all(h$prior_bound_rejections == 0))
})
