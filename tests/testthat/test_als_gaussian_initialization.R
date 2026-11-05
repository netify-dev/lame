skip_on_cran()

test_that('dynamic fits preserve their initializer and certify their final objective', {
	set.seed(11L);n = 8L;Y = lapply(
		1:4,
		function(tt) {a = seq(-1, 1, length.out = n) * sin(tt / 2);b = rev(a) * cos(tt / 3);z = .25 + outer(a, b, '+') + matrix(rnorm(n * n, sd = .2), n);diag(z) = NA_real_;z}
	)
	local_mocked_bindings(
		.ae_gaussian_certify = function(...) stop('reported Gaussian certificate called'),
		.package = 'lame'
	)
	fit = lame(Y, family = 'normal', R = 0L, method = 'als', dynamic_ab = TRUE, als_tol = 1e-6, verbose = FALSE)
	expect_s3_class(fit, 'lame_dynamic_als')
	expect_true(fit$converged)
	expect_true(fit$diagnostics$stationarity$pass)
	# the static fixed-effects request uses the Gaussian certificate; the empirical-Bayes default has its own
	expect_error(
		lame_als(Y, family = 'normal', R = 0L, bias_correct = 'none', verbose = FALSE),
		'reported Gaussian certificate called'
	)
})
