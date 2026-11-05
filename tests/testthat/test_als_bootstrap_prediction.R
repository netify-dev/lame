skip_on_cran()

test_that("bootstrap response generation retains means above one million", {
	draw = get('.ae_boot_draw_response', asNamespace('lame'))
	f = list(family = 'poisson', mode = 'bipartite', link = 'log')
	set.seed(782023L);mu = 2e6;ez = matrix(log(mu), 150, 100)
	y = draw(f, ez, 0, 0)
	expect_true(all(is.finite(y)))
	expect_lt(abs(mean(y) - mu) / sqrt(mu / length(y)), 5)
	# the latent normal variance changes both the mean and the count variance
	ve = .2;mu = 10;ez = matrix(log(mu) - ve / 2, 500, 300)
	y = draw(f, ez, 0, ve)
	variance = mu + mu^2 * expm1(ve)
	expect_lt(abs(mean(y) - mu) / sqrt(variance / length(y)), 5)
})

test_that("bootstrap binary thresholds retain extreme link probabilities", {
	draw = get('.ae_boot_draw_response', asNamespace('lame'))
	local_mocked_bindings(simZ = function(EZ, rho, s2) EZ - 6, .package = 'lame')
	for (link in c('probit', 'logit')) {
		f = list(family = 'binary', mode = 'unipartite', link = link)
		eta = if (link == 'probit') 8 else 40
		expect_equal(draw(f, matrix(eta, 3, 3), 0, 1), matrix(1, 3, 3))
	}
})

test_that("all bootstrap draw routes preserve square bipartite observed cells", {
	set.seed(782024L);y = lapply(1:2, function(t) matrix(rnorm(144), 12, 12))
	y[[1]][2, 3] = NA_real_
	f = lame_als(y, R = 0, mode = 'bipartite', family = 'normal', verbose = FALSE)
	captured = list()
	local_mocked_bindings(ame_als_refit = function(object, Y_new = NULL, Z_new = NULL, ...) {
		captured[[length(captured) + 1L]] <<- if (is.null(Z_new)) Y_new else Z_new
		object
	}, .package = 'lame')
	for (route in c('normal', 'binary', 'poisson', 'transform')) {
		g = f;g$family = if (route == 'transform') 'binary' else route
		g$non_normal_method = if (route == 'transform') 'transform' else 'irls'
		g$link = if (route == 'poisson') 'log' else 'probit'
		# the draw routes are the subject here: a corrected binary fit is bootstrapped
		g$bias_correct = list(method = 'jackknife', applied = TRUE)
		captured = list()
		b = ame_als_bootstrap(g, R = 12, type = 'parametric', seed = 782025L, verbose = FALSE)
		expect_equal(b$n_valid, 12)
		expect_length(captured, 12)
		for (z in captured) expect_identical(is.na(z), is.na(g$Y), ignore_attr = TRUE)
	}
})

test_that("a parametric bootstrap refuses unidentified binary reciprocity", {
	set.seed(782026L);y = matrix(rnorm(144), 12, 12);diag(y) = NA_real_
	f = ame_als(y, R = 0, family = 'normal', verbose = FALSE)
	f$family = 'binary';f$rho_estimation = list(identified = FALSE);f$VC['rho'] = NA_real_
	expect_error(ame_als_bootstrap(f, R = 12, verbose = FALSE), class = 'lame_als_prediction')
})
