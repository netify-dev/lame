skip_on_cran()

test_that("count reciprocity inverts the exact cross-dyad covariance", {
	n = 8L;y = array(NA_real_, c(2, 2, n));y[1, 2, ] = c(0, 0, 0, 0, 2, 2, 2, 2);y[2, 1, ] = c(0, 0, 0, 2, 0, 2, 2, 2)
	eta = array(0, dim(y));r = lame:::.ae_poisson_pair_rho(y, eta, log(2))
	expect_equal(r$rho, log(1.5) / log(2), tolerance = 1e-12);expect_true(r$identified);expect_equal(r$n_pairs, 8)
	for (scale in c(1e12, 1e200)) expect_equal(
		lame:::.ae_poisson_pair_rho(y * scale, eta + log(scale), log(2))$rho,
		r$rho,
		tolerance = 1e-11
	)
	yp = y;yp[1, 2, ] = yp[2, 1, ] = 0
	expect_equal(lame:::.ae_poisson_pair_rho(yp, eta, log(2))$rho, 1)
	yp[1, 2, ] = 2
	expect_equal(lame:::.ae_poisson_pair_rho(yp, eta, log(2))$rho, -1)
	yp[] = NA_real_;expect_false(lame:::.ae_poisson_pair_rho(yp, eta, .4)$identified)
	expect_false(lame:::.ae_poisson_pair_rho(y, eta, 0)$identified)
	expect_error(lame:::.ae_poisson_pair_rho(y, eta, 1000), class = 'lame_als_numerical')
})

test_that("count prediction distinguishes zero variance from missing dependence", {
	set.seed(76019);y = matrix(rnorm(36), 6);diag(y) = NA_real_;f = ame_als(
		y,
		R = 0,
		family = 'normal',
		verbose = FALSE
	)
	f$family = 'poisson';f$link = 'log';f$EZ = lapply(
		f$EZ,
		function(z) {z[is.finite(z)] = log(2);z}
	);f$VC['ve'] = 0;f$VC['rho'] = NA_real_;f$rho_estimation = list(identified = FALSE)
	expect_no_error(simulate(f, nsim = 2, seed = 1))
	f$VC['ve'] = 1
	expect_error(simulate(f, nsim = 2, seed = 1), class = 'lame_als_prediction')
	expect_error(ame_als_bootstrap(f, R = 12, seed = 1, verbose = FALSE), class = 'lame_als_prediction')
})

test_that("dynamic count and binary dependence use reported predictors", {
	set.seed(76021);n = 12L;nt = 4L;Y = lapply(1:nt, function(t) {z = matrix(rpois(n * n, 3), n);diag(z) = NA_real_;z})
	for (family in c('poisson', 'binary')) {
		ys = if (family == 'binary') lapply(Y, function(y) 1 * (y > 2)) else Y
		f = lame(
			ys,
			R = 0,
			method = 'als',
			family = family,
			dynamic_ab = TRUE,
			dynamic_uv = FALSE,
			als_max_iter = 1000,
			als_tol = 1e-6,
			verbose = FALSE
		)
		expect_true(f$converged)
		# equidispersed counts: the latent variance is zero or not distinguishable from it
		expect_true(
			f$rho_estimation$identified || f$VC['ve'] == 0 || identical(
				f$rho_estimation$reason,
				'latent variance not distinguishable from zero'
			)
		)
		expect_length(f$rho_path_estimation, nt);expect_length(f$rho_path_working, nt)
		eta = array(unlist(f$EZ), dim(f$Y))
		if (family == 'binary') {
			# the binary empirical-Bayes fit of a directed network models the dyadic correlation (lame()'s default): the reported value is the model's,
			# and dcor = FALSE is the independent-cell fit, which reports the correlation of the pairs of cells at the fitted predictor
			expect_false(is.null(f$eb$pair))
			expect_equal(unname(f$VC['rho']), f$eb$pair$rho, tolerance = 1e-12)
			f0 = lame(
				ys,
				R = 0,
				method = 'als',
				family = family,
				dynamic_ab = TRUE,
				dynamic_uv = FALSE,
				dcor = FALSE,
				als_max_iter = 1000,
				als_tol = 1e-6,
				verbose = FALSE
			)
			direct = lame:::.ae_binary_pair_rho(f0$Y, array(unlist(f0$EZ), dim(f0$Y)), f0$link)
			expect_equal(unname(f0$VC['rho']), direct$rho, tolerance = 1e-12)
		} else {
			direct = lame:::.ae_poisson_pair_rho(
				f$Y,
				eta + (if (is.null(f$pois_latent_var)) 0 else f$pois_latent_var) / 2,
				f$VC['ve'],
				overdispersed = isTRUE(f$overdispersion$overdispersed)
			)
			expect_equal(unname(f$VC['rho']), direct$rho, tolerance = 1e-12)
		}
	}
})

test_that("count reciprocity is not estimated when the latent variance is not distinguishable from zero", {
	# equidispersed counts with no latent noise or reciprocity: ve is a small positive root of
	# the Pearson identity and the reciprocity moment sits outside [expm1(-ve), expm1(ve)]
	set.seed(6); n = 40
	a = rnorm(n, 0, .4); b = rnorm(n, 0, .4)
	Y = matrix(rpois(n * n, exp(0.2 + outer(a, b, "+"))), n, n); diag(Y) = NA
	# the empirical-Bayes default reports the same: no correlation where the latent variance is not distinguishable from zero
	fe = suppressMessages(ame_als(Y, family = "poisson", verbose = FALSE))
	expect_false(fe$overdispersion$overdispersed)
	expect_true(is.na(fe$VC["rho"]))
	expect_false(fe$rho_estimation$identified)
	expect_identical(fe$rho_estimation$reason, "latent variance not distinguishable from zero")
	f = suppressMessages(ame_als(Y, family = "poisson", bias_correct = "none", verbose = FALSE))
	expect_gt(unname(f$VC["ve"]), 0)
	expect_lt(f$rho_estimation$moment, expm1(-unname(f$VC["ve"])))
	expect_false(f$overdispersion$overdispersed)
	expect_true(is.na(f$VC["rho"]))
	expect_false(f$rho_estimation$identified)
	expect_identical(f$rho_estimation$reason, "latent variance not distinguishable from zero")
	expect_true(is.finite(f$rho_estimation$covariance))
	out = paste(
		utils::capture.output(print(f), type = "message"),
		utils::capture.output(print(f)),
		collapse = "\n"
	)
	expect_match(out, "rho: not estimated", fixed = TRUE)
	out = paste(
		utils::capture.output(print(summary(f)), type = "message"),
		utils::capture.output(print(summary(f))),
		collapse = "\n"
	)
	expect_match(out, "overdispersion test p", fixed = TRUE)
	# simulation treats the counts as conditionally independent
	expect_no_error(simulate(f, nsim = 2, seed = 1))
	# overdispersed reciprocal counts keep their estimate
	set.seed(101); n = 40
	Z = matrix(rnorm(n * n), n, n)
	E = 0.5 * Z
	E[upper.tri(E)] = (0.5 * (0.5 * t(Z) + sqrt(0.75) * matrix(rnorm(n * n), n, n)))[upper.tri(E)]
	Y = matrix(rpois(n * n, exp(0.2 + outer(a, b, "+") + E)), n, n); diag(Y) = NA
	f = suppressMessages(ame_als(Y, family = "poisson", bias_correct = "none", verbose = FALSE))
	expect_true(f$overdispersion$overdispersed)
	expect_true(f$rho_estimation$identified)
	expect_true(is.finite(f$VC["rho"]))
	# the empirical-Bayes default estimates it from the pairs
	fe = suppressMessages(ame_als(Y, family = "poisson", verbose = FALSE))
	expect_true(fe$rho_estimation$identified)
	expect_identical(fe$rho_estimation$method, "integrated pair likelihood")
	expect_true(is.finite(fe$VC["rho"]))
})
