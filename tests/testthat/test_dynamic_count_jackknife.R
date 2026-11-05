skip_on_cran()

test_that("dynamic count correction preserves means and warm-start coordinates", {
	skip_on_cran()
	set.seed(3112001)
	n = 24L
	Y = X = vector("list", 5L)
	a = rnorm(n, 0, .3)
	b = rnorm(n, 0, .3)
	for (t in seq_along(Y)) {
		x = matrix(rnorm(n * n), n)
		eta = .5 + .5 * x + outer(a, b, "+")
		y = matrix(rpois(n * n, exp(eta + rnorm(n * n, 0, sqrt(.4)))), n)
		diag(y) = NA_real_
		dimnames(y) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1L), dimnames = c(dimnames(y), list("x")))
	}
	args = list(
		Y = Y,
		Xdyad = X,
		R = 0,
		family = "poisson",
		method = "als",
		dynamic_ab = TRUE,
		als_max_iter = 200,
		als_tol = 1e-6,
		verbose = FALSE
	)
	rng = .Random.seed
	raw = do.call(lame, c(args, list(bias_correct = "none")))
	fit = do.call(lame, c(args, list(bias_correct = "jackknife")))
	expect_identical(.Random.seed, rng)
	expect_true(fit$converged)
	expect_true(fit$diagnostics$stationarity$pass)
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_true(all(fit$bias_correct$converged))
	expect_identical(fit$raw$VC, raw$VC)
	expect_identical(fit$raw$coefficients, raw$coefficients)
	expect_identical(fit$fitted, raw$fitted)
	expect_equal(predict(fit), predict(raw), tolerance = 1e-10)
	expect_equal(predict(fit, newdata = fit[["X"]]), predict(raw), tolerance = 1e-10)
	expect_equal(
		lame:::.ae_static_seed(fit)$coef_path,
		lame:::.ae_static_seed(raw)$coef_path,
		tolerance = 1e-12
	)
	expect_equal(unname(coef(fit)), unname(rowMeans(fit$coef_path)))
	expect_equal(unname(fit$BETA[1, ]), unname(coef(fit)))
	expect_equal(unname(fit$mu), mean(fit$coef_path["intercept", ]))
	expect_identical(unname(fit$VC["ve"]), fit$pois_latent_var)
	expect_identical(fit$s2, fit$pois_latent_var)
	expect_identical(fit$objective, raw$objective)
	expect_identical(summary(fit)$bias_correct, fit$bias_correct)
	expect_output(print(fit), "count variance and reciprocity")

	# rebuild the moment correction from the two recorded final estimates
	# and check the response-mean intercept identity
	halves = lapply(fit$bias_correct$attempts, function(x) tail(x, 1L)[[1L]])
	fractions = fit$bias_correct$fraction
	vh = vapply(halves, function(x) unname(x$VC["ve"]), numeric(1))
	ch = vapply(halves, function(x) unname(x$VC["ve"] * x$VC["rho"]), numeric(1))
	expected_v = 2 * unname(raw$VC["ve"]) - sum(fractions * vh)
	expected_c = 2 * unname(raw$VC["ve"] * raw$VC["rho"]) - sum(fractions * ch)
	expect_equal(fit$bias_correct$moments$unrestricted$variance, expected_v)
	expect_equal(fit$bias_correct$moments$unrestricted$covariance, expected_c)
	expect_equal(
		unname(fit$mu + fit$VC["ve"] / 2),
		unname(raw$mu + raw$VC["ve"] / 2),
		tolerance = 1e-12
	)
	refit = ame_als_refit(fit)
	expect_identical(refit$bias_correct$method, "jackknife")
	expect_true(refit$converged)
	uncorrected = ame_als_refit(fit, bias_correct = "none")
	expect_identical(uncorrected$bias_correct$method, "none")
	expect_equal(predict(uncorrected), predict(refit), tolerance = 1e-3)
	expect_identical(.Random.seed, rng)
})

test_that("count moment constraints report boundaries and unidentified correlation", {
	mock = function(v, r) list(VC = c(ve = v, rho = r), rho_path = c(r, r))
	z = lame:::.lda_count_moments(
		mock(.1, .4),
		list(mock(.4, .3), mock(.4, .3)),
		c(.5, .5)
	)
	expect_equal(z$unrestricted$variance, -.2)
	expect_equal(z$variance, 0)
	expect_equal(z$covariance, 0)
	expect_true(is.na(z$rho))
	expect_true(all(is.na(z$rho_path)))
	expect_true(z$boundary)
	z = lame:::.lda_count_moments(
		mock(1, .9),
		list(mock(1, -.9), mock(1, -.9)),
		c(.5, .5)
	)
	expect_equal(z$unrestricted$covariance, 2.7)
	expect_equal(z$variance, 1)
	expect_equal(z$rho, 1)
	expect_true(z$boundary)
	# a zero-variance half still contributes zero covariance when its
	# correlation is unidentified; NA times zero must not propagate
	z = lame:::.lda_count_moments(
		mock(.2, .4),
		list(mock(0, NA_real_), mock(.2, .4)),
		c(.5, .5)
	)
	expect_equal(z$variance, .3)
	expect_equal(z$covariance, .12)
	expect_equal(z$rho, .4)
})

test_that("count splits preserve reciprocal histories and ignore outcomes", {
	for (n in c(4L, 17L, 24L, 60L)) {
		y = array(1, c(n, n, 3L))
		for (t in 1:3) diag(y[, , t]) = NA_real_
		a = lame:::.lda_count_split(y)
		y[is.finite(y)] = seq_len(sum(is.finite(y))) %% 11
		b = lame:::.lda_count_split(y)
		expect_identical(a, b)
		expect_identical(a$mask, aperm(a$mask, c(2, 1, 3)))
		expect_true(all(apply(a$mask, c(1, 3), sum) > 0))
		expect_true(all(apply(!a$mask & is.finite(y), c(1, 3), sum) > 0))
		expect_equal(sum(a$fraction), 1)
	}
})

test_that("the count correction is not stopped by warnings the data raise in every fit", {
	# a row covariate that varies within actor across periods warns in the
	# full fit and in each half; the halves are judged on convergence and
	# stationarity
	set.seed(3)
	n = 20L; nt = 6L
	nm = paste0("a", seq_len(n))
	a = rnorm(n, 0, .4)
	Xr = lapply(seq_len(nt), function(t) matrix(rnorm(n), n, 1, dimnames = list(nm, "xr")))
	Y = lapply(seq_len(nt), function(t) {
		m = matrix(
			rpois(
				n * n,
				exp(
					0.3 + 0.4 * Xr[[t]][, 1] + outer(a, a, "+") +
						rnorm(n * n, 0, .3)
				)
			),
			n,
			n
		)
		diag(m) = NA_real_
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(lame(
		Y,
		Xrow = Xr,
		family = "poisson",
		dynamic_ab = TRUE,
		method = "als",
		bias_correct = "jackknife",
		verbose = FALSE
	))
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_true(fit$bias_correct$applied)
	expect_true(all(fit$bias_correct$converged))
	half_warnings = unlist(lapply(fit$bias_correct$attempts, function(h) h[[1]]$warnings))
	expect_true(any(grepl("vary within actor", half_warnings)))
})

test_that("the count jackknife combines unrestricted covariances", {
	est = function(cv) list(covariance = cv)
	full = list(VC = c(ve = 0.04, rho = -1), rho_estimation = est(-0.06))
	halves = list(
		list(VC = c(ve = 0.05, rho = -1), rho_estimation = est(-0.08)),
		list(VC = c(ve = 0.05, rho = NA), rho_estimation = est(-0.02))
	)
	m = lame:::.lda_count_moments(full, halves, c(0.5, 0.5))
	expect_equal(m$unrestricted$covariance, 2 * -0.06 - 0.5 * (-0.08 - 0.02))
	expect_equal(m$unrestricted$variance, 2 * 0.04 - 0.05)
})
