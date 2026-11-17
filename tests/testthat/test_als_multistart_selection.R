skip_on_cran()

test_that("a converged ALS start takes priority over unfinished lower scores", {
	choose = lame:::.ae_multistart_choose
	expect_equal(choose(c(1, 2, 3), c(FALSE, TRUE, TRUE)), 2L)
	expect_equal(choose(c(Inf, 2, 3), c(TRUE, TRUE, TRUE)), 2L)
	expect_equal(choose(c(1, 2, 3), rep(FALSE, 3)), 1L)
	expect_error(choose(c(Inf, NA_real_), c(FALSE, FALSE)), class = 'lame_als_numerical')
})

# stand-ins for .ae_irls_run() results: data deviance, loop objective, the actor-effect covariance at
# convergence, and the fitted predictor with the weights that measure its distance from another start;
# `shift` moves the predictor of every cell and `w` is their weight
fake_start = function(dev, obj, pass = TRUE, sab = diag(c(0.3, 0.2)), shift = 0, w = 0.2) {
	list(
		data_deviance = dev,
		best_obj = obj,
		fit = list(sab = sab),
		eta = matrix(seq(-1, 1, length.out = 12), 3, 4) + shift,
		W = matrix(w, 3, 4),
		stationarity = list(pass = pass)
	)
}
fake_prior = list(kind = "directed", Sab0 = diag(2) * 0.2, eta0 = 8, shift = -3, nu = 1)

test_that("an irls multistart keeps the start with the lowest objective, not the lowest deviance", {
	# the second start fits better (deviance 200) and is rated worse by the
	# objective the loop minimises (348 against 302): inflated scales fit better
	irs = list(fake_start(300, 302), fake_start(200, 348))
	eb = lame:::.ae_irls_select(irs, fake_prior)
	expect_identical(eb$best, 1L)
	expect_identical(eb$criterion, "empirical-Bayes objective")
	# the raw deviance would have kept the other one
	expect_identical(which.min(vapply(irs, `[[`, numeric(1), "data_deviance")), 2L)
	# a converged start still takes priority over an unfinished lower objective
	irs = list(fake_start(300, 302, pass = FALSE), fake_start(310, 320))
	expect_identical(lame:::.ae_irls_select(irs, fake_prior)$best, 2L)
})

test_that("the empirical-Bayes value of a start carries the hyperprior of its covariance", {
	# the same objective of the loop at two covariances: the hyperprior rates
	# the one that is far from the prior covariance worse
	near = fake_start(300, 302, sab = diag(2) * 0.2)
	far = fake_start(300, 302, sab = diag(2) * 3)
	sel = lame:::.ae_irls_select(list(far, near), fake_prior)
	expect_identical(sel$best, 2L)
	expect_gt(sel$value[1L], sel$value[2L])
	# the value is the objective plus c log|Sab| + tr(eta0 Sab0 Sab^-1), c = eta0 + shift
	h = (8 - 3) * log(0.2^2) + 8 * 2
	expect_equal(sel$value[2L], 302 + h)
	# the symmetric and bipartite priors have scalar covariances
	sym = modifyList(fake_prior, list(kind = "symmetric", shift = -2, nu = 0.5))
	expect_equal(
		lame:::.ae_eb_hyper(diag(c(0.4, 0.4)), sym),
		((8 - 2) * log(0.4) + 8 * 0.2 / 0.4) / 0.5
	)
	bip = modifyList(fake_prior, list(kind = "bipartite", shift = -2))
	expect_equal(
		lame:::.ae_eb_hyper(diag(c(0.4, 0.1)), bip),
		sum((8 - 2) * log(c(0.4, 0.1)) + 8 * 0.2 / c(0.4, 0.1))
	)
})

test_that("starts at the same fitted predictor are one optimum, whatever their objectives", {
	# eight starts of a fit whose factors shrink to nothing: one solution,
	# fitted values equal to rounding, objectives 0.02 apart through factor
	# scales that had not settled
	obj = c(-34.2037, -34.2162, -34.1954, -34.2074, -34.2133, -34.1998, -34.2155, -34.2105)
	irs = lapply(obj, function(o) fake_start(178.065, o))
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_identical(sel$best, 2L)
	expect_equal(sel$distance, rep(0, 8))
	expect_false(any(sel$different))
	expect_false(sel$material)
	# however far apart the objectives are: only the fitted values decide
	irs = list(fake_start(178.065, -34.2), fake_start(178.065, 5), fake_start(178.065, 250))
	expect_false(lame:::.ae_irls_select(irs, fake_prior)$material)
})

test_that("a converged start more than a deviance unit away in fitted values is a materially different optimum", {
	# the distance is the weighted squared difference of the predictors, here
	# 12 cells of weight 0.2 shifted by one: a start at a different optimum
	# with an objective 47 higher
	irs = list(fake_start(294, 302), fake_start(198, 349, shift = 1))
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_equal(sel$distance, c(0, 12 * 0.2 * 1^2))
	expect_identical(sel$different, c(FALSE, TRUE))
	expect_true(sel$material)
	expect_equal(sel$excess[2L], 47)
	# the weights of the two starts enter as their mean, so the distance is symmetric
	a = fake_start(294, 302, w = 0.4)
	b = fake_start(198, 349, shift = 1, w = 0.2)
	expect_equal(lame:::.ae_irls_distance(b, a), 12 * 0.3 * 1^2)
	expect_equal(lame:::.ae_irls_distance(a, b), lame:::.ae_irls_distance(b, a))
	# less than a unit apart is one optimum
	irs = list(fake_start(294, 302), fake_start(294, 303, shift = 0.5))
	expect_equal(lame:::.ae_irls_select(irs, fake_prior)$distance[2L], 12 * 0.2 * 0.25)
	expect_false(lame:::.ae_irls_select(irs, fake_prior)$material)
	# the objective does not decide: a distant start rated as well is still a different optimum
	irs = list(fake_start(294, 302), fake_start(250, 302.01, shift = 1))
	expect_true(lame:::.ae_irls_select(irs, fake_prior)$material)
	# a start that has not converged is not an optimum to compare, but its distance is recorded
	irs = list(fake_start(294, 302), fake_start(198, 349, pass = FALSE, shift = 1))
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_false(sel$material)
	expect_equal(sel$distance[2L], 2.4)
	# nothing converged: the kept start is the best that finished and there is nothing to compare
	irs = list(fake_start(294, 302, pass = FALSE), fake_start(198, 349, pass = FALSE, shift = 1))
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_identical(sel$best, 1L)
	expect_false(sel$material)
})

# a directed binary network with a rank-one term and one covariate, fitted at
# rank two: some starts leave the second dimension unused and some use it
ms_network = function(data_seed, n = 30) {
	set.seed(data_seed)
	nm = sprintf("a%03d", seq_len(n))
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	x = matrix(rnorm(n * n), n, dimnames = list(nm, nm))
	P = tcrossprod(matrix(rnorm(n, 0, 0.7), n), matrix(rnorm(n, 0, 0.7), n))
	y = 1 * (-0.8 + 0.7 * x + outer(a, b, "+") + P + matrix(rnorm(n * n), n) > 0)
	diag(y) = NA
	dimnames(y) = list(nm, nm)
	list(Y = y, X = array(x, c(n, n, 1), dimnames = list(nm, nm, "x")))
}

test_that("a multistart fit keeps the coefficients of every start and the distance of their products", {
	d = ms_network(6040, n = 20)
	xr = matrix(seq_len(20) %% 3, ncol = 1, dimnames = list(rownames(d$Y), "xr"))
	xc = matrix(rev(seq_len(20)) %% 4, ncol = 1, dimnames = list(rownames(d$Y), "xc"))
	one = function(...) suppressWarnings(suppressMessages(
		ame_als(
			d$Y,
			d$X,
			xr,
			xc,
			R = 1,
			family = "binary",
			multistart = "cheap",
			seed = 7886,
			verbose = FALSE,
			...
		)
	))
	routes = list(
		eb = one(),
		normal = suppressMessages(ame_als(
			d$Y,
			d$X,
			xr,
			xc,
			R = 1,
			family = "normal",
			multistart = "cheap",
			seed = 7886,
			verbose = FALSE
		))
	)
	for (nm in names(routes)) {
		fit = routes[[nm]]
		sel = fit$multistart_selected
		m = fit$multistart_coefficients
		# a starts by coefficients matrix with the names of coef(fit), node covariates included
		expect_identical(dim(m), c(4L, length(coef(fit))), info = nm)
		expect_identical(colnames(m), names(coef(fit)), info = nm)
		expect_true(all(c("intercept", "x_dyad", "xr_row", "xc_col") %in% colnames(m)), info = nm)
		expect_true(all(is.finite(m)), info = nm)
		# the kept start's row is the vector the fit reports
		expect_identical(m[sel, ], coef(fit), info = nm)
		# the product of the kept start is at zero distance from itself, the others near or far
		dp = fit$multistart_product_distance
		expect_length(dp, 4L)
		expect_identical(dp[sel], 0, info = nm)
		expect_true(all(is.finite(dp) & dp >= 0), info = nm)
	}
	# the summaries leave the kept fit as it is: its reciprocity is estimated as
	# in a single-start fit
	expect_false(is.null(routes$eb$rho_estimation))
	# a fit without a multiplicative term runs no multistart and keeps no summaries
	r0 = suppressWarnings(
		suppressMessages(ame_als(d$Y, d$X, R = 0, family = "binary", multistart = "cheap", verbose = FALSE)),
		classes = "lame_fast_binary_density"
	)
	expect_null(r0$multistart_coefficients)
	expect_null(r0$multistart_product_distance)
})

test_that("the summaries of a multistart follow the assembled starts, and fail softly", {
	prod = function(s) s * matrix(c(1, 2, 3, 4), 2, 2)
	assembled = list(
		list(coefficients = c(a = 1, b = 2), UVPM = prod(1)),
		list(coefficients = c(a = 1.5, b = 2.5), UVPM = prod(1.1)),
		list(coefficients = c(a = 3, b = 4), UVPM = prod(-1))
	)
	sm = lame:::.ae_multistart_summary(function(j) assembled[[j]], 3L, 1L)
	expect_identical(dim(sm$coefficients), c(3L, 2L))
	expect_identical(colnames(sm$coefficients), c("a", "b"))
	expect_identical(sm$coefficients[2L, ], c(a = 1.5, b = 2.5))
	# the Frobenius norm of the difference over the norm of the kept product
	expect_equal(sm$product_distance, c(0, 0.1, 2))
	# the kept start can be any of them, and the distance is measured from it
	expect_equal(
		lame:::.ae_multistart_summary(function(j) assembled[[j]], 3L, 3L)$product_distance,
		c(2, 2.1, 0)
	)
	# without a multiplicative term there is no distance
	none = lapply(assembled, function(z) list(coefficients = z$coefficients, UVPM = NULL))
	sm = lame:::.ae_multistart_summary(function(j) none[[j]], 3L, 2L)
	expect_identical(sm$product_distance, rep(NA_real_, 3L))
	expect_identical(sm$coefficients[3L, ], c(a = 3, b = 4))
	# a start that cannot be assembled, or warns while it is, does not matter
	flaky = function(j) {
		if (j == 2L) stop("no")
		if (j == 3L) warning("careful")
		assembled[[j]]
	}
	expect_no_warning(sm <- lame:::.ae_multistart_summary(flaky, 3L, 1L))
	expect_identical(sm$coefficients[2L, ], c(a = NA_real_, b = NA_real_))
	expect_identical(sm$product_distance[2L], NA_real_)
	expect_equal(sm$product_distance[c(1L, 3L)], c(0, 2))
	# the kept start cannot be missing: nothing is kept then
	expect_null(lame:::.ae_multistart_summary(flaky, 3L, 2L))
	# a kept product of no size is measured against sqrt(eps): equal vanishing products agree,
	# a product of size one against it is a large but finite distance
	tiny = list(
		list(coefficients = c(a = 1), UVPM = prod(1e-30)),
		list(coefficients = c(a = 1), UVPM = prod(2e-30)),
		list(coefficients = c(a = 1), UVPM = prod(1))
	)
	sm = lame:::.ae_multistart_summary(function(j) tiny[[j]], 3L, 1L)
	expect_lt(sm$product_distance[2L], 1e-20)
	expect_true(is.finite(sm$product_distance[3L]) && sm$product_distance[3L] > 1e6)
})

test_that("an empirical-Bayes start with a latent rank is rated by its Laplace marginal likelihood", {
	# the objective of a start moves with the scale of a factor column the
	# data do not use (here by 200 between two starts at one fitted predictor);
	# the marginal likelihood integrates the factors and does not
	with_ev = function(dev, obj, ev, ...) {
		z = fake_start(dev, obj, ...)
		z$fit$U = matrix(0, 3, 2)
		z$fit$eb_vcov = list(evidence = list(value = ev))
		z
	}
	irs = list(with_ev(300, 302, 910.0), with_ev(300, 102, 910.2), with_ev(310, 150, 930, shift = 1))
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_identical(sel$criterion, "Laplace marginal likelihood")
	expect_equal(sel$value, c(910, 910.2, 930))
	expect_identical(sel$best, 1L)
	expect_equal(sel$excess, c(0, 0.2, 20))
	# the objective alone would have kept the second start
	expect_identical(lame:::.ae_irls_select(lapply(irs, function(z) { z$fit$eb_vcov = NULL; z }), fake_prior)$best, 2L)
	# a converged start without a marginal likelihood puts every start back on
	# the objective
	irs[[2]]$fit$eb_vcov = NULL
	sel = lame:::.ae_irls_select(irs, fake_prior)
	expect_identical(sel$criterion, "empirical-Bayes objective")
	expect_identical(sel$best, 2L)
	# a start that has not converged does not need one
	irs = list(with_ev(300, 302, 910.0), with_ev(300, 102, 910.2, pass = FALSE))
	irs[[2]]$fit$eb_vcov = NULL
	expect_identical(lame:::.ae_irls_select(irs, fake_prior)$criterion, "Laplace marginal likelihood")
})

test_that("the starts of a variational fit end at one optimum", {
	skip_on_cran()
	d = ms_network(6040)
	warns = character()
	fit = withCallingHandlers(
		suppressMessages(ame_als(
			d$Y,
			d$X,
			R = 2,
			family = "binary",
			multistart = "cheap",
			seed = 7886,
			verbose = FALSE
		)),
		warning = function(w) {
			warns <<- c(warns, conditionMessage(w))
			invokeRestart("muffleWarning")
		}
	)
	expect_identical(fit$eb$estimator, "variational")
	# the starts are rated by the Laplace marginal likelihood of their modes,
	# which differ; the variational fits that follow do not
	expect_identical(fit$multistart_criterion, "Laplace marginal likelihood")
	ok = fit$multistart_converged
	expect_true(all(ok))
	expect_lt(max(fit$multistart_distance), 1e-2)
	expect_false(any(grepl("materially different local optima", warns)))
	one = suppressMessages(ame_als(d$Y, d$X, R = 2, family = "binary", verbose = FALSE))
	expect_equal(unname(coef(fit)), unname(coef(one)), tolerance = 1e-4)
})

# a sparse directed network without latent structure, fitted at rank two: the
# factors shrink to nothing and every start ends at the same solution
ms_null_network = function(data_seed, n = 40, b0 = -2.8) {
	set.seed(data_seed)
	nm = sprintf("a%03d", seq_len(n))
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	x = matrix(rnorm(n * n), n, dimnames = list(nm, nm))
	y = 1 * (b0 + 0.7 * x + outer(a, b, "+") + matrix(rnorm(n * n), n) > 0)
	diag(y) = NA
	dimnames(y) = list(nm, nm)
	list(Y = y, X = array(x, c(n, n, 1), dimnames = list(nm, nm, "x")))
}

test_that("starts that end at the same fitted values do not warn, whatever their objectives", {
	skip_on_cran()
	# rank 1
	d = ms_null_network(5)
	warns = character()
	fit = withCallingHandlers(
		suppressMessages(ame_als(
			d$Y,
			d$X,
			R = 1,
			family = "binary",
			multistart = "full",
			seed = 7886,
			verbose = FALSE
		)),
		warning = function(w) {
			warns <<- c(warns, conditionMessage(w))
			invokeRestart("muffleWarning")
		}
	)
	expect_true(all(fit$multistart_converged))
	# one solution: the fitted values agree to the tolerance of the iteration
	expect_lt(max(fit$multistart_distance), 1e-4)
	expect_lt(diff(range(fit$multistart_deviance)), 1e-4)
	expect_lt(max(fit$multistart_product_distance), 2e-3)
	expect_lt(max(apply(fit$multistart_coefficients, 2, function(z) diff(range(z)))), 1e-3)
	# the objectives do not: the factor scales still move when the convergence test
	# accepts them, and the objective is not stationary in them
	val = fit$multistart_sse
	# with the objectives equal or apart, the starts are one optimum and do not warn
	expect_lt(diff(range(val)), 0.5)
	expect_false(any(grepl("different local optima", warns)))
})

test_that("the covariance of the actor effects is a minimum of the objective its start is rated by", {
	skip_on_cran()
	d = ms_network(6040, n = 24)
	prep = suppressMessages(lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", FALSE))
	ab = lame:::.ae_eb_prior(prep, FALSE)
	ir = suppressWarnings(lame:::.ae_irls_run(
		prep,
		"binary",
		"probit",
		1L,
		FALSE,
		200,
		1e-6,
		lowrank_method = "als",
		ab_prior = ab
	))
	expect_true(isTRUE(ir$stationarity$pass))
	wsum = rowSums(ir$W, dims = 2L, na.rm = TRUE)
	sab = ir$fit$sab
	J = function(m) lame:::.ae_eb_penalty(ir$fit, sab * m, ab, wsum = wsum) +
		lame:::.ae_eb_hyper(sab * m, ab)
	grid = seq(0.8, 1.25, by = 0.01)
	m_star = grid[which.min(vapply(grid, J, numeric(1)))]
	expect_gt(m_star, 0.9)
	expect_lt(m_star, 1.1)
	# the penalty alone puts the minimum elsewhere: the hyperprior is part of the objective
	Jp = function(m) lame:::.ae_eb_penalty(ir$fit, sab * m, ab, wsum = wsum)
	expect_gt(Jp(1) - min(vapply(grid, Jp, numeric(1))), 1e-3)
})
