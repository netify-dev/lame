skip_on_cran()

library(lame)

# residuals(type = "pearson") scales the response residual by the standard
# deviation the family assigns to a response with that fitted value: the
# residual variance for normal, mu (1 - mu) for binary and cbin,
# mu + mu^2 (exp(ve) - 1) for the Poisson-lognormal count, the variance of the
# category distribution for ordinal; frn has none

quick = list(nscan = 60, burn = 20, odens = 2, verbose = FALSE, plot = FALSE, gof = FALSE)
fit_ame = function(d, family, ...) {
	suppressWarnings(suppressMessages(do.call(
		ame,
		c(list(d$Y, Xdyad = d$Xdyad, family = family, R = 0), quick, list(...))
	)))
}
fit_lame = function(d, family, ...) {
	suppressWarnings(suppressMessages(do.call(
		lame,
		c(list(d$Y, Xdyad = d$Xdyad, family = family, R = 0), quick, list(...))
	)))
}
pearson_variance = function(...) lame:::.pearson_variance(...)
# the observed panel of a lame fit, as a list of period matrices
ylist = function(f) if (is.list(f$Y)) f$Y else lapply(seq_len(dim(f$Y)[3]), function(t) f$Y[, , t])

test_that("the variance functions are those of the model", {
	set.seed(1)
	# normal: the posterior mean of the residual variance
	fit = list(family = "normal", VC = cbind(va = c(1, 1), ve = c(1, 3)))
	expect_equal(pearson_variance(fit, matrix(0, 2, 2)), matrix(2, 2, 2))
	# binary and cbin: the Bernoulli variance of the fitted probability
	p = c(.05, .5, .8)
	expect_equal(pearson_variance(list(family = "binary"), p), p * (1 - p))
	expect_equal(pearson_variance(list(family = "cbin"), p), p * (1 - p))
	# poisson: a count with mean mu and lognormal layer of variance ve
	eta = 0.4; ve = 0.5
	z = rnorm(1e6, eta, sqrt(ve))
	y = rpois(1e6, exp(z))
	mu = exp(eta + ve / 2)
	fit = list(family = "poisson", VC = cbind(va = 1, ve = ve))
	expect_equal(pearson_variance(fit, mu), mu + mu^2 * (exp(ve) - 1))
	expect_equal(mean(y), mu, tolerance = 0.01)
	expect_equal(var(y), pearson_variance(fit, mu), tolerance = 0.03)
	# the Poisson variance exceeds the Poisson's own mean by the lognormal layer
	expect_gt(pearson_variance(fit, mu), mu)
	# ordinal: the variance of the category distribution
	alpha = c(0, 0.7, 1.5); lev = c(1, 2, 4, 5)
	for (e in c(-1, 0.3, 1.2)) {
		category = lev[findInterval(e + rnorm(4e5), alpha, left.open = TRUE) + 1]
		expect_equal(lame:::.ordinal_response_variance(e, alpha, lev), var(category), tolerance = 0.02)
	}
	# frn is refused
	expect_error(pearson_variance(list(family = "frn"), matrix(1, 2, 2)), class = "lame_pearson_residuals")
})

test_that("data-induced cutpoints reproduce the category shares of the observed matrix", {
	set.seed(2)
	eta = rnorm(300, 0, 1.3)
	shares = c(0.2, 0.55, 0.8, 1)
	cuts = lame:::.ordinal_induced_cutpoints(eta, shares)
	expect_length(cuts, 3L)
	expect_equal(vapply(cuts, function(c) mean(pnorm(c - eta)), 0), shares[-4], tolerance = 1e-8)
	expect_true(all(diff(cuts) > 0))
})

test_that("ame() Pearson residuals use each family's variance", {
	skip_on_cran()
	n = 10
	# normal: the residual variance of the fit
	d = simulate_test_network(n = n, family = "normal", seed = 3)
	f = fit_ame(d, "normal")
	ve = mean(f$VC[, "ve"])
	expect_equal(residuals(f, "pearson"), (f$Y - f$YPM) / sqrt(ve))
	# binary: sqrt(p (1 - p)) of the fitted probability
	d = simulate_test_network(n = n, family = "binary", seed = 4)
	f = fit_ame(d, "binary")
	expect_equal(residuals(f, "pearson"), (f$Y - f$YPM) / sqrt(f$YPM * (1 - f$YPM)))
	# cbin: the same marginal Bernoulli variance
	d = simulate_test_network(n = n, family = "cbin", odmax = 3, seed = 5)
	f = fit_ame(d, "cbin", odmax = 3)
	expect_equal(residuals(f, "pearson"), (f$Y - f$YPM) / sqrt(f$YPM * (1 - f$YPM)))
	# poisson: mu + mu^2 (exp(ve) - 1)
	d = simulate_test_network(n = n, family = "poisson", seed = 6)
	f = fit_ame(d, "poisson")
	ve = mean(f$VC[, "ve"])
	mu = f$YPM
	expect_equal(residuals(f, "pearson"), (f$Y - mu) / sqrt(mu + mu^2 * (exp(ve) - 1)))
	# the response residual is unchanged
	expect_equal(residuals(f), f$Y - f$YPM)
	# the diagonal and missing cells stay missing
	expect_true(all(is.na(diag(residuals(f, "pearson")))))
	# frn has no cell-wise variance
	d = simulate_test_network(n = n, family = "frn", odmax = 3, seed = 7)
	f = fit_ame(d, "frn", odmax = 3)
	expect_error(residuals(f, "pearson"), class = "lame_pearson_residuals")
	expect_true(is.matrix(residuals(f)))
})

test_that("ame() ordinal Pearson residuals use the variance of the category distribution", {
	skip_on_cran()
	d = simulate_test_network(n = 12, family = "ordinal", seed = 8)
	# data-induced cutpoints: those that give the observed category shares
	f = fit_ame(d, "ordinal")
	eta = lame:::reconstruct_EZ(f)
	diag(eta) = NA
	y = f$Y
	lev = sort(unique(y[is.finite(y)]))
	share = cumsum(table(factor(y[is.finite(y)], levels = lev))) / sum(is.finite(y))
	cuts = vapply(
		share[-length(share)],
		function(s)
			uniroot(function(c) mean(pnorm(c - eta), na.rm = TRUE) - s, c(-15, 15), tol = 1e-12)$root,
		0
	)
	cdf = c(lapply(cuts, function(c) pnorm(c - eta)), list(eta * 0 + 1))
	prob = lapply(seq_along(cdf), function(j) if (j == 1L) cdf[[1L]] else cdf[[j]] - cdf[[j - 1L]])
	m1 = Reduce(`+`, Map(function(l, p) l * p, lev, prob))
	m2 = Reduce(`+`, Map(function(l, p) l^2 * p, lev, prob))
	expect_equal(residuals(f, "pearson"), (y - f$YPM) / sqrt(m2 - m1^2))
	expect_true(all(is.finite(residuals(f, "pearson")[!is.na(y)])))
	# explicit cutpoints: the posterior-mean cutpoints
	fe = fit_ame(d, "ordinal", ordinal_cutpoints = "explicit")
	alpha = c(0, colMeans(fe$ALPHA))
	eta = lame:::reconstruct_EZ(fe)
	diag(eta) = NA
	lev = fe$ordinal_levels
	cdf = c(lapply(alpha, function(a) pnorm(a - eta)), list(eta * 0 + 1))
	prob = lapply(seq_along(cdf), function(j) if (j == 1L) cdf[[1L]] else cdf[[j]] - cdf[[j - 1L]])
	m1 = Reduce(`+`, Map(function(l, p) l * p, lev, prob))
	m2 = Reduce(`+`, Map(function(l, p) l^2 * p, lev, prob))
	expect_equal(residuals(fe, "pearson"), (fe$Y - fe$YPM) / sqrt(m2 - m1^2))
})

test_that("lame() Pearson residuals use each family's variance in every period", {
	skip_on_cran()
	n = 8
	# binary, with a missing cell
	d = simulate_test_network(n = n, n_time = 3, family = "binary", seed = 9)
	d$Y[[2]][3, 4] = NA
	f = fit_lame(d, "binary")
	r = residuals(f, "pearson")
	expect_length(r, 3L)
	for (t in 1:3) {
		expect_equal(r[[t]], (ylist(f)[[t]] - f$YPM[[t]]) / sqrt(f$YPM[[t]] * (1 - f$YPM[[t]])))
	}
	expect_true(is.na(r[[2]][3, 4]))
	# poisson with period exposures: the fitted mean carries the exposure
	d = simulate_test_network(n = n, n_time = 3, family = "poisson", seed = 10)
	f = fit_lame(d, "poisson", period_exposure = c(1, 2, 0.5))
	ve = mean(f$VC[, "ve"])
	r = residuals(f, "pearson")
	for (t in 1:3) {
		mu = f$YPM[[t]]
		expect_equal(r[[t]], (ylist(f)[[t]] - mu) / sqrt(mu + mu^2 * (exp(ve) - 1)))
	}
	# a period with zero exposure has a zero fitted mean, so no standardised residual
	d$Y[[3]][!is.na(d$Y[[3]])] = 0
	f0 = fit_lame(d, "poisson", period_exposure = c(1, 2, 0))
	r0 = residuals(f0, "pearson")
	expect_true(all(is.na(r0[[3]])))
	expect_true(any(is.finite(r0[[1]])))
	# normal
	d = simulate_test_network(n = n, n_time = 2, family = "normal", seed = 11)
	f = fit_lame(d, "normal")
	ve = mean(f$VC[, "ve"])
	r = residuals(f, "pearson")
	expect_equal(r[[2]], (ylist(f)[[2]] - f$YPM[[2]]) / sqrt(ve))
	# ordinal: finite, and the response residual is unchanged
	d = simulate_test_network(n = n, n_time = 2, family = "ordinal", seed = 12)
	f = fit_lame(d, "ordinal")
	r = residuals(f, "pearson")
	for (t in 1:2) {
		ok = !is.na(ylist(f)[[t]])
		expect_true(all(is.finite(r[[t]][ok])))
		expect_equal(sign(r[[t]][ok]), sign((ylist(f)[[t]] - f$YPM[[t]])[ok]))
	}
	expect_equal(residuals(f)[[1]], ylist(f)[[1]] - f$YPM[[1]])
	# frn is refused
	d = simulate_test_network(n = n, n_time = 2, family = "frn", odmax = 3, seed = 13)
	f = fit_lame(d, "frn", odmax = 3)
	expect_error(residuals(f, "pearson"), class = "lame_pearson_residuals")
})

test_that("bipartite fits have Pearson residuals on the rectangular matrix", {
	d = simulate_test_network(n = 10, nA = 10, nB = 7, mode = "bipartite", family = "binary", seed = 14)
	f = suppressWarnings(suppressMessages(ame(
		d$Y,
		Xdyad = d$Xdyad,
		family = "binary",
		R = 0,
		mode = "bipartite",
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	r = residuals(f, "pearson")
	expect_equal(dim(r), c(10, 7))
	expect_equal(r, (f$Y - f$YPM) / sqrt(f$YPM * (1 - f$YPM)))
	# ordinal: the category variance at the rectangular predictor, every cell of a bipartite network
	d = simulate_test_network(n = 10, nA = 10, nB = 7, mode = "bipartite", family = "ordinal", seed = 17)
	fo = suppressWarnings(suppressMessages(ame(
		d$Y,
		Xdyad = d$Xdyad,
		family = "ordinal",
		R = 0,
		mode = "bipartite",
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	ro = residuals(fo, "pearson")
	expect_equal(dim(ro), c(10, 7))
	expect_true(all(is.finite(ro[!is.na(fo$Y)])))
})

test_that("a cell with no variance has no standardised residual", {
	fit = structure(
		list(family = "binary", Y = matrix(c(1, 0, NA, 1), 2), YPM = matrix(c(1, 0, .5, .4), 2)),
		class = "ame"
	)
	r = residuals(fit, "pearson")
	# fitted probabilities of exactly one and zero, and the missing cell
	expect_true(all(is.na(r[1:3])))
	expect_equal(r[4], (1 - .4) / sqrt(.4 * .6))
	# a fit that stores no residual variance cannot be standardised for the normal and Poisson families
	nofit = structure(list(family = "normal", Y = matrix(1, 2, 2), YPM = matrix(0, 2, 2)), class = "ame")
	expect_error(residuals(nofit, "pearson"), class = "lame_pearson_residuals")
	expect_equal(residuals(nofit), matrix(1, 2, 2))
})

test_that("standardised residuals of fits to data from the model have unit scale", {
	skip_on_cran()
	# binary: sd about 1 on the Bernoulli scale of the fitted probability
	d = simulate_test_network(n = 30, family = "binary", beta_dyad = 0.8, s2 = 1, seed = 15)
	f = suppressWarnings(suppressMessages(ame(
		d$Y,
		Xdyad = d$Xdyad,
		family = "binary",
		R = 0,
		nscan = 400,
		burn = 100,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	r = residuals(f, "pearson")
	expect_gt(sd(r, na.rm = TRUE), 0.85)
	expect_lt(sd(r, na.rm = TRUE), 1.15)
	# poisson-lognormal counts
	set.seed(16)
	n = 30
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	z = 0.6 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.7), n)
	Y = matrix(rpois(n * n, exp(z)), n); diag(Y) = NA
	dimnames(Y) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n))
	f = suppressWarnings(suppressMessages(ame(
		Y,
		family = "poisson",
		R = 0,
		nscan = 600,
		burn = 200,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	r = residuals(f, "pearson")
	expect_gt(sd(r, na.rm = TRUE), 0.8)
	expect_lt(sd(r, na.rm = TRUE), 1.25)
})
