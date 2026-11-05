skip_on_cran()

robust_fixture = function(seed = 741L, rank = 2L, df = 4, burn = 20, nscan = 80, odens = 2) {
	set.seed(14821)
	actors = LETTERS[1:5]
	Y = lapply(1:3, function(t) matrix(rnorm(25), 5, 5, dimnames = list(actors, actors)))
	names(Y) = as.character(2001:2003)
	X = lapply(
		Y,
		function(y) array(
			rnorm(50),
			c(5, 5, 2),
			dimnames = list(actors, actors, c("distance", "lag"))
		)
	)
	list(
		Y = Y,
		X = X,
		fit = lame_robust(
			Y,
			X,
			R = rank,
			df = df,
			burn = burn,
			nscan = nscan,
			odens = odens,
			seed = seed,
			verbose = FALSE
		)
	)
}

test_that("robust inputs align actors, periods and covariates without changing draws", {
	z = robust_fixture()
	order = c(5, 2, 1, 4, 3)
	Y = z$Y
	Y[-1] = lapply(Y[-1], function(y) y[order, rev(order)])
	X = lapply(z$X, function(x) x[order, rev(order), , drop = FALSE])
	X[-1] = lapply(X[-1], function(x) x[, , 2:1, drop = FALSE])
	X = X[3:1]
	f = lame_robust(Y, X, R = 2, burn = 20, nscan = 80, odens = 2, seed = 741, verbose = FALSE)
	expect_identical(f$draws, z$fit$draws)
	new = z$X[1:2]
	permuted = lapply(new, function(x) x[order, rev(order), 2:1, drop = FALSE])
	expect_identical(predict(f, h = 2, newdata = new), predict(f, h = 2, newdata = permuted))
	expect_identical(dim(predict(f, by_draw = TRUE)), c(5L, 5L, 3L, 40L))
	expect_identical(attr(predict(f, by_draw = TRUE), "draw_index"), 1:40)
})

test_that("robust input failures do not silently alter data or priors", {
	z = robust_fixture(rank = 0)
	fit = function(...) lame_robust(z$Y, R = 0, burn = 0, nscan = 4, odens = 1, verbose = FALSE, ...)
	bad = z$X
	bad[[1]][1, 2, 1] = NA_real_
	expect_error(fit(Xdyad = bad), "finite off the diagonal")
	expect_error(fit(df = 2), "df must exceed")
	expect_error(fit(df = 4 + 1i), "df must exceed")
	expect_error(fit(seed = 1 + 1i), "seed must")
	expect_error(fit(prior = list(beta_precision = NULL)), "prior must")
	expect_error(fit(prior = list(ve_shape = 1)), "prior setting")
	expect_error(fit(prior = list(level_scale = matrix(1, 2, 2))), "positive-definite")
	expect_error(fit(prior = list(beta_precision = 1 + 1i)), "beta_precision")
	expect_error(predict(z$fit, h = 1), "one dyadic-covariate array")
	expect_error(predict(z$fit, h = 0, newdata = z$X), "positive forecast horizons")
	expect_error(predict(z$fit, h = 1 + 1i), "h must")
	expect_error(predict(z$fit, unsupported = TRUE), "Unsupported argument")
	expect_error(confint(z$fit, parm = "absent"), "parm must")
	expect_error(confint(z$fit, parm = 0), "parm must")
})

test_that("robust seed scope and rank-zero Gaussian predictions preserve their contracts", {
	z = robust_fixture(rank = 0, df = Inf)
	set.seed(9381)
	old = .Random.seed
	f = lame_robust(z$Y, R = 0, df = Inf, burn = 2, nscan = 4, odens = 1, seed = 11, verbose = FALSE)
	expect_identical(.Random.seed, old)
	expect_true(all(f$draws$lambda == 1))
	expect_identical(dim(predict(f, h = 1)), c(5L, 5L, 1L))
	expect_identical(dim(simulate(f, h = 1, nsim = 1, seed = 12)), c(5L, 5L, 1L, 1L))
	expect_identical(.Random.seed, old)
	expect_identical(nobs(f), 60L)
	expect_equal(fitted(f), predict(f))
	expect_equal(vcov(f), stats::cov(lame:::.lame_robust_coefficient_draws(f)))
	expect_identical(dim(confint(f, parm = names(coef(f))[1])), c(1L, 2L))
	expect_false(any(c("rho_uv", "var_uv") %in% rownames(summary(f)$coefficients)))
	one = lapply(z$X, function(x) unname(x[, , 1]))
	g = lame_robust(z$Y, one, R = 0, burn = 2, nscan = 4, odens = 1, seed = 11, verbose = FALSE)
	expect_identical(dim(predict(g, h = 1, newdata = one[[1]])), c(5L, 5L, 1L))
	# an explicit seed also restores the absence of a caller RNG state
	rm(".Random.seed", envir = .GlobalEnv)
	on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
	invisible(simulate(g, seed = 11))
	expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("robust joint simulation uses each saved persistent pair variance", {
	skip_on_cran()
	z = robust_fixture()
	f = z$fit
	mean = predict(f, by_draw = TRUE)
	N = 5000L
	y = simulate(f, nsim = N, seed = 73201)
	index = attr(y, "draw_index")
	standardized = sapply(
		1:3,
		function(t)
			(y[1, 2, t, ] - mean[1, 2, t, index]) /
			sqrt(f$draws$hyper[index, "ve"] / f$draws$lambda[1, 2, index])
	)
	expect_lt(max(abs(colMeans(standardized))), 6 / sqrt(N))
	expect_lt(max(abs(apply(standardized, 2L, var) - 1)), 6 * sqrt(2 / (N - 1)))
	expect_lt(max(abs(stats::cor(standardized)[upper.tri(diag(3))])), 6 / sqrt(N))
	forecast = predict(f, h = 2, newdata = z$X[1:2])
	future = predict(f, h = 2, newdata = z$X[1:2], by_draw = TRUE, nsim = N, seed = 71391)
	for (t in 1:2) {
		x = future[1, 2, t, ]
		expect_lt(abs(mean(x) - forecast[1, 2, t]) / (sd(x) / sqrt(N)), 6)
	}
})

test_that("robust diagnostics include identified products and reject invalid chains", {
	skip_on_cran()
	skip_if_not_installed("posterior")
	f = robust_fixture(seed = 741, burn = 500, nscan = 2000, odens = 4)$fit
	g = robust_fixture(seed = 742, burn = 500, nscan = 2000, odens = 4)$fit
	d = compute_mcmc_diagnostics(list(f, g))
	expect_equal(sum(grepl("^UV\\[", d$param_names)), 5 * 4 * 3)
	expect_equal(sum(grepl("^precision\\[", d$param_names)), 5 * 4)
	x = cbind(f$draws$mu[, 1], g$draws$mu[, 1])
	W = mean(apply(x, 2, var))
	B = nrow(x) * var(colMeans(x))
	expect_equal(unname(d$rhat_classic["intercept[2001]"]), sqrt(((nrow(x) - 1) * W + B) / (nrow(x) * W)))
	expect_identical(unname(d$rhat["intercept[2001]"]), posterior::rhat(x))
	expect_identical(unname(d$ess_bulk["intercept[2001]"]), posterior::ess_bulk(x))
	expect_error(compute_mcmc_diagnostics(list(f, f)), "different sampler seeds")
	expect_error(compute_mcmc_diagnostics(list(f)), "at least two")
	expect_error(compute_mcmc_diagnostics(list(f, list())), "independently seeded")
	changed = g
	changed$prior$ve_rate = changed$prior$ve_rate * 2
	expect_error(compute_mcmc_diagnostics(list(f, changed)), "same data, model, and priors")
	shifted = g
	shifted$draws$mu[, 1] = shifted$draws$mu[, 1] + 100
	expect_false(compute_mcmc_diagnostics(list(f, shifted))$pass)
	f0 = robust_fixture(seed = 741, rank = 0, df = Inf, burn = 500, nscan = 2000, odens = 4)$fit
	g0 = robust_fixture(seed = 742, rank = 0, df = Inf, burn = 500, nscan = 2000, odens = 4)$fit
	d0 = compute_mcmc_diagnostics(list(f0, g0))
	expect_false(any(grepl("rho_uv|var_uv|UV\\[|precision\\[", d0$param_names)))
})

test_that("actor volatility is explicit, aligned and included in diagnostics", {
	z = robust_fixture(rank = 0, burn = 2, nscan = 12, odens = 1)
	fit = function(actor_df) lame_robust(
		z$Y,
		z$X,
		R = 0,
		burn = 2,
		nscan = 12,
		odens = 1,
		seed = 741,
		verbose = FALSE,
		actor_df = actor_df
	)
	expect_identical(fit(Inf)$draws, z$fit$draws)
	for (bad in list(2, -Inf, NA_real_, c(4, 5), 4 + 1i)) expect_error(fit(bad), "actor_df")
	f = fit(4)
	expect_identical(dim(f$draws$tau_ab), c(2L, 3L, 12L))
	expect_true(all(is.finite(f$draws$tau_ab) & f$draws$tau_ab > 0))
	expect_equal(dimnames(f$draws$tau_ab)[[1]], c("sender", "receiver"))
	expect_equal(dimnames(f$draws$tau_ab)[[2]], names(z$Y))
	d = lame:::.lame_robust_identified_draws(f)
	expect_equal(sum(grepl("^actor_precision", colnames(d))), 6)
	expect_equal(summary(f)$actor_df, 4)
})

test_that("future actor volatility has the independently integrated variance", {
	skip_on_cran()
	z = robust_fixture(rank = 0, burn = 2, nscan = 4, odens = 1)
	f = lame_robust(
		z$Y,
		R = 0,
		df = Inf,
		actor_df = 8,
		burn = 2,
		nscan = 1,
		odens = 1,
		seed = 713,
		verbose = FALSE
	)
	f$draws$hyper[1, c("var_mu", "var_ab", "rho_mu", "rho_ab")] = c(.11, .7, .5, .6)
	N = 10000L
	point = predict(f, h = 2)
	future = predict(f, h = 2, by_draw = TRUE, nsim = N, seed = 81239)
	n = f$data$n
	for (t in 1:2) {
		err = future[1, 2, t, ] - point[1, 2, t]
		variance = .11 * sum(.5^(2 * (0:(t - 1)))) +
			2 * (1 - 1 / n) * .7 * (8 / (8 - 2)) * sum(.6^(2 * (0:(t - 1))))
		expect_lt(abs(mean(err)), 6 * sqrt(variance / N))
		expect_lt(abs(mean(err^2) - variance), 6 * sd(err^2) / sqrt(N))
	}
})
