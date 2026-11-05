# independent checks of the residual-scale extension
skip_on_cran()

test_that("observation scales are explicit and their observed cells enter diagnostics", {
	set.seed(83192)
	nm = LETTERS[1:4]
	Y = setNames(
		lapply(
			1:3,
			function(t) matrix(
				rnorm(16),
				4,
				4,
				dimnames = list(nm, nm)
			)
		),
		as.character(2001:2003)
	)
	Y[[2]][1, 3] = NA_real_
	fit = function(...) lame_robust(
		Y,
		R = 0,
		burn = 2,
		nscan = 12,
		odens = 1,
		seed = 817,
		verbose = FALSE,
		...
	)
	for (bad in list(2, -Inf, NA_real_, c(4, 5), 4 + 1i))
		expect_error(fit(observation_df = bad), "observation_df")
	expect_identical(fit()$draws, fit(observation_df = Inf)$draws)
	for (df in c(4, Inf)) for (actor_df in c(4, Inf)) {
		f = fit(df = df, actor_df = actor_df, observation_df = 4)
		expect_identical(dim(f$draws$omega), c(4L, 4L, 3L, 12L))
		expect_identical(dimnames(f$draws$omega)[1:3], list(nm, nm, names(Y)))
		expect_true(all(is.finite(f$draws$omega) & f$draws$omega > 0))
		expect_true(all(vapply(
			1:12,
			function(s) all(vapply(
				1:3,
				function(t) all(diag(f$draws$omega[, , t, s]) == 1),
				logical(1)
			)),
			logical(1)
		)))
		d = lame:::.lame_robust_identified_draws(f)
		expect_equal(sum(grepl("^observation_precision", colnames(d))), nobs(f))
		expect_true(all(is.finite(d)))
		expect_equal(summary(f)$observation_df, 4)
	}
})

test_that("transient and persistent scale conditionals match the full joint density", {
	set.seed(73191)
	n = 3L; tt = 3L
	Y = lapply(1:tt, function(t) { y = matrix(rnorm(n * n), n, n); diag(y) = NA; y })
	Y[[2]][1, 3] = NA
	D = lame:::.lame_robust_data(Y, NULL, rank = 0L)
	prior = lame:::.lame_robust_prior(0, 0)
	prior$observation_df = 4
	state = list(
		ve = .7,
		lambda = matrix(runif(n * n, .3, 2), n, n),
		omega = array(runif(n * n * tt, .3, 2), c(n, n, tt))
	)
	rr = D$observed * D$Y^2; ss = apply(state$omega * rr, c(1, 2), sum)
	vs = prior$ve_shape + sum(D$observed) / 2; vr = prior$ve_rate + sum(state$lambda * ss) / 2
	keep = rep(D$pairs, tt)
	shapes = c(
		vs,
		(prior$nu + D$counts[D$pairs]) / 2,
		(prior$observation_df + D$observed[keep]) / 2
	)
	N = 3000L
	draws = t(replicate(N, {
		z = lame:::.lame_robust_variance_update(state, D$Y, D, prior)
		lr = (prior$nu + ss[D$pairs] / z$ve) / 2
		or = (prior$observation_df + array(rep(z$lambda, tt), dim(rr)) * rr / z$ve) / 2
		c(vr / z$ve, z$lambda[D$pairs] * lr, z$omega[keep] * or[keep])
	}))
	expect_lt(max(abs((colMeans(draws) - shapes) / sqrt(shapes / N))), 6)
	expect_lt(max(abs((apply(draws, 2, var) - shapes) / sqrt((2 * shapes^2 + 6 * shapes) / (N - 1)))), 6)
	joint = function(v, l, o) {
		a = dgamma(1 / v, prior$ve_shape, prior$ve_rate, log = TRUE) - 2 * log(v) +
			sum(dgamma(l[D$pairs], prior$nu / 2, prior$nu / 2, log = TRUE)) +
			sum(dgamma(o[keep], prior$observation_df / 2, prior$observation_df / 2, log = TRUE))
		for (t in 1:tt) { ok <- D$observed[, , t]; sd <- sqrt(v / (l * o[, , t]));
			a = a + sum(dnorm(D$Y[, , t][ok], 0, sd[ok], log = TRUE)) }
		a
	}
	base = joint(state$ve, state$lambda, state$omega)
	v = state$ve * 1.7; ig = function(v) dgamma(1 / v, vs, vr, log = TRUE) - 2 * log(v)
	errors = joint(v, state$lambda, state$omega) - base - ig(v) + ig(state$ve)
	for (i in which(D$pairs)) {
		l = state$lambda;l[i] = l[i] * 1.4;sh = (prior$nu + D$counts[i]) / 2;rate = (prior$nu + ss[i] / state$ve) / 2
		errors = c(
			errors,
			joint(state$ve, l, state$omega) - base - dgamma(l[i], sh, rate, log = TRUE) + dgamma(
				state$lambda[i],
				sh,
				rate,
				log = TRUE
			)
		)
	}
	rate = (prior$observation_df + array(rep(state$lambda, tt), dim(rr)) * rr / state$ve) / 2
	for (i in which(keep)) {
		o = state$omega;o[i] = o[i] * 1.4;sh = (prior$observation_df + D$observed[i]) / 2
		errors = c(
			errors,
			joint(state$ve, state$lambda, o) - base - dgamma(o[i], sh, rate[i], log = TRUE) + dgamma(
				state$omega[i],
				sh,
				rate[i],
				log = TRUE
			)
		)
	}
	expect_lt(max(abs(errors)), 1e-10)
})

test_that("replication draws fresh temporary shocks with the correct Student law", {
	skip_on_cran()
	set.seed(192)
	Y = lapply(1:2, function(t) matrix(rnorm(16), 4, 4))
	f = lame_robust(
		Y,
		R = 0,
		df = Inf,
		observation_df = 8,
		burn = 2,
		nscan = 1,
		odens = 1,
		seed = 817,
		verbose = FALSE
	)
	f$draws$hyper[1, "ve"] = .7
	# reusing the fitted precisions would give a grossly wrong replicated law
	f$draws$omega[] = 1e-8
	N = 6000L
	p = predict(f, by_draw = TRUE);y = simulate(f, nsim = N, seed = 91328)
	x = cbind(
		(y[1, 2, 1, ] - p[1, 2, 1, 1]) / sqrt(.7),
		(y[1, 2, 2, ] - p[1, 2, 2, 1]) / sqrt(.7),
		(y[2, 1, 1, ] - p[2, 1, 1, 1]) / sqrt(.7)
	)
	expect_lt(max(abs(colMeans(x) / (apply(x, 2, sd) / sqrt(N)))), 6)
	expect_lt(max(abs((colMeans(x^2) - 8 / 6) / (apply(x^2, 2, sd) / sqrt(N)))), 6)
	for (q in c(.01, .05, .5, .95, .99))
		expect_lt(max(abs((colMeans(x <= qt(q, 8)) - q) / sqrt(q * (1 - q) / N))), 6)
	products = cbind(x[, 1] * x[, 2], x[, 1] * x[, 3])
	expect_lt(max(abs(colMeans(products) / (apply(products, 2, sd) / sqrt(N)))), 6)
	# with future mean innovations fixed at zero, the same marginal residual law holds at a positive horizon
	f$draws$hyper[1, c("var_mu", "var_ab")] = 0
	m = predict(f, h = 1);future = simulate(f, nsim = N, h = 1, seed = 17291)
	x = (future[1, 2, 1, ] - m[1, 2, 1]) / sqrt(.7)
	expect_lt(abs(mean(x)) / sqrt(var(x) / N), 6)
	expect_lt(abs(mean(x^2) - 8 / 6) / (sd(x^2) / sqrt(N)), 6)
})
