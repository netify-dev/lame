skip_on_cran()

library(lame)

# a dynamic empirical-Bayes fit reports its effects at mean zero over actors, with the matching intercept
# and the standard error of the identified combination (intercept + mean of the effects)

dce_panel = function(family = "normal", kind = "directed", n = 14, Tt = 4, seed = 8) {
	set.seed(seed)
	bip = kind == "bipartite"; sym = kind == "symmetric"
	m = if (bip) n - 3L else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("b%02d", seq_len(m)) else rn
	al = rnorm(n, 0, 0.6); be = if (sym) al else rnorm(m, 0, 0.6)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * m), n, m)
		if (sym) x[lower.tri(x)] = t(x)[lower.tri(x)]
		eta = 0.4 + outer(al, be, "+") + 0.5 * x + 0.2 * (t - 2)
		y = switch(
			family,
			normal = eta + matrix(rnorm(n * m), n, m),
			binary = 1 * (eta + matrix(rnorm(n * m), n, m) > 0),
			poisson = matrix(rpois(n * m, exp(eta - 0.5)), n, m)
		)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, m, 1), dimnames = list(rn, cn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, bip = bip, sym = sym)
}
dce_fit = function(d, family, ...) {
	args = list(
		d$Y,
		Xdyad = d$X,
		family = family,
		method = "als",
		dynamic_ab = TRUE,
		symmetric = d$sym,
		verbose = FALSE,
		...
	)
	if (d$bip) args$mode = "bipartite"
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

test_that("the reported effects have mean zero over the actors and the fit is unchanged by the centring", {
	for (family in c("normal", "binary", "poisson")) for (kind in c("directed", "symmetric", "bipartite")) {
		d = dce_panel(family, kind)
		fit = dce_fit(d, family)
		info = paste(family, kind)
		expect_lt(abs(mean(fit$a_dynamic)), 1e-10, label = info)
		expect_lt(abs(mean(fit$b_dynamic)), 1e-10, label = info)
		# the linear predictor rebuilds from the reported parts
		expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ, info = info)
		expect_equal(
			unname(fit$coef_path["intercept", ]),
			rep(unname(fit$coefficients[["intercept"]]), 4L),
			info = info
		)
		# the intercept variance is that of the identified combination, not of the intercept alone
		expect_lt(fit$eb$se[["intercept"]], 0.2, label = info)
		expect_true(is.finite(fit$eb$se[["x_dyad"]]), info = info)
	}
})

test_that("the intercept and its standard error agree with the sampler's", {
	skip_on_cran()
	d = dce_panel("normal", "directed", n = 16, Tt = 3, seed = 12)
	fit = dce_fit(d, "normal")
	m = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "normal",
		dynamic_ab = TRUE,
		R = 0,
		burn = 500,
		nscan = 3000,
		odens = 5,
		seed = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	for (nm in c("intercept", "x_dyad")) {
		draws = m$BETA[, nm]
		expect_lt(abs(fit$coefficients[[nm]] - mean(draws)) / sd(draws), 0.5, label = nm)
		expect_gt(fit$eb$se[[nm]] / sd(draws), 0.7)
		expect_lt(fit$eb$se[[nm]] / sd(draws), 1.4)
	}
})

test_that("a moving intercept is centred the same way, with the standard error of each period's combination", {
	d = dce_panel("normal", "directed", n = 16, Tt = 4, seed = 13)
	fit = suppressWarnings(
		suppressMessages(
			lame(
				d$Y,
				Xdyad = d$X,
				family = "normal",
				method = "als",
				dynamic_beta = TRUE,
				dynamic_ab = TRUE,
				verbose = FALSE
			)
		)
	)
	expect_lt(abs(mean(fit$a_dynamic)), 1e-10)
	expect_lt(abs(mean(fit$b_dynamic)), 1e-10)
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	expect_equal(unname(rowMeans(fit$coef_path)), unname(fit$coefficients))
	# the intercept path is combined with the mean of the effects, so its standard errors are small
	expect_true(all(fit$eb$se_path["intercept", ] < 0.5))
	expect_lt(fit$eb$se[["intercept"]], 0.2)
})
