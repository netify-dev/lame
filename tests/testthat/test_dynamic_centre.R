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
		# the reported intercept is the period mean of its path, the slope stays static
		expect_equal(unname(rowMeans(fit$coef_path)), unname(fit$coefficients), info = info)
		expect_equal(
			unname(fit$coef_path["x_dyad", ]),
			rep(unname(fit$coefficients[["x_dyad"]]), 4L),
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
	# each period of the intercept path agrees with the sampler
	for (tt in seq_len(3L)) {
		draws_t = m$BETA[, "intercept", tt]
		expect_lt(
			abs(fit$coef_path["intercept", tt] - mean(draws_t)) / sd(draws_t),
			0.5,
			label = paste("intercept", tt)
		)
	}
	for (nm in c("intercept", "x_dyad")) {
		draws = rowMeans(m$BETA[, nm, ])
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

# the mean of the latent product is identified only jointly with the intercept: a bipartite static product is
# reported at mean zero with the shift in the intercept, and with moving positions the intercept SE carries their common shifts

dcp_panel = function(kind = "directed", n = 18, Tt = 4, seed = 21, R = 1, moving = TRUE) {
	set.seed(seed)
	bip = kind == "bipartite"; sym = kind == "symmetric"
	m = if (bip) n - 4L else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("b%02d", seq_len(m)) else rn
	U = matrix(rnorm(n * R, 0, 0.8), n, R); V = if (sym) U else matrix(rnorm(m * R, 0, 0.8), m, R)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		if (moving) { U = U + matrix(rnorm(n * R, 0, 0.15), n, R); if (sym) V = U else V = V + matrix(rnorm(m * R, 0, 0.15), m, R) }
		x = matrix(rnorm(n * m), n, m)
		if (sym) x[lower.tri(x)] = t(x)[lower.tri(x)]
		y = 0.3 + 0.5 * x + U %*% t(V) + matrix(rnorm(n * m), n, m)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, m, 1), dimnames = list(rn, cn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, bip = bip, sym = sym)
}
dcp_fit = function(d, R = 1, moving = TRUE, ...) {
	args = list(
		d$Y,
		Xdyad = d$X,
		family = "normal",
		method = "als",
		R = R,
		dynamic_uv = moving,
		dynamic_ab = TRUE,
		symmetric = d$sym,
		verbose = FALSE,
		...
	)
	if (d$bip) args$mode = "bipartite"
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

test_that("the centre of a fit is stored with the shift of the product", {
	for (kind in c("directed", "bipartite")) {
		d = dcp_panel(kind, moving = FALSE)
		fit = dcp_fit(d, moving = FALSE)
		expect_equal(names(fit$eb$centre), c("a", "b", "o"), info = kind)
		expect_true(all(is.finite(fit$eb$centre)), info = kind)
	}
})

test_that("a bipartite static product has mean zero and the linear predictor is unchanged by the shift", {
	d = dcp_panel("bipartite", moving = FALSE, seed = 22)
	fit = dcp_fit(d, moving = FALSE)
	pm = Reduce(`+`, fit$UVPM) / length(fit$UVPM)
	if (length(dim(fit$UVPM)) == 3L) pm = apply(fit$UVPM, c(1, 2), mean)
	expect_lt(abs(mean(pm)), 1e-8)
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
})

test_that("the intercept of a fit with moving positions carries the uncertainty of their common shifts", {
	skip_on_cran()
	for (kind in c("symmetric", "directed")) {
		d = dcp_panel(kind, n = 16, Tt = 3, seed = 23)
		fit = dcp_fit(d)
		m = suppressWarnings(suppressMessages(lame(
			d$Y,
			Xdyad = d$X,
			family = "normal",
			symmetric = d$sym,
			R = 1,
			dynamic_uv = TRUE,
			dynamic_ab = TRUE,
			burn = 500,
			nscan = 3000,
			odens = 5,
			seed = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)))
		# the reported intercept and its se are those of the period mean of its path
		draws = rowMeans(m$BETA[, "intercept", ])
		expect_gt(fit$eb$se[["intercept"]] / sd(draws), 0.75, label = kind)
		expect_lt(fit$eb$se[["intercept"]] / sd(draws), 1.35, label = kind)
		expect_lt(abs(fit$coefficients[["intercept"]] - mean(draws)) / sd(draws), 0.7, label = kind)
	}
})
