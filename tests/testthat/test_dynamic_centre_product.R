skip_on_cran()

library(lame)

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
		draws = m$BETA[, "intercept"]
		expect_gt(fit$eb$se[["intercept"]] / sd(draws), 0.75, label = kind)
		expect_lt(fit$eb$se[["intercept"]] / sd(draws), 1.35, label = kind)
		expect_lt(abs(fit$coefficients[["intercept"]] - mean(draws)) / sd(draws), 0.7, label = kind)
	}
})
