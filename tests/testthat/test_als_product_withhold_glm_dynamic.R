# the latent product the data do not determine is withheld in binary and count fits and in the dynamic
# normal estimator (R/ame_als_gauss.R, R/lame_dynamic_als_eb.R)

pw_glm = function(fam, suv, n = 40, seed = 1, kind = "directed") {
	set.seed(seed)
	sym = kind == "symmetric"
	U = matrix(rnorm(n, 0, sqrt(suv)), n); V = if (sym) U else matrix(rnorm(n, 0, sqrt(suv)), n)
	X = matrix(rnorm(n * n), n); if (sym) X[lower.tri(X)] = t(X)[lower.tri(X)]
	E = matrix(rnorm(n * n), n); if (sym) E[lower.tri(E)] = t(E)[lower.tri(E)]
	a = rnorm(n, 0, 0.5); b = if (sym) a else rnorm(n, 0, 0.5)
	z = outer(a, b, "+") + U %*% t(V) + 0.5 * X + E
	Y = if (fam == "binary") 1 * (z > stats::quantile(z, 0.75)) else {
		y = matrix(rpois(n * n, exp(0.6 * (z - mean(z)))), n)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		y
	}
	diag(Y) = NA
	rn = sprintf("a%02d", seq_len(n)); dimnames(Y) = list(rn, rn)
	list(Y = Y, X = array(X, c(n, n, 1), dimnames = list(rn, rn, "x")), sym = sym, fam = fam)
}
pw_fit = function(d, R = 1) suppressWarnings(
	ame_als(d$Y, Xdyad = d$X, R = R, family = d$fam, symmetric = d$sym, verbose = FALSE)
)
pw_print = function(f) paste(utils::capture.output(print(f)), collapse = " ")

test_that("a binary or count product below the information threshold is withheld, and one above it is reported", {
	skip_on_cran()
	for (fam in c("binary", "poisson")) {
		weak = pw_fit(pw_glm(fam, 0.02, seed = 3))
		info = weak$eb$variational$information
		expect_true(
			all(c("events", "per_loading", "signal", "signal_per_loading", "threshold") %in% names(info)),
			label = fam
		)
		expect_equal(info[["threshold"]], 2.1)
		expect_lt(info[["signal_per_loading"]], 2.1)
		expect_true(lame:::.ae_product_is_withheld(weak), label = fam)
		expect_true(all(is.na(weak$U)) && all(is.na(weak$UVPM)), label = fam)
		expect_identical(weak$eb$product_withheld$family, fam)
		expect_match(pw_print(weak), "units of signal per loading", fixed = TRUE)
		expect_error(uv_plot(weak), class = "lame_product_withheld")
		strong = pw_fit(pw_glm(fam, 1.5, seed = 4))
		expect_gte(strong$eb$variational$information[["signal_per_loading"]], 2.1)
		expect_false(lame:::.ae_product_is_withheld(strong), label = fam)
		expect_true(all(is.finite(strong$UVPM)), label = fam)
	}
})

test_that("a withheld binary or count fit has the predictive values, coefficients and standard errors of the fit", {
	skip_on_cran()
	for (fam in c("binary", "poisson")) {
		d = pw_glm(fam, 0.02, seed = 3)
		f = pw_fit(d)
		g = withr::with_options(list(lame.eb_withhold_product = FALSE), pw_fit(d))
		expect_true(lame:::.ae_product_is_withheld(f)); expect_false(lame:::.ae_product_is_withheld(g))
		# predictive probabilities and counts integrate the variational posterior of the product
		expect_identical(fitted(f), fitted(g))
		expect_identical(fitted(f, plug_in = FALSE), fitted(g, plug_in = FALSE))
		expect_identical(predict(f, type = "link"), predict(g, type = "link"))
		expect_identical(residuals(f), residuals(g))
		expect_identical(coef(f), coef(g))
		expect_identical(vcov(f), vcov(g))
		expect_identical(f$eb$product_withheld$means$UVPM, g$UVPM)
		expect_true(all(is.finite(fitted(f)[!is.na(d$Y)])))
	}
})

pw_dyn = function(cfg, suv, n = 20, Tt = 3L, seed = 1) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n))
	U = array(rnorm(n * Tt, 0, sqrt(suv)), c(n, 1, Tt)); V = array(rnorm(n * Tt, 0, sqrt(suv)), c(n, 1, Tt))
	if (cfg == "dab") for (t in 2:Tt) { U[, , t] = U[, , 1]; V[, , t] = V[, , 1] } else for (t in 2:Tt) {
		U[, , t] = 0.9 * U[, , t - 1] + sqrt(1 - 0.81) * U[, , t]; V[, , t] = 0.9 * V[, , t - 1] + sqrt(1 - 0.81) * V[, , t] }
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n)
		y = 0.2 + outer(
			a + rnorm(n, 0, 0.2),
			b + rnorm(n, 0, 0.2),
			"+"
		) + U[, , t] %*% t(V[, , t]) + 0.5 * x + matrix(rnorm(n * n), n)
		diag(y) = NA; dimnames(y) = list(rn, rn)
		Y[[t]] = y; X[[t]] = array(x, c(n, n, 1), dimnames = list(rn, rn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, cfg = cfg)
}
pw_dfit = function(d) suppressWarnings(suppressMessages(lame(
	d$Y,
	Xdyad = d$X,
	R = 1,
	family = "normal",
	method = "als",
	verbose = FALSE,
	dynamic_ab = d$cfg == "dab",
	dynamic_uv = d$cfg == "duv"
)))

test_that("the dynamic normal estimator withholds static factors and moving positions the data do not determine", {
	skip_on_cran()
	for (cfg in c("dab", "duv")) {
		weak = pw_dfit(pw_dyn(cfg, 0.02, seed = 5))
		info = weak$eb$information
		expect_equal(info[["threshold"]], if (cfg == "duv") 1.3 else 2.1, label = cfg)
		expect_lt(info[["per_loading"]], info[["threshold"]])
		expect_true(lame:::.ae_product_is_withheld(weak), label = cfg)
		expect_true(all(is.na(weak$UVPM)) && all(is.na(weak$U)) && all(is.na(weak$Omat)), label = cfg)
		expect_identical(isTRUE(weak$eb$product_withheld$moving), cfg == "duv")
		expect_match(pw_print(weak), "not identified by the data at this rank", fixed = TRUE)
		lp = suppressMessages(latent_positions(weak))
		expect_true(all(is.na(lp$value)))
		expect_error(uv_plot(weak), class = "lame_product_withheld")
		strong = pw_dfit(pw_dyn(cfg, 0.8, seed = 6))
		expect_gte(strong$eb$information[["per_loading"]], strong$eb$information[["threshold"]])
		expect_false(lame:::.ae_product_is_withheld(strong), label = cfg)
	}
})

test_that("a withheld dynamic fit has the fitted values, predictions and coefficients of the fit", {
	skip_on_cran()
	for (cfg in c("dab", "duv")) {
		d = pw_dyn(cfg, 0.02, seed = 5)
		f = pw_dfit(d)
		g = withr::with_options(list(lame.eb_withhold_product = FALSE), pw_dfit(d))
		expect_true(lame:::.ae_product_is_withheld(f)); expect_false(lame:::.ae_product_is_withheld(g))
		expect_identical(fitted(f), fitted(g))
		expect_identical(residuals(f), residuals(g))
		Xa = array(unlist(d$X), c(20, 20, 1, 3))
		expect_identical(predict(f, newdata = Xa), predict(g, newdata = Xa))
		expect_identical(coef(f), coef(g))
		expect_identical(f$eb$product_withheld$means$UVPM, g$UVPM)
		expect_identical(f$eb$product_withheld$means$Omat, g$Omat)
	}
})
