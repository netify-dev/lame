skip_on_cran()

library(lame)

# estimand of fitted(), predict(type = "response") and residuals() for dynamic ALS fits (R/dynamic_als_predictive.R):
# normal: the predictor; count: exp(eta + ve / 2); binary: predictive probability (empirical Bayes) or plug-in (penalised)

dpr_panel = function(family, n = 12, Tt = 3, seed = 5, sym = FALSE) {
	set.seed(seed)
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	al = rnorm(n, 0, 0.6)
	be = if (sym) al else rnorm(n, 0, 0.6)
	w = rnorm(n)
	Y = X = Wr = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x1"))
		if (sym) x[, , 1] = (x[, , 1] + t(x[, , 1])) / 2
		eta = -0.3 + outer(al, be, "+") + 0.5 * x[, , 1] + if (sym) outer(w, w, "+") * 0.3 else 0
		y = switch(
			family,
			binary = 1 * (eta + rnorm(n * n) > 0),
			poisson = matrix(rpois(n * n, exp(eta)), n),
			normal = eta + rnorm(n * n)
		)
		if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
		X[[t]] = x
		Wr[[t]] = matrix(w, n, 1, dimnames = list(rn, "w"))
	}
	names(Y) = names(X) = names(Wr) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X, Xrow = Wr)
}
dpr_fit = function(d, family, bias_correct = NULL, ...) {
	args = c(list(d$Y, d$X, family = family, dynamic_ab = TRUE, verbose = FALSE), list(...))
	if (!is.null(bias_correct)) args$bias_correct = bias_correct
	suppressWarnings(do.call(lame:::lame_dynamic_als, args))
}

test_that("the fitted values follow the rules of the package: normal the predictor, count the plug-in mean, binary the predictive probability", {
	for (family in c("normal", "poisson", "binary")) {
		d = dpr_panel(family)
		for (route in c("eb", "none")) {
			fit = dpr_fit(d, family, if (route == "none") "none")
			info = paste(family, route)
			fv = fitted(fit)
			expect_length(fv, 3L)
			expect_identical(predict(fit), fv, info = info)
			expect_equal(predict(fit, type = "link"), fit$EZ, info = info)
			# stored values stay the plug-in ones, as for the static estimators
			expect_identical(fitted(fit, plug_in = TRUE), fit$fitted, info = info)
			# the penalised count fit has no latent layer, so its noise variance is zero
			ve = if (is.null(fit$pois_latent_var)) 0 else fit$pois_latent_var
			for (t in 1:3) {
				ez = fit$EZ[[t]]
				expected = switch(
					family,
					normal = ez,
					poisson = exp(ez + ve / 2),
					binary = pnorm(ez)
				)
				expect_equal(fit$fitted[[t]], expected, info = info)
				if (family != "binary" || route == "none") expect_equal(fv[[t]], expected, info = info)
			}
			if (family == "poisson") expect_equal(ve > 0, route == "eb")
			# the one fit whose fitted values are not the plug-in ones
			expect_equal(isTRUE(all.equal(fv, fit$fitted)), !(family == "binary" && route == "eb"), info = info)
			# residuals are response minus fitted values, plug-in ones with plug_in = TRUE
			r = residuals(fit)
			expect_equal(r[[2]], fit$Y[, , 2] - fv[[2]], ignore_attr = TRUE, info = info)
			expect_equal(
				residuals(fit, plug_in = TRUE)[[2]],
				fit$Y[, , 2] - fit$fitted[[2]],
				ignore_attr = TRUE,
				info = info
			)
			expect_identical(residuals(fit, type = "working"), fit$residuals, info = info)
		}
	}
})

test_that("an empirical-Bayes binary fit reports the predictive probability, the link integrated over the posterior of the predictor", {
	d = dpr_panel("binary")
	fit = dpr_fit(d, "binary")
	v = lame:::.lde_predictive_var(fit, fit[["X"]])
	pr = fitted(fit)
	for (t in 1:3) {
		expect_true(all(v[[t]][!is.na(fit$EZ[[t]])] > 0))
		expect_equal(pr[[t]], pnorm(fit$EZ[[t]] / sqrt(1 + v[[t]])))
		# integration pulls a probability towards one half but never past it
		ok = !is.na(fit$EZ[[t]])
		expect_true(all(abs(pr[[t]][ok] - 0.5) <= abs(fit$fitted[[t]][ok] - 0.5) + 1e-12))
		expect_true(all(is.na(diag(pr[[t]]))))
	}
	expect_gt(max(abs(unlist(pr) - unlist(fit$fitted)), na.rm = TRUE), 1e-3)
	# the penalised fit has no predictor posterior, so its probability is the plug-in one
	fx = dpr_fit(d, "binary", "none")
	expect_equal(fitted(fx), fx$fitted)
	# the logit link uses its own approximation
	fl = dpr_fit(d, "binary", link = "logit")
	vl = lame:::.lde_predictive_var(fl, fl[["X"]])
	expect_equal(fitted(fl)[[1]], plogis(fl$EZ[[1]] / sqrt(1 + pi * vl[[1]] / 8)))
	# a fit without the predictive record is refused
	fit$eb$predictive = NULL
	expect_error(fitted(fit), class = "lame_predictive_missing")
	expect_s3_class(tryCatch(fitted(fit), error = identity), "lame_predictive_missing")
	expect_identical(fitted(fit, plug_in = TRUE), fit$fitted)
})

test_that("predict(newdata =) follows the same rules, at the variance of the design it is asked at", {
	d = dpr_panel("binary", n = 10, Tt = 3, seed = 8, sym = TRUE)
	fit = suppressWarnings(
		lame:::lame_dynamic_als(
			d$Y,
			d$X,
			Xrow = d$Xrow,
			family = "binary",
			symmetric = TRUE,
			dynamic_ab = TRUE,
			verbose = FALSE
		)
	)
	# the design of the fit reproduces the fitted values, with node covariates and a symmetric network
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	expect_equal(predict(fit, newdata = fit[["X"]]), fitted(fit))
	# at another design the link moves with the coefficient and the probability integrates over that design's variance
	Xn = fit[["X"]] * 3
	lk = predict(fit, newdata = Xn, type = "link")
	pr = predict(fit, newdata = Xn)
	vn = lame:::.lde_predictive_var(fit, Xn)
	vf = lame:::.lde_predictive_var(fit, fit[["X"]])
	for (t in 1:3) {
		expect_equal(pr[[t]], pnorm(lk[[t]] / sqrt(1 + vn[[t]])))
		# x1 has a variance, so a larger covariate gives a larger predictor variance
		expect_gt(mean(vn[[t]] - vf[[t]], na.rm = TRUE), 0)
	}
	# the other families give the plug-in value at the new design
	dc = dpr_panel("poisson", n = 10, Tt = 3, seed = 9)
	fc = dpr_fit(dc, "poisson")
	pc = predict(fc, newdata = fc[["X"]] * 2)
	lc = predict(fc, newdata = fc[["X"]] * 2, type = "link")
	expect_gt(fc$pois_latent_var, 0)
	expect_equal(pc[[2]], exp(lc[[2]] + fc$pois_latent_var / 2))
	dn = dpr_panel("normal", n = 10, Tt = 3, seed = 10)
	fn = dpr_fit(dn, "normal")
	expect_equal(predict(fn, newdata = fn[["X"]] * 2), predict(fn, newdata = fn[["X"]] * 2, type = "link"))
})

test_that("the dynamic route has no forecasts: a horizon is refused by name", {
	d = dpr_panel("binary")
	fit = dpr_fit(d, "binary")
	expect_error(predict(fit, h = 1), "h-step forecasting")
	expect_error(predict(fit, h = 0.5), "non-negative integer")
})

test_that("the cells of an actor with no observed cell carry that actor's prior variance", {
	d = dpr_panel("binary", n = 12, Tt = 3, seed = 5)
	d$Y = lapply(d$Y, function(y) { y[3, ] = NA; y[, 3] = NA; y })
	fit = dpr_fit(d, "binary")
	ve = fit$eb$predictive$var_effects
	prior = fit$eb$sab[1, 1] + fit$sigma_ab^2 / (1 - fit$rho_ab^2)
	expect_gt(min(ve[3, -3, ]), prior)
	expect_true(all(ve[1, 2, ] < prior))
	expect_true(all(is.finite(unlist(fitted(fit))[!is.na(unlist(fit$EZ))])))
})
