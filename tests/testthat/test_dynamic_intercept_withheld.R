skip_on_cran()

library(lame)

# the intercept of a dynamic binary fit with static factors and moving effects is withheld above a signal per loading of 1.5,
# where the shrunken factor posterior makes the intercept compensate; slopes, fitted values, predictions and residuals are kept

diw_panel = function(family = "binary", strength = 1, n = 28, Tt = 4, seed = 31) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n))
	U = matrix(rnorm(n, 0, strength), n, 1); V = matrix(rnorm(n, 0, strength), n, 1)
	al = rnorm(n, 0, 0.5); be = rnorm(n, 0, 0.5)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		al = al + rnorm(n, 0, 0.1); be = be + rnorm(n, 0, 0.1)
		x = matrix(rnorm(n * n), n, n)
		eta = 0.3 + outer(al, be, "+") + 0.6 * x + U %*% t(V)
		y = if (family == "binary") 1 * (eta + matrix(rnorm(n * n), n, n) > 0) else matrix(
			rpois(n * n, exp(0.5 * eta - 0.5)),
			n,
			n
		)
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(rn, rn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}
diw_fit = function(d, family = "binary", ...) {
	suppressWarnings(
		suppressMessages(
			lame(d$Y, Xdyad = d$X, family = family, method = "als", R = 1, dynamic_ab = TRUE, verbose = FALSE, ...)
		)
	)
}

test_that("a small binary panel with a strong static product withholds the intercept, keeps the slope and says so", {
	fit = diw_fit(diw_panel("binary", 1.2, n = 16))
	expect_gt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_true(is.na(fit$coefficients[["intercept"]]))
	expect_true(is.finite(fit$coefficients[["x_dyad"]]) && is.finite(fit$eb$se[["x_dyad"]]))
	expect_match(fit$eb$intercept_withheld$reason, "intercept")
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	expect_match(paste(capture.output(print(fit)), collapse = " "), "intercept is withheld")
})

test_that("a strong static product withholds the intercept of a binary fit and nothing else", {
	skip_on_cran()
	d = diw_panel("binary", 1.1)
	fit = diw_fit(d)
	expect_gt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_true(is.na(fit$coefficients[["intercept"]]))
	expect_true(is.na(fit$eb$se[["intercept"]]))
	expect_true(all(is.na(fit$eb$vcov["intercept", ])) && all(is.na(fit$eb$vcov[, "intercept"])))
	expect_true(all(is.na(fit$coef_path["intercept", ])))
	expect_true(is.finite(fit$coefficients[["x_dyad"]]) && is.finite(fit$eb$se[["x_dyad"]]))
	expect_true(is.list(fit$eb$intercept_withheld))
	expect_match(fit$eb$intercept_withheld$reason, "intercept")
	# withheld, not removed: predictor, fitted values, predictions and residuals match a fit without the rule
	old = options(lame.eb_withhold_intercept = FALSE); on.exit(options(old))
	ref = diw_fit(d)
	options(old)
	expect_null(ref$eb$intercept_withheld)
	expect_true(is.finite(ref$coefficients[["intercept"]]))
	expect_identical(fit$EZ, ref$EZ)
	expect_identical(fit$coefficients[["x_dyad"]], ref$coefficients[["x_dyad"]])
	expect_identical(fit$eb$se[["x_dyad"]], ref$eb$se[["x_dyad"]])
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	expect_equal(predict(fit, newdata = fit[["X"]]), predict(ref, newdata = ref[["X"]]))
	expect_identical(fit$residuals, ref$residuals)
	expect_true(all(is.finite(unlist(fit$YPM)[!is.na(unlist(ref$YPM))])))
})

test_that("print() and summary() state that the intercept is withheld", {
	skip_on_cran()
	fit = diw_fit(diw_panel("binary", 1.1))
	expect_match(paste(capture.output(print(fit)), collapse = " "), "intercept is withheld")
	expect_match(paste(capture.output(print(summary(fit))), collapse = " "), "intercept is withheld")
})

test_that("a weak static product keeps the intercept of a binary fit", {
	skip_on_cran()
	fit = diw_fit(diw_panel("binary", 0.3))
	expect_lt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
})

test_that("a count fit keeps its intercept at any signal", {
	skip_on_cran()
	fit = diw_fit(diw_panel("poisson", 1.1), "poisson")
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
})

test_that("a binary fit with moving positions keeps its intercept", {
	skip_on_cran()
	d = diw_panel("binary", 1.1)
	fit = suppressWarnings(
		suppressMessages(
			lame(d$Y, Xdyad = d$X, family = "binary", method = "als", R = 1, dynamic_uv = TRUE, verbose = FALSE)
		)
	)
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.finite(fit$coefficients[["intercept"]]))
})
