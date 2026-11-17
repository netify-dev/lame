skip_on_cran()

library(lame)

# the static intercept of a dynamic binary fit with static factors is not reproduced above a signal per loading of 1.5,
# where the shrunken factor posterior makes it compensate; the fit then takes the intercept as an rw1 path, reported
# at any signal (options(lame.eb_intercept_path_off = TRUE) keeps the withheld static intercept).
# an intercept path (dynamic_ab or dynamic_beta) is kept at any signal

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
diw_fit = function(d, family = "binary", dynamic_ab = FALSE, dynamic_beta = "dyad", ...) {
	suppressWarnings(
		suppressMessages(
			lame(
				d$Y,
				Xdyad = d$X,
				family = family,
				method = "als",
				R = 1,
				dynamic_ab = dynamic_ab,
				dynamic_beta = dynamic_beta,
				verbose = FALSE,
				...
			)
		)
	)
}

test_that("a strong static product makes the static intercept of a binary fit a path", {
	skip_on_cran()
	d = diw_panel("binary", 1.1)
	fit = diw_fit(d)
	expect_gt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.list(fit$eb$intercept_path))
	expect_identical(fit$eb$intercept_path$kind, "rw1")
	expect_true(isTRUE(fit$beta_dynamic_mask[["intercept"]]))
	expect_true(all(is.finite(fit$coef_path["intercept", ])))
	expect_true(all(is.finite(fit$eb$se_path["intercept", ])))
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
	expect_identical(unname(fit$coef_path["intercept", ]), unname(fit$BETA[1, "intercept", ]))
	expect_equal(predict(fit, newdata = fit[["X"]], type = "link"), fit$EZ)
	# the internal switch keeps the withheld static intercept; the slopes stay close to it
	old = options(lame.eb_intercept_path_off = TRUE); on.exit(options(old))
	ref = diw_fit(d)
	options(old)
	expect_true(is.list(ref$eb$intercept_withheld))
	expect_true(is.na(ref$coefficients[["intercept"]]))
	expect_true(is.finite(ref$coefficients[["x_dyad"]]) && is.finite(ref$eb$se[["x_dyad"]]))
	expect_equal(fit$coefficients[["x_dyad"]], ref$coefficients[["x_dyad"]], tolerance = 0.05)
})

test_that("print() and summary() name the intercept path, not a withheld intercept", {
	skip_on_cran()
	fit = diw_fit(diw_panel("binary", 1.1))
	txt = c(capture.output(print(fit)), capture.output(print(summary(fit))))
	expect_false(any(grepl("intercept is withheld", txt)))
	expect_match(paste(txt, collapse = " "), "intercept is fitted as a path")
})

test_that("a weak static product keeps the intercept of a binary fit", {
	skip_on_cran()
	fit = diw_fit(diw_panel("binary", 0.3))
	expect_lt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
})

test_that("a binary fit with dynamic_ab keeps its rw1 intercept path at any signal", {
	skip_on_cran()
	fit = diw_fit(diw_panel("binary", 1.1), dynamic_ab = TRUE, dynamic_beta = FALSE)
	expect_gt(fit$eb$information[["signal_per_loading"]], 1.5)
	expect_null(fit$eb$intercept_withheld)
	expect_true(is.finite(fit$coefficients[["intercept"]]) && is.finite(fit$eb$se[["intercept"]]))
	expect_true(all(is.finite(fit$coef_path["intercept", ])))
	expect_identical(unname(fit$coef_path["intercept", ]), unname(fit$BETA[1, "intercept", ]))
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
