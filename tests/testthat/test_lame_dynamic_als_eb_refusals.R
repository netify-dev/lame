skip_on_cran()

library(lame)

# the dynamic empirical-Bayes estimator returns a converged fit, with any position scale at an interior minimiser
# of the profile objective, or refuses with a classed error whose `fit` field holds the last iterate

rf_panel = function(family = "normal", n = 12, Tt = 4, seed = 3) {
	set.seed(seed)
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	a = rnorm(n, 0, 0.7)
	b = rnorm(n, 0, 0.7)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = array(rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x1"))
		eta = -0.3 + outer(a, b, "+") + 0.4 * x[, , 1]
		y = switch(
			family,
			normal = eta + rnorm(n * n),
			poisson = matrix(rpois(n * n, exp(eta)), n),
			binary = 1 * (eta + rnorm(n * n) > 0)
		)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		Y[[t]] = y
		X[[t]] = x
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

test_that("an exhausted budget is a classed error with the state, for each family", {
	for (fam in c("normal", "poisson", "binary")) {
		d = rf_panel(fam)
		err = expect_error(
			lame:::lame_dynamic_als(d$Y, d$X, family = fam, dynamic_ab = TRUE, max_iter = 2, verbose = FALSE),
			class = "lame_als_convergence"
		)
		expect_s3_class(err, "lame_als_eb")
		msg = gsub("\\s+", " ", conditionMessage(err))
		expect_match(
			msg,
			"did not reach the empirical-Bayes convergence criterion: it stopped after 2 passes (iteration limit)",
			fixed = TRUE
		)
		expect_match(msg, "largest scaled residual of the empirical-Bayes equations is")
		# the condition is an error, not a warning
		expect_false(inherits(err, "warning"))
		fit = err$fit
		expect_s3_class(fit, "lame_dynamic_als")
		expect_false(fit$converged)
		expect_equal(fit$iterations, 2L)
		expect_equal(fit$family, fam)
		# the same panel at the defaults converges with no condition
		ok = expect_no_warning(lame:::lame_dynamic_als(d$Y, d$X, family = fam, dynamic_ab = TRUE, verbose = FALSE))
		expect_true(ok$converged)
		expect_lt(ok$iterations, 200L)
	}
})

test_that("the budget of the empirical-Bayes iteration is that of lame(method = 'als') unless the caller sets one", {
	d = rf_panel("normal")
	seen = NULL
	testthat::local_mocked_bindings(
		.lde_dynamic_fit = function(prep, family, mode, symmetric, link, prior, max_iter, ...) {
			seen <<- c(seen, max_iter)
			stop("stop here")
		},
		.package = "lame"
	)
	try(lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, verbose = FALSE), silent = TRUE)
	try(
		lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, max_iter = 37, verbose = FALSE),
		silent = TRUE
	)
	try(lame(d$Y, Xdyad = d$X, family = "normal", method = "als", dynamic_ab = TRUE, verbose = FALSE), silent = TRUE)
	expect_equal(seen, c(1000, 37, 1000))
})

test_that("a replay of an empirical-Bayes fit has the budget the fit used", {
	for (fam in c("normal", "poisson", "binary")) {
		d = rf_panel(fam, n = 8, Tt = 3)
		f = lame:::lame_dynamic_als(d$Y, d$X, family = fam, dynamic_ab = TRUE, verbose = FALSE)
		expect_equal(f$als_refit$args$max_iter, 1000, info = fam)
		g = lame:::lame_dynamic_als(d$Y, d$X, family = fam, dynamic_ab = TRUE, max_iter = 300, verbose = FALSE)
		expect_equal(g$als_refit$args$max_iter, 300, info = fam)
	}
})

test_that("the bootstrap of a dynamic empirical-Bayes normal fit is refused, as lame_dynamic_als(bootstrap = ) refuses it", {
	d = rf_panel("normal", n = 8, Tt = 3)
	fit = lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
	expect_identical(fit$bias_correct$method, "eb")
	err = expect_error(ame_als_bootstrap(fit, R = 3, verbose = FALSE), class = "lame_dynamic_eb_unsupported")
	expect_match(
		gsub("\\s+", " ", conditionMessage(err)),
		"A bootstrap is not offered for the dynamic empirical-Bayes estimator"
	)
	err2 = expect_error(
		lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, bootstrap = 3, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
	expect_match(
		gsub("\\s+", " ", conditionMessage(err2)),
		"A bootstrap is not offered for the dynamic empirical-Bayes estimator"
	)
	# the penalised fit of the same panel is not refused
	fit0 = lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "normal",
		dynamic_ab = TRUE,
		verbose = FALSE,
		bias_correct = "none"
	)
	res = tryCatch(suppressWarnings(ame_als_bootstrap(fit0, R = 3, verbose = FALSE)), error = function(e) e)
	expect_false(inherits(res, "lame_dynamic_eb_unsupported"))
})

test_that("a position scale search that ends at an edge or unsettled is refused for normal and count panels as for binary", {
	skip_on_cran()
	real = lame:::.lde_pos_profile
	for (fam in c("normal", "poisson")) {
		d = rf_panel(fam, n = 10, Tt = 3)
		testthat::local_mocked_bindings(
			.lde_pos_profile = function(...) { pr = real(...); pr$at_edge = TRUE; pr$S = pr$range[2L]; pr },
			.package = "lame"
		)
		err = expect_error(
			suppressMessages(
				lame:::lame_dynamic_als(
					d$Y,
					d$X,
					family = fam,
					dynamic_ab = TRUE,
					R = 1,
					dynamic_uv = TRUE,
					verbose = FALSE
				)
			),
			class = "lame_dynamic_uv_scale_edge"
		)
		expect_s3_class(err$fit, "lame_dynamic_als")
		expect_equal(err$fit$family, fam)
		expect_true(is.finite(err$fit$eb$uv_scale))
		testthat::local_mocked_bindings(
			.lde_pos_profile = function(...) { pr = real(...); pr$unsettled = TRUE; pr$n_near = 5L; pr$n_unsettled = 2L; pr },
			.package = "lame"
		)
		err = expect_error(
			suppressMessages(
				lame:::lame_dynamic_als(
					d$Y,
					d$X,
					family = fam,
					dynamic_ab = TRUE,
					R = 1,
					dynamic_uv = TRUE,
					verbose = FALSE
				)
			),
			class = "lame_dynamic_uv_scale_unsettled"
		)
		expect_s3_class(err, "lame_als_convergence")
		expect_s3_class(err$fit, "lame_dynamic_als")
		# the count in the message is the search's own evaluations, not a recount at the refined scale
		expect_match(
			gsub("\\s+", " ", conditionMessage(err)),
			"did not settle at 2 of the 5 position scales nearest the minimum"
		)
	}
})

test_that("a constant normal response is refused: it has no variance to model", {
	d = rf_panel("normal")
	Y1 = lapply(d$Y, function(y) { y[!is.na(y)] = 2; y })
	err = expect_error(
		lame:::lame_dynamic_als(Y1, d$X, family = "normal", dynamic_ab = TRUE, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
	expect_match(gsub("\\s+", " ", conditionMessage(err)), "a constant normal response")
})

test_that("static factors of a normal panel are found at any unit of the response: the start is relative to the weights", {
	skip_on_cran()
	set.seed(7)
	n = 14
	u = matrix(rnorm(n), n, 1)
	v = matrix(rnorm(n), n, 1)
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	Y = lapply(1:3, function(t) {
		y = 0.2 + 1.2 * tcrossprod(u, v) + matrix(rnorm(n * n, 0, 0.5), n)
		dimnames(y) = list(rn, rn)
		diag(y) = NA
		y
	})
	fit = function(s) lame:::lame_dynamic_als(
		lapply(Y, function(y) y * s),
		NULL,
		family = "normal",
		dynamic_ab = TRUE,
		R = 1,
		verbose = FALSE
	)
	f1 = fit(1)
	# the product is found, not the null fixed point where the noise absorbs it (a risk for responses in large units,
	# where the weights 1 / ve are far below one)
	expect_gt(sqrt(mean(f1$UVPM^2, na.rm = TRUE)), 0.5)
	for (s in c(1e-6, 1e4, 1e8)) {
		g = fit(s)
		expect_true(g$converged)
		expect_equal(g$VC[["ve"]], s^2 * f1$VC[["ve"]], tolerance = 2e-2)
		# equal to the scaled product up to the factor prior, whose scale follows the response variance, so a larger
		# unit shrinks the product slightly less
		expect_equal(sqrt(mean(g$UVPM^2, na.rm = TRUE)), s * sqrt(mean(f1$UVPM^2, na.rm = TRUE)), tolerance = 6e-2)
	}
})

test_that("a normal response in any unit has the fit of the unit one: the floors of the variances are relative", {
	d = rf_panel("normal")
	f = lame:::lame_dynamic_als(d$Y, d$X, family = "normal", dynamic_ab = TRUE, verbose = FALSE)
	for (s in c(1e-8, 1e8)) {
		g = lame:::lame_dynamic_als(
			lapply(d$Y, function(y) y * s),
			d$X,
			family = "normal",
			dynamic_ab = TRUE,
			verbose = FALSE
		)
		expect_true(g$converged)
		expect_equal(g$VC[["ve"]], s^2 * f$VC[["ve"]], tolerance = 5e-3)
		expect_equal(g$sigma_ab, s * f$sigma_ab, tolerance = 5e-3)
		expect_equal(unname(g$coefficients[-1]), s * unname(f$coefficients[-1]), tolerance = 5e-3)
	}
})
