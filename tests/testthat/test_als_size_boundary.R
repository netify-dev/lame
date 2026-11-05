library(lame)
skip_on_cran()

# the size boundary of the fixed-effects binary ALS fit: a network the split-sample correction does not reach
# is fit uncorrected at or above the size floor and refused (or bootstrap-corrected) below it; an undirected fit
# with a latent rank takes this default, and bias_correct = "none" opts out at any size

# a directed binary panel with known effects (probit, density about 0.4)
sim_size_panel = function(
	n = 40,
	Tt = 1,
	R = 0,
	seed = 1,
	b0 = -0.3,
	sd_ab = 0.5,
	sd_uv = 0.6
) {
	set.seed(seed)
	nm = sprintf("a%03d", seq_len(n))
	a = rnorm(n, 0, sd_ab); b = rnorm(n, 0, sd_ab)
	U = matrix(rnorm(n * max(R, 1), 0, sd_uv), n)
	V = matrix(rnorm(n * max(R, 1), 0, sd_uv), n)
	P = if (R > 0) U[, seq_len(R), drop = FALSE] %*% t(V[, seq_len(R), drop = FALSE]) else matrix(0, n, n)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, n)
		z = b0 + 0.7 * x + outer(a, b, "+") + P + matrix(rnorm(n * n), n, n)
		y = 1 * (z > 0); diag(y) = NA; diag(x) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(nm, nm, "x1"))
	}
	list(Y = Y, X = X)
}

# the size rule, applied to a fixed-effects fit
size_rule = function(fit, max_iter = 200, tol = 1e-6) {
	lame:::.ae_jackknife(fit, max_iter, tol, explicit = FALSE)
}
fe_default = function(Y, X, ...) {
	size_rule(suppressMessages(ame_als(
		Y,
		X,
		family = "binary",
		verbose = FALSE,
		bias_correct = "none",
		...
	)))
}
lfe_default = function(Y, X, ...) {
	size_rule(suppressMessages(lame_als(
		Y,
		X,
		family = "binary",
		verbose = FALSE,
		bias_correct = "none",
		...
	)))
}

shown = function(expr) {
	msgs = character(0)
	out = utils::capture.output(withCallingHandlers(expr, message = function(m) {
		msgs <<- c(msgs, conditionMessage(m))
		invokeRestart("muffleMessage")
	}))
	paste(c(out, msgs), collapse = "\n")
}

test_that("the size of a binary panel is its cells per actor times the rarer-outcome share", {
	y = matrix(
		c(
			NA, 1, 0, 0,
			1, NA, 0, 1,
			0, 1, NA, 1,
			1, 0, 1, NA
		),
		4,
		4,
		byrow = TRUE,
		dimnames = list(letters[1:4], letters[1:4])
	)
	sz = lame:::.ae_binary_size(array(y, c(4, 4, 1)), FALSE, FALSE, letters[1:4], letters[1:4])
	# 3 observed cells per actor, 7 ties in 12 cells: the rarer outcome is
	# the non-tie, 5 in 12
	expect_equal(sz$cells, 3)
	expect_equal(sz$size, 3 * 5 / 12)
	expect_length(sz$boundary, 0L)
	# a sender without a tie is a boundary actor
	y2 = y; y2["a", ] = c(NA, 0, 0, 0)
	sz2 = lame:::.ae_binary_size(array(y2, c(4, 4, 1)), FALSE, FALSE, letters[1:4], letters[1:4])
	expect_identical(sz2$boundary, "a (sender)")
	# two periods double the cells per actor
	sz3 = lame:::.ae_binary_size(array(c(y, y), c(4, 4, 2)), FALSE, FALSE, letters[1:4], letters[1:4])
	expect_equal(sz3$size, 2 * sz$size)
})

test_that("the size rule corrects a fixed-effects fit below the floor by the parametric bootstrap", {
	d = sim_size_panel(n = 40)
	fit = fe_default(d$Y[[1]], d$X[[1]])
	f0 = ame_als(d$Y[[1]], d$X[[1]], family = "binary", verbose = FALSE, bias_correct = "none")
	bc = fit$bias_correct
	expect_identical(bc$method, "bootstrap")
	expect_true(bc$applied)
	expect_lt(bc$size, bc$floor)
	expect_match(bc$correction_problem, "rarer-outcome cells per actor")
	expect_true(all(bc$valid >= 24))
	expect_identical(fit$raw$coefficients, f0$coefficients)
	# the stretch is pulled in: coefficient and effect variances move down
	expect_lt(abs(fit$beta[["x1_dyad"]]), abs(f0$beta[["x1_dyad"]]))
	expect_lt(fit$VC[["va"]], f0$VC[["va"]])
	expect_lt(fit$VC[["vb"]], f0$VC[["vb"]])
	# the fit is reproducible and leaves the caller's random numbers alone
	set.seed(1); before = runif(1); set.seed(1)
	fit2 = fe_default(d$Y[[1]], d$X[[1]])
	expect_identical(runif(1), before)
	expect_identical(coef(fit2), coef(fit))
	# the linear predictor is rebuilt from the corrected parts
	ez = fit$mu + fit$beta[["x1_dyad"]] * d$X[[1]][, , 1] + outer(fit$a, fit$b, "+")
	diag(ez) = NA
	expect_equal(unname(fit$EZ[[1]]), unname(ez), tolerance = 1e-10)
	expect_match(shown(print(fit)), "parametric bootstrap")
	# the public default is the empirical-Bayes estimator, also through ame()
	fa = suppressWarnings(ame(
		d$Y[[1]],
		d$X[[1]],
		family = "binary",
		method = "als",
		verbose = FALSE
	))
	expect_identical(fa$bias_correct$method, "eb")
	# the parametric bootstrap on request gives the same correction
	fr = ame_als(
		d$Y[[1]],
		d$X[[1]],
		family = "binary",
		verbose = FALSE,
		bias_correct = "bootstrap"
	)
	expect_identical(fr$bias_correct$method, "bootstrap")
	expect_equal(coef(fr), coef(fit), tolerance = 1e-8)
	# an explicit request for the split-sample correction keeps its own error
	expect_error(
		ame_als(
			d$Y[[1]],
			d$X[[1]],
			family = "binary",
			R = 1,
			verbose = FALSE,
			bias_correct = "jackknife"
		),
		class = "lame_jk_split"
	)
	# a warm-start refit and a bootstrap replicate remove the stored correction
	rf = ame_als_refit(fit)
	expect_identical(rf$bias_correct$method, "bootstrap")
	expect_true(rf$bias_correct$stored)
	expect_equal(coef(rf), coef(fit), tolerance = 1e-3)
	bt = ame_als_bootstrap(fit, R = 10, seed = 1, verbose = FALSE)
	expect_equal(bt$n_valid, 10)
	expect_lt(
		abs(mean(bt$coefs[, "x1_dyad"]) - coef(fit)[["x1_dyad"]]),
		abs(mean(bt$coefs[, "x1_dyad"]) - coef(f0)[["x1_dyad"]])
	)
	# the correction can be asked for where another one is the default
	big = sim_size_panel(n = 50, Tt = 2, seed = 4)
	expect_identical(lfe_default(big$Y, big$X)$bias_correct$method, "jackknife")
	fb = lame_als(big$Y, big$X, family = "binary", verbose = FALSE, bias_correct = "bootstrap")
	expect_identical(fb$bias_correct$method, "bootstrap")
	expect_error(
		ame_als(
			big$Y[[1]] + 0.5,
			big$X[[1]],
			family = "normal",
			verbose = FALSE,
			bias_correct = "bootstrap"
		),
		class = "lame_jk_split"
	)
})

test_that("a short panel with a multiplicative term below the floor is corrected", {
	p = sim_size_panel(n = 30, Tt = 3, R = 1, seed = 2)
	fit = suppressWarnings(lfe_default(p$Y, p$X, R = 1))
	expect_identical(fit$bias_correct$method, "bootstrap")
	expect_match(fit$bias_correct$correction_problem, "periods")
	expect_equal(unname(fit$UVPM), unname(fit$U %*% t(fit[["V"]])), tolerance = 1e-10)
	expect_error(
		lame(
			p$Y,
			p$X,
			family = "binary",
			R = 1,
			method = "als",
			dynamic_ab = TRUE,
			bias_correct = "bootstrap",
			verbose = FALSE
		),
		class = "lame_jk_split"
	)
})

test_that("bias_correct = 'none' runs below the floor and records the request", {
	d = sim_size_panel(n = 40)
	expect_no_warning(
		fit <- ame_als(
			d$Y[[1]],
			d$X[[1]],
			family = "binary",
			verbose = FALSE,
			bias_correct = "none"
		)
	)
	expect_identical(fit$bias_correct$method, "none")
	expect_false(fit$bias_correct$applied)
	expect_identical(fit$bias_correct$requested, "none")
	expect_match(shown(print(fit)), "bias_correct = \"none\"", fixed = TRUE)
	# the normal family is not corrected and its default records no request
	fn = ame_als(d$Y[[1]] + rnorm(1600), d$X[[1]], family = "normal", verbose = FALSE)
	expect_null(fn$bias_correct$requested)
})

test_that("a network above the floor that the correction does not reach is fit uncorrected without a flag", {
	# three periods at R = 1: the correction needs five, the size rule and the
	# split-sample check decide
	p = sim_size_panel(n = 150, Tt = 3, R = 1, seed = 3)
	expect_no_warning(fit <- lfe_default(p$Y, p$X, R = 1))
	bc = fit$bias_correct
	expect_identical(bc$method, "none")
	expect_false(bc$applied)
	expect_null(bc$requested)
	expect_gte(bc$size, bc$floor)
	expect_equal(bc$floor, lame:::.ae_size_floor(1L, bc$spread))
	expect_match(bc$basis, "rarer-outcome cells per actor")
	expect_match(bc$basis, "fitted additive-effect spread")
	expect_match(bc$correction_problem, "periods")
	expect_null(fit$raw)
	out = shown(print(fit))
	expect_match(out, "none needed at this size")
	expect_false(grepl("not applied", out))
	expect_match(shown(print(summary(fit))), "none needed at this size")
	# the uncorrected fit is the opt-out fit
	f0 = lame_als(p$Y, p$X, family = "binary", R = 1, verbose = FALSE, bias_correct = "none")
	expect_identical(coef(fit), coef(f0))
	expect_identical(fit$VC, f0$VC)
	# a warm-start refit keeps the record
	rf = suppressWarnings(ame_als_refit(fit))
	expect_identical(rf$bias_correct, fit$bias_correct)
})

test_that("a network whose fitted effect spread is above the measured band is refused", {
	# very heterogeneous sender and receiver effects: no size brings the bias of the uncorrected
	# estimate within the tolerance, so the fit is refused whatever its size
	p = sim_size_panel(n = 150, Tt = 3, R = 1, seed = 1, sd_ab = 3)
	e = tryCatch(lfe_default(p$Y, p$X, R = 1), error = function(e) e)
	expect_s3_class(e, "lame_als_size")
	expect_match(conditionMessage(e), "fitted additive-effect spread")
	expect_match(conditionMessage(e), "ame()", fixed = TRUE)
	f0 = suppressMessages(lame_als(
		p$Y,
		p$X,
		family = "binary",
		R = 1,
		verbose = FALSE,
		bias_correct = "none"
	))
	expect_gt(lame:::.ae_fit_spread(f0), max(lame:::.ae_size_basis$spread))
	expect_identical(lame:::.ae_size_floor(1L, lame:::.ae_fit_spread(f0)), Inf)
})

test_that("in a fixed-effects fit an actor without a tie has its effect on that side not estimated", {
	p = sim_size_panel(n = 150, Tt = 3, R = 1, seed = 3)
	for (t in seq_along(p$Y)) p$Y[[t]]["a007", ] = c(rep(0, 6), NA, rep(0, 143))
	expect_message(
		fit <- lame_als(
			p$Y,
			p$X,
			family = "binary",
			R = 1,
			verbose = FALSE,
			bias_correct = "none"
		),
		"a007 \\(sender, no tie\\)"
	)
	expect_true(is.na(fit$a[["a007"]]))
	expect_true(is.finite(fit$b[["a007"]]))
	expect_identical(fit$not_estimated$actor, "a007")
	expect_identical(fit$not_estimated$side, "sender")
	# its sender cells are left out of the fit and of the fitted values
	expect_true(all(is.na(fit$Y["a007", , ])))
	expect_true(all(is.na(fit$EZ[[1]]["a007", ])))
	expect_true(all(is.finite(fit$EZ[[1]]["a008", -8])))
	# va describes the estimated effects and stays near that of the network
	# without the change
	p0 = sim_size_panel(n = 150, Tt = 3, R = 1, seed = 3)
	f0 = lame_als(p0$Y, p0$X, family = "binary", R = 1, verbose = FALSE, bias_correct = "none")
	expect_equal(fit$VC[["va"]], f0$VC[["va"]], tolerance = 0.05)
	expect_equal(fit$beta[["x1_dyad"]], f0$beta[["x1_dyad"]], tolerance = 0.02)
	# a warm-start refit keeps the actor out
	rf = suppressWarnings(ame_als_refit(fit))
	expect_true(is.na(rf$a[["a007"]]))
})

test_that("binary IRLS reaches the fixed-effects maximum-likelihood fit", {
	d = sim_size_panel(n = 40, seed = 11, b0 = -1.2)
	y = d$Y[[1]]; x = d$X[[1]][, , 1]
	for (link in c("probit", "logit")) {
		f = suppressMessages(ame_als(
			y,
			d$X[[1]],
			family = "binary",
			link = link,
			verbose = FALSE,
			bias_correct = "none"
		))
		keep = !is.na(as.vector(f$Y[, , 1]))
		df = data.frame(
			y = as.vector(y),
			x = as.vector(x),
			s = factor(as.vector(row(y))),
			r = factor(as.vector(col(y)))
		)[keep, ]
		g = suppressWarnings(glm(
			y ~ x + s + r,
			family = binomial(link),
			data = df,
			control = glm.control(epsilon = 1e-12, maxit = 100)
		))
		expect_equal(f$beta[["x1_dyad"]], unname(coef(g)[["x"]]), tolerance = 1e-4)
		expect_equal(f$deviance, deviance(g), tolerance = 1e-6)
	}
})

test_that("a single cross-section at R = 0 is corrected by the checkerboard split from its floor", {
	d = sim_size_panel(n = 120, seed = 4)
	fit = ame_als(
		d$Y[[1]],
		d$X[[1]],
		family = "binary",
		verbose = FALSE,
		bias_correct = "jackknife"
	)
	bc = fit$bias_correct
	expect_identical(bc$method, "jackknife")
	expect_true(bc$applied)
	expect_identical(bc$split, "checkerboard")
	expect_equal(unname(bc$weights), c(2, -0.5))
	f0 = ame_als(d$Y[[1]], d$X[[1]], family = "binary", verbose = FALSE, bias_correct = "none")
	expect_identical(fit$raw$coefficients, f0$coefficients)
	expect_identical(fit$raw$VC, f0$VC)
	# the stretch is pulled in
	expect_lt(fit$coefficients[["x1_dyad"]], f0$coefficients[["x1_dyad"]])
	expect_lt(fit$VC[["va"]], f0$VC[["va"]])
	# the linear predictor is rebuilt from the corrected parts
	ez = fit$mu + fit$beta[["x1_dyad"]] * d$X[[1]][, , 1] + outer(fit$a, fit$b, "+")
	expect_equal(unname(fit$EZ[[1]]), unname(ez), tolerance = 1e-10)
	expect_equal(sum(fit$a), 0, tolerance = 1e-8)
	expect_match(shown(print(fit)), "checkerboard split")
	# on normal data the combination leaves the estimate where it was
	set.seed(5)
	zn = -0.3 + 0.7 * d$X[[1]][, , 1] + outer(rnorm(120, 0, 0.5), rnorm(120, 0, 0.5), "+") +
		matrix(rnorm(120 * 120), 120)
	fn = ame_als(zn, d$X[[1]], family = "normal", verbose = FALSE)
	fj = ame_als(zn, d$X[[1]], family = "normal", verbose = FALSE, bias_correct = "jackknife")
	expect_identical(fj$bias_correct$split, "checkerboard")
	expect_equal(fj$coefficients[["x1_dyad"]], fn$coefficients[["x1_dyad"]], tolerance = 0.01)
	expect_equal(fj$VC[["va"]], fn$VC[["va"]], tolerance = 0.03)
	# below the checkerboard floor the cross-section is corrected by the
	# parametric bootstrap
	d45 = sim_size_panel(n = 45, seed = 4)
	f45 = fe_default(d45$Y[[1]], d45$X[[1]])
	expect_identical(f45$bias_correct$method, "bootstrap")
	expect_match(f45$bias_correct$correction_problem, "checkerboard")
})

test_that("the bootstrap correction refuses an undirected fit with a multiplicative term", {
	# undirected networks with a rank-1 term: the fitted product is flattened (here removed
	# by the factor penalty) and the correction would move the slope away from the truth
	set.seed(150101)
	n = 50
	a = rnorm(n, 0, 0.6); u = rnorm(n, 0, sqrt(0.6))
	x = matrix(rnorm(n * n), n); x[lower.tri(x)] = t(x)[lower.tri(x)]
	e = matrix(rnorm(n * n), n); e[lower.tri(e)] = t(e)[lower.tri(e)]
	y = 1 * (-1 + 0.7 * x + outer(a, a, "+") + outer(u, u) + e > 0); diag(y) = NA
	e = tryCatch(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			R = 1,
			symmetric = TRUE,
			verbose = FALSE,
			bias_correct = "bootstrap"
		)),
		error = function(e) e
	)
	expect_s3_class(e, "lame_als_size")
	expect_match(conditionMessage(e), "undirected")
	expect_match(conditionMessage(e), "R = 0", fixed = TRUE)
	# the opt-out runs, and R = 0 is corrected
	expect_identical(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			R = 1,
			symmetric = TRUE,
			verbose = FALSE,
			bias_correct = "none"
		))$bias_correct$method,
		"none"
	)
	expect_identical(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			symmetric = TRUE,
			verbose = FALSE,
			bias_correct = "bootstrap"
		))$bias_correct$method,
		"bootstrap"
	)
	# the default, with or without a latent rank, is the empirical-Bayes estimator
	expect_identical(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			symmetric = TRUE,
			verbose = FALSE
		))$bias_correct$method,
		"eb"
	)
	expect_identical(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			R = 1,
			symmetric = TRUE,
			verbose = FALSE
		))$bias_correct$method,
		"eb"
	)
	expect_identical(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			R = 1,
			symmetric = TRUE,
			verbose = FALSE,
			bias_correct = "eb"
		))$bias_correct$method,
		"eb"
	)
})
