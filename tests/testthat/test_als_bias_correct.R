library(lame)
skip_on_cran()

# split-sample jackknife correction of the static binary ALS estimator

# a bipartite binary panel with known effects; the identified additive effects
# carry the row / column means of the product, as the estimator reports them
sim_bc_panel = function(
	nA = 24,
	nB = 16,
	Tt = 5,
	R = 1,
	seed = 1,
	family = "binary",
	b0 = -0.5,
	b1 = 0.7,
	sd_ab = 0.5,
	sd_uv = 0.6
) {
	set.seed(seed)
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	a = rnorm(nA, 0, sd_ab); b = rnorm(nB, 0, sd_ab)
	U = matrix(rnorm(nA * R, 0, sd_uv), nA, R)
	V = matrix(rnorm(nB * R, 0, sd_uv), nB, R)
	P = if (R > 0) U %*% t(V) else matrix(0, nA, nB)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(nA * nB), nA, nB)
		z = b0 + b1 * x + outer(a, b, "+") + P + matrix(rnorm(nA * nB), nA, nB)
		y = if (family == "binary") 1 * (z > 0) else z
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(nA, nB, 1), dimnames = list(rn, cn, "x1"))
	}
	list(
		Y = Y,
		X = X,
		b1 = b1,
		a_ident = a + rowMeans(P),
		b_ident = b + colMeans(P)
	)
}

bc_fit = function(d, ..., R = 1, family = "binary") {
	suppressMessages(lame_als(
		d$Y,
		d$X,
		R = R,
		family = family,
		mode = "bipartite",
		verbose = FALSE,
		...
	))
}

# the size rule that sets the default of an undirected binary fit with a latent rank,
# applied to a fixed-effects fit
size_rule = function(fit, max_iter = 200, tol = 1e-6) {
	lame:::.ae_jackknife(fit, max_iter, tol, explicit = FALSE)
}

test_that("the empirical-Bayes estimator is the default for binary and normal fits, the jackknife is explicit, bias_correct = 'none' is uncorrected", {
	d = sim_bc_panel(seed = 3)
	expect_identical(bc_fit(d)$bias_correct$method, "eb")
	fit = bc_fit(d, bias_correct = "jackknife")
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_false(identical(coef(fit), fit$raw$coefficients))

	dn = sim_bc_panel(seed = 3, family = "normal")
	fn = bc_fit(dn, family = "normal")
	fn0 = bc_fit(dn, family = "normal", bias_correct = "none")
	expect_identical(fn$bias_correct$method, "eb")
	expect_identical(fn0$bias_correct$method, "none")
	expect_null(fn$raw); expect_null(fn0$raw)
	expect_false(identical(coef(fn), coef(fn0)))
	expect_equal(unname(coef(fn)), unname(coef(fn0)), tolerance = 0.05)
})

test_that("bias_correct = 'none' reproduces the uncorrected fit", {
	d = sim_bc_panel(seed = 3)
	fit = bc_fit(d, bias_correct = "jackknife")
	fit0 = bc_fit(d, bias_correct = "none")
	expect_identical(fit0$bias_correct$method, "none")
	expect_null(fit0$raw)
	expect_identical(fit0$coefficients, fit$raw$coefficients)
	expect_identical(fit0$a, fit$raw$a)
	expect_identical(fit0$b, fit$raw$b)
	expect_identical(fit0$VC, fit$raw$VC)
	# the product surface the uncorrected fit's linear predictor carries
	O0 = fit0$UVPM
	O0 = O0 - outer(rowMeans(O0), rep(1, ncol(O0))) -
		outer(rep(1, nrow(O0)), colMeans(O0)) + mean(O0)
	expect_equal(unname(fit$raw$O), unname(O0), tolerance = 1e-12)
})

test_that("the correction moves an inflated small-panel estimate toward the truth", {
	d = sim_bc_panel(nA = 30, nB = 20, Tt = 4, R = 0, seed = 1)
	fit = bc_fit(d, R = 0, bias_correct = "jackknife")
	raw = fit$raw
	b1_raw = unname(raw$coefficients["x1_dyad"])
	b1_jk = unname(fit$coefficients["x1_dyad"])
	# the fixed-effects probit stretches the latent scale outward; the
	# jackknife pulls the slope and the effect spread back toward the truth
	expect_gt(b1_raw, d$b1)
	expect_lt(abs(b1_jk - d$b1), abs(b1_raw - d$b1))
	slope = function(est, truth) unname(coef(lm(est ~ truth))[2])
	expect_lt(abs(slope(fit$a, d$a_ident) - 1), abs(slope(raw$a, d$a_ident) - 1))
	va_true = var(d$a_ident)
	expect_lt(abs(fit$VC[["va"]] - va_true), abs(raw$VC[["va"]] - va_true))
})

test_that("the corrected fit is internally consistent and records the correction", {
	d = sim_bc_panel(seed = 5)
	fit = bc_fit(d, bias_correct = "jackknife")
	bc = fit$bias_correct
	expect_setequal(
		names(bc),
		c(
			"method", "applied", "split", "weights", "tried",
			"converged", "iterations", "irls_iterations", "irls_diagnostics",
			"shrink"
		)
	)
	expect_true(bc$applied)
	expect_identical(bc$split, "checkerboard")
	expect_equal(unname(bc$weights), c(2, -0.5))
	# the correction pulls the additive effects in
	expect_lt(bc$shrink[["a"]], 1)
	expect_lt(bc$shrink[["b"]], 1)
	expect_length(bc$converged, 2L)
	expect_type(bc$converged, "logical")
	expect_length(bc$iterations, 2L)
	expect_setequal(names(fit$raw), c("coefficients", "a", "b", "VC", "O", "deviance"))
	# the linear predictor is rebuilt from the corrected components
	for (t in seq_along(d$Y)) {
		ez = fit$mu + fit$beta[["x1_dyad"]] * d$X[[t]][, , 1] +
			outer(fit$a, fit$b, "+") + fit$U %*% t(fit$V)
		expect_equal(unname(fit$EZ[[t]]), unname(ez), tolerance = 1e-10)
		expect_equal(
			unname(fit$fitted[[t]]),
			pnorm(unname(fit$EZ[[t]])),
			tolerance = 1e-12
		)
	}
	expect_equal(unname(fit$UVPM), unname(fit$U %*% t(fit$V)), tolerance = 1e-12)
	r = residuals(fit)
	expect_equal(unname(r[[1]]), unname(d$Y[[1]] - fit$fitted[[1]]), tolerance = 1e-12)
	# the additive effects keep the sum-to-zero gauge and the variances stay nonnegative
	expect_equal(sum(fit$a), 0, tolerance = 1e-8)
	expect_equal(sum(fit$b), 0, tolerance = 1e-8)
	expect_gte(fit$VC[["va"]], 0)
	expect_gte(fit$VC[["vb"]], 0)
})

test_that("a bootstrap on a corrected fit refits each replicate with the correction", {
	d = sim_bc_panel(nA = 24, nB = 16, Tt = 4, R = 0, seed = 3, b0 = 0)
	fit = bc_fit(d, R = 0, bias_correct = "jackknife")
	refit = ame_als_refit(fit)
	expect_identical(refit$bias_correct$method, "jackknife")
	expect_identical(ame_als_refit(fit, bias_correct = "none")$bias_correct$method, "none")
	bt = ame_als_bootstrap(fit, R = 10, seed = 1, verbose = FALSE)
	expect_s3_class(bt, "boot_ame")
	expect_equal(bt$n_valid + bt$n_split_dropped, 10)
	expect_gte(bt$n_valid, 8)
	expect_true(all(is.finite(bt$se)))
	# the replicates scatter around the corrected estimate, not the raw one
	b1_boot = median(bt$coefs[, "x1_dyad"], na.rm = TRUE)
	expect_lt(
		abs(b1_boot - coef(fit)[["x1_dyad"]]),
		abs(b1_boot - fit$raw$coefficients[["x1_dyad"]])
	)
})

test_that("a symmetric panel is corrected and a small single cross-section is bootstrap-corrected", {
	set.seed(4)
	n = 40
	x = matrix(rnorm(n * n), n, n)
	a = rnorm(n, 0, 0.3)
	z = 0.7 * x + outer(a, a, "+") + matrix(rnorm(n * n), n, n)
	y = 1 * (z > 0); diag(y) = NA
	# a 40-actor cross-section is below the size where the checkerboard correction holds and below the size floor
	# of the uncorrected fit: the size rule corrects it by the parametric bootstrap, an explicit split-sample
	# request is refused, and the opt-out runs
	f_none = suppressMessages(ame_als(
		y,
		x,
		family = "binary",
		verbose = FALSE,
		bias_correct = "none"
	))
	fb = size_rule(f_none)
	expect_identical(fb$bias_correct$method, "bootstrap")
	expect_error(
		suppressMessages(ame_als(
			y,
			x,
			family = "binary",
			verbose = FALSE,
			bias_correct = "jackknife"
		)),
		class = "lame_jk_split"
	)
	fit = suppressMessages(ame_als(
		y,
		x,
		family = "binary",
		verbose = FALSE,
		bias_correct = "none"
	))
	expect_identical(fit$bias_correct$method, "none")

	sym_slice = function(seed) {
		set.seed(seed)
		zt = 0.7 * x + outer(a, a, "+") + matrix(rnorm(n * n), n, n)
		zt = (zt + t(zt)) / 2
		yt = 1 * (zt > 0); diag(yt) = NA
		yt
	}
	xs = (x + t(x)) / 2
	fs = suppressMessages(lame_als(
		lapply(101:105, sym_slice),
		rep(list(xs), 5),
		R = 1,
		family = "binary",
		symmetric = TRUE,
		verbose = FALSE,
		bias_correct = "jackknife"
	))
	expect_identical(fs$bias_correct$method, "jackknife")
	# the default of an undirected panel with a latent rank is the
	# empirical-Bayes estimator
	fe = suppressMessages(lame_als(
		lapply(101:105, sym_slice),
		rep(list(xs), 5),
		R = 1,
		family = "binary",
		symmetric = TRUE,
		verbose = FALSE
	))
	expect_identical(fe$bias_correct$method, "eb")
	expect_identical(fs$a, fs$b)
	expect_equal(unname(fs$ULUPM), unname(fs$U %*% (fs$L * t(fs$U))), tolerance = 1e-12)
	expect_equal(unname(fs$EZ[[1]]), t(unname(fs$EZ[[1]])), tolerance = 1e-10)
})

test_that("the correction refuses a panel it cannot split", {
	d = sim_bc_panel(seed = 3)
	expect_error(bc_fit(d, bias_correct = "yes"), "bias_correct")
	expect_error(
		suppressMessages(lame(
			d$Y,
			d$X,
			family = "binary",
			mode = "bipartite",
			method = "als",
			dynamic_ab = TRUE,
			bias_correct = "jackknife",
			verbose = FALSE
		)),
		"not available"
	)
	# an actor with a single tie in the panel keeps none in one half, which
	# cannot estimate its effect
	for (t in seq_along(d$Y)) d$Y[[t]]["r01", ] = 0
	d$Y[[1]]["r01", "c01"] = 1
	expect_error(bc_fit(d, bias_correct = "jackknife"), class = "lame_jk_split")
	# the parametric bootstrap corrects the panel the splits cannot carry, on
	# request and by the size rule, which records the actor the splits lose
	fd = bc_fit(d, bias_correct = "bootstrap")
	expect_identical(fd$bias_correct$method, "bootstrap")
	fd_rule = size_rule(bc_fit(d, bias_correct = "none"))
	expect_identical(fd_rule$bias_correct$method, "bootstrap")
	expect_match(fd_rule$bias_correct$correction_problem, "r01")
	# every admissible split was tried before refusing
	sp = lame:::.ae_jk_split(simplify2array(d$Y), "binary", FALSE, TRUE, 1L)
	expect_setequal(sp$tried$split, c("checkerboard", "random", "periods"))
	# a stored fit that asked for the correction and could not compute it still says so when printed
	fit = bc_fit(d, bias_correct = "none")
	fit$bias_correct = list(
		method = "none", applied = FALSE,
		requested = "jackknife", reason = "r01 without a tie in a half"
	)
	shown = function(expr) {
		msgs = character(0)
		out = utils::capture.output(withCallingHandlers(expr, message = function(m) {
			msgs <<- c(msgs, conditionMessage(m))
			invokeRestart("muffleMessage")
		}))
		paste(c(out, msgs), collapse = "\n")
	}
	expect_match(shown(print(fit)), "not applied")
	expect_match(shown(print(summary(fit))), "not applied")
	# an actor with no tie at all is left out on that side (its effect is at
	# the boundary of the likelihood), and the correction runs on the rest
	for (t in seq_along(d$Y)) d$Y[[t]]["r01", ] = 0
	fj = bc_fit(d, bias_correct = "jackknife")
	expect_identical(fj$bias_correct$method, "jackknife")
	expect_true(is.na(fj$a[["r01"]]))
	expect_identical(fj$not_estimated$actor, "r01")
	# a multiplicative term needs five or more periods for the split; this
	# panel is below the size floor of the uncorrected fit and the size rule
	# corrects it by the parametric bootstrap
	d3 = sim_bc_panel(seed = 3, Tt = 3)
	expect_error(bc_fit(d3, bias_correct = "jackknife"), class = "lame_jk_split")
	f3 = size_rule(bc_fit(d3, bias_correct = "none"))
	expect_identical(f3$bias_correct$method, "bootstrap")
	expect_match(f3$bias_correct$correction_problem, "five|periods")
	# at R = 0 three periods are enough from 25 rarer-outcome cells per actor
	d3b = sim_bc_panel(nA = 40, nB = 30, Tt = 3, seed = 3)
	expect_identical(bc_fit(d3b, R = 0, bias_correct = "jackknife")$bias_correct$method, "jackknife")
	expect_identical(size_rule(bc_fit(d3, R = 0, bias_correct = "none"))$bias_correct$method, "bootstrap")
})

test_that("another admissible split is used when the checkerboard fails", {
	d = sim_bc_panel(seed = 3)
	# two ties of r01, both in cells the checkerboard puts in the same half
	for (t in seq_along(d$Y)) d$Y[[t]]["r01", ] = 0
	d$Y[[1]]["r01", "c01"] = 1
	d$Y[[2]]["r01", "c02"] = 1
	set.seed(99)
	seed_before = .Random.seed
	fit = bc_fit(d, bias_correct = "jackknife")
	# the random split draws from its own seed and leaves the caller's stream
	expect_identical(.Random.seed, seed_before)
	bc = fit$bias_correct
	expect_true(bc$applied)
	expect_true(bc$split %in% c("random", "periods"))
	expect_gt(bc$tried$no_tie[bc$tried$split == "checkerboard"], 0)
	expect_equal(sum(bc$weights[["full"]], 2 * bc$weights[["half"]]), 1)
	# the same split on a refit
	expect_identical(bc_fit(d, bias_correct = "jackknife")$coefficients, fit$coefficients)
})

test_that("an uncorrected fit reports the surface its linear predictor carries", {
	for (sd in 1:3) {
		d = sim_bc_panel(seed = sd, R = 2)
		fit = bc_fit(d, R = 2, bias_correct = "none")
		O = fit$U %*% t(fit$V)
		expect_equal(unname(fit$UVPM), unname(O), tolerance = 1e-10)
		expect_equal(unname(rowMeans(fit$UVPM)), rep(0, nrow(O)), tolerance = 1e-10)
		expect_equal(unname(colMeans(fit$UVPM)), rep(0, ncol(O)), tolerance = 1e-10)
		for (t in seq_along(d$Y)) {
			ez = fit$mu + fit$beta[["x1_dyad"]] * d$X[[t]][, , 1] +
				outer(fit$a, fit$b, "+") + O
			expect_equal(unname(fit$EZ[[t]]), unname(ez), tolerance = 1e-10)
		}
	}
})

test_that("fit$deviance is the deviance of the reported EZ, the optimiser's value fit$objective", {
	dev_ez = function(fit) {
		d = 0
		for (t in seq_along(fit$EZ)) {
			y = fit$Y[, , t]; ez = fit$EZ[[t]]
			ok = is.finite(y) & is.finite(ez)
			if (identical(fit$family, "binary")) {
				p = pmin(pmax(pnorm(ez[ok]), 1e-12), 1 - 1e-12)
				d = d - 2 * sum(y[ok] * log(p) + (1 - y[ok]) * log(1 - p))
			} else {
				d = d + sum((y[ok] - ez[ok])^2)
			}
		}
		d
	}
	d = sim_bc_panel(seed = 3)
	# the split-sample correction rebuilds EZ from corrected components
	fit = bc_fit(d, bias_correct = "jackknife")
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_equal(fit$deviance, dev_ez(fit), tolerance = 1e-10)
	expect_true(is.finite(fit$raw$deviance))
	# uncorrected: the deviance of EZ, and the penalised objective the IRLS
	# loop minimised, which adds the factor penalty
	f0 = bc_fit(d, bias_correct = "none")
	expect_equal(f0$deviance, dev_ez(f0), tolerance = 1e-10)
	expect_true(is.finite(f0$objective))
	expect_false(isTRUE(all.equal(f0$objective, f0$deviance)))
	expect_equal(fit$raw$deviance, f0$deviance, tolerance = 1e-12)
	# the normal family: the residual sum of squares of EZ
	dn = sim_bc_panel(seed = 3, family = "normal")
	fn = bc_fit(dn, family = "normal")
	expect_equal(fn$deviance, dev_ez(fn), tolerance = 1e-10)
	# a warm-start refit reports the deviance of its own EZ
	rf = suppressWarnings(suppressMessages(ame_als_refit(f0)))
	expect_equal(rf$deviance, dev_ez(rf), tolerance = 1e-10)
	# and the summary prints it
	expect_match(paste(capture.output(print(summary(f0))), collapse = " "), "deviance =")
})

test_that("a bootstrap keeps replicates of a fit with a not-estimated actor", {
	# a sender without a tie is left out of a binary fit; its effect is NA in
	# the fit and in every replicate, and the other effects are summarised
	d = sim_bc_panel(nA = 24, nB = 16, Tt = 4, R = 0, seed = 3, b0 = 0)
	for (t in seq_along(d$Y)) d$Y[[t]]["r01", ] = 0
	fit = bc_fit(d, R = 0, bias_correct = "jackknife")
	expect_true(is.na(fit$a[["r01"]]))
	bt = ame_als_bootstrap(fit, R = 10, seed = 1, verbose = FALSE)
	expect_gte(bt$n_valid, 8)
	expect_true(all(is.na(bt$a_coefs[, 1])))
	expect_true(is.na(bt$se_a[1]))
	expect_true(all(is.finite(bt$se_a[-1]) | is.infinite(bt$se_a[-1])))
	ci = confint(bt, which = "a")
	expect_true(all(is.na(ci["a_r01", ])))
	expect_true(all(!is.na(ci[-1, ])))
	# the one-call route carries the correction through the replicates
	fit_b = suppressWarnings(bc_fit(
		d,
		R = 0,
		bias_correct = "jackknife",
		bootstrap = 30,
		bootstrap_seed = 2
	))
	expect_s3_class(fit_b$bootstrap, "boot_ame")
	expect_gte(fit_b$bootstrap$n_valid, 28)
	# the uncorrected fit has no bootstrap, asked for in the call or afterwards
	expect_error(
		bc_fit(d, R = 0, bias_correct = "none", bootstrap = 30),
		class = "lame_als_bootstrap_bias"
	)
	expect_error(
		ame_als_bootstrap(bc_fit(d, R = 0, bias_correct = "none"), R = 10, verbose = FALSE),
		class = "lame_als_bootstrap_bias"
	)
})

test_that("an actor at the boundary in a replicate records an infinite effect", {
	obj = list(a = c(x = NA, y = 0.1, z = -0.2), b = c(x = 0.3, y = 0.1, z = 0.2))
	fb = list(
		symmetric = FALSE,
		row_names = c("x", "y", "z"),
		col_names = c("x", "y", "z"),
		a = c(x = NA, y = NA, z = -0.1),
		b = c(x = NA, y = 0.2, z = 0.1),
		not_estimated = data.frame(
			actor = c("x", "y", "x"),
			side = c("sender", "sender", "receiver"),
			outcome = c("no tie", "no tie", "only ties"),
			stringsAsFactors = FALSE
		)
	)
	e = lame:::.ae_boot_effects(fb, obj)
	expect_true(e$ok)
	expect_identical(unname(e$a), c(NA, -Inf, -0.1))
	expect_identical(unname(e$b), c(Inf, 0.2, 0.1))
	# a non-finite effect the refit does not list is a failed replicate
	fb$a[["z"]] = NA_real_
	expect_false(lame:::.ae_boot_effects(fb, obj)$ok)
	M = rbind(c(NA, -Inf, 1), c(NA, 0.5, 2), c(NA, 0.2, 3))
	expect_identical(lame:::.ae_boot_col_sd(M, rep(TRUE, 3)), c(NA, Inf, 1))
})

test_that("a bootstrap-corrected binary ALS fit reports its correction at default verbosity", {
	skip_on_cran()
	set.seed(3301); n = 30
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Y = 1 * (-1 + outer(a, b, "+") + matrix(rnorm(n * n), n) > 0); diag(Y) = NA
	msgs = character()
	fit = withCallingHandlers(
		ame_als(Y, family = "binary", R = 0, bias_correct = "bootstrap"),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_identical(fit$bias_correct$method, "bootstrap")
	expect_true(any(grepl("two rounds of [0-9]+ simulated refits", msgs)))
})
