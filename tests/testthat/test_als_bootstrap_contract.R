library(lame)
skip_on_cran()

# the ALS bootstrap contract: a minimum number of replicates, no unconverged replicate in an interval, a capped
# share of dropped replicates, no replicate at another local optimum, a point estimate inside its own replicates,
# and no column-wise standard error for factor columns whose singular values are not separated

boot_fixture = function(n = 14, Tt = 3, R = 0, seed = 813114L, family = "normal") {
	set.seed(seed)
	lapply(seq_len(Tt), function(t) {
		z = matrix(rnorm(n * n), n)
		if (R > 0) z = z + tcrossprod(matrix(rnorm(n * R, 0, 1.2), n, R), matrix(rnorm(n * R, 0, 1.2), n, R))
		diag(z) = NA_real_
		z
	})
}

# the minimum number of replicates

test_that("bootstrap = N below 30 is refused before anything is fitted", {
	y = boot_fixture()
	for (n in c(1L, 10L, 29L)) {
		err = tryCatch(
			lame_als(y, R = 0, family = "normal", bootstrap = n, verbose = FALSE),
			error = identity
		)
		expect_s3_class(err, "lame_als_bootstrap_replicates")
		expect_match(conditionMessage(err), "minimum of 30", fixed = TRUE)
	}
	y1 = y[[1]]
	expect_error(
		ame_als(y1, R = 0, family = "normal", bootstrap = 12L, verbose = FALSE),
		class = "lame_als_bootstrap_replicates"
	)
	# the dispatchers inherit it
	expect_error(
		lame(y, R = 0, family = "normal", method = "als", bootstrap = 12L, verbose = FALSE),
		class = "lame_als_bootstrap_replicates"
	)
	# 30 is accepted and bootstrap = 0 is no bootstrap
	fit = lame_als(y, R = 0, family = "normal", bootstrap = 30L, bootstrap_seed = 1, verbose = FALSE)
	expect_s3_class(fit$bootstrap, "boot_ame")
	expect_equal(fit$bootstrap$n_total, 30L)
	expect_null(lame_als(y, R = 0, family = "normal", bootstrap = 0L, verbose = FALSE)$bootstrap)
})

test_that("ame_als_bootstrap() needs ten valid replicates", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	for (n in c(1L, 5L, 9L)) {
		expect_error(
			ame_als_bootstrap(f, R = n, verbose = FALSE),
			class = "lame_als_bootstrap_replicates"
		)
	}
	expect_error(ame_als_bootstrap(f, R = 0, verbose = FALSE), "positive whole number")
	expect_error(ame_als_bootstrap(f, R = 12.5, verbose = FALSE), "positive whole number")
	b = ame_als_bootstrap(f, R = 10, seed = 1, verbose = FALSE)
	expect_equal(b$n_valid, 10L)
})

# replicates that did not converge

test_that("replicates that did not converge are dropped, counted, and kept out of the intervals", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	drop = c(4L, 17L)
	k = 0L
	real = lame:::.ae_boot_refit
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
			k <<- k + 1L
			if (k %in% drop) {
				return(list(
					fit = NULL,
					attempts = 3L,
					iterations = 4L * max_iter,
					start = NA_character_,
					error = NULL,
					last_fit = object
				))
			}
			real(object, Y_b, X_b, Z_b, dyn, max_iter, tol)
		},
		.package = "lame"
	)
	# 2 of 40 is within the 5 percent cap
	b = ame_als_bootstrap(f, R = 40, seed = 3, verbose = FALSE)
	expect_equal(b$n_valid, 38L)
	expect_equal(b$n_failed, 2L)
	expect_equal(which(nzchar(b$fail_reason)), drop)
	expect_true(all(b$fail_reason[drop] == "nonconverged"))
	expect_true(all(is.na(b$coefs[drop, ])))
	expect_true(all(b$boot_converged[b$n_valid > 0 & !nzchar(b$fail_reason)]))
	expect_equal(b$n_nonconverged, 0)
	# the interval is the quantile of the 38 valid replicates
	ok = !nzchar(b$fail_reason)
	expect_equal(unname(b$ci_lo), unname(apply(b$coefs[ok, , drop = FALSE], 2, quantile, 0.025)))
	expect_equal(b$max_fail_share, 0.05)
})

test_that("more than 5 percent of dropped replicates is refused, at the replicate that crosses the cap", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	calls = 0L
	real = lame:::.ae_boot_refit
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
			calls <<- calls + 1L
			if (calls %in% c(3L, 5L, 9L)) {
				return(list(
					fit = NULL,
					attempts = 3L,
					iterations = 4L * max_iter,
					start = NA_character_,
					error = NULL,
					last_fit = object
				))
			}
			real(object, Y_b, X_b, Z_b, dyn, max_iter, tol)
		},
		.package = "lame"
	)
	# R = 40 allows 2: the third failure, at replicate 9, stops the run
	err = tryCatch(ame_als_bootstrap(f, R = 40, seed = 3, verbose = FALSE), error = identity)
	expect_s3_class(err, "lame_als_bootstrap_convergence")
	expect_equal(calls, 9L)
	expect_equal(err$replicate, 3L)
	expect_equal(err$n_failed, 3L)
	expect_equal(err$max_failed, 2L)
	expect_equal(err$n_total, 40L)
	expect_match(conditionMessage(err), "5%", fixed = TRUE)
	expect_match(conditionMessage(err), "max_iter", fixed = TRUE)
	# the one-call route
	expect_error(
		lame_als(
			boot_fixture(),
			R = 0,
			family = "normal",
			bootstrap = 40L,
			bootstrap_seed = 3,
			verbose = FALSE
		),
		NA
	)
})

test_that("a replicate that errors or returns non-finite coefficients counts as dropped", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	k = 0L
	real = lame:::.ae_boot_refit
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
			k <<- k + 1L
			if (k == 2L) stop("refit failed")
			out = real(object, Y_b, X_b, Z_b, dyn, max_iter, tol)
			if (k == 6L) out$fit$coefficients[1] = NA_real_
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 60, seed = 3, verbose = FALSE)
	expect_equal(b$n_failed, 2L)
	expect_equal(b$fail_reason[c(2L, 6L)], c("error", "nonfinite"))
})

# no replicate lands at another local optimum

test_that("a fit with a multiplicative term is also refit from the cold start, and the lowest objective wins", {
	y = boot_fixture(n = 16, Tt = 2, R = 1)
	f = lame_als(y, R = 1, family = "normal", bias_correct = "none", verbose = FALSE)
	starts = character()
	real = lame::ame_als_refit
	# the cold start (no factors in the seed object) reaches a lower objective
	local_mocked_bindings(
		ame_als_refit = function(object, ...) {
			out = real(object, ...)
			cold = is.null(object$U)
			starts <<- c(starts, if (cold) "cold" else "warm")
			if (cold) {
				out$deviance = out$deviance - 5
				out$objective = out$objective - 5
			}
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 12, seed = 2, verbose = FALSE)
	expect_equal(sum(starts == "cold"), 12L)
	expect_equal(sum(starts == "warm"), 12L)
	expect_true(all(b$boot_refit_start == "cold"))
	expect_equal(b$n_other_start, 12L)
	expect_match(b$optimum_check, "cold start", fixed = TRUE)
})

test_that("the warm start is kept when the cold start reaches the same optimum", {
	y = boot_fixture(n = 16, Tt = 2, R = 1)
	f = lame_als(y, R = 1, family = "normal", bias_correct = "none", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 12, seed = 2, verbose = FALSE)
	expect_true(all(b$boot_refit_start == "warm"))
	expect_equal(b$n_other_start, 0L)
	# the cold start was run: two attempts' worth of refits per replicate at least
	expect_true(all(b$boot_refit_attempts >= 2L))
})

test_that("a fit without a multiplicative term has a unique optimum and is refit once", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", bias_correct = "none", verbose = FALSE)
	calls = 0L
	real = lame::ame_als_refit
	local_mocked_bindings(
		ame_als_refit = function(object, ...) {
			calls <<- calls + 1L
			real(object, ...)
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 12, seed = 2, verbose = FALSE)
	expect_equal(calls, 12L)
	expect_match(b$optimum_check, "unique", fixed = TRUE)
})

test_that("a multistart fit is checked against the random starts of its multistart", {
	y = boot_fixture(n = 16, Tt = 2, R = 1)
	f = lame_als(y, R = 1, family = "normal", bias_correct = "none", multistart = "cheap", verbose = FALSE)
	calls = 0L
	real = lame::ame_als_refit
	local_mocked_bindings(
		ame_als_refit = function(object, ...) {
			calls <<- calls + 1L
			real(object, ...)
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 10, seed = 2, verbose = FALSE)
	# warm + cold + three random starts
	expect_equal(calls, 10L * 5L)
	expect_match(b$optimum_check, "3 random starts", fixed = TRUE)
})

# the point estimate lies inside its own replicates

test_that("a point estimate outside its own bootstrap interval is an error that carries the bootstrap", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	real = lame:::.ae_boot_refit
	# every replicate is another optimum: the coefficients sit 50 units away
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
			out = real(object, Y_b, X_b, Z_b, dyn, max_iter, tol)
			out$fit$coefficients = out$fit$coefficients + 50
			out
		},
		.package = "lame"
	)
	err = tryCatch(ame_als_bootstrap(f, R = 20, seed = 1, verbose = FALSE), error = identity)
	expect_s3_class(err, "lame_als_bootstrap_optimum")
	expect_match(conditionMessage(err), "outside the 95% interval", fixed = TRUE)
	expect_match(conditionMessage(err), "intercept", fixed = TRUE)
	expect_s3_class(err$boot, "boot_ame")
	expect_equal(err$boot$n_valid, 20L)
	# through the one-call route as well
	expect_error(
		lame_als(
			boot_fixture(),
			R = 0,
			family = "normal",
			bootstrap = 30L,
			verbose = FALSE
		),
		class = "lame_als_bootstrap_optimum"
	)
})

test_that("a variance component outside its interval gets no interval, and the rest of the bootstrap stands", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", bias_correct = "none", verbose = FALSE)
	real = lame:::.ae_boot_refit
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
			out = real(object, Y_b, X_b, Z_b, dyn, max_iter, tol)
			out$fit$VC[["va"]] = out$fit$VC[["va"]] + 5
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 20, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
	expect_identical(b$vc_uncentred$component, "va")
	expect_equal(b$vc_uncentred$estimate, unname(f$VC[["va"]]))
	expect_gt(b$vc_uncentred$lo, b$vc_uncentred$estimate)
	expect_true(is.na(b$vc_se[["va"]]) && is.na(b$vc_ci_lo[["va"]]) && is.na(b$vc_ci_hi[["va"]]))
	# the other components and the coefficients keep their intervals
	ok = setdiff(names(b$vc_se), c("va", "cab", "rho"))
	expect_true(all(is.finite(b$vc_se[ok])))
	expect_true(all(is.finite(b$ci_lo)))
	ci = confint(b, which = "vc")
	expect_true(all(is.na(ci["va", ])))
	expect_true(all(is.finite(ci["vb", ])))
	txt = paste(cli::ansi_strip(capture.output(print(summary(b)))), collapse = " ")
	expect_match(txt, "No interval for va", fixed = TRUE)
})

test_that("the lognormal variance of a count network without extra-Poisson noise gets no interval, and its coefficients keep theirs", {
	# plain Poisson counts: the variance of the lognormal layer is zero, its estimate is
	# above the boundary, and the replicates drawn at that estimate sit higher still
	set.seed(96)
	n = 22; Tt = 5
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Yl = lapply(seq_len(Tt), function(t) {
		Xd = matrix(rnorm(n * n), n, n)
		m = matrix(rpois(n * n, exp(0.3 + 0.4 * Xd + outer(a, b, "+"))), n, n)
		diag(m) = NA
		list(Y = m, X = Xd)
	})
	fit = lame_als(
		lapply(Yl, `[[`, "Y"),
		Xdyad = lapply(Yl, `[[`, "X"),
		R = 0,
		family = "poisson",
		verbose = FALSE
	)
	expect_identical(fit$bias_correct$method, "eb")
	b = ame_als_bootstrap(fit, R = 30, seed = 1, verbose = FALSE)
	expect_identical(b$vc_uncentred$component, "ve")
	expect_lt(b$vc_uncentred$estimate, b$vc_uncentred$lo)
	expect_true(is.na(b$vc_ci_lo[["ve"]]))
	expect_true(all(is.finite(b$se)))
	expect_true(all(b$point_est >= b$ci_lo & b$point_est <= b$ci_hi))
})

test_that("an ordinary bootstrap keeps every point estimate inside its interval", {
	f = lame_als(boot_fixture(n = 16, Tt = 4), R = 0, family = "normal", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 40, seed = 5, verbose = FALSE)
	expect_true(all(b$point_est >= b$ci_lo & b$point_est <= b$ci_hi))
	ok = is.finite(b$vc_point_est) & is.finite(b$vc_ci_lo)
	expect_true(all(b$vc_point_est[ok] >= b$vc_ci_lo[ok] - 1e-8 & b$vc_point_est[ok] <= b$vc_ci_hi[ok] + 1e-8))
})

# near-equal singular values

# a rank-two term whose two singular values are set directly: orthogonal
# unit-norm factor columns times the singular values
rank2_panel = function(s, n = 30, Tt = 3, seed = 21, noise = 0.5) {
	set.seed(seed)
	qu = qr.Q(qr(matrix(rnorm(n * 2), n, 2)))
	qv = qr.Q(qr(matrix(rnorm(n * 2), n, 2)))
	O = qu %*% diag(s) %*% t(qv)
	lapply(seq_len(Tt), function(t) {
		m = O + matrix(rnorm(n * n, 0, noise), n)
		diag(m) = NA_real_
		m
	})
}

test_that("columns with near-equal singular values get NA column-wise standard errors and a subspace summary", {
	y = rank2_panel(c(14, 13.8))
	f = lame_als(y, R = 2, family = "normal", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 60, seed = 4, verbose = FALSE)
	sp = b$sv_separation
	expect_equal(sp$threshold, lame:::.ae_boot_sv_z)
	expect_length(sp$z, 1L)
	expect_lt(sp$z, sp$threshold)
	expect_identical(sp$resolved, c(FALSE, FALSE))
	# no column-wise standard errors for the columns, raw or aligned
	expect_true(all(is.na(b$se_U)))
	expect_true(all(is.na(b$se_V)))
	expect_true(all(is.na(b$se_U_raw)))
	# the subspace is described instead
	expect_equal(dim(b$subspace$U), c(60L, 2L))
	expect_true(all(is.finite(b$subspace$U)))
	expect_true(all(b$subspace$U >= 0 & b$subspace$U <= 90))
	expect_equal(rownames(b$subspace$largest), c("U", "V"))
	expect_true(all(is.finite(b$subspace$largest)))
	# and the intervals of the unresolved columns are NA
	ci = confint(b, which = "U")
	expect_true(all(is.na(ci)))
	expect_false(anyNA(confint(b, which = "beta")))
	# the summary says so
	txt = paste(cli::ansi_strip(capture.output(print(summary(b)))), collapse = " ")
	expect_match(txt, "not separated", fixed = TRUE)
	expect_match(txt, "principal angle", fixed = TRUE)
})

test_that("separated singular values keep their column-wise standard errors", {
	y = rank2_panel(c(14, 6))
	f = lame_als(y, R = 2, family = "normal", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 60, seed = 4, verbose = FALSE)
	expect_gte(b$sv_separation$z, b$sv_separation$threshold)
	expect_identical(b$sv_separation$resolved, c(TRUE, TRUE))
	expect_true(all(is.finite(b$se_U)))
	expect_true(all(is.finite(b$se_V)))
	expect_false(anyNA(confint(b, which = "U")))
	# the subspace summary is there either way
	expect_true(all(is.finite(b$subspace$largest)))
})

test_that("a rank-one term has no pair to separate and keeps its standard errors", {
	# one period: the product of the fixture is redrawn in each period, which a static fit sees averaged (and, below two
	# noise variances of signal per loading, withholds)
	y = boot_fixture(n = 16, Tt = 1, R = 1)
	f = lame_als(y, R = 1, family = "normal", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 30, seed = 4, verbose = FALSE)
	expect_identical(b$sv_separation$resolved, TRUE)
	expect_true(all(is.finite(b$se_U)))
})

test_that("the principal angles are zero for the same subspace and 90 degrees for an orthogonal one", {
	ang = lame:::.ae_boot_angles
	set.seed(1)
	A = matrix(rnorm(20 * 2), 20, 2)
	expect_equal(ang(A, A), c(0, 0), tolerance = 1e-6)
	rot = matrix(c(cos(0.7), sin(0.7), -sin(0.7), cos(0.7)), 2)
	expect_equal(ang(A %*% rot, A), c(0, 0), tolerance = 1e-6)
	B = qr.Q(qr(cbind(A, rnorm(20), rnorm(20))))[, 3:4]
	expect_equal(ang(B, A), c(90, 90), tolerance = 1e-6)
})
