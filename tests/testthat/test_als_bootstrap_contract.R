library(lame)
skip_on_cran()

# the ALS bootstrap contract: a minimum number of replicates, no unconverged replicate in an interval, a capped
# share of dropped replicates, one refit per replicate, a point estimate inside its own replicates,
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
		.ae_boot_refit = function(object, Y_b, X_b, max_iter, tol) {
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
			real(object, Y_b, X_b, max_iter, tol)
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
		.ae_boot_refit = function(object, Y_b, X_b, max_iter, tol) {
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
			real(object, Y_b, X_b, max_iter, tol)
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
		.ae_boot_refit = function(object, Y_b, X_b, max_iter, tol) {
			k <<- k + 1L
			if (k == 2L) stop("refit failed")
			out = real(object, Y_b, X_b, max_iter, tol)
			if (k == 6L) out$fit$coefficients[1] = NA_real_
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 60, seed = 3, verbose = FALSE)
	expect_equal(b$n_failed, 2L)
	expect_equal(b$fail_reason[c(2L, 6L)], c("error", "nonfinite"))
})

# one refit per replicate

test_that("a replicate is refit once, from the variational state of the fit", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
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
	expect_match(b$optimum_check, "variational state of the fit", fixed = TRUE)
})

# the point estimate lies inside its own replicates

test_that("a coefficient outside its own bootstrap interval takes the analytic error of the fit", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	real = lame:::.ae_boot_refit
	# every replicate is another optimum: the coefficients sit 50 units away
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, max_iter, tol) {
			out = real(object, Y_b, X_b, max_iter, tol)
			out$fit$coefficients = out$fit$coefficients + 50
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 20, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
	expect_equal(b$n_valid, 20L)
	expect_setequal(b$coef_uncentred$parameter, names(coef(f)))
	expect_equal(unname(b$se), unname(f$eb$se[names(coef(f))]), tolerance = 1e-12)
	expect_true(all(names(coef(f)) %in% b$se_source$parameter))
	expect_match(b$se_source$reason[1], "outside the interval of its replicates", fixed = TRUE)
	expect_true(all(b$point_est > b$ci_lo & b$point_est < b$ci_hi))
	# through the one-call route the fit carries the bootstrap and its analytic errors
	g = lame_als(boot_fixture(), R = 0, family = "normal", bootstrap = 30L, verbose = FALSE)
	expect_equal(sqrt(diag(vcov(g))), sqrt(diag(vcov(f))), tolerance = 1e-12)
	expect_true(all(is.finite(confint(g, parm = names(coef(g))))))
})

test_that("a variance component outside its interval gets the analytic one, and the rest of the bootstrap stands", {
	f = lame_als(boot_fixture(), R = 0, family = "normal", verbose = FALSE)
	real = lame:::.ae_boot_refit
	local_mocked_bindings(
		.ae_boot_refit = function(object, Y_b, X_b, max_iter, tol) {
			out = real(object, Y_b, X_b, max_iter, tol)
			out$fit$VC[["va"]] = out$fit$VC[["va"]] + 5
			out
		},
		.package = "lame"
	)
	b = ame_als_bootstrap(f, R = 20, seed = 1, verbose = FALSE)
	expect_s3_class(b, "boot_ame")
	# va is pushed off its replicates (the empirical-Bayes vb may be too, by its own bias)
	un = b$vc_uncentred[b$vc_uncentred$component == "va", ]
	expect_equal(nrow(un), 1L)
	expect_equal(un$estimate, unname(f$VC[["va"]]))
	# outside its interval (for an empirical-Bayes fit the bias-corrected basic interval, which sits below)
	expect_true(un$estimate < un$lo || un$estimate > un$hi)
	expect_equal(b$vc_se[["va"]], f$VC_se[["va"]], tolerance = 1e-12)
	expect_equal(unname(c(b$vc_ci_lo[["va"]], b$vc_ci_hi[["va"]])), unname(f$VC_ci["va", ]), tolerance = 1e-12)
	expect_true("va" %in% b$se_source$parameter)
	# the other components and the coefficients keep their intervals
	ok = setdiff(names(b$vc_se), c("va", "cab", "rho"))
	expect_true(all(is.finite(b$vc_se[ok])))
	expect_true(all(is.finite(b$ci_lo)))
	ci = confint(b, which = "vc")
	expect_equal(unname(ci["va", ]), unname(f$VC_ci["va", ]), tolerance = 1e-12)
	expect_true(all(is.finite(ci["vb", ])))
	txt = paste(cli::ansi_strip(capture.output(print(summary(b)))), collapse = " ")
	expect_match(txt, "Analytic standard error and interval for va", fixed = TRUE)
})

test_that("the lognormal variance of a count network without extra-Poisson noise gets the analytic interval, and its coefficients keep theirs", {
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
	expect_null(fit$bias_correct)
	b = ame_als_bootstrap(fit, R = 30, seed = 1, verbose = FALSE)
	# ve has a standard error and an interval either way: the analytic ones where its replicates
	# sit off the estimate, the replicates' otherwise
	expect_true(is.finite(b$vc_se[["ve"]]) && b$vc_se[["ve"]] > 0)
	expect_true(is.finite(b$vc_ci_lo[["ve"]]) && is.finite(b$vc_ci_hi[["ve"]]))
	if ("ve" %in% b$vc_uncentred$component) {
		expect_lt(b$vc_uncentred$estimate[b$vc_uncentred$component == "ve"], b$vc_uncentred$lo[b$vc_uncentred$component == "ve"])
		expect_equal(b$vc_se[["ve"]], fit$VC_se[["ve"]], tolerance = 1e-12)
		expect_equal(b$vc_ci_lo[["ve"]], fit$VC_ci["ve", 1L], tolerance = 1e-12)
	}
	expect_true(all(is.finite(b$se)))
	expect_true(all(b$point_est >= b$ci_lo & b$point_est <= b$ci_hi))
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

test_that("columns with near-equal singular values take the analytic column-wise standard errors and get a subspace summary", {
	y = rank2_panel(c(14, 13.8))
	f = lame_als(y, R = 2, family = "normal", verbose = FALSE)
	b = ame_als_bootstrap(f, R = 60, seed = 4, verbose = FALSE)
	sp = b$sv_separation
	expect_equal(sp$threshold, lame:::.ae_boot_sv_z)
	expect_length(sp$z, 1L)
	expect_lt(sp$z, sp$threshold)
	expect_identical(sp$resolved, c(FALSE, FALSE))
	# the replicates give the columns no column-wise standard errors, so they keep the analytic ones
	expect_equal(b$se_U, as.vector(f$U_sd), tolerance = 1e-12)
	expect_equal(b$se_V, as.vector(f$V_sd), tolerance = 1e-12)
	expect_equal(b$se_U_raw, as.vector(f$U_sd), tolerance = 1e-12)
	expect_true(all(c("U", "V") %in% b$se_source$parameter))
	# the subspace is described instead
	expect_equal(dim(b$subspace$U), c(60L, 2L))
	expect_true(all(is.finite(b$subspace$U)))
	expect_true(all(b$subspace$U >= 0 & b$subspace$U <= 90))
	expect_equal(rownames(b$subspace$largest), c("U", "V"))
	expect_true(all(is.finite(b$subspace$largest)))
	# and the intervals of the unresolved columns are the analytic Wald ones
	ci = confint(b, which = "U")
	expect_true(all(is.finite(ci)))
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
