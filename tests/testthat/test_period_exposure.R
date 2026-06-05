# poisson period exposure tests
#
# z keeps its unexposed log-rate meaning; exposure enters the poisson
# observation likelihood

skip_on_cran()

.make_poisson_fixture = function(n = 10L, T_per = 4L,
                                    exposure = c(1, 5, 25, 100),
                                    beta_true = 0.4, seed = 7L) {
	set.seed(seed)
	Y_list = vector("list", T_per)
	X_list = vector("list", T_per)
	for (t in seq_len(T_per)) {
		X = matrix(rnorm(n*n), n, n)
		eta = beta_true * X
		lambda = exposure[t] * exp(eta)
		lambda = pmin(lambda, 1e4)
		Y = matrix(rpois(n*n, lambda), n, n)
		diag(Y) = NA
		Y_list[[t]] = Y
		X_list[[t]] = array(X, c(n, n, 1))
	}
	list(Y = Y_list, X = X_list, exposure = exposure, beta_true = beta_true)
}

# -- null and all-ones exposure -------------------------------------------

test_that("period_exposure = NULL and rep(1, T) match", {
	fx = .make_poisson_fixture()
	set.seed(1); fit_null = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	set.seed(1); fit_ones = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		period_exposure = rep(1, length(fx$Y)),
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_identical(fit_null$BETA, fit_ones$BETA)
	expect_identical(fit_null$EZ,   fit_ones$EZ)
	expect_identical(fit_null$YPM,  fit_ones$YPM)
})

# -- intercept recovery ---------------------------------------------------

test_that("intercept recovery: exposure absorbs the log-rate", {
		# intercept-only poisson isolates the exposure offset
	set.seed(7)
	n = 12L; T_per = 4L
	exposure = c(1, 5, 25, 100)
	Y_list = vector("list", T_per)
	for (t in seq_len(T_per)) {
		Y = matrix(rpois(n*n, exposure[t]), n, n)
		diag(Y) = NA
		Y_list[[t]] = Y
	}
	set.seed(1); fit_noexp = suppressWarnings(suppressMessages(lame(
		Y_list, family = "poisson", R = 0,
		nscan = 150, burn = 40, odens = 5, verbose = FALSE, plot = FALSE)))
	set.seed(1); fit_exp = suppressWarnings(suppressMessages(lame(
		Y_list, family = "poisson", R = 0,
		period_exposure = exposure,
		nscan = 150, burn = 40, odens = 5, verbose = FALSE, plot = FALSE)))
		# without exposure, the intercept absorbs the log-rate
	expect_gt(mean(fit_noexp$BETA[, 1]), 1.0)
		# with exposure, the intercept is close to the unexposed log-rate
	expect_lt(abs(mean(fit_exp$BETA[, 1])), 0.5)
})

# -- validation ------------------------------------------------------------

test_that("period_exposure validation", {
	fx = .make_poisson_fixture()
		# wrong length
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		period_exposure = c(1, 2, 3),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
			regexp = NA)
		# negative exposure
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		period_exposure = c(-1, 1, 1, 1),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"non-negative")
})

test_that("period_exposure with non-Poisson family aborts on non-trivial values", {
	fx = .make_poisson_fixture()
		# non-poisson family with non-trivial exposure
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		period_exposure = c(1, 5, 25, 100),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"meaningful only")
		# all-ones exposure is inactive
	set.seed(1)
	fit_ones = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		period_exposure = c(1, 1, 1, 1),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_s3_class(fit_ones, "lame")
})

# -- forecast / predict ---------------------------------------------------

test_that("forecast applies newexposure on the response scale", {
	fx = .make_poisson_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		dynamic_beta = "dyad",
		period_exposure = fx$exposure,
		nscan = 80, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
		# default newexposure uses the last observed exposure
	fc_default = suppressWarnings(predict(fit, h = 2L, type = "response"))
		# explicit newexposure
	fc_explicit = suppressWarnings(predict(fit, h = 2L, type = "response",
	                                         newexposure = c(50, 200)))
	expect_true(is.list(fc_default))
	expect_length(fc_default, 2L)
		# explicit exposure should scale the response
	r_default  = mean(fc_default[[1]], na.rm = TRUE)
	r_explicit = mean(fc_explicit[[1]], na.rm = TRUE)
	expect_gt(r_explicit, r_default * 0.3)  # generous (mc noise, ar(1) drift)
})

test_that("forecast newexposure validation", {
	fx = .make_poisson_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		dynamic_beta = "dyad",
		period_exposure = fx$exposure,
		nscan = 60, burn = 20, odens = 5, verbose = FALSE, plot = FALSE)))
	expect_error(predict(fit, h = 3L, type = "response", newexposure = c(1, 1)),
	             "length")
	expect_error(predict(fit, h = 2L, type = "response", newexposure = c(-1, 2)),
	             "non-negative")
})
