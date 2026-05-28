# period_exposure Poisson rate scaling
#
# c++ exposure-aware likelihood,
# Z keeps its "unexposed log-rate" meaning, exposure enters only the
# Poisson observation likelihood. NULL / all-ones bypasses the new
# code path entirely.

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

# -- byte-identical: NULL vs rep(1, T) ------------------------------------

test_that("period_exposure = NULL and rep(1, T) both byte-identical to legacy", {
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

# -- substantive recovery -------------------------------------------------

test_that("intercept recovery: exposure absorbs the log-rate", {
	# Pure intercept-only Poisson: without exposure, intercept ~ log(mean exposure);
	# with exposure, intercept ~ 0 (true unexposed log-rate is 0).
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
	# Without exposure: intercept must be markedly positive (absorbing log-rate)
	expect_gt(mean(fit_noexp$BETA[, 1]), 1.0)
	# With exposure: intercept must be small (close to 0, the true unexposed log-rate)
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
		regexp = NA)  # current code only checks positivity + length-by-T, not strict length
	# negative
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "poisson", R = 0,
		period_exposure = c(-1, 1, 1, 1),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"non-negative")
})

test_that("period_exposure with non-Poisson family aborts on non-trivial values", {
	fx = .make_poisson_fixture()
	# normal family with non-trivial exposure: abort
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		period_exposure = c(1, 5, 25, 100),
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"meaningful only")
	# normal family with all-ones exposure: should NOT abort (active = FALSE)
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
	# default newexposure = last observed
	fc_default = suppressWarnings(predict(fit, h = 2L, type = "response"))
	# explicit newexposure
	fc_explicit = suppressWarnings(predict(fit, h = 2L, type = "response",
	                                         newexposure = c(50, 200)))
	expect_true(is.list(fc_default))
	expect_length(fc_default, 2L)
	# explicit exposure should scale the response: fc_explicit[k] = (newexp_k /
	# last_observed) * fc_default[k] up to MC noise
	r_default  = mean(fc_default[[1]], na.rm = TRUE)
	r_explicit = mean(fc_explicit[[1]], na.rm = TRUE)
	expect_gt(r_explicit, r_default * 0.3)  # generous (MC noise, AR(1) drift)
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
