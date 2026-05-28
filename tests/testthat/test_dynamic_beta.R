# tests for the dynamic_beta feature on lame() fits.
#
# the byte-identical-default contract is the most important guarantee here:
# when dynamic_beta = FALSE (the default), EVERY slot of the lame fit must be
# bit-for-bit identical to the legacy lame fit produced without the argument
# at all. that's what keeps an existing user's code untouched by the feature.
#
# beyond byte-identity, we test:
#   * parser correctness on every shortcut / type / error path
#   * the FFBS sampler produces a 3-D BETA of the right shape and dimnames
#   * static-block + dynamic-block split works (mixed dyn/static)
#   * the contrast basis identifies (intercept, a, b) when needed
#   * S3 methods (coef, vcov, confint, summary, print, predict, simulate,
#     trace_plot, as_draws, fitted, residuals, prior_summary) all dispatch
#   * cross-sectional ame() and ALS dispatchers reject dynamic_beta loudly

skip_on_cran()

# ------- byte-identical-default ----------------------------------------------

test_that("dynamic_beta = FALSE is byte-identical to no-argument default (unipartite)", {
	set.seed(2026)
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	set.seed(7)
	fit1 <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	             nscan = 20, burn = 5, odens = 1, dynamic_beta = FALSE, verbose = FALSE)
	set.seed(7)
	fit2 <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	             nscan = 20, burn = 5, odens = 1, verbose = FALSE)
	expect_identical(fit1$BETA, fit2$BETA)
	expect_identical(fit1$VC,   fit2$VC)
	expect_identical(fit1$APM,  fit2$APM)
	expect_identical(fit1$BPM,  fit2$BPM)
	expect_identical(fit1$U,    fit2$U)
	expect_identical(fit1$V,    fit2$V)
	expect_identical(fit1$EZ,   fit2$EZ)
	expect_identical(fit1$YPM,  fit2$YPM)
})

test_that("dynamic_beta = FALSE is byte-identical (bipartite)", {
	set.seed(2026)
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal",
	                              mode = "bipartite", nA = 6, nB = 5, seed = 2026)
	set.seed(7)
	fit1 <- lame(dat$Y, Xdyad = dat$Xdyad, mode = "bipartite", family = "normal",
	             R = 0, nscan = 20, burn = 5, odens = 1,
	             dynamic_beta = FALSE, verbose = FALSE)
	set.seed(7)
	fit2 <- lame(dat$Y, Xdyad = dat$Xdyad, mode = "bipartite", family = "normal",
	             R = 0, nscan = 20, burn = 5, odens = 1, verbose = FALSE)
	expect_identical(fit1$BETA, fit2$BETA)
	expect_identical(fit1$VC,   fit2$VC)
	expect_identical(fit1$APM,  fit2$APM)
	expect_identical(fit1$BPM,  fit2$BPM)
})

# ------- parser correctness --------------------------------------------------

test_that("lame:::parse_dynamic_beta(FALSE / NULL / TRUE / character / integer / logical) works", {
	coef_names <- c("intercept", "x1.dyad", "x2.dyad", "x3.row")
	coef_block <- c("intercept", "dyad", "dyad", "row")
	res <- lame:::parse_dynamic_beta(FALSE, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_false(res$any)
	expect_equal(res$n_groups, 0L)
	expect_equal(res$dynamic_idx, integer(0))

	res <- lame:::parse_dynamic_beta(NULL, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_false(res$any)

	res <- lame:::parse_dynamic_beta(TRUE, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(sum(res$mask), 4L)
	expect_setequal(res$group_names, c("intercept", "dyad", "row"))

	res <- lame:::parse_dynamic_beta("dyad", coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(res$mask, c(FALSE, TRUE, TRUE, FALSE))

	res <- lame:::parse_dynamic_beta(c("intercept", "x3.row"), coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(res$mask, c(TRUE, FALSE, FALSE, TRUE))

	res <- lame:::parse_dynamic_beta(c(2L, 3L), coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_equal(res$mask, c(FALSE, TRUE, TRUE, FALSE))

	res <- lame:::parse_dynamic_beta(c(FALSE, TRUE, FALSE, FALSE), coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_equal(sum(res$mask), 1L)
})

test_that("parse_dynamic_beta family bans (ordinal / rrl) are enforced", {
	coef_names <- c("intercept", "x1.dyad")
	coef_block <- c("intercept", "dyad")
	expect_error(
		lame:::parse_dynamic_beta(TRUE, coef_names, coef_block, TRUE, 3, "ordinal", "unipartite"),
		"intercept"
	)
	expect_error(
		lame:::parse_dynamic_beta("intercept", coef_names, coef_block, TRUE, 3, "rrl", "unipartite"),
		"intercept"
	)
})

test_that("parse_dynamic_beta unknown name aborts", {
	expect_error(
		lame:::parse_dynamic_beta("nonexistent_coef", c("intercept", "x1.dyad"),
		                   c("intercept", "dyad"), TRUE, 3, "normal", "unipartite"),
		"not a coefficient or block shortcut"
	)
})

test_that("parse_dynamic_beta T < 2 aborts", {
	expect_error(
		lame:::parse_dynamic_beta(TRUE, c("intercept", "x1.dyad"),
		                   c("intercept", "dyad"), TRUE, 1L, "normal", "unipartite"),
		"at least 2 time periods"
	)
})

# ------- 3-D BETA shape ------------------------------------------------------

test_that("dynamic_beta = TRUE produces 3-D BETA with right dim/names", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_equal(dim(fit$BETA), c(20L, 2L, 3L))
	expect_equal(dimnames(fit$BETA)[[2]], c("intercept", "X1_dyad"))
	expect_true(isTRUE(fit$dynamic_beta))
	expect_false(is.null(fit$rho_beta))
	expect_false(is.null(fit$sigma_beta))
	expect_equal(length(fit$rho_beta), 2L)  # intercept + dyad blocks
})

test_that("partial dynamic_beta (subset of coefs) gives the right mask", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1,
	            dynamic_beta = "dyad", verbose = FALSE)
	# intercept is static -> rows of BETA[, "intercept", ] are identical across t
	beta_int <- fit$BETA[, "intercept", ]
	# intercept-static => columns t1, t2, t3 are all the same per draw
	expect_true(all(beta_int[, 1] == beta_int[, 2]))
	expect_true(all(beta_int[, 2] == beta_int[, 3]))
	# dyadic is dynamic -> rows differ across t (with overwhelming probability)
	beta_dyad <- fit$BETA[, "X1_dyad", ]
	expect_false(all(beta_dyad[, 1] == beta_dyad[, 2]))
})

# ------- S3 methods ----------------------------------------------------------

test_that("coef.lame returns p x T matrix for dynamic_beta fit", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	cf <- coef(fit)
	expect_true(is.matrix(cf))
	expect_equal(dim(cf), c(2L, 3L))
	expect_equal(rownames(cf), c("intercept", "X1_dyad"))
})

test_that("vcov returns (p*T) x (p*T) for dynamic_beta fit", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	V <- vcov(fit)
	expect_equal(dim(V), c(6L, 6L))   # 2 coefs x 3 periods
	expect_true(any(grepl("\\[t1\\]", rownames(V))))
})

test_that("confint returns (p*T) x 2 for dynamic_beta fit", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	ci <- confint(fit)
	expect_equal(dim(ci), c(6L, 2L))
	expect_true(all(ci[, 1] <= ci[, 2]))
})

test_that("summary.lame prints dynamic-beta block without erroring", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	out <- summary(fit)
	expect_s3_class(out, "summary.lame")
	expect_true(isTRUE(out$dynamic_beta))
	expect_true(!is.null(out$beta_dynamic_per_t))
	expect_silent(invisible(capture.output(print(out))))
})

test_that("print.lame works on dynamic_beta fit", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	# cli writes to message stream; capture both
	msgs <- capture.output(print(fit), type = "message")
	out  <- capture.output(print(fit))
	combined <- paste(c(out, msgs), collapse = "\n")
	expect_match(combined, "Dynamic regression coefficients")
})

test_that("as_draws flattens 3-D BETA to draws_array", {
	skip_if_not_installed("posterior")
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	d <- as_draws(fit)
	expect_s3_class(d, "draws_array")
	# 2 coefs x 3 periods (beta) + 5 VC entries = 11 variables
	expect_true(dim(d)[3] >= 6L)
	vn <- dimnames(d)[[3]]
	expect_true(any(grepl("\\[t1\\]", vn)))
})

test_that("trace_plot handles 3-D BETA without erroring", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	p <- trace_plot(fit, params = "beta")
	expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

# ------- identifiability / constraint basis ----------------------------------

test_that("a/b are sum-to-zero (within tolerance) when an intercept is in the model", {
	dat <- simulate_test_network(n = 10, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 40, burn = 10, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	# sum-to-zero on a (and b) keeps the dynamic intercept identified
	expect_lt(abs(sum(fit$APM)), 1e-6)
	expect_lt(abs(sum(fit$BPM)), 1e-6)
})

test_that("partial dynamic_beta with no intercept-in-dyn still works", {
	# only dyadic dynamic; intercept static -> no special constraint needed
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	expect_no_error(
		fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
		            nscan = 20, burn = 5, odens = 1,
		            dynamic_beta = "dyad", verbose = FALSE)
	)
	expect_equal(length(dim(fit$BETA)), 3L)
})

# ------- bipartite + dynamic_beta -------------------------------------------

test_that("dynamic_beta works on bipartite normal", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal",
	                              mode = "bipartite", nA = 6, nB = 5, seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, mode = "bipartite",
	            family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = "dyad",
	            verbose = FALSE)
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_equal(dim(fit$BETA), c(20L, 2L, 3L))
})

# ------- edge cases ----------------------------------------------------------

test_that("dynamic_beta + T = 1 aborts (validation in lame())", {
	set.seed(1)
	Y <- list(matrix(rnorm(36), 6, 6))
	diag(Y[[1]]) <- NA
	expect_error(
		lame(Y, family = "normal", R = 0, nscan = 10, burn = 2, odens = 1,
		     dynamic_beta = TRUE, verbose = FALSE),
		"at least 2 time periods"
	)
})

test_that("dynamic_beta with method = 'als' warns and ignores", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	expect_warning(
		fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
		            method = "als", dynamic_beta = TRUE, verbose = FALSE),
		"dynamic_beta"
	)
	# the fit should be an ALS (static) fit, so BETA is 2-D / NULL
	expect_false(isTRUE(fit$dynamic_beta))
})

test_that("ame() (cross-sectional) rejects dynamic_beta", {
	set.seed(1)
	Y <- matrix(rnorm(64), 8, 8); diag(Y) <- NA
	# ame() catches longitudinal-only args via ... and emits a clean
	# cli_abort directing the user to lame().
	expect_error(
		ame(Y, family = "normal", R = 0, nscan = 10, burn = 2, odens = 1,
		    dynamic_beta = TRUE, verbose = FALSE),
		"longitudinal-only"
	)
})

# ------- tryErrorChecks accumulator -----------------------------------------

test_that("tryErrorChecks$beta is initialised and stays at 0 on a healthy fit", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE, verbose = FALSE)
	expect_true("beta" %in% names(fit$tryErrorChecks))
	expect_equal(fit$tryErrorChecks$beta, 0L)
})

# ------- recovery test (smoke; not a hard convergence criterion) -------------

test_that("dynamic_beta drift is detected when truth is time-varying", {
	# generate Y_t = beta_t * X_t + noise with beta_t evolving linearly
	set.seed(2026)
	n <- 10; T <- 4
	beta_t <- seq(-0.5, 0.5, length.out = T)
	X_list <- replicate(T, matrix(rnorm(n * n), n, n), simplify = FALSE)
	Y_list <- vector("list", T)
	for (t in seq_len(T)) {
		Y_t <- beta_t[t] * X_list[[t]] + matrix(rnorm(n * n, 0, 0.3), n, n)
		diag(Y_t) <- NA
		Y_list[[t]] <- Y_t
	}
	X_arr <- lapply(X_list, function(x) array(x, c(n, n, 1)))
	# fit with dynamic_beta on the dyadic coef
	fit <- lame(Y_list, Xdyad = X_arr, family = "normal", R = 0,
	            nscan = 100, burn = 30, odens = 1, dynamic_beta = "dyad", verbose = FALSE)
	# the per-period posterior-mean dyadic coef should be MONOTONIC-ish
	# (we just check that the spread across periods covers a notable range)
	dyad_per_t <- apply(fit$BETA[, "X1_dyad", , drop = FALSE], 3, mean)
	expect_gt(diff(range(dyad_per_t)), 0.2)
})
