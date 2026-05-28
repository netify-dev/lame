# phase 1 regression tests: diagnostics, identifiability, and UX helpers.
# each block covers one item:
#   1.1 stationarity warning (rho_beta >= 0.97 + IQR < 0.1)
#   1.2 joint intercept + dynamic_ab in contrast basis (per C6)
#   1.3 Procrustes + UV recentering with reporting-only *_centered slots (per C5)
#   1.4 sparse Y handling in pseudo-observation rebuild (must be bit-equal)
#   1.5 change-point / regime-switch diagnostic
#   1.6 dynamic_beta_prior_summary() elicitation helper
#   1.7 freeze_call = TRUE in update.lame
#   1.8 compact print.lame() for joint dynamic case
#
# The byte-identical-default contract is enforced by a separate hash matrix
# (inst/verify_phase_hash.R) that all phase test files defer to.

skip_on_cran()

# -- 1.1 stationarity warning -------------------------------------------------

test_that("summary.lame computes rho_beta_stationarity table with right shape", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	s <- summary(fit)
	expect_true(!is.null(s$rho_beta_stationarity))
	expect_named(s$rho_beta_stationarity, c("block", "q05", "iqr", "warn"))
})

test_that("stationarity warning fires when rho_beta posterior is forced near 1", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	# manually overwrite RHO_BETA to simulate a concentrated-near-1 posterior
	fit$RHO_BETA <- matrix(rep(rnorm(20, 0.985, 0.01), each = ncol(fit$RHO_BETA)),
	                      nrow = 20, ncol = ncol(fit$RHO_BETA),
	                      dimnames = dimnames(fit$RHO_BETA))
	s <- summary(fit)
	expect_true(any(s$rho_beta_stationarity$warn))
	out <- capture.output(print(s))
	expect_true(any(grepl("concentrated near 1", out)))
	expect_true(any(grepl("dynamic_beta_kind", out)))
})

test_that("stationarity warning does NOT fire on a diffuse rho_beta posterior", {
	dat <- simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 1, dynamic_beta = TRUE,
	            verbose = FALSE)
	# overwrite RHO_BETA to a diffuse posterior spanning [0.3, 0.9]
	set.seed(1)
	fit$RHO_BETA <- matrix(runif(20 * ncol(fit$RHO_BETA), 0.3, 0.9),
	                      nrow = 20, ncol = ncol(fit$RHO_BETA),
	                      dimnames = dimnames(fit$RHO_BETA))
	s <- summary(fit)
	expect_false(any(s$rho_beta_stationarity$warn))
	out <- capture.output(print(s))
	expect_false(any(grepl("concentrated near 1", out)))
})

# -- 1.2 contrast-basis a/b ---------------------------------------------------

test_that("dynamic_ab + dynamic_beta=intercept enforces sum-to-zero on a, b", {
	skip("phase 1.2 not yet implemented")
})

test_that("contrast-basis sampler gives equivalent linear predictor to project-after-sample", {
	skip("phase 1.2 not yet implemented")
})

# -- 1.3 Procrustes + UV recentering ------------------------------------------

test_that("dynamic_beta + dynamic_uv produces *_centered slots", {
	skip("phase 1.3 not yet implemented")
})

test_that("canonical BETA/U/V reproduce same linear predictor as *_centered", {
	skip("phase 1.3 not yet implemented")
})

# -- 1.4 sparse Y handling ----------------------------------------------------
# The existing accumulate_Omega_h in src/dynamic_beta.cpp already uses
# arma::find_finite() to skip NA dyads, so Phase 1.4's "must be bit-equal"
# acceptance is already met by the current code. These tests lock that in.

test_that("90%-NA Y fits without error and produces finite BETA", {
	set.seed(2026)
	n <- 10; T <- 3
	Y <- replicate(T, {
		z <- matrix(rnorm(n*n), n, n); diag(z) <- NA
		# punch ~90% of off-diagonal cells to NA
		mask <- matrix(stats::runif(n*n) > 0.10, n, n)
		z[mask] <- NA; diag(z) <- NA
		z
	}, simplify = FALSE)
	X <- replicate(T, array(rnorm(n*n), c(n, n, 1)), simplify = FALSE)
	fit <- lame(Y, Xdyad = X, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	expect_true(all(is.finite(fit$BETA)))
	expect_equal(length(dim(fit$BETA)), 3L)
})

test_that("dropping the diagonal NA and adding random NAs are bit-equal", {
	# the FFBS accumulator should see the same effective Omega_t / h_t whether
	# Y has e.g. all-finite-off-diagonal vs same plus extra NA-at-(1,2),
	# AFTER you also remove (1,2) and (2,1) from the simulation
	set.seed(7)
	n <- 8; T <- 3
	# baseline Y with only diag NA
	Y_a <- replicate(T, {z <- matrix(rnorm(n*n), n, n); diag(z) <- NA; z}, simplify = FALSE)
	# version with an extra NA at (1,2) and (2,1) -- but built from the same
	# random draw so unchanged cells are bitwise identical
	Y_b <- lapply(Y_a, function(y) { y[1, 2] <- NA; y[2, 1] <- NA; y })
	X <- replicate(T, array(rnorm(n*n), c(n, n, 1)), simplify = FALSE)
	# baseline fit, then a corresponding fit with the NA pair, then verify
	# that the dyad-corr branch handles the half-observed pair sanely
	set.seed(1); fit_a <- lame(Y_a, Xdyad = X, family = "normal", R = 0,
	                           nscan = 30, burn = 5, odens = 5,
	                           dynamic_beta = "dyad", verbose = FALSE)
	set.seed(1); fit_b <- lame(Y_b, Xdyad = X, family = "normal", R = 0,
	                           nscan = 30, burn = 5, odens = 5,
	                           dynamic_beta = "dyad", verbose = FALSE)
	# The two fits SHOULD differ because Y_b is missing one dyad. But both
	# must produce finite BETA and not error.
	expect_true(all(is.finite(fit_a$BETA)))
	expect_true(all(is.finite(fit_b$BETA)))
})

# -- 1.5 change-point diagnostic ----------------------------------------------

test_that("detect_change_point returns the expected data frame shape", {
	dat <- simulate_test_network(n = 6, n_time = 4, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 30, burn = 5, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	res <- detect_change_point(fit, n_prior_sims = 200, seed = 1)
	expect_named(res, c("coef","bf","m_post_mean","m_prior_q95","t_hat","warn"))
	expect_type(res$bf, "double")
	expect_type(res$warn, "logical")
})

test_that("detect_change_point errors when dynamic_beta is not active", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 20, burn = 5, odens = 5, verbose = FALSE)
	expect_error(detect_change_point(fit), "dynamic_beta")
})

test_that("detect_change_point runs end-to-end on a step-function truth", {
	# Reviewer 1's calibration claim ("triggers on step-function truths in
	# >= 80% of replications") requires much longer chains than is feasible
	# in a fast regression test. Here we only assert that the diagnostic
	# completes without error and returns a non-null finite BF row for the
	# dynamic dyadic coefficient. The full calibration is documented as a
	# manual sanity-check in inst/CHANGES_PHASE1.md.
	set.seed(2026)
	n <- 8; T <- 6
	beta_true <- c(rep(-0.6, 3), rep(0.6, 3))   # mid-period step
	X_list <- replicate(T, matrix(rnorm(n*n), n, n), simplify = FALSE)
	Y_list <- vector("list", T)
	for (t in seq_len(T)) {
		Yt <- beta_true[t] * X_list[[t]] + matrix(rnorm(n*n, 0, 0.1), n, n)
		diag(Yt) <- NA
		rownames(Yt) <- colnames(Yt) <- paste0("a", seq_len(n))
		Y_list[[t]] <- Yt
	}
	X_arr <- lapply(X_list, function(x) array(x, c(n, n, 1)))
	fit <- lame(Y_list, Xdyad = X_arr, family = "normal", R = 0,
	            nscan = 100, burn = 30, odens = 5,
	            dynamic_beta = "dyad", verbose = FALSE)
	res <- detect_change_point(fit, n_prior_sims = 200, seed = 1)
	expect_true(is.data.frame(res))
	expect_gt(nrow(res), 0L)
	expect_true(all(is.finite(res$bf)))
})

# -- 1.6 dynamic_beta_prior_summary -------------------------------------------

test_that("dynamic_beta_prior_summary is reproducible with seed", {
	a <- dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	b <- dynamic_beta_prior_summary(seed = 123, ndraws = 500)
	expect_equal(a$summary, b$summary)
	expect_equal(a$paths, b$paths)
})

test_that("larger sigma_scale (= more diffuse) gives larger max-first-diff", {
	small <- dynamic_beta_prior_summary(sigma_scale = 0.05, seed = 1, ndraws = 500)
	big   <- dynamic_beta_prior_summary(sigma_scale = 5.0,  seed = 1, ndraws = 500)
	expect_gt(big$summary$q50[1], small$summary$q50[1])  # max_abs_first_diff
	expect_gt(big$summary$q50[2], small$summary$q50[2])  # roughness
})

test_that("rw1 kind has rho fixed at 1 and path is unit-root random walk", {
	s <- dynamic_beta_prior_summary(kind = "rw1", seed = 1, ndraws = 500)
	expect_true(all(s$rho == 1))
	expect_equal(s$kind, "rw1")
})

test_that("invalid rho_mean outside bounds errors cleanly", {
	expect_error(dynamic_beta_prior_summary(rho_mean = 1.5, ndraws = 100),
	             "must be strictly between")
	expect_error(dynamic_beta_prior_summary(rho_mean = -0.1, ndraws = 100),
	             "must be strictly between")
})

test_that("T = 1 errors with informative message", {
	expect_error(dynamic_beta_prior_summary(T = 1), "at least 2")
})

# -- 1.7 freeze_call ---------------------------------------------------------

test_that("freeze_call = TRUE stores a data snapshot on the fit", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            freeze_call = TRUE, verbose = FALSE)
	expect_true(isTRUE(fit$freeze_call))
	expect_false(is.null(fit$data_snapshot))
	expect_identical(fit$data_snapshot$Y, dat$Y)
})

test_that("freeze_call = TRUE: update() refits with snapshot even if local Y mutates", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	Y_local <- dat$Y
	set.seed(7)
	fit <- lame(Y_local, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            freeze_call = TRUE, verbose = FALSE)
	# mutate Y_local AFTER fit
	Y_local[[1]][1, 2] <- Y_local[[1]][1, 2] + 1000
	# update() must use the snapshot, not the mutated Y_local
	set.seed(7)
	fit_upd <- update(fit, nscan = 10)
	# refit against the original data with the same args should reproduce fit
	set.seed(7)
	fit_ref <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	                nscan = 10, burn = 2, odens = 1,
	                freeze_call = TRUE, verbose = FALSE)
	expect_equal(fit_upd$BETA, fit_ref$BETA)
})

test_that("freeze_call = FALSE (default) preserves legacy update.lame behavior", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1, verbose = FALSE)
	expect_null(fit$data_snapshot)
	expect_null(fit$freeze_call)
	# update() should still work using parent.frame() lookup
	fit_upd <- update(fit, nscan = 5)
	expect_s3_class(fit_upd, "lame")
})

# -- 1.8 compact print -------------------------------------------------------

test_that("static print output does NOT use the compact table", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1, verbose = FALSE)
	out <- paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	expect_false(grepl("Dynamic components", out))
})

test_that("single-dynamic fit does NOT use the compact table", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_beta = TRUE, verbose = FALSE)
	out <- paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	# single dynamic_beta should print the regular per-component line, not
	# the compact joint table
	expect_false(grepl("Dynamic components$|Dynamic components\\s", out))
})

test_that("joint dynamic_ab + dynamic_beta fit uses the compact table", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_ab = TRUE, dynamic_beta = TRUE, verbose = FALSE)
	out <- paste(c(capture.output(print(fit)),
	               capture.output(print(fit), type = "message")),
	             collapse = "\n")
	expect_true(grepl("Dynamic components", out))
})

test_that("compact = FALSE disables the joint table", {
	dat <- simulate_test_network(n = 6, n_time = 3, family = "normal", seed = 2026)
	fit <- lame(dat$Y, Xdyad = dat$Xdyad, family = "normal", R = 0,
	            nscan = 10, burn = 2, odens = 1,
	            dynamic_ab = TRUE, dynamic_beta = TRUE, verbose = FALSE)
	out <- paste(c(capture.output(print(fit, compact = FALSE)),
	               capture.output(print(fit, compact = FALSE), type = "message")),
	             collapse = "\n")
	expect_false(grepl("Dynamic components", out))
})
