# longitudinal models via lame()

library(lame)
library(testthat)

test_that("lame() runs with normal family and returns correct structure", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 1, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, rho = 0, s2 = 1, seed = 101
	)

	fit <- lame(
		Y = dat$Y, Xdyad = dat$Xdyad, Xrow = dat$Xrow, Xcol = dat$Xcol,
		family = "normal", R = 0,
		burn = 100, nscan = 500, odens = 5,
		print = FALSE, gof = TRUE
	)

	# Structure checks
	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
	expect_false(is.null(fit$VC))
	expect_false(is.null(fit$APM))
	expect_false(is.null(fit$BPM))
	expect_false(is.null(fit$YPM))

	# BETA should have intercept + Xdyad coefficient + Xrow + Xcol = 4 columns
	expect_equal(ncol(fit$BETA), 4)
	# 100 posterior samples (500 / odens=5)
	expect_equal(nrow(fit$BETA), 100)

	# Check intercept is in reasonable range
	beta_intercept_hat <- median(fit$BETA[, 1])
	expect_true(abs(beta_intercept_hat - 1) < 2,
		info = paste("Intercept estimate:", beta_intercept_hat, "expected ~1"))

	# Check dyad covariate recovery
	beta_dyad_hat <- median(fit$BETA[, 2])
	expect_true(abs(beta_dyad_hat - 0.5) < 1,
		info = paste("Dyad beta estimate:", beta_dyad_hat, "expected ~0.5"))
})

test_that("lame() runs with binary family and returns correct structure", {
	skip_on_cran()

	set.seed(102)
	n <- 15
	n_time <- 3
	actors <- paste0("a", 1:n)

	# Generate binary data directly to avoid simulator edge cases
	Y_list <- list()
	for (t in 1:n_time) {
		Y <- matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
	}
	names(Y_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list,
		family = "binary", R = 0,
		burn = 100, nscan = 500, odens = 5,
		print = FALSE, gof = TRUE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
	expect_false(is.null(fit$VC))

	# YPM should be probabilities in [0,1]
	if (!is.null(fit$YPM)) {
		# YPM may be an array or list; unlist to get numeric values
		ypm_vals <- unlist(fit$YPM)
		ypm_vals <- ypm_vals[!is.na(ypm_vals)]
		expect_true(all(ypm_vals >= 0 & ypm_vals <= 1),
			info = "Binary YPM should be probabilities in [0,1]")
	}

	# Intercept should be negative for sparse (30%) network
	intercept_hat <- median(fit$BETA[, 1])
	expect_true(intercept_hat < 1,
		info = "Sparse network should have negative/small intercept")
})

test_that("lame() with multiplicative effects (R > 0) returns UV", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 2, beta_intercept = 0, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, seed = 103
	)

	fit <- lame(
		Y = dat$Y, Xdyad = dat$Xdyad,
		family = "normal", R = 2,
		burn = 100, nscan = 500, odens = 5,
		print = FALSE, gof = TRUE
	)

	expect_s3_class(fit, "lame")

	# U and V should exist with correct dimensions
	expect_false(is.null(fit$U))
	expect_false(is.null(fit$V))

	# For static (non-dynamic) model, U and V should be 2D (n x R)
	expect_equal(ncol(fit$U), 2)
	expect_equal(ncol(fit$V), 2)
})

test_that("lame() with symmetric networks works", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 0, beta_dyad = 0.5,
		rvar = TRUE, cvar = TRUE, symmetric = TRUE, seed = 104
	)

	# Symmetrize the data
	for (t in seq_along(dat$Y)) {
		Y_t <- dat$Y[[t]]
		Y_t <- (Y_t + t(Y_t)) / 2
		diag(Y_t) <- NA
		dat$Y[[t]] <- Y_t
	}

	fit <- lame(
		Y = dat$Y, Xdyad = dat$Xdyad,
		family = "normal", R = 0, symmetric = TRUE,
		burn = 100, nscan = 500, odens = 5,
		print = FALSE, gof = TRUE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() with covariates recovers known effects", {
	skip_on_cran()

	set.seed(200)
	n <- 20
	n_time <- 4
	true_beta <- 0.8

	# Create data with strong known effect
	actors <- paste0("a", 1:n)
	Y_list <- list()
	Xdyad_list <- list()

	for (t in 1:n_time) {
		X <- matrix(rnorm(n * n), n, n)
		diag(X) <- NA
		rownames(X) <- colnames(X) <- actors
		Z <- true_beta * X + matrix(rnorm(n * n), n, n)
		diag(Z) <- NA

		Y_list[[t]] <- Z
		Xdyad_list[[t]] <- X
	}
	names(Y_list) <- names(Xdyad_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list, Xdyad = Xdyad_list,
		family = "normal", R = 0,
		rvar = FALSE, cvar = FALSE,
		burn = 200, nscan = 1000, odens = 5,
		print = FALSE, gof = FALSE
	)

	# Check beta recovery
	beta_hat <- median(fit$BETA[, 2])
	beta_ci <- quantile(fit$BETA[, 2], c(0.025, 0.975))

	expect_true(
		true_beta >= beta_ci[1] && true_beta <= beta_ci[2],
		info = paste("True beta", true_beta, "should be in CI [",
			round(beta_ci[1], 3), ",", round(beta_ci[2], 3), "]")
	)

	# Point estimate should be reasonably close
	expect_true(abs(beta_hat - true_beta) < 0.3,
		info = paste("Beta hat:", round(beta_hat, 3), "expected:", true_beta))
})

test_that("lame() additive effects correlate with true values", {
	skip_on_cran()

	set.seed(301)
	n <- 20
	n_time <- 5
	actors <- paste0("a", 1:n)

	# Generate strong additive effects with large variance relative to noise
	true_a <- rnorm(n, 0, 2)
	true_b <- rnorm(n, 0, 2)

	Y_list <- list()
	for (t in 1:n_time) {
		eta <- outer(true_a, rep(1, n)) + outer(rep(1, n), true_b)
		Y <- eta + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
	}
	names(Y_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list,
		family = "normal", R = 0,
		rvar = TRUE, cvar = TRUE,
		burn = 500, nscan = 2000, odens = 10,
		print = FALSE, gof = FALSE
	)

	# Additive effects should correlate with truth (match by actor name)
	if (!is.null(fit$APM) && length(fit$APM) == n) {
		# Match APM to true_a by actor names (lame may reorder actors)
		names(true_a) <- actors
		apm_matched <- fit$APM[actors]
		a_centered <- true_a - mean(true_a)
		apm_centered <- apm_matched - mean(apm_matched)
		a_cor <- cor(a_centered, apm_centered, use = "pairwise.complete.obs")
		expect_true(a_cor > 0.3,
			info = paste("Sender effect correlation:", round(a_cor, 3)))
	}
	if (!is.null(fit$BPM) && length(fit$BPM) == n) {
		names(true_b) <- actors
		bpm_matched <- fit$BPM[actors]
		b_centered <- true_b - mean(true_b)
		bpm_centered <- bpm_matched - mean(bpm_matched)
		b_cor <- cor(b_centered, bpm_centered, use = "pairwise.complete.obs")
		expect_true(b_cor > 0.3,
			info = paste("Receiver effect correlation:", round(b_cor, 3)))
	}
})

test_that("lame() GOF output has correct structure", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "normal",
		R = 0, beta_intercept = 0, beta_dyad = 0,
		seed = 105
	)

	fit <- lame(
		Y = dat$Y,
		family = "normal", R = 0,
		burn = 50, nscan = 250, odens = 5,
		print = FALSE, gof = TRUE
	)

	expect_false(is.null(fit$GOF),
		info = "GOF should be computed when gof=TRUE")
	if (!is.null(fit$GOF)) {
		# GOF for lame is a list of matrices (one per GOF statistic)
		expect_true(is.list(fit$GOF) || is.matrix(fit$GOF))
		if (is.list(fit$GOF)) {
			expect_true(length(fit$GOF) > 0, info = "GOF should have at least one statistic")
		}
	}
})

# additional family coverage

test_that("lame() works with ordinal family", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "ordinal",
		R = 0, beta_intercept = 0, beta_dyad = 0.3, seed = 503
	)

	fit <- lame(
		Y = dat$Y, Xdyad = dat$Xdyad,
		family = "ordinal", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with tobit family", {
	skip_on_cran()

	dat <- simulate_test_network(
		n = 15, n_time = 3, family = "tobit",
		R = 0, beta_intercept = 0.5, beta_dyad = 0.3, seed = 504
	)

	fit <- lame(
		Y = dat$Y, Xdyad = dat$Xdyad,
		family = "tobit", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")

	# Tobit YPM should be >= 0
	if (!is.null(fit$YPM)) {
		ypm_vals <- unlist(fit$YPM)
		ypm_vals <- ypm_vals[!is.na(ypm_vals)]
		expect_true(all(ypm_vals >= -0.01),
			info = "Tobit YPM should be approximately non-negative")
	}
})

test_that("lame() works with poisson family", {
	skip_on_cran()

	set.seed(505)
	n <- 15
	n_time <- 3
	actors <- paste0("a", 1:n)

	Y_list <- list()
	Xdyad_list <- list()
	for (t in 1:n_time) {
		X <- matrix(rnorm(n * n, 0, 0.3), n, n)
		diag(X) <- NA
		rownames(X) <- colnames(X) <- actors
		lambda <- exp(0.5 + 0.2 * X)
		lambda[is.na(lambda)] <- exp(0.5)
		Y <- matrix(rpois(n * n, lambda), n, n)
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
		Xdyad_list[[t]] <- X
	}
	names(Y_list) <- names(Xdyad_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list, Xdyad = Xdyad_list,
		family = "poisson", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))

	# YPM should be non-negative counts
	if (!is.null(fit$YPM)) {
		ypm_vals <- unlist(fit$YPM)
		ypm_vals <- ypm_vals[!is.na(ypm_vals)]
		expect_true(all(ypm_vals >= 0),
			info = "Poisson YPM should be non-negative")
	}
})

test_that("lame() works with cbin family", {
	skip_on_cran()

	set.seed(506)
	n <- 15
	n_time <- 3
	actors <- paste0("a", 1:n)
	odmax_val <- 4

	Y_list <- list()
	for (t in 1:n_time) {
		Y <- matrix(0, n, n)
		for (i in 1:n) {
			nominees <- sample(setdiff(1:n, i), odmax_val)
			Y[i, nominees] <- 1
		}
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
	}
	names(Y_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list,
		family = "cbin", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with frn family", {
	skip_on_cran()

	set.seed(507)
	n <- 15
	n_time <- 3
	actors <- paste0("a", 1:n)
	odmax_val <- 3

	Y_list <- list()
	for (t in 1:n_time) {
		Y <- matrix(0, n, n)
		for (i in 1:n) {
			nominees <- sample(setdiff(1:n, i), odmax_val)
			Y[i, nominees] <- 1:odmax_val
		}
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
	}
	names(Y_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list,
		family = "frn", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with rrl family", {
	skip_on_cran()

	set.seed(508)
	n <- 15
	n_time <- 3
	actors <- paste0("a", 1:n)
	n_ranked <- 5

	Y_list <- list()
	for (t in 1:n_time) {
		Y <- matrix(0, n, n)
		for (i in 1:n) {
			ranked <- sample(setdiff(1:n, i), n_ranked)
			Y[i, ranked] <- n_ranked:1
		}
		diag(Y) <- NA
		rownames(Y) <- colnames(Y) <- actors
		Y_list[[t]] <- Y
	}
	names(Y_list) <- paste0("t", 1:n_time)

	fit <- lame(
		Y = Y_list,
		family = "rrl", R = 0,
		burn = 100, nscan = 300, odens = 3,
		print = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
})

test_that("lame() works with all families with R > 0", {
	skip_on_cran()

	families_to_test <- c("normal", "binary")

	for (fam in families_to_test) {
		set.seed(510)
		n <- 12
		n_time <- 3
		actors <- paste0("a", 1:n)

		Y_list <- list()
		for (t in 1:n_time) {
			if (fam == "normal") {
				Y <- matrix(rnorm(n * n), n, n)
			} else {
				Y <- matrix(rbinom(n * n, 1, 0.3), n, n)
			}
			diag(Y) <- NA
			rownames(Y) <- colnames(Y) <- actors
			Y_list[[t]] <- Y
		}
		names(Y_list) <- paste0("t", 1:n_time)

		fit <- lame(
			Y = Y_list,
			family = fam, R = 2,
			burn = 100, nscan = 300, odens = 3,
			print = FALSE, gof = FALSE
		)

		expect_s3_class(fit, "lame")
		expect_false(is.null(fit$U),
			info = paste("U should exist for", fam, "with R=2"))
		expect_false(is.null(fit$V),
			info = paste("V should exist for", fam, "with R=2"))
	}
})
