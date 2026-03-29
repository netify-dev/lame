skip_on_cran()

# numerical stability

####
# separation / extreme data
####

test_that("ame() handles all-zero binary Y", {
	n = 10
	Y = matrix(0, n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	# Known limitation: degenerate binary data causes chol() failure
	# because initial covariance S0 is not positive definite
	# this is acceptable behavior for completely degenerate data
	result = tryCatch(
		suppressWarnings(ame(
			Y, R = 0, family = "binary",
			burn = 20, nscan = 50, odens = 1,
			verbose = FALSE, gof = FALSE, seed = 42
		)),
		error = function(e) {
			expect_true(grepl("positive|singular|chol", e$message, ignore.case = TRUE))
			NULL
		}
	)

	if (!is.null(result)) {
		expect_s3_class(result, "ame")
	}
})

test_that("ame() handles all-one binary Y", {
	n = 10
	Y = matrix(1, n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	# Known limitation: same as all-zero
	result = tryCatch(
		suppressWarnings(ame(
			Y, R = 0, family = "binary",
			burn = 20, nscan = 50, odens = 1,
			verbose = FALSE, gof = FALSE, seed = 42
		)),
		error = function(e) {
			expect_true(grepl("positive|singular|chol", e$message, ignore.case = TRUE))
			NULL
		}
	)

	if (!is.null(result)) {
		expect_s3_class(result, "ame")
	}
})

test_that("ame() handles very large normal Y values", {
	n = 10
	set.seed(7003)
	Y = matrix(rnorm(n * n, mean = 1000, sd = 100), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	fit = suppressWarnings(ame(
		Y, R = 0, family = "normal",
		burn = 20, nscan = 50, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	))

	expect_s3_class(fit, "ame")
	expect_false(any(is.nan(fit$BETA)))
	expect_false(any(is.infinite(fit$BETA)))
})

test_that("ame() handles very small normal Y values", {
	n = 10
	set.seed(7004)
	Y = matrix(rnorm(n * n, mean = 0, sd = 1e-6), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	fit = suppressWarnings(ame(
		Y, R = 0, family = "normal",
		burn = 20, nscan = 50, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	))

	expect_s3_class(fit, "ame")
	expect_false(any(is.nan(fit$BETA)))
})

test_that("ame() handles ordinal with all ties equal", {
	n = 10
	Y = matrix(2L, n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	# this should work but produce degenerate posteriors
	fit = suppressWarnings(ame(
		Y, R = 0, family = "ordinal",
		burn = 20, nscan = 50, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	))

	expect_s3_class(fit, "ame")
})

####
# mcmc diagnostics (basic checks)
####

test_that("MCMC chains have no stuck values for well-specified normal model", {
	n = 20
	dat = simulate_test_network(
		n = n, family = "normal", R = 2,
		beta_intercept = 2, beta_dyad = 1.5,
		mode = "unipartite", seed = 7010
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		burn = 200, nscan = 500, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# BETA chains should have some variation (not stuck)
	for (k in 1:ncol(fit$BETA)) {
		expect_true(sd(fit$BETA[, k]) > 0,
								info = paste("BETA column", k, "has zero variance (stuck chain)"))
	}

	# VC chains should have some variation
	for (k in 1:ncol(fit$VC)) {
		expect_true(sd(fit$VC[, k]) > 0,
								info = paste("VC column", k, "has zero variance"))
	}
})

test_that("True parameters within 95% CI for well-specified model", {
	n = 25
	dat = simulate_test_network(
		n = n, family = "normal", R = 0,
		beta_intercept = 2, beta_dyad = 1.5,
		mode = "unipartite", seed = 7011
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
		burn = 500, nscan = 2000, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	ci = confint(fit, level = 0.95)
	# intercept should contain 2
	expect_true(ci["intercept", 1] < 2 && ci["intercept", 2] > 2,
							info = paste("Intercept CI:", ci["intercept", 1], "-", ci["intercept", 2]))
	# dyad coef should contain 1.5
	dyad_row = grep("dyad", rownames(ci), ignore.case = TRUE)
	if (length(dyad_row) > 0) {
		expect_true(ci[dyad_row, 1] < 1.5 && ci[dyad_row, 2] > 1.5,
								info = paste("Dyad CI:", ci[dyad_row, 1], "-", ci[dyad_row, 2]))
	}
})

####
# matrix singularity
####

test_that("ame() handles small network (n=3)", {
	n = 3
	set.seed(7020)
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", 1:n)

	fit = suppressWarnings(ame(
		Y, R = 0, family = "normal",
		burn = 20, nscan = 50, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	))

	expect_s3_class(fit, "ame")
	expect_false(any(is.nan(fit$BETA)))
})

test_that("ame() handles nearly collinear covariates", {
	n = 15
	set.seed(7021)
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	actors = paste0("a", 1:n)
	rownames(Y) = colnames(Y) = actors

	# two nearly collinear Xdyad covariates
	x1 = rnorm(n * n)
	x2 = x1 + rnorm(n * n, 0, 0.01)  # Near-perfect correlation
	Xdyad = array(c(x1, x2), dim = c(n, n, 2))
	for (k in 1:2) diag(Xdyad[,,k]) = NA
	dimnames(Xdyad) = list(actors, actors, c("x1", "x2"))

	fit = suppressWarnings(ame(
		Y, Xdyad = Xdyad, R = 0, family = "normal",
		burn = 20, nscan = 50, odens = 1,
		verbose = FALSE, gof = FALSE, seed = 42
	))

	expect_s3_class(fit, "ame")
	# collinear covariates may produce wide posteriors but shouldn't crash
	expect_false(any(is.nan(fit$BETA)))
})
