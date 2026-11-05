# symmetric ordinal in ame(): the path runs, returns a symmetric YPM, recovers a dyadic-covariate sign, and uses
# variance 1 (each unordered dyad once) so the dyadic coefficients are not attenuated
skip_on_cran()

skip_if_too_slow = function() {
	if (identical(Sys.getenv("LAME_FAST_TESTS"), "1")) skip("LAME_FAST_TESTS = 1")
}

test_that("ame(symmetric = TRUE, family = 'ordinal') runs end-to-end", {
	skip_if_too_slow()
	set.seed(5201L)
	n = 30
	X1 = matrix(stats::rnorm(n * n), n, n)
	X1 = (X1 + t(X1)) / 2
	diag(X1) = 0
	dimnames(X1) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))

	beta_true = 0.8
	a = stats::rnorm(n, 0, 0.3)
	eta = beta_true * X1 + outer(a, a, "+")
	Zsym = eta + matrix(stats::rnorm(n * n), n, n)
	Zsym = (Zsym + t(Zsym)) / 2

	Y = as.integer(cut(Zsym, c(-Inf, -0.5, 0.3, 1.0, Inf), labels = FALSE))
	dim(Y) = c(n, n)
	dimnames(Y) = dimnames(X1)
	# symmetry from the upper triangle
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA

	Xa = array(X1, c(n, n, 1))
	dimnames(Xa) = list(NULL, NULL, "Xdyad")

	fit = ame(
		Y,
		Xdyad = Xa,
		family = "ordinal",
		symmetric = TRUE,
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5201L
	)

	expect_true(inherits(fit, "ame"))
	# BETA has one column (the dyadic covariate)
	expect_equal(ncol(fit$BETA), 1L)
	expect_true("Xdyad_dyad" %in% colnames(fit$BETA))
	# sign recovery
	expect_gt(mean(fit$BETA[, "Xdyad_dyad"]), 0.3)
})

test_that("symmetric ordinal posterior YPM (and EZ) is symmetric", {
	skip_if_too_slow()
	set.seed(5202L)
	n = 25
	X1 = matrix(stats::rnorm(n * n), n, n)
	X1 = (X1 + t(X1)) / 2
	diag(X1) = 0
	dimnames(X1) = list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))

	beta_true = 0.6
	a = stats::rnorm(n, 0, 0.2)
	eta = beta_true * X1 + outer(a, a, "+")
	Zsym = eta + matrix(stats::rnorm(n * n), n, n)
	Zsym = (Zsym + t(Zsym)) / 2
	Y = as.integer(cut(Zsym, c(-Inf, -0.4, 0.5, Inf), labels = FALSE))
	dim(Y) = c(n, n)
	dimnames(Y) = dimnames(X1)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA

	Xa = array(X1, c(n, n, 1))
	dimnames(Xa) = list(NULL, NULL, "Xdyad")

	fit = ame(
		Y,
		Xdyad = Xa,
		family = "ordinal",
		symmetric = TRUE,
		R = 0,
		nscan = 150,
		burn = 30,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5202L
	)

	# the stored EZ is symmetric to numerical precision
	if (!is.null(fit$EZ)) {
		expect_true(isTRUE(all.equal(
			fit$EZ,
			t(fit$EZ),
			check.attributes = FALSE,
			tolerance = 1e-8
		)))
	}
	# YPM is the posterior-predictive mean of Y; for ordinal a probability-weighted category mean, still symmetric
	ypm = fit$YPM
	if (is.matrix(ypm)) {
		# copy the diagonal so the symmetry check ignores it
		diag(ypm) = diag(t(ypm))
		expect_true(isTRUE(all.equal(
			ypm,
			t(ypm),
			check.attributes = FALSE,
			tolerance = 1e-6
		)))
	}
})

test_that("rZ_ord_sym_fc stays finite under extreme predictors and an inverted start", {
	set.seed(6886)
	n = 10
	Y = matrix(sample(1:3, n * n, replace = TRUE), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	# a block of actors far below every cutpoint saturates the normal CDF
	EZ = matrix(0, n, n)
	EZ[1:3, ] = -30
	EZ[, 1:3] = -30
	# a start whose ranks disagree with Y inverts the truncation window
	Z = matrix(runif(n * n, 2, 3), n, n)
	Z[lower.tri(Z)] = t(Z)[lower.tri(Z)]
	for (k in 1:30) Z = rZ_ord_sym_fc(Z, EZ, Y)
	expect_true(all(is.finite(Z)))
	expect_equal(Z, t(Z))
})
