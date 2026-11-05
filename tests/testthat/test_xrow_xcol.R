skip_on_cran()

# Xrow and Xcol covariates for static gaussian models
library(lame)
library(testthat)

test_that("xrow/xcol work correctly for static unipartite Gaussian models", {
	set.seed(123)
	n = 30

	# row and column covariates
	Xrow = matrix(rnorm(n * 2), n, 2)
	colnames(Xrow) = c("row_cov1", "row_cov2")
	Xcol = matrix(rnorm(n * 2), n, 2)
	colnames(Xcol) = c("col_cov1", "col_cov2")

	beta_row = c(0.5, -0.3)
	beta_col = c(0.8, 0.4)
	intercept = 0.2

	# y from row and column effects
	Y = matrix(intercept, n, n)

	# row effects, constant across columns
	for (i in 1:n) {
		Y[i, ] = Y[i, ] + Xrow[i, 1] * beta_row[1] + Xrow[i, 2] * beta_row[2]
	}

	# column effects, constant across rows
	for (j in 1:n) {
		Y[, j] = Y[, j] + Xcol[j, 1] * beta_col[1] + Xcol[j, 2] * beta_col[2]
	}

	noise = matrix(rnorm(n * n, 0, 0.5), n, n)
	Y = Y + noise
	diag(Y) = NA

	fit = ame(
		Y,
		Xrow = Xrow,
		Xcol = Xcol,
		R = 0,
		family = "normal",
		burn = 200,
		nscan = 800,
		odens = 25,
		verbose = FALSE
	)

	# intercept + 2 row + 2 col = 5 coefficients
	expect_equal(ncol(fit$BETA), 5)

	beta_est = colMeans(fit$BETA)

	if (!is.null(colnames(fit$BETA))) {
		expect_true(
			"intercept" %in% colnames(fit$BETA) ||
				colnames(fit$BETA)[1] == "(Intercept)"
		)
		expect_true(any(grepl("row_cov1", colnames(fit$BETA))))
		expect_true(any(grepl("col_cov1", colnames(fit$BETA))))
	}

	# coefficients close to truth; order is intercept, row covariates, column covariates
	expect_equal(
		unname(beta_est[1]),
		intercept,
		tolerance = 0.5,
		label = "Intercept estimate"
	)
	expect_equal(
		unname(beta_est[2]),
		beta_row[1],
		tolerance = 0.5,
		label = "Row covariate 1 estimate"
	)
	expect_equal(
		unname(beta_est[3]),
		beta_row[2],
		tolerance = 0.5,
		label = "Row covariate 2 estimate"
	)
	expect_equal(
		unname(beta_est[4]),
		beta_col[1],
		tolerance = 0.5,
		label = "Column covariate 1 estimate"
	)
	expect_equal(
		unname(beta_est[5]),
		beta_col[2],
		tolerance = 0.5,
		label = "Column covariate 2 estimate"
	)

	# design array construction
	expect_false(is.null(fit$X))
	expect_equal(dim(fit$X), c(n, n, 5))

	# first layer is the intercept
	expect_true(all(fit$X[, , 1] == 1, na.rm = TRUE))

	# row covariate layers are constant across columns
	for (k in 2:3) {
		for (i in 1:n) {
			row_vals = fit$X[i, , k]
			expect_true(
				length(unique(row_vals[!is.na(row_vals)])) == 1,
				info = paste("Row covariate", k - 1, "not constant across columns")
			)
		}
	}

	# column covariate layers are constant across rows
	for (k in 4:5) {
		for (j in 1:n) {
			col_vals = fit$X[, j, k]
			expect_true(
				length(unique(col_vals[!is.na(col_vals)])) == 1,
				info = paste("Column covariate", k - 3, "not constant across rows")
			)
		}
	}

	pred = predict(fit, type = "link")
	expect_equal(dim(pred), c(n, n))

	# predictions correlate with the true linear predictor
	# the linear predictor that generated Y
	true_eta = Y - noise
	corr = cor(c(true_eta), c(pred), use = "pairwise.complete.obs")
	expect_gt(corr, 0.95, label = "Predictions should correlate with truth")
})

test_that("xrow/xcol work correctly for static bipartite Gaussian models", {
	set.seed(456)
	nA = 20 # row nodes
	nB = 15 # column nodes

	# row covariates (set a)
	Xrow = matrix(rnorm(nA * 2), nA, 2)
	colnames(Xrow) = c("rowA_cov1", "rowA_cov2")

	# column covariates (set b)
	Xcol = matrix(rnorm(nB * 3), nB, 3)
	colnames(Xcol) = c("colB_cov1", "colB_cov2", "colB_cov3")

	beta_row = c(0.6, -0.4)
	beta_col = c(0.3, 0.7, -0.2)
	intercept = -0.1

	Y = matrix(intercept, nA, nB)

	# row effects from set a
	for (i in 1:nA) {
		Y[i, ] = Y[i, ] + Xrow[i, 1] * beta_row[1] + Xrow[i, 2] * beta_row[2]
	}

	# column effects from set b
	for (j in 1:nB) {
		Y[, j] = Y[, j] + Xcol[j, 1] * beta_col[1] +
			Xcol[j, 2] * beta_col[2] + Xcol[j, 3] * beta_col[3]
	}

	noise = matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	Y = Y + noise

	fit = ame(
		Y,
		Xrow = Xrow,
		Xcol = Xcol,
		mode = "bipartite",
		R_row = 0,
		R_col = 0,
		family = "normal",
		burn = 200,
		nscan = 800,
		odens = 25,
		verbose = FALSE
	)

	# intercept + 2 row + 3 col = 6 coefficients
	expect_equal(ncol(fit$BETA), 6)

	# unname to avoid name mismatch in expect_equal
	beta_est = unname(colMeans(fit$BETA))

	# recovery, with wider tolerance for bipartite
	expect_equal(
		beta_est[1],
		intercept,
		tolerance = 0.5,
		label = "Bipartite intercept"
	)
	expect_equal(
		beta_est[2],
		beta_row[1],
		tolerance = 0.5,
		label = "Bipartite row cov 1"
	)
	expect_equal(
		beta_est[3],
		beta_row[2],
		tolerance = 0.5,
		label = "Bipartite row cov 2"
	)
	expect_equal(
		beta_est[4],
		beta_col[1],
		tolerance = 0.5,
		label = "Bipartite col cov 1"
	)
	expect_equal(
		beta_est[5],
		beta_col[2],
		tolerance = 0.5,
		label = "Bipartite col cov 2"
	)
	expect_equal(
		beta_est[6],
		beta_col[3],
		tolerance = 0.5,
		label = "Bipartite col cov 3"
	)

	# design array structure
	expect_false(is.null(fit$X))
	expect_equal(dim(fit$X), c(nA, nB, 6))

	# bipartite-specific attributes
	expect_equal(fit$nA, nA)
	expect_equal(fit$nB, nB)
	expect_equal(fit$mode, "bipartite")

	EZ = reconstruct_EZ(fit)
	expect_equal(dim(EZ), c(nA, nB))

	# correlation with truth
	# the linear predictor that generated Y
	true_eta = Y - noise
	corr = cor(c(true_eta), c(EZ), use = "pairwise.complete.obs")
	expect_gt(corr, 0.95)
})

test_that("Mixed dyadic and nodal covariates work together", {
	set.seed(789)
	n = 25

	# all three covariate types
	Xdyad = matrix(rnorm(n * n), n, n)
	diag(Xdyad) = NA
	Xrow = matrix(rnorm(n), n, 1)
	Xcol = matrix(rnorm(n), n, 1)

	beta_dyad = 1.0
	beta_row = 0.5
	beta_col = -0.5
	intercept = 0.0

	Y = intercept + beta_dyad * Xdyad

	for (i in 1:n) {
		Y[i, ] = Y[i, ] + Xrow[i] * beta_row
	}

	for (j in 1:n) {
		Y[, j] = Y[, j] + Xcol[j] * beta_col
	}

	Y = Y + matrix(rnorm(n * n, 0, 0.3), n, n)
	diag(Y) = NA

	fit = ame(
		Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		R = 0,
		family = "normal",
		burn = 200,
		nscan = 800,
		odens = 25,
		verbose = FALSE
	)

	# intercept, row, column and dyad coefficients
	expect_equal(ncol(fit$BETA), 4)

	beta_est = colMeans(fit$BETA)

	# intercept comes first
	expect_equal(unname(beta_est[1]), intercept, tolerance = 0.5)
	# the remaining coefficients may come in any order, so match each true value to one estimate
	all_true = c(beta_dyad, beta_row, beta_col)
	all_est = unname(beta_est[2:4])

	expect_true(
		any(abs(all_est - beta_dyad) < 0.3),
		info = "Dyadic coefficient not recovered"
	)
	expect_true(
		any(abs(all_est - beta_row) < 0.5),
		info = "Row coefficient not recovered"
	)
	expect_true(
		any(abs(all_est - beta_col) < 0.5),
		info = "Column coefficient not recovered"
	)

	# design array dimensions
	expect_equal(dim(fit$X), c(n, n, 4))
})
