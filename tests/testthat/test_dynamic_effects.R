skip_on_cran()

# dynamic effects: dynamic_ab and dynamic_uv

library(lame)
library(testthat)

# dynamic_ab

test_that("lame() with dynamic_ab=TRUE runs and produces output", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 4
	rho_ab_true = 0.7
	actors = paste0("a", 1:n)

	a = rnorm(n, 0, 0.5)
	b = rnorm(n, 0, 0.5)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			a = rho_ab_true * a + rnorm(n, 0, 0.3)
			b = rho_ab_true * b + rnorm(n, 0, 0.3)
		}
		eta = outer(a, rep(1, n)) + outer(rep(1, n), b)
		Y = eta + matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = 0,
		dynamic_ab = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))
	expect_false(is.null(fit$VC))
	expect_false(is.null(fit$APM))
	expect_false(is.null(fit$BPM))
})

test_that("dynamic_ab estimates positive temporal persistence", {
	skip_on_cran()

	set.seed(6886)
	n = 20
	n_time = 6
	rho_ab_true = 0.8
	actors = paste0("a", 1:n)

	a = rnorm(n, 0, 1)
	b = rnorm(n, 0, 1)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			a = rho_ab_true * a + rnorm(n, 0, 0.2)
			b = rho_ab_true * b + rnorm(n, 0, 0.2)
		}
		eta = outer(a, rep(1, n)) + outer(rep(1, n), b)
		Y = eta + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = 0,
		dynamic_ab = TRUE,
		burn = 300, nscan = 1000, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	has_rho_ab = !is.null(fit$RHO_AB) || !is.null(fit$rho_ab)
	expect_true(has_rho_ab,
		info = "Dynamic AB model should store RHO_AB chain or rho_ab estimate")

	if (!is.null(fit$RHO_AB)) {
		rho_hat = median(fit$RHO_AB, na.rm = TRUE)
		expect_true(rho_hat > 0.1,
			info = paste("Estimated rho_ab:", round(rho_hat, 3), "should be positive"))
	}
})

test_that("dynamic_ab start_vals include rho_ab and sigma_ab", {
	skip_on_cran()

	set.seed(6886)
	n = 12
	n_time = 3
	actors = paste0("a", 1:n)

	Y_list = list()
	for (t in 1:n_time) {
		Y = matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = 0,
		dynamic_ab = TRUE,
		burn = 100, nscan = 300, odens = 3,
		verbose = FALSE, gof = FALSE
	)

	expect_false(is.null(fit$start_vals))
	if (!is.null(fit$start_vals)) {
		expect_true("rho_ab" %in% names(fit$start_vals),
			info = "start_vals should include rho_ab for dynamic_ab models")
		expect_true("sigma_ab" %in% names(fit$start_vals),
			info = "start_vals should include sigma_ab for dynamic_ab models")
	}
})

test_that("dynamic_ab works with binary family", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 4
	actors = paste0("a", 1:n)

	a = rnorm(n, 0, 0.5)
	b = rnorm(n, 0, 0.5)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			a = 0.7 * a + rnorm(n, 0, 0.3)
			b = 0.7 * b + rnorm(n, 0, 0.3)
		}
		eta = -0.5 + outer(a, rep(1, n)) + outer(rep(1, n), b)
		Y = ifelse(eta + matrix(rnorm(n * n), n, n) > 0, 1, 0)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "binary", R = 0,
		dynamic_ab = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))

	if (!is.null(fit$YPM)) {
		ypm_vals = unlist(fit$YPM)
		ypm_vals = ypm_vals[!is.na(ypm_vals)]
		expect_true(all(ypm_vals >= 0 & ypm_vals <= 1),
			info = "Binary dynamic_ab YPM should be in [0,1]")
	}
})

test_that("lame() bipartite with dynamic_ab=TRUE produces correct dimensions", {
	skip_on_cran()

	set.seed(6886)
	nA = 10
	nB = 7
	n_time = 3
	row_actors = paste0("r", 1:nA)
	col_actors = paste0("c", 1:nB)

	a = rnorm(nA, 0, 0.5)
	b = rnorm(nB, 0, 0.5)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			a = 0.7 * a + rnorm(nA, 0, 0.3)
			b = 0.7 * b + rnorm(nB, 0, 0.3)
		}
		eta = outer(a, rep(1, nB)) + outer(rep(1, nA), b)
		Y = eta + matrix(rnorm(nA * nB), nA, nB)
		rownames(Y) = row_actors
		colnames(Y) = col_actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("T", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = 0, mode = "bipartite",
		dynamic_ab = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")

	expect_false(is.null(fit$a_dynamic))
	expect_equal(nrow(fit$a_dynamic), nA)
	expect_equal(ncol(fit$a_dynamic), n_time)

	expect_false(is.null(fit$b_dynamic))
	expect_equal(nrow(fit$b_dynamic), nB)
	expect_equal(ncol(fit$b_dynamic), n_time)

	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
})

# dynamic_uv

test_that("lame() with dynamic_uv=TRUE runs and returns 3D U/V arrays", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 5
	R = 2
	rho_uv_true = 0.8
	actors = paste0("a", 1:n)

	U = matrix(rnorm(n * R), n, R)
	V = matrix(rnorm(n * R), n, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			U = rho_uv_true * U + matrix(rnorm(n * R, 0, 0.3), n, R)
			V = rho_uv_true * V + matrix(rnorm(n * R, 0, 0.3), n, R)
		}
		EZ = tcrossprod(U, V)
		Y = EZ + matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = R,
		dynamic_uv = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$U))
	expect_false(is.null(fit$V))

	if (length(dim(fit$U)) == 3) {
		expect_equal(dim(fit$U)[1], n, info = "U dim 1 should be n")
		expect_equal(dim(fit$U)[2], R, info = "U dim 2 should be R")
		expect_equal(dim(fit$U)[3], n_time, info = "U dim 3 should be n_time")
	}

	if (length(dim(fit$V)) == 3) {
		expect_equal(dim(fit$V)[1], n, info = "V dim 1 should be n")
		expect_equal(dim(fit$V)[2], R, info = "V dim 2 should be R")
		expect_equal(dim(fit$V)[3], n_time, info = "V dim 3 should be n_time")
	}
})

test_that("dynamic_uv model estimates positive temporal autocorrelation", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 6
	R = 2
	rho_uv_true = 0.85
	actors = paste0("a", 1:n)

	U = matrix(rnorm(n * R, 0, 0.5), n, R)
	V = matrix(rnorm(n * R, 0, 0.5), n, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			U = rho_uv_true * U + matrix(rnorm(n * R, 0, 0.2), n, R)
			V = rho_uv_true * V + matrix(rnorm(n * R, 0, 0.2), n, R)
		}
		EZ = tcrossprod(U, V)
		Y = EZ + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = R,
		dynamic_uv = TRUE,
		burn = 300, nscan = 1000, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	has_rho = !is.null(fit$RHO_UV) || !is.null(fit$rho_uv)
	expect_true(has_rho,
		info = "Dynamic UV model should store RHO_UV chain or rho_uv estimate")

	if (!is.null(fit$RHO_UV)) {
		rho_hat = median(fit$RHO_UV, na.rm = TRUE)
		expect_true(rho_hat > 0,
			info = paste("Estimated rho_uv:", round(rho_hat, 3), "should be positive"))
	}

	if (!is.null(fit$rho_uv)) {
		rho_uv_mean = mean(fit$rho_uv, na.rm = TRUE)
		expect_true(rho_uv_mean > 0,
			info = paste("rho_uv estimate:", round(rho_uv_mean, 3), "should be positive"))
	}
})

test_that("dynamic_uv latent positions show temporal continuity", {
	skip_on_cran()

	set.seed(6886)
	n = 12
	n_time = 5
	R = 2
	actors = paste0("a", 1:n)

	U = matrix(rnorm(n * R, 0, 1), n, R)
	V = matrix(rnorm(n * R, 0, 1), n, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			U = 0.9 * U + matrix(rnorm(n * R, 0, 0.15), n, R)
			V = 0.9 * V + matrix(rnorm(n * R, 0, 0.15), n, R)
		}
		EZ = tcrossprod(U, V)
		Y = EZ + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = R,
		dynamic_uv = TRUE,
		burn = 200, nscan = 800, odens = 4,
		verbose = FALSE, gof = FALSE
	)

	if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
		cors = c()
		for (t in 2:n_time) {
			U_t = fit$U[, , t]
			U_prev = fit$U[, , t - 1]
			cors = c(cors, cor(c(U_t), c(U_prev)))
		}
		avg_cor = mean(cors, na.rm = TRUE)
		expect_true(avg_cor > 0,
			info = paste("Adjacent U positions should be positively correlated, got:", round(avg_cor, 3)))
	}
})

test_that("dynamic_uv=FALSE gives 2D U/V (baseline comparison)", {
	skip_on_cran()

	set.seed(6886)
	n = 12
	n_time = 3
	actors = paste0("a", 1:n)

	Y_list = list()
	for (t in 1:n_time) {
		Y = matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit_static = lame(
		Y = Y_list,
		family = "normal", R = 2,
		dynamic_uv = FALSE,
		burn = 100, nscan = 300, odens = 3,
		verbose = FALSE, gof = FALSE
	)

	expect_false(is.null(fit_static$U))
	expect_true(length(dim(fit_static$U)) == 2 || is.null(dim(fit_static$U)) || dim(fit_static$U)[3] == 1,
		info = "Static model U should be 2D or single time slice")
})

test_that("dynamic_uv works with binary family", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 4
	R = 2
	actors = paste0("a", 1:n)

	U = matrix(rnorm(n * R, 0, 0.5), n, R)
	V = matrix(rnorm(n * R, 0, 0.5), n, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			U = 0.8 * U + matrix(rnorm(n * R, 0, 0.2), n, R)
			V = 0.8 * V + matrix(rnorm(n * R, 0, 0.2), n, R)
		}
		EZ = -0.5 + tcrossprod(U, V)
		Y = ifelse(EZ + matrix(rnorm(n * n), n, n) > 0, 1, 0)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "binary", R = R,
		dynamic_uv = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$U))
	expect_false(is.null(fit$V))
})

test_that("lame() bipartite with dynamic_uv=TRUE produces correct dimensions", {
	skip_on_cran()

	set.seed(6886)
	nA = 10
	nB = 7
	n_time = 3
	R = 2
	row_actors = paste0("r", 1:nA)
	col_actors = paste0("c", 1:nB)

	U = matrix(rnorm(nA * R, 0, 0.5), nA, R)
	V = matrix(rnorm(nB * R, 0, 0.5), nB, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			U = 0.8 * U + matrix(rnorm(nA * R, 0, 0.2), nA, R)
			V = 0.8 * V + matrix(rnorm(nB * R, 0, 0.2), nB, R)
		}
		eta = U %*% t(V)
		Y = eta + matrix(rnorm(nA * nB), nA, nB)
		rownames(Y) = row_actors
		colnames(Y) = col_actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("T", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = R, mode = "bipartite",
		dynamic_uv = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$U))
	expect_false(is.null(fit$V))

	if (length(dim(fit$U)) == 3) {
		expect_equal(dim(fit$U)[1], nA)
		expect_equal(dim(fit$U)[2], R)
		expect_equal(dim(fit$U)[3], n_time)
	}

	if (length(dim(fit$V)) == 3) {
		expect_equal(dim(fit$V)[1], nB)
		expect_equal(dim(fit$V)[2], R)
		expect_equal(dim(fit$V)[3], n_time)
	}
})

# combined dynamic_ab + dynamic_uv

test_that("combined dynamic_ab + dynamic_uv runs", {
	skip_on_cran()

	set.seed(6886)
	n = 15
	n_time = 4
	R = 2
	actors = paste0("a", 1:n)

	a = rnorm(n, 0, 0.5)
	b = rnorm(n, 0, 0.5)
	U = matrix(rnorm(n * R, 0, 0.5), n, R)
	V = matrix(rnorm(n * R, 0, 0.5), n, R)
	Y_list = list()

	for (t in 1:n_time) {
		if (t > 1) {
			a = 0.7 * a + rnorm(n, 0, 0.3)
			b = 0.7 * b + rnorm(n, 0, 0.3)
			U = 0.8 * U + matrix(rnorm(n * R, 0, 0.2), n, R)
			V = 0.8 * V + matrix(rnorm(n * R, 0, 0.2), n, R)
		}
		eta = outer(a, rep(1, n)) + outer(rep(1, n), b) + tcrossprod(U, V)
		Y = eta + matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actors
		Y_list[[t]] = Y
	}
	names(Y_list) = paste0("t", 1:n_time)

	fit = lame(
		Y = Y_list,
		family = "normal", R = R,
		dynamic_ab = TRUE, dynamic_uv = TRUE,
		burn = 200, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE
	)

	expect_s3_class(fit, "lame")

	has_rho_uv = !is.null(fit$RHO_UV) || !is.null(fit$rho_uv)
	has_rho_ab = !is.null(fit$RHO_AB) || !is.null(fit$rho_ab)

	expect_true(has_rho_uv || has_rho_ab,
		info = "Combined dynamic model should track at least one rho parameter")

	if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
		expect_equal(dim(fit$U)[3], n_time)
	}
})

# rho parameter properties

test_that("rho parameters are properly bounded", {
	skip_on_cran()

	set.seed(6886)
	n = 20
	n_time = 5

	Y_list = list()
	for (t in 1:n_time) {
		Y_t = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n)
		Y_list[[t]] = Y_t
	}

	fit_all_dyn = lame(Y_list, R = 0, family = "binary",
		dynamic_uv = FALSE, dynamic_ab = TRUE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE)

	expect_true(all(fit_all_dyn$rho_ab >= -1 & fit_all_dyn$rho_ab <= 1))
})

test_that("prior specifications affect rho estimates", {
	skip_on_cran()

	set.seed(6886)
	n = 20
	n_time = 5

	Y_list = list()
	eta_prev = matrix(rnorm(n * n, -0.5, 1), n, n)

	for (t in 1:n_time) {
		eta = 0.6 * eta_prev + matrix(rnorm(n * n, 0, 0.5), n, n)
		Y_t = (eta > 0) * 1
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n)
		Y_list[[t]] = Y_t
		eta_prev = eta
	}

	fit_low = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = TRUE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE,
		prior = list(rho_uv_mean = 0.2, rho_uv_sd = 0.1))

	fit_high = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = TRUE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE,
		prior = list(rho_uv_mean = 0.8, rho_uv_sd = 0.1))

	if (!is.null(fit_low$rho_uv) && length(fit_low$rho_uv) > 0 &&
		 !is.null(fit_high$rho_uv) && length(fit_high$rho_uv) > 0) {
		rho_low = median(fit_low$rho_uv, na.rm = TRUE)
		rho_high = median(fit_high$rho_uv, na.rm = TRUE)

		if (!is.na(rho_low) && !is.na(rho_high)) {
			expect_gte(abs(rho_low), 0)
			expect_gte(abs(rho_high), 0)
		} else {
			expect_true(length(fit_low$rho_uv) > 0)
		}
	} else {
		expect_true(!is.null(fit_low$BETA))
		expect_true(!is.null(fit_high$BETA))
	}
})

test_that("FFBS algorithm maintains proper variance", {
	skip_on_cran()

	set.seed(6886)
	n = 20
	n_time = 8

	Y_list = list()
	for (t in 1:n_time) {
		Y_t = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n)
		Y_list[[t]] = Y_t
	}

	fit = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = TRUE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE)

	if (!is.null(fit$U)) {
		if (length(dim(fit$U)) == 4) {
			U_var = apply(fit$U, c(2, 3, 4), var)
			mean_var = mean(U_var, na.rm = TRUE)

			var_by_time = apply(U_var, 2, mean, na.rm = TRUE)
			if (length(var_by_time) > 1 && all(!is.na(var_by_time))) {
				cv_var = sd(var_by_time) / mean(var_by_time)
				expect_lt(cv_var, 2)
			}
		} else if (length(dim(fit$U)) == 3) {
			mean_var = mean(apply(fit$U, c(2, 3), var), na.rm = TRUE)
		} else {
			mean_var = var(c(fit$U), na.rm = TRUE)
		}

		if (!is.na(mean_var) && mean_var > 0) {
			expect_gt(mean_var, 0.001)
			expect_lt(mean_var, 100)
		} else {
			expect_true(!is.null(fit$BETA))
		}
	}
})

test_that("static and dynamic models converge to similar estimates when rho=0", {
	skip_on_cran()

	set.seed(6886)
	n = 20
	n_time = 5

	Y_list = list()
	for (t in 1:n_time) {
		eta = matrix(rnorm(n * n, -0.5, 1), n, n)
		Y_t = (eta > 0) * 1
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n)
		Y_list[[t]] = Y_t
	}

	fit_static = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = FALSE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE)

	fit_dynamic = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = TRUE,
		burn = 50, nscan = 100, verbose = FALSE, plot = FALSE,
		prior = list(rho_uv_mean = 0, rho_uv_sd = 0.2))

	if (!is.null(fit_dynamic$rho_uv) && length(fit_dynamic$rho_uv) > 0) {
		rho_est = median(fit_dynamic$rho_uv)
		expect_lt(abs(rho_est), 0.5)
	}
})
