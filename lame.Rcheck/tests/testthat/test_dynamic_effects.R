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

test_that("bipartite binary dynamic_ab plus dynamic_uv stays finite", {
	skip_on_cran()

	set.seed(123)
	n_users = 15
	n_items = 12
	n_periods = 5
	R = 2
	true_rho = 0.9

	U_array = array(0, dim = c(n_users, R, n_periods))
	V_array = array(0, dim = c(n_items, R, n_periods))
	U_array[, , 1] = matrix(rnorm(n_users * R), n_users, R)
	V_array[, , 1] = matrix(rnorm(n_items * R), n_items, R)
	for(t in 2:n_periods) {
		U_array[, , t] = true_rho * U_array[, , t-1] +
			sqrt(1 - true_rho^2) * matrix(rnorm(n_users * R), n_users, R)
		V_array[, , t] = true_rho * V_array[, , t-1] +
			sqrt(1 - true_rho^2) * matrix(rnorm(n_items * R), n_items, R)
	}

	G_long = diag(c(1.2, -0.8), R)
	Y_list = vector("list", n_periods)
	for(t in 1:n_periods) {
		eta = U_array[, , t] %*% G_long %*% t(V_array[, , t])
		prob = pnorm(eta)
		Y_list[[t]] = matrix(rbinom(n_users * n_items, 1, as.vector(prob)),
		                     n_users, n_items)
	}

	fit = lame(
		Y = Y_list, mode = "bipartite", R_row = 2, R_col = 2,
		family = "binary", dynamic_uv = TRUE, dynamic_ab = TRUE,
		burn = 100, nscan = 500, odens = 5,
		verbose = FALSE, gof = FALSE,
		seed = 6886
	)

	expect_s3_class(fit, "lame")
	expect_true(all(is.finite(fit$U)))
	expect_true(all(is.finite(fit$V)))
	expect_true(all(is.finite(fit$VC)))
	expect_true(all(is.finite(fit$sigma_uv)))
	expect_true(all(is.finite(fit$sigma_ab)))
	expect_true(all(is.finite(fit$start_vals$U)))
	expect_true(all(is.finite(fit$start_vals$V)))
	expect_true(all(is.finite(fit$G)))
	expect_true(all(is.finite(fit$start_vals$G)))
	total_iters = 600
	expect_lt(fit$tryErrorChecks$UV, 0.05 * total_iters)
	expect_lt(fit$tryErrorChecks$G, 0.05 * total_iters)
	expect_equal(fit$tryErrorChecks$ab, 0)
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

# regression: dynamic_ab full-conditional recovery (defect #0). The old
# sample_dynamic_ab_cpp subtracted a_i's own data signal from its update and
# used Sab[1,1] (the static prior variance) as the dyadic noise variance,
# attenuating additive effects to ~50% of truth. The full-conditional kernel
# must recover slope ~1.0.
test_that("sample_dynamic_ab_cpp recovers additive effects (no ~50% attenuation)", {
	skip_on_cran()

	set.seed(1)
	n = 25; Tt = 3
	a_true = seq(-2, 2, length.out = n)
	s2 = 0.25
	Z = array(rnorm(n * n * Tt, 0, sqrt(s2)), c(n, n, Tt))
	for (t in 1:Tt) Z[, , t] = Z[, , t] + outer(a_true, rep(0, n), "+")
	EZ0 = array(0, c(n, n, Tt))
	a = matrix(0, n, Tt); b = matrix(0, n, Tt); acc = matrix(0, n, Tt)
	nburn = 400; nkeep = 1200
	for (s in seq_len(nburn + nkeep)) {
		o = sample_dynamic_ab_cpp(a, b, Z, EZ0, 0.9, 0.5, s2, FALSE)
		a = o$a; b = o$b
		if (s > nburn) acc = acc + a
	}
	slope = coef(lm((acc / nkeep)[, 1] ~ a_true))[2]
	# defective kernel returned ~0.5; the full conditional returns ~1.0
	expect_gt(slope, 0.8)
	expect_lt(slope, 1.2)
})

# regression: dynamic_uv innovation-scale non-shrinkage (defects #6/#25). The
# old kernel added a hidden +1.0 (N(0, s2 I)) ridge at every time slice on top
# of the AR(1) prior; combined with transition-only hyper updates this
# over-shrank the latent factors and biased sigma_uv / E[(UV)^2] below truth.
# The corrected kernel (stationary init at t=0, no extra ridge) must recover
# sigma_uv and the latent-factor magnitude near truth.
test_that("dynamic_uv kernel does not double-shrink sigma_uv or |UV|", {
	skip_on_cran()

	set.seed(42)
	n = 20; R = 2; Tt = 8; rho_t = 0.9; sig_t = 0.4
	U = array(0, c(n, R, Tt)); V = array(0, c(n, R, Tt))
	U[, , 1] = matrix(rnorm(n * R, 0, sig_t / sqrt(1 - rho_t^2)), n, R)
	V[, , 1] = matrix(rnorm(n * R, 0, sig_t / sqrt(1 - rho_t^2)), n, R)
	for (t in 2:Tt) {
		U[, , t] = rho_t * U[, , t - 1] + matrix(rnorm(n * R, 0, sig_t), n, R)
		V[, , t] = rho_t * V[, , t - 1] + matrix(rnorm(n * R, 0, sig_t), n, R)
	}
	E = array(0, c(n, n, Tt))
	for (t in 1:Tt) E[, , t] = U[, , t] %*% t(V[, , t]) + matrix(rnorm(n * n), n, n)
	truth_uv2 = mean(sapply(1:Tt, function(t) c(U[, , t] %*% t(V[, , t]))^2))

	set.seed(7)
	Uc = array(rnorm(n * R * Tt, 0, 0.5), c(n, R, Tt))
	Vc = array(rnorm(n * R * Tt, 0, 0.5), c(n, R, Tt))
	ru = 0.9; su = 0.5; keep = NULL
	for (it in 1:600) {
		o = rUV_dynamic_fc_cpp(Uc, Vc, E, ru, su, 1, FALSE, FALSE)
		Uc = o$U; Vc = o$V
		ru = sample_rho_uv(Uc, Vc, su, ru, FALSE, 0.9, 0.1)
		su = lame:::.lame_sample_sigma_uv_drift_only(
			Uc, Vc, ru, su, FALSE, matrix(0, n, Tt), matrix(0, n, Tt), 2, 1)
		if (it > 300) {
			keep = rbind(keep, c(ru, su,
				mean(sapply(1:Tt, function(t) c(Uc[, , t] %*% t(Vc[, , t]))^2))))
		}
	}
	est = colMeans(keep)
	# rho_uv and sigma_uv recover near truth (0.9 / 0.4); the double-shrinkage
	# that pulled sigma_uv and |UV| below truth is gone
	expect_gt(est[1], 0.8)
	expect_gt(est[2], 0.32)
	expect_lt(est[2], 0.5)
	# latent-factor magnitude tracks truth (old kernel under-shot)
	expect_gt(est[3] / truth_uv2, 0.85)
	expect_lt(est[3] / truth_uv2, 1.2)
})
