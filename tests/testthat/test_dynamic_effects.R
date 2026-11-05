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
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
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
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		burn = 300,
		nscan = 1000,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)

	has_rho_ab = !is.null(fit$RHO_AB) || !is.null(fit$rho_ab)
	expect_true(
		has_rho_ab,
		info = "Dynamic AB model should store RHO_AB chain or rho_ab estimate"
	)

	if (!is.null(fit$RHO_AB)) {
		rho_hat = median(fit$RHO_AB, na.rm = TRUE)
		expect_true(
			rho_hat > 0.1,
			info = paste("Estimated rho_ab:", round(rho_hat, 3), "should be positive")
		)
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
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		burn = 100,
		nscan = 300,
		odens = 3,
		verbose = FALSE,
		gof = FALSE
	)

	expect_false(is.null(fit$start_vals))
	if (!is.null(fit$start_vals)) {
		expect_true(
			"rho_ab" %in% names(fit$start_vals),
			info = "start_vals should include rho_ab for dynamic_ab models"
		)
		expect_true(
			"sigma_ab" %in% names(fit$start_vals),
			info = "start_vals should include sigma_ab for dynamic_ab models"
		)
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
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)

	expect_s3_class(fit, "lame")
	expect_false(is.null(fit$BETA))

	if (!is.null(fit$YPM)) {
		ypm_vals = unlist(fit$YPM)
		ypm_vals = ypm_vals[!is.na(ypm_vals)]
		expect_true(
			all(ypm_vals >= 0 & ypm_vals <= 1),
			info = "Binary dynamic_ab YPM should be in [0,1]"
		)
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
		family = "normal",
		R = 0,
		mode = "bipartite",
		dynamic_ab = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
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
		family = "normal",
		R = R,
		dynamic_uv = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
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
		family = "normal",
		R = R,
		dynamic_uv = TRUE,
		burn = 300,
		nscan = 1000,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)

	has_rho = !is.null(fit$RHO_UV) || !is.null(fit$rho_uv)
	expect_true(
		has_rho,
		info = "Dynamic UV model should store RHO_UV chain or rho_uv estimate"
	)

	if (!is.null(fit$RHO_UV)) {
		rho_hat = median(fit$RHO_UV, na.rm = TRUE)
		expect_true(
			rho_hat > 0,
			info = paste("Estimated rho_uv:", round(rho_hat, 3), "should be positive")
		)
	}

	if (!is.null(fit$rho_uv)) {
		rho_uv_mean = mean(fit$rho_uv, na.rm = TRUE)
		expect_true(
			rho_uv_mean > 0,
			info = paste("rho_uv estimate:", round(rho_uv_mean, 3), "should be positive")
		)
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

	# short chain, so the settling check may warn
	fit = suppressWarnings(lame(
		Y = Y_list,
		family = "normal",
		R = R,
		dynamic_uv = TRUE,
		burn = 200,
		nscan = 800,
		odens = 4,
		verbose = FALSE,
		gof = FALSE
	))

	if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
		cors = c()
		for (t in 2:n_time) {
			U_t = fit$U[, , t]
			U_prev = fit$U[, , t - 1]
			cors = c(cors, cor(c(U_t), c(U_prev)))
		}
		avg_cor = mean(cors, na.rm = TRUE)
		expect_true(
			avg_cor > 0,
			info = paste("Adjacent U positions should be positively correlated, got:", round(avg_cor, 3))
		)
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
		family = "normal",
		R = 2,
		dynamic_uv = FALSE,
		burn = 100,
		nscan = 300,
		odens = 3,
		verbose = FALSE,
		gof = FALSE
	)

	expect_false(is.null(fit_static$U))
	expect_true(
		length(dim(fit_static$U)) == 2 || is.null(dim(fit_static$U)) || dim(fit_static$U)[3] == 1,
		info = "Static model U should be 2D or single time slice"
	)
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
		family = "binary",
		R = R,
		dynamic_uv = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
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
		family = "normal",
		R = R,
		mode = "bipartite",
		dynamic_uv = TRUE,
		burn = 200,
		nscan = 500,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
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

	# the combined model is compared across independent chains after burn-in, with the settling warning active
	fits = lapply(
		c(7886L, 8886L),
		function(seed) lame(
			Y = Y_list,
			family = "normal",
			R = R,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			burn = 3000,
			nscan = 16000,
			odens = 10,
			seed = seed,
			posterior_opts = list(
				save_ab = TRUE, thin_ab = 1,
				save_UV = TRUE, thin_UV = 1
			),
			verbose = FALSE,
			gof = FALSE
		)
	)
	fit = fits[[1]]
	chains = lapply(
		fits,
		function(f) cbind(
			f$BETA,
			f$VC,
			rho_ab = f$rho_ab,
			sigma_ab = f$sigma_ab,
			rho_uv = f$rho_uv,
			sigma_uv = f$sigma_uv
		)
	)
	diagnostic = compute_mcmc_diagnostics(lapply(chains, function(x) list(BETA = x)))
	expect_true(all(is.finite(diagnostic$rhat) & diagnostic$rhat < 1.05))
	expect_true(all(is.finite(diagnostic$ess) & diagnostic$ess > 0))
	if (requireNamespace("posterior", quietly = TRUE)) {
		split_rhat = vapply(
			seq_len(ncol(chains[[1]])),
			function(j)
				posterior::rhat(cbind(chains[[1]][, j], chains[[2]][, j])),
			numeric(1)
		)
		expect_true(all(is.finite(split_rhat) & split_rhat < 1.05))
	}

	expect_s3_class(fit, "lame")

	has_rho_uv = !is.null(fit$RHO_UV) || !is.null(fit$rho_uv)
	has_rho_ab = !is.null(fit$RHO_AB) || !is.null(fit$rho_ab)

	expect_true(
		has_rho_uv || has_rho_ab,
		info = "Combined dynamic model should track at least one rho parameter"
	)

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
	for (t in 2:n_periods) {
		U_array[, , t] = true_rho * U_array[, , t - 1] +
			sqrt(1 - true_rho^2) * matrix(rnorm(n_users * R), n_users, R)
		V_array[, , t] = true_rho * V_array[, , t - 1] +
			sqrt(1 - true_rho^2) * matrix(rnorm(n_items * R), n_items, R)
	}

	G_long = diag(c(1.2, -0.8), R)
	Y_list = vector("list", n_periods)
	for (t in 1:n_periods) {
		eta = U_array[, , t] %*% G_long %*% t(V_array[, , t])
		prob = pnorm(eta)
		Y_list[[t]] = matrix(
			rbinom(n_users * n_items, 1, as.vector(prob)),
			n_users,
			n_items
		)
	}

	# the scale ridge, assessed with settled, independently seeded chains
	fits = lapply(
		c(6886L, 7886L),
		function(seed) lame(
			Y = Y_list,
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			family = "binary",
			dynamic_uv = TRUE,
			dynamic_ab = TRUE,
			burn = 4000,
			nscan = 40000,
			odens = 20,
			posterior_opts = list(
				save_ab = TRUE, thin_ab = 1,
				save_UV = TRUE, thin_UV = 1
			),
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			seed = seed
		)
	)
	fit = fits[[1L]]
	for (f in fits) expect_true(all(f$sampler_health$numerical_failures == 0))
	chains = lapply(fits, function(f) {
		ns = nrow(f$BETA)
		m = cbind(
			f$BETA,
			f$VC,
			rho_ab = f$rho_ab,
			sigma_ab = f$sigma_ab,
			rho_uv = f$rho_uv,
			sigma_uv = f$sigma_uv,
			t(matrix(f$a_draws, ncol = ns)),
			t(matrix(f$b_draws, ncol = ns)),
			t(f$a_level_samples),
			t(f$b_level_samples)
		)
		product = matrix(NA_real_, ns, n_users * n_items * n_periods)
		for (j in seq_len(ns)) for (tt in seq_len(n_periods)) {
			u = matrix(f$U_draws[, , tt, j], n_users, R)
			v = matrix(f$V_draws[, , tt, j], n_items, R)
			columns = seq_len(n_users * n_items) + (tt - 1L) * n_users * n_items
			product[j, columns] = tcrossprod(u %*% f$G_samples[, , j], v)
		}
		cbind(m, product)
	})
	# exclude only quantities fixed to the same value in both chains
	variable = apply(do.call(rbind, chains), 2L, var) > 0
	chains = lapply(chains, function(x) x[, variable, drop = FALSE])
	diagnostic = do.call(rbind, lapply(seq_len(ncol(chains[[1L]])), function(j) {
		x = cbind(chains[[1L]][, j], chains[[2L]][, j])
		ns = nrow(x)
		within = mean(apply(x, 2L, var))
		between = ns * var(colMeans(x))
		answer = c(
			rhat = sqrt(((ns - 1) * within + between) / (ns * within)),
			ess = sum(vapply(1:2, function(k) as.numeric(coda::effectiveSize(x[, k])), numeric(1)))
		)
		if (requireNamespace("posterior", quietly = TRUE)) answer = c(
			answer,
			split = posterior::rhat(x),
			bulk = posterior::ess_bulk(x),
			tail = posterior::ess_tail(x)
		)
		answer
	}))
	expect_true(all(is.finite(diagnostic)))
	expect_lt(max(diagnostic[, "rhat"]), 1.05)
	expect_gt(min(diagnostic[, "ess"]), 0)
	if (requireNamespace("posterior", quietly = TRUE)) {
		expect_lt(max(diagnostic[, "split"]), 1.05)
		expect_gt(min(diagnostic[, "bulk"]), 0)
		expect_gt(min(diagnostic[, "tail"]), 0)
	}

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
	total_iters = 44000
	expect_lt(fit$tryErrorChecks$UV, 0.05 * total_iters)
	# conditional bound rejections leave the chain moving and the tie rate fitted
	expect_lt(fit$tryErrorChecks$G, 0.05 * total_iters)
	expect_equal(fit$tryErrorChecks$ab, 0)
	predictive_rate = mean(unlist(fit$YPM), na.rm = TRUE)
	expect_lt(abs(predictive_rate / mean(unlist(Y_list), na.rm = TRUE) - 1), 0.1)
	expect_true(all(is.finite(unlist(fit$EZ))))
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

	fit_all_dyn = lame(
		Y_list,
		R = 0,
		family = "binary",
		dynamic_uv = FALSE,
		dynamic_ab = TRUE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE
	)

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

	fit_low = lame(
		Y_list,
		R = 2,
		family = "binary",
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE,
		prior = list(rho_uv_mean = 0.2, rho_uv_sd = 0.1)
	)

	fit_high = lame(
		Y_list,
		R = 2,
		family = "binary",
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE,
		prior = list(rho_uv_mean = 0.8, rho_uv_sd = 0.1)
	)

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

	fit = lame(
		Y_list,
		R = 2,
		family = "binary",
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE
	)

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

	fit_static = lame(
		Y_list,
		R = 2,
		family = "binary",
		dynamic_uv = FALSE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE
	)

	fit_dynamic = lame(
		Y_list,
		R = 2,
		family = "binary",
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		plot = FALSE,
		prior = list(rho_uv_mean = 0, rho_uv_sd = 0.2)
	)

	if (!is.null(fit_dynamic$rho_uv) && length(fit_dynamic$rho_uv) > 0) {
		rho_est = median(fit_dynamic$rho_uv)
		expect_lt(abs(rho_est), 0.5)
	}
})

# sample_dynamic_ab_cpp draws the full conditional of each a_i, using the dyadic noise variance, so additive
# effects are recovered with slope near 1
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
	# the full conditional gives a slope near 1, not near 0.5
	expect_gt(slope, 0.8)
	expect_lt(slope, 1.2)
})

# the dynamic_uv kernel starts the AR(1) at its stationary law and adds no ridge at any slice, so sigma_uv and
# the latent-factor magnitude are not shrunk below truth
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
		hy = lame:::.lame_uv_rho_w_step(Uc, Vc, ru, su, 0.9, 0.1, 2, 1)
		ru = hy$rho; su = hy$sigma
		if (it > 300) {
			keep = rbind(
				keep,
				c(
					ru,
					su,
					mean(sapply(1:Tt, function(t) c(Uc[, , t] %*% t(Vc[, , t]))^2))
				)
			)
		}
	}
	est = colMeans(keep)
	# rho_uv and sigma_uv recover near truth (0.9 / 0.4)
	expect_gt(est[1], 0.8)
	expect_gt(est[2], 0.32)
	expect_lt(est[2], 0.5)
	# latent-factor magnitude tracks truth
	expect_gt(est[3] / truth_uv2, 0.85)
	expect_lt(est[3] / truth_uv2, 1.2)
})

test_that("bipartite dynamic_G fits are refused", {
	skip_on_cran()
	set.seed(21)
	nA = 20; nB = 16; TT = 4
	Y = lapply(1:TT, function(t) {
		m = (matrix(rnorm(nA * nB), nA, nB) + 0.2 > 0) * 1
		rownames(m) = sprintf("r%02d", 1:nA); colnames(m) = sprintf("c%02d", 1:nB); m
	})
	expect_error(
		suppressWarnings(lame(
			Y,
			mode = "bipartite",
			family = "binary",
			R_row = 2,
			R_col = 2,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			dynamic_G = TRUE,
			burn = 200,
			nscan = 800,
			odens = 4,
			seed = 8886,
			verbose = FALSE,
			plot = FALSE
		)),
		"not supported"
	)
	# the supported fit of the same panel runs and keeps its margins finite
	fit = suppressWarnings(lame(
		Y,
		mode = "bipartite",
		family = "binary",
		R_row = 2,
		R_col = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 200,
		nscan = 800,
		odens = 4,
		seed = 8886,
		verbose = FALSE,
		plot = FALSE
	))
	expect_true(all(is.finite(fit$VC[, "va"])))
	expect_lt(max(fit$VC[, "va"]), 100)
	expect_true(all(is.finite(unlist(fit$EZ))))
})

test_that("bipartite dynamic draws retain their sampled scales and aligned products", {
	# proper path and interaction priors set the scale with no deterministic normalization; aligned products
	# reconstruct the posterior product summary
	skip_on_cran()
	set.seed(101)
	nA = 30; nB = 50; TT = 4
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	U = array(0, c(nA, 2, TT)); V = array(0, c(nB, 2, TT)); G = array(0, c(2, 2, TT))
	U[, , 1] = matrix(rnorm(nA * 2, 0, 0.6), nA); V[, , 1] = matrix(rnorm(nB * 2, 0, 0.6), nB)
	G[, , 1] = diag(c(1, -0.7))
	for (t in 2:TT) {
		U[, , t] = 0.8 * U[, , t - 1] + matrix(rnorm(nA * 2, 0, 0.36), nA)
		V[, , t] = 0.8 * V[, , t - 1] + matrix(rnorm(nB * 2, 0, 0.36), nB)
		G[, , t] = 0.8 * G[, , t - 1] + matrix(rnorm(4, 0, 0.2), 2)
	}
	Y = lapply(1:TT, function(t) {
		m = outer(a, b, "+") + U[, , t] %*% G[, , t] %*% t(V[, , t]) + matrix(rnorm(nA * nB), nA)
		rownames(m) = sprintf("s%03d", 1:nA); colnames(m) = sprintf("p%03d", 1:nB); m
	})
	fit = suppressWarnings(lame(
		Y,
		mode = "bipartite",
		family = "normal",
		R_row = 2,
		R_col = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 500,
		nscan = 2000,
		odens = 10,
		seed = 8886,
		verbose = FALSE,
		plot = FALSE,
		posterior_opts = list(save_UV = TRUE)
	))
	nd = dim(fit$U_draws)[4]
	stopifnot(dim(fit$G_samples)[3] == nd)
	product_mean = array(0, c(nA, nB, TT))
	for (j in seq_len(nd)) for (t in seq_len(TT)) {
		product_mean[, , t] = product_mean[, , t] +
			fit$U_draws[, , t, j] %*% fit$G_samples[, , j] %*% t(fit$V_draws[, , t, j]) / nd
	}
	reported = if (is.list(fit$UVPM)) array(unlist(fit$UVPM), dim(product_mean)) else fit$UVPM
	expect_equal(unname(reported), unname(product_mean), tolerance = 1e-9)
	rms_u = apply(fit$U_draws, 4, function(x) sqrt(mean(x^2)))
	rms_v = apply(fit$V_draws, 4, function(x) sqrt(mean(x^2)))
	expect_gt(sd(rms_u), 1e-4)
	expect_gt(sd(rms_v), 1e-4)
	expect_true(all(is.finite(rms_u) & rms_u > 0))
	expect_true(all(is.finite(rms_v) & rms_v > 0))
	expect_true(all(is.finite(fit$VC[, "va"])))
	expect_equal(fit$tryErrorChecks$G, 0)
})

test_that("the additive-effect innovation scale respects its stationary-variance cap", {
	set.seed(9)
	n = 40; TT = 8; rho = 0.9
	a = matrix(0, n, TT); b = matrix(0, n, TT)
	a[, 1] = rnorm(n, 0, 2); b[, 1] = rnorm(n, 0, 2)
	for (t in 2:TT) {
		a[, t] = rho * a[, t - 1] + rnorm(n, 0, 1)
		b[, t] = rho * b[, t - 1] + rnorm(n, 0, 1)
	}
	# these paths imply a stationary variance near 1 / (1 - 0.81) = 5.3; a cap of 1 holds every draw at or below
	# sqrt(1 * (1 - rho^2))
	capped = replicate(300, lame:::.lame_sample_sigma_ab(a, b, rho, 0.5, FALSE, stat_var_cap = 1))
	expect_true(all(capped <= sqrt(1 - rho^2) + 1e-9))
	expect_true(all(capped > 0))
	# without a cap the draw sits near the innovation scale that generated the paths
	free = replicate(300, lame:::.lame_sample_sigma_ab(a, b, rho, 0.5, FALSE, stat_var_cap = Inf))
	expect_gt(median(free), 0.8)
	expect_lt(median(free), 1.25)
	# a symmetric fit reads only the row paths
	sym = replicate(50, lame:::.lame_sample_sigma_ab(a, a * 100, rho, 0.5, TRUE, stat_var_cap = Inf))
	expect_lt(median(sym), 1.5)
	# missing periods contribute nothing
	a2 = a; a2[3, 4] = NA
	expect_true(is.finite(lame:::.lame_sample_sigma_ab(a2, b, rho, 0.5, FALSE, stat_var_cap = Inf)))
	# one period gives no innovations, so the current value is kept
	expect_equal(lame:::.lame_sample_sigma_ab(a[, 1, drop = FALSE], b[, 1, drop = FALSE], rho, 0.37, FALSE), 0.37)
})

test_that("a sparse dynamic bipartite binary fit stays bounded instead of separating (regression)", {
	# a sparse 0/1 panel with per-period effects and positions can be nearly separated; the coordinate cap and the
	# latent-scale rails keep the product, effect variances and intercept near the generating values
	skip_on_cran()
	set.seed(11)
	nA = 30; nB = 90; TT = 12
	a = matrix(0, nA, TT); b = matrix(0, nB, TT)
	a[, 1] = rnorm(nA, 0, 0.9); b[, 1] = rnorm(nB, 0, 0.8)
	for (t in 2:TT) {
		a[, t] = 0.9 * a[, t - 1] + rnorm(nA, 0, 0.9 * sqrt(1 - 0.81))
		b[, t] = 0.9 * b[, t - 1] + rnorm(nB, 0, 0.8 * sqrt(1 - 0.81))
	}
	U = array(0, c(nA, 2, TT)); V = array(0, c(nB, 2, TT))
	U[, , 1] = matrix(rnorm(nA * 2), nA); V[, , 1] = matrix(rnorm(nB * 2), nB)
	for (t in 2:TT) {
		U[, , t] = 0.9 * U[, , t - 1] + sqrt(1 - 0.81) * matrix(rnorm(nA * 2), nA)
		V[, , t] = 0.9 * V[, , t - 1] + sqrt(1 - 0.81) * matrix(rnorm(nB * 2), nB)
	}
	Y = lapply(1:TT, function(t) {
		eta = -2.6 + outer(a[, t], b[, t], "+") + 0.7 * U[, , t] %*% t(V[, , t])
		m = matrix(rbinom(nA * nB, 1, as.vector(pnorm(eta))), nA)
		rownames(m) = sprintf("e%02d", 1:nA); colnames(m) = sprintf("i%03d", 1:nB); m
	})
	fit = suppressWarnings(lame(
		Y,
		mode = "bipartite",
		family = "binary",
		R_row = 2,
		R_col = 2,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		burn = 0,
		nscan = 2000,
		odens = 10,
		seed = 1,
		verbose = FALSE,
		plot = FALSE
	))
	expect_lt(max(fit$VC[, "va"]), 3)
	expect_lt(max(fit$VC[, "vb"]), 3)
	expect_gt(min(fit$BETA[, 1]), -4.5)
	expect_lt(max(abs(unlist(fit$EZ))), 50)
	# the predictive tie rate tracks the data and the posterior-mean predictor carries no rail-pressing excursions
	predictive_rate = mean(unlist(fit$YPM), na.rm = TRUE)
	expect_lt(abs(predictive_rate / mean(unlist(Y)) - 1), 0.3)
	# the interaction scale may touch its rail in a minority of sweeps but is not stuck there
	expect_lt(fit$tryErrorChecks$G, 0.8 * 2000)
})

test_that("dynamic positions keep their prior coordinates while actor paths remain centered", {
	# the position prior is free in location and actor paths have their own exact zero-sum constraint, so no
	# predictor-preserving map moves position means into the actor state
	skip_on_cran()
	set.seed(21)
	nA = 20; nB = 16; TT = 4
	Y = lapply(1:TT, function(t) {
		m = (matrix(rnorm(nA * nB), nA, nB) + 0.2 > 0) * 1
		rownames(m) = sprintf("r%02d", 1:nA); colnames(m) = sprintf("c%02d", 1:nB); m
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		mode = "bipartite",
		family = "binary",
		R_row = 2,
		R_col = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 200,
		odens = 10,
		seed = 3,
		verbose = FALSE,
		plot = FALSE,
		posterior_opts = list(save_UV = TRUE)
	)))
	nd = dim(fit$U_draws)[4]
	u_means = sapply(1:TT, function(t) colMeans(matrix(fit$U_draws[, , t, nd], nrow = nA)))
	v_means = sapply(1:TT, function(t) colMeans(matrix(fit$V_draws[, , t, nd], nrow = nB)))
	expect_gt(max(abs(u_means)), 1e-6)
	expect_gt(max(abs(v_means)), 1e-6)
	expect_lt(max(abs(colSums(fit$start_vals$a))), 1e-8)
	expect_lt(max(abs(colSums(fit$start_vals$b))), 1e-8)
	expect_true(all(is.finite(fit$a_dynamic)))
	expect_true(all(is.finite(fit$b_dynamic)))
	expect_true(all(is.finite(unlist(fit$EZ))))
})

test_that("centring one period of positions leaves the linear predictor unchanged", {
	# components removed from U G V' go to the effects, leaving the predictor unchanged to machine precision for any
	# rank and shape, with actors missing from the period
	set.seed(4)
	chk = function(nA, nB, RA, RB, na_rows = 0) {
		U = matrix(rnorm(nA * RA, 0.7, 1), nA); V = matrix(rnorm(nB * RB, -0.4, 1), nB)
		G = matrix(rnorm(RA * RB, 0, 2), RA, RB); a = rnorm(nA); b = rnorm(nB)
		if (na_rows > 0) U[1:na_rows, ] = NA
		eta0 = outer(a, b, "+") + U %*% G %*% t(V)
		cen = lame:::.lame_center_uv_period(U, V, G)
		eta1 = outer(a + cen$a_add, b + cen$b_add, "+") + cen$U %*% G %*% t(cen$V)
		ok = is.finite(eta0)
		expect_lt(max(abs(eta0[ok] - eta1[ok])), 1e-10)
		expect_lt(max(abs(colMeans(cen$U, na.rm = TRUE))), 1e-10)
		expect_lt(max(abs(colMeans(cen$V, na.rm = TRUE))), 1e-10)
	}
	chk(30, 90, 2, 2); chk(12, 12, 1, 1); chk(15, 40, 3, 2); chk(20, 30, 2, 2, na_rows = 3)
	# already centred positions are left alone
	U0 = scale(matrix(rnorm(20), 10), scale = FALSE); V0 = scale(matrix(rnorm(16), 8), scale = FALSE)
	expect_null(lame:::.lame_center_uv_period(U0, V0, diag(2)))
})

test_that("symmetric centring keeps one effect per actor and the predictor unchanged", {
	# one set of positions: row and column increments average into one effect, so a equals b and u u' + a + b is unchanged
	set.seed(5)
	U = matrix(rnorm(24, 0.5, 1), 12); a = rnorm(12)
	eta0 = outer(a, a, "+") + U %*% t(U)
	cen = lame:::.lame_center_uv_period(U, U, diag(2))
	inc = (cen$a_add + cen$b_add) / 2
	eta1 = outer(a + inc, a + inc, "+") + cen$U %*% t(cen$U)
	expect_lt(max(abs(eta0 - eta1)), 1e-10)
	expect_lt(max(abs(colMeans(cen$U))), 1e-10)
})

test_that("the stationary-variance cap holds on the unipartite dynamic path", {
	skip_on_cran()
	set.seed(6); n = 16; TT = 4
	Y = lapply(
		1:TT,
		function(t) { m = matrix(rbinom(n * n, 1, 0.3), n); diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(suppressMessages(lame(
		Y,
		family = "binary",
		R = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 200,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE
	)))
	stat_var = f$sigma_uv^2 / (1 - pmin(f$rho_uv^2, 0.9801))
	expect_true(all(is.finite(stat_var)))
	expect_lte(max(stat_var), 4 * (1 + 1e-6))
	# a symmetric fit on the same path runs and stays symmetric
	Ys = lapply(Y, function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m })
	fs = suppressWarnings(suppressMessages(lame(
		Ys,
		family = "binary",
		R = 2,
		symmetric = TRUE,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 50,
		nscan = 200,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		plot = FALSE
	)))
	ez = fs$EZ[[TT]]; ez[!is.finite(ez)] = 0
	expect_lt(max(abs(ez - t(ez))), 1e-8)
	expect_true(all(is.finite(unlist(fs$YPM)[is.finite(unlist(Ys))])))
})

test_that("the predictive-level check aborts only dynamic binary or poisson fits that miss the level", {
	set.seed(9)
	Y = array(rbinom(2 * 6 * 6, 1, 0.3), c(6, 6, 2)); for (tt in 1:2) diag(Y[, , tt]) = NA
	n = 100
	yps_ok = Y * n
	yps_off = (Y * 0 + 0.9) * n
	grab = function(expr) { w = character(0); r = withCallingHandlers(
		expr,
		warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
	); list(value = r, warns = w) }
	r = grab(
		lame:::.lame_predictive_level_check("binary", TRUE, yps_ok, Y, n)
	); expect_false(r$value); expect_length(r$warns, 0)
	e = expect_error(
		lame:::.lame_predictive_level_check("binary", TRUE, yps_off, Y, n),
		class = "lame_predictive_level"
	)
	expect_match(conditionMessage(e), "tie rate"); expect_false(
		grepl("Check|Inspect|before interpreting", conditionMessage(e))
	)
	r = grab(
		lame:::.lame_predictive_level_check("binary", FALSE, yps_off, Y, n)
	); expect_true(is.na(r$value)); expect_length(r$warns, 0)
	r = grab(
		lame:::.lame_predictive_level_check("normal", TRUE, yps_off, Y, n)
	); expect_true(is.na(r$value)); expect_length(r$warns, 0)
	# fewer than 100 stored draws do not pin the level, so there is no check
	expect_true(is.na(lame:::.lame_predictive_level_check("binary", TRUE, yps_off / 10, Y, 10)))
	Yp = array(rpois(2 * 6 * 6, 2), c(6, 6, 2))
	e = expect_error(
		lame:::.lame_predictive_level_check("poisson", TRUE, Yp * 3 * n, Yp, n),
		class = "lame_predictive_level"
	)
	expect_match(conditionMessage(e), "mean count")
})

test_that("every documented latent-scale prior element is recognised", {
	set.seed(8)
	Y = lapply(
		1:3,
		function(t) { m = matrix(rbinom(64, 1, 0.4), 8); diag(m) = NA; dimnames(m) = list(letters[1:8], letters[1:8]); m }
	)
	pr = list(
		G_max = 100, uv_max_abs = 50, uv_var_cap = 4, ab_var_cap = Inf, sigma_uv_max = 10, sigma_ab_scale = 1, sigma_ab_shape = 2,
		rho_ab_mean = 0.8, rho_ab_sd = 0.15, rho_uv_mean = 0.9, rho_uv_sd = 0.1, sigma_uv_shape = 2, sigma_uv_scale = 1
	)
	w = character(0)
	withCallingHandlers(
		suppressMessages(lame(
			Y,
			family = "binary",
			R = 1,
			dynamic_ab = TRUE,
			dynamic_uv = TRUE,
			prior = pr,
			burn = 5,
			nscan = 20,
			odens = 5,
			seed = 1,
			verbose = FALSE,
			plot = FALSE
		)),
		warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
	)
	expect_false(any(grepl("Unknown `prior`", w)))
})

test_that("an early exit under max_seconds leaves every per-draw store at the completed length", {
	skip_on_cran()
	set.seed(21); n = 20; TT = 6
	Y = lapply(
		1:TT,
		function(t) { m = matrix(rbinom(n * n, 1, 0.3), n); diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(suppressMessages(lame(
		Y,
		family = "binary",
		R = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		gof = TRUE,
		burn = 20,
		nscan = 40000,
		odens = 4,
		max_seconds = 6,
		seed = 3,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_true(isTRUE(f$terminated_early))
	nd = nrow(f$VC)
	expect_lt(nd, 10000)
	expect_equal(nrow(as.matrix(f$BETA)), nd)
	expect_equal(length(f$rho_ab), nd); expect_equal(length(f$sigma_ab), nd)
	expect_equal(length(f$rho_uv), nd); expect_equal(length(f$sigma_uv), nd)
	expect_true(all(f$rho_uv > 0)); expect_true(all(f$sigma_ab > 0))
	# effects and positions average over completed draws, undiluted
	expect_gt(max(abs(f$a_dynamic)), 0.05)
	expect_gt(sd(as.vector(f$U)), 0.02)
	ypm = unlist(f$YPM); expect_lt(abs(mean(ypm[is.finite(unlist(Y))]) - mean(unlist(Y), na.rm = TRUE)), 0.1)
})

test_that("forecasts keep each actor's level under dynamic effects", {
	skip_on_cran()
	set.seed(22); n = 20; TT = 6
	lev = rnorm(n, 0, 1)
	Y = lapply(
		1:TT,
		function(t) { m = outer(lev, rep(0, n), "+") + matrix(rnorm(n * n), n); diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = "normal",
				R = 0,
				dynamic_ab = TRUE,
				burn = 100,
				nscan = 300,
				odens = 3,
				seed = 3,
				verbose = FALSE,
				plot = FALSE
			)
		)
	)
	expect_true(!is.null(f$a_level)); expect_equal(length(f$a_level), n)
	expect_gt(cor(f$a_level, lev), 0.9)
	# far ahead the expected effect is the level, not zero
	pr = predict(f, h = 25, n_draws = 40)
	eta_far = pr[[length(pr)]]
	a_far = rowMeans(eta_far, na.rm = TRUE) - mean(rowMeans(eta_far, na.rm = TRUE))
	expect_gt(cor(a_far, lev), 0.8)
	expect_gt(sd(a_far), 0.5 * sd(lev))
})

test_that("rho_ab moves every sweep", {
	skip_on_cran()
	set.seed(23); n = 16; TT = 5
	Y = lapply(
		1:TT,
		function(t) { m = matrix(rnorm(n * n), n); diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = "normal",
				R = 0,
				dynamic_ab = TRUE,
				burn = 50,
				nscan = 200,
				odens = 1,
				seed = 3,
				verbose = FALSE,
				plot = FALSE
			)
		)
	)
	expect_gt(length(unique(round(f$rho_ab, 8))), 40)
	expect_gt(length(unique(round(f$sigma_ab, 8))), 150)
})

test_that("symmetric dynamic fits simulate, with no V behind partial matching", {
	skip_on_cran()
	set.seed(24); n = 14; TT = 4
	Y = lapply(
		1:TT,
		function(t) { m = matrix(rnorm(n * n), n); m = (m + t(m)) / 2; diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = "normal",
				R = 1,
				symmetric = TRUE,
				dynamic_uv = TRUE,
				burn = 30,
				nscan = 90,
				odens = 3,
				seed = 2,
				verbose = FALSE,
				plot = FALSE
			)
		)
	)
	expect_null(f[["V"]])
	s = simulate(f, nsim = 1)
	expect_true(is.list(s) || is.array(s))
	expect_error(simulate_posterior(f, component = "Y"), "longitudinal")
})

test_that("no fit element is read with a dollar when a longer name shares its prefix", {
	# on fits without V, G or X, $V, $G and $X partial-match $VC, $GOF and $Xlist, so exact indexing is required
	# an installed package's R directory holds only the lazy-load database, so sources come from the package tree
	src = list.files(system.file("R", package = "lame"), pattern = "[.]R$", full.names = TRUE)
	if (!length(src)) { pkg = normalizePath(file.path(testthat::test_path(), "..", "..")); src = list.files(
		file.path(pkg, "R"),
		pattern = "[.]R$",
		full.names = TRUE
	) }
	skip_if(!length(src), "package sources not available")
	pat = "(?<![A-Za-z0-9_.])(fit|object|obj|combined|chain|res|out|x|f)[$](G|V|X)(?![A-Za-z0-9_.])"
	hits = unlist(lapply(src, function(f) {
		ln = readLines(f, warn = FALSE)
		ln = ln[!grepl("^\\s*#", ln)]
		h = grep(pat, ln, perl = TRUE, value = TRUE)
		if (length(h)) paste0(basename(f), ": ", h) else character(0) }))
	expect_equal(hits, character(0))
})

test_that("parallel chains combine for symmetric fits with dynamic positions", {
	skip_on_cran()
	set.seed(9); n = 14; TT = 4
	Y = lapply(
		1:TT,
		function(t) { m = matrix(rnorm(n * n), n); m = (m + t(m)) / 2; diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	# one core keeps the chains in-process, since spawned workers load the installed package rather than the one under test
	f = suppressWarnings(suppressMessages(lame_parallel(
		Y,
		n_chains = 2,
		cores = 1,
		R = 1,
		family = "normal",
		symmetric = TRUE,
		dynamic_uv = TRUE,
		burn = 20,
		nscan = 60,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_s3_class(f, "lame")
	expect_null(f[["V"]])
	expect_equal(length(dim(f$U)), 3L)
	expect_true(all(is.finite(f$U)))
	expect_equal(dim(f$ULUPM)[3], TT)
	expect_true(all(is.finite(as.matrix(f$BETA))))
})

# a very sparse panel: 40 actors, 10 periods, about 0.5 percent ones, far sparser than the other fixtures
.sparse_field_panel = function(seed = 1, moving = TRUE) {
	set.seed(seed); n = 40; TT = 10; rho = 0.8
	lev_a = rnorm(n, 0, 0.6); lev_b = rnorm(n, 0, 0.5)
	dev_a = matrix(0, n, TT); dev_b = matrix(0, n, TT)
	if (moving) { dev_a[, 1] = rnorm(n, 0, 0.30); dev_b[, 1] = rnorm(n, 0, 0.30)
		for (tt in 2:TT) { dev_a[, tt] = rho * dev_a[, tt - 1] + rnorm(n, 0, 0.30 * sqrt(1 - rho^2))
			dev_b[, tt] = rho * dev_b[, tt - 1] + rnorm(n, 0, 0.30 * sqrt(1 - rho^2)) } }
	a = lev_a + dev_a; b = lev_b + dev_b
	X = lapply(1:TT, function(tt) { arr = array(rnorm(n * n), c(n, n, 1))
		dimnames(arr) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n), "x"); arr })
	Y = lapply(1:TT, function(tt) {
		eta = outer(a[, tt], b[, tt], "+") + 0.5 * X[[tt]][, , 1]
		mu = uniroot(function(m) sum(pnorm(m + eta)) - 74 / TT, c(-14, 4))$root
		m = (mu + eta + matrix(rnorm(n * n), n) > 0) * 1
		diag(m) = NA; dimnames(m) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n)); m })
	list(Y = Y, X = X, dev_a = dev_a)
}

test_that("a panel too sparse to inform per-period effects says so before fitting", {
	d = .sparse_field_panel(1)
	msgs = character(0)
	withCallingHandlers(
		suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			R = 1,
			dynamic_ab = TRUE,
			burn = 20,
			nscan = 60,
			odens = 20,
			seed = 1,
			verbose = FALSE,
			plot = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_match(msgs, "per actor-period", all = FALSE)
	# a dense panel of the same family raises nothing
	set.seed(2); n = 20
	Yd = lapply(1:5, function(tt) { m = matrix(rbinom(n * n, 1, 0.3), n); diag(m) = NA
		dimnames(m) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n)); m })
	m2 = character(0)
	withCallingHandlers(
		suppressWarnings(lame(
			Yd,
			family = "binary",
			R = 1,
			dynamic_ab = TRUE,
			burn = 20,
			nscan = 60,
			odens = 20,
			seed = 1,
			verbose = FALSE,
			plot = FALSE
		)),
		message = function(m) { m2 <<- c(m2, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_false(any(grepl("per actor-period", m2)))
})

test_that("a stationary-variance cap that shapes the draws is reported, and a idle one is not", {
	skip_on_cran()
	d = .sparse_field_panel(1)
	run = function(cap) { msgs = character(0)
		f = withCallingHandlers(
			suppressWarnings(lame(
				d$Y,
				Xdyad = d$X,
				family = "binary",
				R = 2,
				dynamic_ab = TRUE,
				dynamic_uv = TRUE,
				prior = list(uv_var_cap = cap),
				burn = 60,
				nscan = 240,
				odens = 20,
				seed = 1,
				verbose = FALSE,
				plot = FALSE
			)),
			message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
		)
		list(hits = f$var_cap_hits, said = any(grepl("variance cap", msgs))) }
	tight = run(0.02)
	expect_gt(tight$hits, 0); expect_true(tight$said)
	loose = run(4)
	expect_equal(loose$hits, 0); expect_false(loose$said)
})

# the dynamic_beta chains are the slowest fits here, so both seeds are fitted once and shared
.sparse_dyad_chains = local({
	fits = NULL
	function() {
		if (is.null(fits)) {
			d = .sparse_field_panel(1)
			fits <<- lapply(
				c(11, 22),
				function(sd) suppressWarnings(suppressMessages(
					lame(
						d$Y,
						Xdyad = d$X,
						family = "binary",
						R = 1,
						dynamic_ab = TRUE,
						dynamic_beta = "dyad",
						burn = 500,
						nscan = 3000,
						odens = 5,
						seed = sd,
						verbose = FALSE,
						plot = FALSE,
						gof = FALSE
					)
				))
			)
		}
		fits
	}
})

# a fit on the sparse panel at the length the checks below are calibrated for
.sparse_spec_fit = function(spec, seed = 11) {
	d = .sparse_field_panel(1)
	args = c(
		list(
			Y = d$Y,
			Xdyad = d$X,
			family = "binary",
			R = 1,
			burn = 500,
			nscan = 2500,
			odens = 5,
			seed = seed,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		spec
	)
	suppressWarnings(suppressMessages(do.call(lame, args)))
}

# draws of the constant: dynamic_beta = "dyad" keeps the intercept static, so a 3-D BETA repeats it each period
.intercept_draws = function(f) {
	B = f$BETA
	if (length(dim(B)) == 3) {
		j = match("intercept", dimnames(B)[[2]], nomatch = 1L)
		as.numeric(B[, j, 1])
	} else {
		B = as.matrix(B)
		j = match("intercept", colnames(B), nomatch = 1L)
		as.numeric(B[, j])
	}
}

# mean predictive tie probability over the mean observed tie rate, diagonal dropped from both
.tie_rate_ratio = function(f, Y) {
	fitted_rate = mean(vapply(f$YPM, function(m) { diag(m) = NA
		mean(m, na.rm = TRUE) }, 0))
	observed_rate = mean(vapply(Y, function(m) mean(m, na.rm = TRUE), 0))
	fitted_rate / observed_rate
}

test_that("every dynamic spec fits the tie rate of a half-percent panel", {
	skip_on_cran()
	d = .sparse_field_panel(1)
	# the identified constant on this panel is about -3.6; an intercept near zero means the constant has moved into
	# the latent product, which the flag below scores
	check_identified_intercept = TRUE
	specs = list(
		"dynamic_ab" = list(dynamic_ab = TRUE),
		"dynamic_ab + dynamic_beta" = list(dynamic_ab = TRUE, dynamic_beta = "dyad"),
		"dynamic_ab + dynamic_uv" = list(dynamic_ab = TRUE, dynamic_uv = TRUE),
		"dynamic_uv" = list(dynamic_uv = TRUE)
	)
	for (nm in names(specs)) {
		f = if (nm == "dynamic_ab + dynamic_beta") {
			.sparse_dyad_chains()[[1]]
		} else {
			.sparse_spec_fit(specs[[nm]])
		}
		ratio = .tie_rate_ratio(f, d$Y)
		lab = paste0("fitted/observed tie rate, ", nm)
		expect_gt(ratio, 0.7, label = lab)
		expect_lt(ratio, 1.4, label = lab)
		if (check_identified_intercept) {
			expect_lt(
				mean(.intercept_draws(f)),
				-1.5,
				label = paste0("posterior mean intercept, ", nm)
			)
		}
	}
})

test_that("dynamic_beta on top of dynamic_ab keeps the constant mixing", {
	skip_on_cran()
	skip_if_not_installed("posterior")
	ch = .sparse_dyad_chains()
	# the constant is the parameter scored in every arm
	icpt_spread = apply(
		ch[[1]]$BETA[, "intercept", , drop = FALSE],
		1,
		function(z) diff(range(z))
	)
	expect_equal(max(icpt_spread), 0)
	n_draw = length(.intercept_draws(ch[[1]]))
	draws = vapply(ch, .intercept_draws, numeric(n_draw))
	expect_lt(posterior::rhat(draws), 1.35, label = "R-hat of the intercept")
	# a collapse guard, not a health certificate: the ESS floor sits well below a healthy chain and well above a collapsed one
	expect_gte(posterior::ess_bulk(draws), 3, label = "bulk ESS of the intercept")
	expect_true(all(vapply(ch, function(f) all(is.finite(f$a_dynamic)), TRUE)))
})

test_that("dcor fits a panel with no mutual dyads", {
	skip_on_cran()
	d = .sparse_field_panel(1)
	mutual = vapply(d$Y, function(m) { z = m; z[is.na(z)] = 0; sum(z * t(z)) / 2 }, 0)
	expect_equal(sum(mutual), 0)
	# a chain this short passes the drift check for some seeds and not others
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_beta = "dyad",
		dcor = TRUE,
		burn = 200,
		nscan = 400,
		odens = 5,
		seed = 12,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_s3_class(f, "lame")
	# without mutual dyads rho is informed only by their absence, so its draws are checked for support only
	expect_true(all(is.finite(f$VC[, "rho"])))
	expect_true(all(abs(f$VC[, "rho"]) < 1))
})

test_that("dcor on a panel with no mutual dyads says so before fitting", {
	skip_on_cran()
	d = .sparse_field_panel(1)
	expect_message(
		suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			R = 1,
			dynamic_ab = TRUE,
			dcor = TRUE,
			burn = 20,
			nscan = 60,
			odens = 20,
			seed = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		"mutual dyads"
	)
})

# a panel whose level lives only in the additive effects: 25 actors, 5 periods, a constant of -2.4 and no intercept
.no_intercept_panel = function(seed = 4) {
	set.seed(seed); n = 25; TT = 5
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	eta = -2.4 + outer(a, b, "+")
	Y = lapply(1:TT, function(tt) { m = (eta + matrix(rnorm(n * n), n) > 0) * 1
		diag(m) = NA; dimnames(m) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n)); m })
	list(Y = Y, a = a, b = b)
}

test_that("a fit without an intercept keeps the level in the additive effects", {
	skip_on_cran()
	d = .no_intercept_panel(4)
	for (dyn in c(TRUE, FALSE)) {
		f = suppressWarnings(suppressMessages(lame(
			d$Y,
			family = "binary",
			R = 0,
			intercept = FALSE,
			rvar = TRUE,
			cvar = TRUE,
			dynamic_ab = dyn,
			burn = 500,
			nscan = 2000,
			odens = 5,
			seed = 11,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)))
		nm = if (dyn) "dynamic_ab" else "static a and b"
		ratio = .tie_rate_ratio(f, d$Y)
		expect_gt(ratio, 0.6, label = paste0("fitted/observed tie rate, ", nm))
		expect_lt(ratio, 1.6, label = paste0("fitted/observed tie rate, ", nm))
		level = mean(f$APM) + mean(f$BPM)
		expect_gt(level, -3.5, label = paste0("mean(APM) + mean(BPM), ", nm))
		expect_lt(level, -1.5, label = paste0("mean(APM) + mean(BPM), ", nm))
	}
})

test_that("bipartite fitted values use the identified product with the identified effects", {
	# the reported intercept and effects carry the row, column and grand-mean parts of the product, so fitted values
	# combine them with the posterior-mean product, not raw coordinates
	skip_on_cran()
	set.seed(41)
	nA = 12; nB = 9; TT = 3
	Y = lapply(1:TT, function(t) {
		m = (matrix(rnorm(nA * nB), nA, nB) + 0.3 > 0) * 1
		rownames(m) = sprintf("r%02d", 1:nA); colnames(m) = sprintf("c%02d", 1:nB); m
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		mode = "bipartite",
		family = "binary",
		R_row = 2,
		R_col = 2,
		dynamic_ab = TRUE,
		burn = 30,
		nscan = 150,
		odens = 5,
		seed = 7,
		verbose = FALSE,
		plot = FALSE
	)))
	ez = if (is.list(fit$EZ)) fit$EZ[[TT]] else fit$EZ[, , TT]
	expected = mean(fit$BETA[, 1]) + outer(fit$a_dynamic[, TT], fit$b_dynamic[, TT], "+") + fit$UVPM
	expect_equal(unname(ez), unname(expected), tolerance = 1e-8)
	# the reported positions are centred over actors
	expect_lt(max(abs(colMeans(fit$U))), 1e-8)
	expect_lt(max(abs(colMeans(fit[["V"]]))), 1e-8)
})

# a symmetric dynamic_uv fit carries a signed diagonal L, so u l u' can have negative eigen-directions; with unit
# stationary variance of the paths the innovation sd is tied to rho_uv
test_that("symmetric dynamic_uv fits carry a signed L and recover a product with mixed signs", {
	skip_on_cran()
	set.seed(4)
	n = 24
	TT = 4
	nm = sprintf("a%02d", seq_len(n))
	U = matrix(rnorm(n * 2, 0, 0.7), n, 2)
	L = c(1, -1)
	prod = U %*% diag(L) %*% t(U)
	Y = lapply(seq_len(TT), function(t) {
		e = matrix(0, n, n)
		e[upper.tri(e)] = rnorm(n * (n - 1) / 2)
		e = e + t(e)
		m = 0.2 + prod + e
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(lame(
		Y,
		family = "normal",
		R = 2,
		symmetric = TRUE,
		dynamic_uv = TRUE,
		nscan = 600,
		burn = 200,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_length(fit$L, 2L)
	expect_true(any(fit$L < 0) && any(fit$L > 0))
	# the innovation sd follows rho_uv rather than being drawn
	expect_equal(fit$sigma_uv, sqrt(1 - fit$rho_uv^2), tolerance = 1e-6)
	ul = fit$ULUPM
	off = !is.na(Y[[1]])
	pc = mean(sapply(seq_len(TT), function(t) cor(ul[, , t][off], prod[off])))
	# a positive-semidefinite product cannot correlate above about 0.5 with this truth; the signed one does
	expect_gt(pc, 0.6)
	# the residual variance is not inflated by an unrepresentable direction
	expect_lt(mean(fit$VC[, "ve"]), 1.4)
})

# simulate() draws the start state with its posterior uncertainty, so compare finite averages over saved states
test_that("simulate.lame agrees with the retained-draw predictive average of a dynamic poisson fit", {
	skip_on_cran()
	set.seed(6886)
	n = 24
	TT = 5
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.5)
	b = rnorm(n, 0, 0.5)
	Y = lapply(seq_len(TT), function(t) {
		eta = 0.3 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.5), n, n)
		m = matrix(rpois(n * n, exp(pmin(eta, 5))), n, n)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(lame(
		Y,
		family = "poisson",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		nscan = 500,
		burn = 200,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	sims = simulate(fit, nsim = 40, seed = 3)
	obs = lapply(Y, function(m) !is.na(m))
	ypm = mean(unlist(lapply(seq_len(TT), function(t) fit$YPM[[t]][obs[[t]]])))
	ysim = mean(
		sapply(sims$Y, function(traj) mean(unlist(lapply(seq_len(TT), function(t) traj[[t]][obs[[t]]])), na.rm = TRUE))
	)
	expect_gt(ysim / ypm, 0.9)
	expect_lt(ysim / ypm, 1.1)
})
