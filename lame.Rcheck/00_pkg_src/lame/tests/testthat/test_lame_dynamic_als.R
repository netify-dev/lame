library(lame)
skip_on_cran()

test_that("lame method als fits Gaussian dynamic additive effects", {
	set.seed(11)
	n = 8L
	Tt = 4L
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		a = seq(-1, 1, length.out = n) * sin(t / 2)
		b = rev(a) * cos(t / 3)
		Yt = 0.25 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}

	fit = lame(Y, family = "normal", R = 0, method = "als",
	            dynamic_ab = TRUE, verbose = FALSE)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_ab))
	expect_equal(dim(fit$a_dynamic), c(n, Tt))
	expect_equal(dim(fit$b_dynamic), c(n, Tt))
	expect_lt(max(abs(colMeans(fit$a_dynamic))), 1e-8)
	expect_lt(max(abs(colMeans(fit$b_dynamic))), 1e-8)
	expect_equal(length(fitted(fit)), Tt)
	conv_trace = fit$diagnostics$convergence_trace
	expect_s3_class(conv_trace, "data.frame")
	expect_named(
		conv_trace,
		c("iter", "objective", "relative_objective_delta",
		  "max_fitted_delta", "max_parameter_delta")
	)
	expect_equal(nrow(conv_trace), fit$iterations)
	expect_equal(utils::tail(conv_trace$objective, 1L), fit$deviance)
	expect_equal(fit$diagnostics$last_max_fitted_delta,
	             utils::tail(conv_trace$max_fitted_delta, 1L))
	expect_equal(fit$diagnostics$convergence_tolerance, 5e-5)
})

test_that("dynamic ALS fits symmetric dynamic additive effects", {
	set.seed(111)
	n = 7L
	Tt = 4L
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		a = seq(-1, 1, length.out = n) * cos(t / 2)
		Yt = 0.2 + outer(a, a, "+") + matrix(rnorm(n * n, 0, 0.08), n, n)
		Yt = (Yt + t(Yt)) / 2
		diag(Yt) = NA
		Y[[t]] = Yt
	}

	fit = lame(Y, family = "normal", R = 0, method = "als",
	            symmetric = TRUE, dynamic_ab = TRUE,
	            als_max_iter = 30, verbose = FALSE)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$symmetric))
	expect_true(isTRUE(fit$dynamic_ab))
	expect_equal(dim(fit$a_dynamic), c(n, Tt))
	expect_equal(fit$b_dynamic, fit$a_dynamic)
	expect_lt(max(abs(colMeans(fit$a_dynamic))), 1e-8)
	expect_equal(fit$EZ[[2L]], t(fit$EZ[[2L]]), tolerance = 1e-8)
})

test_that("lame method als fits selected dynamic dyadic beta", {
	set.seed(12)
	n = 10L
	Tt = 5L
	beta_true = seq(-1, 1, length.out = Tt)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Xt = matrix(rnorm(n * n), n, n)
		diag(Xt) = 0
		X[[t]] = array(Xt, c(n, n, 1L),
		                dimnames = list(NULL, NULL, "trade"))
		Yt = 0.1 + beta_true[t] * Xt + matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}

	fit = lame(Y, Xdyad = X, family = "normal", R = 0, method = "als",
	            dynamic_beta = "dyad", verbose = FALSE)
	cf = coef(fit)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_beta))
	expect_equal(dim(cf), c(2L, Tt))
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_gt(diff(range(cf["trade_dyad", ])), 1)
	expect_gt(stats::cor(cf["trade_dyad", ], beta_true), 0.95)
})

test_that("dynamic ALS can warm-start dynamic UV MCMC", {
	set.seed(901)
	n = 7L
	Tt = 3L
	R = 1L
	actors = paste0("a", seq_len(n))
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		u = sin(t / 2) * seq(-1, 1, length.out = n)
		v = cos(t / 3) * rev(seq(-1, 1, length.out = n))
		Yt = 0.1 + tcrossprod(matrix(u, ncol = 1), matrix(v, ncol = 1)) +
			matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(Yt) = NA
		rownames(Yt) = colnames(Yt) = actors
		Y[[t]] = Yt
	}
	names(Y) = paste0("t", seq_len(Tt))

	fit_als = lame(Y, family = "normal", R = R, dynamic_uv = TRUE,
	               method = "als", verbose = FALSE)
	start = als_start_vals(fit_als, jitter = 0.01, seed = 902)
	fit_mcmc = lame(Y, family = "normal", R = R, dynamic_uv = TRUE,
	                start_vals = start, nscan = 20, burn = 5, odens = 5,
	                gof = FALSE, verbose = FALSE, seed = 903,
	                prior = list(uv_max_abs = 10))

	expect_s3_class(fit_mcmc, "lame")
	expect_equal(dim(start$U), c(n, R, Tt))
	expect_equal(dim(start$V), c(n, R, Tt))
	expect_equal(dim(fit_mcmc$start_vals$U), c(n, R, Tt))
	expect_equal(dim(fit_mcmc$start_vals$V), c(n, R, Tt))
	expect_true("rho_uv" %in% names(fit_mcmc$start_vals))
	expect_true("sigma_uv" %in% names(fit_mcmc$start_vals))
	expect_lt(abs(mean(fit_mcmc$APM)), 1e-8)
	expect_lt(abs(mean(fit_mcmc$BPM)), 1e-8)
	expect_lte(max(abs(c(fit_mcmc$start_vals$U, fit_mcmc$start_vals$V))), 10)
})

test_that("dynamic ALS reports residual reciprocity for directed normal panels", {
	set.seed(22)
	n = 9L
	Tt = 4L
	rho_true = 0.5
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Yt = matrix(NA_real_, n, n)
		a = seq(-0.7, 0.7, length.out = n) + rnorm(n, 0, 0.05)
		b = rev(a) + rnorm(n, 0, 0.05)
		for (i in seq_len(n - 1L)) {
			for (j in (i + 1L):n) {
				z1 = rnorm(1, 0, 0.25)
				z2 = rnorm(1, 0, 0.25)
				Yt[i, j] = 0.2 + a[i] + b[j] + z1
				Yt[j, i] = 0.2 + a[j] + b[i] +
					rho_true * z1 + sqrt(1 - rho_true^2) * z2
			}
		}
		Y[[t]] = Yt
	}

	fit = lame(Y, family = "normal", R = 0, method = "als",
	           dynamic_ab = TRUE, als_max_iter = 200,
	           als_tol = 1e-4, verbose = FALSE, seed = 23)
	start = als_start_vals(fit)

	expect_true(isTRUE(fit$converged))
	expect_true(is.finite(fit$VC["rho"]))
	expect_gt(fit$VC["rho"], 0.25)
	expect_equal(length(fit$rho_path), Tt)
	expect_equal(unname(start$rho), unname(fit$rho_path))
})

test_that("dynamic rho MCMC stores a per-period reciprocity path", {
	set.seed(30)
	n = 8L
	Tt = 5L
	actors = paste0("n", seq_len(n))
	rho_true = seq(0.1, 0.7, length.out = Tt)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Yt = matrix(NA_real_, n, n, dimnames = list(actors, actors))
		a = rnorm(n, 0, 0.3)
		b = rnorm(n, 0, 0.3)
		for (i in seq_len(n - 1L)) {
			for (j in (i + 1L):n) {
				z1 = rnorm(1, 0, 0.4)
				z2 = rnorm(1, 0, 0.4)
				Yt[i, j] = 0.1 + a[i] + b[j] + z1
				Yt[j, i] = 0.1 + a[j] + b[i] +
					rho_true[t] * z1 + sqrt(1 - rho_true[t]^2) * z2
			}
		}
		Y[[t]] = Yt
	}
	names(Y) = paste0("t", seq_len(Tt))

	fit = lame(Y, R = 1, family = "normal", dynamic_uv = TRUE,
	           dynamic_rho = TRUE, nscan = 40, burn = 20,
	           odens = 5, gof = TRUE, verbose = FALSE, seed = 31)

	expect_true(isTRUE(fit$dynamic_rho))
	expect_equal(dim(fit$RHO), c(8L, Tt))
	expect_equal(length(fit$rho_path), Tt)
	expect_equal(length(fit$start_vals$rho), Tt)
	expect_true(all(is.finite(fit$rho_path)))
})

test_that("dynamic ALS fits symmetric dynamic dyadic beta", {
	set.seed(112)
	n = 8L
	Tt = 4L
	beta_true = seq(-0.8, 0.8, length.out = Tt)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Xraw = matrix(rnorm(n * n), n, n)
		diag(Xraw) = 0
		Xsym = (Xraw + t(Xraw)) / 2
		X[[t]] = array(Xraw, c(n, n, 1L),
		                dimnames = list(NULL, NULL, "trade"))
		Yt = -0.1 + beta_true[t] * Xsym +
			matrix(rnorm(n * n, 0, 0.08), n, n)
		Yt = (Yt + t(Yt)) / 2
		diag(Yt) = NA
		Y[[t]] = Yt
	}

	fit = suppressWarnings(lame(Y, Xdyad = X, family = "normal", R = 0,
	                             method = "als", symmetric = TRUE,
	                             dynamic_beta = "dyad",
	                             als_max_iter = 30, verbose = FALSE))
	cf = coef(fit)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$symmetric))
	expect_true(isTRUE(fit$dynamic_beta))
	expect_equal(dim(cf), c(2L, Tt))
	expect_equal(fit$X[, , 1L, 1L], t(fit$X[, , 1L, 1L]))
	expect_equal(fit$EZ[[3L]], t(fit$EZ[[3L]]), tolerance = 1e-8)
	expect_gt(stats::cor(cf["trade_dyad", ], beta_true), 0.8)
})

test_that("dynamic ALS supports node covariates by orthogonal decomposition", {
	set.seed(113)
	n = 8L
	Tt = 3L
	w_row = seq(-1, 1, length.out = n)
	w_col = rev(w_row)
	Xrow = replicate(Tt, matrix(w_row, n, 1,
	                             dimnames = list(NULL, "sender")), simplify = FALSE)
	Xcol = replicate(Tt, matrix(w_col, n, 1,
	                             dimnames = list(NULL, "receiver")), simplify = FALSE)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Xt = matrix(rnorm(n * n), n, n)
		diag(Xt) = 0
		X[[t]] = array(Xt, c(n, n, 1L),
		                dimnames = list(NULL, NULL, "trade"))
		a = sin(seq_len(n) + t) / 4
		b = cos(seq_len(n) - t) / 5
		Yt = 0.1 + 0.6 * Xt + 1.1 * w_row - 0.7 * rep(w_col, each = n) +
			outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.08), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}

	fit = lame(Y, Xdyad = X, Xrow = Xrow, Xcol = Xcol,
	            family = "normal", R = 0, method = "als",
	            dynamic_ab = TRUE, als_max_iter = 20, verbose = FALSE)
	X4 = array(0, dim = c(n, n, 1L, Tt))
	for (t in seq_len(Tt)) X4[, , 1L, t] = X[[t]][, , 1L]
	pred = predict(fit, newdata = X4, type = "link")

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true("sender_row" %in% names(fit$coefficients))
	expect_true("receiver_col" %in% names(fit$coefficients))
	expect_lt(max(abs(crossprod(fit$W_row, rowMeans(fit$a_dynamic)))), 1e-8)
	expect_lt(max(abs(crossprod(fit$W_col, rowMeans(fit$b_dynamic)))), 1e-8)
	expect_equal(pred, fit$EZ, tolerance = 1e-8)
})

test_that("dynamic ALS aligns named unipartite panels with changing actor composition", {
	set.seed(301)
	actors = paste0("a", 1:5)
	sets = list(t1 = actors[1:4], t2 = actors[2:5],
	             t3 = actors[c(1, 3, 4, 5)])
	Y = X = Xrow = vector("list", length(sets))
	names(Y) = names(X) = names(Xrow) = names(sets)
	for (tt in seq_along(sets)) {
		a = sets[[tt]]
		n = length(a)
		x = matrix(rnorm(n * n), n, n, dimnames = list(a, a))
		diag(x) = 0
		w = seq(-1, 1, length.out = n) + tt / 10
		names(w) = a
		eta = 0.1 + 0.4 * x + 0.6 * matrix(w, n, n) +
			matrix(rnorm(n * n, 0, 0.05), n, n)
		diag(eta) = NA_real_
		dimnames(eta) = list(a, a)
		Y[[tt]] = eta
		X[[tt]] = array(x, c(n, n, 1L),
		                 dimnames = list(a, a, "trade"))
		Xrow[[tt]] = matrix(w, n, 1L, dimnames = list(a, "sender"))
	}

	fit = suppressWarnings(lame(
		Y, Xdyad = X, Xrow = Xrow, family = "normal", R = 0,
		method = "als", dynamic_ab = TRUE,
		dynamic_beta = c("dyad", "row"),
		prior = list(lambda_beta_als = 0),
		als_max_iter = 20, verbose = FALSE, seed = 301))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$changing_composition))
	expect_equal(dim(fit$Y), c(5L, 5L, 3L))
	expect_equal(rownames(fit$actor_presence), actors)
	expect_false(fit$actor_presence["a5", "t1"])
	expect_false(fit$actor_presence["a1", "t2"])
	expect_true(all(is.na(fit$Y["a5", , "t1"])))
	expect_true(all(is.na(fit$Y[, "a5", "t1"])))
	expect_true(all(c("trade_dyad", "sender_row") %in% rownames(fit$coef_path)))
	expect_equal(dim(fit$a_dynamic), c(5L, 3L))
	expect_equal(rownames(fit$a_dynamic), actors)
})

test_that("dynamic ALS aligns named bipartite panels with changing actor composition", {
	set.seed(302)
	row_sets = list(t1 = paste0("r", 1:4), t2 = paste0("r", 2:5),
	                 t3 = paste0("r", c(1, 3, 5)))
	col_sets = list(t1 = paste0("c", 1:3), t2 = paste0("c", 2:4),
	                 t3 = paste0("c", c(1, 3, 4)))
	Y = vector("list", 3L)
	names(Y) = names(row_sets)
	for (tt in seq_along(Y)) {
		rn = row_sets[[tt]]
		cn = col_sets[[tt]]
		nr = length(rn)
		nc = length(cn)
		u = seq(-0.6, 0.6, length.out = nr)
		v = seq(0.5, -0.5, length.out = nc)
		y = 0.2 + outer(u, v, "+") +
			0.4 * tcrossprod(u, v) + matrix(rnorm(nr * nc, 0, 0.05), nr, nc)
		dimnames(y) = list(rn, cn)
		Y[[tt]] = y
	}

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R = 1, R_row = 1, R_col = 1, method = "als",
		dynamic_uv = TRUE, als_max_iter = 6,
		verbose = FALSE, seed = 302))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$changing_composition))
	expect_equal(dim(fit$Y), c(5L, 4L, 3L))
	expect_equal(rownames(fit$row_presence), paste0("r", 1:5))
	expect_equal(rownames(fit$col_presence), paste0("c", 1:4))
	expect_false(fit$row_presence["r5", "t1"])
	expect_false(fit$col_presence["c1", "t2"])
	expect_equal(dim(fit$U), c(5L, 1L, 3L))
	expect_equal(dim(fit$V), c(4L, 1L, 3L))
	expect_true(all(is.na(fit$Y["r5", , "t1"])))
	expect_true(all(is.na(fit$Y[, "c1", "t2"])))
})

test_that("dynamic ALS supports bipartite dynamic_G with changing actor composition", {
	set.seed(303)
	row_sets = list(t1 = paste0("r", 1:4), t2 = paste0("r", 2:5),
	                 t3 = paste0("r", c(1, 3, 5)))
	col_sets = list(t1 = paste0("c", 1:3), t2 = paste0("c", 2:4),
	                 t3 = paste0("c", c(1, 3, 4)))
	Y = vector("list", 3L)
	names(Y) = names(row_sets)
	for (tt in seq_along(Y)) {
		rn = row_sets[[tt]]
		cn = col_sets[[tt]]
		Y[[tt]] = matrix(rnorm(length(rn) * length(cn)),
		                  length(rn), length(cn),
		                  dimnames = list(rn, cn))
	}

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R = 1, R_row = 1, R_col = 1, method = "als",
		dynamic_G = TRUE, als_max_iter = 4,
		verbose = FALSE, seed = 303))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$changing_composition))
	expect_equal(dim(fit$Y), c(5L, 4L, 3L))
	expect_equal(dim(fit$G_cube), c(1L, 1L, 3L))
	expect_false(fit$row_presence["r5", "t1"])
	expect_false(fit$col_presence["c1", "t2"])
})

test_that("dynamic ALS rejects unnamed changing-composition panels", {
	Y = list(
		matrix(rnorm(9), 3, 3),
		matrix(rnorm(16), 4, 4)
	)
	expect_error(
		lame(Y, family = "normal", R = 0, method = "als",
		     dynamic_ab = TRUE, verbose = FALSE),
		"requires row and column names")
})

test_that("dynamic ALS estimates unipartite row and column node beta paths", {
	set.seed(117)
	n = 10L
	Tt = 5L
	w_row = as.numeric(scale(seq(-1, 1, length.out = n), scale = FALSE))
	w_col = as.numeric(scale(rev(w_row), scale = FALSE))
	beta_row = seq(-1.2, 1.2, length.out = Tt)
	beta_col = seq(0.9, -0.9, length.out = Tt)
	Xrow = replicate(Tt, matrix(w_row, n, 1,
	                             dimnames = list(NULL, "sender")), simplify = FALSE)
	Xcol = replicate(Tt, matrix(w_col, n, 1,
	                             dimnames = list(NULL, "receiver")), simplify = FALSE)
	Y = lapply(seq_len(Tt), function(t) {
		y = 0.1 + beta_row[t] * matrix(w_row, n, n) +
			beta_col[t] * matrix(w_col, n, n, byrow = TRUE) +
			matrix(rnorm(n * n, 0, 0.04), n, n)
		diag(y) = NA
		y
	})

	fit = suppressWarnings(lame(
		Y, Xrow = Xrow, Xcol = Xcol, family = "normal", R = 0,
		method = "als", dynamic_beta = c("row", "col"),
		prior = list(lambda_beta_als = 0),
		als_max_iter = 30, verbose = FALSE, seed = 117))
	cf = coef(fit)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_beta))
	expect_true(all(c("sender_row", "receiver_col") %in% rownames(cf)))
	expect_true(fit$beta_dynamic_mask["sender_row"])
	expect_true(fit$beta_dynamic_mask["receiver_col"])
	expect_true(fit$beta_dynamic_mask["intercept"])
	expect_gt(diff(range(cf["sender_row", ])), 1)
	expect_gt(diff(range(cf["receiver_col", ])), 0.8)
	expect_gt(stats::cor(cf["sender_row", ], beta_row), 0.85)
	expect_gt(stats::cor(cf["receiver_col", ], beta_col), 0.85)
})

test_that("dynamic ALS estimates bipartite row and column node beta paths", {
	set.seed(118)
	nr = 8L
	nc = 7L
	Tt = 4L
	w_row = as.numeric(scale(seq(-1, 1, length.out = nr), scale = FALSE))
	w_col = as.numeric(scale(seq(-0.8, 0.8, length.out = nc), scale = FALSE))
	beta_row = seq(-0.8, 0.8, length.out = Tt)
	beta_col = seq(1.1, -0.5, length.out = Tt)
	Xrow = replicate(Tt, matrix(w_row, nr, 1,
	                             dimnames = list(NULL, "exporter")), simplify = FALSE)
	Xcol = replicate(Tt, matrix(w_col, nc, 1,
	                             dimnames = list(NULL, "importer")), simplify = FALSE)
	Y = lapply(seq_len(Tt), function(t) {
		beta_row[t] * matrix(w_row, nr, nc) +
			beta_col[t] * matrix(w_col, nr, nc, byrow = TRUE) +
			matrix(rnorm(nr * nc, 0, 0.04), nr, nc)
	})

	fit = suppressWarnings(lame(
		Y, Xrow = Xrow, Xcol = Xcol, family = "normal", mode = "bipartite",
		R = 0, method = "als", dynamic_beta = c("row", "col"),
		prior = list(lambda_beta_als = 0),
		als_max_iter = 30, verbose = FALSE, seed = 118))
	cf = coef(fit)

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(identical(fit$mode, "bipartite"))
	expect_true(all(c("exporter_row", "importer_col") %in% rownames(cf)))
	expect_true(fit$beta_dynamic_mask["exporter_row"])
	expect_true(fit$beta_dynamic_mask["importer_col"])
	expect_equal(dim(fit$BETA), c(1L, nrow(cf), Tt))
	expect_gt(diff(range(cf["exporter_row", ])), 0.6)
	expect_gt(diff(range(cf["importer_col", ])), 0.6)
	expect_gt(stats::cor(cf["exporter_row", ], beta_row), 0.8)
	expect_gt(stats::cor(cf["importer_col", ], beta_col), 0.8)
})

test_that("dynamic ALS uses period-specific row node values selected by dynamic_beta", {
	set.seed(202)
	n = 9L
	Tt = 5L
	actors = paste0("a", seq_len(n))
	beta_row = seq(-1, 1, length.out = Tt)
	Xrow = vector("list", Tt)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	eta_true = array(NA_real_, c(n, n, Tt))
	for (t in seq_len(Tt)) {
		w = as.numeric(scale(rnorm(n), scale = FALSE))
		Xrow[[t]] = matrix(w, n, 1,
		                     dimnames = list(actors, "sender"))
		xd = matrix(rnorm(n * n), n, n)
		diag(xd) = 0
		X[[t]] = array(xd, c(n, n, 1L),
		                dimnames = list(actors, actors, "noise"))
		eta = 0.2 + beta_row[t] * matrix(w, n, n)
		diag(eta) = NA
		eta_true[, , t] = eta
		y = eta + matrix(rnorm(n * n, 0, 0.03), n, n)
		diag(y) = NA
		dimnames(y) = list(actors, actors)
		Y[[t]] = y
	}

	fit = NULL
	expect_no_warning({
		fit = lame(
			Y, Xdyad = X, Xrow = Xrow, family = "normal", R = 0,
			method = "als", dynamic_beta = "row",
			prior = list(lambda_beta_als = 0),
			als_max_iter = 40, verbose = FALSE, seed = 202)
	})
	cf = coef(fit)
	X4 = array(0, dim = c(n, n, 1L, Tt),
	            dimnames = list(actors, actors, "noise", NULL))
	for (t in seq_len(Tt)) X4[, , 1L, t] = X[[t]][, , 1L]

	expect_true(fit$beta_dynamic_mask["sender_row"])
	expect_equal(dim(fit$W_row_arr), c(n, 1L, Tt))
	expect_equal(as.numeric(fit$W_row_arr[, 1L, 2L]),
	             as.numeric(Xrow[[2L]][, 1L]))
	expect_gt(stats::cor(cf["sender_row", ], beta_row), 0.99)
	expect_gt(stats::cor(unlist(fit$EZ), as.vector(eta_true),
	                     use = "complete.obs"), 0.99)
	expect_equal(predict(fit, newdata = X4, type = "link"),
	             fit$EZ, tolerance = 1e-8)
})

test_that("dynamic ALS warns when varying node values are left static", {
	set.seed(203)
	n = 7L
	Tt = 3L
	Xrow = vector("list", Tt)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		w = as.numeric(scale(rnorm(n), scale = FALSE))
		Xrow[[t]] = matrix(w, n, 1, dimnames = list(NULL, "sender"))
		xd = matrix(rnorm(n * n), n, n)
		diag(xd) = 0
		X[[t]] = array(xd, c(n, n, 1L),
		                dimnames = list(NULL, NULL, "trade"))
		y = 0.1 + (t / 5) * xd + matrix(rnorm(n * n, 0, 0.05), n, n)
		diag(y) = NA
		Y[[t]] = y
	}

	fit = NULL
	expect_warning({
		fit = lame(
			Y, Xdyad = X, Xrow = Xrow, family = "normal", R = 0,
			method = "als", dynamic_beta = "dyad",
			als_max_iter = 10, verbose = FALSE, seed = 203)
	}, "per-actor means for static row coefficient")

	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$beta_dynamic_mask["sender_row"])
	expect_true("sender_row" %in% rownames(fit$coef_path))
})

test_that("dynamic ALS uses period-specific bipartite node values", {
	set.seed(204)
	nr = 8L
	nc = 7L
	Tt = 4L
	beta_row = seq(-0.9, 0.9, length.out = Tt)
	beta_col = seq(0.8, -0.4, length.out = Tt)
	Xrow = vector("list", Tt)
	Xcol = vector("list", Tt)
	X = vector("list", Tt)
	Y = vector("list", Tt)
	eta_true = array(NA_real_, c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		wr = as.numeric(scale(rnorm(nr), scale = FALSE))
		wc = as.numeric(scale(rnorm(nc), scale = FALSE))
		Xrow[[t]] = matrix(wr, nr, 1,
		                     dimnames = list(NULL, "exporter"))
		Xcol[[t]] = matrix(wc, nc, 1,
		                     dimnames = list(NULL, "importer"))
		xd = matrix(rnorm(nr * nc), nr, nc)
		X[[t]] = array(xd, c(nr, nc, 1L),
		                dimnames = list(NULL, NULL, "noise"))
		eta = 0.05 + beta_row[t] * matrix(wr, nr, nc) +
			beta_col[t] * matrix(wc, nr, nc, byrow = TRUE)
		eta_true[, , t] = eta
		Y[[t]] = eta + matrix(rnorm(nr * nc, 0, 0.03), nr, nc)
	}

	fit = NULL
	expect_no_warning({
		fit = lame(
			Y, Xdyad = X, Xrow = Xrow, Xcol = Xcol,
			family = "normal", mode = "bipartite", R = 0,
			method = "als", dynamic_beta = c("row", "col"),
			prior = list(lambda_beta_als = 0),
			als_max_iter = 50, verbose = FALSE, seed = 204)
	})
	cf = coef(fit)

	expect_true(fit$beta_dynamic_mask["exporter_row"])
	expect_true(fit$beta_dynamic_mask["importer_col"])
	expect_equal(dim(fit$W_row_arr), c(nr, 1L, Tt))
	expect_equal(dim(fit$W_col_arr), c(nc, 1L, Tt))
	expect_equal(as.numeric(fit$W_row_arr[, 1L, 3L]),
	             as.numeric(Xrow[[3L]][, 1L]))
	expect_equal(as.numeric(fit$W_col_arr[, 1L, 3L]),
	             as.numeric(Xcol[[3L]][, 1L]))
	expect_gt(stats::cor(cf["exporter_row", ], beta_row), 0.99)
	expect_gt(stats::cor(cf["importer_col", ], beta_col), 0.99)
	expect_gt(stats::cor(unlist(fit$EZ), as.vector(eta_true)), 0.99)
})

test_that("dynamic ALS keeps logical dynamic_beta masks compatible with node covariates", {
	set.seed(119)
	n = 8L
	Tt = 3L
	X = replicate(Tt, {
		x = matrix(rnorm(n * n), n, n)
		diag(x) = 0
		array(x, c(n, n, 1L), dimnames = list(NULL, NULL, "trade"))
	}, simplify = FALSE)
	w = seq(-1, 1, length.out = n)
	Xrow = replicate(Tt, matrix(w, n, 1,
	                             dimnames = list(NULL, "sender")), simplify = FALSE)
	Y = lapply(seq_len(Tt), function(t) {
		y = 0.2 + t / 10 * X[[t]][, , 1L] + 0.5 * matrix(w, n, n) +
			matrix(rnorm(n * n, 0, 0.08), n, n)
		diag(y) = NA
		y
	})

	fit = suppressWarnings(lame(
		Y, Xdyad = X, Xrow = Xrow, family = "normal", R = 0,
		method = "als", dynamic_beta = c(FALSE, TRUE),
		als_max_iter = 20, verbose = FALSE, seed = 119))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(fit$beta_dynamic_mask["trade_dyad"])
	expect_false(fit$beta_dynamic_mask["sender_row"])
	expect_true("sender_row" %in% rownames(fit$coef_path))
})

test_that("dynamic ALS prediction handles dynamic beta newdata", {
	set.seed(13)
	n = 7L
	Tt = 3L
	X = replicate(Tt, {
		x = matrix(rnorm(n * n), n, n)
		diag(x) = 0
		array(x, c(n, n, 1L), dimnames = list(NULL, NULL, "x"))
	}, simplify = FALSE)
	Y = lapply(seq_len(Tt), function(t) {
		y = X[[t]][, , 1] + matrix(rnorm(n * n, 0, 0.1), n, n)
		diag(y) = NA
		y
	})
	fit = lame(Y, Xdyad = X, family = "normal", R = 0, method = "als",
	            dynamic_beta = "dyad", verbose = FALSE)
	X4 = array(0, dim = c(n, n, 1L, Tt))
	for (t in seq_len(Tt)) X4[, , 1L, t] = X[[t]][, , 1L]
	p0 = predict(fit, type = "link")
	p1 = predict(fit, newdata = X4, type = "link")
	expect_equal(p0, p1)
})

test_that("lame method als fits Gaussian AR1 dynamic latent factors", {
	set.seed(14)
	n = 7L
	Tt = 4L
	R = 1L
	U = array(0, dim = c(n, R, Tt))
	V = array(0, dim = c(n, R, Tt))
	U[, , 1L] = rnorm(n, 0, 0.7)
	V[, , 1L] = rnorm(n, 0, 0.7)
	for (t in 2:Tt) {
		U[, , t] = 0.8 * U[, , t - 1L] + rnorm(n, 0, 0.15)
		V[, , t] = 0.8 * V[, , t - 1L] + rnorm(n, 0, 0.15)
	}
	actors = paste0("a", seq_len(n))
	Y = lapply(seq_len(Tt), function(t) {
		y = tcrossprod(U[, , t], V[, , t]) +
			matrix(rnorm(n * n, 0, 0.15), n, n)
		diag(y) = NA
		dimnames(y) = list(actors, actors)
		y
	})

	fit = suppressWarnings(lame(Y, family = "normal", R = R,
	                             method = "als", dynamic_uv = TRUE,
	                             dynamic_uv_kind = "ar1",
	                             dynamic_ab = FALSE, verbose = FALSE))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_uv))
	expect_equal(dim(fit$U), c(n, R, Tt))
	expect_equal(dim(fit$V), c(n, R, Tt))
	expect_equal(dim(fit$Oarr), c(n, n, Tt))
	expect_equal(dimnames(fit$U)[[1L]], actors)
	expect_equal(dimnames(fit$V)[[1L]], actors)
	expect_equal(dimnames(fit$Oarr)[[3L]], paste0("t", seq_len(Tt)))
	expect_equal(length(fitted(fit)), Tt)
	expect_true(all(is.finite(fit$Oarr[is.finite(fit$Y)])))
	expect_gt(max(abs(diff(fit$Oarr[1, 2, ]))), 1e-6)
})

test_that("dynamic ALS estimates bipartite dynamic G with static U and V", {
	set.seed(154)
	nr = 7L
	nc = 6L
	Tt = 4L
	R = 2L
	rows = paste0("r", seq_len(nr))
	cols = paste0("c", seq_len(nc))
	U = matrix(rnorm(nr * R, 0, 0.7), nr, R)
	V = matrix(rnorm(nc * R, 0, 0.7), nc, R)
	G = array(0, dim = c(R, R, Tt))
	for (t in seq_len(Tt)) {
		G[, , t] = matrix(c(1 + 0.25 * t, 0.15 * sin(t),
		                     -0.1 * cos(t), 0.6 - 0.08 * t), R, R)
	}
	Y = lapply(seq_len(Tt), function(t) {
		y = U %*% G[, , t] %*% t(V) +
			matrix(rnorm(nr * nc, 0, 0.08), nr, nc)
		dimnames(y) = list(rows, cols)
		y
	})

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite", R = R,
		method = "als", dynamic_G = TRUE,
		als_max_iter = 12, verbose = FALSE, seed = 154))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_G))
	expect_true(isTRUE(fit$dynamic$G))
	expect_false(isTRUE(fit$dynamic_uv))
	expect_equal(dim(fit$U), c(nr, R))
	expect_equal(dim(fit$V), c(nc, R))
	expect_equal(dim(fit$G), c(R, R))
	expect_equal(dim(fit$G_cube), c(R, R, Tt))
	expect_equal(dim(fit$G_cube_post_mean), c(R, R, Tt))
	expect_equal(fit$G_cube_post_mean, fit$G_cube)
	expect_equal(fit$G_cube_post_sd,
	             array(0, dim = c(R, R, Tt), dimnames = dimnames(fit$G_cube)))
	expect_equal(dim(fit$Oarr), c(nr, nc, Tt))
	expect_equal(dim(fit$UVPM), c(nr, nc, Tt))
	expect_equal(dimnames(fit$Oarr)[[3L]], paste0("t", seq_len(Tt)))
	expect_true(all(is.finite(fit$Oarr)))
	expect_true(is.finite(fit$rho_G))
	expect_true(is.finite(fit$RHO_G))
	expect_true(is.finite(fit$SIGMA_G2))
	expect_true(is.finite(fit$lambda_G_als))

	recon = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fit$Oarr))
	for (t in seq_len(Tt)) {
		recon[, , t] = fit$U %*% fit$G_cube[, , t] %*% t(fit$V)
	}
	expect_equal(fit$Oarr, recon, tolerance = 1e-8)
})

test_that("dynamic ALS estimates bipartite dynamic G with dynamic U and V", {
	set.seed(155)
	nr = 6L
	nc = 5L
	Tt = 3L
	R = 1L
	U = array(rnorm(nr * R * Tt, 0, 0.6), c(nr, R, Tt))
	V = array(rnorm(nc * R * Tt, 0, 0.6), c(nc, R, Tt))
	G = array(seq(0.6, 1.2, length.out = Tt), c(R, R, Tt))
	Y = lapply(seq_len(Tt), function(t) {
		Ut = matrix(U[, , t], nrow = nr, ncol = R)
		Vt = matrix(V[, , t], nrow = nc, ncol = R)
		Gt = matrix(G[, , t], nrow = R, ncol = R)
		Ut %*% Gt %*% t(Vt) +
			matrix(rnorm(nr * nc, 0, 0.12), nr, nc)
	})

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite", R = R,
		method = "als", dynamic_uv = TRUE, dynamic_G = TRUE,
		als_max_iter = 6, verbose = FALSE, seed = 155))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_uv))
	expect_true(isTRUE(fit$dynamic_G))
	expect_equal(dim(fit$U), c(nr, R, Tt))
	expect_equal(dim(fit$V), c(nc, R, Tt))
	expect_equal(dim(fit$G_cube), c(R, R, Tt))
	expect_equal(dim(fit$Oarr), c(nr, nc, Tt))
	expect_true(all(is.finite(fit$Oarr)))

	recon = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fit$Oarr))
	for (t in seq_len(Tt)) {
		Ut = matrix(fit$U[, , t], nrow = nr, ncol = R)
		Vt = matrix(fit$V[, , t], nrow = nc, ncol = R)
		Gt = matrix(fit$G_cube[, , t], nrow = R, ncol = R)
		recon[, , t] = Ut %*% Gt %*% t(Vt)
	}
	expect_equal(fit$Oarr, recon, tolerance = 1e-8)
})

test_that("dynamic ALS honors unequal bipartite R_row and R_col", {
	set.seed(156)
	nr = 6L
	nc = 5L
	Tt = 3L
	RA = 2L
	RB = 1L
	U = matrix(rnorm(nr * RA, 0, 0.6), nr, RA)
	V = matrix(rnorm(nc * RB, 0, 0.6), nc, RB)
	G = array(0, c(RA, RB, Tt))
	for (t in seq_len(Tt)) {
		G[, , t] = matrix(c(1 + 0.1 * t, -0.2 + 0.05 * t), RA, RB)
	}
	Y = lapply(seq_len(Tt), function(t) {
		U %*% G[, , t] %*% t(V) +
			matrix(rnorm(nr * nc, 0, 0.1), nr, nc)
	})

	fg = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_G = TRUE, als_max_iter = 5,
		verbose = FALSE, seed = 156))
	expect_s3_class(fg, "lame_dynamic_als")
	expect_equal(fg$R_row, RA)
	expect_equal(fg$R_col, RB)
	expect_equal(dim(fg$U), c(nr, RA))
	expect_equal(dim(fg$V), c(nc, RB))
	expect_equal(dim(fg$G), c(RA, RB))
	expect_equal(dim(fg$G_cube), c(RA, RB, Tt))
	recon_g = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fg$Oarr))
	for (t in seq_len(Tt)) {
		recon_g[, , t] = fg$U %*% fg$G_cube[, , t] %*% t(fg$V)
	}
	expect_equal(fg$Oarr, recon_g, tolerance = 1e-8)

	fuvg = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_uv = TRUE, dynamic_G = TRUE,
		als_max_iter = 4, verbose = FALSE, seed = 157))
	expect_equal(dim(fuvg$U), c(nr, RA, Tt))
	expect_equal(dim(fuvg$V), c(nc, RB, Tt))
	expect_equal(dim(fuvg$G_cube), c(RA, RB, Tt))
	recon_uvg = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fuvg$Oarr))
	for (t in seq_len(Tt)) {
		Ut = matrix(fuvg$U[, , t], nrow = nr, ncol = RA)
		Vt = matrix(fuvg$V[, , t], nrow = nc, ncol = RB)
		Gt = matrix(fuvg$G_cube[, , t], nrow = RA, ncol = RB)
		recon_uvg[, , t] = Ut %*% Gt %*% t(Vt)
	}
	expect_equal(fuvg$Oarr, recon_uvg, tolerance = 1e-8)

	fuv = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_uv = TRUE, als_max_iter = 4,
		verbose = FALSE, seed = 158))
	expect_equal(dim(fuv$U), c(nr, RA, Tt))
	expect_equal(dim(fuv$V), c(nc, RB, Tt))
	expect_equal(dim(fuv$G), c(RA, RB))
	recon_uv = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fuv$Oarr))
	for (t in seq_len(Tt)) {
		Ut = matrix(fuv$U[, , t], nrow = nr, ncol = RA)
		Vt = matrix(fuv$V[, , t], nrow = nc, ncol = RB)
		recon_uv[, , t] = Ut %*% fuv$G %*% t(Vt)
	}
	expect_equal(fuv$Oarr, recon_uv, tolerance = 1e-8)
})

test_that("dynamic ALS supports bipartite t dynamic UV", {
	set.seed(161)
	nr = 6L
	nc = 4L
	Tt = 3L
	RA = 2L
	RB = 1L
	Y = lapply(seq_len(Tt), function(t) {
		matrix(rnorm(nr * nc, mean = 0.15 * t), nr, nc)
	})

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_uv = TRUE, dynamic_uv_kind = "t",
		als_max_iter = 4, verbose = FALSE, seed = 161))
	expect_s3_class(fit, "lame_dynamic_als")
	expect_equal(fit$dynamic_uv_kind, "t")
	expect_equal(fit$dynamic$uv_kind, "t")
	expect_true(isTRUE(fit$dynamic$t_model_estimated))
	expect_equal(fit$R_row, RA)
	expect_equal(fit$R_col, RB)
	expect_equal(dim(fit$U), c(nr, RA, Tt))
	expect_equal(dim(fit$V), c(nc, RB, Tt))
	expect_equal(dim(fit$G), c(RA, RB))
	expect_equal(dim(fit$lambda_u), c(nr, Tt))
	expect_equal(dim(fit$lambda_v), c(nc, Tt))
	expect_equal(unname(fit$lambda_u[, 1L]), rep(1, nr))
	expect_equal(unname(fit$lambda_v[, 1L]), rep(1, nc))
	expect_equal(rownames(fit$lambda_u), rownames(fit$U))
	expect_equal(rownames(fit$lambda_v), rownames(fit$V))
	expect_true(all(is.finite(fit$lambda_u)))
	expect_true(all(is.finite(fit$lambda_v)))
	expect_true(all(fit$lambda_u > 0))
	expect_true(all(fit$lambda_v > 0))
	recon = array(0, dim = c(nr, nc, Tt), dimnames = dimnames(fit$Oarr))
	for (t in seq_len(Tt)) {
		Ut = matrix(fit$U[, , t], nrow = nr, ncol = RA)
		Vt = matrix(fit$V[, , t], nrow = nc, ncol = RB)
		recon[, , t] = Ut %*% fit$G %*% t(Vt)
	}
	expect_equal(fit$Oarr, recon, tolerance = 1e-8)

	fit_g = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_uv = TRUE, dynamic_uv_kind = "t", dynamic_G = TRUE,
		als_max_iter = 3, verbose = FALSE, seed = 162))
	expect_equal(fit_g$dynamic_uv_kind, "t")
	expect_true(isTRUE(fit_g$dynamic$t_model_estimated))
	expect_true(isTRUE(fit_g$dynamic_G))
	expect_equal(dim(fit_g$G_cube), c(RA, RB, Tt))
	expect_equal(dim(fit_g$lambda_u), c(nr, Tt))
	expect_equal(dim(fit_g$lambda_v), c(nc, Tt))
})

test_that("dynamic ALS expands static bipartite low-rank terms across time", {
	set.seed(159)
	nr = 6L
	nc = 4L
	Tt = 3L
	RA = 2L
	RB = 1L
	Y = lapply(seq_len(Tt), function(t) {
		matrix(rnorm(nr * nc, 0.1 * t, 0.4), nr, nc)
	})

	fit = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = RB, method = "als",
		dynamic_ab = TRUE, als_max_iter = 4,
		verbose = FALSE, seed = 159))
	expect_equal(dim(fit$U), c(nr, RA))
	expect_equal(dim(fit$V), c(nc, RB))
	expect_equal(dim(fit$G), c(RA, RB))
	expect_equal(dim(fit$Oarr), c(nr, nc, Tt))
	recon = fit$U %*% fit$G %*% t(fit$V)
	for (t in seq_len(Tt)) {
		expect_equal(fit$Oarr[, , t], recon, tolerance = 1e-8)
	}

	fit_zero = suppressWarnings(lame(
		Y, family = "normal", mode = "bipartite",
		R_row = RA, R_col = 0, method = "als",
		dynamic_ab = TRUE, als_max_iter = 3,
		verbose = FALSE, seed = 160))
	expect_equal(dim(fit_zero$Oarr), c(nr, nc, Tt))
	expect_null(fit_zero$G)
	expect_true(all(is.finite(fit_zero$Oarr)))
})

test_that("dynamic ALS fits binary and poisson IWLS dynamic additive effects", {
	set.seed(15)
	n = 7L
	Tt = 3L
	x = seq(-1, 1, length.out = n)
	eta = outer(x, -x, "+")
	Y_bin = lapply(seq_len(Tt), function(t) {
		y = matrix(rbinom(n * n, 1, stats::plogis(eta + 0.3 * t)), n, n)
		diag(y) = NA
		y
	})
	Y_pois = lapply(seq_len(Tt), function(t) {
		y = matrix(rpois(n * n, exp(0.2 + eta / 3 + 0.1 * t)), n, n)
		diag(y) = NA
		y
	})

	fb = suppressWarnings(lame(Y_bin, family = "binary", R = 0,
	                            method = "als", dynamic_ab = TRUE,
	                            verbose = FALSE))
	fp = suppressWarnings(lame(Y_pois, family = "poisson", R = 0,
	                            method = "als", dynamic_ab = TRUE,
	                            verbose = FALSE))

	expect_s3_class(fb, "lame_dynamic_als")
	expect_s3_class(fp, "lame_dynamic_als")
	expect_equal(fb$non_normal_method, "irls")
	expect_equal(fp$link, "log")
	expect_true(all(fitted(fb)[[1]][is.finite(Y_bin[[1]])] >= 0))
	expect_true(all(fitted(fb)[[1]][is.finite(Y_bin[[1]])] <= 1))
	expect_true(all(fitted(fp)[[1]][is.finite(Y_pois[[1]])] >= 0))
})

test_that("lame method als exposes t dynamic UV and stability controls", {
	set.seed(151)
	n = 5L
	Tt = 3L
	R = 1L
	actors = paste0("a", seq_len(n))
	U = array(rnorm(n * R * Tt, 0, 0.4), c(n, R, Tt))
	V = array(rnorm(n * R * Tt, 0, 0.4), c(n, R, Tt))
	Y = lapply(seq_len(Tt), function(t) {
		y = tcrossprod(U[, , t], V[, , t]) +
			matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(y) = NA
		dimnames(y) = list(actors, actors)
		y
	})

	fit = suppressWarnings(lame(
		Y, family = "normal", R = R, method = "als",
		dynamic_uv = TRUE, dynamic_uv_kind = "t",
		als_stability = "quick", als_max_iter = 4,
		verbose = FALSE, seed = 151))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_uv))
	expect_true(isTRUE(fit$dynamic$t_model_estimated))
	expect_equal(fit$stability$preset, "quick")
	expect_equal(fit$stability$n_start, 2L)
	expect_equal(dim(fit$U), c(n, R, Tt))
	expect_equal(dim(fit$V), c(n, R, Tt))
	expect_equal(dim(fit$lambda_u), c(n, Tt))
	expect_equal(dim(fit$lambda_v), c(n, Tt))
	expect_true(all(is.finite(fit$lambda_u)))
	expect_true(all(is.finite(fit$lambda_v)))
})

test_that("dynamic ALS keeps remaining unsupported stacks explicit", {
	Y = replicate(3, {
		m = matrix(rnorm(36), 6, 6)
		diag(m) = NA
		m
	}, simplify = FALSE)
	expect_error(
		lame(Y, family = "normal", R = 1, method = "als",
		     dynamic_uv = TRUE, dynamic_uv_kind = "snap",
		     dynamic_ab = TRUE, verbose = FALSE),
		"smooth .*ar1.*t"
	)
	expect_error(
		lame(Y, family = "ordinal", R = 0, method = "als",
		     dynamic_ab = TRUE, intercept = TRUE, verbose = FALSE),
		"normal.*binary.*poisson"
	)
	expect_error(
		lame(Y, family = "normal", R = 1, method = "als",
		     dynamic_G = TRUE, verbose = FALSE),
		"bipartite"
	)
})

test_that("dynamic ALS exposes stability and dynamic bootstrap refits", {
	set.seed(16)
	n = 6L
	Tt = 3L
	Y = lapply(seq_len(Tt), function(t) {
		y = matrix(rnorm(n * n), n, n)
		diag(y) = NA
		y
	})

	fit = suppressWarnings(lame(Y, family = "normal", R = 0,
	                             method = "als", dynamic_ab = TRUE,
	                             als_stability = "quick",
	                             bootstrap = 1,
	                             als_max_iter = 5,
	                             verbose = FALSE, seed = 16))

	expect_s3_class(fit, "lame_dynamic_als")
	expect_equal(fit$stability$preset, "quick")
	expect_equal(fit$stability$n_start, 2L)
	expect_true("max_surface_delta" %in% names(fit$stability))
	expect_equal(fit$bootstrap_dynamic$R, 1L)
	expect_true("coef_paths" %in% names(fit$bootstrap_dynamic))
	expect_true(isTRUE(fit$meta$uncertainty_available))
})

test_that("dynamic bootstrap preserves changing-composition missingness", {
	set.seed(17)
	sets = list(t1 = paste0("a", 1:4), t2 = paste0("a", 2:5),
	             t3 = paste0("a", c(1, 3, 5)))
	Y = lapply(seq_along(sets), function(tt) {
		a = sets[[tt]]
		y = matrix(rnorm(length(a)^2), length(a), length(a),
		            dimnames = list(a, a))
		diag(y) = NA_real_
		y
	})
	names(Y) = names(sets)

	fit = suppressWarnings(lame(
		Y, family = "normal", R = 0, method = "als",
		dynamic_ab = TRUE, bootstrap = 1, als_max_iter = 4,
		verbose = FALSE, seed = 17))
	refit = fit$bootstrap_dynamic$refits[[1L]]

	expect_true(isTRUE(fit$changing_composition))
	expect_true(isTRUE(refit$changing_composition))
	expect_equal(is.finite(refit$Y), is.finite(fit$Y))
	expect_false(refit$actor_presence["a5", "t1"])
	expect_false(refit$actor_presence["a2", "t3"])
})

test_that("dynamic UV level prior bounds a null binary panel (defect 5)", {
	# with no real latent structure the difference penalty alone leaves the
	# common level of the factor path unpenalized, so a binary IRLS fit can
	# drive var(UV') to hundreds. the fixed marginal level prior must keep
	# the factor magnitudes bounded and let the fit converge.
	set.seed(42)
	n = 20L
	Tt = 6L
	actors = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(Tt), function(t) {
		z = -0.5 + matrix(rnorm(n * n), n, n)
		y = 1 * (z > 0)
		diag(y) = NA
		dimnames(y) = list(actors, actors)
		y
	})
	fit = suppressWarnings(suppressMessages(lame_dynamic_als(
		Y, R = 2, family = "binary", dynamic_uv = TRUE,
		verbose = FALSE, seed = 6886)))
	expect_true(isTRUE(fit$converged))
	# pre-fix this panel produced var(UV') in the hundreds/thousands
	expect_lt(max(abs(fit$Oarr[is.finite(fit$Oarr)])), 25)
	expect_lt(stats::var(fit$Oarr[is.finite(fit$Oarr)]), 50)
})

test_that("dynamic UV level prior does not degrade real drift recovery (defect 5)", {
	# a genuinely drifting rank-1 structure must still be recovered: the
	# level prior only bounds the magnitude, it must not flatten real signal.
	set.seed(7)
	n = 24L
	Tt = 6L
	actors = sprintf("b%02d", seq_len(n))
	U0 = matrix(rnorm(n, 0, 0.8), n, 1)
	V0 = matrix(rnorm(n, 0, 0.8), n, 1)
	O_true = array(0, c(n, n, Tt))
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		if (t > 1L) {
			U0 = 0.9 * U0 + matrix(rnorm(n, 0, 0.35), n, 1)
			V0 = 0.9 * V0 + matrix(rnorm(n, 0, 0.35), n, 1)
		}
		O_true[, , t] = tcrossprod(U0, V0)
		z = -0.3 + O_true[, , t] + matrix(rnorm(n * n), n, n)
		y = 1 * (z > 0)
		diag(y) = NA
		dimnames(y) = list(actors, actors)
		Y[[t]] = y
	}
	fit = suppressWarnings(suppressMessages(lame_dynamic_als(
		Y, R = 1, family = "binary", dynamic_uv = TRUE,
		verbose = FALSE, seed = 6886)))
	per_t_cor = vapply(seq_len(Tt), function(t) {
		off = row(O_true[, , t]) != col(O_true[, , t])
		stats::cor(O_true[, , t][off], fit$Oarr[, , t][off])
	}, numeric(1))
	# recovery stays strong at every period (pre-fix mean was ~0.28)
	expect_gt(mean(per_t_cor), 0.7)
	expect_gt(min(per_t_cor), 0.6)
})
