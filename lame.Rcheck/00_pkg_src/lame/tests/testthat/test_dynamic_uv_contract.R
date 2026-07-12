library(lame)

test_that("dynamic_uv_kind cannot request a transition when dynamic_uv is off", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_error(
		lame(Y, family = "normal", R = 1, dynamic_uv_kind = "snap",
		     nscan = 10, burn = 2, odens = 1, gof = FALSE, verbose = FALSE),
		"dynamic_uv = TRUE"
	)
	expect_error(
		lame(Y, family = "normal", R = 1, dynamic_uv_kind = "t",
		     nscan = 10, burn = 2, odens = 1, gof = FALSE, verbose = FALSE),
		"dynamic_uv = TRUE"
	)
})

test_that("dynamic_uv requires positive latent rank", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_error(
		lame(Y, family = "normal", R = 0, dynamic_uv = TRUE,
		     nscan = 10, burn = 2, odens = 1, gof = FALSE, verbose = FALSE),
		"positive multiplicative latent rank"
	)
})

test_that("MCMC bipartite snap/t dynamic UV requests error", {
	Y = replicate(2, matrix(rnorm(24), 6, 4), simplify = FALSE)

	expect_error(
		lame(Y, family = "normal", R = 1, mode = "bipartite",
		     dynamic_uv = TRUE, dynamic_uv_kind = "snap",
		     nscan = 10, burn = 2, odens = 1, gof = FALSE, verbose = FALSE),
		"bipartite MCMC path"
	)
	expect_error(
		lame(Y, family = "normal", R = 1, mode = "bipartite",
		     dynamic_uv = TRUE, dynamic_uv_kind = "t",
		     nscan = 10, burn = 2, odens = 1, gof = FALSE, verbose = FALSE),
		"bipartite MCMC path"
	)
})

test_that("ALS bipartite dynamic UV accepts t transition", {
	Y = replicate(2, matrix(rnorm(24), 6, 4), simplify = FALSE)

	fit = suppressWarnings(lame(
		Y, family = "normal", R = 1, mode = "bipartite",
		method = "als", dynamic_uv = TRUE, dynamic_uv_kind = "t",
		als_max_iter = 3, verbose = FALSE))
	expect_s3_class(fit, "lame_dynamic_als")
	expect_equal(fit$dynamic_uv_kind, "t")
	expect_true(isTRUE(fit$dynamic$t_model_estimated))
	expect_equal(dim(fit$U), c(6L, 1L, 2L))
	expect_equal(dim(fit$V), c(4L, 1L, 2L))
	expect_equal(dim(fit$lambda_u), c(6L, 2L))
	expect_equal(dim(fit$lambda_v), c(4L, 2L))
	expect_true(all(is.finite(fit$lambda_u)))
	expect_true(all(is.finite(fit$lambda_v)))
})

test_that("method = als routes supported smooth dynamic UV transitions", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_s3_class(
		suppressWarnings(lame(Y, family = "normal", R = 1, method = "als",
		                      dynamic_uv = TRUE, nscan = 10, burn = 2,
		                      odens = 1, verbose = FALSE)),
		"lame_dynamic_als"
	)
	expect_s3_class(
		suppressWarnings(lame(Y, family = "normal", R = 1, method = "als",
		                      dynamic_uv = TRUE, dynamic_uv_kind = "t",
		                      nscan = 10, burn = 2, odens = 1,
		                      verbose = FALSE)),
		"lame_dynamic_als"
	)
	expect_s3_class(
		lame(Y, family = "normal", R = 0, method = "als",
		     dynamic_ab = TRUE, verbose = FALSE),
		"lame_dynamic_als"
	)
})

test_that("front-door dynamic ALS handles aliases, ignored controls, and update calls", {
	Y = replicate(2, {
		m = matrix(rbinom(25, 1, 0.35), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_warning(
		expect_warning(
			{
				fit = lame(Y, family = "bin", R = 0, method = "als",
				           dynamic_ab = TRUE, nscan = 10, burn = 2, odens = 1,
				           keep_snap_draws = "draws", time_index = c(1, 3),
				           period_exposure = c(1, 2), checkpoint_path = tempfile(),
				           als_max_iter = 2, verbose = FALSE)
			},
			"did not meet"),
		"ignored MCMC-only arguments")
	expect_s3_class(fit, "lame_dynamic_als")
	expect_equal(fit$family, "binary")
	expect_null(fit$snap_draws)

	fit_static = lame(Y, family = "normal", R = 0, method = "als",
	                   als_max_iter = 2, verbose = FALSE)
	expect_s3_class(suppressWarnings(update(fit_static, R = 1,
	                                        als_max_iter = 2,
	                                        verbose = FALSE)), "ame_als")
	expect_error(
		lame(Y, family = "normal", R = 0, method = "als",
		     intercept = FALSE, verbose = FALSE),
		"intercept = FALSE"
	)
})

test_that("ALS dynamic UV requires positive rank and supports symmetric and bipartite panels", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)
	Ysym = lapply(Y, function(m) {
		ms = (m + t(m)) / 2
		diag(ms) = NA
		ms
	})

	expect_error(
		lame(Y, family = "normal", R = 0, method = "als",
		     dynamic_uv = TRUE, verbose = FALSE),
		"positive multiplicative latent rank"
	)
		fs = suppressWarnings(lame(Ysym, family = "normal", R = 1,
		                            method = "als", symmetric = TRUE,
		                            dynamic_uv = TRUE, verbose = FALSE))
		expect_s3_class(fs, "lame_dynamic_als")
		expect_true(isTRUE(fs$symmetric))
		expect_equal(fs$EZ[[1]], t(fs$EZ[[1]]), tolerance = 1e-8)
		O_from_factors = array(NA_real_, dim = dim(fs$Oarr))
		dimnames(O_from_factors) = dimnames(fs$Oarr)
		for (tt in seq_len(dim(fs$Oarr)[3L])) {
			O_from_factors[, , tt] = fs$U[, , tt] %*% t(fs$V[, , tt])
			expect_lte(sum(svd(fs$Oarr[, , tt])$d > 1e-8), 1L)
		}
		expect_equal(fs$Oarr, O_from_factors, tolerance = 1e-8)

	fb = suppressWarnings(lame(replicate(2, matrix(rnorm(20), 5, 4),
	                                      simplify = FALSE),
	                            family = "normal", R = 1, mode = "bipartite",
	                            method = "als", dynamic_uv = TRUE,
	                            verbose = FALSE))
	expect_s3_class(fb, "lame_dynamic_als")
	expect_equal(dim(fb$U), c(5L, 1L, 2L))
	expect_equal(dim(fb$V), c(4L, 1L, 2L))
	fsa = suppressWarnings(lame(Ysym, family = "normal", R = 1,
	                             method = "als", symmetric = TRUE,
	                             dynamic_uv = TRUE, dynamic_ab = TRUE,
	                             als_max_iter = 8, verbose = FALSE))
	expect_s3_class(fsa, "lame_dynamic_als")
	expect_true(isTRUE(fsa$dynamic_uv))
	expect_true(isTRUE(fsa$dynamic_ab))
	expect_equal(fsa$b_dynamic, fsa$a_dynamic)
	expect_equal(fsa$EZ[[1]], t(fsa$EZ[[1]]), tolerance = 1e-8)
})

test_that("dynamic wrappers expose complete signatures and run", {
	U = array(rnorm(4 * 1 * 2), dim = c(4, 1, 2))
	V = array(rnorm(4 * 1 * 2), dim = c(4, 1, 2))
	ET = array(rnorm(4 * 4 * 2), dim = c(4, 4, 2))

	expect_true("delta_u_current" %in% names(formals(rUV_dynamic_snap_fc)))
	expect_true("delta_v_current" %in% names(formals(rUV_dynamic_snap_fc)))
	out = expect_silent(
		rUV_dynamic_snap_fc(U, V, ET, rho_uv = 0.5, sigma_uv = 1,
		                    s2 = 1, kappa = 2, pi_snap = 0.1,
		                    shrink = FALSE, symmetric = FALSE)
	)
	expect_equal(dim(out$delta_u), c(4L, 2L))
	expect_equal(dim(out$delta_v), c(4L, 2L))
	expect_equal(out$delta_u[, 1], rep(0, 4))

	out_t = expect_silent(
		rUV_dynamic_t_fc(U, V, ET, rho_uv = 0.5, sigma_uv = 1,
		                  s2 = 1, nu = 4, shrink = FALSE,
		                  symmetric = FALSE)
	)
	expect_equal(dim(out_t$lambda_u), c(4L, 2L))
	expect_equal(dim(out_t$lambda_v), c(4L, 2L))
	ET_bad = array(rnorm(3 * 4 * 2), dim = c(3, 4, 2))
	expect_error(
		rUV_dynamic_fc(U, V, ET_bad, rho_uv = 0.5, sigma_uv = 1,
		               s2 = 1, shrink = FALSE, symmetric = FALSE),
		"nrow\\(U\\).*nrow\\(V\\).*n_time")
	expect_error(
		rUV_dynamic_snap_fc(U, V, ET_bad, rho_uv = 0.5, sigma_uv = 1,
		                    s2 = 1, kappa = 2, pi_snap = 0.1,
		                    shrink = FALSE, symmetric = FALSE),
		"nrow\\(U\\).*nrow\\(V\\).*n_time")
	expect_error(
		rUV_dynamic_t_fc(U, V, ET_bad, rho_uv = 0.5, sigma_uv = 1,
		                 s2 = 1, nu = 4, shrink = FALSE,
		                 symmetric = FALSE),
		"nrow\\(U\\).*nrow\\(V\\).*n_time")
	V_bad = array(rnorm(4 * 2 * 2), dim = c(4, 2, 2))
	expect_error(
		rUV_dynamic_fc(U, V_bad, ET, rho_uv = 0.5, sigma_uv = 1,
		               s2 = 1, shrink = FALSE, symmetric = FALSE),
		"same latent rank")
	V_time_bad = array(rnorm(4 * 1 * 3), dim = c(4, 1, 3))
	expect_error(
		rUV_dynamic_fc(U, V_time_bad, ET, rho_uv = 0.5, sigma_uv = 1,
		               s2 = 1, shrink = FALSE, symmetric = FALSE),
		"same number of time periods")
})

test_that("lame_als advertises static dynamic metadata", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)
	fit = lame_als(Y, R = 0, family = "normal", verbose = FALSE)

	expect_false(isTRUE(fit$dynamic$uv))
	expect_false(isTRUE(fit$dynamic$ab))
	expect_false(isTRUE(fit$dynamic$G))
	expect_null(fit$snap_prob)
	expect_null(fit$rho_uv)
})

test_that("dynamic stack does not use bare rbeta", {
	src = test_path("../../R/lame.R")
	skip_if_not(file.exists(src))
	lines = readLines(src, warn = FALSE)
	expect_false(any(grepl("[^:]\\brbeta\\s*\\(", lines)))
})
