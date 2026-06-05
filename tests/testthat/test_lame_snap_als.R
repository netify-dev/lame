library(lame)
skip_on_cran()

simulate_snap_als_panel = function(n = 8, Tt = 4, R = 1, rho = 0.85,
                                    sigma = 0.05, snap = FALSE,
                                    seed = 11) {
	set.seed(seed)
	U = array(0, dim = c(n, R, Tt))
	V = array(0, dim = c(n, R, Tt))
	U[, , 1] = rnorm(n * R)
	V[, , 1] = rnorm(n * R)
	for (t in 2:Tt) {
		U[, , t] = rho * U[, , t - 1L] + rnorm(n * R, 0, sigma)
		V[, , t] = rho * V[, , t - 1L] + rnorm(n * R, 0, sigma)
	}
	if (isTRUE(snap)) {
		U[2, , 3] = -3 * sign(U[2, , 2] + 1e-3)
		V[5, , 3] = -3 * sign(V[5, , 2] + 1e-3)
	}
	lapply(seq_len(Tt), function(t) {
		Ut = U[, , t, drop = FALSE]
		Vt = V[, , t, drop = FALSE]
		dim(Ut) = c(n, R)
		dim(Vt) = c(n, R)
		Y = Ut %*% t(Vt) + matrix(rnorm(n * n, 0, 0.03), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
		Y
	})
}

simulate_broad_snap_panel = function(n = 24, Tt = 6, R = 2, rho = 0.985,
                                      sigma = 0.01, noise = 0.04,
                                      snap_t = Tt, seed = 101) {
	set.seed(seed)
	U = array(0, dim = c(n, R, Tt))
	U[, , 1] = matrix(rnorm(n * R, 0, 0.5), n, R)
	for (t in 2:Tt) {
		U[, , t] = rho * U[, , t - 1L] + matrix(rnorm(n * R, 0, sigma), n, R)
	}
	snap_idx = seq_len(ceiling(0.75 * n))
	U[snap_idx, , snap_t] = matrix(rnorm(length(snap_idx) * R, 0, 1.4),
	                                length(snap_idx), R)
	if (snap_t < Tt) {
		for (t in (snap_t + 1L):Tt) {
			U[, , t] = rho * U[, , t - 1L] + matrix(rnorm(n * R, 0, sigma), n, R)
		}
	}
	Y = lapply(seq_len(Tt), function(t) {
		M = U[, , t] %*% t(U[, , t]) + matrix(rnorm(n * n, 0, noise), n, n)
		M = (M + t(M)) / 2
		diag(M) = NA
		rownames(M) = colnames(M) = paste0("a", seq_len(n))
		M
	})
	names(Y) = paste0("t", seq_len(Tt))
	Y
}

simulate_bipartite_snap_panel = function(nr = 7, nc = 6, Tt = 4,
                                          R_row = 1, R_col = 1,
                                          rho = 0.85, sigma = 0.04,
                                          noise = 0.03, snap = FALSE,
                                          seed = 31) {
	set.seed(seed)
	U = array(0, dim = c(nr, R_row, Tt))
	V = array(0, dim = c(nc, R_col, Tt))
	G = matrix(0.9, R_row, R_col)
	U[, , 1] = matrix(rnorm(nr * R_row), nr, R_row)
	V[, , 1] = matrix(rnorm(nc * R_col), nc, R_col)
	for (t in 2:Tt) {
		U[, , t] = rho * U[, , t - 1L] +
			matrix(rnorm(nr * R_row, 0, sigma), nr, R_row)
		V[, , t] = rho * V[, , t - 1L] +
			matrix(rnorm(nc * R_col, 0, sigma), nc, R_col)
	}
	if (isTRUE(snap)) {
		U[2, , 3] = matrix(rnorm(R_row, 0, 2.5), 1, R_row)
		V[4, , 3] = matrix(rnorm(R_col, 0, 2.5), 1, R_col)
	}
	Y = lapply(seq_len(Tt), function(t) {
		Ut = matrix(U[, , t], nr, R_row)
		Vt = matrix(V[, , t], nc, R_col)
		M = 0.1 + Ut %*% G %*% t(Vt) +
			matrix(rnorm(nr * nc, 0, noise), nr, nc)
		rownames(M) = paste0("r", seq_len(nr))
		colnames(M) = paste0("c", seq_len(nc))
		M
	})
	names(Y) = paste0("t", seq_len(Tt))
	Y
}

test_that("lame_snap_als returns dynamic snap outputs for directed networks", {
	Y = simulate_snap_als_panel(snap = FALSE)
	fit = lame_snap_als(Y, R = 1, max_iter = 20, verbose = FALSE)

	expect_s3_class(fit, "lame_snap_als")
	expect_equal(fit$method, "als_snap")
	expect_equal(dim(fit$U), c(8L, 1L, 4L))
	expect_equal(dim(fit$V), c(8L, 1L, 4L))
	expect_equal(dim(fit$snap_prob), c(8L, 4L))
	expect_equal(dim(fit$snap_prob_v), c(8L, 4L))
	expect_true(all(is.na(fit$snap_prob[, 1])))
	expect_true(all(is.na(fit$snap_prob_v[, 1])))
	expect_true(isTRUE(fit$dynamic$snap_model_estimated))
	expect_true(is.numeric(fit$rho_uv))
	expect_true(is.numeric(fit$sigma_uv))
	expect_true(is.numeric(fit$pi_snap))
	expect_false(is.null(fit$convergence))
})

test_that("lame method=als routes snap dynamic UV to lame_snap_als", {
	Y = simulate_snap_als_panel(snap = FALSE)
	fit = lame(Y, R = 1, family = "normal", method = "als",
	            dynamic_uv = TRUE, dynamic_uv_kind = "snap",
	            als_stability = "quick", als_max_iter = 3,
	            als_tol = 1e-12, verbose = FALSE)

	expect_s3_class(fit, "lame_snap_als")
	expect_equal(fit$method, "als_snap")
	expect_lte(fit$iterations, 3L)
	expect_equal(fit$convergence$tolerance, 1e-12)
	expect_true(all(is.na(fit$snap_prob[, 1])))
	expect_equal(fit$stability$preset, "quick")
	expect_equal(nrow(fit$stability$runs), 2L)

	fit_alias = lame(Y, R = 1, family = "nrm", method = "als",
	                  dynamic_uv = TRUE, dynamic_uv_kind = "snap",
	                  als_max_iter = 2, verbose = FALSE)
	expect_s3_class(fit_alias, "lame_snap_als")
	expect_equal(fit_alias$family, "normal")
	expect_s3_class(update(fit_alias, als_max_iter = 2,
	                       verbose = FALSE), "lame_snap_als")
})

test_that("lame snap ALS route forwards snap prior controls", {
	Y = simulate_snap_als_panel(snap = FALSE)
	fit = lame(Y, R = 1, family = "normal", method = "als",
	            dynamic_uv = TRUE, dynamic_uv_kind = "snap",
	            prior = list(snap_kappa = 4, snap_pi_a = 2,
	                         snap_pi_b = 8, rho_uv_mean = 0.97),
	            verbose = FALSE)

	expect_equal(fit$meta$snap_kappa, 4)
	expect_equal(unname(fit$meta$snap_pi_prior), c(2, 8))
	expect_equal(fit$meta$rho_prior_mean, 0.97)
})

test_that("lame snap ALS route warns on ignored MCMC and bootstrap controls", {
	Y = simulate_snap_als_panel(snap = FALSE)
	expect_warning(
		expect_warning(
			{
				fit = lame(Y, R = 1, family = "normal", method = "als",
				           dynamic_uv = TRUE, dynamic_uv_kind = "snap",
				           nscan = 10, burn = 2, odens = 1, bootstrap = 1,
				           als_max_iter = 3, verbose = FALSE)
			},
			"bootstrap.*ignored"),
		"ignored MCMC-only arguments")

	expect_s3_class(fit, "lame_snap_als")
	expect_null(fit$bootstrap_dynamic)
})

test_that("lame_snap_als supports symmetric normal panels", {
	Y = lapply(simulate_snap_als_panel(snap = FALSE), function(m) {
		out = (m + t(m)) / 2
		diag(out) = NA
		out
	})
	fit = lame_snap_als(Y, R = 1, symmetric = TRUE,
	                     max_iter = 20, verbose = FALSE)

	expect_s3_class(fit, "lame_snap_als")
	expect_equal(dim(fit$U), c(8L, 1L, 4L))
	expect_equal(dim(fit$V), c(8L, 1L, 4L))
	expect_null(fit$snap_prob_v)
	expect_true(all(is.na(fit$snap_prob[, 1])))
})

test_that("snap summary helpers label score fallback and handle one-year ranks", {
	prob = matrix(c(NA, 0.2, 0.8, NA, 0.4, 0.6), nrow = 2, byrow = TRUE,
	               dimnames = list(c("a1", "a2"), c("t1", "t2", "t3")))
	fit = list(snap_prob = prob)
	class(fit) = c("lame_snap_als", "lame_als", "ame_als")

	expect_warning({
		idx = snap_index_draws(fit)
	},
		"ALS score row"
	)
	expect_equal(attr(idx, "snap_draw_source"), "als_score")
	expect_true(is.na(idx$snap_index[idx$year == "t1"]))

	expect_warning({
		rank_one = snap_rank_summary(fit, years = "t2")
	},
		"ALS score row"
	)
	expect_equal(nrow(rank_one), 1L)
	expect_equal(rank_one$year, "t2")
	expect_equal(rank_one$prob_rank_1, 1)
})

test_that("lame_snap_als pads named panels with changing actor composition", {
	full = simulate_snap_als_panel(n = 6, Tt = 3, snap = FALSE)
	Y = list(
		t1 = full[[1]][paste0("a", 1:5), paste0("a", 1:5)],
		t2 = full[[2]][paste0("a", 1:6), paste0("a", 1:6)],
		t3 = full[[3]][paste0("a", 2:6), paste0("a", 2:6)]
	)

	fit = lame_snap_als(Y, R = 1, max_iter = 20, verbose = FALSE)

	expect_equal(dim(fit$U), c(6L, 1L, 3L))
	expect_equal(rownames(fit$snap_prob), paste0("a", 1:6))
	expect_true(isTRUE(fit$balanced_actor_set))
	expect_equal(dim(fit$actor_presence), c(6L, 3L))
	expect_false(fit$actor_presence["a6", "t1"])
	expect_false(fit$actor_presence["a1", "t3"])
	expect_true(is.na(fit$snap_prob["a6", "t1"]))
	expect_true(is.na(fit$snap_prob["a1", "t3"]))
	expect_true(is.na(fit$snap_prob["a6", "t2"]))
})

test_that("lame_snap_als aligns changing-composition Xdyad by names", {
	full = simulate_snap_als_panel(n = 5, Tt = 3, snap = FALSE)
	Y = list(
		t1 = full[[1]][paste0("a", 1:4), paste0("a", 1:4)],
		t2 = full[[2]][paste0("a", 2:5), paste0("a", 2:5)],
		t3 = full[[3]][paste0("a", c(1, 3, 5)), paste0("a", c(1, 3, 5))]
	)
	make_x = function(y, offset = 0) {
		actors = rownames(y)
		val = matrix(seq_along(y) + offset, nrow(y), ncol(y),
		              dimnames = list(rev(actors), rev(actors)))
		array(val, c(nrow(y), ncol(y), 1L),
		      dimnames = list(rev(actors), rev(actors), "trade"))
	}
	X = list(t1 = make_x(Y$t1, 0), t2 = make_x(Y$t2, 100),
	          t3 = make_x(Y$t3, 200))

	fit = lame_snap_als(Y, Xdyad = X, R = 1, max_iter = 4,
	                     verbose = FALSE)

	expect_equal(fit$X["a1", "a4", "trade_dyad", "t1"],
	             X$t1["a1", "a4", "trade"])
	expect_equal(fit$X["a5", "a2", "trade_dyad", "t2"],
	             X$t2["a5", "a2", "trade"])
	expect_equal(fit$X["a3", "a1", "trade_dyad", "t3"],
	             X$t3["a3", "a1", "trade"])
	expect_true(all(fit$X["a2", , "trade_dyad", "t3"] == 0))
})

test_that("classification-only snap convergence requires eligible transitions", {
	full = simulate_snap_als_panel(n = 6, Tt = 3, snap = FALSE)
	Y = list(
		t1 = full[[1]][paste0("a", 1:3), paste0("a", 1:3)],
		t2 = full[[2]][paste0("a", 4:6), paste0("a", 4:6)],
		t3 = full[[3]][paste0("a", 1:3), paste0("a", 1:3)]
	)

	fit = suppressWarnings(lame_snap_als(
		Y, R = 1, max_iter = 5, snap_convergence = "classification",
		verbose = FALSE))

	expect_false(isTRUE(fit$converged))
	expect_false(isTRUE(fit$convergence$snap_stable))
	expect_equal(sum(is.finite(fit$snap_prob)), 0L)
	expect_true(is.na(summary(fit)$mean_snap))
	expect_true(is.na(summary(fit)$max_snap))
})

test_that("lame_snap_als rejects unsupported cases", {
	Y = simulate_snap_als_panel(snap = FALSE)
	expect_error(lame_snap_als(Y, R = 1, family = "binary", verbose = FALSE),
	             "normal")
	expect_error(lame_snap_als(Y, R = 0, verbose = FALSE),
	             "positive integer")
	expect_error(lame_snap_als(Y, R = 1, Xrow = list(matrix(1, 8, 1)),
	                           verbose = FALSE),
	             "node covariates")
})

test_that("lame_snap_als returns bipartite row and column snap outputs", {
	Y = simulate_bipartite_snap_panel(snap = FALSE)
	fit = lame_snap_als(Y, R = 1, mode = "bipartite",
	                     max_iter = 8, verbose = FALSE)

	expect_s3_class(fit, "lame_snap_als")
	expect_equal(fit$mode, "bipartite")
	expect_equal(fit$method, "als_snap")
	expect_equal(fit$R_row, 1L)
	expect_equal(fit$R_col, 1L)
	expect_equal(dim(fit$U), c(7L, 1L, 4L))
	expect_equal(dim(fit$V), c(6L, 1L, 4L))
	expect_equal(dim(fit$G), c(1L, 1L))
	expect_equal(dim(fit$snap_prob), c(7L, 4L))
	expect_equal(dim(fit$snap_prob_v), c(6L, 4L))
	expect_equal(dim(fit$surface_jump), c(7L, 4L))
	expect_equal(dim(fit$surface_jump_v), c(6L, 4L))
	expect_true(all(is.na(fit$snap_prob[, 1])))
	expect_true(all(is.na(fit$snap_prob_v[, 1])))
	expect_equal(names(fit$rho_uv), c("row", "col"))
	expect_equal(names(fit$sigma_uv), c("row", "col"))
	expect_equal(names(fit$pi_snap), c("row", "col"))
	expect_false(is.null(fit$G_diagnostics$singular_values))
	expect_true(isTRUE(fit$dynamic$snap_model_estimated))
})

test_that("snap ALS print, tidy, and list-newdata predict keep dyadic names", {
	Y = simulate_bipartite_snap_panel(snap = FALSE)
	X = lapply(Y, function(y) {
		array(matrix(rnorm(length(y)), nrow(y), ncol(y),
		             dimnames = dimnames(y)),
		      c(nrow(y), ncol(y), 1L),
		      dimnames = list(rownames(y), colnames(y), "trade"))
	})
	fit = lame_snap_als(Y, Xdyad = X, R = 1, mode = "bipartite",
	                     max_iter = 4, verbose = FALSE)

	expect_error(print(fit), NA)
	expect_error(print(summary(fit)), NA)
	td = tidy(fit, conf.int = FALSE)
	expect_true("trade_dyad" %in% td$term)
	expect_false(anyNA(td$estimate[td$term == "trade_dyad"]))
	expect_equal(predict(fit, newdata = X, type = "link"), fit$EZ,
	             tolerance = 1e-8)
})

test_that("lame method=als routes bipartite snap dynamic UV to lame_snap_als", {
	Y = simulate_bipartite_snap_panel(snap = FALSE)
	fit = lame(Y, R = 0, R_row = 1, R_col = 1,
	            family = "normal", mode = "bipartite",
	            method = "als",
	            dynamic_uv = TRUE, dynamic_uv_kind = "snap",
	            als_stability = "quick", als_max_iter = 4,
	            verbose = FALSE)

	expect_s3_class(fit, "lame_snap_als")
	expect_equal(fit$mode, "bipartite")
	expect_equal(fit$R_row, 1L)
	expect_equal(fit$R_col, 1L)
	expect_lte(fit$iterations, 4L)
	expect_true(all(is.na(fit$snap_prob[, 1])))
	expect_true(all(is.na(fit$snap_prob_v[, 1])))
	expect_equal(fit$stability$preset, "quick")
	expect_equal(nrow(fit$stability$runs), 2L)
})

test_that("bipartite snap ALS handles named changing composition panels", {
	full = simulate_bipartite_snap_panel(nr = 6, nc = 5, Tt = 3,
	                                      snap = FALSE)
	Y = list(
		t1 = full[[1]][paste0("r", 1:5), paste0("c", 1:5)],
		t2 = full[[2]][paste0("r", 1:6), paste0("c", 1:5)],
		t3 = full[[3]][paste0("r", 2:6), paste0("c", 2:5)]
	)

	fit = lame_snap_als(Y, R = 1, mode = "bipartite",
	                     max_iter = 8, verbose = FALSE)

	expect_true(isTRUE(fit$changing_composition))
	expect_equal(dim(fit$U), c(6L, 1L, 3L))
	expect_equal(dim(fit$V), c(5L, 1L, 3L))
	expect_equal(rownames(fit$snap_prob), paste0("r", 1:6))
	expect_equal(rownames(fit$snap_prob_v), paste0("c", 1:5))
	expect_false(fit$row_presence["r6", "t1"])
	expect_false(fit$row_presence["r1", "t3"])
	expect_false(fit$col_presence["c1", "t3"])
	expect_true(is.na(fit$snap_prob["r6", "t1"]))
	expect_true(is.na(fit$snap_prob["r6", "t2"]))
	expect_true(is.na(fit$snap_prob["r1", "t3"]))
	expect_true(is.na(fit$snap_prob_v["c1", "t3"]))
})

test_that("bipartite snap ALS records hyperparameter-update controls", {
	Y = simulate_bipartite_snap_panel(snap = FALSE)
	fit = lame_snap_als(
		Y, R = 1, mode = "bipartite", max_iter = 4,
		hyper_update = "em", drift_quantile = 0.2,
		drift_min_transitions = 2, sigma_floor_fraction = 0.25,
		verbose = FALSE)

	expect_equal(fit$meta$hyper_update, "em")
	expect_equal(fit$meta$drift_quantile, 0.2)
	expect_equal(fit$meta$drift_min_transitions, 2)
	expect_equal(fit$meta$sigma_floor_fraction, 0.25)
})

test_that("planted bipartite row and column snaps rank first", {
	Y = simulate_bipartite_snap_panel(nr = 10, nc = 9, Tt = 5,
	                                   snap = TRUE, noise = 0.015,
	                                   seed = 44)
	fit = lame_snap_als(
		Y, R = 1, mode = "bipartite", max_iter = 40,
		estimate_rho_uv = FALSE, estimate_sigma_uv = FALSE,
		rho_uv = c(0.85, 0.85), sigma_uv = c(0.04, 0.04),
		snap_pi_prior = c(1, 99), snap_kappa = 3,
		verbose = FALSE)

	expect_equal(unname(rank(-fit$snap_prob[, 3],
	                         ties.method = "min")[2]), 1)
	expect_equal(unname(rank(-fit$snap_prob_v[, 3],
	                         ties.method = "min")[4]), 1)
	expect_gt(fit$snap_prob[2, 3], 0.9)
	expect_gt(fit$snap_prob_v[4, 3], 0.9)
})

test_that("snap ALS stability preset attaches start-sensitivity diagnostics", {
	Y = simulate_bipartite_snap_panel(snap = FALSE)
	fit = lame_snap_als(Y, R = 1, mode = "bipartite",
	                     max_iter = 4, stability = "quick",
	                     verbose = FALSE)

	expect_equal(fit$stability$preset, "quick")
	expect_equal(nrow(fit$stability$runs), 2L)
	expect_true(all(c("snap_cor", "q95_abs_snap_diff",
	                  "class_change_share", "fitted_rmse",
	                  "top_period_agrees") %in%
	                names(fit$stability$runs)))
	expect_type(fit$convergence$start_stable, "logical")
})

test_that("drift-only snap ALS panel has low snap scores", {
	Y = simulate_snap_als_panel(snap = FALSE)
	fit = lame_snap_als(Y, R = 1, max_iter = 30,
	                     snap_pi_prior = c(1, 999),
	                     snap_kappa = 2, verbose = FALSE)

	expect_lt(mean(fit$snap_prob, na.rm = TRUE), 0.05)
	expect_lt(max(fit$snap_prob, na.rm = TRUE), 0.2)
})

test_that("planted snap actor-years rank above non-planted transitions", {
	Y = simulate_snap_als_panel(snap = TRUE)
	fit = lame_snap_als(Y, R = 1, max_iter = 40,
	                     snap_pi_prior = c(1, 99),
	                     snap_kappa = 3, verbose = FALSE)

	expect_equal(unname(rank(-fit$snap_prob[, 3], ties.method = "min")[2]), 1)
	expect_equal(unname(rank(-fit$snap_prob_v[, 3], ties.method = "min")[5]), 1)
	expect_gt(fit$snap_prob[2, 3], 0.75)
	expect_gt(fit$snap_prob_v[5, 3], 0.75)
})

test_that("broad terminal snap year is not shifted one period earlier", {
	Y = simulate_broad_snap_panel(snap_t = 6)
	fit = lame_snap_als(Y, R = 2, symmetric = TRUE, max_iter = 60,
	                     estimate_rho_uv = FALSE, estimate_sigma_uv = FALSE,
	                     rho_uv = 0.985, sigma_uv = 0.01,
	                     snap_kappa = 2, snap_pi_prior = c(1, 9),
	                     verbose = FALSE)
	idx = colMeans(fit$snap_prob, na.rm = TRUE)

	expect_equal(unname(rank(-idx, ties.method = "min")["t6"]), 1)
	expect_gt(idx["t6"], 0.7)
	expect_lt(idx["t5"], 0.5)
	expect_true(isTRUE(fit$convergence$snap_stable) || !is.na(fit$convergence$final_max_snap_delta))
	expect_false(is.null(fit$transition_diagnostics$profile_delta_obj))
})

test_that("broad non-terminal snap year is ranked at the planted year", {
	Y = simulate_broad_snap_panel(snap_t = 5)
	fit = lame_snap_als(Y, R = 2, symmetric = TRUE, max_iter = 60,
	                     estimate_rho_uv = FALSE, estimate_sigma_uv = FALSE,
	                     rho_uv = 0.985, sigma_uv = 0.01,
	                     snap_kappa = 2, snap_pi_prior = c(1, 9),
	                     verbose = FALSE)
	idx = colMeans(fit$snap_prob, na.rm = TRUE)

	expect_equal(unname(rank(-idx, ties.method = "min")["t5"]), 1)
	expect_gt(idx["t5"], 0.7)
	expect_lt(idx["t6"], 0.5)
})

test_that("estimated robust hyperparameters recover a broad high-density snap", {
	Y = simulate_broad_snap_panel(n = 32, Tt = 8, snap_t = 8,
	                               seed = 108)
	fit = lame_snap_als(Y, R = 2, symmetric = TRUE, max_iter = 50,
	                     verbose = FALSE, seed = 123)
	idx = colMeans(fit$snap_prob, na.rm = TRUE)

	expect_equal(unname(rank(-idx, ties.method = "min")["t8"]), 1)
	expect_gt(idx["t8"], 0.7)
	expect_lt(fit$sigma_uv, 0.05)
	expect_gt(fit$meta$sigma_update_floor, 1e-4)
	expect_true(all(c("snap_delta_q95", "snap_delta_q99",
	                  "snap_delta_conv", "n_snap_class_changed") %in%
	                 names(fit$param_trace)))
	expect_true(all(c("final_snap_delta_q95", "final_snap_delta_q99",
	                  "final_snap_delta_conv",
	                  "final_n_snap_class_changed",
	                  "max_snap_stable") %in%
	                 names(fit$convergence)))
	expect_s3_class(fit$convergence$unstable_transitions, "data.frame")
	expect_true(all(c("actor", "time", "delta", "prev_score",
	                  "snap_score", "class_changed") %in%
	                 names(fit$convergence$unstable_transitions)))
})

test_that("snap ALS scores align with MCMC snap probabilities on a small panel", {
	skip_on_cran()
	Y = simulate_snap_als_panel(n = 6, Tt = 3, snap = TRUE, seed = 21)
	fit_als = lame_snap_als(Y, R = 1, max_iter = 30,
	                         snap_pi_prior = c(1, 99),
	                         snap_kappa = 3, align = "global",
	                         verbose = FALSE)
	fit_mcmc = suppressMessages(lame(
		Y, R = 1, family = "normal", mode = "unipartite",
		dynamic_uv = TRUE, dynamic_uv_kind = "snap",
		nscan = 120, burn = 60, odens = 5,
		keep_snap_draws = "draws",
		gof = FALSE, plot = FALSE,
		verbose = FALSE, seed = 123))

	expect_false(is.null(fit_mcmc$snap_prob))
	expect_false(is.null(fit_mcmc$snap_prob_v))
	expect_equal(dim(fit_mcmc$snap_draws), c(24L, 6L, 3L))
	expect_equal(dim(fit_mcmc$snap_draws_v), c(24L, 6L, 3L))
	expect_true(all(is.na(fit_mcmc$snap_draws[, , 1L])))
	expect_true(all(is.na(fit_mcmc$snap_draws_v[, , 1L])))
	expect_equal(fit_mcmc$snap_prob,
	             apply(fit_mcmc$snap_draws, c(2L, 3L), mean, na.rm = TRUE),
	             tolerance = 1e-12)
	expect_equal(fit_mcmc$snap_prob_v,
	             apply(fit_mcmc$snap_draws_v, c(2L, 3L), mean, na.rm = TRUE),
	             tolerance = 1e-12)
	expect_equal(fit_mcmc$snap_draws_meta$dim_order,
	             c("draw", "actor", "time"))
	expect_false(is.null(fit_mcmc$snap_draws_summary))
	expect_equal(fit_mcmc$snap_draws_summary$mean, fit_mcmc$snap_prob)
	idx = snap_index_draws(fit_mcmc, actors = 1:3, years = 2:3)
	expect_equal(nrow(idx), 24L * 2L)
	expect_true(all(c("draw", "year", "snap_index") %in% names(idx)))
	idx_sum = snap_index_summary(fit_mcmc, actors = 1:3, years = 2:3)
	expect_equal(nrow(idx_sum), 2L)
	expect_true(all(c("mean", "q2_5", "q97_5") %in% names(idx_sum)))
	cat_sum = snap_category_summary(
		fit_mcmc,
		groups = list(focal = c("a1", "a2"), others = c("a3", "a4")),
		years = 2:3)
	expect_equal(sort(unique(cat_sum$group)), c("focal", "others"))
	rank_sum = snap_rank_summary(fit_mcmc, actors = 1:4, years = 2:3)
	expect_true(all(c("prob_rank_1", "prob_top_3", "rank_mean") %in%
	                names(rank_sum)))
	expect_gt(stats::cor(as.vector(fit_als$snap_prob),
	                     as.vector(fit_mcmc$snap_prob),
	                     use = "complete.obs"), 0.5)
	expect_gt(stats::cor(as.vector(fit_als$snap_prob_v),
	                     as.vector(fit_mcmc$snap_prob_v),
	                     use = "complete.obs"), 0.5)
	expect_equal(unname(rank(-fit_als$snap_prob[, 3], ties.method = "min")[2]), 1)
	expect_equal(unname(rank(-fit_mcmc$snap_prob[, 3], ties.method = "min")[2]), 1)
})

test_that("MCMC snap draw summaries respect changing actor composition", {
	Y_full = simulate_snap_als_panel(n = 4, Tt = 3, snap = FALSE, seed = 99)
	Y = list(
		t1 = Y_full[[1]][paste0("a", 1:3), paste0("a", 1:3)],
		t2 = Y_full[[2]][paste0("a", 1:4), paste0("a", 1:4)],
		t3 = Y_full[[3]][paste0("a", 2:4), paste0("a", 2:4)]
	)
	fit = suppressWarnings(suppressMessages(lame(
		Y, R = 1, family = "normal", mode = "unipartite",
		dynamic_uv = TRUE, dynamic_uv_kind = "snap",
		nscan = 30, burn = 10, odens = 5,
		keep_snap_draws = "draws",
		gof = FALSE, plot = FALSE,
		verbose = FALSE, seed = 321)))

	expect_true(all(is.na(fit$snap_prob[, "t1"])))
	expect_true(is.na(fit$snap_prob["a4", "t2"])) # entry transition
	expect_true(is.na(fit$snap_prob["a1", "t3"])) # absent at t3
	expect_true(is.finite(fit$snap_prob["a2", "t2"]))
	expect_equal(fit$snap_draws_summary$n["a4", "t2"], 0)
	expect_equal(fit$snap_draws_summary$n["a1", "t3"], 0)
	expect_true(all(is.na(fit$snap_draws[, "a4", "t2"])))
	expect_true(all(is.na(fit$snap_draws[, "a1", "t3"])))
})
