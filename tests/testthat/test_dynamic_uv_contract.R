skip_on_cran()

library(lame)

test_that("dynamic_uv_kind cannot request a transition when dynamic_uv is off", {
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_error(
		lame(
			Y,
			family = "normal",
			R = 1,
			dynamic_uv_kind = "snap",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		"dynamic_uv = TRUE"
	)
	expect_error(
		lame(
			Y,
			family = "normal",
			R = 1,
			dynamic_uv_kind = "t",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
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
		lame(
			Y,
			family = "normal",
			R = 0,
			dynamic_uv = TRUE,
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		"positive multiplicative latent rank"
	)
})

test_that("MCMC bipartite snap/t dynamic UV requests error", {
	Y = replicate(2, matrix(rnorm(24), 6, 4), simplify = FALSE)

	expect_error(
		lame(
			Y,
			family = "normal",
			R = 1,
			mode = "bipartite",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		"bipartite MCMC path"
	)
	expect_error(
		lame(
			Y,
			family = "normal",
			R = 1,
			mode = "bipartite",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		"bipartite MCMC path"
	)
	# both messages point at the ar1 kind only, since the bipartite snap and t kinds have no ALS route either
	msg_snap = tryCatch(
		lame(
			Y,
			family = "normal",
			R = 1,
			mode = "bipartite",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		error = function(e) gsub("\\s+", " ", conditionMessage(e))
	)
	expect_match(msg_snap, "dynamic_uv_kind = \"ar1\"", fixed = TRUE)
	expect_no_match(msg_snap, "point estimate", fixed = TRUE)
	msg_t = tryCatch(
		lame(
			Y,
			family = "normal",
			R = 1,
			mode = "bipartite",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			nscan = 10,
			burn = 2,
			odens = 1,
			gof = FALSE,
			verbose = FALSE
		),
		error = function(e) gsub("\\s+", " ", conditionMessage(e))
	)
	expect_match(msg_t, "dynamic_uv_kind = \"ar1\"", fixed = TRUE)
	expect_no_match(msg_t, "point estimate", fixed = TRUE)
})

test_that("ALS bipartite dynamic UV refuses the t transition", {
	Y = replicate(2, matrix(rnorm(24), 6, 4), simplify = FALSE)

	expect_error(
		suppressWarnings(lame(
			Y,
			family = "normal",
			R = 1,
			mode = "bipartite",
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			als_max_iter = 3,
			verbose = FALSE
		)),
		"not available for bipartite networks"
	)
})

test_that("method = als routes supported smooth dynamic UV transitions", {
	skip_on_cran()
	Y = replicate(2, {
		m = matrix(rnorm(25), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	expect_s3_class(
		suppressWarnings(lame(
			Y,
			family = "normal",
			R = 1,
			method = "als",
			dynamic_uv = TRUE,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE
		)),
		"lame_dynamic_als"
	)
	expect_s3_class(
		suppressWarnings(lame(
			Y,
			family = "normal",
			R = 1,
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE
		)),
		"lame_dynamic_als"
	)
	expect_s3_class(
		lame(
			Y,
			family = "normal",
			R = 0,
			method = "als",
			dynamic_ab = TRUE,
			verbose = FALSE
		),
		"lame_dynamic_als"
	)
})

test_that("front-door dynamic ALS handles aliases, ignored controls, and update calls", {
	set.seed(9109)
	Y = replicate(2, {
		m = matrix(rbinom(25, 1, 0.35), 5, 5)
		diag(m) = NA
		m
	}, simplify = FALSE)

	# the aliases are honoured: the family alias resolves, the ALS iteration
	# cap reaches the fitter, and MCMC-only controls are reported as ignored
	warns = character()
	err = withCallingHandlers(
		expect_error(
			lame(
				Y,
				family = "bin",
				R = 0,
				method = "als",
				dynamic_ab = TRUE,
				nscan = 10,
				burn = 2,
				odens = 1,
				keep_snap_draws = "draws",
				time_index = c(1, 3),
				period_exposure = c(1, 2),
				checkpoint_path = tempfile(),
				als_max_iter = 2,
				verbose = FALSE
			),
			class = "lame_als_convergence"
		),
		warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
	)
	expect_true(any(grepl("ignored MCMC-only arguments", warns)))
	# two iterations cannot meet the convergence criterion: the fit is refused and its state is in the condition
	fit = err$fit
	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$converged)
	expect_equal(fit$family, "binary")
	expect_equal(fit$als_refit$args$max_iter, 2)
	expect_lte(length(fit$dev_history), 2L)
	expect_null(fit$snap_draws)

	fit_static = lame(
		Y,
		family = "normal",
		R = 0,
		method = "als",
		als_max_iter = 2,
		bias_correct = "none",
		verbose = FALSE
	)
	expect_s3_class(
		suppressWarnings(update(
			fit_static,
			R = 1,
			als_max_iter = 2,
			verbose = FALSE
		)),
		"ame_als"
	)
	expect_error(
		lame(
			Y,
			family = "normal",
			R = 0,
			method = "als",
			intercept = FALSE,
			verbose = FALSE
		),
		"intercept = FALSE"
	)
})

test_that("ALS dynamic UV requires positive rank and supports symmetric and bipartite panels", {
	skip_on_cran()
	set.seed(9201)
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
		lame(
			Y,
			family = "normal",
			R = 0,
			method = "als",
			dynamic_uv = TRUE,
			verbose = FALSE
		),
		"positive multiplicative latent rank"
	)
	fs = suppressWarnings(lame(
		Ysym,
		family = "normal",
		R = 1,
		method = "als",
		symmetric = TRUE,
		dynamic_uv = TRUE,
		verbose = FALSE
	))
	expect_s3_class(fs, "lame_dynamic_als")
	expect_true(isTRUE(fs$symmetric))
	expect_equal(fs$EZ[[1]], t(fs$EZ[[1]]), tolerance = 1e-8)
	O_from_factors = array(NA_real_, dim = dim(fs$Oarr))
	dimnames(O_from_factors) = dimnames(fs$Oarr)
	for (tt in seq_len(dim(fs$Oarr)[3L])) {
		O_from_factors[, , tt] = fs$U[, , tt] %*% t(fs$V[, , tt])
		# the product of an undirected fit carries the shift of the intercept's posterior mean: the reported product is the factors' less that constant
		expect_lte(sum(svd(fs$Oarr[, , tt] + fs$eb$centre[["o"]])$d > 1e-8), 1L)
	}
	expect_equal(fs$Oarr + fs$eb$centre[["o"]], O_from_factors, tolerance = 1e-8)

	fb = suppressWarnings(lame(
		replicate(
			2,
			matrix(rnorm(20), 5, 4),
			simplify = FALSE
		),
		family = "normal",
		R = 1,
		mode = "bipartite",
		method = "als",
		dynamic_uv = TRUE,
		verbose = FALSE
	))
	expect_s3_class(fb, "lame_dynamic_als")
	expect_equal(dim(fb$U), c(5L, 1L, 2L))
	expect_equal(dim(fb$V), c(4L, 1L, 2L))
	fsa = suppressWarnings(lame(
		Ysym,
		family = "normal",
		R = 1,
		method = "als",
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		verbose = FALSE
	))
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
		rUV_dynamic_snap_fc(
			U,
			V,
			ET,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			kappa = 2,
			pi_snap = 0.1,
			shrink = FALSE,
			symmetric = FALSE
		)
	)
	expect_equal(dim(out$delta_u), c(4L, 2L))
	expect_equal(dim(out$delta_v), c(4L, 2L))
	expect_equal(out$delta_u[, 1], rep(0, 4))

	out_t = expect_silent(
		rUV_dynamic_t_fc(
			U,
			V,
			ET,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			nu = 4,
			shrink = FALSE,
			symmetric = FALSE
		)
	)
	expect_equal(dim(out_t$lambda_u), c(4L, 2L))
	expect_equal(dim(out_t$lambda_v), c(4L, 2L))
	ET_bad = array(rnorm(3 * 4 * 2), dim = c(3, 4, 2))
	expect_error(
		rUV_dynamic_fc(
			U,
			V,
			ET_bad,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			shrink = FALSE,
			symmetric = FALSE
		),
		"nrow\\(U\\).*nrow\\(V\\).*n_time"
	)
	expect_error(
		rUV_dynamic_snap_fc(
			U,
			V,
			ET_bad,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			kappa = 2,
			pi_snap = 0.1,
			shrink = FALSE,
			symmetric = FALSE
		),
		"nrow\\(U\\).*nrow\\(V\\).*n_time"
	)
	expect_error(
		rUV_dynamic_t_fc(
			U,
			V,
			ET_bad,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			nu = 4,
			shrink = FALSE,
			symmetric = FALSE
		),
		"nrow\\(U\\).*nrow\\(V\\).*n_time"
	)
	V_bad = array(rnorm(4 * 2 * 2), dim = c(4, 2, 2))
	expect_error(
		rUV_dynamic_fc(
			U,
			V_bad,
			ET,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			shrink = FALSE,
			symmetric = FALSE
		),
		"same latent rank"
	)
	V_time_bad = array(rnorm(4 * 1 * 3), dim = c(4, 1, 3))
	expect_error(
		rUV_dynamic_fc(
			U,
			V_time_bad,
			ET,
			rho_uv = 0.5,
			sigma_uv = 1,
			s2 = 1,
			shrink = FALSE,
			symmetric = FALSE
		),
		"same number of time periods"
	)
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

# snap-shift MCMC draws and their summaries; the planted panel resets actor 2's sender and actor 5's receiver
# positions in period 3
simulate_snap_panel = function(
	n = 8,
	Tt = 4,
	R = 1,
	rho = 0.85,
	sigma = 0.05,
	snap = FALSE,
	seed = 11
) {
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

test_that("the ALS route refuses the snap kind in every shape", {
	# the snap ALS point estimate does not track the per-period product or the planted resets well enough to report,
	# so the message names the routes that do fit
	flat = function(e) gsub("\\s+", " ", conditionMessage(e))
	Yd = simulate_snap_panel(snap = FALSE)
	msg_d = tryCatch(
		lame(
			Yd,
			R = 1,
			family = "normal",
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			verbose = FALSE
		),
		error = flat
	)
	expect_match(msg_d, "is not available", fixed = TRUE)
	expect_match(msg_d, "dynamic_uv_kind = \"snap\"` with `method = \"mcmc\"", fixed = TRUE)
	expect_match(msg_d, "dynamic_uv_kind = \"ar1\"` with `method = \"als\"", fixed = TRUE)
	expect_match(msg_d, "misses resets", fixed = TRUE)
	Ys = lapply(Yd, function(m) { out = (m + t(m)) / 2; diag(out) = NA; out })
	msg_s = tryCatch(
		lame(
			Ys,
			R = 1,
			family = "normal",
			symmetric = TRUE,
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			verbose = FALSE
		),
		error = flat
	)
	expect_match(msg_s, "symmetric snap panels", fixed = TRUE)
	expect_match(msg_s, "dynamic_uv_kind = \"ar1\"", fixed = TRUE)
	expect_no_match(msg_s, "\"snap\"` with `method = \"mcmc\"", fixed = TRUE)
	Yb = replicate(3, matrix(rnorm(24), 6, 4), simplify = FALSE)
	msg_b = tryCatch(
		lame(
			Yb,
			R = 1,
			family = "normal",
			mode = "bipartite",
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			dynamic_ab = TRUE,
			verbose = FALSE
		),
		error = flat
	)
	expect_match(msg_b, "bipartite snap panels", fixed = TRUE)
	expect_error(
		lame:::lame_dynamic_als(
			Yd,
			family = "normal",
			R = 1,
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			verbose = FALSE
		),
		"no ALS route"
	)
	# lame_snap_als() stays exported, as a refusal naming the same routes
	expect_true("lame_snap_als" %in% getNamespaceExports("lame"))
	msg_x = tryCatch(lame_snap_als(Yd, R = 1, verbose = FALSE), error = flat)
	expect_match(msg_x, "was removed in lame 1.4.0", fixed = TRUE)
	expect_match(msg_x, "dynamic_uv_kind = \"snap\"` with `method = \"mcmc\"", fixed = TRUE)
	msg_xb = tryCatch(
		lame_snap_als(Yb, R_row = 1, R_col = 1, mode = "bipartite"),
		error = flat
	)
	expect_match(msg_xb, "bipartite snap panels", fixed = TRUE)
	msg_xs = tryCatch(lame_snap_als(Ys, R = 1, symmetric = TRUE), error = flat)
	expect_match(msg_xs, "symmetric snap panels", fixed = TRUE)
})

test_that("snap summary helpers fall back to the posterior-mean row and handle one-year ranks", {
	prob = matrix(
		c(NA, 0.2, 0.8, NA, 0.4, 0.6),
		nrow = 2,
		byrow = TRUE,
		dimnames = list(c("a1", "a2"), c("t1", "t2", "t3"))
	)
	fit = list(snap_prob = prob)
	class(fit) = "lame"

	expect_warning(
		{
			idx = snap_index_draws(fit)
		},
		"posterior-mean row"
	)
	expect_equal(attr(idx, "snap_draw_source"), "posterior_mean")
	expect_true(is.na(idx$snap_index[idx$year == "t1"]))

	expect_warning(
		{
			rank_one = snap_rank_summary(fit, years = "t2")
		},
		"posterior-mean row"
	)
	expect_equal(nrow(rank_one), 1L)
	expect_equal(rank_one$year, "t2")
	expect_equal(rank_one$prob_rank_1, 1)
})

test_that("MCMC snap draws and their summaries agree on a small planted panel", {
	skip_on_cran()
	Y = simulate_snap_panel(n = 6, Tt = 3, snap = TRUE, seed = 21)
	# both independently seeded chains identify the planted reset actors; the short chain's still-moving warning is not the contract here
	fit_mcmc = suppressWarnings(suppressMessages(lame(
		Y,
		R = 1,
		family = "normal",
		mode = "unipartite",
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		nscan = 600,
		burn = 300,
		odens = 5,
		keep_snap_draws = "draws",
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		seed = 123
	)))
	expect_false(is.null(fit_mcmc$snap_prob))
	expect_false(is.null(fit_mcmc$snap_prob_v))
	expect_equal(dim(fit_mcmc$snap_draws), c(120L, 6L, 3L))
	expect_equal(dim(fit_mcmc$snap_draws_v), c(120L, 6L, 3L))
	expect_true(all(is.na(fit_mcmc$snap_draws[, , 1L])))
	expect_true(all(is.na(fit_mcmc$snap_draws_v[, , 1L])))
	expect_equal(
		fit_mcmc$snap_prob,
		apply(fit_mcmc$snap_draws, c(2L, 3L), mean, na.rm = TRUE),
		tolerance = 1e-12
	)
	expect_equal(
		fit_mcmc$snap_prob_v,
		apply(fit_mcmc$snap_draws_v, c(2L, 3L), mean, na.rm = TRUE),
		tolerance = 1e-12
	)
	expect_equal(
		fit_mcmc$snap_draws_meta$dim_order,
		c("draw", "actor", "time")
	)
	expect_false(is.null(fit_mcmc$snap_draws_summary))
	expect_equal(fit_mcmc$snap_draws_summary$mean, fit_mcmc$snap_prob)
	idx = snap_index_draws(fit_mcmc, actors = 1:3, years = 2:3)
	expect_equal(nrow(idx), 120L * 2L)
	expect_true(all(c("draw", "year", "snap_index") %in% names(idx)))
	idx_sum = snap_index_summary(fit_mcmc, actors = 1:3, years = 2:3)
	expect_equal(nrow(idx_sum), 2L)
	expect_true(all(c("mean", "q2_5", "q97_5") %in% names(idx_sum)))
	cat_sum = snap_category_summary(
		fit_mcmc,
		groups = list(focal = c("a1", "a2"), others = c("a3", "a4")),
		years = 2:3
	)
	expect_equal(sort(unique(cat_sum$group)), c("focal", "others"))
	rank_sum = snap_rank_summary(fit_mcmc, actors = 1:4, years = 2:3)
	expect_true(all(c("prob_rank_1", "prob_top_3", "rank_mean") %in%
			names(rank_sum)))
	expect_equal(unname(rank(-fit_mcmc$snap_prob[, 3], ties.method = "min")[2]), 1)
	fit_second = suppressMessages(lame(
		Y,
		R = 1,
		family = "normal",
		mode = "unipartite",
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		nscan = 600,
		burn = 300,
		odens = 5,
		keep_snap_draws = "draws",
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		seed = 124
	))
	for (fit in list(fit_mcmc, fit_second)) {
		expect_equal(unname(rank(-fit$snap_prob[, 3], ties.method = "min")[2]), 1)
		expect_equal(unname(rank(-fit$snap_prob_v[, 3], ties.method = "min")[5]), 1)
		expect_gt(fit$snap_prob[2, 3], .8)
		expect_gt(fit$snap_prob_v[5, 3], .8)
	}
})

test_that("MCMC snap draw summaries respect changing actor composition", {
	Y_full = simulate_snap_panel(n = 4, Tt = 3, snap = FALSE, seed = 99)
	Y = list(
		t1 = Y_full[[1]][paste0("a", 1:3), paste0("a", 1:3)],
		t2 = Y_full[[2]][paste0("a", 1:4), paste0("a", 1:4)],
		t3 = Y_full[[3]][paste0("a", 2:4), paste0("a", 2:4)]
	)
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		R = 1,
		family = "normal",
		mode = "unipartite",
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		nscan = 30,
		burn = 10,
		odens = 5,
		keep_snap_draws = "draws",
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE,
		seed = 321
	)))

	expect_true(all(is.na(fit$snap_prob[, "t1"])))
	expect_true(is.na(fit$snap_prob["a4", "t2"])) # entry transition
	expect_true(is.na(fit$snap_prob["a1", "t3"])) # absent at t3
	expect_true(is.finite(fit$snap_prob["a2", "t2"]))
	expect_equal(fit$snap_draws_summary$n["a4", "t2"], 0)
	expect_equal(fit$snap_draws_summary$n["a1", "t3"], 0)
	expect_true(all(is.na(fit$snap_draws[, "a4", "t2"])))
	expect_true(all(is.na(fit$snap_draws[, "a1", "t3"])))
})

test_that("an exhausted dynamic ALS budget is an error that reports its unfinished score", {
	set.seed(9110)
	Y = replicate(3, {
		m = matrix(rbinom(64, 1, 0.35), 8, 8)
		diag(m) = NA
		m
	}, simplify = FALSE)
	warnings = character()
	err = withCallingHandlers(
		expect_error(
			lame(
				Y,
				family = "binary",
				R = 0,
				method = "als",
				dynamic_ab = TRUE,
				als_max_iter = 1,
				bias_correct = "none",
				verbose = FALSE
			),
			class = "lame_als_convergence"
		),
		warning = function(w) {
			warnings <<- c(warnings, conditionMessage(w))
			invokeRestart("muffleWarning")
		}
	)
	# nothing precedes the error: the throwaway static pass that seeds the paths is silent (its stop is recorded)
	expect_length(warnings, 0L)
	expect_match(
		gsub("\\s+", " ", conditionMessage(err)),
		"did not reach its convergence criterion within `max_iter` = 1\\."
	)
	fit = err$fit
	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$converged)
	expect_false(fit$diagnostics$stationarity$pass)
	expect_gt(
		fit$diagnostics$stationarity$max_score,
		fit$diagnostics$stationarity$threshold
	)
	expect_false(fit$diagnostics$static_seed$converged)
	expect_identical(err$kind, "bcd")
	expect_false(isTRUE(err$retry))
})
