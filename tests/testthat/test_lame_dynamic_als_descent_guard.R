skip_on_cran()

library(lame)

# the penalized dynamic ALS route (bias_correct = "none"): a sweep of exact block
# minimizations cannot raise the objective, so a block update that does is a
# solver fault. the descent guard undoes the sweep and redoes it with every block
# checked: a block that raises the objective is stepped back along its own
# displacement, or kept at its start. the faults below are injected into the
# block solvers

gd_gaussian = function(n = 8L, Tt = 5L, seed = 11) {
	set.seed(seed)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		a = seq(-1, 1, length.out = n) * sin(t / 2)
		b = rev(a) * cos(t / 3)
		Yt = 0.25 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	Y
}

gd_bilinear = function(n = 9L, Tt = 5L, seed = 5) {
	set.seed(seed)
	U0 = matrix(rnorm(n), n, 1)
	V0 = matrix(rnorm(n), n, 1)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		Yt = 0.2 + (0.6 + 0.1 * t) * tcrossprod(U0, V0) +
			matrix(rnorm(n * n, 0, 0.3), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	Y
}

gd_counts = function(n = 10L, Tt = 4L, seed = 21) {
	set.seed(seed)
	Y = vector("list", Tt)
	for (t in seq_len(Tt)) {
		a = seq(-0.6, 0.6, length.out = n) * cos(t / 2)
		b = rev(a) * 0.5
		Yt = matrix(rpois(n * n, exp(-0.2 + outer(a, b, "+"))), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	Y
}

gd_bipartite = function(nr = 8L, nc = 10L, Tt = 4L, seed = 4) {
	set.seed(seed)
	Y = vector("list", Tt)
	u = rnorm(nr)
	v = rnorm(nc)
	for (t in seq_len(Tt)) {
		a = seq(-1, 1, length.out = nr) * sin(t / 2)
		b = seq(-1, 1, length.out = nc) * cos(t / 3)
		Yt = 0.25 + outer(a, b, "+") + 0.8 * tcrossprod(u, v) +
			matrix(rnorm(nr * nc, 0, 0.3), nr, nc)
		dimnames(Yt) = list(paste0("r", seq_len(nr)), paste0("c", seq_len(nc)))
		Y[[t]] = Yt
	}
	Y
}

gd_fit = function(Y, ...) {
	lame(Y, bias_correct = "none", method = "als", verbose = FALSE, ...)
}

gd_non_increasing = function(fit) {
	tr = fit$diagnostics$objective
	all(diff(tr) <= 1e-10 * (abs(utils::head(tr, -1L)) + 1))
}

# the additive block's update pushed k times past its own displacement: beyond
# twice the step it overshoots the minimum and raises the objective
gd_overshoot_ab = function(k) {
	orig = get(".lda_update_dynamic_ab", envir = asNamespace("lame"))
	function(Z, W, Xb, Oarr, a_mat, b_mat, ...) {
		r = orig(Z, W, Xb, Oarr, a_mat, b_mat, ...)
		r$a = a_mat + k * (r$a - a_mat)
		r$b = b_mat + k * (r$b - b_mat)
		r
	}
}

test_that("the descent guard does not touch a fit whose blocks descend", {
	fit = expect_no_warning(gd_fit(
		gd_gaussian(),
		family = "normal",
		R = 0,
		dynamic_ab = TRUE
	))
	g = fit$diagnostics$descent_guard
	expect_false(g$active)
	expect_identical(c(g$sweeps_redone, g$blocks_halved, g$blocks_rolled_back), c(0L, 0L, 0L))
	expect_false(g$stalled)
	expect_true(fit$converged)
	expect_identical(fit$diagnostics$objective_increases, 0L)
	expect_equal(nrow(fit$diagnostics$convergence_trace), fit$iterations)
})

test_that("a panel with actors entering and leaving does not set the guard off: the first sweep puts the start on the gauge", {
	set.seed(17)
	sets = list(t1 = paste0("a", 1:4), t2 = paste0("a", 2:5), t3 = paste0("a", c(1, 3, 5)))
	Y = lapply(seq_along(sets), function(tt) {
		a = sets[[tt]]
		y = matrix(rnorm(length(a)^2), length(a), length(a), dimnames = list(a, a))
		diag(y) = NA_real_
		y
	})
	names(Y) = names(sets)
	# the paths enter centred over all actors and the block solves hold the sum to zero over the actors present; the sweep
	# that moves the start onto that gauge raises the objective with every block exact, which is not a solver fault
	fit = suppressWarnings(gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE, als_max_iter = 30))
	expect_false(fit$diagnostics$descent_guard$active)
	expect_identical(fit$diagnostics$descent_guard$sweeps_redone, 0L)
	tr = fit$diagnostics$objective
	expect_gt(tr[1L], 21)
	expect_true(all(diff(tr) <= 1e-10 * (abs(utils::head(tr, -1L)) + 1)))
	expect_identical(fit$diagnostics$objective_increases, 0L)
	# the score of the additive paths is projected onto the same zero-sum gauge, so a panel with actors entering and leaving converges
	expect_true(fit$converged)
	expect_lt(fit$iterations, 10L)
})

test_that("an additive block that overshoots is stepped back and the fit converges", {
	Y = gd_gaussian()
	ref = gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE)
	testthat::local_mocked_bindings(
		.lda_update_dynamic_ab = gd_overshoot_ab(3),
		.package = "lame"
	)
	# no warning of any kind: a fit that needed the guard is a converged fit
	fit = expect_no_warning(gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE))
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_identical(g$sweeps_redone, 1L)
	expect_gt(g$blocks_halved, 0L)
	expect_identical(g$blocks_rolled_back, 0L)
	expect_identical(g$blocks, "ab")
	expect_false(g$stalled)
	expect_true(fit$converged)
	expect_true(gd_non_increasing(fit))
	expect_identical(fit$diagnostics$objective_increases, 0L)
	# the solution is the unmocked one, to the tolerance
	expect_lt(max(abs(fit$a_dynamic - ref$a_dynamic)), 5e-3)
	expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 5e-3)
	expect_equal(fit$deviance, ref$deviance, tolerance = 1e-4)
	# the monitored sequence is the accepted sweeps only
	expect_lt(nrow(fit$diagnostics$convergence_trace), fit$iterations)
})

test_that("a block that cannot descend is kept at its start and the fit is refused with its state", {
	Y = gd_gaussian()
	orig = get(".lda_update_dynamic_ab", envir = asNamespace("lame"))
	# the displacement reversed: every point of the segment is worse than its start
	testthat::local_mocked_bindings(
		.lda_update_dynamic_ab = function(Z, W, Xb, Oarr, a_mat, b_mat, ...) {
			r = orig(Z, W, Xb, Oarr, a_mat, b_mat, ...)
			r$a = a_mat - (r$a - a_mat)
			r$b = b_mat - (r$b - b_mat)
			r
		},
		.package = "lame"
	)
	w = NULL
	err = withCallingHandlers(
		expect_error(
			gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE),
			class = "lame_als_convergence"
		),
		warning = function(cnd) {
			w <<- c(w, list(cnd))
			invokeRestart("muffleWarning")
		}
	)
	# one condition, an error: no warning of any kind precedes it
	expect_length(w, 0L)
	msg = gsub("\\s+", " ", conditionMessage(err))
	expect_match(msg, "stopped without converging: a block update could not lower the penalized objective")
	expect_match(msg, "ab block update raised the objective at every step the descent guard tried")
	expect_false(grepl("caution", msg))
	fit = err$fit
	expect_s3_class(fit, "lame_dynamic_als")
	expect_false(fit$converged)
	g = fit$diagnostics$descent_guard
	expect_true(g$stalled)
	expect_gt(g$blocks_rolled_back, 0L)
	expect_identical(g$blocks_halved, 0L)
	expect_identical(g$blocks, "ab")
	expect_true(gd_non_increasing(fit))
	expect_identical(fit$diagnostics$objective_increases, 0L)
	expect_match(fit$diagnostics$message, "could not lower the penalized objective")
	expect_true(all(is.finite(fit$a_dynamic)))
	# the state is the last one that lowered the objective: the score at it is reported
	expect_false(fit$diagnostics$stationarity$pass)
})

test_that("a block that returns non-finite values is kept at its start", {
	skip_on_cran()
	Y = gd_gaussian()
	orig = get(".lda_update_beta_path", envir = asNamespace("lame"))
	testthat::local_mocked_bindings(
		.lda_update_beta_path = function(Z, W, prep, coef_path_old, ...) {
			r = orig(Z, W, prep, coef_path_old, ...)
			r[] = NA_real_
			r
		},
		.package = "lame"
	)
	fit = suppressWarnings(gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE))
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_gt(g$blocks_rolled_back, 0L)
	expect_true("beta" %in% g$blocks)
	expect_true(all(is.finite(coef(fit))))
	expect_true(all(is.finite(unlist(fitted(fit)[[1L]])[!is.na(Y[[1L]])])))
	expect_true(gd_non_increasing(fit))
})

test_that("a dynamic factor block that overshoots is stepped back", {
	skip_on_cran()
	Y = gd_bilinear()
	ref = gd_fit(Y, family = "normal", R = 1, dynamic_uv = TRUE)
	expect_false(ref$diagnostics$descent_guard$active)
	orig = get(".lda_update_dynamic_uv_directed", envir = asNamespace("lame"))
	to_o = get(".lda_uv_to_o_array", envir = asNamespace("lame"))
	testthat::local_mocked_bindings(
		.lda_update_dynamic_uv_directed = function(
			Z,
			W,
			Xb,
			a_mat,
			b_mat,
			U_cube,
			V_cube,
			...
		) {
			r = orig(Z, W, Xb, a_mat, b_mat, U_cube, V_cube, ...)
			r$U = U_cube + 3 * (r$U - U_cube)
			r$V = V_cube + 3 * (r$V - V_cube)
			r$Oarr = to_o(r$U, r$V)
			r
		},
		.package = "lame"
	)
	fit = expect_no_warning(gd_fit(Y, family = "normal", R = 1, dynamic_uv = TRUE))
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_gt(g$blocks_halved, 0L)
	expect_identical(g$blocks, "mult")
	expect_true(fit$converged)
	expect_true(gd_non_increasing(fit))
	expect_equal(fit$deviance, ref$deviance, tolerance = 1e-4)
	expect_identical(fit$diagnostics$objective_increases, 0L)
})

test_that("the static low-rank block, which cannot be stepped back, is rolled back", {
	skip_on_cran()
	Y = gd_bilinear()
	orig = get(".lda_update_static_lowrank", envir = asNamespace("lame"))
	testthat::local_mocked_bindings(
		.lda_update_static_lowrank = function(...) {
			r = orig(...)
			r$Omat = -r$Omat
			r
		},
		.package = "lame"
	)
	err = expect_error(
		gd_fit(Y, family = "normal", R = 1, dynamic_ab = TRUE),
		class = "lame_als_convergence"
	)
	fit = err$fit
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_identical(g$blocks_halved, 0L)
	expect_gt(g$blocks_rolled_back, 0L)
	expect_identical(g$blocks, "mult")
	expect_true(gd_non_increasing(fit))
	expect_false(fit$converged)
})

test_that("the guard also covers the working problems of the IRLS families", {
	skip_on_cran()
	Y = gd_counts()
	ref = gd_fit(Y, family = "poisson", R = 0, dynamic_ab = TRUE)
	expect_false(ref$diagnostics$descent_guard$active)
	testthat::local_mocked_bindings(
		.lda_update_dynamic_ab = gd_overshoot_ab(3),
		.package = "lame"
	)
	fit = expect_no_warning(gd_fit(Y, family = "poisson", R = 0, dynamic_ab = TRUE))
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_gt(g$blocks_halved, 0L)
	expect_true(fit$converged)
	expect_identical(fit$diagnostics$objective_increases, 0L)
	expect_equal(fit$deviance, ref$deviance, tolerance = 1e-5)
	expect_lt(max(abs(fit$a_dynamic - ref$a_dynamic)), 5e-3)
})

test_that("a stall that the certified joint refinement of a bipartite panel resolves is a converged fit, not a refusal", {
	Y = gd_bipartite()
	orig = get(".lda_update_dynamic_uv_bip_g", envir = asNamespace("lame"))
	to_o = get(".lda_uvg_to_o_array", envir = asNamespace("lame"))
	# the displacement of the factor block reversed: the block cannot descend, the sweeps stall, and the joint refinement
	# that follows a stalled normal bipartite fit certifies the state
	testthat::local_mocked_bindings(
		.lda_update_dynamic_uv_bip_g = function(Z, W, Xb, a_mat, b_mat, U_cube, V_cube, G, ...) {
			r = orig(Z, W, Xb, a_mat, b_mat, U_cube, V_cube, G, ...)
			r$U = U_cube - (r$U - U_cube)
			r$V = V_cube - (r$V - V_cube)
			r$Oarr = to_o(r$U, r$V, r$G)
			r
		},
		.package = "lame"
	)
	fit = expect_no_warning(expect_no_error(
		gd_fit(Y, family = "normal", mode = "bipartite", R = 1, dynamic_uv = TRUE)
	))
	expect_true(fit$converged)
	expect_true(fit$diagnostics$joint_refinement$accepted)
	expect_true(fit$diagnostics$stationarity$pass)
	g = fit$diagnostics$descent_guard
	expect_true(g$active)
	expect_false(g$stalled)
	expect_gt(g$blocks_rolled_back, 0L)
	expect_identical(g$blocks, "mult")
})

test_that("the tolerance of the guard follows the unit of the response", {
	orig = get(".lda_update_dynamic_ab", envir = asNamespace("lame"))
	testthat::local_mocked_bindings(
		.lda_update_dynamic_ab = gd_overshoot_ab(2.6),
		.package = "lame"
	)
	for (unit in c(1e-6, 1, 1e3)) {
		Y = lapply(gd_gaussian(), function(y) unit * y)
		fit = suppressWarnings(gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE))
		g = fit$diagnostics$descent_guard
		info = paste("unit", unit)
		# the overshoot is caught at every unit of the response, a tiny objective included
		expect_true(g$active, info = info)
		expect_gt(g$blocks_halved + g$blocks_rolled_back, 0L, label = info)
		tr = fit$diagnostics$objective
		expect_true(all(diff(tr) <= 2e-10 * abs(utils::head(tr, -1L))), info = info)
	}
})

test_that("the summary reports a guarded fit and says nothing of an unguarded one", {
	skip_on_cran()
	Y = gd_gaussian()
	plain = gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE)
	out = utils::capture.output(print(summary(plain)))
	expect_false(any(grepl("Descent guard", out)))
	testthat::local_mocked_bindings(
		.lda_update_dynamic_ab = gd_overshoot_ab(3),
		.package = "lame"
	)
	fit = gd_fit(Y, family = "normal", R = 0, dynamic_ab = TRUE)
	out = utils::capture.output(print(summary(fit)))
	expect_true(any(grepl("Descent guard from iteration", out)))
	expect_true(any(grepl("stepped back", out)))
})
