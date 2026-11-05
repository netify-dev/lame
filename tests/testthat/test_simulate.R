skip_on_cran()

# simulate.ame and simulate.lame

library(testthat)
library(lame)

test_that("simulate.ame works for binary networks", {
	# data and a fit
	data(YX_bin)
	set.seed(6886)

	# a simple binary ame model
	fit = ame(
		YX_bin$Y,
		YX_bin$X[, , 1:2],
		burn = 20,
		nscan = 30,
		odens = 10,
		family = "binary",
		verbose = FALSE
	)

	# basic simulation
	sims = simulate(
		fit,
		nsim = 5,
		seed = 456,
		newdata = list(Xdyad = YX_bin$X[, , 1:2])
	)

	expect_s3_class(sims, "ame.sim")
	expect_equal(length(sims$Y), 5)
	expect_equal(sims$family, "binary")
	expect_equal(sims$mode, "unipartite")

	# dimensions of the simulated networks
	expect_equal(dim(sims$Y[[1]]), dim(YX_bin$Y))

	# binary networks hold only 0, 1 and NA
	vals = unique(as.vector(sims$Y[[1]]))
	vals = vals[!is.na(vals)]
	expect_true(all(vals %in% c(0, 1)))

	# the diagonal is NA (unipartite)
	expect_true(all(is.na(diag(sims$Y[[1]]))))
})

test_that("simulate.ame captures uncertainty", {
	# data and a fit
	data(YX_bin)
	set.seed(6886)

	fit = ame(
		YX_bin$Y,
		YX_bin$X[, , 1:2],
		burn = 20,
		nscan = 30,
		odens = 10,
		family = "binary",
		verbose = FALSE
	)

	# simulate multiple networks
	sims = simulate(
		fit,
		nsim = 20,
		seed = 789,
		newdata = list(Xdyad = YX_bin$X[, , 1:2])
	)

	# calculate densities
	densities = sapply(sims$Y, function(y) mean(y == 1, na.rm = TRUE))

	# different draws give different networks
	expect_true(sd(densities) > 0)

	# the mean density is reasonable
	orig_density = mean(YX_bin$Y == 1, na.rm = TRUE)
	expect_true(abs(mean(densities) - orig_density) < 0.1)

	# latent Z matrices can be returned
	sims_z = simulate(
		fit,
		nsim = 3,
		return_latent = TRUE,
		seed = 111,
		newdata = list(Xdyad = YX_bin$X[, , 1:2])
	)
	expect_false(is.null(sims_z$Z))
	expect_equal(length(sims_z$Z), 3)
})

test_that("simulate.lame works for longitudinal binary networks", {
	# data
	data(YX_bin_list)
	set.seed(6886)

	# a subset for speed
	Y_sub = YX_bin_list$Y[1:3]
	X_sub = YX_bin_list$X[1:3]

	# a longitudinal fit (y is continuous; family = "binary" thresholds it)
	fit = suppressWarnings(lame(
		Y_sub,
		X_sub,
		burn = 20,
		nscan = 30,
		odens = 10,
		family = "binary",
		verbose = FALSE,
		plot = FALSE
	))

	# an extended horizon is the conditional approximation, which has to be asked for
	expect_error(simulate(fit, nsim = 5, n_time = 4, seed = 456), class = "lame_simulate_periods")
	expect_error(simulate(fit, nsim = 5, n_time = 4, seed = 456), "approximate")

	# basic simulation
	sims = simulate(
		fit,
		nsim = 5,
		n_time = 4,
		seed = 456,
		type = "approximate",
		newdata = list(Xdyad = X_sub)
	)

	expect_s3_class(sims, "lame.sim")
	expect_equal(length(sims$Y), 5) # 5 trajectories
	expect_equal(sims$n_time, 4) # 4 time points each
	expect_equal(sims$family, "binary")

	# each trajectory has the right structure
	expect_equal(length(sims$Y[[1]]), 4) # 4 time points

	# network dimensions at each time
	for (t in 1:4) {
		expect_equal(dim(sims$Y[[1]][[t]]), dim(Y_sub[[1]]))
	}

	# different starting points
	expect_error(
		simulate(fit, nsim = 3, n_time = 3, start_from = "random", seed = 222),
		class = "lame_simulate_periods"
	)
	sims_random = simulate(
		fit,
		nsim = 3,
		n_time = 3,
		type = "approximate",
		start_from = "random",
		seed = 222
	)
	expect_equal(length(sims_random$Y), 3)
})

test_that("simulate.lame captures temporal dynamics", {
	# data
	data(YX_bin_list)
	set.seed(6886)

	# a subset
	Y_sub = YX_bin_list$Y[1:3]
	X_sub = YX_bin_list$X[1:3]

	# a fit with dynamic effects
	fit = suppressWarnings(lame(
		Y_sub,
		X_sub,
		burn = 20,
		nscan = 30,
		odens = 10,
		family = "binary",
		verbose = FALSE,
		plot = FALSE
	))

	# a longer horizon
	sims = simulate(fit, nsim = 10, n_time = 10, seed = 333, type = "approximate")

	# the trajectories have the requested length
	expect_equal(length(sims$Y[[1]]), 10)

	# densities over time for each trajectory
	density_matrix = matrix(NA, 10, 10) # 10 sims x 10 time points

	for (i in 1:10) {
		for (t in 1:10) {
			density_matrix[i, t] = mean(sims$Y[[i]][[t]] == 1, na.rm = TRUE)
		}
	}

	# densities vary over time and across draws (temporal dynamics and uncertainty)
	time_variation = apply(density_matrix, 1, sd) # sd across time for each sim
	sim_variation = apply(density_matrix, 2, sd) # sd across sims for each time

	expect_true(mean(time_variation) > 0) # variation over time
	expect_true(mean(sim_variation) > 0) # variation across simulations

	# print and summary methods
	expect_no_error(print(sims))
	expect_no_error(summary(sims))
})

test_that("simulate.lame mirrors every draw of a symmetric fit", {
	set.seed(6886)
	n = 12
	Tn = 2
	nm = paste0("a", seq_len(n))
	u = rnorm(n)
	Y = lapply(seq_len(Tn), function(t) {
		z = 0.8 * outer(u, u) + matrix(rnorm(n * n), n, n)
		z[lower.tri(z)] = t(z)[lower.tri(z)]
		m = 1 * (z > 0)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(lame(
		Y,
		family = "binary",
		R = 1,
		symmetric = TRUE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	sims = simulate(fit, nsim = 3, n_time = Tn, seed = 1)
	expect_s3_class(sims, "lame.sim")
	for (traj in sims$Y) {
		for (m in traj) {
			m[is.na(m)] = 0
			expect_equal(m, t(m))
		}
	}
})

test_that("simulate.ame applies the symmetric L once and restores the effects' spread", {
	set.seed(6886)
	n = 24
	nm = sprintf("a%02d", seq_len(n))
	U = matrix(rnorm(n * 2, 0, 0.8), n, 2)
	prod = U %*% diag(c(1, -1)) %*% t(U)
	a = rnorm(n, 0, 0.5)
	E = matrix(0, n, n)
	E[upper.tri(E)] = rnorm(n * (n - 1) / 2)
	E = E + t(E)
	Y = matrix(rpois(n * n, exp(pmin(0.3 + outer(a, a, "+") + prod + 0.3 * E, 5))), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	dimnames(Y) = list(nm, nm)
	fit = suppressWarnings(ame(
		Y,
		family = "poisson",
		R = 2,
		symmetric = TRUE,
		nscan = 300,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	sims = simulate(fit, nsim = 40, seed = 1)
	ratio = mean(sapply(sims$Y, mean, na.rm = TRUE)) / mean(fit$YPM, na.rm = TRUE)
	# L is applied once, so the simulated mean stays near YPM
	expect_gt(ratio, 0.75)
	expect_lt(ratio, 1.35)
	# a sparse directed binary fit: the simulated mean tracks YPM
	set.seed(3)
	Yb = 1 * ((-1.7 + outer(rnorm(n, 0, 0.6), rnorm(n, 0, 0.6), "+") + matrix(rnorm(n * n), n, n)) > 0)
	diag(Yb) = NA
	dimnames(Yb) = list(nm, nm)
	fb = suppressWarnings(ame(
		Yb,
		family = "binary",
		R = 0,
		nscan = 400,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	sb = simulate(fb, nsim = 200, seed = 2)
	rb = mean(sapply(sb$Y, mean, na.rm = TRUE)) / mean(fb$YPM, na.rm = TRUE)
	expect_gt(rb, 0.9)
	expect_lt(rb, 1.12)
})

test_that("simulate() on a fit without saved draws stops unless the approximation is requested", {
	set.seed(21)
	n = 10
	nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	fit = suppressWarnings(ame(
		Y,
		family = "normal",
		R = 0,
		burn = 5,
		nscan = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		posterior_opts = posterior_options()
	))
	for (ty in c("posterior", "auto")) {
		e = expect_error(simulate(fit, nsim = 2, type = ty), class = "lame_simulate_draws")
		expect_match(conditionMessage(e), "approximate")
		expect_match(conditionMessage(e), "posterior_options")
	}
	expect_error(simulate(fit, nsim = 2), class = "lame_simulate_draws")
	expect_error(simulate(fit, nsim = 2), class = "lame_predictive_draws")
	sims = simulate(fit, nsim = 2, type = "approximate", seed = 1)
	expect_identical(sims$method, "conditional approximation")
	expect_length(sims$Y, 2L)
	# a fit that kept its draws is simulated from them by default, and "auto" is the same route
	fit2 = suppressWarnings(ame(
		Y,
		family = "normal",
		R = 0,
		burn = 5,
		nscan = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE
	))
	a = simulate(fit2, nsim = 3, seed = 4)
	expect_identical(a$method, "joint posterior")
	expect_identical(simulate(fit2, nsim = 3, seed = 4, type = "auto")$Y, a$Y)
	# a panel: an extended n_time and start_from belong to the approximation
	Yl = list(Y, Y + rnorm(n * n, 0, .1)); Yl = lapply(Yl, function(y) { diag(y) = NA; y })
	fl = suppressWarnings(lame(
		Yl,
		family = "normal",
		R = 0,
		burn = 5,
		nscan = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_error(simulate(fl, nsim = 2, n_time = 3), class = "lame_simulate_periods")
	expect_error(simulate(fl, nsim = 2, start_from = "random"), class = "lame_simulate_periods")
	ok = simulate(fl, nsim = 2, n_time = 3, type = "approximate", seed = 2)
	expect_length(ok$Y[[1]], 3L)
})

test_that("the approximation follows a dynamic_beta intercept per period and stops past the fitted periods", {
	set.seed(31)
	n = 10
	Tn = 4
	nm = paste0("a", seq_len(n))
	# the intercept rises by one per period
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(t + rnorm(n * n, 0, 0.5), n, n, dimnames = list(nm, nm))
		diag(y) = NA
		y
	})
	names(Yl) = paste0("t", seq_len(Tn))
	fit = suppressWarnings(lame(
		Yl,
		family = "normal",
		R = 0,
		dynamic_beta = "intercept",
		burn = 20,
		nscan = 60,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		chains = 1
	))
	expect_length(dim(fit$BETA), 3L)
	# a request beyond the fitted periods has no coefficient path to follow
	expect_error(
		simulate(fit, nsim = 2, type = "approximate", n_time = Tn + 2, seed = 1),
		class = "lame_simulate_periods"
	)
	# a two-dimensional covariate slice falls back to the intercept, which is the period's own
	nd = list(Xdyad = lapply(seq_len(Tn), function(t) matrix(0, n, n)))
	sims = simulate(fit, nsim = 20, type = "approximate", newdata = nd, seed = 2)
	period_mean = sapply(seq_len(Tn), function(t) {
		mean(sapply(sims$Y, function(y) mean(y[[t]], na.rm = TRUE)))
	})
	expect_gt(period_mean[Tn] - period_mean[1], 2)
})

test_that("latent draws of a square bipartite fit carry no reciprocity", {
	set.seed(32)
	n = 10
	Y = matrix(
		rnorm(n * n),
		n,
		n,
		dimnames = list(paste0("r", seq_len(n)), paste0("c", seq_len(n)))
	)
	fit = suppressWarnings(ame(
		Y,
		family = "normal",
		mode = "bipartite",
		R = 0,
		burn = 10,
		nscan = 40,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	))
	expect_true(all(fit$VC[, "rho"] == 0))
	# a nonzero rho column must not couple cell (i, j) with cell (j, i) of two different actor sets
	fit$VC[, "rho"] = 0.9
	sims = simulate(fit, nsim = 200, type = "approximate", return_latent = TRUE, seed = 3)
	pair_cor = mean(sapply(sims$Z, function(z) {
		z = z - mean(z)
		cor(z[upper.tri(z)], t(z)[upper.tri(z)])
	}))
	expect_lt(abs(pair_cor), 0.15)
	# the panel simulator draws its latent cells the same way
	Yl = list(t1 = Y, t2 = Y + rnorm(n * n, 0, 0.1))
	fl = suppressWarnings(lame(
		Yl,
		family = "normal",
		mode = "bipartite",
		R = 0,
		burn = 10,
		nscan = 40,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		chains = 1
	))
	fl$VC[, "rho"] = 0.9
	sl = simulate(fl, nsim = 100, type = "approximate", return_latent = TRUE, seed = 4)
	pair_cor_l = mean(sapply(sl$Z, function(zt) {
		z = zt[[1]] - mean(zt[[1]])
		cor(z[upper.tri(z)], t(z)[upper.tri(z)])
	}))
	expect_lt(abs(pair_cor_l), 0.15)
})

test_that("approximate simulation builds newdata designs with the intercept and covariates in coefficient order", {
	skip_on_cran()
	set.seed(2)
	n = 15
	nm = sprintf("a%02d", 1:n)
	X = lapply(1:3, function(t) array(rnorm(n * n), c(n, n, 1), dimnames = list(nm, nm, "x")))
	Y = lapply(1:3, function(t) {
		m = 1 + 2 * X[[t]][, , 1] + matrix(rnorm(n * n), n, n, dimnames = list(nm, nm))
		diag(m) = NA
		m
	})
	f = suppressMessages(suppressWarnings(lame(
		Y,
		Xdyad = X,
		family = "normal",
		nscan = 300,
		burn = 100,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	sim_mean = function(s) mean(sapply(s$Y, function(yl) mean(unlist(yl), na.rm = TRUE)))
	s_fit = simulate(f, nsim = 40, seed = 1, type = "approximate")
	s_new = simulate(f, nsim = 40, seed = 1, type = "approximate", newdata = list(Xdyad = X))
	s_zero = simulate(
		f,
		nsim = 40,
		seed = 1,
		type = "approximate",
		newdata = list(Xdyad = lapply(X, function(a) a * 0))
	)
	expect_equal(sim_mean(s_new), sim_mean(s_fit), tolerance = 0.05)
	expect_gt(cor(as.vector(s_new$Y[[1]][[1]]), as.vector(X[[1]][, , 1]), use = "complete.obs"), 0.8)
	# with the covariate at zero the intercept is still in the predictor
	expect_gt(sim_mean(s_zero), 0.5)
})
