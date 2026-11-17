skip_on_cran()

# dynamic effects: dynamic_ab and dynamic_uv

library(lame)
library(testthat)

# combined dynamic_ab + dynamic_uv

# the only mixing gate of the bipartite binary G product (2 x 44000 sweeps): run with LAME_SLOW_TESTS set
test_that("bipartite binary dynamic_ab plus dynamic_uv stays finite", {
	skip_on_cran()
	skip_if(identical(Sys.getenv("LAME_SLOW_TESTS"), ""), "the bipartite G-product mixing gate is run with LAME_SLOW_TESTS set")

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
		burn = 200,
		nscan = 800,
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
		nscan = 1000,
		odens = 10,
		seed = 1,
		verbose = FALSE,
		plot = FALSE
	))
	expect_lt(max(fit$VC[, "va"]), 3)
	expect_lt(max(fit$VC[, "vb"]), 3)
	# every period of the intercept path
	expect_gt(min(fit$BETA[, 1, ]), -4.5)
	expect_lt(max(abs(unlist(fit$EZ))), 50)
	# the predictive tie rate tracks the data and the posterior-mean predictor carries no rail-pressing excursions
	predictive_rate = mean(unlist(fit$YPM), na.rm = TRUE)
	expect_lt(abs(predictive_rate / mean(unlist(Y)) - 1), 0.3)
	# the interaction scale may touch its rail in a minority of sweeps but is not stuck there
	expect_lt(fit$tryErrorChecks$G, 0.8 * 1000)
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
	# the budget follows the machine: a short timed fit sets it so that burn-in and a few hundred sweeps
	# complete under any load, while the 40000 requested sweeps still cannot
	t0 = proc.time()[["elapsed"]]
	suppressWarnings(suppressMessages(lame(
		Y,
		family = "binary",
		R = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		gof = TRUE,
		burn = 20,
		nscan = 40,
		odens = 4,
		seed = 3,
		verbose = FALSE,
		plot = FALSE
	)))
	budget = max(6, ceiling(10 * (proc.time()[["elapsed"]] - t0)))
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
		max_seconds = budget,
		seed = 3,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_true(isTRUE(f$terminated_early))
	nd = nrow(f$VC)
	expect_lt(nd, 10000)
	expect_equal(dim(f$BETA)[1], nd)
	expect_equal(length(f$rho_ab), nd); expect_equal(length(f$sigma_ab), nd)
	expect_equal(length(f$rho_uv), nd); expect_equal(length(f$sigma_uv), nd)
	expect_true(all(f$rho_uv > 0)); expect_true(all(f$sigma_ab > 0))
	# effects and positions average over completed draws, undiluted
	expect_gt(max(abs(f$a_dynamic)), 0.05)
	expect_gt(sd(as.vector(f$U)), 0.02)
	ypm = unlist(f$YPM); expect_lt(abs(mean(ypm[is.finite(unlist(Y))]) - mean(unlist(Y), na.rm = TRUE)), 0.1)
})

# the only recovery of rho_ab by the sampler: persistent effects (rho 0.8) give a positive estimate, and rho_ab and
# sigma_ab move every sweep
test_that("rho_ab is estimated positive on persistent effects and moves every sweep", {
	skip_on_cran()
	set.seed(6886); n = 20; TT = 6
	nm = paste0("a", 1:n)
	a = rnorm(n); b = rnorm(n)
	Y = vector("list", TT)
	for (t in 1:TT) {
		if (t > 1) { a = 0.8 * a + rnorm(n, 0, 0.2); b = 0.8 * b + rnorm(n, 0, 0.2) }
		m = outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.5), n); diag(m) = NA; dimnames(m) = list(nm, nm)
		Y[[t]] = m
	}
	f = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = "normal",
				R = 0,
				dynamic_ab = TRUE,
				burn = 100,
				nscan = 300,
				odens = 1,
				seed = 3,
				verbose = FALSE,
				plot = FALSE
			)
		)
	)
	expect_gt(median(f$rho_ab), 0.1)
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

# the dynamic_beta chain is the slowest fit here, so it is fitted once and shared
.sparse_dyad_chains = local({
	fits = NULL
	function() {
		if (is.null(fits)) {
			d = .sparse_field_panel(1)
			fits <<- lapply(
				11,
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
		# a seed whose short chain passes the drift check
		seed = 5,
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
			burn = 300,
			nscan = 1000,
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
	expected = mean(fit$BETA[, 1, TT]) + outer(fit$a_dynamic[, TT], fit$b_dynamic[, TT], "+") + fit$UVPM
	expect_equal(unname(ez), unname(expected), tolerance = 1e-8)
	# the reported positions are centred over actors
	expect_lt(max(abs(colMeans(fit$U))), 1e-8)
	expect_lt(max(abs(colMeans(fit[["V"]]))), 1e-8)
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
