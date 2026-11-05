skip_on_cran()

forecast_joint_panel = function(bip = FALSE, symmetric = FALSE) {
	set.seed(205)
	nr = 8L; nc = if (bip) 7L else nr
	rn = paste0("a", seq_len(nr)); cn = if (bip) paste0("b", seq_len(nc)) else rn
	Y = lapply(1:3, function(tt) {
		m = matrix(rnorm(nr * nc), nr, nc, dimnames = list(rn, cn))
		if (symmetric) m = (m + t(m)) / 2
		if (!bip) diag(m) = NA
		m
	})
	list(
		Y = Y,
		family = "normal",
		mode = if (bip) "bipartite" else "unipartite",
		symmetric = symmetric,
		R = 1L,
		burn = 20L,
		nscan = 60L,
		odens = 2L,
		seed = 205L,
		gof = FALSE,
		verbose = FALSE,
		plot = FALSE,
		posterior_opts = posterior_options()
	)
}

test_that("terminal posterior states preserve uncertainty and reproduce fitted predictors", {
	skip_on_cran()
	for (mode in c("directed", "bipartite", "symmetric")) for (ab in c(FALSE, TRUE)) for (uv in c(FALSE, TRUE)) {
		a = forecast_joint_panel(mode == "bipartite", mode == "symmetric")
		a$dynamic_ab = ab; a$dynamic_uv = uv
		f = suppressWarnings(suppressMessages(do.call(lame, a)))
		state = f$forecast_state
		expect_identical(state$draw_index, 1:30)
		expect_null(f$a_draws)
		expect_null(f$U_draws)
		if (ab) { f$rho_ab[] = 1; f$sigma_ab[] = 0 }
		if (uv) { f$rho_uv[] = 1; f$sigma_uv[] = 0 }
		fc = suppressWarnings(predict(f, h = 1, type = "link", by_draw = TRUE, seed = 31))
		set.seed(31); ord = sample.int(30, 30)
		for (j in c(1L, 10L, 30L)) {
			i = ord[j]; U = state$U[, , i, drop = FALSE]; dim(U) = dim(state$U)[1:2]
			V = state$V[, , i, drop = FALSE]; dim(V) = dim(state$V)[1:2]
			prod = if (mode == "bipartite") U %*% state$G[, , i] %*% t(V) else U %*% t(V)
			expected = f$BETA[i, 1] + outer(state$a[, i], state$b[, i], "+") + prod
			expect_equal(unname(fc[, , 1, j]), unname(expected), tolerance = 1e-10)
		}
		average = apply(fc[, , 1, , drop = FALSE], c(1, 2), mean)
		target = f$EZ[[3]]; keep = is.finite(target)
		expect_equal(unname(average[keep]), unname(target[keep]), tolerance = 1e-10)
		if (mode == "symmetric" && ab) {
			f$sigma_ab[] = .2
			fresh = predict(f, h = 1, by_draw = TRUE, seed = 47)
			expect_equal(fresh[, , 1, 1], t(fresh[, , 1, 1]), tolerance = 1e-12)
		}
	}
})

test_that("forecasts retain cancellation between correlated posterior blocks", {
	beta = c(-3, -1, 1, 3); a = matrix(rep(-beta / 2, each = 2), 2)
	f = structure(
		list(
			BETA = matrix(beta, 4, 1, dimnames = list(NULL, "intercept")),
			VC = matrix(1, 4, 1, dimnames = list(NULL, "ve")),
			family = "normal",
			R = 0,
			APM = c(a = 0, b = 0),
			BPM = c(a = 0, b = 0),
			Xlist = list(array(1, c(2, 2, 1))),
			forecast_state = list(a = a, b = a, draw_index = 1:4)
		),
		class = "lame"
	)
	fc = suppressWarnings(predict(f, h = 1, type = "link", by_draw = TRUE, seed = 1))
	expect_equal(as.numeric(fc), rep(0, 16), tolerance = 0)
	f$forecast_state = NULL
	expect_error(suppressWarnings(predict(f, h = 1)), class = "lame_forecast_draws")
})

test_that("legacy saved trajectories supply matched terminal draws", {
	a = forecast_joint_panel(TRUE); a$dynamic_ab = TRUE; a$dynamic_uv = TRUE
	a$posterior_opts = list(save_ab = TRUE, thin_ab = 1, save_UV = TRUE)
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	legacy = f; legacy$forecast_state = NULL
	expect_equal(
		predict(f, h = 2, by_draw = TRUE, seed = 19),
		predict(legacy, h = 2, by_draw = TRUE, seed = 19),
		tolerance = 0
	)
	a$posterior_opts$thin_ab = 3
	thin = suppressWarnings(suppressMessages(do.call(lame, a))); thin$forecast_state = NULL
	state = lame:::.lame_forecast_state(thin)
	expect_identical(as.integer(state$draw_index), seq(3L, 30L, 3L))
	expect_equal(unname(state$a), unname(f$forecast_state$a[, state$draw_index]), tolerance = 0)
	expect_equal(unname(state$U), unname(f$forecast_state$U[, , state$draw_index, drop = FALSE]), tolerance = 0)
})

test_that("snap and heavy-tailed forecasts follow their fitted innovation laws", {
	set.seed(322)
	n = 100000L; z = matrix(2, n, 2)
	x = lame:::.lame_forecast_uv_step(z, .7, .3, "snap", pi = .25, kappa = 1.2)
	mu = .75 * .7 * 2
	variance = .75 * .3^2 + .25 * 1.2^2 + .25 * .75 * (.7 * 2)^2
	expect_equal(mean(x), mu, tolerance = .01)
	expect_equal(var(x[, 1]), variance, tolerance = .02)
	x = lame:::.lame_forecast_uv_step(z, .7, .3, "t", nu = 6) - 1.4
	expect_equal(mean(x), 0, tolerance = .005)
	expect_equal(var(x[, 1]), .3^2 * 6 / 4, tolerance = .006)
	expect_gt(cor(x[, 1]^2, x[, 2]^2), .1)
})

test_that("pooled dynamic chains keep joint states and stored product means", {
	a = forecast_joint_panel(TRUE); a$dynamic_ab = TRUE; a$dynamic_uv = TRUE
	a$posterior_opts = list(save_ab = TRUE, thin_ab = 1, save_UV = TRUE)
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	a$seed = 206; a$nscan = 40
	g = suppressWarnings(suppressMessages(do.call(lame, a)))
	f$log_lik = matrix(seq_len(60), 30, 2); g$log_lik = matrix(seq_len(40), 20, 2)
	p = combine_ame_chains(list(f, g), diagnostics = FALSE)
	expect_equal(p$UVPM, (30 * f$UVPM + 20 * g$UVPM) / 50, tolerance = 1e-12)
	expect_equal(p$a_dynamic, (30 * f$a_dynamic + 20 * g$a_dynamic) / 50, tolerance = 1e-12)
	expect_identical(p$rho_ab, c(f$rho_ab, g$rho_ab))
	expect_equal(p$ab_draw_index, 1:50)
	expect_equal(p$forecast_state$draw_index, 1:50)
	expect_equal(unname(p$forecast_state$a[, 31:50]), unname(g$forecast_state$a), tolerance = 0)
	expect_equal(dim(p$U_draws)[4], 50)
	expect_equal(dim(p$a_draws)[3], 50)
	expect_equal(dim(p$G_samples)[3], 50)
	expect_equal(ncol(p$GOF[[1]]), 51)
	expect_equal(unname(p$log_lik), unname(rbind(f$log_lik, g$log_lik)), tolerance = 0)
	expected_sd = sqrt(
		(30 * (f$a_dynamic_sd^2 + f$a_dynamic^2) +
				20 * (g$a_dynamic_sd^2 + g$a_dynamic^2)) / 50 - p$a_dynamic^2
	)
	expect_equal(p$a_dynamic_sd, expected_sd, tolerance = 1e-12)
	p$rho_ab[] = 1; p$sigma_ab[] = 0; p$rho_uv[] = 1; p$sigma_uv[] = 0
	fc = predict(p, h = 1, type = "link", seed = 6)[[1]]
	expect_equal(unname(fc), unname(p$EZ[[3]]), tolerance = 1e-10)
})


test_that("forecast covariates respect actor order and merged symmetric nodal terms", {
	X = array(
		0,
		c(2, 2, 3),
		dimnames = list(
			c("a", "b"),
			c("a", "b"),
			c("intercept", "size_row", "size_col")
		)
	)
	X[, , 1] = 1; X[, , 2] = c(1, 3); X[, , 3] = t(X[, , 2])
	Z = lame:::.lame_forecast_design(
		X[2:1, 2:1, ],
		c("intercept", "size_node"),
		c("a", "b"),
		c("a", "b"),
		2,
		2
	)
	expect_equal(unname(Z[, , 2]), outer(c(1, 3), c(1, 3), "+"))
	expect_equal(unname(Z[, , 1]), matrix(1, 2, 2))
})

test_that("centered actor forecasts follow the constrained Gaussian transition", {
	set.seed(7123)
	N = 4000L
	initial = c(1, -1, 0)
	rho = .6
	sigma = .7
	draws = t(replicate(
		N,
		lame:::.propagate_ar1_vec(
			initial,
			rho,
			sigma,
			3L,
			center = TRUE
		)[3L, ]
	))
	variance = sigma^2 * sum(rho^(2 * (0:2)))
	V = variance * (diag(3) - matrix(1 / 3, 3, 3))
	mu = rho^3 * initial
	expect_lt(max(abs(rowSums(draws))), 1e-12)
	expect_lt(max(abs(colMeans(draws) - mu) / sqrt(diag(V) / N)), 6)
	expect_lt(max(abs(cov(draws) - V) / sqrt((V^2 + outer(diag(V), diag(V))) / (N - 1))), 6)
})

test_that("forecasting keeps disabled actor margins absent", {
	S = 2000L
	f = structure(
		list(
			BETA = matrix(0, S, 1, dimnames = list(NULL, "intercept")),
			VC = matrix(1, S, 1, dimnames = list(NULL, "ve")),
			family = "normal",
			mode = "bipartite",
			R = 0,
			intercept = TRUE,
			dynamic_ab = TRUE,
			rho_ab = rep(.6, S),
			sigma_ab = rep(.5, S),
			APM = setNames(numeric(3), paste0("a", 1:3)),
			BPM = setNames(numeric(4), paste0("b", 1:4)),
			Xlist = list(array(1, c(3, 4, 1))),
			forecast_state = list(
				a = matrix(0, 3, S),
				b = matrix(0, 4, S),
				a_level = matrix(0, 3, S),
				b_level = matrix(0, 4, S),
				draw_index = seq_len(S)
			)
		),
		class = "lame"
	)
	for (margin in c("sender", "receiver", "neither")) {
		f$rvar = margin == "sender"
		f$cvar = margin == "receiver"
		x = predict(f, h = 1, type = "link", by_draw = TRUE, seed = 51)
		expect_lt(max(abs(apply(x, 4, mean))), 1e-12)
		if (margin == "neither") expect_equal(as.numeric(x), rep(0, length(x)))
		if (margin == "sender") {
			expect_equal(x[, 1, , ], x[, 2, , ], tolerance = 0)
			expect_equal(var(as.numeric(x[1, 1, , ])), .25 * (1 - 1 / 3), tolerance = .015)
		}
		if (margin == "receiver") {
			expect_equal(x[1, , , ], x[2, , , ], tolerance = 0)
			expect_equal(var(as.numeric(x[1, 1, , ])), .25 * (1 - 1 / 4), tolerance = .015)
		}
	}
})
