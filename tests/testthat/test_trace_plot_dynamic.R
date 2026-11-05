skip_on_cran()

# trace_plot() draws the dynamic hyperparameters, leaves out the quantities a
# model fixes, and applies burn-in and thinning within each chain

tp_data = function(p) if (inherits(p, "patchwork")) p[[1]]$data else p$data

test_that("trace_plot draws the persistence and innovation of dynamic effects and positions", {
	dat = simulate_test_network(n = 10, n_time = 3, family = "normal", seed = 11)
	fit = lame(
		dat$Y,
		R = 2,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 10,
		nscan = 60,
		odens = 2,
		verbose = FALSE,
		seed = 11
	)
	d = tp_data(trace_plot(fit))
	dyn = c(
		"Effect persistence (rho_ab)", "Effect innovation (sigma_ab)",
		"Position persistence (rho_uv)", "Position innovation (sigma_uv)"
	)
	expect_true(all(dyn %in% d$parameter))
	expect_equal(d$value[d$parameter == "Position persistence (rho_uv)"], as.numeric(fit$rho_uv))
	d_dyn = tp_data(trace_plot(fit, params = "dynamic"))
	expect_setequal(unique(d_dyn$parameter), dyn)
	# raw names and labels both filter
	d_raw = tp_data(trace_plot(fit, include = c("rho_uv", "sigma_ab")))
	expect_setequal(unique(d_raw$parameter), c("Position persistence (rho_uv)", "Effect innovation (sigma_ab)"))
	# display labels filter regardless of case
	d_ex = tp_data(trace_plot(fit, params = "dynamic", exclude = "Position Innovation (sigma_uv)"))
	expect_false("Position innovation (sigma_uv)" %in% d_ex$parameter)
	expect_true("Position persistence (rho_uv)" %in% d_ex$parameter)
})

test_that("trace_plot draws the per-group hyperparameters of dynamic coefficients", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 12)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		R = 0,
		dynamic_beta = "intercept",
		burn = 5,
		nscan = 40,
		odens = 2,
		verbose = FALSE,
		seed = 12
	)
	d = tp_data(trace_plot(fit, params = "dynamic"))
	expect_true(any(grepl("^Coefficient persistence \\(rho_beta\\)", d$parameter)))
	expect_true(any(grepl("^Coefficient innovation \\(sigma_beta\\)", d$parameter)))
	d_inc = tp_data(trace_plot(fit, include = "rho_beta"))
	expect_true(all(grepl("^Coefficient persistence", d_inc$parameter)))
})

test_that("trace_plot leaves out the quantities a binary bipartite model fixes", {
	set.seed(13)
	Y = lapply(1:3, function(t) {
		m = matrix(
			rbinom(8 * 11, 1, 0.3),
			8,
			11,
			dimnames = list(paste0("r", 1:8), paste0("c", 1:11))
		)
		m
	})
	fit = lame(
		Y,
		R = 0,
		family = "binary",
		mode = "bipartite",
		burn = 10,
		nscan = 60,
		odens = 2,
		verbose = FALSE,
		seed = 13
	)
	expect_true(all(fit$VC[, "ve"] == 1))
	d = tp_data(trace_plot(fit, params = "variance"))
	expect_false(any(c("Error variance", "Dyadic correlation", "Sender-receiver covariance") %in% d$parameter))
	expect_true(all(c("Sender variance", "Receiver variance") %in% d$parameter))
})

test_that("trace_plot keeps an estimated quantity whose chain is stuck", {
	dat = simulate_test_network(n = 10, n_time = 1, family = "normal", seed = 14)
	fit = ame(dat$Y, R = 0, burn = 10, nscan = 60, odens = 2, verbose = FALSE, seed = 14)
	fit$VC[, "va"] = fit$VC[1, "va"]
	fit$VC[, "rho"] = 0
	d = tp_data(trace_plot(fit, params = "variance"))
	expect_true(all(c("Sender variance", "Dyadic correlation") %in% d$parameter))
})

test_that("trace_plot applies burn-in and thinning within each chain", {
	dat = simulate_test_network(n = 10, n_time = 1, family = "normal", seed = 15)
	fit = ame(
		dat$Y,
		R = 0,
		burn = 10,
		nscan = 120,
		odens = 2,
		verbose = FALSE,
		seed = 15,
		chains = 2
	)
	n_per = table(fit$chain_indicator)
	d = tp_data(trace_plot(fit, params = "beta", burn.in = 5, thin = 2))
	one = d[d$parameter == d$parameter[1], ]
	kept = table(one$chain)
	expect_equal(as.numeric(kept), as.numeric(ceiling((n_per - 5) / 2)))
	expect_equal(
		one$value[one$chain == "2"][1],
		unname(fit$BETA[fit$chain_indicator == 2, 1][6])
	)
})

test_that("trace_plot leaves out what summary() reports as not estimated, with a message", {
	dat = simulate_test_network(n = 10, n_time = 3, family = "normal", seed = 16)
	fit = lame(dat$Y, R = 0, burn = 10, nscan = 40, odens = 2, verbose = FALSE, seed = 16)
	local_mocked_bindings(
		.lame_ab_var_unident_reason = function(fit) "sparse",
		.lame_ab_coef_reason = function(fit) "sparse",
		.package = "lame"
	)
	expect_message(p <- trace_plot(fit), "not estimated")
	d = tp_data(p)
	expect_false(any(c("Sender variance", "Receiver variance", "intercept") %in% d$parameter))
})
