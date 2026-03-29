skip_on_cran()

# plotting functions

####
# setup: fit representative models
####

fit_ame_uni = local({
	n = 15
	dat = simulate_test_network(
		n = n, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.3,
		mode = "unipartite", seed = 5001
	)
	ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		burn = 50, nscan = 100, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_ame_bip = local({
	nA = 12; nB = 8
	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.5,
		mode = "bipartite", seed = 5002
	)
	ame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		mode = "bipartite", burn = 50, nscan = 100, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_lame_uni = local({
	n = 15
	dat = simulate_test_network(
		n = n, n_time = 3, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.3,
		mode = "unipartite", seed = 5003
	)
	lame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		burn = 50, nscan = 100, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_lame_bip = local({
	nA = 12; nB = 8
	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 3, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.3,
		mode = "bipartite", seed = 5004
	)
	Xdyad_list = lapply(dat$Xdyad, function(x) x[, , 1])
	lame(
		dat$Y, Xdyad = Xdyad_list, R = 2, family = "normal",
		mode = "bipartite", burn = 50, nscan = 100, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

fit_lame_dyn = local({
	n = 15
	dat = simulate_test_network(
		n = n, n_time = 3, family = "normal", R = 2,
		beta_intercept = 1, beta_dyad = 0.3,
		rho_ab = 0.7, rho_uv = 0.5,
		mode = "unipartite", seed = 5005
	)
	lame(
		dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
		dynamic_ab = TRUE, dynamic_uv = TRUE,
		burn = 50, nscan = 100, odens = 1,
		verbose = FALSE, gof = TRUE, seed = 42
	)
})

####
# trace_plot()
####

test_that("trace_plot works for ame unipartite", {
	p = trace_plot(fit_ame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("trace_plot works for ame bipartite", {
	p = trace_plot(fit_ame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("trace_plot works for lame unipartite", {
	p = trace_plot(fit_lame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("trace_plot works for lame bipartite", {
	p = trace_plot(fit_lame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("trace_plot works for lame dynamic", {
	p = trace_plot(fit_lame_dyn)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

####
# gof_plot()
####

test_that("gof_plot works for ame unipartite", {
	p = gof_plot(fit_ame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("gof_plot works for ame bipartite", {
	p = gof_plot(fit_ame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("gof_plot works for lame unipartite", {
	p = gof_plot(fit_lame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

test_that("gof_plot works for lame bipartite", {
	p = gof_plot(fit_lame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot") || is.list(p))
})

####
# ab_plot()
####

test_that("ab_plot works for ame unipartite", {
	p = ab_plot(fit_ame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ab_plot works for ame bipartite", {
	p = ab_plot(fit_ame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ab_plot works for lame unipartite", {
	p = ab_plot(fit_lame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ab_plot works for lame dynamic", {
	# default mode for dynamic should be trajectory
	p = suppressWarnings(ab_plot(fit_lame_dyn))
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ab_plot trajectory plot_type for lame dynamic", {
	p = suppressWarnings(ab_plot(fit_lame_dyn, plot_type = "trajectory"))
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ab_plot snapshot plot_type for lame dynamic", {
	p = suppressWarnings(ab_plot(fit_lame_dyn, plot_type = "snapshot"))
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

####
# uv_plot()
####

test_that("uv_plot works for ame unipartite", {
	p = uv_plot(fit_ame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("uv_plot works for ame bipartite", {
	p = uv_plot(fit_ame_bip)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("uv_plot works for lame unipartite", {
	p = uv_plot(fit_lame_uni)
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("uv_plot circle layout for ame unipartite", {
	p = uv_plot(fit_ame_uni, layout = "circle")
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("uv_plot biplot layout for ame unipartite", {
	p = uv_plot(fit_ame_uni, layout = "biplot")
	expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

####
# plot.ame() and plot.lame()
####

test_that("plot.ame works for unipartite", {
	expect_no_error(suppressWarnings(plot(fit_ame_uni)))
})

test_that("plot.ame works for bipartite", {
	expect_no_error(suppressWarnings(plot(fit_ame_bip)))
})

test_that("plot.lame works for unipartite", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni)))
})

test_that("plot.lame works for bipartite", {
	expect_no_error(suppressWarnings(plot(fit_lame_bip)))
})

test_that("plot.lame works for dynamic", {
	expect_no_error(suppressWarnings(plot(fit_lame_dyn)))
})

####
# plot.lame which= options
####

test_that("plot.lame which='trace' works", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni, which = "trace")))
})

test_that("plot.lame which='density' works", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni, which = "density")))
})

test_that("plot.lame which='gof' works", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni, which = "gof")))
})

test_that("plot.lame which='effects' works for unipartite", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni, which = "effects")))
})

test_that("plot.lame which='effects' works for bipartite", {
	expect_no_error(suppressWarnings(plot(fit_lame_bip, which = "effects")))
})

test_that("plot.lame which='network' works", {
	expect_no_error(suppressWarnings(plot(fit_lame_uni, which = "network")))
})
