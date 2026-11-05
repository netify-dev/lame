.health_static_bipartite_args = function(family = "binary", rank = 2L) {
	n = 6L; m = 7L
	Y = matrix(
		rep(c(0, 1, 0, 0, 1, 1, 0), length.out = n * m),
		n,
		m,
		dimnames = list(paste0("sender", seq_len(n)), paste0("receiver", seq_len(m)))
	)
	if (family == "normal") Y = Y + matrix(seq_len(n * m) / 100, n, m)
	# failures are injected into the Gaussian conditional of the interaction matrix, so that sampler runs alone
	# (the collapsed actor update and the scale-ridge draw also move the positions)
	list(
		Y = Y,
		mode = "bipartite",
		family = family,
		R_row = rank,
		R_col = rank,
		burn = 4L,
		nscan = 16L,
		odens = 4L,
		seed = 47337821L,
		gof = FALSE,
		prior = list(actor_move = FALSE, uv_g_scale_move = FALSE),
		verbose = FALSE,
		posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1L, thin_UV = 1L)
	)
}

.health_static_bipartite_fitter = function(mode = "clean") {
	ns = asNamespace("lame"); internal = get("ame_bipartite", ns)
	e = new.env(parent = environment(internal))
	if (mode %in% c("exception", "nonfinite", "invalid_dimensions")) {
		e$sample_G_bip_cpp = switch(
			mode,
			exception = function(...) stop("injected Gaussian draw failure"),
			nonfinite = function(...) matrix(NaN, 2L, 2L),
			invalid_dimensions = function(...) matrix(0, 2L, 3L)
		)
	}
	environment(internal) = e; fitter = get("ame", ns)
	wrapper = new.env(parent = environment(fitter)); wrapper$ame_bipartite = internal
	environment(fitter) = wrapper; fitter
}

.health_static_bipartite_start = function(Y) {
	list(
		Z = 2 * Y - 1,
		beta = 0,
		a = rep(0, nrow(Y)),
		b = rep(0, ncol(Y)),
		Sab = diag(2),
		s2 = 1,
		U = matrix(seq_len(nrow(Y) * 2L) / 100, nrow(Y), 2L),
		V = matrix(seq_len(ncol(Y) * 2L) / 100, ncol(Y), 2L),
		G = diag(2)
	)
}

.health_static_bipartite_capture = function(expr) {
	warnings = list()
	value = withCallingHandlers(expr, warning = function(w) {
		warnings[[length(warnings) + 1L]] <<- list(message = conditionMessage(w), class = class(w))
		invokeRestart("muffleWarning")
	})
	# optional test-only capture keeps the full injected-control outputs
	sink = getOption("lame.test.static_bipartite_health_capture")
	if (is.function(sink)) sink(value, warnings)
	list(value = value, warnings = warnings)
}

test_that("static bipartite ordinary transitions report measured zero numerical failures", {
	for (family in c("binary", "normal")) for (rank in c(0L, 2L)) {
		a = .health_static_bipartite_args(family, rank); set.seed(58337821L); caller = .Random.seed
		z = .health_static_bipartite_capture(do.call(.health_static_bipartite_fitter(), a))
		f = z$value; expect_length(z$warnings, 0L); expect_identical(.Random.seed, caller)
		expect_identical(f$tryErrorChecks, list(U = 0L, V = 0L, G = 0L))
		expect_identical(f$mh_counters, f$tryErrorChecks)
		expect_identical(f$sampler_health_iterations, 20L)
		expect_identical(f$sampler_health$step, c("U", "V", "G"))
		expect_true(all(f$sampler_health$numerical_failures == 0))
		expect_true(all(f$sampler_health$prior_bound_rejections == 0))
	}
})

for (.health_G_mode in c("exception", "nonfinite", "invalid_dimensions")) local({
	mode = .health_G_mode
	test_that(paste("static bipartite Gaussian", mode, "is counted and aborts the fit"), {
		a = .health_static_bipartite_args(); a$start_vals = .health_static_bipartite_start(a$Y)
		set.seed(58337821L); caller = .Random.seed
		e = expect_error(do.call(.health_static_bipartite_fitter(mode), a), class = "lame_sampler_failure")
		expect_identical(.Random.seed, caller)
		expect_identical(e$health$step, c("U", "V", "G"))
		expect_identical(e$health$numerical_failures, c(0, 0, 20))
		expect_true(all(e$health$prior_bound_rejections == 0))
	})
})
