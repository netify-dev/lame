skip_on_cran()

chains_panel = function(n = 12L, nt = 3L, seed = 41) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(nt), function(t) { m = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	X = lapply(seq_len(nt), function(t) array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x")))
	list(Y = Y, X = X)
}
cross_section = function(n = 12L, seed = 42) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	Y
}
quiet = function(expr) suppressMessages(suppressWarnings(expr))
fit_ame = function(Y, ...) quiet(ame(Y, R = 1, burn = 20, nscan = 100, odens = 5, seed = 7, verbose = FALSE, ...))

test_that("the default is two chains and chains = 1 is the single-chain fit", {
	Y = cross_section()
	withr::local_options(lame.chains = NULL)
	f2 = fit_ame(Y)
	expect_identical(f2$chains, 2L)
	expect_identical(f2$chain_seeds, c(7, 1007))
	expect_equal(nrow(f2$BETA), 40L)
	expect_identical(unique(f2$chain_indicator), 1:2)
	f1 = fit_ame(Y, chains = 1)
	expect_null(f1$chains)
	expect_null(f1$chain_summary)
	expect_equal(nrow(f1$BETA), 20L)
	# chain 1 of the pooled fit is the single-chain fit of the same seed
	expect_identical(f2$BETA[1:20, , drop = FALSE], f1$BETA)
	expect_identical(f2$VC[1:20, , drop = FALSE], f1$VC)
	expect_identical(as.list(f2$call)$seed, 7)
	expect_identical(as.list(f2$call)$Y, quote(Y))
})

test_that("the pooled fit carries the chain record and works with the methods", {
	skip_on_cran()
	d = chains_panel()
	f = quiet(lame(
		d$Y,
		Xdyad = d$X,
		R = 1,
		dynamic_ab = TRUE,
		burn = 100,
		nscan = 1500,
		odens = 5,
		seed = 7,
		chains = 2,
		prior = list(chain_extensions = 0)
	))
	expect_s3_class(f, "lame")
	expect_identical(f$chains, 2L)
	expect_equal(dim(f$BETA)[1], 600L)
	expect_s3_class(f$rhat, "data.frame")
	# under dynamic_ab the intercept is a path, monitored per period
	expect_true(all(c("intercept[t1]", "sigma_ab") %in% f$rhat$parameter))
	expect_identical(f$auto_extend$chains, 2L)
	expect_length(f$chain_start_vals, 2L)
	expect_match(paste(capture.output(print(f)), collapse = " "), "2 chains")
	expect_no_error(quiet(capture.output(summary(f))))
	expect_no_error(quiet(predict(f)))
	expect_no_error(quiet(simulate(f, nsim = 2)))
	skip_if_not_installed("posterior")
	dr = posterior::as_draws(f)
	expect_equal(posterior::nchains(dr), 2L)
	expect_equal(posterior::niterations(dr), 300L)
})

test_that("the results do not depend on the number of cores", {
	skip_on_cran()
	skip_on_os("windows")
	Y = cross_section()
	a = fit_ame(Y, chains = 2, cores = 1)
	b = fit_ame(Y, chains = 2, cores = 2)
	expect_identical(a$BETA, b$BETA)
	expect_identical(a$VC, b$VC)
})

test_that("Windows runs the chains one after the other", {
	Y = cross_section()
	local_mocked_bindings(.lame_chains_os = function() "windows", .package = "lame")
	expect_identical(lame:::.lame_chains_cores(8L, 4L), 1L)
	a = fit_ame(Y, chains = 2, cores = 2)
	local_mocked_bindings(.lame_chains_os = function() "unix", .package = "lame")
	expect_identical(lame:::.lame_chains_cores(8L, 4L), 4L)
	b = fit_ame(Y, chains = 2, cores = 1)
	expect_identical(a$BETA, b$BETA)
})

test_that("ame_parallel() and lame_parallel() are the chains of ame() and lame()", {
	Y = cross_section()
	a = quiet(
		ame_parallel(Y, n_chains = 2, cores = 1, R = 1, burn = 20, nscan = 100, odens = 5, seed = 7, verbose = FALSE)
	)
	b = fit_ame(Y, chains = 2)
	expect_identical(a$BETA, b$BETA)
	expect_identical(a$chain_seeds, b$chain_seeds)
	d = chains_panel()
	l1 = quiet(
		lame_parallel(
			d$Y,
			n_chains = 2,
			cores = 1,
			Xdyad = d$X,
			R = 0,
			burn = 20,
			nscan = 100,
			odens = 5,
			seed = 7,
			verbose = FALSE
		)
	)
	l2 = quiet(lame(d$Y, Xdyad = d$X, R = 0, burn = 20, nscan = 100, odens = 5, seed = 7, verbose = FALSE, chains = 2))
	expect_identical(l1$BETA, l2$BETA)
	lst = quiet(
		ame_parallel(
			Y,
			n_chains = 2,
			cores = 1,
			combine_method = "list",
			R = 1,
			burn = 20,
			nscan = 100,
			odens = 5,
			seed = 7,
			verbose = FALSE
		)
	)
	expect_length(lst, 2L)
	expect_identical(lst[[2]]$BETA, b$BETA[21:40, , drop = FALSE])
})

test_that("the chains argument is checked", {
	Y = cross_section()
	expect_error(fit_ame(Y, chains = 0), class = "lame_chains_arg")
	expect_error(fit_ame(Y, chains = 1.5), class = "lame_chains_arg")
	d = chains_panel()
	expect_error(quiet(lame(d$Y, R = 0, burn = 5, nscan = 20, odens = 5, chains = "a")), class = "lame_chains_arg")
	expect_error(
		quiet(lame(d$Y, R = 0, burn = 5, nscan = 20, odens = 5, chains = 0, sampler = "tempered")),
		class = "lame_chains_arg"
	)
})

test_that("a split rank-normalised R-hat agrees with the posterior package", {
	skip_if_not_installed("posterior")
	set.seed(5)
	for (i in 1:20) {
		n = sample(c(60, 200), 1); m = sample(2:4, 1)
		x = matrix(rnorm(n * m), n, m) + rep(rnorm(m, 0, runif(1)), each = n)
		expect_equal(lame:::.lame_rhat(x), as.numeric(posterior::rhat(x)))
	}
})

test_that("the gate compares the chains on the monitored set and not on the products", {
	set.seed(6)
	mk = function(shift) list(
		BETA = cbind(intercept = rnorm(100) + shift),
		VC = cbind(va = rnorm(100, 1, .1)),
		family = "normal",
		mode = "unipartite"
	)
	tab = lame:::.lame_gate_table(list(mk(0), mk(0)))
	expect_true(all(c("intercept", "va") %in% tab$parameter))
	expect_named(tab, c("parameter", "rhat", "ess", "ess_chain", "threshold", "informative", "ess_needed"))
	expect_true(all(tab$informative))
	expect_true(all(tab$threshold >= 1.05))
	expect_equal(nrow(lame:::.lame_gate_flagged(tab)), 0L)
	tab = lame:::.lame_gate_table(list(mk(0), mk(2)))
	flagged = lame:::.lame_gate_flagged(tab)
	expect_identical(flagged$parameter, "intercept")
	expect_gt(flagged$rhat, flagged$threshold)
	# short runs are not compared
	short = function() list(BETA = cbind(intercept = rnorm(30)), family = "normal", mode = "unipartite")
	expect_equal(nrow(lame:::.lame_gate_table(list(short(), short()))), 0L)
})

test_that("chains too short to compare get a longer block and are returned, not refused", {
	Y = cross_section()
	short = data.frame(
		parameter = "va",
		rhat = 1.01,
		ess = 10,
		ess_chain = 10,
		threshold = 1.2,
		informative = FALSE,
		ess_needed = 35,
		stringsAsFactors = FALSE
	)
	calls = 0L
	local_mocked_bindings(.lame_gate_table = function(fits) {
		calls <<- calls + 1L
		if (calls == 1L) short else short[0, ]
	}, .package = "lame")
	f = fit_ame(Y, chains = 2)
	# 100 sweeps held 10 of 35 effective draws: 100 * 1.25 * 3.5 = 440 sweeps, capped at 4 * nscan
	expect_identical(f$auto_extend$extensions, 1L)
	expect_identical(f$auto_extend$block_nscan, 400L)
	expect_identical(f$auto_extend$discarded_sweeps, 100L)
	expect_equal(nrow(f$BETA), 2L * 80L)
	expect_equal(f$auto_extend$next_seeds, c(7, 1007) + 20 + 100 + 400)
	# still too short after the last extension: returned with a warning that names the series
	local_mocked_bindings(.lame_gate_table = function(fits) short, .package = "lame")
	expect_warning(
		g <- suppressMessages(ame(Y, R = 1, burn = 20, nscan = 100, odens = 5, seed = 7, verbose = FALSE, chains = 2)),
		class = "lame_chain_not_compared"
	)
	expect_s3_class(g, "ame")
	expect_identical(g$auto_extend$extensions, 3L)
	expect_identical(g$auto_extend$block_nscan, 400L)
	expect_false(isTRUE(g$auto_extend$refused))
})

test_that("chains that disagree are extended together and then refused with every chain", {
	Y = cross_section()
	bad = data.frame(parameter = "intercept", rhat = 1.4, stringsAsFactors = FALSE)
	local_mocked_bindings(.lame_gate_table = function(fits) bad, .package = "lame")
	e = tryCatch(fit_ame(Y, chains = 2), error = function(e) e)
	expect_s3_class(e, "lame_chain_unsettled")
	expect_identical(e$extensions, 3L)
	expect_length(e$fits, 2L)
	expect_equal(e$rhat$rhat, 1.4)
	expect_identical(e$fit$chains, 2L)
	expect_match(conditionMessage(e), "intercept")
	expect_match(conditionMessage(e), "R-hat 1.400 above 1.050", fixed = TRUE)
	expect_match(conditionMessage(e), "lame_resume", fixed = TRUE)
	expect_identical(e$fit$auto_extend$alpha, lame:::.lame_chain_alpha)
	# no extension at all
	e0 = tryCatch(fit_ame(Y, chains = 2, prior = list(chain_extensions = 0)), error = function(e) e)
	expect_identical(e0$extensions, 0L)
	# a refused fit continues: every chain, or one
	local_mocked_bindings(
		.lame_gate_table = function(fits) data.frame(parameter = character(0), rhat = numeric(0)),
		.package = "lame"
	)
	r = quiet(lame_resume(e$fit, nscan = 50))
	expect_identical(r$chains, 2L)
	expect_equal(nrow(r$BETA), 20L)
	r1 = quiet(lame_resume(e$fits[[2]], nscan = 50))
	expect_null(r1$chains)
	expect_equal(nrow(r1$BETA), 10L)
})

test_that("a pooled fit that settled resumes every chain from its own state", {
	Y = cross_section()
	f = fit_ame(Y, chains = 2)
	r = quiet(lame_resume(f, nscan = 50))
	expect_identical(r$chains, 2L)
	expect_equal(nrow(r$BETA), 20L)
	expect_identical(r$chain_seeds, f$auto_extend$next_seeds)
	expect_false(identical(r$BETA[1:10, ], f$BETA[1:10, ]))
})

test_that("a resumed two-chain lame() fit holds nscan_more / odens draws per chain, and a fit stopped by max_seconds keeps its time limit unless lifted", {
	skip_on_cran()
	P = chains_panel()
	fl = quiet(
		lame(
			P$Y,
			P$X,
			family = "normal",
			R = 0,
			burn = 20,
			nscan = 100,
			odens = 5,
			seed = 7,
			verbose = FALSE,
			chains = 2
		)
	)
	rl = quiet(lame_resume(fl, nscan_more = 200))
	expect_equal(nrow(rl$BETA), 80L)
	cut = quiet(
		lame(
			P$Y,
			P$X,
			family = "normal",
			R = 0,
			burn = 0,
			nscan = 100000,
			odens = 5,
			seed = 7,
			verbose = FALSE,
			chains = 2,
			max_seconds = 0.2
		)
	)
	expect_true(cut$terminated_early)
	# the resume inherits the time limit, so more scans than fit in it force the stop
	again = quiet(lame_resume(cut, nscan_more = 100000))
	expect_true(again$terminated_early)
	expect_lt(nrow(again$BETA), nrow(cut$BETA) + 20000L)
	full = quiet(lame_resume(cut, nscan_more = 200, max_seconds = Inf))
	expect_false(isTRUE(full$terminated_early))
	expect_equal(nrow(full$BETA), 80L)
	full0 = quiet(lame_resume(cut, nscan_more = 200, max_seconds = NULL))
	expect_equal(nrow(full0$BETA), 80L)
})

test_that("a drifting chain extends every chain by a block", {
	Y = cross_section()
	calls = 0L
	local_mocked_bindings(.lame_chain_moving = function(fit, ...) {
		calls <<- calls + 1L
		if (calls == 2L) data.frame(
			parameter = "va", mean_first = 0, mean_last = 1, z = 9, shift = 3, text = "va (0 to 1)", stringsAsFactors = FALSE
		)
		else data.frame(parameter = character(0), text = character(0))
	}, .package = "lame")
	f = fit_ame(Y, chains = 2)
	expect_identical(f$auto_extend$extensions, 1L)
	expect_match(f$auto_extend$triggers[[1]], "chain 2: va")
	expect_equal(nrow(f$BETA), 40L)
})

test_that("a multi-chain fit leaves the global random stream where it found it", {
	Y = cross_section()
	for (cores in c(1, 2)) {
		set.seed(99); r0 = runif(1); set.seed(99)
		fit_ame(Y, chains = 2, cores = cores)
		expect_identical(runif(1), r0)
	}
})

test_that("a single-period lame() fit pools its chains and plots its goodness of fit", {
	Y = cross_section()
	f = quiet(lame(list(Y), R = 0, burn = 20, nscan = 100, odens = 5, seed = 7, verbose = FALSE, chains = 2))
	expect_identical(unique(f$chain_indicator), 1:2)
	# one period row; the observed column plus 20 draws from each chain
	expect_true(all(vapply(f$GOF, function(g) identical(dim(g), c(1L, 41L)), logical(1))))
	expect_s3_class(quiet(gof_plot(f)), "ggplot")
})

test_that("chunked likelihood storage keeps a directory per chain and reads back the pooled draws", {
	d = chains_panel()
	ll_dir = withr::local_tempdir()
	args = list(
		Y = d$Y, Xdyad = d$X, R = 1, dcor = FALSE, burn = 50, nscan = 400, odens = 5, seed = 7,
		verbose = FALSE, chains = 2
	)
	mem = quiet(do.call(lame, c(args, list(save_log_lik = TRUE))))
	chunk = quiet(do.call(lame, c(args, list(
		save_log_lik = "chunked", log_lik_path = file.path(ll_dir, "ll"), log_lik_chunk_size = 100L
	))))
	expect_setequal(list.files(ll_dir), c("ll_chain1", "ll_chain2"))
	expect_null(chunk$log_lik)
	expect_identical(chunk$BETA, mem$BETA)
	expect_identical(nrow(read_log_lik(chunk)), nrow(chunk$BETA))
	expect_identical(unname(read_log_lik(chunk)), unname(mem$log_lik))
	if (requireNamespace("loo", quietly = TRUE)) {
		expect_equal(suppressWarnings(loo(chunk))$estimates, suppressWarnings(loo(mem))$estimates)
	}
	# a resumed block has its own directories, and the fit it continues still reads its own draws
	more = quiet(lame_resume(chunk, nscan_more = 100))
	more_mem = quiet(lame_resume(mem, nscan_more = 100))
	expect_true(all(c("ll_resume1_chain1", "ll_resume1_chain2") %in% list.files(ll_dir)))
	expect_identical(unname(read_log_lik(more)), unname(more_mem$log_lik))
	expect_identical(unname(read_log_lik(chunk)), unname(mem$log_lik))
	# chains pooled by hand, with either storage
	one = function(seed, ...) quiet(lame(
		d$Y, Xdyad = d$X, R = 1, dcor = FALSE, burn = 20, nscan = 100, odens = 5, seed = seed,
		verbose = FALSE, chains = 1, ...
	))
	ref = rbind(one(7, save_log_lik = TRUE)$log_lik, one(1007, save_log_lik = TRUE)$log_lik)
	both = quiet(combine_ame_chains(list(one(7, save_log_lik = "chunked"), one(1007, save_log_lik = "chunked"))))
	expect_identical(read_log_lik(both), ref)
	mixed = quiet(combine_ame_chains(list(one(7, save_log_lik = "chunked"), one(1007, save_log_lik = TRUE))))
	expect_identical(unname(read_log_lik(mixed)), ref)
})

test_that("a time limit leaves chunked likelihood storage with the draws that were kept", {
	d = chains_panel()
	stopped = quiet(lame(
		d$Y, Xdyad = d$X, R = 1, dcor = FALSE, burn = 20, nscan = 50000, odens = 5, seed = 7,
		verbose = FALSE, chains = 1, save_log_lik = "chunked", max_seconds = 4
	))
	ll = read_log_lik(stopped)
	expect_equal(stopped$log_lik_chunks$n_stored, nrow(stopped$BETA))
	expect_identical(nrow(ll), nrow(stopped$BETA))
	expect_true(all(is.finite(ll)))
})

test_that("a 2-d dyadic covariate is named after its argument on every route", {
	Y = cross_section()
	set.seed(3)
	Xd = matrix(rnorm(length(Y)), nrow(Y), dimnames = dimnames(Y))
	nm = c("intercept", "Xd_dyad")
	expect_identical(names(coef(fit_ame(Y, Xdyad = Xd))), nm)
	expect_identical(names(coef(fit_ame(Y, Xdyad = Xd, chains = 1))), nm)
	expect_identical(names(coef(quiet(ame(Y, Xdyad = Xd, method = "als", verbose = FALSE)))), nm)
	expect_identical(names(coef(quiet(ame_als(Y, Xdyad = Xd, verbose = FALSE)))), nm)
})
