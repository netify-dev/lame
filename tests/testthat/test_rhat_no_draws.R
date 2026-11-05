skip_on_cran()

# every reported R-hat is the gate's statistic on the gate's draws, and a fit
# that stored no draw has empty draw arrays and methods that refuse it

rnd_section = function(n = 12L, seed = 42) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	X = array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x"))
	list(Y = Y, X = X)
}
rnd_panel = function(n = 10L, nt = 3L, seed = 41) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(nt), function(t) { m = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	names(Y) = paste0("t", seq_len(nt))
	X = lapply(seq_len(nt), function(t) array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x")))
	names(X) = names(Y)
	list(Y = Y, X = X)
}
quiet_rh = function(expr) suppressMessages(suppressWarnings(expr))
# split rank-normalised R-hat of one column of the pooled draws, per chain
rhat_pooled = function(draws, chain) lame:::.lame_rhat(do.call(cbind, split(draws, chain)))

test_that("summary's Rhat is fit$rhat for every monitored quantity of a two-chain fit", {
	d = rnd_section()
	f = quiet_rh(ame(d$Y, Xdyad = d$X, R = 1, burn = 20, nscan = 300, odens = 5, seed = 7, verbose = FALSE, chains = 2))
	expect_identical(f$chains, 2L)
	expect_gt(nrow(f$rhat), 0L)
	s = summary(f)
	tab = rbind(s$beta[, "Rhat", drop = FALSE], s$variance[, "Rhat", drop = FALSE])
	hit = intersect(rownames(tab), f$rhat$parameter)
	expect_true(all(c(colnames(f$BETA), "va", "ve") %in% hit))
	expect_identical(unname(tab[hit, "Rhat"]), f$rhat$rhat[match(hit, f$rhat$parameter)])
	# compute_mcmc_diagnostics() gives the same numbers without the override
	expect_identical(unname(f$diagnostics$rhat[hit]), f$rhat$rhat[match(hit, f$rhat$parameter)])
	# and they are the split rank-normalised R-hat of the draws the fit holds
	for (v in colnames(f$BETA)) expect_equal(unname(s$beta[v, "Rhat"]), rhat_pooled(f$BETA[, v], f$chain_indicator))
	# a settled fit shows no monitored Rhat at or above the gate's threshold
	expect_true(all(tab[hit, "Rhat"] < 1.05, na.rm = TRUE))
	# the classic statistic stays available under its own name
	expect_named(f$diagnostics$rhat_classic, names(f$diagnostics$rhat))
})

test_that("after an extension the gate, summary() and the returned draws are the last block", {
	d = rnd_section()
	calls = 0L
	local_mocked_bindings(.lame_chain_moving = function(fit, ...) {
		calls <<- calls + 1L
		if (calls == 1L) data.frame(
			parameter = "va", mean_first = 0, mean_last = 1, z = 9, shift = 3, text = "va (0 to 1)", stringsAsFactors = FALSE
		)
		else data.frame(parameter = character(0), text = character(0))
	}, .package = "lame")
	f = quiet_rh(ame(d$Y, R = 0, burn = 20, nscan = 300, odens = 5, seed = 7, verbose = FALSE, chains = 2))
	expect_identical(f$auto_extend$extensions, 1L)
	expect_equal(nrow(f$BETA), 2L * 60L)
	expect_equal(f$rhat$rhat[f$rhat$parameter == "intercept"], rhat_pooled(f$BETA[, "intercept"], f$chain_indicator))
	expect_identical(unname(summary(f)$beta["intercept", "Rhat"]), f$rhat$rhat[f$rhat$parameter == "intercept"])
})

test_that("the two-chain lame() fit reports the gate's R-hat in its diagnostics", {
	skip_on_cran()
	p = rnd_panel()
	f = quiet_rh(lame(p$Y, p$X, R = 0, burn = 20, nscan = 300, odens = 5, seed = 7, verbose = FALSE, chains = 2))
	hit = intersect(names(f$diagnostics$rhat), f$rhat$parameter)
	expect_true(length(hit) >= 3L)
	expect_identical(unname(f$diagnostics$rhat[hit]), f$rhat$rhat[match(hit, f$rhat$parameter)])
})

test_that("a one-chain fit reports no R-hat, as before", {
	d = rnd_section()
	f = quiet_rh(ame(d$Y, R = 0, burn = 10, nscan = 100, odens = 5, seed = 7, verbose = FALSE, chains = 1))
	expect_null(f$rhat)
	expect_null(f$diagnostics)
	s = summary(f)
	expect_false("Rhat" %in% colnames(s$beta))
	expect_false("Rhat" %in% colnames(s$variance))
	expect_null(s$diagnostics)
	expect_warning(out <- compute_mcmc_diagnostics(list(f)), "at least 2 chains")
	expect_null(out)
})

test_that("compute_mcmc_diagnostics() and rhat_dynamic_beta() use the gate's statistic", {
	set.seed(9)
	ch = lapply(
		1:2,
		function(k) list(
			BETA = cbind(b1 = rnorm(80) + k / 4, b2 = rnorm(80)),
			VC = cbind(va = rnorm(80, 1), ve = rep(1, 80))
		)
	)
	dg = quiet_rh(compute_mcmc_diagnostics(ch))
	for (v in c("b1", "b2")) {
		x = cbind(ch[[1]]$BETA[, v], ch[[2]]$BETA[, v])
		expect_identical(unname(dg$rhat[v]), lame:::.lame_rhat(x))
		W = mean(apply(x, 2, var)); B = 80 * var(colMeans(x))
		expect_equal(unname(dg$rhat_classic[v]), sqrt(((80 - 1) * W + B) / 80 / W))
	}
	# a quantity fixed in every chain has no R-hat
	expect_true(is.na(dg$rhat[["ve"]]))
	expect_identical(lame:::.lame_rhat_chains(matrix(1, 50, 2)), NA_real_)
	expect_identical(lame:::.lame_rhat_chains(matrix(rnorm(10), 5, 2)), NA_real_)
	# unequal chains are compared on their common length, as the gate does
	ch[[2]]$BETA = rbind(ch[[2]]$BETA, ch[[2]]$BETA[1:10, ])
	dg2 = quiet_rh(compute_mcmc_diagnostics(lapply(ch, function(z) z["BETA"])))
	expect_identical(unname(dg2$rhat["b1"]), unname(dg$rhat["b1"]))
	# the univariate column of rhat_dynamic_beta()
	mk = function(shift) structure(
		list(
			BETA = array(
				rnorm(60 * 2 * 3) + shift,
				c(60, 2, 3),
				dimnames = list(NULL, c("intercept", "x"), c("t1", "t2", "t3"))
			)
		),
		class = c("lame", "ame")
	)
	fl = list(mk(0), mk(.3))
	rd = rhat_dynamic_beta(fl)
	for (k in 1:2) {
		per = vapply(1:3, function(t) lame:::.lame_rhat(cbind(fl[[1]]$BETA[, k, t], fl[[2]]$BETA[, k, t])), numeric(1))
		expect_equal(rd$rhat_max_univariate[k], max(per))
	}
})

test_that("a lame() fit stopped before its first stored draw has empty draw arrays and refusing methods", {
	p = rnd_panel()
	ref = quiet_rh(lame(p$Y, p$X, R = 1, burn = 10, nscan = 50, odens = 5, seed = 3, verbose = FALSE, chains = 1))
	expect_warning(
		f <- suppressMessages(lame(
			p$Y,
			p$X,
			R = 1,
			burn = 10,
			nscan = 1000,
			odens = 5,
			seed = 3,
			verbose = FALSE,
			chains = 1,
			max_seconds = 1e-12
		)),
		"before any draw was stored"
	)
	expect_true(f$terminated_early)
	expect_true(f$no_draws_stored)
	expect_identical(dim(f$BETA), c(0L, ncol(ref$BETA)))
	expect_identical(colnames(f$BETA), colnames(ref$BETA))
	expect_identical(dim(f$VC), c(0L, ncol(ref$VC)))
	expect_identical(colnames(f$VC), colnames(ref$VC))
	expect_identical(dim(f$U_samples)[3], 0L)
	expect_true(all(is.na(f$APM)))
	expect_true(all(is.na(f$UVPM)))
	expect_identical(dimnames(f$UVPM), dimnames(ref$UVPM))
	for (m in list(coef, vcov, confint, summary, predict, fitted, residuals, tidy, as_draws)) {
		expect_error(m(f), class = "lame_no_draws")
	}
	expect_error(simulate(f, nsim = 1), class = "lame_no_draws")
	expect_error(trace_plot(f), class = "lame_no_draws")
	out = paste(capture.output(print(f)), collapse = " ")
	expect_match(out, "No posterior draw was stored")
	# the start values are kept, so the fit continues once the time limit is lifted
	r = quiet_rh(lame_resume(f, nscan_more = 50, max_seconds = Inf))
	expect_false(isTRUE(r$no_draws_stored))
	expect_equal(nrow(r$BETA), 10L)
	expect_true(all(is.finite(coef(r))))
})

test_that("the dynamic draw arrays of a fit with no stored draw are empty along the draw axis", {
	skip_on_cran()
	p = rnd_panel()
	f = quiet_rh(lame(
		p$Y,
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		burn = 10,
		nscan = 1000,
		odens = 5,
		seed = 3,
		verbose = FALSE,
		chains = 1,
		max_seconds = 1e-12,
		posterior_opts = list(save_UV = TRUE, save_ab = TRUE, thin_UV = 1, thin_ab = 1)
	))
	expect_true(f$no_draws_stored)
	expect_length(f$rho_ab, 0L)
	expect_length(f$sigma_uv, 0L)
	expect_identical(dim(f$a_draws)[3], 0L)
	expect_identical(dim(f$U_draws)[4], 0L)
	expect_identical(dim(f$forecast_state$a)[2], 0L)
	expect_true(all(vapply(f$GOF, ncol, integer(1)) == 1L))
	expect_error(summary(f), class = "lame_no_draws")
})

test_that("two chains stopped before their first stored draw pool to an empty fit", {
	p = rnd_panel()
	f = quiet_rh(
		lame(p$Y, R = 0, burn = 10, nscan = 1000, odens = 5, seed = 3, verbose = FALSE, chains = 2, max_seconds = 1e-12)
	)
	expect_identical(f$chains, 2L)
	expect_true(f$no_draws_stored)
	expect_identical(nrow(f$BETA), 0L)
	expect_identical(colnames(f$BETA), "intercept")
	expect_length(f$chain_indicator, 0L)
	expect_equal(nrow(f$rhat), 0L)
	expect_error(coef(f), class = "lame_no_draws")
	expect_match(paste(capture.output(print(f)), collapse = " "), "No posterior draw was stored")
	r = quiet_rh(lame_resume(f, nscan_more = 50, max_seconds = Inf))
	expect_identical(r$chains, 2L)
	expect_equal(nrow(r$BETA), 20L)
})

test_that("a chain with no stored draw does not enter the pooled averages", {
	p = rnd_panel()
	a = quiet_rh(lame(p$Y, R = 1, burn = 10, nscan = 50, odens = 5, seed = 3, verbose = FALSE, chains = 1))
	b = lame:::.lame_empty_draws(
		quiet_rh(lame(p$Y, R = 1, burn = 10, nscan = 50, odens = 5, seed = 4, verbose = FALSE, chains = 1))
	)
	pooled = quiet_rh(combine_ame_chains(list(a, b), diagnostics = FALSE))
	expect_identical(nrow(pooled$BETA), nrow(a$BETA))
	expect_equal(pooled$APM, a$APM)
	expect_equal(pooled$UVPM, a$UVPM)
	expect_identical(as.integer(table(factor(pooled$chain_indicator, levels = 1:2))), c(nrow(a$BETA), 0L))
})
