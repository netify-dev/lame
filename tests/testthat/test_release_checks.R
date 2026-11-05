skip_on_cran()

# argument validation, methods on unusual fits, and diagnostics

quiet = function(expr) suppressMessages(suppressWarnings(expr))
warns_of = function(expr) {
	w = character()
	withCallingHandlers(
		suppressMessages(expr),
		warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
	)
	w
}
msgs_of = function(expr) {
	m = character()
	withCallingHandlers(
		suppressWarnings(expr),
		message = function(x) { m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage") }
	)
	m
}
sym_panel = function(n = 12, TT = 3, seed = 3) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	x = rnorm(n)
	Y = X = Xr = vector("list", TT)
	for (t in seq_len(TT)) {
		d = matrix(rnorm(n * n), n, n); d = (d + t(d)) / 2
		e = matrix(rnorm(n * n, 0, 0.5), n, n); e = (e + t(e)) / 2
		y = 0.3 + 0.8 * outer(x, x, "+") + 0.5 * d + e
		diag(y) = NA; dimnames(y) = list(nm, nm)
		Y[[t]] = y
		X[[t]] = array(d, c(n, n, 1), dimnames = list(nm, nm, "d"))
		Xr[[t]] = matrix(x, n, 1, dimnames = list(nm, "x"))
	}
	list(Y = Y, X = X, Xr = Xr, nm = nm)
}
bin_panel = function(n = 12, TT = 3, seed = 6886) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n); diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

# a symmetric fit stores one node coefficient per nodal covariate while the sampler carries row and column slices;
# dynamic_beta indexes the sampler's layout so that every coefficient is updated
test_that("symmetric dynamic_beta with nodal and dyadic covariates updates every coefficient", {
	d = sym_panel()
	fit = quiet(lame(
		d$Y,
		Xdyad = d$X,
		Xrow = d$Xr,
		family = "normal",
		R = 0,
		symmetric = TRUE,
		dynamic_beta = TRUE,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(colnames(fit$BETA), c("intercept", "x_node", "d_dyad"))
	sds = apply(fit$BETA, c(2, 3), sd)
	expect_true(all(sds > 1e-8))
	expect_length(fit$beta_dynamic_mask, 3L)
	fit_d = quiet(lame(
		d$Y,
		Xdyad = d$X,
		Xrow = d$Xr,
		family = "normal",
		R = 0,
		symmetric = TRUE,
		dynamic_beta = "dyad",
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_true(all(apply(fit_d$BETA, c(2, 3), sd) > 1e-8))
	expect_equal(sum(fit_d$beta_dynamic_mask), 1L)
})

test_that("a fit without regression coefficients has empty confint and tidy tables", {
	n = 12; nm = sprintf("a%02d", seq_len(n))
	set.seed(6886)
	Y = lapply(
		1:2,
		function(t) { m = matrix(sample(0:3, n * n, TRUE), n, n); diag(m) = NA; dimnames(m) = list(nm, nm); m }
	)
	fit = quiet(
		lame(Y, family = "ordinal", R = 0, nscan = 20, burn = 5, odens = 1, verbose = FALSE, plot = FALSE, gof = FALSE)
	)
	expect_equal(ncol(fit$BETA), 0L)
	ci = confint(fit)
	expect_true(is.matrix(ci) && nrow(ci) == 0L)
	td = tidy(fit)
	expect_true(is.data.frame(td) && nrow(td) == 0L)
	expect_error(ggplot2::autoplot(fit, which = "beta"), "no regression coefficients")
})

test_that("odmax must be a single value or one value per row", {
	set.seed(6886)
	n = 12; nm = sprintf("a%02d", seq_len(n))
	Y = matrix(0, n, n)
	for (i in seq_len(n)) Y[i, sample(setdiff(seq_len(n), i), 3)] = 1:3
	diag(Y) = NA; dimnames(Y) = list(nm, nm)
	expect_error(
		quiet(ame(
			Y,
			family = "frn",
			odmax = rep(3, 5),
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		"one value per row"
	)
	expect_error(
		quiet(lame(
			list(1 * (Y > 0), 1 * (Y > 0)),
			family = "cbin",
			odmax = rep(3, 5),
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		"one value per row"
	)
})

test_that("ame() drops a one-sided bipartite rank with a warning and uv_plot() says so", {
	set.seed(6886)
	Yb = matrix(rbinom(70, 1, 0.3), 10, 7, dimnames = list(sprintf("r%02d", 1:10), sprintf("c%02d", 1:7)))
	w = warns_of(
		fb <- ame(
			Yb,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 0,
			nscan = 20,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_true(any(grepl("Only one of", w)))
	expect_true(is.null(fb$U) || NCOL(fb$U) == 0L || all(fb$UVPM == 0))
	expect_error(uv_plot(fb), "no latent positions")
})

test_that("predict(newdata) works on a bipartite fit with nodal covariates", {
	set.seed(3)
	Yb = matrix(rbinom(70, 1, 0.3), 10, 7, dimnames = list(sprintf("r%02d", 1:10), sprintf("c%02d", 1:7)))
	Xr = matrix(rnorm(10), 10, 1, dimnames = list(rownames(Yb), "w"))
	Xd = array(rnorm(70), c(10, 7, 1), dimnames = list(rownames(Yb), colnames(Yb), "d"))
	fit = quiet(ame(
		Yb,
		Xdyad = Xd,
		Xrow = Xr,
		family = "binary",
		mode = "bipartite",
		R = 0,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(dimnames(fit$X)[[3]], colnames(fit$BETA))
	p1 = quiet(predict(fit, newdata = Xd, type = "link"))
	p2 = quiet(predict(fit, newdata = Xd + 1, type = "link"))
	expect_equal(dim(p1), c(10L, 7L))
	expect_true(all(is.finite(p1)))
	expect_equal(mean(p2 - p1), mean(fit$BETA[, "d_dyad"]), tolerance = 1e-6)
	# a fit without coefficients simulates without a covariate warning
	f0 = quiet(ame(
		Yb,
		family = "binary",
		mode = "bipartite",
		R = 0,
		intercept = FALSE,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_false(any(grepl("Original covariates", warns_of(simulate(f0, nsim = 2)))))
})

test_that("the settling check names a trending parameter and stays quiet on a settled chain", {
	set.seed(1)
	trending = list(
		BETA = cbind(intercept = seq(0, 1, length.out = 300) + rnorm(300, 0, 0.02)),
		VC = cbind(va = rnorm(300, 1, 0.1))
	)
	moving = lame:::.lame_chain_moving(trending)$text
	expect_equal(sub(" .*", "", moving), "intercept")
	expect_equal(
		nrow(lame:::.lame_chain_moving(list(
			BETA = cbind(intercept = rnorm(300)),
			VC = cbind(va = rnorm(300, 1, 0.1))
		))),
		0L
	)
})

test_that("the settling check leaves out the quantities a fit reports as not estimated", {
	set.seed(2)
	trend = seq(0, 1, length.out = 300) + rnorm(300, 0, 0.02)
	flat = rnorm(300, 0.5, 0.1)
	# a one-period saturated nomination cap withholds the intercept: its
	# drift is not named, a drifting va or dyadic coefficient still is
	sat = list(
		family = "cbin",
		mode = "unipartite",
		odmax_saturated_share = 1,
		n_time = 1,
		BETA = cbind(intercept = trend, x_dyad = flat),
		VC = cbind(va = flat, vb = flat, cab = flat)
	)
	expect_equal(nrow(lame:::.lame_chain_moving(sat)), 0L)
	sat$BETA[, "x_dyad"] = trend
	expect_equal(lame:::.lame_chain_moving(sat)$parameter, "x_dyad")
	sat$VC[, "va"] = trend
	expect_setequal(lame:::.lame_chain_moving(sat)$parameter, c("x_dyad", "va"))
	# the same panel over four periods is checked in full
	open = sat; open$n_time = 4; open$odmax_saturated_share = 0.3; open$BETA[, "x_dyad"] = flat
	moving = lame:::.lame_chain_moving(open)$text
	expect_setequal(sub(" .*", "", moving), c("intercept", "va"))
	# a sparse dynamic_ab panel withholds va and vb and reports rho_ab and sigma_ab
	sparse = list(
		family = "binary",
		mode = "unipartite",
		events_per_actor_period = 0.4,
		dynamic_ab = TRUE,
		BETA = cbind(intercept = flat),
		VC = cbind(va = trend, vb = trend, cab = flat),
		rho_ab = trend,
		sigma_ab = trend
	)
	expect_setequal(
		lame:::.lame_chain_moving(sparse)$parameter,
		c("rho_ab", "sigma_ab")
	)
	# the snap and t kinds withhold rho_uv and sigma_uv (and snap pi_snap) at
	# rank two or more
	snap = list(
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		BETA = cbind(intercept = flat),
		rho_uv = trend,
		sigma_uv = trend,
		pi_snap = trend
	)
	expect_equal(nrow(lame:::.lame_chain_moving(snap)), 0L)
	snap$R = 2
	expect_equal(nrow(lame:::.lame_chain_moving(snap)), 0L)
	# at rank one they are reported, and checked
	snap$R = 1
	expect_setequal(
		lame:::.lame_chain_moving(snap)$parameter,
		c("rho_uv", "sigma_uv", "pi_snap")
	)
	snap$R = NULL; snap$pi_snap = NULL
	snap$dynamic_uv_kind = "ar1"
	expect_setequal(lame:::.lame_chain_moving(snap)$parameter, c("rho_uv", "sigma_uv"))
})

test_that("print() reports a single dynamic component and the symmetric gauge", {
	Y = bin_panel()
	fab = quiet(lame(
		Y,
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	out = paste(c(capture.output(print(fab)), capture.output(print(fab), type = "message")), collapse = "\n")
	expect_true(grepl("Additive effects AR\\(1\\): rho_ab", out))
	expect_false(grepl("Dynamic components", out))
	Ys = lapply(Y, function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m })
	fs = quiet(lame(
		Ys,
		family = "binary",
		R = 1,
		symmetric = TRUE,
		dynamic_uv = TRUE,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	out = paste(c(capture.output(print(fs)), capture.output(print(fs), type = "message")), collapse = "\n")
	expect_true(grepl("gauge", out))
})

test_that("plot() skips the GOF panel with a message when no draws are stored", {
	skip_on_cran()
	Y = bin_panel()
	fit = quiet(lame(
		Y,
		family = "binary",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	pdf(NULL)
	on.exit(dev.off(), add = TRUE)
	w = warns_of(m <- msgs_of(plot(fit, which = 1:4)))
	expect_true(any(grepl("GOF\\s+panel is skipped", m)))
	expect_false(any(grepl("geom_line", w)))
})

test_that("posterior_quantiles(component = 'UV') needs the saved position draws", {
	Y = bin_panel()
	fit = quiet(lame(
		Y,
		family = "binary",
		R = 2,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		posterior_opts = posterior_options()
	))
	expect_error(posterior_quantiles(fit, component = "UV"), "save_UV")
	fs = quiet(lame(
		Y,
		family = "binary",
		R = 2,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		posterior_opts = posterior_options(save_UV = TRUE)
	))
	pq = quiet(posterior_quantiles(fs, component = "UV"))
	expect_equal(dim(pq)[1:2], c(12L, 12L))
	f0 = quiet(lame(
		Y,
		family = "binary",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_error(posterior_quantiles(f0, component = "UV"), "R = 0")
	expect_message(latent_positions(f0), "no latent positions")
})

test_that("procrustes_align() keeps U V' on a directed fit", {
	Y = bin_panel()
	fit = quiet(lame(
		Y,
		family = "binary",
		R = 2,
		dynamic_uv = TRUE,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	al = quiet(procrustes_align(fit))
	before = sapply(1:3, function(t) fit$U[, , t] %*% t(fit[["V"]][, , t]))
	after = sapply(1:3, function(t) al$U[, , t] %*% t(al[["V"]][, , t]))
	expect_lt(max(abs(before - after)), 1e-8)
})

test_that("ALS Poisson draws do not warn on the network diagonal", {
	set.seed(6886)
	n = 12; nm = sprintf("a%02d", seq_len(n))
	Y = matrix(rpois(n * n, 1.5), n, n); diag(Y) = NA; dimnames(Y) = list(nm, nm)
	fit = quiet(ame_als(Y, family = "poisson", R = 1))
	expect_length(warns_of(simulate(fit, nsim = 2)), 0L)
})

test_that("dynamic_beta_prior_summary() rejects a non-numeric n_periods", {
	expect_error(dynamic_beta_prior_summary(list(a = 1)), "single number")
})

test_that("a handful of mutual dyads in total gets the weak-identification message", {
	set.seed(6886)
	n = 30; nm = sprintf("b%02d", seq_len(n))
	Y = lapply(1:4, function(t) {
		m = matrix(0, n, n)
		ev = cbind(sample(n, 25, TRUE), sample(n, 25, TRUE)); ev = ev[ev[, 1] != ev[, 2], , drop = FALSE]
		m[ev] = 1
		m[ev[1:3, 2:1, drop = FALSE]] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	m = msgs_of(lame(
		Y,
		family = "binary",
		R = 0,
		dcor = TRUE,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_true(any(grepl("mutual dyads in total", m)))
	# the absence of mutual ties informs rho, so the note does not say its posterior stays at the prior
	expect_false(any(grepl("stays near its prior", m)))
	expect_true(any(grepl("move well away from it", m)))
})

test_that("the pre-fit Not-estimated note reads the rank and names the supplied prior", {
	# one event per actor-period: under the intercept rate of a rank-one fit (1.6), above that of a rank-zero fit (0.6)
	set.seed(6886)
	n = 30; nm = sprintf("b%02d", seq_len(n))
	Y = lapply(1:4, function(t) {
		m = matrix(0, n, n)
		m[cbind(seq_len(n), c(2:n, 1))] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	run = function(R, prior = list()) msgs_of(try(
		lame(
			Y,
			family = "binary",
			R = R,
			prior = prior,
			dcor = FALSE,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		silent = TRUE
	))
	m1 = run(1, list(eta0 = 40))
	expect_true(any(grepl("The intercept is not estimated.*below 1.6 for a fit with a multiplicative term", m1)))
	expect_true(any(grepl("reported under the supplied prior (eta0 = 40)", m1, fixed = TRUE)))
	expect_false(any(grepl("(the supplied prior)", m1, fixed = TRUE)))
	m0 = run(0)
	expect_false(any(grepl("The intercept is not estimated", m0)))
})

test_that("bipartite G is a posterior mean and V positions carry a posterior sd", {
	set.seed(6886)
	rn = sprintf("r%02d", 1:12); cn = sprintf("c%02d", 1:9)
	Yb = lapply(1:3, function(t) matrix(rbinom(108, 1, 0.3), 12, 9, dimnames = list(rn, cn)))
	fs = quiet(lame(
		Yb,
		family = "binary",
		mode = "bipartite",
		R_row = 2,
		R_col = 2,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(dim(fs[["G"]]), c(2L, 2L))
	expect_true(all(is.finite(fs[["G"]])))
	fd = quiet(lame(
		Yb,
		family = "binary",
		mode = "bipartite",
		R_row = 2,
		R_col = 2,
		dynamic_uv = TRUE,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_equal(dim(fd[["G"]]), c(2L, 2L))
	expect_true(all(is.finite(fd[["G"]])))
	lp = latent_positions(fd)
	expect_true(all(is.finite(lp$posterior_sd)))
	expect_true(all(c("U", "V") %in% lp$type))
	# a per-period interaction matrix is refused, and nothing in the reporting
	# offers one
	expect_error(
		lame(
			Yb,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			dynamic_uv = TRUE,
			dynamic_G = TRUE,
			nscan = 40,
			burn = 10,
			odens = 2,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		"not supported"
	)
	out = paste(c(capture.output(print(fd)), capture.output(print(fd), type = "message")), collapse = "\n")
	expect_false(grepl("G_t|rho_G", out))
})

test_that("frequent rejected factor conditionals are stated as a property in print and summary", {
	skip_on_cran()
	Y = bin_panel(n = 16)
	fit = quiet(lame(
		Y,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		prior = list(uv_max_abs = 1e-4),
		burn = 20,
		nscan = 60,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_gt(fit$refusal_share[["UV"]], 0.5)
	out = paste(c(capture.output(print(fit)), capture.output(print(fit), type = "message")), collapse = "\n")
	expect_true(grepl("prior box", out))
	expect_false(grepl("sensitivity to the prior bounds|before interpreting", out))
	so = paste(capture.output(print(summary(fit))), collapse = "\n")
	expect_true(grepl("prior box", so))
	# a fit that rejects few proposals states nothing
	expect_null(lame:::.lame_bound_share_line(list(
		sampler_health = data.frame(step = "UV", numerical_failures = 0, prior_bound_rejections = 1),
		sampler_health_iterations = 1000
	)))
})

test_that("a 4-D Xdyad array is accepted as the per-period list it stands for", {
	Y = bin_panel()
	set.seed(3)
	Xl = lapply(Y, function(m) array(rnorm(length(m)), c(dim(m), 1), dimnames = list(rownames(m), colnames(m), "d")))
	X4 = array(unlist(Xl), c(dim(Y[[1]]), 1, length(Y)), dimnames = list(rownames(Y[[1]]), colnames(Y[[1]]), "d", NULL))
	f1 = quiet(
		lame(
			Y,
			Xdyad = Xl,
			family = "binary",
			R = 0,
			nscan = 20,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	f2 = quiet(
		lame(
			Y,
			Xdyad = X4,
			family = "binary",
			R = 0,
			nscan = 20,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_equal(f1$BETA, f2$BETA)
})

test_that("large covariates fit on the MCMC route with a note, and the ALS route refuses them", {
	Y = bin_panel()
	set.seed(4)
	Xbig = lapply(
		Y,
		function(m) array(rnorm(length(m), 0, 1e5), c(dim(m), 1), dimnames = list(rownames(m), colnames(m), "gdp"))
	)
	m = msgs_of(
		f <- lame(
			Y,
			Xdyad = Xbig,
			family = "binary",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_true(any(grepl("does not depend on the units", m)))
	expect_true(all(is.finite(f$BETA)))
	Ya = Y[[1]]
	expect_true(
		all(
			is.finite(
				quiet(
					ame(
						Ya,
						Xdyad = Xbig[[1]],
						family = "binary",
						R = 0,
						nscan = 10,
						burn = 2,
						odens = 1,
						verbose = FALSE,
						plot = FALSE,
						gof = FALSE
					)
				)$BETA
			)
		)
	)
	expect_error(
		quiet(lame(Y, Xdyad = Xbig, family = "binary", R = 0, method = "als", verbose = FALSE)),
		class = "lame_als_covariate_scale"
	)
})

test_that("the bundled addhealthc9 covariates align with its network", {
	data(addhealthc9)
	expect_equal(nrow(addhealthc9$X), nrow(addhealthc9$Y))
	expect_equal(rownames(addhealthc9$X), rownames(addhealthc9$Y))
})

test_that("UVPM marks actors absent from a period as NA under changing composition", {
	set.seed(6886)
	n = 10; TT = 3; nm = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		keep = nm[-t]
		m = matrix(rnorm(length(keep)^2), length(keep), length(keep), dimnames = list(keep, keep)); diag(m) = NA; m
	})
	fit = quiet(
		lame(
			Y,
			family = "normal",
			R = 2,
			dynamic_uv = TRUE,
			nscan = 40,
			burn = 10,
			odens = 2,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)
	)
	expect_equal(dim(fit$UVPM), c(n, n, TT))
	for (t in seq_len(TT)) {
		expect_true(all(is.na(fit$UVPM[nm[t], , t])))
		expect_true(all(is.na(fit$UVPM[, nm[t], t])))
		present = nm[-t]
		expect_true(all(is.finite(fit$UVPM[present, present, t])))
		expect_equal(rownames(fit$YPM[[t]]), present)
	}
})
