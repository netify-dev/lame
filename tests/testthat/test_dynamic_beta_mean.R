# long-run means of ar1 coefficient paths (R/dynamic_beta_mean.R): an ar1 path reverts to its own mean mu_k, with a
# stationary first state and the static coefficient prior on mu_k. checked here: the exactness of every conditional
# under the centred prior, recovery of mu_k, forecast reversion, the reports, and the ALS counterpart
skip_on_cran()

# exactness

# one step from the joint law of (mu, sigma, rho, path, data): mu_k ~ N(0, vmu), sigma_g^2 ~ IG(shape, scale), rho_g
# truncated normal, each path stationary around its mean, and data Z_t = X_t beta_t + N(0, s2) on an n x n grid.
# returns the z scores of the paired moment differences before and after the step
mean_geweke = function(step, ndraw, seed, gaps = rep(1, 4), pool = FALSE) {
	set.seed(seed)
	ns = asNamespace("lame")
	n = 4L; p = 2L; nt = length(gaps) + 1L
	group_id = c(1L, 1L); lam = c(0.8, 1.4); vmu = 2; s2 = 0.7
	shape = 3; scale = 0.4; rho_mean = 0.6; rho_sd = 0.25; lo = 0; hi = 0.999; pmu = 0.3; ptau = 0.8
	X = lapply(seq_len(nt), function(t) matrix(rnorm(n * n * p), n * n, p))
	rtn = function() repeat { x = rnorm(1, rho_mean, rho_sd); if (x > lo && x < hi) return(x) }
	stat_of = function(s) c(
		s$mu, s$mu^2, s$sig, s$sig^2, s$rho, s$rho^2, s$rho * s$sig,
		as.vector(s$path[c(1, nt), ]), as.vector(s$path[c(1, nt), ])^2,
		s$mu[1] * s$path[nt, 1], s$rho * s$path[1, 2],
		sum(vapply(seq_len(nt), function(t) sum((as.vector(s$z[[t]]) - X[[t]] %*% s$path[t, ])^2, na.rm = TRUE), 0))
	)
	D = vector("list", ndraw)
	for (k in seq_len(ndraw)) {
		mu = rnorm(p, 0, sqrt(vmu))
		sig = sqrt(1 / rgamma(1, shape, scale))
		rho = if (pool) lo + (hi - lo) * plogis(rnorm(1, pmu, ptau)) else rtn()
		path = matrix(0, nt, p)
		z0 = ns$.lame_beta_gap_transition(rho, gaps)
		for (j in seq_len(p)) {
			d = rnorm(1, 0, sig * sqrt(lam[j] / (1 - rho^2)))
			path[1, j] = d
			for (t in 2:nt) path[t, j] = z0$f[t - 1] * path[t - 1, j] + sig * sqrt(lam[j] * z0$q[t - 1]) * rnorm(1)
		}
		path = sweep(path, 2, mu, "+")
		z = lapply(seq_len(nt), function(t) {
			zt = matrix(X[[t]] %*% path[t, ] + rnorm(n * n, 0, sqrt(s2)), n, n)
			diag(zt) = NA
			zt
		})
		s = list(mu = mu, sig = sig, rho = rho, path = path, z = z)
		ctx = list(
			X = X, n = n, p = p, nt = nt, gaps = gaps, group_id = group_id, lam = lam, vmu = vmu, s2 = s2,
			shape = shape, scale = scale, rho_mean = rho_mean, rho_sd = rho_sd, lo = lo, hi = hi, pmu = pmu, ptau = ptau
		)
		s1 = step(s, ctx)
		D[[k]] = stat_of(s1) - stat_of(s)
	}
	D = do.call(rbind, D)
	sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, colMeans(D) / (sdv / sqrt(nrow(D))), 0)
	list(maxz = max(abs(z)), moved = mean(rowSums(abs(D) > 1e-12) > 0))
}

# the observation statistics of the test model, with the diagonal left out as in the sampler
mean_stats = function(s, ctx) {
	ns = asNamespace("lame")
	h = matrix(0, ctx$nt, ctx$p)
	Om = array(0, c(ctx$p, ctx$p, ctx$nt))
	for (t in seq_len(ctx$nt)) {
		zz = ns$.lame_beta_observation_stats(ctx$X[[t]], s$z[[t]], ctx$s2, 0, FALSE, FALSE)
		Om[, , t] = zz$Q
		h[t, ] = zz$h
	}
	list(Omega = Om, h = h)
}

test_that("the conditional of the long-run means keeps the joint law", {
	step = function(s, ctx) {
		s$mu = lame:::.lame_beta_mean_draw(
			s$path, s$mu, c(TRUE, TRUE), rep(s$rho, 2), rep(s$sig, 2), ctx$lam, rep(ctx$vmu, 2), ctx$gaps
		)
		s
	}
	r = mean_geweke(step, 3000, 31)
	expect_gt(r$moved, 0.99)
	expect_lt(r$maxz, 4.5)
	r = mean_geweke(step, 3000, 32, gaps = c(1, 2, 1, 3))
	expect_lt(r$maxz, 4.5)
})

test_that("the persistence and scale updates with a stationary first state keep the joint law", {
	rho_step = function(s, ctx) {
		dev = sweep(s$path, 2, s$mu)
		s$rho = lame:::.lame_sample_rho_beta_stationary(
			dev, ctx$group_id, 1L, ctx$lam, rep(s$sig, 2), s$rho, ctx$rho_mean, ctx$rho_sd, ctx$lo, ctx$hi, TRUE
		)
		s
	}
	r = mean_geweke(rho_step, 3000, 41)
	expect_gt(r$moved, 0.3)
	expect_lt(r$maxz, 4.5)
	gap_step = function(s, ctx) {
		dev = sweep(s$path, 2, s$mu)
		s$rho = lame:::.lame_sample_rho_beta_gaps(
			dev, ctx$group_id, 1L, diag(1 / ctx$lam), rep(s$sig, 2), s$rho, ctx$rho_mean, ctx$rho_sd, ctx$lo, ctx$hi,
			c(0, cumsum(ctx$gaps)), stationary = TRUE
		)
		s
	}
	r = mean_geweke(gap_step, 3000, 42, gaps = c(1, 2, 1, 3))
	expect_lt(r$maxz, 4.5)
	pool_step = function(s, ctx) {
		dev = sweep(s$path, 2, s$mu)
		s$rho = lame:::.pool_sample_rho(
			dev, ctx$group_id, 1L, diag(ctx$lam), rep(s$sig, 2), s$rho, ctx$pmu, ctx$ptau, ctx$lo, ctx$hi,
			c(0, cumsum(ctx$gaps)), stationary = TRUE
		)
		s
	}
	r = mean_geweke(pool_step, 3000, 43, pool = TRUE)
	expect_lt(r$maxz, 4.5)
	sig_step = function(s, ctx) {
		dev = sweep(s$path, 2, s$mu)
		s$sig = lame:::.sample_sigma_beta_alt(
			dev, "ar1", prior_shape = ctx$shape, prior_scale = ctx$scale, Lambda = diag(ctx$lam), group_id = ctx$group_id,
			n_groups = 1L, time_positions = c(0, cumsum(ctx$gaps)), rho_by_group = s$rho, stationary = TRUE
		)
		s
	}
	r = mean_geweke(sig_step, 3000, 44)
	expect_lt(r$maxz, 4.5)
	r = mean_geweke(sig_step, 3000, 45, gaps = c(1, 2, 1, 3))
	expect_lt(r$maxz, 4.5)
})

test_that("the exactness check fails for a persistence step without its first-state correction", {
	src = deparse(lame:::.lame_sample_rho_beta_stationary)
	src = sub("log(stats::runif(1L)) < lr", "TRUE", src, fixed = TRUE)
	broken = eval(parse(text = paste(src, collapse = "\n")), envir = asNamespace("lame"))
	step = function(s, ctx) {
		dev = sweep(s$path, 2, s$mu)
		s$rho = broken(dev, ctx$group_id, 1L, ctx$lam, rep(s$sig, 2), s$rho, ctx$rho_mean, ctx$rho_sd, ctx$lo, ctx$hi, TRUE)
		s
	}
	expect_gt(mean_geweke(step, 3000, 41)$maxz, 6)
})

test_that("the path draw given the long-run means keeps the joint law", {
	step = function(s, ctx) {
		off = lapply(seq_len(ctx$nt), function(t) matrix(0, ctx$n, ctx$n))
		mu_off = lame:::.lame_beta_mean_offsets(off, ctx$X, s$mu)
		init = lame:::.lame_beta_init_cov(diag(10, 2), c(TRUE, TRUE), rep(s$rho, 2), rep(s$sig, 2), ctx$lam)
		out = lame:::sample_beta_dynamic_cpp(
			Xdyn_list = ctx$X,
			Xstat_list = lapply(ctx$X, function(x) matrix(0, nrow(x), 0)),
			Z_list = s$z,
			offset_list = mu_off,
			beta_static = numeric(0),
			rho_by_coef = rep(s$rho, 2),
			sigma_by_coef = rep(s$sig, 2),
			Lambda = diag(ctx$lam),
			beta0_mean = c(0, 0),
			beta0_cov = init,
			s2 = ctx$s2,
			dyad_rho = 0,
			bipartite = FALSE,
			symmetric = FALSE,
			use_dyad_rho = FALSE
		)
		s$path = sweep(out$path, 2, s$mu, "+")
		s
	}
	r = mean_geweke(step, 2500, 51)
	expect_lt(r$maxz, 4.5)
	# the precomputed statistics, moved by Omega mu, give the same draw as the offsets
	set.seed(52)
	s = list(z = lapply(1:5, function(t) matrix(rnorm(16), 4)), mu = c(0.7, -0.4))
	X = lapply(1:5, function(t) matrix(rnorm(32), 16, 2))
	ctx = list(X = X, nt = 5L, p = 2L, s2 = 0.7)
	off = lapply(1:5, function(t) matrix(0, 4, 4))
	base = list(
		Xdyn_list = X, Xstat_list = lapply(X, function(x) matrix(0, 16, 0)), Z_list = s$z, beta_static = numeric(0),
		rho_by_coef = c(0.6, 0.6), sigma_by_coef = c(0.3, 0.3), Lambda = diag(2), beta0_mean = c(0, 0),
		beta0_cov = diag(0.25, 2), s2 = 0.7, dyad_rho = 0, bipartite = FALSE, symmetric = FALSE, use_dyad_rho = FALSE
	)
	st = lame:::.lame_beta_mean_stats(mean_stats(s, ctx), s$mu)
	set.seed(9)
	a = do.call(lame:::sample_beta_dynamic_cpp, c(base, list(offset_list = lame:::.lame_beta_mean_offsets(off, X, s$mu))))
	set.seed(9)
	b = do.call(lame:::sample_beta_dynamic_cpp, c(base, list(offset_list = off, Omega_in = st$Omega, h_in = st$h)))
	expect_equal(a$path, b$path, tolerance = 1e-10)
})

test_that("the interweaving step with long-run means keeps the joint law", {
	step = function(s, ctx) {
		st = mean_stats(s, ctx)
		out = lame:::.lame_beta_nc_step(
			s$path, st, ctx$group_id, 1L, "ar1", ctx$lam, s$sig, s$rho, ctx$shape, ctx$scale,
			ctx$rho_mean, ctx$rho_sd, ctx$lo, ctx$hi, gaps = ctx$gaps, mu = s$mu
		)
		s$path = out$path; s$sig = out$sigma; s$rho = out$rho
		s
	}
	r = mean_geweke(step, 3000, 61)
	expect_gt(r$moved, 0.05)
	expect_lt(r$maxz, 4.5)
	# the map of an ar1 deviation path with a stationary first state inverts itself
	x = matrix(rnorm(12), 6, 2)
	moved = lame:::.lame_beta_nc_path(x, "ar1", 0.3, 0.7, 0.5, 0.4, c(1, 2), rep(1, 5), TRUE)
	back = lame:::.lame_beta_nc_path(moved, "ar1", 0.5, 0.4, 0.3, 0.7, c(1, 2), rep(1, 5), TRUE)
	expect_equal(back, x, tolerance = 1e-10)
	expect_equal(moved[1, ] / x[1, ], rep(0.5 / 0.3 * sqrt((1 - 0.49) / (1 - 0.16)), 2), tolerance = 1e-12)
})

# fits

# a directed panel whose dyadic slope moves around 0.8 and whose second slope is constant at 0.5
mean_panel = function(seed, n = 15, nt = 6, family = "normal", rho = 0.7, sd = 0.2) {
	set.seed(seed)
	d = numeric(nt)
	d[1] = rnorm(1, 0, sd / sqrt(1 - rho^2))
	for (t in 2:nt) d[t] = rho * d[t - 1] + rnorm(1, 0, sd)
	slope = 0.8 + d
	a = rnorm(n, 0, 0.3)
	rn = sprintf("a%02d", seq_len(n))
	Y = X = vector("list", nt)
	for (t in seq_len(nt)) {
		x1 = matrix(rnorm(n * n), n, n)
		x2 = matrix(rnorm(n * n), n, n)
		eta = 0.2 + slope[t] * x1 + 0.5 * x2 + outer(a, a, "+")
		y = if (family == "binary") 1 * (eta + rnorm(n * n) > 0) else eta + matrix(rnorm(n * n), n, n)
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		Y[[t]] = y
		X[[t]] = array(c(x1, x2), c(n, n, 2), dimnames = list(rn, rn, c("x1", "x2")))
	}
	names(Y) = names(X) = paste0("t", seq_len(nt))
	list(Y = Y, X = X, slope = slope)
}

mean_fit = function(d, nt = 4, ...) {
	suppressWarnings(suppressMessages(lame(
		d$Y[seq_len(nt)],
		Xdyad = d$X[seq_len(nt)],
		R = 0,
		dynamic_beta = "dyad",
		burn = 200,
		nscan = 1200,
		odens = 5,
		chains = 1,
		verbose = FALSE,
		...
	)))
}

test_that("a slope path is recovered around its long-run mean, which is reported and monitored", {
	d = mean_panel(71)
	fit = mean_fit(d, family = "normal")
	expect_true(is.matrix(fit$MU_BETA))
	expect_identical(colnames(fit$MU_BETA), c("x1_dyad", "x2_dyad"))
	expect_equal(nrow(fit$MU_BETA), nrow(fit$BETA))
	expect_equal(fit$mu_beta, colMeans(fit$MU_BETA))
	# the constant slope's long-run mean sits at its value, the moving slope's path at its truth
	q = quantile(fit$MU_BETA[, "x2_dyad"], c(0.025, 0.975))
	expect_true(q[1] < 0.5 && q[2] > 0.5)
	path = colMeans(fit$BETA[, "x1_dyad", ])
	expect_lt(max(abs(path - d$slope[1:4]) / apply(fit$BETA[, "x1_dyad", ], 2, sd)), 3)
	# reports
	out = paste(capture.output(print(fit)), collapse = "\n")
	expect_match(out, "Long-run mean mu_beta[x1_dyad]", fixed = TRUE)
	sm = summary(fit)
	expect_identical(rownames(sm$mu_beta), c("x1_dyad", "x2_dyad"))
	expect_match(paste(capture.output(print(sm)), collapse = "\n"), "mu_beta = x1_dyad", fixed = TRUE)
	ci = confint(fit, parm = "mu_beta")
	expect_identical(rownames(ci), c("mu_beta[x1_dyad]", "mu_beta[x2_dyad]"))
	expect_equal(unname(ci[2, ]), unname(quantile(fit$MU_BETA[, "x2_dyad"], c(0.025, 0.975))))
	expect_identical(rownames(confint(fit, parm = c("mu_beta[x2_dyad]", "ve"))), c("mu_beta[x2_dyad]", "ve"))
	expect_true(all(c("mu_beta[x1_dyad]", "mu_beta[x2_dyad]") %in% colnames(lame:::.lame_gate_series(fit))))
	expect_true(all(c("mu_beta[x1_dyad]", "mu_beta[x2_dyad]") %in% names(lame:::.lame_chain_series(fit))))
	tp = trace_plot(fit, params = "dynamic")
	tp_data = if (inherits(tp, "patchwork")) tp[[1]]$data else tp$data
	expect_true(any(grepl("long-run mean", unique(as.character(tp_data$parameter)), fixed = TRUE)))
})

test_that("coefficient forecasts revert to the long-run mean", {
	d = mean_panel(72)
	fit = mean_fit(d, family = "normal")
	# with no innovations the forecast is the AR(1) mean recursion around mu
	fit$SIGMA_BETA[] = 0
	n = nrow(d$Y[[1]])
	rn = rownames(d$Y[[1]])
	nd = function(v1) lapply(1:3, function(k) array(c(rep(v1, n * n), rep(0, n * n)), c(n, n, 2), dimnames = list(rn, rn, c("x1", "x2"))))
	p1 = suppressWarnings(predict(fit, h = 3, type = "link", newdata = nd(1), seed = 1))
	p0 = suppressWarnings(predict(fit, h = 3, type = "link", newdata = nd(0), seed = 1))
	fc = sapply(1:3, function(k) mean(p1[[k]] - p0[[k]], na.rm = TRUE))
	mu = fit$MU_BETA[, "x1_dyad"]
	rho = fit$RHO_BETA[, "dyad"]
	last = fit$BETA[, "x1_dyad", 4]
	expect_equal(fc, sapply(1:3, function(h) mean(mu + rho^h * (last - mu))), tolerance = 1e-8)
	# the forecast of a long horizon is the long-run mean, not zero
	p_far = suppressWarnings(predict(fit, h = 40, type = "link", newdata = c(nd(1), rep(nd(1)[1], 37)), seed = 1))
	p_far0 = suppressWarnings(predict(fit, h = 40, type = "link", newdata = c(nd(0), rep(nd(0)[1], 37)), seed = 1))
	expect_equal(mean(p_far[[40]] - p_far0[[40]], na.rm = TRUE), mean(mu + rho^40 * (last - mu)), tolerance = 1e-8)
	expect_gt(mean(p_far[[40]] - p_far0[[40]], na.rm = TRUE), 0.3)
})

test_that("random-walk and static fits carry no long-run means", {
	d = mean_panel(73)
	f_rw = mean_fit(d, family = "normal", dynamic_beta_kind = "rw1")
	expect_null(f_rw$MU_BETA)
	expect_null(f_rw$mu_beta)
	expect_identical(nrow(confint(f_rw, parm = "mu_beta")), 0L)
	f_ab = suppressWarnings(suppressMessages(lame(
		d$Y[1:4], Xdyad = d$X[1:4], R = 0, dynamic_ab = TRUE, burn = 50, nscan = 200, odens = 5, chains = 1, verbose = FALSE
	)))
	expect_identical(unname(f_ab$dynamic_beta_kind_by_block), "rw1")
	expect_null(f_ab$MU_BETA)
	# a fit with a dynamic_ab intercept (rw1) and ar1 slopes centres the slopes only
	f_mix = mean_fit(d, family = "normal", dynamic_ab = TRUE)
	expect_identical(colnames(f_mix$MU_BETA), c("x1_dyad", "x2_dyad"))
	expect_true(all(f_mix$RHO_BETA[, "intercept"] == 1))
})

test_that("long-run means are sampled in every mode and family that offers ar1 paths", {
	run = function(...) suppressWarnings(suppressMessages(lame(
		..., burn = 30, nscan = 120, odens = 5, chains = 1, verbose = FALSE
	)))
	d = mean_panel(74, n = 10, nt = 4)
	f_bin = run(d$Y, Xdyad = d$X, R = 1, family = "binary", dynamic_beta = TRUE)
	expect_identical(colnames(f_bin$MU_BETA), c("intercept", "x1_dyad", "x2_dyad"))
	expect_true(all(is.finite(f_bin$MU_BETA)))
	f_pois = run(lapply(d$Y, function(y) round(exp(pmin(y, 3)))), Xdyad = d$X, R = 0, family = "poisson", dynamic_beta = "dyad")
	expect_true(all(is.finite(f_pois$MU_BETA)))
	f_ord = run(lapply(d$Y, function(y) { z = y; z[] = findInterval(y, c(-1, 0, 1)); z }), Xdyad = d$X, R = 0, family = "ordinal",
		dynamic_beta = "dyad")
	expect_true(all(is.finite(f_ord$MU_BETA)))
	# symmetric binary with a tied nodal covariate
	sim = simulate_test_network(n = 10, n_time = 4, family = "binary", symmetric = TRUE, seed = 75)
	xr = lapply(seq_along(sim$Y), function(t) matrix(rnorm(10), 10, 1, dimnames = list(rownames(sim$Y[[t]]), "xn")))
	names(xr) = names(sim$Y)
	f_sym = run(sim$Y, Xrow = xr, R = 0, family = "binary", symmetric = TRUE, dynamic_beta = "row")
	expect_true(is.matrix(f_sym$MU_BETA) && all(is.finite(f_sym$MU_BETA)))
	# bipartite normal with a row covariate path
	bp = simulate_test_network(n = 8, n_time = 4, family = "normal", mode = "bipartite", nA = 8, nB = 6, seed = 76)
	f_bip = run(bp$Y, Xdyad = bp$Xdyad, mode = "bipartite", R_row = 1, R_col = 1, family = "normal", dynamic_beta = "dyad")
	expect_true(is.matrix(f_bip$MU_BETA) && all(is.finite(f_bip$MU_BETA)))
	# unequal time gaps and a changing actor set
	f_gap = run(d$Y, Xdyad = d$X, R = 0, family = "normal", dynamic_beta = "dyad", time_index = c(1, 2, 4, 5))
	expect_true(all(is.finite(f_gap$MU_BETA)))
	yc = d$Y; xc = d$X
	yc[[2]] = yc[[2]][-1, -1]; xc[[2]] = xc[[2]][-1, -1, , drop = FALSE]
	f_chg = run(yc, Xdyad = xc, R = 0, family = "normal", dynamic_beta = "dyad")
	expect_true(all(is.finite(f_chg$MU_BETA)))
	# pooled persistences across blocks
	f_pool = run(d$Y, Xdyad = d$X, R = 0, family = "normal", dynamic_beta = TRUE, dynamic_beta_pool = "both")
	expect_true(all(is.finite(f_pool$MU_BETA)))
})

test_that("the ALS estimator fits the same long-run means", {
	d = mean_panel(77, n = 20)
	fit = mean_fit(d, family = "normal")
	als = suppressWarnings(suppressMessages(lame(
		d$Y[1:4], Xdyad = d$X[1:4], R = 0, dynamic_beta = "dyad", method = "als", verbose = FALSE
	)))
	expect_identical(names(als$mu_beta), c("x1_dyad", "x2_dyad"))
	expect_true(all(is.finite(als$eb$mu_beta_se)))
	gap = (als$coef_path["x1_dyad", ] - colMeans(fit$BETA[, "x1_dyad", ])) / apply(fit$BETA[, "x1_dyad", ], 2, sd)
	expect_lt(max(abs(gap)), 1)
	expect_lt(abs(als$mu_beta[["x2_dyad"]] - fit$mu_beta[["x2_dyad"]]) / sd(fit$MU_BETA[, "x2_dyad"]), 1)
	# the penalised fit, as the sampler's warm start asks for it, profiles the mean out of the path penalty
	pen = suppressWarnings(suppressMessages(lame(
		d$Y[1:4], Xdyad = d$X[1:4], R = 0, dynamic_beta = "dyad", method = "als",
		prior = structure(list(), lame_als_estimator = "none"), verbose = FALSE
	)))
	expect_identical(names(pen$mu_beta), c("x1_dyad", "x2_dyad"))
	expect_lt(abs(pen$mu_beta[["x2_dyad"]] - 0.5), 0.15)
	K = lame:::.lda_beta_path_penalty(5, 2, "ar1", 0.6)
	expect_equal(drop(K %*% rep(1, 5)), rep(0, 5), tolerance = 1e-12)
})

# forecasts of rw2 and matern32 paths

test_that("rw2 and matern32 forecast maps are the conditionals of their priors", {
	nt = 5L
	h = 3L
	fake = function(kind) list(
		BETA = array(0, c(2, 1, nt), dimnames = list(NULL, "x1_dyad", paste0("t", seq_len(nt)))),
		beta_dynamic_mask = TRUE,
		beta_dynamic_groups = "dyad",
		RHO_BETA = matrix(1, 2, 1, dimnames = list(NULL, "dyad")),
		dynamic_beta_kind = kind,
		prior = list()
	)
	m = lame:::.lame_beta_smooth_forecast_map(fake("rw2"), h)$rw2
	# second differences: the mean extends the last slope, the variance grows as the sum of squares
	expect_equal(unname(m$A[, (nt - 1):nt]), cbind(-(1:h), 2:(h + 1)), tolerance = 1e-10)
	expect_equal(sum(abs(m$A[, 1:(nt - 2)])), 0, tolerance = 1e-10)
	expect_equal(diag(tcrossprod(m$L)), cumsum((1:h)^2), tolerance = 1e-10)
	# matern32: the change since the first period of a Matern 3/2 process, conditioned on the fitted changes
	m = lame:::.lame_beta_smooth_forecast_map(fake("matern32"), h)$matern32
	pos = seq_len(nt + h)
	ls = max(2, nt / 4)
	a = sqrt(3) * abs(outer(pos, pos, "-")) / ls
	S = (1 + a) * exp(-a)
	M = cbind(-1, diag(nt + h - 1))
	Sc = M %*% S %*% t(M)
	p_c = 1:(nt - 1)
	f_c = nt - 1 + 1:h
	Ac = Sc[f_c, p_c] %*% solve(Sc[p_c, p_c])
	# the forecast of beta is beta_1 plus the forecast change; in the path's coordinates
	A = cbind(1 - rowSums(Ac), Ac)
	expect_equal(unname(m$A), unname(A), tolerance = 1e-8)
	expect_equal(unname(tcrossprod(m$L)), unname(Sc[f_c, f_c] - Ac %*% Sc[p_c, f_c]), tolerance = 1e-8)
	# the forecast keeps the level: a constant path forecasts itself
	expect_equal(drop(m$A %*% rep(2.5, nt)), rep(2.5, h), tolerance = 1e-10)
	expect_null(lame:::.lame_beta_smooth_forecast_map(fake("rw1"), h))
})

test_that("rw2 and matern32 fits with dynamic_ab forecast at their conditional mean and spread", {
	d = mean_panel(78, n = 10, nt = 5)
	for (kind in c("rw2", "matern32")) {
		fit = suppressWarnings(suppressMessages(lame(
			d$Y, Xdyad = d$X, R = 0, dynamic_ab = TRUE, dynamic_beta = "dyad", dynamic_beta_kind = kind,
			burn = 30, nscan = 400, odens = 1, chains = 1, verbose = FALSE
		)))
		expect_identical(unname(fit$dynamic_beta_kind_by_block[c("intercept", "dyad")]), c(kind, kind))
		# one fitted path and scale in every draw, so the forecast draws are the conditional law
		fit$BETA[] = rep(fit$BETA[1, , ], each = dim(fit$BETA)[1])
		fit$SIGMA_BETA[] = 0.2
		n = nrow(d$Y[[1]])
		rn = rownames(d$Y[[1]])
		nd = function(v1) lapply(1:2, function(k) array(c(rep(v1, n * n), rep(0, n * n)), c(n, n, 2), dimnames = list(rn, rn, c("x1", "x2"))))
		p1 = suppressWarnings(predict(fit, h = 2, type = "link", newdata = nd(1), seed = 3, by_draw = TRUE))
		p0 = suppressWarnings(predict(fit, h = 2, type = "link", newdata = nd(0), seed = 3, by_draw = TRUE))
		fc = p1[2, 1, , ] - p0[2, 1, , ]
		m = lame:::.lame_beta_smooth_forecast_map(fit, 2)[[kind]]
		sdk = lame:::.lame_beta_transition_scale(fit, c("ar1", "rw1", "rw2", "matern32"))$sd[["x1_dyad"]]
		mean_cf = drop(m$A %*% fit$BETA[1, "x1_dyad", ])
		sd_cf = 0.2 * sdk * sqrt(diag(tcrossprod(m$L)))
		S = ncol(fc)
		expect_lt(max(abs(rowMeans(fc) - mean_cf) / (sd_cf / sqrt(S))), 4)
		expect_lt(max(abs(apply(fc, 1, sd) / sd_cf - 1)), 4 / sqrt(2 * S))
	}
})
