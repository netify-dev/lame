library(lame)
skip_on_cran()

# panels with known additive levels, an optional AR(1) deviation, and a gaussian or lognormal-poisson response
sim_als_vc_panel = function(
	n = 30,
	Tt = 5,
	family = "normal",
	seed = 1,
	sd_lev = 0.5,
	sd_dev = 0,
	rho_dev = 0.7,
	s2 = if (family == "poisson") 0.3 else 1,
	b0 = if (family == "poisson") 0.3 else 0.2,
	b1 = 0.7
) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	a_lev = rnorm(n, 0, sd_lev); b_lev = rnorm(n, 0, sd_lev)
	ar_path = function(lev) {
		out = matrix(lev, n, Tt)
		if (sd_dev > 0) {
			d = matrix(0, n, Tt)
			d[, 1] = rnorm(n, 0, sd_dev)
			if (Tt > 1) for (t in 2:Tt) {
				d[, t] = rho_dev * d[, t - 1] + rnorm(n, 0, sd_dev * sqrt(1 - rho_dev^2))
			}
			out = out + d
		}
		out
	}
	a_path = ar_path(a_lev); b_path = ar_path(b_lev)
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, n); diag(x) = NA
		eta = b0 + b1 * x + outer(a_path[, t], b_path[, t], "+")
		z = eta + matrix(rnorm(n * n, 0, sqrt(s2)), n, n)
		y = if (family == "poisson") {
			matrix(rpois(n * n, exp(pmin(ifelse(is.na(z), 0, z), 6))), n, n)
		} else z
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(nm, nm, "x1"))
	}
	list(
		Y = Y,
		X = X,
		a_lev = a_lev,
		b_lev = b_lev,
		s2 = s2,
		b0 = b0,
		b1 = b1,
		a_path = a_path,
		b_path = b_path
	)
}

test_that("poisson ALS reports the level and ve on the latent scale", {
	# the IRLS fit solves the marginal log-mean E[y] = exp(eta + ve / 2) of the lognormal-poisson model; the reported
	# level subtracts the Jensen term and ve is the latent variance from the count overdispersion
	icpt = ve = slope = numeric(0)
	for (s in 1:4) {
		d = sim_als_vc_panel(family = "poisson", seed = 100 + s)
		f = suppressWarnings(lame_als(
			d$Y,
			Xdyad = d$X,
			R = 0,
			family = "poisson",
			verbose = FALSE
		))
		icpt = c(icpt, as.numeric(f$coefficients["intercept"]))
		ve = c(ve, as.numeric(f$VC["ve"]))
		slope = c(slope, as.numeric(f$coefficients["x1_dyad"]))
		# the response scale adds the Jensen term back, so fitted values are the marginal count mean
		expect_equal(f$fitted[[1]], exp(f$EZ[[1]] + f$pois_latent_var / 2))
		expect_true(is.finite(f$ve_working) && f$ve_working > 0)
	}
	# the generating level is b0 plus a small identified additive/mean shift
	expect_lt(abs(mean(icpt) - 0.3), 0.075)
	expect_lt(abs(mean(ve) / 0.3 - 1), 0.3)
	expect_lt(abs(mean(slope) - 0.7), 0.035)
})

test_that("equidispersed counts leave the poisson level unshifted", {
	set.seed(7)
	n = 25
	x = matrix(rnorm(n * n), n, n)
	y = matrix(rpois(n * n, exp(0.3 + 0.5 * x)), n, n)
	diag(x) = NA; diag(y) = NA
	f = suppressWarnings(ame_als(
		y,
		Xdyad = array(x, c(n, n, 1)),
		R = 0,
		family = "poisson",
		bias_correct = "none",
		verbose = FALSE
	))
	# with no latent layer ve sits at or near zero and the level keeps the marginal log-mean
	expect_lt(as.numeric(f$VC["ve"]), 0.05)
	expect_lt(abs(as.numeric(f$coefficients["intercept"]) - 0.3), 0.15)
	# the empirical-Bayes default floors ve at the posterior mean of its inverse gamma prior, as the sampler does,
	# and reports the latent level
	fe = suppressWarnings(ame_als(y, Xdyad = array(x, c(n, n, 1)), R = 0, family = "poisson", verbose = FALSE))
	expect_lt(as.numeric(fe$VC["ve"]), 0.15)
	expect_lt(abs(as.numeric(fe$coefficients["intercept"]) - 0.19), 0.1)
})

test_that("static ALS va / vb subtract the additive sampling variance", {
	d = sim_als_vc_panel(seed = 3)
	f = suppressWarnings(lame_als(
		d$Y,
		Xdyad = d$X,
		R = 0,
		family = "normal",
		bias_correct = "none",
		verbose = FALSE
	))
	expect_true(all(is.finite(f$vc_raw)))
	expect_true(all(is.finite(f$vc_noise)))
	expect_equal(unname(f$VC["va"]), max(unname(f$vc_raw["va"] - f$vc_noise["va"]), 0))
	expect_equal(unname(f$VC["vb"]), max(unname(f$vc_raw["vb"] - f$vc_noise["vb"]), 0))
	# the correction moves the report toward the variance of the true levels
	expect_lt(
		abs(unname(f$VC["va"]) - var(d$a_lev)),
		abs(unname(f$vc_raw["va"]) - var(d$a_lev))
	)
	# and is small relative to the variance it corrects
	expect_lt(unname(f$vc_noise["va"]), 0.5 * unname(f$vc_raw["va"]))
})

test_that("bipartite ALS corrects the short side more than the long side", {
	# a row effect from 20 columns carries more error than a column effect from 30 rows, so its correction is larger
	set.seed(5)
	nA = 30; nB = 20; Tt = 4
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(nA * nB), nA, nB)
		y = 0.2 + 0.7 * x + outer(a, b, "+") + matrix(rnorm(nA * nB), nA, nB)
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(nA, nB, 1), dimnames = list(rn, cn, "x1"))
	}
	f = suppressWarnings(lame_als(
		Y,
		Xdyad = X,
		R = 0,
		family = "normal",
		bias_correct = "none",
		mode = "bipartite",
		verbose = FALSE
	))
	expect_gt(unname(f$vc_noise["va"]), unname(f$vc_noise["vb"]))
	expect_lt(abs(unname(f$VC["va"]) - var(a)), abs(unname(f$vc_raw["va"]) - var(a)))
})

test_that("dynamic_ab ALS reports the level variance, not the per-period one", {
	d = sim_als_vc_panel(seed = 11, sd_dev = 0.3)
	f = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		R = 0,
		family = "normal",
		dynamic_ab = TRUE,
		method = "als",
		bias_correct = "none",
		verbose = FALSE,
		plot = FALSE
	))
	expect_true(all(c("vc_period", "vc_raw") %in% names(f)))
	# the per-period variance includes the AR(1) deviation, so it is larger, and the reported va targets the level
	expect_gt(unname(f$vc_period["va"]), unname(f$VC["va"]))
	expect_lt(
		abs(unname(f$VC["va"]) - var(d$a_lev)),
		abs(unname(f$vc_period["va"]) - var(d$a_lev))
	)
	expect_gte(unname(f$VC["va"]), 0)
})

test_that("normal dynamic ALS reports ve on the residual degrees of freedom", {
	# penalized paths count by their effective df, so ve matches the realised noise variance
	ratio = numeric(0)
	for (s in 1:3) {
		d = sim_als_vc_panel(seed = 40 + s, sd_dev = 0.3)
		f = suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			R = 0,
			family = "normal",
			dynamic_ab = TRUE,
			method = "als",
			verbose = FALSE,
			plot = FALSE
		))
		ratio = c(ratio, unname(f$VC["ve"]) / d$s2)
	}
	expect_lt(abs(mean(ratio) - 1), 0.06)
})

test_that("t-kind dynamic ALS counts its position paths in the ve degrees of freedom", {
	# heavy-tailed paths enter the residual df with their fitted weights and the curvature of the t penalty
	cell = grid_cell(
		"tve", "normal", "lame_dynamic_als", R = 2,
		dynamic_uv = TRUE, dynamic_uv_kind = "t", T = 5
	)
	ratio = numeric(0)
	for (s in 1:3) {
		sim = grid_sim(cell, s)
		f = suppressWarnings(lame(
			sim$Y,
			sim$X,
			R = 2,
			family = "normal",
			method = "als",
			dynamic_uv = TRUE,
			dynamic_uv_kind = "t",
			verbose = FALSE,
			seed = s + 11L
		))
		e = unlist(lapply(seq_along(sim$E), function(t) {
			sim$E[[t]][is.finite(sim$Y[[t]]) & is.finite(sim$E[[t]])]
		}))
		ratio = c(ratio, unname(f$VC["ve"]) / var(e))
	}
	expect_gt(mean(ratio), 0.95)
	expect_lt(mean(ratio), 1.03)
})

test_that("symmetric dynamic ALS counts its position paths in the ve degrees of freedom", {
	# one path per actor, counted by the curvature of the symmetric solve
	for (kind in c("ar1", "t")) {
		cell = grid_cell(
			"sve",
			"normal",
			"lame_dynamic_als",
			symmetric = TRUE,
			R = 2,
			dynamic_uv = TRUE,
			dynamic_uv_kind = kind,
			T = 5
		)
		ratio = numeric(0)
		for (s in 1:3) {
			sim = grid_sim(cell, s)
			f = suppressWarnings(lame(
				sim$Y,
				sim$X,
				R = 2,
				family = "normal",
				symmetric = TRUE,
				method = "als",
				dynamic_uv = TRUE,
				dynamic_uv_kind = kind,
				verbose = FALSE,
				seed = s + 11L
			))
			e = unlist(lapply(seq_along(sim$E), function(t) {
				sim$E[[t]][is.finite(sim$Y[[t]]) & is.finite(sim$E[[t]])]
			}))
			ratio = c(ratio, unname(f$VC["ve"]) / var(e))
		}
		expect_gt(mean(ratio), if (kind == "ar1") 0.96 else 0.95)
		expect_lt(mean(ratio), 1.03)
	}
})

test_that("symmetric dynamic_ab and bipartite position paths enter the ve df in full", {
	# the symmetric additive path counts once on the unique-pair scale, and bipartite position paths count by their
	# joint smoother with the residual curvature of the bilinear term
	realised = function(sim) {
		var(unlist(lapply(seq_along(sim$E), function(t) {
			sim$E[[t]][is.finite(sim$Y[[t]]) & is.finite(sim$E[[t]])]
		})))
	}
	cell_s = grid_cell(
		"sab", "normal", "lame_dynamic_als", symmetric = TRUE,
		R = 1, dynamic_ab = TRUE, T = 5
	)
	ratio_s = vapply(1:3, function(s) {
		sim = grid_sim(cell_s, s)
		f = suppressWarnings(lame(
			sim$Y,
			sim$X,
			R = 1,
			family = "normal",
			symmetric = TRUE,
			method = "als",
			dynamic_ab = TRUE,
			verbose = FALSE,
			seed = s + 11L
		))
		unname(f$VC["ve"]) / realised(sim)
	}, numeric(1))
	expect_gt(mean(ratio_s), 0.99)
	expect_lt(mean(ratio_s), 1.02)
	cell_b = grid_cell(
		"bip", "normal", "lame_dynamic_als", mode = "bipartite",
		nA = 30, nB = 20, R_row = 2, R_col = 2,
		dynamic_uv = TRUE, T = 5
	)
	ratio_b = vapply(1:2, function(s) {
		sim = grid_sim(cell_b, s)
		f = suppressWarnings(lame(
			sim$Y,
			sim$X,
			family = "normal",
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			method = "als",
			dynamic_uv = TRUE,
			verbose = FALSE,
			seed = s + 11L
		))
		unname(f$VC["ve"]) / realised(sim)
	}, numeric(1))
	expect_gt(mean(ratio_b), 0.975)
	expect_lt(mean(ratio_b), 1.04)
})
