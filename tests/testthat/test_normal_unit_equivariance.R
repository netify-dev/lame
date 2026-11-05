# the normal family has no unit of its own: default priors and starting values are in the unit of the response
# (R/normal_scale.R, R/uv_prior.R), so a fit of c * y is the fit of y in other units; same-seed sampler draws scale
# exactly and the fast estimators scale to their stopping tolerance
skip_on_cran()

eq_net = function(mode = c("dir", "sym", "bip"), n = 12, nB = 9, Tt = 1, R = 2, seed = 11, miss = 0.05) {
	mode = match.arg(mode)
	set.seed(seed)
	m = if (mode == "bip") nB else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (mode == "bip") sprintf("c%02d", seq_len(m)) else rn
	a = rnorm(n, 0, 0.5); b = if (mode == "sym") a else rnorm(m, 0, 0.5)
	U = matrix(rnorm(n * 2), n, 2) * 0.7; V = if (mode == "sym") U else matrix(rnorm(m * 2), m, 2) * 0.7
	xr = matrix(
		rnorm(n),
		n,
		1,
		dimnames = list(rn, "xr")
	); xc = if (mode == "sym") xr else matrix(rnorm(m), m, 1, dimnames = list(cn, "xc"))
	Y = Xd = vector("list", Tt)
	for (t in seq_len(Tt)) {
		xd = matrix(rnorm(n * m), n, m)
		if (mode == "sym") xd = (xd + t(xd)) / 2
		mu = 0.5 + 0.8 * xd + 0.35 * matrix(xr, n, m) + 0.35 * matrix(xc, n, m, byrow = TRUE) + outer(a, b, "+") + U %*% diag(c(1.2, -0.8)) %*% t(V)
		e = matrix(rnorm(n * m), n, m)
		if (mode == "sym") e = (e + t(e)) / sqrt(2)
		y = mu + e
		if (mode != "bip") diag(y) = NA
		if (miss > 0) { k = sample(
			which(!is.na(y)),
			ceiling(miss * sum(!is.na(y)))
		); y[k] = NA; if (mode == "sym") y[t(is.na(y))] = NA }
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		Xd[[t]] = array(xd, c(n, m, 1), dimnames = list(rn, cn, "xd"))
	}
	list(Y = Y, Xdyad = Xd, Xrow = xr, Xcol = if (mode == "sym") NULL else xc, mode = mode, R = R)
}

# largest deviation of the draws of c * y from c^k times the draws of y, relative to the largest draw
eq_dev = function(f1, fc, cc) {
	rel = function(a, b) max(abs(a - b)) / max(abs(b), .Machine$double.xmin)
	vc1 = as.matrix(f1$VC); vcc = as.matrix(fc$VC)
	k = rep(2, ncol(vc1)); k[colnames(vc1) == "rho"] = 0
	out = c(
		beta = rel(as.matrix(fc$BETA), as.matrix(f1$BETA) * cc),
		vc = max(vapply(seq_len(ncol(vc1)), function(j) rel(vcc[, j], vc1[, j] * cc^k[j]), 0))
	)
	if (!is.null(f1$U_samples)) {
		prod = function(f) vapply(seq_len(dim(f$U_samples)[3]), function(s) {
			Um = matrix(
				f$U_samples[, , s],
				dim(f$U_samples)[1],
				dim(f$U_samples)[2]
			); Vm = matrix(f$V_samples[, , s], dim(f$V_samples)[1], dim(f$V_samples)[2])
			if (!is.null(f$G_samples)) Um %*% matrix(
				f$G_samples[, , s],
				dim(f$G_samples)[1],
				dim(f$G_samples)[2]
			) %*% t(Vm) else Um %*% t(Vm)
		}, matrix(0, dim(f$U_samples)[1], dim(f$V_samples)[1]))
		out["product"] = rel(prod(fc), prod(f1) * cc)
	}
	out
}

eq_fit = function(d, scale, fun, nscan = 25, burn = 8, ...) {
	args = list(
		family = "normal",
		nscan = nscan,
		burn = burn,
		odens = 1,
		seed = 7,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		...
	)
	if (d$mode == "bip") { args$mode = "bipartite"; args$R_row = d$R; args$R_col = d$R } else args$R = d$R
	if (d$mode == "sym") args$symmetric = TRUE
	if (fun == "ame") {
		args = c(
			list(
				Y = d$Y[[1]] * scale,
				Xdyad = d$Xdyad[[1]],
				Xrow = d$Xrow,
				Xcol = d$Xcol,
				posterior_opts = list(save_UV = TRUE)
			),
			args
		)
		suppressWarnings(suppressMessages(do.call(ame, args)))
	} else {
		args = c(
			list(
				Y = lapply(d$Y, function(y) y * scale),
				Xdyad = d$Xdyad,
				Xrow = rep(list(d$Xrow), length(d$Y)),
				Xcol = if (is.null(d$Xcol)) NULL else rep(list(d$Xcol), length(d$Y)),
				posterior_opts = list(save_UV_draws = TRUE)
			),
			args
		)
		suppressWarnings(suppressMessages(do.call(lame, args)))
	}
}

test_that("the unit of a normal response scales with the response and is one where the link fixes the scale", {
	set.seed(3)
	n = 14
	Y = matrix(rnorm(n * n, 2, 1.5), n) + outer(rnorm(n), rnorm(n), "+"); diag(Y) = NA
	for (fam in c("binary", "ordinal", "cbin", "frn", "poisson")) {
		expect_identical(lame:::.lame_unit_var(Y, fam), 1)
		expect_identical(lame:::.lame_unit_sd(Y, fam), 1)
		expect_identical(lame:::.lame_weight_prior_var(Y, fam), 1)
		expect_identical(lame:::.lame_sigma_uv_scale_default(Y, fam), 1)
		expect_identical(lame:::.lame_beta0_var_default(Y, fam), 10)
		expect_identical(lame:::.lame_rail_units(Y, fam), list(uv = 1, G = 1))
	}
	expect_equal(lame:::.lame_unit_var(Y, "normal"), var(c(Y), na.rm = TRUE))
	for (cc in c(1e-3, 1e3)) {
		expect_equal(
			lame:::.lame_unit_var(cc * Y, "normal"),
			cc^2 * lame:::.lame_unit_var(Y, "normal"),
			tolerance = 1e-12
		)
		# the additive scale carries the unit squared and the factor scale the unit, since u_i'v_j is in the unit of the response
		expect_equal(lame:::.lame_ab_scale(cc * Y), cc^2 * lame:::.lame_ab_scale(Y), tolerance = 1e-12)
		expect_equal(
			lame:::.lame_uv_prior_scale(cc * Y, "normal"),
			cc * lame:::.lame_uv_prior_scale(Y, "normal"),
			tolerance = 1e-12
		)
		expect_equal(
			lame:::.lame_beta0_var_default(cc * Y, "normal"),
			cc^2 * 10 * lame:::.lame_unit_var(Y, "normal"),
			tolerance = 1e-12
		)
		expect_equal(
			lame:::.lame_residual_prior_scale(cc * Y, "normal"),
			cc^2 * lame:::.lame_residual_prior_scale(Y, "normal"),
			tolerance = 1e-12
		)
		expect_equal(
			lame:::.lame_sigma_uv_scale_default(cc * Y, "normal"),
			cc * lame:::.lame_sigma_uv_scale_default(Y, "normal"),
			tolerance = 1e-12
		)
		# bipartite positions and symmetric coordinates are unit free
		expect_identical(lame:::.lame_sigma_uv_scale_default(cc * Y, "normal", bip = TRUE), 1)
		expect_identical(lame:::.lame_sigma_uv_scale_default(cc * Y, "normal", symmetric = TRUE), 1)
		ru = lame:::.lame_rail_units(cc * Y, "normal")
		expect_equal(ru$G, cc * lame:::.lame_unit_sd(Y, "normal"), tolerance = 1e-12)
		expect_equal(ru$uv, sqrt(ru$G), tolerance = 1e-12)
		expect_equal(lame:::.lame_rail_units(cc * Y, "normal", bip = TRUE)$uv, 1)
	}
	# the factor-row prior follows the unit: psi = (2R + 2) s
	pr1 = lame:::.lame_uv_prior(Y, "normal", 2L); pr3 = lame:::.lame_uv_prior(1e3 * Y, "normal", 2L)
	expect_equal(pr3$scale, 1e3 * pr1$scale, tolerance = 1e-12)
	expect_equal(pr3$psi, 1e3 * pr1$psi, tolerance = 1e-12)
	# floors of the dynamic scales are relative to the scale of the effects and the residual
	S0 = diag(2) * 1e-9
	expect_equal(
		lame:::.lame_sigma_ab_prior_scale(S0, 6, 0.8, TRUE, TRUE),
		1e-9 * lame:::.lame_sigma_ab_prior_scale(diag(2), 6, 0.8, TRUE, TRUE),
		tolerance = 1e-12
	)
	expect_equal(
		lame:::.lame_sigma_beta_prior_scale(1e-6 * Y, "normal"),
		1e-12 * lame:::.lame_sigma_beta_prior_scale(Y, "normal"),
		tolerance = 1e-10
	)
})

test_that("the samplers of ame() fit c * y as y: same-seed draws scale exactly, directed, symmetric and bipartite, ranks 0 and 2", {
	skip_on_cran()
	for (mode in c("dir", "sym", "bip")) for (R in c(0L, 2L)) {
		d = eq_net(mode, R = R)
		f1 = eq_fit(d, 1, "ame")
		for (cc in c(1e-3, 1e3)) {
			dv = eq_dev(f1, eq_fit(d, cc, "ame"), cc)
			# a symmetric fit solves a near-singular reciprocal system, so its rounding error is larger
			tol = if (mode == "sym") 1e-6 else 1e-9
			expect_lt(
				max(dv),
				tol,
				label = sprintf("ame %s R=%d c=%g (%s)", mode, R, cc, paste(names(dv), signif(dv, 2), collapse = " "))
			)
		}
	}
})

test_that("the samplers of lame() fit c * y as y: same-seed draws scale exactly, directed, symmetric and bipartite, ranks 0 and 2", {
	skip_on_cran()
	for (mode in c("dir", "sym", "bip")) for (R in c(0L, 2L)) {
		d = eq_net(mode, R = R, Tt = 2)
		f1 = eq_fit(d, 1, "lame")
		for (cc in c(1e-3, 1e3)) {
			dv = eq_dev(f1, eq_fit(d, cc, "lame"), cc)
			tol = if (mode == "sym") 1e-6 else 1e-9
			expect_lt(
				max(dv),
				tol,
				label = sprintf("lame %s R=%d c=%g (%s)", mode, R, cc, paste(names(dv), signif(dv, 2), collapse = " "))
			)
		}
	}
	# a longer burn-in reaches the scale split of the directed factors with their covariance held at n s I (R/uv_scale_split.R)
	d = eq_net("dir", R = 2L, Tt = 2)
	f1 = eq_fit(d, 1, "lame", burn = 60)
	dv = eq_dev(f1, eq_fit(d, 1e3, "lame", burn = 60), 1e3)
	expect_lt(
		max(dv),
		1e-9,
		label = sprintf("lame dir R=2 burn 60 (%s)", paste(names(dv), signif(dv, 2), collapse = " "))
	)
})

test_that("the dynamic samplers of lame() are equivariant too: moving effects, coefficients and positions", {
	skip_on_cran()
	d = eq_net("dir", n = 10, R = 2, Tt = 3)
	for (spec in list(
		list(dynamic_ab = TRUE),
		list(dynamic_beta = TRUE),
		list(dynamic_uv = TRUE, prior = list(als_start = FALSE))
	)) {
		f1 = do.call(eq_fit, c(list(d, 1, "lame", nscan = 20, burn = 6), spec))
		for (cc in c(1e-3, 1e3)) {
			fc = do.call(eq_fit, c(list(d, cc, "lame", nscan = 20, burn = 6), spec))
			expect_lt(max(eq_dev(f1, fc, cc)), 1e-8, label = paste(names(spec)[1], cc))
			if (!is.null(f1$sigma_ab)) expect_equal(fc$sigma_ab, cc * f1$sigma_ab, tolerance = 1e-9)
			if (!is.null(f1$sigma_beta)) expect_equal(fc$sigma_beta, cc * f1$sigma_beta, tolerance = 1e-9)
			# the innovation sd of unipartite positions carries the square root of the unit
			if (!is.null(f1$sigma_uv)) expect_equal(fc$sigma_uv, sqrt(cc) * f1$sigma_uv, tolerance = 1e-9)
		}
	}
})

test_that("the bipartite interaction no longer collapses at a large unit and the symmetric weights follow the unit", {
	d = eq_net("bip", R = 2, n = 14, nB = 10)
	# at sd(y) of order 1e3 the interaction prior variance is in the unit of the response, so it is not shrunk to zero
	f = eq_fit(d, 1e4, "ame", nscan = 60, burn = 20)
	g = eq_fit(d, 1, "ame", nscan = 60, burn = 20)
	expect_equal(mean(abs(f$G)) / 1e4, mean(abs(g$G)), tolerance = 1e-8)
	expect_gt(mean(abs(f$G)), 1e3)
	# symmetric kernel: the prior variance of an eigen-weight enters its conditional precision
	set.seed(1); n = 8
	U = matrix(rnorm(n * 2), n); E = U %*% diag(c(2, -1)) %*% t(U); E = (E + t(E)) / 2
	ids = as.integer(rep(sample(n), 4) - 1L)
	set.seed(5); a = lame:::rUV_sym_fc_cpp(E * 50, U, U %*% diag(c(50, -50)), 2500, TRUE, ids, 2500)
	set.seed(5); b = lame:::rUV_sym_fc_cpp(E, U, U %*% diag(c(1, -1)), 1, TRUE, ids, 1)
	expect_equal(a$U, b$U, tolerance = 1e-9)
	expect_equal(a[["V"]], 50 * b[["V"]], tolerance = 1e-9)
	expect_error(lame:::rUV_sym_fc_cpp(E, U, U, 1, TRUE, ids, -1), "positive")
})

test_that("lame_robust() fits c * y as y: same-seed draws of the coefficients and the residual variance scale exactly", {
	d = eq_net("dir", n = 10, R = 1, Tt = 3, miss = 0)
	fit = function(scale) {
		Y = lapply(d$Y, function(y) y * scale); names(Y) = paste0("t", seq_along(Y))
		X = d$Xdyad; names(X) = names(Y)
		suppressWarnings(
			suppressMessages(
				lame_robust(Y, Xdyad = X, R = 1L, burn = 20L, nscan = 40L, odens = 2L, seed = 3L, verbose = FALSE)
			)
		)
	}
	f1 = fit(1)
	for (cc in c(1e-3, 1e3)) {
		fc = fit(cc)
		rel = function(a, b) max(abs(a - b)) / max(abs(b))
		expect_lt(rel(fc$draws$beta, f1$draws$beta * cc), 1e-8, label = sprintf("beta c=%g", cc))
		expect_lt(rel(fc$draws$hyper[, "ve"], f1$draws$hyper[, "ve"] * cc^2), 1e-8, label = sprintf("ve c=%g", cc))
	}
})

test_that("ame_als() and lame_als() fit c * y as y to their stopping tolerance", {
	skip_on_cran()
	for (mode in c("dir", "sym", "bip")) for (R in c(0L, 2L)) {
		d = eq_net(mode, R = R, Tt = 2, n = 14, nB = 10)
		fit = function(scale, fun) {
			a = list(family = "normal", R = R, verbose = FALSE)
			if (mode == "bip") a$mode = "bipartite"
			if (mode == "sym") a$symmetric = TRUE
			if (fun == "ame_als") {
				a = c(list(d$Y[[1]] * scale, Xdyad = d$Xdyad[[1]], Xrow = d$Xrow, Xcol = d$Xcol), a)
				suppressWarnings(suppressMessages(do.call(ame_als, a)))
			} else {
				a = c(
					list(
						lapply(d$Y, function(y) y * scale),
						Xdyad = d$Xdyad,
						Xrow = rep(list(d$Xrow), 2),
						Xcol = if (is.null(d$Xcol)) NULL else rep(list(d$Xcol), 2)
					),
					a
				)
				suppressWarnings(suppressMessages(do.call(lame_als, a)))
			}
		}
		for (fun in c("ame_als", "lame_als")) {
			f1 = fit(1, fun)
			for (cc in c(1e-3, 1e3)) {
				fc = fit(cc, fun)
				rel = function(x, y) max(abs(x - y)) / max(abs(y))
				lab = sprintf("%s %s R=%d c=%g", fun, mode, R, cc)
				# the fits stop at a tolerance of 1e-6 in the unit of each fit, so they agree to about that
				expect_lt(rel(unname(fc$coefficients), cc * unname(f1$coefficients)), 3e-4, label = lab)
				expect_lt(rel(fc$VC[["ve"]], cc^2 * f1$VC[["ve"]]), 3e-4, label = lab)
				if (R > 0L) {
					# the product withholding decision (.ae_gauss_withhold) and the product are unit free
					expect_identical(
						lame:::.ae_product_is_withheld(fc),
						lame:::.ae_product_is_withheld(f1),
						label = lab
					)
					prod = function(f) if (lame:::.ae_product_is_withheld(f)) f$eb$product_withheld$means$UVPM else f$UVPM
					expect_lt(rel(prod(fc), cc * prod(f1)), 3e-4, label = lab)
				}
				expect_lt(rel(unlist(fc$a), cc * unlist(f1$a)), 3e-4, label = lab)
			}
		}
	}
})

test_that("the chain check judges a series in any unit as in its own", {
	set.seed(3)
	s = list(intercept = cumsum(rnorm(300)) / 30 + rnorm(300), ve = 2 + rnorm(300, 0, .1), rho = rep(0.4, 300))
	st1 = lame:::.lame_chain_stats(s)
	for (cc in c(1e-12, 1e-6, 1e6)) {
		stc = lame:::.lame_chain_stats(lapply(s, function(x) cc * x))
		expect_identical(stc$parameter, st1$parameter)
		expect_equal(stc$z, st1$z, tolerance = 1e-10)
		expect_equal(stc$shift, st1$shift, tolerance = 1e-10)
	}
	# a constant series is left out
	expect_false("rho" %in% st1$parameter)
})
