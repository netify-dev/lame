# the product of these small panels is withheld in the default fit, so the validation reads it with the withholding off
withr::local_options(list(lame.eb_withhold_product = FALSE))

library(lame)

# recovery and equivalence checks of the dynamic empirical-Bayes estimator (R/lame_dynamic_als_eb.R)
# simulations draw from the dynamic_ab model: a level plus a zero-mean AR(1) deviation per actor, static factors, probit links
# the fast contract tests are in test_lame_dynamic_als_eb.R

vsim = function(
	n = 30,
	nb = n,
	Tt = 5,
	mode = "directed",
	seed = 1,
	mu = -1.2,
	beta = c(0.6, -0.5),
	va = 0.6,
	vb = 0.6,
	rho = 0.7,
	sigma = 0.35,
	R = 0,
	suv = 0.5,
	beta_sd = 0,
	quiet = 0,
	sender = TRUE,
	dyn_uv = FALSE,
	rho_uv = 0.9,
	sigma_uv = 0.3,
	dcor = 0
) {
	set.seed(seed)
	sym = identical(mode, "symmetric")
	bip = identical(mode, "bipartite")
	if (!bip) nb = n
	ar1 = function(m) {
		out = matrix(0, m, Tt)
		out[, 1] = rnorm(m, 0, sigma / sqrt(1 - rho^2))
		for (t in 2:Tt) out[, t] = rho * out[, t - 1] + rnorm(m, 0, sigma)
		out
	}
	alpha = rnorm(n, 0, sqrt(va))
	gam = if (sym) alpha else rnorm(nb, 0, sqrt(vb))
	if (quiet > 0) {
		alpha[seq_len(quiet)] = alpha[seq_len(quiet)] - 4
		gam[seq_len(quiet)] = gam[seq_len(quiet)] - 4
		if (sym) gam = alpha
	}
	d = ar1(n)
	e = if (sym) d else ar1(nb)
	# a panel without sender effects: neither levels nor paths
	if (!sender) { alpha[] = 0; d[] = 0 }
	bpath = matrix(beta, length(beta), Tt)
	if (beta_sd > 0) for (t in 2:Tt) bpath[1, t] = bpath[1, t - 1] + rnorm(1, 0, beta_sd)
	if (R > 0 && !dyn_uv) {
		U = matrix(rnorm(n * R, 0, sqrt(suv)), n, R)
		V = if (sym) U else matrix(rnorm(nb * R, 0, sqrt(suv)), nb, R)
	}
	# moving positions: every coordinate a zero-mean stationary AR(1)
	if (R > 0 && dyn_uv) {
		ar1uv = function(m) {
			out = matrix(0, m, Tt)
			out[, 1] = rnorm(m, 0, sigma_uv / sqrt(1 - rho_uv^2))
			for (t in 2:Tt) out[, t] = rho_uv * out[, t - 1] + rnorm(m, 0, sigma_uv)
			out
		}
		Ut = lapply(seq_len(R), function(r) ar1uv(n))
		Vt = if (sym) Ut else lapply(seq_len(R), function(r) ar1uv(nb))
	}
	Otrue = array(0, c(n, nb, Tt))
	rn = sprintf("a%03d", seq_len(n))
	cn = if (bip) sprintf("b%03d", seq_len(nb)) else rn
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		xt = array(
			rnorm(n * nb * length(beta)),
			c(n, nb, length(beta)),
			dimnames = list(rn, cn, paste0("x", seq_along(beta)))
		)
		if (sym) for (k in seq_along(beta)) {
			m = xt[, , k]
			m[lower.tri(m)] = t(m)[lower.tri(m)]
			xt[, , k] = m
		}
		eta = mu + outer(alpha + d[, t], gam + e[, t], "+")
		for (k in seq_along(beta)) eta = eta + bpath[k, t] * xt[, , k]
		if (R > 0 && !dyn_uv) eta = eta + tcrossprod(U, V)
		if (R > 0 && dyn_uv) {
			Ot = Reduce(`+`, lapply(seq_len(R), function(r) tcrossprod(Ut[[r]][, t], Vt[[r]][, t])))
			eta = eta + Ot
			Otrue[, , t] = Ot
		}
		err = matrix(rnorm(n * nb), n, nb)
		# the two errors of a dyad are correlated (directed network; the draws exist only when dcor is set)
		if (dcor != 0 && !sym && !bip) {
			lo = lower.tri(err)
			err[lo] = dcor * t(err)[lo] + sqrt(1 - dcor^2) * matrix(rnorm(n * nb), n, nb)[lo]
		}
		z = eta + err
		if (sym) z[lower.tri(z)] = t(z)[lower.tri(z)]
		y = 1 * (z > 0)
		dimnames(y) = list(rn, cn)
		if (!bip) diag(y) = NA
		Y[[t]] = y
		X[[t]] = xt
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(
		Y = Y,
		X = X,
		truth = list(
			mu = mu,
			beta = beta,
			bpath = bpath,
			alpha = alpha,
			gam = gam,
			d = d,
			e = e,
			O = Otrue
		)
	)
}

eb_run = function(d, mode = "directed", bias_correct = NULL, dynamic_ab = TRUE, ...) {
	md = if (identical(mode, "bipartite")) "bipartite" else "unipartite"
	suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		method = "als",
		mode = md,
		symmetric = identical(mode, "symmetric"),
		dynamic_ab = dynamic_ab,
		bias_correct = bias_correct,
		verbose = FALSE,
		...
	)))
}

# mean relative bias of the two slopes over replicates
slope_bias = function(fits, truth = c(0.6, -0.5)) {
	B = do.call(rbind, lapply(fits, function(f) f$coefficients[c("x1_dyad", "x2_dyad")]))
	100 * (colMeans(B) / truth - 1)
}

test_that("with static effects the dynamic engine is the static empirical-Bayes fit", {
	skip_on_cran()
	for (R in c(0L, 1L, 2L)) {
		d = vsim(
			n = 30,
			Tt = 4,
			seed = 7,
			mu = -1.2,
			va = 0.5,
			vb = 0.5,
			sigma = 0.001,
			R = R,
			suv = 0.6
		)
		prep = lame:::.ae_prepare(
			d$Y,
			d$X,
			NULL,
			NULL,
			"unipartite",
			longitudinal = TRUE,
			node_time_policy = "quiet"
		)
		ctx = lame:::.lde_context(
			prep,
			FALSE,
			"probit",
			coef_prior = "flat",
			dynamic_ab = FALSE,
			R = R
		)
		fit = lame:::.lde_fit(ctx, max_iter = 300, tol = 1e-6)
		expect_true(fit$converged)
		th = fit$st$theta
		th[1] = th[1] - sum(ctx$centre * th)
		# covariance of the coefficients, with factors marginal over them and over the scale of their prior (as the static estimator reports it)
		Vc = if (R > 0L) lame:::.lde_fac_vcov(fit$ctx, fit$st, fit$ev)$vcov else fit$ev$fac$Szi[ctx$q + seq_len(ctx$mth), ctx$q + seq_len(ctx$mth)]
		se = sqrt(diag(Vc))
		run = function(...) suppressWarnings(suppressMessages(lame:::.ae_run(
			d$Y,
			d$X,
			NULL,
			NULL,
			R,
			"binary",
			"unipartite",
			FALSE,
			300,
			1e-8,
			FALSE,
			6886,
			longitudinal = TRUE,
			call = quote(lame_als()),
			lowrank_method = "mm",
			non_normal_method = "irls",
			link = "probit",
			linear_solver = "eigen",
			multistart = "none",
			bias_correct = NULL,
			...
		)))
		# the dynamic engine's factors are the static estimator's variational posterior, so the reference is the static default fit
		# the two treat the actor effects and the working weights differently, so the coefficients agree to about a percent
		st = run()
		expect_equal(unname(th), unname(st$coefficients), tolerance = 1e-2)
		expect_equal(unname(se), unname(st$eb$se), tolerance = 6e-2)
		expect_equal(as.vector(fit$st$hyp$sab), as.vector(st$eb$sab), tolerance = 0.05)
		if (R > 0L) {
			expect_identical(st$eb$estimator, "variational")
			# the factor prior scale follows the static estimator's cavity rule (an isotropic prior of one scale)
			expect_identical(st$eb$variational$scale$rule, "cavity")
			expect_equal(fit$st$hyp$fac$Suv[1L, 1L], st$eb$variational$scale$s, tolerance = 0.05)
			md = run(eb_posterior = "mode")
			expect_identical(md$eb$estimator, "laplace")
			expect_lt(max(abs(th - st$coefficients)), max(abs(th - md$coefficients)))
		}
	}
})

test_that("the empirical-Bayes slopes are nearly unbiased where the fixed-effects ones are not", {
	skip_on_cran()
	fits = lapply(1:16, function(r) {
		d = vsim(n = 30, Tt = 5, seed = 1000 + r, mu = -1.3)
		list(eb = eb_run(d), none = eb_run(d, bias_correct = "none"))
	})
	eb = lapply(fits, `[[`, "eb")
	none = lapply(fits, `[[`, "none")
	expect_true(all(vapply(eb, function(f) identical(f$bias_correct$method, "eb"), logical(1))))
	expect_true(all(vapply(eb, function(f) isTRUE(f$converged), logical(1))))
	expect_lt(max(abs(slope_bias(eb))), 5)
	expect_gt(min(abs(slope_bias(none))), 6)
	# the uncorrected penalised fit shrinks both slopes towards zero, the
	# corrected one does not
	expect_true(all(slope_bias(none) < 0))
	# the levels and the deviation scale are recovered to within the sampling noise
	hv = rowMeans(sapply(eb, function(f) c(f$eb$sab[1, 1], f$eb$sab[2, 2], f$rho_ab, f$sigma_ab)))
	expect_equal(unname(hv), c(0.6, 0.6, 0.7, 0.35), tolerance = 0.25)
	# the Laplace standard errors are of the order of the sampling spread (a coarse check)
	B = do.call(rbind, lapply(eb, function(f) f$coefficients[c("x1_dyad", "x2_dyad")]))
	S = do.call(rbind, lapply(eb, function(f) f$eb$se[c("x1_dyad", "x2_dyad")]))
	expect_gt(min(colMeans(S) / apply(B, 2, stats::sd)), 0.5)
	expect_lt(max(colMeans(S) / apply(B, 2, stats::sd)), 2)
})

test_that("undirected and bipartite panels are corrected as well", {
	skip_on_cran()
	for (mode in c("symmetric", "bipartite")) {
		fits = lapply(1:12, function(r) {
			d = vsim(n = 36, nb = 26, Tt = 5, mode = mode, seed = 200 + r, mu = -1.2)
			list(eb = eb_run(d, mode = mode), none = eb_run(d, mode = mode, bias_correct = "none"))
		})
		eb = lapply(fits, `[[`, "eb")
		none = lapply(fits, `[[`, "none")
		expect_true(all(vapply(eb, function(f) isTRUE(f$converged), logical(1))))
		expect_lt(max(abs(slope_bias(eb))), 7)
		expect_gt(mean(abs(slope_bias(none))), mean(abs(slope_bias(eb))))
	}
})

test_that("static latent factors are fitted jointly on directed, undirected and bipartite panels", {
	skip_on_cran()
	options_w = options(lame.eb_withhold_intercept = FALSE)
	on.exit(options(options_w), add = TRUE)
	for (mode in c("directed", "symmetric", "bipartite")) {
		fits = lapply(1:6, function(r) {
			d = vsim(n = 40, nb = 30, Tt = 4, mode = mode, seed = 300 + r, mu = -1.2, R = 1, suv = 0.6)
			eb_run(d, mode = mode, R = 1)
		})
		expect_true(all(vapply(fits, function(f) isTRUE(f$converged), logical(1))))
		expect_true(all(vapply(fits, function(f) all(is.finite(f$eb$se)), logical(1))))
		expect_lt(max(abs(slope_bias(fits))), 8)
		expect_equal(fits[[1]]$R, 1)
	}
})

test_that("moving coefficients are recovered with their own persistence and scale", {
	skip_on_cran()
	fits = lapply(1:8, function(r) {
		d = vsim(n = 40, Tt = 8, seed = 400 + r, mu = -1.2, beta_sd = 0.15)
		f = eb_run(d, dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw1")
		list(f = f, truth = d$truth$bpath[1, ])
	})
	expect_true(all(vapply(fits, function(z) isTRUE(z$f$converged), logical(1))))
	cors = vapply(fits, function(z) stats::cor(z$f$coef_path["x1_dyad", ], z$truth), numeric(1))
	expect_gt(stats::median(cors), 0.7)
	cover = vapply(fits, function(z) {
		mean(abs(z$f$coef_path["x1_dyad", ] - z$truth) < 1.96 * z$f$eb$se_path["x1_dyad", ])
	}, numeric(1))
	expect_gt(mean(cover), 0.8)
	# the static coefficient stays a constant path
	expect_true(all(vapply(fits, function(z) {
		diff(range(z$f$coef_path["x2_dyad", ])) < 1e-8
	}, logical(1))))
})

test_that("an actor without events has a finite, negative effect", {
	skip_on_cran()
	d = vsim(n = 30, Tt = 5, seed = 11, mu = -1.2, quiet = 3)
	# the three quiet senders send (almost) nothing in any period
	sent = sum(vapply(d$Y, function(y) sum(y[1:3, ], na.rm = TRUE), numeric(1)))
	expect_lt(sent, 15)
	f = eb_run(d)
	expect_true(isTRUE(f$converged))
	expect_true(all(is.finite(f$a_level)) && all(is.finite(f$b_level)))
	expect_true(all(f$a_level[1:3] < 0))
	expect_true(all(is.finite(f$eb$se)))
	# and not minus infinity
	expect_gt(min(f$a_dynamic[1:3, ]), -6)
})

test_that("node covariates, missing cells and a changing composition are handled", {
	skip_on_cran()
	d = vsim(n = 24, Tt = 5, seed = 21, mu = -1.2)
	rn = rownames(d$Y[[1]])
	Xr = lapply(1:5, function(t) matrix(rnorm(24), 24, 1, dimnames = list(rn, "wr")))
	Xc = lapply(1:5, function(t) matrix(rnorm(24), 24, 1, dimnames = list(rn, "wc")))
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	)))
	expect_true(isTRUE(f$converged))
	expect_true(all(c("wr_row", "wc_col") %in% names(f$coefficients)))
	expect_true(all(is.finite(f$eb$se)))
	# missing cells
	Y2 = d$Y
	set.seed(5)
	for (t in 1:5) Y2[[t]][sample(24 * 24, 40)] = NA
	diag_na = lapply(Y2, function(y) { diag(y) = NA; y })
	f2 = suppressWarnings(suppressMessages(lame(
		diag_na,
		Xdyad = d$X,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	)))
	expect_true(isTRUE(f2$converged))
	# a changing composition: actors leave and enter
	keep = lapply(1:5, function(t) setdiff(rn, rn[((t - 1) * 3 + 1):((t - 1) * 3 + 3)]))
	Y3 = lapply(1:5, function(t) d$Y[[t]][keep[[t]], keep[[t]]])
	X3 = lapply(1:5, function(t) d$X[[t]][keep[[t]], keep[[t]], , drop = FALSE])
	names(Y3) = names(X3) = names(d$Y)
	f3 = suppressWarnings(suppressMessages(lame(
		Y3,
		Xdyad = X3,
		family = "binary",
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	)))
	expect_true(isTRUE(f3$converged))
	expect_true(isTRUE(f3$changing_composition))
	expect_true(all(is.finite(f3$eb$se)))
})

test_that("the generic methods work on an empirical-Bayes dynamic fit", {
	skip_on_cran()
	d = vsim(n = 20, Tt = 4, seed = 31, mu = -1.2, R = 1)
	f = eb_run(d, R = 1)
	expect_s3_class(f, "lame_dynamic_als")
	msg = capture.output(print(f))
	expect_true(any(grepl("empirical-Bayes", msg)))
	msg = capture.output(print(summary(f)))
	expect_true(any(grepl("Laplace", msg)))
	expect_equal(names(coef(f)), c("intercept", "x1_dyad", "x2_dyad"))
	V = vcov(f)
	expect_equal(dim(V), c(3L, 3L))
	expect_true(all(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 0))
	ci = confint(f)
	expect_equal(dim(ci), c(3L, 2L))
	expect_true(all(ci[, 1] < coef(f) & coef(f) < ci[, 2]))
	pr = predict(f)
	expect_true(length(pr) > 0)
	expect_equal(length(fitted(f)), 4L)
	expect_equal(length(residuals(f)), 4L)
	sv = als_start_vals(f)
	expect_true(all(c("rho_ab", "sigma_ab") %in% names(sv)))
	expect_true(is.finite(sv$rho_ab) && is.finite(sv$sigma_ab))
	# a warm restart from the fit stays at the fixed point
	g = suppressWarnings(suppressMessages(ame_als_refit(f)))
	expect_equal(unname(coef(g)), unname(coef(f)), tolerance = 1e-3)
})

test_that("a warm restart of a fit with static factors stays at the fixed point", {
	skip_on_cran()
	# the factors, their scales and the product are part of the start, so a converged fit is its own fixed point
	d = vsim(n = 14, Tt = 4, seed = 71, mu = -1, R = 1)
	f = eb_run(d, R = 1)
	g = suppressWarnings(suppressMessages(ame_als_refit(f)))
	expect_equal(g$iterations, 0)
	expect_equal(unname(coef(g)), unname(coef(f)), tolerance = 1e-8)
	expect_equal(g$Omat, f$Omat, tolerance = 1e-8)
	# a rank-two directed fit, a bipartite fit (the larger factor padded with zero columns) and an undirected one
	f2 = eb_run(d, R = 2)
	g2 = suppressWarnings(suppressMessages(ame_als_refit(f2)))
	expect_equal(g2$iterations, 0)
	expect_equal(unname(coef(g2)), unname(coef(f2)), tolerance = 1e-8)
	db = vsim(n = 16, nb = 11, Tt = 4, mode = "bipartite", seed = 72, mu = -1, R = 1)
	fb = suppressMessages(eb_run(db, mode = "bipartite", R_row = 2, R_col = 1))
	gb = suppressWarnings(suppressMessages(ame_als_refit(fb)))
	expect_equal(gb$iterations, 0)
	expect_equal(unname(coef(gb)), unname(coef(fb)), tolerance = 1e-8)
	ds = vsim(n = 16, Tt = 4, mode = "symmetric", seed = 73, mu = -1, R = 1)
	fs = eb_run(ds, mode = "symmetric", R = 1)
	gs = suppressWarnings(suppressMessages(ame_als_refit(fs)))
	expect_equal(gs$iterations, 0)
	expect_equal(unname(coef(gs)), unname(coef(fs)), tolerance = 1e-8)
	# moving coefficients: the period-0 state and the hyperparameters of the blocks are part of the record
	for (cfg in list(
		list(dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw1", dynamic_ab = FALSE),
		list(dynamic_beta = c("intercept", "x1_dyad"), dynamic_beta_kind = "ar1"),
		list(dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw2", R = 1),
		list(dynamic_beta = c("intercept", "x1_dyad"), dynamic_beta_kind = "rw1", dynamic_beta_pool = "sigma"),
		list(rvar = FALSE)
	)) {
		fm = do.call(eb_run, c(list(d), cfg))
		gm = suppressWarnings(suppressMessages(ame_als_refit(fm)))
		expect_equal(gm$iterations, 0)
		expect_equal(unname(coef(gm)), unname(coef(fm)), tolerance = 1e-8)
		if (!is.null(fm$coef_path)) expect_equal(gm$coef_path, fm$coef_path, tolerance = 1e-8)
	}
	# a start of another rank or with other actors does not start the factors, so the refit runs
	# the restart record wins over the reported fields, so the edit drops it
	seed = lame:::.ae_static_seed(f)
	seed$eb$restart = NULL
	seed$U = seed$U[1:5, , drop = FALSE]
	h = suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
		d$Y,
		d$X,
		family = "binary",
		R = 1,
		dynamic_ab = TRUE,
		start = seed,
		verbose = FALSE
	)))
	expect_true(h$converged)
	expect_gt(h$iterations, 0)
	expect_equal(unname(coef(h)), unname(coef(f)), tolerance = 5e-3)
})

test_that("a warm restart of a fit with moving positions stays at the fixed point", {
	skip_on_cran()
	d = vsim(n = 14, Tt = 4, seed = 74, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	for (spec in list(
		list(mode = "directed", d = d, args = list(R = 1)),
		list(
			mode = "bipartite",
			d = vsim(
				n = 16, nb = 11, Tt = 4, mode = "bipartite", seed = 75, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4
			),
			args = list(R = 1)
		),
		list(
			mode = "symmetric",
			d = vsim(n = 16, Tt = 4, mode = "symmetric", seed = 76, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4),
			args = list(R = 1)
		),
		list(mode = "directed", d = d, args = list(R = 1, dynamic_beta = "x1_dyad", dynamic_beta_kind = "rw1"))
	)) {
		f = do.call(eb_run, c(list(spec$d, mode = spec$mode, dynamic_uv = TRUE), spec$args))
		# a refit of a converged fit gives no warning
		expect_no_warning(g <- suppressMessages(ame_als_refit(f)))
		expect_equal(g$iterations, 0)
		expect_equal(unname(coef(g)), unname(coef(f)), tolerance = 1e-8)
		# the product moves by one sweep of the position update at most
		expect_equal(g$Oarr, f$Oarr, tolerance = 5e-3)
		expect_equal(g$eb$uv_scale, f$eb$uv_scale)
	}
})

test_that("a refit of moving positions on changed data searches the scale afresh", {
	skip_on_cran()
	d = vsim(n = 14, Tt = 4, seed = 74, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	f = eb_run(d, R = 1, dynamic_uv = TRUE)
	# flipped cells move the fixed point, so the refit searches the scale again and agrees with a cold fit of the new data
	Yn = f$Y
	set.seed(5)
	for (t in seq_len(dim(Yn)[3L])) {
		cells = which(is.finite(Yn[, , t]))
		pick = sample(cells, 6L)
		Yn[, , t][pick] = 1 - Yn[, , t][pick]
	}
	g = suppressWarnings(suppressMessages(ame_als_refit(f, Y_new = Yn)))
	expect_true(g$converged)
	expect_gt(g$iterations, 0)
	expect_false(isTRUE(all.equal(g$eb$uv_scale, f$eb$uv_scale, tolerance = 1e-6)))
	Yl = lapply(seq_len(dim(Yn)[3L]), function(t) Yn[, , t])
	names(Yl) = names(d$Y)
	cold = eb_run(list(Y = Yl, X = d$X), R = 1, dynamic_uv = TRUE)
	expect_equal(unname(coef(g)), unname(coef(cold)), tolerance = 2e-2)
	expect_equal(g$eb$uv_scale, cold$eb$uv_scale, tolerance = 0.15)
	expect_s3_class(g$eb$uv_search, "data.frame")
	# the same data with a stricter tolerance is a search again, and lands where the fit was
	h = suppressWarnings(suppressMessages(ame_als_refit(f, tol = 1e-8, max_iter = 200)))
	expect_equal(unname(coef(h)), unname(coef(f)), tolerance = 5e-3)
	expect_equal(h$eb$uv_scale, f$eb$uv_scale, tolerance = 0.1)
})

test_that("a converged moving-position fit is within the tolerance of its fixed point on a slowly contracting panel", {
	skip_on_cran()
	for (seed in c(152, 153)) {
		d = vsim(n = 18, Tt = 5, seed = seed, mu = -0.5, R = 1, dyn_uv = TRUE, rho_uv = 0.95, sigma_uv = 1.0)
		prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
		ctx = lame:::.lde_context(prep, FALSE, "probit", R = 1L, dynamic_uv = TRUE)
		fit = lame:::.lde_fit(ctx, lame:::.lde_init(ctx), max_iter = 600L, tol = 5e-5)
		expect_true(fit$converged)
		# at a tolerance a million times tighter the objective and the coefficients barely move
		tight = lame:::.lde_fit_core_dyn(fit$ctx, fit$st, max_iter = 400L, tol = 1e-9)
		expect_lt(abs(tight$ev$objective - fit$ev$objective), 2e-4)
		expect_lt(max(abs(tight$st$theta - fit$st$theta)), 1e-4)
	}
})

test_that("moving positions: an actor without an observed cell has the prior's sd and zero means", {
	skip_on_cran()
	d = vsim(n = 12, Tt = 4, seed = 83, mu = -0.8, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	for (t in seq_along(d$Y)) { d$Y[[t]][3, ] = NA; d$Y[[t]][, 3] = NA }
	f = eb_run(d, R = 1, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_true(all(f$U[3, , ] == 0) && all(f$V[3, , ] == 0))
	expect_equal(unique(as.vector(f$U_sd[3, , ])), sqrt(f$eb$uv_scale))
	expect_equal(unique(as.vector(f$V_sd[3, , ])), sqrt(f$eb$uv_scale))
	expect_true(all(f$U_sd[-3, , ] < sqrt(f$eb$uv_scale)))
})

test_that("moving positions: an undirected panel with a changing composition is fitted and restarts at the fixed point", {
	skip_on_cran()
	d = vsim(n = 16, Tt = 4, mode = "symmetric", seed = 86, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	for (t in 2:3) { d$Y[[t]][5, ] = NA; d$Y[[t]][, 5] = NA }
	f = eb_run(d, mode = "symmetric", R = 1, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_true(all(is.finite(f$U)) && all(is.finite(f$U_sd)))
	expect_no_warning(g <- suppressMessages(ame_als_refit(f)))
	expect_equal(g$iterations, 0)
	expect_equal(unname(coef(g)), unname(coef(f)), tolerance = 1e-8)
})

test_that("moving positions: two periods, and static additive effects, are fitted and restart at the fixed point", {
	skip_on_cran()
	d2 = vsim(n = 14, Tt = 2, seed = 84, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	f2 = eb_run(d2, R = 1, dynamic_uv = TRUE)
	expect_true(f2$converged)
	expect_equal(dim(f2$U), c(14L, 1L, 2L))
	expect_no_warning(g2 <- suppressMessages(ame_als_refit(f2)))
	expect_equal(g2$iterations, 0)
	expect_equal(unname(coef(g2)), unname(coef(f2)), tolerance = 1e-8)
	d = vsim(n = 14, Tt = 4, seed = 85, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	fs = eb_run(d, R = 1, dynamic_uv = TRUE, dynamic_ab = FALSE)
	expect_true(fs$converged)
	expect_null(fs$rho_ab)
	expect_no_warning(gs <- suppressMessages(ame_als_refit(fs)))
	expect_equal(gs$iterations, 0)
	expect_equal(unname(coef(gs)), unname(coef(fs)), tolerance = 1e-8)
})

test_that("moving positions: an undirected panel is fitted, and ignoring the positions biases the slopes", {
	skip_on_cran()
	fits = lapply(61:63, function(seed) {
		d = vsim(n = 30, Tt = 5, mode = "symmetric", seed = seed, mu = -1.5, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
		list(d = d, uv = eb_run(d, mode = "symmetric", R = 1, dynamic_uv = TRUE), r0 = eb_run(d, mode = "symmetric"))
	})
	expect_true(all(vapply(fits, function(f) isTRUE(f$uv$converged), logical(1))))
	expect_lt(max(abs(slope_bias(lapply(fits, `[[`, "uv")))), 15)
	expect_gt(min(abs(slope_bias(lapply(fits, `[[`, "r0")))), 8)
	f = fits[[1]]$uv
	expect_equal(dim(f$U), c(30L, 1L, 5L))
	expect_equal(unname(f$L), 1)
	for (t in 1:5) {
		expect_equal(f$Oarr[, , t], t(f$Oarr[, , t]), tolerance = 1e-10)
		expect_equal(
			f$Oarr[, , t] + f$eb$centre[["o"]],
			f$U[, , t] %*% diag(f$L, 1) %*% t(f$U[, , t]),
			ignore_attr = TRUE
		)
	}
	for (fit in fits) {
		ob = array(FALSE, dim(fit$d$truth$O))
		for (t in seq_along(fit$d$Y)) ob[, , t] = is.finite(fit$d$Y[[t]])
		expect_gt(stats::cor(fit$uv$Oarr[ob], fit$d$truth$O[ob]), 0.5)
	}
	# a negative direction (heterophily) keeps its sign: the sign pattern of the start is part of the fit
	fh = eb_run(fits[[1]]$d, mode = "symmetric", R = 2, dynamic_uv = TRUE)
	expect_equal(length(fh$L), 2L)
	expect_true(all(abs(fh$L) == 1))
})

test_that("moving positions are an empirical-Bayes fit by default and an uncorrected one on request", {
	skip_on_cran()
	d = vsim(n = 12, Tt = 4, seed = 41, mu = -1, R = 1, dyn_uv = TRUE)
	fn = eb_run(d, R = 1, dynamic_uv = TRUE, bias_correct = "none")
	expect_identical(fn$bias_correct$method, "none")
	expect_true(isTRUE(fn$dynamic_uv))
	expect_null(fn$eb)
	fe = eb_run(d, R = 1, dynamic_uv = TRUE)
	expect_identical(fe$bias_correct$method, "eb")
	expect_true(isTRUE(fe$dynamic_uv))
	expect_true(is.list(fe$eb) && is.finite(fe$eb$rho_uv) && fe$eb$sigma_uv > 0)
	# static latent positions are fitted
	fs = eb_run(d, R = 1)
	expect_identical(fs$bias_correct$method, "eb")
	expect_false(isTRUE(fs$dynamic_uv))
})

test_that("moving positions: the slopes are nearly unbiased where ignoring the positions biases them", {
	skip_on_cran()
	fits = lapply(51:53, function(seed) {
		d = vsim(n = 30, Tt = 5, seed = seed, mu = -1.2, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
		list(d = d, uv = eb_run(d, R = 1, dynamic_uv = TRUE), r0 = eb_run(d))
	})
	expect_true(all(vapply(fits, function(f) isTRUE(f$uv$converged), logical(1))))
	b_uv = slope_bias(lapply(fits, `[[`, "uv"))
	b_0 = slope_bias(lapply(fits, `[[`, "r0"))
	# the product here carries a quarter of the slopes' signal: leaving it out biases the slopes, fitting it removes the bias
	expect_lt(max(abs(b_uv)), 12)
	expect_gt(min(abs(b_0)), 12)
	for (f in fits) {
		ob = array(FALSE, dim(f$d$truth$O))
		for (t in seq_along(f$d$Y)) ob[, , t] = is.finite(f$d$Y[[t]])
		expect_gt(stats::cor(f$uv$Oarr[ob], f$d$truth$O[ob]), 0.6)
		expect_gt(f$uv$rho_uv, 0.6)
		expect_lt(f$uv$rho_uv, 0.99)
	}
})

test_that("moving positions: the reported object is complete on directed and bipartite panels", {
	skip_on_cran()
	d = vsim(n = 20, Tt = 4, seed = 61, mu = -1, R = 2, dyn_uv = TRUE, sigma_uv = 0.4)
	f = eb_run(d, R = 2, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_equal(dim(f$U), c(20L, 2L, 4L))
	expect_equal(dim(f$V), c(20L, 2L, 4L))
	expect_equal(dim(f$U_sd), dim(f$U))
	expect_true(all(f$U_sd > 0))
	for (t in 1:4) expect_equal(f$Oarr[, , t], tcrossprod(f$U[, , t], f$V[, , t]), ignore_attr = TRUE)
	expect_equal(f$UVPM, f$Oarr)
	expect_identical(f$dynamic_uv_kind, "ar1")
	expect_true("positions" %in% names(f$diagnostics$stationarity$by_block))
	expect_true(f$diagnostics$stationarity$pass)
	expect_true(is.finite(f$eb$uv_scale) && f$eb$uv_scale > 0.08 && f$eb$uv_scale < 4)
	expect_true(is.data.frame(f$eb$uv_search) && nrow(f$eb$uv_search) >= 4L)
	# the linear predictor carries the product, the fitted probabilities follow from it
	expect_equal(
		f$EZ[[2]][3, 7],
		f$Oarr[3, 7, 2] + f$a_dynamic[3, 2] + f$b_dynamic[7, 2] +
			sum(f$coef_path[, 2] * c(1, d$X[[2]][3, 7, ])),
		tolerance = 1e-8
	)
	expect_equal(f$fitted[[2]][3, 7], stats::pnorm(f$EZ[[2]][3, 7]), tolerance = 1e-10)
	# the generic methods read it
	msgs = capture.output(print(f))
	expect_true(any(grepl("Moving latent positions", msgs)))
	expect_s3_class(summary(f), "summary.lame_dynamic_als")
	expect_equal(length(coef(f)), 3L)
	expect_equal(dim(vcov(f)), c(3L, 3L))
	expect_equal(length(residuals(f)), 4L)
	lp = latent_positions(f)
	expect_s3_class(lp, "data.frame")
	expect_equal(nrow(lp), 2L * 20L * 2L * 4L)
	expect_true(all(lp$posterior_sd > 0))
	# a bipartite panel; the interaction has rank min(R_row, R_col) and the larger factor carries zero columns
	db = vsim(n = 24, nb = 16, Tt = 4, mode = "bipartite", seed = 62, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	fb = eb_run(db, mode = "bipartite", R = 1, dynamic_uv = TRUE)
	expect_true(fb$converged)
	expect_equal(dim(fb$U), c(24L, 1L, 4L))
	expect_equal(dim(fb$V), c(16L, 1L, 4L))
	for (t in 1:4) expect_equal(fb$Oarr[, , t], tcrossprod(fb$U[, , t], fb$V[, , t]), ignore_attr = TRUE)
	f2 = suppressMessages(eb_run(db, mode = "bipartite", R_row = 2, R_col = 1, dynamic_uv = TRUE))
	expect_equal(dim(f2$U), c(24L, 2L, 4L))
	expect_equal(dim(f2$V), c(16L, 1L, 4L))
	expect_true(all(f2$U[, 2, ] == 0))
	expect_equal(f2$coefficients, fb$coefficients, tolerance = 1e-6)
})

test_that("moving positions: missing cells and a changing composition are handled", {
	skip_on_cran()
	d = vsim(n = 18, Tt = 5, seed = 63, mu = -1, R = 1, dyn_uv = TRUE, sigma_uv = 0.4)
	# an actor absent in two periods, scattered missing cells
	for (t in c(2, 3)) { d$Y[[t]][4, ] = NA; d$Y[[t]][, 4] = NA }
	set.seed(1)
	for (t in 1:5) d$Y[[t]][sample(which(is.finite(d$Y[[t]])), 12)] = NA
	f = eb_run(d, R = 1, dynamic_uv = TRUE)
	expect_true(f$converged)
	expect_true(all(is.finite(f$U)) && all(is.finite(f$V)))
	expect_true(all(is.finite(f$coefficients)))
	expect_equal(dim(f$Oarr), c(18L, 18L, 5L))
})

test_that("a sparse 100-actor panel on which the full fixed-point step cycles converges", {
	# the full fixed-point step overshoots on this panel (an eigenvalue below -1) and the half step cycles
	# the relaxation halves the step while successive moves reverse
	skip_on_cran()
	d = vsim(n = 100, Tt = 10, seed = 1002, mu = -3)
	f = eb_run(d)
	expect_true(isTRUE(f$converged))
	expect_lt(f$iterations, 120)
	expect_equal(unname(f$coefficients[c("x1_dyad", "x2_dyad")]), c(0.6, -0.5), tolerance = 0.1)
	# the trace shows the relaxation at work
	tr = f$diagnostics$convergence_trace
	expect_lt(min(tr$step, na.rm = TRUE), 0.5)
	expect_gt(utils::tail(tr$step, 1), 0.5)
})

test_that("strongly identified small panels settle: the quadrature follows the variance of the product, and the relaxation keeps its cap", {
	# nearly deterministic cells: the product of the positions has a large variance, which the quadrature has to follow
	# the cells are independent (dcor = FALSE) since the subject is the quadrature and the step controller
	skip_on_cran()
	pos_sim = function(seed, n = 10, Tt = 4, mu = -0.5, rho_uv = 0.95, sigma_uv = 1) {
		set.seed(seed)
		ar1 = function(m) {
			out = matrix(0, m, Tt)
			out[, 1] = rnorm(m, 0, sigma_uv / sqrt(1 - rho_uv^2))
			for (t in 2:Tt) out[, t] = rho_uv * out[, t - 1] + rnorm(m, 0, sigma_uv)
			out
		}
		u = ar1(n); v = ar1(n)
		rn = paste0("a", sprintf("%02d", seq_len(n)))
		Y = lapply(seq_len(Tt), function(t) {
			y = 1 * (mu + tcrossprod(u[, t], v[, t]) + matrix(rnorm(n * n), n) > 0)
			dimnames(y) = list(rn, rn)
			diag(y) = NA
			y
		})
		names(Y) = paste0("t", seq_len(Tt))
		Y
	}
	for (seed in c(3, 8, 11)) {
		Y = pos_sim(seed)
		expect_no_warning(
			f <- suppressMessages(lame(
				Y,
				family = "binary",
				method = "als",
				R = 1,
				dynamic_ab = TRUE,
				dynamic_uv = TRUE,
				dcor = FALSE,
				verbose = FALSE
			))
		)
		expect_true(isTRUE(f$converged), info = paste("seed", seed))
		# the positions are large and the stationary variance of the order of the generating one (10)
		expect_gt(f$eb$uv_scale, 3)
		expect_lt(f$eb$uv_scale, 30)
	}
})

test_that("first-order persistences pool through a logit-normal law and a random walk has none to pool", {
	skip_on_cran()
	d = vsim(n = 30, Tt = 8, seed = 51, mu = -1.2, beta_sd = 0.15)
	args = list(dynamic_beta = c("intercept", "x1_dyad"))
	fit = function(...) {
		suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			dynamic_ab = TRUE,
			verbose = FALSE,
			...
		)))
	}
	fa = do.call(fit, c(args, list(dynamic_beta_kind = "ar1", dynamic_beta_pool = "both")))
	expect_true(fa$converged)
	expect_true(is.finite(fa$eb$beta_pool$rho_mu) && fa$eb$beta_pool$rho_tau > 0)
	expect_true(is.finite(fa$eb$beta_pool$sigma_share) && fa$eb$beta_pool$sigma_share > 0)
	expect_true(all(fa$rho_beta > 0 & fa$rho_beta < 0.999))
	# the persistence-only pool leaves the scales on their own priors
	fr = do.call(fit, c(args, list(dynamic_beta_kind = "ar1", dynamic_beta_pool = "rho")))
	expect_null(fr$eb$beta_pool$sigma_share)
	expect_true(is.finite(fr$eb$beta_pool$rho_mu))
	# a random walk has no persistence: "both" is "sigma"
	fs = do.call(fit, c(args, list(dynamic_beta_kind = "rw1", dynamic_beta_pool = "sigma")))
	fb = do.call(fit, c(args, list(dynamic_beta_kind = "rw1", dynamic_beta_pool = "both")))
	expect_equal(fb$sigma_beta, fs$sigma_beta)
	expect_null(fb$eb$beta_pool$rho_mu)
})

test_that("a disabled margin recovers a panel generated without those effects", {
	skip_on_cran()
	# the fitter called directly, which takes the switch whether or not lame() forwards it
	eb_direct = function(d, ...) {
		suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
			d$Y,
			d$X,
			family = "binary",
			dynamic_ab = TRUE,
			verbose = FALSE,
			...
		)))
	}
	fits = lapply(1:8, function(r) {
		d = vsim(n = 30, Tt = 5, seed = 600 + r, mu = -1.2, sender = FALSE)
		list(off = eb_direct(d, rvar = FALSE), on = eb_direct(d))
	})
	expect_true(all(vapply(fits, function(z) isTRUE(z$off$converged), logical(1))))
	expect_true(all(vapply(fits, function(z) all(z$off$a_dynamic == 0), logical(1))))
	off = lapply(fits, `[[`, "off"); on = lapply(fits, `[[`, "on")
	expect_lt(max(abs(slope_bias(off))), 5)
	# the switch is information: the standard errors are no larger than with the
	# sender effects estimated
	se_off = rowMeans(sapply(off, function(f) f$eb$se[c("x1_dyad", "x2_dyad")]))
	se_on = rowMeans(sapply(on, function(f) f$eb$se[c("x1_dyad", "x2_dyad")]))
	expect_true(all(se_off <= se_on * 1.001))
})

test_that("the front door forwards the sampler switches where lame() has been taught to", {
	skip_on_cran()
	# lame() hands rvar, cvar, nvar, dcor and dynamic_beta_pool to the dynamic route when its dispatcher forwards them
	forwards = grepl(
		'"rvar", "cvar", "nvar", "dcor", "dynamic_beta_pool"',
		paste(deparse(body(lame), width.cutoff = 500L), collapse = " "),
		fixed = TRUE
	)
	skip_if_not(forwards, "lame() does not forward the sampler switches yet")
	d = vsim(n = 14, Tt = 4, seed = 61, mu = -1.2)
	f1 = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		method = "als",
		dynamic_ab = TRUE,
		rvar = FALSE,
		verbose = FALSE
	))
	expect_true(all(f1$a_dynamic == 0))
	f2 = eb_run(d, rvar = FALSE)
	expect_equal(f1$coefficients, f2$coefficients)
	# dcor is honoured by the empirical-Bayes route (no warning, a record of the pair likelihood, or none when it is off) and
	# named as ignored by the penalised route
	catch = function(expr) {
		w = character()
		val = withCallingHandlers(
			expr,
			warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
		)
		list(val = val, w = w)
	}
	f3 = catch(
		lame(d$Y, Xdyad = d$X, family = "binary", method = "als", dynamic_ab = TRUE, dcor = FALSE, verbose = FALSE)
	)
	expect_false(any(grepl("ignored MCMC-only", f3$w)))
	expect_null(f3$val$eb$pair)
	f4 = catch(
		lame(d$Y, Xdyad = d$X, family = "binary", method = "als", dynamic_ab = TRUE, dcor = TRUE, verbose = FALSE)
	)
	expect_false(any(grepl("ignored MCMC-only", f4$w)))
	expect_false(is.null(f4$val$eb$pair))
	f5 = catch(
		lame(
			d$Y,
			Xdyad = d$X,
			family = "binary",
			method = "als",
			dynamic_ab = TRUE,
			dcor = FALSE,
			bias_correct = "none",
			verbose = FALSE
		)
	)
	expect_true(any(grepl("ignored MCMC-only arguments.*dcor", f5$w)))
})


# the dyadic correlation: with static effects the dynamic engine takes the pair likelihood by its own route (R/ame_als_pair.R)
# and agrees with the static estimator on the correlation, coefficients and standard errors
test_that("with static effects the dynamic engine with the dyadic correlation is the static estimator's pair fit", {
	skip_on_cran()
	for (cf in list(
		list(n = 40, Tt = 3, seed = 61, dc = 0.6),
		list(n = 40, Tt = 3, seed = 62, dc = -0.4),
		list(n = 60, Tt = 2, seed = 63, dc = 0.3)
	)) {
		d = vsim(n = cf$n, Tt = cf$Tt, seed = cf$seed, mu = -1.2, va = 0.5, vb = 0.5, sigma = 0.001, dcor = cf$dc)
		prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
		ctx = lame:::.lde_context(prep, FALSE, "probit", coef_prior = "flat", dynamic_ab = FALSE, R = 0L, dcor = TRUE)
		fit = lame:::.lde_fit(ctx, max_iter = 600, tol = 1e-6)
		info = paste("truth", cf$dc)
		expect_true(fit$converged, info = info)
		th = fit$st$theta
		th[1] = th[1] - sum(ctx$centre * th)
		q = fit$ctx$q; m = fit$ctx$mth
		se = sqrt(diag(lame:::.lde_report_vcov(fit$ctx, fit$ev$fac$Szi[q + seq_len(m), q + seq_len(m), drop = FALSE])))
		st = suppressWarnings(suppressMessages(lame:::.ae_run(
			d$Y,
			d$X,
			NULL,
			NULL,
			0L,
			"binary",
			"unipartite",
			FALSE,
			600,
			1e-8,
			FALSE,
			6886,
			longitudinal = TRUE,
			call = quote(lame_als()),
			lowrank_method = "mm",
			non_normal_method = "irls",
			link = "probit",
			linear_solver = "eigen",
			multistart = "none",
			bias_correct = NULL,
			dcor = TRUE,
			dcor_given = TRUE
		)))
		expect_false(is.null(st$eb$pair))
		expect_lt(abs(fit$st$hyp$rho_d - unname(st$VC["rho"])), 0.02, label = paste("rho", info))
		expect_lt(max(abs(th - st$coefficients) / st$eb$se), 0.15, label = paste("z", info))
		expect_true(all(se / st$eb$se > 0.96 & se / st$eb$se < 1.04), info = info)
		# and the correlation is the panel's: within three posterior standard deviations of the truth
		expect_lt(abs(fit$st$hyp$rho_d - cf$dc), 3 * st$eb$pair$sd, label = paste("truth", info))
	}
})

# static effects and factors against the static estimator on panels where the factors are weakly informed
# the cavity rule for the scale of the factor prior decides these fits, and the fixed point does not depend on the start
fac_panels = list(
	list(mode = "directed", n = 40, nb = 40, Tt = 2, R = 2, suv = 0.25, mu = -2, seed = 103),
	list(mode = "directed", n = 24, nb = 24, Tt = 3, R = 1, suv = 0.15, mu = -1.5, seed = 105),
	list(mode = "symmetric", n = 40, nb = 40, Tt = 2, R = 2, suv = 0.25, mu = -2, seed = 107),
	list(mode = "bipartite", n = 30, nb = 18, Tt = 4, R = 1, suv = 0.6, mu = -1.2, seed = 109),
	list(mode = "bipartite", n = 24, nb = 16, Tt = 3, R = 1, suv = 0.15, mu = -1.5, seed = 111)
)

test_that("weakly identified static factors: the dynamic engine tracks the static estimator, and the start does not matter", {
	skip_on_cran()
	for (cf in fac_panels) {
		d = vsim(
			n = cf$n,
			nb = cf$nb,
			Tt = cf$Tt,
			mode = cf$mode,
			seed = cf$seed,
			mu = cf$mu,
			va = 0.5,
			vb = 0.5,
			sigma = 0.001,
			R = cf$R,
			suv = cf$suv
		)
		sym = identical(cf$mode, "symmetric")
		md = if (identical(cf$mode, "bipartite")) "bipartite" else "unipartite"
		prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, md, longitudinal = TRUE, node_time_policy = "quiet")
		ctx = lame:::.lde_context(prep, sym, "probit", coef_prior = "flat", dynamic_ab = FALSE, R = as.integer(cf$R))
		fit = lame:::.lde_fit(ctx, max_iter = 600, tol = 1e-6)
		expect_true(fit$converged, info = paste(cf$mode, cf$seed))
		th = fit$st$theta
		th[1] = th[1] - sum(ctx$centre * th)
		se = sqrt(diag(lame:::.lde_fac_vcov(fit$ctx, fit$st, fit$ev)$vcov))
		st = suppressWarnings(suppressMessages(lame:::.ae_run(
			d$Y,
			d$X,
			NULL,
			NULL,
			as.integer(cf$R),
			"binary",
			md,
			sym,
			600,
			1e-8,
			FALSE,
			6886,
			longitudinal = TRUE,
			call = quote(lame_als()),
			lowrank_method = "mm",
			non_normal_method = "irls",
			link = "probit",
			linear_solver = "eigen",
			multistart = "none",
			bias_correct = NULL
		)))
		info = paste(cf$mode, cf$seed)
		# coefficients within a third of a posterior standard deviation, standard errors within 12 percent
		expect_lt(max(abs(th - st$coefficients) / st$eb$se), 0.35, label = paste("z", info))
		expect_true(all(se / st$eb$se > 0.88 & se / st$eb$se < 1.12), info = info)
		# the start of the iteration does not matter over a range of sizes (0.3 to 1 of the least-squares fit)
		for (sc in c(0.35, 0.9)) {
			ctx2 = ctx
			ctx2$fac_start = sc
			fit2 = lame:::.lde_fit(ctx2, max_iter = 600, tol = 1e-6)
			expect_true(fit2$converged, info = paste(info, sc))
			expect_equal(unname(fit2$st$theta), unname(fit$st$theta), tolerance = 2e-3, info = paste(info, sc))
			expect_equal(
				mean(lame:::.lde_fac_scales(fit2$st$hyp$fac)),
				mean(lame:::.lde_fac_scales(fit$st$hyp$fac)),
				tolerance = 0.02,
				info = paste(info, sc)
			)
		}
	}
})

test_that("static factors with the logit link are the static estimator's logit fit", {
	skip_on_cran()
	d = vsim(n = 30, Tt = 4, seed = 7, mu = -1.2, va = 0.5, vb = 0.5, sigma = 0.001, R = 1, suv = 0.6)
	prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", longitudinal = TRUE, node_time_policy = "quiet")
	ctx = lame:::.lde_context(prep, FALSE, "logit", coef_prior = "flat", dynamic_ab = FALSE, R = 1L)
	fit = lame:::.lde_fit(ctx, max_iter = 600, tol = 1e-6)
	expect_true(fit$converged)
	th = fit$st$theta
	th[1] = th[1] - sum(ctx$centre * th)
	st = suppressWarnings(suppressMessages(lame:::.ae_run(
		d$Y,
		d$X,
		NULL,
		NULL,
		1L,
		"binary",
		"unipartite",
		FALSE,
		600,
		1e-8,
		FALSE,
		6886,
		longitudinal = TRUE,
		call = quote(lame_als()),
		lowrank_method = "mm",
		non_normal_method = "irls",
		link = "logit",
		linear_solver = "eigen",
		multistart = "none",
		bias_correct = NULL
	)))
	expect_lt(max(abs(th - st$coefficients) / st$eb$se), 0.5)
})
