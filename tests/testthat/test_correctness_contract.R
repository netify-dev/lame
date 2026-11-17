# correctness contracts that check the numbers, not only that a fit runs, across method, mode and family
# fits are small for speed; the assertions are exact where the math is exact
skip_on_cran()

# helpers

.mk_uni_ame = function(seed = 1, n = 16, nodal = TRUE, family = "binary") {
	set.seed(seed)
	Xd = array(rnorm(n * n * 2), c(n, n, 2), dimnames = list(NULL, NULL, c("d1", "d2")))
	Y = if (family == "normal") {
		0.5 + 0.8 * Xd[, , 1] - 0.6 * Xd[, , 2] + matrix(rnorm(n * n, 0, 0.5), n, n)
	} else {
		eta = -0.3 + 0.9 * Xd[, , 1] - 0.7 * Xd[, , 2]
		matrix(rbinom(n * n, 1, pnorm(eta)), n, n)
	}
	diag(Y) = NA; rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	args = list(
		Y = Y,
		Xdyad = Xd,
		family = family,
		R = 1,
		burn = 60,
		nscan = 400,
		odens = 5,
		verbose = FALSE
	)
	if (nodal) {
		args$Xrow = matrix(rnorm(n), n, 1, dimnames = list(NULL, "r1"))
		args$Xcol = matrix(rnorm(n), n, 1, dimnames = list(NULL, "c1"))
	}
	list(fit = do.call(ame, args), Y = Y, Xd = Xd)
}

# coef and confint are the posterior summaries

test_that("C1/C2: coef == posterior mean and confint == posterior quantiles", {
	obj = .mk_uni_ame(); fit = obj$fit
	expect_equal(unname(coef(fit)), unname(colMeans(fit$BETA)), tolerance = 1e-10)
	ci = confint(fit)
	q = t(apply(fit$BETA, 2, quantile, c(0.025, 0.975)))
	expect_equal(unname(as.matrix(ci)), unname(q), tolerance = 1e-8)
})

# fitted, predict and residuals are mutually consistent

test_that("C3/C4: fitted == predict(response) and residuals == Y - fitted", {
	for (fam in c("binary", "normal")) {
		obj = .mk_uni_ame(family = fam); fit = obj$fit
		fv = fitted(fit)
		pr = predict(fit, type = "response")
		expect_equal(
			as.numeric(fv),
			as.numeric(pr),
			tolerance = 1e-8,
			info = paste("fitted==predict", fam)
		)
		rs = residuals(fit)
		# residuals defined as observed - fitted on observed cells
		obs = obj$Y
		expect_equal(
			as.numeric(rs[!is.na(obs)]),
			as.numeric((obs - fv)[!is.na(obs)]),
			tolerance = 1e-6,
			info = paste("resid==Y-fitted", fam)
		)
	}
})

# predict on the fit's own covariates reproduces the in-sample linear predictor

test_that("C5: predict(newdata = fit's own Xdyad) reproduces the in-sample EZ", {
	# with nodal covariates present
	obj = .mk_uni_ame(nodal = TRUE); fit = obj$fit
	ez_in = predict(fit, type = "link")
	ez_nd = predict(fit, type = "link", newdata = obj$Xd)
	expect_lt(max(abs(ez_in - ez_nd), na.rm = TRUE), 1e-8)
	# bipartite too
	set.seed(3); nA = 12; nB = 9
	Yb = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
	rownames(Yb) = paste0("a", seq_len(nA)); colnames(Yb) = paste0("b", seq_len(nB))
	Xb = array(rnorm(nA * nB * 2), c(nA, nB, 2), dimnames = list(NULL, NULL, c("d1", "d2")))
	fb = ame(
		Yb,
		Xdyad = Xb,
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		family = "binary",
		burn = 50,
		nscan = 300,
		odens = 5,
		verbose = FALSE
	)
	expect_lt(
		max(
			abs(predict(fb, type = "link") -
					predict(fb, type = "link", newdata = Xb)),
			na.rm = TRUE
		),
		1e-8
	)
})

# a shift of delta in a covariate moves the link by delta times its own coefficient

test_that("C6: +delta on a covariate moves the link by delta * its own coef", {
	obj = .mk_uni_ame(nodal = TRUE); fit = obj$fit
	bm = colMeans(fit$BETA)
	base = predict(fit, type = "link", newdata = obj$Xd)
	for (slice in 1:2) {
		nm = paste0("d", slice, "_dyad")
		Xup = obj$Xd; Xup[, , slice] = Xup[, , slice] + 2.5
		up = predict(fit, type = "link", newdata = Xup)
		expect_equal(
			mean(up - base, na.rm = TRUE),
			unname(2.5 * bm[nm]),
			tolerance = 1e-6,
			info = nm
		)
	}
})

# known-dgp recovery: the truth lies inside the 95% ci

test_that("C7: a strong dyadic effect is recovered inside its 95% CI", {
	set.seed(11); n = 30
	Xd = array(rnorm(n * n), c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	eta = -0.5 + 1.0 * Xd[, , 1] + outer(a, b, "+")
	Y = matrix(rbinom(n * n, 1, pnorm(eta)), n, n); diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	fit = ame(
		Y,
		Xdyad = Xd,
		family = "binary",
		R = 0,
		burn = 200,
		nscan = 1500,
		odens = 5,
		verbose = FALSE
	)
	ci = confint(fit)["x1_dyad", ]
	expect_true(ci[1] <= 1.0 && 1.0 <= ci[2]) # truth 1.0 inside the ci
	expect_gt(coef(fit)["x1_dyad"], 0.5) # right sign and ballpark
})

# ame() and lame() agree on the default g-prior

test_that("C8: ame() and lame() use the same default g for the same data", {
	# normal family with a large-magnitude outcome: both scale g by the variance and recover the slope
	set.seed(5); n = 20
	Xd = matrix(rnorm(n * n), n, n)
	Y = 300 * (0.7 * Xd) + matrix(rnorm(n * n, 0, 150), n, n)
	diag(Y) = NA; rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	Xa = array(Xd, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	f_ame = ame(
		Y,
		Xdyad = Xa,
		family = "normal",
		R = 0,
		burn = 100,
		nscan = 600,
		odens = 10,
		verbose = FALSE
	)
	# lame with a single-period list reduces to the same cross-section
	f_lame = lame(
		list(t1 = Y),
		Xdyad = list(t1 = Xa),
		family = "normal",
		R = 0,
		burn = 100,
		nscan = 600,
		odens = 10,
		verbose = FALSE
	)
	truth = 300 * 0.7
	# a single-period lame() fit is static, so coef() is a vector (a dynamic_beta fit gives a p x t matrix);
	# index by name only
	cl = coef(f_lame)
	lame_dyad = if (is.matrix(cl)) cl["x1_dyad", 1] else cl["x1_dyad"]
	expect_gt(coef(f_ame)["x1_dyad"], 0.5 * truth) # ame recovers
	expect_gt(lame_dyad, 0.5 * truth) # lame recovers
})

# covariates are stored and reused by every method

test_that("a binary fit keeps its covariates and its link: predict(link) == reconstruct_EZ and YPM ~ pnorm(EZ)", {
	set.seed(6886)
	n = 15

	X = matrix(rnorm(n * n), n, n)
	diag(X) = NA

	# binary probit data with a dyadic covariate
	eta_true = -0.5 + 0.8 * X
	p_true = pnorm(eta_true)
	Y = matrix(NA, n, n)
	non_na_idx = which(!is.na(p_true))
	Y[non_na_idx] = rbinom(length(non_na_idx), 1, p_true[non_na_idx])
	diag(Y) = NA

	fit = ame(
		Y,
		Xdyad = X,
		R = 0,
		family = "binary",
		burn = 200,
		nscan = 500,
		verbose = FALSE,
		seed = 456
	)

	expect_false(is.null(fit$X), label = "X should be stored in fit object")
	expect_equal(dim(fit$X), c(n, n, 2))

	EZ = reconstruct_EZ(fit)
	expect_false(is.null(EZ))

	# reconstructed link tracks the data-generating link
	corr = cor(c(eta_true[!is.na(Y)]), c(EZ[!is.na(Y)]))
	expect_gt(corr, 0.5)

	pred_link = predict(fit, type = "link")
	expect_equal(
		pred_link,
		EZ,
		tolerance = 1e-10,
		info = "predict(type='link') should match reconstruct_EZ"
	)

	pred_resp = predict(fit, type = "response")
	expect_true(
		all(pred_resp >= 0 & pred_resp <= 1, na.rm = TRUE),
		info = "Binary predictions should be probabilities"
	)

	expected_probs = pnorm(EZ)
	mean_diff = mean(abs(fit$YPM - expected_probs), na.rm = TRUE)
	expect_lt(mean_diff, 0.1, label = "YPM should approximate Φ(EZ)")
})

# the contracts on bipartite and longitudinal fits and the cross-function checks, on small normal fits
# built on first use

.cc_cache = new.env()
.cc_build = list(
	ame_uni = function() {
		n = 15
		dat = simulate_test_network(
			n = n,
			family = "normal",
			R = 2,
			beta_intercept = 1,
			beta_dyad = 0.5,
			mode = "unipartite",
			seed = 8001
		)
		ame(
			dat$Y,
			Xdyad = dat$Xdyad,
			R = 2,
			family = "normal",
			burn = 100,
			nscan = 200,
			odens = 1,
			verbose = FALSE,
			gof = TRUE,
			seed = 42
		)
	},
	ame_bip = function() {
		nA = 12; nB = 8
		dat = simulate_test_network(
			n = nA,
			nA = nA,
			nB = nB,
			family = "normal",
			R = 2,
			beta_intercept = 1,
			beta_dyad = 0.5,
			mode = "bipartite",
			seed = 8002
		)
		ame(
			dat$Y,
			Xdyad = dat$Xdyad,
			R = 2,
			family = "normal",
			mode = "bipartite",
			burn = 100,
			nscan = 200,
			odens = 1,
			verbose = FALSE,
			gof = TRUE,
			seed = 42
		)
	},
	lame_uni = function() {
		n = 15
		dat = simulate_test_network(
			n = n,
			n_time = 3,
			family = "normal",
			R = 2,
			beta_intercept = 1,
			beta_dyad = 0.5,
			mode = "unipartite",
			seed = 8003
		)
		lame(
			dat$Y,
			Xdyad = dat$Xdyad,
			R = 2,
			family = "normal",
			burn = 100,
			nscan = 200,
			odens = 1,
			verbose = FALSE,
			gof = TRUE,
			seed = 42,
			posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)
		)
	},
	lame_bip = function() {
		nA = 12; nB = 8
		dat = simulate_test_network(
			n = nA,
			nA = nA,
			nB = nB,
			n_time = 3,
			family = "normal",
			R = 2,
			beta_intercept = 1,
			beta_dyad = 0.5,
			mode = "bipartite",
			seed = 8004
		)
		Xdyad_list = lapply(dat$Xdyad, function(x) x[, , 1])
		lame(
			dat$Y,
			Xdyad = Xdyad_list,
			R = 2,
			family = "normal",
			mode = "bipartite",
			burn = 100,
			nscan = 200,
			odens = 1,
			verbose = FALSE,
			gof = TRUE,
			seed = 42
		)
	}
)
.cc_fit = function(kind) {
	if (is.null(.cc_cache[[kind]])) assign(kind, .cc_build[[kind]](), envir = .cc_cache)
	.cc_cache[[kind]]
}

test_that("C1/C2 for bipartite and longitudinal fits: coef == posterior mean and confint == posterior quantiles", {
	for (kind in c("ame_bip", "lame_uni")) {
		fit = .cc_fit(kind)
		expect_equal(unname(coef(fit)), unname(colMeans(fit$BETA)), tolerance = 1e-10, info = kind)
		q = t(apply(fit$BETA, 2, quantile, c(0.025, 0.975)))
		expect_equal(unname(as.matrix(confint(fit))), unname(q), tolerance = 1e-8, info = kind)
	}
})

test_that("C3/C4 for bipartite and longitudinal fits: fitted == predict(response) and residuals == Y - fitted", {
	fit = .cc_fit("ame_bip")
	expect_equal(as.numeric(predict(fit, type = "response")), as.numeric(fitted(fit)), tolerance = 1e-8)
	expect_equal(residuals(fit, type = "response"), fit$Y - fitted(fit))
	for (kind in c("lame_uni", "lame_bip")) {
		fit = .cc_fit(kind)
		pr = predict(fit, type = "response")
		fv = fitted(fit)
		expect_equal(length(pr), length(fv))
		for (t in seq_along(pr)) expect_equal(as.numeric(pr[[t]]), as.numeric(fv[[t]]), tolerance = 1e-8, info = kind)
	}
	fit = .cc_fit("lame_uni")
	r = suppressWarnings(residuals(fit, type = "response"))
	fv = fitted(fit)
	for (t in seq_along(r)) expect_equal(r[[t]], fit$Y[, , t] - fv[[t]])
})

test_that("diag(vcov) consistent with confint width for ame unipartite", {
	v = vcov(.cc_fit("ame_uni"))
	ci = confint(.cc_fit("ame_uni"), level = 0.95)
	# interval width is about 2 * 1.96 * sqrt(diag(vcov)) for a near-normal posterior
	ci_width = ci[, 2] - ci[, 1]
	approx_width = 2 * 1.96 * sqrt(diag(v))
	# roughly proportional, within a factor of 2
	ratio = ci_width / approx_width
	expect_true(
		all(ratio > 0.3 & ratio < 3),
		info = paste("Width ratios:", paste(round(ratio, 2), collapse = ", "))
	)
})

# simulate mean vs predict

test_that("simulate mean approximates predict for ame normal", {
	# for the normal family the mean of many simulations approaches predict
	sim_obj = simulate(.cc_fit("ame_uni"), nsim = 50)
	sims = sim_obj$Y
	sim_mean = Reduce("+", sims) / length(sims)
	pred = predict(.cc_fit("ame_uni"), type = "response")
	diff = abs(sim_mean - pred)
	diff = diff[!is.na(diff)]
	# mean absolute difference is moderate (generous tolerance for MCMC)
	expect_true(mean(diff) < 2.0)
})

# reconstruct_EZ vs stored EZ (lame only)

test_that("reconstruct_EZ matches stored EZ for lame unipartite", {
	skip_if(is.null(.cc_fit("lame_uni")$EZ), "EZ not stored in fit object")
	recon = reconstruct_EZ(.cc_fit("lame_uni"))
	stored = .cc_fit("lame_uni")$EZ
	expect_equal(length(recon), length(stored))
	for (t in seq_along(recon)) {
		diff = abs(recon[[t]] - stored[[t]])
		diff = diff[!is.na(diff)]
		expect_true(
			max(diff) < 0.1,
			info = paste("Max EZ diff at time", t, ":", max(diff))
		)
	}
})

test_that("reconstruct_EZ matches stored EZ for lame bipartite", {
	skip_if(is.null(.cc_fit("lame_bip")$EZ), "EZ not stored in fit object")
	recon = reconstruct_EZ(.cc_fit("lame_bip"))
	stored = .cc_fit("lame_bip")$EZ
	expect_equal(length(recon), length(stored))
	for (t in seq_along(recon)) {
		diff = abs(recon[[t]] - stored[[t]])
		diff = diff[!is.na(diff)]
		expect_true(max(diff) < 0.1)
	}
})

# gof() on lame objects

test_that("gof() works on lame objects", {
	result = gof(.cc_fit("lame_uni"), nsim = 5, verbose = FALSE)
	expect_true(is.matrix(result))
	expect_equal(nrow(result), 6) # 1 obs + 5 sims
	expect_true("sd.rowmean" %in% colnames(result))
	expect_true(all(is.finite(result[1, ])))
})
