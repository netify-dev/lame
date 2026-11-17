# h-step forecasts propagate the multiplicative term as the posterior mean of U V' decayed by rho_uv^(2h),
# not as a product of posterior-mean coordinates
skip_on_cran()

# freeze the AR(1) hyperparameters (rho = 1, sigma = 0) so the forecast carries the last-period state forward unchanged
.freeze_dynamics = function(fit, rho_uv = 1, rho_ab = 1) {
	if (!is.null(fit$rho_uv)) fit$rho_uv = rep(rho_uv, length(fit$rho_uv))
	if (!is.null(fit$sigma_uv)) fit$sigma_uv = rep(0, length(fit$sigma_uv))
	if (!is.null(fit$rho_ab)) fit$rho_ab = rep(rho_ab, length(fit$rho_ab))
	if (!is.null(fit$sigma_ab)) fit$sigma_ab = rep(0, length(fit$sigma_ab))
	# freeze the coefficient paths too, the dynamic_ab intercept among them
	if (!is.null(fit$RHO_BETA)) fit$RHO_BETA[] = 1
	if (!is.null(fit$SIGMA_BETA)) fit$SIGMA_BETA[] = 0
	fit
}

# the stored EZ has an NA diagonal and the forecast fills every cell, so comparisons drop the diagonal
.offdiag = function(m) m[row(m) != col(m)]

.fc_panel = function(
	n = 12,
	TT = 4,
	seed = 11,
	symmetric = FALSE,
	nA = NULL,
	nB = NULL
) {
	set.seed(seed)
	if (!is.null(nA)) {
		rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
		X = lapply(seq_len(TT), function(t) {
			arr = array(rnorm(nA * nB), c(nA, nB, 1))
			dimnames(arr) = list(rn, cn, "x1"); arr })
		Y = lapply(seq_len(TT), function(t) {
			m = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
			dimnames(m) = list(rn, cn); m })
		return(list(Y = Y, X = X, TT = TT))
	}
	nm = sprintf("a%02d", seq_len(n))
	X = lapply(seq_len(TT), function(t) {
		x = matrix(rnorm(n * n), n, n)
		# a symmetric fit needs a symmetric covariate
		if (symmetric) x = (x + t(x)) / 2
		arr = array(x, c(n, n, 1))
		dimnames(arr) = list(nm, nm, "x1"); arr })
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n)
		if (symmetric) m[lower.tri(m)] = t(m)[lower.tri(m)]
		diag(m) = NA; dimnames(m) = list(nm, nm); m })
	list(Y = Y, X = X, TT = TT, n = n)
}

.fc_fit_dynamic_uv = function() {
	d = .fc_panel()
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 2,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 5L
	))
	list(fit = fit, d = d)
}

test_that("the h = 1 multiplicative term is rho_uv^2 times the posterior-mean product", {
	skip_on_cran()
	obj = .fc_fit_dynamic_uv(); fit = obj$fit; TT = obj$d$TT
	rho = 0.6
	frozen = .freeze_dynamics(fit, rho_uv = rho)
	fc = suppressWarnings(predict(frozen, h = 1L, type = "link"))[[1L]]
	# all but the multiplicative term match the period-T linear predictor, isolating the uv contribution
	expect_equal(
		.offdiag(fc) - .offdiag(fit$EZ[[TT]]),
		.offdiag((rho^2 - 1) * fit$UVPM[, , TT]),
		tolerance = 1e-10
	)
	# at h = 2 the decay is rho^4
	fc2 = suppressWarnings(predict(frozen, h = 2L, type = "link"))[[2L]]
	expect_equal(
		.offdiag(fc2) - .offdiag(fit$EZ[[TT]]),
		.offdiag((rho^4 - 1) * fit$UVPM[, , TT]),
		tolerance = 1e-10
	)
})

test_that("symmetric forecasts carry the ULU' term at its posterior mean", {
	skip_on_cran()
	d = .fc_panel(symmetric = TRUE, seed = 3)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 2,
		symmetric = TRUE,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 5L
	))
	# symmetric fits store ULUPM and no V
	expect_null(fit[["V"]])
	expect_false(is.null(fit$ULUPM))
	frozen = .freeze_dynamics(fit)
	fc = suppressWarnings(predict(frozen, h = 1L, type = "link"))[[1L]]
	bm = if (length(dim(fit$BETA)) == 3L) {
		colMeans(fit$BETA[, , dim(fit$BETA)[3]])
	} else {
		colMeans(fit$BETA)
	}
	# a symmetric fit stores one effect path for both margins
	a_T = fit$a_dynamic[, d$TT]
	# each unordered dyad uses the mean of the two directed design entries
	x_T = d$X[[d$TT]][, , 1L]
	expected = bm[["intercept"]] + bm[["x1_dyad"]] * (x_T + t(x_T)) / 2 +
		outer(a_T, a_T, "+") + fit$ULUPM[, , d$TT]
	expect_equal(fc, expected, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("bipartite forecasts carry U G V' at its posterior mean", {
	skip_on_cran()
	d = .fc_panel(nA = 10, nB = 8, seed = 9)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 2,
		mode = "bipartite",
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 5L
	))
	rho = 0.7
	frozen = .freeze_dynamics(fit, rho_uv = rho)
	fc = suppressWarnings(predict(frozen, h = 1L, type = "link"))[[1L]]
	bm = if (length(dim(fit$BETA)) == 3L) {
		colMeans(fit$BETA[, , dim(fit$BETA)[3]])
	} else {
		colMeans(fit$BETA)
	}
	expected = bm[["intercept"]] + bm[["x1_dyad"]] * d$X[[d$TT]][, , 1L] +
		outer(fit$a_dynamic[, d$TT], fit$b_dynamic[, d$TT], "+") +
		rho^2 * fit$UVPM[, , d$TT]
	expect_equal(fc, expected, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("the by_draw forecast averages to the point forecast", {
	skip_on_cran()
	obj = .fc_fit_dynamic_uv(); fit = obj$fit
	pt = suppressWarnings(predict(fit, h = 2L, type = "link", seed = 99L))
	bd = suppressWarnings(predict(
		fit,
		h = 2L,
		type = "link",
		by_draw = TRUE,
		seed = 99L
	))
	for (k in seq_len(2L)) {
		expect_equal(
			apply(bd[, , k, , drop = FALSE], c(1L, 2L), mean),
			pt[[k]],
			tolerance = 1e-12,
			ignore_attr = TRUE
		)
	}
})

test_that("R = 0 fits forecast without a multiplicative term", {
	skip_on_cran()
	d = .fc_panel()
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 5L
	))
	frozen = .freeze_dynamics(fit)
	fc = suppressWarnings(predict(frozen, h = 1L, type = "link"))[[1L]]
	bm = if (length(dim(fit$BETA)) == 3L) {
		colMeans(fit$BETA[, , dim(fit$BETA)[3]])
	} else {
		colMeans(fit$BETA)
	}
	expected = bm[["intercept"]] + bm[["x1_dyad"]] * d$X[[d$TT]][, , 1L] +
		outer(fit$a_dynamic[, d$TT], fit$b_dynamic[, d$TT], "+")
	expect_equal(fc, expected, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("bipartite forecasts with static positions keep their own V", {
	skip_on_cran()
	d = .fc_panel(nA = 10, nB = 8, seed = 13)
	fit = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 2,
		mode = "bipartite",
		dynamic_ab = TRUE,
		dynamic_uv = FALSE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 5L
	)))
	frozen = .freeze_dynamics(fit)
	fc = suppressWarnings(predict(frozen, h = 1L, type = "link"))[[1L]]
	expect_equal(dim(fc), c(10L, 8L))
	# static positions: the multiplicative term is the posterior-mean product
	bm = if (length(dim(fit$BETA)) == 3L) {
		colMeans(fit$BETA[, , dim(fit$BETA)[3]])
	} else {
		colMeans(fit$BETA)
	}
	expected = bm[["intercept"]] + bm[["x1_dyad"]] * d$X[[d$TT]][, , 1L] +
		outer(fit$a_dynamic[, d$TT], fit$b_dynamic[, d$TT], "+") + fit$UVPM
	expect_equal(fc, expected, tolerance = 1e-8, ignore_attr = TRUE)
})

# given one Poisson state the mean is exp(eta + s2 / 2), so response forecasts carry the overdispersion factor
test_that("poisson response forecasts carry the overdispersion factor", {
	skip_on_cran()
	set.seed(6886)
	n = 20
	TT = 5
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.4)
	b = rnorm(n, 0, 0.4)
	Y = lapply(seq_len(TT), function(t) {
		z = 0.5 + outer(a, b, "+") + matrix(rnorm(n * n, 0, sqrt(0.5)), n, n)
		m = matrix(rpois(n * n, exp(pmin(z, 5))), n, n)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(lame(
		Y[1:4],
		family = "poisson",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 400,
		burn = 150,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	expect_no_warning({ fc = predict(fit, h = 1L, type = "response") })
	obs = !is.na(Y[[5]])
	ratio = mean(fc[[1]][obs]) / mean(Y[[5]][obs])
	# with s2 near 0.5 the bare exp(eta) is about 0.78 of the mean
	expect_gt(ratio, 0.85)
	expect_lt(ratio, 1.2)
})

# interval = "credible" covers the conditional mean; "prediction" adds the outcome noise
.fc_normal_fit = function(n = 14, TT = 5, seed = 21, nscan = 400) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	X = lapply(seq_len(TT), function(t) {
		arr = array(rnorm(n * n), c(n, n, 1))
		dimnames(arr) = list(nm, nm, "x1"); arr })
	Y = lapply(seq_len(TT), function(t) {
		m = 0.2 + 0.7 * X[[t]][, , 1] + outer(a, b, "+") +
			matrix(rnorm(n * n, 0, 1), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m })
	fit = suppressWarnings(lame(
		Y[seq_len(TT - 1L)],
		Xdyad = X[seq_len(TT - 1L)],
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		nscan = nscan,
		burn = 150,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 4L
	))
	list(fit = fit, Y = Y, X = X, TT = TT, n = n)
}

test_that("a normal prediction interval is the credible interval widened by ve", {
	skip_on_cran()
	o = .fc_normal_fit()
	nd = list(o$X[[o$TT]])
	cr = predict(o$fit, h = 1L, newdata = nd, interval = "credible", seed = 7L)[[1]]
	pr = predict(o$fit, h = 1L, newdata = nd, interval = "prediction", seed = 7L)[[1]]
	w_cr = mean(.offdiag(cr$upper - cr$lower))
	w_pr = mean(.offdiag(pr$upper - pr$lower))
	expect_gt(w_pr, w_cr)
	# the variances add: the prediction width is the credible width and 2 * 1.96 * sqrt(ve) in quadrature
	ve = mean(o$fit$VC[, "ve"])
	expect_equal(
		w_pr,
		sqrt(w_cr^2 + (2 * qnorm(0.975) * sqrt(ve))^2),
		tolerance = 0.12
	)
	# the medians agree, so the point forecast is untouched
	expect_equal(
		mean(.offdiag(pr$median)),
		mean(.offdiag(cr$median)),
		tolerance = 0.1
	)
})

test_that("interval = \"credible\" still quantiles the by_draw mean path", {
	skip_on_cran()
	o = .fc_normal_fit()
	nd = list(o$X[[o$TT]])
	cr = predict(o$fit, h = 1L, newdata = nd, interval = "credible", seed = 7L)[[1]]
	bd = predict(o$fit, h = 1L, newdata = nd, by_draw = TRUE, seed = 7L)
	slab = array(bd[, , 1L, ], c(dim(bd)[1L], dim(bd)[2L], dim(bd)[4L]))
	expect_equal(
		cr$lower,
		apply(slab, c(1L, 2L), quantile, probs = 0.025),
		tolerance = 1e-12,
		ignore_attr = TRUE
	)
	expect_equal(
		cr$upper,
		apply(slab, c(1L, 2L), quantile, probs = 0.975),
		tolerance = 1e-12,
		ignore_attr = TRUE
	)
})

test_that("a 0/1 outcome keeps the interval on the probability and says so", {
	skip_on_cran()
	obj = .fc_fit_dynamic_uv(); fit = obj$fit
	expect_message(
		pr <- suppressWarnings(predict(fit, h = 1L, interval = "prediction", seed = 3L)),
		"interval on the tie probability"
	)
	cr = suppressWarnings(predict(fit, h = 1L, interval = "credible", seed = 3L))
	expect_equal(pr[[1]]$lower, cr[[1]]$lower, tolerance = 1e-12)
	expect_equal(pr[[1]]$upper, cr[[1]]$upper, tolerance = 1e-12)
})

test_that("a poisson prediction interval holds counts and needs the response scale", {
	skip_on_cran()
	set.seed(6886)
	n = 14; TT = 4
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.3); b = rnorm(n, 0, 0.3)
	Y = lapply(seq_len(TT), function(t) {
		z = 0.5 + outer(a, b, "+") + matrix(rnorm(n * n, 0, sqrt(0.4)), n, n)
		m = matrix(rpois(n * n, exp(pmin(z, 5))), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m })
	fit = suppressWarnings(lame(
		Y,
		family = "poisson",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 400,
		burn = 150,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	pr = predict(fit, h = 1L, interval = "prediction", seed = 5L)[[1]]
	cr = predict(fit, h = 1L, interval = "credible", seed = 5L)[[1]]
	# bounds are quantiles of drawn counts: non-negative and wider than the interval of the mean
	expect_true(all(is.finite(.offdiag(pr$lower)) & .offdiag(pr$lower) >= 0))
	expect_true(all(is.na(diag(pr$lower))))
	expect_gt(
		mean(.offdiag(pr$upper - pr$lower)),
		2 * mean(.offdiag(cr$upper - cr$lower))
	)
	# a count has no link-scale draw
	expect_error(
		predict(fit, h = 1L, type = "link", interval = "prediction"),
		"type = \"response\""
	)
})

test_that("an interval request with h = 0 warns for either kind", {
	skip_on_cran()
	obj = .fc_fit_dynamic_uv(); fit = obj$fit
	expect_warning(predict(fit, interval = "credible"), "only has an effect")
	expect_warning(predict(fit, interval = "prediction"), "only has an effect")
})


test_that("dynamic coefficient forecasts preserve covariate rescaling", {
	dat = simulate_test_network(n = 5, n_time = 3, family = "normal", seed = 2726)
	f = suppressWarnings(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		burn = 5,
		nscan = 30,
		odens = 2,
		seed = 91,
		verbose = FALSE
	))
	g = f;nm = dimnames(f$BETA)[[2]][f$beta_dynamic_mask][1]
	g$Xlist = lapply(g$Xlist, function(x) {x[, , nm] = x[, , nm] * 10; x})
	g$BETA[, nm, ] = g$BETA[, nm, ] / 10
	# the long-run mean of an ar1 path is in the coefficient's units as well
	if (!is.null(g$MU_BETA)) g$MU_BETA[, nm] = g$MU_BETA[, nm] / 10
	a = suppressWarnings(predict(f, h = 2, type = "link", n_draws = 15, by_draw = TRUE, seed = 19))
	b = suppressWarnings(predict(g, h = 2, type = "link", n_draws = 15, by_draw = TRUE, seed = 19))
	expect_equal(a, b, tolerance = 1e-10)
	g = f;g$SIGMA_BETA = g$SIGMA_BETA[-1, , drop = FALSE]
	expect_error(suppressWarnings(predict(g, h = 1)), 'aligned')
	g = f;g$RHO_BETA = NULL
	expect_error(suppressWarnings(predict(g, h = 1)), 'aligned')
})

test_that("a dynamic coefficient block without transition draws is refused", {
	dat = simulate_test_network(n = 5, n_time = 3, family = "normal", seed = 2726)
	f = suppressWarnings(lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		burn = 5,
		nscan = 30,
		odens = 2,
		seed = 91,
		verbose = FALSE
	))
	# the transition check runs first and already refuses a misaligned fit; bypass it to reach the forecast
	scale = lame:::.lame_beta_transition_scale(f)
	local_mocked_bindings(.lame_beta_transition_scale = function(object, ...) scale, .package = "lame")
	g = f
	colnames(g$RHO_BETA) = colnames(g$SIGMA_BETA) = paste0(colnames(f$RHO_BETA), "_other")
	expect_error(
		suppressWarnings(predict(g, h = 1, type = "link")),
		class = "lame_forecast_block_missing"
	)
	expect_no_error(suppressWarnings(predict(f, h = 1, type = "link")))
})

# counterfactual and newdata forecasts

build_fc_fit = function(seed = 2026, nscan = 300, burn = 80) {
	set.seed(seed)
	n = 25; Tn = 5
	rho_true = 0.7; beta_bar = 0.5
	beta_t = numeric(Tn); beta_t[1] = 0.9
	for (t in 2:Tn) beta_t[t] = beta_bar + rho_true * (beta_t[t - 1] - beta_bar) +
		rnorm(1, 0, 0.15)
	Xdyad = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n * n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "trade"))
	})
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Y = lapply(seq_len(Tn), function(t) {
		eta = -0.5 + beta_t[t] * Xdyad[[t]][, , 1] + outer(a, b, "+")
		Yt = matrix(rbinom(n * n, 1, pnorm(eta)), n, n); diag(Yt) = NA
		rownames(Yt) = colnames(Yt) = sprintf("a%02d", seq_len(n)); Yt
	})
	names(Y) = paste0("t", seq_len(Tn))
	list(
		fit = lame(
			Y,
			Xdyad = Xdyad,
			family = "binary",
			R = 0,
			dynamic_beta = "dyad",
			nscan = nscan,
			burn = burn,
			odens = 5,
			verbose = FALSE
		),
		Xdyad = Xdyad,
		Tn = Tn
	)
}

test_that("counterfactual raising a positive-coef covariate raises the link predictor", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit; Xdyad = obj$Xdyad; Tn = obj$Tn
	expect_gt(mean(coef(fit)["trade_dyad", ]), 0.2)

	last_X = Xdyad[[Tn]]
	X_future = list(last_X, last_X, last_X)
	X_up = lapply(X_future, function(x) { x[, , 1] = x[, , 1] + 1; x })

	lk_cf = suppressWarnings(predict(fit, h = 3, type = "link", newdata = X_future))
	lk_up = suppressWarnings(predict(fit, h = 3, type = "link", newdata = X_up))

	# the linear predictor increases
	delta_link = mean(lk_up[[3]] - lk_cf[[3]], na.rm = TRUE)
	expect_gt(delta_link, 0)

	# response-scale mean shift is also positive
	rc_cf = suppressWarnings(predict(fit, h = 3, type = "response", newdata = X_future))
	rc_up = suppressWarnings(predict(fit, h = 3, type = "response", newdata = X_up))
	expect_gt(mean(rc_up[[3]] - rc_cf[[3]], na.rm = TRUE), 0)
})

test_that("forecast newdata with the wrong slice count errors clearly", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit
	bad = array(0, c(25, 25, 4)) # 4 slices, model has 1 dyadic coefficient
	expect_error(
		predict(fit, h = 1, type = "response", newdata = list(bad)),
		"slice"
	)
})

test_that("forecasts keep each actor's level under dynamic effects", {
	skip_on_cran()
	set.seed(22); n = 20; TT = 6
	lev = rnorm(n, 0, 1)
	Y = lapply(
		1:TT,
		function(t) { m = outer(lev, rep(0, n), "+") + matrix(rnorm(n * n), n); diag(m) = NA; dimnames(m) = list(
			sprintf("a%02d", 1:n),
			sprintf("a%02d", 1:n)
		); m }
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				family = "normal",
				R = 0,
				dynamic_ab = TRUE,
				burn = 100,
				nscan = 300,
				odens = 3,
				seed = 3,
				verbose = FALSE,
				plot = FALSE
			)
		)
	)
	expect_true(!is.null(f$a_level)); expect_equal(length(f$a_level), n)
	expect_gt(cor(f$a_level, lev), 0.9)
	# far ahead the expected effect is the level, not zero
	pr = predict(f, h = 25, n_draws = 40)
	eta_far = pr[[length(pr)]]
	a_far = rowMeans(eta_far, na.rm = TRUE) - mean(rowMeans(eta_far, na.rm = TRUE))
	expect_gt(cor(a_far, lev), 0.8)
	expect_gt(sd(a_far), 0.5 * sd(lev))
})
