# rolling-origin leave-future-out: the one-step forecast uses the held-out period's own covariates
# and each dyad is scored by its posterior predictive density, log mean_s p(y | eta_s)

.lfo_panel = function(n = 20, TT = 5, seed = 42, family = "binary") {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	X = lapply(seq_len(TT), function(t) {
		arr = array(rnorm(n * n), c(n, n, 1))
		dimnames(arr) = list(nm, nm, "x1"); arr })
	Y = lapply(seq_len(TT), function(t) {
		eta = -0.4 + 1.6 * X[[t]][, , 1L] + outer(a, b, "+")
		m = if (identical(family, "binary")) {
			matrix(rbinom(n * n, 1, pnorm(eta)), n, n)
		} else if (identical(family, "poisson")) {
			matrix(rpois(n * n, exp(pmin(eta, 3))), n, n)
		} else {
			eta + matrix(rnorm(n * n, 0, 0.5), n, n)
		}
		diag(m) = NA; dimnames(m) = list(nm, nm); m })
	list(Y = Y, X = X, TT = TT, n = n)
}

# the stored design has an intercept slice, which lfo() strips before refitting and forecasting; the tests do the same
.lfo_stripped_X = function(fit) lapply(fit$Xlist, function(A) A[, , -1L, drop = FALSE])

# posterior predictive elpd from a [n, n, h, draw] link-scale forecast
.lfo_pp_elpd = function(eta_draws, Y_test, family = "binary", s2 = 1) {
	dm = dim(eta_draws)
	dim(eta_draws) = c(dm[1L] * dm[2L], dm[4L])
	obs = which(!is.na(Y_test))
	y = Y_test[obs]
	E = eta_draws[obs, , drop = FALSE]
	ym = matrix(y, nrow = length(y), ncol = ncol(E))
	lp = switch(
		family,
		binary = {
			p = pmin(pmax(pnorm(E), 1e-12), 1 - 1e-12)
			ym * log(p) + (1 - ym) * log(1 - p)
		},
		poisson = vapply(
			seq_len(ncol(E)),
			function(j)
				lame:::.lame_pois_marginal_loglik(y, E[, j], s2[j]),
			numeric(length(y))
		),
		dnorm(ym, mean = E, sd = matrix(rep(sqrt(s2), each = nrow(E)), nrow(E)), log = TRUE)
	)
	mx = apply(lp, 1L, max)
	sum(mx + log(rowMeans(exp(lp - mx))))
}

# plug-in elpd: the family density at the posterior-mean linear predictor
.lfo_plugin_elpd = function(eta_draws, Y_test) {
	dm = dim(eta_draws)
	dim(eta_draws) = c(dm[1L] * dm[2L], dm[4L])
	obs = which(!is.na(Y_test))
	y = Y_test[obs]
	p = pmin(pmax(pnorm(rowMeans(eta_draws)[obs]), 1e-12), 1 - 1e-12)
	sum(y * log(p) + (1 - y) * log(1 - p))
}

# run lfo() keeping the refit it builds, so the same object can be scored independently
.lfo_capture = function(fit, ...) {
	cap = new.env()
	real_lame = lame
	testthat::local_mocked_bindings(lame = function(...) {
		cap$fit = do.call(real_lame, list(...))
		cap$fit
	}, .package = "lame", .env = parent.frame())
	cap$res = suppressWarnings(lfo(fit, ...))
	cap
}

test_that("lfo scores the leave-out period's own covariates", {
	skip_on_cran()
	d = .lfo_panel()
	# a static fit makes the one-step forecast a deterministic function of the stored draws, so the recomputation is exact
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 0,
		nscan = 100,
		burn = 25,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	cap = .lfo_capture(
		fit,
		periods = 5L,
		refit = TRUE,
		nscan = 60,
		burn = 15,
		odens = 5
	)
	Xs = .lfo_stripped_X(fit)
	Y_test = d$Y[[5L]]

	eta_right = suppressWarnings(predict(
		cap$fit,
		h = 1L,
		type = "link",
		newdata = list(Xs[[5L]]),
		by_draw = TRUE
	))
	expect_equal(
		cap$res$elpd_lfo,
		.lfo_pp_elpd(eta_right, Y_test),
		tolerance = 1e-8
	)

	# without newdata the forecaster uses the last training period's covariates, a worse score when the dyadic covariate moves
	eta_stale = suppressWarnings(predict(
		cap$fit,
		h = 1L,
		type = "link",
		by_draw = TRUE
	))
	expect_gt(cap$res$elpd_lfo, .lfo_pp_elpd(eta_stale, Y_test) + 50)
})

test_that("lfo averages the family density over draws rather than scoring the mean", {
	skip_on_cran()
	d = .lfo_panel(seed = 13)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 2,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		nscan = 120,
		burn = 30,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	Xs = .lfo_stripped_X(fit)
	Y_test = d$Y[[d$TT]]
	# widen the AR(1) innovations so the one-step draws spread and the two scoring rules differ
	fit$sigma_ab = rep(0.5, length(fit$sigma_ab))
	fit$sigma_uv = rep(0.5, length(fit$sigma_uv))
	# refit = FALSE leaves the RNG stream untouched before the forecast, so equal seeds give the same draws
	set.seed(404)
	res = suppressWarnings(lfo(fit, periods = d$TT, refit = FALSE))
	set.seed(404)
	eta = suppressWarnings(predict(
		fit,
		h = 1L,
		type = "link",
		newdata = list(Xs[[d$TT]]),
		by_draw = TRUE
	))

	expect_equal(res$elpd_lfo, .lfo_pp_elpd(eta, Y_test), tolerance = 1e-8)
	# the plug-in score at the posterior-mean predictor is a different number
	expect_gt(abs(res$elpd_lfo - .lfo_plugin_elpd(eta, Y_test)), 2)
})

test_that("lfo keeps its pointwise and per-period structure", {
	skip_on_cran()
	d = .lfo_panel(n = 14, TT = 4, seed = 21)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 0,
		nscan = 60,
		burn = 15,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	res = suppressWarnings(lfo(
		fit,
		periods = c(3L, 4L),
		refit = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5
	))
	expect_s3_class(res, "lfo_lame")
	expect_named(res$pointwise, c("3", "4"))
	expect_true(all(vapply(res$pointwise, is.numeric, logical(1L))))
	expect_equal(
		vapply(res$pointwise, length, integer(1L)),
		c("3" = res$per_period$n_obs[1L], "4" = res$per_period$n_obs[2L])
	)
	expect_equal(
		vapply(res$pointwise, sum, numeric(1L)),
		setNames(res$per_period$elpd, c("3", "4"))
	)
	expect_equal(res$elpd_lfo, sum(unlist(res$pointwise)))
	expect_equal(res$per_period$n_obs, rep(14L * 13L, 2L))
	expect_true(all(is.finite(unlist(res$pointwise))))
})

test_that("lfo scores the normal family on the posterior predictive scale", {
	skip_on_cran()
	d = .lfo_panel(n = 14, TT = 4, seed = 33, family = "normal")
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "normal",
		R = 0,
		nscan = 60,
		burn = 15,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	cap = .lfo_capture(
		fit,
		periods = 4L,
		refit = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5
	)
	Xs = .lfo_stripped_X(fit)
	eta = suppressWarnings(predict(
		cap$fit,
		h = 1L,
		type = "link",
		newdata = list(Xs[[4L]]),
		by_draw = TRUE
	))
	s2 = cap$fit$VC[attr(eta, "draw_index"), "ve"]
	expect_equal(
		cap$res$elpd_lfo,
		.lfo_pp_elpd(eta, d$Y[[4L]], family = "normal", s2 = s2),
		tolerance = 1e-8
	)
})

test_that("lfo scores the poisson family on the posterior predictive scale", {
	skip_on_cran()
	d = .lfo_panel(n = 12, TT = 4, seed = 51, family = "poisson")
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 0,
		nscan = 60,
		burn = 15,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	cap = .lfo_capture(
		fit,
		periods = 4L,
		refit = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5
	)
	Xs = .lfo_stripped_X(fit)
	eta = suppressWarnings(predict(
		cap$fit,
		h = 1L,
		type = "link",
		newdata = list(Xs[[4L]]),
		by_draw = TRUE
	))
	expect_equal(
		cap$res$elpd_lfo,
		.lfo_pp_elpd(
			eta,
			d$Y[[4L]],
			family = "poisson",
			s2 = cap$fit$VC[attr(eta, "draw_index"), "ve"]
		),
		tolerance = 1e-8
	)
	expect_true(is.finite(cap$res$elpd_lfo))
})

test_that("lfo still rejects a leave-out period before the second", {
	skip_on_cran()
	d = .lfo_panel(n = 10, TT = 3, seed = 61)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "binary",
		R = 0,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		seed = 7L
	))
	expect_error(
		suppressWarnings(lfo(fit, periods = 1L, refit = FALSE)),
		"leave-out period"
	)
})

test_that("the posterior predictive average is taken on the log scale", {
	# row-wise log mean exp over draws, so an underflowing draw cannot take the average with it
	lp = matrix(c(-1, -2, -1000, -1001), nrow = 2, byrow = TRUE)
	expect_equal(
		lame:::.lfo_log_mean_exp_rows(lp),
		log(rowMeans(exp(lp - c(-1, -1000)))) + c(-1, -1000)
	)
	expect_equal(
		lame:::.lfo_log_mean_exp_rows(lp)[1L],
		log(mean(exp(c(-1, -2))))
	)
	expect_true(all(is.finite(lame:::.lfo_log_mean_exp_rows(lp))))
})

test_that("a held-out period is aligned to the refit's forecast grid by name", {
	nm = sprintf("a%02d", 1:4)
	Y_t = matrix(1, 4, 4, dimnames = list(nm, nm))
	# same grid: returned untouched
	expect_identical(lame:::.lfo_align_heldout(Y_t, nm, nm, 4L, 4L, 2L), Y_t)
	# refit grid covers three of the four actors: the rest is dropped
	al = lame:::.lfo_align_heldout(Y_t, nm[1:3], nm[1:3], 3L, 3L, 2L)
	expect_equal(dim(al), c(3L, 3L))
	expect_equal(dimnames(al), list(nm[1:3], nm[1:3]))
	expect_true(all(al == 1))
	# actors the held-out period does not carry stay NA
	al2 = lame:::.lfo_align_heldout(Y_t[1:3, 1:3, drop = FALSE], nm, nm, 4L, 4L, 2L)
	expect_true(all(is.na(al2[4L, ])) && all(is.na(al2[, 4L])))
	expect_true(all(al2[1:3, 1:3] == 1))
	# without names there is nothing to align on
	Y_bare = matrix(1, 4, 4)
	expect_error(
		lame:::.lfo_align_heldout(Y_bare, nm[1:3], nm[1:3], 3L, 3L, 2L),
		"Actor names"
	)
})

test_that("lfo refits the same model shape as the fit it scores", {
	skip_on_cran()
	set.seed(6886)
	TT = 3
	# bipartite: the refit keeps the mode and both ranks
	nA = 10
	nB = 7
	rn = paste0("r", seq_len(nA))
	cn = paste0("c", seq_len(nB))
	Yb = lapply(seq_len(TT), function(t) {
		m = matrix(rnorm(nA * nB), nA, nB)
		dimnames(m) = list(rn, cn)
		m
	})
	fit_b = suppressWarnings(lame(
		Yb,
		family = "normal",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	cap_b = .lfo_capture(fit_b, periods = 3L, nscan = 20, burn = 10, odens = 5)
	expect_identical(cap_b$fit$mode, "bipartite")
	expect_true(is.finite(cap_b$res$elpd_lfo))
	# symmetric: the refit stays symmetric
	n = 10
	nm = paste0("a", seq_len(n))
	Ys = lapply(seq_len(TT), function(t) {
		m = matrix(rnorm(n * n), n, n)
		m = (m + t(m)) / 2
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit_s = suppressWarnings(lame(
		Ys,
		family = "normal",
		R = 0,
		symmetric = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	))
	cap_s = .lfo_capture(fit_s, periods = 3L, nscan = 20, burn = 10, odens = 5)
	expect_true(isTRUE(cap_s$fit$symmetric))
	expect_true(is.finite(cap_s$res$elpd_lfo))
})


test_that("held-out actor permutations do not silently change the scoring targets", {
	nm = letters[1:4]; Y = matrix(seq_len(16), 4, dimnames = list(nm, nm))
	perm = c(3, 1, 4, 2)
	expect_equal(lame:::.lfo_align_heldout(Y[perm, rev(perm)], nm, nm, 4L, 4L, 3L), Y)
})

test_that("LFO preserves nodal roles, selected coefficients, transitions and user priors", {
	skip_on_cran()
	d = .lfo_panel(n = 10, TT = 3, family = "normal")
	r = lapply(
		seq_len(3),
		function(t) matrix(
			seq_len(10) / 10 + t,
			10,
			1,
			dimnames = list(rownames(d$Y[[t]]), "size")
		)
	)
	prior = list(rho_ab_mean = .6, rho_ab_sd = .2)
	f = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		Xrow = r,
		family = "normal",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		dynamic_uv_kind = "t",
		dynamic_beta = "x1_dyad",
		dynamic_beta_pool = "both",
		dcor = FALSE,
		prior = prior,
		nscan = 30,
		burn = 10,
		odens = 2,
		gof = FALSE
	))
	cap = .lfo_capture(f, periods = 3L, nscan = 20, burn = 5, odens = 2)
	expect_identical(cap$fit$refit_options$dynamic_beta, "x1_dyad")
	expect_identical(cap$fit$refit_options$dynamic_beta_pool, "both")
	expect_identical(cap$fit$dynamic_uv_kind, "t")
	expect_false(cap$fit$dcor)
	expect_identical(cap$fit$refit_options$prior, prior)
	expect_identical(cap$fit$refit_data$Xrow, r[1:2])
	expect_identical(dimnames(cap$fit$BETA)[[2]], dimnames(f$BETA)[[2]])
	expect_true(is.finite(cap$res$elpd_lfo))
	expect_true(is.na(cap$res$p_lfo))
	old = f; old$refit_options = NULL
	expect_error(lfo(old, periods = 3), class = "lame_lfo_refit")
})

test_that("Poisson LFO integrates lognormal noise and held-out exposure", {
	# an independent numerical integral, with deliberately different variances
	eta = c(-.5, .7); ve = c(.2, 1.1); observed = 2; exposure = 3
	y = matrix(observed, 1, 1, dimnames = list("r", "c"))
	f = structure(
		list(
			family = "poisson",
			mode = "bipartite",
			R = 0,
			Y_obs = list(y, y),
			Xlist = list(array(1, c(1, 1, 1)), array(1, c(1, 1, 1))),
			BETA = matrix(eta, 2, 1, dimnames = list(NULL, "intercept")),
			VC = matrix(ve, 2, 1, dimnames = list(NULL, "ve")),
			APM = c(r = 0),
			BPM = c(c = 0),
			period_exposure = c(1, exposure),
			forecast_state = list(a = matrix(0, 1, 2), b = matrix(0, 1, 2), draw_index = 1:2)
		),
		class = c("lame", "ame")
	)
	truth = log(mean(vapply(
		1:2,
		function(j) integrate(
			function(z)
				dpois(observed, exp(z) * exposure) * dnorm(z, eta[j], sqrt(ve[j])),
			-12,
			12,
			rel.tol = 1e-11
		)$value,
		numeric(1)
	)))
	res = suppressWarnings(lfo(f, periods = 2, refit = FALSE))
	expect_equal(res$elpd_lfo, truth, tolerance = 2e-6)
	f$family = "cbin"
	expect_error(lfo(f, periods = 2, refit = FALSE), class = "lame_lfo_family")
})

test_that("symmetric LFO scores each unordered dyad once", {
	skip_on_cran()
	d = .lfo_panel(n = 10, TT = 3, family = "normal")
	Y = lapply(d$Y, function(y) (y + t(y)) / 2)
	f = suppressWarnings(lame(
		Y,
		family = "normal",
		symmetric = TRUE,
		nscan = 20,
		burn = 5,
		odens = 2,
		gof = FALSE
	))
	res = suppressWarnings(lfo(f, periods = 3, nscan = 20, burn = 5, odens = 2))
	expect_equal(res$per_period$n_obs, 45)
})
