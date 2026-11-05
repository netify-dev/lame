skip_on_cran()

ordinal_prediction_fixture = function(long = TRUE) {
	nm = c("a", "b", "c")
	X = array(1, c(3, 3, 1), dimnames = list(nm, nm, "x_dyad"))
	Y = matrix(0, 3, 3, dimnames = list(nm, nm))
	f = list(
		BETA = matrix(c(-2, .5, 2), 3, 1, dimnames = list(NULL, "x_dyad")),
		VC = cbind(ve = rep(1, 3), rho = rep(.3, 3)),
		family = "ordinal",
		ordinal_cutpoints = "explicit",
		ordinal_levels = c(10, 20, 50),
		ALPHA = matrix(c(.2, 1, 3), 3, 1),
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		mode = "unipartite",
		symmetric = FALSE,
		n = 3,
		APM = setNames(numeric(3), nm),
		BPM = setNames(numeric(3), nm),
		a_draws = matrix(0, 3, 3),
		b_draws = matrix(0, 3, 3),
		X = X,
		Xlist = list(X),
		YPM = if (long) list(Y) else Y,
		EZ = if (long) list(Y) else Y,
		forecast_state = list(a = matrix(0, 3, 3), b = matrix(0, 3, 3), draw_index = 1:3)
	)
	class(f) = if (long) c("lame", "ame") else "ame"
	f
}

test_that("explicit ordinal means use category probabilities and matched forecast cuts", {
	f = ordinal_prediction_fixture()
	set.seed(47328311)
	link = suppressWarnings(lame:::.forecast_lame(f, h = 3, type = "link", by_draw = TRUE, n_draws = 2))
	rng = .Random.seed
	idx = attr(link, "draw_index")
	set.seed(47328311)
	response = suppressWarnings(lame:::.forecast_lame(f, h = 3, type = "response", by_draw = TRUE, n_draws = 2))
	expect_identical(.Random.seed, rng)
	expect_identical(attr(response, "draw_index"), idx)
	for (s in seq_along(idx)) {
		eta = link[, , , s]
		alpha = f$ALPHA[idx[s], 1]
		expected = 10 * pnorm(-eta) + 20 * (pnorm(alpha - eta) - pnorm(-eta)) + 50 * pnorm(eta - alpha)
		expect_equal(response[, , , s], expected, tolerance = 1e-12)
	}
	outcomes = suppressWarnings(lame:::.forecast_lame(f, h = 4, type = "response", by_draw = TRUE, predictive = TRUE))
	expect_true(all(outcomes[is.finite(outcomes)] %in% f$ordinal_levels))
	expect_true(all(is.na(outcomes[cbind(1:3, 1:3, 1, 1)])))
})

test_that("explicit ordinal newdata returns joint posterior response means", {
	for (long in c(FALSE, TRUE)) {
		f = ordinal_prediction_fixture(long)
		X = f$X; X[] = 2
		eta = 2 * f$BETA[, 1]; alpha = f$ALPHA[, 1]
		expected = mean(10 * pnorm(-eta) + 20 * (pnorm(alpha - eta) - pnorm(-eta)) + 50 * pnorm(eta - alpha))
		y = if (long) predict(f, newdata = list(X))[[1]] else predict(f, newdata = X)
		expect_equal(unname(y), matrix(expected, 3, 3), tolerance = 1e-12)
		if (long) {
			y = predict(f, newdata = list(X[c(1, 3), c(1, 3), , drop = FALSE]))[[1]]
			expect_equal(unname(y), matrix(expected, 2, 2), tolerance = 1e-12)
			expect_identical(rownames(y), c("a", "c"))
		}
	}
})

test_that("explicit ordinal simulations retain labels and reject unmatched cuts", {
	f = ordinal_prediction_fixture()
	eta = matrix(0, 3, 3)
	for (s in 1:3) {
		z = lame:::.predictive_observation(eta, f, s, 4)
		expected = z$Z
		keep = is.finite(expected)
		expected[keep] = f$ordinal_levels[findInterval(expected[keep], c(0, f$ALPHA[s, 1]), left.open = TRUE) + 1L]
		expect_identical(z$Y, expected)
	}
	set.seed(47328312)
	y = lame:::.simY_ord_explicit(eta, .9999, c(0, 1), f$ordinal_levels, TRUE)
	expect_identical(y, t(y))
	expect_true(all(is.na(diag(y))))
	f$ALPHA = f$ALPHA[-1, , drop = FALSE]
	before = .Random.seed
	expect_error(lame:::.ordinal_saved_cutpoints(f, 1), "aligned")
	expect_identical(.Random.seed, before)
	expect_error(lame:::.predictive_observation(eta, f, 1), "aligned")
	expect_identical(.Random.seed, before)
})

test_that("fitter ordinal predictive draws use each stored cutpoint row", {
	original = lame:::.simY_ord_explicit
	records = list()
	testthat::local_mocked_bindings(.simY_ord_explicit = function(EZ, rho, alpha, levels, symmetric = FALSE) {
		y = original(EZ, rho, alpha, levels, symmetric)
		records[[length(records) + 1L]] <<- list(alpha = alpha, y = y)
		y
	}, .package = "lame")
	Y = matrix(rep(c(10, 20, 50), 12), 6, 6); diag(Y) = NA
	for (long in c(FALSE, TRUE)) {
		records = list()
		args = list(
			Y = if (long) list(Y, Y) else Y,
			family = "ordinal",
			ordinal_cutpoints = "explicit",
			dcor = FALSE,
			R = 0,
			nvar = FALSE,
			rvar = FALSE,
			cvar = FALSE,
			intercept = TRUE,
			nscan = 4,
			burn = 1,
			odens = 1,
			gof = FALSE,
			verbose = FALSE,
			seed = 47328313
		)
		fit = do.call(if (long) lame else ame, args)
		nt = if (long) 2L else 1L
		expect_length(records, nrow(fit$BETA) * nt)
		for (t in seq_len(nt)) {
			ys = lapply(seq_len(nrow(fit$BETA)), function(s) {
				r = records[[(s - 1L) * nt + t]]
				expect_equal(r$alpha, c(0, unname(fit$ALPHA[s, ])), tolerance = 0)
				r$y
			})
			expect_equal(
				unname(if (long) fit$YPM[[t]] else fit$YPM),
				unname(Reduce(`+`, ys) / length(ys)),
				tolerance = 0
			)
		}
	}
})


test_that("symmetric ordinal predictions average asymmetric designs with inactive actors", {
	f = ordinal_prediction_fixture()
	f$symmetric = TRUE
	f$a_draws = f$b_draws = NULL
	state = lame:::.joint_predictive_state(f)
	expect_identical(state$index, 1:3)
	X = f$X
	X[, , 1] = matrix(seq(-3, 5), 3, 3)
	average = (X[, , 1] + t(X[, , 1])) / 2
	means = lapply(1:3, function(s) {
		eta = average * f$BETA[s, 1]
		alpha = f$ALPHA[s, 1]
		10 * pnorm(-eta) + 20 * (pnorm(alpha - eta) - pnorm(-eta)) + 50 * pnorm(eta - alpha)
	})
	prediction = predict(f, newdata = list(X))[[1]]
	expect_equal(unname(prediction), unname(Reduce(`+`, means) / 3), tolerance = 1e-12)
	forecast = suppressWarnings(lame:::.forecast_lame(f, h = 1, newdesign = list(X), type = "response", by_draw = TRUE))
	idx = attr(forecast, "draw_index")
	for (s in seq_along(idx)) expect_equal(unname(forecast[, , 1, s]), unname(means[[idx[s]]]), tolerance = 1e-12)
	set.seed(47328314)
	raw = lame:::.predictive_observation(X[, , 1], f, 1)
	set.seed(47328314)
	averaged = lame:::.predictive_observation(average, f, 1)
	expect_identical(raw, averaged)
	expect_identical(raw$Y, t(raw$Y))
})


test_that("intercept warnings distinguish explicit and data-induced ordinal laws", {
	Y = matrix(rep(c(10, 20, 50), 12), 6, 6); diag(Y) = NA
	for (long in c(FALSE, TRUE)) {
		args = list(
			Y = if (long) list(Y, Y) else Y,
			family = "ordinal",
			dcor = FALSE,
			R = 0,
			nvar = FALSE,
			rvar = FALSE,
			cvar = FALSE,
			intercept = TRUE,
			nscan = 3,
			burn = 1,
			odens = 1,
			gof = FALSE,
			verbose = FALSE,
			seed = 47328315
		)
		fun = if (long) lame else ame
		expect_no_warning(do.call(fun, c(args, list(ordinal_cutpoints = "explicit"))))
		expect_warning(
			do.call(fun, c(args, list(ordinal_cutpoints = "data_induced"))),
			"An intercept is not estimable"
		)
	}
})


test_that("symmetric predictions ignore antisymmetric design perturbations across families", {
	f = ordinal_prediction_fixture(); f$symmetric = TRUE
	f$a_draws = f$b_draws = NULL
	x = matrix(seq(-.4, .4, length.out = 9), 3, 3)
	avg = (x + t(x)) / 2
	X = f$X; X[, , 1] = x
	average_design = X; average_design[, , 1] = avg
	for (family in c("normal", "binary", "poisson", "ordinal", "cbin", "frn")) {
		f$family = family
		set.seed(47328316)
		perturbed = lame:::.predictive_observation(x, f, 1)
		set.seed(47328316)
		unperturbed = lame:::.predictive_observation(avg, f, 1)
		expect_identical(perturbed, unperturbed)
		expect_equal(perturbed$Y, t(perturbed$Y), tolerance = 0)
		a = suppressWarnings(
			lame:::.forecast_lame(f, h = 1, newdesign = list(X), type = "link", by_draw = TRUE, seed = 47328317)
		)
		b = suppressWarnings(
			lame:::.forecast_lame(
				f,
				h = 1,
				newdesign = list(average_design),
				type = "link",
				by_draw = TRUE,
				seed = 47328317
			)
		)
		expect_equal(a, b, tolerance = 1e-12)
	}
})


test_that("explicit ordinal simulation retains structural missing predictors", {
	eta = matrix(seq(-.4, .4, length.out = 9), 3, 3)
	missing_diagonal = eta; diag(missing_diagonal) = NA_real_
	set.seed(47328318)
	finite = lame:::.simY_ord_explicit(eta, .3, c(0, 1), c(10, 20, 50))
	seed_after = .Random.seed
	set.seed(47328318)
	missing = lame:::.simY_ord_explicit(missing_diagonal, .3, c(0, 1), c(10, 20, 50))
	expect_identical(missing, finite)
	expect_identical(.Random.seed, seed_after)
	missing_diagonal[1, 2] = NA_real_
	y = lame:::.simY_ord_explicit(missing_diagonal, .3, c(0, 1), c(10, 20, 50))
	expect_true(is.na(y[1, 2]))
	missing_diagonal[1, 2] = Inf
	before = .Random.seed
	expect_error(lame:::.simY_ord_explicit(missing_diagonal, .3, c(0, 1), c(10, 20, 50)), "finite or missing")
	expect_identical(.Random.seed, before)
	for (long in c(FALSE, TRUE)) {
		f = ordinal_prediction_fixture(long)
		X = f$X; diag(X[, , 1]) = NA_real_
		design = if (long) list(X) else X
		z = simulate(f, nsim = 1, type = "approximate", newdata = list(Xdyad = design), seed = 47328319)
		y = if (long) z$Y[[1]][[1]] else z$Y[[1]]
		expect_true(all(is.na(diag(y))))
		expect_true(all(y[is.finite(y)] %in% f$ordinal_levels))
	}
})
