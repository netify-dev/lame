skip_on_cran()

pit_distribution_fit = function(
	beta = c(-2, -2, 2, 2),
	ve = rep(.01, length(beta)),
	family = "normal",
	bip = FALSE,
	symmetric = FALSE
) {
	S = length(beta); nm = c("a", "b", "c")
	f = structure(
		list(
			BETA = matrix(beta, S, 1, dimnames = list(NULL, "intercept")),
			VC = cbind(ve = ve, rho = rep(.2, S), sentinel = rep(1000, S)),
			family = family,
			R = 0,
			mode = if (bip) "bipartite" else "unipartite",
			symmetric = symmetric,
			APM = setNames(numeric(3), nm),
			BPM = setNames(numeric(3), nm),
			Xlist = list(array(1, c(3, 3, 1), dimnames = list(nm, nm, "intercept"))),
			forecast_state = list(a = matrix(0, 3, S), b = matrix(0, 3, S), draw_index = seq_len(S))
		),
		class = "lame"
	)
	f
}
pit_distribution_future = function(value = -1) matrix(value, 3, 3, dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
pit_distribution_link = function(fit, n_draws = NULL, h = 1L, newdesign = NULL) {
	get(
		".forecast_lame",
		envir = asNamespace("lame")
	)(fit, h = h, type = "link", by_draw = TRUE, n_draws = n_draws, newdesign = newdesign)
}
pit_distribution_private = function(name) get(name, envir = environment(forecast_pit), inherits = TRUE)

test_that("normal PIT averages matched Gaussian laws and permits one draw", {
	f = pit_distribution_fit(); Y = pit_distribution_future()
	set.seed(691); eta = suppressWarnings(pit_distribution_link(f)); expected_rng = .Random.seed
	idx = attr(eta, "draw_index")
	expected = mean(pnorm(-1, mean = f$BETA[idx, 1], sd = sqrt(f$VC[idx, "ve"])))
	set.seed(691); result = suppressWarnings(forecast_pit(f, Y))
	expect_equal(result$pit, rep(expected, 6), tolerance = 1e-14)
	expect_identical(.Random.seed, expected_rng)
	expect_identical(result$draw_index, idx)
	expect_gt(abs(expected - pnorm(-1, mean(f$BETA), sqrt(var(f$BETA) + mean(f$VC[, "ve"])))), .1)
	f = pit_distribution_fit(c(0, 0, 0), c(1, 1000, .01))
	set.seed(692); eta = suppressWarnings(pit_distribution_link(f, 1)); idx = attr(eta, "draw_index")
	set.seed(692); result = suppressWarnings(forecast_pit(f, pit_distribution_future(1), n_draws = 1))
	expect_equal(result$pit, rep(pnorm(1, sd = sqrt(f$VC[idx, "ve"])), 6), tolerance = 1e-14)
	expect_equal(result$cdf_lower, result$cdf_upper, tolerance = 0)
})

test_that("normal PIT preserves joint actor cancellation and future covariates", {
	f = pit_distribution_fit(c(-3, -1, 1, 3), rep(1, 4))
	f$forecast_state$a[, ] = f$forecast_state$b[, ] = rep(-f$BETA[, 1] / 2, each = 3)
	set.seed(693); result = suppressWarnings(forecast_pit(f, pit_distribution_future(0)))
	expect_equal(result$pit, rep(.5, 6), tolerance = 0)
	f$forecast_state$a[] = f$forecast_state$b[] = 0
	X = f$Xlist[[1]]; X[] = 2
	set.seed(694); result = suppressWarnings(forecast_pit(f, pit_distribution_future(1), newdesign = list(X)))
	expect_equal(result$pit, rep(mean(pnorm(1, mean = 2 * f$BETA[, 1])), 6), tolerance = 1e-14)
})

test_that("binary PIT keeps exact masses and randomized RNG order", {
	f = pit_distribution_fit(c(-1, 1), c(1, 1), "binary"); Y = pit_distribution_future(0)
	Y[1, 2] = 1
	set.seed(695); eta = suppressWarnings(pit_distribution_link(f)); uniform = runif(6); expected_rng = .Random.seed
	y = Y[row(Y) != col(Y)]; lo = ifelse(y == 0, 0, .5); hi = ifelse(y == 0, .5, 1)
	set.seed(695); result = suppressWarnings(forecast_pit(f, Y))
	expect_equal(result$pit, lo + uniform * (hi - lo), tolerance = 1e-14)
	expect_equal(result$cdf_lower, lo, tolerance = 1e-14)
	expect_equal(result$cdf_upper, hi, tolerance = 1e-14)
	expect_identical(.Random.seed, expected_rng)
	f$BETA[] = -40; result = suppressWarnings(forecast_pit(f, pit_distribution_future(0)))
	expect_equal(result$cdf_upper, rep(1, 6), tolerance = 0)
	expect_error(suppressWarnings(forecast_pit(f, pit_distribution_future(2))), "zero or one")
})

test_that("Poisson PIT integrates rates and lognormal effects with exposure", {
	f = pit_distribution_fit(log(c(.1, 9.9)), c(0, 0), "poisson")
	set.seed(696); result = suppressWarnings(forecast_pit(f, pit_distribution_future(0)))
	expected = mean(exp(-c(.1, 9.9)))
	expect_equal(result$cdf_lower, rep(0, 6), tolerance = 0)
	expect_equal(result$cdf_upper, rep(expected, 6), tolerance = 1e-14)
	expect_gt(abs(expected - exp(-5)), .4)
	f = pit_distribution_fit(0, 4, "poisson")
	cdf = integrate(
		function(z) exp(-2 * exp(2 * z)) * dnorm(z),
		-Inf,
		Inf,
		rel.tol = 1e-11,
		subdivisions = 1000
	)$value
	result = suppressWarnings(forecast_pit(f, pit_distribution_future(0), newexposure = 2))
	expect_equal(result$cdf_upper, rep(cdf, 6), tolerance = 1e-9)
	fn = pit_distribution_private(".lame_pois_cdf_mass")
	for (v in c(1e-12, .5, 4, 25)) {
		oracle = integrate(
			function(z) ppois(2, exp(-.3 + sqrt(v) * z)) * dnorm(z),
			-Inf,
			Inf,
			rel.tol = 1e-11,
			subdivisions = 1000
		)$value
		ans = fn(3, -.3, v)
		expect_equal(unname(ans[1, "lower"]), oracle, tolerance = 1e-9)
	}
	expect_equal(unname(fn(0, 1000, 0)), matrix(c(0, 0), 1), tolerance = 0)
	expect_equal(unname(fn(0, -1000, 0)), matrix(c(0, 1), 1), tolerance = 0)
})

test_that("zero exposure has the degenerate Poisson law and masks stay structural", {
	f = pit_distribution_fit(0, 1, "poisson"); Y = pit_distribution_future(0)
	Y[1, 2] = 1; Y[3, 2] = NA
	result = suppressWarnings(forecast_pit(f, Y, newexposure = 0))
	expect_equal(result$n_obs, 5)
	expect_equal(result$cdf_lower, c(0, 0, 1, 0, 0), tolerance = 0)
	expect_equal(result$cdf_upper, rep(1, 5), tolerance = 0)
	expect_equal(result$pit[3], 1, tolerance = 0)
	expect_error(suppressWarnings(forecast_pit(f, Y, newexposure = -1)), "nonnegative")
	expect_error(suppressWarnings(forecast_pit(f, pit_distribution_future(.5))), "integer counts")
})

test_that("ordinal PIT uses matched posterior cutpoints and original labels", {
	f = pit_distribution_fit(c(.5, .5), c(1, 1), "ordinal")
	f$ordinal_cutpoints = "explicit"; f$ordinal_levels = c(10, 20, 30)
	f$ALPHA = matrix(c(1, 2), 2, 1)
	Y = pit_distribution_future(20); Y[2, 1] = 10; Y[3, 1] = 30
	set.seed(697); eta = suppressWarnings(pit_distribution_link(f, 1)); idx = attr(eta, "draw_index"); u = runif(6)
	labels = match(Y[row(Y) != col(Y)], f$ordinal_levels); cuts = c(-Inf, 0, f$ALPHA[idx, 1], Inf)
	lo = pnorm(cuts[labels] - .5); hi = pnorm(cuts[labels + 1] - .5)
	set.seed(697); result = suppressWarnings(forecast_pit(f, Y, n_draws = 1))
	expect_equal(result$cdf_lower, lo, tolerance = 1e-14)
	expect_equal(result$cdf_upper, hi, tolerance = 1e-14)
	expect_equal(result$pit, lo + u * (hi - lo), tolerance = 1e-14)
	f$ordinal_cutpoints = "data_induced"
	expect_error(forecast_pit(f, Y), class = "lame_forecast_ordinal_model")
	f$ordinal_cutpoints = "explicit"; f$ALPHA[] = -1
	expect_error(suppressWarnings(forecast_pit(f, Y)), "strictly increasing")
})

test_that("PIT aligns named actors and excludes only unipartite self cells", {
	f = pit_distribution_fit(0, 1); Y = pit_distribution_future(0)
	Y[2, 1] = 1; Y[3, 2] = NA
	a = suppressWarnings(forecast_pit(f, Y)); b = suppressWarnings(forecast_pit(f, Y[3:1, 3:1]))
	expect_equal(a$pit, b$pit, tolerance = 0); expect_equal(a$cells, b$cells)
	expect_equal(a$n_obs, 5)
	f$symmetric = TRUE
	expect_equal(suppressWarnings(forecast_pit(f, Y))$n_obs, 3)
	f$symmetric = FALSE; f$mode = "bipartite"
	expect_equal(suppressWarnings(forecast_pit(f, Y))$n_obs, 8)
	Y[, ] = NA
	empty = suppressWarnings(forecast_pit(f, Y)); expect_equal(empty$n_obs, 0)
	expect_true(is.na(empty$ks_stat)); expect_true(is.na(empty$cover_95))
	expect_match(a$ks_inference, "not calibrated")
	print_body = paste(deparse(body(pit_distribution_private("print.forecast_pit"))), collapse = " ")
	expect_match(print_body, "descriptive"); expect_false(grepl("KS rejects", print_body, fixed = TRUE))
})

pit_distribution_robust_fit = function(df = 6) {
	n = 3L; tt = 2L; S = 3L; rank = 1L; actors = c("a", "b", "c")
	hyper = cbind(
		ve = c(.2, 2, 20),
		rho_mu = .6,
		rho_ab = .7,
		rho_uv = .4,
		var_mu = .3,
		var_ab = .2,
		var_uv = .1
	)
	structure(
		list(
			data = list(
				n = n,
				rank = rank,
				tt = tt,
				p = 0L,
				actor_names = actors,
				years = c("t1", "t2"),
				X = matrix(numeric(), n * n * tt, 0)
			),
			R = rank,
			family = "persistent and transient scale mixtures",
			actor_df = 5,
			observation_df = df,
			draws = list(
				hyper = hyper,
				mu = matrix(seq_len(S * tt) / 10, S, tt),
				A = array(0, c(n, tt, S)),
				B = array(0, c(n, tt, S)),
				a0 = matrix(0, n, S),
				b0 = matrix(0, n, S),
				U = array(.2, c(n, rank, tt, S)),
				V = array(.3, c(n, rank, tt, S)),
				beta = matrix(numeric(), S, 0),
				lambda = array(rep(c(.25, 1, 4), each = n * n), c(n, n, S))
			)
		),
		class = "lame_robust"
	)
}

test_that("robust PIT uses genuine future innovations and fresh transient Student law", {
	for (df in c(6, Inf)) {
		f = pit_distribution_robust_fit(df); Y = list(pit_distribution_future(1), pit_distribution_future(-.3))
		Y[[2]][3, 2] = NA
		set.seed(698); eta = predict(
			f,
			h = 2,
			by_draw = TRUE,
			nsim = 2
		); idx = attr(eta, "draw_index"); expected_rng = .Random.seed
		oracle = numeric()
		for (k in 1:2) {
			keep = which(is.finite(Y[[k]]) & row(Y[[k]]) != col(Y[[k]]))
			E = matrix(eta[, , k, , drop = FALSE], 9, 2)[keep, , drop = FALSE]
			L = matrix(f$draws$lambda[, , idx, drop = FALSE], 9, 2)[keep, , drop = FALSE]
			scale = sqrt(matrix(rep(f$draws$hyper[idx, "ve"], each = length(keep)), length(keep)) / L)
			x = (Y[[k]][keep] - E) / scale
			oracle = c(oracle, rowMeans(if (is.finite(df)) pt(x, df) else pnorm(x)))
		}
		set.seed(698); result = forecast_pit(f, Y, n_draws = 2)
		expect_equal(result$pit, oracle, tolerance = 1e-14)
		expect_identical(result$draw_index, idx); expect_identical(.Random.seed, expected_rng)
		expect_equal(result$n_obs, 11)
		set.seed(698); static = get(".lame_robust_prediction", asNamespace("lame"))(f,
			get(".lame_robust_design", asNamespace("lame"))(f, 2, NULL), idx, 2, innovations = FALSE, noise = FALSE)
		expect_gt(max(abs(eta - static), na.rm = TRUE), .01)
		f$draws$lambda[1, 2, ] = 0
		expect_error(forecast_pit(f, Y), "positive pair[[:space:]]+precisions")
	}
})

test_that("LFO and shared Poisson log likelihood preserve valid tiny variances", {
	f = pit_distribution_fit(0, 1e-12); Y = pit_distribution_future(0)
	f$Y = list(Y, Y); f$Xlist = rep(f$Xlist, 2)
	result = suppressWarnings(lfo(f, periods = 2, refit = FALSE))
	expect_equal(result$pointwise[[1]], rep(dnorm(0, sd = 1e-6, log = TRUE), 6), tolerance = 1e-12)
	f$VC[, "ve"] = 0
	expect_error(suppressWarnings(lfo(f, periods = 2, refit = FALSE)), "positive")
	fn = get(".lame_pois_marginal_loglik", envir = environment(lfo), inherits = TRUE)
	value = fn(0, log(2), 1e-12)
	expect_equal(value, -2, tolerance = 2e-10)
	expect_gt(abs(fn(0, log(2), 1e-8) - value), 1e-9)
	expect_error(fn(0, 0, 0), "positive")
	expect_error(fn(0, 0, Inf), "positive")
})


test_that("robust transient CDF matches an independent Gamma-scale integral", {
	x = .8; df = 6
	oracle = integrate(
		function(w) pnorm(x * sqrt(w)) * dgamma(w, df / 2, rate = df / 2),
		0,
		Inf,
		rel.tol = 1e-12
	)$value
	expect_equal(pt(x, df), oracle, tolerance = 1e-11)
	f = pit_distribution_robust_fit(); f$draws$lambda = NULL
	expect_error(forecast_pit(f, pit_distribution_future()), class = "lame_forecast_robust_draws")
})

test_that("post-hoc AME Poisson scoring reuses the named actual-variance law", {
	f = pit_distribution_fit(log(2), 1e-20, "poisson")
	fn = get(".ame_compute_log_lik_post_hoc", envir = environment(forecast_pit), inherits = TRUE)
	value = fn(f, pit_distribution_future(0), "poisson")
	expect_equal(as.numeric(value), rep(-2, 6), tolerance = 1e-12)
	f$VC[, "ve"] = 0
	expect_error(fn(f, pit_distribution_future(0), "poisson"), "positive")
})


test_that("standardized Poisson quadrature retains tiny physical perturbations", {
	fn = get(".lame_pois_marginal_loglik", envir = environment(lfo), inherits = TRUE)
	for (v in c(1e-20, 1e-16, 1e-12, 1e-8)) {
		oracle = log(
			integrate(
				function(z) exp(-2 * exp(sqrt(v) * z)) * dnorm(z),
				-Inf,
				Inf,
				rel.tol = 1e-12
			)$value
		)
		expect_equal(fn(0, log(2), v), oracle, tolerance = 1e-12)
	}
	for (v in c(.01, .5, 4)) {
		oracle = log(
			integrate(
				function(z) dpois(2, exp(.3 + sqrt(v) * z)) * dnorm(z),
				-Inf,
				Inf,
				rel.tol = 1e-12,
				subdivisions = 1000
			)$value
		)
		expect_equal(fn(2, .3, v), oracle, tolerance = 1e-9)
	}
})


test_that("named future subsets score only shared actors and reject ambiguous names", {
	f = pit_distribution_fit(0, 1)
	Y = matrix(0, 4, 4, dimnames = list(c("d", "c", "a", "b"), c("d", "c", "a", "b")))
	result = suppressWarnings(forecast_pit(f, Y))
	expect_equal(result$n_obs, 6); expect_equal(result$pit, rep(.5, 6), tolerance = 0)
	expect_false(any(result$cells$row_actor == "d" | result$cells$col_actor == "d"))
	rownames(Y)[1] = "a"
	expect_error(forecast_pit(f, Y), "unique")
})

test_that("Poisson CDF and mass match independent adaptive integration across scales", {
	fn = pit_distribution_private(".lame_pois_cdf_mass")
	grid = expand.grid(y = c(0, 1, 20), eta = c(-4, .3, 4), ve = c(.0001, .1, 1, 4, 16, 100))
	for (i in seq_len(nrow(grid))) {
		x = grid[i, ]
		ans = fn(x$y, x$eta, x$ve)
		for (mass in c(FALSE, TRUE)) {
			integrand = function(z) {
				lambda = exp(x$eta + sqrt(x$ve) * z)
				(if (mass) dpois(x$y, lambda) else ppois(x$y - 1, lambda)) * dnorm(z)
			}
			oracle = integrate(integrand, -Inf, Inf, rel.tol = 1e-10, abs.tol = 1e-11, subdivisions = 1000)$value
			# the integration bounds absolute probability error, including near-zero masses
			expect_lt(abs(unname(ans[1, if (mass) "mass" else "lower"]) - oracle), 1e-9)
		}
	}
})


test_that("large-count Poisson transitions cannot alias a quadrature rule", {
	fn = pit_distribution_private(".lame_pois_cdf_mass")
	for (count in c(1e6, 1e12)) for (offset in c(-.1, .2)) {
		eta = log(count) + offset
		# oracle from P(Poisson(lambda) < y) = P(Gamma(y) > lambda)
		oracle = integrate(
			function(u) pnorm(log(qgamma(u, count)) - eta),
			0,
			1,
			rel.tol = 1e-11,
			abs.tol = 1e-12,
			subdivisions = 1000
		)$value
		ans = fn(count, eta, 1)
		expect_equal(unname(ans[1, "lower"]), oracle, tolerance = 1e-9)
		if (count == 1e12 && offset == -.1) expect_gt(ans[1, "lower"], .53)
		mass_oracle = integrate(
			function(u) dlnorm(qgamma(u, count + 1), eta, 1),
			0,
			1,
			rel.tol = 1e-9,
			abs.tol = 1e-18,
			subdivisions = 1000
		)$value
		expect_equal(unname(ans[1, "mass"]), mass_oracle, tolerance = 1e-9)
	}
	result = suppressWarnings(forecast_pit(pit_distribution_fit(), pit_distribution_future()))
	expect_true(is.na(result$ks_p))
})


test_that("tiny-variance Poisson CDF resolves the concentrated prior mass", {
	fn = pit_distribution_private(".lame_pois_cdf_mass")
	for (v in c(1e-20, 1e-12, .0001)) {
		ans = fn(1, log(.5), v)
		oracle = integrate(
			function(z) ppois(0, exp(log(.5) + sqrt(v) * z)) * dnorm(z),
			-Inf,
			Inf,
			rel.tol = 1e-11
		)$value
		expect_equal(unname(ans[1, "lower"]), oracle, tolerance = 1e-9)
	}
})
