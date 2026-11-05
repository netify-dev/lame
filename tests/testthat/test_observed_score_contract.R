skip_on_cran()

score_contract_private = function(name) get(name, envir = environment(lame::lame), inherits = TRUE)

test_that("observed exact rank and cap score requests refuse the wrong likelihood", {
	for (family in c("cbin", "frn", "ordinal")) {
		expect_error(
			score_contract_private(".lame_validate_score_request")(family, "data_induced", "observed_exact", TRUE),
			class = "lame_observed_score_unsupported"
		)
		expect_invisible(
			score_contract_private(".lame_validate_score_request")(family, "data_induced", "augmented", TRUE)
		)
		expect_invisible(
			score_contract_private(".lame_validate_score_request")(family, "data_induced", "observed_exact", FALSE)
		)
	}
	expect_invisible(
		score_contract_private(".lame_validate_score_request")("ordinal", "explicit", "observed_exact", TRUE)
	)
})

test_that("observed score masks distinguish unordered dyads from bipartite cells", {
	Y = array(1, c(3, 3, 2)); Y[1, 3, 2] = NA
	i = score_contract_private(".lame_score_index")(Y, FALSE, TRUE)
	expect_equal(nrow(i), 5L); expect_true(all(i[, 1] < i[, 2]))
	expect_equal(nrow(score_contract_private(".lame_score_index")(Y, FALSE, FALSE)), 11L)
	expect_equal(nrow(score_contract_private(".lame_score_index")(Y, TRUE, FALSE)), 17L)
})

test_that("explicit ordinal scores use actual cutpoints and original category labels", {
	y = c(-8, 2, 19, 2); eta = c(1, -.5, 40, -40); alpha = c(0, 1.7)
	actual = score_contract_private(".lame_explicit_ordinal_loglik")(y, eta, alpha, c(-8, 2, 19))
	expect_equal(actual[1], pnorm(-1, log.p = TRUE), tolerance = 1e-14)
	expect_equal(actual[2], log(pnorm(2.2) - pnorm(.5)), tolerance = 1e-14)
	expect_equal(actual[3], pnorm(1.7 - 40, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-14)
	# tail interval reference integrated around its dominant edge
	ref = dnorm(40, log = TRUE) + log(integrate(function(t) exp(-40 * t - t * t / 2), 0, 1.7, rel.tol = 1e-12)$value)
	expect_equal(actual[4], ref, tolerance = 1e-12)
	expect_equal(
		score_contract_private(".lame_explicit_ordinal_loglik")(c(7, 9), c(-.2, .8), 0, c(7, 9)),
		c(pnorm(.2, log.p = TRUE), pnorm(.8, log.p = TRUE)),
		tolerance = 1e-14
	)
	expect_error(score_contract_private(".lame_explicit_ordinal_loglik")(3, 0, alpha, c(-8, 2, 19)), "categories")
})

test_that("the complete saved predictor includes all contemporaneous components without RNG", {
	set.seed(611); n = 3; m = 2; nt = 2; p = 2
	X = lapply(seq_len(nt), function(t) array(rnorm(n * m * p), c(n, m, p)))
	beta = matrix(rnorm(nt * p), nt, p); a = matrix(rnorm(n * nt), n, nt); b = matrix(rnorm(m * nt), m, nt)
	U = array(rnorm(n * 2 * nt), c(n, 2, nt)); V = array(rnorm(m * nt), c(m, 1, nt))
	G = array(c(.2, .5, 1.2, -.3), c(2, 1, nt)); offset = array(rnorm(n * m * nt), c(n, m, nt))
	before = .Random.seed
	actual = score_contract_private(".lame_saved_predictor")(X, beta, a, b, U, V, G, offset = offset)
	expect_identical(.Random.seed, before)
	for (t in seq_len(nt)) for (i in seq_len(n)) for (j in seq_len(m)) {
		ref = sum(X[[t]][i, j, ] * beta[t, ]) + a[i, t] + b[j, t] + sum(U[i, , t] * G[, 1, t]) * V[j, 1, t] + offset[i, j, t]
		expect_equal(actual[i, j, t], ref, tolerance = 1e-14)
	}
	Xs = lapply(seq_len(nt), function(t) array(rnorm(n * n * p), c(n, n, p)))
	bs = matrix(rnorm(n * nt), n, nt); Vs = array(rnorm(n * 2 * nt), c(n, 2, nt))
	for (active in c(FALSE, TRUE)) {
		raw = score_contract_private(".lame_saved_predictor")(Xs, beta, a, bs, U, Vs, rank_active = active)
		sym = score_contract_private(".lame_saved_predictor")(Xs, beta, a, bs, U, Vs, rank_active = active, symmetric = TRUE)
		for (t in seq_len(nt)) expect_equal(sym[, , t], (raw[, , t] + t(raw[, , t])) / 2, tolerance = 0)
	}
})

test_that("AME posthoc scores match each active draw and its bipartite interaction", {
	Y = matrix(c(1, 2, NA, 4), 2, 2)
	fit = list(
		mode = "bipartite",
		symmetric = FALSE,
		R = 1,
		rvar = TRUE,
		cvar = FALSE,
		BETA = matrix(c(-.3, .8), 2, 1, dimnames = list(NULL, "intercept")),
		VC = matrix(c(.5, 1.2), 2, 1, dimnames = list(NULL, "ve")),
		a_samples = matrix(c(.4, -.1, .7, .2), 2, 2),
		ab_draw_index = c(2L, 1L),
		U_samples = array(c(1, 2, 3, 4), c(2, 1, 2)),
		V_samples = array(c(.2, .5, .7, -.3), c(2, 1, 2)),
		G_samples = array(c(3, -2), c(1, 1, 2)),
		uv_draw_index = c(2L, 1L)
	)
	ll = score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "normal")
	expect_equal(ncol(ll), 3L)
	for (s in 1:2) {
		k = match(s, c(2L, 1L));eta = matrix(fit$BETA[s, 1], 2, 2) + fit$a_samples[, k]
		eta = eta + outer(fit$U_samples[, 1, k], fit$V_samples[, 1, k]) * fit$G_samples[1, 1, k]
		expect_equal(
			ll[s, ],
			dnorm(Y[is.finite(Y)], eta[is.finite(Y)], sqrt(fit$VC[s, 1]), log = TRUE),
			tolerance = 1e-14
		)
	}
	fit$uv_draw_index = c(4L, 2L)
	expect_error(score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "normal"), "matched")
	fit$uv_draw_index = c(2L, 1L);fit$G_samples = NULL
	expect_error(score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "normal"), "interaction")
})

test_that("AME explicit scores match ALPHA rows and symmetric unique observations", {
	Y = matrix(c(NA, 3, 9, 3, NA, 9, 9, 9, NA), 3, 3)
	fit = list(
		mode = "unipartite",
		symmetric = TRUE,
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		BETA = matrix(c(-.2, .8), 2, 1, dimnames = list(NULL, "intercept")),
		VC = matrix(1, 2, 1, dimnames = list(NULL, "ve")),
		ordinal_cutpoints = "explicit",
		ordinal_levels = c(0, 3, 9),
		ALPHA = matrix(c(.5, 2), 2, 1)
	)
	ll = score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "ordinal")
	expect_equal(dim(ll), c(2L, 3L))
	for (s in 1:2) {
		y = Y[upper.tri(Y)];code = match(y, fit$ordinal_levels);cuts = c(-Inf, 0, fit$ALPHA[s, 1], Inf)
		ref = log(pnorm(cuts[code + 1] - fit$BETA[s, 1]) - pnorm(cuts[code] - fit$BETA[s, 1]))
		expect_equal(ll[s, ], ref, tolerance = 1e-14)
	}
	fit$ALPHA = fit$ALPHA[1, , drop = FALSE]
	expect_error(score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "ordinal"), "ALPHA rows")
})

test_that("LOO guards latent, approximate, correlated and legacy scores explicitly", {
	idx = matrix(c(1L, 2L, 1L), 1, 3)
	contract = function(
		family = "binary",
		method = "observed_exact",
		cut = "explicit",
		bip = FALSE,
		sym = FALSE,
		dcor = FALSE,
		rho = NULL
	)
		score_contract_private(".lame_score_contract")(family, method, cut, bip, sym, dcor, idx, rho)
	fit = structure(list(log_lik = matrix(-1, 10, 1)), class = c("lame", "ame"))
	expect_error(loo(fit), class = "lame_loo_contract_missing")
	expect_error(waic(fit), class = "lame_loo_contract_missing")
	bad = list(
		contract(method = "augmented"),
		contract(family = "cbin", method = "observed_ghk"),
		contract(family = "ordinal", cut = "data_induced", method = "observed_ghk"),
		contract(dcor = TRUE),
		contract(rho = .3)
	)
	for (c in bad) {
		fit$log_lik_contract = c
		expect_error(loo(fit), class = "lame_loo_factorization_unsupported")
		expect_error(waic(fit), class = "lame_loo_factorization_unsupported")
		expect_identical(read_log_lik(fit), fit$log_lik)
	}
	for (c in list(contract(), contract(bip = TRUE, dcor = TRUE), contract(sym = TRUE, dcor = TRUE))) {
		expect_true(c$factorization$valid_for_loo)
		fit$log_lik_contract = c
		expect_invisible(score_contract_private(".lame_require_loo_contract")(fit))
	}
})

test_that("saved LAME likelihoods match the current parameter draw and chunked storage", {
	set.seed(771);n = 6;nt = 2
	X = lapply(1:nt, function(t) array(rnorm(n * n), c(n, n, 1)))
	Y = lapply(1:nt, function(t) { y = matrix(rbinom(n * n, 1, .4), n, n);diag(y) = NA;y })
	args = list(
		Y = Y,
		Xdyad = X,
		family = "binary",
		R = 0,
		dcor = FALSE,
		nscan = 8,
		burn = 4,
		odens = 1,
		seed = 991,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		posterior_opts = list(save_ab = TRUE)
	)
	fit = do.call(lame, c(args, list(save_log_lik = TRUE)));rng = .Random.seed
	for (s in seq_len(nrow(fit$BETA))) for (t in 1:nt) {
		eta = matrix(fit$BETA[s, 1], n, n) + fit$BETA[s, 2] * X[[t]][, , 1] + outer(
			fit$a_samples[, s],
			fit$b_samples[, s],
			"+"
		)
		index = fit$log_lik_index;cols = which(index[, 3] == t);pairs = index[cols, 1:2, drop = FALSE]
		ref = pnorm(ifelse(Y[[t]][pairs] > 0, eta[pairs], -eta[pairs]), log.p = TRUE)
		expect_equal(fit$log_lik[s, cols], ref, tolerance = 1e-12)
	}
	expect_true(fit$log_lik_contract$factorization$valid_for_loo)
	plain = do.call(lame, c(args, list(save_log_lik = FALSE)))
	expect_identical(plain$BETA, fit$BETA);expect_identical(plain$VC, fit$VC);expect_identical(.Random.seed, rng)
	chunk = do.call(lame, c(args, list(save_log_lik = "chunked", log_lik_chunk_size = 10L)))
	expect_identical(read_log_lik(chunk), fit$log_lik);expect_identical(.Random.seed, rng)
	if (requireNamespace("loo", quietly = TRUE)) {
		expect_s3_class(suppressWarnings(loo(chunk)), "loo")
		expect_s3_class(suppressWarnings(waic(chunk)), "waic")
	}
})


test_that("beta-only AME scores do not invent a rank or actor block", {
	fit = list(
		BETA = matrix(c(-.3, .7), 2, 1, dimnames = list(NULL, "intercept")),
		VC = matrix(1, 2, 1, dimnames = list(NULL, "ve"))
	)
	Y = matrix(c(0, 1, 1, 0, 0, 1), 2, 3)
	ll = score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "binary")
	for (s in 1:2) expect_equal(
		ll[s, ],
		pnorm(ifelse(Y > 0, fit$BETA[s, 1], -fit$BETA[s, 1]), log.p = TRUE)[seq_along(Y)],
		tolerance = 0
	)
	fit$R = 1L;fit$R_row = 0L;fit$R_col = 1L;fit$mode = "bipartite"
	expect_equal(score_contract_private(".ame_compute_log_lik_post_hoc")(fit, Y, "binary"), ll, tolerance = 0)
})

test_that("short ordinal intervals retain their original width in scores and Cowles target", {
	fn = score_contract_private(".lame_explicit_ordinal_loglik")
	marginal = score_contract_private(".ord_marginal_loglik")
	for (w in c(1e-4, 1e-12, 1e-20, 1e-100)) for (eta in c(0, -40, 40)) {
		# independent change of variables z = w*t on the unshifted cut interval
		ref = dnorm(-eta, log = TRUE) + log(w) + log(
			integrate(
				function(t) exp(eta * w * t - (w * t)^2 / 2),
				0,
				1,
				rel.tol = 1e-12,
				abs.tol = 0
			)$value
		)
		actual = fn(2, eta, c(0, w), c(1, 2, 3))
		expect_equal(actual, ref, tolerance = 1e-12)
		expect_equal(marginal(c(-Inf, 0, w, Inf), eta, 2), ref, tolerance = 1e-12)
	}
	expect_equal(fn(2, 1e6, c(0, 1e-20), c(1, 2, 3)), dnorm(-1e6, log = TRUE) + log(1e-20), tolerance = 1e-14)
	expect_identical(score_contract_private("log_phi_diff")(Inf, -Inf), 0)
	# the Cowles acceptance log ratio is finite for these valid proposals
	before = marginal(c(-Inf, 0, 1e-20, Inf), c(0, -40, 40), rep(2, 3))
	after = marginal(c(-Inf, 0, 2e-20, Inf), c(0, -40, 40), rep(2, 3))
	expect_equal(after - before, 3 * log(2), tolerance = 1e-12)
})


test_that("public rank scoring requests fail with actionable typed contracts", {
	Y = matrix(rep(c(0, 1), 18), 6, 6);diag(Y) = NA
	for (family in c("cbin", "frn", "ordinal")) {
		expect_error(
			suppressWarnings(ame(
				Y,
				family = family,
				R = 0,
				nscan = 4,
				burn = 2,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				save_log_lik = TRUE
			)),
			class = "lame_observed_score_unsupported"
		)
		expect_error(
			suppressWarnings(lame(
				list(Y, Y),
				family = family,
				R = 0,
				nscan = 4,
				burn = 2,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				save_log_lik = TRUE
			)),
			class = "lame_observed_score_unsupported"
		)
	}
})


test_that("explicit approximate AME diagnostics never become observed LOO", {
	Y = matrix(rep(c(0, 1), 18), 6, 6);diag(Y) = NA
	fit = suppressWarnings(ame(
		Y,
		family = "cbin",
		R = 0,
		nscan = 4,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		save_log_lik = TRUE,
		log_lik_method = "observed_approx"
	))
	expect_identical(fit$log_lik_method, "observed_approx")
	expect_identical(fit$log_lik_contract$response_scale, "approximate_proxy")
	expect_false(fit$log_lik_contract$factorization$valid_for_loo)
	expect_error(loo(fit), class = "lame_loo_factorization_unsupported")
	expect_true(is.matrix(read_log_lik(fit)))
})
