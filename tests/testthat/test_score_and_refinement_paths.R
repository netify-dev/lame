# the empirical-cutpoint ordinal pointwise log-likelihood (log_lik_method = "observed_ghk" for ordinal fits)
# and the joint refinement of a directed Poisson dynamic_uv penalised ALS fit that stops short of its criterion

test_that("the ordinal pointwise log-likelihood is the cumulative-probit cell probability at the empirical cutpoints", {
	y = c(1, 1, 2, 3, 3, 3); eta = c(-.5, .2, 0, 1, -1, 2.5)
	cuts = c(-Inf, qnorm(2 / 6), qnorm(3 / 6), Inf)
	k = y - min(y) + 1
	expected = log(pnorm(cuts[k + 1] - eta) - pnorm(cuts[k] - eta))
	expect_equal(lame:::.ordinal_pointwise_loglik(y, eta), expected, tolerance = 1e-12)
	# the labels only order the categories
	expect_equal(lame:::.ordinal_pointwise_loglik(y + 10, eta), expected, tolerance = 1e-12)
	# at the empirical cutpoints and eta = 0 a cell's probability is its category's share
	expect_equal(exp(lame:::.ordinal_pointwise_loglik(y, rep(0, 6))), c(2, 2, 1, 3, 3, 3) / 6, tolerance = 1e-12)
	# far from its category the probability is floored, so the score stays finite
	expect_equal(lame:::.ordinal_pointwise_loglik(c(1, 2), c(40, 0))[1], log(1e-12))
})

test_that("an ordinal fit stores finite observed_ghk scores, one per stored draw and observed cell", {
	set.seed(5); n = 8; nm = sprintf("a%02d", 1:n)
	Y = lapply(
		1:3,
		function(t) { m = matrix(sample(1:3, n * n, TRUE), n, n, dimnames = list(nm, nm)); diag(m) = NA; m }
	)
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "ordinal",
		R = 0,
		burn = 5,
		nscan = 20,
		odens = 2,
		seed = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		save_log_lik = TRUE,
		log_lik_method = "observed_ghk"
	)))
	ll = read_log_lik(fit)
	expect_equal(dim(ll), c(10L, 3L * n * (n - 1L)))
	expect_true(all(is.finite(ll) & ll <= 0))
	# the scores are approximate diagnostics, which loo() refuses
	skip_if_not_installed("loo")
	expect_error(loo::loo(fit), class = "lame_loo_factorization_unsupported")
})

test_that("a directed Poisson dynamic_uv penalised fit that stops short is finished by the joint refinement", {
	skip_on_cran()
	set.seed(6); n = 12; nm = sprintf("a%02d", 1:n)
	Y = lapply(1:4, function(t) { m = matrix(rpois(n * n, 1.5), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	fit = function(max_iter) suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
		Y,
		family = "poisson",
		R = 1,
		dynamic_uv = TRUE,
		bias_correct = "none",
		max_iter = max_iter,
		verbose = FALSE
	)))
	short = fit(5)
	r = short$diagnostics$joint_refinement
	expect_true(isTRUE(short$converged))
	expect_true(isTRUE(r$accepted)); expect_identical(r$termination, "model score reached")
	expect_lte(r$objective_after, r$objective_before + 1e-8 * abs(r$objective_before))
	expect_lte(r$max_score, r$threshold)
	expect_true(isTRUE(short$diagnostics$stationarity$pass))
	# a fit that reaches its criterion on its own sweeps is not refined, and both land on the same optimum
	long = fit(50)
	expect_null(long$diagnostics$joint_refinement)
	expect_true(isTRUE(long$converged))
	expect_equal(unlist(short$coefficients), unlist(long$coefficients), tolerance = 1e-3)
})
