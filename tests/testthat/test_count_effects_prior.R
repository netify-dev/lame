# default effects prior of the count sampler (.ame_pois_effects_prior): one actor's weight at the
# fast route's effects covariance, the log-rate scale where that fit fails

sim_count = function(n, nc = n, seed = 1, sym = FALSE) {
	set.seed(seed)
	rn = sprintf("a%02d", 1:n); cn = if (nc == n) rn else sprintf("c%02d", 1:nc)
	a = rnorm(n, 0, 0.5); b = if (sym) a else rnorm(nc, 0, 0.5)
	eta = -0.3 + outer(a, b, "+")
	y = matrix(rpois(n * nc, exp(eta)), n, nc)
	if (sym) y[lower.tri(y)] = t(y)[lower.tri(y)]
	if (nc == n) diag(y) = NA
	dimnames(y) = list(rn, cn)
	y
}
pos_def = function(S) all(is.finite(S)) && all(eigen(S, symmetric = TRUE, only.values = TRUE)$values > 0)

test_that("the count prior has weight one and a positive definite centre in every mode", {
	skip_on_cran()
	y = sim_count(15)
	p = lame:::.ame_pois_effects_prior(NULL, y, NULL, NULL, NULL, 0, "unipartite", FALSE, FALSE)
	expect_identical(p$eta0, 1)
	expect_true(pos_def(p$Sab0))
	ys = sim_count(15, sym = TRUE)
	p = lame:::.ame_pois_effects_prior(NULL, ys, NULL, NULL, NULL, 0, "unipartite", TRUE, FALSE)
	expect_identical(p$eta0, 1)
	expect_true(pos_def(p$Sab0))
	expect_equal(p$Sab0[1, 2], 0)
	yb = sim_count(15, 12)
	p = lame:::.ame_pois_effects_prior(NULL, yb, NULL, NULL, NULL, 0, "bipartite", FALSE, FALSE)
	expect_identical(p$eta0, 1)
	expect_true(pos_def(p$Sab0))
	expect_equal(p$Sab0[1, 2], 0)
	yp = list(t1 = sim_count(15, seed = 2), t2 = sim_count(15, seed = 3))
	p = lame:::.ame_pois_effects_prior(NULL, yp, NULL, NULL, NULL, 0, "unipartite", FALSE, TRUE)
	expect_identical(p$eta0, 1)
	expect_true(pos_def(p$Sab0))
})

test_that("a supplied effects prior is kept as given", {
	skip_on_cran()
	y = sim_count(12)
	expect_false(lame:::.lame_ab_prior_given(NULL))
	expect_true(lame:::.lame_ab_prior_given(list(Sab0 = diag(2))))
	expect_true(lame:::.lame_ab_prior_given(list(eta0 = 3)))
	fit = ame(y, R = 0, family = "poisson", prior = list(Sab0 = diag(2) * 2, eta0 = 3), nscan = 100, burn = 20,
		odens = 5, verbose = FALSE, gof = FALSE)
	expect_equal(fit$prior$Sab0, diag(2) * 2)
	expect_identical(fit$prior$eta0, 3)
})

test_that("the fallback is the log-rate scale at weight one", {
	skip_on_cran()
	y = sim_count(12)
	# a design the fast route cannot take makes the inner fit fail
	p = lame:::.ame_pois_effects_prior(NULL, y, "not a design", NULL, NULL, 0, "unipartite", FALSE, FALSE)
	expect_identical(p$eta0, 1)
	expect_equal(p$Sab0, diag(2) * lame:::.ame_pois_effect_scale(y))
})

test_that("the count prior leaves the global random number state unchanged", {
	skip_on_cran()
	y = sim_count(12)
	set.seed(9); r0 = runif(3)
	set.seed(9)
	p = lame:::.ame_pois_effects_prior(NULL, y, NULL, NULL, NULL, 0, "unipartite", FALSE, FALSE)
	expect_identical(runif(3), r0)
})

test_that("a small count sampler fit returns finite effect variances", {
	skip_on_cran()
	y = sim_count(12)
	fit = ame(y, R = 0, family = "poisson", nscan = 100, burn = 20, odens = 5, verbose = FALSE, gof = FALSE)
	expect_true(all(is.finite(fit$VC[, "va"])) && all(fit$VC[, "va"] > 0))
	expect_true(all(is.finite(fit$VC[, "vb"])) && all(fit$VC[, "vb"] > 0))
	expect_identical(fit$prior$eta0, 1)
})
