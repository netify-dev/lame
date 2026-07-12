# in-loop per-actor slope sampler
#
# with pairwise contrast ffbs. zero-sum is preserved per period; ar(1)
# prior on each actor's path; hyperparameters are separate from the
# per-block dynamic_beta hyperparameters.

skip_on_cran()

.make_per_actor_fixture = function(n = 20L, T_per = 5L,
                                      beta_pop = 0.3, sigma_noise = 0.15,
                                      seed = 7L) {
	set.seed(seed)
	actor_nms = sprintf("a%02d", seq_len(n))
	X_list = lapply(seq_len(T_per), function(t) {
		X = matrix(rnorm(n*n), n, n)
		rownames(X) = colnames(X) = actor_nms
		array(X, c(n, n, 1), dimnames = list(actor_nms, actor_nms, "X1"))
	})
	theta_true = matrix(rnorm(n*T_per, 0, 0.5), n, T_per,
	                     dimnames = list(actor_nms, paste0("t", 1:T_per)))
	theta_true = sweep(theta_true, 2L, colMeans(theta_true), "-")
	Y_list = vector("list", T_per)
	for (t in seq_len(T_per)) {
		X_t = X_list[[t]][, , 1]
		eta = beta_pop * X_t + theta_true[, t] * X_t
		Y = eta + matrix(rnorm(n*n, 0, sigma_noise), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actor_nms
		Y_list[[t]] = Y
	}
	list(Y = Y_list, X = X_list, theta_true = theta_true,
	     beta_pop = beta_pop, actor_nms = actor_nms)
}

# -- structural --------------------------------------------------------

test_that("dynamic_beta_per_actor records on fit + dim checks", {
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "draws",
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_identical(fit$dynamic_beta_per_actor, "row")
	expect_true(!is.null(fit$THETA_ACTOR))
	expect_equal(dim(fit$THETA_ACTOR)[2L], 20L)  # n_actors
	expect_equal(dim(fit$THETA_ACTOR)[3L], 5L)   # t_per
	expect_true(all(is.finite(fit$THETA_ACTOR)))
})

# -- zero-sum constraint per period -----------------------------------

test_that("pairwise contrast preserves sum_i theta_{i,t} = 0 exactly", {
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "draws",
		nscan = 80, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	# sum over actors should be ~0 for every draw and every period
	zero_sum_check = apply(fit$THETA_ACTOR, c(1, 3), sum)
	expect_lt(max(abs(zero_sum_check)), 1e-8)
})

# -- substantive recovery ---------------------------------------------

test_that("per-actor sampler recovers truth with high correlation", {
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "draws",
		nscan = 600, burn = 200, odens = 10,
		verbose = FALSE, plot = FALSE)))
	rho = cor(as.vector(fit$theta_actor_mean), as.vector(fx$theta_true))
	expect_gt(rho, 0.85)  # 0.99 in longer runs; relax for fast tests
	expect_lt(sqrt(mean((fit$theta_actor_mean - fx$theta_true)^2)), 0.30)
})

# -- byte-identical default -------------------------------------------

test_that("dynamic_beta_per_actor = NULL leaves the base sampler unchanged", {
	# null keeps the per-actor block inactive
	fx = .make_per_actor_fixture()
	set.seed(1); fit_null = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		nscan = 30, burn = 10, odens = 5,
		verbose = FALSE, plot = FALSE)))
	set.seed(1); fit_explicit_null = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = NULL,
		nscan = 30, burn = 10, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_identical(fit_null$BETA, fit_explicit_null$BETA)
	expect_identical(fit_null$EZ,   fit_explicit_null$EZ)
})

# -- storage modes ----------------------------------------------------

test_that("keep_per_actor = 'summary' stores streaming mean + sd, no full draws", {
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "summary",
		nscan = 60, burn = 20, odens = 5,
		verbose = FALSE, plot = FALSE)))
	expect_null(fit$THETA_ACTOR)
	expect_true(!is.null(fit$theta_actor_mean))
	expect_true(!is.null(fit$theta_actor_sd))
	expect_equal(dim(fit$theta_actor_mean), c(20L, 5L))
})

# -- validation -------------------------------------------------------

test_that("invalid dynamic_beta_per_actor aborts", {
	fx = .make_per_actor_fixture()
	expect_error(suppressWarnings(suppressMessages(lame(
		fx$Y, fx$X, family = "normal", R = 0,
		dynamic_beta_per_actor = "magic",
		nscan = 30, burn = 10, odens = 5, verbose = FALSE, plot = FALSE))),
		"must be one of")
})
