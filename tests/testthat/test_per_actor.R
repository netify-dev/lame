# per-actor dynamic beta with exact Gaussian conditioning on the sum-to-zero constraint:
# the default "center" contract, exact zero sums under "exact_center"
# and recovery on a known dgp

test_that("per_actor_identifiability = 'center' is byte-identical to default", {
	skip_on_cran()
	set.seed(2026)
	n = 8; Tn = 3
	Xdyad_list = lapply(seq_len(Tn), function(t) {
		x = array(rnorm(n * n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(Tn), function(t) {
		mat = matrix(rnorm(n * n), n, n); diag(mat) = NA; mat
	})
	f1 = lame(
		Y_list,
		Xdyad = Xdyad_list,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		per_actor_covariate_idx = 1L,
		keep_per_actor = "draws",
		nscan = 20,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		seed = 7
	)
	f2 = lame(
		Y_list,
		Xdyad = Xdyad_list,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		per_actor_covariate_idx = 1L,
		per_actor_identifiability = "center",
		keep_per_actor = "draws",
		nscan = 20,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		seed = 7
	)
	expect_identical(f1$THETA_ACTOR, f2$THETA_ACTOR)
	expect_identical(f1$RHO_ACTOR, f2$RHO_ACTOR)
})

test_that("exact_center preserves zero-sum exactly across all draws", {
	skip_on_cran()
	set.seed(2027)
	n = 10; Tn = 4
	Xdyad_list = lapply(seq_len(Tn), function(t) {
		x = array(rnorm(n * n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(Tn), function(t) {
		mat = matrix(rnorm(n * n), n, n); diag(mat) = NA; mat
	})
	fit = lame(
		Y_list,
		Xdyad = Xdyad_list,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		per_actor_covariate_idx = 1L,
		per_actor_identifiability = "exact_center",
		keep_per_actor = "draws",
		nscan = 30,
		burn = 15,
		odens = 5,
		verbose = FALSE,
		seed = 13
	)
	expect_equal(fit$per_actor_identifiability, "exact_center")
	# zero-sum to machine precision across every draw and period
	per_draw_per_period_sums = apply(fit$THETA_ACTOR, c(1, 3), sum)
	expect_lt(max(abs(per_draw_per_period_sums)), 1e-10)
})

# the in-loop per-actor slope sampler of lame(dynamic_beta_per_actor = )

.make_per_actor_fixture = function(
	n = 20L,
	T_per = 5L,
	beta_pop = 0.3,
	sigma_noise = 0.15,
	seed = 7L
) {
	set.seed(seed)
	actor_nms = sprintf("a%02d", seq_len(n))
	X_list = lapply(seq_len(T_per), function(t) {
		X = matrix(rnorm(n * n), n, n)
		rownames(X) = colnames(X) = actor_nms
		array(X, c(n, n, 1), dimnames = list(actor_nms, actor_nms, "X1"))
	})
	theta_true = matrix(
		rnorm(n * T_per, 0, 0.5),
		n,
		T_per,
		dimnames = list(actor_nms, paste0("t", 1:T_per))
	)
	theta_true = sweep(theta_true, 2L, colMeans(theta_true), "-")
	Y_list = vector("list", T_per)
	for (t in seq_len(T_per)) {
		X_t = X_list[[t]][, , 1]
		eta = beta_pop * X_t + theta_true[, t] * X_t
		Y = eta + matrix(rnorm(n * n, 0, sigma_noise), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = actor_nms
		Y_list[[t]] = Y
	}
	list(
		Y = Y_list,
		X = X_list,
		theta_true = theta_true,
		beta_pop = beta_pop,
		actor_nms = actor_nms
	)
}

test_that("pairwise contrast preserves sum_i theta_{i,t} = 0 exactly", {
	skip_on_cran()
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y,
		fx$X,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "draws",
		nscan = 80,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	# the sum over actors is about 0 for every draw and period
	zero_sum_check = apply(fit$THETA_ACTOR, c(1, 3), sum)
	expect_lt(max(abs(zero_sum_check)), 1e-8)
})

test_that("per-actor sampler recovers truth with high correlation", {
	skip_on_cran()
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y,
		fx$X,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "draws",
		nscan = 600,
		burn = 200,
		odens = 10,
		verbose = FALSE,
		plot = FALSE
	)))
	rho = cor(as.vector(fit$theta_actor_mean), as.vector(fx$theta_true))
	expect_gt(rho, 0.85) # about 0.99 in longer runs
	expect_lt(sqrt(mean((fit$theta_actor_mean - fx$theta_true)^2)), 0.30)
})

test_that("keep_per_actor = 'summary' stores streaming mean + sd, no full draws", {
	skip_on_cran()
	fx = .make_per_actor_fixture()
	set.seed(1); fit = suppressWarnings(suppressMessages(lame(
		fx$Y,
		fx$X,
		family = "normal",
		R = 0,
		dynamic_beta_per_actor = "row",
		keep_per_actor = "summary",
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_null(fit$THETA_ACTOR)
	expect_true(!is.null(fit$theta_actor_mean))
	expect_true(!is.null(fit$theta_actor_sd))
	expect_equal(dim(fit$theta_actor_mean), c(20L, 5L))
})

test_that("invalid dynamic_beta_per_actor aborts", {
	skip_on_cran()
	fx = .make_per_actor_fixture()
	expect_error(
		suppressWarnings(suppressMessages(lame(
			fx$Y,
			fx$X,
			family = "normal",
			R = 0,
			dynamic_beta_per_actor = "magic",
			nscan = 30,
			burn = 10,
			odens = 5,
			verbose = FALSE,
			plot = FALSE
		))),
		"must be one of"
	)
})
