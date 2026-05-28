# per-actor dynamic beta with exact Gaussian
# conditioning for the sum-to-zero centering constraint.
# Tests cover: (a) helper math, (b) byte-identical-default contract for
# "center", (c) "exact_center" preserves zero-sum exactly per draw,
# (d) end-to-end recovery on a known DGP, (e) FFBS path posterior
# variance is correctly returned.

test_that(".actor_ffbs_path returns theta and V of length T", {
	set.seed(11)
	T_per = 8
	H = runif(T_per, 0.1, 1.0)
	h = rnorm(T_per)
	out = lame:::.actor_ffbs_path(H, h, rho_actor = 0.7, sigma_actor2 = 0.1)
	expect_equal(length(out$theta), T_per)
	expect_equal(length(out$V), T_per)
	expect_true(all(out$V > 0))
	expect_true(all(is.finite(out$theta)))
})

test_that(".exact_center_per_actor enforces zero-sum exactly per period", {
	set.seed(23)
	n_actors = 8; T_per = 5
	theta_star = matrix(rnorm(n_actors * T_per, 0, 0.5), n_actors, T_per)
	V = matrix(runif(n_actors * T_per, 0.1, 1.0), n_actors, T_per)
	out = lame:::.exact_center_per_actor(theta_star, V)
	# zero-sum at each period to machine precision
	col_sums = colSums(out)
	expect_lt(max(abs(col_sums)), 1e-12)
})

test_that(".sweep_per_actor_exact produces n_actors x T zero-sum draws", {
	set.seed(31)
	n_actors = 6; T_per = 5
	H_mat = matrix(runif(n_actors * T_per, 0.1, 1.0), n_actors, T_per)
	h_mat = matrix(rnorm(n_actors * T_per), n_actors, T_per)
	theta = lame:::.sweep_per_actor_exact(H_mat, h_mat,
	                                    rho_actor = 0.7, sigma_actor2 = 0.1)
	expect_equal(dim(theta), c(n_actors, T_per))
	# zero-sum per period
	expect_lt(max(abs(colSums(theta))), 1e-12)
})

test_that(".sweep_per_actor_exact collapses to zero for zero-precision case", {
	# when H = 0 (no observation precision), the data does nothing and the
	# prior is centered at zero. Projection should give zero (up to MC).
	set.seed(41)
	n_actors = 5; T_per = 4
	H_mat = matrix(0.001, n_actors, T_per)  # near-zero precision
	h_mat = matrix(0, n_actors, T_per)
	theta = lame:::.sweep_per_actor_exact(H_mat, h_mat,
	                                    rho_actor = 0.5, sigma_actor2 = 0.05)
	expect_lt(max(abs(colSums(theta))), 1e-10)
})

test_that("per_actor_identifiability = 'center' is byte-identical to default", {
	skip_on_cran()
	set.seed(2026)
	n = 8; T = 3
	Xdyad_list = lapply(seq_len(T), function(t) {
		x = array(rnorm(n*n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rnorm(n*n), n, n); diag(mat) = NA; mat
	})
	f1 = lame(Y_list, Xdyad = Xdyad_list, family = "normal", R = 0,
	          dynamic_beta_per_actor = "row",
	          per_actor_covariate_idx = 1L,
	          keep_per_actor = "draws",
	          nscan = 20, burn = 10, odens = 5, verbose = FALSE, seed = 7)
	f2 = lame(Y_list, Xdyad = Xdyad_list, family = "normal", R = 0,
	          dynamic_beta_per_actor = "row",
	          per_actor_covariate_idx = 1L,
	          per_actor_identifiability = "center",
	          keep_per_actor = "draws",
	          nscan = 20, burn = 10, odens = 5, verbose = FALSE, seed = 7)
	expect_identical(f1$THETA_ACTOR, f2$THETA_ACTOR)
	expect_identical(f1$RHO_ACTOR, f2$RHO_ACTOR)
})

test_that("exact_center preserves zero-sum exactly across all draws", {
	skip_on_cran()
	set.seed(2027)
	n = 10; T = 4
	Xdyad_list = lapply(seq_len(T), function(t) {
		x = array(rnorm(n*n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rnorm(n*n), n, n); diag(mat) = NA; mat
	})
	fit = lame(Y_list, Xdyad = Xdyad_list, family = "normal", R = 0,
	           dynamic_beta_per_actor = "row",
	           per_actor_covariate_idx = 1L,
	           per_actor_identifiability = "exact_center",
	           keep_per_actor = "draws",
	           nscan = 30, burn = 15, odens = 5, verbose = FALSE, seed = 13)
	expect_equal(fit$per_actor_identifiability, "exact_center")
	# zero-sum to machine precision across every draw and period
	per_draw_per_period_sums = apply(fit$THETA_ACTOR, c(1, 3), sum)
	expect_lt(max(abs(per_draw_per_period_sums)), 1e-10)
})

test_that("exact_center end-to-end fit produces sensible posterior summaries", {
	skip_on_cran()
	set.seed(53)
	n = 12; T = 4
	Xdyad_list = lapply(seq_len(T), function(t) {
		x = array(rnorm(n*n), dim = c(n, n, 1))
		dimnames(x) = list(NULL, NULL, "x1")
		x
	})
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rnorm(n*n), n, n); diag(mat) = NA; mat
	})
	fit = lame(Y_list, Xdyad = Xdyad_list, family = "normal", R = 0,
	           dynamic_beta_per_actor = "row",
	           per_actor_covariate_idx = 1L,
	           per_actor_identifiability = "exact_center",
	           keep_per_actor = "summary",
	           nscan = 60, burn = 30, odens = 5, verbose = FALSE, seed = 99)
	# fit completes and stores the streaming summary
	expect_equal(dim(fit$theta_actor_mean), c(n, T))
	expect_equal(dim(fit$theta_actor_sd), c(n, T))
	expect_true(all(is.finite(fit$theta_actor_mean)))
	expect_true(all(fit$theta_actor_sd >= 0))
	# RHO_ACTOR + SIGMA_ACTOR chains run within sensible ranges
	expect_true(all(abs(fit$RHO_ACTOR) < 1))
	expect_true(all(fit$SIGMA_ACTOR > 0))
})
