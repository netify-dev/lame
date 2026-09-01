skip_on_cran()

# posterior utility functions
library(lame)

test_that("posterior_options creates correct structure", {
	# test default options
	opts1 = posterior_options()
	expect_equal(opts1$save_UV, FALSE)
	expect_equal(opts1$save_ab, FALSE)
	expect_equal(opts1$thin_UV, 10)
	expect_equal(opts1$thin_ab, 10)
	
	# test custom options
	opts2 = posterior_options(save_UV = TRUE, save_ab = TRUE, 
														thin_UV = 5, thin_ab = 20)
	expect_equal(opts2$save_UV, TRUE)
	expect_equal(opts2$save_ab, TRUE)
	expect_equal(opts2$thin_UV, 5)
	expect_equal(opts2$thin_ab, 20)
})

test_that("simulate_posterior works for beta coefficients", {
	set.seed(6886)
	n = 20
	
	# generate simple network
	X = array(rnorm(n * n * 2), c(n, n, 2))
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	
	# fit model
	fit = ame(Y, Xdyad = X, R = 0, 
						burn = 100, nscan = 500, verbose = FALSE)
	
	# simulate posterior for beta (may warn about fewer samples than requested)
	beta_post = suppressWarnings(
		simulate_posterior(fit, "beta", n_samples = 100)
	)
	
	# check dimensions (may have fewer samples than requested)
	expect_lte(nrow(beta_post), 100)
	expect_gt(nrow(beta_post), 0)
	expect_equal(ncol(beta_post), ncol(fit$BETA))
	
	# check values are numeric
	expect_true(all(is.finite(beta_post)))
	
	# check mean is close to mcmc mean
	beta_mcmc_mean = colMeans(fit$BETA)
	beta_sim_mean = colMeans(beta_post)
	expect_equal(beta_sim_mean, beta_mcmc_mean, tolerance = 0.5)
})

test_that("simulate_posterior errors for ab without saved samples", {
	set.seed(6886)
	n = 20

	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA

	# fit model without saving ab samples
	fit = ame(Y, R = 0, rvar = TRUE, cvar = TRUE,
						burn = 100, nscan = 500, verbose = FALSE)

	# should error — no saved samples available
	expect_error(simulate_posterior(fit, "ab", n_samples = 50),
		"saved MCMC samples")
})

test_that("simulate_posterior errors for UV without saved samples", {
	set.seed(6886)
	n = 20
	R = 2

	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA

	# fit model without saving uv samples
	fit = ame(Y, R = R, burn = 100, nscan = 500, verbose = FALSE)

	# should error — no saved samples available
	expect_error(simulate_posterior(fit, "UV", n_samples = 10),
		"saved MCMC samples")
})

test_that("posterior_quantiles computes correct intervals", {
	set.seed(6886)
	n = 20
	
	X = array(rnorm(n * n * 2), c(n, n, 2))
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	
	# fit model
	fit = ame(Y, Xdyad = X, R = 0,
						burn = 100, nscan = 500, verbose = FALSE)
	
	# get quantiles for beta
	beta_quants = posterior_quantiles(fit, "beta", probs = c(0.025, 0.5, 0.975))
	
	# check dimensions
	expect_equal(nrow(beta_quants), 3)  # 3 quantiles
	expect_equal(ncol(beta_quants), ncol(fit$BETA))  # number of parameters
	
	# check ordering (lower < median < upper)
	expect_true(all(beta_quants[1,] <= beta_quants[2,]))
	expect_true(all(beta_quants[2,] <= beta_quants[3,]))
	
	# check median is close to mean for approximately normal posteriors
	beta_means = colMeans(fit$BETA)
	beta_medians = beta_quants[2,]
	expect_equal(beta_medians, beta_means, tolerance = 0.5)
})

test_that("posterior samples are saved when requested", {
	set.seed(6886)
	n = 15
	R = 2
	
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	
	# fit model with posterior saving
	opts = posterior_options(save_UV = TRUE, save_ab = TRUE, 
													 thin_UV = 5, thin_ab = 5)
	
	fit = ame(Y, R = R, rvar = TRUE, cvar = TRUE,
						burn = 50, nscan = 200, odens = 5,
						posterior_opts = opts, verbose = FALSE)
	
	# check that samples were saved
	if(!is.null(opts) && opts$save_UV) {
		expect_false(is.null(fit$U_samples))
		expect_false(is.null(fit$V_samples))
		
		# check dimensions
		n_expected = floor(200 / (5 * 5))  # nscan / (odens * thin)
		expect_equal(dim(fit$U_samples)[1], n)
		expect_equal(dim(fit$U_samples)[2], R)
		expect_lte(dim(fit$U_samples)[3], n_expected + 1)
	}
	
	if(!is.null(opts) && opts$save_ab) {
		expect_false(is.null(fit$a_samples))
		expect_false(is.null(fit$b_samples))
		
		# check dimensions
		expect_equal(nrow(fit$a_samples), n)
		expect_equal(nrow(fit$b_samples), n)
	}
})

test_that("simulate_posterior works with saved samples", {
	set.seed(6886)
	n = 15
	R = 2
	
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	
	# fit model with posterior saving
	opts = posterior_options(save_UV = TRUE, thin_UV = 10)
	
	fit = ame(Y, R = R, burn = 50, nscan = 200, odens = 5,
						posterior_opts = opts, verbose = FALSE)
	
	# simulate from saved samples
	if(!is.null(fit$U_samples) && !is.null(fit$V_samples)) {
		UV_post = suppressWarnings(
			simulate_posterior(fit, "UV", n_samples = 5)
		)
		
		# should use saved samples
		expect_equal(dim(UV_post)[1:2], c(n, n))
		expect_lte(dim(UV_post)[3], 5)
	}
})

test_that("posterior functions work for bipartite networks", {
	set.seed(6886)
	nA = 10
	nB = 12
	
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# fit bipartite model
	fit = ame(Y, mode = "bipartite", R_row = 1, R_col = 1,
						burn = 50, nscan = 200, verbose = FALSE)
	
	# check we can simulate posteriors
	beta_post = suppressWarnings(
		simulate_posterior(fit, "beta", n_samples = 20)
	)
	expect_lte(nrow(beta_post), 20)  # may have fewer if not enough mcmc samples
	expect_gt(nrow(beta_post), 0)    # but should have some
	
	# check quantiles work
	beta_quants = posterior_quantiles(fit, "beta")
	expect_equal(nrow(beta_quants), 3)  # default 3 quantiles
})

test_that("simulate_Y_posterior generates reasonable predictions", {
	set.seed(6886)
	n = 20
	
	# simple model
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	
	fit = ame(Y, R = 0, burn = 100, nscan = 300, verbose = FALSE)
	
	# simulate posterior predictive
	Y_post = simulate_posterior(fit, "Y", n_samples = 10)
	
	# check dimensions
	expect_equal(dim(Y_post), c(n, n, 10))
	
	# check missing values are preserved
	for(i in 1:10) {
		expect_true(all(is.na(diag(Y_post[,,i]))))
	}
	
	# check values are reasonable (similar scale to original)
	Y_mean = mean(Y, na.rm = TRUE)
	Y_sd = sd(Y, na.rm = TRUE)
	Y_post_mean = mean(Y_post, na.rm = TRUE)
	Y_post_sd = sd(Y_post, na.rm = TRUE)
	
	expect_equal(Y_post_mean, Y_mean, tolerance = 1)
	expect_equal(Y_post_sd, Y_sd, tolerance = 1)
})

test_that("simulate_posterior uses every stored draw by default", {
	skip_on_cran()
	set.seed(11); n = 14
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	fit = ame(Y, R = 2, symmetric = TRUE, family = "binary", burn = 20,
	          nscan = 120, odens = 1, verbose = FALSE, gof = FALSE, seed = 5,
	          posterior_opts = list(save_UV = TRUE, thin_UV = 1))
	n_saved = dim(fit$U_samples)[3]

	# default returns the full posterior sample, not a 100-draw subsample
	expect_equal(dim(simulate_posterior(fit, "UV"))[3], n_saved)
	# an explicit smaller request still subsamples
	expect_equal(dim(simulate_posterior(fit, "UV", n_samples = 25))[3], 25)
	# asking for more than exist warns and returns what is stored
	expect_warning(uv <- simulate_posterior(fit, "UV", n_samples = n_saved + 50))
	expect_equal(dim(uv)[3], n_saved)
	# the draw mean reproduces the stored posterior-mean product
	expect_equal(apply(simulate_posterior(fit, "UV"), c(1, 2), mean),
	             unname(fit$ULUPM), tolerance = 1e-8)
})

# --- regression: dynamic-fit posterior extraction (1.3.5) --------------------
# dynamic fits store latent draws in the 4-D U_draws/V_draws (not U_samples)
# and a 3-D BETA; simulate_posterior() returns draws for both and the credible
# forecasts keep their actor names.

.mk_sq_panel = function(n, T = 3, sym = FALSE, fam = "binary") {
	lapply(seq_len(T), function(t) {
		Y = matrix(rnorm(n * n), n, n)
		if (sym) Y = (Y + t(Y)) / 2
		if (fam == "binary") Y = (Y > 0) * 1
		diag(Y) = NA
		rownames(Y) = colnames(Y) = sprintf("a%02d", seq_len(n))
		Y
	})
}

test_that("simulate_posterior('UV') works on directed dynamic_uv fits", {
	set.seed(6886)
	f = lame(.mk_sq_panel(8), family = "binary", R = 2, dynamic_uv = TRUE,
	         burn = 20, nscan = 60, odens = 2, verbose = FALSE, plot = FALSE,
	         posterior_opts = list(save_UV = TRUE))
	uv = simulate_posterior(f, "UV", n_samples = 6)
	expect_length(dim(uv), 4L)          # actor x actor x time x draw
	expect_equal(dim(uv)[1:2], c(8L, 8L))
	expect_equal(dim(uv)[4], 6L)
	expect_true(all(is.finite(uv)))
})

test_that("simulate_posterior('UV') is symmetric per period for symmetric dynamic fits", {
	set.seed(1)
	f = lame(.mk_sq_panel(8, sym = TRUE, fam = "normal"), family = "normal",
	         symmetric = TRUE, R = 2, dynamic_uv = TRUE, burn = 20, nscan = 60,
	         odens = 2, verbose = FALSE, plot = FALSE,
	         posterior_opts = list(save_UV = TRUE))
	uv = simulate_posterior(f, "UV", n_samples = 6)
	m = apply(uv[, , 1, ], c(1, 2), mean)
	expect_lt(max(abs(m - t(m)), na.rm = TRUE), 1e-6)
})

test_that("simulate_posterior('UV') on bipartite dynamic fits is nA x nB x T x draw", {
	set.seed(2)
	Y = lapply(1:3, function(t) {
		m = matrix(rnorm(8 * 6), 8, 6)
		rownames(m) = sprintf("r%02d", 1:8); colnames(m) = sprintf("c%02d", 1:6); m
	})
	f = lame(Y, mode = "bipartite", family = "normal", R = 2, dynamic_uv = TRUE,
	         burn = 20, nscan = 60, odens = 2, verbose = FALSE, plot = FALSE,
	         posterior_opts = list(save_UV = TRUE))
	uv = simulate_posterior(f, "UV", n_samples = 5)
	expect_equal(dim(uv)[1:2], c(8L, 6L))
	expect_equal(dim(uv)[3], 3L)
})

test_that("simulate_posterior('beta') handles a 3-D dynamic_beta array", {
	set.seed(3)
	Y = .mk_sq_panel(8)
	X = lapply(1:3, function(t) { a = array(rnorm(64), c(8, 8, 1)); dimnames(a)[[3]] = "x"; a })
	f = lame(Y, Xdyad = X, family = "binary", R = 2, dynamic_beta = "dyad",
	         burn = 20, nscan = 60, odens = 2, verbose = FALSE, plot = FALSE)
	b = simulate_posterior(f, "beta", n_samples = 12)
	expect_length(dim(b), 3L)           # draws x coefs x time
	expect_equal(dim(b)[1], 12L)
	expect_equal(dim(b)[3], 3L)
})

test_that("credible-interval forecasts carry actor dimnames", {
	set.seed(4)
	f = lame(.mk_sq_panel(8), family = "binary", R = 2, dynamic_uv = TRUE,
	         burn = 20, nscan = 60, odens = 2, verbose = FALSE, plot = FALSE)
	r = predict(f, h = 2, type = "response", interval = "credible")
	expect_false(is.null(rownames(r[[1]]$median)))
	expect_false(is.null(colnames(r[[1]]$lower)))
	expect_false(is.null(rownames(r[[2]]$upper)))
})

test_that("bipartite simulate_posterior('UV') includes G: draw-mean reproduces UVPM", {
	set.seed(11)
	mk = function() {
		m = matrix(rnorm(8 * 6), 8, 6)
		rownames(m) = sprintf("r%02d", 1:8); colnames(m) = sprintf("c%02d", 1:6); m
	}
	po = list(save_UV = TRUE)
	# static ame()
	fa = ame(mk(), mode = "bipartite", family = "normal", R = 2, burn = 20,
	         nscan = 60, odens = 2, verbose = FALSE, plot = FALSE, posterior_opts = po)
	expect_false(is.null(fa$G_samples))
	ma = apply(simulate_posterior(fa, "UV"), c(1, 2), mean)
	expect_equal(unname(ma), unname(fa$UVPM), tolerance = 1e-8)
	# static lame()
	fl = lame(replicate(3, mk(), simplify = FALSE), mode = "bipartite",
	          family = "normal", R = 2, burn = 20, nscan = 60, odens = 2,
	          verbose = FALSE, plot = FALSE, posterior_opts = po)
	ml = apply(simulate_posterior(fl, "UV"), c(1, 2), mean)
	expect_equal(unname(ml), unname(fl$UVPM), tolerance = 1e-8)
	# dynamic_uv lame(): per-period product matches the stored per-period mean
	fd = lame(replicate(3, mk(), simplify = FALSE), mode = "bipartite",
	          family = "normal", R = 2, dynamic_uv = TRUE, burn = 20, nscan = 60,
	          odens = 2, verbose = FALSE, plot = FALSE, posterior_opts = po)
	ud = simulate_posterior(fd, "UV")
	expect_equal(dim(ud)[1:3], c(8L, 6L, 3L))
	md = apply(ud[, , 1, ], c(1, 2), mean)
	ref = if (length(dim(fd$UVPM)) == 3L) fd$UVPM[, , 1] else fd$UVPM[[1]]
	expect_gt(cor(c(md), c(ref)), 0.999)
})
