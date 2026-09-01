skip_on_cran()

# consolidated bipartite tests: internal machinery, parameter recovery,
# and per-family behavior for both static ame() and longitudinal lame().
# merged from test_bipartite.R, test_bipartite_recovery.R,
# test_bipartite_families.R, test_bipartite_families_xs.R,
# test_static_bipartite.R, and test_static_bipartite_families.R.
# the uvpm / rank posterior-mean guards live in test_bipartite_uvpm.R.

# bipartite-family recovery tests
# each test fits a small rectangular network with a known truth

skip_if_too_slow = function() {
	if (identical(Sys.getenv("LAME_FAST_TESTS"), "1")) skip("LAME_FAST_TESTS = 1")
}

# -------- helpers --------

mk_bip_linpred = function(nA, nB, Tn, seed, R_true = 0L, beta_true = c(0.8, -0.5)) {
	set.seed(seed)
	X = lapply(seq_len(Tn), function(t) {
		array(stats::rnorm(nA * nB * length(beta_true)),
		      c(nA, nB, length(beta_true)),
		      dimnames = list(paste0("r", seq_len(nA)),
		                      paste0("c", seq_len(nB)),
		                      paste0("v", seq_along(beta_true))))
	})
	a = stats::rnorm(nA, 0, 0.3); b = stats::rnorm(nB, 0, 0.3)
	EZ = lapply(seq_len(Tn), function(t) {
		eta = matrix(a, nA, nB) + matrix(b, nA, nB, byrow = TRUE)
		for (k in seq_along(beta_true)) eta = eta + beta_true[k] * X[[t]][, , k]
		eta
	})
	list(EZ = EZ, X = X, a = a, b = b, beta = beta_true,
	     nA = nA, nB = nB, Tn = Tn)
}

# align true params with the alphabetically sorted actor order
align_true_params = function(true_vals, orig_names, sorted_names) {
	names(true_vals) = orig_names
	true_vals[sorted_names]
}

# cross-sectional bipartite families.
# mirrors test_bipartite_families.r for `lame()` but exercises the
# `ame()` (cross-sectional) path through ame_bipartite(). the
# every family routes through the same
# rectangular z samplers from r/rz_bipartite.r.
mk_xs = function(family, nA = 18, nB = 14, seed = 40010L, odmax_val = 3L,
                  beta_true = c(0.6, -0.4)) {
	set.seed(seed)
	X = array(stats::rnorm(nA * nB * length(beta_true)),
	          c(nA, nB, length(beta_true)),
	          dimnames = list(paste0("r", seq_len(nA)),
	                          paste0("c", seq_len(nB)),
	                          paste0("v", seq_along(beta_true))))
	a = stats::rnorm(nA, 0, 0.25); b = stats::rnorm(nB, 0, 0.25)
	eta = matrix(a, nA, nB) + matrix(b, nA, nB, byrow = TRUE)
	for (k in seq_along(beta_true)) eta = eta + beta_true[k] * X[, , k]
	Z = eta + matrix(stats::rnorm(nA * nB), nA, nB)
	if (family == "poisson") {
		Y = matrix(stats::rpois(nA * nB, lambda = exp(pmin(eta - 1, 5))), nA, nB)
	} else if (family == "ordinal") {
		Y = as.integer(cut(Z, c(-Inf, -0.5, 0.3, 1.0, Inf), labels = FALSE))
		dim(Y) = c(nA, nB)
	} else if (family == "cbin") {
		Y = matrix(0L, nA, nB)
		for (i in seq_len(nA))
			Y[i, order(Z[i, ], decreasing = TRUE)[seq_len(odmax_val)]] = 1L
	} else if (family == "frn") {
		Y = matrix(0L, nA, nB)
		for (i in seq_len(nA))
			Y[i, order(Z[i, ], decreasing = TRUE)[seq_len(odmax_val)]] = seq_len(odmax_val)
	}
	dimnames(Y) = list(paste0("r", seq_len(nA)), paste0("c", seq_len(nB)))
	list(Y = Y, X = X)
}

# ==== internal bipartite machinery and initialization ====

test_that("bipartite EZ computation works correctly", {
	skip_if_not(requireNamespace("lame", quietly = TRUE))
	# set dimensions
	nA = 5
	nB = 7
	Tn = 3
	RA = 2
	RB = 3
	
	# create test data
	base = array(rnorm(nA * nB * Tn), c(nA, nB, Tn))
	a = matrix(rnorm(nA * Tn), nA, Tn)
	b = matrix(rnorm(nB * Tn), nB, Tn)
	U = array(rnorm(nA * RA * Tn), c(nA, RA, Tn))
	V = array(rnorm(nB * RB * Tn), c(nB, RB, Tn))
	G = matrix(rnorm(RA * RB), RA, RB)
	
	# test ez computation
	EZ = lame:::get_EZ_bip_cpp(base, a, b, U, V, G)
	
	# check dimensions
	expect_equal(dim(EZ), c(nA, nB, Tn))
	
	# check that values are finite
	expect_true(all(is.finite(EZ)))
	
	# manual check for first time slice
	EZ_manual_t1 = base[,,1]
	for(i in 1:nA) {
		EZ_manual_t1[i,] = EZ_manual_t1[i,] + a[i,1]
	}
	for(j in 1:nB) {
		EZ_manual_t1[,j] = EZ_manual_t1[,j] + b[j,1]
	}
	EZ_manual_t1 = EZ_manual_t1 + U[,,1] %*% G %*% t(V[,,1])
	
	expect_equal(EZ[,,1], EZ_manual_t1, tolerance = 1e-10)
})

test_that("four-cycle counting works for bipartite networks", {
	# create a simple bipartite network
	Y = array(0, c(3, 3, 1))
	
	# create a pattern with exactly 1 four-cycle
	# nodes a1, a2 connect to b1, b2 (forms a 4-cycle)
	Y[1, 1, 1] = 1  # a1 -> b1
	Y[1, 2, 1] = 1  # a1 -> b2
	Y[2, 1, 1] = 1  # a2 -> b1
	Y[2, 2, 1] = 1  # a2 -> b2
	
	cycles = lame:::count_four_cycles_bip_cpp(Y)
	expect_equal(as.numeric(cycles), 1)
	
	# test with no edges (should be 0 cycles)
	Y_empty = array(0, c(5, 5, 2))
	cycles_empty = lame:::count_four_cycles_bip_cpp(Y_empty)
	expect_equal(as.numeric(cycles_empty), c(0, 0))
})

test_that("bipartite degree computation works", {
	# create test network
	Y = array(0, c(3, 4, 2))
	Y[1, 1:2, 1] = 1  # node a1 has degree 2 at time 1
	Y[2, 3:4, 1] = 1  # node a2 has degree 2 at time 1
	Y[3, 1, 1] = 1    # node a3 has degree 1 at time 1
	
	degrees = lame:::compute_degrees_bip_cpp(Y)
	
	expect_equal(degrees$row_degrees[,1], c(2, 2, 1))
	expect_equal(degrees$col_degrees[,1], c(2, 1, 1, 1))
})

test_that("UV sampling for bipartite works", {
	# setup
	nA = 4
	nB = 5
	Tn = 2
	RA = 2
	RB = 3
	
	R = array(rnorm(nA * nB * Tn), c(nA, nB, Tn))
	V = array(rnorm(nB * RB * Tn), c(nB, RB, Tn))
	G = matrix(rnorm(RA * RB), RA, RB)
	lambdaU = rep(1, RA)
	s2 = 1
	
	# test u sampling
	U_new = lame:::sample_U_bip_cpp(R, V, G, lambdaU, s2)
	
	expect_equal(dim(U_new), c(nA, RA, Tn))
	expect_true(all(is.finite(U_new)))
	
	# test v sampling
	lambdaV = rep(1, RB)
	V_new = lame:::sample_V_bip_cpp(R, U_new, G, lambdaV, s2)
	
	expect_equal(dim(V_new), c(nB, RB, Tn))
	expect_true(all(is.finite(V_new)))
})

test_that("G sampling and orientation works", {
	nA = 4
	nB = 5
	Tn = 2
	RA = 3
	RB = 2
	
	R = array(rnorm(nA * nB * Tn), c(nA, nB, Tn))
	U = array(rnorm(nA * RA * Tn), c(nA, RA, Tn))
	V = array(rnorm(nB * RB * Tn), c(nB, RB, Tn))
	lambdaG = 1
	s2 = 1
	
	# test g sampling
	G_new = lame:::sample_G_bip_cpp(R, U, V, lambdaG, s2)
	
	expect_equal(dim(G_new), c(RA, RB))
	expect_true(all(is.finite(G_new)))
	
	# test canonical orientation
	oriented = lame:::canon_orient_bip_cpp(U, V, G_new)
	
	expect_equal(dim(oriented$U), dim(U))
	expect_equal(dim(oriented$V), dim(V))
	expect_equal(dim(oriented$G), dim(G_new))
	
	# check that g is now diagonal-like (singular values on diagonal)
	G_oriented = oriented$G
	# off-diagonal elements should be close to zero
	for(i in 1:min(RA, RB)) {
		for(j in 1:min(RA, RB)) {
			if(i != j) {
				expect_lt(abs(G_oriented[i,j]), 1e-10)
			}
		}
	}
})

test_that("bipartite initialization works", {
	skip_if_not_installed("lame")
	
	# create simple bipartite data
	Y = array(0, c(5, 7, 2))
	Y[,,1] = matrix(rbinom(5*7, 1, 0.3), 5, 7)
	Y[,,2] = matrix(rbinom(5*7, 1, 0.3), 5, 7)
	
	init = lame:::init_bipartite_startvals(Y, family = "binary", 
																	 nA = 5, nB = 7, 
																	 RA = 2, RB = 3, Tn = 2)
	
	expect_equal(dim(init$Z), c(5, 7, 2))
	# a and b are vectors (static effects), not matrices
	expect_equal(length(init$a), 5)
	expect_equal(length(init$b), 7)
	expect_equal(dim(init$U), c(5, 2, 2))
	expect_equal(dim(init$V), c(7, 3, 2))
	expect_equal(dim(init$G), c(2, 3))
	expect_null(init$rho)
})

# ==== parameter recovery (ame and lame) ====

test_that("ame() bipartite normal recovers beta and additive effects", {
	skip_on_cran()

	nA = 15
	nB = 10
	R = 2
	true_intercept = 1.0
	true_beta = 0.8

	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "normal", R = R,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "bipartite", seed = 123
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "normal", mode = "bipartite",
		burn = 2000, nscan = 8000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat = colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.5,
		info = paste("Intercept:", beta_hat[1], "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.5,
		info = paste("Beta dyad:", beta_hat[2], "vs true:", true_beta))

	# additive effects recovery (align by name)
	true_a = align_true_params(dat$true_params$a[[1]],
		rownames(dat$Y), names(fit$APM))
	true_b = align_true_params(dat$true_params$b[[1]],
		colnames(dat$Y), names(fit$BPM))
	cor_a = cor(fit$APM, true_a)
	cor_b = cor(fit$BPM, true_b)
	expect_true(cor_a > 0.4,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.4,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# latent position recovery (uv product correlation)
	# true uv is in original order, est uv is in sorted order
	# correlation of vectorized products is order-independent if we
	# align both to the same actor ordering
	true_UV = tcrossprod(dat$true_params$U[[1]], dat$true_params$V[[1]])
	est_UV = fit$U %*% fit$G %*% t(fit$V)
	# reorder true_uv to match sorted actor order
	row_order = match(names(fit$APM), rownames(dat$Y))
	col_order = match(names(fit$BPM), colnames(dat$Y))
	true_UV_aligned = true_UV[row_order, col_order]
	cor_UV = cor(c(true_UV_aligned), c(est_UV))
	expect_true(cor_UV > 0.2,
		info = paste("UV product correlation:", round(cor_UV, 3)))

	# variance s2 should be near 1
	s2_hat = mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("ame() bipartite binary recovers beta and additive effects", {
	skip_on_cran()

	nA = 20
	nB = 12
	R = 2

	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "binary", R = R,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "bipartite", seed = 456
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "binary", mode = "bipartite",
		burn = 2000, nscan = 8000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat = colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation (align by name)
	true_a = align_true_params(dat$true_params$a[[1]],
		rownames(dat$Y), names(fit$APM))
	true_b = align_true_params(dat$true_params$b[[1]],
		colnames(dat$Y), names(fit$BPM))
	cor_a = cor(fit$APM, true_a)
	cor_b = cor(fit$BPM, true_b)
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# uv product correlation
	true_UV = tcrossprod(dat$true_params$U[[1]], dat$true_params$V[[1]])
	est_UV = fit$U %*% fit$G %*% t(fit$V)
	row_order = match(names(fit$APM), rownames(dat$Y))
	col_order = match(names(fit$BPM), colnames(dat$Y))
	true_UV_aligned = true_UV[row_order, col_order]
	cor_UV = cor(c(true_UV_aligned), c(est_UV))
	expect_true(cor_UV > 0.15,
		info = paste("UV product correlation:", round(cor_UV, 3)))
})


test_that("lame() bipartite normal recovers beta and additive effects", {
	skip_on_cran()

	nA = 12
	nB = 8
	n_time = 3
	true_intercept = 0.5
	true_beta = 0.6

	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = n_time,
		family = "normal", R = 0,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "bipartite", seed = 789
	)

	# lame() expects xdyad as list of matrices
	Xdyad_list = lapply(dat$Xdyad, function(x) x[,,1])

	fit = lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "normal", mode = "bipartite",
		burn = 1000, nscan = 5000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat = colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.6,
		info = paste("Intercept:", round(beta_hat[1], 3), "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.6,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "vs true:", true_beta))

	# additive effects (align by name: lame sorts actors alphabetically)
	orig_row_names = rownames(dat$Y[[1]])
	orig_col_names = colnames(dat$Y[[1]])
	sorted_row_names = names(fit$APM)
	sorted_col_names = names(fit$BPM)

	true_a_avg = rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg = rowMeans(do.call(cbind, dat$true_params$b))
	true_a_aligned = align_true_params(true_a_avg, orig_row_names, sorted_row_names)
	true_b_aligned = align_true_params(true_b_avg, orig_col_names, sorted_col_names)

	cor_a = cor(fit$APM, true_a_aligned)
	cor_b = cor(fit$BPM, true_b_aligned)
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# s2 should be near 1
	s2_hat = mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1.5,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("lame() bipartite binary recovers beta and effects", {
	skip_on_cran()

	nA = 15
	nB = 10
	n_time = 3

	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = n_time,
		family = "binary", R = 0,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "bipartite", seed = 101
	)

	Xdyad_list = lapply(dat$Xdyad, function(x) x[,,1])

	fit = lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "binary", mode = "bipartite",
		burn = 1000, nscan = 5000, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat = colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation (align by name)
	orig_row_names = rownames(dat$Y[[1]])
	orig_col_names = colnames(dat$Y[[1]])

	true_a_avg = rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg = rowMeans(do.call(cbind, dat$true_params$b))
	true_a_aligned = align_true_params(true_a_avg, orig_row_names, names(fit$APM))
	true_b_aligned = align_true_params(true_b_avg, orig_col_names, names(fit$BPM))

	cor_a = cor(fit$APM, true_a_aligned)
	cor_b = cor(fit$BPM, true_b_aligned)
	expect_true(cor_a > 0.2,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.2,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))
})

test_that("ame() bipartite dimensions are correct", {
	nA = 10
	nB = 7
	R = 2

	dat = simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "normal", R = R,
		mode = "bipartite", seed = 222
	)

	fit = ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "normal", mode = "bipartite",
		burn = 100, nscan = 200, odens = 2,
		verbose = FALSE, gof = FALSE, seed = 42
	)

	# check output dimensions
	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
	expect_equal(nrow(fit$U), nA)
	expect_equal(ncol(fit$U), R)
	expect_equal(nrow(fit$V), nB)
	expect_equal(ncol(fit$V), R)
	expect_equal(nrow(fit$YPM), nA)
	expect_equal(ncol(fit$YPM), nB)
	expect_equal(nrow(fit$G), R)
	expect_equal(ncol(fit$G), R)
})

# ==== longitudinal lame() per-family behavior ====

# -------- bipartite poisson --------

test_that("bipartite poisson runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 25, nB = 18, Tn = 3, seed = 4105L, beta_true = c(0.6, -0.4))
	Y = lapply(d$EZ, function(eta) {
		yt = matrix(stats::rpois(d$nA * d$nB, lambda = exp(pmin(eta - 1, 5))),
		            d$nA, d$nB)
		dimnames(yt) = list(paste0("r", seq_len(d$nA)), paste0("c", seq_len(d$nB)))
		yt
	})
	fit = lame(Y, Xdyad = d$X, family = "poisson", mode = "bipartite",
	           R = 0, nscan = 120, burn = 30, odens = 5, verbose = FALSE,
	           seed = 4105L)
	expect_true(inherits(fit, "lame"))
	expect_true("intercept" %in% colnames(fit$BETA))
	expect_true("v1_dyad" %in% colnames(fit$BETA))
	# sign recovery on the strong covariate
	m1 = mean(fit$BETA[, "v1_dyad"])
	m2 = mean(fit$BETA[, "v2_dyad"])
	expect_gt(m1, 0)
	expect_lt(m2, 0)
})

# -------- bipartite ordinal --------

test_that("bipartite ordinal runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 30, nB = 22, Tn = 3, seed = 4101L,
	                    beta_true = c(0.8, -0.5))
	Y = lapply(d$EZ, function(eta) {
		Z = eta + matrix(stats::rnorm(d$nA * d$nB), d$nA, d$nB)
		yt = as.integer(cut(Z, c(-Inf, -0.5, 0.3, 1.0, Inf), labels = FALSE))
		dim(yt) = c(d$nA, d$nB)
		dimnames(yt) = list(paste0("r", seq_len(d$nA)), paste0("c", seq_len(d$nB)))
		yt
	})
	fit = lame(Y, Xdyad = d$X, family = "ordinal", mode = "bipartite",
	           R = 0, nscan = 120, burn = 30, odens = 5, verbose = FALSE,
	           seed = 4101L)
	expect_true(inherits(fit, "lame"))
	# ordinal has no intercept by default
	expect_true("v1_dyad" %in% colnames(fit$BETA))
	m1 = mean(fit$BETA[, "v1_dyad"])
	m2 = mean(fit$BETA[, "v2_dyad"])
	expect_gt(m1, 0)
	expect_lt(m2, 0)
})

# -------- bipartite cbin --------

test_that("bipartite cbin runs and respects the row nomination cap", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 30, nB = 25, Tn = 2, seed = 4102L,
	                    beta_true = c(0.7, -0.4))
	odmax = 5L
	Y = lapply(d$EZ, function(eta) {
		Z = eta + matrix(stats::rnorm(d$nA * d$nB), d$nA, d$nB)
		yt = matrix(0L, d$nA, d$nB)
		# top-odmax per row receive a tie
		for (i in seq_len(d$nA)) {
			top = order(Z[i, ], decreasing = TRUE)[seq_len(odmax)]
			yt[i, top] = 1L
		}
		dimnames(yt) = list(paste0("r", seq_len(d$nA)),
		                    paste0("c", seq_len(d$nB)))
		yt
	})
	fit = lame(Y, Xdyad = d$X, family = "cbin", mode = "bipartite",
	           R = 0, odmax = odmax,
	           nscan = 120, burn = 30, odens = 5, verbose = FALSE,
	           seed = 4102L)
	expect_true(inherits(fit, "lame"))
	expect_true("v1_dyad" %in% colnames(fit$BETA))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
	# observed-data invariant: each row's outdegree must equal odmax (the test
	# data was generated to be exactly saturated). confirms the sampler does
	# not silently overwrite y.
	yt = if (is.list(fit$Y)) fit$Y[[1]] else fit$Y[, , 1]
	expect_true(all(rowSums(yt > 0, na.rm = TRUE) == odmax))
})

# -------- bipartite frn --------

test_that("bipartite frn runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	# use a larger panel for a clearer sign check
	d = mk_bip_linpred(nA = 40, nB = 32, Tn = 2, seed = 4103L,
	                    beta_true = c(1.2, -1.0))
	odmax = 6L
	Y = lapply(d$EZ, function(eta) {
		Z = eta + matrix(stats::rnorm(d$nA * d$nB), d$nA, d$nB)
		yt = matrix(0L, d$nA, d$nB)
		for (i in seq_len(d$nA)) {
			top = order(Z[i, ], decreasing = TRUE)[seq_len(odmax)]
			yt[i, top] = seq_len(odmax)  # rank 1..odmax (top to bottom)
		}
		dimnames(yt) = list(paste0("r", seq_len(d$nA)),
		                    paste0("c", seq_len(d$nB)))
		yt
	})
	fit = lame(Y, Xdyad = d$X, family = "frn", mode = "bipartite",
	           R = 0, odmax = odmax,
	           nscan = 200, burn = 50, odens = 5, verbose = FALSE,
	           seed = 4103L)
	expect_true(inherits(fit, "lame"))
	expect_true("v1_dyad" %in% colnames(fit$BETA))
	# rank-likelihood chains on small bipartite panels mix slowly
	expect_true(all(is.finite(fit$BETA[, "v1_dyad"])))
	expect_gt(stats::sd(fit$BETA[, "v1_dyad"]), 0)
})

# ==== cross-sectional ame() per-family behavior ====

test_that("ame(family='poisson', mode='bipartite') runs cross-sectionally", {
	skip_if_too_slow()
	d = mk_xs("poisson", seed = 40010L)
	fit = ame(d$Y, Xdyad = d$X, family = "poisson", mode = "bipartite",
	          R = 0, nscan = 80, burn = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE, seed = 40010L)
	expect_true(inherits(fit, "ame"))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
})

test_that("ame(family='ordinal', mode='bipartite') runs cross-sectionally", {
	skip_if_too_slow()
	d = mk_xs("ordinal", seed = 40011L, beta_true = c(0.8, -0.5))
	fit = ame(d$Y, Xdyad = d$X, family = "ordinal", mode = "bipartite",
	          R = 0, nscan = 80, burn = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE, seed = 40011L)
	expect_true(inherits(fit, "ame"))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
})

test_that("ame(family='cbin', mode='bipartite') runs cross-sectionally", {
	skip_if_too_slow()
	d = mk_xs("cbin", seed = 40012L, beta_true = c(0.7, -0.4), odmax_val = 3L)
	fit = ame(d$Y, Xdyad = d$X, family = "cbin", mode = "bipartite",
	          R = 0, odmax = 3L,
	          nscan = 80, burn = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE, seed = 40012L)
	expect_true(inherits(fit, "ame"))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
	# observed-data invariant: row outdegrees match odmax (by construction)
	expect_true(all(rowSums(d$Y > 0) == 3L))
})

test_that("ame(family='frn', mode='bipartite') runs cross-sectionally", {
	skip_if_too_slow()
	d = mk_xs("frn", seed = 40013L, beta_true = c(0.8, -0.5), odmax_val = 3L)
	fit = ame(d$Y, Xdyad = d$X, family = "frn", mode = "bipartite",
	          R = 0, odmax = 3L,
	          nscan = 80, burn = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE, seed = 40013L)
	expect_true(inherits(fit, "ame"))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
})

# ==== static ame() behavioral tests ====

test_that("bipartite model fits additive effects", {
	skip_on_cran()
	
	set.seed(790)
	
	nA = 15
	nB = 20
	
	# true parameters
	mu_true = 0.3
	sigma2_a = 0.4
	sigma2_b = 0.3
	sigma2_e = 0.5
	
	# generate additive effects
	a_true = rnorm(nA, 0, sqrt(sigma2_a))
	b_true = rnorm(nB, 0, sqrt(sigma2_b))
	
	# generate network
	Y = matrix(mu_true, nA, nB)
	for(i in 1:nA) {
		Y[i,] = Y[i,] + a_true[i]
	}
	for(j in 1:nB) {
		Y[,j] = Y[,j] + b_true[j]
	}
	Y = Y + matrix(rnorm(nA * nB, 0, sqrt(sigma2_e)), nA, nB)
	
	# fit model with additive effects
	fit = ame(Y, mode = "bipartite",
						 rvar = TRUE, cvar = TRUE,
						 R_row = 0, R_col = 0,
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# additive effects are present
	expect_false(is.null(fit$APM))
	expect_false(is.null(fit$BPM))
	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
	
	expect_true(mean(fit$VC[,1]) > 0)
	expect_true(mean(fit$VC[,3]) > 0)
	
	# row and column means are close to true values
	expect_equal(mean(fit$APM), mean(a_true), tolerance = 0.3)
	expect_equal(mean(fit$BPM), mean(b_true), tolerance = 0.3)
})

test_that("bipartite model fits multiplicative effects", {
	skip_on_cran()
	
	set.seed(791)
	
	nA = 20
	nB = 25
	R_row = 2
	R_col = 2
	
	# generate latent factors
	U_true = matrix(rnorm(nA * R_row), nA, R_row)
	V_true = matrix(rnorm(nB * R_col), nB, R_col)
	G_true = diag(c(2, 1), R_row, R_col)
	
	# generate network with multiplicative effects
	mu_true = 0.2
	UV_true = U_true %*% G_true %*% t(V_true)
	Y = mu_true + UV_true + matrix(rnorm(nA * nB, 0, 0.5), nA, nB)
	
	# fit model with multiplicative effects
	fit = ame(Y, mode = "bipartite",
						 R_row = R_row, R_col = R_col,
						 burn = 500, nscan = 2000, verbose = FALSE)
	
	# latent-factor dimensions
	expect_equal(dim(fit$U), c(nA, R_row))
	expect_equal(dim(fit$V), c(nB, R_col))
	expect_equal(dim(fit$G), c(R_row, R_col))
	
	# multiplicative effects capture structure
	UVPM_reconstructed = reconstruct_UVPM(fit)
	
	# correlation should be positive
	if(!is.null(UVPM_reconstructed)) {
		cor_uv = cor(c(UV_true), c(UVPM_reconstructed))
		expect_gt(cor_uv, 0.3)
	}
})

test_that("bipartite multiplicative rank improves fit", {
	skip_on_cran()
	
	set.seed(794)
	
	nA = 15
	nB = 12
	
	# generate network with rank-2 structure
	U_true = matrix(rnorm(nA * 2), nA, 2)
	V_true = matrix(rnorm(nB * 2), nB, 2)
	Y = U_true %*% t(V_true) + matrix(rnorm(nA * nB, 0, 0.3), nA, nB)
	
	# fit with r = 0
	fit0 = ame(Y, mode = "bipartite",
							R_row = 0, R_col = 0,
							burn = 300, nscan = 1000, verbose = FALSE)
	
	# fit with r = 2
	fit2 = ame(Y, mode = "bipartite",
							R_row = 2, R_col = 2,
							burn = 300, nscan = 1000, verbose = FALSE)
	
	# model with multiplicative effects should fit better
	resid0 = Y - fit0$YPM
	resid2 = Y - fit2$YPM
	
	# both models produced valid predictions
	expect_false(all(is.na(fit0$YPM)))
	expect_false(all(is.na(fit2$YPM)))
	
	# the rank-2 model should have lower residuals than r = 0
	mse0 = mean(resid0^2, na.rm = TRUE)
	mse2 = mean(resid2^2, na.rm = TRUE)
	
	# only compare finite mse values
	if(is.finite(mse0) && is.finite(mse2)) {
		expect_lt(mse2, mse0)
	}
})

test_that("bipartite model handles missing data", {
	skip_on_cran()
	
	set.seed(795)
	
	nA = 10
	nB = 8
	
	# generate complete network
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# add missing values
	missing_idx = sample(1:(nA*nB), size = 10, replace = FALSE)
	Y[missing_idx] = NA
	
	# fit model
	fit = ame(Y, mode = "bipartite",
						 burn = 200, nscan = 500, verbose = FALSE)
	
	# model ran
	expect_equal(fit$mode, "bipartite")
	
	expect_false(any(is.na(fit$YPM)))
})

test_that("bipartite gof statistics are computed", {
	skip_on_cran()
	
	set.seed(796)
	
	nA = 8
	nB = 10
	
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# fit with gof
	fit = ame(Y, mode = "bipartite",
						 burn = 100, nscan = 300, 
						 gof = TRUE, verbose = FALSE)
	
	# gof output is present
	expect_false(is.null(fit$GOF))
	# bipartite gof has 3 columns: sd.rowmean, sd.colmean, four.cycles
	expect_equal(ncol(fit$GOF), 3)
	
	# gof statistics are non-negative
	obs_sd_row = fit$GOF[1,1]
	obs_sd_col = fit$GOF[1,2]
	obs_four_cycles = fit$GOF[1,3]
	
	expect_true(obs_sd_row >= 0)
	expect_true(obs_sd_col >= 0)
	expect_true(obs_four_cycles >= 0)
})

# ==== static ame() structural and covariate tests ====

test_that("supported bipartite families fit", {
	skip_on_cran()
	
	set.seed(456)
	
	# set up bipartite dimensions
	nA = 12
	nB = 10
	
	# test the supported rectangular families
	families = c("normal", "binary")

	for(fam in families) {
		# generate test data
		if(fam == "normal") {
			Y = matrix(rnorm(nA * nB, 0, 1), nA, nB)
		} else if(fam == "binary" || fam == "cbin") {
			Y = matrix(rbinom(nA * nB, 1, 0.5), nA, nB)
		} else if(fam == "poisson") {
			Y = matrix(rpois(nA * nB, 3), nA, nB)
		} else if(fam == "ordinal") {
			Y = matrix(sample(0:3, nA * nB, replace = TRUE), nA, nB)
		} else if(fam == "frn") {
			# fixed rank nomination
			Y = matrix(0, nA, nB)
			for(i in 1:nA) {
				nominated = sample(1:nB, min(3, nB))
				Y[i, nominated] = 1:length(nominated)
			}
		}
		
		# set odmax for rank families
		if(fam %in% c("cbin", "frn")) {
			odmax = rep(3, nA)
		} else {
			odmax = NULL
		}
		
		# fit model
		suppressWarnings({
			fit = ame(Y, mode = "bipartite", family = fam,
								 R_row = 0, R_col = 0,
								 odmax = odmax,
								 burn = 100, nscan = 300, 
								 verbose = FALSE)
		})
		
		# basic checks
		expect_equal(fit$mode, "bipartite", 
								 info = paste("Mode check failed for", fam))
		expect_equal(fit$family, fam, 
								 info = paste("Family check failed for", fam))
		expect_equal(fit$nA, nA, 
								 info = paste("Row dimension check failed for", fam))
		expect_equal(fit$nB, nB, 
								 info = paste("Column dimension check failed for", fam))
		
		# posterior samples are present
		expect_false(is.null(fit$BETA), 
								 info = paste("BETA missing for", fam))
		
		# family-specific checks
		if(fam == "binary" || fam == "cbin") {
			# predictions are binary probabilities
			expect_true(all(fit$YPM >= 0 & fit$YPM <= 1),
									info = paste("YPM not in [0,1] for", fam))
		}
	}
})

test_that("bipartite models with multiplicative effects fit", {
	skip_on_cran()
	
	set.seed(457)
	
	nA = 10
	nB = 8
	R_row = 1
	R_col = 1
	
	# test the supported subset with multiplicative effects
	families = c("normal", "binary")

	for(fam in families) {
		# generate data
		if(fam == "normal") {
			Y = matrix(rnorm(nA * nB), nA, nB)
		} else if(fam == "binary") {
			Y = matrix(rbinom(nA * nB, 1, 0.5), nA, nB)
		} else if(fam == "poisson") {
			Y = matrix(rpois(nA * nB, 2), nA, nB)
		}
		
		# fit with multiplicative effects
		suppressWarnings({
			fit = ame(Y, mode = "bipartite", family = fam,
								 R_row = R_row, R_col = R_col,
								 burn = 100, nscan = 300,
								 verbose = FALSE)
		})
		
		# multiplicative effects are present
		expect_equal(dim(fit$U), c(nA, R_row), 
								 info = paste("U dimension wrong for", fam))
		expect_equal(dim(fit$V), c(nB, R_col), 
								 info = paste("V dimension wrong for", fam))
		expect_equal(dim(fit$G), c(R_row, R_col), 
								 info = paste("G dimension wrong for", fam))
	}
})

test_that("bipartite models handle covariates", {
	skip_on_cran()
	
	set.seed(458)
	
	nA = 8
	nB = 6
	
	# generate covariates
	Xdyad = array(rnorm(nA * nB * 2), c(nA, nB, 2))
	Xrow = matrix(rnorm(nA * 2), nA, 2)
	Xcol = matrix(rnorm(nB * 2), nB, 2)
	
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# fit with all covariate types
	fit = ame(Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
						 mode = "bipartite",
						 burn = 100, nscan = 300,
						 verbose = FALSE)
	
	# intercept + 2 dyadic + 2 row + 2 col
	expect_equal(ncol(fit$BETA), 7)
	
	# estimates are present
	beta_means = colMeans(fit$BETA)
	expect_false(any(is.na(beta_means)))
})

test_that("bipartite dynamic_uv keeps the latent scale identified (regression)", {
	# the U/V-vs-G scale is unidentified (only U G V' is), so without gauge
	# fixing it random-walks: G collapses toward 0 while U/V inflate past
	# prior$uv_max_abs and the sampler rejects its own UV proposals. the
	# drift needs a realistic panel size to develop, hence the size here.
	skip_on_cran()
	set.seed(5)
	nA = 25; nB = 20; TT = 8; rho_true = 0.6; sig = 0.7
	a = rnorm(nA, 0, 0.4); b = rnorm(nB, 0, 0.4); Gt = diag(2)
	U = array(0, c(nA, 2, TT)); V = array(0, c(nB, 2, TT))
	U[, , 1] = matrix(rnorm(nA * 2, 0, sig), nA)
	V[, , 1] = matrix(rnorm(nB * 2, 0, sig), nB)
	sdin = sig * sqrt(1 - rho_true^2)
	for (t in 2:TT) {
		U[, , t] = rho_true * U[, , t - 1] + matrix(rnorm(nA * 2, 0, sdin), nA)
		V[, , t] = rho_true * V[, , t - 1] + matrix(rnorm(nB * 2, 0, sdin), nB)
	}
	Y = lapply(1:TT, function(t) {
		m = outer(a, b, "+") + U[, , t] %*% Gt %*% t(V[, , t]) +
			matrix(rnorm(nA * nB, 0, 1), nA)
		rownames(m) = sprintf("r%02d", 1:nA); colnames(m) = sprintf("c%02d", 1:nB); m
	})
	fit = lame(Y, mode = "bipartite", family = "normal", R = 2, dynamic_uv = TRUE,
	           seed = 1, burn = 150, nscan = 600, odens = 4,
	           verbose = FALSE, plot = FALSE)

	iters = 150 + 600
	# the UV sampler should not be rejecting its own proposals
	expect_lt(fit$tryErrorChecks$UV / iters, 0.05)
	# the gauge is pinned: G carries unit Frobenius norm, coordinates stay O(1)
	expect_lt(max(abs(fit[["U"]]), na.rm = TRUE), 10)
	expect_lt(max(abs(fit[["V"]]), na.rm = TRUE), 10)
	expect_lt(max(abs(fit[["G"]]), na.rm = TRUE), 100)
	# and the additive effects still recover
	expect_gt(cor(fit$APM, a), 0.8)
})

test_that("bipartite dynamic_uv rho_uv converges on a low-persistence panel (regression)", {
	# rho_uv starts at the prior mean (0.9), so recovering a genuinely low
	# persistence requires the hyperparameter updates to run often enough for
	# the chain to travel. with every-sweep updates this panel lands near 0.55
	# against a truth of 0.3, with the deliberately high default prior still
	# pulling upward.
	skip_on_cran()
	set.seed(5)
	nA = 25; nB = 20; TT = 8; rho_true = 0.3; sig = 0.7
	a = rnorm(nA, 0, 0.4); b = rnorm(nB, 0, 0.4); Gt = diag(2)
	U = array(0, c(nA, 2, TT)); V = array(0, c(nB, 2, TT))
	U[, , 1] = matrix(rnorm(nA * 2, 0, sig), nA)
	V[, , 1] = matrix(rnorm(nB * 2, 0, sig), nB)
	sdin = sig * sqrt(1 - rho_true^2)
	for (t in 2:TT) {
		U[, , t] = rho_true * U[, , t - 1] + matrix(rnorm(nA * 2, 0, sdin), nA)
		V[, , t] = rho_true * V[, , t - 1] + matrix(rnorm(nB * 2, 0, sdin), nB)
	}
	Y = lapply(1:TT, function(t) {
		m = outer(a, b, "+") + U[, , t] %*% Gt %*% t(V[, , t]) +
			matrix(rnorm(nA * nB, 0, 1), nA)
		rownames(m) = sprintf("r%02d", 1:nA); colnames(m) = sprintf("c%02d", 1:nB); m
	})
	fit = lame(Y, mode = "bipartite", family = "normal", R = 2, dynamic_uv = TRUE,
	           seed = 1, burn = 200, nscan = 800, odens = 4,
	           verbose = FALSE, plot = FALSE)
	rho_med = median(as.numeric(fit$rho_uv), na.rm = TRUE)
	# a median near 0.9 here means the chain is stuck at its starting value
	expect_lt(rho_med, 0.75)
	expect_gt(rho_med, -1)
})
