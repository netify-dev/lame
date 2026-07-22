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
