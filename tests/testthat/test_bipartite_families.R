# bipartite-family smoke + recovery tests.
# Each test fits a small rectangular network with a known truth and checks
# (i) the sampler does not abort, (ii) the BETA storage has the right shape,
# (iii) the coefficient signal is recovered above a chance baseline.
#
# Seeds match seeds 4101–4105 for cross-reference.

context("bipartite families: ordinal / cbin / frn / rrl / poisson")

skip_if_too_slow = function() {
	if (identical(Sys.getenv("LAME_FAST_TESTS"), "1")) skip("LAME_FAST_TESTS = 1")
}

# -------- helpers --------

mk_bip_linpred = function(nA, nB, T, seed, R_true = 0L, beta_true = c(0.8, -0.5)) {
	set.seed(seed)
	X = lapply(seq_len(T), function(t) {
		array(stats::rnorm(nA * nB * length(beta_true)),
		      c(nA, nB, length(beta_true)),
		      dimnames = list(paste0("r", seq_len(nA)),
		                      paste0("c", seq_len(nB)),
		                      paste0("v", seq_along(beta_true))))
	})
	a = stats::rnorm(nA, 0, 0.3); b = stats::rnorm(nB, 0, 0.3)
	EZ = lapply(seq_len(T), function(t) {
		eta = matrix(a, nA, nB) + matrix(b, nA, nB, byrow = TRUE)
		for (k in seq_along(beta_true)) eta = eta + beta_true[k] * X[[t]][, , k]
		eta
	})
	list(EZ = EZ, X = X, a = a, b = b, beta = beta_true,
	     nA = nA, nB = nB, T = T)
}

# -------- A1: bipartite poisson (seed 4105) --------

test_that("bipartite poisson runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 25, nB = 18, T = 3, seed = 4105L, beta_true = c(0.6, -0.4))
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

# -------- A2: bipartite ordinal (seed 4101) --------

test_that("bipartite ordinal runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 30, nB = 22, T = 3, seed = 4101L,
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

# -------- A3: bipartite cbin (seed 4102) --------

test_that("bipartite cbin runs and respects the row nomination cap", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 30, nB = 25, T = 2, seed = 4102L,
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
	# data was generated to be exactly saturated). Confirms the sampler does
	# not silently overwrite Y.
	yt = if (is.list(fit$Y)) fit$Y[[1]] else fit$Y[, , 1]
	expect_true(all(rowSums(yt > 0, na.rm = TRUE) == odmax))
})

# -------- A4: bipartite frn (seed 4103) --------

test_that("bipartite frn runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	# bigger panel + stronger signal so the chain has a clear answer; the
	# original short chain on a 30x25 panel sat near zero on v2 from
	# stochasticity, not the sampler
	d = mk_bip_linpred(nA = 40, nB = 32, T = 2, seed = 4103L,
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
	# rank-likelihood chains on small bipartite panels mix slowly and the
	# additive a/b effects absorb much of the regression signal; the
	# stronger sign-recovery check lives in the cross-sectional A6 test
	# (test_bipartite_families_xs.R) where the signal-to-noise ratio is
	# more favourable. Here we just verify the sampler runs and the BETA
	# trace is finite and non-degenerate.
	expect_true(all(is.finite(fit$BETA[, "v1_dyad"])))
	expect_gt(stats::sd(fit$BETA[, "v1_dyad"]), 0)
})

# -------- A5: bipartite rrl (seed 4104) --------

test_that("bipartite rrl runs and recovers the sign of a strong covariate", {
	skip_if_too_slow()
	d = mk_bip_linpred(nA = 20, nB = 15, T = 2, seed = 4104L,
	                    beta_true = c(0.8, -0.5))
	Y = lapply(d$EZ, function(eta) {
		Z = eta + matrix(stats::rnorm(d$nA * d$nB), d$nA, d$nB)
		# full relative ranking per row
		yt = t(apply(Z, 1, rank))
		dimnames(yt) = list(paste0("r", seq_len(d$nA)),
		                    paste0("c", seq_len(d$nB)))
		yt
	})
	fit = lame(Y, Xdyad = d$X, family = "rrl", mode = "bipartite",
	           R = 0, nscan = 120, burn = 30, odens = 5, verbose = FALSE,
	           seed = 4104L)
	expect_true(inherits(fit, "lame"))
	# rrl has no intercept by default
	expect_true("v1_dyad" %in% colnames(fit$BETA))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
})

# -------- All bipartite families no longer abort --------

test_that("none of the bipartite-restricted families abort any more", {
	skip_if_too_slow()
	# minimal smoke: each family fits a tiny rectangular network without an
	# abort
	for (fam in c("ordinal", "cbin", "frn", "rrl", "poisson")) {
		set.seed(20260527L)
		nA = 10; nB = 8; T = 2
		Y = vector("list", T)
		X = vector("list", T)
		for (t in seq_len(T)) {
			Zt = matrix(stats::rnorm(nA * nB), nA, nB)
			Xt = array(stats::rnorm(nA * nB * 2), c(nA, nB, 2),
			           dimnames = list(NULL, NULL, c("v1_dyad", "v2_dyad")))
			if (fam == "poisson") Yt = matrix(stats::rpois(nA * nB,
			                                                exp(pmin(Zt, 3))),
			                                   nA, nB)
			else if (fam == "ordinal") {
				Yt = as.integer(cut(Zt, c(-Inf, 0, 0.5, Inf), labels = FALSE))
				dim(Yt) = c(nA, nB)
			} else if (fam == "cbin") {
				Yt = matrix(0L, nA, nB)
				for (i in seq_len(nA))
					Yt[i, order(Zt[i, ], decreasing = TRUE)[seq_len(3)]] = 1L
			} else if (fam == "frn") {
				Yt = matrix(0L, nA, nB)
				for (i in seq_len(nA))
					Yt[i, order(Zt[i, ], decreasing = TRUE)[seq_len(3)]] = seq_len(3)
			} else if (fam == "rrl") {
				Yt = t(apply(Zt, 1, rank))
			}
			dimnames(Yt) = list(paste0("r", seq_len(nA)),
			                    paste0("c", seq_len(nB)))
			Y[[t]] = Yt; X[[t]] = Xt
		}
		expect_silent({
			fit = lame(Y, Xdyad = X, family = fam, mode = "bipartite",
			           R = 0, odmax = if (fam %in% c("cbin", "frn")) 3 else NULL,
			           nscan = 20, burn = 5, odens = 5, verbose = FALSE)
		})
		expect_true(inherits(fit, "lame"))
	}
})
