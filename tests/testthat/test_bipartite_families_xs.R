# cross-sectional bipartite families.
# mirrors test_bipartite_families.r for `lame()` but exercises the
# `ame()` (cross-sectional) path through ame_bipartite(). the five
# families that used to abort or warn now route through the same
# rectangular z samplers from r/rz_bipartite.r.

context("cross-sectional bipartite families (ame): poisson / ordinal / cbin / frn / rrl")

skip_if_too_slow = function() {
	if (identical(Sys.getenv("LAME_FAST_TESTS"), "1")) skip("LAME_FAST_TESTS = 1")
}

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
	} else if (family == "rrl") {
		Y = t(apply(Z, 1, rank))
	}
	dimnames(Y) = list(paste0("r", seq_len(nA)), paste0("c", seq_len(nB)))
	list(Y = Y, X = X)
}

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

test_that("ame(family='rrl', mode='bipartite') runs cross-sectionally", {
	skip_if_too_slow()
	d = mk_xs("rrl", seed = 40014L, beta_true = c(0.8, -0.5))
	fit = ame(d$Y, Xdyad = d$X, family = "rrl", mode = "bipartite",
	          R = 0, nscan = 80, burn = 20, odens = 5,
	          verbose = FALSE, plot = FALSE, gof = FALSE, seed = 40014L)
	expect_true(inherits(fit, "ame"))
	expect_gt(mean(fit$BETA[, "v1_dyad"]), 0)
	expect_lt(mean(fit$BETA[, "v2_dyad"]), 0)
})
