skip_on_cran()

# one cheap parameter-recovery check per supported family on the
# cross-sectional ame() path. family-specific edge cases live in the
# per-family test_static_*.r files.

library(lame)
library(testthat)

# generate a synthetic ame network with the requested family. eta is mu +
# beta * x plus standard normal noise; family map applies the appropriate
# link or thresholding without collapsing the matrix shape.
make_net = function(family, n = 20, mu = 0.5, beta = 0.6, seed = 6886) {
	set.seed(seed)
	# iid symmetric dyadic covariate. a rank-1 tcrossprod(u) covariate is
	# partially confounded with the additive row/col effects (u_i u_j
	# correlates with a_i + b_j), which attenuates beta by ~13% on any
	# single realization and makes the coverage assertion below a test of
	# the confound, not of coefficient recovery.
	X = matrix(rnorm(n * n), n, n)
	X = (X + t(X)) / 2
	diag(X) = NA
	eta = mu + beta * X
	noise = matrix(rnorm(n * n, 0, 1), n, n)
	rn = paste0("a", sprintf("%02d", seq_len(n)))

	Y = switch(family,
		normal  = eta + noise,
		binary  = (eta + noise > 0) * 1,
		poisson = matrix(rpois(n * n, pmin(exp(eta), 1e6)), n, n),
		ordinal = matrix(pmin(pmax(round(eta + noise), 0), 3), n, n),
		cbin    = (eta + noise > 0) * 1,
		stop("family not handled: ", family)
	)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = rn
	list(Y = Y, X = X, mu = mu, beta = beta)
}

# 95% ci on the dyadic coefficient should cover the truth
recovery_test = function(
	family, n = 20, mu = 0.5, beta = 0.6,
	burn = 200, nscan = 600, R = 0
) {
	d = make_net(family, n = n, mu = mu, beta = beta)
	fit = ame(
		Y = d$Y, Xdyad = d$X, R = R, family = family,
		burn = burn, nscan = nscan, verbose = FALSE, plot = FALSE
	)
	# locate the dyadic coefficient column by name
	beta_name = grep("\\.dyad$", colnames(fit$BETA), value = TRUE)[1]
	if (is.na(beta_name)) {
		beta_post = fit$BETA[, ncol(fit$BETA)]
	} else {
		beta_post = fit$BETA[, beta_name]
	}
	ci = quantile(beta_post, c(0.025, 0.975), na.rm = TRUE)
	expect_true(
		ci[1] <= d$beta && d$beta <= ci[2],
		info = sprintf(
			"%s: beta=%.2f outside [%.2f, %.2f]",
			family, d$beta, ci[1], ci[2]
		)
	)
}

test_that("normal AME recovers the dyadic coefficient", {
	recovery_test("normal")
})

test_that("binary AME recovers the dyadic coefficient", {
	recovery_test("binary")
})

test_that("poisson AME recovers the dyadic coefficient", {
	recovery_test("poisson", beta = 0.4)
})

test_that("ordinal AME recovers the dyadic coefficient", {
	recovery_test("ordinal")
})

# bipartite headline check: the rectangular ame() path with a dyadic covariate
test_that("normal bipartite AME runs and produces a fit object", {
	set.seed(6886)
	nr = 12; nc = 10
	Y = matrix(rnorm(nr * nc), nr, nc)
	rownames(Y) = paste0("r", seq_len(nr))
	colnames(Y) = paste0("c", seq_len(nc))
	X = array(rnorm(nr * nc), c(nr, nc, 1))
	dimnames(X)[[3]] = "x"
	fit = ame(Y, Xdyad = X, mode = "bipartite", R = 0,
	          family = "normal", burn = 50, nscan = 200,
	          verbose = FALSE, plot = FALSE)
	expect_s3_class(fit, "ame")
	cf = coef(fit)
	dyad_name = grep("dyad$", names(cf), value = TRUE)[1]
	expect_true(!is.na(dyad_name) && is.finite(cf[dyad_name]))
})

test_that("binary bipartite AME runs", {
	set.seed(6886)
	nr = 12; nc = 10
	Y = matrix(rbinom(nr * nc, 1, 0.3), nr, nc)
	rownames(Y) = paste0("r", seq_len(nr))
	colnames(Y) = paste0("c", seq_len(nc))
	fit = ame(Y, mode = "bipartite", R = 0, family = "binary",
	          burn = 50, nscan = 200, verbose = FALSE, plot = FALSE)
	expect_s3_class(fit, "ame")
})
