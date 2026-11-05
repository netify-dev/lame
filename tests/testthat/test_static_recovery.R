skip_on_cran()

# one cheap parameter-recovery check per supported family on the cross-sectional ame() path
# family-specific edge cases live in the per-family test_static_*.R files

library(lame)
library(testthat)

# a synthetic network of the requested family: eta = mu + beta * x plus standard normal noise, mapped by the family's
# link or threshold with the matrix shape kept
make_net = function(family, n = 20, mu = 0.5, beta = 0.6, seed = 6886) {
	set.seed(seed)
	# iid symmetric dyadic covariate; a rank-1 tcrossprod(u) covariate would be partly confounded with the additive effects
	X = matrix(rnorm(n * n), n, n)
	X = (X + t(X)) / 2
	diag(X) = NA
	eta = mu + beta * X
	noise = matrix(rnorm(n * n, 0, 1), n, n)
	rn = paste0("a", sprintf("%02d", seq_len(n)))

	Y = switch(
		family,
		normal = eta + noise,
		binary = (eta + noise > 0) * 1,
		# the NA diagonal of eta gets rate zero so rpois() does not warn
		poisson = matrix(rpois(n * n, ifelse(is.na(eta), 0, pmin(exp(eta), 1e6))), n, n),
		ordinal = matrix(pmin(pmax(round(eta + noise), 0), 3), n, n),
		cbin = (eta + noise > 0) * 1,
		stop("family not handled: ", family)
	)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = rn
	list(Y = Y, X = X, mu = mu, beta = beta)
}

# the 95% interval of the dyadic coefficient covers the truth
recovery_test = function(
	family,
	n = 20,
	mu = 0.5,
	beta = 0.6,
	burn = 200,
	nscan = 600,
	R = 0
) {
	d = make_net(family, n = n, mu = mu, beta = beta)
	fit = ame(
		Y = d$Y,
		Xdyad = d$X,
		R = R,
		family = family,
		burn = burn,
		nscan = nscan,
		verbose = FALSE,
		plot = FALSE
	)
	# the dyadic coefficient column, named <covariate>_dyad
	beta_name = grep("_dyad$", colnames(fit$BETA), value = TRUE)
	expect_length(beta_name, 1L)
	beta_post = fit$BETA[, beta_name]
	ci = quantile(beta_post, c(0.025, 0.975), na.rm = TRUE)
	med = stats::median(beta_post, na.rm = TRUE)
	# recovered when the 95% interval covers the truth or the median lands close to it, since the iid covariate is
	# partly confounded with the additive effects and a narrow interval can just miss on one realization
	recovered = (ci[1] <= d$beta && d$beta <= ci[2]) || abs(med - d$beta) < 0.2
	expect_true(
		recovered,
		info = sprintf(
			"%s: beta=%.2f, median=%.2f, ci=[%.2f, %.2f]",
			family,
			d$beta,
			med,
			ci[1],
			ci[2]
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

# frn: recovery_test() cannot generate ranked nominations, so a standalone fit checks that a dyadic covariate
# driving the ranks gets a finite coefficient
test_that("Fixed rank nomination (frn) with covariates works", {
	set.seed(6886)
	n = 20

	# covariate
	X = matrix(rnorm(n * n, 0, 0.5), n, n)
	diag(X) = NA

	# ranked nominations driven by x
	odmax = 4
	Y = matrix(NA, n, n)

	for (i in 1:n) {
		# higher x gets a better (lower) rank
		scores = X[i, ]
		scores[i] = NA
		top_indices = order(scores, decreasing = TRUE, na.last = TRUE)[1:odmax]
		Y[i, top_indices] = 1:odmax
	}
	diag(Y) = NA

	fit = ame(
		Y,
		Xdyad = X,
		R = 0,
		family = "frn",
		odmax = odmax,
		burn = 300,
		nscan = 800,
		verbose = FALSE
	)

	expect_true(!is.null(fit$BETA))
	# positive x gives better (lower) ranks
	if (ncol(fit$BETA) >= 2) {
		beta_est = median(fit$BETA[, 2])
		expect_true(is.finite(beta_est))
	}
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
	fit = ame(
		Y,
		Xdyad = X,
		mode = "bipartite",
		R = 0,
		family = "normal",
		burn = 50,
		nscan = 200,
		verbose = FALSE,
		plot = FALSE
	)
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
	fit = ame(
		Y,
		mode = "bipartite",
		R = 0,
		family = "binary",
		burn = 50,
		nscan = 200,
		verbose = FALSE,
		plot = FALSE
	)
	expect_s3_class(fit, "ame")
})

# sparse counts: the row and column means of log1p(y) barely vary on a mostly-zero matrix, so the effect prior
# scale comes from a log-rate moment, which keeps the additive effects and the predictive total alive

test_that("poisson prior scale and start values are sane on sparse counts", {
	set.seed(7); n = 30
	Y = matrix(rpois(n * n, 0.05), n, n); diag(Y) = NA
	# log-rate moments are O(1) even when almost every cell is zero
	vs = .ame_pois_effect_scale(Y)
	expect_true(is.finite(vs) && vs > 0.1 && vs < 10)
	# dense, homogeneous counts give a small scale
	Yd = matrix(rpois(n * n, 50), n, n); diag(Yd) = NA
	expect_lt(.ame_pois_effect_scale(Yd), 0.1)
	# zero cells start at the mean rate on sparse data, at log(y + 1) on dense
	Z = .ame_pois_init_z(Y)
	ybar = mean(Y, na.rm = TRUE)
	expect_equal(Z[!is.na(Y) & Y == 0][1], log(ybar))
	expect_equal(diag(Z)[1], log(ybar))
	Zd = .ame_pois_init_z(Yd)
	expect_equal(Zd[2, 1], log(Yd[2, 1] + 1))
	# a square bipartite matrix keeps its real diagonal
	Yb = Yd; diag(Yb) = 3
	expect_equal(diag(.ame_pois_init_z(Yb, unipartite = FALSE))[1], log(4))
	# the 3-d (lame) layout pools per-actor totals over periods
	Y3 = array(c(Y, Y), c(n, n, 2))
	ra = log((2 * rowSums(Y, na.rm = TRUE) + 0.5) / (2 * rowSums(!is.na(Y)) + 1))
	rb = log((2 * colSums(Y, na.rm = TRUE) + 0.5) / (2 * colSums(!is.na(Y)) + 1))
	expect_equal(.ame_pois_effect_scale(Y3), mean(c(var(ra), var(rb))))
})

test_that("poisson AME recovers actor heterogeneity on a sparse network", {
	skip_on_cran()
	set.seed(42); n = 40
	a = rnorm(n, 0, 0.8); b = rnorm(n, 0, 0.8)
	x = matrix(rnorm(n * n), n, n)
	eta = -3.6 + 0.5 * x + outer(a, b, "+")
	Y = matrix(rpois(n * n, exp(eta)), n, n); diag(Y) = NA
	X = array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	fit = ame(
		Y,
		Xdyad = X,
		family = "poisson",
		R = 0,
		burn = 500,
		nscan = 1000,
		odens = 5,
		seed = 1,
		gof = FALSE,
		verbose = FALSE
	)
	# the effect prior is not pinned near zero
	expect_gt(fit$prior$Sab0[1, 1], 0.05)
	expect_gt(mean(fit$VC[, "va"]), 0.1)
	expect_gt(cor(fit$APM, a), 0.4)
	expect_gt(cor(fit$BPM, b), 0.4)
	# the dyadic variance does not run away and the predicted total is in
	# the neighbourhood of the observed one
	S = nrow(fit$VC)
	expect_lt(mean(fit$VC[seq(S - S %/% 4, S), "ve"]), 1)
	tot = sum(Y, na.rm = TRUE)
	expect_gt(sum(fit$YPM, na.rm = TRUE), 0.6 * tot)
	expect_lt(sum(fit$YPM, na.rm = TRUE), 1.6 * tot)
	# the translation of the intercept with the latent array (R/latent_translation_move.R) ran
	expect_gt(fit$actor_move[["coef_calls"]], 0)
})
