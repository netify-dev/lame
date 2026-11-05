skip_on_cran()

test_that("residual prior scale is fixed in observed units and validates overrides", {
	scale = lame:::.lame_residual_prior_scale
	Y = matrix(c(-2, 0.3, 0.8, 1.7, NA, 0.2), 2, 3)
	expect_equal(scale(Y, "normal"), var(c(Y), na.rm = TRUE))
	expect_equal(scale(100 * Y, "normal"), 10000 * scale(Y, "normal"))
	expect_equal(scale(Y, "normal", 0.7), 0.7)
	expect_equal(scale(Y, "poisson"), 1)
	expect_equal(scale(matrix(0, 2, 3), "normal"), 1)
	expect_equal(scale(matrix(2, 2, 3), "normal"), 4)
	expect_null(scale(Y, "binary"))
	for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1")) {
		expect_error(scale(Y, "normal", bad), "s20.*finite positive scalar")
	}
})

test_that("variance kernels match a fixed-prior dyadic joint-density oracle", {
	set.seed(441)
	E = array(rnorm(18), c(3, 3, 2))
	rho = c(-0.4, 0.6)
	s20 = 0.7
	# likelihood quadratic recovered from independent bivariate Gaussian densities, with normalizers, each dyad once
	loglik = function(v, rs) sum(vapply(1:2, function(t) {
		S = v * matrix(c(1, rs[t], rs[t], 1), 2)
		z = E[, , t]; i = which(upper.tri(z), arr.ind = TRUE)
		sum(vapply(seq_len(nrow(i)), function(k) {
			x = c(z[i[k, 1], i[k, 2]], z[i[k, 2], i[k, 1]])
			-log(2 * pi) - log(det(S)) / 2 - drop(t(x) %*% solve(S, x)) / 2
		}, numeric(1)))
	}, numeric(1)))
	count = 12
	oracle = function(rs) {
		v1 = 0.4; v2 = 1.3
		Q = -2 * (loglik(v1, rs) - loglik(v2, rs) + count / 2 * log(v1 / v2)) /
			(1 / v1 - 1 / v2)
		set.seed(572)
		1 / rgamma(20, (count + 1) / 2, (Q + s20) / 2)
	}
	truth = oracle(rep(rho[1], 2))
	set.seed(572)
	x = replicate(20, lame:::rs2_rep_fc_cpp(E, solve(matrix(c(1, rho[1], rho[1], 1), 2)), s20))
	expect_equal(x, truth, tolerance = 1e-12)
	set.seed(572)
	expect_equal(replicate(20, rs2_rep_fc(E, rho[1], s20)), truth, tolerance = 1e-12)
	truth_dynamic = oracle(rho)
	set.seed(572)
	expect_equal(
		replicate(20, lame:::.lame_sample_s2_dynamic_rho(E, rho, s20 = s20)),
		truth_dynamic,
		tolerance = 1e-12
	)
	# a fixed default also gives a proper conditional when residuals are zero
	expect_true(is.finite(lame:::rs2_rep_fc_cpp(E * 0, diag(2))))
	expect_error(lame:::rs2_rep_fc_cpp(E, diag(2), 0), "finite positive scalar")
})

test_that("model paths report and honor a residual prior set before sampling", {
	set.seed(171)
	Y = matrix(
		rnorm(30, sd = 0.3),
		5,
		6,
		dimnames = list(paste0("r", 1:5), paste0("c", 1:6))
	)
	common = list(
		family = "normal", R = 0, rvar = FALSE, cvar = FALSE,
		nscan = 40, burn = 10, odens = 2, gof = FALSE, verbose = FALSE
	)
	for (panel in c(FALSE, TRUE)) for (mode in c("unipartite", "bipartite")) {
		y = if (mode == "bipartite") Y else Y[, 1:5]
		if (mode == "unipartite") {
			colnames(y) = rownames(y)
			diag(y) = NA
		}
		if (panel) y = list(y, y + 0.1)
		fn = if (panel) lame else ame
		a = c(list(Y = y, mode = mode, dcor = FALSE), common)
		f = do.call(fn, a)
		expect_equal(f$s2_prior_scale, var(as.numeric(unlist(y)), na.rm = TRUE))
		a$prior = list(s20 = 0.6)
		f = do.call(fn, a)
		expect_equal(f$s2_prior_scale, 0.6)
		expect_true(all(is.finite(f$VC[, "ve"]) & f$VC[, "ve"] > 0))
		a$prior = list(s20 = -1)
		expect_error(do.call(fn, a), "s20.*finite positive scalar")
	}
	S = Y[, 1:5]; colnames(S) = rownames(S); S = (S + t(S)) / 2; diag(S) = NA
	for (panel in c(FALSE, TRUE)) {
		y = if (panel) list(S, S + 0.1) else S
		f = do.call(
			if (panel) lame else ame,
			c(list(Y = y, symmetric = TRUE, dcor = FALSE, nvar = FALSE, prior = list(s20 = 0.4)), common)
		)
		expect_equal(f$s2_prior_scale, 0.4)
		expect_true(all(is.finite(f$VC[, "ve"]) & f$VC[, "ve"] > 0))
	}
})

test_that("variance-scaled coefficient priors contribute their density", {
	beta = c(0.4, -0.7); g = c(0.5, 2)
	terms = lame:::.lame_s2_beta_terms(beta, g)
	E = array(seq(-0.8, 0.6, length.out = 18), c(3, 3, 2)); s20 = 0.6
	data = unlist(lapply(1:2, function(t) E[, , t][row(E[, , t]) != col(E[, , t])]))
	log_joint = function(v) sum(dnorm(data, 0, sqrt(v), log = TRUE)) +
		sum(dnorm(beta, 0, sqrt(g * v), log = TRUE))
	count = length(data) + length(beta); v = c(0.3, 1.4)
	Q = -2 * (log_joint(v[1]) - log_joint(v[2]) + count / 2 * log(v[1] / v[2])) /
		(1 / v[1] - 1 / v[2])
	set.seed(6682); truth = 1 / rgamma(20, (count + 1) / 2, (Q + s20) / 2)
	set.seed(6682)
	expect_equal(
		replicate(20, lame:::rs2_rep_fc_cpp(E, diag(2), s20, terms["df"], terms["ss"])),
		truth,
		tolerance = 1e-12
	)
	set.seed(6682)
	expect_equal(
		replicate(
			20,
			lame:::.lame_sample_s2_dynamic_rho(
				E,
				c(0, 0),
				s20 = s20,
				beta_df = terms["df"],
				beta_ss = terms["ss"]
			)
		),
		truth,
		tolerance = 1e-12
	)
	expect_equal(lame:::.lame_static_g_scale(4, 2L), 4)
	expect_equal(lame:::.lame_static_g_scale(c(2, 4, 6), c(2L, 3L)), c(4, 6))
	expect_equal(lame:::.lame_s2_beta_terms(numeric(0), numeric(0)), c(df = 0, ss = 0))
})

test_that("collapsed diagonal refresh uses new covariance and preserves real cells", {
	refresh = lame:::.lame_refresh_diagonal
	set.seed(871)
	Z = array(rnorm(18), c(3, 3, 2)); means = matrix(seq(-0.5, 0.5, length.out = 6), 3, 2)
	rho = c(-0.4, 0.7); s2 = 0.6
	z = refresh(Z, means, rho, s2)
	for (t in 1:2) {
		mask = row(Z[, , t]) != col(Z[, , t])
		expect_identical(z[, , t][mask], Z[, , t][mask])
	}
	x = replicate(3000, as.numeric(lame:::.lame_diagonal(refresh(Z, means, rho, s2))))
	truth_sd = rep(sqrt((1 + rho) * s2), each = 3)
	expect_true(all(abs(rowMeans(x) - as.numeric(means)) < 6 * truth_sd / sqrt(ncol(x))))
	expect_true(all(abs(apply(x, 1, sd) - truth_sd) < 6 * truth_sd / sqrt(2 * ncol(x))))
	for (nt in c(1L, 2L)) {
		z = array(0, c(1, 1, nt)); m = lame:::.lame_diagonal(z)
		expect_identical(dim(m), c(1L, nt))
		expect_identical(dim(refresh(z, m, 0, 1)), dim(z))
	}
})

test_that("one static coefficient with scalar g is sampled in dynamic-beta fits", {
	set.seed(932)
	Y = lapply(
		1:3,
		function(t) {z = matrix(rnorm(36), 6); diag(z) = NA; dimnames(z) = list(letters[1:6], letters[1:6]); z}
	)
	X = lapply(
		1:3,
		function(t) array(
			rnorm(36),
			c(6, 6, 1),
			dimnames = list(letters[1:6], letters[1:6], "x")
		)
	)
	f = lame(
		Y,
		Xdyad = X,
		R = 0,
		rvar = FALSE,
		cvar = FALSE,
		dcor = FALSE,
		dynamic_beta = "intercept",
		g = 5,
		prior = list(s20 = 0.7),
		burn = 10,
		nscan = 60,
		odens = 2,
		gof = FALSE,
		verbose = FALSE
	)
	expect_true(all(is.finite(f$BETA)))
	expect_length(dim(f$BETA), 3L)
	expect_gt(sd(f$BETA[, "x_dyad", 1L]), 0)
	for (t in 2:3) expect_identical(f$BETA[, "x_dyad", t], f$BETA[, "x_dyad", 1L])
	expect_equal(f$tryErrorChecks$beta_static, 0)
})
