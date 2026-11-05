skip_on_cran()

library(lame)

# compiled kernels of the binary probit empirical-Bayes estimator (src/ae_kernels.cpp) against their R
# references (suffix _r); they agree to rounding on realistic predictors, and the tolerance is looser where the
# references are noisy themselves: at the hand-over between the series and the log-domain integral, and far in the tails

kn = function(f) get(f, envir = asNamespace("lame"))
rel_diff = function(x, y, floor = 1e-8) max(abs(x - y) / (abs(y) + floor * max(abs(y[is.finite(y)]))))

test_that("the Gauss-Hermite rule of the kernels is the rule of the R reference to the last bit", {
	for (k in c(1L, 2L, 3L, 5L, 9L, 15L, 21L)) {
		g0 = lame:::.ae_mf_gh(k); g1 = kn("ae_gh_rule_cpp")(k)
		expect_identical(g1$z, g0$z)
		expect_identical(g1$w, g0$w)
	}
})

test_that("the log normal cdf kernel agrees with pnorm over the whole line", {
	x = c(seq(-60, 12, by = 0.173), -1e3, -1e5, 38, 100, 0)
	l0 = pnorm(x, log.p = TRUE); l1 = kn("ae_log_pnorm_cpp")(x)
	expect_lt(max(abs(l1 - l0) / pmax(abs(l0), 1e-300)), 1e-13)
	expect_identical(kn("ae_log_pnorm_cpp")(c(-Inf, Inf, NaN)), c(-Inf, 0, NaN))
})

test_that("the bivariate normal probability kernel reproduces pbivnorm", {
	skip_if_not_installed("pbivnorm")
	g = expand.grid(
		a = seq(-8, 8, by = 0.53),
		b = seq(-8, 8, by = 0.47),
		r = c(
			-0.99, -0.95, -0.93, -0.92, -0.75, -0.7, -0.3, -0.29, 0, 0.29, 0.3, 0.5, 0.75, 0.9, 0.924, 0.926, 0.97, 0.99
		)
	)
	p0 = pbivnorm::pbivnorm(g$a, g$b, g$r)
	p1 = kn("ae_pbivnorm_cpp")(g$a, g$b, g$r)
	expect_lt(max(abs(p1 - p0)), 1e-15)
	# infinite limits are the exact limits, a scalar correlation is recycled
	expect_equal(kn("ae_pbivnorm_cpp")(c(-Inf, Inf, Inf, 0.3), c(1, 0.5, Inf, Inf), 0.4)[1:3], c(0, pnorm(0.5), 1))
	expect_equal(kn("ae_pbivnorm_cpp")(0.3, Inf, 0.4), pnorm(0.3))
	# the result does not depend on the number of threads
	expect_identical(kn("ae_pbivnorm_cpp")(g$a, g$b, g$r, 4L), p1)
})

test_that("the log bivariate probability and its derivatives agree with the R reference", {
	skip_if_not_installed("pbivnorm")
	set.seed(5)
	n = 3000
	a = rnorm(n, -1, 2); b = rnorm(n, -1, 2); r = runif(n, -0.97, 0.97)
	l0 = lame:::.ae_pair_logphi2_r(a, b, r); l1 = kn("ae_pair_logphi2_cpp")(a, b, r)
	expect_lt(max(abs(l1 - l0)), 1e-10)
	d0 = lame:::.ae_pair_logp_r(a, b, r); d1 = kn("ae_pair_logp_cpp")(a, b, r)
	for (f in names(d0)) expect_lt(rel_diff(d1[[f]], d0[[f]]), 1e-8, label = f)
	# scalar correlation: the same as the vector of copies
	d2 = kn("ae_pair_logp_cpp")(a, b, 0.45)
	d3 = kn("ae_pair_logp_cpp")(a, b, rep(0.45, n))
	expect_identical(d2, d3)
	# a large and a very small threshold: the log-domain integral runs under the smaller limit
	# (a window under the larger one misses the mass of the integrand)
	a2 = c(31, 25, 40, 20, 14, -5.56); b2 = c(-5.56, -9, -12, -8, -6, 31); r2 = c(0, 0.7, -0.3, -0.5, 0.2, 0)
	truth = vapply(seq_along(a2), function(k) {
		lo = min(a2[k], b2[k]); hi = max(a2[k], b2[k]); kap = sqrt(1 - r2[k]^2)
		lf = function(x) dnorm(x, log = TRUE) + pnorm((hi - r2[k] * x) / kap, log.p = TRUE)
		m = max(vapply(seq(lo - 40, lo, length.out = 2001), lf, 0))
		m + log(integrate(function(x) exp(lf(x) - m), -Inf, lo, rel.tol = 1e-12, subdivisions = 1000L)$value)
	}, 0)
	expect_equal(truth[1], pnorm(-5.56, log.p = TRUE) + pnorm(31, log.p = TRUE), tolerance = 1e-12)
	expect_lt(max(abs(kn("ae_pair_logphi2_cpp")(a2, b2, r2) - truth)), 1e-8)
	expect_lt(max(abs(lame:::.ae_pair_logphi2_r(a2, b2, r2) - truth)), 1e-8)
	# very small probabilities take the log-domain integral, the same value as the reference
	ab = rbind(c(-8, -8, -0.5), c(-10, -5, -0.9), c(-6, -6, -0.3), c(-8, -8, 0.8), c(-3, -8, 0.8), c(-4, -4, 0.9))
	expect_lt(
		max(
			abs(
				kn("ae_pair_logphi2_cpp")(ab[, 1], ab[, 2], ab[, 3]) - lame:::.ae_pair_logphi2_r(ab[, 1], ab[, 2], ab[, 3])
			)
		),
		1e-10
	)
})

test_that("the cell quadrature kernels agree with the R references, observed cells only", {
	set.seed(6)
	nc = 6000
	y = array(sample(c(0, 1, NA), nc, TRUE, prob = c(.45, .45, .1)), c(nc, 1L, 1L))
	m = array(rnorm(nc, -1, 2), dim(y)); v = array(runif(nc, 0, 2), dim(y))
	for (link in c("probit", "logit")) for (k in c(3L, 9L)) {
		w0 = lame:::.ae_mf_cells_r(y, m, v, link, k); w1 = kn("ae_cells_cpp")(y, m, v, link == "probit", k)
		for (f in c("z", "w", "ell")) expect_lt(
			rel_diff(as.vector(w1[[f]]), as.vector(w0[[f]])),
			1e-10,
			label = paste(link, k, f)
		)
		expect_true(all(as.vector(w1$w)[is.na(y)] == 0))
		i0 = lame:::.ae_mf_info_w_r(y, m, v, link, k); i1 = kn("ae_info_w_cpp")(y, m, v, link == "probit", k)
		expect_lt(rel_diff(i1, as.vector(i0)), 1e-9, label = paste(link, k, "info_w"))
	}
	# predictors tens of standard deviations out: the cells are conditioned badly and both engines are noisy
	m[1:20] = runif(20, -40, -20); v[1:20] = 5
	w0 = lame:::.ae_mf_cells_r(y, m, v, "probit", 9L); w1 = kn("ae_cells_cpp")(y, m, v, TRUE, 9L)
	expect_lt(rel_diff(as.vector(w1$w), as.vector(w0$w)), 1e-6)
	# the wrapper takes the engine from the option, and the arrays keep their shape
	old = options(lame.als_engine = "r")
	on.exit(options(old), add = TRUE)
	r_wrap = lame:::.ae_mf_cells(y, m, v, "probit", 9L)
	options(lame.als_engine = "cpp")
	c_wrap = lame:::.ae_mf_cells(y, m, v, "probit", 9L)
	expect_identical(dim(c_wrap$z), dim(y))
	expect_identical(r_wrap, w0)
	expect_lt(rel_diff(as.vector(c_wrap$w), as.vector(r_wrap$w)), 1e-6)
})

pair_inputs = function(np = 800, seed = 7, mu = -1.2) {
	set.seed(seed)
	m1 = rnorm(np, mu, 1.5); m2 = rnorm(np, mu, 1.5)
	v1 = runif(np, 0, 1.5); v2 = runif(np, 0, 1.5)
	cv = runif(np, -0.9, 0.9) * sqrt(v1 * v2)
	p1 = pnorm(mu)
	list(
		m1 = m1,
		m2 = m2,
		v1 = v1,
		v2 = v2,
		cv = cv,
		s1 = ifelse(runif(np) < p1, 1, -1),
		s2 = ifelse(runif(np) < p1, 1, -1)
	)
}

test_that("the pair quadrature kernel (observed information) agrees with the R reference", {
	skip_if_not_installed("pbivnorm")
	d = pair_inputs()
	for (rho in c(-0.8, -0.3, 0, 0.5, 0.9)) for (k in c(3L, 5L)) {
		q0 = lame:::.ae_pair_quad_r(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, rho, k)
		q1 = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, rho, k, 0L)
		expect_identical(names(q1), names(q0))
		scale = max(abs(unlist(q0[c("psi1", "psi2", "w11", "w22")])))
		for (f in names(q0)) {
			# the cross information is zero at rho = 0: compare in absolute terms to the scale of the other columns
			expect_lt(max(abs(q1[[f]] - q0[[f]])) / max(abs(q0[[f]]), scale * 1e-3), 1e-8, label = paste(f, rho, k))
		}
	}
	# a node count per pair is the pair's own rule, and the result does not depend on the threads
	k = rep(c(3L, 5L), length.out = length(d$m1))
	q = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, 0.4, k, 0L)
	for (kk in c(3L, 5L)) {
		ix = which(k == kk)
		ref = kn("ae_pair_quad_cpp")(d$m1[ix], d$m2[ix], d$v1[ix], d$v2[ix], d$cv[ix], d$s1[ix], d$s2[ix], 0.4, kk, 0L)
		for (f in names(q)) expect_identical(q[[f]][ix], ref[[f]])
	}
	expect_identical(kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, 0.4, k, 0L, FALSE, FALSE, 3L), q)
})

test_that("the pair quadrature kernel (Fisher information) agrees with the R reference", {
	skip_if_not_installed("pbivnorm")
	d = pair_inputs(np = 500, seed = 8)
	k = rep(c(1L, 5L, 9L), length.out = length(d$m1))
	for (rho in c(-0.6, 0.4)) {
		q0 = lame:::.ae_pair_quad_fisher_r(
			d$m1,
			d$m2,
			d$v1,
			d$v2,
			0,
			d$s1,
			d$s2,
			rho,
			k,
			deriv = TRUE,
			rho_deriv = TRUE
		)
		q1 = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, rho, k, 1L, TRUE, TRUE)
		expect_identical(names(q1), names(q0))
		for (f in names(q0)) expect_lt(rel_diff(q1[[f]], q0[[f]], floor = 1e-6), 1e-7, label = paste(f, rho))
	}
	# the optional columns are the prefixes of the full list
	q2 = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.2, k, 1L, TRUE, FALSE)
	q3 = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.2, k, 1L, FALSE, FALSE)
	q4 = kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.2, k, 1L, TRUE, TRUE)
	expect_identical(names(q3), c("ell", "psi1", "psi2", "I11", "I22", "I12"))
	expect_identical(q2[names(q3)], q3)
	expect_identical(q4[names(q2)], q2)
	expect_identical(kn("ae_pair_quad_cpp")(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.2, k, 1L, FALSE, TRUE), q3)
	# at rho = 0 the cells are independent: no cross information, and the diagonal is the Fisher weight of a cell
	z = kn("ae_pair_quad_cpp")(c(-0.5, 0.7), c(0.3, -1), 0, 0, 0, c(1, -1), c(-1, -1), 0, 1L, 1L)
	expect_equal(z$I12, c(0, 0), tolerance = 1e-12)
	expect_equal(z$I11, dnorm(c(-0.5, 0.7))^2 / (pnorm(c(-0.5, 0.7)) * pnorm(-c(-0.5, 0.7))), tolerance = 1e-10)
	# the derivatives in the means and the correlation match finite differences of the expectations
	h = 1e-5
	f = function(m1, m2, rho, s1, s2) kn("ae_pair_quad_cpp")(m1, m2, 0.3, 0.4, 0, s1, s2, rho, 9L, 1L, TRUE, TRUE)
	for (cs in list(c(0.3, -0.8, 0.35, 1, -1), c(-1.2, -2.5, -0.4, -1, -1))) {
		q = f(cs[1], cs[2], cs[3], cs[4], cs[5])
		d_rho = function(nm) (f(
			cs[1],
			cs[2],
			cs[3] + h,
			cs[4],
			cs[5]
		)[[nm]] - f(cs[1], cs[2], cs[3] - h, cs[4], cs[5])[[nm]]) / (2 * h)
		d_m1 = function(nm) (f(
			cs[1] + h,
			cs[2],
			cs[3],
			cs[4],
			cs[5]
		)[[nm]] - f(cs[1] - h, cs[2], cs[3], cs[4], cs[5])[[nm]]) / (2 * h)
		d_m2 = function(nm) (f(
			cs[1],
			cs[2] + h,
			cs[3],
			cs[4],
			cs[5]
		)[[nm]] - f(cs[1], cs[2] - h, cs[3], cs[4], cs[5])[[nm]]) / (2 * h)
		expect_equal(q$s0, d_rho("ell"), tolerance = 1e-5)
		expect_equal(q$psi1r, d_rho("psi1"), tolerance = 1e-5)
		expect_equal(q$psi2r, d_rho("psi2"), tolerance = 1e-5)
		expect_equal(q$I11r, d_rho("I11"), tolerance = 1e-5)
		expect_equal(q$I22r, d_rho("I22"), tolerance = 1e-5)
		expect_equal(q$I12r, d_rho("I12"), tolerance = 1e-5)
		expect_equal(q$rc0, d_rho("s0"), tolerance = 1e-5)
		expect_equal(q$psi1, d_m1("ell"), tolerance = 1e-5)
		expect_equal(q$psi2, d_m2("ell"), tolerance = 1e-5)
		expect_equal(q$I11a, d_m1("I11"), tolerance = 1e-5)
		expect_equal(q$I11b, d_m2("I11"), tolerance = 1e-5)
		expect_equal(q$I22a, d_m1("I22"), tolerance = 1e-5)
		expect_equal(q$I22b, d_m2("I22"), tolerance = 1e-5)
		expect_equal(q$I12a, d_m1("I12"), tolerance = 1e-5)
		expect_equal(q$I12b, d_m2("I12"), tolerance = 1e-5)
	}
})

test_that("the pair quadrature wrappers select the engine with the option", {
	skip_if_not_installed("pbivnorm")
	d = pair_inputs(np = 200, seed = 9)
	old = options(lame.als_engine = "r", lame.als_threads = 1L)
	on.exit(options(old), add = TRUE)
	a = lame:::.ae_pair_quad(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, 0.3, 5L)
	f_a = lame:::.ae_pair_quad_fisher(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.3, 5L, TRUE, TRUE)
	options(lame.als_engine = "cpp")
	b = lame:::.ae_pair_quad(d$m1, d$m2, d$v1, d$v2, d$cv, d$s1, d$s2, 0.3, 5L)
	f_b = lame:::.ae_pair_quad_fisher(d$m1, d$m2, d$v1, d$v2, 0, d$s1, d$s2, 0.3, 5L, TRUE, TRUE)
	for (f in names(a)) expect_lt(rel_diff(b[[f]], a[[f]]), 1e-8, label = f)
	for (f in names(f_a)) expect_lt(rel_diff(f_b[[f]], f_a[[f]], floor = 1e-6), 1e-7, label = f)
	options(lame.als_engine = "nonsense")
	expect_error(lame:::.ae_engine(), "must be")
	options(lame.als_engine = "cpp", lame.als_threads = 0)
	expect_error(lame:::.ae_threads(), "positive")
})

test_that("the sums of the pair engine and the signed weighted sums equal their references", {
	set.seed(10)
	n = 25L; Tt = 3L; p = 2L
	X = array(rnorm(n * n * p * Tt), c(n, n, p, Tt)); X[sample(length(X), 40)] = NA
	naX = which(is.na(X), arr.ind = TRUE)
	cw = array(rnorm(n * n * Tt), c(n, n, Tt)); cw[sample(length(cw), 20)] = NA
	cw[cbind(naX[, 1], naX[, 2], naX[, 4])] = 0
	cx = list(nr = n, Tt = Tt, p = p, X = X)
	a = lame:::.ae_pair_cross_r(cx, cw); b = lame:::.ae_pair_cross(cx, cw)
	expect_identical(names(a), names(b))
	for (f in c("Csum", "gmc", "xx", "d", "th_a", "th_b")) expect_equal(b[[f]], a[[f]], tolerance = 1e-12, label = f)
	for (k in seq_along(a$Mk)) expect_equal(b$Mk[[k]], a$Mk[[k]], tolerance = 1e-12)
	hp = array(rnorm(n * n * Tt), c(n, n, Tt)); hp[hp > 1] = 0; hp[sample(length(hp), 15)] = NA
	hp[cbind(naX[, 1], naX[, 2], naX[, 4])] = 0
	expect_equal(
		lame:::.ae_eb_hprime_sums(hp, X, n, n, Tt, p),
		lame:::.ae_eb_hprime_sums_r(hp, X, n, n, Tt, p),
		tolerance = 1e-12
	)
	# the terms of the probit working response and the derivative of its weight
	eta = array(rnorm(n * n * Tt, 0, 3), c(n, n, Tt))
	t0 = lame:::.ae_probit_terms_r(eta, which(rbinom(length(eta), 1, .2) == 1)); t1 = lame:::.ae_probit_terms(eta, NULL)
	for (f in c("ld", "lp1", "lp0")) expect_equal(t1[[f]], t0[[f]], tolerance = 1e-13, label = f)
	expect_identical(dim(t1$lp1), dim(eta))
	h0 = lame:::.ae_eb_hprime(eta, "probit", t0)
	old = options(lame.als_engine = "r")
	h1 = lame:::.ae_eb_hprime(eta, "probit", t0)
	options(old)
	expect_equal(h0, h1, tolerance = 1e-12)
	# the log bivariate cdf of the correlation estimate is that of the sampler's function
	a = rnorm(500, 0, 1.5); b = rnorm(500, 0, 1.5); r = runif(500, -0.95, 0.95)
	expect_equal(lame:::.ae_bvn_log(a, b, r), lame:::.lame_binary_dynamic_bvn_log(a, b, r), tolerance = 1e-12)
	# the data sums of the bcd do not depend on the number of threads
	W = array(runif(n * n * Tt), c(n, n, Tt)); W[sample(length(W), 60)] = 0
	Z = array(rnorm(n * n * Tt), c(n, n, Tt))
	X0 = X; X0[is.na(X0)] = 0
	expect_identical(
		kn("ae_bcd_stats_cpp")(W, Z, X0, n, n, Tt, p, 1L),
		kn("ae_bcd_stats_cpp")(W, Z, X0, n, n, Tt, p, 4L)
	)
})

test_that("the row solves of the factor sweeps agree with the one-at-a-time reference, indefinite rows included", {
	set.seed(11)
	for (R in c(3L, 4L, 6L)) {
		n = 60L
		Gflat = t(
			vapply(
				seq_len(n),
				function(i) as.vector(crossprod(matrix(rnorm(R * R), R)) * runif(1, 0.1, 3)),
				numeric(R * R)
			)
		)
		Lm = crossprod(matrix(rnorm(R * R), R)) / 3 + diag(R)
		rhs = matrix(rnorm(n * R), n)
		a = lame:::.ae_mf_row_solve_r(Gflat, Lm, rhs); b = lame:::.ae_mf_row_solve(Gflat, Lm, rhs)
		expect_equal(b$x, a$x, tolerance = 1e-12)
		expect_equal(b$Sig, a$Sig, tolerance = 1e-12)
	}
	R = 3L; n = 6L
	Gflat = t(vapply(seq_len(n), function(i) as.vector(crossprod(matrix(rnorm(R * R), R))), numeric(R * R)))
	Gflat[2, ] = as.vector(diag(c(1, -5, 1)))
	Lm = diag(R) * 0.1; rhs = matrix(rnorm(n * R), n)
	a = lame:::.ae_mf_row_solve_r(Gflat, Lm, rhs); b = lame:::.ae_mf_row_solve(Gflat, Lm, rhs)
	expect_equal(b$x, a$x, tolerance = 1e-12)
	expect_equal(b$Sig, a$Sig, tolerance = 1e-12)
})

test_that("the effect step forms the inverse of the information and of its actor block from one factor", {
	set.seed(12)
	for (kind in c("directed", "bipartite", "symmetric")) for (pgn in c(0L, 2L)) {
		nr = 14L; nc = if (kind == "bipartite") 10L else nr
		Y = array(rbinom(nr * nc, 1, .3), c(nr, nc, 1)); if (kind != "bipartite") diag(Y[, , 1]) = NA
		p = 2L; X = array(rnorm(nr * nc * p), c(nr, nc, p, 1))
		prep = list(
			Yarr = Y,
			bip = kind == "bipartite",
			nr = nr,
			nc = nc,
			Tt = 1L,
			p_dyad = p,
			X = X,
			W_row = if (pgn) matrix(rnorm(nr * pgn), nr) else NULL,
			W_col = if (pgn && kind != "symmetric") matrix(rnorm(nc * pgn), nc) else NULL
		)
		ab = lame:::.ae_eb_prior(prep, kind == "symmetric")
		W = array(runif(nr * nc, .02, .2), c(nr, nc, 1)); Z = array(rnorm(nr * nc), c(nr, nc, 1))
		if (kind != "bipartite") { W[cbind(1:nr, 1:nr, 1)] = 0; Z[cbind(1:nr, 1:nr, 1)] = 0 }
		st = kn("ae_bcd_stats_cpp")(W, Z, X, nr, nc, 1L, p)
		ctx = lame:::.ae_eb_context(ab, st, rowSums(W, dims = 2L), p)
		sab = matrix(c(.9, -.3, -.3, .6), 2)
		Lam = lame:::.ae_eb_precision(sab, ab)
		# the penalty added by index is the penalty added as a matrix
		Lbare = Lam; attr(Lbare, "ae_prior") = NULL
		H = lame:::.ae_eb_penalise(ctx, Lam)
		expect_identical(H, lame:::.ae_eb_penalise(ctx, Lbare))
		s = ctx$i_s
		inv = lame:::.ae_eb_inverse_pair(H, s)
		expect_equal(inv$S, lame:::.ae_eb_inverse(H), tolerance = 1e-10)
		expect_equal(inv$Sc, lame:::.ae_eb_inverse(H[s, s]), tolerance = 1e-10)
	}
	# a singular information (a collinear design) falls back to the separate inverses and their pseudo-inverse route
	H = diag(6); H[1:2, 1:2] = 1; s = 3:6
	inv = lame:::.ae_eb_inverse_pair(H, s)
	expect_equal(inv$S, lame:::.ae_eb_inverse(H), tolerance = 1e-12)
})

als_engine_panel = function(n = 14L, Tt = 2L, seed = 13) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n))
	Y = X = vector("list", Tt)
	a = rnorm(n, 0, .5); b = rnorm(n, 0, .5); U = rnorm(n, 0, .7); V = rnorm(n, 0, .7)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n)
		e = matrix(rnorm(n * n), n); z2 = matrix(rnorm(n * n), n)
		e[lower.tri(e)] = (0.5 * t(e) + sqrt(0.75) * z2)[lower.tri(e)]
		y = 1 * (-0.7 + 0.6 * x + outer(a, b, "+") + tcrossprod(U, V) + e > 0); diag(y) = NA
		dimnames(y) = list(rn, rn); Y[[t]] = y
		X[[t]] = array(x, c(n, n, 1), dimnames = list(rn, rn, "x"))
	}
	names(Y) = names(X) = paste0("t", seq_len(Tt))
	list(Y = Y, X = X)
}

als_engine_fit = function(d, engine, threads, dcor, R) {
	old = options(lame.als_engine = engine, lame.als_threads = threads)
	on.exit(options(old), add = TRUE)
	suppressWarnings(suppressMessages(lame_als(d$Y, d$X, R = R, family = "binary", verbose = FALSE, dcor = dcor)))
}

als_engine_same = function(f_c, f_r, R) {
	expect_equal(unname(f_c$coefficients), unname(f_r$coefficients), tolerance = 1e-7)
	expect_equal(f_c$VC[["rho"]], f_r$VC[["rho"]], tolerance = 1e-6)
	expect_equal(f_c$eb$variational$iterations, f_r$eb$variational$iterations)
	if (R > 0L) expect_equal(f_c$UVPM, f_r$UVPM, tolerance = 1e-6)
	expect_equal(unname(f_c$eb$vcov), unname(f_r$eb$vcov), tolerance = 1e-5)
}

test_that("the estimates are the same on either engine with the dyadic correlation, and for any thread count", {
	skip_if_not_installed("pbivnorm")
	d = als_engine_panel()
	f_r = als_engine_fit(d, "r", 1L, TRUE, 1L)
	f_c = als_engine_fit(d, "cpp", 1L, TRUE, 1L)
	als_engine_same(f_c, f_r, 1L)
	f3 = als_engine_fit(d, "cpp", 3L, TRUE, 1L)
	expect_identical(f3$coefficients, f_c$coefficients)
	expect_identical(f3$UVPM, f_c$UVPM)
})

test_that("the estimates are the same on either engine for independent cells and for rank zero", {
	skip_on_cran()
	skip_if_not_installed("pbivnorm")
	d = als_engine_panel()
	for (cfg in list(c(FALSE, 0L), c(FALSE, 1L), c(TRUE, 0L))) {
		f_r = als_engine_fit(d, "r", 1L, as.logical(cfg[1]), as.integer(cfg[2]))
		f_c = als_engine_fit(d, "cpp", 1L, as.logical(cfg[1]), as.integer(cfg[2]))
		als_engine_same(f_c, f_r, as.integer(cfg[2]))
	}
})
