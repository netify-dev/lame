skip_on_cran()

library(lame)

# the empirical-Bayes Poisson estimator (ame_als(), lame_als() with family = "poisson"): the Poisson-lognormal cell
# (src/ae_pln_kernels.cpp, R/ame_als_pln.R), its quadrature against an independent integral, the compiled engine
# against the R reference, and the fits it produces (zero-count actors, every kind of network, prediction)

kn = function(f) get(f, envir = asNamespace("lame"))
rule = function(K = 24L) lame:::.ae_pln_rule(K)
# the largest difference relative to the size of the reference, with a floor for references near zero
agree = function(x, y, tol, floor = 1e-6) {
	x = as.vector(x); y = as.vector(y)
	ok = max(abs(x - y) / pmax(abs(y), floor)) < tol
	expect_true(ok)
}

# the cell by numerical integration of the posterior of the noise, independent of the quadrature rule of the kernel
pln_integral = function(y, eta, ve) {
	f = function(e) exp(y * e - exp(eta + e) - e^2 / (2 * ve))
	zeta = lame:::.ae_pln_node_r(y, eta, ve)
	e0 = zeta$me
	lo = e0 - 25 * sqrt(ve); hi = e0 + 25 * sqrt(ve)
	br = sort(
		unique(c(lo, e0 - 6 * sqrt(ve), e0 - 2 * sqrt(ve), e0, e0 + 2 * sqrt(ve), min(hi, e0 + 6 * sqrt(ve)), hi))
	)
	I = function(g) sum(
		vapply(
			seq_len(length(br) - 1L),
			function(k) integrate(g, br[k], br[k + 1L], rel.tol = 1e-12, subdivisions = 1000L)$value,
			0
		)
	)
	I0 = I(f); I1 = I(function(e) e * f(e)); I2 = I(function(e) e^2 * f(e)); I3 = I(function(e) e^3 * f(e))
	m1 = I1 / I0; m2 = I2 / I0; k3 = I3 / I0 - 3 * m1 * m2 + 2 * m1^3
	list(
		l = y * eta - lgamma(y + 1) - 0.5 * log(2 * pi * ve) + log(I0),
		s = m1 / ve,
		w = (ve - (m2 - m1^2)) / ve^2,
		hp = -k3 / ve^3,
		ee = m2
	)
}

test_that("the Poisson-lognormal cell agrees with direct integration of the posterior of the noise", {
	g = expand.grid(y = c(0, 1, 4, 30), eta = c(-4, -1, 0.5, 2.5), ve = c(0.05, 0.4, 1.5))
	nd = lame:::.ae_pln_node_r(g$y, g$eta, g$ve)
	for (i in seq_len(nrow(g))) {
		ref = pln_integral(g$y[i], g$eta[i], g$ve[i])
		expect_equal(nd$l[i], ref$l, tolerance = 1e-6, label = paste("l", i))
		expect_equal(nd$s[i], ref$s, tolerance = 1e-6, label = paste("s", i))
		expect_equal(nd$w[i], ref$w, tolerance = 1e-4, label = paste("w", i))
		expect_equal(nd$ee[i], ref$ee, tolerance = 1e-5, label = paste("ee", i))
		expect_lt(abs(nd$hp[i] - ref$hp) / max(abs(ref$hp), 1e-3), 1e-3, label = paste("hp", i))
	}
})

test_that("the score, the information and its derivative are the derivatives of the log likelihood", {
	g = expand.grid(y = c(0, 2, 12), eta = c(-3, 0, 1.5), ve = c(0.1, 0.8))
	h = 1e-4
	at = function(d) lame:::.ae_pln_node_r(g$y, g$eta + d, g$ve)
	n0 = at(0); np = at(h); nm = at(-h)
	expect_equal((np$l - nm$l) / (2 * h), n0$s, tolerance = 1e-5)
	expect_equal(-(np$s - nm$s) / (2 * h), n0$w, tolerance = 1e-4)
	agree((np$w - nm$w) / (2 * h), n0$hp, 1e-3, floor = 1)
	# the derivative in the variance is the second moment of the noise: dl/dve = (E[e^2] - ve) / (2 ve^2)
	hv = 1e-4
	lv = function(f) lame:::.ae_pln_node_r(g$y, g$eta, g$ve * (1 + f))$l
	expect_equal((lv(hv) - lv(-hv)) / (2 * g$ve * hv), (n0$ee - g$ve) / (2 * g$ve^2), tolerance = 1e-4)
	# the log likelihood is concave in the predictor
	expect_true(all(n0$w > 0))
})

test_that("a vanishing noise variance gives the Poisson cell", {
	y = c(0, 2, 7, 25); eta = c(-1, 0.5, 2, 3)
	nd = lame:::.ae_pln_node_r(y, eta, 1e-7)
	expect_equal(nd$l, dpois(y, exp(eta), log = TRUE), tolerance = 1e-5)
	expect_equal(nd$s, y - exp(eta), tolerance = 1e-5)
	expect_equal(nd$w, exp(eta), tolerance = 1e-5)
	expect_equal(nd$hp, exp(eta), tolerance = 1e-4)
})

test_that("the cell stays finite at extreme predictors, counts and variances", {
	g = expand.grid(y = c(0, 1, 50, 1e5), eta = c(-300, -30, 0, 30, 300), ve = c(1e-6, 1, 100))
	nd = lame:::.ae_pln_node_r(g$y, g$eta, g$ve, 24L)
	for (f in c("l", "s", "w", "hp", "ee", "me")) expect_true(all(is.finite(nd[[f]])), label = f)
	expect_true(all(nd$w >= 0))
	cp = kn("ae_pln_node_cpp")(g$y, g$eta, g$ve, lgamma(g$y + 1), rule()$z, rule()$lw, 1L)
	for (f in c("l", "s", "w", "ee")) agree(cp[[f]], nd[[f]], 1e-9, floor = 1e-3)
})

test_that("the compiled cell kernels agree with the R references, observed cells only", {
	set.seed(7)
	n = 2500
	y = rpois(n, 2) * (runif(n) < 0.8) + rpois(n, 40) * (runif(n) < 0.05)
	y[sample(n, 30)] = NA
	m = rnorm(n, 0.3, 1.3); v = runif(n, 0, 1.5)
	r = rule()
	r28 = rule(lame:::.ae_pln_kin(0.6))
	nd0 = lame:::.ae_pln_node_r(y[!is.na(y)], m[!is.na(y)], 0.6, 24L)
	nd1 = kn("ae_pln_node_cpp")(y, m, 0.6, lgamma(y + 1), r$z, r$lw, 1L)
	for (f in names(nd0)) agree(nd1[[f]][!is.na(y)], nd0[[f]], 1e-10)
	expect_true(all(nd1$w[is.na(y)] == 0))
	# the outer rules by the variance of the predictor (3, 5, 7 or 9 nodes), and one rule of any size
	for (K in c(9L, 5L)) {
		og = lame:::.ae_pln_outer(K)
		c0 = lame:::.ae_pln_cells_r(y, m, v, 0.6, K)
		c1 = kn("ae_pln_cells_cpp")(y, lgamma(y + 1), m, v, 0.6, og$z, og$w, og$vcut, r28$z, r28$lw, 1L)
		for (f in c("z", "w", "ell", "ee")) agree(c1[[f]], c0[[f]], 1e-10)
		expect_true(all(c1$w[is.na(y)] == 0))
		# the result does not depend on the number of threads
		expect_identical(
			kn("ae_pln_cells_cpp")(y, lgamma(y + 1), m, v, 0.6, og$z, og$w, og$vcut, r28$z, r28$lw, 3L),
			c1
		)
	}
	og = lame:::.ae_pln_outer(9L)
	expect_length(og$z, 4L)
	expect_equal(lengths(og$z), c(3L, 5L, 7L, 9L))
	# the rules that depend on the variance of the cell are within 1e-6 of the rule of forty-one nodes
	vt = c(1e-4, 0.003, 0.02, 0.1, 0.2, 0.5)
	yt = rep(c(0, 1, 4), length.out = length(vt)); mt = c(-1, 0.3, 1, 0.2, -0.4, 0.8)
	a41 = lame:::.ae_pln_cells_r(yt, mt, vt, 0.5, 41L); aa = lame:::.ae_pln_cells_r(yt, mt, vt, 0.5, 9L)
	agree(
		aa$ell,
		a41$ell,
		1e-6,
		floor = 1
	); agree(aa$w, a41$w, 1e-5, floor = 1e-3); agree(aa$z, a41$z, 1e-5, floor = 1e-3)
	w0 = lame:::.ae_pln_info_w_r(y, m, v, 0.6) # the rule of the largest total variance
	w1 = kn("ae_pln_info_w_cpp")(y, lgamma(y + 1), m, v, 0.6, lame:::.ae_pln_rule(
		lame:::.ae_pln_kin(max(0.6 + v, na.rm = TRUE))
	)$z,
		lame:::.ae_pln_rule(lame:::.ae_pln_kin(max(0.6 + v, na.rm = TRUE)))$lw, 1L)
	agree(w1, w0, 1e-10)
	# at a vanishing predictor variance the expected cell is the cell
	c2 = lame:::.ae_pln_cells_r(y, m, 0, 0.6, 9L, Kin = 24L)
	expect_equal(
		as.vector(c2$ell)[!is.na(y)],
		lame:::.ae_pln_node_r(y[!is.na(y)], m[!is.na(y)], 0.6, 24L)$l,
		tolerance = 1e-12
	)
})

test_that("the cavity weight is the information of the cell at the total variance", {
	# the count averaged over the predictor and the noise is a Poisson-lognormal count with variance v + ve
	y = c(0, 1, 3, 9); m = c(0.2, -0.5, 1, 1.8); v = c(0.3, 1, 0.1, 0.6); ve = 0.5
	w = lame:::.ae_pln_info_w_r(y, m, v, ve)
	nd = lame:::.ae_pln_node_r(y, m, ve + v)
	expect_equal(as.vector(w), nd$w)
	# and it is minus the second derivative in the mean of the log of the likelihood averaged over the predictor
	Z = function(mm) {
		og = lame:::.ae_mf_gh(41L)
		vapply(
			seq_along(y),
			function(i) sum(og$w * exp(lame:::.ae_pln_node_r(rep(y[i], 41L), mm[i] + sqrt(v[i]) * og$z, ve)$l)),
			0
		)
	}
	h = 1e-3
	expect_equal(-(log(Z(m + h)) - 2 * log(Z(m)) + log(Z(m - h))) / h^2, as.vector(w), tolerance = 1e-4)
})

test_that("the variance of the noise has the em image of the sampler's full conditional", {
	# (nu0 s20 + SS) / (N + nu0 - 2) with nu0 = s20 = 1: the posterior mean of the inverse gamma that rs2_fc() draws from
	expect_equal(lame:::.ae_pln_ve_image(40, 100), (1 + 40) / 99)
	expect_equal(lame:::.ae_pln_ve_image(40, 100, 0.5), (1 + 20) / 49)
	# with no data the image is the floor, and the hyperprior term is the log of the density it comes from
	expect_equal(lame:::.ae_pln_ve_image(0, 0), 1e-8)
	h = 1e-6
	N = 80; SS = 30
	obj = function(ve) -(N / 2) * log(ve) - SS / (2 * ve) + lame:::.ae_pln_ve_logprior(ve)
	ve_hat = lame:::.ae_pln_ve_image(SS, N)
	expect_equal((obj(ve_hat + h) - obj(ve_hat - h)) / (2 * h), 0, tolerance = 1e-5)
	# the geometric tail of a slowly converging sequence is extrapolated, and a short or erratic one is not
	x = log(0.2) + 0.4 * 0.985^(0:40)
	expect_equal(lame:::.ae_pln_aitken(x[1:30]), 0.5 * (log(0.2) - x[30]), tolerance = 0.02)
	expect_true(is.na(lame:::.ae_pln_aitken(x[1:8])))
	expect_true(is.na(lame:::.ae_pln_aitken(c(x[1:20], rev(x[1:11])))))
})

# the pair of cells of a dyad

pair_brute = function(y1, y2, m1, m2, s11, s22, s12) {
	Sig = matrix(c(s11, s12, s12, s22), 2); P = solve(Sig); dS = det(Sig)
	f = function(e1, e2) exp(
		y1 * (m1 + e1) - exp(m1 + e1) + y2 * (m2 + e2) - exp(m2 + e2) - lgamma(y1 + 1) - lgamma(y2 + 1) -
			0.5 * (P[1, 1] * e1^2 + 2 * P[1, 2] * e1 * e2 + P[2, 2] * e2^2)
	) / (2 * pi * sqrt(dS))
	lim = 9 * sqrt(max(s11, s22)) + 4
	integrate(
		function(e1) vapply(
			e1,
			function(x) integrate(function(e2) f(x, e2), -lim, lim, rel.tol = 1e-10, subdivisions = 400L)$value,
			0
		),
		-lim,
		lim,
		rel.tol = 1e-10,
		subdivisions = 400L
	)$value
}

test_that("the pair of cells of a dyad agrees with direct integration, and is two cells when the noises are independent", {
	cases = list(c(0, 1, -0.5, 0.2, 0.4, 0.5, 0.25), c(3, 2, 0.5, 1, 0.6, 0.5, 0.4), c(5, 1, 1.5, 0, 0.3, 0.3, -0.15))
	for (cs in cases) {
		nd = lame:::.ae_pln2_node_r(cs[1], cs[2], cs[3], cs[4], cs[5], cs[6], cs[7], K = 16L)
		expect_equal(nd$L, log(pair_brute(cs[1], cs[2], cs[3], cs[4], cs[5], cs[6], cs[7])), tolerance = 1e-8)
		# the score and the information are the first two derivatives of the log probability in the means
		h = 1e-4
		L = function(d1, d2) lame:::.ae_pln2_node_r(
			cs[1],
			cs[2],
			cs[3] + d1,
			cs[4] + d2,
			cs[5],
			cs[6],
			cs[7],
			K = 16L
		)$L
		expect_equal(nd$L1, (L(h, 0) - L(-h, 0)) / (2 * h), tolerance = 1e-6)
		expect_equal(nd$L2, (L(0, h) - L(0, -h)) / (2 * h), tolerance = 1e-6)
		expect_equal(nd$w11, -(L(h, 0) - 2 * nd$L + L(-h, 0)) / h^2, tolerance = 1e-4)
		expect_equal(nd$w22, -(L(0, h) - 2 * nd$L + L(0, -h)) / h^2, tolerance = 1e-4)
		expect_equal(nd$w12, -(L(h, h) - L(h, -h) - L(-h, h) + L(-h, -h)) / (4 * h^2), tolerance = 1e-4)
	}
	# independent noises: the log probabilities add, the cells' informations are theirs and the cross information is zero
	y1 = c(0, 2, 5, 1); y2 = c(1, 0, 3, 7); m1 = c(-1, 0.5, 1, 0); m2 = c(0.2, -0.5, 0.8, 1.5); s = c(0.4, 0.4, 0.8, 0.3)
	nd = lame:::.ae_pln2_node_r(y1, y2, m1, m2, s, s, 0, K = 16L)
	n1 = lame:::.ae_pln_node_r(y1, m1, s, 28L); n2 = lame:::.ae_pln_node_r(y2, m2, s, 28L)
	expect_equal(nd$L, n1$l + n2$l, tolerance = 1e-7)
	expect_equal(nd$w11, n1$w, tolerance = 1e-5); expect_equal(nd$w22, n2$w, tolerance = 1e-5)
	expect_lt(max(abs(nd$w12)), 1e-12)
	# positively correlated noises make the cross information negative, negatively correlated ones positive
	ps = lame:::.ae_pln2_node_r(c(2, 2), c(2, 2), c(0.3, 0.3), c(0.3, 0.3), 0.5, 0.5, c(0.3, -0.3), K = 12L)
	expect_lt(ps$w12[1], 0); expect_gt(ps$w12[2], 0)
})

test_that("the compiled pair kernel agrees with the R reference", {
	set.seed(8)
	n = 400
	y1 = rpois(n, 2); y2 = rpois(n, 2); y1[1:3] = c(0, 40, 1e3); m1 = rnorm(n, 0.3, 1); m2 = rnorm(n, 0.3, 1)
	s11 = runif(n, 0.1, 1.2); s22 = runif(n, 0.1, 1.2); s12 = runif(n, -0.8, 0.8) * sqrt(s11 * s22)
	for (K in c(8L, 12L)) {
		a = lame:::.ae_pln2_node_r(y1, y2, m1, m2, s11, s22, s12, K)
		r = lame:::.ae_pln_rule(K)
		b = kn("ae_pln2_node_cpp")(y1, y2, lgamma(y1 + 1), lgamma(y2 + 1), m1, m2, s11, s22, s12, r$z, r$lw, 1L)
		for (f in names(a)) agree(b[[f]], a[[f]], 1e-8, floor = 1e-3)
		expect_identical(
			kn("ae_pln2_node_cpp")(y1, y2, lgamma(y1 + 1), lgamma(y2 + 1), m1, m2, s11, s22, s12, r$z, r$lw, 3L),
			b
		)
	}
})

# the fits

sim_pln_net = function(n = 24, R = 0, mu = -0.3, beta = 0.5, ve = 0.4, seed = 1, zero = 0L, kind = "directed", m = n) {
	set.seed(seed)
	if (kind == "symmetric") m = n
	x = matrix(rnorm(n * m), n, m)
	if (kind == "symmetric") x = (x + t(x)) / sqrt(2)
	a = rnorm(n, 0, 0.7); b = if (kind == "symmetric") a else rnorm(m, 0, 0.7)
	uv = if (R > 0) tcrossprod(matrix(rnorm(n * R, 0, 0.7), n, R), matrix(rnorm(m * R, 0, 0.7), m, R)) else 0
	if (kind == "symmetric" && R > 0) { u = matrix(rnorm(n * R, 0, 0.7), n, R); uv = tcrossprod(u) }
	eta = mu + beta * x + outer(a, b, "+") + uv
	e = matrix(rnorm(n * m, 0, sqrt(ve)), n, m)
	if (kind == "symmetric") e[lower.tri(e)] = t(e)[lower.tri(e)]
	Y = matrix(rpois(n * m, exp(eta + e)), n, m)
	if (kind == "symmetric") Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]
	if (zero > 0) { Y[seq_len(zero), ] = 0; Y[, seq_len(zero) + (if (kind == "symmetric") 0L else 3L)] = 0 }
	if (kind != "bipartite") diag(Y) = NA
	rn = sprintf("r%02d", seq_len(n)); cn = if (kind == "bipartite") sprintf("c%02d", seq_len(m)) else rn
	dimnames(Y) = list(rn, cn)
	list(Y = Y, X = array(x, c(n, m, 1L), dimnames = list(rn, cn, "x")))
}

test_that("an actor without a count gets a finite, shrunk effect, as in the sampler", {
	d = sim_pln_net(24, zero = 1L, seed = 3)
	# the fixed-effects fit sends the effect of a zero-count actor to minus infinity
	fx = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, bias_correct = "none", verbose = FALSE))
	expect_lt(fx$a[1], -10)
	expect_identical(fx$bias_correct$method, "none")
	expect_null(fx$eb)
	fit = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE))
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(isTRUE(fit$bias_correct$default))
	expect_identical(fit$eb$estimator, "variational")
	expect_true(fit$converged)
	expect_true(all(is.finite(c(fit$a, fit$b))))
	expect_gt(fit$a[1], -5)
	expect_lt(fit$a[1], min(fit$a[-1]))
	expect_true(is.finite(fit$VC[["ve"]]) && fit$VC[["ve"]] > 0)
	expect_equal(unname(fit$s2), unname(fit$VC[["ve"]]))
	expect_equal(unname(fit$eb$ve), unname(fit$VC[["ve"]]))
	# the reported predictor is on the latent scale and the fitted counts add the lognormal mean
	expect_equal(unname(fit$fitted[[1]][2, 3]), unname(exp(fit$EZ[[1]][2, 3] + fit$VC[["ve"]] / 2)))
	# requesting the empirical-Bayes estimator by name is the same fit
	fit2 = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, bias_correct = "eb", verbose = FALSE))
	expect_equal(fit2$coefficients, fit$coefficients)
	expect_false(isTRUE(fit2$bias_correct$default))
	# the transformed route of a count fit has no such estimator (a normal fit has it, test_als_normal_eb.R)
	expect_error(
		ame_als(d$Y, family = "poisson", non_normal_method = "transform", bias_correct = "eb", verbose = FALSE),
		class = "lame_als_eb"
	)
})

test_that("the prior of the actor effects of a count network is the sampler's", {
	d = sim_pln_net(20, seed = 5)
	prep = lame:::.ae_prepare(d$Y, d$X, NULL, NULL, "unipartite", FALSE, node_time_policy = "mean", symmetric = FALSE)
	ab = lame:::.ae_eb_prior(prep, FALSE, "poisson")
	expect_equal(ab$Sab0, diag(2) * lame:::.ame_pois_effect_scale(d$Y))
	expect_equal(ab$eta0, round(4 + 3 * 20 / 100))
	expect_identical(ab$family, "poisson")
	# the binary prior is untouched
	expect_null(lame:::.ae_eb_prior(prep, FALSE)$family)
})

test_that("a count network with no finite estimate is refused, not returned unconverged", {
	d = sim_pln_net(16, seed = 6)
	Y0 = d$Y; Y0[!is.na(Y0)] = 0
	expect_error(ame_als(Y0, Xdyad = d$X, family = "poisson", verbose = FALSE), class = "lame_als_separation")
	# a dyadic covariate of one sign that is non-zero only where every count is zero
	blk = matrix(0, 16, 16); blk[1:4, 5:8] = 1
	Yb = d$Y; Yb[blk == 1] = 0
	Xb = array(c(d$X[, , 1], blk), c(16, 16, 2), dimnames = list(rownames(Yb), colnames(Yb), c("x", "blk")))
	expect_error(
		ame_als(Yb, Xdyad = Xb, family = "poisson", verbose = FALSE),
		class = "lame_als_separation",
		regexp = "blk"
	)
	# the fixed-effects fit is not refused here (it returns what it can)
	expect_no_error(
		suppressWarnings(ame_als(Yb, Xdyad = Xb, family = "poisson", bias_correct = "none", verbose = FALSE))
	)
})

test_that("the compiled and the R engines give the same Poisson fit", {
	skip_on_cran()
	d = sim_pln_net(16, R = 1, seed = 4)
	op = options(lame.als_engine = "r"); on.exit(options(op))
	f_r = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 1, verbose = FALSE))
	options(lame.als_engine = "cpp")
	f_c = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 1, verbose = FALSE))
	expect_equal(f_c$coefficients, f_r$coefficients, tolerance = 1e-6)
	expect_equal(f_c$VC[["ve"]], f_r$VC[["ve"]], tolerance = 1e-6)
	expect_equal(f_c$a, f_r$a, tolerance = 1e-6)
	expect_equal(f_c$eb$se, f_r$eb$se, tolerance = 1e-5)
})

test_that("every kind of network, with and without a latent rank, converges to a variational fit", {
	skip_on_cran()
	cases = list(
		list(label = "directed R0", d = sim_pln_net(22, seed = 5), R = 0),
		list(label = "directed R1", d = sim_pln_net(22, R = 1, seed = 6), R = 1),
		list(
			label = "bipartite R1",
			d = sim_pln_net(22, m = 16, R = 1, seed = 7, kind = "bipartite"),
			R = 1,
			mode = "bipartite"
		),
		list(label = "symmetric R1", d = sim_pln_net(22, R = 1, seed = 8, kind = "symmetric"), R = 1, symmetric = TRUE)
	)
	for (cs in cases) {
		args = list(cs$d$Y, Xdyad = cs$d$X, R = cs$R, family = "poisson", verbose = FALSE)
		if (!is.null(cs$mode)) args$mode = cs$mode
		if (!is.null(cs$symmetric)) args$symmetric = cs$symmetric
		fit = suppressWarnings(do.call(ame_als, args))
		expect_true(fit$converged, label = cs$label)
		expect_identical(fit$eb$estimator, "variational", label = cs$label)
		expect_true(all(is.finite(fit$eb$se)), label = cs$label)
		expect_true(all(is.finite(fit$coefficients)), label = cs$label)
		expect_gt(fit$VC[["ve"]], 0)
		expect_equal(unname(fit$coefficients["x_dyad"]), 0.5, tolerance = 0.5, label = cs$label)
	}
})

test_that("a longitudinal count panel is fitted by the empirical-Bayes estimator", {
	skip_on_cran()
	set.seed(9)
	n = 18
	Ys = lapply(1:3, function(t) { d = sim_pln_net(n, seed = 20 + t); d$Y })
	Xs = lapply(1:3, function(t) { d = sim_pln_net(n, seed = 20 + t); d$X })
	fit = suppressWarnings(lame_als(Ys, Xdyad = Xs, R = 0, family = "poisson", verbose = FALSE))
	expect_identical(fit$bias_correct$method, "eb")
	expect_true(fit$converged)
	expect_length(fit$fitted, 3L)
})

test_that("predictions, fitted counts, simulation and covariances of a Poisson fit follow the model", {
	skip_on_cran()
	d = sim_pln_net(20, seed = 11)
	fit = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE))
	ve = fit$VC[["ve"]]
	# the default of a count is the plug-in mean exp(EZ + ve / 2), as the samplers report it: the posterior predictive mean of
	# a count is not finite in this model; the mean that integrates the predictor is plug_in = FALSE
	pl = predict(fit)
	pr = predict(fit, plug_in = FALSE)
	lk = predict(fit, type = "link")
	expect_equal(pl, predict(fit, plug_in = TRUE))
	# the integrated mean adds half the variance of the predictor to the plug-in mean
	expect_true(all(pr[!is.na(pr)] > pl[!is.na(pl)] - 1e-12))
	expect_equal(unname(pl[2, 3]), unname(exp(lk[2, 3] + ve / 2)))
	expect_equal(fitted(fit), pl)
	expect_equal(fitted(fit, plug_in = TRUE), pl)
	expect_equal(fitted(fit, plug_in = FALSE), pr)
	off = row(pl) != col(pl)
	expect_equal(unname(pl[off]), unname(fit$fitted[[1]][off]), tolerance = 1e-12)
	expect_equal(unname(pl[off]), unname(fit$YPM[[1]][off]), tolerance = 1e-12)
	expect_equal(unname(residuals(fit)[off]), unname((d$Y - pl)[off]), tolerance = 1e-12)
	expect_equal(unname(residuals(fit, plug_in = FALSE)[off]), unname((d$Y - pr)[off]), tolerance = 1e-12)
	expect_error(fitted(fit, plug_in = NA), "plug_in")
	expect_error(predict(fit, plug_in = "yes"), "plug_in")
	V = vcov(fit)
	expect_true(all(is.finite(V)))
	expect_equal(sqrt(diag(V)), fit$eb$se, tolerance = 1e-8)
	ci = confint(fit)
	expect_true(all(ci[, 1] < ci[, 2]))
	sm = simulate(fit, nsim = 2, seed = 1)
	expect_length(sm$Y, 2L)
	expect_true(all(sm$Y[[1]][!is.na(sm$Y[[1]])] >= 0))
})

test_that("a refit on new counts uses the empirical-Bayes estimator", {
	skip_on_cran()
	d = sim_pln_net(18, seed = 12)
	fit = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE))
	d2 = sim_pln_net(18, seed = 13)
	rf = suppressWarnings(
		ame_als_refit(
			fit,
			Y_new = d2$Y,
			X_new = array(d2$X, c(dim(d2$X), 1L), dimnames = c(dimnames(d2$X), list(NULL)))
		)
	)
	expect_identical(rf$bias_correct$method, "eb")
	expect_true(all(is.finite(rf$eb$se)))
})

test_that("the dyadic correlation of a count network is estimated from the pairs and enters the standard errors", {
	skip_on_cran()
	set.seed(12)
	n = 36
	x = matrix(rnorm(n * n), n, n); x = (x + t(x)) / sqrt(2) # a covariate that is the same in both cells of a dyad
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	e1 = matrix(rnorm(n * n, 0, sqrt(0.6)), n, n); e2 = matrix(rnorm(n * n, 0, sqrt(0.6)), n, n)
	e = e1; up = upper.tri(e); et = t(e1); et[up] = 0.8 * e1[up] + sqrt(1 - 0.64) * e2[up]; e[lower.tri(e)] = t(et)[lower.tri(e)]
	Y = matrix(rpois(n * n, exp(2 + 0.5 * x + outer(a, b, "+") + e)), n, n); diag(Y) = NA
	rn = sprintf("r%02d", seq_len(n)); dimnames(Y) = list(rn, rn)
	X = array(x, c(n, n, 1), dimnames = list(rn, rn, "x"))
	f1 = suppressWarnings(ame_als(Y, Xdyad = X, family = "poisson", R = 0, verbose = FALSE))
	f0 = suppressWarnings(ame_als(Y, Xdyad = X, family = "poisson", R = 0, dcor = FALSE, verbose = FALSE))
	expect_identical(f1$eb$pair$method, "integrated pair likelihood")
	expect_null(f0$eb$pair)
	# the correlation is the pair estimate, near the generating 0.8; the post hoc moment estimate is kept beside it
	expect_equal(f1$VC[["rho"]], f1$eb$pair$rho)
	expect_gt(f1$VC[["rho"]], 0.55); expect_lt(f1$VC[["rho"]], 0.97)
	expect_true(is.finite(f1$rho_estimation$sd))
	# the cross information of the two cells of a dyad (the covariate is the same in both) raises the standard error
	expect_gt(f1$eb$se[[2]], f0$eb$se[[2]] * 1.1)
	expect_equal(f1$coefficients, f0$coefficients, tolerance = 1e-8)
	# the two engines give the same pair stage
	op = options(lame.als_engine = "r"); on.exit(options(op))
	f1r = suppressWarnings(ame_als(Y, Xdyad = X, family = "poisson", R = 0, verbose = FALSE))
	expect_equal(f1r$VC[["rho"]], f1$VC[["rho"]], tolerance = 1e-6)
	expect_equal(f1r$eb$se, f1$eb$se, tolerance = 1e-5)
})

test_that("a count network without reciprocity keeps its correlation near zero, and the pair stage is skipped where it does not apply", {
	skip_on_cran()
	d = sim_pln_net(24, seed = 14)
	f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE))
	expect_lt(abs(f$VC[["rho"]]), 0.4)
	# undirected and bipartite networks have no pair of cells
	fb = suppressWarnings(
		ame_als(
			sim_pln_net(20, m = 14, seed = 15, kind = "bipartite")$Y,
			family = "poisson",
			mode = "bipartite",
			R = 0,
			verbose = FALSE
		)
	)
	expect_null(fb$eb$pair)
	fs = suppressWarnings(
		ame_als(
			sim_pln_net(20, seed = 16, kind = "symmetric")$Y,
			family = "poisson",
			symmetric = TRUE,
			R = 0,
			verbose = FALSE
		)
	)
	expect_null(fs$eb$pair)
})

test_that("a pair stage that fails is reported and the fit keeps independent cells", {
	skip_on_cran()
	d = sim_pln_net(24, seed = 14)
	testthat::local_mocked_bindings(.ae_pln_pair_stage = function(...) stop("no pair"), .package = "lame")
	expect_warning(f <- ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE), class = "lame_als_pair")
	expect_null(f$eb$pair)
	expect_true(f$converged)
	f0 = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, dcor = FALSE, verbose = FALSE))
	expect_equal(f$coefficients, f0$coefficients, tolerance = 1e-10)
	expect_equal(f$eb$se, f0$eb$se, tolerance = 1e-10)
})

test_that("the replicates of a bootstrap form the pair stage, so that their correlation is the one of the fit", {
	skip_on_cran()
	d = sim_pln_net(24, seed = 14)
	f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, verbose = FALSE))
	orig = lame:::.ae_pln_pair_stage
	calls = 0L
	testthat::local_mocked_bindings(
		.ae_pln_pair_stage = function(...) { calls <<- calls + 1L; orig(...) },
		.package = "lame"
	)
	b = suppressWarnings(ame_als_bootstrap(f, R = 10, seed = 3, verbose = FALSE))
	expect_gte(calls, 10L)
	expect_equal(dim(b$vc_coefs), c(10L, 5L))
})

# a weak-information network of the validation against the sampler: a directed network of n actors with a latent product of
# rank R in which most cells are zero (a few nonzero counts for every loading of the factors)
sim_weak = function(n, R, mu, seed, vuv = 0.6, ve = 0.4) {
	set.seed(seed)
	S = diag(0.5, 2); ab = matrix(rnorm(2 * n), n, 2) %*% chol(S); a = ab[, 1]; b = ab[, 2]
	U = matrix(rnorm(n * R, 0, sqrt(vuv)), n, R); V = matrix(rnorm(n * R, 0, sqrt(vuv)), n, R)
	x = matrix(rnorm(n * n), n, n)
	eta = mu + 0.5 * x + outer(a, b, "+") + tcrossprod(U, V)
	E = matrix(rnorm(n * n, 0, sqrt(ve)), n, n)
	Y = matrix(rpois(n * n, exp(eta + E)), n, n); diag(Y) = NA
	rn = sprintf("r%03d", seq_len(n)); dimnames(Y) = list(rn, rn)
	list(Y = Y, X = array(x, c(n, n, 1), dimnames = list(rn, rn, "x")))
}
# the mode of a fitted object at a noise variance held at a value, cold
mode_at = function(f, ve, from = NULL) {
	st = lame:::.ae_mf_object_state(f)
	lame:::.ae_irls_run(
		st$prep,
		"poisson",
		"log",
		f$R,
		FALSE,
		200L,
		1e-6,
		lowrank_method = "als",
		linear_solver = "eigen",
		ab_prior = st$ab,
		init = if (is.null(from)) NULL else lame:::.ae_pln_mode_init(from$fit, ve),
		warm = !is.null(from),
		ve_fixed = ve,
		warm_retry = FALSE
	)
}

test_that("the mode holds a noise variance that it is given, and the search finds the maximum of the Laplace marginal likelihood", {
	skip_on_cran()
	d = sim_pln_net(30, R = 1, mu = -0.3, ve = 0.4, seed = 21)
	f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 1, dcor = FALSE, verbose = FALSE))
	m3 = mode_at(f, 0.3)
	expect_equal(m3$fit$ve, 0.3)
	expect_true(m3$stationarity$pass)
	# a warm mode at a nearby value is the cold mode there (a live column: the same evidence)
	m35 = mode_at(f, 0.35); m35w = mode_at(f, 0.35, from = m3)
	expect_equal(m35w$fit$eb_vcov$evidence$log_marginal, m35$fit$eb_vcov$evidence$log_marginal, tolerance = 1e-5)
	# the search: the reported variance is the one of the best point, and no cold mode at a value 25 percent away does better
	vs = f$eb$ve_search
	expect_equal(f$VC[["ve"]], vs$ve)
	expect_equal(f$eb$ve, vs$ve)
	expect_equal(vs$objective, max(vs$curve$objective))
	J = function(ve) { ir = mode_at(f, ve); ir$fit$eb_vcov$evidence$log_marginal + lame:::.ae_pln_ve_logprior(ve) }
	expect_lte(J(vs$ve * 0.8), vs$objective + 0.05)
	expect_lte(J(vs$ve * 1.25), vs$objective + 0.05)
})

test_that("a noise variance of a latent product is not overstated where the counts are few: the sparse rank-3 network of the validation", {
	skip_on_cran()
	d = sim_weak(50, 3, -1.2, 402)
	f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 3, dcor = FALSE, verbose = FALSE))
	expect_true(f$converged)
	# the reported ve lies within the range of the sampler's posterior on this network
	expect_gt(f$VC[["ve"]], 0.2); expect_lt(f$VC[["ve"]], 0.43)
	expect_gt(svd(f$UVPM)$d[3], 5)
})

test_that("a column of the factors that the mode collapsed is restarted by the variational stage", {
	skip_on_cran()
	d = sim_weak(40, 2, -2.5, 81)
	f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 2, dcor = FALSE, verbose = FALSE))
	expect_gt(svd(lame:::.ae_product_restore(f)$UVPM)$d[2], 1)
	# a mode whose second column is zero: the variational start restarts it from the score of the cells,
	# since the updates of a zero column keep it at zero
	st = lame:::.ae_mf_object_state(f)
	ir = mode_at(f, f$VC[["ve"]])
	ir$fit$U[, 2] = 0; ir$fit[["V"]][, 2] = 0
	cx = lame:::.ae_mf_context(st$prep, st$ab, 2L, "log")
	S = lame:::.ae_mf_init(cx, ir$fit, ir$W)
	nu = sqrt(colSums(S$U^2)); nv = sqrt(colSums(S[["V"]]^2))
	expect_gt(nu[2], 0.1 * nu[1]); expect_gt(nv[2], 0.1 * nv[1])
})

test_that("a rank above the generating rank converges to the maximum of the marginal likelihood, with the noise variance of the model", {
	skip_on_cran()
	for (seed in c(11, 13)) {
		d = sim_pln_net(30, R = 0, mu = 0, ve = 0.5, seed = seed)
		f = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 2, dcor = FALSE, verbose = FALSE))
		expect_true(f$converged)
		vs = f$eb$ve_search
		# the mode at the reported variance is one that solved its equations
		expect_true(vs$curve$solved[which.max(vs$curve$objective)])
		expect_gt(f$VC[["ve"]], 0.3); expect_lt(f$VC[["ve"]], 0.75)
	}
})

# the priors the estimator shares with the samplers: the scale of the prior of the noise variance (prior$s20) and the degrees
# of freedom of the additive effects of a bipartite panel (the number of row actors, as lame() counts them)
test_that("the scale of the prior of the noise variance is the sampler's prior$s20", {
	skip_on_cran()
	d = sim_pln_net(30, R = 0, mu = 0, ve = 0.5, seed = 31)
	f1 = suppressWarnings(ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, dcor = FALSE, verbose = FALSE))
	f2 = suppressWarnings(
		ame_als(d$Y, Xdyad = d$X, family = "poisson", R = 0, dcor = FALSE, verbose = FALSE, prior = list(s20 = 4))
	)
	expect_equal(f1$eb$ve_prior[["s20"]], 1)
	expect_equal(f2$eb$ve_prior[["s20"]], 4)
	# the larger scale of the prior raises the variance a little
	expect_gt(f2$VC[["ve"]], f1$VC[["ve"]])
	expect_lt(f2$VC[["ve"]] / f1$VC[["ve"]], 1.08)
	# the EM image and the log density read the prior they are given
	expect_equal(lame:::.ae_pln_ve_image(10, 100, 1, c(nu0 = 1, s20 = 4)), (4 + 10) / (100 - 1))
	expect_equal(lame:::.ae_pln_ve_logprior(0.5, c(nu0 = 1, s20 = 4)), 0.5 * log(0.5) - 0.5 * 4 / 0.5)
	# a refit keeps the prior of its fit
	r = suppressWarnings(ame_als_refit(f2, Y_new = d$Y, verbose = FALSE))
	expect_equal(r$eb$ve_prior[["s20"]], 4)
	expect_equal(r$VC[["ve"]], f2$VC[["ve"]], tolerance = 1e-3)
	# ame() and lame() hand the prior over to the static route
	g = suppressWarnings(
		ame(
			d$Y,
			Xdyad = d$X,
			family = "poisson",
			R = 0,
			method = "als",
			dcor = FALSE,
			verbose = FALSE,
			prior = list(s20 = 4)
		)
	)
	expect_equal(g$eb$ve_prior[["s20"]], 4)
	# the keys the estimator does not read, and a scale for a family without a noise variance, are reported
	expect_warning(
		ame_als(
			d$Y,
			Xdyad = d$X,
			family = "poisson",
			R = 0,
			dcor = FALSE,
			verbose = FALSE,
			prior = list(Sab0 = diag(2))
		),
		"no counterpart"
	)
	expect_warning(
		ame_als(1 * (d$Y > 0), family = "binary", R = 0, dcor = FALSE, verbose = FALSE, prior = list(s20 = 2)),
		"Poisson and normal estimators only"
	)
	expect_error(ame_als(d$Y, family = "poisson", R = 0, verbose = FALSE, prior = list(s20 = -1)), "positive scalar")
})

test_that("the additive effects of a bipartite count panel have the degrees of freedom lame() gives them", {
	for (dims in list(c(30L, 50L), c(50L, 30L), c(40L, 40L))) {
		Y = array(1, c(dims, 1L))
		expect_equal(lame:::.ae_eb_hyperprior_pois(Y)$eta0, round(4 + 3 * dims[1L] / 100))
	}
})
