skip_on_cran()

# a small network from the model: effects, a rank-2 product, one dyadic and one row covariate
sim_se_net = function(family, n = 16, seed = 1, kind = "directed") {
	set.seed(seed)
	nc = if (kind == "bipartite") 12 else n
	rn = paste0("r", 1:n); cn = if (kind == "bipartite") paste0("c", 1:nc) else rn
	a = rnorm(n, 0, 0.5); b = if (kind == "symmetric") a else rnorm(nc, 0, 0.5)
	u = matrix(rnorm(2 * n), n); v = if (kind == "symmetric") u else matrix(rnorm(2 * nc), nc)
	x = matrix(rnorm(n * nc), n, nc); e = matrix(rnorm(n * nc), n, nc)
	if (kind == "symmetric") { x = (x + t(x)) / 2; e = (e + t(e)) / sqrt(2) }
	w = rnorm(n)
	eta = -0.3 + 0.5 * x + 0.3 * outer(w, if (kind == "symmetric") w else rep(0, nc), "+") + outer(a, b, "+") + tcrossprod(u, v)
	y = switch(family, normal = eta + e, binary = 1 * (eta + e > 0), poisson = matrix(rpois(n * nc, exp(pmin(eta + 0.5 * e, 4))), n))
	if (kind != "bipartite") diag(y) = NA
	dimnames(y) = list(rn, cn)
	list(y = y, x = array(x, c(n, nc, 1), dimnames = list(rn, cn, "x")), xr = matrix(w, n, 1, dimnames = list(rn, "w")))
}
finite_pos = function(x) length(x) > 0 && all(is.finite(x)) && all(x > 0)

test_that("a directed fit reports sds of the effects, positions and product cells", {
	d = sim_se_net("binary", seed = 11)
	f = ame_als(d$y, Xdyad = d$x, Xrow = d$xr, R = 2, family = "binary", verbose = FALSE)
	expect_true(finite_pos(f$a_sd))
	expect_true(finite_pos(f$b_sd))
	expect_identical(names(f$a_sd), names(f$a))
	if (!lame:::.ae_product_is_withheld(f)) {
		expect_identical(dimnames(f$U_sd), dimnames(f$U))
		expect_identical(dimnames(f$UVPM_sd), dimnames(f$UVPM))
		expect_true(finite_pos(f$U_sd))
		expect_true(finite_pos(f$V_sd))
		expect_true(finite_pos(f$UVPM_sd[row(f$UVPM_sd) != col(f$UVPM_sd)]))
		lp = suppressMessages(latent_positions(f))
		expect_true(finite_pos(lp$posterior_sd))
	} else {
		# a withheld product withholds its sds with it
		expect_true(all(is.na(f$UVPM_sd)) && all(is.na(f$U_sd)))
	}
})

test_that("variance components carry a standard error and an interval inside the parameter space", {
	d = sim_se_net("normal", seed = 12)
	f = ame_als(d$y, Xdyad = d$x, Xrow = d$xr, R = 2, family = "normal", verbose = FALSE)
	for (k in c("va", "cab", "vb", "rho", "ve")) {
		expect_true(finite_pos(f$VC_se[[k]]), info = k)
		expect_true(f$VC_ci[k, 1] < f$VC[[k]] && f$VC[[k]] < f$VC_ci[k, 2], info = k)
	}
	expect_true(all(f$VC_ci[c("va", "vb", "ve"), 1] > 0))
	expect_true(abs(f$VC_ci["rho", 1]) < 1 && abs(f$VC_ci["rho", 2]) < 1)
	bound = sqrt(f$VC[["va"]] * f$VC[["vb"]])
	expect_true(all(abs(f$VC_ci["cab", ]) <= bound))
})

test_that("a symmetric fit reports eigenvector and eigenvalue sds and no sender-receiver terms", {
	d = sim_se_net("binary", n = 18, seed = 13, kind = "symmetric")
	f = ame_als(d$y, Xdyad = d$x, Xrow = d$xr, R = 1, family = "binary", symmetric = TRUE, verbose = FALSE)
	expect_true(finite_pos(f$a_sd))
	expect_identical(f$a_sd, f$b_sd)
	expect_true(finite_pos(f$VC_se[["va"]]))
	expect_true(all(is.na(f$VC_se[c("cab", "vb", "rho", "ve")])))
	if (!lame:::.ae_product_is_withheld(f)) {
		expect_true(finite_pos(f$U_sd))
		expect_true(finite_pos(f$L_sd))
		expect_identical(f$ULUPM_sd, f$UVPM_sd)
	}
})

test_that("bipartite and panel fits size their sds by side, without a product at rank zero", {
	d = sim_se_net("normal", seed = 14, kind = "bipartite")
	f = ame_als(d$y, Xdyad = d$x, Xrow = d$xr, R = 2, family = "normal", mode = "bipartite", verbose = FALSE)
	expect_length(f$a_sd, nrow(d$y))
	expect_length(f$b_sd, ncol(d$y))
	expect_identical(dim(f$V_sd), dim(f[["V"]]))
	expect_true(all(is.na(f$VC_se[c("cab", "rho")])))
	p = lapply(1:3, function(t) sim_se_net("poisson", seed = 20 + t))
	g = lame_als(lapply(p, `[[`, "y"), Xdyad = lapply(p, `[[`, "x"), R = 0, family = "poisson", verbose = FALSE)
	expect_true(finite_pos(g$a_sd))
	expect_null(g$U_sd)
	expect_true(finite_pos(g$VC_se[c("va", "vb", "ve")]))
})

test_that("the effect sds are the delta method under the joint laplace covariance", {
	d = sim_se_net("normal", seed = 15)
	f = ame_als(d$y, R = 0, family = "normal", verbose = FALSE)
	st = lame:::.ae_eb_laplace_state(f)
	# with no product and no covariates the reported a is the stacked effect itself
	S = st$jp$S
	expect_equal(unname(f$a_sd), sqrt(diag(S)[st$jp$ctx$i_s[seq_len(nrow(d$y))]]), tolerance = 0.05)
})

test_that("summary and print show the standard errors of the variance components", {
	d = sim_se_net("binary", seed = 16)
	f = ame_als(d$y, Xdyad = d$x, R = 1, family = "binary", verbose = FALSE)
	s = summary(f)
	expect_true("Std.Error" %in% colnames(s$VC_table))
	txt = paste(capture.output(print(s), type = "message"), capture.output(print(s)), collapse = " ")
	txt = paste(txt, paste(capture.output(print(f)), collapse = " "))
	expect_false(grepl("Point estimate only|No standard errors", txt))
	expect_true(grepl("a_sd", txt))
})

test_that("the information inverse equilibrates before it gives up a direction", {
	H = diag(c(1e10, 1, 1e-6))
	H[1, 2] = H[2, 1] = 10
	S = lame:::.ae_eb_inverse(H)
	expect_null(attr(S, "ae_singular"))
	expect_equal(S %*% H, diag(3), tolerance = 1e-8)
	Hs = matrix(1, 3, 3)
	expect_type(attr(lame:::.ae_eb_inverse(Hs), "ae_singular"), "character")
})

test_that("a normal fit adds the noise a weak factor column takes up to its noise variance", {
	excess = lame:::.ae_gauss_bulk_excess
	set.seed(12); n = 50
	uv = tcrossprod(rnorm(n), rnorm(n)); uv = uv / sqrt(sum(uv^2))
	E = array(rnorm(n * n * 2), c(n, n, 2))
	one = function(P, y) list(cx = list(nr = n, nc = n, obs = is.finite(y), y = y), ev = list(eta = array(P, dim(y)), Pm = P))
	# a weak column fitted at the data's own top direction takes up noise; a strong one next to none
	weak = sqrt(1.5 * n) * uv
	y = E + array(weak, dim(E))
	sv = svd(rowMeans(y, dims = 2L) * sqrt(2))
	P = sv$d[1] * tcrossprod(sv$u[, 1], sv$v[, 1]) / sqrt(2)
	a = one(P, y)
	ex_weak = excess(a$cx, a$ev, 1L)
	expect_gt(ex_weak, 0.2 * n)
	strong = 10 * n * uv
	y2 = E + array(strong, dim(E))
	a2 = one(strong, y2)
	ex_strong = excess(a2$cx, a2$ev, 1L)
	expect_lt(abs(ex_strong), 1)
	# the noise variance moves to the EM image with the excess (here far from the start), and stops once it settles
	cx = c(a$cx, list(R = 1L, n_obs = sum(a$cx$obs), nu = 1, ab = list(family = "normal", ve_prior = c(nu0 = 1, s20 = 1))))
	ee = array(0.9, dim(y))
	res = list(cx = cx, ev = c(a$ev, list(wk = list(ee = ee))), S = list(ve = 1), converged = TRUE)
	calls = 0
	solve = function(S1, stat) { calls <<- calls + 1; r = res; r$S = S1; r }
	out = lame:::.ae_gauss_bulk_ve(res, solve)
	expect_equal(out$S$ve, (1 + sum(ee) + ex_weak) / (cx$n_obs - 1))
	expect_identical(calls, 1)
	# nothing moves without a product
	none = res; none$ev$Pm = 0 * P; none$ev$eta = array(0, dim(y))
	expect_identical(lame:::.ae_gauss_bulk_ve(none, solve)$S$ve, 1)
})

# a network whose product (columns of scale `sig`) is below the signal threshold
sim_weak_net = function(family, n, nc, seed, kind, sig = 0.35) {
	set.seed(seed)
	rn = sprintf("a%02d", 1:n); cn = if (kind == "bipartite") sprintf("c%02d", 1:nc) else rn
	a = rnorm(n, 0, 0.5); b = if (kind == "symmetric") a else rnorm(nc, 0, 0.5)
	u = matrix(rnorm(2 * n), n) * sig; v = if (kind == "symmetric") u else matrix(rnorm(2 * nc), nc) * sig
	x = matrix(rnorm(n * nc), n, nc); e = matrix(rnorm(n * nc), n, nc)
	if (kind == "symmetric") { x = (x + t(x)) / 2; e = (e + t(e)) / sqrt(2) }
	eta = -0.3 + 0.5 * x + outer(a, b, "+") + tcrossprod(u, v)
	y = switch(family, normal = eta + e, poisson = matrix(rpois(n * nc, exp(eta + 0.5 * e)), n))
	if (kind == "symmetric") y[lower.tri(y)] = t(y)[lower.tri(y)]
	if (kind != "bipartite") diag(y) = NA
	dimnames(y) = list(rn, cn)
	list(y = y, x = array(x, c(n, nc, 1), dimnames = list(rn, cn, "x")))
}

test_that("a weak undirected fit reports eigenvalue intervals on both signs and the uninformed direction sds", {
	d = sim_weak_net("normal", 30, 30, 4, "symmetric")
	f = ame_als(d$y, Xdyad = d$x, R = 2, family = "normal", symmetric = TRUE, verbose = FALSE)
	info = f$eb$variational$information
	expect_lt(lame:::.ae_product_statistic(f$family, info), info[["threshold"]])
	ok = is.finite(f$L) & f$L != 0
	expect_true(all(f$L_ci[ok, 1] < 0 & f$L_ci[ok, 2] > 0))
	expect_equal(-f$L_ci[ok, 1], f$L_ci[ok, 2])
	n = nrow(f$U); ec = exp(lgamma((n - 1) / 2) - lgamma(n / 2)) / sqrt(pi)
	floor_sd = sqrt(f$U[, ok]^2 * (1 - 2 * ec) + 1 / n)
	expect_true(all(f$U_sd[, ok] >= floor_sd - 1e-12))
})

test_that("a count product held at zero does not charge its prior width to the coefficients", {
	d = sim_weak_net("poisson", 30, 24, 1, "bipartite")
	f = ame_als(d$y, Xdyad = d$x, R = 2, family = "poisson", mode = "bipartite", verbose = FALSE)
	f0 = ame_als(d$y, Xdyad = d$x, R = 0, family = "poisson", mode = "bipartite", verbose = FALSE)
	expect_true(any(colSums(f$U^2) == 0 | colSums(f[["V"]]^2) == 0))
	se = function(g) sqrt(diag(vcov(g)))
	# the rank-0 fit of the same counts is the reference: the held columns add at most their
	# below-threshold cell variance
	expect_lt(max(se(f) / se(f0)), 1.3)
	expect_true(finite_pos(f$UVPM_sd))
})

test_that("a weak directed fit adds to each row the spread it takes from the other side", {
	d = sim_weak_net("normal", 40, 40, 2, "directed")
	f = ame_als(d$y, Xdyad = d$x, R = 2, family = "normal", verbose = FALSE)
	info = f$eb$variational$information
	expect_lt(lame:::.ae_product_statistic(f$family, info), info[["threshold"]])
	# the variational row sds condition on the other side at its mean: the reported sds exceed them
	fp = f$eb$factor_posterior
	row_sd = function(A) t(sqrt(apply(A, 3, diag)))
	expect_gt(mean(f$U_sd), 1.05 * mean(row_sd(fp$Sig)))
	expect_gt(mean(f$V_sd), 1.05 * mean(row_sd(fp$Om)))
	expect_true(finite_pos(f$U_sd) && finite_pos(f$V_sd))
})

test_that("a weak undirected dynamic fit reports eigenvalue intervals on both signs", {
	withr::local_options(lame.eb_withhold_product = FALSE)
	set.seed(5)
	n = 30; rn = sprintf("a%02d", 1:n)
	fac = cbind(rnorm(n, 0, 0.35), rnorm(n, 0, 0.21))
	P = tcrossprod(fac[, 1]) - tcrossprod(fac[, 2])
	a = rnorm(n, 0, 0.5)
	one = function() {
		x = matrix(rnorm(n * n), n); x = (x + t(x)) / 2
		e = matrix(rnorm(n * n), n); e = (e + t(e)) / sqrt(2)
		y = -0.3 + 0.5 * x + outer(a, a, "+") + P + e
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		list(y = y, x = array(x, c(n, n, 1), dimnames = list(rn, rn, "x")))
	}
	d = lapply(1:4, function(t) one())
	f = lame(lapply(d, `[[`, "y"), Xdyad = lapply(d, `[[`, "x"), R = 2, family = "normal", symmetric = TRUE,
		method = "als", dynamic_ab = TRUE, verbose = FALSE)
	info = f$eb$information
	expect_lt(lame:::.ae_product_statistic(f$family, info), info[["threshold"]])
	expect_equal(-f$L_ci[, 1], f$L_ci[, 2])
	expect_true(all(f$L_ci[, 2] > abs(f$L)))
	expect_true(finite_pos(f$U_sd))
})

test_that("an undirected column inside the sign margin takes the two-sign interval", {
	set.seed(3)
	n = 20
	U = qr.Q(qr(matrix(rnorm(n), n)))
	out = list(symmetric = TRUE, family = "binary", U = U, L = 4, L_ci = cbind(lower = 3, upper = 6), L_sd = 0.8, U_sd = matrix(0.05, n, 1))
	info = c(signal = 20, loadings = n, signal_per_loading = 1, threshold = 0.7)
	Sig = array(0.05, c(1, 1, n))
	# 1 / 0.7 = 1.4 times the threshold: determined at margin 1, undetermined at the dynamic glm margin
	same = lame:::.ae_eb_weak_symmetric(out, info, Sig, centred = FALSE)
	expect_identical(same$L_ci, out$L_ci)
	two = lame:::.ae_eb_weak_symmetric(out, info, Sig, centred = FALSE, margin = lame:::.lde_glm_sign_margin)
	expect_equal(-unname(two$L_ci[1, 1]), unname(two$L_ci[1, 2]))
	expect_gt(two$L_ci[1, 2], 4)
	expect_true(all(two$U_sd > out$U_sd))
})

test_that("an undirected column held at zero reports its zero positions with the uniform-direction sd", {
	d = sim_weak_net("poisson", 25, 25, 3, "symmetric")
	f = ame_als(d$y, Xdyad = d$x, R = 2, family = "poisson", symmetric = TRUE, verbose = FALSE)
	zero = f$L == 0
	expect_true(any(zero))
	expect_true(all(f$U[, zero] == 0))
	expect_equal(unname(f$U_sd[, zero]), matrix(1 / sqrt(nrow(f$U)), nrow(f$U), sum(zero)))
	lp = latent_positions(f)
	expect_true(all(is.finite(lp$value)) && finite_pos(lp$posterior_sd))
})

test_that("the product cells of a weak undirected fit carry its undetermined direction and sign", {
	d = sim_weak_net("normal", 30, 30, 4, "symmetric")
	f = ame_als(d$y, Xdyad = d$x, R = 2, family = "normal", symmetric = TRUE, verbose = FALSE)
	info = f$eb$variational$information
	expect_lt(lame:::.ae_product_statistic(f$family, info), info[["threshold"]])
	# the columns below the threshold by their share of the statistic
	st = lame:::.ae_product_statistic(f$family, info) * 2 * f$L^2 / sum(f$L^2)
	weak = which(st < info[["threshold"]])
	expect_true(length(weak) > 0)
	n = nrow(f$U); k = (n - 1) * (n + 1)
	cf = Reduce(`+`, lapply(weak, function(r) (f$L[r] * tcrossprod(f$U[, r]))^2 + f$L[r]^2 / k))
	off = row(cf) != col(cf)
	expect_true(all(f$UVPM_sd[off] >= sqrt(cf[off]) - 1e-12))
	expect_true(finite_pos(f$UVPM_sd[off]))
})

test_that("the uncentred weak row variances (dynamic report) are the gauge-projected diagonal", {
	set.seed(7)
	nr = 5; nc = 4; R = 2
	Fu = matrix(rnorm(nr * R), nr); Fv = matrix(rnorm(nc * R), nc)
	spd = function() { A = matrix(rnorm(R * R), R); crossprod(A) + diag(R) }
	Sig = array(vapply(1:nr, function(i) spd(), matrix(0, R, R)), c(R, R, nr))
	Om = array(vapply(1:nc, function(j) spd(), matrix(0, R, R)), c(R, R, nc))
	W = matrix(runif(nr * nc), nr, nc)
	# dense reference: rows with the move from the other side, block diagonal, gauge projected
	D = matrix(0, (nr + nc) * R, (nr + nc) * R)
	for (i in 1:nr) { M = Reduce(`+`, lapply(1:nc, function(j) W[i, j] * Om[, , j])); ix = (0:(R - 1)) * nr + i; D[ix, ix] = Sig[, , i] + Sig[, , i] %*% M %*% Sig[, , i] }
	for (j in 1:nc) { M = Reduce(`+`, lapply(1:nr, function(i) W[i, j] * Sig[, , i])); ix = R * nr + (0:(R - 1)) * nc + j; D[ix, ix] = Om[, , j] + Om[, , j] %*% M %*% Om[, , j] }
	Q = qr.Q(qr(lame:::.ae_eb_gauge_dirs(Fu, Fv, rep(1, R), FALSE)))
	P = diag(nrow(D)) - tcrossprod(Q)
	expect_equal(lame:::.ae_eb_weak_row_var(Fu, Fv, Sig, Om, W, centre = FALSE), diag(P %*% D %*% P))
})
