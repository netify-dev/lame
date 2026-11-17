skip_on_cran()

library(lame)

# empirical-Bayes ALS fits with a latent product report the intercept, effects and positions as the samplers
# store them: positions centred over the actors, the means of the product in the effects and the intercept

# a network whose positions have a mean, so that the report moves the intercept and the effects
rf_net = function(kind = "directed", fam = "normal", n = 30, R = 2, seed = 1, uv_mean = 0.5, uv_sd = 0.9) {
	set.seed(seed)
	sym = kind == "symmetric"; bip = kind == "bipartite"
	nc = if (bip) round(0.7 * n) else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("c%02d", seq_len(nc)) else rn
	a = rnorm(n, 0, 0.5); b = if (sym) a else rnorm(nc, 0, 0.5)
	U = matrix(rnorm(n * R, uv_mean, uv_sd), n)
	V = if (sym) U else matrix(rnorm(nc * R, -uv_mean, uv_sd), nc)
	x = matrix(rnorm(n * nc), n, nc); e = matrix(rnorm(n * nc), n, nc)
	if (sym) { x[lower.tri(x)] = t(x)[lower.tri(x)]; e[lower.tri(e)] = t(e)[lower.tri(e)] }
	z = 0.3 + 0.5 * x + outer(a, b, "+") + tcrossprod(U, V) + e
	y = if (fam == "binary") 1 * (z > 0) else z
	if (!bip) diag(y) = NA
	dimnames(y) = list(rn, cn)
	list(Y = y, X = array(x, c(n, nc, 1), dimnames = list(rn, cn, "x")), kind = kind, fam = fam, R = R)
}
rf_fit = function(d, ...) suppressWarnings(ame_als(
	d$Y, Xdyad = d$X, R = d$R, family = d$fam, symmetric = d$kind == "symmetric",
	mode = if (d$kind == "bipartite") "bipartite" else "unipartite", verbose = FALSE, ...
))
# the part of the linear predictor the report rearranges
rf_part = function(f) unname(coef(f)[[1]] + outer(f$a, f$b, "+") + f$UVPM)
rf_grid = expand.grid(
	kind = c("directed", "symmetric", "bipartite"), fam = c("normal", "binary"), R = 1:2,
	stringsAsFactors = FALSE
)
rf_fits = lapply(seq_len(nrow(rf_grid)), function(i) {
	rf_fit(rf_net(rf_grid$kind[i], rf_grid$fam[i], n = if (rf_grid$fam[i] == "binary") 60 else 30, R = rf_grid$R[i], seed = i))
})
rf_label = paste(rf_grid$kind, rf_grid$fam, rf_grid$R)

test_that("the report centres the positions and leaves the linear predictor where the fit put it", {
	for (i in seq_along(rf_fits)) {
		f = rf_fits[[i]]; lab = rf_label[i]
		rf = f$eb$report_fold
		expect_false(is.null(rf), label = lab)
		expect_identical(rf$view, "report", label = lab)
		expect_lt(max(abs(colMeans(f$U))), 1e-10, label = lab)
		expect_lt(max(abs(colMeans(f[["V"]]))), 1e-10, label = lab)
		expect_gt(abs(rf$const), 1e-3, label = lab)
		s = lame:::.ae_product_restore(f)
		expect_equal(rf_part(f), rf_part(s), tolerance = 1e-10, label = lab)
		expect_equal(coef(f)[[1]] - f$eb$population$intercept, rf$const, tolerance = 1e-10, label = lab)
		expect_equal(coef(f)[-1], coef(s)[-1], label = lab)
		if (rf_grid$kind[i] == "symmetric") {
			expect_equal(unname(crossprod(f$U)), diag(1, rf_grid$R[i]), tolerance = 1e-8, label = lab)
			expect_equal(unname(f$U %*% (f$L * t(f$U))), unname(f$UVPM), tolerance = 1e-8, label = lab)
			expect_equal(f$a, f$b, label = lab)
		} else {
			expect_equal(unname(tcrossprod(f$U, f[["V"]])), unname(f$UVPM), tolerance = 1e-8, label = lab)
		}
	}
})

test_that("a normal fit's fitted values are the reported intercept, slope, effects and product", {
	for (i in which(rf_grid$fam == "normal")) {
		f = rf_fits[[i]]
		d = rf_net(rf_grid$kind[i], "normal", n = 30, R = rf_grid$R[i], seed = i)
		eta = rf_part(f) + coef(f)[[2]] * d$X[, , 1]
		ok = !is.na(d$Y)
		expect_equal(unname(fitted(f))[ok], unname(eta)[ok], tolerance = 1e-8, label = rf_label[i])
	}
})

test_that("the fitted state comes back exactly, once or twice, and the report after it", {
	for (i in seq_along(rf_fits)) {
		f = rf_fits[[i]]; lab = rf_label[i]
		rf = f$eb$report_fold
		s = lame:::.ae_product_restore(f)
		expect_identical(s$eb$report_fold$view, "state", label = lab)
		expect_identical(s[names(rf$state)], rf$state, label = lab)
		expect_identical(lame:::.ae_product_restore(s), s, label = lab)
		expect_identical(lame:::.ae_report_view(s), f, label = lab)
		expect_identical(lame:::.ae_product_restore(f, state = FALSE), f, label = lab)
	}
})

test_that("the intercept's variance carries the grand mean of the product, the slopes' is untouched", {
	for (i in seq_along(rf_fits)) {
		f = rf_fits[[i]]; lab = rf_label[i]
		vc = vcov(f); pop = f$eb$population$vcov
		expect_equal(dim(vc), dim(pop), label = lab)
		expect_equal(unname(vc[-1, -1]), unname(pop[-1, -1]), label = lab)
		expect_true(all(is.finite(vc)) && all(diag(vc) > 0), label = lab)
		expect_equal(vc, t(vc), label = lab)
		ev = eigen(vc, symmetric = TRUE, only.values = TRUE)$values
		expect_gte(min(ev), -1e-10 * max(ev), label = lab)
	}
})

test_that("the gradient of the product's grand mean matches a finite difference", {
	set.seed(11)
	n = 9; nc = 6; R = 2; h = 1e-6
	fd = function(g, x) vapply(seq_along(x), function(k) {
		up = x; dn = x; up[k] = up[k] + h; dn[k] = dn[k] - h
		(g(up) - g(dn)) / (2 * h)
	}, 0)
	U = matrix(rnorm(n * R, 0.4), n); V = matrix(rnorm(nc * R, -0.3), nc)
	g_dir = function(x) mean(tcrossprod(matrix(x[seq_len(n * R)], n), matrix(x[-seq_len(n * R)], nc)))
	expect_equal(
		lame:::.ae_eb_fold_grad(list(kind = "bipartite"), list(U = U, V = V)), fd(g_dir, c(U, V)),
		tolerance = 1e-7
	)
	# symmetric coordinates are the scaled factors, a sign per column
	Q = qr.Q(qr(matrix(rnorm(n * R, 0.4), n))); L = c(1.7, -0.6)
	Fm = sweep(Q, 2L, sqrt(abs(L)), "*"); sg = sign(L)
	g_sym = function(x) { Fx = matrix(x, n); mean(Fx %*% (sg * t(Fx))) }
	expect_equal(
		lame:::.ae_eb_fold_grad(list(kind = "symmetric"), list(U = Q, L = L)), fd(g_sym, c(Fm)),
		tolerance = 1e-7
	)
})

test_that("a weak product is reported and restores to the same fitted state however often it is asked for", {
	w = rf_fit(rf_net("directed", "binary", n = 40, R = 1, seed = 3, uv_mean = 0, uv_sd = 0.14))
	expect_false(lame:::.ae_product_is_withheld(w))
	expect_true(all(is.finite(w$U)) && all(is.finite(w$UVPM)))
	s = lame:::.ae_product_restore(w)
	expect_identical(lame:::.ae_product_restore(s), s)
	r = lame:::.ae_product_restore(w, state = FALSE)
	expect_identical(r[c("coefficients", "a", "b", "U", "V", "UVPM")], w[c("coefficients", "a", "b", "U", "V", "UVPM")])
})

test_that("bootstrap point estimates are the reported ones", {
	for (i in which(rf_grid$fam == "normal" & rf_grid$R == 2 & rf_grid$kind != "bipartite")) {
		f = rf_fits[[i]]
		b = suppressWarnings(ame_als_bootstrap(f, R = 10, seed = 1, verbose = FALSE))
		expect_equal(b$point_est, coef(f), label = rf_label[i])
		expect_equal(b$point_a, f$a, label = rf_label[i])
		expect_equal(b$point_U, f$U, label = rf_label[i])
	}
})

test_that("a fit without a product has nothing to rearrange", {
	d0 = rf_net("directed", "normal", n = 30, R = 2, seed = 1)
	d0$R = 0
	f0 = rf_fit(d0)
	expect_null(f0$eb$report_fold)
	expect_equal(coef(f0)[[1]], f0$eb$population$intercept)
})
