skip_on_cran()

# a symmetric network with node covariates x and z (truth 0.5 and 0.3) and a dyadic covariate
sym_node_data = function(n = 40, seed = 11) {
	set.seed(seed)
	x = rnorm(n); z = rnorm(n); a0 = rnorm(n, 0, 0.3)
	Xd = matrix(rnorm(n * n), n); Xd = (Xd + t(Xd)) / sqrt(2)
	e = matrix(rnorm(n * n), n); e = (e + t(e)) / sqrt(2)
	Y = 0.2 + outer(a0 + 0.5 * x + 0.3 * z, a0 + 0.5 * x + 0.3 * z, "+") + 0.5 * Xd + e
	diag(Y) = NA
	nm = sprintf("a%02d", seq_len(n)); dimnames(Y) = list(nm, nm)
	list(
		Y = Y,
		Xd = array(Xd, c(n, n, 1), list(nm, nm, "xd")),
		X = matrix(x, n, 1, dimnames = list(nm, "x")),
		Z = matrix(z, n, 1, dimnames = list(nm, "z"))
	)
}
sym_fit = function(d, ...) ame_als(d$Y, R = 0, symmetric = TRUE, verbose = FALSE, ...)

test_that("a covariate given in Xrow and Xcol enters once", {
	d = sym_node_data()
	f = sym_fit(d, Xdyad = d$Xd, Xrow = d$X, Xcol = d$X)
	g = sym_fit(d, Xdyad = d$Xd, Xrow = d$X)
	expect_setequal(names(coef(f)), c("intercept", "xd_dyad", "x_node"))
	expect_true(all(is.finite(sqrt(diag(vcov(f))))))
	expect_equal(coef(f), coef(g))
	expect_equal(vcov(f), vcov(g))
})

test_that("a covariate given only in Xcol is a node covariate", {
	d = sym_node_data()
	f = sym_fit(d, Xcol = d$X)
	expect_setequal(names(coef(f)), c("intercept", "x_node"))
	expect_equal(coef(f), coef(sym_fit(d, Xrow = d$X)))
})

test_that("different covariates in Xrow and Xcol are both estimated", {
	d = sym_node_data()
	f = sym_fit(d, Xrow = d$X, Xcol = d$Z)
	expect_setequal(names(coef(f)), c("intercept", "x_node", "z_node"))
	ci = confint(f)
	expect_true(ci["x_node", 1] < 0.5 && ci["x_node", 2] > 0.5)
	expect_true(ci["z_node", 1] < 0.3 && ci["z_node", 2] > 0.3)
	expect_equal(coef(f), coef(sym_fit(d, Xrow = cbind(d$X, d$Z)))[names(coef(f))])
	expect_true("z_node" %in% rownames(summary(f)$coefficients))
	expect_output(print(f), "z_node")
})

test_that("a name clash between Xrow and Xcol keeps both covariates under unique names", {
	d = sym_node_data()
	Zx = d$Z; colnames(Zx) = "x"
	f = sym_fit(d, Xrow = d$X, Xcol = Zx)
	expect_identical(sort(names(coef(f))), c("intercept", "x.1_node", "x_node"))
	Z1 = d$Z; colnames(Z1) = "Xcol1"; X1 = d$X; colnames(X1) = "Xrow1"
	expect_false(anyDuplicated(names(coef(sym_fit(d, Xrow = X1, Xcol = Z1)))) > 0)
})

test_that("unnamed node covariates take the samplers' default names", {
	d = sym_node_data()
	f = sym_fit(d, Xrow = unname(cbind(d$X, d$Z)))
	expect_setequal(names(coef(f)), c("intercept", "Xnode1_node", "Xnode2_node"))
	g = sym_fit(d, Xrow = unname(d$X), Xcol = unname(d$Z))
	expect_setequal(names(coef(g)), c("intercept", "Xnode1_node", "Xnode2_node"))
})

test_that("S3 consumers of a symmetric node fit work", {
	d = sym_node_data(n = 30)
	f = sym_fit(d, Xdyad = d$Xd, Xrow = d$X, Xcol = d$Z)
	p = predict(f)
	expect_equal(p, t(p))
	expect_equal(fitted(f), p)
	expect_equal(predict(f, newdata = d$Xd), p)
	expect_s3_class(simulate(f, nsim = 2, seed = 1), "ame_als.sim")
	r = ame_als_refit(f)
	expect_equal(coef(r), coef(f), tolerance = 1e-3)
	b = suppressWarnings(ame_als_bootstrap(f, R = 10, seed = 1, verbose = FALSE))
	expect_true(all(c("x_node", "z_node") %in% b$param_names))
})

test_that("a symmetric binary fit with the covariate in Xrow and Xcol fits", {
	d = sym_node_data(seed = 5)
	f = ame_als(1 * (d$Y > 0), Xrow = d$X, Xcol = d$X, R = 1, family = "binary", symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(f)), c("intercept", "x_node"))
	expect_true(all(is.finite(sqrt(diag(vcov(f))))))
})

test_that("symmetric panels use each distinct node covariate once", {
	d1 = sym_node_data(n = 25, seed = 1); d2 = sym_node_data(n = 25, seed = 2)
	XR = list(d1$X, d1$X); XC = list(d1$Z, d1$Z)
	f = lame_als(list(d1$Y, d2$Y), Xrow = XR, Xcol = XR, R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(f)), c("intercept", "x_node"))
	expect_equal(coef(f), coef(lame_als(list(d1$Y, d2$Y), Xrow = XR, R = 0, symmetric = TRUE, verbose = FALSE)))
	g = lame_als(list(d1$Y, d2$Y), Xrow = XR, Xcol = XC, R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(g)), c("intercept", "x_node", "z_node"))
	expect_length(predict(g), 2L)
	# a column that matches Xrow in one period only is a different covariate
	XC2 = list(d1$X, d1$Z); colnames(XC2[[2]]) = "x"
	h = lame_als(list(d1$Y, d2$Y), Xrow = XR, Xcol = XC2, R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(h)), c("intercept", "x_node", "x.1_node"))
	k = lame(list(d1$Y, d2$Y), Xrow = XR, Xcol = XC, R = 0, symmetric = TRUE, dynamic_ab = TRUE,
		method = "als", verbose = FALSE)
	expect_true(all(c("x_node", "z_node") %in% rownames(k$coef_path)))
})

test_that("both routes name the node covariates alike", {
	d = sym_node_data(n = 30)
	for (a in list(list(Xrow = d$X), list(Xcol = d$X), list(Xrow = d$X, Xcol = d$X), list(Xrow = d$X, Xcol = d$Z))) {
		f = do.call(sym_fit, c(list(d), a))
		m = do.call(ame, c(list(d$Y, R = 0, symmetric = TRUE, burn = 20, nscan = 60, verbose = FALSE,
			plot = FALSE, gof = FALSE), a))
		expect_setequal(names(coef(f)), colnames(m$BETA))
	}
})

test_that("the fast route fits Xrow = cbind(x, z) and Xrow = x, Xcol = z alike", {
	d = sym_node_data(n = 30)
	XZ = cbind(d$X, d$Z)
	expect_identical(coef(sym_fit(d, Xrow = XZ)), coef(sym_fit(d, Xrow = d$X, Xcol = d$Z)))
	expect_identical(fitted(sym_fit(d, Xrow = XZ)), fitted(sym_fit(d, Xrow = d$X, Xcol = d$Z)))
	G = function(...) lame_als(list(d$Y, d$Y), R = 0, symmetric = TRUE, verbose = FALSE, ...)
	a = G(Xrow = list(XZ, XZ)); b = G(Xrow = list(d$X, d$X), Xcol = list(d$Z, d$Z))
	expect_identical(coef(a), coef(b))
	expect_identical(fitted(a), fitted(b))
})
