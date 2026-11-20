skip_on_cran()

# a symmetric network with node covariates `growth` and `x` and one dyadic covariate
sym_names_data = function(n = 25, seed = 4) {
	set.seed(seed)
	g = rnorm(n); x = rnorm(n); nm = sprintf("a%02d", seq_len(n))
	Xd = matrix(rnorm(n * n), n); Xd = (Xd + t(Xd)) / 2
	e = matrix(rnorm(n * n), n); e = (e + t(e)) / 2
	Y = outer(0.5 * g + 0.2 * x, 0.5 * g + 0.2 * x, "+") + 0.5 * Xd + e
	diag(Y) = NA; dimnames(Y) = list(nm, nm)
	Xn = cbind(growth = g, x = x); rownames(Xn) = nm
	Xu = unname(Xn); rownames(Xu) = nm
	list(Y = Y, Xn = Xn, Xu = Xu, Xd = array(Xd, c(n, n, 1), list(nm, nm, "xd")))
}

test_that("the generated row marker alone becomes the node marker", {
	expect_identical(
		.lame_node_name(c("growth_row", "x_row", "arrow_col", "Xrow1_row", "Xcol2_col")),
		c("growth_node", "x_node", "arrow_node", "Xnode1_node", "Xnode2_node")
	)
})

test_that("symmetric sampler fits name node covariates by their own names", {
	d = sym_names_data()
	m = ame(d$Y, Xdyad = d$Xd, Xrow = d$Xn, R = 0, symmetric = TRUE, burn = 20, nscan = 100,
		verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_identical(colnames(m$BETA), c("intercept", "growth_node", "x_node", "xd_dyad"))
	# newdata keeps the node term: the same dyadic data reproduce the fit, a new slice moves it
	off = function(p) { diag(p) = NA; p }
	expect_equal(off(predict(m, newdata = d$Xd)), predict(m))
	X0 = m$X; X0[, , c("growth_row", "growth_col")] = 0
	b = colMeans(m$BETA)
	expect_gt(max(abs(.xbeta_newdata(b, d$Xd, m$X, 25, 25) - .xbeta_newdata(b, d$Xd, X0, 25, 25))), 0.1)
	mu = ame(d$Y, Xdyad = d$Xd, Xrow = d$Xu, R = 0, symmetric = TRUE, burn = 20, nscan = 100,
		verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_identical(colnames(mu$BETA), c("intercept", "Xnode1_node", "Xnode2_node", "xd_dyad"))
	expect_equal(off(predict(mu, newdata = d$Xd)), predict(mu))
	l = lame(list(d$Y, d$Y), Xdyad = list(d$Xd, d$Xd), Xrow = list(d$Xn, d$Xn), R = 0, symmetric = TRUE,
		burn = 20, nscan = 100, verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_identical(colnames(l$BETA), c("intercept", "growth_node", "x_node", "xd_dyad"))
	expect_no_error(predict(l, newdata = list(d$Xd, d$Xd)))
})

test_that("symmetric fast-route fits name node covariates as the sampler does", {
	d = sym_names_data()
	f = ame_als(d$Y, Xdyad = d$Xd, Xrow = d$Xn, R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(f)), c("intercept", "xd_dyad", "growth_node", "x_node"))
	expect_equal(predict(f, newdata = d$Xd), predict(f))
	expect_false(isTRUE(all.equal(predict(f, newdata = d$Xd * 0), predict(f))))
	fu = ame_als(d$Y, Xrow = d$Xu, R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(fu)), c("intercept", "Xnode1_node", "Xnode2_node"))
	fl = lame_als(list(d$Y, d$Y), Xrow = list(d$Xu, d$Xu), R = 0, symmetric = TRUE, verbose = FALSE)
	expect_setequal(names(coef(fl)), c("intercept", "Xnode1_node", "Xnode2_node"))
	r = ame_als_refit(f)
	expect_setequal(names(coef(r)), names(coef(f)))
})

test_that("dynamic coefficient paths of node covariates named with row or col are matched", {
	set.seed(7); n = 20; nn = sprintf("a%02d", 1:n)
	X = cbind(growth = rnorm(n), color = rnorm(n)); rownames(X) = nn
	Yl = lapply(1:3, function(t) {
		Y = outer(X[, 1], X[, 1], "+") * 0.4 + matrix(rnorm(n * n), n); Y = (Y + t(Y)) / 2
		diag(Y) = NA; dimnames(Y) = list(nn, nn); Y
	})
	fit = lame(Yl, Xrow = rep(list(X), 3), R = 0, symmetric = TRUE, dynamic_beta = TRUE, burn = 20,
		nscan = 100, odens = 5, verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_true(all(c("growth_node", "color_node") %in% dimnames(fit$BETA)[[2L]]))
	expect_no_error(detect_change_point(fit))
	expect_no_error(predict(fit, h = 1))
	expect_no_error(simulate(fit, nsim = 1))
	expect_identical(.lame_node_match(c("growth_row", "growth_col", "color_col"), c("gnodeth_node", "color_node")), c(1L, 1L, 2L))
})

test_that("the sampler fits Xrow = cbind(x, z) and Xrow = x, Xcol = z alike", {
	d = sym_names_data()
	X = d$Xn[, "growth", drop = FALSE]; Z = d$Xn[, "x", drop = FALSE]
	A = function(...) ame(d$Y, R = 0, symmetric = TRUE, burn = 20, nscan = 60, verbose = FALSE, plot = FALSE, gof = FALSE, ...)
	a = A(Xrow = d$Xn); b = A(Xrow = X, Xcol = Z)
	expect_identical(a$BETA, b$BETA)
	expect_identical(a$YPM, b$YPM)
	L = function(...) lame(list(d$Y, d$Y), R = 0, symmetric = TRUE, burn = 20, nscan = 100, verbose = FALSE, plot = FALSE, gof = FALSE, ...)
	a = L(Xrow = list(d$Xn, d$Xn)); b = L(Xrow = list(X, X), Xcol = list(Z, Z))
	expect_identical(a$BETA, b$BETA)
	expect_identical(a$YPM, b$YPM)
})
