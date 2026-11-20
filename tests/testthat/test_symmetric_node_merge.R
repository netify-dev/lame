skip_on_cran()

# how a symmetric model takes node covariates from Xrow and Xcol: by actor name, once per distinct
# covariate, under unique names, with the type checks on the argument the user passed
sym_merge_data = function(n = 20, seed = 5) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n)); x = rnorm(n); z = rnorm(n)
	e = matrix(rnorm(n * n), n)
	Y = outer(0.5 * x + 0.3 * z, 0.5 * x + 0.3 * z, "+") + (e + t(e)) / 2
	diag(Y) = NA; dimnames(Y) = list(nm, nm)
	d = matrix(rnorm(n * n), n); d = (d + t(d)) / 2
	M = function(...) { m = cbind(...); rownames(m) = nm; m }
	list(Y = Y, x = x, z = z, nm = nm, M = M, Xd = array(d, c(n, n, 1), list(nm, nm, "d")))
}
fit_ame = function(...) ame(..., symmetric = TRUE, R = 0, burn = 20, nscan = 40, odens = 1,
	verbose = FALSE, plot = FALSE, gof = FALSE)
off = function(p) { diag(p) = NA; p }

test_that("each node coefficient takes only its own design slices, whatever the names", {
	d = sym_merge_data()
	m = suppressMessages(fit_ame(d$Y, Xdyad = d$Xd, Xrow = d$M(brow = d$x, bnode = d$z)))
	expect_identical(colnames(m$BETA), c("intercept", "brow_node", "bnode_node", "d_dyad"))
	expect_equal(off(predict(m, newdata = d$Xd)), predict(m))
})

test_that("Xcol of a symmetric model is checked and aligned by actor name", {
	d = sym_merge_data()
	X = d$M(x = d$x)
	Zbad = d$M(z = d$z); rownames(Zbad) = sprintf("b%02d", seq_along(d$nm))
	expect_error(suppressMessages(fit_ame(d$Y, Xrow = X, Xcol = Zbad)), "carries actor names that do not match")
	expect_error(
		suppressMessages(lame(list(d$Y, d$Y), Xrow = list(X, X), Xcol = list(Zbad, Zbad), symmetric = TRUE,
			R = 0, burn = 10, nscan = 20, odens = 1, verbose = FALSE, plot = FALSE, gof = FALSE)),
		"carries actor names that do not match"
	)
	Zrev = d$M(z = d$z)[rev(d$nm), , drop = FALSE]
	a = suppressMessages(fit_ame(d$Y, Xrow = d$M(x = d$x, z = d$z)))
	b = suppressMessages(fit_ame(d$Y, Xrow = X, Xcol = Zrev))
	expect_identical(a$BETA, b$BETA)
})

test_that("an Xcol column equal to an Xrow column up to rounding is the same covariate", {
	d = sym_merge_data()
	X = d$M(x = d$x)
	expect_false(.lame_symmetric_node_covariates(list(X), list(X * (1 + 1e-15)))$changed)
	expect_true(.lame_symmetric_node_covariates(list(X), list(X * (1 + 1e-6)))$changed)
	a = ame_als(d$Y, Xrow = X, symmetric = TRUE, R = 0, verbose = FALSE)
	b = ame_als(d$Y, Xrow = X, Xcol = d$M(x = d$x * (1 + 1e-15)), symmetric = TRUE, R = 0, verbose = FALSE)
	expect_identical(coef(b), coef(a))
})

test_that("a non-numeric Xcol is reported as Xcol, without a coercion warning", {
	d = sym_merge_data()
	df = data.frame(g = factor(rep(c("u", "v"), length.out = length(d$nm))), row.names = d$nm)
	first = function(e) tryCatch(
		withCallingHandlers(e, warning = function(w) stop("warning: ", conditionMessage(w))),
		error = conditionMessage
	)
	expect_match(first(fit_ame(d$Y, Xrow = d$M(x = d$x), Xcol = df)), "`Xcol` must be numeric")
	expect_match(first(ame_als(d$Y, Xrow = d$M(x = d$x), Xcol = df, symmetric = TRUE, R = 0, verbose = FALSE)), "`Xcol` must be numeric")
})

test_that("node coefficient names are unique when Xrow names collide", {
	d = sym_merge_data()
	m = fit_ame(d$Y, Xrow = d$M(x = d$x, x_row = d$z))
	expect_identical(colnames(m$BETA), c("intercept", "x_node", "x.1_node"))
	f = ame_als(d$Y, Xrow = d$M(Xrow1 = d$x, Xcol1 = d$z), symmetric = TRUE, R = 0, verbose = FALSE)
	expect_setequal(names(coef(f)), c("intercept", "Xnode1_node", "Xnode1.1_node"))
	l = lame(list(d$Y, d$Y), Xrow = rep(list(d$M(x = d$x, x_row = d$z)), 2), symmetric = TRUE, R = 0,
		burn = 10, nscan = 20, odens = 1, verbose = FALSE, plot = FALSE, gof = FALSE)
	expect_identical(colnames(l$BETA), c("intercept", "x_node", "x.1_node"))
})
