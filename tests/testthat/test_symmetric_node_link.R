skip_on_cran()

# a symmetric lame fit stores one <v>_node coefficient for the <v>_row and <v>_col slices of its
# design: its stored link, its predictions and its simulations take the term through both actors
sym_link_data = function(n = 12, seed = 5) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n)); x = rnorm(n)
	e = matrix(rnorm(n * n), n)
	Y = 0.5 + outer(x, x, "+") + (e + t(e)) / 2
	diag(Y) = NA; dimnames(Y) = list(nm, nm)
	d = matrix(rnorm(n * n), n); d = (d + t(d)) / 2
	Xd = array(d, c(n, n, 1), list(nm, nm, "d"))
	G = matrix(x, n, 1, dimnames = list(nm, "growth"))
	list(YL = list(t1 = Y, t2 = Y, t3 = Y), XdL = list(Xd, Xd, Xd), GL = list(G, G, G), x = x)
}
fit_sym = function(d, ...) suppressWarnings(suppressMessages(lame(d$YL, Xdyad = d$XdL, Xrow = d$GL,
	symmetric = TRUE, R = 0, burn = 20, nscan = 40, odens = 1, verbose = FALSE, plot = FALSE,
	gof = FALSE, ...)))
link_gap = function(a, b) max(unlist(Map(function(p, q) max(abs(p - q)[row(p) != col(p)]), a, b)))

test_that("a static symmetric lame fit with a node covariate has one link with and without newdata", {
	d = sym_link_data()
	f = fit_sym(d)
	p = predict(f, newdata = d$XdL, type = "link")
	expect_lt(link_gap(f$EZ, p), 1e-10)
	expect_lt(link_gap(predict(f, type = "link"), p), 1e-10)
	# the node term enters through both actors: beta * (x_i + x_j)
	b = colMeans(f$BETA)
	own = b[["intercept"]] + b[["growth_node"]] * outer(d$x, d$x, "+") + b[["d_dyad"]] * d$XdL[[1]][, , 1] +
		outer(f$APM, f$APM, "+")
	expect_lt(max(abs(own - f$EZ[[1]])[row(own) != col(own)]), 1e-10)
	# simulated networks centre on the fitted values
	s = simulate(f, nsim = 200, type = "approximate", seed = 1)
	m = Reduce(`+`, lapply(s$Y, function(y) y[[1]])) / 200
	expect_lt(mean(abs(m - fitted(f)[[1]])[row(m) != col(m)]), 0.15)
})

test_that("a dynamic_beta symmetric lame fit with a node covariate has one link with and without newdata", {
	d = sym_link_data()
	f = fit_sym(d, dynamic_beta = TRUE)
	p = predict(f, newdata = d$XdL, type = "link")
	expect_lt(link_gap(f$EZ, p), 1e-10)
	expect_lt(link_gap(predict(f, type = "link"), p), 1e-10)
})

test_that("a symmetric Poisson lame fit with a node covariate forms its fitted counts from that link", {
	d = sym_link_data()
	d$YL = lapply(d$YL, function(y) round(exp(y / 3)))
	f = fit_sym(d, family = "poisson")
	p = predict(f, newdata = d$XdL, type = "link")
	ve = mean(f$VC[, "ve"])
	expect_lt(link_gap(fitted(f), lapply(p, function(e) exp(e + ve / 2))), 1e-8)
})
