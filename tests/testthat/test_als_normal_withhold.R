skip_on_cran()

# the latent product of a static normal empirical-Bayes fit that the data barely determine is reported with the
# standard errors of its cells and positions (R/ame_als_gauss.R, R/ame_als_report_se.R)

wh_net = function(n = 30, kind = "directed", R = 2, suv = 0.5, seed = 1) {
	set.seed(seed)
	sym = kind == "symmetric"; bip = kind == "bipartite"
	m = if (bip) n + 6 else n
	rn = sprintf("a%02d", seq_len(n)); cn = if (bip) sprintf("b%02d", seq_len(m)) else rn
	a = rnorm(n, 0, 0.7); b = if (sym) a else rnorm(m, 0, 0.7)
	O = if (sym) { F = matrix(rnorm(n * R, 0, sqrt(suv)), n); F %*% t(F) } else tcrossprod(
		matrix(rnorm(n * R, 0, sqrt(suv)), n),
		matrix(rnorm(m * R, 0, sqrt(suv)), m)
	)
	X = array(rnorm(n * m), c(n, m, 1), dimnames = list(rn, cn, "x"))
	if (sym) X[, , 1][lower.tri(X[, , 1])] = t(X[, , 1])[lower.tri(X[, , 1])]
	E = matrix(rnorm(n * m), n, m)
	if (sym) E[lower.tri(E)] = t(E)[lower.tri(E)]
	Y = 0.2 + outer(a, b, "+") + O + 0.5 * X[, , 1] + E
	if (!bip) diag(Y) = NA
	dimnames(Y) = list(rn, cn)
	list(Y = Y, X = X, R = R, kind = kind)
}
wh_fit = function(d, ...) suppressWarnings(ame_als(
	d$Y,
	Xdyad = d$X,
	R = d$R,
	family = "normal",
	mode = if (d$kind == "bipartite") "bipartite" else "unipartite",
	symmetric = d$kind == "symmetric",
	verbose = FALSE,
	...
))
wh_print = function(f) paste(utils::capture.output(print(f)), collapse = " ")

# the rule at one kind of network: a weak product is reported with its standard errors and the information note,
# a strong one without the note
wh_threshold_check = function(kind) {
	weak = wh_fit(wh_net(kind = kind, suv = 0.05, seed = 3))
	info = weak$eb$variational$information
	expect_equal(info[["threshold"]], 2.1, label = kind)
	expect_lt(info[["per_loading"]], 2.1, label = kind)
	expect_false(lame:::.ae_product_is_withheld(weak), label = kind)
	expect_true(all(is.finite(weak$U)) && all(is.finite(weak$UVPM)), label = kind)
	if (kind != "symmetric") expect_true(
		all(is.finite(weak[["V"]])) && all(weak$V_sd > 0),
		label = kind
	) else expect_true(all(is.finite(weak$L)) && all(weak$L_sd > 0), label = kind)
	expect_true(all(weak$U_sd > 0) && all(is.finite(weak$U_sd)), label = kind)
	expect_identical(dimnames(weak$UVPM_sd), dimnames(weak$UVPM))
	expect_true(all(weak$UVPM_sd[!is.na(weak$Y[, , 1])] > 0), label = kind)
	expect_match(wh_print(weak), "barely determine it at this rank", fixed = TRUE)
	expect_false(grepl("not identified", wh_print(weak), fixed = TRUE))
	strong = wh_fit(wh_net(kind = kind, suv = 0.8, seed = 4))
	expect_gte(strong$eb$variational$information[["per_loading"]], 2.1)
	expect_false(lame:::.ae_product_is_withheld(strong), label = kind)
	expect_true(all(is.finite(strong$U)), label = kind)
	expect_false(grepl("barely determine", wh_print(strong), fixed = TRUE))
}

test_that("below 2.1 noise variances of signal per loading the product is reported with its standard errors and a note: a directed network", {
	wh_threshold_check("directed")
})

test_that("below 2.1 noise variances of signal per loading the product is reported with its standard errors and a note: undirected and bipartite networks", {
	skip_on_cran()
	for (kind in c("symmetric", "bipartite")) wh_threshold_check(kind)
})

test_that("a weak normal product is the fit the withholding switch leaves, above and below the threshold", {
	skip_on_cran()
	for (suv in c(0.05, 0.8)) {
		d = wh_net(suv = suv, seed = if (suv < 0.5) 3 else 4)
		f = wh_fit(d)
		g = withr::with_options(list(lame.eb_withhold_product = FALSE), wh_fit(d))
		f$call = g$call = NULL
		expect_identical(f, g)
	}
})

test_that("the positions of a weak normal product carry sds, are plotted, and a refit reports them as the fit does", {
	f = wh_fit(wh_net(suv = 0.05, seed = 3))
	lp = suppressMessages(latent_positions(f))
	expect_gt(nrow(lp), 0L)
	expect_true(all(is.finite(lp$value)) && all(lp$posterior_sd > 0))
	expect_s3_class(suppressMessages(uv_plot(f)), "ggplot")
	r = ame_als_refit(f, verbose = FALSE)
	expect_false(lame:::.ae_product_is_withheld(r))
	expect_true(r$converged)
	sv = als_start_vals(f)
	expect_true(all(is.finite(sv$U)))
})

test_that("the bootstrap of a weak normal product gives finite coefficient intervals", {
	skip_on_cran()
	f = wh_fit(wh_net(suv = 0.05, seed = 3))
	b = ame_als_bootstrap(f, R = 30, seed = 2, verbose = FALSE)
	expect_true(all(is.finite(b$ci_lo)) && all(is.finite(b$ci_hi)))
})
