skip_on_cran()

# the latent product of a normal empirical-Bayes fit that the data do not determine is withheld (R/ame_als_gauss.R)

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

# the withholding rule at one kind of network: a weak product is withheld with its rank evidence, a strong one is reported
wh_threshold_check = function(kind) {
	weak = wh_fit(wh_net(kind = kind, suv = 0.05, seed = 3))
	info = weak$eb$variational$information
	expect_equal(info[["threshold"]], 2.1, label = kind)
	expect_lt(info[["per_loading"]], 2.1, label = kind)
	expect_true(lame:::.ae_product_is_withheld(weak), label = kind)
	pw = weak$eb$product_withheld
	expect_true(all(is.na(weak$U)) && all(is.na(weak$UVPM)), label = kind)
	if (kind != "symmetric") expect_true(
		all(is.na(weak[["V"]])),
		label = kind
	) else expect_true(all(is.na(weak$L)), label = kind)
	# the slots keep their dimensions and names; the means are kept for the computations that start from them
	expect_identical(dim(weak$U), dim(pw$means$U)); expect_identical(dimnames(weak$UVPM), dimnames(pw$means$UVPM))
	expect_true(all(is.finite(pw$means$UVPM[!is.na(weak$Y[, , 1])])), label = kind)
	expect_identical(pw$rank, 2L)
	expect_length(pw$log_marginal, 3L)
	expect_true(all(is.finite(pw$log_marginal)), label = kind)
	expect_identical(pw$rank_chosen, as.integer(which.max(pw$log_marginal) - 1L))
	expect_match(wh_print(weak), "not identified by the data at this rank", fixed = TRUE)
	expect_match(wh_print(weak), sprintf("chooses rank %d", pw$rank_chosen), fixed = TRUE)
	expect_match(
		paste(utils::capture.output(print(summary(weak))), collapse = " "),
		"not identified by the data",
		fixed = TRUE
	)
	strong = wh_fit(wh_net(kind = kind, suv = 0.8, seed = 4))
	expect_gte(strong$eb$variational$information[["per_loading"]], 2.1)
	expect_false(lame:::.ae_product_is_withheld(strong), label = kind)
	expect_true(all(is.finite(strong$U)), label = kind)
	expect_false(grepl("not identified", wh_print(strong), fixed = TRUE))
}

test_that("below 2.1 noise variances of signal per loading the product is withheld, with the rank evidence; above it is reported: a directed network", {
	wh_threshold_check("directed")
})

test_that("below 2.1 noise variances of signal per loading the product is withheld, with the rank evidence; above it is reported: undirected and bipartite networks", {
	skip_on_cran()
	for (kind in c("symmetric", "bipartite")) wh_threshold_check(kind)
})

test_that("a withheld fit has the fitted values, predictions, residuals, coefficients and standard errors of the fit, either way", {
	skip_on_cran()
	d = wh_net(suv = 0.05, seed = 3)
	f = wh_fit(d)
	expect_true(lame:::.ae_product_is_withheld(f))
	g = withr::with_options(list(lame.eb_withhold_product = FALSE), wh_fit(d))
	expect_false(lame:::.ae_product_is_withheld(g))
	expect_identical(fitted(f), fitted(g))
	expect_identical(predict(f), predict(g))
	expect_identical(predict(f, newdata = d$X * 0.5), predict(g, newdata = d$X * 0.5))
	expect_identical(residuals(f), residuals(g))
	expect_identical(coef(f), coef(g))
	expect_identical(vcov(f), vcov(g))
	expect_identical(f$VC, g$VC)
	# fitted values of a normal fit use the variational predictor, whose product is the withheld mean
	expect_identical(f$eb$product_withheld$means$UVPM, g$UVPM)
	expect_identical(f$eb$product_withheld$means$U, g$U)
	# above the threshold the switch changes nothing
	s = wh_net(suv = 0.8, seed = 4)
	f2 = wh_fit(s); g2 = withr::with_options(list(lame.eb_withhold_product = FALSE), wh_fit(s))
	f2$call = g2$call = NULL
	expect_identical(f2, g2)
})

test_that("the positions of a withheld fit are NA, its positions are not plotted, and a refit withholds as the fit does", {
	f = wh_fit(wh_net(suv = 0.05, seed = 3))
	lp = suppressMessages(latent_positions(f))
	expect_gt(nrow(lp), 0L)
	expect_true(all(is.na(lp$value)))
	expect_message(latent_positions(f), "withheld")
	expect_error(uv_plot(f), class = "lame_product_withheld")
	r = ame_als_refit(f, verbose = FALSE)
	expect_true(lame:::.ae_product_is_withheld(r))
	expect_true(r$converged)
	sv = als_start_vals(f)
	expect_true(all(is.finite(sv$U)))
})

test_that("the bootstrap of a withheld fit draws from the fit and reports no interval for the positions", {
	skip_on_cran()
	f = wh_fit(wh_net(suv = 0.05, seed = 3))
	b = ame_als_bootstrap(f, R = 30, seed = 2, verbose = FALSE)
	expect_true(all(c("U", "V") %in% b$withheld$parameter))
	expect_null(b$se_U)
	expect_true(all(is.finite(b$ci_lo)) && all(is.finite(b$ci_hi)))
})
