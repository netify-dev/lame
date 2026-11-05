skip_on_cran()

# symmetric = TRUE needs a symmetric regression term: an asymmetric Xdyad slice is refused by name

symmetric_inputs = function(n = 8L, seed = 1) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n); Y = (Y + t(Y)) / 2
	diag(Y) = NA; dimnames(Y) = list(nm, nm)
	even = matrix(rnorm(n * n), n, n); even = (even + t(even)) / 2
	odd = matrix(rnorm(n * n), n, n)
	list(
		Y = Y,
		nm = nm,
		n = n,
		X = array(c(even, odd), c(n, n, 2), dimnames = list(nm, nm, c("x_even", "x_odd")))
	)
}

test_that("ame() refuses an asymmetric Xdyad slice under symmetric = TRUE and names it", {
	d = symmetric_inputs()
	e = expect_error(
		ame(
			d$Y,
			Xdyad = d$X,
			symmetric = TRUE,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		class = "lame_symmetric_xdyad"
	)
	msg = conditionMessage(e)
	expect_match(msg, "x_odd")
	expect_false(grepl("x_even", msg))
	expect_match(msg, "(X + t(X)) / 2", fixed = TRUE)
	# a slice without a name is named by its position
	X2 = d$X; dimnames(X2)[[3]] = NULL
	e2 = expect_error(
		ame(
			d$Y,
			Xdyad = X2,
			symmetric = TRUE,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		class = "lame_symmetric_xdyad"
	)
	expect_match(conditionMessage(e2), "k=2")
	# two bad slices are both named
	X3 = d$X; X3[, , 1] = matrix(rnorm(d$n^2), d$n, d$n)
	e3 = expect_error(
		ame(
			d$Y,
			Xdyad = X3,
			symmetric = TRUE,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		class = "lame_symmetric_xdyad"
	)
	expect_match(conditionMessage(e3), "x_even")
	expect_match(conditionMessage(e3), "x_odd")
})

test_that("the symmetrised slice, and a directed fit with the asymmetric one, are accepted", {
	d = symmetric_inputs()
	X = d$X; X[, , 2] = (X[, , 2] + t(X[, , 2])) / 2
	fit = suppressWarnings(ame(
		d$Y,
		Xdyad = X,
		symmetric = TRUE,
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	))
	expect_s3_class(fit, "ame")
	expect_true(isSymmetric(unname(fit$YPM)))
	# the covariate is only required to be symmetric when the fit is
	Yd = d$Y; Yd[upper.tri(Yd)] = Yd[upper.tri(Yd)] + rnorm(sum(upper.tri(Yd)), 0, .5)
	fitd = suppressWarnings(ame(
		Yd,
		Xdyad = d$X,
		symmetric = FALSE,
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	))
	expect_s3_class(fitd, "ame")
	# a single matrix covariate is checked as a one-slice array
	e = expect_error(
		ame(
			d$Y,
			Xdyad = d$X[, , 2],
			symmetric = TRUE,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		),
		class = "lame_symmetric_xdyad"
	)
})

# the same checks for a panel: lame() refuses an asymmetric Xdyad slice under
# symmetric = TRUE, by slice name and first period, as ame() does
symmetric_panel = function(n = 8L, Tt = 3L, seed = 2) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	sym = function() { m = matrix(rnorm(n * n), n, n); (m + t(m)) / 2 }
	Y = lapply(seq_len(Tt), function(t) { y = sym(); diag(y) = NA; dimnames(y) = list(nm, nm); y })
	names(Y) = paste0("t", seq_len(Tt))
	X = lapply(
		seq_len(Tt),
		function(t)
			array(c(sym(), sym()), c(n, n, 2), dimnames = list(nm, nm, c("x_even", "x_odd")))
	)
	names(X) = names(Y)
	list(Y = Y, X = X, nm = nm, n = n, Tt = Tt)
}
fit_panel = function(p, X = p$X, ..., symmetric = TRUE) {
	lame(
		p$Y,
		Xdyad = X,
		symmetric = symmetric,
		family = "normal",
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		...
	)
}

test_that("lame() refuses an asymmetric Xdyad slice under symmetric = TRUE and names it and its period", {
	p = symmetric_panel()
	# asymmetric in the second period only
	p$X[[2]][, , 2] = matrix(rnorm(p$n^2), p$n, p$n)
	e = expect_error(fit_panel(p), class = "lame_symmetric_xdyad")
	msg = conditionMessage(e)
	expect_match(msg, "x_odd")
	expect_false(grepl("x_even", msg))
	expect_match(msg, "period t2")
	expect_match(msg, "(X + t(X)) / 2", fixed = TRUE)
	# a slice without a name is named by its position (X2: the name lame() gives it)
	X2 = lapply(p$X, function(x) { dimnames(x)[[3]] = NULL; x })
	e2 = expect_error(fit_panel(p, X = X2), class = "lame_symmetric_xdyad")
	expect_match(conditionMessage(e2), "X2")
	# two bad slices are both named
	X3 = p$X; X3[[1]][, , 1] = matrix(rnorm(p$n^2), p$n, p$n)
	e3 = expect_error(fit_panel(p, X = X3), class = "lame_symmetric_xdyad")
	expect_match(conditionMessage(e3), "x_even")
	expect_match(conditionMessage(e3), "x_odd")
	expect_match(conditionMessage(e3), "period t1")
	# a single-matrix slice per period is a one-slice array
	X4 = lapply(p$X, function(x) x[, , 2])
	expect_error(fit_panel(p, X = X4), class = "lame_symmetric_xdyad")
})

test_that("lame() accepts the symmetrised slice, and a directed fit with the asymmetric one", {
	p = symmetric_panel()
	p$X[[2]][, , 2] = matrix(rnorm(p$n^2), p$n, p$n)
	X = lapply(p$X, function(x) { x[, , 2] = (x[, , 2] + t(x[, , 2])) / 2; x })
	fit = suppressWarnings(fit_panel(p, X = X))
	expect_s3_class(fit, "lame")
	# the covariate is only required to be symmetric when the fit is
	pd = p
	pd$Y = lapply(p$Y, function(y) { y[upper.tri(y)] = y[upper.tri(y)] + rnorm(sum(upper.tri(y)), 0, .5); y })
	fitd = suppressWarnings(fit_panel(pd, symmetric = FALSE))
	expect_s3_class(fitd, "lame")
})

test_that("a hole in one triangle of a symmetric covariate is not an asymmetry", {
	p = symmetric_panel()
	p$X[[1]][2, 5, 1] = NA
	p$X[[3]][4, 1, 2] = NA
	fit = suppressWarnings(suppressMessages(fit_panel(p)))
	expect_s3_class(fit, "lame")
	expect_false(anyNA(fit$BETA))
	y1 = suppressWarnings(suppressMessages(ame(
		p$Y[[1]],
		Xdyad = p$X[[1]],
		symmetric = TRUE,
		R = 0,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	)))
	expect_s3_class(y1, "ame")
	expect_false(is.na(y1$Y[2, 5]))
})

# exactly symmetric Y with symmetric = FALSE: a directed normal model with a dyadic correlation has an improper posterior
# (rho runs to one) once the pairs outnumber what the mean absorbs; lame() refuses it as ame() does, other directed fits run and state it
test_that("lame() refuses a directed normal model with a dyadic correlation on exactly symmetric Y", {
	p = symmetric_panel()
	flat = function(e) gsub("[[:space:]]+", " ", conditionMessage(e))
	args = list(family = "normal", nscan = 10, burn = 2, odens = 1, verbose = FALSE, gof = FALSE, plot = FALSE)
	# static or dynamic positions, or none: the pairs exceed the bound in every period
	for (cfg in list(list(R = 0), list(R = 1), list(R = 1, dynamic_uv = TRUE))) {
		e = expect_error(do.call(lame, c(list(p$Y), cfg, args)), class = "lame_symmetric_y")
		expect_s3_class(e, "lame_improper_normal_reciprocity")
		expect_match(flat(e), "exactly symmetric")
		expect_match(flat(e), "symmetric = TRUE")
	}
	# a refused fit gets the error alone, not the note about the posterior of rho
	expect_no_message(
		expect_error(do.call(lame, c(list(p$Y, R = 0), args)), class = "lame_symmetric_y"),
		class = "lame_symmetric_y_note"
	)
	# the undirected model fits the same panel, and the directed model without the dyadic correlation too
	fit = suppressMessages(
		suppressWarnings(do.call(lame, c(list(p$Y, R = 1, dynamic_uv = TRUE, symmetric = TRUE), args)))
	)
	expect_s3_class(fit, "lame")
	note = "lame_symmetric_y_note"
	expect_message(
		suppressWarnings(do.call(lame, c(list(p$Y, R = 1, dynamic_uv = TRUE, dcor = FALSE), args))),
		class = note
	)
	# one asymmetric dyad lifts the refusal
	pd = p
	pd$Y[[1]][2, 1] = pd$Y[[1]][2, 1] + 1
	fitd = suppressMessages(suppressWarnings(do.call(lame, c(list(pd$Y, R = 1, dynamic_uv = TRUE), args))))
	expect_s3_class(fitd, "lame")
})

test_that("a small normal panel that the pair-count bound leaves proper is fitted, with the note", {
	# four actors: six pairs against an actor rank of three
	Y4 = matrix(c(NA, .4, .8, 1.1, .4, NA, 1.2, .3, .8, 1.2, NA, .6, 1.1, .3, .6, NA), 4L, 4L)
	dimnames(Y4) = list(letters[1:4], letters[1:4])
	expect_message({fit = suppressWarnings(lame(
		list(t1 = Y4),
		family = "normal",
		R = 0,
		g = 4,
		nscan = 30,
		burn = 5,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))}, class = "lame_symmetric_y_note")
	expect_s3_class(fit, "lame")
})

# the propriety bound of a panel is pooled: the periods share rho, so the excess of every period adds to the
# singular power of the density of rho at one, and parameters the periods share absorb their contrasts once
test_that("lame() applies the propriety bound to the pooled panel", {
	Y4 = matrix(c(NA, .4, .8, 1.1, .4, NA, 1.2, .3, .8, 1.2, NA, .6, 1.1, .3, .6, NA), 4L, 4L)
	dimnames(Y4) = list(letters[1:4], letters[1:4])
	flat = function(e) gsub("[[:space:]]+", " ", conditionMessage(e))
	args = list(
		family = "normal", R = 0, g = 4, nscan = 20, burn = 4, odens = 1, verbose = FALSE, gof = FALSE, plot = FALSE
	)
	# each period alone is proper (six pairs, actor rank three), two periods are not: twelve pairs, the same three directions
	e = expect_error(do.call(lame, c(list(list(t1 = Y4, t2 = Y4)), args)), class = "lame_symmetric_y")
	expect_s3_class(e, "lame_improper_normal_reciprocity")
	expect_match(flat(e), "all 12 reciprocal pairs")
	expect_match(flat(e), "at most 3 of their contrasts")
	# actor effects of their own in every period absorb three per period; the periods' excess still adds (12 - 6 >= 4)
	expect_error(do.call(lame, c(list(list(t1 = Y4, t2 = Y4), dynamic_ab = TRUE), args)), class = "lame_symmetric_y")
	# a panel whose pairs stay under the pooled bound by four is fitted: two periods of three observed pairs each
	Ya = Y4; Ya[1, 4] = Ya[4, 1] = Ya[2, 4] = Ya[4, 2] = Ya[1, 3] = Ya[3, 1] = NA
	Yb = Y4; Yb[1, 2] = Yb[2, 1] = Yb[2, 3] = Yb[3, 2] = Yb[3, 4] = Yb[4, 3] = NA
	expect_message(
		{fit = suppressWarnings(do.call(lame, c(list(list(t1 = Ya, t2 = Yb)), args)))}, class = "lame_symmetric_y_note")
	expect_s3_class(fit, "lame")
	# the pooled count agrees with the single-network certificate on one period
	b = lame:::.lame_panel_reciprocity_bound(array(Y4, c(4, 4, 1)), list(array(1, c(4, 4, 1))), TRUE, TRUE, 0L)
	expect_equal(c(b$pairs, b$rank), c(6, 3))
	expect_false(b$improper)
	b2 = lame:::.lame_panel_reciprocity_bound(array(Y4, c(4, 4, 2)), rep(list(array(1, c(4, 4, 1))), 2), TRUE, TRUE, 0L)
	expect_equal(c(b2$pairs, b2$rank), c(12, 3))
	expect_true(b2$improper)
	# moving positions and moving coefficients each count per period
	b3 = lame:::.lame_panel_reciprocity_bound(
		array(Y4, c(4, 4, 2)),
		rep(list(array(1, c(4, 4, 1))), 2),
		TRUE,
		TRUE,
		1L,
		dynamic_uv = TRUE
	)
	expect_equal(b3$factor, 8L)
	# a covariate that is asymmetric at a pair adds a direction, once if it is static and per period if it moves
	X = array(1, c(4, 4, 2)); X[1, 2, 2] = 3
	b4 = lame:::.lame_panel_reciprocity_bound(array(Y4, c(4, 4, 2)), rep(list(X), 2), TRUE, TRUE, 0L)
	expect_equal(b4$coef, 1L)
	b5 = lame:::.lame_panel_reciprocity_bound(
		array(Y4, c(4, 4, 2)),
		rep(list(X), 2),
		TRUE,
		TRUE,
		0L,
		dynamic_beta = TRUE
	)
	expect_equal(b5$coef, 2L)
})

test_that("exactly symmetric binary and Poisson panels fit under the directed dyadic-correlation model with dynamic positions", {
	skip_on_cran()
	p = symmetric_panel(n = 10, Tt = 3)
	note = "lame_symmetric_y_note"
	args = list(R = 1, dynamic_uv = TRUE, nscan = 20, burn = 5, odens = 5, verbose = FALSE, gof = FALSE, plot = FALSE)
	yb = lapply(p$Y, function(y) { z = 1 * (y > 0); diag(z) = NA; z })
	e = expect_message({fit = suppressWarnings(do.call(lame, c(list(yb, family = "binary"), args)))}, class = note)
	expect_s3_class(fit, "lame")
	expect_match(gsub("[[:space:]]+", " ", conditionMessage(e)), "concentrates near one")
	yp = lapply(p$Y, function(y) { z = round(exp(y) * 3); diag(z) = NA; z })
	expect_message({fit = suppressWarnings(do.call(lame, c(list(yp, family = "poisson"), args)))}, class = note)
	expect_s3_class(fit, "lame")
})

# the fast estimators refuse an asymmetric dyadic covariate under symmetric = TRUE, as the samplers do
test_that("the ALS routes refuse an asymmetric Xdyad slice under symmetric = TRUE and name it", {
	d = symmetric_inputs()
	e = expect_error(
		ame_als(d$Y, Xdyad = d$X, symmetric = TRUE, R = 0, verbose = FALSE),
		class = "lame_symmetric_xdyad"
	)
	expect_match(conditionMessage(e), "x_odd")
	expect_false(grepl("x_even", conditionMessage(e)))
	expect_error(
		suppressMessages(ame(d$Y, Xdyad = d$X, symmetric = TRUE, R = 0, method = "als", verbose = FALSE)),
		class = "lame_symmetric_xdyad"
	)
	p = symmetric_panel()
	p$X[[2]][, , 2] = matrix(rnorm(p$n^2), p$n, p$n)
	e = expect_error(
		lame_als(p$Y, Xdyad = p$X, symmetric = TRUE, R = 0, verbose = FALSE),
		class = "lame_symmetric_xdyad"
	)
	expect_match(conditionMessage(e), "x_odd")
	expect_match(conditionMessage(e), "period t2")
	expect_error(
		suppressMessages(
			lame(p$Y, Xdyad = p$X, symmetric = TRUE, R = 0, method = "als", family = "normal", verbose = FALSE)
		),
		class = "lame_symmetric_xdyad"
	)
	expect_error(
		suppressMessages(lame(
			p$Y,
			Xdyad = p$X,
			symmetric = TRUE,
			R = 0,
			method = "als",
			family = "normal",
			dynamic_beta = "dyad",
			verbose = FALSE
		)),
		class = "lame_symmetric_xdyad"
	)
	# the symmetrised slice is fitted, and a directed fit keeps the asymmetric one
	X = lapply(p$X, function(x) { x[, , 2] = (x[, , 2] + t(x[, , 2])) / 2; x })
	expect_s3_class(suppressMessages(lame_als(p$Y, Xdyad = X, symmetric = TRUE, R = 0, verbose = FALSE)), "lame_als")
	Yd = d$Y; Yd[upper.tri(Yd)] = Yd[upper.tri(Yd)] + rnorm(sum(upper.tri(Yd)), 0, .5)
	expect_s3_class(suppressMessages(ame_als(Yd, Xdyad = d$X, symmetric = FALSE, R = 0, verbose = FALSE)), "ame_als")
})

# the fast normal estimator's refusal counts the absorbed rank with the helper of the samplers (static effects, one count
# for what the periods share)
test_that("the ALS Gaussian reciprocity check shares the pooled rank count", {
	Y4 = matrix(c(NA, .4, .8, 1.1, .4, NA, 1.2, .3, .8, 1.2, NA, .6, 1.1, .3, .6, NA), 4L, 4L)
	prep = function(Tt, X = NULL) list(
		Yarr = array(Y4, c(4, 4, Tt)),
		nr = 4L,
		Tt = Tt,
		p_dyad = if (is.null(X)) 0L else dim(X)[3L],
		X = X
	)
	chk = function(...) lame:::.ae_gauss_check_reciprocity(...)
	expect_null(chk(prep(1L), 0L))
	e = expect_error(chk(prep(2L), 0L), class = "lame_improper_normal_reciprocity")
	expect_match(gsub("[[:space:]]+", " ", conditionMessage(e)), "All 12 observed reciprocal pairs.*bound is 3")
	# a factor of rank one absorbs n of the contrasts per fit, not per period
	expect_null(chk(prep(1L), 1L))
	X = array(1, c(4, 4, 1, 2)); X[1, 2, 1, 2] = 3
	e = expect_error(chk(prep(2L, X), 0L), class = "lame_improper_normal_reciprocity")
	expect_match(gsub("[[:space:]]+", " ", conditionMessage(e)), "bound is 4")
})
