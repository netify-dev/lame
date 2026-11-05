skip_on_cran()

# a regression design that is exactly rank-deficient over the observed cells has
# coefficients the data cannot separate: ame() and lame() refuse it, in every
# mode and with dynamic_beta, with a classed error that names the aliased columns

rank_panel = function(kind = "dup", mode = "unipartite", Tn = 1L, family = "normal", n = 9L, nB = 6L, seed = 1) {
	set.seed(seed)
	bip = mode == "bipartite"; sym = mode == "symmetric"
	nB = if (bip) nB else n
	ra = paste0("r", seq_len(n)); ca = if (bip) paste0("c", seq_len(nB)) else ra
	Y = X = vector("list", Tn)
	for (t in seq_len(Tn)) {
		x1 = matrix(rnorm(n * nB), n, nB); x2 = matrix(rnorm(n * nB), n, nB)
		if (sym) { x1 = (x1 + t(x1)) / 2; x2 = (x2 + t(x2)) / 2 }
		x3 = switch(
			kind,
			dup = x1,
			combo = x1 + 2 * x2,
			indep = {
				x = matrix(rnorm(n * nB), n, nB); if (sym) (x + t(x)) / 2 else x },
			near = x1 + 1e-5 * matrix(rnorm(n * nB), n, nB)
		)
		z = 0.3 * x1 + matrix(rnorm(n * nB, 0, .7), n, nB)
		if (sym) z = (z + t(z)) / 2
		y = switch(family, normal = z, binary = 1 * (z > 0), poisson = matrix(rpois(n * nB, exp(.3 * x1 + .2)), n, nB))
		if (sym && family == "poisson") y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (!bip) diag(y) = NA
		dimnames(y) = list(ra, ca); Y[[t]] = y
		X[[t]] = array(c(x1, x2, x3), c(n, nB, 3), dimnames = list(ra, ca, c("x1", "x2", "x3")))
	}
	list(Y = Y, X = X, mode = mode, bip = bip, sym = sym)
}
quick = list(nscan = 12, burn = 4, odens = 2, verbose = FALSE, gof = FALSE)
fit_rank = function(d, panel = FALSE, family = "normal", ...) {
	mode = if (d$bip) "bipartite" else "unipartite"
	if (panel) {
		args = c(
			list(d$Y, Xdyad = d$X, family = family, R = 0, mode = mode, symmetric = d$sym, plot = FALSE),
			quick,
			list(...)
		)
	} else {
		args = c(
			list(d$Y[[1]], Xdyad = d$X[[1]], family = family, R = 0, mode = mode, symmetric = d$sym),
			quick,
			list(...)
		)
	}
	if (d$bip) { args$R = NULL; args$R_row = 0; args$R_col = 0 }
	suppressWarnings(suppressMessages(do.call(if (panel) lame else ame, args)))
}

test_that("the check names the columns aliased over the observed cells", {
	check = lame:::.lame_check_design_rank
	set.seed(3); n = 7L
	x1 = matrix(rnorm(n * n), n, n); x2 = matrix(rnorm(n * n), n, n)
	y = matrix(rnorm(n * n), n, n); diag(y) = NA
	arr = function(...) { l = list(...); array(unlist(l), c(n, n, length(l)), dimnames = list(NULL, NULL, names(l))) }
	expect_null(check(arr(a = x1, b = x2), y))
	# a duplicate: the later column is the aliased one, as lm() reports it
	e = expect_error(check(arr(a = x1, b = x2, c = x1), y), class = "lame_design_rank")
	expect_identical(e$aliased, "c"); expect_identical(e$rank, 2L); expect_identical(e$n_col, 3L)
	expect_match(conditionMessage(e), "`c` is aliased")
	expect_match(conditionMessage(e), "c = 1 a", fixed = TRUE)
	expect_match(conditionMessage(e), "On the 42 observed cells of 1 period")
	# a combination with weights
	e = expect_error(check(arr(a = x1, b = x2, c = 2 * x1 - 3 * x2), y), class = "lame_design_rank")
	expect_match(conditionMessage(e), "c = -3 b + 2 a", fixed = TRUE)
	# two aliased columns are both named
	e = expect_error(check(arr(a = x1, b = x2, c = x1 + x2, d = x2), y), class = "lame_design_rank")
	expect_identical(e$aliased, c("c", "d"))
	# a column that is zero on every observed cell; a single such column too
	e = expect_error(check(arr(a = x1, z = 0 * x1), y), class = "lame_design_rank")
	expect_identical(e$aliased, "z"); expect_match(conditionMessage(e), "z is zero on every observed cell")
	expect_error(check(arr(z = 0 * x1), y), class = "lame_design_rank")
	# unnamed columns are named by position; braces in a name are not markup
	unnamed = arr(a = x1, b = x1); dimnames(unnamed) = NULL
	e = expect_error(check(unnamed, y), class = "lame_design_rank")
	expect_identical(e$aliased, "column 2")
	e = expect_error(check(arr(`a{1}` = x1, `b}2{` = x1), y), class = "lame_design_rank")
	expect_identical(e$aliased, "b}2{")
	expect_match(conditionMessage(e), "b}2{ = 1 a{1}", fixed = TRUE)
	# near-collinear designs are not refused, and neither is a design without columns
	expect_null(check(arr(a = x1, b = x1 + 1e-6 * rnorm(n * n)), y))
	expect_null(check(array(0, c(n, n, 0)), y))
	expect_null(check(NULL, y))
})

test_that("only the cells with an observed response and a finite design count", {
	check = lame:::.lame_check_design_rank
	set.seed(4); n = 7L
	x1 = matrix(rnorm(n * n), n, n); x3 = x1
	y = matrix(rnorm(n * n), n, n); diag(y) = NA
	# the columns differ only where the response is missing: aliased on what is observed
	hide = matrix(FALSE, n, n); hide[2, 3] = hide[5, 1] = TRUE
	x3[hide] = x3[hide] + 3; y[hide] = NA
	arr = array(c(x1, x3), c(n, n, 2), dimnames = list(NULL, NULL, c("a", "b")))
	e = expect_error(check(arr, y), class = "lame_design_rank")
	expect_identical(e$aliased, "b"); expect_match(conditionMessage(e), "On the 40 observed cells")
	# where they differ on an observed cell they are separate
	x3[2, 4] = x3[2, 4] + 3
	arr = array(c(x1, x3), c(n, n, 2), dimnames = list(NULL, NULL, c("a", "b")))
	expect_null(check(arr, y))
	# a cell with a non-finite design value is out, whatever its response
	arr[2, 4, 2] = NA
	expect_error(check(arr, y), class = "lame_design_rank")
})

test_that("the pooled design is what counts, and each period is reduced exactly", {
	check = lame:::.lame_check_design_rank
	set.seed(5); n = 6L; Tn = 4L
	y = array(rnorm(n * n * Tn), c(n, n, Tn)); for (t in 1:Tn) diag(y[, , t]) = NA
	mk = function(dup_in) lapply(1:Tn, function(t) {
		x1 = matrix(rnorm(n * n), n, n)
		array(
			c(x1, if (t %in% dup_in) x1 else matrix(rnorm(n * n), n, n)),
			c(n, n, 2),
			dimnames = list(NULL, NULL, c("p", "q"))
		)
	})
	# aliased in every period: refused; in some periods only: the pooled design is full rank
	e = expect_error(check(mk(1:Tn), y), class = "lame_design_rank")
	expect_match(conditionMessage(e), "of 4 periods")
	expect_null(check(mk(1L), y))
	expect_null(check(mk(1:3), y))
	# a set of period indicators is full rank pooled; with a constant column it is not
	ind = lapply(
		1:Tn,
		function(t) { a = array(0, c(n, n, Tn), dimnames = list(NULL, NULL, paste0("t", 1:Tn))); a[, , t] = 1; a }
	)
	expect_null(check(ind, y))
	withint = lapply(
		ind,
		function(a) array(
			c(rep(1, n * n), a),
			c(n, n, Tn + 1L),
			dimnames = list(NULL, NULL, c("int", paste0("t", 1:Tn)))
		)
	)
	e = expect_error(check(withint, y), class = "lame_design_rank")
	expect_identical(e$rank, Tn)
	# the rank of the reduced blocks is the rank of the stacked design
	xs = mk(2:3)
	stacked = do.call(rbind, lapply(1:Tn, function(t) {
		m = matrix(xs[[t]], ncol = 2); m[!is.na(c(y[, , t])), , drop = FALSE] }))
	expect_identical(qr(stacked)$rank, 2L)
})

test_that("ame() refuses an exactly aliased design in every mode and family", {
	for (mode in c(
		"unipartite", "symmetric", "bipartite"
	)) for (kind in c("dup", "combo")) for (family in c("normal", "binary")) {
		d = rank_panel(kind, mode, family = family)
		e = expect_error(fit_rank(d, family = family), class = "lame_design_rank", info = paste(mode, kind, family))
		expect_identical(e$aliased, "x3_dyad", info = paste(mode, kind, family))
		expect_match(
			conditionMessage(e),
			if (kind == "dup") "x3_dyad = 1 x1_dyad" else "x3_dyad = 2 x2_dyad + 1 x1_dyad",
			fixed = TRUE
		)
	}
	# counts: the sampler never ran
	d = rank_panel("dup", "unipartite", family = "poisson")
	expect_error(fit_rank(d, family = "poisson"), class = "lame_design_rank")
})

test_that("lame() refuses an exactly aliased panel design, static or with dynamic_beta, in every mode", {
	for (mode in c("unipartite", "symmetric", "bipartite")) for (kind in c("dup", "combo")) {
		d = rank_panel(kind, mode, Tn = 3L)
		for (dyn in list(FALSE, TRUE, "dyad")) {
			info = paste(mode, kind, format(dyn))
			e = expect_error(fit_rank(d, panel = TRUE, dynamic_beta = dyn), class = "lame_design_rank", info = info)
			expect_identical(e$aliased, "x3_dyad", info = info)
			expect_match(conditionMessage(e), "of 3 periods", info = info)
		}
	}
	# the aliased columns are the ones the coefficients index, node covariates included
	n = 9L; set.seed(6); nm = paste0("r", 1:n)
	y = lapply(1:2, function(t) { m = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	r1 = rnorm(n); r2 = rnorm(n)
	e = expect_error(
		suppressWarnings(suppressMessages(lame(
			y,
			Xrow = lapply(1:2, function(t) cbind(a = r1, b = r2, c = r1 + r2)),
			R = 0,
			plot = FALSE,
			nscan = 6,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		class = "lame_design_rank"
	)
	expect_identical(e$aliased, "c_row")
})

test_that("node covariates and a dyadic covariate built from them are aliased, in ame()", {
	n = 9L; set.seed(7); nm = paste0("r", 1:n)
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	r1 = rnorm(n); r2 = rnorm(n)
	fit = function(...) suppressWarnings(
		suppressMessages(ame(Y, R = 0, nscan = 6, burn = 2, odens = 1, verbose = FALSE, gof = FALSE, ...))
	)
	e = expect_error(fit(Xrow = cbind(a = r1, b = r1)), class = "lame_design_rank")
	expect_identical(e$aliased, "b_row")
	e = expect_error(fit(Xcol = cbind(a = r1, b = r1)), class = "lame_design_rank")
	expect_identical(e$aliased, "b_col")
	sumd = array(outer(r1, r2, "+"), c(n, n, 1), dimnames = list(nm, nm, "s"))
	e = expect_error(fit(Xrow = cbind(r = r1), Xcol = cbind(c = r2), Xdyad = sumd), class = "lame_design_rank")
	expect_identical(e$aliased, "s_dyad")
	expect_match(conditionMessage(e), "s_dyad = 1 r_row + 1 c_col", fixed = TRUE)
	e = expect_error(fit(Xdyad = array(0, c(n, n, 1), dimnames = list(nm, nm, "z"))), class = "lame_design_rank")
	expect_match(conditionMessage(e), "z_dyad is zero on every observed cell")
	# independent covariates fit
	expect_s3_class(fit(Xrow = cbind(a = r1, b = r2), Xcol = cbind(c = rnorm(n))), "ame")
})

test_that("near-collinear designs are fitted, and so are designs aliased in some periods only", {
	d = rank_panel("near", "unipartite", family = "normal")
	expect_s3_class(fit_rank(d), "ame")
	d = rank_panel("indep", "symmetric", family = "binary")
	expect_s3_class(fit_rank(d, family = "binary"), "ame")
	d = rank_panel("near", "bipartite", Tn = 2L)
	expect_s3_class(fit_rank(d, panel = TRUE), "lame")
	# x3 repeats x1 in the first period only: the pooled design has full rank
	d = rank_panel("indep", "unipartite", Tn = 3L)
	d$X[[1]][, , 3] = d$X[[1]][, , 1]
	expect_s3_class(fit_rank(d, panel = TRUE), "lame")
	expect_s3_class(fit_rank(d, panel = TRUE, dynamic_beta = TRUE), "lame")
})

test_that("a constant covariate is aliased with the intercept where the intercept is kept", {
	set.seed(8); nA = 8L; nB = 6L
	Y = matrix(rnorm(nA * nB), nA, nB, dimnames = list(paste0("r", 1:nA), paste0("c", 1:nB)))
	e = expect_error(
		suppressWarnings(suppressMessages(ame(
			Y,
			Xcol = cbind(k = rep(1, nB)),
			mode = "bipartite",
			R_row = 0,
			R_col = 0,
			nscan = 6,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		class = "lame_design_rank"
	)
	expect_identical(e$aliased, "k_col")
	expect_match(conditionMessage(e), "k_col = 1 intercept", fixed = TRUE)
	# the directed design drops the redundant intercept itself
	n = 8L; nm = paste0("a", 1:n)
	Yd = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Yd) = NA
	fit = suppressWarnings(suppressMessages(ame(
		Yd,
		Xrow = cbind(k = rep(2, n)),
		R = 0,
		nscan = 6,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	)))
	expect_false("intercept" %in% colnames(fit$BETA))
})
