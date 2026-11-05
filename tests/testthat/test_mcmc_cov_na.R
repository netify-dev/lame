skip_on_cran()

library(lame)

# a cell whose dyadic, sender or receiver covariate is missing in its period is
# unobserved in ame() and lame() for unipartite, bipartite and symmetric
# networks, as it is in the ALS routes (test_als_cov_na.R): one rule, one note

watch = function(expr) {
	msgs = character(); warns = character()
	val = withCallingHandlers(
		expr,
		message = function(m) {
			msgs <<- c(msgs, conditionMessage(m))
			invokeRestart("muffleMessage")
		},
		warning = function(w) {
			warns <<- c(warns, conditionMessage(w))
			invokeRestart("muffleWarning")
		}
	)
	list(val = val, msgs = msgs, warns = warns)
}
cov_note = "had a missing covariate and (is|are) treated as unobserved"
fill_note = "missing in one triangle (was|were) filled from the other triangle"
fill_count = function(msgs) {
	hit = grep(fill_note, msgs, value = TRUE)
	if (!length(hit)) return(NA_integer_)
	as.integer(sub("^.*?([0-9]+) `Xdyad` value.*$", "\\1", gsub("\n", " ", hit[1L])))
}
note_count = function(msgs) {
	hit = grep(cov_note, msgs, value = TRUE)
	if (!length(hit)) return(NA_integer_)
	as.integer(sub("^.*?([0-9]+) `Y` cell.*$", "\\1", gsub("\n", " ", hit[1L])))
}
na_cells = function(Yarr) {
	if (length(dim(Yarr)) == 2L) return(list(unname(is.na(Yarr))))
	lapply(seq_len(dim(Yarr)[3]), function(t) unname(is.na(Yarr[, , t])))
}
quick = list(nscan = 20, burn = 5, odens = 5, verbose = FALSE, plot = FALSE, gof = FALSE)

# small unipartite, bipartite and symmetric panels with the three kinds of
# covariate; the node covariates are constant over time
cov_data = function(kind = c("directed", "bipartite", "symmetric"), Tt = 3, seed = 1) {
	kind = match.arg(kind)
	set.seed(seed)
	bip = kind == "bipartite"
	nr = if (bip) 10L else 12L
	nc = if (bip) 8L else 12L
	rn = sprintf("r%02d", seq_len(nr))
	cn = if (bip) sprintf("c%02d", seq_len(nc)) else rn
	wr = rnorm(nr); wc = rnorm(nc)
	Y = Xd = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(nr * nc), nr, nc, dimnames = list(rn, cn))
		e = matrix(rnorm(nr * nc), nr, nc)
		if (kind == "symmetric") { x = (x + t(x)) / 2; e = (e + t(e)) / sqrt(2) }
		y = 1 * (-0.4 + 0.7 * x + e > 0)
		if (!bip) diag(y) = NA
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		Xd[[t]] = array(x, c(nr, nc, 1), dimnames = list(rn, cn, "xd"))
	}
	Xr = lapply(seq_len(Tt), function(t) matrix(wr, nr, 1, dimnames = list(rn, "xr")))
	Xc = if (kind == "symmetric") NULL else
		lapply(seq_len(Tt), function(t) matrix(wc, nc, 1, dimnames = list(cn, "xc")))
	list(Y = Y, Xd = Xd, Xr = Xr, Xc = Xc, kind = kind, bip = bip, symmetric = kind == "symmetric")
}
# missing values in each kind of covariate (one triangle only for the symmetric
# dyadic covariate)
cov_hole = function(d, nodes = TRUE) {
	Tt = length(d$Y)
	d$Xd[[1]][2, 5, 1] = NA
	d$Xd[[min(2L, Tt)]][7, 3, 1] = NA
	if (nodes) {
		d$Xr[[min(2L, Tt)]][4, 1] = NA
		if (!is.null(d$Xc)) d$Xc[[min(3L, Tt)]][6, 1] = NA
	}
	d
}
# the cells that are out when a cell with a missing covariate is unobserved; a
# symmetric dyad has one covariate value, taken from the other triangle when
# one cell is missing, so it is unobserved only when both cells are
oracle = function(d) {
	lapply(seq_along(d$Y), function(t) {
		m = is.na(d$Y[[t]])
		cm = apply(is.na(d$Xd[[t]]), c(1, 2), any)
		if (d$symmetric) cm = cm & t(cm)
		m = m | cm
		m[apply(is.na(d$Xr[[t]]), 1, any), ] = TRUE
		if (d$symmetric) m[, apply(is.na(d$Xr[[t]]), 1, any)] = TRUE
		if (!is.null(d$Xc)) m[, apply(is.na(d$Xc[[t]]), 1, any)] = TRUE
		unname(m)
	})
}
# the data with those cells hidden and the missing covariate values filled in;
# a symmetric covariate takes the value of the other triangle first
hidden = function(d) {
	o = oracle(d)
	h = d
	for (t in seq_along(d$Y)) {
		y = d$Y[[t]]; y[o[[t]]] = NA; h$Y[[t]] = y
		if (d$symmetric) {
			x = h$Xd[[t]]
			flip = aperm(x, c(2, 1, 3))
			x[is.na(x)] = flip[is.na(x)]
			h$Xd[[t]] = x
		}
		h$Xd[[t]][is.na(h$Xd[[t]])] = 0
		h$Xr[[t]][is.na(h$Xr[[t]])] = 0
		if (!is.null(h$Xc)) h$Xc[[t]][is.na(h$Xc[[t]])] = 0
	}
	h
}
fit_ame = function(d, ...) {
	args = c(
		list(
			d$Y[[1]],
			Xdyad = d$Xd[[1]],
			Xrow = d$Xr[[1]],
			family = "binary",
			R = 0,
			g = 100,
			mode = if (d$bip) "bipartite" else "unipartite",
			symmetric = d$symmetric
		),
		quick,
		list(...)
	)
	if (!is.null(d$Xc)) args$Xcol = d$Xc[[1]]
	do.call(ame, args)
}
fit_lame = function(d, ...) {
	args = c(
		list(
			d$Y,
			Xdyad = d$Xd,
			Xrow = d$Xr,
			family = "binary",
			R = 0,
			mode = if (d$bip) "bipartite" else "unipartite",
			symmetric = d$symmetric
		),
		quick,
		list(...)
	)
	if (!is.null(d$Xc)) args$Xcol = d$Xc
	do.call(lame, args)
}

test_that("the propagation helper handles a rectangular network and a symmetric dyad", {
	set.seed(1)
	y = matrix(rbinom(6 * 4, 1, 0.5), 6, 4)
	xd = matrix(rnorm(6 * 4), 6, 4)
	xd[2, 3] = NA
	# a 2-D dyadic covariate takes the shape of Y
	expect_identical(which(is.na(lame:::.ame_propagate_cov_na(y, Xdyad = xd))), 14L)
	xr = matrix(c(1, NA, 3, 4, 5, 6), 6, 1)
	xc = matrix(c(1, 2, NA, 4), 4, 1)
	out = lame:::.ame_propagate_cov_na(y, Xrow = xr, Xcol = xc)
	expect_true(all(is.na(out[2, ])) && all(is.na(out[, 3])))
	expect_equal(sum(is.na(out)), 4 + 6 - 1)
	# a square matrix: only the cell with the missing covariate drops
	ys = matrix(1, 5, 5); diag(ys) = NA
	xs = matrix(0, 5, 5); xs[2, 4] = NA
	one = lame:::.ame_propagate_cov_na(ys, Xdyad = xs)
	expect_true(is.na(one[2, 4]) && !is.na(one[4, 2]))
	# a symmetric dyad has one covariate value
	both = lame:::.ame_propagate_cov_na(ys, Xdyad = xs, symmetric = TRUE)
	expect_true(is.na(both[2, 4]) && is.na(both[4, 2]))
	expect_equal(sum(is.na(both)), 5 + 2)
	# complete covariates leave Y as it is
	expect_identical(lame:::.ame_propagate_cov_na(ys, Xdyad = matrix(0, 5, 5), Xrow = matrix(1, 5, 1)), ys)
})

test_that("a panel design hides the whole symmetric dyad", {
	n = 5; Tt = 2
	nm = letters[seq_len(n)]
	Y = array(1, c(n, n, Tt), dimnames = list(nm, nm, c("t1", "t2")))
	for (t in 1:Tt) diag(Y[, , t]) = NA
	X = array(rnorm(n * n * Tt), c(n, n, 1, Tt), dimnames = list(nm, nm, "xd", c("t1", "t2")))
	X[2, 4, 1, 1] = NA
	one = get_design_rep(Y, X, NULL, NULL, nm, TRUE, n, Tt, 0L, 0L, 1L)
	expect_true(is.na(one$Y[2, 4, 1]) && !is.na(one$Y[4, 2, 1]))
	both = get_design_rep(Y, X, NULL, NULL, nm, TRUE, n, Tt, 0L, 0L, 1L, symmetric = TRUE)
	expect_true(is.na(both$Y[2, 4, 1]) && is.na(both$Y[4, 2, 1]))
	expect_equal(sum(is.na(both$Y[, , 2])), n)
	# the covariate of an unobserved cell is 0
	expect_equal(both$X[2, 4, 2, 1], 0)
	expect_false(anyNA(both$X))
})

test_that("the note gives the number of cells once, in one wording", {
	expect_message(lame:::.lame_inform_cov_na(3L), "3 `Y` cells had a missing covariate and are treated as unobserved")
	expect_message(lame:::.lame_inform_cov_na(1L), "1 `Y` cell had a missing covariate and is treated as unobserved")
	expect_no_message(lame:::.lame_inform_cov_na(0L))
	op = options(lame.warn.na = FALSE)
	on.exit(options(op))
	expect_no_message(lame:::.lame_inform_cov_na(5L))
})

test_that("ame() on a bipartite network drops the cells with a missing covariate", {
	d = cov_hole(cov_data("bipartite", Tt = 1))
	r = watch(fit_ame(d))
	o = oracle(d)
	expect_identical(na_cells(r$val$Y), o)
	expect_identical(note_count(r$msgs), as.integer(sum(o[[1]])))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	# no data-augmentation warning
	expect_false(any(grepl("data-augmented", r$warns)))
	expect_false(anyNA(r$val$X))
	# the fit is that of the data with those cells hidden and the values filled
	# in (dyadic holes only: ame() stops on a bipartite row or column that the
	# user leaves with no observed cell)
	d1 = cov_hole(cov_data("bipartite", Tt = 1), nodes = FALSE)
	same = suppressMessages(fit_ame(d1))
	hand = suppressMessages(fit_ame(hidden(d1)))
	expect_equal(same$BETA, hand$BETA)
	expect_equal(same$VC, hand$VC)
})

test_that("lame() gives the note once, with the count, for unipartite and bipartite panels", {
	for (kind in c("directed", "bipartite")) {
		d = cov_hole(cov_data(kind, Tt = 3))
		r = watch(fit_lame(d))
		o = oracle(d)
		expect_identical(na_cells(r$val$Y), o)
		expect_identical(
			note_count(r$msgs),
			as.integer(sum(vapply(o, sum, numeric(1))) -
					sum(vapply(d$Y, function(y) sum(is.na(y)), numeric(1))))
		)
		expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
		expect_false(any(grepl("data-augmented", r$warns)))
	}
	# complete covariates give no note
	r = watch(fit_lame(cov_data("directed", Tt = 2)))
	expect_equal(sum(grepl(cov_note, r$msgs)), 0L)
})

test_that("a symmetric fit fills a hole in one triangle and hides a dyad missing in both", {
	d = cov_hole(cov_data("symmetric", Tt = 3))
	o = oracle(d)
	# a hole in one triangle of the dyadic covariate takes the value of the other
	# triangle: the dyad stays observed. ame() on the first period
	r = watch(fit_ame(d))
	expect_identical(na_cells(r$val$Y), o[1])
	expect_false(is.na(r$val$Y[2, 5]) || is.na(r$val$Y[5, 2]))
	expect_equal(fill_count(r$msgs), 1L)
	expect_equal(sum(grepl(cov_note, r$msgs)), 0L)
	# the fit is the fit of the data with the hole filled by hand
	hand = suppressWarnings(suppressMessages(fit_ame(hidden(d))))
	expect_equal(r$val$BETA, hand$BETA)
	expect_equal(r$val$VC, hand$VC)
	# lame() on the panel: one note for the two filled values, one for the actor
	rl = watch(fit_lame(d))
	expect_identical(na_cells(rl$val$Y), o)
	expect_equal(fill_count(rl$msgs), 2L)
	expect_equal(sum(grepl(fill_note, rl$msgs)), 1L)
	expect_equal(sum(grepl(cov_note, rl$msgs)), 1L)
	expect_false(is.na(rl$val$Y[2, 5, 1]) || is.na(rl$val$Y[7, 3, 2]))
	# the cells in the row and the column of the actor with a missing covariate
	expect_true(all(is.na(rl$val$Y[4, -4, 2])) && all(is.na(rl$val$Y[-4, 4, 2])))
	handl = suppressWarnings(suppressMessages(fit_lame(hidden(d))))
	expect_equal(rl$val$BETA, handl$BETA)
	# a dyad missing in both triangles is unobserved, and nothing is filled
	ds = d
	ds$Xd[[1]][5, 2, 1] = NA
	ds$Xd[[2]][3, 7, 1] = NA
	os = oracle(ds)
	expect_false(identical(os, o))
	rs = watch(fit_ame(ds))
	expect_identical(na_cells(rs$val$Y), os[1])
	expect_true(is.na(rs$val$Y[2, 5]) && is.na(rs$val$Y[5, 2]))
	expect_equal(sum(grepl(fill_note, rs$msgs)), 0L)
	expect_equal(sum(grepl(cov_note, rs$msgs)), 1L)
	hand = suppressWarnings(suppressMessages(fit_ame(hidden(ds))))
	expect_equal(rs$val$BETA, hand$BETA)
	rsl = watch(fit_lame(ds))
	expect_identical(na_cells(rsl$val$Y), os)
	expect_equal(sum(grepl(fill_note, rsl$msgs)), 0L)
	handl = suppressWarnings(suppressMessages(fit_lame(hidden(ds))))
	expect_equal(rsl$val$BETA, handl$BETA)
	# a slice with both kinds of hole: the dyad missing in both is out, the other is filled
	dm = d
	dm$Xd[[1]][5, 2, 1] = NA
	dm$Xd[[1]][1, 9, 1] = NA
	rmix = watch(fit_ame(dm))
	expect_identical(na_cells(rmix$val$Y), oracle(dm)[1])
	expect_equal(fill_count(rmix$msgs), 1L)
	expect_true(is.na(rmix$val$Y[2, 5]) && !is.na(rmix$val$Y[1, 9]) && !is.na(rmix$val$Y[9, 1]))
	# the options switch silences the note, not the fill
	op = options(lame.warn.na = FALSE)
	on.exit(options(op))
	quiet = watch(fit_ame(d))
	expect_equal(length(quiet$msgs), 0L)
	expect_equal(quiet$val$BETA, r$val$BETA)
	expect_identical(na_cells(quiet$val$Y), o[1])
})

test_that("the fill helper takes a hole from the other triangle and leaves the rest as it is", {
	set.seed(3)
	n = 6
	x = array(rnorm(n * n * 2), c(n, n, 2), dimnames = list(letters[1:n], letters[1:n], c("u", "v")))
	x[, , 1] = (x[, , 1] + t(x[, , 1])) / 2
	x[, , 2] = (x[, , 2] + t(x[, , 2])) / 2
	full = x
	x[2, 5, 1] = NA
	x[4, 1, 2] = NA; x[1, 4, 2] = NA
	x[3, 3, 2] = NA
	op = options(lame.warn.na = FALSE)
	on.exit(options(op))
	out = lame:::.lame_fill_symmetric_dyadic(x)
	expect_equal(out[2, 5, 1], full[5, 2, 1])
	expect_equal(out[, , 1], t(out[, , 1]))
	# missing in both triangles: unchanged; a missing diagonal value: unchanged
	expect_true(is.na(out[4, 1, 2]) && is.na(out[1, 4, 2]) && is.na(out[3, 3, 2]))
	expect_identical(sum(is.na(out)), 3L)
	expect_identical(dimnames(out), dimnames(x))
	# a panel is filled period by period
	pan = array(NA_real_, c(n, n, 2, 3))
	for (tt in 1:3) pan[, , , tt] = full * tt
	pan[2, 5, 1, 2] = NA
	pan[6, 1, 2, 3] = NA
	pan_out = lame:::.lame_fill_symmetric_dyadic(pan)
	expect_equal(pan_out[2, 5, 1, 2], full[5, 2, 1] * 2)
	expect_equal(pan_out[6, 1, 2, 3], full[1, 6, 2] * 3)
	expect_false(anyNA(pan_out))
	# a covariate without a square leading block is returned as it is
	expect_identical(
		lame:::.lame_fill_symmetric_dyadic(matrix(c(1, NA, 3, 4, 5, 6), 2)),
		matrix(c(1, NA, 3, 4, 5, 6), 2)
	)
	# the note counts the values filled, once
	options(lame.warn.na = TRUE)
	expect_message(
		lame:::.lame_fill_symmetric_dyadic(x),
		"1 `Xdyad` value missing in one triangle was filled from the other triangle",
		class = "lame_symmetric_xdyad_fill"
	)
	expect_no_message(lame:::.lame_fill_symmetric_dyadic(full))
})

test_that("ame(), lame() and the ALS routes leave the same cells unobserved", {
	skip_on_cran()
	for (kind in c("directed", "bipartite", "symmetric")) {
		d = cov_hole(cov_data(kind, Tt = 3))
		mode = if (d$bip) "bipartite" else "unipartite"
		# a cross-section
		mc = watch(fit_ame(d))
		al_args = list(
			d$Y[[1]],
			Xdyad = d$Xd[[1]],
			Xrow = d$Xr[[1]],
			family = "binary",
			R = 0,
			mode = mode,
			symmetric = d$symmetric,
			verbose = FALSE
		)
		if (!is.null(d$Xc)) al_args$Xcol = d$Xc[[1]]
		al = watch(do.call(ame_als, al_args))
		expect_identical(na_cells(al$val$Y), na_cells(mc$val$Y), info = kind)
		expect_identical(note_count(al$msgs), note_count(mc$msgs), info = kind)
		expect_equal(length(mc$warns[grepl("covariate", mc$warns)]), 0L, info = kind)
		# a panel
		mc = watch(fit_lame(d))
		al_args = list(
			d$Y,
			Xdyad = d$Xd,
			Xrow = d$Xr,
			family = "binary",
			R = 0,
			mode = mode,
			symmetric = d$symmetric,
			verbose = FALSE
		)
		if (!is.null(d$Xc)) al_args$Xcol = d$Xc
		al = watch(do.call(lame_als, al_args))
		expect_identical(na_cells(al$val$Y), na_cells(mc$val$Y), info = kind)
		expect_identical(note_count(al$msgs), note_count(mc$msgs), info = kind)
	}
})

test_that("the slope of a bipartite ame() fit is not biased by a missing covariate", {
	skip_on_cran()
	nA = 20; nB = 15
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	set.seed(2)
	x = matrix(rnorm(nA * nB, 3, 1), nA, nB, dimnames = list(rn, cn))
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	y = matrix(
		1 * (-2.4 + 0.7 * x + outer(a, b, "+") + rnorm(nA * nB) > 0),
		nA,
		nB,
		dimnames = list(rn, cn)
	)
	hide = matrix(runif(nA * nB) < 0.2, nA)
	xm = x; xm[hide] = NA
	fit = function(X) suppressMessages(ame(
		y,
		Xdyad = array(X, c(nA, nB, 1), dimnames = list(rn, cn, "x")),
		family = "binary",
		mode = "bipartite",
		R = 0,
		nscan = 1500,
		burn = 300,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	fit_complete = fit(x)
	fit_missing = fit(xm)
	slope = function(f) mean(f$BETA[, "x_dyad"])
	complete = slope(fit_complete)
	missing = slope(fit_missing)
	# dropping the cells keeps the slope near the complete-data fit (true slope 0.7)
	expect_lt(abs(missing - complete), 0.2)
	expect_gt(missing, 0.4)
	# both fits recover the true slope: it lies inside their 95% credible intervals
	for (f in list(fit_complete, fit_missing)) {
		ci = quantile(f$BETA[, "x_dyad"], c(0.025, 0.975))
		expect_gt(0.7, ci[[1]])
		expect_lt(0.7, ci[[2]])
	}
})

test_that("the log-likelihood scoring skips the cells a missing covariate removed", {
	skip_on_cran()
	for (kind in c("directed", "bipartite")) {
		d = cov_hole(cov_data(kind, Tt = 1))
		fit = suppressMessages(fit_ame(d, save_log_lik = TRUE))
		expect_equal(ncol(fit$log_lik), sum(!is.na(fit$Y)))
		expect_equal(nrow(fit$log_lik_index), sum(!is.na(fit$Y)))
	}
})
