skip_on_cran()

library(lame)

# a cell whose dyadic, sender or receiver covariate is missing in its period is
# unobserved in every ALS route, as in ame() and lame(): the covariate
# coefficients are estimated from the complete cells, and a missing value is
# not modelled as 0

# a directed panel with a dyadic, a sender and a receiver covariate, complete.
# the node covariates are constant over time unless varying = TRUE
cov_panel = function(n = 14, Tt = 3, seed = 1, varying = FALSE) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	wr = rnorm(n); wc = rnorm(n)
	Y = Xd = Xr = Xc = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(n * n), n, dimnames = list(nm, nm))
		y = 1 * (-0.5 + 0.7 * x + matrix(rnorm(n * n), n) > 0)
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		Y[[t]] = y
		Xd[[t]] = array(x, c(n, n, 1), dimnames = list(nm, nm, "xd"))
		Xr[[t]] = matrix(
			if (varying) wr + rnorm(n, 0, 0.3) else wr,
			n,
			1,
			dimnames = list(nm, "xr")
		)
		Xc[[t]] = matrix(
			if (varying) wc + rnorm(n, 0, 0.3) else wc,
			n,
			1,
			dimnames = list(nm, "xc")
		)
	}
	list(Y = Y, Xd = Xd, Xr = Xr, Xc = Xc, nm = nm, n = n, Tt = Tt)
}

# the cells the MCMC rule leaves unobserved, written out independently of the
# package: a cell is out when its outcome, its dyadic covariate, its sender's
# covariate or its receiver's covariate is missing
cov_na_oracle = function(Y, Xd, Xr, Xc) {
	lapply(seq_along(Y), function(t) {
		m = is.na(Y[[t]])
		if (!is.null(Xd)) m = m | apply(is.na(Xd[[t]]), c(1, 2), any)
		if (!is.null(Xr)) m[apply(is.na(Xr[[t]]), 1, any), ] = TRUE
		if (!is.null(Xc)) m[, apply(is.na(Xc[[t]]), 1, any)] = TRUE
		unname(m)
	})
}
na_cells = function(Yarr) {
	lapply(seq_len(dim(Yarr)[3]), function(t) unname(is.na(Yarr[, , t])))
}
n_out = function(oracle, Y) {
	as.integer(sum(vapply(oracle, sum, numeric(1))) - sum(vapply(Y, function(y) sum(is.na(y)), numeric(1))))
}

# the messages and warnings of an expression, and its value
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
note_count = function(msgs) {
	hit = grep(cov_note, msgs, value = TRUE)
	if (!length(hit)) return(NA_integer_)
	as.integer(sub("^.*?([0-9]+) `Y` cell.*$", "\\1", gsub("\n", " ", hit[1L])))
}

# static data with missing values in all three kinds of covariate
cov_static = function(seed = 2) {
	d = cov_panel(n = 14, Tt = 1, seed = seed)
	xd = d$Xd[[1]]; xr = d$Xr[[1]]; xc = d$Xc[[1]]
	xd[2, 5, 1] = NA; xd[7, 3, 1] = NA
	xr[4, 1] = NA
	xc[9, 1] = NA
	list(y = d$Y[[1]], xd = xd, xr = xr, xc = xc)
}

test_that("a cross-section: a cell with a missing covariate is unobserved", {
	d = cov_static()
	r = watch(ame_als(
		d$y,
		Xdyad = d$xd,
		Xrow = d$xr,
		Xcol = d$xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	# the cells ame() sets to NA (.ame_propagate_cov_na), and the note it gives
	expected = is.na(lame:::.ame_propagate_cov_na(d$y, d$xr, d$xc, d$xd))
	expect_identical(unname(is.na(r$val$Y[, , 1])), unname(expected))
	expect_identical(note_count(r$msgs), as.integer(sum(expected) - sum(is.na(d$y))))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	# no zero-fill warning
	expect_length(r$warns, 0L)
	expect_true(r$val$converged)
	expect_true(all(is.finite(coef(r$val))))
	# the fit is that of the data with those cells hidden and the values filled in
	xz = d$xd; xz[is.na(xz)] = 0
	xrz = d$xr; xrz[is.na(xrz)] = 0
	xcz = d$xc; xcz[is.na(xcz)] = 0
	yh = d$y; yh[expected] = NA
	hand = suppressMessages(ame_als(
		yh,
		Xdyad = xz,
		Xrow = xrz,
		Xcol = xcz,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_equal(coef(r$val), coef(hand), tolerance = 1e-8)
	# a call that names the argument gets the same
	expect_message(
		ame_als(d$y, Xdyad = d$xd, family = "binary", R = 0, verbose = FALSE),
		cov_note
	)
})

test_that("a panel: the cells of a period with a missing covariate are unobserved", {
	d = cov_panel(n = 14, Tt = 3, seed = 3)
	Xd = d$Xd; Xr = d$Xr; Xc = d$Xc
	Xd[[1]][2, 5, 1] = NA; Xd[[2]][7, 3, 1] = NA
	Xr[[2]][4, 1] = NA
	Xc[[3]][9, 1] = NA
	oracle = cov_na_oracle(d$Y, Xd, Xr, Xc)
	r = watch(lame_als(
		d$Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_identical(na_cells(r$val$Y), oracle)
	expect_identical(note_count(r$msgs), n_out(oracle, d$Y))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	# the node covariates vary over time here, which is a different warning
	expect_false(any(grepl("replaced with 0", r$warns)))
	# the same route through lame()
	r2 = watch(lame(
		d$Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		method = "als",
		verbose = FALSE
	))
	expect_identical(na_cells(r2$val$Y), oracle)
	expect_equal(sum(grepl(cov_note, r2$msgs)), 1L)
})

test_that("the per-actor mean of a node covariate ignores the periods in which it is missing", {
	d = cov_panel(n = 14, Tt = 3, seed = 4, varying = TRUE)
	Xr = d$Xr; Xc = d$Xc
	Xr[[2]][4, 1] = NA
	Xc[[1]][9, 1] = NA; Xc[[3]][9, 1] = NA
	for (t in 1:3) Xr[[t]][6, 1] = NA # actor 6 is missing in every period
	prep = lame:::.ae_prepare(d$Y, d$Xd, Xr, Xc, "unipartite", TRUE, "quiet")
	rows = cbind(Xr[[1]][, 1], Xr[[2]][, 1], Xr[[3]][, 1])
	expect_equal(prep$W_row[4, 1], mean(rows[4, ], na.rm = TRUE))
	expect_equal(unname(prep$W_row[-6, 1]), unname(rowMeans(rows, na.rm = TRUE)[-6]))
	expect_equal(prep$W_col[9, 1], d$Xc[[2]][9, 1])
	# the periods in which a value is missing hold the mean, so they add no variation
	expect_equal(
		unname(prep$W_row_arr[4, 1, ]),
		c(Xr[[1]][4, 1], prep$W_row[4, 1], Xr[[3]][4, 1])
	)
	# the cells of those periods are unobserved
	expect_true(all(is.na(prep$Yarr[4, , 2])))
	expect_true(all(is.na(prep$Yarr[, 9, c(1, 3)])))
	expect_true(all(is.na(prep$Yarr[6, , ])))
	expect_false(all(is.na(prep$Yarr[4, , c(1, 3)])))
	# an actor missing in every period has no sender cell left
	expect_identical(prep$cov_na$row, 6L)
	expect_identical(prep$cov_na$col, integer(0))
})

test_that("an actor without a covariate in any period has no effect to estimate", {
	d = cov_panel(n = 16, Tt = 3, seed = 5)
	Xr = d$Xr; Xc = d$Xc
	for (t in 1:3) { Xr[[t]][6, 1] = NA; Xc[[t]][6, 1] = NA }
	# the empirical-Bayes prior holds its effects at the mean, with no note on
	# the actor and no disconnected-graph warning
	r = watch(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_identical(r$val$bias_correct$method, "eb")
	expect_true(all(is.na(r$val$Y[6, , ])) && all(is.na(r$val$Y[, 6, ])))
	expect_identical(unname(r$val$a[6]), 0)
	expect_identical(unname(r$val$b[6]), 0)
	expect_null(r$val$not_estimated)
	expect_false(any(grepl("disconnected", r$warns)))
	expect_true(all(is.finite(r$val$a)) && all(is.finite(r$val$b)))
	# a fixed-effects fit reports the effects as not estimated
	r = watch(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE,
		bias_correct = "none"
	))
	fe = r$val
	expect_true(is.na(fe$a[["a06"]]) && is.na(fe$b[["a06"]]))
	expect_true(all(is.finite(fe$a[-6])) && all(is.finite(fe$b[-6])))
	expect_identical(sort(fe$not_estimated$side), c("receiver", "sender"))
	expect_true(all(fe$not_estimated$actor == "a06"))
	expect_true(all(fe$not_estimated$outcome == "no observed cell"))
	expect_true(any(grepl("a06 \\(sender, no observed cell\\)", r$msgs)))
	expect_false(any(grepl("disconnected", r$warns)))
	expect_true(all(is.na(fe$fitted[[1]][6, ])) && all(is.na(fe$fitted[[1]][, 6])))
	# the variance components describe the actors that are estimated
	expect_true(all(is.finite(fe$VC[c("va", "vb")])))
	# a normal fit does the same: the empirical-Bayes default holds the effects at the mean, the fixed-effects fit reports them as not estimated
	r = watch(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "normal",
		R = 0,
		verbose = FALSE
	))
	expect_identical(r$val$bias_correct$method, "eb")
	expect_identical(unname(r$val$a[6]), 0)
	expect_identical(unname(r$val$b[6]), 0)
	expect_null(r$val$not_estimated)
	r = watch(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "normal",
		R = 0,
		verbose = FALSE,
		bias_correct = "none"
	))
	expect_true(is.na(r$val$a[["a06"]]) && is.na(r$val$b[["a06"]]))
	# only the sender side missing: the receiver effect stays
	Xr1 = d$Xr; for (t in 1:3) Xr1[[t]][6, 1] = NA
	fe1 = suppressMessages(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr1,
		Xcol = d$Xc,
		family = "binary",
		R = 0,
		verbose = FALSE,
		bias_correct = "none"
	))
	expect_true(is.na(fe1$a[["a06"]]))
	expect_true(is.finite(fe1$b[["a06"]]))
	expect_identical(fe1$not_estimated$side, "sender")
})

test_that("a bootstrap of a fit with unobserved actors keeps them out", {
	skip_on_cran()
	d = cov_panel(n = 16, Tt = 2, seed = 6)
	Xr = d$Xr; Xc = d$Xc
	for (t in 1:2) { Xr[[t]][6, 1] = NA; Xc[[t]][6, 1] = NA }
	fe = suppressMessages(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE,
		bias_correct = "none"
	))
	# an uncorrected binary fixed-effects fit has no bootstrap
	expect_error(
		ame_als_bootstrap(fe, R = 10, type = "parametric", seed = 1, verbose = FALSE),
		class = "lame_als_bootstrap_bias"
	)
	refit = ame_als_refit(fe)
	expect_true(is.na(refit$a[["a06"]]))
	expect_identical(refit$not_estimated$outcome, c("no observed cell", "no observed cell"))
	eb = suppressMessages(lame_als(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	bt = suppressWarnings(ame_als_bootstrap(
		eb,
		R = 10,
		type = "parametric",
		seed = 1,
		verbose = FALSE
	))
	expect_gte(bt$n_valid, 10L)
})

test_that("bipartite: a cell with a missing covariate is unobserved", {
	nA = 10; nB = 8; Tt = 3
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	set.seed(7)
	wr = rnorm(nA); wc = rnorm(nB)
	Y = Xd = Xr = Xc = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(nA * nB), nA, nB)
		Xd[[t]] = array(x, c(nA, nB, 1), dimnames = list(rn, cn, "xd"))
		Xr[[t]] = matrix(wr, nA, 1, dimnames = list(rn, "xr"))
		Xc[[t]] = matrix(wc, nB, 1, dimnames = list(cn, "xc"))
		Y[[t]] = matrix(
			1 * (-0.3 + 0.6 * x + rnorm(nA * nB) > 0),
			nA,
			nB,
			dimnames = list(rn, cn)
		)
	}
	Xd[[1]][2, 5, 1] = NA; Xd[[2]][7, 3, 1] = NA
	Xr[[2]][4, 1] = NA
	Xc[[3]][6, 1] = NA
	oracle = cov_na_oracle(Y, Xd, Xr, Xc)
	r = watch(lame_als(
		Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		mode = "bipartite",
		R = 0,
		verbose = FALSE
	))
	expect_identical(na_cells(r$val$Y), oracle)
	expect_identical(note_count(r$msgs), n_out(oracle, Y))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	expect_false(any(grepl("replaced with 0", r$warns)))
	# a rectangular cross-section
	r = watch(ame_als(
		Y[[1]],
		Xdyad = Xd[[1]],
		Xrow = Xr[[1]],
		Xcol = Xc[[1]],
		family = "binary",
		mode = "bipartite",
		R = 0,
		verbose = FALSE
	))
	expect_identical(unname(is.na(r$val$Y[, , 1])), oracle[[1]])
	# a square bipartite network has no diagonal: a missing value on it is missing
	nS = 9
	rs = sprintf("r%02d", seq_len(nS)); cs = sprintf("c%02d", seq_len(nS))
	set.seed(8)
	ys = matrix(1 * (rnorm(nS * nS) > 0), nS, dimnames = list(rs, cs))
	xs = matrix(rnorm(nS * nS), nS, dimnames = list(rs, cs))
	xs[3, 3] = NA; xs[5, 5] = NA
	r = watch(ame_als(
		ys,
		Xdyad = array(xs, c(nS, nS, 1), dimnames = list(rs, cs, "x")),
		family = "binary",
		mode = "bipartite",
		R = 0,
		verbose = FALSE
	))
	expect_identical(note_count(r$msgs), 2L)
	expect_true(is.na(r$val$Y[3, 3, 1]) && is.na(r$val$Y[5, 5, 1]))
	expect_equal(sum(is.na(r$val$Y)), 2L)
})

test_that("a missing value on the diagonal of a dyadic covariate is not a missing value", {
	d = cov_panel(n = 14, Tt = 2, seed = 9)
	Xd = d$Xd
	for (t in 1:2) { x = Xd[[t]][, , 1]; diag(x) = NA; Xd[[t]][, , 1] = x }
	r = watch(lame_als(d$Y, Xdyad = Xd, family = "binary", R = 0, verbose = FALSE))
	expect_length(r$msgs, 0L)
	expect_length(r$warns, 0L)
	expect_identical(na_cells(r$val$Y), lapply(d$Y, function(y) unname(is.na(y))))
	# the fit is that of the covariate with its diagonal filled
	Xz = Xd
	for (t in 1:2) { x = Xz[[t]][, , 1]; diag(x) = 0; Xz[[t]][, , 1] = x }
	r0 = lame_als(d$Y, Xdyad = Xz, family = "binary", R = 0, verbose = FALSE)
	expect_equal(coef(r$val), coef(r0), tolerance = 1e-10)
})

test_that("a covariate that is missing everywhere stops the fit, as in ame() and lame()", {
	d = cov_panel(n = 12, Tt = 3, seed = 10)
	# a cross-section, as in ame()
	x1 = matrix(NA_real_, 12, 12, dimnames = list(d$nm, d$nm))
	expect_error(
		ame_als(d$Y[[1]], Xdyad = x1, family = "binary", R = 0, verbose = FALSE),
		"`Xdyad` covariate 1 is entirely NA\\."
	)
	xr = matrix(c(rnorm(12), rep(NA, 12)), 12, 2, dimnames = list(d$nm, c("a", "b")))
	expect_error(
		ame_als(d$Y[[1]], Xrow = xr, family = "binary", R = 0, verbose = FALSE),
		"`Xrow` covariate 2 is entirely NA\\."
	)
	expect_error(
		ame_als(d$Y[[1]], Xcol = xr, family = "binary", R = 0, verbose = FALSE),
		"`Xcol` covariate 2 is entirely NA\\."
	)
	expect_error(
		ame_als(d$Y[[1]], Xrow = xr, family = "binary", R = 0, verbose = FALSE),
		"cannot be identified"
	)
	# a panel, as in lame(): missing in every period
	Xd = d$Xd
	for (t in 1:3) Xd[[t]][, , 1] = NA
	expect_error(
		lame_als(d$Y, Xdyad = Xd, family = "binary", R = 0, verbose = FALSE),
		"`Xdyad` covariate 1 is entirely NA in every period\\."
	)
	Xc = d$Xc
	for (t in 1:3) Xc[[t]][, 1] = NA
	expect_error(
		lame_als(d$Y, Xcol = Xc, family = "binary", R = 0, verbose = FALSE),
		"`Xcol` covariate 1 is entirely NA in every period\\."
	)
	expect_error(
		lame(d$Y, Xcol = Xc, family = "binary", R = 0, method = "als", verbose = FALSE),
		"entirely NA in every period"
	)
	# missing in one period only is an ordinary missing value
	Xc2 = d$Xc
	Xc2[[2]][, 1] = NA
	r = watch(lame_als(d$Y, Xcol = Xc2, family = "binary", R = 0, verbose = FALSE))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	expect_true(all(is.na(r$val$Y[, , 2])))
})

test_that("a symmetric fit fills a hole in one triangle, drops a dyad missing in both, and drops the whole actor", {
	n = 14
	nm = sprintf("s%02d", seq_len(n))
	set.seed(11)
	ys = matrix(0, n, n)
	ys[upper.tri(ys)] = rbinom(n * (n - 1) / 2, 1, 0.4)
	ys = ys + t(ys)
	diag(ys) = NA
	dimnames(ys) = list(nm, nm)
	xs = matrix(rnorm(n * n), n, dimnames = list(nm, nm))
	xs = (xs + t(xs)) / 2
	arr = function(x) array(x, c(n, n, 1), dimnames = list(nm, nm, "xd"))
	full = suppressMessages(ame_als(
		ys,
		Xdyad = arr(xs),
		family = "binary",
		R = 0,
		symmetric = TRUE,
		verbose = FALSE
	))
	# a value missing in one triangle only: the dyad has the covariate of the
	# other triangle, so it stays observed and the fit is the fit of the complete
	# covariate
	x1 = xs; x1[2, 5] = NA
	r = watch(ame_als(
		ys,
		Xdyad = arr(x1),
		family = "binary",
		R = 0,
		symmetric = TRUE,
		verbose = FALSE
	))
	expect_false(is.na(r$val$Y[2, 5, 1]) || is.na(r$val$Y[5, 2, 1]))
	expect_identical(note_count(r$msgs), NA_integer_)
	expect_equal(
		sum(grepl(
			"1 `Xdyad` value missing in one triangle was filled from the other triangle",
			r$msgs
		)),
		1L
	)
	expect_equal(coef(r$val), coef(full), tolerance = 1e-10)
	# a panel is filled period by period, in one note
	pan = list(t1 = ys, t2 = ys)
	xp = list(t1 = arr(xs), t2 = arr(xs * 2))
	xph = xp; xph$t1[2, 5, 1] = NA; xph$t2[7, 3, 1] = NA
	fullp = suppressMessages(lame_als(pan, Xdyad = xp, family = "binary", R = 0, symmetric = TRUE, verbose = FALSE))
	rp = watch(lame_als(pan, Xdyad = xph, family = "binary", R = 0, symmetric = TRUE, verbose = FALSE))
	expect_equal(
		sum(grepl(
			"2 `Xdyad` values missing in one triangle were filled from the other triangle",
			rp$msgs
		)),
		1L
	)
	expect_identical(sum(is.na(rp$val$Y)), sum(is.na(fullp$Y)))
	expect_equal(coef(rp$val), coef(fullp), tolerance = 1e-10)
	# in both triangles: the dyad has no covariate
	x2 = xs; x2[2, 5] = NA; x2[5, 2] = NA
	r2 = watch(ame_als(
		ys,
		Xdyad = arr(x2),
		family = "binary",
		R = 0,
		symmetric = TRUE,
		verbose = FALSE
	))
	expect_true(is.na(r2$val$Y[2, 5, 1]) && is.na(r2$val$Y[5, 2, 1]))
	expect_identical(note_count(r2$msgs), 2L)
	expect_equal(sum(grepl("filled from the other triangle", r2$msgs)), 0L)
	# an actor missing a covariate is out of every dyad
	xr = matrix(rnorm(n), n, 1, dimnames = list(nm, "xr"))
	xr[4, 1] = NA
	r3 = watch(ame_als(ys, Xrow = xr, family = "binary", R = 0, symmetric = TRUE, verbose = FALSE))
	expect_true(all(is.na(r3$val$Y[4, -4, 1])) && all(is.na(r3$val$Y[-4, 4, 1])))
	expect_identical(note_count(r3$msgs), as.integer(2 * (n - 1)))
	expect_identical(unname(r3$val$a[4]), 0)
})

test_that("a dynamic ALS fit informs once and drops the same cells", {
	d = cov_panel(n = 14, Tt = 4, seed = 12)
	Xd = d$Xd
	Xd[[2]][3, 7, 1] = NA; Xd[[3]][8, 2, 1] = NA
	oracle = cov_na_oracle(d$Y, Xd, NULL, NULL)
	r = watch(lame(
		d$Y,
		Xdyad = Xd,
		family = "binary",
		R = 0,
		method = "als",
		dynamic_ab = TRUE,
		verbose = FALSE
	))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	expect_identical(note_count(r$msgs), n_out(oracle, d$Y))
	expect_identical(na_cells(r$val$Y), oracle)
	expect_false(any(grepl("replaced with 0", r$warns)))
})

test_that("a dynamic fit with an actor that has no observed cell reaches its certificate", {
	skip_on_cran()
	d = cov_panel(n = 18, Tt = 5, seed = 17)
	Xr = d$Xr; Xc = d$Xc
	for (t in 1:5) { Xr[[t]][6, 1] = NA; Xc[[t]][6, 1] = NA }
	# the static additive block and the dynamic coefficient or positions: the actor has no score and no information
	r = watch(lame(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		method = "als",
		dynamic_beta = c(intercept = TRUE),
		verbose = FALSE
	))
	expect_true(r$val$converged)
	expect_true(all(is.finite(r$val$diagnostics$stationarity$by_block)))
	expect_equal(sum(grepl(cov_note, r$msgs)), 1L)
	expect_false(any(grepl("not estimated", r$msgs)))
	r = watch(lame(
		d$Y,
		Xdyad = d$Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 1,
		method = "als",
		dynamic_uv = TRUE,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_true(r$val$converged)
	# the fit with the actor blanked in Y is the same model
	Yb = d$Y
	for (t in 1:5) { Yb[[t]][6, ] = NA; Yb[[t]][, 6] = NA }
	rb = watch(lame(
		Yb,
		Xdyad = d$Xd,
		Xrow = d$Xr,
		Xcol = d$Xc,
		family = "binary",
		R = 1,
		method = "als",
		dynamic_uv = TRUE,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_equal(unname(coef(r$val)), unname(coef(rb$val)), tolerance = 1e-6)
})

test_that("a changing actor set is not a missing covariate", {
	d = cov_panel(n = 14, Tt = 3, seed = 13)
	keep = list(1:14, c(1:12), c(1:11, 13:14))
	Y = Xd = Xr = Xc = vector("list", 3)
	for (t in 1:3) {
		k = keep[[t]]
		Y[[t]] = d$Y[[t]][k, k]
		Xd[[t]] = d$Xd[[t]][k, k, , drop = FALSE]
		Xr[[t]] = d$Xr[[t]][k, , drop = FALSE]
		Xc[[t]] = d$Xc[[t]][k, , drop = FALSE]
	}
	r = watch(lame_als(
		Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_false(any(grepl(cov_note, r$msgs)))
	expect_false(any(grepl("replaced with 0", r$warns)))
	# a missing value for a present actor is a missing covariate
	Xr[[2]][3, 1] = NA
	r = watch(lame_als(
		Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_identical(note_count(r$msgs), 11L)
})

test_that("zero-filling a covariate that is not centred would bias its coefficient", {
	n = 40
	nm = sprintf("a%02d", seq_len(n))
	set.seed(1)
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	x = matrix(rnorm(n * n, 3, 1), n, dimnames = list(nm, nm))
	y = 1 * (-2.6 + 0.7 * x + outer(a, b, "+") + matrix(rnorm(n * n), n) > 0)
	diag(y) = NA
	dimnames(y) = list(nm, nm)
	hide = matrix(runif(n * n) < 0.15, n)
	diag(hide) = FALSE
	xm = x; xm[hide] = NA
	arr = function(m) array(m, c(n, n, 1), dimnames = list(nm, nm, "x"))
	fit = function(...) suppressMessages(ame_als(..., family = "binary", R = 0, verbose = FALSE))
	slope = function(f) unname(coef(f)["x_dyad"])
	b_new = slope(fit(y, Xdyad = arr(xm)))
	b_full = slope(fit(y, Xdyad = arr(x)))
	# the missing values modelled as 0 with the outcomes kept
	xz = xm; xz[is.na(xz)] = 0
	b_zero = slope(fit(y, Xdyad = arr(xz)))
	expect_lt(b_zero, 0.5 * b_full)
	expect_lt(abs(b_new - b_full), 0.1)
	expect_lt(abs(b_new - 0.7), 0.15)
	# the fit is that of the complete cells
	yh = y; yh[hide] = NA
	expect_equal(b_new, slope(fit(yh, Xdyad = arr(xz))), tolerance = 1e-8)
})

test_that("ALS and MCMC leave the same cells unobserved and say the same", {
	skip_on_cran()
	# a cross-section: ame()
	d = cov_static(seed = 14)
	mc = watch(ame(
		d$y,
		Xdyad = d$xd,
		Xrow = d$xr,
		Xcol = d$xc,
		family = "binary",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	al = watch(ame_als(
		d$y,
		Xdyad = d$xd,
		Xrow = d$xr,
		Xcol = d$xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_identical(unname(is.na(al$val$Y[, , 1])), unname(is.na(mc$val$Y)))
	expect_equal(sum(grepl(cov_note, mc$msgs)), 1L)
	expect_equal(sum(grepl(cov_note, al$msgs)), 1L)
	expect_identical(note_count(al$msgs), note_count(mc$msgs))
	# a panel: lame()
	p = cov_panel(n = 12, Tt = 3, seed = 15)
	Xd = p$Xd; Xr = p$Xr; Xc = p$Xc
	Xd[[1]][2, 5, 1] = NA; Xd[[2]][7, 3, 1] = NA
	Xr[[2]][4, 1] = NA
	Xc[[3]][9, 1] = NA
	mc = suppressMessages(lame(
		p$Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	))
	al = suppressMessages(lame_als(
		p$Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		R = 0,
		verbose = FALSE
	))
	expect_identical(na_cells(al$Y), na_cells(mc$Y))
	# a bipartite panel: lame()
	nA = 10; nB = 8
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	set.seed(16)
	wr = rnorm(nA); wc = rnorm(nB)
	Y = Xd = Xr = Xc = vector("list", 3)
	for (t in 1:3) {
		x = matrix(rnorm(nA * nB), nA, nB)
		Xd[[t]] = array(x, c(nA, nB, 1), dimnames = list(rn, cn, "xd"))
		Xr[[t]] = matrix(wr, nA, 1, dimnames = list(rn, "xr"))
		Xc[[t]] = matrix(wc, nB, 1, dimnames = list(cn, "xc"))
		Y[[t]] = matrix(1 * (-0.3 + 0.6 * x + rnorm(nA * nB) > 0), nA, nB, dimnames = list(rn, cn))
	}
	Xd[[1]][2, 5, 1] = NA; Xd[[2]][7, 3, 1] = NA
	Xr[[2]][4, 1] = NA
	Xc[[3]][6, 1] = NA
	mc = suppressWarnings(suppressMessages(lame(
		Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		mode = "bipartite",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)))
	al = suppressMessages(lame_als(
		Y,
		Xdyad = Xd,
		Xrow = Xr,
		Xcol = Xc,
		family = "binary",
		mode = "bipartite",
		R = 0,
		verbose = FALSE
	))
	expect_identical(na_cells(al$Y), na_cells(mc$Y))
})
