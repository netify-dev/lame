skip_on_cran()

library(lame)

test_that("the rank-family prior scale is the binary probit-moment scale", {
	set.seed(3)
	n = 20
	z = matrix(rnorm(n * n), n, n)
	diag(z) = NA
	yb = 1 * (z > 0.8)
	sc = lame:::.lame_probit_ab_scale(yb, "ordinal")
	ybar = mean(yb, na.rm = TRUE)
	E = (yb - ybar) / dnorm(qnorm(ybar))
	diag(E) = NA
	ref = mean(c(var(rowMeans(E, na.rm = TRUE)), var(colMeans(E, na.rm = TRUE))))
	# two categories: the multi-threshold slope reduces to the binary one
	expect_equal(sc$vscale, ref)
	expect_equal(sc$vdfmlt, 1)
	# cbin / frn read nominations as ties and inflate the prior df as binary does
	sc_c = lame:::.lame_probit_ab_scale(yb, "cbin")
	expect_equal(sc_c$vscale, ref)
	expect_gte(sc_c$vdfmlt, 1)
})

test_that("the probit-moment scale tracks the latent effect variance", {
	# rank-normal scores estimate E[z | y] and shrink the actor effects; the
	# probit-moment scale keeps them on the latent scale up to the spread of
	# the other terms, as the binary family's scale does
	set.seed(11)
	n = 60
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	z = outer(a, b, "+") + matrix(rnorm(n * n), n, n)
	diag(z) = NA
	y = matrix(cut(z, c(-Inf, -0.8, 0, 0.9, Inf), labels = FALSE), n, n)
	sc = lame:::.lame_probit_ab_scale(y, "ordinal")
	zs = array(lame:::zscores(y), dim(y))
	rank_scale = mean(c(var(rowMeans(zs, na.rm = TRUE)), var(colMeans(zs, na.rm = TRUE))))
	latent = mean(c(var(rowMeans(z, na.rm = TRUE)), var(colMeans(z, na.rm = TRUE))))
	expect_gt(sc$vscale, 1.4 * rank_scale)
	expect_gt(sc$vscale, 0.5 * latent)
	expect_lt(sc$vscale, 1.2 * latent)
	# a 3-D array and a list of slices give the same scale
	arr = array(c(y, y), c(n, n, 2))
	expect_equal(
		lame:::.lame_probit_ab_scale(arr, "ordinal")$vscale,
		lame:::.lame_probit_ab_scale(list(y, y), "ordinal")$vscale
	)
	# degenerate input falls back to the unit scale
	expect_equal(lame:::.lame_probit_ab_scale(matrix(1, 3, 3), "ordinal")$vscale, 1)
})

# a ranked-nomination panel from the model: each row nominates its top `cap`
# cells among those with a positive latent value
.capped_nomination_panel = function(n = 30, TT = 4, cap = 5, mu = -1.3, sd_ab = 0.5, seed = 22) {
	set.seed(seed)
	a = rnorm(n, 0, sd_ab); b = rnorm(n, 0, sd_ab)
	lapply(seq_len(TT), function(t) {
		z = mu + outer(a, b, "+") + matrix(rnorm(n * n, 0, sqrt(1.49)), n, n)
		diag(z) = -Inf
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			top = top[z[i, top] > 0]
			y[i, top] = 1
		}
		diag(y) = NA
		y
	})
}

test_that("the cbin / frn prior scale is read off the rows below the cap", {
	Y = .capped_nomination_panel()
	# about half the row-periods at the cap
	share = lame:::.lame_odmax_saturated_share(Y, 5, "cbin")
	expect_gt(share, 0.45)
	expect_lt(share, 0.6)
	all_rows = lame:::.lame_probit_ab_scale(Y, "cbin")
	below = lame:::.lame_probit_ab_scale(Y, "cbin", odmax = 5)
	# the rows at the cap share one censored row mean and pull the scale down
	expect_gte(below$vscale, all_rows$vscale)
	# a per-row cap and a 3-d array read the same rows
	expect_equal(lame:::.lame_probit_ab_scale(Y, "cbin", odmax = rep(5, 30)), below)
	expect_equal(lame:::.lame_probit_ab_scale(array(unlist(Y), c(30, 30, 4)), "cbin", odmax = 5), below)
	# frn reads nominations as cbin does
	expect_equal(lame:::.lame_probit_ab_scale(Y, "frn", odmax = 5), below)
	# ordinal has no cap: the argument is ignored
	expect_equal(
		lame:::.lame_probit_ab_scale(Y, "ordinal", odmax = 5),
		lame:::.lame_probit_ab_scale(Y, "ordinal")
	)
})

test_that("a panel with no row at the cap gives the same prior scale either way", {
	Y = .capped_nomination_panel()
	all_rows = lame:::.lame_probit_ab_scale(Y, "cbin")
	# no outdegree reaches a cap of 20
	expect_equal(lame:::.lame_odmax_saturated_share(Y, 20, "cbin"), 0)
	expect_equal(lame:::.lame_probit_ab_scale(Y, "cbin", odmax = 20), all_rows)
	# nor with an absent actor
	Ya = lapply(Y, function(m) { m[1, ] = NA; m[, 1] = NA; m })
	expect_equal(
		lame:::.lame_probit_ab_scale(Ya, "cbin", odmax = 20),
		lame:::.lame_probit_ab_scale(Ya, "cbin")
	)
	# with fewer than five row actors below the cap every row is used (and,
	# the rows at the cap being past the reporting threshold, one pooled
	# scale serves both sides)
	Ys = .capped_nomination_panel(mu = 1)
	expect_gt(lame:::.lame_odmax_saturated_share(Ys, 5, "cbin"), 0.9)
	expect_null(lame:::.lame_prior_rows_below_cap(Ys, 5))
	sc_cap = lame:::.lame_probit_ab_scale(Ys, "cbin", odmax = 5)
	sc_all = lame:::.lame_probit_ab_scale(Ys, "cbin")
	expect_equal(sc_cap$vscale, sc_all$vscale)
	expect_equal(sc_cap$vdfmlt, sc_all$vdfmlt)
	expect_equal(sc_cap$vscale_row, sc_cap$vscale)
})

test_that("rows with hidden cells are read against the cap thinned at their observed share", {
	Y = .capped_nomination_panel(n = 10, TT = 1, mu = 0, seed = 5)[[1]]
	# row 1 is complete with 5 nominations: at a cap of 5
	# row 2 has 5 of its 9 cells hidden and 2 nominations showing: under the
	# cap thinned at its observed share, 5 * 4 / 9
	# row 3 has 5 of its 9 cells hidden and 3 nominations showing: not under it
	Y[1:3, ] = 0; diag(Y) = NA
	Y[1, 2:6] = 1
	Y[2, c(1, 3)] = 1; Y[2, 4:8] = NA
	Y[3, c(1, 2, 4)] = 1; Y[3, 5:9] = NA
	keep = lame:::.lame_prior_rows_below_cap(list(Y), 5)[[1]]
	expect_false(keep[1])
	expect_true(keep[2])
	expect_false(keep[3])
})

test_that("cbin / frn take one prior scale per side", {
	# senders spread more than receivers; half the rows at the cap
	set.seed(31)
	n = 30; cap = 5
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.3)
	Y = lapply(1:4, function(t) {
		z = -1.3 + outer(a, b, "+") + matrix(rnorm(n * n, 0, sqrt(1.49)), n, n)
		diag(z) = -Inf
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			top = top[z[i, top] > 0]
			y[i, top] = 1
		}
		diag(y) = NA
		y
	})
	sc = lame:::.lame_probit_ab_scale(Y, "cbin", odmax = cap)
	expect_true(all(c("vscale_row", "vscale_col") %in% names(sc)))
	# the row side reads the larger sender spread although the rows at the
	# cap carry none of it in their means
	expect_gt(sc$vscale_row, sc$vscale_col)
	expect_gt(sc$vscale_row, 0.4 * var(a))
	# and stays within a factor of four of the column side
	expect_lte(sc$vscale_row, 4 * sc$vscale_col + 1e-12)
	expect_gte(sc$vscale_row, sc$vscale_col / 4 - 1e-12)
	# a single period, or a panel whose rows mostly sit at the cap, keeps the
	# pooled scale for both sides
	s1 = lame:::.lame_probit_ab_scale(Y[[1]], "cbin", odmax = cap)
	expect_equal(s1$vscale_row, s1$vscale)
	expect_equal(s1$vscale_col, s1$vscale)
	Ysat = .capped_nomination_panel(mu = 0.8)
	expect_gte(lame:::.lame_odmax_saturated_share(Ysat, cap, "cbin"), lame:::.lame_odmax_sat_threshold)
	ss = lame:::.lame_probit_ab_scale(Ysat, "cbin", odmax = cap)
	expect_equal(ss$vscale_row, ss$vscale)
	expect_equal(ss$vscale_col, ss$vscale)
	# ordinal keeps one pooled scale for both sides
	so = lame:::.lame_probit_ab_scale(Y, "ordinal")
	expect_equal(so$vscale_row, so$vscale)
	expect_equal(so$vscale_col, so$vscale)
	# lame() puts the two sides on the diagonal of Sab0
	nm = paste0("a", seq_len(n))
	Yn = lapply(Y, function(m) { dimnames(m) = list(nm, nm); m })
	fit = suppressMessages(suppressWarnings(lame(
		Yn,
		family = "cbin",
		odmax = cap,
		R = 0,
		nscan = 20,
		burn = 10,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(unname(diag(fit$prior$Sab0)), c(sc$vscale_row, sc$vscale_col), tolerance = 1e-8)
	expect_equal(fit$prior$Sab0[1, 2], 0)
})

test_that("the censored row scale reads the sender spread with most rows at the cap", {
	set.seed(32)
	n = 30; cap = 5
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Y = lapply(1:4, function(t) {
		z = -0.9 + outer(a, b, "+") + matrix(rnorm(n * n, 0, sqrt(1.49)), n, n)
		diag(z) = -Inf
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			y[i, top[z[i, top] > 0]] = 1
		}
		diag(y) = NA
		y
	})
	share = lame:::.lame_odmax_saturated_share(Y, cap, "cbin")
	expect_gt(share, 0.6)
	s2 = lame:::.lame_censored_row_scale(Y, cap)
	# the probit-scale spread with the receivers' spread absorbed in the link
	expect_gt(s2, 0.3 * var(a))
	expect_lt(s2, 2 * var(a))
	# the rows below the cap alone understate it
	rows = lame:::.lame_prior_rows_below_cap(Y, cap)
	expect_false(is.null(rows))
	# no cap binding, or every row at the cap: nothing to fit
	expect_true(is.na(lame:::.lame_censored_row_scale(Y, 29)))
})

test_that("a single cbin / frn slice reads the prior scale over every row", {
	Y = .capped_nomination_panel()
	y1 = Y[[1]]
	# about half of the slice's rows are at the cap and at least five below it
	expect_false(is.null(lame:::.lame_prior_rows_below_cap(list(y1), 5)))
	for (fam in c("cbin", "frn")) {
		one = lame:::.lame_probit_ab_scale(y1, fam, odmax = 5)
		expect_equal(one, lame:::.lame_probit_ab_scale(y1, fam))
		expect_equal(lame:::.lame_probit_ab_scale(list(y1), fam, odmax = 5), one)
		expect_equal(one$vscale_row, one$vscale)
		expect_equal(one$vscale_col, one$vscale)
	}
	# two or more periods still read the rows below the cap
	expect_false(isTRUE(all.equal(
		lame:::.lame_probit_ab_scale(Y[1:2], "cbin", odmax = 5)$vscale,
		lame:::.lame_probit_ab_scale(Y[1:2], "cbin")$vscale
	)))
})
