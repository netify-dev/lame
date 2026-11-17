# a cbin or frn panel with many rows at the cap and va in [0.5, 1.5) reports the intercept and
# row-covariate coefficients as not estimated; va is reported

.shrunk_fit = function(share = 0.7, va = 1, n_time = 4, family = "cbin", mode = "unipartite") {
	set.seed(1)
	VC = cbind(va = rnorm(50, va, 0.05), cab = rnorm(50, 0, 0.05), vb = rnorm(50, 1, 0.05), rho = 0, ve = 1)
	list(family = family, mode = mode, n_time = n_time, odmax_saturated_share = share, VC = VC)
}

test_that("the rule reads the share, the estimated va and the number of periods", {
	f = lame:::.lame_cap_shrunk_unreported
	expect_equal(lame:::.lame_cap_shrunk_share, 0.5)
	expect_equal(lame:::.lame_cap_shrunk_va, c(0.5, 1.5))
	expect_true(f(.shrunk_fit()))
	expect_true(f(.shrunk_fit(share = 0.5, va = 0.55)))
	expect_true(f(.shrunk_fit(family = "frn")))
	# outside the share, the band, the periods and the designs covered
	expect_false(f(.shrunk_fit(share = 0.45)))
	expect_false(f(.shrunk_fit(va = 0.3)))
	expect_false(f(.shrunk_fit(va = 1.6)))
	expect_false(f(.shrunk_fit(va = 1.5)))
	expect_false(f(.shrunk_fit(n_time = 1)))
	expect_false(f(.shrunk_fit(family = "binary")))
	expect_false(f(.shrunk_fit(mode = "bipartite")))
	expect_false(f(list(family = "cbin", n_time = 4, odmax_saturated_share = 0.7)))
	expect_false(f(NULL))
})

test_that("the reasons, the hidden rows and the statements follow", {
	fit = .shrunk_fit()
	# va is reported: no variance reason, no hidden variance row
	expect_true(is.na(lame:::.lame_ab_data_reason(fit)))
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fit)))
	expect_identical(lame:::.lame_ab_coef_reason(fit), "cap_shrunk")
	expect_length(lame:::.lame_ab_var_hidden_rows("cap_shrunk"), 0L)
	expect_identical(
		lame:::.lame_ab_coef_hidden_rows("cap_shrunk", c("intercept", "x_row", "x_col", "x_dyad")),
		c("intercept", "x_row")
	)
	st = lame:::.lame_ab_statements(fit)
	expect_length(st, 1L)
	expect_match(st[1], "The intercept and any row-covariate coefficient are not estimated")
	expect_match(st[1], "va is estimated")
	expect_match(st[1], "0\\.5 and 1\\.5")
	expect_false(grepl("conservative|pending|caution|sensitiv", st[1]))
	# a supplied effects prior does not identify the level of a row
	fp = modifyList(fit, list(prior_supplied = TRUE, prior = list(Sab0 = diag(2), eta0 = 4)))
	expect_identical(lame:::.lame_ab_coef_reason(fp), "cap_shrunk")
	# the one-period level rule comes first
	one = .shrunk_fit(share = 1, n_time = 1)
	expect_identical(lame:::.lame_ab_coef_reason(one), "saturated")
})

test_that("a fit in the band reports the variances and withholds the intercept", {
	skip_on_cran()
	n = 16; nm = paste0("a", seq_len(n))
	set.seed(3)
	Y = lapply(1:3, function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), 3)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fit$odmax_saturated_share, 1)
	# put the posterior mean of va in the band without refitting
	fit$VC[, "va"] = 1
	s = summary(fit)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab"), ])))
	expect_true(all(is.na(s$coefficients["intercept", ])))
	expect_true(all(is.na(confint(fit, parm = "intercept"))))
	expect_false(any(is.na(confint(fit, parm = c("va", "vb")))))
	expect_true(is.na(s$ab_data_reason))
	expect_identical(s$ab_coef_reason, "cap_shrunk")
	expect_match(
		paste(capture.output(print(s)), collapse = " "),
		"The intercept and any row-covariate coefficient are not estimated"
	)
	td = tidy(fit)
	expect_true(is.na(td$estimate[td$term == "intercept"]))
	# outside the band everything is reported
	fit$VC[, "va"] = 2.5
	s2 = summary(fit)
	expect_false(any(is.na(s2$variance[c("va", "vb", "cab"), ])))
	expect_false(any(is.na(s2$coefficients["intercept", ])))
})

# saturated caps: the share of rows at the cap, the reporting levels and their named constants

# a cbin or frn nomination cap sitting at (almost) every row's outdegree erases
# the cross-row spread the additive effects rest on
.sat_cbin_panel = function(n = 20, TT = 3, k = 3, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) {
			m[i, sample(setdiff(seq_len(n), i), k)] = 1
		}
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

test_that("a saturating odmax in a multi-period fit reports every parameter", {
	skip_on_cran()
	Y = .sat_cbin_panel()
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "cbin",
			odmax = 3,
			R = 0,
			nscan = 200,
			burn = 100,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	# every row nominates exactly the cap
	expect_equal(fit$odmax_saturated_share, 1)
	expect_gt(fit$events_per_actor_period, 1)
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fit)))
	expect_true(is.na(lame:::.lame_ab_coef_reason(fit)))
	expect_false(grepl("not estimated", gsub("\\s+", " ", paste(msgs, collapse = " "))))
	s = summary(fit)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab", "ve"), ])))
	expect_false(any(is.na(s$coefficients["intercept", ])))
	expect_equal(s$odmax_saturated_share, 1)
	expect_false(any(is.na(confint(fit, parm = c("intercept", "va", "vb", "ve")))))
	td = tidy(fit)
	expect_false(is.na(td$estimate[td$term == "intercept"]))
	expect_false(grepl("not estimated", gsub("\\s+", " ", paste(capture.output(print(s)), collapse = " "))))
})

test_that("hidden cells do not hide a saturated cap", {
	skip_on_cran()
	# every row nominates the full cap; a fifth of the cells is then hidden,
	# nominations included, so the observed outdegrees fall below the cap
	Y = .sat_cbin_panel()
	set.seed(9)
	Ym = lapply(Y, function(m) {
		off = which(!is.na(m))
		m[sample(off, round(0.2 * length(off)))] = NA
		m
	})
	full = lame:::.lame_odmax_saturated_share(Y, 3, "cbin")
	hid = lame:::.lame_odmax_saturated_share(Ym, 3, "cbin")
	expect_equal(full, 1)
	# the share is estimated from the thinned outdegrees: it reaches the
	# level at which a one-period design withholds its intercept
	expect_gte(as.numeric(hid), lame:::.lame_cap_level_share_hidden)
	one = list(family = "cbin", mode = "unipartite", odmax_saturated_share = hid, cells_hidden = TRUE, n_time = 1)
	expect_identical(lame:::.lame_ab_coef_reason(one), "saturated")
	expect_true(is.na(lame:::.lame_ab_coef_reason(modifyList(one, list(cells_hidden = FALSE)))))
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(one)))
	expect_match(lame:::.lame_sat_phrase(hid), "nominate the full odmax")
	# a cap far above the outdegrees stays unsaturated with the same holes
	far = lame:::.lame_odmax_saturated_share(Ym, 15, "cbin")
	expect_lt(far, 0.2)
	expect_false(isTRUE(attr(far, "not_ruled_out")))
	# outdegrees spread below a small cap are not read as saturated
	set.seed(12)
	Ys = lapply(1:4, function(t) {
		m = matrix(0, 30, 30)
		for (i in 1:30) m[i, sample(setdiff(1:30, i), sample(1:5, 1))] = 1
		diag(m) = NA
		off = which(!is.na(m)); m[sample(off, round(0.2 * length(off)))] = NA
		m
	})
	spread = lame:::.lame_odmax_saturated_share(Ys, 5, "cbin")
	expect_lt(spread, 0.5)
	expect_false(isTRUE(attr(spread, "not_ruled_out")))
	# an absent actor (a whole row and column missing) is not a hidden cell
	Ya = lapply(Y, function(m) { m[1, ] = NA; m[, 1] = NA; m })
	expect_equal(
		lame:::.lame_odmax_saturated_share(Ya, 3, "cbin"),
		lame:::.lame_odmax_saturated_share(lapply(Y, function(m) m[-1, -1]), 3, "cbin")
	)
})

test_that("a cap above the typical outdegree reports everything", {
	skip_on_cran()
	set.seed(4)
	n = 20; TT = 3
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) {
			k = sample(0:5, 1)
			if (k > 0) m[i, sample(setdiff(seq_len(n), i), k)] = 1
		}
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		family = "cbin",
		odmax = 9,
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_lt(fit$odmax_saturated_share, lame:::.lame_odmax_sat_threshold)
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fit)))
	expect_false(any(is.na(summary(fit)$variance[c("va", "vb"), ])))
	expect_false(any(is.na(confint(fit, parm = "va"))))
	# a one-period fit at a saturated share withholds its intercept
	one = fit; one$n_time = 1; one$n.periods = 1; one$odmax_saturated_share = 1
	expect_identical(lame:::.lame_ab_coef_reason(one), "saturated")
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(one)))
	# a fit without the recorded share reports its intercept
	one$odmax_saturated_share = NULL
	expect_true(is.na(lame:::.lame_ab_coef_reason(one)))
})

test_that("the saturation share is computed only where a cap applies", {
	skip_on_cran()
	Y = .sat_cbin_panel(n = 12, TT = 2, k = 2, seed = 3)
	f = lame:::.lame_odmax_saturated_share
	expect_equal(f(Y, 2, "cbin"), 1)
	expect_equal(f(Y, 2, "frn"), 1)
	# a cap above every row's outdegree saturates nothing
	expect_equal(f(Y, 5, "cbin"), 0)
	# families without a cap carry no share at all
	expect_true(is.na(f(Y, 2, "binary")))
	expect_true(is.na(f(Y, 2, "normal")))
	# a missing cap resolves to the largest observed outdegree, as lame() does
	expect_equal(f(Y, NULL, "cbin"), 1)
	# a per-row cap is read row by row
	expect_equal(f(Y, c(rep(2, 6), rep(5, 6)), "cbin"), 0.5)
	# a 3-d array is the same panel as a list of slices
	arr = array(unlist(Y), c(12, 12, 2))
	expect_equal(f(arr, 2, "cbin"), 1)
})

test_that("the named constants set the saturated-cap setting and the reporting levels", {
	skip_on_cran()
	thr = lame:::.lame_odmax_sat_threshold
	expect_true(is.numeric(thr) && length(thr) == 1L && thr > 0.5 && thr < 1)
	expect_equal(lame:::.lame_cap_level_periods, 1L)
	expect_true(lame:::.lame_cap_level_share_hidden < lame:::.lame_cap_level_share)
	reason = function(s, ...) lame:::.lame_ab_coef_reason(list(
		family = "cbin",
		mode = "unipartite",
		odmax_saturated_share = s,
		n_time = 1,
		...
	))
	expect_identical(reason(1), "saturated")
	expect_true(is.na(reason(lame:::.lame_cap_level_share - 0.01)))
	expect_identical(reason(lame:::.lame_cap_level_share_hidden, cells_hidden = TRUE), "saturated")
	expect_true(is.na(reason(lame:::.lame_cap_level_share_hidden - 0.01, cells_hidden = TRUE)))
	expect_match(lame:::.lame_sat_phrase(0.97), "97 percent of the observed rows nominate the full odmax")
	# the share is computed on the observed outdegrees
	n = 20
	k_at = ceiling(thr * n - 1e-9)
	Y = lapply(1:2, function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, setdiff(seq_len(n), i)[seq_len(if (i <= k_at) 3 else 2)]] = 1
		diag(m) = NA; m
	})
	sh = lame:::.lame_odmax_saturated_share(Y, 3, "cbin")
	expect_equal(as.numeric(sh), k_at / n)
	# the pre-fit message of a one-period fit fires at a full cap and not under it
	msg_of = function(Yp) {
		msgs = character()
		withCallingHandlers(
			suppressWarnings(ame(
				Yp,
				family = "cbin",
				odmax = 3,
				R = 0,
				nscan = 20,
				burn = 10,
				odens = 1,
				verbose = FALSE,
				gof = FALSE
			)),
			message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
		)
		gsub("\\s+", " ", paste(msgs, collapse = " "))
	}
	nm = paste0("a", seq_len(n))
	full = .sat_cbin_panel(n = n, TT = 1, k = 3, seed = 5)[[1]]
	under = full; under[1:5, ] = 0; under[cbind(1:5, 6:10)] = 1; under[cbind(1:5, 11:15)] = 1; diag(under) = NA
	expect_match(msg_of(full), "nominate the full odmax")
	expect_false(grepl("nominate the full odmax", msg_of(under)))
	rd = tools::Rd_db("lame")
	for (f in c("lame.Rd", "ame.Rd")) {
		txt = gsub(
			"\\s+",
			" ",
			paste(
				capture.output(tools::Rd2txt(
					rd[[f]],
					options = list(underline_titles = FALSE)
				)),
				collapse = " "
			)
		)
		expect_match(txt, paste0(round(100 * thr), " percent of the observed rows"))
	}
})

test_that("a single period is judged against its own saturated-cap setting", {
	skip_on_cran()
	thr1 = lame:::.lame_odmax_sat_threshold_single
	thr = lame:::.lame_odmax_sat_threshold
	expect_true(is.numeric(thr1) && length(thr1) == 1L && thr1 > 0.3 && thr1 < thr)
	expect_equal(lame:::.lame_odmax_threshold_for(1L), thr1)
	expect_equal(lame:::.lame_odmax_threshold_for(4L), thr)
	reason = function(s, cls = NULL, ...) lame:::.lame_ab_coef_reason(structure(
		list(
			family = "cbin",
			mode = "unipartite",
			odmax_saturated_share = s,
			...
		),
		class = cls
	))
	# an ame() fit has one period, a lame() fit the periods it carries, and an unclassed list the panel rule
	expect_identical(reason(1, "ame"), "saturated")
	expect_true(is.na(reason(1, c("lame", "ame"), n_time = 4L)))
	expect_identical(reason(1, c("lame", "ame"), n_time = 1L), "saturated")
	expect_true(is.na(reason(1)))
	# a slice with about seven in ten rows at the cap reports everything
	n = 20; k_at = 14
	nm = paste0("a", seq_len(n))
	mk = function(seed) {
		set.seed(seed)
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), if (i <= k_at) 3 else 1)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	}
	Y1 = mk(1)
	expect_equal(as.numeric(lame:::.lame_odmax_saturated_share(Y1, 3, "cbin")), k_at / n)
	fa = suppressMessages(suppressWarnings(ame(
		Y1,
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)))
	expect_equal(fa$odmax_sat_threshold, thr1)
	expect_true(is.na(lame:::.lame_ab_coef_reason(fa)))
	expect_false(any(is.na(summary(fa)$variance["va", c("Estimate", "StdError")])))
	fl = suppressMessages(suppressWarnings(lame(
		list(Y1, mk(2)),
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fl$odmax_sat_threshold, thr)
	# the one-period rule does not apply to two periods
	expect_false(identical(lame:::.lame_ab_coef_reason(fl), "saturated"))
	f1 = suppressMessages(suppressWarnings(lame(
		list(Y1),
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(f1$odmax_sat_threshold, thr1)
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(f1)))
	rd = tools::Rd_db("lame")
	for (f in c("lame.Rd", "ame.Rd")) {
		txt = gsub(
			"\\s+",
			" ",
			paste(
				capture.output(tools::Rd2txt(
					rd[[f]],
					options = list(underline_titles = FALSE)
				)),
				collapse = " "
			)
		)
		expect_match(txt, paste0(round(100 * thr1), " percent in a single period"))
	}
})
