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
