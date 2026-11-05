# a sparse directed binary panel reports va, vb and its intercept as not estimated (see ?lame);
# a saturated nomination cap withholds only the intercept of a one-period design

.sparse_panel = function(n, TT, rate, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, rate / (n - 1)), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}
.flat = function(x) gsub("\\s+", " ", paste(x, collapse = " "))
# what a print method writes to the console and as messages
.printed = function(x) {
	msg = character()
	out = capture.output(msg <- capture.output(print(x), type = "message"))
	.flat(c(out, msg))
}

# every actor nominates exactly k targets per period: rate k per period, k * TT events per actor
.fixed_panel = function(n, TT, k = 1, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), k)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

test_that("the reporting rules read the measured rates", {
	pp = lame:::.lame_sparse_period_rate
	expect_equal(pp, 1.2)
	expect_true(lame:::.lame_ab_sparse(pp - 0.01))
	expect_false(lame:::.lame_ab_sparse(pp))
	expect_false(lame:::.lame_ab_sparse(NA_real_))
	expect_false(lame:::.lame_ab_sparse(NULL))
	f = list(family = "binary", mode = "unipartite", events_per_actor_period = 1, R = 0)
	expect_identical(lame:::.lame_ab_var_unident_reason(f), "sparse")
	# the intercept threshold is lower: 0.6, and 1.6 with a multiplicative term
	expect_true(is.na(lame:::.lame_ab_coef_reason(f)))
	expect_identical(lame:::.lame_ab_coef_reason(modifyList(f, list(events_per_actor_period = 0.5))), "sparse")
	expect_true(is.na(lame:::.lame_ab_coef_reason(modifyList(f, list(R = 1, events_per_actor_period = 1.7)))))
	expect_identical(lame:::.lame_ab_coef_reason(modifyList(f, list(R = 1, events_per_actor_period = 1.5))), "sparse")
	# a supplied effects prior reports va and vb, not the intercept
	fp = modifyList(f, list(events_per_actor_period = 0.5, prior_supplied = TRUE))
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fp)))
	expect_identical(lame:::.lame_ab_coef_reason(fp), "sparse")
	# count families and nominations report at any rate
	for (fam in c("cbin", "frn", "poisson")) {
		expect_true(is.na(lame:::.lame_ab_var_unident_reason(list(
			family = fam,
			mode = "unipartite",
			events_per_actor_period = 0.3
		))))
	}
	st = lame:::.lame_ab_statements(modifyList(f, list(events_per_actor_period = 0.5)))
	expect_length(st, 2L)
	expect_match(st[1], "va and vb are not estimated")
	expect_match(st[1], "too few for their intervals")
	expect_false(grepl("conservative|pending|sensitiv", paste(st, collapse = " ")))
})

test_that("a sparse binary fit withholds va, vb and the intercept and keeps the draws", {
	skip_on_cran()
	pp = lame:::.lame_sparse_period_rate
	# twelve periods at half an event per actor-period: six events per actor
	Y = .sparse_panel(24, 12, 0.5)
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "binary",
			R = 0,
			nscan = 100,
			burn = 50,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_lt(fit$events_per_actor_period, pp)
	expect_equal(fit$events_per_actor, sum(sapply(Y, sum, na.rm = TRUE)) / 24)
	expect_match(.flat(msgs), "va and vb are not estimated")
	expect_true(all(c("va", "vb") %in% colnames(fit$VC)))
	s = summary(fit)
	expect_true(all(is.na(s$variance[c("va", "vb"), ])))
	expect_true(all(is.na(s$coefficients["intercept", ])))
	expect_true(all(is.na(confint(fit, parm = c("va", "vb")))))
	expect_match(.printed(fit), "not estimated")
	# a supplied effects prior reports va and vb
	fp = fit; fp$prior_supplied = TRUE
	expect_false(any(is.na(summary(fp)$variance[c("va", "vb"), ])))
	# the total counts every period's events over actors present in any period, so an actor present in
	# fewer periods counts at its own share
	Y_drop = Y[1:3]
	Y_drop[[3]] = Y_drop[[3]][-1, -1]
	fd = suppressMessages(suppressWarnings(lame(
		Y_drop,
		family = "binary",
		R = 0,
		nscan = 20,
		burn = 10,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fd$events_per_actor, sum(sapply(Y_drop, sum, na.rm = TRUE)) / 24)
	expect_equal(
		fd$events_per_actor_period,
		mean(c(
			sum(Y_drop[[1]], na.rm = TRUE) / 24,
			sum(Y_drop[[2]], na.rm = TRUE) / 24,
			sum(Y_drop[[3]], na.rm = TRUE) / 23
		))
	)
})

test_that("a sparse dynamic_ab fit reports rho_ab and sigma_ab and withholds va and vb", {
	skip_on_cran()
	Y = .sparse_panel(24, 4, 0.6)
	f = suppressMessages(suppressWarnings(lame(
		Y,
		family = "binary",
		R = 0,
		dynamic_ab = TRUE,
		nscan = 100,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	s = summary(f)
	expect_true(is.numeric(s$rho_ab))
	expect_true(is.numeric(s$sigma_ab))
	expect_true(all(is.na(s$variance[c("va", "vb"), ])))
	expect_false(any(is.na(s$rho_ab)))
})

test_that("a supplied effects prior is recorded", {
	expect_true(lame:::.lame_ab_prior_given(list(eta0 = 4)))
	expect_true(lame:::.lame_ab_prior_given(list(etaab = 4)))
	expect_true(lame:::.lame_ab_prior_given(list(Sab0 = diag(2))))
	expect_false(lame:::.lame_ab_prior_given(list(rho_uv_mean = 0.5)))
	expect_false(lame:::.lame_ab_prior_given(list()))
	expect_false(lame:::.lame_ab_prior_given(NULL))
})

.sat_panel = function(n = 20, TT = 3, k = 3, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), k)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

.sat_panel = function(n = 20, TT = 3, k = 3, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), k)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

test_that("a one-period saturated cap withholds the intercept and reports va, vb and cab", {
	skip_on_cran()
	Y = .sat_panel(n = 20, TT = 1, k = 3, seed = 7)[[1]]
	msgs = character()
	f = withCallingHandlers(
		suppressWarnings(ame(
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
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_match(.flat(msgs), "The intercept and any row-covariate coefficient are not estimated")
	expect_identical(lame:::.lame_ab_coef_reason(f), "saturated")
	expect_true(is.na(lame:::.lame_ab_data_reason(f)))
	s = summary(f)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab"), c("Estimate", "StdError")])))
	expect_true(all(is.na(s$coefficients["intercept", ])))
	expect_true(all(is.na(confint(f, parm = "intercept"))))
	expect_false(any(is.na(confint(f, parm = "va"))))
	td = tidy(f)
	expect_true(is.na(td$estimate[td$term == "intercept"]))
	expect_match(.flat(capture.output(print(s))), "nothing bounds the level of a row above")
	expect_match(.printed(f), "not estimated")
})

test_that("a multi-period saturated cap reports the intercept and the variances", {
	skip_on_cran()
	set.seed(11)
	n = 20; nm = paste0("a", seq_len(n))
	Y = lapply(1:3, function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) if (runif(1) < 0.9) m[i, sample(setdiff(seq_len(n), i), 1)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		family = "cbin",
		odmax = 1,
		R = 0,
		nscan = 100,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_gte(fit$odmax_saturated_share, lame:::.lame_odmax_sat_threshold)
	expect_true(is.na(lame:::.lame_ab_coef_reason(fit)))
	s = summary(fit)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab"), ])))
	expect_false(any(is.na(s$coefficients["intercept", ])))
})
