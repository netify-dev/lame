# the response scale the dynamic-coefficient innovation prior is set from, and
# the reporting rule for the position hyperparameters of the two
# heavy-transition kinds
skip_on_cran()

test_that("the response-scale unit is one for the families whose latent scale is fixed", {
	set.seed(6886)
	n = 12; TT = 3
	Yb = array(rbinom(n * n * TT, 1, 0.3), c(n, n, TT))
	for (fam in c("binary", "ordinal", "cbin", "frn")) {
		expect_equal(lame:::.lame_response_resid_var(Yb, fam), 1)
	}
	# the normal family reads the variance of the two-way centred outcome
	Yn = array(rnorm(n * n * TT), c(n, n, TT))
	tw = function(M) {
		r = rowMeans(M, na.rm = TRUE); cc = colMeans(M, na.rm = TRUE)
		M - outer(r, cc, "+") + mean(M, na.rm = TRUE)
	}
	want = var(unlist(lapply(seq_len(TT), function(t) tw(Yn[, , t]))))
	expect_equal(lame:::.lame_response_resid_var(Yn, "normal"), want)
	# a scaled outcome scales the unit by the square of the factor
	expect_equal(lame:::.lame_response_resid_var(Yn * 10, "normal"), want * 100)
	# Poisson reads it on the log scale
	Yp = array(rpois(n * n * TT, 3), c(n, n, TT))
	wantp = var(unlist(lapply(seq_len(TT), function(t) tw(log1p(Yp[, , t])))))
	expect_equal(lame:::.lame_response_resid_var(Yp, "poisson"), wantp)
	# a degenerate outcome falls back to one rather than to zero
	expect_equal(lame:::.lame_response_resid_var(array(1, c(n, n, TT)), "normal"), 1)
})

test_that("the sigma_beta prior scale is a fifth of a residual sd, squared", {
	set.seed(6886)
	n = 12; TT = 3
	Yb = array(rbinom(n * n * TT, 1, 0.3), c(n, n, TT))
	# shape 2: the prior mean of sigma_beta^2 is the scale itself
	expect_equal(lame:::.lame_sigma_beta_prior_scale(Yb, "binary"), 0.04)
	# the prior mean is held as the shape moves
	expect_equal(lame:::.lame_sigma_beta_prior_scale(Yb, "binary", shape = 3), 0.08)
	expect_equal(
		lame:::.lame_sigma_beta_prior_scale(Yb, "binary", shape = 3) / 2,
		lame:::.lame_sigma_beta_prior_scale(Yb, "binary", shape = 2)
	)
	# and it follows the response scale for a normal outcome
	Yn = array(rnorm(n * n * TT), c(n, n, TT))
	expect_equal(
		lame:::.lame_sigma_beta_prior_scale(Yn * 10, "normal"),
		100 * lame:::.lame_sigma_beta_prior_scale(Yn, "normal")
	)
	expect_gt(lame:::.lame_sigma_beta_prior_scale(Yn, "normal"), 0)
})

test_that("a dynamic_beta fit reports the scale it used and a user value overrides it", {
	skip_on_cran()
	set.seed(6886)
	n = 14; TT = 4
	nm = paste0("a", seq_len(n))
	X = lapply(seq_len(TT), function(t) {
		x = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm))
		array(x, c(n, n, 1), dimnames = list(nm, nm, "x1"))
	})
	Y = lapply(seq_len(TT), function(t) {
		m = 0.5 * X[[t]][, , 1] + matrix(rnorm(n * n), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		X,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	Yarr = array(unlist(Y), c(n, n, TT))
	expect_equal(
		fit$sigma_beta_prior_scale,
		lame:::.lame_sigma_beta_prior_scale(Yarr, "normal", 2)
	)
	expect_equal(fit$prior$sigma_beta_scale, fit$sigma_beta_prior_scale)
	expect_true(all(is.finite(fit$SIGMA_BETA[, "dyad"])))
	# a user value is kept as given
	fit2 = suppressMessages(suppressWarnings(lame(
		Y,
		X,
		family = "normal",
		R = 0,
		dynamic_beta = "dyad",
		prior = list(sigma_beta_scale = 0.5),
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fit2$prior$sigma_beta_scale, 0.5)
	expect_equal(fit2$sigma_beta_prior_scale, 0.5)
	# a fit without a dynamic coefficient block neither reports nor changes it
	fit3 = suppressMessages(suppressWarnings(lame(
		Y,
		X,
		family = "normal",
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_null(fit3$sigma_beta_prior_scale)
	expect_equal(fit3$prior$sigma_beta_scale, 1)
})

test_that("the position hyperparameters are reported for ar1, and for snap or t at rank one", {
	f = list(dynamic_uv = TRUE, dynamic_uv_kind = "ar1")
	expect_false(lame:::.lame_uv_hyper_unidentified(f))
	f$dynamic_uv_kind = "snap"
	expect_true(lame:::.lame_uv_hyper_unidentified(f))
	expect_true(lame:::.lame_snap_rate_unidentified(f))
	f$dynamic_uv_kind = "t"
	expect_true(lame:::.lame_uv_hyper_unidentified(f))
	expect_false(lame:::.lame_snap_rate_unidentified(f))
	# the rank decides: rank one is reported, rank two or more is withheld
	for (k in c("snap", "t")) {
		expect_false(lame:::.lame_uv_hyper_unidentified(list(
			dynamic_uv = TRUE,
			dynamic_uv_kind = k,
			R = 1
		)))
		expect_true(lame:::.lame_uv_hyper_unidentified(list(
			dynamic_uv = TRUE,
			dynamic_uv_kind = k,
			R = 2
		)))
		expect_true(lame:::.lame_uv_hyper_unidentified(list(
			dynamic_uv = TRUE,
			dynamic_uv_kind = k,
			R = 3
		)))
	}
	expect_false(lame:::.lame_snap_rate_unidentified(list(
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap", R = 1
	)))
	expect_true(lame:::.lame_snap_rate_unidentified(list(
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap", R = 2
	)))
	# a fit without R is read from the position array
	expect_false(lame:::.lame_uv_hyper_unidentified(list(
		dynamic_uv = TRUE,
		dynamic_uv_kind = "t",
		U = array(0, c(5, 1, 3))
	)))
	expect_true(lame:::.lame_uv_hyper_unidentified(list(
		dynamic_uv = TRUE,
		dynamic_uv_kind = "t",
		U = array(0, c(5, 2, 3))
	)))
	# a static-position fit is never flagged, whatever the kind slot says
	expect_false(lame:::.lame_uv_hyper_unidentified(list(
		dynamic_uv = FALSE,
		dynamic_uv_kind = "snap"
	)))
	# a fit that records the kind only in the dynamic list is read too
	expect_true(lame:::.lame_uv_hyper_unidentified(list(
		dynamic_uv = TRUE,
		dynamic = list(uv_kind = "t")
	)))
	# and one with no kind at all keeps the default reporting
	expect_false(lame:::.lame_uv_hyper_unidentified(list(dynamic_uv = TRUE)))
})

test_that("a rank-two snap fit says why rho_uv, sigma_uv and pi_snap are not reported", {
	set.seed(6886)
	n = 16; TT = 4
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "binary",
			R = 2,
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			nscan = 200,
			burn = 100,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_true(lame:::.lame_uv_hyper_unidentified(fit))
	expect_true(lame:::.lame_snap_rate_unidentified(fit))
	msg = gsub("\\s+", " ", paste(msgs, collapse = " "))
	expect_match(msg, "at rank two or more")
	expect_match(msg, "extra resets stand in for drift")
	expect_false(grepl("not yet validated", msg))
	s = summary(fit)
	expect_null(s$rho_uv)
	expect_null(s$sigma_uv)
	expect_null(s$pi_snap)
	expect_true(isTRUE(s$uv_hyper_unidentified))
	expect_true(isTRUE(s$pi_snap_unidentified))
	out_s = gsub("\\s+", " ", paste(capture.output(print(s)), collapse = " "))
	expect_match(out_s, "rho_uv and sigma_uv are not reported")
	expect_match(out_s, "pi_snap is not reported")
	out = gsub("\\s+", " ", paste(capture.output(print(fit)), collapse = " "))
	expect_match(out, "not reported")
	# the draws are kept
	expect_true(all(is.finite(fit$rho_uv)))
	expect_true(all(is.finite(fit$sigma_uv)))
	expect_true(all(is.finite(fit$pi_snap)))
})

test_that("rank-one snap and t fits report rho_uv, sigma_uv and pi_snap", {
	skip_on_cran()
	set.seed(6886)
	n = 16; TT = 4
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	for (k in c("snap", "t")) {
		msgs = character()
		fit = withCallingHandlers(
			suppressWarnings(lame(
				Y,
				family = "binary",
				R = 1,
				dynamic_uv = TRUE,
				dynamic_uv_kind = k,
				nscan = 200,
				burn = 100,
				odens = 5,
				verbose = FALSE,
				plot = FALSE,
				gof = FALSE
			)),
			message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
		)
		expect_false(lame:::.lame_uv_hyper_unidentified(fit))
		expect_false(grepl("not report", paste(msgs, collapse = " ")))
		s = summary(fit)
		expect_equal(s$rho_uv, mean(fit$rho_uv))
		expect_equal(s$sigma_uv, mean(fit$sigma_uv))
		expect_false(isTRUE(s$uv_hyper_unidentified))
		out_s = gsub("\\s+", " ", paste(capture.output(print(s)), collapse = " "))
		expect_false(grepl("not reported", out_s))
		out = gsub("\\s+", " ", paste(capture.output(print(fit)), collapse = " "))
		expect_false(grepl("not reported", out))
		if (k == "snap") {
			expect_false(isTRUE(s$pi_snap_unidentified))
			expect_equal(s$pi_snap, mean(fit$pi_snap))
			expect_match(out_s, "pi_snap =")
		} else {
			expect_null(s$pi_snap)
		}
	}
})

test_that("an ar1 dynamic_uv fit keeps reporting rho_uv and sigma_uv", {
	skip_on_cran()
	set.seed(6886)
	n = 16; TT = 4
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		family = "binary",
		R = 1,
		dynamic_uv = TRUE,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_false(lame:::.lame_uv_hyper_unidentified(fit))
	s = summary(fit)
	expect_true(is.numeric(s$rho_uv))
	expect_true(is.numeric(s$sigma_uv))
	expect_false(isTRUE(s$uv_hyper_unidentified))
	out = gsub("\\s+", " ", paste(capture.output(print(fit)), collapse = " "))
	expect_false(grepl("not reported", out))
})
