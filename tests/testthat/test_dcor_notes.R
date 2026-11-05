skip_on_cran()

test_that("a poisson fit without overdispersion says its dyadic correlation is not identified", {
	# this panel's latent residual has almost no variance; the note reads that variance, which the mutual-dyad count cannot see
	skip_on_cran()
	set.seed(5)
	n = 30; TT = 6
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Y = lapply(1:TT, function(t) {
		eta = 0.3 + outer(a, b, "+")
		m = matrix(rpois(n * n, exp(eta)), n); diag(m) = NA; m
	})
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "poisson",
			R = 1,
			dcor = TRUE,
			dynamic_ab = TRUE,
			burn = 100,
			nscan = 400,
			odens = 5,
			seed = 2,
			verbose = FALSE,
			plot = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_lt(mean(fit$VC[, "ve"]), 0.1)
	expect_true(any(grepl("residual variance of this Poisson fit", msgs)))
	# an overdispersed panel identifies rho and stays silent
	Y2 = lapply(1:TT, function(t) {
		eta = 0.3 + outer(a, b, "+") + matrix(rnorm(n * n, 0, 0.7), n)
		m = matrix(rpois(n * n, exp(eta)), n); diag(m) = NA; m
	})
	msgs2 = character()
	fit2 = withCallingHandlers(
		suppressWarnings(lame(
			Y2,
			family = "poisson",
			R = 1,
			dcor = TRUE,
			dynamic_ab = TRUE,
			burn = 100,
			nscan = 400,
			odens = 5,
			seed = 2,
			verbose = FALSE,
			plot = FALSE
		)),
		message = function(m) { msgs2 <<- c(msgs2, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_gt(mean(fit2$VC[, "ve"]), 0.1)
	expect_false(any(grepl("residual variance of this Poisson fit", msgs2)))
})

# a directed fit with dynamic positions keeps rho: the position kernel whitens each dyad pair by rho, as the static kernel does
test_that("a dcor fit with dynamic positions keeps rho", {
	skip_on_cran()
	set.seed(6886)
	n = 24; TT = 3
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n); diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "binary",
			R = 1,
			dcor = TRUE,
			dynamic_uv = TRUE,
			burn = 20,
			nscan = 60,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_false(any(grepl("made without the dyadic correlation", msgs)))
	expect_true(fit$dcor)
	expect_true(any(fit$VC[, "rho"] != 0))
	expect_equal(fit$tryErrorChecks$UV, 0)
	msgs2 = character()
	fit2 = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "binary",
			R = 1,
			dcor = TRUE,
			burn = 20,
			nscan = 60,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs2 <<- c(msgs2, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_false(any(grepl("made without", msgs2)))
	expect_true(fit2$dcor)
	expect_true(any(fit2$VC[, "rho"] != 0))
	# a period-varying correlation needs static positions
	expect_error(
		suppressMessages(lame(
			Y,
			family = "binary",
			R = 1,
			dcor = TRUE,
			dynamic_uv = TRUE,
			dynamic_rho = TRUE,
			burn = 20,
			nscan = 60,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		"not available with"
	)
	# a sparse panel with a multiplicative term says its intercept is not estimated
	sparse = lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		ev = cbind(sample(n, 4, TRUE), sample(n, 4, TRUE)); ev = ev[ev[, 1] != ev[, 2], , drop = FALSE]
		m[ev] = 1
		m[ev[1, 2:1, drop = FALSE]] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	msgs3 = character()
	fit3 = withCallingHandlers(
		suppressWarnings(lame(
			sparse,
			family = "binary",
			R = 1,
			dcor = TRUE,
			dynamic_uv = TRUE,
			burn = 20,
			nscan = 60,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs3 <<- c(msgs3, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_true(any(grepl("The intercept is not estimated", msgs3)))
})

# a coordinate proposal outside the prior support leaves that conditional unchanged; the other actor updates continue
test_that("refused coordinate draws are counted as prior-bound rejections and stated by the summary on the unipartite path", {
	set.seed(6886)
	n = 16; TT = 3
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.3), n, n); diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	w = messages = character()
	fit = withCallingHandlers(
		lame(
			Y,
			family = "binary",
			R = 1,
			dynamic_uv = TRUE,
			prior = list(uv_max_abs = 1e-4),
			burn = 20,
			nscan = 60,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
		message = function(x) { messages <<- c(messages, conditionMessage(x)); invokeRestart("muffleMessage") }
	)
	expect_gt(fit$tryErrorChecks$UV, 40)
	h = fit$sampler_health
	expect_equal(h$prior_bound_rejections[h$step == "UV"], fit$tryErrorChecks$UV)
	expect_equal(h$numerical_failures[h$step == "UV"], 0)
	# a valid Metropolis rejection is not announced while fitting; the summary states its share of sweeps
	expect_false(any(grepl("uv_max_abs|individual conditionals", messages)))
	expect_length(w, 0L)
	line = lame:::.lame_bound_share_line(fit)
	expect_match(line, "prior\\$uv_max_abs")
	expect_match(line, "kept its state")
	expect_false(grepl("sensitivity|Check|Inspect", line))
})
