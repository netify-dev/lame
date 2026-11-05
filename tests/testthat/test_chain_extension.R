skip_on_cran()

chain_panel = function(n = 12L, nt = 3L, seed = 31) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(nt), function(t) { m = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(m) = NA; m })
	X = lapply(seq_len(nt), function(t) array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x")))
	list(Y = Y, X = X)
}

# the check as the fit sees it: flagged until `settle_after` calls have been made
flag_until = function(settle_after) {
	calls = 0L
	function(fit, ...) {
		calls <<- calls + 1L
		if (calls <= settle_after) {
			data.frame(
				parameter = "intercept", mean_first = 0, mean_last = 1, z = 9, shift = 3, text = "intercept (0 to 1)",
				stringsAsFactors = FALSE
			)
		} else data.frame(parameter = character(0), text = character(0))
	}
}

lame_args = function(d, ...) list(
	d$Y,
	Xdyad = d$X,
	family = "normal",
	R = 0,
	dynamic_ab = TRUE,
	burn = 20,
	nscan = 100,
	odens = 5,
	seed = 11,
	verbose = FALSE,
	...
)

test_that("the check names a drifting series and not a stationary one", {
	set.seed(5)
	drift = list(
		BETA = cbind(intercept = seq(0, 1, length.out = 300) + rnorm(300, 0, 0.05)),
		VC = cbind(va = rnorm(300, 1, 0.1))
	)
	mv = lame:::.lame_chain_moving(drift)
	expect_equal(mv$parameter, "intercept")
	expect_match(mv$text, "^intercept \\(")
	still = list(BETA = cbind(intercept = rnorm(300)), VC = cbind(va = rnorm(300, 1, 0.1), rho = rnorm(300)))
	expect_equal(nrow(lame:::.lame_chain_moving(still)), 0L)
	# fewer than sixty stored draws are not checked
	short = list(BETA = cbind(intercept = seq(0, 1, length.out = 40)))
	expect_equal(nrow(lame:::.lame_chain_moving(short)), 0L)
	# rho of ame() is monitored
	rho_drift = list(
		BETA = cbind(intercept = rnorm(300)),
		VC = cbind(rho = seq(0, 1, length.out = 300) + rnorm(300, 0, 0.05))
	)
	expect_equal(lame:::.lame_chain_moving(rho_drift)$parameter, "rho")
})

test_that("a settled lame fit carries a record with no extension", {
	d = chain_panel()
	fit = suppressMessages(suppressWarnings(do.call(lame, lame_args(d))))
	expect_identical(fit$auto_extend$extensions, 0L)
	expect_identical(fit$auto_extend$discarded_sweeps, 0L)
	expect_identical(fit$auto_extend$max_extensions, 3L)
	expect_length(fit$auto_extend$triggers, 0L)
	expect_equal(nrow(fit$BETA), 20L)
})

test_that("a flagged chain is continued as a block with the same number of draws, and the record says so", {
	skip_on_cran()
	d = chain_panel()
	args = lame_args(d)
	ref = suppressMessages(suppressWarnings(do.call(lame, args)))
	local_mocked_bindings(.lame_chain_moving = flag_until(1L), .package = "lame")
	fit = suppressMessages(suppressWarnings(do.call(lame, args)))
	expect_equal(nrow(fit$BETA), nrow(ref$BETA))
	expect_identical(fit$auto_extend$extensions, 1L)
	expect_identical(fit$auto_extend$discarded_sweeps, 100L)
	expect_equal(fit$auto_extend$triggers[[1]], "intercept (0 to 1)")
	# the call and the refit options describe the user's call, not the block
	expect_identical(fit$call, ref$call)
	expect_identical(fit$refit_options, ref$refit_options)
	# the block is the chain continued from the final state of the first run with the advanced seed
	cont_args = lame_args(d, start_vals = ref$start_vals)
	cont_args$burn = 0L; cont_args$seed = 11 + 20 + 100
	cont = suppressMessages(suppressWarnings(do.call(lame, cont_args)))
	expect_identical(fit$BETA, cont$BETA)
	expect_identical(fit$VC, cont$VC)
	expect_identical(fit$start_vals, cont$start_vals)
})

test_that("a chain still moving after the extensions is refused with the last block", {
	skip_on_cran()
	d = chain_panel()
	local_mocked_bindings(.lame_chain_moving = flag_until(Inf), .package = "lame")
	e = tryCatch(suppressMessages(suppressWarnings(do.call(lame, lame_args(d)))), error = function(e) e)
	expect_s3_class(e, "lame_chain_unsettled")
	expect_s3_class(e$fit, "lame")
	expect_equal(nrow(e$fit$BETA), 20L)
	expect_identical(e$fit$auto_extend$extensions, 3L)
	expect_true(e$fit$auto_extend$refused)
	expect_equal(e$moving, "intercept (0 to 1)")
	expect_match(conditionMessage(e), "intercept")
	expect_match(conditionMessage(e), "lame_resume(e$fit, nscan = 100)", fixed = TRUE)
	# the refused fit is continued with lame_resume(), which checks the new block again
	local_mocked_bindings(.lame_chain_moving = flag_until(0L), .package = "lame")
	r = suppressMessages(suppressWarnings(lame_resume(e$fit, nscan = 60)))
	expect_s3_class(r, "lame")
	expect_equal(nrow(r$BETA), 12L)
	expect_identical(r$auto_extend$extensions, 0L)
})

test_that("prior$chain_extensions sets the number of extensions, 0 refuses at the first flag", {
	skip_on_cran()
	d = chain_panel()
	local_mocked_bindings(.lame_chain_moving = flag_until(Inf), .package = "lame")
	e0 = tryCatch(
		suppressMessages(suppressWarnings(do.call(lame, lame_args(d, prior = list(chain_extensions = 0))))),
		error = function(e) e
	)
	expect_s3_class(e0, "lame_chain_unsettled")
	expect_identical(e0$extensions, 0L)
	expect_match(conditionMessage(e0), "chain_extensions")
	e1 = tryCatch(
		suppressMessages(suppressWarnings(do.call(lame, lame_args(d, prior = list(chain_extensions = 1))))),
		error = function(e) e
	)
	expect_identical(e1$extensions, 1L)
	expect_error(do.call(lame, lame_args(d, prior = list(chain_extensions = -1))), "non-negative whole number")
	expect_error(do.call(lame, lame_args(d, prior = list(chain_extensions = 1.5))), "non-negative whole number")
})

test_that("an extended fit leaves the global random stream where it found it", {
	d = chain_panel()
	local_mocked_bindings(.lame_chain_moving = flag_until(2L), .package = "lame")
	set.seed(99); r0 = runif(1); set.seed(99)
	suppressMessages(suppressWarnings(do.call(lame, lame_args(d))))
	expect_identical(runif(1), r0)
})

test_that("ame() extends, refuses and resumes the same way", {
	skip_on_cran()
	set.seed(3)
	n = 14; nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	Yb = matrix(rnorm(12 * 9), 12, 9, dimnames = list(paste0("r", 1:12), paste0("c", 1:9)))
	for (case in list(list(Y = Y, R = 1), list(Y = Yb, R = 1, mode = "bipartite"))) {
		run = function(...) suppressMessages(
			suppressWarnings(
				do.call(ame, c(case, list(burn = 20, nscan = 100, odens = 5, seed = 3, verbose = FALSE), list(...)))
			)
		)
		ref = run()
		expect_identical(ref$auto_extend$extensions, 0L)
		expect_false(is.null(ref$start_vals))
		local_mocked_bindings(.lame_chain_moving = flag_until(1L), .package = "lame")
		fit = run()
		expect_identical(fit$auto_extend$extensions, 1L)
		expect_equal(nrow(fit$BETA), nrow(ref$BETA))
		expect_identical(fit$call, ref$call)
		cont = suppressMessages(suppressWarnings(do.call(
			ame,
			c(
				case,
				list(
					burn = 0L,
					nscan = 100,
					odens = 5,
					seed = 3 + 20 + 100,
					verbose = FALSE,
					start_vals = ref$start_vals
				)
			)
		)))
		expect_identical(fit$BETA, cont$BETA)
		local_mocked_bindings(.lame_chain_moving = flag_until(Inf), .package = "lame")
		e = tryCatch(run(), error = function(e) e)
		expect_s3_class(e, "lame_chain_unsettled")
		expect_equal(nrow(e$fit$BETA), 20L)
		local_mocked_bindings(.lame_chain_moving = flag_until(0L), .package = "lame")
		r = suppressMessages(suppressWarnings(lame_resume(e$fit, nscan = 50)))
		expect_s3_class(r, "ame")
		expect_equal(nrow(r$BETA), 10L)
	}
})

test_that("the chains of a multi-chain fit extend together", {
	skip_on_cran()
	set.seed(4)
	n = 12; nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	local_mocked_bindings(.lame_chain_moving = flag_until(1L), .package = "lame")
	res = suppressMessages(
		suppressWarnings(
			ame_parallel(
				Y,
				n_chains = 2,
				cores = 1,
				combine_method = "list",
				burn = 20,
				nscan = 100,
				odens = 5,
				verbose = FALSE
			)
		)
	)
	expect_identical(vapply(res, function(f) f$auto_extend$extensions, integer(1)), c(1L, 1L))
	expect_identical(vapply(res, function(f) nrow(f$BETA), integer(1)), c(20L, 20L))
})

test_that("print() reports an extension", {
	rec = lame:::.lame_chain_record(3L, 2L, 100L, list("sigma_ab (0.4 to 0.3)", "sigma_ab (0.35 to 0.3)"))
	expect_match(lame:::.lame_chain_print_line(rec), "extended 2 times .* sigma_ab")
	expect_null(lame:::.lame_chain_print_line(lame:::.lame_chain_record(3L)))
})

test_that("a multi-chain call that stays unsettled refuses with every chain attached", {
	skip_on_cran()
	set.seed(4)
	n = 12; nm = paste0("a", seq_len(n))
	Y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(Y) = NA
	local_mocked_bindings(.lame_chain_moving = flag_until(Inf), .package = "lame")
	e = tryCatch(
		suppressMessages(
			suppressWarnings(
				ame_parallel(Y, n_chains = 2, cores = 1, burn = 20, nscan = 100, odens = 5, verbose = FALSE)
			)
		),
		error = function(e) e
	)
	expect_s3_class(e, "lame_chain_unsettled")
	expect_identical(e$extensions, 3L)
	expect_length(e$fits, 2L)
	expect_s3_class(e$fits[[1]], "ame")
	expect_s3_class(e$fit, "ame")
	expect_identical(e$fit$chains, 2L)
})
