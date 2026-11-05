# a factor column the data do not support is retired, not refused: the fit is that of the next rank down, and says so
# (R/ame_als_retire.R); a stalled five-node pair rule of a weakly informed fit is refined to seven nodes instead
# (tested here and with the null network in test_als_multistart_selection.R)

# a directed binary panel of 50 actors with a rank-R product: sender and receiver effects (sd 0.6, correlation 0.3),
# one dyadic covariate (slope 0.7) and a product of variance 0.6 per factor
ret_sim_binary = function(n, R, mu, seed, sd_uv = sqrt(0.6)) {
	set.seed(seed)
	nm = sprintf("a%03d", seq_len(n))
	ab = matrix(rnorm(2 * n), n, 2) %*% chol(0.36 * matrix(c(1, 0.3, 0.3, 1), 2))
	U = matrix(rnorm(n * R, 0, sd_uv), n); V = matrix(rnorm(n * R, 0, sd_uv), n)
	x = matrix(rnorm(n * n), n)
	y = 1 * (mu + 0.7 * x + outer(ab[, 1], ab[, 2], "+") + tcrossprod(U, V) + matrix(rnorm(n * n), n) > 0)
	diag(y) = NA
	dimnames(y) = list(nm, nm)
	list(Y = y, X = array(x, c(n, n, 1), dimnames = list(nm, nm, "x_dyad")))
}

test_that("the retirement loop retires the weakest column of an iteration that does not settle, down to the rank that does", {
	mk = function(Rr, settled) {
		U = matrix(1, 6, Rr) %*% diag(seq(Rr, 1), Rr)
		list(
			fit = list(
				U = U,
				V = U,
				sig2_u = rep(1, Rr),
				sig2_v = rep(1, Rr),
				eb_mf = list(
					record = list(converged = settled, termination = if (settled) "converged" else "no improving step")
				)
			),
			stationarity = list(
				pass = settled,
				termination = if (settled) "converged" else "no improving conditional step"
			),
			retained = Rr
		)
	}
	seen = list()
	run = function(init, Rr) { seen[[length(seen) + 1L]] <<- list(init = init, rank = Rr); mk(Rr, settled = Rr <= 2L) }
	init = list(U = matrix(seq_len(24), 6, 4), V = matrix(seq_len(24), 6, 4) + 100)
	ir = lame:::.ae_retire_columns(run, init, 4L, FALSE, TRUE)
	expect_identical(vapply(seen, function(z) z$rank, 0L), c(4L, 3L, 2L))
	expect_identical(ir$retained, 2L)
	ret = ir$retirement
	expect_identical(ret$requested_rank, 4L); expect_identical(ret$effective_rank, 2L); expect_true(ret$settled)
	expect_length(ret$steps, 2L)
	# the weakest column is the last one (the product norms fall with the column), the start of a lower rank drops it
	expect_identical(vapply(ret$steps, function(z) z$column, 0L), c(4L, 3L))
	expect_identical(vapply(ret$steps, function(z) z$rank, 0L), c(4L, 3L))
	expect_equal(ret$steps[[1L]]$share, 1 / 16)
	expect_identical(dim(seen[[3L]]$init$U), c(6L, 2L))
	expect_equal(seen[[3L]]$init$U, init$U[, 1:2])
	# nothing retired when the iteration settles, and a mode-only fit (no variational stage) is not retired
	ok = lame:::.ae_retire_columns(function(init, Rr) mk(Rr, TRUE), NULL, 3L, FALSE, TRUE)
	expect_null(ok$retirement)
	bad = lame:::.ae_retire_columns(function(init, Rr) mk(Rr, FALSE), NULL, 3L, FALSE, FALSE)
	expect_null(bad$retirement)
	# an iteration that settles at no rank is returned unsettled at rank zero, which the caller refuses
	none = lame:::.ae_retire_columns(function(init, Rr) mk(Rr, FALSE), NULL, 3L, FALSE, TRUE)
	expect_identical(none$retirement$effective_rank, 0L); expect_false(none$retirement$settled)
})

test_that("an iteration that only ran out of a shortened budget is refused at its rank, not reduced in rank", {
	mk = function(Rr, termination) {
		U = matrix(1, 6, Rr) %*% diag(seq(Rr, 1), Rr)
		list(
			fit = list(U = U, V = U, eb_mf = list(record = list(converged = FALSE, termination = "converged"))),
			stationarity = list(pass = FALSE, termination = termination),
			retained = Rr
		)
	}
	# below the default budget a stop at the limit is the caller's bound: the iteration may only need more steps
	short = lame:::.ae_retire_columns(
		function(init, Rr) mk(Rr, "iteration limit"),
		NULL,
		3L,
		FALSE,
		TRUE,
		max_iter = 50L
	)
	expect_null(short$retirement); expect_identical(short$retained, 3L)
	# a stall (no step that improves) is a property of the problem at any budget
	stall = lame:::.ae_retire_columns(
		function(init, Rr) mk(Rr, "no improving conditional step"),
		NULL,
		3L,
		FALSE,
		TRUE,
		max_iter = 50L
	)
	expect_identical(stall$retirement$effective_rank, 0L)
	# and at the default budget the limit stop is retired too
	full = lame:::.ae_retire_columns(
		function(init, Rr) mk(Rr, "iteration limit"),
		NULL,
		3L,
		FALSE,
		TRUE,
		max_iter = 200L
	)
	expect_identical(full$retirement$effective_rank, 0L)
	# a fit with a shortened budget is refused at the rank asked for
	d = ret_sim_binary(30, 1, -1, 51)
	err = tryCatch(
		suppressMessages(
			suppressWarnings(ame_als(d$Y, Xdyad = d$X, R = 1, family = "binary", max_iter = 3, verbose = FALSE))
		),
		error = function(e) e
	)
	expect_s3_class(err, "lame_als_convergence")
	expect_match(conditionMessage(err), "(rank 1)", fixed = TRUE)
})

test_that("the weakest column of an undirected fit is the one with the smallest eigenvalue", {
	fit = list(U = diag(5)[, 1:3], L = c(-4, 0.5, 2))
	expect_equal(lame:::.ae_column_strength(fit, TRUE), c(4, 0.5, 2))
	expect_equal(
		lame:::.ae_column_strength(list(U = cbind(1:3, 1), V = cbind(1, c(2, 2, 2))), FALSE),
		sqrt(c(14 * 3, 3 * 12))
	)
	expect_identical(lame:::.ae_fit_rank(list()), 0L)
	# a random start whose last column is retired is the cold start, as the fit at rank 0 of every other start
	expect_null(lame:::.ae_init_keep(list(U = matrix(1, 3, 1), V = matrix(2, 3, 1)), integer(0)))
	expect_identical(dim(lame:::.ae_init_keep(list(U = matrix(1, 3, 2), V = matrix(2, 3, 2)), 2L)$U), c(3L, 1L))
})

test_that("a sparse rank-2 network whose posterior mode does not settle retires its weakest column and reports it", {
	skip_on_cran()
	# this rank-2 network does not settle at rank 2 under the data-scaled factor prior (R/uv_prior.R)
	d = ret_sim_binary(50, 2, -3, 309)
	said = list()
	fit = suppressMessages(withCallingHandlers(
		suppressWarnings(ame_als(d$Y, Xdyad = d$X, R = 2, family = "binary", verbose = FALSE)),
		lame_als_retirement = function(m) said[[length(said) + 1L]] <<- m
	))
	expect_true(fit$converged)
	expect_identical(fit$R, 1L)
	expect_identical(fit$R_requested, 2L)
	# the fit says at once that it is the rank-1 fit, not only when it is printed
	expect_length(said, 1L)
	expect_match(conditionMessage(said[[1L]]), "Rank 2 was requested; the data support 1")
	ret = fit$eb$retirement
	expect_identical(ret$requested_rank, 2L); expect_identical(ret$effective_rank, 1L); expect_true(ret$settled)
	expect_identical(ret$steps[[1L]]$stage, "posterior mode")
	expect_match(ret$steps[[1L]]$termination, "no improving conditional step")
	expect_lt(ret$steps[[1L]]$share, 1)
	# it is the fit of rank 1, with its evidence and its standard errors
	f1 = suppressMessages(suppressWarnings(ame_als(d$Y, Xdyad = d$X, R = 1, family = "binary", verbose = FALSE)))
	expect_equal(unname(coef(fit)), unname(coef(f1)), tolerance = 1e-8)
	expect_equal(fit$UVPM, f1$UVPM, tolerance = 1e-8)
	expect_equal(fit$eb$log_marginal, f1$eb$log_marginal, tolerance = 1e-8)
	expect_equal(unname(fit$eb$se), unname(f1$eb$se), tolerance = 1e-8)
	expect_equal(dim(fit$U), c(50L, 1L))
	expect_identical(ncol(fit$eb$factor_posterior$Sig) %||% 1L, 1L)
	# print and summary report the retirement; the methods work on the fit of the effective rank
	expect_match(
		paste(capture.output(print(fit)), collapse = "\n"),
		"Rank 2 was requested; the data support 1",
		fixed = TRUE
	)
	sm = summary(fit)
	expect_identical(sm$retirement$effective_rank, 1L)
	expect_true(all(is.finite(vcov(fit))))
	expect_equal(dim(fitted(fit)), c(50L, 50L))
	rf = suppressMessages(suppressWarnings(ame_als_refit(fit)))
	expect_identical(rf$R, 1L)
	# the product of this weak network is withheld from the report; the variational means are the fit's state
	expect_true(all(is.finite(lame:::.ae_product_restore(fit)$UVPM)))
})

test_that("a five-node pair rule that stalls is refined to seven nodes and the fit converges without retiring a column", {
	skip_on_cran()
	d = ret_sim_binary(50, 2, -3, 201)
	# the unit factor prior IW(2R + 2, (2R + 2) I) makes the five-node rule stall on this network, so it is set here
	testthat::local_mocked_bindings(
		.lame_uv_prior_moments = function(Y, family, odmax = NULL) list(vscale = 1, vdfmlt = 1),
		.package = "lame"
	)
	fit = suppressMessages(suppressWarnings(ame_als(d$Y, Xdyad = d$X, R = 2, family = "binary", verbose = FALSE)))
	expect_true(fit$converged)
	expect_identical(fit$R, 2L)
	expect_null(fit$eb$retirement)
	expect_identical(fit$eb$pair$quad, 7L)
	expect_identical(fit$eb$pair$refined, 7L)
	# a refit of the fit starts at the rule it settled on
	rf = suppressMessages(suppressWarnings(ame_als_refit(fit)))
	expect_identical(rf$eb$pair$quad, 7L)
})

test_that("an example network at rank 6 returns its rank-5 fit, and a rank that settles is not retired", {
	skip_on_cran()
	# the second Lazega network
	data("lazegalaw", package = "lame", envir = environment())
	Y = lazegalaw$Y[, , 2L]
	diag(Y) = NA
	fit = suppressMessages(suppressWarnings(ame_als(Y, R = 6, family = "binary", verbose = FALSE)))
	expect_true(fit$converged)
	expect_identical(fit$R, 5L)
	expect_identical(fit$R_requested, 6L)
	f5 = suppressMessages(suppressWarnings(ame_als(Y, R = 5, family = "binary", verbose = FALSE)))
	expect_null(f5$eb$retirement)
	expect_equal(unname(coef(fit)), unname(coef(f5)), tolerance = 1e-8)
})

test_that("a refinement that still does not settle after retiring every column is refused with the state and the requested rank", {
	rec = list(termination = "no improving step", iterations = 47L, residual = 4.9e-3, tolerance = 1e-6)
	err = tryCatch(lame:::.ae_mf_refuse(rec, 0L, 3L), error = function(e) e)
	expect_s3_class(err, "lame_als_convergence")
	expect_match(conditionMessage(err), "rank 3 was requested and 3 columns were retired")
	expect_match(conditionMessage(err), "rank 0")
	# without retirement the message does not mention it
	expect_false(grepl("retired", conditionMessage(tryCatch(lame:::.ae_mf_refuse(rec, 2L), error = function(e) e))))
})
