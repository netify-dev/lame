library(lame)
skip_on_cran()

# an unconverged fixed-effects or normal ALS fit is refused (class "lame_als_convergence") after one rerun with
# four times the iteration limit; the additive block is solved exactly or the fit stops; an unconverged
# half-sample refit of the split-sample jackknife is an error ("lame_jk_nonconvergence")

flat_binary = function(seed, n = 25, p = 0.3) {
	set.seed(seed)
	Y = matrix(rbinom(n * n, 1, p), n, n)
	diag(Y) = NA
	Y
}

test_that("an IRLS fit that stops on its iteration limit is refused, not returned at its last iterate", {
	Y = flat_binary(7, n = 30)
	err = tryCatch(
		ame_als(
			Y,
			R = 1,
			family = "binary",
			bias_correct = "none",
			max_iter = 1,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_convergence")
	expect_identical(err$kind, "irls")
	# the state and the remedy are in the message
	expect_match(conditionMessage(err), "max_iter", fixed = TRUE)
	expect_match(conditionMessage(err), "iteration limit", fixed = TRUE)
	expect_match(conditionMessage(err), "largest scaled likelihood score", fixed = TRUE)
	expect_true(is.finite(err$state$max_score))
	expect_gt(err$state$max_score, err$state$threshold)
	# the budget is raised once, from 1 to 4, before the refusal
	expect_match(conditionMessage(err), "raised from 1 to 4", fixed = TRUE)
	expect_identical(err$max_iter, 4L)
	# the same through lame_als and the lame() dispatcher
	Yl = replicate(2, flat_binary(8, n = 20), simplify = FALSE)
	expect_error(
		lame_als(
			Yl,
			R = 1,
			family = "binary",
			bias_correct = "none",
			max_iter = 1,
			verbose = FALSE
		),
		class = "lame_als_convergence"
	)
	expect_error(
		lame(
			Yl,
			R = 1,
			family = "binary",
			method = "als",
			bias_correct = "none",
			als_max_iter = 1,
			verbose = FALSE
		),
		class = "lame_als_convergence"
	)
})

test_that("a fit that needs more than its iteration limit but less than four times it is returned converged", {
	Y = flat_binary(1)
	f0 = ame_als(Y, R = 1, family = "binary", bias_correct = "none", verbose = FALSE)
	expect_true(f0$converged)
	f1 = suppressMessages(ame_als(
		Y,
		R = 1,
		family = "binary",
		bias_correct = "none",
		max_iter = 5,
		verbose = FALSE
	))
	expect_true(f1$converged)
	expect_true(f1$irls_diagnostics$pass)
	# the escalated run takes more reweighting steps than its given limit
	expect_gt(f1$irls_iterations, 5L)
	expect_equal(unname(coef(f1)), unname(coef(f0)), tolerance = 1e-4)
	# the call recorded on the fit is the caller's
	expect_identical(f1$call$max_iter, 5)
	expect_identical(f1$call[[1L]], quote(ame_als))
})

test_that("a warning of the first run that the escalated run repeats is shown once, not twice", {
	# two actor blocks with no observed dyad between them: the graph warning is raised before the fit,
	# in the refused run and again in the escalated one
	set.seed(1)
	n = 24
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA
	Y[1:12, 13:24] = NA
	Y[13:24, 1:12] = NA
	w = character()
	fit = withCallingHandlers(
		suppressMessages(ame_als(
			Y,
			R = 1,
			family = "binary",
			bias_correct = "none",
			max_iter = 5,
			verbose = FALSE
		)),
		warning = function(x) {
			w <<- c(w, conditionMessage(x))
			invokeRestart("muffleWarning")
		}
	)
	expect_true(fit$converged)
	expect_gt(fit$irls_iterations, 5L)
	expect_equal(sum(grepl("disconnected components", w)), 1L)
})

test_that("a two-cycle of the factor-scale updates is damped and the fit converges", {
	# a rank-one term on pure noise at density 0.1, where the undamped penalised deviance cycles between two values
	Y = flat_binary(4, p = 0.1)
	fit = suppressMessages(ame_als(Y, R = 1, family = "binary", bias_correct = "none", verbose = FALSE))
	expect_true(fit$converged)
	expect_identical(fit$irls_diagnostics$termination, "converged")
	expect_lt(fit$irls_iterations, 200L)
})

test_that("a Gaussian fit whose first-order score stays above its threshold is refused with the score", {
	set.seed(5)
	n = 20
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	err = tryCatch(
		ame_als(
			Y,
			R = 2,
			family = "normal",
			bias_correct = "none",
			max_iter = 1,
			tol = 1e-12,
			verbose = FALSE
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_convergence")
	expect_identical(err$kind, "bcd")
	expect_match(conditionMessage(err), "Block coordinate descent did not converge", fixed = TRUE)
	expect_match(conditionMessage(err), "max_iter", fixed = TRUE)
	expect_true(is.finite(err$state$max_score))
	expect_true(is.finite(err$state$threshold))
})

test_that("the throwaway static pass that seeds a dynamic fit is silent: the dynamic fit has its own criterion", {
	set.seed(9110)
	Y = replicate(3, {
		m = matrix(rbinom(64, 1, 0.35), 8, 8)
		diag(m) = NA
		m
	}, simplify = FALSE)
	warns = character()
	record = function(expr) withCallingHandlers(expr, warning = function(w) {
		warns <<- c(warns, conditionMessage(w))
		invokeRestart("muffleWarning")
	})
	# an unfinished dynamic fit is refused, with the last iterate in the condition and nothing shown before it
	err = record(expect_error(
		lame(
			Y,
			family = "binary",
			R = 0,
			method = "als",
			dynamic_ab = TRUE,
			als_max_iter = 1,
			bias_correct = "none",
			verbose = FALSE
		),
		class = "lame_als_convergence"
	))
	expect_length(warns, 0L)
	expect_false(err$fit$converged)
	expect_false(err$fit$diagnostics$static_seed$converged)
	# a dynamic fit that finishes is returned whatever its static seed did, with its stop recorded, not shown
	Yn = lapply(Y, function(m) { m[!is.na(m)] = rnorm(sum(!is.na(m))); m })
	fit = record(lame(
		Yn,
		family = "normal",
		R = 0,
		method = "als",
		dynamic_ab = TRUE,
		als_max_iter = 1,
		bias_correct = "none",
		verbose = FALSE
	))
	expect_length(warns, 0L)
	expect_true(fit$converged)
	expect_false(fit$diagnostics$static_seed$converged)
	expect_gt(fit$diagnostics$static_seed$iterations, 0L)
})

# the additive block

test_that("a weakly connected observed-dyad graph gets the exact additive effects, with no message about approximation", {
	# two dense blocks of senders and receivers joined by a single observed cell
	set.seed(31)
	n = 24
	a = rnorm(n, 0, 0.8); b = rnorm(n, 0, 0.8)
	Y = outer(a, b, "+") + 0.4 + matrix(rnorm(n * n, 0, 0.3), n, n)
	rows1 = 1:12; rows2 = 13:24
	Y[rows1, rows2] = NA
	Y[rows2, rows1] = NA
	Y[3, 15] = 0.4 + a[3] + b[15]
	diag(Y) = NA
	expect_no_warning(fit <- ame_als(Y, R = 0, family = "normal", bias_correct = "none", verbose = FALSE))
	expect_true(fit$converged)
	# oracle: two-way additive least squares on the observed cells
	df = data.frame(y = c(Y), i = factor(c(row(Y))), j = factor(c(col(Y))))
	df = df[is.finite(df$y), ]
	oracle = lm(y ~ i + j, data = df)
	fv = unlist(fit$EZ)[is.finite(c(Y))]
	expect_equal(unname(fv), unname(fitted(oracle)), tolerance = 1e-6)
	expect_false(any(grepl("approximate", unlist(lapply(fit$meta, as.character)))))
})

test_that("the exact additive solve reaches the normal equations on unbalanced and weakly connected weights", {
	solver = lame:::.ae_additive_exact_solver
	resid = function(w, rs, cs, a, b) {
		r1 = rowSums(w) * a + as.vector(w %*% b) - rs
		r2 = colSums(w) * b + as.vector(crossprod(w, a)) - cs
		max(abs(c(r1[rowSums(w) > 0], r2[colSums(w) > 0])))
	}
	set.seed(2)
	cases = list(
		unbalanced = {w = matrix(exp(rnorm(30 * 30, 0, 4)), 30, 30); diag(w) = 0; w},
		bridge = {w = matrix(0, 20, 20); w[1:10, 1:10] = 1; w[11:20, 11:20] = 1
			w[5, 15] = 1e-9; w[15, 5] = 1e-9; diag(w) = 0; w},
		sparse = {w = matrix(rbinom(25 * 25, 1, 0.2), 25, 25); diag(w) = 0; w}
	)
	for (nm in names(cases)) {
		w = cases[[nm]]
		a0 = rnorm(nrow(w)); b0 = rnorm(ncol(w))
		a0[rowSums(w) == 0] = 0; b0[colSums(w) == 0] = 0
		rs = rowSums(w) * a0 + as.vector(w %*% b0)
		cs = colSums(w) * b0 + as.vector(crossprod(w, a0))
		sol = solver(w, rowSums(w), colSums(w), FALSE)(rs, cs)
		expect_lt(resid(w, rs, cs, sol$a, sol$b) / max(1, abs(rs), abs(cs)), 1e-9)
	}
})

test_that("an additive solve that cannot reach its normal equations stops the fit with a numerical error", {
	set.seed(3)
	n = 12
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	local_mocked_bindings(
		ae_additive_cpp = function(wsum, base_rs, base_cs, a0, b0, cnt_row, cnt_col, symmetric, direct) {
			list(a = a0, b = b0, resid = 1e3, scale = 1)
		},
		ae_additive_resid_cpp = function(...) 1e3,
		.package = "lame"
	)
	err = tryCatch(ame_als(Y, R = 0, family = "normal", verbose = FALSE), error = identity)
	expect_s3_class(err, "lame_als_numerical")
	expect_match(conditionMessage(err), "could not be solved exactly", fixed = TRUE)
})

# the split-sample jackknife

sim_jk_panel = function(nA = 24, nB = 16, Tt = 5, seed = 3) {
	set.seed(seed)
	rn = sprintf("r%02d", seq_len(nA)); cn = sprintf("c%02d", seq_len(nB))
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	P = tcrossprod(matrix(rnorm(nA, 0, 0.6), nA, 1), matrix(rnorm(nB, 0, 0.6), nB, 1))
	Y = X = vector("list", Tt)
	for (t in seq_len(Tt)) {
		x = matrix(rnorm(nA * nB), nA, nB)
		z = -0.5 + 0.7 * x + outer(a, b, "+") + P + matrix(rnorm(nA * nB), nA, nB)
		y = 1 * (z > 0)
		dimnames(y) = list(rn, cn)
		Y[[t]] = y
		X[[t]] = array(x, c(nA, nB, 1), dimnames = list(rn, cn, "x1"))
	}
	list(Y = Y, X = X)
}

jk_fit = function(d, ...) {
	suppressMessages(lame_als(
		d$Y,
		d$X,
		R = 1,
		family = "binary",
		mode = "bipartite",
		bias_correct = "jackknife",
		verbose = FALSE,
		...
	))
}

test_that("every half of a split-sample jackknife has converged, and says so", {
	fit = jk_fit(sim_jk_panel())
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_true(all(fit$bias_correct$converged))
	expect_true(fit$converged)
})

test_that("a half-sample refit that does not converge is an error of class lame_jk_nonconvergence", {
	d = sim_jk_panel()
	orig = lame:::.ae_jk_half
	calls = 0L
	local_mocked_bindings(
		.ae_jk_half = function(...) {
			calls <<- calls + 1L
			out = orig(...)
			out$converged = FALSE
			out
		},
		.package = "lame"
	)
	err = tryCatch(jk_fit(d), error = identity)
	expect_s3_class(err, "lame_jk_nonconvergence")
	expect_identical(err$half, 1L)
	# the half is continued once, from its own last iterate, with five times the budget
	expect_length(err$attempts, 2L)
	expect_identical(err$attempts[[2L]]$max_iter, 5L * err$attempts[[1L]]$max_iter)
	expect_equal(calls, 2L)
	expect_match(conditionMessage(err), "max_iter", fixed = TRUE)
	# no fit comes back
	expect_error(jk_fit(d), class = "lame_jk_nonconvergence")
})

test_that("a half that converges on its second attempt is used, and the attempts are recorded", {
	d = sim_jk_panel()
	orig = lame:::.ae_jk_half
	n_call = 0L
	local_mocked_bindings(
		.ae_jk_half = function(...) {
			n_call <<- n_call + 1L
			out = orig(...)
			# the first attempt of the first half reports as unfinished
			if (n_call == 1L) out$converged = FALSE
			out
		},
		.package = "lame"
	)
	fit = jk_fit(d)
	expect_identical(fit$bias_correct$method, "jackknife")
	expect_true(all(fit$bias_correct$converged))
	expect_equal(n_call, 3L)
})
