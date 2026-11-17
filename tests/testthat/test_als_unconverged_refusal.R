library(lame)
skip_on_cran()

# the throwaway static pass that seeds a dynamic fit is silent; the additive block is solved exactly or the
# fit stops

test_that("the throwaway static pass that seeds a dynamic fit is silent: the dynamic fit has its own criterion", {
	# the penalised dynamic fit, as the sampler's warm start asks for it
	pen = structure(list(), lame_als_estimator = "none")
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
			prior = pen,
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
		prior = pen,
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
	# the gaussian fit to the response, which seeds the penalised dynamic fit
	expect_no_warning(fit <- lame:::.ae_run(
		Y,
		NULL,
		NULL,
		NULL,
		0,
		"normal",
		"unipartite",
		FALSE,
		200,
		1e-6,
		FALSE,
		6886,
		FALSE,
		NULL,
		estimator = "transform"
	))
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
