library(lame)
skip_on_cran()

# bootstrap of a penalised normal dynamic ALS fit under the static contract: at least 30 replicates,
# an unconverged replicate is retried at four times the budget and otherwise dropped,
# and more than 5 percent dropped is an error

dyn_panel = function(n = 6L, Tt = 4L, seed = 16L) {
	set.seed(seed)
	lapply(seq_len(Tt), function(t) {
		y = matrix(rnorm(n * n), n, n)
		diag(y) = NA
		y
	})
}

dyn_fit = function(Y, ...) {
	suppressWarnings(lame(
		Y,
		family = "normal",
		bias_correct = "none",
		R = 0,
		method = "als",
		dynamic_ab = TRUE,
		verbose = FALSE,
		seed = 16,
		...
	))
}

test_that("a dynamic bootstrap below 30 replicates is refused before anything is fitted", {
	Y = dyn_panel()
	for (n in c(1L, 10L, 29L)) {
		expect_error(dyn_fit(Y, bootstrap = n), class = "lame_als_bootstrap_replicates")
	}
	expect_error(
		lame:::lame_dynamic_als(
			Y,
			family = "normal",
			dynamic_ab = TRUE,
			bias_correct = "none",
			bootstrap = 5L,
			verbose = FALSE
		),
		class = "lame_als_bootstrap_replicates"
	)
})

test_that("the refusals of the routes without a bootstrap come before the count", {
	Y = replicate(3, { m = matrix(rbinom(36, 1, 0.4), 6, 6); diag(m) = NA; m }, simplify = FALSE)
	expect_error(
		lame:::lame_dynamic_als(Y, family = "binary", dynamic_ab = TRUE, bootstrap = 5L, verbose = FALSE),
		class = "lame_dynamic_eb_unsupported"
	)
})

test_that("a converged penalised fit bootstraps with every kept replicate converged", {
	Y = dyn_panel()
	fit = dyn_fit(Y, bootstrap = 30)
	expect_true(fit$converged)
	b = fit$bootstrap_dynamic
	expect_equal(b$R, 30L)
	expect_equal(b$n_failed, 0L)
	expect_equal(b$n_success, 30L)
	expect_true(all(b$converged))
	expect_length(b$refits, 30L)
	expect_true(all(vapply(b$refits, function(r) isTRUE(r$converged), TRUE)))
	expect_equal(dim(b$coef_paths)[3L], 30L)
	expect_equal(b$max_fail_share, 0.05)
	expect_true(isTRUE(fit$meta$uncertainty_available))
})

test_that("a replicate that has not converged is retried at four times the budget, and kept only if that converges", {
	Y = dyn_panel()
	fit = dyn_fit(Y)
	expect_true(fit$converged)
	real = lame:::lame_dynamic_als
	bud = fit$als_refit$args$max_iter
	budgets = numeric()
	k = 0L
	local_mocked_bindings(
		lame_dynamic_als = function(...) {
			a = list(...)
			out = real(...)
			# a replicate refit is the call that carries a start
			if (!is.null(a$start)) {
				budgets <<- c(budgets, a$max_iter)
				k <<- k + 1L
				# the first attempt of every replicate is reported as unfinished
				if (a$max_iter == bud) out$converged = FALSE
			}
			out
		},
		.package = "lame"
	)
	b = lame:::.lda_dynamic_bootstrap(
		fit,
		R_boot = 30L,
		type = "parametric",
		block_length = 1L,
		seed = 3,
		call_args = c(fit$als_refit$args, list(Y = Y))
	)
	expect_equal(b$n_failed, 0L)
	expect_equal(b$n_success, 30L)
	expect_equal(unique(budgets), c(bud, 4 * bud))
	expect_equal(k, 60L)
})

test_that("replicates that never converge are dropped and counted, never kept", {
	Y = dyn_panel()
	fit = dyn_fit(Y)
	real = lame:::lame_dynamic_als
	bud = fit$als_refit$args$max_iter
	rep_idx = 0L
	local_mocked_bindings(
		lame_dynamic_als = function(...) {
			a = list(...)
			out = real(...)
			if (!is.null(a$start)) {
				# replicates 5 and 20 never converge (both budgets), replicate 11 errors
				if (a$max_iter == bud) rep_idx <<- rep_idx + 1L
				if (rep_idx %in% c(5L, 20L)) out$converged = FALSE
				if (rep_idx == 11L) stop("refit failed")
			}
			out
		},
		.package = "lame"
	)
	# 30 replicates allow one dropped; 3 are dropped, so the run stops at the second
	err = tryCatch(
		lame:::.lda_dynamic_bootstrap(
			fit,
			R_boot = 30L,
			type = "parametric",
			block_length = 1L,
			seed = 3,
			call_args = c(fit$als_refit$args, list(Y = Y))
		),
		error = identity
	)
	expect_s3_class(err, "lame_als_bootstrap_convergence")
	expect_equal(err$replicate, 5L)
	expect_equal(err$n_failed, 2L)
	expect_equal(err$max_failed, 1L)
	expect_equal(err$n_total, 30L)
	# 100 replicates allow five: three dropped replicates are within the cap and are not in the summaries
	rep_idx = 0L
	b = lame:::.lda_dynamic_bootstrap(
		fit,
		R_boot = 100L,
		type = "parametric",
		block_length = 1L,
		seed = 3,
		call_args = c(fit$als_refit$args, list(Y = Y))
	)
	expect_equal(b$n_failed, 3L)
	expect_equal(b$n_success, 97L)
	expect_equal(which(nzchar(b$fail_reason)), c(5L, 11L, 20L))
	expect_equal(b$fail_reason[c(5L, 11L, 20L)], c("nonconverged", "error", "nonconverged"))
	expect_length(b$refits, 97L)
	expect_equal(dim(b$coef_paths)[3L], 97L)
	expect_true(all(b$converged))
	expect_length(b$failures, 3L)
})

test_that("a bootstrap of a dynamic fit that has not converged is refused", {
	Y = dyn_panel()
	fit = dyn_fit(Y)
	fit$converged = FALSE
	expect_error(
		lame:::.lda_dynamic_bootstrap(
			fit,
			R_boot = 30L,
			type = "parametric",
			block_length = 1L,
			seed = 3,
			call_args = c(fit$als_refit$args, list(Y = Y))
		),
		class = "lame_als_bootstrap_convergence"
	)
})

test_that("the changing-composition bootstrap keeps the missingness of the panel", {
	set.seed(17)
	sets = list(t1 = paste0("a", 1:4), t2 = paste0("a", 2:5), t3 = paste0("a", c(1, 3, 5)))
	Y = lapply(seq_along(sets), function(tt) {
		a = sets[[tt]]
		y = matrix(rnorm(length(a)^2), length(a), length(a), dimnames = list(a, a))
		diag(y) = NA_real_
		y
	})
	names(Y) = names(sets)
	fit = dyn_fit(Y, bootstrap = 30)
	expect_true(fit$converged)
	refit = fit$bootstrap_dynamic$refits[[1L]]
	expect_true(isTRUE(refit$changing_composition))
	expect_equal(is.finite(refit$Y), is.finite(fit$Y))
})

test_that("a bootstrap of a penalised non-normal dynamic fit is refused with the class of the post-hoc route", {
	set.seed(21)
	n = 8L
	Y = lapply(1:3, function(t) { m = matrix(rbinom(n * n, 1, 0.4), n, n); diag(m) = NA; m })
	for (fam in c("binary", "poisson")) {
		Yf = if (fam == "poisson") lapply(Y, function(m) m * 2) else Y
		err = tryCatch(
			lame(
				Yf,
				family = fam,
				method = "als",
				dynamic_ab = TRUE,
				bias_correct = "none",
				bootstrap = 30,
				verbose = FALSE
			),
			error = identity
		)
		expect_s3_class(err, "lame_als_bootstrap_unsupported")
	}
})

test_that("the dynamic bootstrap redraws the additive paths as the post-hoc route does", {
	Y = dyn_panel(n = 10L, Tt = 4L, seed = 5L)
	# actor levels with a spread the residual alone does not have
	set.seed(6)
	lev = rnorm(10, 0, 1)
	Y = lapply(Y, function(m) { m = m + outer(lev, lev, "+"); diag(m) = NA; m })
	fit = dyn_fit(Y)
	expect_true(fit$converged)
	s1 = lame:::.lda_boot_sims(fit, 60L, seed = 4)
	s2 = lame:::.lda_boot_sims(fit, 60L, seed = 4)
	expect_identical(s1, s2)
	expect_length(s1$Y, 60L)
	expect_equal(dim(s1$Y[[1L]]), c(10L, 10L, 4L))
	expect_true(all(is.na(s1$Y[[1L]][1, 1, ])))
	# the row mean of an actor over the replicates varies by the actor-level spread, far more than the residual alone gives
	rm1 = vapply(s1$Y, function(a) mean(a[1, -1, ]), 1)
	expect_gt(stats::sd(rm1), 3 * sqrt(fit$s2 / (9 * 4)))
	# the global random stream is left as it was
	set.seed(9); u = runif(1); set.seed(9)
	invisible(lame:::.lda_boot_sims(fit, 3L, seed = 4))
	expect_identical(runif(1), u)
})
