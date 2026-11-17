# numerical failures inside the samplers: a draw that cannot be made keeps its state only when that is
# a valid move (counted in fit$exact_fallbacks and reported); otherwise it is counted and stops the fit

.fallback_bip_panel = function() {
	set.seed(11)
	lapply(1:3, function(t) {
		matrix(
			rnorm(8 * 6), 8, 6,
			dimnames = list(paste0('r', 1:8), paste0('c', 1:6))
		)
	})
}

.fallback_ord_panel = function() {
	set.seed(12)
	lapply(1:3, function(t) {
		y = matrix(sample(0:2, 64, TRUE), 8, 8, dimnames = list(paste0('a', 1:8), paste0('a', 1:8)))
		diag(y) = NA
		y
	})
}

test_that('a failed bipartite coefficient and effect draw is counted and stops the fit', {
	calls = 0L
	real = lame:::rbeta_ab_bip_gibbs_cpp
	local_mocked_bindings(
		rbeta_ab_bip_gibbs_cpp = function(...) {
			calls <<- calls + 1L
			if (calls == 5L) stop('injected bipartite draw failure')
			real(...)
		},
		.package = 'lame'
	)
	e = expect_error(
		suppressWarnings(lame(
			.fallback_bip_panel(), R = 1, family = 'normal', mode = 'bipartite',
			burn = 10, nscan = 20, odens = 2, verbose = FALSE, gof = FALSE, seed = 5
		)),
		class = 'lame_sampler_failure'
	)
	expect_gt(calls, 5L)
	expect_equal(e$health$numerical_failures[e$health$step == 'betaAB'], 1)
})

test_that('a non-finite dynamic coefficient hyperparameter draw is counted, not skipped', {
	skip_on_cran()
	local_mocked_bindings(
		sample_sigma_beta_cpp = function(...) NA_real_,
		.sample_sigma_beta_alt = function(...) NA_real_,
		.package = 'lame'
	)
	set.seed(13)
	Y = lapply(1:4, function(t) {
		y = matrix(rnorm(64), 8, 8, dimnames = list(paste0('a', 1:8), paste0('a', 1:8)))
		diag(y) = NA
		y
	})
	e = expect_error(
		suppressWarnings(lame(
			Y, R = 0, family = 'normal', dynamic_beta = TRUE,
			burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, seed = 5
		)),
		class = 'lame_sampler_failure'
	)
	expect_gt(e$health$numerical_failures[e$health$step == 'rho_sigma_beta'], 0)
})

test_that('a factor prior covariance too close to singular for chol is inverted exactly', {
	# a correlation of rank deficiency 1e-17: chol fails, the eigenvalues stay positive
	set.seed(14)
	Q = qr.Q(qr(matrix(rnorm(9), 3)))
	found = FALSE
	for (k in 15:19) {
		C = Q %*% diag(c(10^-k, 1, 2)) %*% t(Q)
		C = (C + t(C)) / 2
		b = C / sqrt(outer(diag(C), diag(C)))
		ev = eigen((b + t(b)) / 2, symmetric = TRUE, only.values = TRUE)$values
		if (inherits(try(chol((b + t(b)) / 2), silent = TRUE), 'try-error') && min(ev) > 0) {
			found = TRUE
			break
		}
	}
	skip_if_not(found, 'no test matrix with a failed chol and positive eigenvalues on this platform')
	P = lame:::.lame_factor_prior_precision(C)
	# the eigenvalue floor of the old fallback capped the precision near 1 / sqrt(eps)
	expect_gt(max(eigen(P, symmetric = TRUE, only.values = TRUE)$values), 1e12)
	expect_true(all(is.finite(P)))
})

test_that('a factor prior covariance that is not positive definite stops instead of being floored', {
	C = matrix(c(1, 1.5, 1.5, 1), 2)
	expect_error(lame:::.lame_factor_prior_precision(C), 'not numerically positive definite')
})

test_that('a failed observed custom GOF statistic is NA, not zero', {
	set.seed(15)
	Y = lapply(1:2, function(t) {
		y = matrix(rnorm(36), 6, 6, dimnames = list(paste0('a', 1:6), paste0('a', 1:6)))
		diag(y) = NA
		y
	})
	fit = lame(Y, R = 0, family = 'normal', burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, seed = 3)
	g = suppressWarnings(gof(
		fit, nsim = 3, verbose = FALSE,
		custom_gof = list(bad = function(Y) stop('boom'), fine = function(Y) 1)
	))
	G = if (is.list(g) && !is.null(g$GOF)) g$GOF else g
	if (is.list(G) && !is.matrix(G)) G = G[[1L]]
	expect_true(is.na(G['obs', 'bad']))
	expect_equal(unname(G['obs', 'fine']), 1)
})

test_that('a near-singular Matern-3/2 prior path keeps its covariance', {
	s = dynamic_beta_prior_summary(
		n_periods = 10, ndraws = 200, kind = 'matern32', matern32_length_scale = 1e6,
		sigma_shape = 50, sigma_scale = 50 * 1e8, seed = 1
	)
	tab = if (is.data.frame(s)) s else s$summary
	row = tab[grepl('first_diff', tab[[1L]]), , drop = FALSE]
	# white noise of sd about 1e4 gives first differences of that order; the kernel gives about 0.02
	expect_lt(max(abs(unlist(row[vapply(row, is.numeric, logical(1))]))), 10)
})

test_that('a failed translation step keeps its state, is counted and is reported', {
	skip_on_cran()
	real = lame:::.lame_checked_truncated_gaussian
	local_mocked_bindings(
		.lame_checked_truncated_gaussian = function(mean, sd, lower, upper) {
			in_ltm = any(vapply(sys.calls(), function(cl) {
				f = cl[[1L]]
				is.name(f) && as.character(f) %in% c('.lame_ltm_visit', '.lame_ltm_location')
			}, logical(1)))
			if (in_ltm) cli::cli_abort('forced', class = 'lame_truncation_numerical_error')
			real(mean, sd, lower, upper)
		},
		.package = 'lame'
	)
	w = NULL
	fit = withCallingHandlers(
		lame(.fallback_ord_panel(), R = 0, family = 'ordinal', burn = 10, nscan = 20, odens = 2, verbose = FALSE, gof = FALSE, seed = 4),
		lame_exact_fallbacks = function(cnd) {
			w <<- cnd
			invokeRestart('muffleWarning')
		},
		warning = function(cnd) invokeRestart('muffleWarning')
	)
	expect_s3_class(fit, 'lame')
	expect_gt(fit$exact_fallbacks[['translation']], 0)
	expect_false(is.null(w))
	expect_match(lame:::.lame_fallback_share_line(fit), 'translation')
	out = capture.output(print(summary(fit)))
	expect_true(any(grepl('kept its state', out)))
})

test_that('a fit with no numerical failure records no fallback', {
	fit = lame(.fallback_ord_panel(), R = 0, family = 'ordinal', burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, seed = 4)
	expect_length(fit$exact_fallbacks, 0L)
	expect_null(lame:::.lame_fallback_share_line(fit))
})

test_that('the fallback report warns only above a tenth of the sweeps', {
	expect_no_warning(lame:::.lame_report_fallbacks(c(translation = 10), 100))
	expect_warning(lame:::.lame_report_fallbacks(c(translation = 11), 100), class = 'lame_exact_fallbacks')
	expect_match(
		lame:::.lame_fallback_share_line(list(exact_fallbacks = c(ab_shift = 2), sampler_health_iterations = 50)),
		'ab_shift 2 times in 50 sweeps'
	)
})

test_that('a checkpoint that cannot be written warns once and the run continues', {
	skip_on_cran()
	set.seed(16)
	Y = lapply(1:2, function(t) {
		y = matrix(rnorm(36), 6, 6, dimnames = list(paste0('a', 1:6), paste0('a', 1:6)))
		diag(y) = NA
		y
	})
	# a directory in place of the file: every write fails
	path = tempfile()
	dir.create(path)
	n_warn = 0L
	fit = withCallingHandlers(
		lame(
			Y, R = 0, family = 'normal', burn = 5, nscan = 40, odens = 1, verbose = FALSE, gof = FALSE, seed = 3,
			checkpoint_path = path, checkpoint_every = 5
		),
		warning = function(cnd) {
			if (grepl('checkpoint could not be written', conditionMessage(cnd))) n_warn <<- n_warn + 1L
			invokeRestart('muffleWarning')
		}
	)
	expect_s3_class(fit, 'lame')
	expect_equal(n_warn, 1L)
})
