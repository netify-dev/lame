library(lame)

# each numerical failure or fallback is computed exactly, counted, or surfaced with its message;
# every test forces its path

rr_panel = function(n = 8, Tt = 3, binary = FALSE, seed = 1) {
	set.seed(seed)
	nm = paste0('a', seq_len(n))
	lapply(seq_len(Tt), function(t) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm))
		if (binary) y = 1 * (y > 0)
		diag(y) = NA
		y
	})
}

test_that('pooled diagnostics that cannot be computed warn with the error', {
	skip_on_cran()
	local_mocked_bindings(compute_mcmc_diagnostics = function(...) stop('injected diagnostics failure'), .package = 'lame')
	Y = rr_panel()[[1]]
	w = expect_warning(
		fit <- ame(Y, R = 0, family = 'normal', burn = 10, nscan = 20, odens = 2, verbose = FALSE, gof = FALSE, chains = 2),
		class = 'lame_diagnostics_unavailable'
	)
	expect_match(conditionMessage(w), 'injected diagnostics failure')
	expect_null(fit$diagnostics)
})

test_that('a separation check that cannot run says so', {
	local_mocked_bindings(reconstruct_EZ = function(...) stop('injected predictor failure'), .package = 'lame')
	w = expect_warning(
		lame:::.lame_separation_check(list(family = 'binary')),
		class = 'lame_separation_unchecked'
	)
	expect_match(conditionMessage(w), 'injected predictor failure')
})

test_that('a bipartite linear predictor that cannot be formed warns instead of vanishing', {
	local_mocked_bindings(get_EZ_bip_cpp = function(...) stop('injected assembly failure'), .package = 'lame')
	nA = 3L; nB = 2L; Tt = 2L
	rn = paste0('r', 1:nA); cn = paste0('c', 1:nB)
	X = lapply(1:Tt, function(t) array(1, c(nA, nB, 1), dimnames = list(rn, cn, 'intercept')))
	args = list(
		APS = setNames(numeric(nA), rn), BPS = setNames(numeric(nB), cn),
		UVPS = matrix(NA_real_, nA, nB), YPS = array(0, c(nA, nB, Tt)),
		BETA = matrix(0.1, 4, 1, dimnames = list(NULL, 'intercept')),
		VC = matrix(1, 4, 3, dimnames = list(NULL, c('va', 'vb', 've'))), GOF = NULL, Xlist = X,
		actorByYr = rep(list(rn), Tt), colActorByYr = rep(list(cn), Tt), start_vals = NULL,
		symmetric = FALSE, tryErrorChecks = list(), bip = TRUE, family = 'normal',
		nA = nA, nB = nB, n_time = Tt, G = diag(1)
	)
	w = expect_warning(fit <- do.call(lame:::get_fit_object, args), class = 'lame_ez_unavailable')
	expect_match(conditionMessage(w), 'injected assembly failure')
	expect_null(fit$EZ)
})

test_that('a failed bipartite dynamic position draw aborts with its message and sweep', {
	local_mocked_bindings(rUV_dynamic_bip_fc_cpp = function(...) stop('injected position failure'), .package = 'lame')
	set.seed(3)
	Y = lapply(1:3, function(t) matrix(rnorm(6 * 5), 6, 5, dimnames = list(paste0('r', 1:6), paste0('c', 1:5))))
	e = expect_error(
		suppressWarnings(lame(
			Y, family = 'normal', mode = 'bipartite', R_row = 1, R_col = 1, dynamic_uv = TRUE,
			burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1
		)),
		class = 'lame_sampler_failure'
	)
	msg = gsub('[[:space:]]+', ' ', conditionMessage(e))
	expect_match(msg, 'injected position failure')
	expect_match(msg, 'UV, first failure at sweep 1')
	expect_identical(e$first_errors$UV$sweep, 1L)
})

test_that('a non-finite additive level draw is a counted failure, never a silent zero', {
	skip_on_cran()
	expect_error(lame:::.lame_ab_level_joint(c(NaN, 1), c(1, 1), c(1, 1), c(1, 1), diag(2), symmetric = TRUE), 'Nonfinite additive level draw')
	expect_error(lame:::.lame_ab_level_joint(c(NaN, 1, 2), c(1, 1, 1), c(1, 1), c(1, 1), diag(2)), 'Nonfinite additive level draw')
	local_mocked_bindings(.lame_ab_level_draw = function(...) stop('injected level failure'), .package = 'lame')
	e = expect_error(
		suppressWarnings(lame(rr_panel(), R = 0, family = 'normal', dynamic_ab = TRUE, burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1)),
		class = 'lame_sampler_failure'
	)
	expect_match(gsub('[[:space:]]+', ' ', conditionMessage(e)), 'ab, first failure at sweep 1: injected level failure')
})

test_that('the compiled binary latent kernel falls back to the exact R kernel, counted in ame() and lame()', {
	skip_on_cran()
	local_mocked_bindings(rZ_bin_fused_cpp = function(...) stop('injected kernel failure'), .package = 'lame')
	Yb = rr_panel(binary = TRUE)
	fa = suppressWarnings(ame(Yb[[1]], R = 0, family = 'binary', burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1))
	expect_equal(unname(fa$exact_fallbacks['Z_compiled']), 15)
	fl = suppressWarnings(lame(Yb, R = 0, family = 'binary', burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1))
	expect_equal(unname(fl$exact_fallbacks['Z_compiled']), 15 * 3)
	expect_true(all(fl$sampler_health$numerical_failures == 0))
})

test_that('lame_multi issues each panel\'s warnings once, under the panel name', {
	real_lame = lame
	local_mocked_bindings(lame = function(...) {
		warning('panel issue {braces}')
		warning('panel issue {braces}')
		suppressWarnings(real_lame(...))
	}, .package = 'lame')
	Y = rr_panel(n = 6, Tt = 2)
	X = lapply(Y, function(y) array(rnorm(length(y)), c(dim(y), 1), dimnames = c(dimnames(y), list('x'))))
	ws = list()
	withCallingHandlers(
		lame_multi(list(east = Y, west = Y), list(X, X), R = 0, family = 'normal', burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1),
		warning = function(w) {
			ws[[length(ws) + 1L]] <<- w
			invokeRestart('muffleWarning')
		}
	)
	panel_ws = Filter(function(w) inherits(w, 'lame_multi_panel_warning'), ws)
	expect_length(panel_ws, 2L)
	expect_match(conditionMessage(panel_ws[[1]]), 'east')
	expect_match(conditionMessage(panel_ws[[2]]), 'west')
	expect_match(conditionMessage(panel_ws[[1]]), 'panel issue {braces}', fixed = TRUE)
	expect_length(panel_ws[[1]]$warnings, 2L)
})

test_that('an underflowed proposal of the joint factor HMC is counted as underflow, not as a bound', {
	calls = 0L
	local_mocked_bindings(
		.lame_binary_dynamic_ar_hmc_mass = function(cx, h) rep(1, 3),
		.lame_binary_dynamic_ar_hmc_evaluate = function(q, cx, h) {
			calls <<- calls + 1L
			if (calls <= 2L) list(logp = 0, gradient = rep(0, 3), bound = FALSE)
			else list(logp = -Inf, gradient = NULL, bound = FALSE)
		},
		.package = 'lame'
	)
	z = lame:::.lame_binary_dynamic_ar_hmc_step(c(0, 0, 0), list(), list(), 0.1, L = 3L)
	expect_true(z$underflow)
	expect_false(z$bound)
	expect_false(z$nonfinite)
	expect_false(z$accepted)
})

test_that('dynamic_ab with an intercept dropped as redundant fits without the rule-added path', {
	set.seed(11); n = 8; Tt = 3; rn = paste0('a', 1:n)
	x1 = lapply(1:Tt, function(t) matrix(rnorm(n * n), n, n))
	X = lapply(1:Tt, function(t) array(c(x1[[t]], rep(1, n * n)), c(n, n, 2), dimnames = list(rn, rn, c('x1', 'const'))))
	Y = lapply(1:Tt, function(t) { y = 0.5 * x1[[t]] + matrix(rnorm(n * n), n, n); diag(y) = NA; dimnames(y) = list(rn, rn); y })
	fit = suppressWarnings(lame(Y, Xdyad = X, family = 'normal', R = 0, dynamic_ab = TRUE, burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE, chains = 1))
	expect_s3_class(fit, 'lame')
	expect_false(isTRUE(fit$dynamic_intercept_from_ab))
	expect_false('intercept' %in% colnames(fit$BETA))
})

test_that('a design collinear beyond what the cross products resolve is refused, a nearly collinear one is fitted', {
	check = lame:::.lame_check_design_rank
	set.seed(2); n = 10
	x1 = matrix(rnorm(n * n), n, n); y = matrix(rnorm(n * n), n, n); diag(y) = NA
	arr = function(e) array(c(x1, x1 + e * matrix(rnorm(n * n), n, n)), c(n, n, 2), dimnames = list(NULL, NULL, c('a', 'b')))
	e = expect_error(check(arr(1e-8), y), class = 'lame_design_rank')
	expect_match(conditionMessage(e), 'numerically rank-deficient')
	expect_identical(e$aliased, 'b')
	e = expect_error(check(arr(0), y), class = 'lame_design_rank')
	expect_match(conditionMessage(e), 'exactly rank-deficient')
	expect_null(check(arr(1e-5), y))
})

test_that('the ALS solvers abort on a system they cannot decompose, with the error', {
	e = expect_error(lame:::.ae_safe_solve(matrix(NaN, 2, 2), c(1, 1)), class = 'lame_als_numerical')
	e = expect_error(lame:::.lda_safe_pinv(matrix(NaN, 2, 2)), class = 'lame_als_numerical')
	# an empty system has the empty solution
	expect_identical(lame:::.ae_safe_solve(matrix(0, 0, 0), numeric(0)), numeric(0))
	expect_identical(dim(lame:::.lda_safe_pinv(matrix(0, 0, 0))), c(0L, 0L))
})

test_that('a part of an ALS fit that cannot be computed is named in a warning and on the fit', {
	local_mocked_bindings(.ae_eb_evidence_core = function(...) stop('injected evidence failure'), .package = 'lame')
	Y = rr_panel(n = 12, Tt = 1)[[1]]
	w = expect_warning(fit <- ame_als(Y, R = 0, family = 'normal', verbose = FALSE), class = 'lame_als_numerical_fallback')
	expect_match(conditionMessage(w), 'injected evidence failure')
	expect_true('eb_evidence' %in% names(fit$numerical_failures))
	# an exact fallback is counted on the fit of the outermost call
	expect_true(lame:::.ae_notes_open())
	lame:::.ae_note_fallback('pair_rho_start')
	inner = lame:::.ae_notes_isolated({
		top = lame:::.ae_notes_open()
		lame:::.ae_note_fallback('inner_site')
		lame:::.ae_notes_attach(list(), top)
	})
	expect_equal(inner$exact_fallbacks, c(inner_site = 1))
	out = lame:::.ae_notes_attach(list(), TRUE)
	expect_equal(out$exact_fallbacks, c(pair_rho_start = 1))
})

test_that('a moving ALS intercept takes each period\'s mean of the effects, as the samplers do', {
	skip_on_cran()
	set.seed(5); n = 12; Tt = 4; rn = sprintf('r%02d', 1:n)
	al = rnorm(n, 0, 0.5)
	Y = lapply(1:Tt, function(t) {
		al <<- 0.8 * al + rnorm(n, 0, 0.3)
		y = -0.3 + outer(al, al, '+') + matrix(rnorm(n * n), n, n); diag(y) = NA; dimnames(y) = list(rn, rn); y
	})
	fit = suppressWarnings(lame(Y, family = 'normal', R = 0, dynamic_ab = TRUE, method = 'als', verbose = FALSE))
	expect_lt(max(abs(colMeans(fit$a_dynamic))), 1e-10)
	expect_lt(max(abs(colMeans(fit$b_dynamic))), 1e-10)
	expect_equal(dim(fit$eb$centre_path), c(2L, Tt))
	expect_true(all(is.finite(fit$eb$se_path['intercept', ])))
	expect_equal(predict(fit, type = 'link'), fit$EZ)
})

test_that('a static lame_als intercept is that of the centred effects, an ame_als one the population intercept', {
	set.seed(6); n = 15
	a = rnorm(n, 0, 0.7); b = rnorm(n, 0, 0.7)
	Y = lapply(1:3, function(t) { y = 0.2 + outer(a, b, '+') + matrix(rnorm(n * n), n, n); diag(y) = NA; y })
	fl = suppressWarnings(lame_als(Y, R = 0, family = 'normal', verbose = FALSE))
	expect_lt(abs(mean(fl$a)), 1e-10)
	expect_lt(abs(mean(fl$b)), 1e-10)
	fa = suppressWarnings(ame_als(Y[[1]], R = 0, family = 'normal', verbose = FALSE))
	# the population intercept carries the variance of the mean actor effects, the centred one does not
	pop = sqrt((fa$VC[['va']] + fa$VC[['vb']]) / n)
	expect_gt(fa$eb$se[['intercept']], 0.7 * pop)
	expect_lt(fl$eb$se[['intercept']], 0.5 * pop)
	# the predictor is unchanged by the centring
	ez = fl$coefficients[['intercept']] + outer(fl$a, fl$b, '+')
	expect_equal(unname(ez[row(ez) != col(ez)]), unname(fl$EZ[[1]][row(ez) != col(ez)]), tolerance = 1e-8)
})
