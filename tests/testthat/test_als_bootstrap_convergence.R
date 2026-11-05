.bootstrap_convergence_fixture = function() {
	set.seed(813114L);y = lapply(1:3, function(t) {z = matrix(rnorm(144), 12);diag(z) = NA_real_;z})
	lame_als(y, R = 0, family = 'normal', verbose = FALSE)
}

test_that("bootstrap retries the same data and returns only converged refits", {
	f = .bootstrap_convergence_fixture();data = list();budgets = integer()
	local_mocked_bindings(ame_als_refit = function(object, Y_new, max_iter, ...) {
		data[[length(data) + 1L]] <<- Y_new;budgets <<- c(budgets, max_iter)
		object$converged = max_iter >= 40L;object
	}, .package = 'lame')
	set.seed(461L);before = .Random.seed
	b = ame_als_bootstrap(f, R = 12, max_iter = 10L, seed = 92511L, verbose = FALSE)
	expect_identical(.Random.seed, before)
	expect_true(all(b$boot_converged))
	expect_equal(b$n_nonconverged, 0)
	expect_equal(b$boot_refit_attempts, rep(2L, 12))
	expect_equal(b$boot_refit_iterations, rep(40L, 12))
	expect_equal(budgets, rep(c(10L, 40L), 12))
	for (k in 1:12) expect_identical(data[[2 * k - 1]], data[[2 * k]])
})

test_that("a bootstrap cannot turn unfinished refits into intervals", {
	f = .bootstrap_convergence_fixture();calls = 0L
	local_mocked_bindings(ame_als_refit = function(object, ...) {
		calls <<- calls + 1L;object$converged = FALSE;object
	}, .package = 'lame')
	error = tryCatch(ame_als_bootstrap(f, R = 12, max_iter = 2L, seed = 411, verbose = FALSE), error = identity)
	expect_s3_class(error, 'lame_als_bootstrap_convergence')
	# the warm start and the cold start are each tried at the budget and at four times it
	expect_equal(error$replicate, 1L);expect_false(error$fit$converged);expect_equal(calls, 4L)
	# a share of 12 replicates above the 5 percent cap stops the run at the first one that is dropped
	expect_equal(error$n_failed, 1L);expect_equal(error$max_failed, 0L);expect_equal(error$n_total, 12L)
	f$converged = FALSE
	expect_error(ame_als_bootstrap(f, R = 12, verbose = FALSE), class = 'lame_als_bootstrap_convergence')
	expect_equal(calls, 4L)
})

test_that("the bootstrap alias honors explicit refit controls", {
	f = .bootstrap_convergence_fixture();budgets = integer();tolerances = numeric()
	local_mocked_bindings(ame_als_refit = function(object, max_iter, tol, ...) {
		budgets <<- c(budgets, max_iter);tolerances <<- c(tolerances, tol);object
	}, .package = 'lame')
	b = boot_ame(f, R = 12, max_iter = 17L, tol = 1e-8, seed = 411, verbose = FALSE)
	expect_equal(budgets, rep(17L, 12));expect_equal(tolerances, rep(1e-8, 12))
	expect_equal(b$boot_refit_attempts, rep(1L, 12))
})
