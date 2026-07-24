# dynamic_ab fits carry per-time posterior SD of the additive effects
# (a_dynamic_sd / b_dynamic_sd), and ab_plot(plot_type="ribbon") uses them.

test_that("dynamic_ab fit stores calibrated per-time additive SD", {
	skip_on_cran()
	set.seed(11)
	n = 15; N = 5
	a_true = matrix(0, n, N); b_true = matrix(0, n, N)
	a_true[, 1] = rnorm(n, 0, 0.8); b_true[, 1] = rnorm(n, 0, 0.8)
	for (t in 2:N) {
		a_true[, t] = 0.8 * a_true[, t - 1] + rnorm(n, 0, 0.4)
		b_true[, t] = 0.8 * b_true[, t - 1] + rnorm(n, 0, 0.4)
	}
	Y = vector("list", N)
	for (t in 1:N) {
		eta = 0.2 + outer(a_true[, t], b_true[, t], "+") + matrix(rnorm(n * n, 0, 0.5), n, n)
		Yt = 1 * (eta > 0); diag(Yt) = NA
		dimnames(Yt) = list(sprintf("v%02d", 1:n), sprintf("v%02d", 1:n))
		Y[[t]] = Yt
	}
	names(Y) = paste0("t", 1:N)

	fit = lame(Y, family = "binary", dynamic_ab = TRUE, R = 0,
	           burn = 100, nscan = 600, odens = 5, verbose = FALSE, gof = FALSE)

	expect_false(is.null(fit$a_dynamic_sd))
	expect_false(is.null(fit$b_dynamic_sd))
	expect_equal(dim(fit$a_dynamic_sd), dim(fit$a_dynamic))
	expect_identical(dimnames(fit$a_dynamic_sd), dimnames(fit$a_dynamic))
	expect_true(all(fit$a_dynamic_sd > 0))
	# the SD genuinely varies over time (not a constant band)
	per_actor_range = apply(fit$a_dynamic_sd, 1, function(x) diff(range(x)))
	expect_gt(mean(per_actor_range > 1e-6), 0.5)
})

test_that("ab_plot ribbon draws a per-time-varying band from the stored SD", {
	skip_on_cran()
	nm = sprintf("v%02d", 1:6)
	fit = list(
		a_dynamic    = matrix(rnorm(6 * 4), 6, 4, dimnames = list(nm, paste0("t", 1:4))),
		a_dynamic_sd = matrix(runif(6 * 4, 0.1, 0.6), 6, 4, dimnames = list(nm, paste0("t", 1:4))),
		symmetric = FALSE)
	class(fit) = "lame"
	p = ab_plot(fit, "sender", plot_type = "ribbon", show_actors = nm[1])
	expect_s3_class(p, "ggplot")
	w = p$data
	expect_true(all(w$upper > w$lower))
	# band width tracks the per-time SD, so it is not constant
	expect_gt(diff(range(w$upper - w$lower)), 1e-6)

	# a fit without the SD slot aborts with clear guidance, not a wrong plot
	fit$a_dynamic_sd = NULL
	expect_error(ab_plot(fit, "sender", plot_type = "ribbon"),
	             "per-time uncertainty")
})
