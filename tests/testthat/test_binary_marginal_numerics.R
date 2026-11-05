skip_on_cran()

test_that('dynamic actor trajectories remain finite and reversible with a large factor force', {
	a = readRDS(test_path('fixtures', 'binary_dynamic_hmc_large_force.rds'))
	# a stored context without the prior fields was drawn under the unit prior: 2R + 2 degrees of freedom, scale (2R + 2) I
	a$context$uv_df = a$context$uv_psi = 2 * a$context$rank + 2
	forward = lame:::.lame_binary_hmc_integrate(a$q, a$p, a$context, a$hyper, a$mass, a$step, a$L)
	expect_false(is.null(forward))
	expect_true(all(is.finite(c(forward$q, forward$p, forward$evaluation$logp, forward$evaluation$gradient))))
	backward = lame:::.lame_binary_hmc_integrate(forward$q, -forward$p, a$context, a$hyper, a$mass, a$step, a$L)
	expect_false(is.null(backward))
	expect_lt(max(abs(backward$q - a$q), abs(backward$p + a$p)), 1e-7)
})
