test_that("bivariate quadrature resolves near-perfect correlation", {
	quad = get('.lame_binary_dynamic_bvn_log_quad', asNamespace('lame'))
	for (rho in c(-.999999, -.99, -.5, 0, .5, .99, .999999)) {
		reference = log(.25 + asin(rho) / (2 * pi))
		expect_equal(quad(0, 0, rho), reference, tolerance = 1e-10)
	}
})

test_that("extreme bivariate tails agree with analytically bounded marginals", {
	quad = get('.lame_binary_dynamic_bvn_log_quad', asNamespace('lame'))
	cases = list(c(5000, -100, -.9), c(240.0123, 321.1275, -.9311163), c(20, -25, .99), c(100, -110, -.1))
	for (z in cases) {
		a = z[1];b = z[2];rho = z[3]
		if (a == 100) {
			cut = -120
			bound = exp(pnorm(cut, log.p = TRUE) - pnorm(b, log.p = TRUE)) + pnorm(
				(a - rho * cut) / sqrt(1 - rho^2),
				lower.tail = FALSE
			)
		} else if (a == 20) {
			bound = pnorm((a - rho * b) / sqrt(1 - rho^2), lower.tail = FALSE)
		} else bound = exp(pnorm(a, lower.tail = FALSE, log.p = TRUE) - pnorm(b, log.p = TRUE))
		expect_lt(bound, 1e-14)
		expect_equal(quad(a, b, rho), pnorm(b, log.p = TRUE), tolerance = 1e-8)
	}
})


test_that("the retained full-network actor proposal has positive probability", {
	quad = get('.lame_binary_dynamic_bvn_log_quad', asNamespace('lame'))
	a = 6.124086;b = -15.34894;rho = .9636428
	bound = pnorm((a - rho * b) / sqrt((1 - rho) * (1 + rho)), lower.tail = FALSE)
	expect_lt(bound, 1e-14)
	expect_equal(quad(a, b, rho), pnorm(b, log.p = TRUE), tolerance = 1e-10)
	expect_equal(quad(b, a, rho), pnorm(b, log.p = TRUE), tolerance = 1e-10)
})
