skip_on_cran()

test_that("observed binary scores retain both probit tails and infinity limits", {
	y = rep(c(0, 1), each = 5)
	eta = rep(c(-40, -10, 0, 10, 40), 2)
	# independent high-accuracy normal-tail reference, no probability clipping
	expected = c(
		-0, -7.619853024160473e-24, -0.6931471805599453,
		-53.23128515051248, -804.6084420137539,
		-804.6084420137539, -53.23128515051248, -0.6931471805599453,
		-7.619853024160473e-24, -0
	)
	expect_equal(lame:::.binary_pointwise_loglik(y, eta), expected, tolerance = 1e-12)
	idx = matrix(1, length(y), 3)
	expect_equal(
		lame:::.pointwise_loglik_observed_ghk(
			y,
			eta,
			eta,
			idx,
			"binary",
			1
		),
		expected,
		tolerance = 1e-12
	)
	expect_identical(
		lame:::.binary_pointwise_loglik(
			c(0, 1, 0, 1),
			c(-Inf, Inf, Inf, -Inf)
		),
		c(0, 0, -Inf, -Inf)
	)
})

test_that("AME saved scores use stable tails and the actual normal variance", {
	Y = matrix(rep(c(0, 1), 3), 2, 3)
	eta = c(-Inf, -40, -10, 0, 10, 40, Inf)
	fit = list(
		BETA = matrix(eta, ncol = 1, dimnames = list(NULL, "intercept")),
		VC = matrix(1, length(eta), 1, dimnames = list(NULL, "ve")),
		APM = numeric(2),
		BPM = numeric(3)
	)
	ll = lame:::.ame_compute_log_lik_post_hoc(fit, Y, "binary")
	expected = sapply(
		as.numeric(Y),
		function(y)
			pnorm(if (y) eta else -eta, log.p = TRUE)
	)
	expect_equal(ll, expected, tolerance = 1e-12)
	expect_false(anyNA(ll))
	fit$BETA = matrix(0, 1, 1, dimnames = list(NULL, "intercept"))
	for (v in c(1e-20, 1e-100, 1e-300)) {
		fit$VC = matrix(v, 1, 1, dimnames = list(NULL, "ve"))
		for (z in c(0, 1, 3)) {
			ll = lame:::.ame_compute_log_lik_post_hoc(
				fit,
				matrix(z * sqrt(v), 2, 3),
				"normal"
			)
			expect_equal(
				as.numeric(ll),
				rep(-.5 * log(2 * pi) - .5 * log(v) - .5 * z^2, 6),
				tolerance = 1e-12
			)
		}
	}
	fit$VC[, ] = 0
	expect_error(lame:::.ame_compute_log_lik_post_hoc(fit, Y, "normal"), "positive residual variance")
})

test_that("the saved LAME binary block routes to the same stable score", {
	find = function(x) {
		if (!is.call(x)) return(NULL)
		if (identical(x[[1]], as.name("if")) &&
				identical(
					paste(deparse(x[[2]]), collapse = " "),
					'family == "binary" || family == "cbin"'
				) &&
				grepl("ll_row", paste(deparse(x[[3]]), collapse = " "), fixed = TRUE)) return(x[[3]])
		for (y in as.list(x)) {
			z = tryCatch(find(y), error = function(e) NULL)
			if (!is.null(z)) return(z)
		}
		NULL
	}
	expr = find(body(lame::lame))
	expect_false(is.null(expr))
	e = new.env(parent = asNamespace("lame"))
	e$y_obs = c(1, 0, 1, 0); e$ez_obs = c(-40, 40, Inf, -Inf)
	eval(expr, e)
	expect_equal(e$ll_row, c(-804.6084420137539, -804.6084420137539, 0, 0), tolerance = 1e-12)
})
