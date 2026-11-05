test_that("explicit ordinal Gaussian truncation retains positive tail probabilities", {
	mu = matrix(c(40, 40, -40), 1, 3)
	y = matrix(1:3, 1, 3)
	set.seed(47327711)
	x = replicate(1000, lame:::rZ_ord_explicit_fc(matrix(0, 1, 3), mu, y, c(0, 1)))
	x = matrix(x, nrow = 3)
	expect_true(all(is.finite(x)))
	expect_true(all(x[1, ] < 0 & x[2, ] > 0 & x[2, ] < 1 & x[3, ] > 1))
	# independent conditional CDF in its lower-tail and upper-tail forms
	u1 = exp(pnorm(x[1, ] - 40, log.p = TRUE) - pnorm(-40, log.p = TRUE))
	u3 = 1 - exp(pnorm(x[3, ] + 40, lower.tail = FALSE, log.p = TRUE) -
			pnorm(41, lower.tail = FALSE, log.p = TRUE))
	expect_lt(abs(mean(u1) - .5), 5 / sqrt(12 * ncol(x)))
	expect_lt(abs(mean(u3) - .5), 5 / sqrt(12 * ncol(x)))
	expect_identical(lame:::.ord_marginal_loglik(c(-Inf, 0, 0, Inf), 0, 2), -Inf)
})

test_that("symmetric array cutpoints count each unordered pair at its averaged predictor", {
	y = matrix(c(NA, 1, 2, 3, 1, NA, 3, 2, 2, 3, NA, 1, 3, 2, 1, NA), 4)
	y = array(rep(y, 2), c(4, 4, 2))
	eta = array(seq_len(32) / 80, dim(y))
	mask = array(upper.tri(matrix(0, 4, 4)), dim(y))
	em = (eta + aperm(eta, c(2, 1, 3))) / 2
	cat = y[mask]
	loglik = function(a) {
		cuts = c(-Inf, a, Inf)
		sum(log(pnorm(cuts[cat + 1] - em[mask]) - pnorm(cuts[cat] - em[mask])))
	}
	for (seed in 47327712:47327731) {
		set.seed(seed)
		dp = .6 * rnorm(1)
		ap = c(0, exp(dp))
		accept = log(runif(1)) < loglik(ap) - loglik(c(0, 1)) + dp
		set.seed(seed)
		a = lame:::sample_alpha_cowles(c(0, 1), y, eta, .6, TRUE)
		expect_identical(a$accept, accept)
		expect_equal(a$alpha, if (accept) ap else c(0, 1), tolerance = 1e-13)
	}
})

test_that("explicit ordinal supported correlation preserves independent model restrictions", {
	set.seed(47327732)
	before = .Random.seed
	expect_silent(lame:::.lame_explicit_ordinal_support(
		"ordinal", "explicit", FALSE,
		FALSE, TRUE
	))
	expect_error(
		lame:::.lame_explicit_ordinal_support(
			"ordinal", "explicit", FALSE,
			FALSE, FALSE, .2
		),
		class = "lame_explicit_ordinal_correlation"
	)
	expect_error(
		lame:::.lame_explicit_ordinal_support(
			"ordinal", "explicit", TRUE,
			FALSE, FALSE
		),
		class = "lame_explicit_ordinal_mode"
	)
	expect_silent(lame:::.lame_explicit_ordinal_support(
		"ordinal", "explicit", FALSE,
		TRUE, TRUE
	))
	expect_identical(.Random.seed, before)
})
