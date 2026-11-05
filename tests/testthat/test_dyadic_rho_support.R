skip_on_cran()

test_that("dyadic correlation helpers preserve the full documented support", {
	rho = c(-.9998, .9997)
	expect_identical(lame:::.lame_expand_rho(rho, 2L), rho)
	expect_identical(lame:::.lame_expand_rho(.9998, 2L), rep(.9998, 2L))
	expect_identical(lame:::.lame_rho_scalar(.9998), .9998)
	set.seed(3291)
	before = .Random.seed
	for (bad in list(NA_real_, NaN, Inf, -1, 1, numeric(0), c(.2, .3, .4))) {
		expect_error(lame:::.lame_expand_rho(bad, 2L), "Dyadic correlations")
		expect_identical(.Random.seed, before)
	}
})

test_that("near-boundary correlation moves match a covariance-density oracle", {
	skip_on_cran()
	oracle = function(z, rho, prior) {
		u = upper.tri(z)
		a = z[u]; b = t(z)[u]
		proposal_sd = 2 * (1 - cor(a, b)^2) / sqrt(length(a))
		proposal = rho + proposal_sd * qnorm(runif(
			1,
			pnorm((-1 - rho) / proposal_sd),
			pnorm((1 - rho) / proposal_sd)
		))
		if (!is.finite(proposal) || abs(proposal) >= 1) return(rho)
		log_density = function(r) {
			covariance = matrix(c(1, r, r, 1), 2)
			factor = chol(covariance)
			sum(vapply(seq_along(a), function(k) {
				-sum(log(diag(factor))) -
					sum(backsolve(factor, c(a[k], b[k]), transpose = TRUE)^2) / 2
			}, numeric(1))) + prior * log(det(covariance))
		}
		log_normalizer = function(r) log(pnorm((1 - r) / proposal_sd) -
				pnorm((-1 - r) / proposal_sd))
		acceptance = log_density(proposal) - log_density(rho) +
			log_normalizer(rho) - log_normalizer(proposal)
		if (log(runif(1)) < acceptance) proposal else rho
	}
	z = matrix(c(NA, .2, .4, .21, NA, -.3, .39, -.28, NA), 3)
	# ame() and lame() share one prior: the default of rrho_mh() and of the compiled lame() kernel
	package_prior = lame:::.lame_rho_prior_power()
	expect_identical(package_prior, 1)
	for (rho in c(-.9998, .9998)) for (seed in 1:12) {
		for (prior in c(-.5, 0, 1)) {
			set.seed(47332900 + seed)
			expected = oracle(z, rho, prior)
			end = .Random.seed
			# rrho_mh(): the package prior by default, the arc-sine (asp = TRUE) and uniform (asp = FALSE) priors on request
			set.seed(47332900 + seed)
			actual = rrho_mh(z, rho, asp = switch(as.character(prior), `-0.5` = TRUE, `0` = FALSE, `1` = NULL))
			expect_equal(actual, expected, tolerance = 1e-12)
			expect_identical(.Random.seed, end)
			# the compiled kernel of lame(), given the same exponent
			set.seed(47332900 + seed)
			actual = lame:::rrho_mh_rep_cpp(array(z, c(3, 3, 1)), rho, 1, 0, prior)
			expect_equal(actual, expected, tolerance = 1e-12)
			expect_identical(.Random.seed, end)
		}
		# and without the argument it is the package prior
		set.seed(47332900 + seed)
		default = rrho_mh(z, rho)
		set.seed(47332900 + seed)
		compiled = lame:::rrho_mh_rep_cpp(array(z, c(3, 3, 1)), rho, 1, 0)
		set.seed(47332900 + seed)
		expect_identical(compiled, lame:::rrho_mh_rep_cpp(array(z, c(3, 3, 1)), rho, 1, 0, package_prior))
		set.seed(47332900 + seed)
		expect_equal(default, oracle(z, rho, package_prior), tolerance = 1e-12)
	}
})

test_that("collapsed AME reciprocity proposes from the prior it shares with lame()", {
	cx = lame:::.lame_binary_dynamic_collapsed_rho_context(array(NA_real_, c(2, 2, 1)))
	set.seed(47332953)
	expected = replicate(100, {
		r = 2 * rbeta(1, 2, 2) - 1
		runif(1)
		r
	})
	end = .Random.seed
	set.seed(47332953)
	actual = replicate(
		100,
		lame:::.ame_binary_collapsed_rho_mh(
			.9998,
			array(0, c(2, 2, 1)),
			cx,
			tries = 1L
		)$rho
	)
	expect_identical(actual, expected)
	expect_identical(.Random.seed, end)
	# the same draws, step for step, as the collapsed update of lame()
	set.seed(47332953)
	shared = replicate(
		100,
		lame:::.lame_binary_dynamic_collapsed_rho_mh(
			.9998,
			array(0, c(2, 2, 1)),
			cx,
			tries = 1L
		)$rho
	)
	expect_identical(shared, expected)
	# the law is 2 Beta(2, 2) - 1: no mass piled at the ends as the arc-sine prior has
	set.seed(47332954)
	big = replicate(20000, lame:::.lame_rho_prior_draw())
	expect_lt(abs(mean(big)), 0.02)
	expect_lt(abs(var(big) - 0.2), 0.01)
	expect_lt(mean(abs(big) > .95), 0.01)
})

test_that("dyadic variance arithmetic preserves representable extreme scales", {
	extreme = array(c(NA, 1e154, 1e154, NA), c(2, 2, 1))
	set.seed(47332961)
	gamma = rgamma(1, 1.5)
	set.seed(47332961)
	draw = lame:::.lame_sample_s2_dynamic_rho(extreme, .999, s20 = .7)
	expect_equal(draw / (1e308 / 1.999 / gamma), 1, tolerance = 2e-15)
	near = 1 - .Machine$double.eps / 2
	extreme = array(c(NA, 1e-168, -1e-168, NA), c(2, 2, 1))
	rate = (1e-168 / sqrt(1 - near))^2 + 1e-320 / 2
	set.seed(47332962)
	gamma = rgamma(1, 1.5)
	set.seed(47332962)
	draw = lame:::.lame_sample_s2_dynamic_rho(extreme, near, s20 = 1e-320)
	expect_identical(draw, rate / gamma)
	expect_gt(draw, 0)
	reflected = extreme
	reflected[1, 2, 1] = -reflected[1, 2, 1]
	set.seed(47332962)
	expect_identical(lame:::.lame_sample_s2_dynamic_rho(reflected, -near, s20 = 1e-320), draw)
	set.seed(3295)
	prior_draw = .7 / 2 / rgamma(1, .5)
	set.seed(3295)
	expect_identical(
		lame:::.lame_sample_s2_dynamic_rho(
			array(NA_real_, c(2, 2, 1)),
			.9998,
			s20 = .7
		),
		prior_draw
	)
	before = .Random.seed
	expect_error(
		lame:::.lame_sample_s2_dynamic_rho(
			array(NA_real_, c(2, 2, 1)),
			.9998,
			s20 = .Machine$double.xmin * .Machine$double.eps
		),
		"nonpositive rate"
	)
	expect_identical(.Random.seed, before)
	for (sign in c(-1, 1)) {
		huge = matrix(c(NA, sign * .999e154, 1e154, NA), 2)
		expect_true(is.finite(rrho_mh(huge, sign * .999)))
		expect_true(is.finite(lame:::rrho_mh_rep_cpp(array(huge, c(2, 2, 1)), sign * .999, 1, 0)))
	}
	bad = extreme
	bad[1, 2, 1] = Inf
	set.seed(3293)
	before = .Random.seed
	expect_error(lame:::.lame_sample_s2_dynamic_rho(bad, near), "non-finite")
	expect_identical(.Random.seed, before)
	expect_error(rrho_mh(bad[, , 1], .2), "non-finite")
	expect_identical(.Random.seed, before)
	expect_error(lame:::rrho_mh_rep_cpp(bad, .2, 1), "non-finite")
	expect_identical(.Random.seed, before)
})


test_that("improper exact-boundary conditionals fail before drawing", {
	set.seed(47332964)
	before = .Random.seed
	# the prior (1 - rho^2)^power leaves the conditional improper from 2 + 2 power exactly degenerate dyads
	expect_identical(lame:::.lame_rho_improper_dyads(), 4L)
	expect_identical(lame:::.lame_rho_improper_dyads(-0.5), 1L)
	expect_identical(lame:::.lame_rho_improper_dyads(0), 2L)
	for (sign in c(-1, 1)) {
		two = matrix(c(NA, sign * .2, .2, NA), 2)
		three = matrix(0, 3, 3)
		three[upper.tri(three)] = .2
		three[lower.tri(three)] = sign * .2
		big = matrix(0, 4, 4)
		big[upper.tri(big)] = .2
		big[lower.tri(big)] = sign * .2
		# one dyad: improper under the arc-sine prior only, three dyads: proper under the package prior and the uniform one
		expect_error(rrho_mh(two, .2, asp = TRUE), "improper")
		expect_error(rrho_mh(three, .2, asp = FALSE), "improper")
		expect_error(rrho_mh(three, .2, asp = TRUE), "improper")
		expect_identical(.Random.seed, before)
		expect_true(is.finite(rrho_mh(two, .2, asp = FALSE)))
		expect_true(is.finite(rrho_mh(two, .2)))
		expect_true(is.finite(rrho_mh(three, .2)))
		set.seed(47332964)
		expect_error(rrho_mh(big, .2), "improper")
		expect_identical(.Random.seed, before)
		# the compiled kernel agrees, at each exponent
		expect_error(lame:::rrho_mh_rep_cpp(array(big, c(4, 4, 1)), .2, 1), "improper")
		expect_error(lame:::rrho_mh_rep_cpp(array(two, c(2, 2, 1)), .2, 1, 0.1, -0.5), "improper")
		expect_error(lame:::rrho_mh_rep_cpp(array(three, c(3, 3, 1)), .2, 1, 0.1, 0), "improper")
		expect_identical(.Random.seed, before)
		expect_true(is.finite(lame:::rrho_mh_rep_cpp(array(three, c(3, 3, 1)), .2, 1)))
		set.seed(47332964)
	}
	# an exponent that leaves the prior non-integrable is refused
	expect_error(lame:::rrho_mh_rep_cpp(array(0.1, c(3, 3, 1)), .2, 1, 0.1, -1), "greater than -1")
	expect_error(rrho_mh(matrix(0.1, 3, 3), .2, asp = NA), "NULL or a single logical")
})


test_that("lost standardized residuals are diagnosed before drawing", {
	z = matrix(c(NA, 1e-300, 2e-300, NA), 2)
	set.seed(3296)
	before = .Random.seed
	expect_error(rrho_mh(z, .2, 1e300), "cannot be represented")
	expect_identical(.Random.seed, before)
	expect_error(lame:::rrho_mh_rep_cpp(array(z, c(2, 2, 1)), .2, 1e300), "cannot be represented")
	expect_identical(.Random.seed, before)
})
