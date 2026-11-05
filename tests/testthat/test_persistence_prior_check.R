skip_on_cran()

# summary() measures how far the data narrow the persistences rho_uv and
# rho_ab, as it does for rho_beta: the posterior sd against the sd of the
# truncated normal prior the sampler drew them from

persistence_fit = function(seed = 11) {
	set.seed(seed); n = 8L; Tn = 3L
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(Tn), function(t) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm)); diag(y) = NA; y })
	suppressWarnings(suppressMessages(lame(
		Y,
		family = "normal",
		R = 1,
		dynamic_ab = TRUE,
		dynamic_uv = TRUE,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = seed
	)))
}

# draws from the normal prior truncated to the sampler's support
prior_draws = function(m, s, k) {
	x = rnorm(20 * k, m, s)
	x[abs(x) < .99][seq_len(k)]
}

test_that("the sd of the truncated normal prior is exact", {
	f = lame:::.lame_truncnorm_sd
	set.seed(1)
	for (pr in list(c(.9, .1), c(.8, .15), c(.5, .2), c(0, 1), c(.95, .1))) {
		x = rnorm(6e5, pr[1], pr[2]); x = x[abs(x) < .99]
		expect_equal(f(pr[1], pr[2], -.99, .99), sd(x), tolerance = .01)
	}
	# a prior far from its bounds is not truncated
	expect_equal(f(0, 1e-3, -.99, .99), 1e-3, tolerance = 1e-9)
	# the default rho_uv prior is narrower than its nominal sd: it ends at the bound
	expect_lt(f(.9, .1, -.99, .99), .09)
	expect_true(is.na(f(5, .01, -.99, .99)))
})

test_that("summary() reports the posterior sd against the prior sd of rho_uv and rho_ab", {
	fit = persistence_fit()
	sm = summary(fit)
	pc = sm$persistence_prior_check
	expect_s3_class(pc, "data.frame")
	expect_identical(pc$parameter, c("rho_uv", "rho_ab"))
	expect_named(pc, c("parameter", "post_sd", "prior_sd", "sd_ratio", "weak"))
	# the default priors of dynamic fits: rho_uv ~ N(.9, .1^2), rho_ab ~ N(.8, .15^2), both on |rho| < .99
	expect_equal(
		pc$prior_sd,
		c(
			lame:::.lame_truncnorm_sd(.9, .1, -.99, .99),
			lame:::.lame_truncnorm_sd(.8, .15, -.99, .99)
		)
	)
	expect_equal(pc$post_sd, c(sd(fit$rho_uv), sd(fit$rho_ab)))
	expect_equal(pc$sd_ratio, pc$post_sd / pc$prior_sd)
	expect_identical(pc$weak, pc$sd_ratio > .6)
	# no persistence, no check
	static = fit; static$dynamic_uv = static$dynamic_ab = FALSE
	expect_null(summary(static)$persistence_prior_check)
	# a fit without a recorded prior is summarized without the check
	old = fit; old$prior = NULL
	expect_null(summary(old)$persistence_prior_check)
	expect_no_error(capture.output(print(summary(old))))
})

test_that("a persistence the data leave as wide as its prior is said so, in plain terms", {
	fit = persistence_fit()
	set.seed(2)
	k = 3000L
	# posterior = prior: the data add nothing
	fit$rho_ab = prior_draws(.8, .15, k)
	# a posterior about as wide as 0.7 of the prior
	fit$rho_uv = rnorm(k, .9, .055)
	out = capture.output(print(summary(fit)))
	# cli wraps long lines, so the sentences are matched on the joined text
	txt = gsub("\\s+", " ", paste(out, collapse = " "))
	pc = summary(fit)$persistence_prior_check
	expect_gt(pc$sd_ratio[pc$parameter == "rho_ab"], .9)
	expect_length(grep("rho_ab: posterior sd", out), 1L)
	expect_match(
		txt,
		"rho_ab: posterior sd 0\\.[0-9]{3} is [01]\\.[0-9]{2} of the prior sd 0\\.[0-9]{3}: the data barely inform it\\."
	)
	r_uv = pc$sd_ratio[pc$parameter == "rho_uv"]
	expect_gt(r_uv, .6); expect_lt(r_uv, .8)
	expect_length(grep("rho_uv: posterior sd", out), 1L)
	expect_match(
		txt,
		"rho_uv: posterior sd 0\\.[0-9]{3} is [01]\\.[0-9]{2} of the prior sd 0\\.[0-9]{3}: the data inform it only weakly\\."
	)
	# the property is stated; no refit is asked for
	expect_false(grepl("refit", txt))
	expect_false(grepl("different prior", txt))
	# the lines sit under the dynamic component they describe
	i_ab = grep("Dynamic additive effects: enabled", out)
	expect_identical(which(grepl("rho_ab: posterior sd", out)), i_ab + 1L)
	i_uv = grep("Dynamic latent positions: enabled", out)
	expect_true(any(which(grepl("rho_uv: posterior sd", out)) > i_uv))
})

test_that("a persistence the data narrow well is not flagged", {
	fit = persistence_fit()
	set.seed(3)
	fit$rho_ab = rnorm(3000, .8, .02)
	fit$rho_uv = rnorm(3000, .9, .01)
	pc = summary(fit)$persistence_prior_check
	expect_true(all(pc$sd_ratio < .3))
	expect_false(any(pc$weak))
	out = capture.output(print(summary(fit)))
	expect_false(any(grepl("posterior sd .* of the prior sd", out)))
	# a chain that never moved has no posterior sd to compare
	fit$rho_ab = rep(.8, 3000)
	pc = summary(fit)$persistence_prior_check
	expect_identical(pc$parameter, "rho_uv")
})

test_that("a supplied prior sets the sd the posterior is compared with", {
	fit = persistence_fit()
	fit$prior$rho_ab_mean = .5; fit$prior$rho_ab_sd = .3
	set.seed(4)
	fit$rho_ab = prior_draws(.5, .3, 4000)
	pc = summary(fit)$persistence_prior_check
	row = pc[pc$parameter == "rho_ab", ]
	expect_equal(row$prior_sd, lame:::.lame_truncnorm_sd(.5, .3, -.99, .99))
	expect_equal(row$sd_ratio, 1, tolerance = .05)
	expect_true(row$weak)
})
