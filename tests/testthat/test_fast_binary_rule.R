# the fast route's density warning for binary networks (R/fast_binary_rule.R): a user-facing fast call on
# a network with min(d, 1 - d) < 0.20 warns once and returns the fast fit; internal fits do not warn

rule_net = function(n = 12, p = 0.08, seed = 1, nc = n) {
	set.seed(seed)
	rn = sprintf("a%02d", seq_len(n)); cn = if (nc == n) rn else sprintf("c%02d", seq_len(nc))
	y = matrix(rbinom(n * nc, 1, p), n, nc, dimnames = list(rn, cn))
	if (nc == n) diag(y) = NA
	list(y = y, x = array(rnorm(n * nc), c(n, nc, 1), dimnames = list(rn, cn, "x")))
}

# the number of density warnings an expression gives, and its value; other conditions are muffled
count_density_warnings = function(expr) {
	n = 0
	value = suppressMessages(suppressWarnings(withCallingHandlers(
		expr,
		lame_fast_binary_density = function(w) {
			n <<- n + 1
			invokeRestart("muffleWarning")
		}
	)))
	list(n = n, value = value)
}

finite_se = function(fit) {
	se = sqrt(diag(vcov(fit)))
	length(se) > 0 && all(is.finite(se))
}

test_that("the density pools the observed cells, leaves out the diagonal of a one-mode network and reads both sides", {
	y = matrix(c(NA, 1, 0, 0, NA, 0, 1, 0, NA), 3, 3)
	expect_equal(lame:::.lame_fast_binary_density(y), 2 / 6)
	expect_equal(lame:::.lame_fast_binary_density(y, "bipartite"), 2 / 6)
	b = matrix(c(1, 0, 0, 0, 1, NA), 2, 3)
	expect_equal(lame:::.lame_fast_binary_density(b, "bipartite"), 2 / 5)
	# a panel with changing actors pools its periods
	expect_equal(lame:::.lame_fast_binary_density(list(y, matrix(c(NA, 1, 1, NA), 2, 2))), 4 / 8)
	expect_equal(lame:::.lame_fast_binary_density(array(c(y, y), c(3, 3, 2))), 2 / 6)
	expect_true(is.na(lame:::.lame_fast_binary_density(matrix("a", 2, 2))))
	expect_true(is.na(lame:::.lame_fast_binary_density(matrix(NA_real_, 2, 2))))
})

test_that("the warning fires in both tails, only for binary networks and only outside the band", {
	expect_equal(lame:::.lame_fast_binary_min_density, 0.20)
	sparse = rule_net(p = 0.08)$y
	dense = 1 - sparse
	mid = rule_net(p = 0.3, seed = 2)$y
	expect_gt(min(lame:::.lame_fast_binary_density(mid), 1 - lame:::.lame_fast_binary_density(mid)), 0.2)
	expect_warning(lame:::.lame_fast_binary_warn(sparse, "binary", "unipartite"), class = "lame_fast_binary_density")
	expect_warning(lame:::.lame_fast_binary_warn(dense, "binary", "unipartite"), class = "lame_fast_binary_density")
	expect_warning(lame:::.lame_fast_binary_warn(sparse, "binary", "unipartite", "lame"), "lame")
	expect_no_warning(lame:::.lame_fast_binary_warn(mid, "binary", "unipartite"))
	expect_no_warning(lame:::.lame_fast_binary_warn(sparse, "normal", "unipartite"))
	expect_no_warning(lame:::.lame_fast_binary_warn(sparse, "poisson", "unipartite"))
	# a density that cannot be computed stays silent
	expect_no_warning(lame:::.lame_fast_binary_warn(matrix("a", 3, 3), "binary", "unipartite"))
	expect_no_warning(lame:::.lame_fast_binary_warn(NULL, "binary", "unipartite"))
})

test_that("each user-facing fast call on a sparse or dense binary network warns once and returns the fast fit", {
	skip_on_cran()
	d = rule_net(n = 20)
	w = count_density_warnings(ame_als(d$y, Xdyad = d$x, R = 1, family = "binary", verbose = FALSE))
	expect_equal(w$n, 1)
	expect_s3_class(w$value, "ame_als")
	expect_true(finite_se(w$value))
	w = count_density_warnings(ame(1 - d$y, Xdyad = d$x, R = 0, family = "binary", method = "als", verbose = FALSE))
	expect_equal(w$n, 1)
	expect_s3_class(w$value, "ame_als")
	expect_true(finite_se(w$value))
	b = rule_net(n = 14, nc = 10, seed = 3)
	w = count_density_warnings(ame_als(b$y, R = 1, family = "binary", mode = "bipartite", verbose = FALSE))
	expect_equal(w$n, 1)
	expect_true(finite_se(w$value))
	ys = lapply(1:3, function(t) rule_net(n = 12, seed = 10 + t)$y)
	w = count_density_warnings(lame_als(ys, R = 1, family = "binary", verbose = FALSE))
	expect_equal(w$n, 1)
	expect_s3_class(w$value, "lame_als")
	expect_true(finite_se(w$value))
	w = count_density_warnings(lame(ys, R = 1, family = "binary", method = "als", verbose = FALSE))
	expect_equal(w$n, 1)
	expect_s3_class(w$value, "lame_als")
	w = count_density_warnings(lame(ys, R = 1, family = "binary", method = "als", dynamic_uv = TRUE, verbose = FALSE))
	expect_equal(w$n, 1)
	expect_s3_class(w$value, "lame_als")
	expect_true(finite_se(w$value))
})

test_that("a binary network inside the band and other families do not warn", {
	skip_on_cran()
	d = rule_net(n = 16, p = 0.4, seed = 4)
	expect_equal(count_density_warnings(ame_als(d$y, R = 0, family = "binary", verbose = FALSE))$n, 0)
	set.seed(5)
	yn = matrix(rnorm(256) - 3, 16); diag(yn) = NA
	yp = matrix(rpois(256, 0.05), 16); diag(yp) = NA
	expect_equal(count_density_warnings(ame_als(yn, R = 0, family = "normal", verbose = FALSE))$n, 0)
	expect_equal(count_density_warnings(ame_als(yp, R = 0, family = "poisson", verbose = FALSE))$n, 0)
})

test_that("bootstrap replicates, refits and update() of an existing fit do not warn", {
	skip_on_cran()
	d = rule_net(n = 16, seed = 6)
	w = count_density_warnings(ame_als(d$y, R = 0, family = "binary", bootstrap = 30, verbose = FALSE))
	expect_equal(w$n, 1)
	fit = w$value
	expect_true(finite_se(fit))
	expect_equal(count_density_warnings(ame_als_bootstrap(fit, R = 30, verbose = FALSE))$n, 0)
	expect_equal(count_density_warnings(ame_als_refit(fit, verbose = FALSE))$n, 0)
	expect_equal(count_density_warnings(update(fit, max_iter = 50))$n, 0)
})

test_that("the sampler does not warn, including the fast warm start of dynamic positions", {
	skip_on_cran()
	# the warm start runs under suppressWarnings(), so count the calls that would warn
	orig = lame:::.lame_fast_binary_warn
	n = 0
	testthat::local_mocked_bindings(
		.lame_fast_binary_warn = function(...) {
			d = suppressWarnings(orig(...))
			if (!is.null(d)) n <<- n + 1
			invisible(d)
		},
		.package = "lame"
	)
	ys = lapply(1:3, function(t) rule_net(n = 10, seed = 20 + t)$y)
	short = list(nscan = 100, burn = 20, odens = 5, chains = 1, verbose = FALSE, gof = FALSE)
	w = count_density_warnings(do.call(lame, c(list(ys, R = 1, family = "binary", dynamic_uv = TRUE), short)))
	expect_s3_class(w$value, "lame")
	expect_true(isTRUE(w$value$als_start))
	d = rule_net(n = 12, seed = 7)
	w = count_density_warnings(do.call(ame, c(list(d$y, R = 0, family = "binary"), short)))
	expect_s3_class(w$value, "ame")
	expect_equal(n, 0)
	# the counter sees a fast call
	count_density_warnings(ame_als(d$y, R = 0, family = "binary", verbose = FALSE))
	expect_equal(n, 1)
})
