# essential tests for cran — runs in under 2 minutes
# this file does NOT use skip_on_cran() so it always runs
# full suite (1900+ tests) runs via devtools::test() locally
# to run only cran tests: devtools::check(env_vars = c(NOT_CRAN = "false"))

# minimal mcmc settings for speed
burn = 5
nscan = 25
odens = 5

# -------------------------------------------------------
# ame: unipartite binary
# -------------------------------------------------------

test_that("ame fits a unipartite binary model", {
	set.seed(6886)
	n = 15
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	expect_s3_class(fit, "ame")
	expect_equal(nrow(fit$BETA), nscan / odens)
	expect_equal(ncol(fit$YPM), n)
	expect_true(all(is.finite(fit$BETA)))
})

# -------------------------------------------------------
# ame: unipartite normal with covariates
# -------------------------------------------------------

test_that("ame fits a normal model with covariates", {
	set.seed(6886)
	n = 15
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	Xd = array(rnorm(n * n), dim = c(n, n, 1))
	dimnames(Xd)[[3]] = "covariate"
	diag(Xd[,,1]) = NA

	fit = ame(Y, Xdyad = Xd, R = 0, family = "normal",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	expect_s3_class(fit, "ame")
	expect_equal(ncol(fit$BETA), 2)  # intercept + covariate
	expect_true(all(is.finite(coef(fit))))
})

# -------------------------------------------------------
# ame: bipartite binary
# -------------------------------------------------------

test_that("ame fits a bipartite binary model", {
	set.seed(6886)
	nA = 12; nB = 8
	Y = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
	rownames(Y) = paste0("R", 1:nA)
	colnames(Y) = paste0("C", 1:nB)

	fit = ame(Y, mode = "bipartite", R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	expect_s3_class(fit, "ame")
	expect_equal(fit$mode, "bipartite")
	expect_equal(nrow(fit$U), nA)
	expect_equal(nrow(fit$V), nB)
	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
})

# -------------------------------------------------------
# lame: longitudinal binary
# -------------------------------------------------------

test_that("lame fits a longitudinal binary model", {
	set.seed(6886)
	n = 10; n_time = 3
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rbinom(n * n, 1, 0.2), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = paste0("n", 1:n)
		Y
	})

	fit = lame(Y_list, R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	expect_s3_class(fit, "lame")
	expect_equal(nrow(fit$BETA), nscan / odens)
	expect_true(is.list(fit$EZ))
	expect_equal(length(fit$EZ), n_time)
})

# -------------------------------------------------------
# lame: dynamic effects
# -------------------------------------------------------

test_that("lame fits with dynamic_uv and dynamic_ab", {
	set.seed(6886)
	n = 10; n_time = 3
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rbinom(n * n, 1, 0.2), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = paste0("n", 1:n)
		Y
	})

	fit = lame(Y_list, R = 2, family = "binary",
		dynamic_uv = TRUE, dynamic_ab = TRUE,
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	expect_s3_class(fit, "lame")
	expect_true(isTRUE(fit$dynamic_uv))
	expect_true(isTRUE(fit$dynamic_ab))
	# u should be 3D for dynamic models
	expect_equal(length(dim(fit$U)), 3)
})

# -------------------------------------------------------
# lame: bipartite longitudinal
# -------------------------------------------------------

test_that("lame fits a bipartite longitudinal model", {
	set.seed(6886)
	nA = 8; nB = 6; n_time = 2
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(Y) = paste0("R", 1:nA)
		colnames(Y) = paste0("C", 1:nB)
		Y
	})

	fit = lame(Y_list, mode = "bipartite", R = 2,
		family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	expect_s3_class(fit, "lame")
	expect_equal(fit$mode, "bipartite")
})

# -------------------------------------------------------
# s3 methods: coef, vcov, confint, summary, predict, fitted, residuals
# -------------------------------------------------------

test_that("s3 methods work on ame objects", {
	set.seed(6886)
	n = 12
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "normal",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	# coef
	b = coef(fit)
	expect_true(is.numeric(b))
	expect_true(length(b) > 0)

	# vcov
	v = vcov(fit)
	expect_true(is.matrix(v))
	expect_equal(nrow(v), length(b))

	# confint
	ci = confint(fit)
	expect_equal(nrow(ci), length(b))
	expect_true(all(ci[, 2] >= ci[, 1]))

	# summary
	s = summary(fit)
	expect_s3_class(s, "summary.ame")

	# predict
	pred = predict(fit, type = "response")
	expect_equal(dim(pred), dim(Y))

	# fitted
	f = fitted(fit)
	expect_equal(dim(f), dim(Y))

	# residuals
	r = residuals(fit)
	expect_equal(dim(r), dim(Y))
})

test_that("s3 methods work on lame objects", {
	set.seed(6886)
	n = 10; n_time = 2
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rnorm(n * n), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = paste0("n", 1:n)
		Y
	})

	fit = lame(Y_list, R = 0, family = "normal",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	b = coef(fit)
	expect_true(is.numeric(b))

	ci = confint(fit)
	expect_equal(nrow(ci), length(b))

	s = summary(fit)
	expect_s3_class(s, "summary.lame")

	pred = predict(fit, type = "response")
	expect_true(is.list(pred))
	expect_equal(length(pred), n_time)
})

# -------------------------------------------------------
# simulate
# -------------------------------------------------------

test_that("simulate works for ame and lame", {
	set.seed(6886)
	n = 10
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	sims = simulate(fit, nsim = 3)
	expect_true(is.list(sims$Y))
	expect_equal(length(sims$Y), 3)
	expect_true(all(sims$Y[[1]] %in% c(0, 1, NA)))
})

# -------------------------------------------------------
# plotting functions (smoke tests)
# -------------------------------------------------------

test_that("plotting functions produce ggplot objects", {
	set.seed(6886)
	n = 10
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "normal",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	p1 = trace_plot(fit)
	expect_true(inherits(p1, "gg") || inherits(p1, "patchwork"))

	p2 = gof_plot(fit)
	expect_true(inherits(p2, "gg") || inherits(p2, "patchwork"))

	p3 = uv_plot(fit)
	expect_true(inherits(p3, "gg"))

	p4 = ab_plot(fit, effect = "sender")
	expect_true(inherits(p4, "gg"))
})

# -------------------------------------------------------
# latent_positions and procrustes_align
# -------------------------------------------------------

test_that("latent_positions returns tidy data frame", {
	set.seed(6886)
	n = 10
	Y = matrix(rnorm(n * n), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "normal",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	lp = latent_positions(fit)
	expect_true(is.data.frame(lp))
	expect_true("actor" %in% names(lp))
	expect_true("dimension" %in% names(lp))
	expect_true("value" %in% names(lp))
	expect_true("type" %in% names(lp))
})

# -------------------------------------------------------
# input validation
# -------------------------------------------------------

test_that("ame rejects non-square Y for unipartite", {
	Y = matrix(1, 5, 3)
	expect_error(ame(Y, family = "normal", burn = 1, nscan = 1, odens = 1))
})

test_that("check_format validates lame inputs", {
	# y must be a list
	expect_error(check_format(Y = matrix(1, 5, 5)))
	# covariate length must match
	Y_list = list(matrix(1, 5, 5), matrix(1, 5, 5))
	expect_error(check_format(Y = Y_list, Xdyad = list(matrix(1, 5, 5))))
})

# -------------------------------------------------------
# gof and gof_stats
# -------------------------------------------------------

test_that("gof_stats returns named vector", {
	set.seed(6886)
	Y = matrix(rbinom(100, 1, 0.3), 10, 10)
	diag(Y) = NA
	stats = gof_stats(Y)
	expect_true(is.numeric(stats))
	expect_true(length(stats) >= 4)
	expect_true("sd.rowmean" %in% names(stats))
})

test_that("gof_stats_bipartite returns named vector", {
	set.seed(6886)
	Y = matrix(rbinom(80, 1, 0.3), 10, 8)
	stats = gof_stats_bipartite(Y)
	expect_true(is.numeric(stats))
	expect_true("sd.rowmean" %in% names(stats))
	expect_true("four.cycles" %in% names(stats))
})

test_that("bipartite ame GOF first row is observed", {
	set.seed(6886)
	nA = 12; nB = 8
	Y = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
	rownames(Y) = paste0("R", 1:nA)
	colnames(Y) = paste0("C", 1:nB)

	fit = ame(Y, mode = "bipartite", R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	obs = gof_stats(Y, mode = "bipartite")
	expect_equal(unname(fit$GOF[1, ]), unname(obs))
})

# -------------------------------------------------------
# utility functions
# -------------------------------------------------------

test_that("list_to_array converts Y list to 3D array", {
	set.seed(6886)
	n = 8
	nms = paste0("n", 1:n)
	Y_list = list(
		matrix(rnorm(n * n), n, n, dimnames = list(nms, nms)),
		matrix(rnorm(n * n), n, n, dimnames = list(nms, nms))
	)
	diag(Y_list[[1]]) = diag(Y_list[[2]]) = NA

	arr = list_to_array(actors = nms, Y = Y_list,
		Xdyad = NULL, Xrow = NULL, Xcol = NULL)
	expect_equal(length(dim(arr$Y)), 3)
	expect_equal(dim(arr$Y)[3], 2)
})

test_that("el2sm and sm2el round-trip", {
	set.seed(6886)
	n = 6
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = 0
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	el = sm2el(Y)
	expect_true(is.data.frame(el) || is.matrix(el))

	Y2 = el2sm(el)
	expect_equal(dim(Y2), dim(Y))
})

# -------------------------------------------------------
# gof with custom_gof
# -------------------------------------------------------

test_that("gof with custom_gof works for ame", {
	set.seed(6886)
	n = 12
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("n", 1:n)

	fit = ame(Y, R = 2, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE)

	custom_fn = function(Y) {
		c(density = mean(Y, na.rm = TRUE))
	}

	result = gof(fit, custom_gof = custom_fn, nsim = 3, verbose = FALSE)
	expect_true(is.matrix(result))
	expect_true("density" %in% colnames(result))
	expect_true(nrow(result) >= 2)
})

test_that("gof with custom_gof works for lame", {
	set.seed(6886)
	n = 10; n_time = 2
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rbinom(n * n, 1, 0.2), n, n)
		diag(Y) = NA
		rownames(Y) = colnames(Y) = paste0("n", 1:n)
		Y
	})

	fit = lame(Y_list, R = 0, family = "binary",
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	custom_fn = function(Y) {
		c(density = mean(Y, na.rm = TRUE))
	}

	result = gof(fit, custom_gof = custom_fn, nsim = 3, verbose = FALSE)
	expect_true(is.matrix(result))
	expect_true("density" %in% colnames(result))
})

# -------------------------------------------------------
# procrustes_align with bipartite dynamic model
# -------------------------------------------------------

test_that("procrustes_align works for bipartite dynamic model", {
	set.seed(6886)
	nA = 8; nB = 6; n_time = 3
	Y_list = lapply(1:n_time, function(t) {
		Y = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(Y) = paste0("R", 1:nA)
		colnames(Y) = paste0("C", 1:nB)
		Y
	})

	fit = lame(Y_list, mode = "bipartite", R = 2,
		family = "binary", dynamic_uv = TRUE,
		burn = burn, nscan = nscan, odens = odens,
		verbose = FALSE, plot = FALSE)

	aligned = procrustes_align(fit)
	expect_true(!is.null(aligned$U))
	expect_true(!is.null(aligned$V))
	expect_true(!is.null(aligned$G))
	expect_equal(length(dim(aligned$U)), 3)
	expect_equal(length(dim(aligned$V)), 3)
	expect_equal(dim(aligned$U)[1], nA)
	expect_equal(dim(aligned$V)[1], nB)
})
