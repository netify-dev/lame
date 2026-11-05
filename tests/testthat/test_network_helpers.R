# direct checks of the exported helpers for covariates, input conversion, held-out scoring, memory
# reporting and binary latent draws; all fast enough for CRAN
skip_on_cran()

test_that("absdiff, nodematch and nodefactor build the dyadic arrays they document", {
	x = c(1, 4, 6)
	expect_equal(absdiff(x), matrix(c(NA, 3, 5, 3, NA, 2, 5, 2, NA), 3))
	expect_equal(diag(absdiff(x, na_diag = FALSE)), rep(0, 3))
	expect_error(absdiff(c("a", "b")), "numeric")
	g = c("a", "b", "a")
	m = nodematch(g)
	expect_identical(storage.mode(m), "double")
	expect_equal(m[1, 3], 1); expect_equal(m[1, 2], 0); expect_true(all(is.na(diag(m))))
	expect_equal(diag(nodematch(g, na_diag = FALSE)), rep(1, 3))
	expect_error(nodematch("a"), "length")
	expect_equal(nodefactor(x), outer(x, x, "+") + diag(NA, 3))
	f = nodefactor(factor(g, levels = c("b", "a")))
	expect_named(f, c("b", "a"))
	# the "a" slice counts how many of a dyad's two actors are in "a"
	expect_equal(f$a[1, 3], 2); expect_equal(f$a[1, 2], 1); expect_equal(f$b[1, 3], 0)
	expect_named(nodefactor(c("y", "x")), c("x", "y"))
	expect_error(nodefactor(list(1, 2)), "numeric, factor, or character")
})

test_that("as_lame_y converts matrices and numeric data frames and refuses edgelists", {
	m = matrix(c(0, 1, 1, 0), 2)
	y = as_lame_y(m)
	expect_identical(storage.mode(y), "double")
	expect_true(all(is.na(diag(y)))); expect_equal(y[1, 2], 1)
	expect_equal(diag(as_lame_y(m, na_diag = FALSE)), c(0, 0))
	r = matrix(1, 2, 3)
	expect_false(anyNA(as_lame_y(r)))
	expect_equal(as_lame_y(as.data.frame(m)), y, ignore_attr = TRUE)
	expect_error(as_lame_y(data.frame(from = "a", to = "b")), "edgelist")
	expect_error(as_lame_y(matrix(c("a", "b"), 1)), "edgelist")
	expect_error(as_lame_y(list(1)), "must be an")
})

test_that("as_lame_y reads igraph objects, including bipartite ones", {
	skip_if_not_installed("igraph")
	g = igraph::make_graph(c(1, 2, 2, 3), directed = TRUE)
	y = as_lame_y(g)
	expect_equal(dim(y), c(3L, 3L)); expect_equal(y[1, 2], 1); expect_equal(y[2, 1], 0)
	b = igraph::make_bipartite_graph(c(FALSE, FALSE, TRUE, TRUE, TRUE), c(1, 3, 2, 4))
	expect_message(yb <- as_lame_y(b), "incidence")
	expect_equal(dim(yb), c(2L, 3L))
})

test_that("ame_options returns its arguments as a list", {
	expect_identical(ame_options(), list(parallel_chains = 1, verbose = TRUE, odens = 25, use_sparse_matrices = FALSE))
	expect_identical(ame_options(parallel_chains = 2, odens = 5)$odens, 5)
})

test_that("evaluate_heldout scores binary, normal and Poisson cells on the mask only", {
	y = matrix(c(1, 0, 1, 0), 2); p = matrix(c(.8, .3, .6, .9), 2)
	mask = matrix(c(TRUE, TRUE, TRUE, FALSE), 2)
	r = suppressMessages(evaluate_heldout(y, p, mask, family = "binary"))
	keep = c(1, 2, 3)
	expect_equal(r$n_eval, 3L)
	expect_equal(r$brier, mean((p[keep] - y[keep])^2))
	expect_equal(r$logloss, -mean(y[keep] * log(p[keep]) + (1 - y[keep]) * log(1 - p[keep])))
	# lists of periods are flattened
	r2 = suppressMessages(evaluate_heldout(list(y, y), list(p, p), list(mask, mask)))
	expect_equal(r2$n_eval, 6L); expect_equal(r2$brier, r$brier)
	yn = c(1, 2, 3, 4); pn = c(1.5, 2, 2, 4)
	rn = evaluate_heldout(yn, pn, rep(TRUE, 4), family = "normal")
	expect_equal(rn$rmse, sqrt(mean((yn - pn)^2))); expect_equal(rn$mae, mean(abs(yn - pn)))
	rp = evaluate_heldout(c(0, 2, 5), c(1, 2, 4), rep(TRUE, 3), family = "poisson")
	expect_equal(rp$mean_logdens, mean(dpois(c(0, 2, 5), c(1, 2, 4), log = TRUE)))
	expect_error(evaluate_heldout(1:3, 1:2, rep(TRUE, 3)), "same shape")
	expect_error(evaluate_heldout(1:3, 1:3, rep(FALSE, 3)), "selects no")
	expect_warning(r4 <- evaluate_heldout(1:3, 1:3, rep(TRUE, 3), family = "ordinal"), "not handled")
	expect_identical(r4$n_eval, 3L)
})

test_that("ame_memory_settings and ame_memory_usage report sizes", {
	s = suppressMessages(ame_memory_settings(200, R = 2))
	expect_equal(s$matrix_mb, 8 * 200^2 / 1024^2)
	expect_equal(s$factor_mb, 8 * 200 * 2 * 2 / 1024^2)
	expect_equal(s$total_mb, 3 * s$matrix_mb + s$factor_mb)
	expect_equal(s$optimized_mb, s$matrix_mb + s$factor_mb)
	expect_equal(
		suppressMessages(ame_memory_settings(50))$optimized_mb,
		suppressMessages(ame_memory_settings(50))$total_mb
	)
	expect_match(paste(capture_messages(ame_memory_settings(1200)), collapse = ""), "n > 1000")
	# the printed optimized figure is the returned one on both sides of the n > 100 rule
	for (n in c(50, 200)) {
		msg = gsub("\\s+", " ", paste(capture_messages(s_n <- ame_memory_settings(n)), collapse = ""))
		expect_match(msg, paste0("Optimized [^:]*: ", round(s_n$optimized_mb, 1), " MB"))
	}
	expect_match(paste(capture_messages(ame_memory_settings(50)), collapse = ""), "all matrices kept")
	fit = structure(list(BETA = matrix(0, 10, 2), YPM = matrix(0, 30, 30), APM = 1:3), class = "ame")
	u = ame_memory_usage(fit, detailed = FALSE)
	expect_setequal(u$Component, c("BETA", "YPM", "APM"))
	expect_identical(u$Component[1], "YPM")
	expect_equal(u$Bytes[u$Component == "YPM"], as.numeric(object.size(fit$YPM)))
	expect_match(paste(capture_messages(ame_memory_usage(fit)), collapse = ""), "Predictions")
	expect_error(ame_memory_usage(), "required")
})

test_that("compact_ame drops empty slots, keeps the class and sparsifies sparse posterior means", {
	fit = structure(
		list(BETA = matrix(0, 2, 1), U = NULL, YPM = diag(5), EZ = list(matrix(1, 3, 3)), U_samples = array(1, c(3, 1, 2))),
		class = c("lame", "list")
	)
	fit["U"] = list(NULL)
	cf = compact_ame(fit)
	expect_false("U" %in% names(cf)); expect_identical(class(cf), class(fit)); expect_identical(cf$YPM, fit$YPM)
	skip_if_not_installed("Matrix")
	cs = suppressMessages(compact_ame(fit, use_sparse_matrices = TRUE))
	expect_s4_class(cs$YPM, "Matrix")
	expect_true(is.matrix(cs$EZ[[1]]))
	expect_identical(cs$U_samples, fit$U_samples)
})

test_that("rZ_bin_fc draws the latent array on the half-line the ties imply", {
	set.seed(31)
	n = 60
	Y = matrix(rbinom(n * n, 1, .5), n); Y[cbind(1:5, 6:10)] = NA
	EZ = matrix(0, n, n); Z = matrix(0, n, n)
	z = rZ_bin_fc(Z, EZ, 0, Y)
	off = row(Y) != col(Y)
	expect_equal(dim(z), dim(Y))
	expect_true(all(z[off & Y %in% 1] > 0)); expect_true(all(z[off & Y %in% 0] < 0))
	expect_true(all(is.finite(z[is.na(Y)])))
	# with EZ = 0 and no correlation a tie's latent value is half-normal with mean sqrt(2 / pi)
	ones = z[off & Y %in% 1]
	expect_lt(abs(mean(ones) - sqrt(2 / pi)) / (sqrt(1 - 2 / pi) / sqrt(length(ones))), 4)
	zr = rZ_bin_fc(z, EZ, .5, Y)
	expect_true(all(zr[off & Y %in% 1] > 0)); expect_true(all(zr[off & Y %in% 0] < 0))
})

test_that("evaluate_heldout's normal density uses the fit's residual variance", {
	y = c(1, 2, 3, 4); p = c(1.5, 2, 2, 4); m = rep(TRUE, 4)
	# without the fit the density is withheld; error summaries are reported
	r0 = evaluate_heldout(y, p, m, family = "normal")
	expect_true(is.na(r0$mean_logdens)); expect_equal(r0$rmse, sqrt(mean((y - p)^2)))
	# an MCMC fit uses the posterior mean of the ve draws (here 1)
	mc = structure(list(family = "normal", VC = cbind(va = c(.1, .2), ve = c(.5, 1.5))), class = "ame")
	expect_equal(
		evaluate_heldout(y, p, m, family = "normal", fit = mc)$mean_logdens,
		mean(dnorm(y, p, 1, log = TRUE)),
		tolerance = 1e-14
	)
	expect_equal(
		evaluate_heldout(y, p, m, family = "normal", fit = mc)$mean_logdens,
		-0.5 * log(2 * pi) - mean((y - p)^2) / 2,
		tolerance = 1e-14
	)
	# an ALS fit uses its fitted ve (here 2)
	als = structure(list(family = "normal", VC = c(va = .1, vb = .2, ve = 2)), class = c("ame_als", "list"))
	expect_equal(
		evaluate_heldout(y, p, m, family = "normal", fit = als)$mean_logdens,
		mean(dnorm(y, p, sqrt(2), log = TRUE)),
		tolerance = 1e-14
	)
	set.seed(3); n = 8
	Y = matrix(rnorm(n * n), n); diag(Y) = NA
	fit = suppressWarnings(suppressMessages(ame_als(Y, family = "normal", R = 0, verbose = FALSE)))
	pred = Y; pred[] = 0
	expect_equal(
		evaluate_heldout(Y, pred, !is.na(Y), family = "normal", fit = fit)$mean_logdens,
		mean(dnorm(Y[!is.na(Y)], 0, sqrt(fit$VC[["ve"]]), log = TRUE)),
		tolerance = 1e-12
	)
	# the fit must be normal with a residual variance
	expect_error(
		evaluate_heldout(
			y,
			p,
			m,
			family = "normal",
			fit = structure(list(family = "binary", VC = c(ve = 1)), class = "ame")
		),
		class = "lame_heldout_fit"
	)
	expect_error(evaluate_heldout(y, p, m, family = "normal", fit = list(VC = c(ve = 1))), class = "lame_heldout_fit")
	expect_error(
		evaluate_heldout(
			y,
			p,
			m,
			family = "normal",
			fit = structure(list(family = "normal", VC = c(va = 1)), class = "ame")
		),
		class = "lame_heldout_fit"
	)
	# other families ignore it
	expect_equal(
		evaluate_heldout(c(0, 2), c(1, 2), c(TRUE, TRUE), family = "poisson", fit = mc)$mean_logdens,
		mean(dpois(c(0, 2), c(1, 2), log = TRUE))
	)
})
