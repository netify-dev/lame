# regression guards for silent-wrong-answer and crash paths in predict,
# fit storage, plotting, and the bipartite samplers.

# predict.ame(newdata=) applies each covariate's own coefficient: feeding the
# training xdyad back reproduces the in-sample linear predictor, and a
# covariate shift moves the predictor by that covariate's own coefficient
# (guards against mixing nodal and dyadic coefficients).
test_that("predict.ame newdata reproduces in-sample EZ with nodal covariates", {
	skip_on_cran()
	set.seed(1); n = 15
	Y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	Xd = array(rnorm(n*n*2), c(n, n, 2), dimnames = list(NULL, NULL, c("d1","d2")))
	Xr = matrix(rnorm(n*2), n, 2, dimnames = list(NULL, c("r1","r2")))
	Xc = matrix(rnorm(n*2), n, 2, dimnames = list(NULL, c("c1","c2")))
	fit = ame(Y, Xdyad = Xd, Xrow = Xr, Xcol = Xc, family = "binary", R = 1,
	           burn = 50, nscan = 300, odens = 5, verbose = FALSE)
	ez_in  = predict(fit, type = "link")
	ez_new = predict(fit, type = "link", newdata = Xd)
	expect_lt(max(abs(ez_in - ez_new), na.rm = TRUE), 1e-8)
	# counterfactual: +2 on d1 must move the link by exactly 2 * beta_d1_dyad
	bm = colMeans(fit$BETA)
	Xd_up = Xd; Xd_up[, , 1] = Xd_up[, , 1] + 2
	ez_up = predict(fit, type = "link", newdata = Xd_up)
	expect_equal(mean(ez_up - ez_new, na.rm = TRUE),
	             unname(2 * bm["d1_dyad"]), tolerance = 1e-6)
})

test_that("predict.ame newdata stays correct for a dyad-only model", {
	skip_on_cran()
	set.seed(2); n = 15
	Y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(Y) = NA
	rownames(Y) = colnames(Y) = paste0("a", seq_len(n))
	Xd = array(rnorm(n*n), c(n, n, 1), dimnames = list(NULL, NULL, "d1"))
	fit = ame(Y, Xdyad = Xd, family = "binary", R = 1,
	           burn = 50, nscan = 300, odens = 5, verbose = FALSE)
	ez_in  = predict(fit, type = "link")
	ez_new = predict(fit, type = "link", newdata = Xd)
	expect_lt(max(abs(ez_in - ez_new), na.rm = TRUE), 1e-8)
})

# dynamic_uv / dynamic_ab fits store the innovation-sd chains fit$sigma_uv /
# fit$sigma_ab, and print.lame() populates the sigma_innov column.
test_that("dynamic_uv / dynamic_ab fits store innovation-SD chains", {
	skip_on_cran()
	set.seed(6886); n = 18; Tn = 4
	Xl = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	fit = lame(Yl, Xdyad = Xl, family = "binary", R = 1,
	            dynamic_uv = TRUE, dynamic_ab = TRUE,
	            nscan = 200, burn = 40, odens = 5, verbose = FALSE)
	expect_false(is.null(fit$sigma_uv))
	expect_false(is.null(fit$sigma_ab))
	expect_equal(length(fit$sigma_uv), 200 / 5)
	expect_true(all(fit$sigma_uv >= 0))
	expect_true(all(fit$sigma_ab >= 0))
	# the documented round-trip works
	expect_true(is.finite(mean(fit$sigma_uv)^2))
	# print.lame() populates the sigma_innov column (no blank)
	out = capture.output(print(fit))
	uv_row = grep("U/V", out, value = TRUE)
	expect_true(length(uv_row) == 1L && grepl("[0-9]", uv_row))
})

# autoplot(fit, which = "uv") on an r = 0 fit gives a clean message rather
# than crashing (u/v are 0-column matrices, not null).
test_that("uv_plot / autoplot give a clean message on an R = 0 fit", {
	skip_on_cran()
	set.seed(1); n = 12; Tn = 3
	Xl = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	fit0 = lame(Yl, Xdyad = Xl, family = "binary", R = 0,
	             dynamic_beta = "dyad", nscan = 120, burn = 30, odens = 5,
	             verbose = FALSE)
	expect_error(autoplot(fit0, which = "uv"), "no multiplicative effects|R = 0")
	expect_error(uv_plot(fit0), "no multiplicative effects|R = 0")
})

test_that("gof_plot ignores missing longitudinal predictive draws", {
	skip_if_not_installed("ggplot2")

	gof_list = list(
		sd.rowmean = matrix(c(1, 2, 3, 1.1, NA, 3.2, 0.9, NaN, 2.8),
		                    nrow = 3, byrow = FALSE),
		sd.colmean = matrix(c(1, 2, 3, 1.2, 1.8, NA, 0.8, 2.1, NaN),
		                    nrow = 3, byrow = FALSE),
		four.cycles = matrix(c(10, 12, 14, 11, NA, 15, 9, 13, NaN),
		                     nrow = 3, byrow = FALSE)
	)
	fit = list(GOF = gof_list, mode = "bipartite")
	class(fit) = c("lame", "ame")

	p = gof_plot(fit)
	expect_s3_class(p, "ggplot")
})

# gof_plot(statistics=) selects panels instead of aborting. two defects are
# guarded: the warning path raised cli's "Multiple quantities for
# pluralization" error because one message interpolated two vectors, and the
# internal gof column names reported by names(fit$GOF) were rejected as
# unknown even though only the display aliases were ever accepted.
test_that("gof_plot statistics= accepts aliases, internal names, and warns on bad", {
	skip_if_not_installed("ggplot2")

	set.seed(2)
	mk = function(nr, nc) matrix(runif(nr * nc, 1, 2), nr, nc)
	gof_list = list(
		sd.rowmean = mk(3, 9), sd.colmean = mk(3, 9),
		dyad.dep   = mk(3, 9), cycle.dep  = mk(3, 9),
		trans.dep  = mk(3, 9)
	)
	fit_l = structure(list(GOF = gof_list, mode = "unipartite"),
	                  class = c("lame", "ame"))

	n_facet = function(p) length(unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL))

	# documented alias selects exactly one panel
	p1 = gof_plot(fit_l, statistics = "sd.row")
	expect_s3_class(p1, "ggplot")
	expect_identical(n_facet(p1), 1L)

	# the internal column name reported by names(fit$GOF) is equally valid
	p2 = expect_no_warning(gof_plot(fit_l, statistics = "sd.rowmean"))
	expect_identical(n_facet(p2), 1L)
	expect_no_warning(gof_plot(fit_l, statistics = "cycle.dep"))

	# multiple statistics, mixing both spellings
	p3 = expect_no_warning(gof_plot(fit_l, statistics = c("sd.rowmean", "trans.dep")))
	expect_identical(n_facet(p3), 2L)

	# an alias and its internal name name the same panel, not two
	expect_identical(n_facet(gof_plot(fit_l, statistics = c("sd.row", "sd.rowmean"))), 1L)

	# a genuinely unknown name warns cleanly and is dropped
	expect_warning(gof_plot(fit_l, statistics = c("sd.row", "not.a.stat")),
	               "Unknown gof statistic")
	expect_s3_class(
		suppressWarnings(gof_plot(fit_l, statistics = "not.a.stat")), "ggplot")

	# same contract on the cross-sectional path
	gof_mat = matrix(runif(9 * 5, 1, 2), 9, 5,
	                 dimnames = list(NULL, names(gof_list)))
	fit_a = structure(list(GOF = gof_mat, mode = "unipartite"), class = "ame")

	expect_s3_class(expect_no_warning(gof_plot(fit_a, statistics = "sd.rowmean")), "ggplot")
	expect_no_warning(gof_plot(fit_a, statistics = c("sd.row", "sd.col")))
	expect_warning(gof_plot(fit_a, statistics = "bogus"), "Unknown gof statistic")
})

# near-separable sparse binary data stays finite and bounded under
# dynamic_beta. (This test previously ALSO asserted the divergence-safeguard
# warning fired; after the additive-effects update was routed through the
# joint rbeta_ab sampler the fit stays bounded on its own -- max |coef| ~2
# instead of hitting the clip ceiling -- so the safeguard correctly stays
# silent. The substantive guarantees are the boundedness checks below.)
test_that("dynamic_beta stays finite and bounded on sparse binary", {
	skip_on_cran()
	set.seed(7); n = 15; Tn = 6
	Xl = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(rbinom(n*n, 1, 0.02), n, n); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	fit = suppressWarnings(
		lame(Yl, Xdyad = Xl, family = "binary", R = 0,
		     dynamic_beta = "dyad", nscan = 300, burn = 80, odens = 5,
		     verbose = FALSE))
	expect_true(all(is.finite(coef(fit))))
	expect_lt(max(abs(coef(fit))), 1e2)   # bounded on its own, no clip needed
})

test_that("dynamic_beta well-identified normal fit is not clipped", {
	skip_on_cran()
	set.seed(1); n = 15; Tn = 5
	bt = seq(-0.5, 1.0, length.out = Tn)
	Xl = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(Tn), function(t) {
		y = bt[t] * Xl[[t]][, , 1] + matrix(rnorm(n*n, 0, 0.4), n, n)
		diag(y) = NA; rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	# no divergence warning on a clean fit; recovers the ramp
	fit = lame(Yl, Xdyad = Xl, family = "normal", R = 0,
	            dynamic_beta = "dyad", nscan = 300, burn = 50, odens = 5,
	            verbose = FALSE)
	cp = coef(fit)["x1_dyad", ]
	expect_lt(max(abs(cp - bt)), 0.4)
})

# predict.lame(newdata=) in-sample (h = 0) includes the nodal (xrow / xcol)
# contribution, so a fit with nodal covariates reproduces its own ez when fed
# its training xdyad back.
#
# r = 0 isolates the nodal term: with no uv term the stored ez is the running
# mean of (x beta + a + b), which equals x beta_pm + apm + bpm to machine
# precision, so a dropped nodal contribution (order ~1) fails loudly while the
# correct path reproduces to < 1e-8. (r >= 1 adds a mean(uv') vs
# mean(u)mean(v)' reconstruction gap unrelated to the nodal term.)
# names are zero-padded (a01..) so they survive lame's internal alphabetical
# actor reordering.
test_that("predict.lame newdata reproduces in-sample EZ with nodal covariates", {
	skip_on_cran()
	set.seed(11); n = 14; Tn = 4
	pad = sprintf("a%02d", seq_len(n))
	Xl = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(pad, pad, "x1"))
	})
	Xr = lapply(seq_len(Tn), function(t)
		matrix(rnorm(n), n, 1, dimnames = list(pad, "r1")))
	Xc = lapply(seq_len(Tn), function(t)
		matrix(rnorm(n), n, 1, dimnames = list(pad, "c1")))
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(y) = NA
		rownames(y) = colnames(y) = pad; y
	})
	fit = lame(Yl, Xdyad = Xl, Xrow = Xr, Xcol = Xc, family = "binary", R = 0,
	            nscan = 200, burn = 40, odens = 5, verbose = FALSE)
	ez_in  = predict(fit, type = "link")
	ez_new = predict(fit, type = "link", newdata = Xl)
	md = max(mapply(function(a, b) max(abs(a - b), na.rm = TRUE),
	                 ez_in, ez_new))
	expect_lt(md, 1e-8)
	# counterfactual: +2 on the dyad covariate moves the link by 2 * its own
	# coefficient at every period -- guards against re-introducing the nodal mix.
	bm = coef(fit)
	bdyad = if (is.matrix(bm)) bm["x1_dyad", ] else rep(bm["x1_dyad"], Tn)
	Xup = lapply(Xl, function(X) { X[, , 1] = X[, , 1] + 2; X })
	ez_up = predict(fit, type = "link", newdata = Xup)
	for (t in seq_len(Tn)) {
		expect_equal(mean(ez_up[[t]] - ez_new[[t]], na.rm = TRUE),
		             unname(2 * bdyad[t]), tolerance = 1e-6,
		             info = paste("period", t))
	}
})

# lame() sorts actors alphabetically on ingest, so a fit with names a1..a14
# stores effects in order a1, a10, a11, ..., a2. predict.lame realigns newdata
# to the fit's order by name, so newdata supplied in the original input order
# still lines up with the stored a / b / uv effects.
test_that("predict.lame newdata realigns actors by name to the fit's order", {
	skip_on_cran()
	set.seed(11); n = 14; Tn = 3
	nm = paste0("a", seq_len(n))                  # a1..a14 sort non-naturally
	stopifnot(!identical(nm, sort(nm)))            # guard: order really differs
	Xl = lapply(seq_len(Tn), function(t)
		array(matrix(rnorm(n*n), n, n), c(n, n, 1), dimnames = list(nm, nm, "x1")))
	Yl = lapply(seq_len(Tn), function(t) {
		y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(y) = NA
		rownames(y) = colnames(y) = nm; y
	})
	fit = lame(Yl, Xdyad = Xl, family = "binary", R = 0,
	            nscan = 150, burn = 40, odens = 5, verbose = FALSE)
	ez_in  = predict(fit, type = "link")
	# newdata supplied in the original (unsorted) actor order
	ez_new = predict(fit, type = "link", newdata = Xl)
	md = max(mapply(function(a, b) max(abs(a - b), na.rm = TRUE), ez_in, ez_new))
	expect_lt(md, 1e-8)
	# the returned matrices are labelled in the fit's (sorted) order
	expect_identical(rownames(ez_new[[1]]), sort(nm))
})
