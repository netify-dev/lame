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
	set.seed(6886); n = 18; T = 4
	Xl = lapply(seq_len(T), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(T), function(t) {
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
	set.seed(1); n = 12; T = 3
	Xl = lapply(seq_len(T), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(T), function(t) {
		y = matrix(rbinom(n*n, 1, 0.3), n, n); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	fit0 = lame(Yl, Xdyad = Xl, family = "binary", R = 0,
	             dynamic_beta = "dyad", nscan = 120, burn = 30, odens = 5,
	             verbose = FALSE)
	expect_error(autoplot(fit0, which = "uv"), "no multiplicative effects|R = 0")
	expect_error(uv_plot(fit0), "no multiplicative effects|R = 0")
})

# the dynamic_beta divergence safeguard keeps coef() finite and warns on
# near-separable sparse binary data, and leaves a well-identified normal fit
# unaffected (the clip never binds).
test_that("dynamic_beta stays finite on sparse binary and warns", {
	skip_on_cran()
	set.seed(7); n = 15; T = 6
	Xl = lapply(seq_len(T), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(T), function(t) {
		y = matrix(rbinom(n*n, 1, 0.02), n, n); diag(y) = NA
		rownames(y) = colnames(y) = paste0("a", seq_len(n)); y
	})
	expect_warning({
		fit = lame(Yl, Xdyad = Xl, family = "binary", R = 0,
		           dynamic_beta = "dyad", nscan = 300, burn = 80, odens = 5,
		           verbose = FALSE)
	},
		"divergence safeguard")
	expect_true(all(is.finite(coef(fit))))
	expect_lt(max(abs(coef(fit))), 1e4)   # bounded, not 1e25
})

test_that("dynamic_beta well-identified normal fit is not clipped", {
	skip_on_cran()
	set.seed(1); n = 15; T = 5
	bt = seq(-0.5, 1.0, length.out = T)
	Xl = lapply(seq_len(T), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	})
	Yl = lapply(seq_len(T), function(t) {
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
	set.seed(11); n = 14; T = 4
	pad = sprintf("a%02d", seq_len(n))
	Xl = lapply(seq_len(T), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(pad, pad, "x1"))
	})
	Xr = lapply(seq_len(T), function(t)
		matrix(rnorm(n), n, 1, dimnames = list(pad, "r1")))
	Xc = lapply(seq_len(T), function(t)
		matrix(rnorm(n), n, 1, dimnames = list(pad, "c1")))
	Yl = lapply(seq_len(T), function(t) {
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
	bdyad = if (is.matrix(bm)) bm["x1_dyad", ] else rep(bm["x1_dyad"], T)
	Xup = lapply(Xl, function(X) { X[, , 1] = X[, , 1] + 2; X })
	ez_up = predict(fit, type = "link", newdata = Xup)
	for (t in seq_len(T)) {
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
	set.seed(11); n = 14; T = 3
	nm = paste0("a", seq_len(n))                  # a1..a14 sort non-naturally
	stopifnot(!identical(nm, sort(nm)))            # guard: order really differs
	Xl = lapply(seq_len(T), function(t)
		array(matrix(rnorm(n*n), n, n), c(n, n, 1), dimnames = list(nm, nm, "x1")))
	Yl = lapply(seq_len(T), function(t) {
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

# the bipartite cbin / frn / rrl z-samplers draw from a truncated normal via
# qnorm(runif(., pnorm(lb), pnorm(ub))). on wide rows (nb >= 50) a finite bound
# can sit far enough in the tail that both pnorm() values saturate to 0/1 and
# qnorm() returns +/-inf. clamping the truncation probabilities into the open
# interval keeps z finite so the fit runs.
test_that("bipartite ame() rrl runs on wide rows (nB >= 50) without NaN crash", {
	skip_on_cran()
	for (sd in c(7, 42)) {
		set.seed(sd); nA = 15; nB = 50
		Y = matrix(0L, nA, nB)
		for (i in seq_len(nA)) {
			k = sample(3:8, 1); top = sample(nB, k)
			Y[i, top] = sample(seq_len(k))      # ranks 1..k, higher = preferred
		}
		rownames(Y) = paste0("r", seq_len(nA)); colnames(Y) = paste0("c", seq_len(nB))
		Xb = array(rnorm(nA*nB), c(nA, nB, 1),
		            dimnames = list(rownames(Y), colnames(Y), "x1"))
		f = expect_error(
			ame(Y, Xdyad = Xb, mode = "bipartite", family = "rrl",
			    R_row = 1, R_col = 1, nscan = 200, burn = 40, odens = 5,
			    seed = sd, verbose = FALSE, gof = FALSE),
			NA)                                   # expect no error
		expect_true(all(is.finite(f$BETA)))
		expect_true(all(is.finite(unlist(f$APM))) && all(is.finite(unlist(f$BPM))))
	}
})
