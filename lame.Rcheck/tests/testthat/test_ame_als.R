skip_on_cran()

# tests for the fast (mcmc-free) ame estimator: ame_als() / lame_als(),
# ame_als_refit(), ame_als_bootstrap(), and their s3 methods.

# --------------------------------------------------------------------------
# helpers: simulate ame data with known parameters
# --------------------------------------------------------------------------

sim_ame = function(n = 25, Tt = 1, R = 0, mu = 0.3, beta = 0.6,
                    s2 = 0.4, mode = "unipartite", nc = NULL,
                    symmetric = FALSE, seed = 1) {
	set.seed(seed)
	nr = n
	if (is.null(nc)) nc = if (mode == "bipartite") n - 7L else n
	a = rnorm(nr, 0, 0.7); a = a - mean(a)
	b = rnorm(nc, 0, 0.7); b = b - mean(b)
	if (symmetric) { b = a; nc = nr }
	U = V = NULL
	if (R > 0) {
		U = matrix(rnorm(nr * R), nr, R)
		V = if (symmetric) U else matrix(rnorm(nc * R), nc, R)
	}
	Ylist = vector("list", Tt); Xlist = vector("list", Tt)
	for (t in seq_len(Tt)) {
		X = matrix(rnorm(nr * nc), nr, nc)
		EZ = mu + beta * X + outer(a, b, "+")
		if (R > 0) EZ = EZ + tcrossprod(U, V)
		if (symmetric) EZ = (EZ + t(EZ)) / 2
		Y = EZ + matrix(rnorm(nr * nc, 0, sqrt(s2)), nr, nc)
		if (symmetric) Y = (Y + t(Y)) / 2
		if (mode != "bipartite") diag(Y) = NA
		Ylist[[t]] = Y; Xlist[[t]] = X
	}
	list(Ylist = Ylist, Xlist = Xlist, a = a, b = b, U = U, V = V,
	     mu = mu, beta = beta, n_row = nr, n_col = nc)
}

# ==========================================================================
# point estimator: structure and basic behaviour
# ==========================================================================

test_that("ame_als fits a cross-sectional normal network", {
	d = sim_ame(n = 25, Tt = 1, R = 1, seed = 11)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 1,
	                   family = "normal", verbose = FALSE)
	expect_s3_class(fit, "ame_als")
	expect_true(fit$converged)
	expect_length(coef(fit), 2L)                  # intercept + 1 dyad covariate
	expect_named(coef(fit), c("intercept", "dyad1_dyad"))
	expect_true(all(is.finite(coef(fit))))
	expect_true(is.matrix(fitted(fit)))
	expect_equal(dim(fitted(fit)), c(25L, 25L))
	expect_equal(fit$n_time, 1L)
	expect_false(fit$meta$uncertainty_available)
})

test_that("lame_als fits a longitudinal normal network", {
	d = sim_ame(n = 20, Tt = 6, R = 1, seed = 12)
	fit = lame_als(d$Ylist, Xdyad = d$Xlist, R = 1,
	                    family = "normal", verbose = FALSE)
	expect_s3_class(fit, "ame_als")
	expect_true(fit$converged)
	expect_equal(fit$n_time, 6L)
	expect_type(fitted(fit), "list")
	expect_length(fitted(fit), 6L)
})

test_that("lame_als accepts a 3D array as input", {
	d = sim_ame(n = 18, Tt = 4, R = 0, seed = 13)
	Yarr = array(unlist(d$Ylist), dim = c(18, 18, 4))
	fit = lame_als(Yarr, R = 0, family = "normal", verbose = FALSE)
	expect_s3_class(fit, "ame_als")
	expect_equal(fit$n_time, 4L)
})

test_that("ame_als handles bipartite (rectangular) networks", {
	# r = 0 fit: check beta and additive effects are actually recovered
	d = sim_ame(n = 30, nc = 20, Tt = 1, R = 0, beta = 0.7,
	             mode = "bipartite", seed = 14)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 0,
	                   family = "normal", mode = "bipartite", verbose = FALSE)
	expect_s3_class(fit, "ame_als")
	expect_equal(fit$mode, "bipartite")
	expect_length(fit$a, 30L)
	expect_length(fit$b, 20L)
	expect_true(is.na(fit$VC["cab"]))             # no cross-covariance bipartite
	expect_equal(unname(coef(fit)["dyad1_dyad"]), 0.7, tolerance = 0.12)
	expect_gt(cor(fit$a, d$a), 0.85)
	expect_gt(cor(fit$b, d$b), 0.85)
	# r > 0 fit: multiplicative factors take the correct rectangular shapes
	d2 = sim_ame(n = 22, nc = 14, Tt = 1, R = 2, mode = "bipartite", seed = 15)
	fit2 = ame_als(d2$Ylist[[1]], R = 2, family = "normal",
	                    mode = "bipartite", verbose = FALSE)
	expect_equal(dim(fit2$U), c(22L, 2L))
	expect_equal(dim(fit2$V), c(14L, 2L))
})

test_that("ame_als handles symmetric networks", {
	d = sim_ame(n = 24, Tt = 1, R = 1, symmetric = TRUE, seed = 15)
	fit = ame_als(d$Ylist[[1]], R = 1, family = "normal",
	                   symmetric = TRUE, verbose = FALSE)
	expect_true(fit$symmetric)
	expect_equal(fit$a, fit$b, tolerance = 1e-8)  # sender == receiver effects
	expect_false(is.null(fit$L))                  # eigenmodel: l is reported
})

test_that("ame_als runs for binary and poisson families", {
	# als supports normal / binary / poisson; ordinal inference
	# goes through the mcmc path.
	set.seed(16)
	Yb = matrix(rbinom(625, 1, 0.3), 25, 25); diag(Yb) = NA
	Yp = matrix(rpois(625, 2), 25, 25); diag(Yp) = NA
	for (fam in c("binary", "poisson")) {
		Y = switch(fam, binary = Yb, poisson = Yp)
		fit = ame_als(Y, R = 1, family = fam, verbose = FALSE)
		expect_s3_class(fit, "ame_als")
		expect_true(all(is.finite(coef(fit))))
		# both default to irls
		expect_match(fit$meta$approximation_note, "reweighted least squares")
	}
})

test_that("non-normal families recover the coefficient sign", {
	set.seed(41)
	n = 40
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	# binary: a strong positive dyadic effect must come back positive
	EZb = -0.3 + 1.0 * Xd + outer(a, b, "+")
	Yb = 1 * (EZb + matrix(rnorm(n * n), n, n) > 0); diag(Yb) = NA
	fb = ame_als(Yb, Xdyad = Xd, R = 0, family = "binary", verbose = FALSE)
	expect_gt(coef(fb)["dyad1_dyad"], 0)
	# poisson: a positive dyadic effect must come back positive
	EZp = 0.4 + 0.3 * Xd + outer(a, b, "+")
	Yp = matrix(rpois(n * n, exp(EZp)), n, n); diag(Yp) = NA
	fp = ame_als(Yp, Xdyad = Xd, R = 0, family = "poisson", verbose = FALSE)
	expect_gt(coef(fp)["dyad1_dyad"], 0)
})

# ==========================================================================
# input validation / clear errors
# ==========================================================================

test_that("unsupported families raise informative errors", {
	# ordinal inference goes through mcmc; cbin / frn have no
	# als likelihood used by this estimator.
	d = sim_ame(n = 15, Tt = 1, seed = 17)
	for (fam in c("ordinal", "cbin", "frn")) {
		expect_error(
			ame_als(d$Ylist[[1]], family = fam, verbose = FALSE),
			"not supported")
	}
})

test_that("mode and dimension mismatches are caught", {
	set.seed(18)
	Yrect = matrix(rnorm(20 * 14), 20, 14)
	expect_error(
		ame_als(Yrect, mode = "unipartite", verbose = FALSE),
		"rectangular")
	d = sim_ame(n = 12, Tt = 1, seed = 19)
	expect_error(
		ame_als(d$Ylist[[1]], R = 20, family = "normal", verbose = FALSE),
		"must be smaller")
})

test_that("malformed max_iter, tol and R inputs raise clean errors", {
	d = sim_ame(n = 14, Tt = 1, seed = 90)
	Y = d$Ylist[[1]]
	expect_error(ame_als(Y, max_iter = 0, verbose = FALSE), "max_iter")
	expect_error(ame_als(Y, max_iter = -5, verbose = FALSE), "max_iter")
	expect_error(ame_als(Y, tol = c(1e-6, 1e-6), verbose = FALSE), "tol")
	expect_error(ame_als(Y, tol = 0, verbose = FALSE), "tol")
	expect_error(ame_als(Y, R = c(1, 2), verbose = FALSE),
	             "single non-negative")
})

# ==========================================================================
# parameter recovery
# ==========================================================================

test_that("ame_als recovers known beta and additive effects", {
	# r = 0 keeps the intercept identifiable (no multiplicative confounding)
	d = sim_ame(n = 30, Tt = 8, R = 0, mu = 0.4, beta = 0.6,
	             s2 = 0.5, seed = 21)
	fit = lame_als(d$Ylist, Xdyad = d$Xlist, R = 0,
	                    family = "normal", verbose = FALSE)
	# point estimates close to truth
	expect_equal(unname(fit$mu), 0.4, tolerance = 0.12)
	expect_equal(unname(fit$beta), 0.6, tolerance = 0.08)
	# additive effects recovered up to noise
	expect_gt(cor(fit$a, d$a), 0.9)
	expect_gt(cor(fit$b, d$b), 0.9)
	# true beta lies inside the bootstrap ci
	bt = ame_als_bootstrap(fit, R = 150, type = "block",
	                            seed = 5, verbose = FALSE)
	ci = confint(bt)
	expect_lte(ci["dyad1_dyad", 1], 0.6)
	expect_gte(ci["dyad1_dyad", 2], 0.6)
})

test_that("inconsistent covariates and actor names are rejected", {
	set.seed(42)
	Y = lapply(1:3, function(t) { m = matrix(rnorm(144), 12, 12)
	                               rownames(m) = colnames(m) = paste0("a", 1:12); m })
	# xdyad slices with differing covariate counts
	Xbad = list(array(rnorm(144), c(12, 12, 1)),
	             array(rnorm(288), c(12, 12, 2)),
	             array(rnorm(144), c(12, 12, 1)))
	expect_error(lame_als(Y, Xdyad = Xbad, family = "normal", verbose = FALSE),
	             "same number of covariates")
	Xname_bad = list(
		array(rnorm(288), c(12, 12, 2),
		      dimnames = list(paste0("a", 1:12), paste0("a", 1:12),
		                      c("x", "z"))),
		array(rnorm(288), c(12, 12, 2),
		      dimnames = list(paste0("a", 1:12), paste0("a", 1:12),
		                      c("z", "x"))),
		array(rnorm(288), c(12, 12, 2),
		      dimnames = list(paste0("a", 1:12), paste0("a", 1:12),
		                      c("x", "z"))))
	expect_error(lame_als(Y, Xdyad = Xname_bad, family = "normal",
	                      verbose = FALSE),
	             "same covariate names")
	# named slices with changing actor sets are aligned to the union panel
	Ybad = Y
	rownames(Ybad[[2]]) = colnames(Ybad[[2]]) = paste0("z", 1:12)
	fit = suppressWarnings(lame_als(Ybad, family = "normal", verbose = FALSE))
	expect_s3_class(fit, "lame_als")
	expect_true(isTRUE(fit$changing_composition))
	expect_equal(dim(fit$Y), c(24L, 24L, 3L))

	Yunnamed = list(matrix(rnorm(9), 3, 3), matrix(rnorm(16), 4, 4))
	expect_error(lame_als(Yunnamed, family = "normal", verbose = FALSE),
	             "requires row and column names")
})

test_that("named unipartite slices are aligned within slice before fitting", {
	actors = paste0("p", 1:4)
	Y1 = matrix(seq_len(16), 4, 4, dimnames = list(actors, actors))
	Y2 = matrix(100 + seq_len(16), 4, 4,
	            dimnames = list(actors, rev(actors)))
	diag(Y1) = NA_real_
	Y2[cbind(seq_len(4), match(actors, rev(actors)))] = NA_real_

	fit = suppressWarnings(lame_als(list(t1 = Y1, t2 = Y2),
	                                R = 0, family = "normal",
	                                verbose = FALSE))

	expect_s3_class(fit, "lame_als")
	expect_false(isTRUE(fit$changing_composition))
	expect_true(all(fit$actor_presence))
	expect_equal(dimnames(fit$Y)[[1]], actors)
	expect_equal(dimnames(fit$Y)[[2]], actors)
	expect_equal(fit$Y["p1", "p4", "t2"], Y2["p1", "p4"])
	expect_equal(fit$Y["p4", "p1", "t2"], Y2["p4", "p1"])

	Ybad = Y1
	colnames(Ybad) = paste0("q", 1:4)
	expect_error(lame_als(list(Y1, Ybad), R = 0, family = "normal",
	                      verbose = FALSE),
	             "row and column")
})

test_that("fixed-layout named covariates are aligned to Y names", {
	rows = paste0("r", 1:3)
	cols = paste0("c", 1:4)
	Y = matrix(rnorm(12), 3, 4, dimnames = list(rows, cols))
	Ylist = list(t1 = Y, t2 = Y + 0.1)
	Xraw = array(seq_len(12), c(3, 4, 1),
	             dimnames = list(rev(rows), rev(cols), "trade"))
	Xlist = list(t1 = Xraw, t2 = Xraw + 100)
	Wr = matrix(seq_along(rows), 3, 1,
	            dimnames = list(rev(rows), "sender"))
	Wc = matrix(seq_along(cols), 4, 1,
	            dimnames = list(rev(cols), "receiver"))

	fit = suppressWarnings(lame_als(
		Ylist, Xdyad = Xlist, Xrow = list(Wr, Wr), Xcol = list(Wc, Wc),
		mode = "bipartite", R = 0, family = "normal", verbose = FALSE))

	expect_equal(fit$X["r1", "c1", "trade_dyad", "t1"],
	             Xraw["r1", "c1", "trade"])
	expect_equal(fit$X["r3", "c4", "trade_dyad", "t2"],
	             (Xraw + 100)["r3", "c4", "trade"])
	expect_equal(fit$W_row["r1", "sender_row"], Wr["r1", "sender"])
	expect_equal(fit$W_col["c1", "receiver_col"], Wc["c1", "receiver"])

	Ydup = Y
	rownames(Ydup) = c("r1", "r1", "r3")
	expect_error(
		lame_als(list(t1 = Ydup, t2 = Ydup), mode = "bipartite",
		         R = 0, family = "normal", verbose = FALSE),
		"invalid row names")
})

# ==========================================================================
# bootstrap
# ==========================================================================

test_that("parametric bootstrap returns a well-formed boot_ame object", {
	d = sim_ame(n = 20, Tt = 1, R = 1, seed = 22)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 1,
	                   family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                            seed = 3, verbose = FALSE)
	expect_s3_class(bt, "boot_ame")
	expect_equal(bt$n_total, 40L)
	expect_lte(bt$n_valid, bt$n_total)
	expect_gt(bt$n_valid, 0L)
	expect_true(all(is.finite(bt$se)))
	expect_true(all(bt$se >= 0))
	expect_true(all(bt$ci_lo <= bt$ci_hi))
	# expanded confint to include vc/a/b/u/v rows; restrict to beta
	# for the old beta-only assertion.
	ci = confint(bt, which = "beta")
	expect_equal(nrow(ci), length(coef(fit)))
	expect_equal(ncol(ci), 2L)
	# default confint(bt) now also covers vc + actor-level uncertainty
	ci_all = confint(bt)
	expect_gt(nrow(ci_all), nrow(ci))
})

test_that("block bootstrap works for longitudinal and errors cross-sectional", {
	dL = sim_ame(n = 18, Tt = 7, R = 1, seed = 23)
	fitL = lame_als(dL$Ylist, R = 1, family = "normal", verbose = FALSE)
	btL = ame_als_bootstrap(fitL, R = 30, type = "block",
	                             seed = 4, verbose = FALSE)
	expect_s3_class(btL, "boot_ame")
	expect_equal(btL$type, "block")
	expect_gt(btL$n_valid, 10L)

	dC = sim_ame(n = 18, Tt = 1, R = 1, seed = 24)
	fitC = ame_als(dC$Ylist[[1]], R = 1, family = "normal", verbose = FALSE)
	expect_error(
		ame_als_bootstrap(fitC, R = 10, type = "block", verbose = FALSE),
		"Block bootstrap")
})

test_that("block and parametric bootstrap give comparable SEs", {
	d = sim_ame(n = 20, Tt = 8, R = 0, beta = 0.6, seed = 25)
	fit = lame_als(d$Ylist, Xdyad = d$Xlist, R = 0,
	                    family = "normal", verbose = FALSE)
	b_blk = ame_als_bootstrap(fit, R = 80, type = "block",
	                               seed = 6, verbose = FALSE)
	b_par = ame_als_bootstrap(fit, R = 80, type = "parametric",
	                               seed = 6, verbose = FALSE)
	ratio = b_blk$se["dyad1_dyad"] / b_par$se["dyad1_dyad"]
	expect_gt(ratio, 0.2)
	expect_lt(ratio, 5)
})

test_that("bootstrap accounts for replicates and warm-starts each refit", {
	d = sim_ame(n = 18, Tt = 6, R = 1, seed = 26)
	fit = lame_als(d$Ylist, R = 1, family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 50, type = "parametric",
	                            seed = 7, verbose = FALSE)
	# n_valid / n_total are reported and consistent
	expect_true(bt$n_valid <= bt$n_total)
	expect_equal(nrow(bt$coefs), bt$n_total)
	# warm start: a refit from the point estimate on the same data is stable
	refit = ame_als_refit(fit, verbose = FALSE)
	expect_s3_class(refit, "ame_als")
	# warm start lands at essentially the same optimum (small numerical drift ok)
	expect_equal(coef(refit), coef(fit), tolerance = 0.02)
	expect_lte(refit$iterations, 30L)
})

test_that("symmetric models bootstrap correctly (block and parametric)", {
	# symmetric longitudinal network
	set.seed(30)
	n = 20; Tt = 6
	a = rnorm(n, 0, 0.7); a = a - mean(a)
	Yl = lapply(seq_len(Tt), function(t) {
		X = matrix(rnorm(n * n), n, n)
		E = 0.3 + 0.5 * X + outer(a, a, "+")
		Y = E + matrix(rnorm(n * n, 0, sqrt(0.5)), n, n)
		Y = (Y + t(Y)) / 2; diag(Y) = NA
		list(Y = Y, X = (X + t(X)) / 2)
	})
	Y = lapply(Yl, `[[`, "Y")
	Xd = lapply(Yl, `[[`, "X")
	fit = lame_als(Y, Xdyad = Xd, R = 0, family = "normal",
	                    symmetric = TRUE, verbose = FALSE)
	expect_true(fit$symmetric)

	b_par = ame_als_bootstrap(fit, R = 60, type = "parametric",
	                               seed = 1, verbose = FALSE)
	b_blk = ame_als_bootstrap(fit, R = 60, type = "block",
	                               seed = 1, verbose = FALSE)
	expect_true(all(is.finite(b_par$se)) && all(b_par$se > 0))
	expect_true(all(is.finite(b_blk$se)) && all(b_blk$se > 0))
	# the parametric ses must not be deflated relative to the block ses:
	# with the symmetric-mirror fix they are the same order of magnitude
	ratio = b_par$se["dyad1_dyad"] / b_blk$se["dyad1_dyad"]
	expect_gt(ratio, 0.4)
	expect_lt(ratio, 2.5)
})

# ==========================================================================
# identifiability: procrustes alignment of u/v
# ==========================================================================

test_that("Procrustes alignment reduces U/V bootstrap replicate spread", {
	d = sim_ame(n = 24, Tt = 7, R = 2, seed = 27)
	fit = lame_als(d$Ylist, R = 2, family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 100, type = "block",
	                            seed = 8, verbose = FALSE)
	# aligned standard errors must be markedly smaller than the raw ones
	expect_lt(mean(bt$se_U), mean(bt$se_U_raw))
	expect_lt(mean(bt$se_V), mean(bt$se_V_raw))
	expect_true(all(is.finite(bt$se_U)))
	expect_true(all(is.finite(bt$se_V)))
})

# ==========================================================================
# identifiability: additive / multiplicative gauge (double-centered o)
# ==========================================================================

test_that("the multiplicative term O is double-centered for R > 0", {
	d = sim_ame(n = 26, Tt = 5, R = 2, seed = 70)
	fit = lame_als(d$Ylist, R = 2, family = "normal", verbose = FALSE)
	O = tcrossprod(fit$U, fit$V)
	# the gauge: o carries no broadcast (row/column-constant) structure
	expect_lt(max(abs(rowMeans(O))), 1e-7)
	expect_lt(max(abs(colMeans(O))), 1e-7)
})

test_that("the gauge is fit-preserving: EZ = mu + Xb + a + b + UV'", {
	d = sim_ame(n = 24, Tt = 4, R = 2, beta = 0.5, seed = 71)
	fit = lame_als(d$Ylist, Xdyad = d$Xlist, R = 2,
	                    family = "normal", verbose = FALSE)
	O = tcrossprod(fit$U, fit$V)
	off = !diag(TRUE, fit$dims$n_row)
	for (t in seq_len(fit$n_time)) {
		recon = fit$mu + unname(fit$beta["dyad1_dyad"]) * d$Xlist[[t]] +
			outer(fit$a, fit$b, "+") + O
		expect_equal(fit$EZ[[t]][off], recon[off], tolerance = 1e-7)
	}
})

test_that("symmetric O is double-centered and reported as U L U'", {
	d = sim_ame(n = 22, Tt = 1, R = 2, symmetric = TRUE, seed = 72)
	fit = ame_als(d$Ylist[[1]], R = 2, family = "normal",
	                   symmetric = TRUE, verbose = FALSE)
	O = fit$U %*% (fit$L * t(fit$U))
	expect_lt(max(abs(rowMeans(O))), 1e-7)
	expect_lt(max(abs(colMeans(O))), 1e-7)
	expect_equal(fit$a, fit$b, tolerance = 1e-8)
})

test_that("the gauge makes the additive effect recover the full sender effect", {
	# under the double-centered gauge the broadcast part of uv' belongs to the
	# additive effect, so fit$a estimates a + rowmeans(uv') (re-centred)
	d = sim_ame(n = 70, Tt = 6, R = 1, s2 = 0.25, seed = 73)
	fit = lame_als(d$Ylist, R = 1, family = "normal", verbose = FALSE)
	bc = rowMeans(tcrossprod(d$U, d$V))
	a_true_gauged = d$a + bc - mean(bc)
	expect_gt(cor(fit$a, a_true_gauged), 0.9)
})

test_that("the gauged fit is reproduced by a warm-start refit", {
	# the double-centered gauge is deterministic, so a warm-start refit must
	# return the same a, b and multiplicative term
	d = sim_ame(n = 28, Tt = 5, R = 2, seed = 74)
	f1 = lame_als(d$Ylist, R = 2, family = "normal", verbose = FALSE)
	f2 = ame_als_refit(f1, verbose = FALSE)
	expect_equal(f1$a, f2$a, tolerance = 1e-3)
	expect_equal(f1$b, f2$b, tolerance = 1e-3)
	expect_equal(tcrossprod(f1$U, f1$V), tcrossprod(f2$U, f2$V),
	             tolerance = 1e-3)
})

test_that("the bootstrap handles a fit with node covariates", {
	set.seed(80)
	n = 22; Tt = 6
	w_row = rnorm(n); w_col = rnorm(n)
	a = rnorm(n, 0, 0.5); a = a - mean(a)
	b = rnorm(n, 0, 0.5); b = b - mean(b)
	parts = lapply(seq_len(Tt), function(t) {
		Xd = matrix(rnorm(n * n), n, n)
		EZ = 0.2 + 0.5 * Xd + matrix(w_row, n, n) +
		      t(matrix(w_col, n, n)) + outer(a, b, "+")
		m = EZ + matrix(rnorm(n * n, 0, 0.6), n, n)
		diag(m) = NA
		list(Y = m, X = Xd)
	})
	Yl = lapply(parts, `[[`, "Y")
	Xd = lapply(parts, `[[`, "X")
	Xrow_l = replicate(Tt, matrix(w_row, n, 1), simplify = FALSE)
	Xcol_l = replicate(Tt, matrix(w_col, n, 1), simplify = FALSE)
	fit = lame_als(Yl, Xdyad = Xd, Xrow = Xrow_l, Xcol = Xcol_l,
	                    R = 0, family = "normal", verbose = FALSE)
	expect_named(coef(fit), c("intercept", "dyad1_dyad", "row1_row", "col1_col"))
	# bootstrap must produce ses/cis for the node coefficients too
	bt = ame_als_bootstrap(fit, R = 40, type = "block",
	                            seed = 2, verbose = FALSE)
	expect_length(bt$se, 4L)
	expect_true(all(is.finite(bt$se)))
	expect_true(all(bt$ci_lo <= bt$ci_hi))
	# default confint() returns all uncertainty channels; restrict
	# to the regression block for the beta-name assertion.
	ci = confint(bt, which = "beta")
	expect_equal(rownames(ci), c("intercept", "dyad1_dyad", "row1_row", "col1_col"))
	# parametric bootstrap also works with node covariates
	btp = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                             seed = 3, verbose = FALSE)
	expect_length(btp$se, 4L)
	expect_true(all(is.finite(btp$se)))
})

# ==========================================================================
# singular / rank-deficient designs: ginv fallback
# ==========================================================================

test_that("rank-deficient designs fall back to a pseudoinverse, not a crash", {
	set.seed(28)
	Y = matrix(rnorm(625), 25, 25); diag(Y) = NA
	# two identical dyadic covariates make x'x exactly singular
	# (the rank-deficiency warning is expected and tested separately)
	x1 = matrix(rnorm(625), 25, 25)
	Xdup = array(c(x1, x1), dim = c(25, 25, 2))
	fit = suppressWarnings(ame_als(Y, Xdyad = Xdup, R = 0,
	                        family = "normal", verbose = FALSE))
	expect_s3_class(fit, "ame_als")
	expect_true(all(is.finite(coef(fit))))        # ginv kept it finite

	# sparse binary network with a constant covariate
	Yb = matrix(rbinom(625, 1, 0.05), 25, 25); diag(Yb) = NA
	fitb = suppressWarnings(ame_als(Yb, Xdyad = matrix(1, 25, 25),
	                         R = 1, family = "binary", verbose = FALSE))
	expect_true(all(is.finite(coef(fitb))))
})

# ==========================================================================
# regression tests for the review-driven fixes
# ==========================================================================

test_that("weighted Block 4 keeps the descent monotone on an unbalanced panel", {
	# dyads observed at unequal numbers of time points -> weighted low-rank
	set.seed(50)
	n = 22; Tt = 8
	U0 = matrix(rnorm(n * 2), n, 2); V0 = matrix(rnorm(n * 2), n, 2)
	a = rnorm(n, 0, 0.6); b = rnorm(n, 0, 0.6)
	Yl = lapply(seq_len(Tt), function(t) {
		Y = 0.3 + outer(a, b, "+") + tcrossprod(U0, V0) +
			matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		miss = sample(which(!is.na(Y)), round(0.25 * n * n))
		Y[miss] = NA                       # unbalanced missingness
		Y
	})
	fit = lame_als(Yl, R = 2, family = "normal", verbose = FALSE)
	dh = fit$dev_history
	# residual sse must be (numerically) monotone non-increasing
	expect_true(all(diff(dh) <= 1e-6))
	expect_true(fit$converged)
})

test_that("symmetric input supplied as a single triangle is handled", {
	set.seed(51)
	n = 18
	g = rnorm(n, 0, 0.7)
	Yfull = outer(g, g, "+") + matrix(rnorm(n * n, 0, 0.5), n, n)
	Yfull = (Yfull + t(Yfull)) / 2; diag(Yfull) = NA
	# supply only the upper triangle; lower triangle missing
	Ytri = Yfull; Ytri[lower.tri(Ytri)] = NA
	fit_tri  = ame_als(Ytri,  R = 0, family = "normal",
	                        symmetric = TRUE, verbose = FALSE)
	fit_full = ame_als(Yfull, R = 0, family = "normal",
	                        symmetric = TRUE, verbose = FALSE)
	expect_true(all(is.finite(coef(fit_tri))))
	# one-triangle and full input give the same fit (pairwise symmetrisation)
	expect_equal(unname(coef(fit_tri)), unname(coef(fit_full)), tolerance = 1e-6)
})

test_that("a network with no observed cells raises an error", {
	Yempty = matrix(NA_real_, 12, 12)
	# (an empty network is also trivially disconnected -> that warning is
	# expected; the error is the behaviour under test)
	suppressWarnings(
		expect_error(ame_als(Yempty, family = "normal", verbose = FALSE),
		             "No observed cells"))
})

test_that("EZ and fitted values have NA diagonals for unipartite networks", {
	d = sim_ame(n = 20, Tt = 3, R = 1, seed = 52)
	fit = lame_als(d$Ylist, R = 1, family = "normal", verbose = FALSE)
	for (t in seq_len(fit$n_time)) {
		expect_true(all(is.na(diag(fit$EZ[[t]]))))
		expect_true(all(is.na(diag(fit$fitted[[t]]))))
	}
})

test_that("node covariates satisfy the orthogonality identification constraint", {
	set.seed(53)
	n = 40
	w_row = rnorm(n); w_col = rnorm(n)
	a = rnorm(n, 0, 0.5); a = a - mean(a)
	b = rnorm(n, 0, 0.5); b = b - mean(b)
	Xd = matrix(rnorm(n * n), n, n)
	EZ = 0.2 + 0.6 * Xd + matrix(w_row, n, n) * 1.5 +
	      t(matrix(w_col, n, n)) * (-1.0) + outer(a, b, "+")
	Y = EZ + matrix(rnorm(n * n, 0, 0.6), n, n); diag(Y) = NA
	fit = ame_als(Y, Xdyad = Xd, Xrow = matrix(w_row, n, 1),
	                   Xcol = matrix(w_col, n, 1), R = 0,
	                   family = "normal", verbose = FALSE)
	# the additive effects are exactly orthogonal to the node covariates
	expect_lt(max(abs(crossprod(fit$W_row, fit$a))), 1e-8)
	expect_lt(max(abs(crossprod(fit$W_col, fit$b))), 1e-8)
	expect_lt(abs(sum(fit$a)), 1e-8)
	expect_named(coef(fit), c("intercept", "dyad1_dyad", "row1_row", "col1_col"))
	# the dyadic coefficient and the between-actor node coefficients are
	# recovered when the additive effects are independent of the covariates
	expect_equal(unname(coef(fit)["dyad1_dyad"]),     0.6, tolerance = 0.1)
	expect_equal(unname(coef(fit)["row1_row"]),  1.5, tolerance = 0.25)
	expect_equal(unname(coef(fit)["col1_col"]), -1.0, tolerance = 0.25)
})

test_that("ame_als_refit reports non-convergence; bootstrap accounts for it", {
	d = sim_ame(n = 18, Tt = 5, R = 1, seed = 54)
	fit = lame_als(d$Ylist, R = 1, family = "normal", verbose = FALSE)
	# a one-iteration refit cannot converge -> the flag must say so
	refit1 = ame_als_refit(fit, max_iter = 1, verbose = FALSE)
	expect_false(refit1$converged)
	# the bootstrap reports the non-converged count and it is consistent
	bt = ame_als_bootstrap(fit, R = 30, type = "block",
	                            seed = 1, verbose = FALSE)
	expect_true(!is.null(bt$n_nonconverged))
	expect_gte(bt$n_nonconverged, 0L)
	expect_lte(bt$n_nonconverged, bt$n_valid)
})

test_that("symmetric R>0 bootstrap reports eigenvalue (L) uncertainty", {
	d = sim_ame(n = 22, Tt = 1, R = 1, symmetric = TRUE, seed = 55)
	fit = ame_als(d$Ylist[[1]], R = 1, family = "normal",
	                   symmetric = TRUE, verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                            seed = 1, verbose = FALSE)
	expect_false(is.null(bt$se_L))
	expect_length(bt$se_L, fit$R)
	expect_true(all(is.finite(bt$se_L)))
})

# ==========================================================================
# regression tests for the tier-2/3 hardening
# ==========================================================================

test_that("a constant covariate triggers a rank-deficiency warning", {
	set.seed(60)
	Y = matrix(rnorm(625), 25, 25); diag(Y) = NA
	expect_warning({
		fit = ame_als(Y, Xdyad = matrix(1, 25, 25), R = 0,
		              family = "normal", verbose = FALSE)
	},
		"rank-deficient")
	expect_true(all(is.finite(coef(fit))))      # ginv kept it finite
})

test_that("off-diagonal covariate NAs trigger a warning and become 0", {
	set.seed(61)
	Y = matrix(rnorm(625), 25, 25); diag(Y) = NA
	X = matrix(rnorm(625), 25, 25); X[4, 9] = NA
	expect_warning(
		ame_als(Y, Xdyad = X, R = 0, family = "normal", verbose = FALSE),
		"replaced with 0")
})

test_that("a disconnected observed-dyad graph triggers a warning", {
	Yd = matrix(NA_real_, 20, 20)
	set.seed(62)
	Yd[1:10, 1:10]   = rnorm(100)
	Yd[11:20, 11:20] = rnorm(100)
	diag(Yd) = NA
	expect_warning(
		ame_als(Yd, R = 0, family = "normal", verbose = FALSE),
		"disconnected components")
})

test_that("the deterministic estimator does not mutate the global RNG", {
	d = sim_ame(n = 18, Tt = 1, R = 1, seed = 63)
	set.seed(7); x1 = runif(1)
	invisible(ame_als(d$Ylist[[1]], R = 1, family = "normal",
	                      verbose = FALSE, seed = 999))
	x2 = runif(1)
	set.seed(7); y1 = runif(1); y2 = runif(1)
	expect_equal(c(x1, x2), c(y1, y2))
})

test_that("symmetric mode symmetrizes covariates so EZ is symmetric", {
	set.seed(64)
	n = 22
	Ys = matrix(rnorm(n * n), n, n); Ys = (Ys + t(Ys)) / 2; diag(Ys) = NA
	X_asym = matrix(rnorm(n * n), n, n)        # deliberately asymmetric
	fit = ame_als(Ys, Xdyad = X_asym, R = 0, family = "normal",
	                   symmetric = TRUE, verbose = FALSE)
	ez = fit$EZ[[1]]
	expect_lt(max(abs(ez - t(ez)), na.rm = TRUE), 1e-8)
})

test_that("confint supports basic and percentile interval types", {
	d = sim_ame(n = 20, Tt = 1, R = 1, seed = 65)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 1,
	                   family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                            seed = 1, verbose = FALSE)
	ci_p = confint(bt, ci_type = "percentile")
	ci_b = confint(bt, ci_type = "basic")
	expect_equal(dim(ci_p), dim(ci_b))
	expect_true(all(ci_b[, 1] <= ci_b[, 2]))
})

test_that("a time-varying node covariate triggers a within-variation warning", {
	set.seed(75)
	n = 16; Tt = 4
	Y = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
		# xrow genuinely varies across time slices
	Xrow_tv = lapply(seq_len(Tt), function(t) matrix(rnorm(n), n, 1))
	expect_warning(
		lame_als(Y, Xrow = Xrow_tv, R = 0, family = "normal",
		             verbose = FALSE),
		"varies within actor")
		# a node covariate constant across time should not warn
	w = matrix(rnorm(n), n, 1)
	Xrow_const = replicate(Tt, w, simplify = FALSE)
	expect_no_warning(
		lame_als(Y, Xrow = Xrow_const, R = 0, family = "normal",
		             verbose = FALSE))
})

test_that("an all-missing time slice triggers a warning", {
	set.seed(76)
	n = 14; Tt = 3
	Y = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
	Y[[2]][] = NA                              # slice 2 entirely missing
	expect_warning(
		lame_als(Y, R = 0, family = "normal", verbose = FALSE),
		"no observed entries")
})

test_that("a degenerate binary slice triggers a warning", {
	set.seed(83)
	n = 16; Tt = 3
	Y = lapply(seq_len(Tt), function(t) {
		m = matrix(rbinom(n * n, 1, 0.4), n, n); diag(m) = NA; m
	})
	Y[[2]][!is.na(Y[[2]])] = 0                 # slice 2 a single distinct value
	expect_warning(
		lame_als(Y, R = 0, family = "binary", verbose = FALSE),
		"single distinct observed")
})

test_that("the parametric bootstrap preserves the data's missingness pattern", {
	set.seed(81)
	n = 24; Tt = 4
	base = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
	holed = lapply(base, function(m) {
		off = which(!is.na(m))
		m[sample(off, round(0.5 * length(off)))] = NA
		m
	})
	fit_full = lame_als(base, R = 0, family = "normal", verbose = FALSE)
	fit_miss = suppressWarnings(
		lame_als(holed, R = 0, family = "normal", verbose = FALSE))
	se_full = ame_als_bootstrap(fit_full, R = 80, type = "parametric",
	                                 seed = 1, verbose = FALSE)$se[["intercept"]]
	se_miss = ame_als_bootstrap(fit_miss, R = 80, type = "parametric",
	                                 seed = 1, verbose = FALSE)$se[["intercept"]]
	# replicate networks keep the holes, so the sparser fit has a larger
	# bootstrap se; it would be about equal if the holes were silently filled
	expect_gt(se_miss, se_full * 1.15)
})

test_that("disconnected graphs get a reproducible per-component additive gauge", {
	set.seed(82)
	n = 24
	Y = matrix(NA_real_, n, n)
	Y[1:12, 1:12]   = rnorm(144)
	Y[13:24, 13:24] = rnorm(144)
	diag(Y) = NA
	fit = suppressWarnings(
		ame_als(Y, R = 0, family = "normal", verbose = FALSE))
	# the gauge is fit-preserving and deterministic: a warm-start refit
	# reproduces the additive effects exactly
	refit = suppressWarnings(ame_als_refit(fit, verbose = FALSE))
	expect_equal(unname(fit$a), unname(refit$a), tolerance = 1e-5)
	expect_equal(unname(fit$b), unname(refit$b), tolerance = 1e-5)
	# each component carries the precision-weighted minimum-norm gauge:
	# sum(cnt_row * a) == sum(cnt_col * b) within the component
	cr = rowSums(!is.na(Y)); cc = colSums(!is.na(Y))
	for (idx in list(1:12, 13:24)) {
		expect_lt(abs(sum(cr[idx] * fit$a[idx]) -
		              sum(cc[idx] * fit$b[idx])), 1e-6)
	}
})

test_that("near-collinear designs are solved stably (eigen-based solver)", {
	set.seed(84)
	n = 30
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	x1 = matrix(rnorm(n * n), n, n)
	x2 = x1 + matrix(rnorm(n * n, 0, 1e-6), n, n)   # near-collinear, not equal
	Xd = array(c(x1, x2), dim = c(n, n, 2))
	fit = suppressWarnings(ame_als(Y, Xdyad = Xd, R = 0,
	                        family = "normal", verbose = FALSE))
	expect_true(all(is.finite(coef(fit))))
	# the minimum-norm eigen solve keeps coefficients sane; a bare solve() on
	# the near-singular normal equations could blow them up
	expect_lt(max(abs(coef(fit))), 5)
})

# ==========================================================================
# opt-in mode: als / hybrid low-rank solver
# ==========================================================================

test_that("ALS and hybrid low-rank solvers reach the same optimum as MM", {
	set.seed(90)
	n = 30; Tt = 6
	U0 = matrix(rnorm(n * 2), n, 2); V0 = matrix(rnorm(n * 2), n, 2)
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Yl = lapply(seq_len(Tt), function(t) {
		Y = 0.3 + outer(a, b, "+") + tcrossprod(U0, V0) +
			matrix(rnorm(n * n, 0, 0.5), n, n)
		diag(Y) = NA
		if (t > 1) {                                  # unbalanced panel
			m = sample(which(!is.na(Y)), round(0.5 * n * n))
			Y[m] = NA
		}
		Y
	})
	fmm = lame_als(Yl, R = 2, family = "normal",
	                    lowrank_method = "mm", verbose = FALSE)
	fal = lame_als(Yl, R = 2, family = "normal",
	                    lowrank_method = "als", verbose = FALSE)
	fhy = lame_als(Yl, R = 2, family = "normal",
	                    lowrank_method = "hybrid", verbose = FALSE)
	expect_true(fal$converged && fhy$converged)
	# all three minimise the same objective -> same deviance and fitted values
	expect_equal(fal$deviance, fmm$deviance, tolerance = 1e-4)
	expect_equal(fhy$deviance, fmm$deviance, tolerance = 1e-4)
	expect_equal(cor(c(fal$EZ[[1]]), c(fmm$EZ[[1]]), use = "complete.obs"),
	             1, tolerance = 1e-4)
	expect_identical(fal$lowrank_method, "als")
})

test_that("ALS works for bipartite networks", {
	set.seed(91)
	Yb = lapply(seq_len(4), function(t) matrix(rnorm(28 * 20), 28, 20))
	fit = lame_als(Yb, R = 2, family = "normal", mode = "bipartite",
	                    lowrank_method = "als", verbose = FALSE)
	expect_true(fit$converged)
	expect_equal(dim(fit$U), c(28L, 2L))
	expect_equal(dim(fit$V), c(20L, 2L))
})

test_that("ALS is declined for symmetric models with a warning", {
	d = sim_ame(n = 22, Tt = 1, R = 1, symmetric = TRUE, seed = 92)
	expect_warning({
		fit = ame_als(d$Ylist[[1]], R = 1, family = "normal",
		              symmetric = TRUE, lowrank_method = "als",
		              verbose = FALSE)
	},
		"not available for symmetric")
	expect_identical(fit$lowrank_method, "mm")   # fell back
	expect_false(is.null(fit$L))
})

# ==========================================================================
# opt-in mode: irls for non-normal families
# ==========================================================================

test_that("IRLS recovers the Poisson coefficient better than the transform", {
	set.seed(93)
	n = 40
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	EZ = 0.4 + 0.5 * Xd + outer(a, b, "+")
	Y = matrix(rpois(n * n, exp(EZ)), n, n); diag(Y) = NA
	ft = ame_als(Y, Xdyad = Xd, R = 0, family = "poisson",
	                  non_normal_method = "transform", verbose = FALSE)
	fi = ame_als(Y, Xdyad = Xd, R = 0, family = "poisson",
	                  non_normal_method = "irls", verbose = FALSE)
	expect_true(fi$converged)
	expect_identical(fi$non_normal_method, "irls")
	expect_identical(fi$link, "log")
	# log(y+1) attenuates the slope; irls is closer to the true poisson coef
	expect_lt(abs(coef(fi)["dyad1_dyad"] - 0.5), abs(coef(ft)["dyad1_dyad"] - 0.5))
	expect_equal(unname(coef(fi)["dyad1_dyad"]), 0.5, tolerance = 0.15)
})

test_that("IRLS runs for binary with logit and probit links", {
	set.seed(94)
	n = 38
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	EZ = -0.3 + 1.0 * Xd + outer(a, b, "+")
	Y = 1 * (EZ + matrix(rnorm(n * n), n, n) > 0); diag(Y) = NA
	for (lk in c("logit", "probit")) {
		fi = ame_als(Y, Xdyad = Xd, R = 0, family = "binary",
		                  non_normal_method = "irls", link = lk,
		                  verbose = FALSE)
		expect_true(fi$converged)
		expect_identical(fi$link, lk)
		expect_gt(coef(fi)["dyad1_dyad"], 0)            # positive effect recovered
		expect_true(all(fitted(fi) >= 0 & fitted(fi) <= 1, na.rm = TRUE))
	}
})

test_that("the moving-block bootstrap draws contiguous blocks", {
	set.seed(101)
	n = 18; Tt = 9
	Yl = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
	fit = lame_als(Yl, R = 1, family = "normal", verbose = FALSE)
	b1 = ame_als_bootstrap(fit, R = 40, type = "block",
	                            block_length = 1, seed = 1, verbose = FALSE)
	b3 = ame_als_bootstrap(fit, R = 40, type = "block",
	                            block_length = 3, seed = 1, verbose = FALSE)
	expect_identical(b1$block_length, 1L)
	expect_identical(b3$block_length, 3L)
	expect_true(all(is.finite(b1$se)) && all(is.finite(b3$se)))
	expect_true(all(b3$ci_lo <= b3$ci_hi))
	# a malformed block length is rejected
	expect_error(
		ame_als_bootstrap(fit, R = 5, type = "block",
		                      block_length = 0, verbose = FALSE),
		"block_length")
})

test_that("vcov gives a sandwich covariance for the regression coefficients", {
	set.seed(100)
	n = 40; Tt = 4
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	parts = lapply(seq_len(Tt), function(t) {
		Xd = matrix(rnorm(n * n), n, n)
		Y = 0.3 + 0.6 * Xd + outer(a, b, "+") +
			matrix(rnorm(n * n, 0, 0.7), n, n)
		diag(Y) = NA
		list(Y = Y, X = Xd)
	})
	fit = lame_als(lapply(parts, `[[`, "Y"),
	                    Xdyad = lapply(parts, `[[`, "X"),
	                    R = 0, family = "normal", verbose = FALSE)
	V = vcov(fit)
	expect_equal(dim(V), c(2L, 2L))
	expect_equal(rownames(V), c("intercept", "dyad1_dyad"))
	expect_true(isSymmetric(V))
	expect_true(all(diag(V) > 0))
	# hc0 and dyad-clustered both run and are sane
	Vh = vcov(fit, cluster = "none")
	expect_true(all(diag(Vh) > 0))
	# the sandwich se is in the same ballpark as the bootstrap se
	bt = ame_als_bootstrap(fit, R = 100, type = "block",
	                            seed = 1, verbose = FALSE)
	ratio = sqrt(diag(V))[["dyad1_dyad"]] / bt$se[["dyad1_dyad"]]
	expect_gt(ratio, 0.3)
	expect_lt(ratio, 3)
})

test_that("multi-start is reproducible, RNG-safe, and never worse than one start", {
	set.seed(99)
	n = 26; Tt = 5
	U0 = matrix(rnorm(n * 2), n, 2); V0 = matrix(rnorm(n * 2), n, 2)
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Yl = lapply(seq_len(Tt), function(t) {
		Y = 0.3 + outer(a, b, "+") + tcrossprod(U0, V0) +
			matrix(rnorm(n * n, 0, 0.8), n, n)
		diag(Y) = NA
		Y
	})
	f0 = lame_als(Yl, R = 2, family = "normal",
	                   multistart = "none", verbose = FALSE)
	fc = lame_als(Yl, R = 2, family = "normal",
	                   multistart = "cheap", verbose = FALSE, seed = 123)
	# multi-start keeps the lowest-sse fit, so it cannot be worse
	expect_lte(fc$deviance, f0$deviance + 1e-6)
	expect_length(fc$multistart_sse, 4L)
	# reproducible given the seed, and the caller's rng stream is untouched
	set.seed(7); before = runif(2)
	fc2 = lame_als(Yl, R = 2, family = "normal",
	                    multistart = "cheap", verbose = FALSE, seed = 123)
	after = runif(2)
	set.seed(7); ref = runif(4)
	expect_equal(fc$deviance, fc2$deviance, tolerance = 1e-8)
	expect_equal(c(before, after), ref)
})

test_that("the QR and auto linear solvers agree with eigen", {
	set.seed(97)
	n = 30
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	x1 = matrix(rnorm(n * n), n, n)
	x2 = matrix(rnorm(n * n), n, n)
	Xd = array(c(x1, 0.6 * x1 + 0.8 * x2), dim = c(n, n, 2))
	fe = ame_als(Y, Xdyad = Xd, R = 0, family = "normal",
	                  linear_solver = "eigen", verbose = FALSE)
	fq = ame_als(Y, Xdyad = Xd, R = 0, family = "normal",
	                  linear_solver = "qr", verbose = FALSE)
	fa = ame_als(Y, Xdyad = Xd, R = 0, family = "normal",
	                  linear_solver = "auto", verbose = FALSE)
	# all solvers minimise the same least-squares problem
	expect_equal(unname(coef(fq)), unname(coef(fe)), tolerance = 1e-6)
	expect_equal(unname(coef(fa)), unname(coef(fe)), tolerance = 1e-6)
	expect_identical(fq$linear_solver, "qr")
})

test_that("the QR solver handles a rank-deficient design", {
	set.seed(98)
	n = 25
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	x1 = matrix(rnorm(n * n), n, n)
	Xdup = array(c(x1, x1), dim = c(n, n, 2))      # exactly collinear
	fit = suppressWarnings(ame_als(Y, Xdyad = Xdup, R = 0,
	                        family = "normal", linear_solver = "qr",
	                        verbose = FALSE))
	expect_true(all(is.finite(coef(fit))))
})

test_that("an IRLS fit bootstraps (refit re-runs the reweighting loop)", {
	set.seed(96)
	n = 22; Tt = 5
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Yl = lapply(seq_len(Tt), function(t) {
		Xd = matrix(rnorm(n * n), n, n)
		m = matrix(rpois(n * n, exp(0.3 + 0.4 * Xd + outer(a, b, "+"))), n, n)
		diag(m) = NA
		list(Y = m, X = Xd)
	})
	fit = lame_als(lapply(Yl, `[[`, "Y"), Xdyad = lapply(Yl, `[[`, "X"),
	                    R = 0, family = "poisson", non_normal_method = "irls",
	                    verbose = FALSE)
	rf = ame_als_refit(fit, verbose = FALSE)
	expect_identical(rf$non_normal_method, "irls")
	expect_equal(coef(rf), coef(fit), tolerance = 0.05)
	bt = ame_als_bootstrap(fit, R = 30, type = "block",
	                            seed = 1, verbose = FALSE)
	expect_s3_class(bt, "boot_ame")
	expect_true(all(is.finite(bt$se)))
})

# ==========================================================================
# regression tests for the review-driven hardening
# ==========================================================================

test_that("the weighted BCD core is invariant to a global weight rescaling", {
	set.seed(110)
	n = 20; Tt = 3
	Z = array(rnorm(n * n * Tt), dim = c(n, n, Tt))
	for (t in seq_len(Tt)) diag(Z[, , t]) = NA
	X = array(rnorm(n * n * Tt), dim = c(n, n, 1, Tt))
	W = array(runif(n * n * Tt, 0.2, 2), dim = c(n, n, Tt))
	f1 = lame:::.ae_bcd_fit(Z, X, 0, FALSE, 100, 1e-8, weights = W)
	f2 = lame:::.ae_bcd_fit(Z, X, 0, FALSE, 100, 1e-8, weights = 7 * W)
	# scaling every weight by a constant does not change the wls solution
	expect_equal(f1$mu, f2$mu, tolerance = 1e-6)
	expect_equal(f1$beta, f2$beta, tolerance = 1e-6)
	expect_equal(f1$a, f2$a, tolerance = 1e-6)
})

test_that("a zero observation weight is equivalent to a missing cell", {
	set.seed(111)
	n = 18; Tt = 2
	Z = array(rnorm(n * n * Tt), dim = c(n, n, Tt))
	for (t in seq_len(Tt)) diag(Z[, , t]) = NA
	X = array(0, dim = c(n, n, 0, Tt))
	W = array(1, dim = c(n, n, Tt)); W[!is.finite(Z)] = 0
	drop = cbind(sample(2:n, 5), sample(2:n, 5), sample(Tt, 5, TRUE))
	Wz = W; Zna = Z
	for (r in seq_len(nrow(drop))) {
		Wz[drop[r, 1], drop[r, 2], drop[r, 3]] = 0
		Zna[drop[r, 1], drop[r, 2], drop[r, 3]] = NA
	}
	f_zero = lame:::.ae_bcd_fit(Z,   X, 0, FALSE, 100, 1e-8, weights = Wz)
	f_miss = lame:::.ae_bcd_fit(Zna, X, 0, FALSE, 100, 1e-8, weights = W)
	expect_equal(f_zero$mu, f_miss$mu, tolerance = 1e-6)
	expect_equal(f_zero$a,  f_miss$a,  tolerance = 1e-6)
})

test_that("IRLS stays finite on a separated binary network", {
	set.seed(112)
	n = 30
	Xd = matrix(rnorm(n * n), n, n)
	Y = 1 * (Xd > 0); diag(Y) = NA            # perfectly separated
	for (lk in c("logit", "probit")) {
		fit = suppressWarnings(
			ame_als(Y, Xdyad = Xd, R = 0, family = "binary",
			            non_normal_method = "irls", link = lk,
			            verbose = FALSE))
		expect_true(all(is.finite(coef(fit))))
		expect_true(all(is.finite(fit$a)) && all(is.finite(fit$b)))
		ez = unlist(fit$EZ)
		expect_true(all(is.finite(ez) | is.na(ez)))
	}
})

test_that("vcov uses the IRLS weights for an IRLS fit", {
	set.seed(113)
	n = 34; Tt = 3
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	parts = lapply(seq_len(Tt), function(t) {
		Xd = matrix(rnorm(n * n), n, n)
		m = matrix(rpois(n * n, exp(0.3 + 0.4 * Xd + outer(a, b, "+"))), n, n)
		diag(m) = NA
		list(Y = m, X = Xd)
	})
	fit = lame_als(lapply(parts, `[[`, "Y"),
	                   Xdyad = lapply(parts, `[[`, "X"),
	                   R = 0, family = "poisson", non_normal_method = "irls",
	                   verbose = FALSE)
	expect_false(is.null(fit$obs_weights))    # final irls weights stored
	expect_warning({
		V = vcov(fit)
	}, "conditional sandwich")
	expect_equal(dim(V), c(2L, 2L))
	expect_true(isSymmetric(V))
	expect_true(all(diag(V) > 0))
})

test_that("a degenerate moving-block length warns", {
	set.seed(114)
	n = 14; Tt = 6
	Yl = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
	fit = lame_als(Yl, R = 0, family = "normal", verbose = FALSE)
	expect_warning(
		ame_als_bootstrap(fit, R = 10, type = "block", block_length = 5,
		                      seed = 1, verbose = FALSE),
		"fewer than 3 blocks")
})

# ==========================================================================
# regression tests for the user-feedback-driven fixes
# ==========================================================================

test_that("IRLS converges for R > 0 (penalized ALS factor solver)", {
	set.seed(120)
	n = 30
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	U0 = matrix(rnorm(n), n, 1); V0 = matrix(rnorm(n), n, 1)
	EZ = -0.3 + 1.0 * Xd + outer(a, b, "+") + tcrossprod(U0, V0)
	Yb = 1 * (EZ + matrix(rnorm(n * n), n, n) > 0); diag(Yb) = NA
	fb = ame_als(Yb, Xdyad = Xd, R = 1, family = "binary",
	                 non_normal_method = "irls", verbose = FALSE)
	expect_true(fb$converged)
	# a directed R > 0 IRLS fit solves the penalized (MAP) factor sub-problem
	# via the ridged ALS row solver (was the unpenalized mm->hybrid switch
	# before the penalty landed; that solver overfit binary outcomes)
	expect_identical(fb$lowrank_method, "als")
	Yp = matrix(rpois(n * n, exp(0.3 + 0.4 * Xd + outer(a, b, "+") +
	                             tcrossprod(U0, V0))), n, n); diag(Yp) = NA
	fp = ame_als(Yp, Xdyad = Xd, R = 1, family = "poisson",
	                 non_normal_method = "irls", verbose = FALSE)
	expect_true(fp$converged)
})

# --------------------------------------------------------------------------
# regression: penalized (MAP) factor block for the IRLS families
#
# before the penalty, the unpenalized rank-R binary/poisson factor block was
# an unregularized rank-R GLM MLE, which diverges by quasi-separation: on a
# null-latent DGP the fitted var(UV') ran away (10+ on binary, ~1-3 on
# poisson) and dragged the binary slope up by +37%. the penalty gives each
# factor row the ridge implied by the AME prior with an EM-updated scale, so
# the null latent variance shrinks and real structure is still recovered.
# --------------------------------------------------------------------------

test_that("binary R>0 IRLS does not run away on a null-latent DGP", {
	# no latent structure in the DGP: the penalized factor block must shrink
	# var(UV') toward zero and keep the slope near truth, instead of the
	# unpenalized runaway (var(UV') ~ 10, slope +37%) it replaced
	n = 30
	slopes = numeric(3); vars = numeric(3)
	for (s in 1:3) {
		set.seed(s)
		X = matrix(rnorm(n * n), n, n); diag(X) = NA
		Y = 1 * (-1.1 + 0.5 * X + matrix(rnorm(n * n), n, n) > 0); diag(Y) = NA
		f = suppressWarnings(suppressMessages(
			ame_als(Y, Xdyad = array(X, c(n, n, 1)), R = 2, family = "binary",
			        verbose = FALSE)))
		O = tcrossprod(f$U, f$V)
		slopes[s] = unname(coef(f)[2]); vars[s] = var(O[row(O) != col(O)])
	}
	# the latent variance no longer diverges (old unpenalized: ~10)
	expect_true(all(vars < 2))
	# the slope stays bounded near the truth of 0.5 (old unpenalized: ~0.7)
	expect_true(mean(slopes) > 0.3 && mean(slopes) < 0.75)
})

test_that("binary R>0 IRLS recovers real latent structure without collapse", {
	# strong latent structure (u,v ~ N(0, 1), R=2): the penalized fit must
	# still recover it -- var(UV') stays well above 0 (no MAP collapse) and
	# the fitted UV' correlates with the truth
	n = 30
	cors = numeric(3); vars = numeric(3)
	for (s in 1:3) {
		set.seed(100 + s)
		X = matrix(rnorm(n * n), n, n); diag(X) = NA
		U0 = matrix(rnorm(n * 2, 0, 1), n, 2)
		V0 = matrix(rnorm(n * 2, 0, 1), n, 2)
		O0 = tcrossprod(U0, V0)
		Y = 1 * (-1.1 + 0.5 * X + O0 + matrix(rnorm(n * n), n, n) > 0)
		diag(Y) = NA
		f = suppressWarnings(suppressMessages(
			ame_als(Y, Xdyad = array(X, c(n, n, 1)), R = 2, family = "binary",
			        verbose = FALSE)))
		O = tcrossprod(f$U, f$V)
		off = row(O) != col(O)
		cors[s] = cor(O[off], O0[off]); vars[s] = var(O[off])
	}
	expect_true(all(vars > 0.3))       # not collapsed to zero
	expect_true(mean(cors) > 0.6)      # real structure recovered
})

test_that("poisson R>0 IRLS shrinks the null latent variance", {
	# poisson shares the penalized factor block: on a null-latent count DGP
	# var(UV') shrinks toward zero while the slope stays unbiased
	n = 30
	slopes = numeric(3); vars = numeric(3)
	for (s in 1:3) {
		set.seed(200 + s)
		X = matrix(rnorm(n * n), n, n)
		lambda = exp(-0.5 + 0.3 * X)
		Y = matrix(rpois(n * n, c(lambda)), n, n)   # rate finite (no NA diag)
		diag(X) = NA; diag(Y) = NA                   # self-ties undefined
		f = suppressWarnings(suppressMessages(
			ame_als(Y, Xdyad = array(X, c(n, n, 1)), R = 2, family = "poisson",
			        verbose = FALSE)))
		O = tcrossprod(f$U, f$V)
		slopes[s] = unname(coef(f)[2]); vars[s] = var(O[row(O) != col(O)])
	}
	expect_true(all(vars < 0.5))       # old unpenalized: ~1-3
	expect_true(mean(slopes) > 0.2 && mean(slopes) < 0.4)
})

test_that("normal-family ALS is unchanged by the IRLS penalty (R=0 exact LS)", {
	# the penalty is confined to the IRLS families; the normal path must stay
	# the exact least-squares solution at R = 0
	set.seed(11)
	n = 30
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	X = array(rnorm(n * n), c(n, n, 1))
	Y = X[, , 1] + outer(a, b, "+") + matrix(rnorm(n * n), n, n); diag(Y) = NA
	nm = sprintf("a%02d", 1:n); rownames(Y) = colnames(Y) = nm
	dimnames(X) = list(nm, nm, "x")
	als = suppressWarnings(ame_als(Y, Xdyad = X, R = 0, family = "normal",
	                               verbose = FALSE))
	idx = which(row(Y) != col(Y))
	ols = lm(y ~ x + r + c, data = data.frame(
		y = Y[idx], x = X[, , 1][idx],
		r = factor(row(Y)[idx]), c = factor(col(Y)[idx])))
	# coefficients match OLS with actor fixed effects
	expect_equal(unname(coef(als)["x_dyad"]), unname(coef(ols)["x"]),
	             tolerance = 1e-4)
	# ve is now df-corrected and matches the df-corrected OLS residual variance
	s2_ols = sum(residuals(ols)^2) / (length(idx) - length(coef(ols)))
	expect_equal(als$s2, s2_ols, tolerance = 1e-4)
})

test_that("binary ALS reports latent-scale ve = 1 with a working-scale diagnostic", {
	# the probit model fixes the latent error variance at 1 by identification;
	# ALS reports that (matching the MCMC VC) and keeps the working-scale GLM
	# dispersion separately as ve_working
	set.seed(1)
	n = 30
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	X = array(rnorm(n * n), c(n, n, 1))
	Yb = 1 * ((-0.5 + X[, , 1] + outer(a, b, "+") +
	           matrix(rnorm(n * n), n, n)) > 0); diag(Yb) = NA
	f = suppressWarnings(ame_als(Yb, Xdyad = X, R = 0, family = "binary",
	                             verbose = FALSE))
	expect_equal(unname(f$VC["ve"]), 1)
	expect_false(is.null(f$ve_working))
	expect_true(is.finite(f$ve_working) && f$ve_working > 0)
})

test_that("multistart is honored for binary IRLS fits", {
	# before, the use_irls branch ran a single deterministic start and
	# silently ignored multistart; now it runs one IRLS trajectory per start
	set.seed(11)
	n = 30
	X1 = matrix(rnorm(n * n), n, n)
	Y = 1 * (-0.5 + X1 + matrix(rnorm(n * n), n, n) > 0); diag(Y) = NA
	f = suppressWarnings(suppressMessages(
		ame_als(Y, Xdyad = X1, R = 2, family = "binary", multistart = "full",
		        verbose = FALSE)))
	expect_identical(f$multistart, "full")
	expect_false(is.null(f$multistart_sse))
	expect_length(f$multistart_sse, 8L)
})

test_that("confint.ame_als returns Wald intervals with a conditional note", {
	d = sim_ame(n = 25, Tt = 1, R = 0, seed = 121)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 0,
	                  family = "normal", verbose = FALSE)
	expect_message({
		ci = confint(fit)
	}, "anti-conservative")
	expect_true(is.matrix(ci))
	expect_equal(rownames(ci), c("intercept", "dyad1_dyad"))
	expect_true(all(ci[, 1] <= ci[, 2]))
})

test_that("vcov.boot_ame returns the bootstrap coefficient covariance", {
	d = sim_ame(n = 20, Tt = 1, R = 1, seed = 122)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 1,
	                  family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                           seed = 1, verbose = FALSE)
	V = vcov(bt)
	expect_equal(dim(V), c(2L, 2L))
	expect_true(isSymmetric(V))
	expect_equal(unname(sqrt(diag(V))), unname(bt$se), tolerance = 1e-8)
})

test_that("logLik / AIC are not defined for a fast AME fit", {
	d = sim_ame(n = 15, Tt = 1, seed = 123)
	fit = ame_als(d$Ylist[[1]], R = 0, family = "normal", verbose = FALSE)
	expect_error(logLik(fit), "not defined")
	expect_error(AIC(fit), "not defined")
})

test_that("multiplicative factors carry dim column names", {
	d = sim_ame(n = 22, Tt = 1, R = 2, seed = 124)
	fit = ame_als(d$Ylist[[1]], R = 2, family = "normal", verbose = FALSE)
	expect_equal(colnames(fit$U), c("dim1", "dim2"))
	expect_equal(colnames(fit$V), c("dim1", "dim2"))
})

test_that("negative counts with the Poisson family raise an error", {
	set.seed(125)
	Y = matrix(rpois(400, 2), 20, 20); Y[3, 5] = -1; diag(Y) = NA
	expect_error(
		ame_als(Y, R = 0, family = "poisson", verbose = FALSE),
		"non-negative")
})

test_that("non-0/1 data with the binary family warns", {
	set.seed(126)
	Y = matrix(sample(0:2, 400, TRUE), 20, 20); diag(Y) = NA
	expect_warning(
		ame_als(Y, R = 0, family = "binary", verbose = FALSE),
		"expects 0/1")
})

test_that("a 3D Xdyad array for lame_als raises an informative error", {
	set.seed(127)
	n = 16; Tt = 4
	Y = lapply(seq_len(Tt), function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA; m
	})
	Xbad = array(rnorm(n * n * Tt), dim = c(n, n, Tt))   # ambiguous 3d
	expect_error(
		lame_als(Y, Xdyad = Xbad, family = "normal", verbose = FALSE),
		"ambiguous")
})

test_that("ame_als_refit accepts a plain matrix for Y_new", {
	d = sim_ame(n = 20, Tt = 1, R = 1, seed = 128)
	fit = ame_als(d$Ylist[[1]], R = 1, family = "normal", verbose = FALSE)
	rf = ame_als_refit(fit, Y_new = d$Ylist[[1]], verbose = FALSE)
	expect_s3_class(rf, "ame_als")
	expect_equal(coef(rf), coef(fit), tolerance = 0.02)
})

# ==========================================================================
# s3 methods
# ==========================================================================

test_that("S3 methods for ame_als and boot_ame run cleanly", {
	d = sim_ame(n = 20, Tt = 5, R = 1, seed = 29)
	fit = lame_als(d$Ylist, Xdyad = d$Xlist, R = 1,
	                    family = "normal", verbose = FALSE)
	expect_output(print(fit))
	expect_s3_class(summary(fit), "summary.ame_als")
	expect_output(print(summary(fit)))
	expect_invisible(sampler_describe(fit, verbose = FALSE))

	bt = ame_als_bootstrap(fit, R = 30, type = "block",
	                            seed = 10, verbose = FALSE)
	expect_output(print(bt))
	expect_s3_class(summary(bt), "summary.boot_ame")
	expect_output(print(summary(bt)))
	expect_invisible(sampler_describe(bt, verbose = FALSE))
	expect_true(inherits(bt, "boot_ame"))
})

# ==========================================================================
# regression tests: user-feedback fixes
# ==========================================================================

test_that("parametric bootstrap CIs bracket the point estimate (non-normal)", {
	set.seed(220)
	n = 35
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Xd = matrix(rnorm(n * n), n, n)
	Yb = 1 * (-0.3 + 0.8 * Xd + outer(a, b, "+") +
	          matrix(rnorm(n * n), n, n) > 0); diag(Yb) = NA
	for (nm in c("irls", "transform")) {
		fb = ame_als(Yb, Xdyad = Xd, R = 0, family = "binary",
		                 non_normal_method = nm, verbose = FALSE)
		bt = ame_als_bootstrap(fb, R = 80, type = "parametric",
		                           seed = 1, verbose = FALSE)
		expect_true(all(bt$point_est >= bt$ci_lo & bt$point_est <= bt$ci_hi),
		            info = nm)
	}
})

test_that("binary and poisson default to the IRLS method with a probit link", {
	d = sim_ame(n = 20, Tt = 1, seed = 232)
	Yb = 1 * (d$Ylist[[1]] > 0); diag(Yb) = NA
	fb = ame_als(Yb, R = 0, family = "binary", verbose = FALSE)
	expect_identical(fb$non_normal_method, "irls")
	expect_identical(fb$link, "probit")
})

test_that("ame() and lame() reject an unrecognised family", {
	d = sim_ame(n = 15, Tt = 1, seed = 221)
	expect_error(ame(d$Ylist[[1]], family = "Binary", nscan = 10, burn = 5,
	                 odens = 5, verbose = FALSE), "not a recognised family")
	dl = sim_ame(n = 12, Tt = 3, seed = 222)
	expect_error(lame(dl$Ylist, family = "banana", nscan = 10, burn = 5,
	                  odens = 5, verbose = FALSE), "not a recognised family")
})

test_that("lame() accepts Y matrices without actor names", {
	set.seed(223)
	Yl = lapply(1:4, function(t) {
		m = matrix(rnorm(144), 12, 12); diag(m) = NA; m
	})
	fit = lame(Yl, family = "normal", nscan = 20, burn = 10, odens = 5,
	           verbose = FALSE, gof = FALSE)
	expect_s3_class(fit, "lame")
	expect_length(fit$APM, 12)
})

test_that("ame() rejects degenerate MCMC settings and constant continuous Y", {
	d = sim_ame(n = 15, Tt = 1, seed = 233)
	expect_error(ame(d$Ylist[[1]], family = "normal", nscan = 0, burn = 5,
	                 odens = 1, verbose = FALSE), "nscan")
	expect_error(ame(matrix(5, 12, 12), family = "normal", nscan = 20,
	                 burn = 5, odens = 5, verbose = FALSE), "constant")
})

test_that("vcov.ame_als names omitted node-covariate coefficients", {
	d = sim_ame(n = 22, Tt = 1, seed = 224)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]],
	                  Xrow = matrix(rnorm(22), 22, 1), R = 0,
	                  family = "normal", verbose = FALSE)
	expect_message({
		V = vcov(fit)
	}, "Node-covariate")
	expect_lt(nrow(V), length(coef(fit)))
})

test_that("linear_solver = 'qr' is stable on a near-collinear design", {
	set.seed(225); n = 60
	x1 = matrix(rnorm(n * n), n, n)
	x2 = x1 + matrix(rnorm(n * n, 0, 1e-5), n, n)
	X = array(c(x1, x2), dim = c(n, n, 2))
	Y = 0.5 + 0.4 * x1 + matrix(rnorm(n * n), n, n); diag(Y) = NA
	cf = coef(suppressWarnings(
		ame_als(Y, Xdyad = X, R = 0, linear_solver = "qr",
		            verbose = FALSE)))
	expect_true(all(abs(cf) < 50))
})

test_that("predict / residuals / plot methods work for ame_als", {
	d = sim_ame(n = 20, Tt = 1, R = 1, seed = 226)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 1,
	                  family = "normal", verbose = FALSE)
	pr = predict(fit)
	expect_true(is.matrix(pr))
	rs = residuals(fit)
	expect_true(is.matrix(rs))
	expect_equal(unname(rs), unname(fit$Y[, , 1] - fitted(fit)),
	             tolerance = 1e-8)
	pf = tempfile(fileext = ".pdf"); grDevices::pdf(pf)
	expect_invisible(plot(fit))
	grDevices::dev.off()
})

test_that("coef() works on a boot_ame object", {
	d = sim_ame(n = 20, Tt = 1, R = 0, seed = 227)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 0,
	                  family = "normal", verbose = FALSE)
	bt = ame_als_bootstrap(fit, R = 40, type = "parametric",
	                           seed = 1, verbose = FALSE)
	cf = coef(bt)
	expect_equal(cf, bt$point_est)
	expect_named(cf, bt$param_names)
})

test_that("confint rejects an out-of-range level", {
	d = sim_ame(n = 20, Tt = 1, R = 0, seed = 228)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 0,
	                  family = "normal", verbose = FALSE)
	expect_error(suppressMessages(confint(fit, level = -0.2)), "level")
	expect_error(suppressMessages(confint(fit, level = 1.5)), "level")
})

test_that("ame_als rejects a non-integer R", {
	d = sim_ame(n = 18, Tt = 1, seed = 229)
	expect_error(ame_als(d$Ylist[[1]], R = 1.7, family = "normal",
	                         verbose = FALSE), "non-negative integer")
})

test_that("ame_als_bootstrap does not mutate the global RNG", {
	d = sim_ame(n = 18, Tt = 1, R = 0, seed = 230)
	fit = ame_als(d$Ylist[[1]], Xdyad = d$Xlist[[1]], R = 0,
	                  family = "normal", verbose = FALSE)
	set.seed(99); before = runif(1)
	set.seed(99)
	invisible(ame_als_bootstrap(fit, R = 20, type = "parametric",
	                                seed = 7, verbose = FALSE))
	after = runif(1)
	expect_equal(before, after)
})

test_that("uv_plot and latent_positions accept an ame_als fit", {
	d = sim_ame(n = 20, Tt = 1, R = 2, seed = 231)
	fit = ame_als(d$Ylist[[1]], R = 2, family = "normal", verbose = FALSE)
	lp = latent_positions(fit)
	expect_s3_class(lp, "data.frame")
	expect_gt(nrow(lp), 0)
	p = uv_plot(fit = fit)
	expect_true(inherits(p, "ggplot") || inherits(p, "gg") ||
	            inherits(p, "patchwork"))
})

# --------------------------------------------------------------------------
# summary.ame_als integration with bootstrap cis
# --------------------------------------------------------------------------

test_that("summary.ame_als merges bootstrap CIs into the regression coef table", {
	set.seed(7)
	nA = 12L; nB = 8L
	Y = matrix(rbinom(nA*nB, 1, 0.3), nA, nB)
	X = array(rnorm(nA*nB*2), c(nA, nB, 2))
	fit = suppressWarnings(ame_als(Y, Xdyad = X, R = 1, family = "binary",
	                                mode = "bipartite", bootstrap = 50L,
	                                verbose = FALSE))
	s = summary(fit)
	# the coefficient table should now have 4 columns when a bootstrap is
	# attached: estimate / std.error / ci 2.5% / ci 97.5%
	expect_true("Std.Error" %in% colnames(s$coefficients))
	expect_true("CI 2.5%"   %in% colnames(s$coefficients))
	expect_true("CI 97.5%"  %in% colnames(s$coefficients))
	# every row should have finite estimate and std.error
	expect_true(all(is.finite(s$coefficients[, "Estimate"])))
	expect_true(all(is.finite(s$coefficients[, "Std.Error"])))
})

test_that("summary.ame_als without bootstrap shows only Estimate column", {
	set.seed(7)
	nA = 8L; nB = 6L
	Y = matrix(rnorm(nA*nB), nA, nB)
	X = array(rnorm(nA*nB*2), c(nA, nB, 2))
	fit = suppressWarnings(ame_als(Y, Xdyad = X, R = 0, family = "normal",
	                                mode = "bipartite", verbose = FALSE))
	s = summary(fit)
	expect_equal(colnames(s$coefficients), "Estimate")
})
