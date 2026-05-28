# correctness-contract tests: assert the numbers, not just that a fit runs.
# the contracts that must hold for any fit, across method x mode x family:
#   C1  coef() == posterior mean of BETA
#   C2  confint() == posterior quantiles of BETA
#   C3  fitted() == predict(type = "response")
#   C4  residuals() == observed - fitted (where defined)
#   C5  predict(newdata = the fit's own covariates) reproduces the in-sample
#       linear predictor (guards the positional-indexing bug class)
#   C6  a +delta shift on a covariate moves the link predictor by exactly
#       delta * that covariate's OWN coefficient (counterfactual sign/scale)
#   C7  known-DGP parameter recovery: the truth lies in the 95% CI
#   C8  ame() and lame() agree on the default g-prior for the same data
#
# Fits are small for speed; the assertions are exact where the math is exact.

# ---- helpers ---------------------------------------------------------------

.mk_uni_ame <- function(seed = 1, n = 16, nodal = TRUE, family = "binary") {
	set.seed(seed)
	Xd <- array(rnorm(n*n*2), c(n, n, 2), dimnames = list(NULL, NULL, c("d1","d2")))
	Y <- if (family == "normal") {
		0.5 + 0.8*Xd[,,1] - 0.6*Xd[,,2] + matrix(rnorm(n*n, 0, 0.5), n, n)
	} else {
		eta <- -0.3 + 0.9*Xd[,,1] - 0.7*Xd[,,2]
		matrix(rbinom(n*n, 1, pnorm(eta)), n, n)
	}
	diag(Y) <- NA; rownames(Y) <- colnames(Y) <- paste0("a", seq_len(n))
	args <- list(Y = Y, Xdyad = Xd, family = family, R = 1,
	             burn = 60, nscan = 400, odens = 5, verbose = FALSE)
	if (nodal) {
		args$Xrow <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "r1"))
		args$Xcol <- matrix(rnorm(n), n, 1, dimnames = list(NULL, "c1"))
	}
	list(fit = do.call(ame, args), Y = Y, Xd = Xd)
}

# ---- C1 / C2: coef and confint are the posterior summaries -----------------

test_that("C1/C2: coef == posterior mean and confint == posterior quantiles", {
	skip_on_cran()
	obj <- .mk_uni_ame(); fit <- obj$fit
	expect_equal(unname(coef(fit)), unname(colMeans(fit$BETA)), tolerance = 1e-10)
	ci <- confint(fit)
	q <- t(apply(fit$BETA, 2, quantile, c(0.025, 0.975)))
	expect_equal(unname(as.matrix(ci)), unname(q), tolerance = 1e-8)
})

# ---- C3 / C4: fitted, predict, residuals are mutually consistent -----------

test_that("C3/C4: fitted == predict(response) and residuals == Y - fitted", {
	skip_on_cran()
	for (fam in c("binary", "normal")) {
		obj <- .mk_uni_ame(family = fam); fit <- obj$fit
		fv <- fitted(fit)
		pr <- predict(fit, type = "response")
		expect_equal(as.numeric(fv), as.numeric(pr), tolerance = 1e-8,
		             info = paste("fitted==predict", fam))
		rs <- residuals(fit)
		# residuals defined as observed - fitted on observed cells
		obs <- obj$Y
		expect_equal(as.numeric(rs[!is.na(obs)]),
		             as.numeric((obs - fv)[!is.na(obs)]), tolerance = 1e-6,
		             info = paste("resid==Y-fitted", fam))
	}
})

# ---- C5: in-sample reproduction (the positional-indexing guard) ------------

test_that("C5: predict(newdata = fit's own Xdyad) reproduces the in-sample EZ", {
	skip_on_cran()
	# with nodal covariates present -- the case that exposed the bug
	obj <- .mk_uni_ame(nodal = TRUE); fit <- obj$fit
	ez_in <- predict(fit, type = "link")
	ez_nd <- predict(fit, type = "link", newdata = obj$Xd)
	expect_lt(max(abs(ez_in - ez_nd), na.rm = TRUE), 1e-8)
	# bipartite too
	set.seed(3); nA <- 12; nB <- 9
	Yb <- matrix(rbinom(nA*nB, 1, 0.3), nA, nB)
	rownames(Yb) <- paste0("a", seq_len(nA)); colnames(Yb) <- paste0("b", seq_len(nB))
	Xb <- array(rnorm(nA*nB*2), c(nA, nB, 2), dimnames = list(NULL, NULL, c("d1","d2")))
	fb <- ame(Yb, Xdyad = Xb, mode = "bipartite", R_row = 1, R_col = 1,
	          family = "binary", burn = 50, nscan = 300, odens = 5, verbose = FALSE)
	expect_lt(max(abs(predict(fb, type = "link") -
	                  predict(fb, type = "link", newdata = Xb)), na.rm = TRUE), 1e-8)
})

# ---- C6: counterfactual sign and scale -------------------------------------

test_that("C6: +delta on a covariate moves the link by delta * its own coef", {
	skip_on_cran()
	obj <- .mk_uni_ame(nodal = TRUE); fit <- obj$fit
	bm <- colMeans(fit$BETA)
	base <- predict(fit, type = "link", newdata = obj$Xd)
	for (slice in 1:2) {
		nm <- paste0("d", slice, "_dyad")
		Xup <- obj$Xd; Xup[,,slice] <- Xup[,,slice] + 2.5
		up <- predict(fit, type = "link", newdata = Xup)
		expect_equal(mean(up - base, na.rm = TRUE), unname(2.5 * bm[nm]),
		             tolerance = 1e-6, info = nm)
	}
})

# ---- C7: known-DGP parameter recovery (truth inside the 95% CI) ------------

test_that("C7: a strong dyadic effect is recovered inside its 95% CI", {
	skip_on_cran()
	set.seed(11); n <- 30
	Xd <- array(rnorm(n*n), c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	a <- rnorm(n, 0, 0.4); b <- rnorm(n, 0, 0.4)
	eta <- -0.5 + 1.0 * Xd[,,1] + outer(a, b, "+")
	Y <- matrix(rbinom(n*n, 1, pnorm(eta)), n, n); diag(Y) <- NA
	rownames(Y) <- colnames(Y) <- paste0("a", seq_len(n))
	fit <- ame(Y, Xdyad = Xd, family = "binary", R = 0,
	           burn = 200, nscan = 1500, odens = 5, verbose = FALSE)
	ci <- confint(fit)["x1_dyad", ]
	expect_true(ci[1] <= 1.0 && 1.0 <= ci[2])   # truth 1.0 inside the CI
	expect_gt(coef(fit)["x1_dyad"], 0.5)        # right sign and ballpark
})

# ---- C8: ame() and lame() agree on the default g-prior ---------------------

test_that("C8: ame() and lame() use the same default g for the same data", {
	skip_on_cran()
	# normal family, large-magnitude outcome: both must variance-scale g and
	# both must recover the slope (the bug was lame using a count-only g and
	# silently shrinking beta to ~0).
	set.seed(5); n <- 20
	Xd <- matrix(rnorm(n*n), n, n)
	Y <- 300 * (0.7 * Xd) + matrix(rnorm(n*n, 0, 150), n, n)
	diag(Y) <- NA; rownames(Y) <- colnames(Y) <- paste0("a", seq_len(n))
	Xa <- array(Xd, c(n, n, 1), dimnames = list(NULL, NULL, "x1"))
	f_ame <- ame(Y, Xdyad = Xa, family = "normal", R = 0,
	             burn = 100, nscan = 600, odens = 10, verbose = FALSE)
	# lame with a single-period list reduces to the same cross-section
	f_lame <- lame(list(t1 = Y), Xdyad = list(t1 = Xa), family = "normal", R = 0,
	               burn = 100, nscan = 600, odens = 10, verbose = FALSE)
	truth <- 300 * 0.7
	# coef() shape differs by fit: a length-p vector for a static fit, a
	# p x T matrix for a dynamic_beta fit. A single-period lame() fit is
	# static, so coef() is a vector -- index by name only.
	cl <- coef(f_lame)
	lame_dyad <- if (is.matrix(cl)) cl["x1_dyad", 1] else cl["x1_dyad"]
	expect_gt(coef(f_ame)["x1_dyad"], 0.5 * truth)   # ame recovers
	expect_gt(lame_dyad, 0.5 * truth)                # lame recovers (was ~0)
})
