# Regression tests for the "fix batch D" MCMC-correctness defects.
# Each guards against a previously-confirmed bug re-appearing:
#   #1  symmetric ordinal Z kernel used variance 1/2 -> coefs attenuated ~1/sqrt(2)
#   #11 lame() symmetric double-counted every dyad -> dyadic-coef sd too small by sqrt(2)
#   #12 ame() duplicated the posterior_opts storage block -> every draw stored twice
#   #19 ame() never sampled Sab for poisson -> additive variance frozen at start
# Kept fast (small n, short chains); assertions have wide tolerances so they
# fail only if the underlying defect is genuinely back.

test_that("#1 symmetric ordinal dyadic coef is not attenuated (~1/sqrt(2))", {
	skip_on_cran()
	set.seed(3); n = 30
	Xd = matrix(rnorm(n * n), n, n); Xd = (Xd + t(Xd)) / sqrt(2)
	Em = matrix(rnorm(n * n), n, n); Em = (Em + t(Em)) / sqrt(2)
	Zl = 1.0 * Xd + Em
	Y = matrix(cut(Zl, quantile(Zl, seq(0, 1, length = 7), na.rm = TRUE),
	                labels = FALSE, include.lowest = TRUE), n, n)
	diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	Xa = array(Xd, c(n, n, 1), dimnames = list(rownames(Y), colnames(Y), "x1"))
	fs = ame(Y, Xdyad = Xa, R = 0, symmetric = TRUE, family = "ordinal",
	          burn = 300, nscan = 1500, odens = 5, verbose = FALSE, gof = FALSE, seed = 1)
	# the symmetric = FALSE fit on exactly-symmetric Y warns by design
	fa = suppressWarnings(
		ame(Y, Xdyad = Xa, R = 0, symmetric = FALSE, family = "ordinal",
		    burn = 300, nscan = 1500, odens = 5, verbose = FALSE, gof = FALSE, seed = 1))
	ratio = mean(fs$BETA[, "x1_dyad"]) / mean(fa$BETA[, "x1_dyad"])
	# variance-1/2 bug drives this to ~0.71; variance-1 keeps it near 1.
	expect_gt(ratio, 0.85)
	expect_lt(ratio, 1.2)
})

test_that("#11 lame() symmetric dyadic-coef sd matches ame() (not sqrt(2) smaller)", {
	skip_on_cran()
	set.seed(7); n <- 35
	Xd = matrix(rnorm(n * n), n, n); Xd <- (Xd + t(Xd)) / 2
	E = matrix(rnorm(n * n), n, n); E <- (E + t(E)) / sqrt(2)
	av = rnorm(n, 0, 0.4)
	Y = 0.3 + 0.8 * Xd + outer(av, av, "+") + E
	Y = (Y + t(Y)) / 2; diag(Y) <- NA
	dimnames(Y) <- list(paste0("a", 1:n), paste0("a", 1:n))
	Xa = array(Xd, c(n, n, 1), dimnames = list(rownames(Y), colnames(Y), "x1"))
	fa = ame(Y, Xdyad = Xa, R = 0, symmetric = TRUE, family = "normal",
	          burn = 400, nscan = 3000, odens = 2, verbose = FALSE, gof = FALSE, seed = 5)
	fl = lame(list(t1 = Y), Xdyad = list(t1 = Xa), R = 0, symmetric = TRUE, family = "normal",
	           burn = 400, nscan = 3000, odens = 2, verbose = FALSE, gof = FALSE, seed = 5)
	sd_ratio = sd(fa$BETA[, "x1_dyad"]) / sd(fl$BETA[, "x1_dyad"])
	# double-count bug gives ~1.41; correct treatment gives ~1.0.
	expect_gt(sd_ratio, 0.8)
	expect_lt(sd_ratio, 1.25)
	# additive variance should also agree closely between the two entry points
	expect_equal(mean(fl$VC[, "va"]), mean(fa$VC[, "va"]), tolerance = 0.4)
})

test_that("#12 ame() stores each posterior draw once (no duplication, no crash)", {
	skip_on_cran()
	set.seed(1); n <- 18
	Y = matrix(rnorm(n), n) %*% t(matrix(rnorm(n), n)) + matrix(rnorm(n * n), n, n)
	diag(Y) <- NA
	dimnames(Y) <- list(paste0("a", 1:n), paste0("a", 1:n))
	f = ame(Y, R = 1, family = "normal", burn = 60, nscan = 240, odens = 4,
	         verbose = FALSE, gof = FALSE, seed = 2,
	         posterior_opts = list(save_UV = TRUE, save_ab = TRUE, thin_UV = 1, thin_ab = 1))
	n_eff = nrow(f$BETA)
	k = dim(f$U_samples)[3]
	expect_equal(k, n_eff)
	# consecutive pairs must not be identical draws (they were, pre-fix)
	frac_dup = mean(sapply(seq(1, k - 1, 2),
	                        function(i) isTRUE(all.equal(f$U_samples[, , i], f$U_samples[, , i + 1]))))
	expect_lt(frac_dup, 0.1)

	# symmetric + partially-specified posterior_opts must not crash, and
	# V_samples must stay NULL for a symmetric fit
	Ys = Y; Ys[lower.tri(Ys)] <- t(Ys)[lower.tri(Ys)]
	fsym = expect_no_error(
		ame(Ys, R = 1, symmetric = TRUE, family = "normal", burn = 40, nscan = 80,
		    odens = 2, verbose = FALSE, gof = FALSE, seed = 3,
		    posterior_opts = list(save_UV = TRUE, thin_UV = 1)))
	expect_null(fsym$V_samples)
})

test_that("#19 ame() samples Sab for poisson (additive variance not frozen)", {
	skip_on_cran()
	set.seed(3); n <- 28
	a = rnorm(n, 0, 0.5); b <- rnorm(n, 0, 0.5)
	X = array(rnorm(n * n), c(n, n, 1))
	eta = 0.2 + 0.5 * X[, , 1] + outer(a, b, "+")
	Y = matrix(rpois(n * n, exp(pmin(eta, 5))), n, n); diag(Y) <- NA
	nm = sprintf("a%02d", 1:n); rownames(Y) <- colnames(Y) <- nm
	dimnames(X) <- list(nm, nm, "x")
	f = ame(Y, Xdyad = X, R = 0, family = "poisson", nscan = 1200, burn = 400,
	         odens = 5, verbose = FALSE, gof = FALSE, seed = 7)
	# frozen Sab => constant va column (sd == 0); a live sampler moves it.
	expect_gt(sd(f$VC[, "va"]), 1e-4)
	# and va should be nowhere near the frozen ~0.05 start value
	expect_gt(mean(f$VC[, "va"]), 0.1)
})
