skip_on_cran()

# static bipartite ame checks

test_that("supported bipartite families fit", {
	skip_on_cran()
	
	set.seed(456)
	
	# set up bipartite dimensions
	nA = 12
	nB = 10
	
	# test the supported rectangular families
	families = c("normal", "binary", "tobit", "rrl")

	for(fam in families) {
		# generate test data
		if(fam == "normal") {
			Y = matrix(rnorm(nA * nB, 0, 1), nA, nB)
		} else if(fam == "tobit") {
			Y = matrix(rnorm(nA * nB, 0, 1), nA, nB)
			Y[Y < 0] = 0
		} else if(fam == "binary" || fam == "cbin") {
			Y = matrix(rbinom(nA * nB, 1, 0.5), nA, nB)
		} else if(fam == "poisson") {
			Y = matrix(rpois(nA * nB, 3), nA, nB)
		} else if(fam == "ordinal") {
			Y = matrix(sample(0:3, nA * nB, replace = TRUE), nA, nB)
		} else if(fam == "frn") {
			# fixed rank nomination
			Y = matrix(0, nA, nB)
			for(i in 1:nA) {
				nominated = sample(1:nB, min(3, nB))
				Y[i, nominated] = 1:length(nominated)
			}
		} else if(fam == "rrl") {
			# row rank likelihood
			Y = matrix(0, nA, nB)
			for(i in 1:nA) {
				Y[i,] = sample(1:nB)
			}
		}
		
		# set odmax for rank families
		if(fam %in% c("cbin", "frn")) {
			odmax = rep(3, nA)
		} else {
			odmax = NULL
		}
		
		# fit model
		suppressWarnings({
			fit = ame(Y, mode = "bipartite", family = fam,
								 R_row = 0, R_col = 0,
								 odmax = odmax,
								 burn = 100, nscan = 300, 
								 verbose = FALSE)
		})
		
		# basic checks
		expect_equal(fit$mode, "bipartite", 
								 info = paste("Mode check failed for", fam))
		expect_equal(fit$family, fam, 
								 info = paste("Family check failed for", fam))
		expect_equal(fit$nA, nA, 
								 info = paste("Row dimension check failed for", fam))
		expect_equal(fit$nB, nB, 
								 info = paste("Column dimension check failed for", fam))
		
		# posterior samples are present
		expect_false(is.null(fit$BETA), 
								 info = paste("BETA missing for", fam))
		
		# family-specific checks
		if(fam == "binary" || fam == "cbin") {
			# predictions are binary probabilities
			expect_true(all(fit$YPM >= 0 & fit$YPM <= 1), 
									info = paste("YPM not in [0,1] for", fam))
		}
		
		if(fam == "tobit") {
			# censoring at zero
			expect_true(all(fit$YPM >= 0), 
									info = paste("Negative predictions for tobit"))
		}
	}
})

test_that("bipartite models with multiplicative effects fit", {
	skip_on_cran()
	
	set.seed(457)
	
	nA = 10
	nB = 8
	R_row = 1
	R_col = 1
	
	# test the supported subset with multiplicative effects
	families = c("normal", "binary")

	for(fam in families) {
		# generate data
		if(fam == "normal") {
			Y = matrix(rnorm(nA * nB), nA, nB)
		} else if(fam == "binary") {
			Y = matrix(rbinom(nA * nB, 1, 0.5), nA, nB)
		} else if(fam == "poisson") {
			Y = matrix(rpois(nA * nB, 2), nA, nB)
		}
		
		# fit with multiplicative effects
		suppressWarnings({
			fit = ame(Y, mode = "bipartite", family = fam,
								 R_row = R_row, R_col = R_col,
								 burn = 100, nscan = 300,
								 verbose = FALSE)
		})
		
		# multiplicative effects are present
		expect_equal(dim(fit$U), c(nA, R_row), 
								 info = paste("U dimension wrong for", fam))
		expect_equal(dim(fit$V), c(nB, R_col), 
								 info = paste("V dimension wrong for", fam))
		expect_equal(dim(fit$G), c(R_row, R_col), 
								 info = paste("G dimension wrong for", fam))
	}
})

test_that("bipartite models handle covariates", {
	skip_on_cran()
	
	set.seed(458)
	
	nA = 8
	nB = 6
	
	# generate covariates
	Xdyad = array(rnorm(nA * nB * 2), c(nA, nB, 2))
	Xrow = matrix(rnorm(nA * 2), nA, 2)
	Xcol = matrix(rnorm(nB * 2), nB, 2)
	
	Y = matrix(rnorm(nA * nB), nA, nB)
	
	# fit with all covariate types
	fit = ame(Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
						 mode = "bipartite",
						 burn = 100, nscan = 300,
						 verbose = FALSE)
	
	# intercept + 2 dyadic + 2 row + 2 col
	expect_equal(ncol(fit$BETA), 7)
	
	# estimates are present
	beta_means = colMeans(fit$BETA)
	expect_false(any(is.na(beta_means)))
})

test_that("rank families handle their constraints", {
	skip_on_cran()
	
	set.seed(459)
	
	nA = 6
	nB = 8
	
	# test frn with odmax constraint
	Y_frn = matrix(0, nA, nB)
	odmax = rep(3, nA)
	for(i in 1:nA) {
		nominated = sample(1:nB, odmax[i])
		Y_frn[i, nominated] = 1:odmax[i]
	}
	
	# frn / ordinal / rrl use rectangular bipartite samplers
	# rectangular samplers in r/rz_bipartite.r back these families
	fit_frn = ame(Y_frn, mode = "bipartite", family = "frn", odmax = odmax,
	              burn = 5, nscan = 10, odens = 5, verbose = FALSE,
	              plot = FALSE, gof = FALSE)
	expect_equal(fit_frn$family, "frn")

	Y_rrl = matrix(0, nA, nB)
	for(i in 1:nA) Y_rrl[i, ] = sample(1:nB)
	fit_rrl = ame(Y_rrl, mode = "bipartite", family = "rrl",
	              burn = 5, nscan = 10, odens = 5, verbose = FALSE,
	              plot = FALSE, gof = FALSE)
	expect_equal(fit_rrl$family, "rrl")

	Y_ord = matrix(sample(0:4, nA * nB, replace = TRUE), nA, nB)
	fit_ord = ame(Y_ord, mode = "bipartite", family = "ordinal",
	              burn = 5, nscan = 10, odens = 5, verbose = FALSE,
	              plot = FALSE, gof = FALSE)
	expect_equal(fit_ord$family, "ordinal")
})
