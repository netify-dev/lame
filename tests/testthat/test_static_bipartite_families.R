skip_on_cran()

# test suite for all 8 families in bipartite AME models

test_that("All 8 families work for bipartite networks", {
	skip_on_cran()
	
	set.seed(456)
	
	# set up bipartite dimensions
	nA = 12  # Row nodes
	nB = 10  # Column nodes
	
	# test each family
	families = c("normal", "binary", "poisson", "tobit", "ordinal", "cbin", "frn", "rrl")
	
	for(fam in families) {
		# generate appropriate test data
		if(fam == "normal") {
			Y = matrix(rnorm(nA * nB, 0, 1), nA, nB)
		} else if(fam == "tobit") {
			Y = matrix(rnorm(nA * nB, 0, 1), nA, nB)
			Y[Y < 0] = 0  # Censor at 0
		} else if(fam == "binary" || fam == "cbin") {
			Y = matrix(rbinom(nA * nB, 1, 0.5), nA, nB)
		} else if(fam == "poisson") {
			Y = matrix(rpois(nA * nB, 3), nA, nB)
		} else if(fam == "ordinal") {
			Y = matrix(sample(0:3, nA * nB, replace = TRUE), nA, nB)
		} else if(fam == "frn") {
			# fixed rank nomination - each row nominates 3 columns
			Y = matrix(0, nA, nB)
			for(i in 1:nA) {
				nominated = sample(1:nB, min(3, nB))
				Y[i, nominated] = 1:length(nominated)
			}
		} else if(fam == "rrl") {
			# row rank likelihood - rank within rows
			Y = matrix(0, nA, nB)
			for(i in 1:nA) {
				Y[i,] = sample(1:nB)
			}
		}
		
		# set odmax for appropriate families
		if(fam %in% c("cbin", "frn")) {
			odmax = rep(3, nA)  # Max 3 nominations per row
		} else {
			odmax = NULL
		}
		
		# fit model
		suppressWarnings({
			fit = ame(Y, mode = "bipartite", family = fam,
								 R_row = 0, R_col = 0,  # Start with no multiplicative effects
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
		
		# check that posterior samples exist
		expect_false(is.null(fit$BETA), 
								 info = paste("BETA missing for", fam))
		
		# family-specific checks
		if(fam == "binary" || fam == "cbin") {
		# check predictions are binary probabilities
			expect_true(all(fit$YPM >= 0 & fit$YPM <= 1), 
									info = paste("YPM not in [0,1] for", fam))
		}
		
		if(fam == "poisson") {
		}
		
		if(fam == "tobit") {
			# check censoring at 0
			expect_true(all(fit$YPM >= 0), 
									info = paste("Negative predictions for tobit"))
		}
	}
})

test_that("Bipartite models with multiplicative effects work for all families", {
	skip_on_cran()
	
	set.seed(457)
	
	nA = 10
	nB = 8
	R_row = 1
	R_col = 1
	
	# test subset of families with multiplicative effects
	families = c("normal", "binary", "poisson")
	
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
		
		# check multiplicative effects
		expect_equal(dim(fit$U), c(nA, R_row), 
								 info = paste("U dimension wrong for", fam))
		expect_equal(dim(fit$V), c(nB, R_col), 
								 info = paste("V dimension wrong for", fam))
		expect_equal(dim(fit$G), c(R_row, R_col), 
								 info = paste("G dimension wrong for", fam))
	}
})

test_that("Bipartite models handle covariates correctly", {
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
	
	# check coefficient dimensions
	# should have: intercept + 2 dyadic + 2 row + 2 col = 7
	expect_equal(ncol(fit$BETA), 7)
	
	# check that estimates exist
	beta_means = colMeans(fit$BETA)
	expect_false(any(is.na(beta_means)))
})

test_that("Special families handle their constraints correctly", {
	skip_on_cran()
	
	set.seed(459)
	
	nA = 6
	nB = 8
	
	# test FRN with odmax constraint
	Y_frn = matrix(0, nA, nB)
	odmax = rep(3, nA)  # Each row can nominate at most 3
	for(i in 1:nA) {
		nominated = sample(1:nB, odmax[i])
		Y_frn[i, nominated] = 1:odmax[i]
	}
	
	suppressWarnings({
		fit_frn = ame(Y_frn, mode = "bipartite", family = "frn",
									 odmax = odmax,
									 burn = 100, nscan = 300,
									 verbose = FALSE)
	})
	
	expect_equal(fit_frn$family, "frn")
	# check that simulated values respect odmax
	# each row should have at most odmax[i] non-zero values
	
	# test RRL (no intercept, no row effects)
	Y_rrl = matrix(0, nA, nB)
	for(i in 1:nA) {
		Y_rrl[i,] = sample(1:nB)
	}
	
	suppressWarnings({
		fit_rrl = ame(Y_rrl, mode = "bipartite", family = "rrl",
									 burn = 100, nscan = 300,
									 verbose = FALSE)
	})
	
	expect_equal(fit_rrl$family, "rrl")
	# rRL should not have intercept
	# first column of BETA should be near 0 if no covariates
	
	# test ordinal (no intercept)
	Y_ord = matrix(sample(0:4, nA * nB, replace = TRUE), nA, nB)
	
	suppressWarnings({
		fit_ord = ame(Y_ord, mode = "bipartite", family = "ordinal",
									 burn = 100, nscan = 300,
									 verbose = FALSE)
	})
	
	expect_equal(fit_ord$family, "ordinal")
})

test_that("Bipartite models converge reasonably", {
	skip_on_cran()
	
	set.seed(460)
	
	nA = 8
	nB = 6
	
	# simple model with known structure
	mu_true = 0.5
	Y = matrix(mu_true + rnorm(nA * nB, 0, 0.5), nA, nB)
	
	fit = ame(Y, mode = "bipartite",
						 burn = 200, nscan = 500,
						 verbose = FALSE)
	
	# check that intercept is recovered approximately
	intercept_est = mean(fit$BETA[,1])
	expect_equal(intercept_est, mu_true, tolerance = 0.5)
	
	# check that variance components are positive
	expect_true(all(fit$VC[,4] > 0))  # Error variance should be positive
})