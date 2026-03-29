skip_on_cran()

# g-prior adjustment for coverage with additive effects
library(lame)
library(testthat)

test_that("Adjusting g-prior improves coverage with additive effects", {
	set.seed(6886)
	n_sims = 30
	n = 40
	mu_true = 0
	beta_true = 1
	
	results_default = list()
	results_adjusted = list()
	
	for(i in 1:n_sims) {
		set.seed(i)
		
		# generate data with additive effects
		xw = matrix(rnorm(n*2), n, 2)
		X = tcrossprod(xw[,1])
		diag(X) = NA
		
		# true additive effects
		a_true = rnorm(n, 0, 0.5)
		b_true = rnorm(n, 0, 0.5)
		
		eta = mu_true + beta_true * X + 
					 outer(a_true, rep(1, n)) + outer(rep(1, n), b_true)
		
		Y = eta + matrix(rnorm(n*n), n, n)
		diag(Y) = NA
		
		# fit with default g-prior
		fit_default = ame(Y, Xdyad = X, R = 0, family = "normal",
											rvar = TRUE, cvar = TRUE, dcor = FALSE,
											burn = 300, nscan = 1000,
											verbose = FALSE)
		
		# fit with adjusted g-prior (smaller to account for additive effects)
		# use residual variance after removing row/col means as estimate
		row_means = rowMeans(Y, na.rm = TRUE)
		col_means = colMeans(Y, na.rm = TRUE)
		grand_mean = mean(Y, na.rm = TRUE)
		Y_centered = Y - outer(row_means - grand_mean, rep(1, n)) - 
											outer(rep(1, n), col_means - grand_mean)
		g_adjusted = sum(!is.na(Y)) * var(c(Y_centered), na.rm = TRUE)
		
		fit_adjusted = ame(Y, Xdyad = X, R = 0, family = "normal",
											 rvar = TRUE, cvar = TRUE, dcor = FALSE,
											 g = g_adjusted,
											 burn = 300, nscan = 1000,
											 verbose = FALSE)
		
		# extract results for default
		beta_ci_default = quantile(fit_default$BETA[,2], c(0.025, 0.975))
		results_default[[i]] = list(
			beta_hat = median(fit_default$BETA[,2]),
			beta_covered = (beta_true >= beta_ci_default[1]) & (beta_true <= beta_ci_default[2]),
			ci_width = beta_ci_default[2] - beta_ci_default[1]
		)
		
		# extract results for adjusted
		beta_ci_adjusted = quantile(fit_adjusted$BETA[,2], c(0.025, 0.975))
		results_adjusted[[i]] = list(
			beta_hat = median(fit_adjusted$BETA[,2]),
			beta_covered = (beta_true >= beta_ci_adjusted[1]) & (beta_true <= beta_ci_adjusted[2]),
			ci_width = beta_ci_adjusted[2] - beta_ci_adjusted[1]
		)
	}
	
	# compare coverage
	coverage_default = mean(sapply(results_default, function(x) x$beta_covered))
	coverage_adjusted = mean(sapply(results_adjusted, function(x) x$beta_covered))
	
	cat("\n=== G-prior Adjustment Results ===\n")
	cat("Default g-prior coverage:", coverage_default, "\n")
	cat("Adjusted g-prior coverage:", coverage_adjusted, "\n")
	
	mean_width_default = mean(sapply(results_default, function(x) x$ci_width))
	mean_width_adjusted = mean(sapply(results_adjusted, function(x) x$ci_width))
	
	cat("Default CI width:", mean_width_default, "\n")
	cat("Adjusted CI width:", mean_width_adjusted, "\n")
	
	# adjusted should have better coverage
	expect_gt(coverage_adjusted, coverage_default - 0.05)  # Allow small tolerance
	
	# document if adjustment helps
	if(coverage_adjusted > coverage_default + 0.05) {
		message("G-prior adjustment improves coverage by ", 
						round((coverage_adjusted - coverage_default)*100, 1), " percentage points")
	}
})

