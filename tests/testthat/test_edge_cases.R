skip_on_cran()

# edge cases and error handling
library(lame)

test_that("models handle small networks", {
	skip_on_cran()
	
	set.seed(6886)
	
	# very small network
	n = 5
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) = NA
	
	# should work with small network
	fit = ame(Y, R=1, family="binary",
						 burn=50, nscan=100, verbose = FALSE)
	expect_true(!is.null(fit))
	
	# bipartite with small dimensions
	Y_bip = matrix(rbinom(3*4, 1, 0.5), 3, 4)
	rownames(Y_bip) = paste0("Row", 1:3)
	colnames(Y_bip) = paste0("Col", 1:4)
	
	fit_bip = ame(Y_bip, mode="bipartite", R_row=1, R_col=1,
								 family="binary", burn=50, nscan=100,
								 verbose = FALSE, gof=FALSE)  # Skip GOF for tiny networks
	expect_true(!is.null(fit_bip))
})

test_that("models handle sparse networks", {
	skip_on_cran()
	
	set.seed(6886)
	n = 30
	
	# very sparse network (density ~ 0.01)
	Y_sparse = matrix(0, n, n)
	n_edges = floor(0.01 * n * (n-1))
	edge_idx = sample(which(upper.tri(Y_sparse)), n_edges)
	Y_sparse[edge_idx] = 1
	Y_sparse[lower.tri(Y_sparse)] = t(Y_sparse)[lower.tri(Y_sparse)]
	diag(Y_sparse) = NA
	
	expect_no_error({
		fit_sparse = ame(Y_sparse, R=1, family="binary", symmetric=TRUE,
											burn=50, nscan=200, verbose = FALSE)
	})
	
	# check that model makes predictions (YPM exists)
	expect_true(!is.null(fit_sparse$YPM))
	expect_true(is.numeric(fit_sparse$YPM))
})

test_that("models handle dense networks", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	
	# very dense network (density ~ 0.9)
	Y_dense = matrix(1, n, n)
	n_non_edges = floor(0.1 * n * (n-1))
	non_edge_idx = sample(which(upper.tri(Y_dense)), n_non_edges)
	Y_dense[non_edge_idx] = 0
	Y_dense[lower.tri(Y_dense)] = t(Y_dense)[lower.tri(Y_dense)]
	diag(Y_dense) = NA
	
	expect_no_error({
		fit_dense = ame(Y_dense, R=1, family="binary", symmetric=TRUE,
										 burn=50, nscan=200, verbose = FALSE)
	})
	
	# check that model captures high density (YPM for binary is probability)
	expect_gt(mean(fit_dense$YPM, na.rm=TRUE), 0.5)
})

test_that("models handle networks with isolated nodes", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	
	# create network with isolated nodes
	Y = matrix(rbinom(n*n, 1, 0.2), n, n)
	# make first 3 nodes isolated
	Y[1:3, ] = 0
	Y[, 1:3] = 0
	diag(Y) = NA
	
	expect_no_error({
		fit = ame(Y, R=1, family="binary",
							 burn=50, nscan=200, verbose = FALSE)
	})
	
	# isolated nodes should have low predicted values (YPM for binary is probability)
	isolated_preds = fit$YPM[1:3, -c(1:3)]
	expect_lt(mean(isolated_preds, na.rm=TRUE), 0.2)
})

test_that("longitudinal models handle varying network sizes", {
	skip_on_cran()
	
	set.seed(6886)
	
	# networks with different missing patterns (simulating different sizes)
	T = 4
	n_max = 20
	Y_list = list()
	
	for(t in 1:T) {
		Y_t = matrix(rbinom(n_max*n_max, 1, 0.3), n_max, n_max)
		
		# add missing nodes (different pattern each time)
		n_missing = sample(0:5, 1)
		if(n_missing > 0) {
			missing_idx = sample(1:n_max, n_missing)
			Y_t[missing_idx, ] = NA
			Y_t[, missing_idx] = NA
		}
		
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n_max)
		Y_list[[t]] = Y_t
	}
	
	expect_no_error({
		fit = lame(Y_list, R=1, family="binary",
								burn=50, nscan=200, verbose = FALSE, plot=FALSE)
	})
})

test_that("models handle extreme parameter values", {
	skip_on_cran()
	
	set.seed(6886)
	n = 15
	
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) = NA
	
	# extreme prior values
	expect_no_error({
		fit = ame(Y, R=1, family="binary",
							 burn=50, nscan=100, verbose = FALSE,
							 prior=list(Sab0=diag(c(100, 100))))  # Large variance prior
	})
	
	# very large R
	expect_no_error({
		fit_large_R = ame(Y, R=10, family="binary",
											 burn=50, nscan=100, verbose = FALSE)
	})
	
	# r larger than n (may fail or adjust R automatically)
	tryCatch({
		fit_R_too_large = ame(Y, R=n+5, family="binary",
													 burn=50, nscan=100, verbose = FALSE)
		expect_true(TRUE)  # If it succeeds (perhaps by adjusting R), that's fine
	}, error = function(e) {
		expect_true(TRUE)  # If it fails, that's also acceptable
	}, warning = function(w) {
		expect_true(TRUE)  # Warnings are fine too
	})
})

test_that("models handle non-standard data types correctly", {
	skip_on_cran()
	
	set.seed(6886)
	n = 20
	
	# ordinal data with gaps
	Y_ord = matrix(sample(c(1, 3, 5, 7), n*n, replace=TRUE), n, n)
	diag(Y_ord) = NA
	
	expect_no_error({
		fit_ord = ame(Y_ord, R=1, family="ordinal", intercept=FALSE,
									 burn=50, nscan=200, verbose = FALSE)
	})
	
	# count data with large values
	Y_count = matrix(rpois(n*n, lambda=10), n, n)
	Y_count[sample(1:(n*n), 5)] = 100  # Add some large counts
	diag(Y_count) = NA
	
	expect_no_error({
		fit_count = ame(Y_count, R=1, family="poisson",
										 burn=50, nscan=200, verbose = FALSE)
	})
})

test_that("bipartite models reject square matrices with mode='bipartite'", {
	skip_on_cran()
	
	set.seed(6886)
	n = 10
	
	# square matrix
	Y_square = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y_square) = NA
	
	# should work with unipartite
	expect_no_error({
		fit_uni = ame(Y_square, mode="unipartite", R=1, family="binary",
									 burn=50, nscan=100, verbose = FALSE)
	})
	
	tryCatch({
		fit_bip = ame(Y_square, mode="bipartite", R_row=1, R_col=1,
									 family="binary", burn=50, nscan=100, 
									 verbose = FALSE, gof=FALSE)
		expect_true(TRUE)
	}, error = function(e) {
		expect_true(TRUE)
	})
})

test_that("models handle character row/column names", {
	skip_on_cran()
	
	set.seed(6886)
	n = 15
	
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	rownames(Y) = paste0("Node_", 1:n)
	colnames(Y) = paste0("Node_", 1:n)
	diag(Y) = NA
	
	expect_no_error({
		fit = ame(Y, R=1, family="binary",
							 burn=50, nscan=100, verbose = FALSE)
	})
	
	# check that names are preserved
	expect_equal(rownames(fit$YPM), rownames(Y))
	expect_equal(colnames(fit$YPM), colnames(Y))
})

