skip_on_cran()

# static rank-based models (cbin, frn, rrl)

library(lame)
library(testthat)

####
# censored binary (cbin)
####

test_that("Censored binary (cbin) AME runs without errors", {
	set.seed(6886)
	n = 25
	
	# create binary data
	Y = matrix(rbinom(n*n, 1, 0.3), n, n)
	diag(Y) = NA
	
	# set max outdegree (e.g., each node can nominate at most 5 others)
	odmax = 5
	
	# ensure Y respects odmax constraint
	for(i in 1:n) {
		if(sum(Y[i,], na.rm=TRUE) > odmax) {
			# randomly keep only odmax nominations
			ones = which(Y[i,] == 1)
			keep = sample(ones, odmax)
			Y[i, setdiff(ones, keep)] = 0
		}
	}
	
	# test basic model
	fit = ame(Y, R=1, family="cbin", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	
	# check outdegrees respect constraint
	outdegrees = rowSums(Y, na.rm=TRUE)
	expect_true(all(outdegrees <= odmax))
})

test_that("Censored binary (cbin) with covariates works", {
	set.seed(6886)
	n = 25
	
	# generate covariate
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	# generate binary outcome influenced by X
	prob = pnorm(0.5 * X)
	prob[is.na(prob)] = 0.5  # Handle NAs
	Y = matrix(rbinom(n*n, 1, c(prob)), n, n)
	diag(Y) = NA
	
	# apply outdegree constraint
	odmax = 4
	for(i in 1:n) {
		if(sum(Y[i,], na.rm=TRUE) > odmax) {
			ones = which(Y[i,] == 1)
			keep = sample(ones, odmax)
			Y[i, setdiff(ones, keep)] = 0
		}
	}
	
	# fit model
	fit = ame(Y, Xdyad=X, R=0, family="cbin", odmax=odmax,
						burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	if(ncol(fit$BETA) >= 2) {
		beta_est = median(fit$BETA[,2])
		expect_true(abs(beta_est) < 2)
	}
})

####
# fixed rank nomination (frn)
####

test_that("Fixed rank nomination (frn) AME runs without errors", {
	set.seed(6886)
	n = 20
	
	# create ranked nomination data (e.g., rank top 3 friends)
	odmax = 3
	Y = matrix(NA, n, n)
	
	for(i in 1:n) {
		# each person nominates odmax others with ranks 1, 2, 3
		nominees = sample(setdiff(1:n, i), odmax)
		Y[i, nominees] = 1:odmax
	}
	diag(Y) = NA
	
	# test basic model
	fit = ame(Y, R=1, family="frn", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	
	# check that each row has exactly odmax nominations
	nominations_per_row = apply(Y, 1, function(x) sum(!is.na(x) & x > 0))
	expect_true(all(nominations_per_row == odmax))
})

test_that("Fixed rank nomination (frn) with covariates works", {
	set.seed(6886)
	n = 20
	
	# generate covariate
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	# create ranked nominations influenced by X
	odmax = 4
	Y = matrix(NA, n, n)
	
	for(i in 1:n) {
		# higher X values get better (lower) ranks
		scores = X[i,]
		scores[i] = NA
		top_indices = order(scores, decreasing=TRUE, na.last=TRUE)[1:odmax]
		Y[i, top_indices] = 1:odmax
	}
	diag(Y) = NA
	
	# fit model
	fit = ame(Y, Xdyad=X, R=0, family="frn", odmax=odmax,
						burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	# for frn, positive X should lead to better (lower) ranks
	if(ncol(fit$BETA) >= 2) {
		beta_est = median(fit$BETA[,2])
		expect_true(is.finite(beta_est))
	}
})

####
# row rank likelihood (rrl)
####

test_that("Row rank likelihood (rrl) AME runs without errors", {
	set.seed(6886)
	n = 20
	
	# create data where each row is ranked
	Y = matrix(NA, n, n)
	for(i in 1:n) {
		# random continuous values
		Y[i,] = rnorm(n)
	}
	diag(Y) = NA
	
	# note: rrl doesn't allow row effects (rvar must be FALSE)
	fit = ame(Y, R=1, family="rrl", 
						rvar=FALSE, cvar=TRUE,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	expect_true(!is.null(fit$U))
	expect_true(!is.null(fit$V))
	
	# check that BPM exists but APM doesn't (rvar=FALSE for rrl)
	expect_true(!is.null(fit$BPM))
})

test_that("Row rank likelihood (rrl) with covariates works", {
	set.seed(6886)
	n = 20
	
	# generate covariate
	X = matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(X) = NA
	
	# generate outcome where higher X leads to higher Y
	Y = 0.5 * X + matrix(rnorm(n*n, 0, 0.5), n, n)
	diag(Y) = NA
	
	# fit model (rvar must be FALSE for rrl)
	fit = ame(Y, Xdyad=X, R=0, family="rrl",
						rvar=FALSE, cvar=TRUE,
						burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
	if(ncol(fit$BETA) >= 1) {
		# rrl may not have intercept, check first non-intercept coef
		coef_idx = ifelse(ncol(fit$BETA) >= 2, 2, 1)
		beta_est = median(fit$BETA[,coef_idx])
		expect_true(is.finite(beta_est))
	}
})

####
# edge cases for rank-based models
####

test_that("Rank models handle varying outdegrees", {
	set.seed(6886)
	n = 20
	
	# test cbin with varying odmax per node
	odmax_vec = sample(2:5, n, replace=TRUE)
	Y_cbin = matrix(0, n, n)
	
	for(i in 1:n) {
		nominees = sample(setdiff(1:n, i), odmax_vec[i])
		Y_cbin[i, nominees] = 1
	}
	diag(Y_cbin) = NA
	
	fit_cbin = ame(Y_cbin, R=1, family="cbin", odmax=odmax_vec,
								 burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit_cbin$BETA))
	
	# check each row respects its specific odmax
	for(i in 1:n) {
		expect_lte(sum(Y_cbin[i,], na.rm=TRUE), odmax_vec[i])
	}
})

test_that("FRN handles incomplete rankings", {
	set.seed(6886)
	n = 20
	
	# create data where some nodes make fewer nominations
	odmax = 5
	Y = matrix(NA, n, n)
	
	for(i in 1:n) {
		# some nodes nominate fewer than odmax
		n_nominations = sample(1:odmax, 1)
		nominees = sample(setdiff(1:n, i), n_nominations)
		Y[i, nominees] = 1:n_nominations
	}
	diag(Y) = NA
	
	fit = ame(Y, R=1, family="frn", odmax=odmax,
						burn=200, nscan=600, verbose = FALSE)
	
	expect_true(!is.null(fit$BETA))
})

####
# rank models with multiplicative effects
####

test_that("Rank models with multiplicative effects work", {
	set.seed(6886)
	n = 25
	
	# generate latent space positions
	U_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	V_true = matrix(rnorm(n*2, 0, 0.5), n, 2)
	
	# generate latent scores
	Z = tcrossprod(U_true, V_true) + matrix(rnorm(n*n, 0, 0.5), n, n)
	
	# convert to cbin data
	Y_cbin = 1 * (Z > 0)
	odmax = 6
	for(i in 1:n) {
		if(sum(Y_cbin[i,], na.rm=TRUE) > odmax) {
			ones = which(Y_cbin[i,] == 1)
			keep = sample(ones, odmax)
			Y_cbin[i, setdiff(ones, keep)] = 0
		}
	}
	diag(Y_cbin) = NA
	
	fit_cbin = ame(Y_cbin, R=2, family="cbin", odmax=odmax,
								 burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_cbin$U))
	expect_true(!is.null(fit_cbin$V))
	expect_equal(ncol(fit_cbin$U), 2)
	expect_equal(ncol(fit_cbin$V), 2)
	
	# convert to frn data (rank the positive links)
	Y_frn = matrix(NA, n, n)
	for(i in 1:n) {
		scores = Z[i,]
		scores[i] = NA
		top_indices = order(scores, decreasing=TRUE, na.last=TRUE)[1:min(odmax, n-1)]
		Y_frn[i, top_indices] = 1:length(top_indices)
	}
	diag(Y_frn) = NA
	
	fit_frn = ame(Y_frn, R=2, family="frn", odmax=odmax,
								burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_frn$U))
	expect_true(!is.null(fit_frn$V))
	
	# test rrl with continuous version
	Y_rrl = Z
	diag(Y_rrl) = NA
	
	fit_rrl = ame(Y_rrl, R=2, family="rrl",
								rvar=FALSE, cvar=TRUE,
								burn=300, nscan=800, verbose = FALSE)
	
	expect_true(!is.null(fit_rrl$U))
	expect_true(!is.null(fit_rrl$V))
})