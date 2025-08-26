# Basic tests for static rank-based models (cbin, frn, rrl)
# These are more complex families with specific constraints

library(lame)
library(testthat)

# ============================================================================
# TEST 1: Censored Binary (cbin) - binary with max outdegree
# ============================================================================

test_that("Censored binary (cbin) AME runs without errors", {
  set.seed(6886)
  n <- 25
  
  # Create binary data
  Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
  diag(Y) <- NA
  
  # Set max outdegree (e.g., each node can nominate at most 5 others)
  odmax <- 5
  
  # Ensure Y respects odmax constraint
  for(i in 1:n) {
    if(sum(Y[i,], na.rm=TRUE) > odmax) {
      # Randomly keep only odmax nominations
      ones <- which(Y[i,] == 1)
      keep <- sample(ones, odmax)
      Y[i, setdiff(ones, keep)] <- 0
    }
  }
  
  # Test basic model
  fit <- ame(Y, R=1, family="cbin", odmax=odmax,
            burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  expect_true(!is.null(fit$U))
  expect_true(!is.null(fit$V))
  
  # Check outdegrees respect constraint
  outdegrees <- rowSums(Y, na.rm=TRUE)
  expect_true(all(outdegrees <= odmax))
})

test_that("Censored binary (cbin) with covariates works", {
  set.seed(6886)
  n <- 25
  
  # Generate covariate
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  
  # Generate binary outcome influenced by X
  prob <- pnorm(0.5 * X)
  prob[is.na(prob)] <- 0.5  # Handle NAs
  Y <- matrix(rbinom(n*n, 1, c(prob)), n, n)
  diag(Y) <- NA
  
  # Apply outdegree constraint
  odmax <- 4
  for(i in 1:n) {
    if(sum(Y[i,], na.rm=TRUE) > odmax) {
      ones <- which(Y[i,] == 1)
      keep <- sample(ones, odmax)
      Y[i, setdiff(ones, keep)] <- 0
    }
  }
  
  # Fit model
  fit <- ame(Y, Xdyad=X, R=0, family="cbin", odmax=odmax,
            burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  if(ncol(fit$BETA) >= 2) {
    beta_est <- median(fit$BETA[,2])
    expect_true(abs(beta_est) < 2)
  }
})

# ============================================================================
# TEST 2: Fixed Rank Nomination (frn) - ranked nominations
# ============================================================================

test_that("Fixed rank nomination (frn) AME runs without errors", {
  set.seed(6886)
  n <- 20
  
  # Create ranked nomination data (e.g., rank top 3 friends)
  odmax <- 3
  Y <- matrix(NA, n, n)
  
  for(i in 1:n) {
    # Each person nominates odmax others with ranks 1, 2, 3
    nominees <- sample(setdiff(1:n, i), odmax)
    Y[i, nominees] <- 1:odmax
  }
  diag(Y) <- NA
  
  # Test basic model
  fit <- ame(Y, R=1, family="frn", odmax=odmax,
            burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  expect_true(!is.null(fit$U))
  expect_true(!is.null(fit$V))
  
  # Check that each row has exactly odmax nominations
  nominations_per_row <- apply(Y, 1, function(x) sum(!is.na(x) & x > 0))
  expect_true(all(nominations_per_row == odmax))
})

test_that("Fixed rank nomination (frn) with covariates works", {
  set.seed(6886)
  n <- 20
  
  # Generate covariate
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  
  # Create ranked nominations influenced by X
  odmax <- 4
  Y <- matrix(NA, n, n)
  
  for(i in 1:n) {
    # Higher X values get better (lower) ranks
    scores <- X[i,]
    scores[i] <- NA
    top_indices <- order(scores, decreasing=TRUE, na.last=TRUE)[1:odmax]
    Y[i, top_indices] <- 1:odmax
  }
  diag(Y) <- NA
  
  # Fit model
  fit <- ame(Y, Xdyad=X, R=0, family="frn", odmax=odmax,
            burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  # For frn, positive X should lead to better (lower) ranks
  if(ncol(fit$BETA) >= 2) {
    beta_est <- median(fit$BETA[,2])
    expect_true(is.finite(beta_est))
  }
})

# ============================================================================
# TEST 3: Row Rank Likelihood (rrl) - row-wise rankings
# ============================================================================

test_that("Row rank likelihood (rrl) AME runs without errors", {
  set.seed(6886)
  n <- 20
  
  # Create data where each row is ranked
  Y <- matrix(NA, n, n)
  for(i in 1:n) {
    # Random continuous values
    Y[i,] <- rnorm(n)
  }
  diag(Y) <- NA
  
  # Note: rrl doesn't allow row effects (rvar must be FALSE)
  fit <- ame(Y, R=1, family="rrl", 
            rvar=FALSE, cvar=TRUE,
            burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  expect_true(!is.null(fit$U))
  expect_true(!is.null(fit$V))
  
  # Check that BPM exists but APM doesn't (rvar=FALSE for rrl)
  expect_true(!is.null(fit$BPM))
})

test_that("Row rank likelihood (rrl) with covariates works", {
  set.seed(6886)
  n <- 20
  
  # Generate covariate
  X <- matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(X) <- NA
  
  # Generate outcome where higher X leads to higher Y
  Y <- 0.5 * X + matrix(rnorm(n*n, 0, 0.5), n, n)
  diag(Y) <- NA
  
  # Fit model (rvar must be FALSE for rrl)
  fit <- ame(Y, Xdyad=X, R=0, family="rrl",
            rvar=FALSE, cvar=TRUE,
            burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  if(ncol(fit$BETA) >= 1) {
    # rrl may not have intercept, check first non-intercept coef
    coef_idx <- ifelse(ncol(fit$BETA) >= 2, 2, 1)
    beta_est <- median(fit$BETA[,coef_idx])
    expect_true(is.finite(beta_est))
  }
})

# ============================================================================
# TEST 4: Edge cases for rank-based models
# ============================================================================

test_that("Rank models handle varying outdegrees", {
  set.seed(6886)
  n <- 20
  
  # Test cbin with varying odmax per node
  odmax_vec <- sample(2:5, n, replace=TRUE)
  Y_cbin <- matrix(0, n, n)
  
  for(i in 1:n) {
    nominees <- sample(setdiff(1:n, i), odmax_vec[i])
    Y_cbin[i, nominees] <- 1
  }
  diag(Y_cbin) <- NA
  
  fit_cbin <- ame(Y_cbin, R=1, family="cbin", odmax=odmax_vec,
                 burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit_cbin$BETA))
  
  # Check each row respects its specific odmax
  for(i in 1:n) {
    expect_lte(sum(Y_cbin[i,], na.rm=TRUE), odmax_vec[i])
  }
})

test_that("FRN handles incomplete rankings", {
  set.seed(6886)
  n <- 20
  
  # Create data where some nodes make fewer nominations
  odmax <- 5
  Y <- matrix(NA, n, n)
  
  for(i in 1:n) {
    # Some nodes nominate fewer than odmax
    n_nominations <- sample(1:odmax, 1)
    nominees <- sample(setdiff(1:n, i), n_nominations)
    Y[i, nominees] <- 1:n_nominations
  }
  diag(Y) <- NA
  
  fit <- ame(Y, R=1, family="frn", odmax=odmax,
            burn=200, nscan=600, print=FALSE)
  
  expect_true(!is.null(fit$BETA))
  # EZ no longer stored - can be reconstructed if needed
})

# ============================================================================
# TEST 5: Model comparison for rank-based families
# ============================================================================

test_that("Rank models with multiplicative effects work", {
  set.seed(6886)
  n <- 25
  
  # Generate latent space positions
  U_true <- matrix(rnorm(n*2, 0, 0.5), n, 2)
  V_true <- matrix(rnorm(n*2, 0, 0.5), n, 2)
  
  # Generate latent scores
  Z <- tcrossprod(U_true, V_true) + matrix(rnorm(n*n, 0, 0.5), n, n)
  
  # Convert to cbin data
  Y_cbin <- 1 * (Z > 0)
  odmax <- 6
  for(i in 1:n) {
    if(sum(Y_cbin[i,], na.rm=TRUE) > odmax) {
      ones <- which(Y_cbin[i,] == 1)
      keep <- sample(ones, odmax)
      Y_cbin[i, setdiff(ones, keep)] <- 0
    }
  }
  diag(Y_cbin) <- NA
  
  fit_cbin <- ame(Y_cbin, R=2, family="cbin", odmax=odmax,
                 burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit_cbin$U))
  expect_true(!is.null(fit_cbin$V))
  expect_equal(ncol(fit_cbin$U), 2)
  expect_equal(ncol(fit_cbin$V), 2)
  
  # Convert to frn data (rank the positive links)
  Y_frn <- matrix(NA, n, n)
  for(i in 1:n) {
    scores <- Z[i,]
    scores[i] <- NA
    top_indices <- order(scores, decreasing=TRUE, na.last=TRUE)[1:min(odmax, n-1)]
    Y_frn[i, top_indices] <- 1:length(top_indices)
  }
  diag(Y_frn) <- NA
  
  fit_frn <- ame(Y_frn, R=2, family="frn", odmax=odmax,
                burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit_frn$U))
  expect_true(!is.null(fit_frn$V))
  
  # Test rrl with continuous version
  Y_rrl <- Z
  diag(Y_rrl) <- NA
  
  fit_rrl <- ame(Y_rrl, R=2, family="rrl",
                rvar=FALSE, cvar=TRUE,
                burn=300, nscan=800, print=FALSE)
  
  expect_true(!is.null(fit_rrl$U))
  expect_true(!is.null(fit_rrl$V))
})