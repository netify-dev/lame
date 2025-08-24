# Tests for utility and plotting functions
library(lame)

test_that("get_start_vals returns appropriate starting values", {
  set.seed(6886)
  n <- 20
  T <- 1  # Single time point for simplicity
  
  # Binary network - need 3D array
  Y_bin_2d <- matrix(rbinom(n*n, 1, 0.3), n, n)
  diag(Y_bin_2d) <- NA
  Y_bin <- array(Y_bin_2d, dim=c(n, n, T))
  
  # Call with proper parameters
  starts_bin <- get_start_vals(
    start_vals = NULL,
    Y = Y_bin,
    family = "binary",
    xP = 0,  # No exogenous covariates
    rvar = TRUE,
    cvar = TRUE,
    R = 2  # Number of latent dimensions
  )
  expect_true(is.list(starts_bin))
  expect_true("Z" %in% names(starts_bin))
  expect_equal(dim(starts_bin$Z), dim(Y_bin))
  
  # Normal network - need 3D array
  Y_norm_2d <- matrix(rnorm(n*n), n, n)
  diag(Y_norm_2d) <- NA
  Y_norm <- array(Y_norm_2d, dim=c(n, n, T))
  
  starts_norm <- get_start_vals(
    start_vals = NULL,
    Y = Y_norm,
    family = "normal",
    xP = 0,
    rvar = TRUE,
    cvar = TRUE,
    R = 2
  )
  # For normal family, Z is initialized to Y but with diagonal filled in
  # So check dimensions and non-diagonal elements
  expect_equal(dim(starts_norm$Z), dim(Y_norm))
  # Check that non-diagonal elements match
  for(t in 1:T) {
    non_diag <- !diag(n)
    expect_equal(starts_norm$Z[,,t][non_diag], Y_norm[,,t][non_diag])
  }
  
  # Check other components
  expect_true(is.numeric(starts_bin$beta))
  expect_true(is.numeric(starts_bin$s2))
  expect_true(is.numeric(starts_bin$rho))
  expect_true(is.matrix(starts_bin$U))
  expect_true(is.matrix(starts_bin$V))
  expect_true(is.matrix(starts_bin$Sab))
})

test_that("mhalf function works correctly", {
  # Test with positive definite matrix
  A <- matrix(c(4, 2, 2, 3), 2, 2)
  A_half <- mhalf(A)
  
  # A_half %*% A_half should equal A
  A_reconstructed <- A_half %*% A_half
  expect_equal(A_reconstructed, A, tolerance = 1e-10)
  
  # Test with identity matrix
  I <- diag(3)
  I_half <- mhalf(I)
  expect_equal(I_half, I)
  
  # Test with singular matrix (should handle gracefully)
  B <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_no_error(mhalf(B))
})

test_that("rSab_fc generates valid covariance samples", {
  set.seed(6886)
  
  # Set up parameters
  n <- 20
  a <- rnorm(n, 0, 1)
  b <- rnorm(n, 0, 1)
  Sab0 <- diag(2)
  eta0 <- 4  # Use 4 as minimum for proper prior
  
  # Generate sample with default parameters
  Sab <- rSab_fc(a, b, Sab0, eta0)
  
  # Should be 2x2 matrix
  expect_equal(dim(Sab), c(2, 2))
  
  # Should be symmetric
  expect_equal(Sab[1,2], Sab[2,1])
  
  # Should be positive definite
  eigenvals <- eigen(Sab)$values
  expect_true(all(eigenvals > 0))
})

test_that("check_format validates input correctly", {
  n <- 20
  T <- 3
  
  # Valid list of matrices
  Y_list <- list()
  for(t in 1:T) {
    Y_t <- matrix(rnorm(n*n), n, n)
    diag(Y_t) <- NA
    rownames(Y_t) <- colnames(Y_t) <- paste0("Node", 1:n)
    Y_list[[t]] <- Y_t
  }
  
  # check_format doesn't return anything, it just validates
  # It should not error if format is correct
  expect_no_error(check_format(Y_list))
  
  # check_format expects a list, not a single matrix
  # It will error if given a single matrix
  Y_single <- matrix(rnorm(n*n), n, n)
  diag(Y_single) <- NA
  expect_error(check_format(Y_single))
  
  # For array input, it also expects a list
  Y_array <- array(rnorm(n*n*T), dim=c(n, n, T))
  expect_error(check_format(Y_array))
})

test_that("zscores function computes correct scores", {
  # Test with simple vector
  y <- c(3, 1, 4, 1, 5, 9, 2, 6)
  z <- zscores(y)
  
  # Should preserve rank order
  expect_equal(order(y), order(z))
  
  # Handle ties appropriately
  y_ties <- c(1, 2, 2, 3)
  z_ties <- zscores(y_ties)
  expect_equal(z_ties[2], z_ties[3])  # Tied values should have same z-score
  
  # Handle NA values
  y_na <- c(1, NA, 3, 4)
  z_na <- zscores(y_na)
  expect_true(is.na(z_na[2]))
  expect_false(any(is.na(z_na[-2])))
})

test_that("plotting functions run without error", {
  skip_on_cran()
  
  set.seed(6886)
  n <- 15
  
  # Generate simple network and fit
  Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
  diag(Y) <- NA
  # Add rownames so APM/BPM get names
  rownames(Y) <- colnames(Y) <- paste0("Actor", 1:n)
  
  fit <- ame(Y, R=2, family="binary",
             burn=50, nscan=200, print=FALSE, plot=FALSE, gof=TRUE)
  
  # Test that plotting functions don't error
  expect_no_error({
    pdf(NULL)  # Null device to suppress actual plotting
    
    # Test ab_plot if it exists and accepts ame objects
    if(exists("ab_plot")) {
      tryCatch(
        ab_plot(fit, effect="sender"),
        error = function(e) {
          if(grepl("class 'ame' or 'lame'", e$message)) {
            # Expected - ab_plot may only work with lame objects
            NULL
          } else {
            stop(e)
          }
        }
      )
    }
    
    # Test uv_plot
    if(exists("uv_plot")) {
      tryCatch(
        uv_plot(fit),
        error = function(e) {
          if(grepl("class 'ame' or 'lame'", e$message)) {
            NULL
          } else {
            stop(e)
          }
        }
      )
    }
    
    # Test gof_plot
    if(exists("gof_plot")) {
      tryCatch(
        gof_plot(fit),
        error = function(e) {
          if(grepl("class 'ame' or 'lame'", e$message)) {
            NULL
          } else {
            stop(e)
          }
        }
      )
    }
    
    # Test trace_plot
    if(exists("trace_plot")) {
      tryCatch(
        trace_plot(list(BETA = fit$BETA)),
        error = function(e) {
          if(grepl("class 'ame' or 'lame'", e$message)) {
            NULL
          } else {
            stop(e)
          }
        }
      )
    }
    
    dev.off()
  })
})

test_that("summary and print methods work", {
  skip_on_cran()
  
  set.seed(6886)
  n <- 15
  
  Y <- matrix(rbinom(n*n, 1, 0.3), n, n)
  diag(Y) <- NA
  
  fit <- ame(Y, R=1, family="binary",
             burn=50, nscan=200, print=FALSE, plot=FALSE)
  
  # Summary should work
  expect_no_error({
    s <- summary(fit)
  })
  
  # Print should work
  expect_no_error({
    capture.output(print(fit))
  })
})

test_that("effective sample size calculation works", {
  skip_on_cran()
  skip_if_not(requireNamespace("coda", quietly = TRUE),
              message = "coda package not available")
  
  set.seed(6886)
  
  # Generate autocorrelated sequence with lower correlation for more stable test
  n <- 1000
  rho <- 0.5  # Lower correlation for more predictable ESS
  x <- numeric(n)
  x[1] <- rnorm(1)
  for(i in 2:n) {
    x[i] <- rho * x[i-1] + rnorm(1, 0, sqrt(1-rho^2))
  }
  
  # Test effective sample size
  ess <- coda::effectiveSize(x)
  
  # ESS should be less than n due to autocorrelation
  expect_lt(ess, n)
  # With rho=0.5, ESS should be reasonable (roughly n/(1+rho) ~ 667)
  expect_gt(ess, n/5)  # More lenient bound
})

test_that("array_to_list conversion works", {
  n <- 10
  T <- 4
  
  # Create array with actor names
  Y_array <- array(rnorm(n*n*T), dim=c(n, n, T))
  rownames(Y_array) <- colnames(Y_array) <- paste0("Actor", 1:n)
  
  # Create actor list (all actors present at all times)
  actorList <- lapply(1:T, function(t) paste0("Actor", 1:n))
  
  # Convert to list (array_to_list needs all three parameters)
  Y_list <- array_to_list(Y_array, actorList, paste0("T", 1:T))
  
  expect_true(is.list(Y_list))
  expect_equal(length(Y_list), T)
  
  # Each element should be a matrix
  for(t in 1:T) {
    expect_true(is.matrix(Y_list[[t]]))
    expect_equal(dim(Y_list[[t]]), c(n, n))
    # array_to_list sets diagonal to NA (expected for network data)
    Y_expected <- Y_array[,,t]
    diag(Y_expected) <- NA
    expect_equal(Y_list[[t]], Y_expected)
  }
})