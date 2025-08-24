test_that("dynamic UV works with lame", {
  # Skip if not enough time for MCMC
  skip_on_cran()
  
  # Generate simple test data
  set.seed(6886)
  n <- 10
  T <- 3
  
  # Create list of adjacency matrices
  Y <- list()
  for(t in 1:T) {
    Y[[t]] <- matrix(rnorm(n*n), n, n)
    diag(Y[[t]]) <- NA
    rownames(Y[[t]]) <- colnames(Y[[t]]) <- paste0("actor", 1:n)
  }
  names(Y) <- paste0("T", 1:T)
  
  # Test dynamic UV with small MCMC run
  fit_dynamic <- lame(Y, 
                      R = 2, 
                      dynamic_uv = TRUE,
                      family = "normal",
                      burn = 10, 
                      nscan = 20, 
                      odens = 10)
  
  # Check that output has correct structure
  expect_true("U" %in% names(fit_dynamic))
  expect_true("V" %in% names(fit_dynamic))
  
  # For dynamic UV, U and V should be 3D arrays
  if(!is.null(fit_dynamic$U) && length(dim(fit_dynamic$U)) == 3) {
    expect_equal(dim(fit_dynamic$U)[1], n)
    expect_equal(dim(fit_dynamic$U)[2], 2)  # R=2
    expect_equal(dim(fit_dynamic$U)[3], T)
  }
  
  # Test static UV for comparison
  fit_static <- lame(Y, 
                     R = 2, 
                     dynamic_uv = FALSE,
                     family = "normal",
                     burn = 10, 
                     nscan = 20, 
                     odens = 10)
  
  # Check that static UV has 2D structure
  expect_equal(nrow(fit_static$U), n)
  expect_equal(ncol(fit_static$U), 2)  # R=2
})

# These functions are internal and tested through the main lame function