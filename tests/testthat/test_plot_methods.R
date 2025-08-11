test_that("plot.lame method works", {
  
  cnt <- 10  
  time <- 5  
  
  # 
  set.seed(6886)
  
  # 
  Y <- list()
  for(t in 1:time) {
    Y[[t]] <- matrix(rbinom(cnt*cnt, 1, 0.3), cnt, cnt)  # Lower density
    diag(Y[[t]]) <- NA
    rownames(Y[[t]]) <- colnames(Y[[t]]) <- paste0("Cnt", 1:cnt)
  }
  names(Y) <- paste0("Time", 1:time)
  
  # Simpler dyadic covariates (just one to reduce complexity)
  X <- list()
  for(t in 1:time) {
    X[[t]] <- array(rnorm(cnt*cnt*1), dim = c(cnt, cnt, 1))
    dimnames(X[[t]]) <- list(paste0("Cnt", 1:cnt), paste0("Cnt", 1:cnt), c("cov1"))
  }
  names(X) <- paste0("Time", 1:time)  # Match Y names
  
  # lame with GOF - reduced R for stability in tests
  fit <- lame(
    Y = Y, 
    Xdyad = X,
    R = 0,  # No multiplicative effects for test stability
    family = "binary",
    burn = 2,     
    nscan = 4,     
    odens = 2,
    gof = TRUE,     
    print = FALSE, 
    plot = FALSE   
  )
  
  # Test default plot
  expect_no_error({
    plot(fit)
  })
  
  # Test which argument
  expect_no_error({
    plot(fit, which = c(1, 2))
  })
  
  expect_no_error({
    plot(fit, which = c("trace", "density"))
  })
})
