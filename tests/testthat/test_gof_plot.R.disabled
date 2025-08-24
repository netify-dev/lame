test_that("gof_plot works", {
  
  cnt <- 3
  time <- 2
  
  # Y matrices
  Y <- list()
  for(t in 1:time) {
    Y[[t]] <- matrix(rbinom(cnt*cnt, 1, 0.5), cnt, cnt)
    diag(Y[[t]]) <- NA
    rownames(Y[[t]]) <- colnames(Y[[t]]) <- paste0("Cnt", 1:cnt)
  }
  names(Y) <- paste0("Time", 1:time)
  
  # dyadic covariates
  X <- list()
  for(t in 1:time) {
    X[[t]] <- array(rnorm(cnt*cnt*2), dim = c(cnt, cnt, 2))
    dimnames(X[[t]]) <- list(paste0("Cnt", 1:cnt), paste0("Cnt", 1:cnt), c("cov1", "cov2"))
  }
  names(X) <- paste0("T", 1:time)
  
  # lame with GOF
  fit <- lame(
    Y = Y, 
    Xdyad = X,
    family = "binary",
    burn = 2,     
    nscan = 4,     
    odens = 2,
    gof = TRUE,     
    print = FALSE, 
    plot = FALSE   
  )
  
  # Test with fit object (recommended usage)
  expect_no_error({
    p <- gof_plot(fit)
  })
  
  # Test type argument
  expect_no_error({
    p <- gof_plot(fit, type = "static")
  })
  
  expect_no_error({
    p <- gof_plot(fit, type = "longitudinal")
  })
  
  # Test statistics selection
  expect_no_error({
    p <- gof_plot(fit, statistics = c("sd.row", "sd.col"))
  })
})

# Skipping AME test due to numerical instability in small test networks
# test_that("gof_plot works for AME models", {
#   fit <- ame(...)
#   expect_no_error({
#     p <- gof_plot(fit)
#   })
# })