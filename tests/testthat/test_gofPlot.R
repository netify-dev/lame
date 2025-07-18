test_that("gofPlot works", {
  
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
  
  # lame
  fit <- lame(
    Y = Y, 
    Xdyad = X,
    family = "bin",
    burn = 2,     
    nscan = 4,     
    odens = 2,     
    print = FALSE, 
    plot = FALSE   
  )
  
  expect_no_error({
    p <- gofPlot(fit$GOF, symmetric = FALSE)
  })
})
