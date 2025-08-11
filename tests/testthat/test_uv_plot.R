test_that("uv_plot works", {

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
    R = 2,
    family = "bin",
    burn = 2,     
    nscan = 4,     
    odens = 2,     
    print = FALSE, 
    plot = FALSE   
  )

  expect_no_error({
    p <- uv_plot(fit)
  })
  
  # Test with explicit Y matrix
  expect_no_error({
    p <- uv_plot(fit, Y = Y[[1]])
  })
  
  # Test with just U and V
  expect_no_error({
    p <- uv_plot(Y = Y[[1]], U = fit$U, V = fit$V)
  })
})