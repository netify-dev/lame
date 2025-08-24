test_that("bipartite EZ computation works correctly", {
  skip_if_not(requireNamespace("lame", quietly = TRUE))
  # Set dimensions
  nA <- 5
  nB <- 7
  T <- 3
  RA <- 2
  RB <- 3
  
  # Create test data
  base <- array(rnorm(nA * nB * T), c(nA, nB, T))
  a <- matrix(rnorm(nA * T), nA, T)
  b <- matrix(rnorm(nB * T), nB, T)
  U <- array(rnorm(nA * RA * T), c(nA, RA, T))
  V <- array(rnorm(nB * RB * T), c(nB, RB, T))
  G <- matrix(rnorm(RA * RB), RA, RB)
  
  # Test EZ computation
  EZ <- lame:::get_EZ_bip_cpp(base, a, b, U, V, G)
  
  # Check dimensions
  expect_equal(dim(EZ), c(nA, nB, T))
  
  # Check that values are finite
  expect_true(all(is.finite(EZ)))
  
  # Manual check for first time slice
  EZ_manual_t1 <- base[,,1]
  for(i in 1:nA) {
    EZ_manual_t1[i,] <- EZ_manual_t1[i,] + a[i,1]
  }
  for(j in 1:nB) {
    EZ_manual_t1[,j] <- EZ_manual_t1[,j] + b[j,1]
  }
  EZ_manual_t1 <- EZ_manual_t1 + U[,,1] %*% G %*% t(V[,,1])
  
  expect_equal(EZ[,,1], EZ_manual_t1, tolerance = 1e-10)
})

test_that("four-cycle counting works for bipartite networks", {
  # Create a simple bipartite network
  Y <- array(0, c(3, 3, 1))
  
  # Create a pattern with exactly 1 four-cycle
  # Nodes A1, A2 connect to B1, B2 (forms a 4-cycle)
  Y[1, 1, 1] <- 1  # A1 -> B1
  Y[1, 2, 1] <- 1  # A1 -> B2
  Y[2, 1, 1] <- 1  # A2 -> B1
  Y[2, 2, 1] <- 1  # A2 -> B2
  
  cycles <- lame:::count_four_cycles_bip_cpp(Y)
  expect_equal(as.numeric(cycles), 1)
  
  # Test with no edges (should be 0 cycles)
  Y_empty <- array(0, c(5, 5, 2))
  cycles_empty <- lame:::count_four_cycles_bip_cpp(Y_empty)
  expect_equal(as.numeric(cycles_empty), c(0, 0))
})

test_that("bipartite degree computation works", {
  # Create test network
  Y <- array(0, c(3, 4, 2))
  Y[1, 1:2, 1] <- 1  # Node A1 has degree 2 at time 1
  Y[2, 3:4, 1] <- 1  # Node A2 has degree 2 at time 1
  Y[3, 1, 1] <- 1    # Node A3 has degree 1 at time 1
  
  degrees <- lame:::compute_degrees_bip_cpp(Y)
  
  expect_equal(degrees$row_degrees[,1], c(2, 2, 1))
  expect_equal(degrees$col_degrees[,1], c(2, 1, 1, 1))
})

test_that("UV sampling for bipartite works", {
  # Setup
  nA <- 4
  nB <- 5
  T <- 2
  RA <- 2
  RB <- 3
  
  R <- array(rnorm(nA * nB * T), c(nA, nB, T))
  V <- array(rnorm(nB * RB * T), c(nB, RB, T))
  G <- matrix(rnorm(RA * RB), RA, RB)
  lambdaU <- rep(1, RA)
  s2 <- 1
  
  # Test U sampling
  U_new <- lame:::sample_U_bip_cpp(R, V, G, lambdaU, s2)
  
  expect_equal(dim(U_new), c(nA, RA, T))
  expect_true(all(is.finite(U_new)))
  
  # Test V sampling
  lambdaV <- rep(1, RB)
  V_new <- lame:::sample_V_bip_cpp(R, U_new, G, lambdaV, s2)
  
  expect_equal(dim(V_new), c(nB, RB, T))
  expect_true(all(is.finite(V_new)))
})

test_that("G sampling and orientation works", {
  nA <- 4
  nB <- 5
  T <- 2
  RA <- 3
  RB <- 2
  
  R <- array(rnorm(nA * nB * T), c(nA, nB, T))
  U <- array(rnorm(nA * RA * T), c(nA, RA, T))
  V <- array(rnorm(nB * RB * T), c(nB, RB, T))
  lambdaG <- 1
  s2 <- 1
  
  # Test G sampling
  G_new <- lame:::sample_G_bip_cpp(R, U, V, lambdaG, s2)
  
  expect_equal(dim(G_new), c(RA, RB))
  expect_true(all(is.finite(G_new)))
  
  # Test canonical orientation
  oriented <- lame:::canon_orient_bip_cpp(U, V, G_new)
  
  expect_equal(dim(oriented$U), dim(U))
  expect_equal(dim(oriented$V), dim(V))
  expect_equal(dim(oriented$G), dim(G_new))
  
  # Check that G is now diagonal-like (singular values on diagonal)
  G_oriented <- oriented$G
  # Off-diagonal elements should be close to zero
  for(i in 1:min(RA, RB)) {
    for(j in 1:min(RA, RB)) {
      if(i != j) {
        expect_lt(abs(G_oriented[i,j]), 1e-10)
      }
    }
  }
})

test_that("bipartite initialization works", {
  skip_if_not_installed("lame")
  
  # Create simple bipartite data
  Y <- array(0, c(5, 7, 2))
  Y[,,1] <- matrix(rbinom(5*7, 1, 0.3), 5, 7)
  Y[,,2] <- matrix(rbinom(5*7, 1, 0.3), 5, 7)
  
  init <- lame:::init_bipartite_startvals(Y, family = "binary", 
                                   nA = 5, nB = 7, 
                                   RA = 2, RB = 3, T = 2)
  
  expect_equal(dim(init$Z), c(5, 7, 2))
  expect_equal(dim(init$a), c(5, 2))
  expect_equal(dim(init$b), c(7, 2))
  expect_equal(dim(init$U), c(5, 2, 2))
  expect_equal(dim(init$V), c(7, 3, 2))
  expect_equal(dim(init$G), c(2, 3))
  expect_null(init$rho)
})