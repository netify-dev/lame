# edge cases and input validation

####
# missing data
####

test_that("ame() handles 10% random NAs in Y", {
  n <- 15
  set.seed(6001)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  # Add 10% random NAs
  off_diag <- which(!is.na(Y))
  na_idx <- sample(off_diag, size = floor(0.1 * length(off_diag)))
  Y[na_idx] <- NA
  rownames(Y) <- colnames(Y) <- paste0("a", 1:n)

  fit <- ame(
    Y, R = 0, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  # YPM should have values where Y had NAs (imputed)
  expect_false(all(is.na(fit$YPM[na_idx])))
})

test_that("ame() handles 50% random NAs in Y", {
  n <- 15
  set.seed(6002)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  off_diag <- which(!is.na(Y))
  na_idx <- sample(off_diag, size = floor(0.5 * length(off_diag)))
  Y[na_idx] <- NA
  rownames(Y) <- colnames(Y) <- paste0("a", 1:n)

  fit <- ame(
    Y, R = 0, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  expect_false(any(is.nan(fit$BETA)))
})

####
# extreme dimensions
####

test_that("ame() works with R=0 (no latent factors)", {
  n <- 15
  dat <- simulate_test_network(
    n = n, family = "normal", R = 0,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "unipartite", seed = 6010
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  # U and V should be NULL or degenerate
  expect_true(is.null(fit$U) || nrow(fit$U) == 0 || ncol(fit$U) == 0)
})

test_that("ame() works with R=1", {
  n <- 15
  dat <- simulate_test_network(
    n = n, family = "normal", R = 1,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "unipartite", seed = 6011
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 1, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  expect_equal(ncol(fit$U), 1)
})

test_that("ame() works with small network (n=5)", {
  n <- 5
  dat <- simulate_test_network(
    n = n, family = "normal", R = 1,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "unipartite", seed = 6012
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 1, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  expect_equal(length(fit$APM), n)
})

test_that("lame() works with N=1 (single time period)", {
  n <- 15
  set.seed(6013)
  actors <- paste0("a", 1:n)
  Y_list <- list(t1 = {
    Y <- matrix(rnorm(n * n), n, n)
    diag(Y) <- NA
    rownames(Y) <- colnames(Y) <- actors
    Y
  })

  fit <- lame(
    Y_list, R = 0, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "lame")
  expect_equal(fit$n_time, 1)
})

test_that("ame() works with extreme bipartite imbalance (nA=4, nB=20)", {
  nA <- 4; nB <- 20
  set.seed(6014)
  Y <- matrix(rnorm(nA * nB), nA, nB)
  rownames(Y) <- paste0("r", 1:nA)
  colnames(Y) <- paste0("c", 1:nB)

  fit <- ame(
    Y, R = 1, family = "normal", mode = "bipartite",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "ame")
  expect_equal(length(fit$APM), nA)
  expect_equal(length(fit$BPM), nB)
})

####
# changing actor compositions (lame only)
####

test_that("lame() handles actors entering/leaving across time", {
  # Period 1: actors A,B,C,D
  # Period 2: actors A,B,C,E (D leaves, E joins)
  # Period 3: actors A,C,E,F (B leaves, F joins)
  set.seed(6020)
  Y1 <- matrix(rnorm(16), 4, 4)
  diag(Y1) <- NA
  rownames(Y1) <- colnames(Y1) <- c("A", "B", "C", "D")

  Y2 <- matrix(rnorm(16), 4, 4)
  diag(Y2) <- NA
  rownames(Y2) <- colnames(Y2) <- c("A", "B", "C", "E")

  Y3 <- matrix(rnorm(16), 4, 4)
  diag(Y3) <- NA
  rownames(Y3) <- colnames(Y3) <- c("A", "C", "E", "F")

  Y_list <- list(t1 = Y1, t2 = Y2, t3 = Y3)

  fit <- lame(
    Y_list, R = 0, family = "normal",
    burn = 50, nscan = 100, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_s3_class(fit, "lame")
  # Full union: A,B,C,D,E,F = 6 actors
  expect_equal(length(fit$APM), 6)
})

####
# input validation
####

test_that("ame() errors on non-square Y for unipartite", {
  Y <- matrix(rnorm(12), 3, 4)
  expect_error(ame(Y, family = "normal", print = FALSE))
})

test_that("lame() errors on non-list Y", {
  Y <- matrix(rnorm(25), 5, 5)
  expect_error(lame(Y, family = "normal", print = FALSE))
})

test_that("ame() errors on unsupported family", {
  n <- 5
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  rownames(Y) <- colnames(Y) <- paste0("a", 1:n)
  expect_error(ame(Y, family = "unknown", print = FALSE))
})
