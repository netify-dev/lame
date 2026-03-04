# covariate correctness across modes

# mcmc settings
burn <- 50
nscan <- 100
odens <- 1

####
# ame() covariate correctness
####

####
# Xdyad only
####

test_that("ame() unipartite: Xdyad only", {
  n <- 15
  dat <- simulate_test_network(
    n = n, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "unipartite", seed = 1001
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 1 Xdyad = 2 columns
  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("dyad", beta_names, ignore.case = TRUE)))
})

test_that("ame() bipartite: Xdyad only", {
  nA <- 12; nB <- 8
  dat <- simulate_test_network(
    n = nA, nA = nA, nB = nB, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "bipartite", seed = 1002
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
    mode = "bipartite", burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
})

####
# Xrow only
####

test_that("ame() unipartite: Xrow only", {
  n <- 15
  set.seed(1003)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  Xrow <- matrix(rnorm(n), n, 1)
  colnames(Xrow) <- "xr1"
  rownames(Y) <- rownames(Xrow) <- paste0("a", 1:n)
  colnames(Y) <- paste0("a", 1:n)

  fit <- ame(
    Y, Xrow = Xrow, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 1 Xrow = 2 columns
  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("row", beta_names, ignore.case = TRUE)))
})

test_that("ame() bipartite: Xrow only", {
  nA <- 12; nB <- 8
  set.seed(1004)
  Y <- matrix(rnorm(nA * nB), nA, nB)
  Xrow <- matrix(rnorm(nA), nA, 1)
  colnames(Xrow) <- "xr1"
  rownames(Y) <- rownames(Xrow) <- paste0("r", 1:nA)
  colnames(Y) <- paste0("c", 1:nB)

  fit <- ame(
    Y, Xrow = Xrow, R = 0, family = "normal",
    mode = "bipartite", burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("row", beta_names, ignore.case = TRUE)))
})

####
# Xcol only
####

test_that("ame() unipartite: Xcol only", {
  n <- 15
  set.seed(1005)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  Xcol <- matrix(rnorm(n), n, 1)
  colnames(Xcol) <- "xc1"
  rownames(Y) <- colnames(Y) <- paste0("a", 1:n)
  rownames(Xcol) <- paste0("a", 1:n)

  fit <- ame(
    Y, Xcol = Xcol, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("col", beta_names, ignore.case = TRUE)))
})

test_that("ame() bipartite: Xcol only", {
  nA <- 12; nB <- 8
  set.seed(1006)
  Y <- matrix(rnorm(nA * nB), nA, nB)
  Xcol <- matrix(rnorm(nB), nB, 1)
  colnames(Xcol) <- "xc1"
  rownames(Y) <- paste0("r", 1:nA)
  colnames(Y) <- rownames(Xcol) <- paste0("c", 1:nB)

  fit <- ame(
    Y, Xcol = Xcol, R = 0, family = "normal",
    mode = "bipartite", burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("col", beta_names, ignore.case = TRUE)))
})

####
# Xdyad + Xrow + Xcol
####

test_that("ame() unipartite: Xdyad + Xrow + Xcol", {
  n <- 15
  set.seed(1007)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  Xdyad <- array(rnorm(n * n), dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  Xrow <- matrix(rnorm(n), n, 1)
  Xcol <- matrix(rnorm(n), n, 1)
  colnames(Xrow) <- "xr1"
  colnames(Xcol) <- "xc1"
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors
  dimnames(Xdyad) <- list(actors, actors, "xd1")
  rownames(Xrow) <- rownames(Xcol) <- actors

  fit <- ame(
    Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
    R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 1 Xrow + 1 Xcol + 1 Xdyad = 4 columns
  expect_equal(ncol(fit$BETA), 4)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("row", beta_names, ignore.case = TRUE)))
  expect_true(any(grepl("col", beta_names, ignore.case = TRUE)))
  expect_true(any(grepl("dyad", beta_names, ignore.case = TRUE)))
})

test_that("ame() bipartite: Xdyad + Xrow + Xcol", {
  nA <- 12; nB <- 8
  set.seed(1008)
  Y <- matrix(rnorm(nA * nB), nA, nB)
  Xdyad <- array(rnorm(nA * nB), dim = c(nA, nB, 1))
  Xrow <- matrix(rnorm(nA), nA, 1)
  Xcol <- matrix(rnorm(nB), nB, 1)
  colnames(Xrow) <- "xr1"
  colnames(Xcol) <- "xc1"
  row_actors <- paste0("r", 1:nA)
  col_actors <- paste0("c", 1:nB)
  rownames(Y) <- row_actors; colnames(Y) <- col_actors
  dimnames(Xdyad) <- list(row_actors, col_actors, "xd1")
  rownames(Xrow) <- row_actors
  rownames(Xcol) <- col_actors

  fit <- ame(
    Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
    R = 0, family = "normal", mode = "bipartite",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 4)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
})

####
# No covariates (intercept only)
####

test_that("ame() unipartite: no covariates (intercept only)", {
  n <- 15
  set.seed(1009)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  rownames(Y) <- colnames(Y) <- paste0("a", 1:n)

  fit <- ame(
    Y, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept only = 1 column
  expect_equal(ncol(fit$BETA), 1)
  expect_true("intercept" %in% colnames(fit$BETA))
})

test_that("ame() bipartite: no covariates (intercept only)", {
  nA <- 12; nB <- 8
  set.seed(1010)
  Y <- matrix(rnorm(nA * nB), nA, nB)
  rownames(Y) <- paste0("r", 1:nA)
  colnames(Y) <- paste0("c", 1:nB)

  fit <- ame(
    Y, R = 0, family = "normal", mode = "bipartite",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 1)
  expect_true("intercept" %in% colnames(fit$BETA))
})

####
# intercept=FALSE
####

test_that("ame() unipartite: intercept=FALSE with Xdyad", {
  n <- 15
  dat <- simulate_test_network(
    n = n, family = "normal", R = 0,
    beta_intercept = 0, beta_dyad = 0.5,
    mode = "unipartite", seed = 1011
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
    intercept = FALSE,
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # No intercept, 1 Xdyad = 1 column
  expect_equal(ncol(fit$BETA), 1)
  beta_names <- colnames(fit$BETA)
  expect_false("intercept" %in% beta_names)
})

####
# Coefficient recovery for strong signal
####

test_that("ame() recovers strong Xdyad signal", {
  n <- 20
  dat <- simulate_test_network(
    n = n, family = "normal", R = 0,
    beta_intercept = 2, beta_dyad = 3.0,
    mode = "unipartite", seed = 1012
  )

  fit <- ame(
    dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
    burn = 200, nscan = 500, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  beta_mean <- colMeans(fit$BETA)
  # Intercept should be in neighborhood of 2
  expect_true(abs(beta_mean["intercept"] - 2) < 2)
  # Dyad coefficient should be positive and in neighborhood of 3
  dyad_col <- grep("dyad", names(beta_mean), ignore.case = TRUE)
  expect_true(length(dyad_col) > 0)
  expect_true(beta_mean[dyad_col] > 0)
  expect_true(abs(beta_mean[dyad_col] - 3) < 3)
})

####
# lame() covariate correctness
####

n_time <- 3

####
# lame() Xdyad only
####

test_that("lame() unipartite: Xdyad only", {
  n <- 15
  dat <- simulate_test_network(
    n = n, n_time = n_time, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "unipartite", seed = 1101
  )

  fit <- lame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("dyad", beta_names, ignore.case = TRUE)))
})

test_that("lame() bipartite: Xdyad only", {
  nA <- 12; nB <- 8
  dat <- simulate_test_network(
    n = nA, nA = nA, nB = nB, n_time = n_time, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "bipartite", seed = 1102
  )

  # lame() bipartite expects Xdyad as list of matrices
  Xdyad_list <- lapply(dat$Xdyad, function(x) x[, , 1])

  fit <- lame(
    dat$Y, Xdyad = Xdyad_list, R = 2, family = "normal",
    mode = "bipartite", burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
})

####
# lame() Xrow only
####

test_that("lame() unipartite: Xrow only", {
  n <- 15
  set.seed(1103)
  actors <- paste0("a", 1:n)

  Y_list <- lapply(1:n_time, function(t) {
    Y <- matrix(rnorm(n * n), n, n)
    diag(Y) <- NA
    rownames(Y) <- colnames(Y) <- actors
    Y
  })
  names(Y_list) <- paste0("t", 1:n_time)

  Xrow_list <- lapply(1:n_time, function(t) {
    xr <- matrix(rnorm(n), n, 1)
    colnames(xr) <- "xr1"
    rownames(xr) <- actors
    xr
  })
  names(Xrow_list) <- paste0("t", 1:n_time)

  fit <- lame(
    Y_list, Xrow = Xrow_list, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("row", beta_names, ignore.case = TRUE)))
})

####
# lame() Xcol only
####

test_that("lame() unipartite: Xcol only", {
  n <- 15
  set.seed(1104)
  actors <- paste0("a", 1:n)

  Y_list <- lapply(1:n_time, function(t) {
    Y <- matrix(rnorm(n * n), n, n)
    diag(Y) <- NA
    rownames(Y) <- colnames(Y) <- actors
    Y
  })
  names(Y_list) <- paste0("t", 1:n_time)

  Xcol_list <- lapply(1:n_time, function(t) {
    xc <- matrix(rnorm(n), n, 1)
    colnames(xc) <- "xc1"
    rownames(xc) <- actors
    xc
  })
  names(Xcol_list) <- paste0("t", 1:n_time)

  fit <- lame(
    Y_list, Xcol = Xcol_list, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 2)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("col", beta_names, ignore.case = TRUE)))
})

####
# lame() Xdyad + Xrow + Xcol
####

test_that("lame() unipartite: Xdyad + Xrow + Xcol", {
  n <- 15
  set.seed(1105)
  actors <- paste0("a", 1:n)

  Y_list <- lapply(1:n_time, function(t) {
    Y <- matrix(rnorm(n * n), n, n)
    diag(Y) <- NA
    rownames(Y) <- colnames(Y) <- actors
    Y
  })
  names(Y_list) <- paste0("t", 1:n_time)

  Xdyad_list <- lapply(1:n_time, function(t) {
    X <- array(rnorm(n * n), dim = c(n, n, 1))
    diag(X[,,1]) <- NA
    dimnames(X) <- list(actors, actors, "xd1")
    X
  })
  names(Xdyad_list) <- paste0("t", 1:n_time)

  Xrow_list <- lapply(1:n_time, function(t) {
    xr <- matrix(rnorm(n), n, 1)
    colnames(xr) <- "xr1"
    rownames(xr) <- actors
    xr
  })
  names(Xrow_list) <- paste0("t", 1:n_time)

  Xcol_list <- lapply(1:n_time, function(t) {
    xc <- matrix(rnorm(n), n, 1)
    colnames(xc) <- "xc1"
    rownames(xc) <- actors
    xc
  })
  names(Xcol_list) <- paste0("t", 1:n_time)

  fit <- lame(
    Y_list, Xdyad = Xdyad_list, Xrow = Xrow_list, Xcol = Xcol_list,
    R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 1 Xrow + 1 Xcol + 1 Xdyad = 4
  expect_equal(ncol(fit$BETA), 4)
  beta_names <- colnames(fit$BETA)
  expect_true("intercept" %in% beta_names)
  expect_true(any(grepl("row", beta_names, ignore.case = TRUE)))
  expect_true(any(grepl("col", beta_names, ignore.case = TRUE)))
  expect_true(any(grepl("dyad", beta_names, ignore.case = TRUE)))
})

####
# lame() no covariates (intercept only)
####

test_that("lame() unipartite: no covariates (intercept only)", {
  n <- 15
  set.seed(1106)
  actors <- paste0("a", 1:n)

  Y_list <- lapply(1:n_time, function(t) {
    Y <- matrix(rnorm(n * n), n, n)
    diag(Y) <- NA
    rownames(Y) <- colnames(Y) <- actors
    Y
  })
  names(Y_list) <- paste0("t", 1:n_time)

  fit <- lame(
    Y_list, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 1)
  expect_true("intercept" %in% colnames(fit$BETA))
})

####
# lame() intercept=FALSE
####

test_that("lame() unipartite: intercept=FALSE with Xdyad", {
  n <- 15
  dat <- simulate_test_network(
    n = n, n_time = n_time, family = "normal", R = 0,
    beta_intercept = 0, beta_dyad = 0.5,
    mode = "unipartite", seed = 1107
  )

  fit <- lame(
    dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
    intercept = FALSE,
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 1)
  beta_names <- colnames(fit$BETA)
  expect_false("intercept" %in% beta_names)
})

####
# covariate edge cases
####

####
# Multiple dyadic covariates
####

test_that("ame() unipartite: multiple Xdyad covariates", {
  n <- 15
  set.seed(1201)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  # 3 dyadic covariates
  Xdyad <- array(rnorm(n * n * 3), dim = c(n, n, 3))
  for (k in 1:3) diag(Xdyad[,,k]) <- NA
  dimnames(Xdyad) <- list(actors, actors, paste0("xd", 1:3))

  fit <- ame(
    Y, Xdyad = Xdyad, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 3 Xdyad = 4 columns
  expect_equal(ncol(fit$BETA), 4)
})

####
# Multiple Xrow + Xcol covariates
####

test_that("ame() unipartite: multiple Xrow and Xcol", {
  n <- 15
  set.seed(1202)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  Xrow <- matrix(rnorm(n * 2), n, 2)
  colnames(Xrow) <- c("xr1", "xr2")
  rownames(Xrow) <- actors

  Xcol <- matrix(rnorm(n * 2), n, 2)
  colnames(Xcol) <- c("xc1", "xc2")
  rownames(Xcol) <- actors

  fit <- ame(
    Y, Xrow = Xrow, Xcol = Xcol, R = 0, family = "normal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # intercept + 2 Xrow + 2 Xcol = 5 columns
  expect_equal(ncol(fit$BETA), 5)
  beta_names <- colnames(fit$BETA)
  expect_true(sum(grepl("row", beta_names, ignore.case = TRUE)) == 2)
  expect_true(sum(grepl("col", beta_names, ignore.case = TRUE)) == 2)
})

####
# Covariates with missing values
####

test_that("ame() handles NAs in Xdyad covariates", {
  n <- 15
  set.seed(1203)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  Xdyad <- array(rnorm(n * n), dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  # Add 10% additional NAs
  na_idx <- sample(which(!is.na(Xdyad[,,1])), size = floor(0.1 * n * (n - 1)))
  Xdyad[,,1][na_idx] <- NA
  dimnames(Xdyad) <- list(actors, actors, "xd1")

  # Should run with a warning about NAs in design matrix
  expect_no_error(
    fit <- suppressWarnings(ame(
      Y, Xdyad = Xdyad, R = 0, family = "normal",
      burn = burn, nscan = nscan, odens = odens,
      print = FALSE, gof = FALSE, seed = 42
    ))
  )

  expect_s3_class(fit, "ame")
  expect_equal(ncol(fit$BETA), 2)
})

####
# Covariates with zero variance
####

test_that("ame() handles constant covariate column", {
  n <- 15
  set.seed(1204)
  Y <- matrix(rnorm(n * n), n, n)
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  # Constant Xdyad (zero variance)
  Xdyad <- array(1, dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  dimnames(Xdyad) <- list(actors, actors, "xd_const")

  # Constant covariate is collinear with intercept
  # design_array should handle this (it checks var==0 and may not add intercept)
  expect_no_error(
    fit <- suppressWarnings(ame(
      Y, Xdyad = Xdyad, R = 0, family = "normal",
      burn = burn, nscan = nscan, odens = odens,
      print = FALSE, gof = FALSE, seed = 42
    ))
  )

  expect_s3_class(fit, "ame")
  # When a covariate is constant (var==0), design_array skips adding
  # intercept since it would be collinear. Check that BETA exists.
  expect_true(nrow(fit$BETA) > 0)
})

####
# Binary family with covariates
####

test_that("ame() binary with Xdyad + Xrow", {
  n <- 15
  set.seed(1205)
  Z <- matrix(rnorm(n * n, 0.5, 1), n, n)
  Y <- ifelse(Z > 0, 1, 0)
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  Xdyad <- array(rnorm(n * n), dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  dimnames(Xdyad) <- list(actors, actors, "xd1")

  Xrow <- matrix(rnorm(n), n, 1)
  colnames(Xrow) <- "xr1"
  rownames(Xrow) <- actors

  fit <- ame(
    Y, Xdyad = Xdyad, Xrow = Xrow, R = 0, family = "binary",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  expect_equal(ncol(fit$BETA), 3)  # intercept + Xrow + Xdyad
  expect_false(any(is.nan(fit$BETA)))
  expect_false(any(is.infinite(fit$BETA)))
})

####
# ordinal family: intercept=FALSE by default
####

test_that("ame() ordinal has no intercept by default", {
  n <- 15
  set.seed(1206)
  Z <- matrix(rnorm(n * n), n, n)
  breaks <- quantile(Z, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  Y <- matrix(0L, n, n)
  Y[Z > breaks[1]] <- 1L
  Y[Z > breaks[2]] <- 2L
  Y[Z > breaks[3]] <- 3L
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  Xdyad <- array(rnorm(n * n), dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  dimnames(Xdyad) <- list(actors, actors, "xd1")

  fit <- ame(
    Y, Xdyad = Xdyad, R = 0, family = "ordinal",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # ordinal: intercept=FALSE by default, so only 1 Xdyad col
  expect_equal(ncol(fit$BETA), 1)
  expect_false("intercept" %in% colnames(fit$BETA))
})

####
# rrl family: intercept=FALSE by default
####

test_that("ame() rrl has no intercept by default", {
  n <- 15
  set.seed(1207)
  Z <- matrix(rnorm(n * n), n, n)
  diag(Z) <- -Inf
  Y <- matrix(0, n, n)
  for (i in 1:n) {
    zi <- Z[i, ]
    valid <- which(is.finite(zi))
    top <- valid[order(zi[valid], decreasing = TRUE)[1:min(3, length(valid))]]
    Y[i, top] <- seq(length(top), 1)
  }
  diag(Y) <- NA
  actors <- paste0("a", 1:n)
  rownames(Y) <- colnames(Y) <- actors

  Xdyad <- array(rnorm(n * n), dim = c(n, n, 1))
  diag(Xdyad[,,1]) <- NA
  dimnames(Xdyad) <- list(actors, actors, "xd1")

  fit <- ame(
    Y, Xdyad = Xdyad, R = 0, family = "rrl",
    burn = burn, nscan = nscan, odens = odens,
    print = FALSE, gof = FALSE, seed = 42
  )

  # rrl: intercept=FALSE, rvar=FALSE by default
  expect_equal(ncol(fit$BETA), 1)
  expect_false("intercept" %in% colnames(fit$BETA))
})

####
# lame() coefficient recovery
####

test_that("lame() recovers strong Xdyad signal", {
  n <- 20
  dat <- simulate_test_network(
    n = n, n_time = 3, family = "normal", R = 0,
    beta_intercept = 2, beta_dyad = 3.0,
    mode = "unipartite", seed = 1208
  )

  fit <- lame(
    dat$Y, Xdyad = dat$Xdyad, R = 0, family = "normal",
    burn = 200, nscan = 500, odens = 1,
    print = FALSE, gof = FALSE, seed = 42
  )

  beta_mean <- colMeans(fit$BETA)
  # Intercept should be near 2
  expect_true(abs(beta_mean["intercept"] - 2) < 2)
  # Dyad coefficient should be positive
  dyad_col <- grep("dyad", names(beta_mean), ignore.case = TRUE)
  expect_true(length(dyad_col) > 0)
  expect_true(beta_mean[dyad_col] > 0)
})
