# dynamic effects smoke tests

n_time <- 3

####
# unipartite dynamic effects
####

for (fam in c("normal", "binary")) {
  # dynamic_ab only
  test_that(paste0("lame() unipartite dynamic_ab: ", fam), {
    n <- 15
    dat <- simulate_test_network(
      n = n, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_ab = 0.7, sigma_ab = 0.3,
      mode = "unipartite", seed = 600
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      dynamic_ab = TRUE, dynamic_uv = FALSE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_true(isTRUE(fit$dynamic_ab))
    expect_false(isTRUE(fit$dynamic_uv))

    # a_dynamic should be n × T matrix
    if (!is.null(fit$a_dynamic)) {
      expect_equal(nrow(fit$a_dynamic), n)
      expect_equal(ncol(fit$a_dynamic), n_time)
    }
    # rho_ab should be within (-1, 1)
    if (!is.null(fit$rho_ab)) {
      expect_true(all(fit$rho_ab > -1 & fit$rho_ab < 1))
    }
    # APM still a vector of length n
    expect_equal(length(fit$APM), n)
  })

  # dynamic_uv only
  test_that(paste0("lame() unipartite dynamic_uv: ", fam), {
    n <- 15
    dat <- simulate_test_network(
      n = n, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_uv = 0.7, sigma_uv = 0.3,
      mode = "unipartite", seed = 610
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      dynamic_ab = FALSE, dynamic_uv = TRUE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_false(isTRUE(fit$dynamic_ab))
    expect_true(isTRUE(fit$dynamic_uv))

    # U should be 3D: n × R × T
    if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
      expect_equal(dim(fit$U)[1], n)
      expect_equal(dim(fit$U)[3], n_time)
    }
    # rho_uv should be within (-1, 1)
    if (!is.null(fit$rho_uv)) {
      expect_true(all(fit$rho_uv > -1 & fit$rho_uv < 1))
    }
  })

  # both dynamic
  test_that(paste0("lame() unipartite dynamic_ab + dynamic_uv: ", fam), {
    n <- 15
    dat <- simulate_test_network(
      n = n, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_ab = 0.7, rho_uv = 0.5,
      mode = "unipartite", seed = 620
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      dynamic_ab = TRUE, dynamic_uv = TRUE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_true(isTRUE(fit$dynamic_ab))
    expect_true(isTRUE(fit$dynamic_uv))
  })
}

####
# Unipartite symmetric dynamic
####

test_that("lame() unipartite symmetric dynamic_ab: normal", {
  n <- 15
  dat <- simulate_test_network(
    n = n, n_time = n_time, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.3,
    rho_ab = 0.7, symmetric = TRUE,
    mode = "unipartite", seed = 630
  )

  # Ensure symmetric
  Y_sym <- lapply(dat$Y, function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  })
  names(Y_sym) <- names(dat$Y)

  fit <- lame(
    Y_sym, Xdyad = dat$Xdyad, R = 2, family = "normal",
    symmetric = TRUE, dynamic_ab = TRUE,
    burn = 5, nscan = 10, odens = 1,
    verbose = FALSE, gof = TRUE, seed = 42
  )

  expect_s3_class(fit, "lame")
  expect_true(isTRUE(fit$dynamic_ab))
})

####
# bipartite dynamic effects
####

for (fam in c("normal", "binary")) {
  test_that(paste0("lame() bipartite dynamic_ab: ", fam), {
    nA <- 12
    nB <- 8
    dat <- simulate_test_network(
      n = nA, nA = nA, nB = nB, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_ab = 0.7,
      mode = "bipartite", seed = 700
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      mode = "bipartite", dynamic_ab = TRUE, dynamic_uv = FALSE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_true(isTRUE(fit$dynamic_ab))

    # a_dynamic: nA × T, b_dynamic: nB × T
    if (!is.null(fit$a_dynamic)) {
      expect_equal(nrow(fit$a_dynamic), nA)
      expect_equal(ncol(fit$a_dynamic), n_time)
    }
    if (!is.null(fit$b_dynamic)) {
      expect_equal(nrow(fit$b_dynamic), nB)
      expect_equal(ncol(fit$b_dynamic), n_time)
    }
    # APM length nA, BPM length nB
    expect_equal(length(fit$APM), nA)
    expect_equal(length(fit$BPM), nB)
  })

  test_that(paste0("lame() bipartite dynamic_uv: ", fam), {
    nA <- 12
    nB <- 8
    dat <- simulate_test_network(
      n = nA, nA = nA, nB = nB, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_uv = 0.7,
      mode = "bipartite", seed = 710
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      mode = "bipartite", dynamic_ab = FALSE, dynamic_uv = TRUE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_true(isTRUE(fit$dynamic_uv))

    # U: nA × R × T, V: nB × R × T (different first dims)
    if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
      expect_equal(dim(fit$U)[1], nA)
      expect_equal(dim(fit$U)[3], n_time)
    }
    if (!is.null(fit$V) && length(dim(fit$V)) == 3) {
      expect_equal(dim(fit$V)[1], nB)
      expect_equal(dim(fit$V)[3], n_time)
    }
    # G should exist
    expect_true(!is.null(fit$G))
  })

  test_that(paste0("lame() bipartite dynamic_ab + dynamic_uv: ", fam), {
    nA <- 12
    nB <- 8
    dat <- simulate_test_network(
      n = nA, nA = nA, nB = nB, n_time = n_time, family = fam, R = 2,
      beta_intercept = if (fam == "binary") 0.5 else 1,
      beta_dyad = 0.3, rho_ab = 0.7, rho_uv = 0.5,
      mode = "bipartite", seed = 720
    )

    fit <- lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = fam,
      mode = "bipartite", dynamic_ab = TRUE, dynamic_uv = TRUE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = TRUE, seed = 42
    )

    expect_s3_class(fit, "lame")
    expect_true(isTRUE(fit$dynamic_ab))
    expect_true(isTRUE(fit$dynamic_uv))
  })
}

####
# dynamic_G warning
####

test_that("lame() warns about dynamic_G=TRUE not implemented", {
  n <- 10
  dat <- simulate_test_network(
    n = n, nA = n, nB = 8, n_time = 3, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.3, mode = "bipartite", seed = 800
  )

  expect_warning(
    lame(
      dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
      mode = "bipartite", dynamic_G = TRUE,
      burn = 5, nscan = 10, odens = 1,
      verbose = FALSE, gof = FALSE, seed = 42
    ),
    "dynamic_G.*not yet implemented"
  )
})
