# Comprehensive audit tests for C++ performance optimization
# Verifies that optimized C++ functions produce correct results

test_that("rZ_nrm_batch_cpp produces valid Z and residuals", {
  set.seed(42)
  n <- 15; TT <- 5
  Z <- array(rnorm(n*n*TT), dim = c(n, n, TT))
  EZ <- array(rnorm(n*n*TT), dim = c(n, n, TT))
  Y <- array(rnorm(n*n*TT), dim = c(n, n, TT))
  # Set diagonal to NA (standard for unipartite)
  for(t in 1:TT) { diag(Y[,,t]) <- NA; diag(Z[,,t]) <- NA }
  rho <- 0.3; s2 <- 1.5

  result <- rZ_nrm_batch_cpp(Z, EZ, rho, s2, Y)

  expect_true(is.list(result))
  expect_true("Z" %in% names(result))
  expect_true("E_nrm" %in% names(result))
  expect_equal(dim(result$Z), c(n, n, TT))
  expect_equal(dim(result$E_nrm), c(n, n, TT))

  # Residuals should equal Z - EZ
  for(t in 1:TT) {
    expect_equal(result$E_nrm[,,t], result$Z[,,t] - EZ[,,t],
                 tolerance = 1e-10)
  }

  # Non-missing Y values should be preserved in Z
  for(t in 1:TT) {
    obs_idx <- which(!is.na(Y[,,t]))
    expect_equal(result$Z[,,t][obs_idx], Y[,,t][obs_idx])
  }

  # No NaN or Inf in output
  expect_false(any(is.nan(result$Z)))
  expect_false(any(is.infinite(result$Z)))
})

test_that("rZ_nrm_batch_cpp works with rho=0 (bipartite-like)", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3
  Z <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  EZ <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  Y <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  Y[sample(length(Y), 20)] <- NA

  result <- rZ_nrm_batch_cpp(Z, EZ, rho = 0, s2 = 1.0, Y)

  expect_equal(dim(result$Z), c(nA, nB, TT))
  obs_idx <- which(!is.na(Y))
  expect_equal(result$Z[obs_idx], Y[obs_idx])
  # Missing should be sampled (not NA)
  miss_idx <- which(is.na(Y))
  expect_false(any(is.na(result$Z[miss_idx])))
})

test_that("rZ_bin_bip_batch_cpp produces valid binary latent Z", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3
  Z <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  EZ <- array(rnorm(nA*nB*TT, sd = 0.5), dim = c(nA, nB, TT))
  Y <- array(rbinom(nA*nB*TT, 1, 0.5), dim = c(nA, nB, TT))
  Y[sample(length(Y), 10)] <- NA

  result <- rZ_bin_bip_batch_cpp(Z, EZ, Y)

  expect_equal(dim(result), c(nA, nB, TT))

  # Where Y=1, Z should be positive
  for(t in 1:TT) {
    idx_1 <- which(Y[,,t] == 1)
    if(length(idx_1) > 0) {
      expect_true(all(result[,,t][idx_1] > 0),
                  info = paste("Y=1 entries should have Z>0 at time", t))
    }
    # Where Y=0, Z should be negative
    idx_0 <- which(Y[,,t] == 0)
    if(length(idx_0) > 0) {
      expect_true(all(result[,,t][idx_0] < 0),
                  info = paste("Y=0 entries should have Z<0 at time", t))
    }
  }

  # No NaN or Inf
  expect_false(any(is.nan(result)))
  expect_false(any(is.infinite(result)))
})

test_that("rZ_tob_bip_batch_cpp produces valid tobit latent Z", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3
  Z <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  EZ <- array(rnorm(nA*nB*TT, mean = 1), dim = c(nA, nB, TT))
  Y <- array(pmax(rnorm(nA*nB*TT, mean = 1), 0), dim = c(nA, nB, TT))
  Y[sample(length(Y), 10)] <- NA

  result <- rZ_tob_bip_batch_cpp(Z, EZ, s2 = 1.0, Y)

  expect_equal(dim(result), c(nA, nB, TT))

  # Where Y>0, Z should equal Y
  for(t in 1:TT) {
    idx_pos <- which(Y[,,t] > 0)
    if(length(idx_pos) > 0) {
      expect_equal(result[,,t][idx_pos], Y[,,t][idx_pos])
    }
    # Where Y=0, Z should be <= 0
    idx_zero <- which(Y[,,t] == 0)
    if(length(idx_zero) > 0) {
      expect_true(all(result[,,t][idx_zero] <= 0),
                  info = paste("Y=0 entries should have Z<=0 at time", t))
    }
  }
})

test_that("compute_XtX_Xty_bip_cpp matches R implementation", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3; p <- 3

  # Create covariate list
  Xlist <- lapply(1:TT, function(t) array(rnorm(nA*nB*p), dim = c(nA, nB, p)))
  resid <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))

  # C++ version
  cpp_result <- compute_XtX_Xty_bip_cpp(Xlist, resid, p)

  # R reference implementation
  XtX_r <- matrix(0, p, p)
  Xty_r <- rep(0, p)
  for(t in 1:TT) {
    for(k1 in 1:p) {
      for(k2 in k1:p) {
        val <- sum(Xlist[[t]][,,k1] * Xlist[[t]][,,k2], na.rm = TRUE)
        XtX_r[k1, k2] <- XtX_r[k1, k2] + val
        if(k1 != k2) XtX_r[k2, k1] <- XtX_r[k2, k1] + val
      }
      Xty_r[k1] <- Xty_r[k1] + sum(Xlist[[t]][,,k1] * resid[,,t], na.rm = TRUE)
    }
  }

  expect_equal(cpp_result$XtX, XtX_r, tolerance = 1e-10)
  expect_equal(as.vector(cpp_result$Xty), Xty_r, tolerance = 1e-10)
})

test_that("compute_XtX_Xty_bip_cpp handles NAs correctly", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 2; p <- 2

  Xlist <- lapply(1:TT, function(t) array(rnorm(nA*nB*p), dim = c(nA, nB, p)))
  resid <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  # Inject NAs
  resid[1:3, 1:2, 1] <- NA

  cpp_result <- compute_XtX_Xty_bip_cpp(Xlist, resid, p)

  # R reference with na.rm
  Xty_r <- rep(0, p)
  for(t in 1:TT) {
    for(k1 in 1:p) {
      Xty_r[k1] <- Xty_r[k1] + sum(Xlist[[t]][,,k1] * resid[,,t], na.rm = TRUE)
    }
  }

  expect_equal(as.vector(cpp_result$Xty), Xty_r, tolerance = 1e-10)
})

test_that("rbeta_ab_bip_gibbs_cpp produces finite output", {
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3; p <- 2

  Z <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  Xlist <- lapply(1:TT, function(t) array(rnorm(nA*nB*p), dim = c(nA, nB, p)))
  UV_eff <- matrix(rnorm(nA*nB, sd = 0.5), nA, nB)
  a <- rnorm(nA, sd = 0.3)
  b <- rnorm(nB, sd = 0.3)

  result <- rbeta_ab_bip_gibbs_cpp(Z, Xlist, UV_eff, a, b,
                                    s2 = 1.0, g_prior = 100,
                                    va = 1.0, vb = 1.0,
                                    rvar = TRUE, cvar = TRUE)

  expect_true(is.list(result))
  expect_equal(length(result$beta), p)
  expect_equal(length(result$a), nA)
  expect_equal(length(result$b), nB)
  expect_false(any(is.nan(result$beta)))
  expect_false(any(is.nan(result$a)))
  expect_false(any(is.nan(result$b)))
  expect_false(any(is.infinite(result$beta)))
  expect_false(any(is.infinite(result$a)))
  expect_false(any(is.infinite(result$b)))
})

test_that("rUV_dynamic_bip_fc_cpp produces valid output", {
  set.seed(42)
  nA <- 10; nB <- 8; RA <- 2; RB <- 2; TT <- 5

  U_cube <- array(rnorm(nA*RA*TT, sd = 0.5), dim = c(nA, RA, TT))
  V_cube <- array(rnorm(nB*RB*TT, sd = 0.5), dim = c(nB, RB, TT))
  E <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  G <- matrix(c(1, 0.3, 0.3, 1), RA, RB)

  result <- rUV_dynamic_bip_fc_cpp(U_cube, V_cube, E, G,
                                    rho_uv = 0.9, sigma_uv = 0.1, s2 = 1.0)

  expect_true(is.list(result))
  expect_equal(dim(result$U), c(nA, RA, TT))
  expect_equal(dim(result$V), c(nB, RB, TT))
  expect_false(any(is.nan(result$U)))
  expect_false(any(is.nan(result$V)))
  expect_false(any(is.infinite(result$U)))
  expect_false(any(is.infinite(result$V)))
})

test_that("rUV_dynamic_bip_fc_cpp with R=1 works", {
  set.seed(42)
  nA <- 10; nB <- 8; RA <- 1; RB <- 1; TT <- 3

  U_cube <- array(rnorm(nA*RA*TT), dim = c(nA, RA, TT))
  V_cube <- array(rnorm(nB*RB*TT), dim = c(nB, RB, TT))
  E <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  G <- matrix(1, RA, RB)

  result <- rUV_dynamic_bip_fc_cpp(U_cube, V_cube, E, G,
                                    rho_uv = 0.95, sigma_uv = 0.05, s2 = 1.0)

  expect_equal(dim(result$U), c(nA, RA, TT))
  expect_equal(dim(result$V), c(nB, RB, TT))
  expect_false(any(is.nan(result$U)))
  expect_false(any(is.nan(result$V)))
})

test_that("rUV_dynamic_bip_fc_cpp with NAs in E works", {
  set.seed(42)
  nA <- 10; nB <- 8; RA <- 2; RB <- 2; TT <- 3

  U_cube <- array(rnorm(nA*RA*TT, sd = 0.5), dim = c(nA, RA, TT))
  V_cube <- array(rnorm(nB*RB*TT, sd = 0.5), dim = c(nB, RB, TT))
  E <- array(rnorm(nA*nB*TT), dim = c(nA, nB, TT))
  E[sample(length(E), 20)] <- NA
  G <- matrix(c(1, 0, 0, 1), RA, RB)

  result <- rUV_dynamic_bip_fc_cpp(U_cube, V_cube, E, G,
                                    rho_uv = 0.9, sigma_uv = 0.1, s2 = 1.0)

  expect_false(any(is.nan(result$U)))
  expect_false(any(is.nan(result$V)))
})

test_that("get_EZ_cpp optimized version matches expected structure", {
  set.seed(42)
  n <- 15; TT <- 3; p <- 2

  Xlist <- lapply(1:TT, function(t) array(rnorm(n*n*p), dim = c(n, n, p)))
  beta <- c(0.5, -0.3)
  ab <- outer(rnorm(n, sd = 0.3), rnorm(n, sd = 0.3), "+")
  U <- matrix(rnorm(n*2), n, 2)
  V <- matrix(rnorm(n*2), n, 2)

  EZ <- get_EZ_cpp(Xlist, beta, ab, U, V)

  expect_equal(dim(EZ), c(n, n, TT))
  expect_false(any(is.nan(EZ)))

  # Verify structure: EZ[,,t] should equal Xbeta + ab + U*V'
  UV <- U %*% t(V)
  for(t in 1:TT) {
    Xbeta <- beta[1] * Xlist[[t]][,,1] + beta[2] * Xlist[[t]][,,2]
    expected <- Xbeta + ab + UV
    expect_equal(EZ[,,t], expected, tolerance = 1e-10)
  }
})

test_that("latent_positions() returns correct structure", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")
  fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
             burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE)

  lp <- latent_positions(fit)

  expect_true(is.data.frame(lp))
  expect_true(all(c("actor", "dimension", "time", "value", "posterior_sd", "type") %in% names(lp)))
  expect_true(all(lp$type %in% c("U", "V")))
  n <- nrow(YX_nrm$Y)
  expect_equal(nrow(lp), n * 2 * 2)  # n actors * 2 dims * 2 types (U+V)
  # Without posterior_options, posterior_sd should be NA
  expect_true(all(is.na(lp$posterior_sd)))
})

test_that("latent_positions() includes posterior_sd when samples available", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")
  fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
             burn = 5, nscan = 20, odens = 1, verbose = FALSE, gof = FALSE,
             posterior_opts = posterior_options(save_UV = TRUE))

  lp <- latent_positions(fit)

  expect_true("posterior_sd" %in% names(lp))
  # With samples, posterior_sd should be finite (not NA)
  expect_true(all(!is.na(lp$posterior_sd)))
  expect_true(all(lp$posterior_sd >= 0))
  expect_true(all(is.finite(lp$posterior_sd)))
})

test_that("latent_positions() returns empty df for R=0", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")
  fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0,
             burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE)

  lp <- latent_positions(fit)

  expect_true(is.data.frame(lp))
  expect_equal(nrow(lp), 0)
  expect_true(all(c("actor", "dimension", "time", "value", "posterior_sd", "type") %in% names(lp)))
})

test_that("procrustes_align() works on fitted object", {
  skip_if_not_installed("lame")
  data(YX_bin_list, package = "lame")
  fit <- lame(YX_bin_list$Y, Xdyad = YX_bin_list$X, R = 2,
              family = "binary", dynamic_uv = TRUE,
              burn = 5, nscan = 10, odens = 1,
              verbose = FALSE, gof = FALSE)

  aligned <- procrustes_align(fit)

  expect_true(is.list(aligned))
  expect_equal(dim(aligned$U), dim(fit$U))
  expect_false(any(is.nan(aligned$U)))

  # return_fit should preserve class
  fit_al <- procrustes_align(fit, return_fit = TRUE)
  expect_s3_class(fit_al, "lame")
  expect_equal(dim(fit_al$U), dim(fit$U))
})

test_that("procrustes_align() on static model returns unchanged", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")
  fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
             burn = 5, nscan = 10, odens = 1, verbose = FALSE, gof = FALSE)

  result <- procrustes_align(fit)

  expect_equal(result$U, fit$U)
})

test_that("verbose deprecation warning fires for print=FALSE", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")

  expect_warning(
    fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0,
               burn = 5, nscan = 10, odens = 1, print = FALSE, gof = FALSE),
    "deprecated"
  )
})

test_that("lame() accepts 3D array Y", {
  skip_if_not_installed("lame")
  data(YX_bin_list, package = "lame")

  # Convert list Y to 3D array
  Y_list <- YX_bin_list$Y
  n <- nrow(Y_list[[1]])
  TT <- length(Y_list)
  Y_array <- array(NA, dim = c(n, n, TT),
                   dimnames = list(rownames(Y_list[[1]]),
                                   colnames(Y_list[[1]]),
                                   names(Y_list)))
  for(t in 1:TT) Y_array[,,t] <- Y_list[[t]]

  fit_array <- lame(Y_array, R = 0, family = "binary",
                     burn = 5, nscan = 10, odens = 1,
                     verbose = FALSE, gof = FALSE)

  expect_s3_class(fit_array, "lame")
  expect_true(length(fit_array$BETA) > 0)
})

test_that("full ame normal model produces valid output", {
  skip_if_not_installed("lame")
  data(YX_nrm, package = "lame")
  fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
             burn = 50, nscan = 100, odens = 1, verbose = FALSE, gof = TRUE)

  # Check all outputs are finite
  expect_false(any(is.nan(fit$BETA)))
  expect_false(any(is.nan(fit$VC)))
  expect_false(any(is.nan(fit$APM)))
  expect_false(any(is.nan(fit$BPM)))
  expect_false(any(is.nan(fit$YPM)))
  expect_false(any(is.nan(fit$U)))
  expect_false(any(is.nan(fit$V)))

  # GOF should have right structure
  expect_true(!is.null(fit$GOF))
  expect_false(any(is.nan(fit$GOF)))
})

test_that("full lame normal unipartite model produces valid output", {
  skip_if_not_installed("lame")
  data(YX_bin_list, package = "lame")
  fit <- lame(YX_bin_list$Y, Xdyad = YX_bin_list$X, R = 2,
              family = "binary",
              burn = 50, nscan = 100, odens = 1, verbose = FALSE, gof = TRUE)

  expect_false(any(is.nan(fit$BETA)))
  expect_false(any(is.nan(fit$VC)))
  expect_true(!is.null(fit$YPM))
})

test_that("full lame binary bipartite model produces valid output", {
  skip_if_not_installed("lame")
  set.seed(42)
  nA <- 12; nB <- 8; TT <- 3
  Y_list <- lapply(1:TT, function(t) {
    m <- matrix(rbinom(nA*nB, 1, 0.4), nA, nB,
                dimnames = list(paste0("r", 1:nA), paste0("c", 1:nB)))
    m
  })
  names(Y_list) <- paste0("t", 1:TT)
  X_list <- lapply(1:TT, function(t) {
    array(rnorm(nA*nB), dim = c(nA, nB, 1),
          dimnames = list(paste0("r", 1:nA), paste0("c", 1:nB), "x1"))
  })

  fit <- lame(Y_list, Xdyad = X_list, R = 2, family = "binary",
              mode = "bipartite", burn = 50, nscan = 100, odens = 1,
              verbose = FALSE, gof = TRUE)

  expect_false(any(is.nan(fit$BETA)))
  expect_false(any(is.nan(fit$VC)))
})

test_that("full lame with dynamic_uv bipartite produces valid output", {
  skip_if_not_installed("lame")
  set.seed(42)
  nA <- 10; nB <- 8; TT <- 3
  Y_list <- lapply(1:TT, function(t) {
    matrix(rnorm(nA*nB), nA, nB,
           dimnames = list(paste0("r", 1:nA), paste0("c", 1:nB)))
  })
  names(Y_list) <- paste0("t", 1:TT)
  X_list <- lapply(1:TT, function(t) {
    array(rnorm(nA*nB), dim = c(nA, nB, 1),
          dimnames = list(paste0("r", 1:nA), paste0("c", 1:nB), "x1"))
  })

  fit <- lame(Y_list, Xdyad = X_list, R = 2, family = "normal",
              mode = "bipartite", dynamic_uv = TRUE,
              burn = 20, nscan = 50, odens = 1,
              verbose = FALSE, gof = TRUE)

  expect_s3_class(fit, "lame")
  expect_false(any(is.nan(fit$BETA)))
  expect_false(any(is.nan(fit$VC)))
  # U and V should be 3D for dynamic
  expect_equal(length(dim(fit$U)), 3)
  expect_equal(length(dim(fit$V)), 3)
})
