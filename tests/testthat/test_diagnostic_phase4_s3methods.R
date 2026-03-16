# s3 methods for ame and lame objects

####
# setup: fit representative models
####

# ame() unipartite binary
fit_ame_uni <- local({
  n <- 15
  dat <- simulate_test_network(
    n = n, family = "binary", R = 2,
    beta_intercept = 0.5, beta_dyad = 0.3,
    mode = "unipartite", seed = 2001
  )
  ame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "binary",
    burn = 50, nscan = 100, odens = 1,
    verbose = FALSE, gof = TRUE, seed = 42
  )
})

# ame() bipartite normal
fit_ame_bip <- local({
  nA <- 12; nB <- 8
  dat <- simulate_test_network(
    n = nA, nA = nA, nB = nB, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.5,
    mode = "bipartite", seed = 2002
  )
  ame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "normal",
    mode = "bipartite", burn = 50, nscan = 100, odens = 1,
    verbose = FALSE, gof = TRUE, seed = 42
  )
})

# lame() unipartite binary (static)
fit_lame_uni <- local({
  n <- 15
  dat <- simulate_test_network(
    n = n, n_time = 3, family = "binary", R = 2,
    beta_intercept = 0.5, beta_dyad = 0.3,
    mode = "unipartite", seed = 2003
  )
  lame(
    dat$Y, Xdyad = dat$Xdyad, R = 2, family = "binary",
    burn = 50, nscan = 100, odens = 1,
    verbose = FALSE, gof = TRUE, seed = 42
  )
})

# lame() bipartite normal (dynamic_ab + dynamic_uv)
fit_lame_bip_dyn <- local({
  nA <- 12; nB <- 8
  dat <- simulate_test_network(
    n = nA, nA = nA, nB = nB, n_time = 3, family = "normal", R = 2,
    beta_intercept = 1, beta_dyad = 0.3,
    rho_ab = 0.7, rho_uv = 0.5,
    mode = "bipartite", seed = 2004
  )
  Xdyad_list <- lapply(dat$Xdyad, function(x) x[, , 1])
  lame(
    dat$Y, Xdyad = Xdyad_list, R = 2, family = "normal",
    mode = "bipartite", dynamic_ab = TRUE, dynamic_uv = TRUE,
    burn = 50, nscan = 100, odens = 1,
    verbose = FALSE, gof = TRUE, seed = 42
  )
})

####
# print()
####

test_that("print.ame works for unipartite binary", {
  # cli output may not be captured by expect_output; just verify no error
  expect_no_error(print(fit_ame_uni))
})

test_that("print.ame works for bipartite normal", {
  expect_no_error(print(fit_ame_bip))
})

test_that("print.lame works for unipartite binary", {
  expect_no_error(print(fit_lame_uni))
})

test_that("print.lame works for bipartite dynamic", {
  expect_no_error(print(fit_lame_bip_dyn))
})

####
# summary()
####

test_that("summary.ame works for unipartite binary", {
  s <- summary(fit_ame_uni)
  expect_s3_class(s, "summary.ame")
  # summary uses $beta not $coefficients
  expect_true(!is.null(s$beta))
  expect_true(is.data.frame(s$beta) || is.matrix(s$beta))
})

test_that("summary.ame works for bipartite normal", {
  s <- summary(fit_ame_bip)
  expect_s3_class(s, "summary.ame")
  expect_true(!is.null(s$beta))
})

test_that("summary.lame works for unipartite binary", {
  s <- summary(fit_lame_uni)
  expect_s3_class(s, "summary.lame")
  expect_true(!is.null(s$beta))
})

test_that("summary.lame works for bipartite dynamic", {
  s <- summary(fit_lame_bip_dyn)
  expect_s3_class(s, "summary.lame")
  expect_true(!is.null(s$beta))
})

####
# coef()
####

test_that("coef.ame returns named vector for unipartite binary", {
  b <- coef(fit_ame_uni)
  expect_true(is.numeric(b))
  expect_true(length(b) > 0)
  expect_true(!is.null(names(b)))
  expect_equal(length(b), ncol(fit_ame_uni$BETA))
})

test_that("coef.ame returns named vector for bipartite normal", {
  b <- coef(fit_ame_bip)
  expect_true(is.numeric(b))
  expect_true(length(b) > 0)
  expect_true(!is.null(names(b)))
})

test_that("coef.lame works for unipartite binary", {
  b <- coef(fit_lame_uni)
  expect_true(is.numeric(b))
  expect_equal(length(b), ncol(fit_lame_uni$BETA))
})

test_that("coef.lame works for bipartite dynamic", {
  b <- coef(fit_lame_bip_dyn)
  expect_true(is.numeric(b))
  expect_equal(length(b), ncol(fit_lame_bip_dyn$BETA))
})

####
# vcov()
####

test_that("vcov.ame returns matrix for unipartite binary", {
  v <- vcov(fit_ame_uni)
  p <- ncol(fit_ame_uni$BETA)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), p)
  expect_equal(ncol(v), p)
  # Covariance matrix should be symmetric
  expect_equal(v, t(v))
  # Diagonal should be positive
  expect_true(all(diag(v) > 0))
})

test_that("vcov.ame returns matrix for bipartite normal", {
  v <- vcov(fit_ame_bip)
  expect_true(is.matrix(v))
  expect_true(all(diag(v) > 0))
})

test_that("vcov.lame works for unipartite binary", {
  v <- vcov(fit_lame_uni)
  expect_true(is.matrix(v))
  expect_true(all(diag(v) > 0))
})

test_that("vcov.lame works for bipartite dynamic", {
  v <- vcov(fit_lame_bip_dyn)
  expect_true(is.matrix(v))
  expect_true(all(diag(v) > 0))
})

####
# confint()
####

test_that("confint.ame returns intervals for unipartite binary", {
  ci <- confint(fit_ame_uni)
  p <- ncol(fit_ame_uni$BETA)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), p)
  expect_equal(ncol(ci), 2)
  # Upper > lower
  expect_true(all(ci[, 2] > ci[, 1]))
})

test_that("confint.ame with custom level", {
  ci90 <- confint(fit_ame_uni, level = 0.9)
  ci95 <- confint(fit_ame_uni, level = 0.95)
  # 95% intervals should be wider than 90%
  width90 <- ci90[, 2] - ci90[, 1]
  width95 <- ci95[, 2] - ci95[, 1]
  expect_true(all(width95 >= width90))
})

test_that("confint.ame with parm selection", {
  ci <- confint(fit_ame_uni, parm = "intercept")
  expect_equal(nrow(ci), 1)
})

test_that("confint.lame works for unipartite binary", {
  ci <- confint(fit_lame_uni)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), ncol(fit_lame_uni$BETA))
  expect_true(all(ci[, 2] > ci[, 1]))
})

test_that("confint.lame works for bipartite dynamic", {
  ci <- confint(fit_lame_bip_dyn)
  expect_true(is.matrix(ci))
  expect_true(all(ci[, 2] > ci[, 1]))
})

####
# predict()
####

test_that("predict.ame type=response for unipartite binary", {
  pred <- predict(fit_ame_uni, type = "response")
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 15)
  expect_equal(ncol(pred), 15)
  # For binary, probabilities should be in [0,1]
  non_na <- pred[!is.na(pred)]
  expect_true(all(non_na >= 0 & non_na <= 1))
})

test_that("predict.ame type=link for unipartite binary", {
  pred <- predict(fit_ame_uni, type = "link")
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 15)
  expect_equal(ncol(pred), 15)
  # Link scale can be any real number
})

test_that("predict.ame type=response for bipartite normal", {
  pred <- predict(fit_ame_bip, type = "response")
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 12)
  expect_equal(ncol(pred), 8)
})

test_that("predict.lame type=response for unipartite binary", {
  pred <- predict(fit_lame_uni, type = "response")
  expect_true(is.list(pred))
  expect_equal(length(pred), 3)  # 3 time points
  # Each element should be a matrix
  for (t in seq_along(pred)) {
    expect_true(is.matrix(pred[[t]]))
  }
})

test_that("predict.lame type=link for bipartite dynamic", {
  pred <- predict(fit_lame_bip_dyn, type = "link")
  expect_true(is.list(pred))
  expect_equal(length(pred), 3)
})

####
# fitted()
####

test_that("fitted.ame returns matrix for unipartite binary", {
  f <- fitted(fit_ame_uni)
  expect_true(is.matrix(f))
  expect_equal(nrow(f), 15)
  expect_equal(ncol(f), 15)
})

test_that("fitted.ame returns matrix for bipartite normal", {
  f <- fitted(fit_ame_bip)
  expect_true(is.matrix(f))
  expect_equal(nrow(f), 12)
  expect_equal(ncol(f), 8)
})

test_that("fitted.lame returns list for unipartite binary", {
  f <- fitted(fit_lame_uni)
  expect_true(is.list(f))
  expect_equal(length(f), 3)
  for (t in seq_along(f)) {
    expect_true(is.matrix(f[[t]]))
  }
})

test_that("fitted.lame returns list for bipartite dynamic", {
  f <- fitted(fit_lame_bip_dyn)
  expect_true(is.list(f))
  expect_equal(length(f), 3)
  for (t in seq_along(f)) {
    expect_true(is.matrix(f[[t]]))
    expect_equal(nrow(f[[t]]), 12)
    expect_equal(ncol(f[[t]]), 8)
  }
})

####
# residuals()
####

test_that("residuals.ame type=response for unipartite binary", {
  r <- residuals(fit_ame_uni, type = "response")
  expect_true(is.matrix(r))
  expect_equal(nrow(r), 15)
  expect_equal(ncol(r), 15)
})

test_that("residuals.ame type=pearson for bipartite normal", {
  r <- residuals(fit_ame_bip, type = "pearson")
  expect_true(is.matrix(r))
  expect_equal(nrow(r), 12)
  expect_equal(ncol(r), 8)
})

test_that("residuals.lame type=response for unipartite binary", {
  r <- suppressWarnings(residuals(fit_lame_uni, type = "response"))
  expect_true(is.list(r))
  expect_equal(length(r), 3)
})

test_that("residuals.lame type=response for bipartite dynamic", {
  r <- suppressWarnings(residuals(fit_lame_bip_dyn, type = "response"))
  expect_true(is.list(r))
  expect_equal(length(r), 3)
  for (t in seq_along(r)) {
    expect_true(is.matrix(r[[t]]))
  }
})

####
# simulate()
####

test_that("simulate.ame works for unipartite binary", {
  sim_obj <- simulate(fit_ame_uni, nsim = 3)
  expect_s3_class(sim_obj, "ame.sim")
  sims <- sim_obj$Y
  expect_true(is.list(sims))
  expect_equal(length(sims), 3)
  for (s in sims) {
    expect_true(is.matrix(s))
    expect_equal(nrow(s), 15)
    expect_equal(ncol(s), 15)
    # Binary: should be 0/1
    non_na <- s[!is.na(s)]
    expect_true(all(non_na %in% c(0, 1)))
  }
})

test_that("simulate.ame works for bipartite normal", {
  sim_obj <- simulate(fit_ame_bip, nsim = 3)
  expect_s3_class(sim_obj, "ame.sim")
  sims <- sim_obj$Y
  expect_true(is.list(sims))
  expect_equal(length(sims), 3)
  for (s in sims) {
    expect_true(is.matrix(s))
    expect_equal(nrow(s), 12)
    expect_equal(ncol(s), 8)
  }
})

test_that("simulate.lame works for unipartite binary", {
  sim_obj <- simulate(fit_lame_uni, nsim = 3)
  expect_s3_class(sim_obj, "lame.sim")
  sims <- sim_obj$Y
  expect_true(is.list(sims))
  expect_equal(length(sims), 3)
  # Each sim should be a list of time-point matrices
  for (s in sims) {
    expect_true(is.list(s))
    expect_equal(length(s), 3)  # 3 time points
    for (mat in s) {
      expect_true(is.matrix(mat))
    }
  }
})

test_that("simulate.lame works for bipartite dynamic", {
  sim_obj <- simulate(fit_lame_bip_dyn, nsim = 3)
  expect_s3_class(sim_obj, "lame.sim")
  sims <- sim_obj$Y
  expect_true(is.list(sims))
  expect_equal(length(sims), 3)
  for (s in sims) {
    expect_true(is.list(s))
    expect_equal(length(s), 3)
    for (mat in s) {
      expect_true(is.matrix(mat))
      expect_equal(nrow(mat), 12)
      expect_equal(ncol(mat), 8)
    }
  }
})

####
# plot()
####

test_that("plot.ame produces valid output for unipartite binary", {
  expect_no_error(suppressWarnings(plot(fit_ame_uni)))
})

test_that("plot.ame produces valid output for bipartite normal", {
  expect_no_error(suppressWarnings(plot(fit_ame_bip)))
})

test_that("plot.lame produces valid output for unipartite binary", {
  expect_no_error(suppressWarnings(plot(fit_lame_uni)))
})

test_that("plot.lame produces valid output for bipartite dynamic", {
  expect_no_error(suppressWarnings(plot(fit_lame_bip_dyn)))
})
