# dynamic beta als regression tests

skip_on_cran()

test_that("penalised ALS with lambda=0 equals per-period least squares", {
	set.seed(1)
	n = 12L; T_per = 4L
	X = replicate(T_per, array(rnorm(n*n*2), c(n, n, 2)), simplify = FALSE)
	beta_true = rbind(seq(-1, 1, length.out = T_per),
	                   seq(0.5, -0.5, length.out = T_per))
	Y = vector("list", T_per)
	for (t in seq_len(T_per)) {
		Yt = X[[t]][, , 1] * beta_true[1, t] +
		      X[[t]][, , 2] * beta_true[2, t] +
		      matrix(rnorm(n*n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	fit = als_dynamic_beta(Y, X, lambda = 0, intercept = FALSE)
	expect_s3_class(fit, "als_dynamic_beta")
	expect_equal(dim(fit$beta), c(2L, T_per))
	# per-period least squares by hand
	beta_ols = matrix(NA_real_, 2L, T_per)
	for (t in seq_len(T_per)) {
		ok = !is.na(Y[[t]])
		yobs = Y[[t]][ok]
		Xobs = cbind(X[[t]][, , 1][ok], X[[t]][, , 2][ok])
		beta_ols[, t] = as.numeric(solve(crossprod(Xobs), crossprod(Xobs, yobs)))
	}
	expect_lt(max(abs(fit$beta - beta_ols)), 1e-6)
})

test_that("penalised ALS with large lambda gives near-constant beta path", {
	set.seed(2)
	n = 12L; T_per = 5L
	X = replicate(T_per, array(rnorm(n*n*1), c(n, n, 1)), simplify = FALSE)
	beta_true = seq(-1, 1, length.out = T_per)
	Y = vector("list", T_per)
	for (t in seq_len(T_per)) {
		Yt = X[[t]][, , 1] * beta_true[t] +
		      matrix(rnorm(n*n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	fit_smooth = als_dynamic_beta(Y, X, lambda = 1e6, intercept = FALSE)
	expect_lt(diff(range(fit_smooth$beta)), 1e-3)
})

test_that("penalised ALS coef() returns the beta matrix", {
	set.seed(3)
	n = 8L; T_per = 3L
	X = replicate(T_per, array(rnorm(n*n*1), c(n, n, 1)), simplify = FALSE)
	Y = lapply(X, function(xt) {
		Yt = xt[, , 1] + matrix(rnorm(n*n, 0, 0.5), n, n)
		diag(Yt) = NA; Yt
	})
	fit = als_dynamic_beta(Y, X, lambda = 1, intercept = TRUE)
	cf = coef(fit)
	expect_identical(cf, fit$beta)
	expect_equal(nrow(cf), 2L)
})

test_that("penalised ALS errors on negative lambda", {
	set.seed(4)
	n = 5L
	X = list(array(rnorm(n*n), c(n, n, 1)), array(rnorm(n*n), c(n, n, 1)))
	Y = lapply(X, function(xt) xt[, , 1] + rnorm(n*n, 0, 0.1))
	expect_error(als_dynamic_beta(Y, X, lambda = -1), "non-negative")
})
