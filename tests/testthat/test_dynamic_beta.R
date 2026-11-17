# dynamic_beta on lame() fits: fitting, recovery, validation guards, the in-loop per-actor slope sampler, alternative
# transition kinds, lame_multi and the penalised-ALS path; with dynamic_beta = FALSE every slot equals that of a fit without the argument

skip_on_cran()

# the dynamic_beta fit the method checks share, fitted once
.db_s3_fit = local({
	fit = NULL
	function() {
		if (is.null(fit)) {
			dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
			fit <<- lame(
				dat$Y,
				Xdyad = dat$Xdyad,
				family = "normal",
				R = 0,
				nscan = 30,
				burn = 5,
				odens = 1,
				dynamic_beta = TRUE,
				verbose = FALSE
			)
		}
		fit
	}
})

# byte-identical default

test_that("dynamic_beta = FALSE is byte-identical to no-argument default (unipartite)", {
	set.seed(2026)
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	set.seed(7)
	fit1 = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = FALSE,
		verbose = FALSE
	)
	set.seed(7)
	fit2 = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE
	)
	expect_identical(fit1$BETA, fit2$BETA)
	expect_identical(fit1$VC, fit2$VC)
	expect_identical(fit1$APM, fit2$APM)
	expect_identical(fit1$BPM, fit2$BPM)
	expect_identical(fit1$U, fit2$U)
	expect_identical(fit1$V, fit2$V)
	expect_identical(fit1$EZ, fit2$EZ)
	expect_identical(fit1$YPM, fit2$YPM)
})

test_that("dynamic_beta = FALSE is byte-identical (bipartite)", {
	set.seed(2026)
	dat = simulate_test_network(
		n = 8, n_time = 3, family = "normal",
		mode = "bipartite", nA = 6, nB = 5, seed = 2026
	)
	set.seed(7)
	fit1 = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		mode = "bipartite",
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = FALSE,
		verbose = FALSE
	)
	set.seed(7)
	fit2 = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		mode = "bipartite",
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		verbose = FALSE
	)
	expect_identical(fit1$BETA, fit2$BETA)
	expect_identical(fit1$VC, fit2$VC)
	expect_identical(fit1$APM, fit2$APM)
	expect_identical(fit1$BPM, fit2$BPM)
})

# parser correctness

test_that("lame:::parse_dynamic_beta(FALSE / NULL / TRUE / character / integer / logical) works", {
	coef_names = c("intercept", "x1.dyad", "x2.dyad", "x3.row")
	coef_block = c("intercept", "dyad", "dyad", "row")
	res = lame:::parse_dynamic_beta(FALSE, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_false(res$any)
	expect_equal(res$n_groups, 0L)
	expect_equal(res$dynamic_idx, integer(0))

	res = lame:::parse_dynamic_beta(NULL, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_false(res$any)

	res = lame:::parse_dynamic_beta(TRUE, coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(sum(res$mask), 4L)
	expect_setequal(res$group_names, c("intercept", "dyad", "row"))

	res = lame:::parse_dynamic_beta("dyad", coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(res$mask, c(FALSE, TRUE, TRUE, FALSE))

	res = lame:::parse_dynamic_beta(c("intercept", "x3.row"), coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_true(res$any)
	expect_equal(res$mask, c(TRUE, FALSE, FALSE, TRUE))

	res = lame:::parse_dynamic_beta(c(2L, 3L), coef_names, coef_block, TRUE, 3, "normal", "unipartite")
	expect_equal(res$mask, c(FALSE, TRUE, TRUE, FALSE))

	res = lame:::parse_dynamic_beta(
		c(FALSE, TRUE, FALSE, FALSE),
		coef_names,
		coef_block,
		TRUE,
		3,
		"normal",
		"unipartite"
	)
	expect_equal(sum(res$mask), 1L)
})

test_that("parse_dynamic_beta family ban (ordinal) is enforced", {
	coef_names = c("intercept", "x1.dyad")
	coef_block = c("intercept", "dyad")
	expect_error(
		lame:::parse_dynamic_beta(TRUE, coef_names, coef_block, TRUE, 3, "ordinal", "unipartite"),
		"intercept"
	)
})

test_that("parse_dynamic_beta unknown name aborts", {
	expect_error(
		lame:::parse_dynamic_beta(
			"nonexistent_coef",
			c("intercept", "x1.dyad"),
			c("intercept", "dyad"),
			TRUE,
			3,
			"normal",
			"unipartite"
		),
		"not a coefficient or block shortcut"
	)
})

test_that("parse_dynamic_beta Tn < 2 aborts", {
	expect_error(
		lame:::parse_dynamic_beta(
			TRUE,
			c("intercept", "x1.dyad"),
			c("intercept", "dyad"),
			TRUE,
			1L,
			"normal",
			"unipartite"
		),
		"at least 2 time periods"
	)
})

# 3-D BETA shape

test_that("dynamic_beta = TRUE produces 3-D BETA with right dim/names", {
	fit = .db_s3_fit()
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_equal(dim(fit$BETA), c(30L, 2L, 3L))
	expect_equal(dimnames(fit$BETA)[[2]], c("intercept", "X1_dyad"))
	expect_true(isTRUE(fit$dynamic_beta))
	expect_false(is.null(fit$rho_beta))
	expect_false(is.null(fit$sigma_beta))
	expect_equal(length(fit$rho_beta), 2L) # intercept + dyad blocks
})

test_that("partial dynamic_beta (subset of coefs) gives the right mask", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = "dyad",
		verbose = FALSE
	)
	# the static intercept is the same in every period of a draw
	beta_int = fit$BETA[, "intercept", ]
	expect_true(all(beta_int[, 1] == beta_int[, 2]))
	expect_true(all(beta_int[, 2] == beta_int[, 3]))
	# the dynamic dyadic coefficient differs across periods (with overwhelming probability)
	beta_dyad = fit$BETA[, "X1_dyad", ]
	expect_false(all(beta_dyad[, 1] == beta_dyad[, 2]))
})

# S3 methods

test_that("coef.lame returns p x Tn matrix for dynamic_beta fit", {
	fit = .db_s3_fit()
	cf = coef(fit)
	expect_true(is.matrix(cf))
	expect_equal(dim(cf), c(2L, 3L))
	expect_equal(rownames(cf), c("intercept", "X1_dyad"))
})

test_that("vcov returns (p*Tn) x (p*Tn) for dynamic_beta fit", {
	fit = .db_s3_fit()
	V = vcov(fit)
	expect_equal(dim(V), c(6L, 6L)) # 2 coefs x 3 periods
	expect_true(any(grepl("\\[t1\\]", rownames(V))))
})

test_that("confint returns (p*Tn) x 2 for dynamic_beta fit", {
	fit = .db_s3_fit()
	ci = confint(fit)
	expect_equal(dim(ci), c(6L, 2L))
	expect_true(all(ci[, 1] <= ci[, 2]))
})

test_that("summary.lame prints dynamic-beta block without erroring", {
	fit = .db_s3_fit()
	out = summary(fit)
	expect_s3_class(out, "summary.lame")
	expect_true(isTRUE(out$dynamic_beta))
	expect_true(!is.null(out$beta_dynamic_per_t))
	expect_silent(invisible(capture.output(print(out))))
})

test_that("print.lame works on dynamic_beta fit", {
	fit = .db_s3_fit()
	# cli writes to the message stream; capture it
	msgs = capture.output(print(fit), type = "message")
	out = capture.output(print(fit))
	combined = paste(c(out, msgs), collapse = "\n")
	expect_match(combined, "Dynamic regression coefficients")
})

test_that("as_draws flattens 3-D BETA to draws_array", {
	skip_if_not_installed("posterior")
	fit = .db_s3_fit()
	d = as_draws(fit)
	expect_s3_class(d, "draws_array")
	# 2 coefficients x 3 periods plus the five variance components (va, cab, vb, rho, ve)
	expect_equal(dim(d)[3], 11L)
	vn = dimnames(d)[[3]]
	expect_true(any(grepl("\\[t1\\]", vn)))
})

# identifiability and the constraint basis

test_that("a/b are sum-to-zero (within tolerance) when an intercept is in the model", {
	dat = simulate_test_network(n = 10, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		nscan = 40,
		burn = 10,
		odens = 1,
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	# sum-to-zero on a (and b) keeps the dynamic intercept identified
	expect_lt(abs(sum(fit$APM)), 1e-6)
	expect_lt(abs(sum(fit$BPM)), 1e-6)
})

# bipartite with dynamic_beta

test_that("dynamic_beta works on bipartite normal", {
	dat = simulate_test_network(
		n = 8, n_time = 3, family = "normal",
		mode = "bipartite", nA = 6, nB = 5, seed = 2026
	)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		mode = "bipartite",
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 1,
		dynamic_beta = "dyad",
		verbose = FALSE
	)
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_equal(dim(fit$BETA), c(20L, 2L, 3L))
})

# edge cases

test_that("dynamic_beta + Tn = 1 aborts (validation in lame())", {
	set.seed(1)
	Y = list(matrix(rnorm(36), 6, 6))
	diag(Y[[1]]) = NA
	expect_error(
		lame(
			Y,
			family = "normal",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			dynamic_beta = TRUE,
			verbose = FALSE
		),
		"at least 2 time periods"
	)
})

test_that("dynamic_beta with method = 'als' returns dynamic ALS fit", {
	dat = simulate_test_network(n = 8, n_time = 3, family = "normal", seed = 2026)
	fit = lame(
		dat$Y,
		Xdyad = dat$Xdyad,
		family = "normal",
		R = 0,
		method = "als",
		dynamic_beta = TRUE,
		verbose = FALSE
	)
	expect_s3_class(fit, "lame_dynamic_als")
	expect_true(isTRUE(fit$dynamic_beta))
	expect_equal(length(dim(fit$BETA)), 3L)
	expect_equal(dim(fit$BETA)[3L], length(dat$Y))
})

test_that("ame() (cross-sectional) rejects dynamic_beta", {
	set.seed(1)
	Y = matrix(rnorm(64), 8, 8); diag(Y) = NA
	# ame() catches longitudinal-only arguments in ... and points the user to lame()
	expect_error(
		ame(
			Y,
			family = "normal",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			dynamic_beta = TRUE,
			verbose = FALSE
		),
		"longitudinal-only"
	)
})

# recovery check

test_that("dynamic_beta drift is detected when truth is time-varying", {
	# y_t = beta_t * x_t + noise with beta_t moving linearly
	set.seed(2026)
	n = 10; Tn = 4
	beta_t = seq(-0.5, 0.5, length.out = Tn)
	X_list = replicate(Tn, matrix(rnorm(n * n), n, n), simplify = FALSE)
	Y_list = vector("list", Tn)
	for (t in seq_len(Tn)) {
		Y_t = beta_t[t] * X_list[[t]] + matrix(rnorm(n * n, 0, 0.3), n, n)
		diag(Y_t) = NA
		Y_list[[t]] = Y_t
	}
	X_arr = lapply(X_list, function(x) array(x, c(n, n, 1)))
	# dynamic_beta on the dyadic coefficient
	fit = lame(
		Y_list,
		Xdyad = X_arr,
		family = "normal",
		R = 0,
		nscan = 100,
		burn = 30,
		odens = 1,
		dynamic_beta = "dyad",
		verbose = FALSE
	)
	# the per-period posterior-mean dyadic coefficients span a notable range
	dyad_per_t = apply(fit$BETA[, "X1_dyad", , drop = FALSE], 3, mean)
	expect_gt(diff(range(dyad_per_t)), 0.2)
})

# alternative transition kinds, log_lik_method, per_actor_slopes and multiple panels

.fit_dynamic_beta_alt = function(
	kind = "ar1",
	ll_method = "observed_exact",
	suppress_messages = TRUE,
	...
) {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	expr = quote(suppressWarnings(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		dynamic_beta = "dyad",
		dynamic_beta_kind = kind,
		log_lik_method = ll_method,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		...
	)))
	if (suppress_messages) suppressMessages(eval(expr)) else eval(expr)
}

test_that("invalid dynamic_beta_kind aborts", {
	data(YX_bin_list, envir = environment())
	expect_error(
		suppressWarnings(suppressMessages(lame(
			YX_bin_list$Y,
			YX_bin_list$X,
			family = "binary",
			R = 0,
			dynamic_beta = "dyad",
			dynamic_beta_kind = "magic",
			nscan = 30,
			burn = 10,
			odens = 5,
			verbose = FALSE,
			plot = FALSE
		))),
		"dynamic_beta_kind"
	)
})

test_that("log_lik_method = 'observed_ghk' emits the bias note", {
	expect_message(
		.fit_dynamic_beta_alt(ll_method = "observed_ghk", suppress_messages = FALSE),
		"biased downward"
	)
})

test_that("observed_ghk produces a finite log-lik matrix on binary", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		log_lik_method = "observed_ghk",
		save_log_lik = TRUE,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_true(!is.null(fit[["log_lik"]]))
	expect_true(all(is.finite(fit[["log_lik"]])))
})

test_that("log_lik_method = 'augmented' is accepted silently", {
	expect_silent(.fit_dynamic_beta_alt(ll_method = "augmented"))
})

test_that("per_actor_slopes returns the expected shape and label", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	pas = per_actor_slopes(fit, kind = "row", covariate_idx = 1L, lambda = 1)
	expect_s3_class(pas, "per_actor_slopes")
	expect_equal(dim(pas$slopes), c(50L, 4L))
	expect_true(all(is.finite(pas$slopes)))
})

test_that("per_actor_slopes errors on bad lambda and bad kind", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fit = suppressWarnings(suppressMessages(lame(
		YX_bin_list$Y,
		YX_bin_list$X,
		family = "binary",
		R = 0,
		nscan = 30,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_error(per_actor_slopes(fit, lambda = -1), "non-negative")
	expect_error(per_actor_slopes(fit, kind = "bogus"))
})

test_that("lame_multi with K=1 falls through to standard lame()", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X),
		family = "binary",
		R = 0,
		nscan = 30,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_s3_class(fm, "lame_multi")
	expect_equal(fm$K, 1L)
})

test_that("lame_multi with K=2 identical panels has zero deviations", {
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y, YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X, YX_bin_list$X),
		family = "binary",
		R = 0,
		nscan = 30,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_equal(fm$K, 2L)
	expect_length(fm$beta_deviations, 2L)
	# identical panels have near-zero deviations
	for (k in seq_len(fm$K)) {
		expect_lt(max(abs(fm$beta_deviations[[k]])), 1e-10)
	}
})

test_that("lame_multi gives proper 1/K variance shrinkage with K identical panels", {
	skip_on_cran()
	data(YX_bin_list, envir = environment())
	set.seed(7)
	fm = suppressWarnings(suppressMessages(lame_multi(
		Y_list = list(YX_bin_list$Y, YX_bin_list$Y, YX_bin_list$Y),
		Xdyad_list = list(YX_bin_list$X, YX_bin_list$X, YX_bin_list$X),
		family = "binary",
		R = 0,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE
	)))
	expect_equal(fm$K, 3L)
	expect_lt(abs(fm$var_shrinkage - 1 / 3), 0.05)
	expect_true(!is.null(fm$BETA_joint))
	expect_true(!is.null(fm$beta_shared_se))
})

# penalised-ALS dynamic_beta (als_dynamic_beta)

test_that("penalised ALS with lambda=0 equals per-period least squares", {
	set.seed(1)
	n = 12L; T_per = 4L
	X = replicate(T_per, array(rnorm(n * n * 2), c(n, n, 2)), simplify = FALSE)
	beta_true = rbind(
		seq(-1, 1, length.out = T_per),
		seq(0.5, -0.5, length.out = T_per)
	)
	Y = vector("list", T_per)
	for (t in seq_len(T_per)) {
		Yt = X[[t]][, , 1] * beta_true[1, t] +
			X[[t]][, , 2] * beta_true[2, t] +
			matrix(rnorm(n * n, 0, 0.2), n, n)
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
	X = replicate(T_per, array(rnorm(n * n * 1), c(n, n, 1)), simplify = FALSE)
	beta_true = seq(-1, 1, length.out = T_per)
	Y = vector("list", T_per)
	for (t in seq_len(T_per)) {
		Yt = X[[t]][, , 1] * beta_true[t] +
			matrix(rnorm(n * n, 0, 0.2), n, n)
		diag(Yt) = NA
		Y[[t]] = Yt
	}
	fit_smooth = als_dynamic_beta(Y, X, lambda = 1e6, intercept = FALSE)
	expect_lt(diff(range(fit_smooth$beta)), 1e-3)
})

test_that("penalised ALS coef() returns the beta matrix", {
	set.seed(3)
	n = 8L; T_per = 3L
	X = replicate(T_per, array(rnorm(n * n * 1), c(n, n, 1)), simplify = FALSE)
	Y = lapply(X, function(xt) {
		Yt = xt[, , 1] + matrix(rnorm(n * n, 0, 0.5), n, n)
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
	X = list(array(rnorm(n * n), c(n, n, 1)), array(rnorm(n * n), c(n, n, 1)))
	Y = lapply(X, function(xt) xt[, , 1] + rnorm(n * n, 0, 0.1))
	expect_error(als_dynamic_beta(Y, X, lambda = -1), "non-negative")
})

test_that("dynamic_beta + dynamic_ab samples the AR(1) hyperparameters (regression)", {
	# the combination draws per-period additive effects under the AR(1) prior and samples rho_ab / sigma_ab;
	# values stuck at their initialization (0.8 / 0.1) mean the AR(1) draws are not reached
	skip_on_cran()
	set.seed(4)
	n = 12; TT = 5
	Y = lapply(1:TT, function(t) {
		m = matrix(rnorm(n * n), n, n); diag(m) = NA
		rownames(m) = colnames(m) = sprintf("a%02d", 1:n); m
	})
	X = lapply(1:TT, function(t) {
		a = array(rnorm(n * n), c(n, n, 1)); dimnames(a)[[3]] = "x"; a
	})
	fit = suppressWarnings(lame(
		Y,
		Xdyad = X,
		family = "normal",
		R = 0,
		dynamic_ab = TRUE,
		dynamic_beta = "dyad",
		burn = 100,
		nscan = 400,
		odens = 2,
		seed = 1,
		verbose = FALSE,
		plot = FALSE
	))
	expect_gt(sd(as.numeric(fit$rho_ab)), 1e-6)
	expect_gt(sd(as.numeric(fit$sigma_ab)), 1e-6)
	expect_false(is.null(fit$a_dynamic))
	expect_true(all(is.finite(fit$a_dynamic)))
})

test_that("symmetric dynamic_beta + dynamic_ab fits and merges coefficients (regression)", {
	# the symmetric coefficient merge indexes beta with the nodal-covariate index vectors, which must survive the AR(1) additive effect draws
	skip_on_cran()
	set.seed(3)
	n = 12; TT = 4
	Y = lapply(1:TT, function(t) {
		m = matrix(rnorm(n * n), n, n); m = (m + t(m)) / 2; diag(m) = NA
		rownames(m) = colnames(m) = sprintf("a%02d", 1:n); m
	})
	X = lapply(1:TT, function(t) {
		a = matrix(rnorm(n * n), n, n); a = (a + t(a)) / 2
		a = array(a, c(n, n, 1)); dimnames(a)[[3]] = "x"; a
	})
	fit = suppressWarnings(lame(
		Y,
		Xdyad = X,
		family = "normal",
		symmetric = TRUE,
		R = 0,
		dynamic_ab = TRUE,
		dynamic_beta = "dyad",
		burn = 50,
		nscan = 150,
		odens = 5,
		seed = 1,
		verbose = FALSE,
		plot = FALSE
	))
	expect_s3_class(fit, "lame")
	expect_true(all(is.finite(fit$a_dynamic)))
	expect_gt(sd(as.numeric(fit$rho_ab)), 1e-6)
})
