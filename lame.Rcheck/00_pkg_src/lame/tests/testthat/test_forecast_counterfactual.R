# forecast and counterfactual newdata tests

build_fc_fit = function(seed = 2026, nscan = 300, burn = 80) {
	set.seed(seed)
	n = 25; Tn = 5
	rho_true = 0.7; beta_bar = 0.5
	beta_t = numeric(Tn); beta_t[1] = 0.9
	for (t in 2:Tn) beta_t[t] = beta_bar + rho_true * (beta_t[t-1] - beta_bar) +
		rnorm(1, 0, 0.15)
	Xdyad = lapply(seq_len(Tn), function(t) {
		x = matrix(rnorm(n*n), n, n)
		array(x, c(n, n, 1), dimnames = list(NULL, NULL, "trade"))
	})
	a = rnorm(n, 0, 0.4); b = rnorm(n, 0, 0.4)
	Y = lapply(seq_len(Tn), function(t) {
		eta = -0.5 + beta_t[t] * Xdyad[[t]][, , 1] + outer(a, b, "+")
		Yt = matrix(rbinom(n*n, 1, pnorm(eta)), n, n); diag(Yt) = NA
		rownames(Yt) = colnames(Yt) = sprintf("a%02d", seq_len(n)); Yt
	})
	names(Y) = paste0("t", seq_len(Tn))
	list(fit = lame(Y, Xdyad = Xdyad, family = "binary", R = 0,
	                dynamic_beta = "dyad", nscan = nscan, burn = burn,
	                odens = 5, verbose = FALSE),
	     Xdyad = Xdyad, Tn = Tn)
}

test_that("counterfactual raising a positive-coef covariate raises the link predictor", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit; Xdyad = obj$Xdyad; Tn = obj$Tn
		# coefficient is positive
	expect_gt(mean(coef(fit)["trade_dyad", ]), 0.2)

	last_X = Xdyad[[Tn]]
	X_future = list(last_X, last_X, last_X)
	X_up = lapply(X_future, function(x) { x[, , 1] = x[, , 1] + 1; x })

	lk_cf = predict(fit, h = 3, type = "link", newdata = X_future)
	lk_up = predict(fit, h = 3, type = "link", newdata = X_up)

		# the linear predictor should increase
	delta_link = mean(lk_up[[3]] - lk_cf[[3]], na.rm = TRUE)
	expect_gt(delta_link, 0)

	# response-scale mean shift is also positive
	rc_cf = predict(fit, h = 3, type = "response", newdata = X_future)
	rc_up = predict(fit, h = 3, type = "response", newdata = X_up)
	expect_gt(mean(rc_up[[3]] - rc_cf[[3]], na.rm = TRUE), 0)
})

test_that("forecast newdata accepts covariates without an intercept slice", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit; Xdyad = obj$Xdyad
		# user passes the substantive covariate without an intercept slice
	expect_no_error(
		predict(fit, h = 2, type = "response", newdata = list(Xdyad[[1]], Xdyad[[2]])))
})

test_that("forecast newdata with the wrong slice count errors clearly", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit
	bad = array(0, c(25, 25, 4))   # 4 slices, model has 1 dyadic coefficient
	expect_error(
		predict(fit, h = 1, type = "response", newdata = list(bad)),
		"slice")
})

test_that("newdata=NULL forecast is unaffected by the alignment fix", {
	skip_on_cran()
	obj = build_fc_fit()
	fit = obj$fit
		# default extrapolation uses the model's own design
	fc = predict(fit, h = 3, type = "response")
	expect_length(fc, 3L)
	expect_true(all(vapply(fc, function(m) all(m >= 0 & m <= 1, na.rm = TRUE),
	                       logical(1))))
})
