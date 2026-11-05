# the fitted value of a Poisson MCMC fit is the mean count at the posterior means of the log mean and the residual
# variance, E_t exp(EZ + ve / 2), the form the fast estimators report; an average of replicated counts has no
# finite expectation under the inverse-gamma prior on the residual variance
skip_on_cran()

poisson_panel = function(n = 9L, Tn = 3L, bip = FALSE, nB = 6L, expo = NULL, seed = 1) {
	set.seed(seed)
	nB = if (bip) nB else n
	ra = paste0("r", seq_len(n)); ca = if (bip) paste0("c", seq_len(nB)) else ra
	a = rnorm(n, 0, .4); b = rnorm(nB, 0, .4)
	Y = X = vector("list", Tn)
	for (t in seq_len(Tn)) {
		x = matrix(rnorm(n * nB), n, nB)
		z = .2 + .3 * x + outer(a, b, "+") + matrix(rnorm(n * nB, 0, .6), n, nB)
		y = matrix(rpois(n * nB, (if (is.null(expo)) 1 else expo[t]) * exp(z)), n, nB)
		if (!bip) diag(y) = NA
		dimnames(y) = list(ra, ca); Y[[t]] = y
		X[[t]] = array(x, c(n, nB, 1), dimnames = list(ra, ca, "x"))
	}
	list(Y = Y, X = X)
}

# mean count of the fit at its posterior means, built from the pieces
expected_mean = function(fit, ez, t = 1L, expo = NULL, diag_na = TRUE) {
	m = exp(ez + mean(fit$VC[, "ve"]) / 2) * (if (is.null(expo)) 1 else expo[t])
	if (diag_na) diag(m) = NA
	m
}

test_that("ame() Poisson YPM is exp(EZ + ve / 2) at the posterior means, finite, with an empty diagonal", {
	d = poisson_panel(Tn = 1L)
	fit = suppressWarnings(ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		family = "poisson",
		R = 1,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	))
	ez = reconstruct_EZ(fit)
	expect_equal(unname(fit$YPM), unname(expected_mean(fit, ez)), tolerance = 1e-12)
	expect_identical(dimnames(fit$YPM), dimnames(d$Y[[1]]))
	expect_true(all(is.na(diag(fit$YPM))))
	expect_true(all(is.finite(fit$YPM[!is.na(fit$YPM)])) && all(fit$YPM[!is.na(fit$YPM)] > 0))
	# every accessor reads the same value, silently
	expect_silent({ f = fitted(fit) }); expect_silent({ p = predict(fit) })
	expect_identical(f, fit$YPM); expect_identical(p, fit$YPM)
	expect_equal(residuals(fit), d$Y[[1]] - fit$YPM, tolerance = 1e-12)
	# at a new design, the same form at the new linear predictor
	pn = predict(fit, newdata = d$X[[1]], type = "response")
	expect_equal(
		unname(pn),
		unname(exp(predict(fit, newdata = d$X[[1]], type = "link") + mean(fit$VC[, "ve"]) / 2)),
		tolerance = 1e-12
	)
	# independent of the replicates drawn for the goodness-of-fit statistics
	fit2 = suppressWarnings(ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		family = "poisson",
		R = 1,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = TRUE,
		seed = 2
	))
	expect_identical(fit2$YPM, fit$YPM)
})

test_that("a symmetric ame() Poisson fit has a symmetric mean count", {
	set.seed(3); n = 9L; nm = paste0("a", 1:n)
	z = outer(rnorm(n, 0, .4), rnorm(n, 0, .4), "+") + .2
	Y = matrix(rpois(n * n, exp(z + matrix(rnorm(n * n, 0, .5), n, n))), n, n)
	Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]; diag(Y) = NA; dimnames(Y) = list(nm, nm)
	fit = suppressWarnings(ame(
		Y,
		family = "poisson",
		symmetric = TRUE,
		R = 0,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	))
	expect_equal(unname(fit$YPM), unname(expected_mean(fit, reconstruct_EZ(fit))), tolerance = 1e-12)
	expect_true(isSymmetric(unname(fit$YPM)))
	expect_true(all(is.na(diag(fit$YPM))))
})

test_that("a bipartite ame() Poisson fit keeps every cell of its mean count", {
	d = poisson_panel(n = 10L, Tn = 1L, bip = TRUE, nB = 7L)
	fit = suppressWarnings(ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		mode = "bipartite",
		family = "poisson",
		R_row = 1,
		R_col = 1,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	))
	expect_equal(unname(fit$YPM), unname(expected_mean(fit, reconstruct_EZ(fit), diag_na = FALSE)), tolerance = 1e-12)
	expect_identical(dimnames(fit$YPM), dimnames(d$Y[[1]]))
	expect_false(anyNA(fit$YPM))
	# a square bipartite grid keeps its diagonal
	d2 = poisson_panel(n = 8L, Tn = 1L, bip = TRUE, nB = 8L)
	fit2 = suppressWarnings(ame(
		d2$Y[[1]],
		Xdyad = d2$X[[1]],
		mode = "bipartite",
		family = "poisson",
		R_row = 0,
		R_col = 0,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	))
	expect_false(anyNA(fit2$YPM))
})

test_that("lame() Poisson YPM carries the period exposure, static and dynamic, uni- and bipartite", {
	skip_on_cran()
	expo = c(1, 2, .5)
	d = poisson_panel(expo = expo)
	for (dyn in c(FALSE, TRUE)) {
		fit = suppressWarnings(lame(
			d$Y,
			Xdyad = d$X,
			family = "poisson",
			R = 1,
			period_exposure = expo,
			dynamic_ab = dyn,
			dynamic_uv = dyn,
			nscan = 60,
			burn = 20,
			odens = 2,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			seed = 4
		))
		expect_length(fit$YPM, 3L)
		for (t in 1:3) {
			expect_equal(fit$YPM[[t]], expected_mean(fit, fit$EZ[[t]], t, expo), tolerance = 1e-12)
			expect_true(all(is.na(diag(fit$YPM[[t]]))))
		}
		expect_silent({ f = fitted(fit) }); expect_silent({ p = predict(fit) })
		expect_identical(f, fit$YPM); expect_identical(p, fit$YPM)
		r = residuals(fit)
		for (t in 1:3) expect_equal(unname(r[[t]]), unname(d$Y[[t]] - fit$YPM[[t]]), tolerance = 1e-12)
	}
	# without an exposure the mean count is exp(EZ + ve / 2) in every period
	fit0 = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 0,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 4
	))
	for (t in 1:3) expect_equal(fit0$YPM[[t]], expected_mean(fit0, fit0$EZ[[t]], t), tolerance = 1e-12)
	db = poisson_panel(n = 8L, Tn = 3L, bip = TRUE, nB = 5L, expo = expo)
	fitb = suppressWarnings(lame(
		db$Y,
		Xdyad = db$X,
		mode = "bipartite",
		family = "poisson",
		R_row = 1,
		R_col = 1,
		period_exposure = expo,
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 4
	))
	for (t in 1:3) {
		expect_equal(fitb$YPM[[t]], expected_mean(fitb, fitb$EZ[[t]], t, expo, diag_na = FALSE), tolerance = 1e-12)
		expect_false(anyNA(fitb$YPM[[t]]))
	}
})

test_that("a zero-exposure period has a structural zero mean count", {
	expo = c(1, 0, 2)
	d = poisson_panel(expo = expo)
	d$Y[[2]][] = ifelse(is.na(d$Y[[2]]), NA, 0)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 0,
		period_exposure = expo,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 4
	))
	expect_true(all(fit$YPM[[2]][!is.na(fit$YPM[[2]])] == 0))
	expect_true(all(fit$YPM[[1]][!is.na(fit$YPM[[1]])] > 0))
})

test_that("pooled chains and per-draw response draws use the same mean count", {
	skip_on_cran()
	d = poisson_panel(Tn = 1L)
	one = function(s) suppressWarnings(ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		family = "poisson",
		R = 1,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = s
	))
	f1 = one(1); f2 = one(2)
	pooled = suppressWarnings(combine_ame_chains(list(f1, f2), diagnostics = FALSE))
	# formed from the pooled posterior means, not averaged over the chains' fitted values
	expect_equal(unname(pooled$YPM), unname(expected_mean(pooled, reconstruct_EZ(pooled))), tolerance = 1e-12)
	expect_false(isTRUE(all.equal(pooled$YPM, (f1$YPM + f2$YPM) / 2, tolerance = 1e-12)))
	# the same for panels
	dl = poisson_panel(Tn = 2L)
	lone = function(s) suppressWarnings(lame(
		dl$Y,
		Xdyad = dl$X,
		family = "poisson",
		R = 0,
		period_exposure = c(1, 3),
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = s
	))
	l1 = lone(1); l2 = lone(2)
	lp = suppressWarnings(combine_ame_chains(list(l1, l2), diagnostics = FALSE))
	for (t in 1:2) expect_equal(lp$YPM[[t]], expected_mean(lp, lp$EZ[[t]], t, c(1, 3)), tolerance = 1e-12)
	# a draw of the response is the conditional mean of that draw
	resp = prediction_draws_long(l1, type = "response", n_draws = 4, seed = 9)
	link = prediction_draws_long(l1, type = "link", n_draws = 4, seed = 9)
	key = c(".draw", "period", "i", "j")
	both = merge(resp, link, by = key)
	expect_equal(
		both$.value.x,
		c(1, 3)[both$period] * exp(both$.value.y + l1$VC[both$.iteration.x, "ve"] / 2),
		tolerance = 1e-12
	)
})

test_that("the Poisson fitted value tracks the observed level of the counts it was fitted to", {
	skip_on_cran()
	d = poisson_panel(n = 14L, Tn = 3L, seed = 11)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 1,
		nscan = 400,
		burn = 150,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5
	))
	# the fit sorts actors by name, so counts are read from the fit
	y = as.vector(fit$Y); m = unlist(fit$YPM); ok = is.finite(y) & is.finite(m)
	expect_gt(cor(y[ok], m[ok]), .4)
	expect_equal(mean(m[ok]) / mean(y[ok]), 1, tolerance = .12)
})

test_that("the Poisson fitted value is exp of the posterior mean predictor taken draw by draw", {
	skip_on_cran()
	d = poisson_panel(Tn = 1L)
	post = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 1, thin_UV = 1)
	fit = suppressWarnings(ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		family = "poisson",
		R = 2,
		nscan = 80,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 3,
		posterior_opts = post
	))
	S = nrow(fit$BETA)
	eta = Reduce(`+`, lapply(seq_len(S), function(s) {
		xb = 0 * fit$Y
		for (k in seq_len(dim(fit$X)[3])) xb = xb + fit$BETA[s, k] * fit$X[, , k]
		xb + outer(fit$a_samples[, s], fit$b_samples[, s], "+") +
			tcrossprod(fit$U_samples[, , s], fit$V_samples[, , s])
	})) / S
	want = exp(eta + mean(fit$VC[, "ve"]) / 2)
	diag(want) = NA
	expect_equal(unname(fit$YPM), unname(want), tolerance = 1e-10)
	# the same through the stored per-draw pieces of a bipartite fit
	db = poisson_panel(n = 8L, Tn = 1L, bip = TRUE, nB = 5L)
	fb = suppressWarnings(ame(
		db$Y[[1]],
		Xdyad = db$X[[1]],
		mode = "bipartite",
		family = "poisson",
		R_row = 2,
		R_col = 1,
		nscan = 80,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		gof = FALSE,
		seed = 3,
		posterior_opts = post
	))
	mat = function(a, s) matrix(a[, , s], dim(a)[1L], dim(a)[2L])
	etab = Reduce(`+`, lapply(seq_len(nrow(fb$BETA)), function(s) {
		xb = 0 * db$Y[[1]]
		for (k in seq_len(dim(fb$X)[3])) xb = xb + fb$BETA[s, k] * fb$X[, , k]
		xb + outer(fb$a_samples[, s], fb$b_samples[, s], "+") +
			mat(fb$U_samples, s) %*% mat(fb$G_samples, s) %*% t(mat(fb$V_samples, s))
	})) / nrow(fb$BETA)
	expect_equal(unname(fb$YPM), unname(exp(etab + mean(fb$VC[, "ve"]) / 2)), tolerance = 1e-10)
})

test_that("a Poisson panel with moving coefficients reports the mean count of each period's own predictor", {
	expo = c(1, 2, .5)
	d = poisson_panel(Tn = 3L, expo = expo, seed = 8)
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 0,
		period_exposure = expo,
		dynamic_beta = "dyad",
		nscan = 60,
		burn = 20,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 7
	))
	expect_identical(length(dim(fit$BETA)), 3L)
	for (t in 1:3) expect_equal(fit$YPM[[t]], expected_mean(fit, fit$EZ[[t]], t, expo), tolerance = 1e-12)
	# each period's own coefficient enters, so the predictor differs across periods
	expect_false(isTRUE(all.equal(fit$EZ[[1]], fit$EZ[[3]])))
})

test_that("a checkpoint written during a Poisson run carries mean counts, not replicate sums", {
	expo = c(1, 2, .5)
	d = poisson_panel(expo = expo)
	out = tempfile(fileext = ".rda")
	fit = suppressWarnings(lame(
		d$Y,
		Xdyad = d$X,
		family = "poisson",
		R = 0,
		period_exposure = expo,
		nscan = 40,
		burn = 10,
		odens = 2,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 6,
		periodic_save = TRUE,
		out_file = out,
		save_interval = .5
	))
	expect_true(file.exists(out))
	e = new.env()
	load(out, envir = e)
	saved = e$fit
	expect_length(saved$YPM, 3L)
	expect_identical(lapply(saved$YPM, dim), lapply(fit$YPM, dim))
	expect_true(all(vapply(saved$YPM, function(m) all(is.na(diag(m))), logical(1))))
})
