test_that("all observed Poisson scores use period exposure and saved parameters", {
	skip_on_cran()
	set.seed(29011)
	Y = lapply(1:3, function(t) {y = matrix(rpois(36, 2), 6, 6); diag(y) = NA; y})
	Y[[1]][is.finite(Y[[1]])] = 0
	Y[[1]][1, 2] = NA
	for (method in c("observed_exact", "observed_ghk")) {
		fit = lame(
			Y,
			family = "poisson",
			R = 0,
			rvar = FALSE,
			cvar = FALSE,
			dcor = FALSE,
			nscan = 5,
			burn = 3,
			odens = 1,
			seed = 473290,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE,
			period_exposure = c(0, .5, 2),
			save_log_lik = TRUE,
			log_lik_method = method
		)
		index = fit$log_lik_index
		expect_true(all(is.finite(fit$log_lik)))
		expect_identical(unname(fit$log_lik[, index[, 3] == 1, drop = FALSE]), matrix(0, 5, sum(index[, 3] == 1)))
		expect_true(all(fit$YPM[[1]][is.finite(fit$YPM[[1]])] == 0))
		for (s in 1:5) for (t in 2:3) {
			cols = which(index[, 3] == t)
			ij = index[cols, 1:2, drop = FALSE]
			eta = rep(fit$BETA[s, 1] + log(c(0, .5, 2)[t]), length(cols))
			ref = lame:::.lame_pois_marginal_loglik(Y[[t]][ij], eta, fit$VC[s, "ve"])
			expect_equal(fit$log_lik[s, cols], ref, tolerance = 1e-12)
		}
	}
	Y[[1]][2, 1] = 1
	expect_error(
		lame(
			Y,
			family = "poisson",
			R = 0,
			nscan = 5,
			burn = 3,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			period_exposure = c(0, .5, 2)
		),
		"strictly positive"
	)
})

test_that("zero exposure leaves the native directed latent Gaussian target unchanged", {
	Z = matrix(seq(-1.7, 1.8, length.out = 16), 4, 4)
	eta = matrix(seq(.1, .4, length.out = 16), 4, 4)
	Y = matrix(0, 4, 4); diag(Y) = NA; Y[1, 3] = NA
	for (rho in c(-.6, 0, .6)) {
		v = 2; e = Z - eta
		ref = -(e * e + t(e) * t(e) - 2 * rho * e * t(e)) / (2 * v * (1 - rho * rho))
		diag(ref) = -diag(e)^2 / (2 * v * (1 + rho))
		actual = lame:::ldZgbme_pois_exposure_cpp(Z, Y, eta, rho, v, -Inf)
		expect_equal(actual, ref, tolerance = 1e-12)
	}
})

test_that("zero exposure is consistent for symmetric and rectangular fitted periods", {
	for (bip in c(FALSE, TRUE)) {
		n = 6; m = if (bip) 4 else 6
		zero = matrix(0, n, m)
		counts = matrix(rep(c(0, 1, 2), length.out = n * m), n, m)
		if (!bip) {
			counts[lower.tri(counts)] = t(counts)[lower.tri(counts)]
			diag(counts) = NA
			# a structural diagonal is not an observed count
			diag(zero) = 123
		}
		fit = suppressMessages(lame(
			list(zero, counts),
			family = "poisson",
			mode = if (bip) "bipartite" else "unipartite",
			symmetric = !bip,
			R = 0,
			nvar = FALSE,
			rvar = FALSE,
			cvar = FALSE,
			dcor = FALSE,
			nscan = 4,
			burn = 2,
			odens = 1,
			seed = 473291,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			period_exposure = c(0, 2),
			save_log_lik = TRUE
		))
		idx = fit$log_lik_index
		expect_equal(sum(idx[, 3] == 1), if (bip) n * m else n * (n - 1) / 2)
		expect_true(all(fit$log_lik[, idx[, 3] == 1] == 0))
		expect_true(all(fit$YPM[[1]][is.finite(fit$YPM[[1]])] == 0))
		cols = which(idx[, 3] == 2); pairs = idx[cols, 1:2, drop = FALSE]
		for (s in 1:4) {
			ref = lame:::.lame_pois_marginal_loglik(
				counts[pairs],
				rep(fit$BETA[s, 1] + log(2), length(cols)),
				fit$VC[s, "ve"]
			)
			expect_equal(fit$log_lik[s, cols], ref, tolerance = 1e-12)
		}
	}
})
