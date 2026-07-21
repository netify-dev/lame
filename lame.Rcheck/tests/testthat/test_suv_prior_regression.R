# Regression guards for the multiplicative-effects covariance prior.
# ame() once divided prior$Suv0 by kappa0 before passing it to rSuv_fc,
# which multiplies by kappa0 internally -- deflating the effective
# inverse-Wishart scale by a factor of kappa0 and collapsing U/V to zero
# on roughly half of sparse binary datasets (with beta attenuated too).
# These seeds collapse (sd(U) ~ 0.02) under that defect and recover
# (sd(U) ~ 0.75-0.92) with the correct prior, so the margin is wide.

test_that("binary R=1 latent factors do not collapse across seeds", {
	skip_on_cran()
	sim = function(seed, n = 50, mu = -2, beta = 1, gamma = 1) {
		set.seed(seed)
		xw = matrix(rnorm(n * 2), n, 2)
		X = tcrossprod(xw[, 1])
		W = tcrossprod(xw[, 2])
		Y = 1 * (mu + beta * X + gamma * W + matrix(rnorm(n * n), n, n) > 0)
		list(Y = Y, X = X, W = W)
	}
	for (s in c(3, 5, 9)) {
		d = sim(s)
		fit = suppressMessages(
			ame(d$Y, d$X, R = 1, rvar = FALSE, cvar = FALSE, dcor = FALSE,
			    family = "binary", verbose = FALSE, plot = FALSE,
			    nscan = 800, burn = 200, odens = 10))
		expect_gt(sd(as.matrix(fit$U)), 0.3)
		uv = tcrossprod(as.matrix(fit$U), as.matrix(fit$V))
		expect_gt(cor(c(d$W), c(uv)), 0.5)
	}
})

test_that("rUV_rep_fc fallback draws Suv at the same df as the C++ path", {
	skip_on_cran()
	set.seed(8)
	n = 50; R = 1
	U = matrix(rnorm(n * R), n, R)
	V = matrix(rnorm(n * R), n, R)
	UV = cbind(U, V)
	# reference: the amen / rUV_rep_fc_cpp draw, IW(I + UV'UV, n + R + 2)
	ref = replicate(2000,
		solve(lame::rwish(solve(diag(2 * R) + t(UV) %*% UV), n + R + 2))[1, 1])
	# the fallback's call: scale kappa0 * Suv0 = I, df n + kappa0 = n + R + 2
	fb = replicate(2000,
		rSuv_fc(U, V, Suv0 = diag(2 * R) / (R + 2), kappa0 = R + 2)[1, 1])
	expect_lt(abs(mean(fb) / mean(ref) - 1), 0.1)
})

test_that("rUV_rep_fc passes the intended kappa0 to rSuv_fc", {
	# rSuv_fc adds nrow(U) to kappa0 internally, so the fallback must pass
	# R + 2 (not n + R + 2) for the draw to sit at df n + R + 2
	src = gsub("\\s+", " ", paste(deparse(body(rUV_rep_fc)), collapse = " "))
	expect_true(grepl("kappa0 = R + 2", src, fixed = TRUE))
	expect_false(grepl("kappa0 = n + R + 2", src, fixed = TRUE))
})
