skip_on_cran()

.per_actor_target_fixture = function(coupled = FALSE) {
	n = 3L; nt = 3L; rho = .8; variance = .4
	K = variance / (1 - rho * rho) * rho^abs(outer(1:nt, 1:nt, '-'))
	H = rbind(c(.1, 5, .2), c(4, .2, 3), c(.5, 2, .1))
	h = 10 * rbind(c(.1, -.2, .5), c(.3, .1, -.6), c(-.1, .4, .2))
	P = lapply(1:nt, function(t) {
		if (!coupled) return(diag(H[, t]))
		x = matrix(c(1, .2, -.3, .5, 1, .4, -.1, .6, 1), n)
		crossprod(x) + diag(H[, t])
	})
	Q = kronecker(diag(n), solve(K))
	for (t in 1:nt) {i = t + (0:(n - 1)) * nt;Q[i, i] = Q[i, i] + P[[t]]}
	B = kronecker(qr.Q(qr(contr.helmert(n))), diag(nt))
	V = solve(crossprod(B, Q %*% B)); m = drop(V %*% crossprod(B, c(t(h))))
	list(
		n = n,
		nt = nt,
		H = H,
		h = h,
		P = if (coupled) P else NULL,
		rho = rho,
		variance = variance,
		prior_precision = solve(K),
		basis = B,
		mean = m,
		covariance = V,
		target_mean = drop(B %*% m),
		target_covariance = B %*% V %*% t(B)
	)
}

test_that("actor FFBS returns the full smoothed covariance including small scales and one period", {
	for (nt in c(1L, 3L)) for (v in c(.2, 1e-16)) {
		K = v / (1 - .7^2) * .7^abs(outer(seq_len(nt), seq_len(nt), '-'))
		z = lame:::.actor_ffbs_path(rep(.5, nt), rep(.2, nt), .7, v)
		target = solve(solve(K) + diag(.5, nt))
		expect_lt(max(abs(z$covariance / target - 1)), 1e-10)
		expect_equal(z$V, diag(z$covariance))
	}
})

test_that("joint actor centering matches independent constrained Gaussian posteriors", {
	skip_on_cran()
	set.seed(816572);N = 6000L
	for (coupled in c(FALSE, TRUE)) {
		f = .per_actor_target_fixture(coupled)
		z = replicate(N, c(t(lame:::.sweep_per_actor_exact(f$H, f$h, f$rho, f$variance, f$P))))
		mean_error = (rowMeans(z) - f$target_mean) / sqrt(diag(f$target_covariance) / N)
		covariance_se = sqrt((outer(diag(f$target_covariance), diag(f$target_covariance)) + f$target_covariance^2) / N)
		expect_lt(max(abs(mean_error)), 6)
		expect_lt(max(abs((cov(t(z)) - f$target_covariance) / covariance_se)), 6)
		expect_lt(max(abs(apply(array(z, c(f$nt, f$n, N)), c(1, 3), sum))), 1e-10)
	}
})

test_that("pairwise actor updates preserve the same Gaussian posterior with reciprocity", {
	skip_on_cran()
	set.seed(173586)
	for (coupled in c(FALSE, TRUE)) {
		f = .per_actor_target_fixture(coupled)
		d = replicate(4000L, {
			q = f$mean + drop(t(chol(f$covariance)) %*% rnorm(length(f$mean)))
			x = t(matrix(drop(f$basis %*% q), f$nt, f$n))
			y = lame:::.sweep_per_actor(x, f$H, f$h, f$prior_precision, f$P)
			c(t(y - x), t(y * y - x * x))
		})
		mcse = rowMeans(d) / (apply(d, 1, sd) / sqrt(ncol(d)))
		expect_true(all(is.finite(mcse)))
		expect_lt(max(abs(mcse)), 6)
	}
})

test_that("centered actor-slope hyperparameters preserve their stationary joint prior", {
	skip_on_cran()
	set.seed(952631);B = qr.Q(qr(contr.helmert(4)))
	for (nt in c(1L, 3L)) {
		d = replicate(5000L, {
			repeat{rho = 2 * rbeta(1, 2, 2) - 1;if (abs(rho) < .99) break}
			variance = 1 / rgamma(1, 3, rate = .4)
			K = variance / (1 - rho * rho) * rho^abs(outer(seq_len(nt), seq_len(nt), '-'))
			theta = B %*% (matrix(rnorm(3 * nt), 3, nt) %*% chol(K))
			v = lame:::.sample_sigma_actor(theta, rho, 3, .4)
			r = lame:::.sample_rho_actor(theta, v, rho)
			c(
				log(v) - log(variance),
				atan(v) - atan(variance),
				r - rho,
				r * r - rho * rho,
				r * log(v) - rho * log(variance)
			)
		})
		mcse = rowMeans(d) / (apply(d, 1, sd) / sqrt(ncol(d)))
		expect_lt(max(abs(mcse)), 6)
	}
})

test_that("reciprocal actor likelihood statistics match an independent Gaussian design", {
	set.seed(187356)
	for (kind in c('row', 'col')) for (rho in c(-.7, 0, .6)) {
		n = 4L;R = matrix(rnorm(n * n), n);X = matrix(rnorm(n * n), n);s2 = .4
		z = lame:::.lame_per_actor_suffstats(R, X, kind, s2, rho, FALSE)
		P = matrix(0, n, n);h = numeric(n)
		for (i in 1:n) {
			D = matrix(0, 1, n);D[1, i] = X[i, i]
			P = P + crossprod(D) / (s2 * (1 + rho));h = h + drop(t(D) * R[i, i] / (s2 * (1 + rho)))
			if (i < n) for (j in (i + 1):n) {
				D = matrix(0, 2, n);D[1, if (kind == 'row') i else j] = X[i, j];D[2, if (kind == 'row') j else i] = X[j, i]
				W = solve(s2 * matrix(c(1, rho, rho, 1), 2))
				P = P + crossprod(D, W %*% D);h = h + drop(crossprod(D, W %*% c(R[i, j], R[j, i])))
			}
		}
		expect_equal(if (is.null(z$precision)) diag(z$H, n) else z$precision, P, tolerance = 1e-10)
		expect_equal(z$h, h, tolerance = 1e-10)
	}
})

test_that("stored predictors include actor slopes on the selected dyadic covariate", {
	set.seed(725613)
	for (mode in c('unipartite', 'bipartite')) {
		na = 4L;nb = if (mode == 'bipartite') 5L else na;nt = 3L
		X = lapply(1:nt, function(t) array(rnorm(na * nb), c(na, nb, 1), dimnames = list(NULL, NULL, 'dx')))
		xr = lapply(1:nt, function(t) matrix(rnorm(na), na, 1, dimnames = list(NULL, 'rx')))
		xc = lapply(1:nt, function(t) matrix(rnorm(nb), nb, 1, dimnames = list(NULL, 'cx')))
		Y = lapply(1:nt, function(t) {y = matrix(rnorm(na * nb), na, nb);if (mode == 'unipartite') diag(y) = NA;y})
		f = lame(
			Y,
			Xdyad = X,
			Xrow = xr,
			Xcol = xc,
			mode = mode,
			R = 0,
			family = 'normal',
			dynamic_beta_per_actor = 'row',
			per_actor_identifiability = 'exact_center',
			keep_per_actor = 'draws',
			burn = 20,
			nscan = 80,
			odens = 2,
			seed = 73516,
			verbose = FALSE,
			gof = FALSE
		)
		expect_identical(f$per_actor_model_version, 2L)
		expect_equal(f$per_actor_design_idx, if (mode == 'bipartite') 2 else 4)
		for (t in 1:nt) {
			base = matrix(matrix(f$Xlist[[t]], na * nb, 4) %*% colMeans(f$BETA), na, nb) + outer(f$APM, f$BPM, '+')
			expected = base + sweep(f$Xlist[[t]][, , f$per_actor_design_idx], 1L, f$theta_actor_mean[, t], '*')
			ok = is.finite(f$EZ[[t]])
			expect_equal(f$EZ[[t]][ok], expected[ok], tolerance = 1e-10)
		}
		expect_true(all(unlist(f$tryErrorChecks) == 0))
	}
})
