# dense Gaussian and prior-predictive targets for the R-level coefficient samplers,
# including missing reciprocal cells and nonunit gaps
test_that("alternative coefficient paths and shared scales target the declared joint law", {
	skip_on_cran()
	new = asNamespace("lame")
	set.seed(45123)
	n = 3L
	nt = 3L
	p = 2L
	ve = .6
	rho = .5
	sig = c(.6, .9)
	lam = c(.4, 1.3)
	times = c(1, 2, 5)
	X = lapply(1:nt, function(t) {
		x = rnorm(n * n)
		cbind(x, .97 * x + .2 * rnorm(n * n))
	})
	Z = lapply(1:nt, function(t) {
		z = matrix(rnorm(n * n), n)
		diag(z) = 999
		z
	})
	Z[[2]][1, 2] = NA
	zero = lapply(Z, function(z) matrix(0, n, n))
	empty = lapply(Z, function(z) matrix(0, n * n, 0))
	# independent dense observation design and reciprocal covariance
	A = NULL
	y = NULL
	W = NULL
	for (t in 1:nt) {
		ij = which(is.finite(Z[[t]]) & row(Z[[t]]) != col(Z[[t]]))
		xy = arrayInd(ij, c(n, n))
		V = diag(ve, length(ij))
		for (a in seq_along(ij)) for (b in seq_along(ij)) if (all(xy[a, ] == rev(xy[b, ]))) V[a, b] = ve * rho
		at = matrix(0, length(ij), nt * p)
		at[, t + (0:(p - 1)) * nt] = X[[t]][ij, ]
		A = rbind(A, at)
		y = c(y, Z[[t]][ij])
		W = if (is.null(W)) solve(V) else as.matrix(Matrix::bdiag(W, solve(V)))
	}
	results = list()
	for (kind in c("matern32", "rw2", "ar1", "rw1")) {
		Q = matrix(0, nt * p, nt * p)
		h = numeric(nt * p)
		for (k in 1:p) {
			idx = 1:nt + (k - 1) * nt
			if (kind == "matern32") {
				dst = abs(outer(times, times, "-"))
				a = sqrt(3) * dst / 2
				K = solve((1 + a) * exp(-a) + diag(1e-8, nt))
			}
			if (kind == "rw2") {
				D = matrix(c(1, -4 / 3, 1 / 3) / sqrt(2), 1, 3)
				K = crossprod(D) + diag(1e-6, nt)
			}
			if (kind %in% c("ar1", "rw1")) {
				ar = if (kind == "rw1") 1 else .6
				gp = diff(times)
				D = matrix(0, nt - 1, nt)
				for (t in 1:(nt - 1)) {
					fac = ar^gp[t]
					vr = if (ar == 1) gp[t] else sum(ar^(2 * (0:(gp[t] - 1))))
					D[t, t] = -fac / sqrt(vr)
					D[t, t + 1] = 1 / sqrt(vr)
				}
				K = crossprod(D)
			}
			Q[idx, idx] = K / (sig[k]^2 * lam[k])
			if (kind %in% c("ar1", "rw1")) {
				Q[idx[1], idx[1]] = Q[idx[1], idx[1]] + 1 / 3
				h[idx[1]] = .7 / 3
			}
		}
		Q = Q + crossprod(A, W %*% A)
		h = h + crossprod(A, W %*% y)
		V = solve(Q)
		mu = drop(V %*% h)
		args = list(
			Xdyn_list = X,
			Xstat_list = empty,
			Z_list = Z,
			offset_list = zero,
			beta_static = numeric(0),
			sigma_by_coef = sig,
			Lambda = diag(lam),
			s2 = ve,
			kind = kind,
			length_scale = 2,
			time_positions = times,
			rho_by_coef = rep(.6, p),
			dyad_rho = rho,
			use_dyad_rho = TRUE,
			beta0_mean = rep(.7, p),
			beta0_cov = diag(3, p)
		)
		N = 1800L
		draw = t(replicate(N, as.vector(do.call(new$.sample_beta_dynamic_alt, args)$path)))
		zm = abs(colMeans(draw) - mu) / sqrt(diag(V) / N)
		zc = abs(cov(draw) - V) / sqrt((V^2 + outer(diag(V), diag(V))) / (N - 1))
		results[[kind]] = list(pass = max(zm, zc) < 6, new_mean_mcse = max(zm), new_cov_mcse = max(zc))
		expect_true(max(zm, zc) < 6, info = kind)
	}
	# stationarity of the shared scale prior: two coefficients share one SD
	set.seed(78932)
	N = 2500L
	for (kind in c("rw2", "matern32", "ar1", "rw1")) {
		dif = t(replicate(N, {
			v = 1 / rgamma(1, 4, 2)
			path = matrix(0, nt, p)
			for (k in 1:p) {
				if (kind == "matern32") {
					dst = abs(outer(times, times, "-"))
					a = sqrt(3) * dst / 2
					path[, k] = drop(t(chol((1 + a) * exp(-a) + diag(1e-8, nt))) %*% rnorm(nt)) * sqrt(v * lam[k])
				}
				if (kind == "rw2") {
					K = crossprod(matrix(c(1, -4 / 3, 1 / 3) / sqrt(2), 1, 3)) + diag(1e-6, nt)
					path[, k] = backsolve(chol(K), rnorm(nt)) * sqrt(v * lam[k])
				}
				if (kind %in% c("ar1", "rw1")) {
					ar = if (kind == "rw1") 1 else .6
					path[1, k] = rnorm(1, .7, sqrt(3))
					for (t in 2:nt) {
						gap = times[t] - times[t - 1]
						q = if (ar == 1) gap else sum(ar^(2 * (0:(gap - 1))))
						path[t, k] = ar^gap * path[t - 1, k] + rnorm(1, sd = sqrt(v * lam[k] * q))
					}
				}
			}
			vn = new$.sample_sigma_beta_alt(path, kind, 2, 4, 2, diag(lam), c(1, 1), 1, times, .6)^2
			c(vn - v, log(vn) - log(v))
		}))
		z = colMeans(dif) / (apply(dif, 2, sd) / sqrt(N))
		results[[paste0(kind, "_scale")]] = list(pass = max(abs(z)) < 6, max_mcse = max(abs(z)))
		expect_true(max(abs(z)) < 6, info = paste(kind, "scale"))
	}
})

test_that("pooled persistence preserves its joint hierarchy and path prior", {
	skip_on_cran()
	e = asNamespace("lame")
	set.seed(725199)
	N = 3000L
	p = 3L
	nt = 4L
	times = c(1, 2, 4, 5)
	sg = c(.5, .6, .8)
	lam = c(.6, 1.2, .9)
	feature = function(mu, tau, rho) c(
		mu = mu,
		mu2 = mu^2,
		tau2 = tau^2,
		logtau = log(tau),
		rho = mean(rho),
		rho2 = mean(rho^2),
		mu_rho = mu * mean(rho)
	)
	d = t(replicate(N, {
		mu = rnorm(1, sd = 1)
		tau = sqrt(1 / rgamma(1, 3, 1))
		rho = .999 * plogis(rnorm(p, mu, tau))
		path = matrix(0, nt, p)
		for (k in 1:p) {
			path[1, k] = rnorm(1, sd = sqrt(2))
			for (t in 2:nt) {
				gap = times[t] - times[t - 1]
				q = sum(rho[k]^(2 * (0:(gap - 1))))
				path[t, k] = rho[k]^gap * path[t - 1, k] + rnorm(1, sd = sg[k] * sqrt(lam[k] * q))
			}
		}
		before = feature(mu, tau, rho)
		h = e$.pool_update_rho_hyper(rho, prior_mu_sd = 1, prior_tau_a = 3, prior_tau_b = 1, current_tau = tau)
		rho2 = e$.pool_sample_rho(path, 1:p, p, diag(lam), sg, rho, h$mu, h$tau, 0, .999, times)
		feature(h$mu, h$tau, rho2) - before
	}))
	z = colMeans(d) / (apply(d, 2, sd) / sqrt(N))
	expect_lt(max(abs(z)), 6)
})


test_that("nodal dynamics without an intercept preserve unconstrained actor levels", {
	set.seed(9931)
	Y = lapply(1:3, function(t) {
		y = matrix(rnorm(36), 6, 6)
		diag(y) = NA_real_
		y
	})
	Xr = lapply(
		Y,
		function(y) matrix(
			rnorm(6),
			6,
			1,
			dimnames = list(NULL, "x")
		)
	)
	f = suppressMessages(lame(
		Y,
		Xrow = Xr,
		family = "normal",
		R = 0,
		intercept = FALSE,
		dynamic_beta = "row",
		dcor = FALSE,
		nscan = 20,
		burn = 10,
		odens = 1,
		seed = 678,
		posterior_opts = list(save_ab = TRUE, thin_ab = 1),
		verbose = FALSE,
		gof = FALSE
	))
	s = f$start_vals$checkpoint_state$values
	expect_gt(abs(sum(s$a)), 1e-6)
	expect_gt(abs(sum(s$b)), 1e-6)
	expect_true(all(is.finite(f$BETA)))
})
