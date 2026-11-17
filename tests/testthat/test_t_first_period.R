# the first period of the t kind: a Student-t of scale sigma_uv / sqrt(1 - rho_uv^2), drawn as a gaussian
# with its own weight (src/dynamic_latent_factors.cpp, R/uv_rho_w_step.R, R/uv_split_dynamic.R,
# R/uv_product_scale.R, R/dynamic_position_move.R)
skip_on_cran()

# z of a mean against its expectation, the Monte Carlo error from the draws
tfp_z = function(x, m) (mean(x) - m) / (sd(x) / sqrt(length(x)))

test_that('the t kernel draws the first period with its own weight', {
	n = 6L; R = 2L; N = 3L
	rho_uv = 0.8; sigma_uv = 0.5; s2 = 4; nu = 4
	set.seed(31)
	U = array(rnorm(n * R * N, 0, 0.6), c(n, R, N))
	V = array(rnorm(n * R * N, 0, 0.6), c(n, R, N))
	E = array(rnorm(n * n * N, 0, 1), c(n, n, N))
	lu = matrix(rgamma(n * N, 2, 2), n, N); lv = matrix(rgamma(n * N, 2, 2), n, N)
	lu[1, 1] = 0.3; lu[1, 2] = 1.7
	# conditional of the first row of the first period: data V_1' V_1 / s2, the first period's weight and
	# the transition into period 2
	cond = function(l1) {
		P = crossprod(V[, , 1]) / s2 + diag(l1 * (1 - rho_uv^2) / sigma_uv^2 + lu[1, 2] * rho_uv^2 / sigma_uv^2, R)
		h = drop(E[1, , 1] %*% V[, , 1]) / s2 + lu[1, 2] * rho_uv * U[1, , 2] / sigma_uv^2
		S = solve(P)
		list(mean = drop(S %*% h), cov = S)
	}
	draws = function(lam_u) {
		t(vapply(seq_len(4000), function(k) {
			set.seed(5e5 + k)
			o = lame:::rUV_dynamic_t_fc_cpp(U, V, E, rho_uv, sigma_uv, s2, nu, lam_u, lv, FALSE, FALSE)
			u = o$U[1, , 1]
			c(u, o$lambda_u[1, 1] * (nu + sum(u^2) * (1 - rho_uv^2) / sigma_uv^2) / (nu + R))
		}, numeric(R + 1L)))
	}
	X = draws(lu)
	cg = cond(lu[1, 1])
	zm = (colMeans(X[, 1:R]) - cg$mean) / sqrt(diag(cg$cov) / nrow(X))
	zv = (apply(X[, 1:R], 2, var) - diag(cg$cov)) / (diag(cg$cov) * sqrt(2 / nrow(X)))
	expect_lt(max(abs(c(zm, zv))), 4.5)
	# the weight given the drawn position is Gamma((nu + R) / 2, (nu + q) / 2): the scaled weight has mean one
	expect_lt(abs(tfp_z(X[, R + 1L], 1)), 4.5)
	# a kernel with the gaussian first period (weight one) fails the position check, and the weight kept as it
	# came in fails the weight check
	X0 = draws(replace(lu, 1L, 1))
	zm0 = (colMeans(X0[, 1:R]) - cg$mean) / sqrt(diag(cg$cov) / nrow(X0))
	zv0 = (apply(X0[, 1:R], 2, var) - diag(cg$cov)) / (diag(cg$cov) * sqrt(2 / nrow(X0)))
	expect_gt(max(abs(c(zm0, zv0))), 5)
	kept = lu[1, 1] * (nu + rowSums(X[, 1:R]^2) * (1 - rho_uv^2) / sigma_uv^2) / (nu + R)
	expect_gt(abs(tfp_z(kept, 1)), 5)
})

# paths of the t kind and their weights: lam[, 1] the first period's, lam[, t] the transition into t
tfp_paths = function(n, R, N, rho, sigma, nu = 4) {
	lam = matrix(rgamma(n * N, nu / 2, nu / 2), n, N)
	X = array(0, c(n, R, N))
	X[, , 1] = matrix(rnorm(n * R), n, R) * sigma / sqrt((1 - rho^2) * lam[, 1])
	for (t in 2:N) X[, , t] = rho * X[, , t - 1] + matrix(rnorm(n * R), n, R) * sigma / sqrt(lam[, t])
	list(X = X, lam = lam)
}

test_that('the persistence and stationary variance read the first period with its weight', {
	set.seed(41)
	N = 4L; rho = 0.95; sig = 0.4
	a = tfp_paths(25L, 1L, N, rho, sig); b = tfp_paths(25L, 1L, N, rho, sig)
	p = list(U = a$X, V = b$X)
	# exact conditional of rho by a grid with w integrated out, the prior N(0.9, 0.1) and IG(2, 1) on sigma^2
	grid_mean = function(lu, lv) {
		S1 = Stt = Sxy = Syy = 0; n1 = ntr = 0
		for (k in 1:2) {
			X = list(p$U, p$V)[[k]]; lam = list(lu, lv)[[k]]
			S1 = S1 + sum(lam[, 1] * X[, 1, 1]^2); n1 = n1 + dim(X)[1]
			for (t in 2:N) {
				Stt = Stt + sum(lam[, t] * X[, 1, t]^2); Sxy = Sxy + sum(lam[, t] * X[, 1, t] * X[, 1, t - 1])
				Syy = Syy + sum(lam[, t] * X[, 1, t - 1]^2); ntr = ntr + dim(X)[1]
			}
		}
		rg = seq(-0.989, 0.989, length.out = 3001)
		sw = (n1 + ntr) / 2 + 2
		lw = vapply(rg, function(r) {
			om = 1 - r^2
			rate = (S1 + (Stt - 2 * r * Sxy + r^2 * Syy + 2) / om) / 2
			-(ntr / 2 + 2) * log(om) + dnorm(r, 0.9, 0.1, log = TRUE) - sw * log(rate)
		}, 0)
		wt = exp(lw - max(lw)); wt = wt / sum(wt)
		c(mean = sum(wt * rg), sd = sqrt(sum(wt * rg^2) - sum(wt * rg)^2))
	}
	g = grid_mean(a$lam, b$lam)
	run = function(lu, lv, n = 12000L) {
		set.seed(7)
		s = c(0.9, 0.3); x = numeric(n)
		for (i in seq_len(n)) {
			o = lame:::.lame_uv_rho_w_step(p$U, p$V, s[1], s[2], 0.9, 0.1, 2, 1, lambda_u = lu, lambda_v = lv)
			s = c(o$rho, o$sigma); x[i] = s[1]
		}
		x[-(1:500)]
	}
	ess = function(x) { a = acf(x, lag.max = 200, plot = FALSE)$acf[-1]; k = which(a < 0.05)[1]; length(x) / (1 + 2 * sum(a[seq_len(k - 1L)])) }
	x = run(a$lam, b$lam)
	z = (mean(x) - g[["mean"]]) / (sd(x) / sqrt(ess(x)))
	expect_lt(abs(z), 4)
	expect_lt(abs(sd(x) / g[["sd"]] - 1), 0.12)
	# the step that reads the first period without its weight targets another conditional
	x0 = run(replace(a$lam, cbind(seq_len(25L), 1L), 1), replace(b$lam, cbind(seq_len(25L), 1L), 1))
	expect_gt(abs(mean(x0) - g[["mean"]]) / (sd(x0) / sqrt(ess(x0))), 5)
})

test_that('the kind weights carry the first period, and the split and gauge keep the t prior', {
	lu = matrix(c(0.3, 2, 1.1, 0.7, 1.4, 0.9), 2, 3); lv = lu[2:1, ]
	kw = lame:::.lame_uv_dyn_kind_weights('t', lu, lv)
	expect_identical(kw$w_u, lu)
	expect_identical(kw$w_v, lv)
	expect_null(kw$k_u)
	# the paths stay draws of the t prior after the split and the gauge when the weights carry the first period
	law = function(first_one, ndraw = 4000L, seed = 51L) {
		set.seed(seed)
		n = 3L; R = 2L; N = 2L; rho = 0.9; sigma = 0.3
		stat = function(U, V) {
			x = c(log(sum(U[, , 1]^2)), log(sum(V[, , 1]^2)), log(sum(U[, , 2]^2)), log(sum(U[, 1, ]^2) / sum(U[, 2, ]^2)))
			c(x, x^2)
		}
		out = matrix(0, ndraw, 8L)
		for (r in seq_len(ndraw)) {
			a = tfp_paths(n, R, N, rho, sigma); b = tfp_paths(n, R, N, rho, sigma)
			kw = lame:::.lame_uv_dyn_kind_weights('t', a$lam, b$lam)
			if (first_one) { kw$w_u[, 1] = 1; kw$w_v[, 1] = 1 }
			U = a$X; V = b$X
			for (k in 1:2) {
				o = lame:::.lame_uv_split_dynamic_draw(U, V, rho, sigma, kw$w_u, kw$w_v)
				o = lame:::.lame_uv_gauge_dynamic_draw(o$U, o$V, rho, sigma, kw$w_u, kw$w_v)
				U = o$U; V = o$V
			}
			out[r, ] = stat(U, V) - stat(a$X, b$X)
		}
		max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))))
	}
	expect_lt(law(FALSE), 4.5)
	expect_gt(law(TRUE, seed = 52L), 5)
})

# one joint draw of a directed t-kind state with its panel, as in test_uv_product_scale.R
tfp_state = function(family, n = 6L, R = 2L, N = 3L, rho_uv = 0.85, shape = 3, scale = 2, s2 = 0.3, nu = 4) {
	sigma = sqrt(1 / rgamma(1, shape, rate = scale))
	a = tfp_paths(n, R, N, rho_uv, sigma, nu); b = tfp_paths(n, R, N, rho_uv, sigma, nu)
	base = array(-0.3 + rnorm(n * n * N, 0, 0.3), c(n, n, N))
	EZ = base + lame:::.lame_uvps_product(a$X, b$X)
	if (family == 'poisson') {
		Z = EZ + rnorm(length(EZ), 0, sqrt(s2))
		Y = array(rpois(length(Z), exp(Z)), dim(Z))
	} else {
		Z = EZ + rnorm(length(EZ))
		Y = 1 * (Z > 0)
	}
	for (t in 1:N) diag(Y[, , t]) = NA
	list(Y = Y, Z = Z, base = base, U = a$X, V = b$X, sigma = sigma, lam_u = a$lam, lam_v = b$lam)
}

# the size visit keeps the posterior: two visits from a joint draw leave the law of the state unchanged
tfp_size_law = function(ndraw, seed, family) {
	set.seed(seed)
	stat = function(st) {
		P = lame:::.lame_uvps_product(st$U, st$V)
		x = c(log(apply(P^2, 3, sum)), log(apply(st$U^2, 3, sum) / apply(st$V^2, 3, sum)), log(st$sigma),
			log(sum(st$lam_u[, 1])), log(sum(st$lam_v[, 1])))
		c(x, x^2)
	}
	out = NULL
	for (r in seq_len(ndraw)) {
		st = tfp_state(family)
		cx = lame:::.lame_uvps_context(st$Y, family, 'dir')
		o = st
		for (k in 1:2) {
			v = lame:::.lame_uvps_visit(
				cx, o$Z, o$base + lame:::.lame_uvps_product(o$U, o$V), o$Y, o$U, o$V, 0.85, o$sigma, 3, 2,
				kind = 't', lambda_u = o$lam_u, lambda_v = o$lam_v, nu = 4
			)
			o$Z = v$Z; o$U = v$U; o$V = v$V; o$sigma = v$sigma; o$lam_u = v$lambda_u; o$lam_v = v$lambda_v
		}
		s1 = stat(o) - stat(st)
		if (is.null(out)) out = matrix(0, ndraw, length(s1))
		out[r, ] = s1
	}
	max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(ndraw))))
}

test_that('the size visit integrates the first-period weight with the transition weights', {
	expect_lt(tfp_size_law(1500, 9301, 'binary'), 4.5)
	expect_lt(tfp_size_law(1500, 9302, 'poisson'), 4.5)
	# the visit with the gaussian first period fails the check
	terms = lame:::.lame_uvps_t_terms
	gauss_first = function(X, t, rho, s2, nu) {
		f = terms(X, t, rho, s2, nu)
		if (t != 1L) return(f)
		sq = rowSums(matrix(X[, , 1L], dim(X)[1L])^2)
		tt = function(s) -(nu + dim(X)[2L]) / 2 * sum(log1p((1 - rho^2) * s^2 * sq / (nu * s2)))
		function(s) f(s) - tt(s) - (1 - rho^2) * s^2 * sum(sq) / (2 * s2)
	}
	local_mocked_bindings(.lame_uvps_t_terms = gauss_first, .package = 'lame')
	expect_gt(tfp_size_law(1500, 9303, 'binary'), 5)
})

test_that('the path covariance of the position move carries the first-period weight', {
	nt = 4; rho = 0.9; sg = 0.5
	lam = c(0.25, 1, 1, 1)
	cv = lame:::.lame_dpm_cov(nt, rho, sg, lam)
	expect_equal(cv[1, 1], sg^2 / ((1 - rho^2) * 0.25))
	expect_equal(cv[1, 4], rho^3 * cv[1, 1])
	expect_equal(lame:::.lame_dpm_cov(nt, rho, sg, rep(1, nt)), lame:::.lame_dpm_cov(nt, rho, sg))
})

test_that('lame() draws the first-period weights of a t fit and carries them across a resume', {
	set.seed(61)
	n = 10L; N = 4L
	nm = sprintf('a%02d', 1:n)
	Y = lapply(1:N, function(t) {
		y = matrix(rbinom(n * n, 1, 0.3), n, n, dimnames = list(nm, nm))
		diag(y) = NA
		y
	})
	names(Y) = paste0('t', 1:N)
	a = list(Y = Y, family = 'binary', R = 2, dynamic_uv = TRUE, dynamic_uv_kind = 't', burn = 20, nscan = 40,
		odens = 2, chains = 1, verbose = FALSE, plot = FALSE, prior = list(als_start = FALSE))
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	lam = f$start_vals$checkpoint_state$values$lambda_u
	expect_equal(dim(lam), c(n, N))
	expect_true(all(is.finite(lam) & lam > 0))
	expect_false(all(lam[, 1] == 1))
	expect_true(all(is.finite(f$rho_uv)) && all(is.finite(f$sigma_uv)))
	a2 = a; a2$start_vals = f$start_vals; a2$burn = 0; a2$nscan = 10; a2$odens = 1
	f2 = suppressWarnings(suppressMessages(do.call(lame, a2)))
	expect_true(all(is.finite(f2$rho_uv)))
})
