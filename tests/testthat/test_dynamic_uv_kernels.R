# dynamic position kernels against an independently assembled full conditional: for i != j the pair
# (e_ij, e_ji) is bivariate normal (variance s2, correlation rho), e_ii has variance s2 (1 + rho), e = E - U_t V_t',
# and positions follow a stationary AR(1) with transitions N(rho_uv x_{t-1}, sigma_uv^2 / lambda I)
# each kernel is called with shrink = FALSE, so the returned row is the draw
skip_on_cran()

kern_pars = list(rho_uv = 0.7, sigma_uv = 0.5, s2 = 1.3, kappa = 2, pi_snap = 0.2, nu = 4)

kern_state = function(n, R, TT, seed, sym = FALSE) {
	set.seed(seed)
	U = array(rnorm(n * R * TT, 0, 0.8), c(n, R, TT))
	V = if (sym) U else array(rnorm(n * R * TT, 0, 0.8), c(n, R, TT))
	for (t in seq_len(TT)) {
		U[, 2, t] = 0.6 * U[, 1, t] + 0.8 * U[, 2, t]
		V[, 2, t] = if (sym) U[, 2, t] else 0.5 * V[, 1, t] + 0.85 * V[, 2, t]
	}
	E = array(0, c(n, n, TT))
	for (t in seq_len(TT)) {
		z1 = matrix(rnorm(n * n), n, n)
		z2 = matrix(rnorm(n * n), n, n)
		e = z1
		e[lower.tri(e)] = (0.4 * t(z1) + sqrt(1 - 0.16) * z2)[lower.tri(e)]
		m = U[, , t] %*% t(V[, , t]) + 0.3 + e
		if (sym) m[lower.tri(m)] = t(m)[lower.tri(m)]
		E[, , t] = m
	}
	list(U = U, V = V, E = E)
}

# data part of row i of U at period t, one term per pair and one for the self cell
kern_data_u = function(i, t, U, V, E, s2, rho) {
	n = dim(U)[1]; R = dim(U)[2]
	Ut = U[, , t]; Vt = V[, , t]; Et = E[, , t]
	P = matrix(0, R, R); h = numeric(R)
	cp = s2 * (1 - rho^2)
	for (j in seq_len(n)[-i]) {
		y = Et[j, i] - sum(Ut[j, ] * Vt[i, ])
		P = P + tcrossprod(Vt[j, ]) / cp
		h = h + (Et[i, j] - rho * y) * Vt[j, ] / cp
	}
	P = P + tcrossprod(Vt[i, ]) / (s2 * (1 + rho))
	h = h + Et[i, i] * Vt[i, ] / (s2 * (1 + rho))
	list(P = P, h = h)
}

kern_data_v = function(j, t, U, V, E, s2, rho) {
	n = dim(U)[1]; R = dim(U)[2]
	Ut = U[, , t]; Vt = V[, , t]; Et = E[, , t]
	P = matrix(0, R, R); h = numeric(R)
	cp = s2 * (1 - rho^2)
	for (i in seq_len(n)[-j]) {
		y = Et[j, i] - sum(Ut[j, ] * Vt[i, ])
		P = P + tcrossprod(Ut[i, ]) / cp
		h = h + (Et[i, j] - rho * y) * Ut[i, ] / cp
	}
	P = P + tcrossprod(Ut[j, ]) / (s2 * (1 + rho))
	h = h + Et[j, j] * Ut[j, ] / (s2 * (1 + rho))
	list(P = P, h = h)
}

kern_prior = function(i, t, X, lam = NULL, fwd = TRUE) {
	TT = dim(X)[3]; R = dim(X)[2]
	ru = kern_pars$rho_uv; s2u = kern_pars$sigma_uv^2
	lw = function(tt) if (is.null(lam)) 1 else lam[i, tt]
	P = matrix(0, R, R); h = numeric(R)
	if (t == 1) {
		P = P + diag((1 - ru^2) / s2u, R)
	} else {
		P = P + diag(lw(t) / s2u, R)
		h = h + lw(t) * ru * X[i, , t - 1] / s2u
	}
	if (t < TT && fwd) {
		P = P + diag(lw(t + 1) * ru^2 / s2u, R)
		h = h + lw(t + 1) * ru * X[i, , t + 1] / s2u
	}
	list(P = P, h = h)
}

kern_cond = function(a, b) {
	P = a$P + b$P
	list(P = P, mean = drop(solve(P, a$h + b$h)), cov = solve(P))
}

# kernel draws of u_1 at period 1 against the analytic mean and covariance
kern_check_first = function(draw_fun, cg, N, seed0) {
	X = t(vapply(seq_len(N), function(k) { set.seed(seed0 + k); draw_fun() }, numeric(length(cg$mean))))
	S = cg$cov
	z_mean = (colMeans(X) - cg$mean) / sqrt(diag(S) / N)
	z_cov = (cov(X) - S) / sqrt((outer(diag(S), diag(S)) + S^2) / N)
	expect_lt(max(abs(z_mean)), 4.5)
	expect_lt(max(abs(z_cov)), 4.5)
	invisible(list(z_mean = z_mean, z_cov = z_cov))
}

test_that("ar1 kernel draws the first row from the exact pair-model conditional", {
	N = 4000
	for (rho in c(0, 0.4, 0.7)) {
		st = kern_state(6, 2, 3, seed = 11)
		dp = kern_data_u(1, 1, st$U, st$V, st$E, kern_pars$s2, rho)
		cg = kern_cond(dp, kern_prior(1, 1, st$U))
		# the assembly is the log density's quadratic: central differences agree
		ld = function(u) {
			U = st$U; U[1, , 1] = u
			res = st$E[, , 1] - U[, , 1] %*% t(st$V[, , 1])
			S = kern_pars$s2 * matrix(c(1, rho, rho, 1), 2)
			ll = sum(vapply(2:6, function(j) {
				q = c(res[1, j], res[j, 1]); -0.5 * sum(q * solve(S, q))
			}, 0)) - 0.5 * res[1, 1]^2 / (kern_pars$s2 * (1 + rho))
			ll + sum(dnorm(u, 0, kern_pars$sigma_uv / sqrt(1 - kern_pars$rho_uv^2), log = TRUE)) +
				sum(dnorm(U[1, , 2], kern_pars$rho_uv * u, kern_pars$sigma_uv, log = TRUE))
		}
		eps = 1e-3
		g = vapply(1:2, function(a) (ld(eps * (1:2 == a)) - ld(-eps * (1:2 == a))) / (2 * eps), 0)
		expect_equal(g, dp$h + kern_prior(1, 1, st$U)$h, tolerance = 1e-6)
		if (rho == 0) {
			# no reciprocal term: the independent-dyad regression
			expect_equal(dp$P, crossprod(st$V[, , 1]) / kern_pars$s2, tolerance = 1e-12)
			expect_equal(dp$h, drop(st$E[1, , 1] %*% st$V[, , 1]) / kern_pars$s2, tolerance = 1e-12)
		}
		kern_check_first(
			function() rUV_dynamic_fc(
				st$U,
				st$V,
				st$E,
				kern_pars$rho_uv,
				kern_pars$sigma_uv,
				kern_pars$s2,
				shrink = FALSE,
				symmetric = FALSE,
				rho = rho
			)$U[1, , 1],
			cg,
			N,
			1e6 * (1 + 10 * rho)
		)
	}
})

test_that("ar1 kernel's correlated path reduces to the independent-dyad path at rho = 0", {
	st = kern_state(6, 2, 3, seed = 12)
	set.seed(3)
	a = rUV_dynamic_fc(st$U, st$V, st$E, 0.7, 0.5, 1.3, shrink = FALSE, rho = 0)
	set.seed(3)
	b = rUV_dynamic_fc(st$U, st$V, st$E, 0.7, 0.5, 1.3, shrink = FALSE, rho = 1e-9)
	expect_lt(max(abs(a$U - b$U), abs(a$V - b$V)), 1e-6)
})

test_that("every conditional of an ar1 sweep with rho = 0.7 is exact", {
	skip_on_cran()
	# each drawn row standardised by its own conditional given the state the
	# kernel had reached (earlier rows and periods at their new values)
	N = 1200; rho = 0.7; n = 6; TT = 3; R = 2
	st = kern_state(n, R, TT, seed = 13)
	Z = matrix(0, N, 2 * n * TT * R)
	for (k in seq_len(N)) {
		set.seed(2e6 + k)
		fin = rUV_dynamic_fc(
			st$U,
			st$V,
			st$E,
			kern_pars$rho_uv,
			kern_pars$sigma_uv,
			kern_pars$s2,
			shrink = FALSE,
			rho = rho
		)
		col = 0
		for (t in seq_len(TT)) {
			U = st$U; V = st$V
			if (t > 1) { U[, , 1:(t - 1)] = fin$U[, , 1:(t - 1)]; V[, , 1:(t - 1)] = fin$V[, , 1:(t - 1)] }
			for (i in seq_len(n)) {
				if (i > 1) U[i - 1, , t] = fin$U[i - 1, , t]
				cg = kern_cond(kern_data_u(i, t, U, V, st$E, kern_pars$s2, rho), kern_prior(i, t, U))
				Z[k, col + 1:R] = drop(chol(cg$P) %*% (fin$U[i, , t] - cg$mean)); col = col + R
			}
			U[n, , t] = fin$U[n, , t]
			for (j in seq_len(n)) {
				if (j > 1) V[j - 1, , t] = fin$V[j - 1, , t]
				cg = kern_cond(kern_data_v(j, t, U, V, st$E, kern_pars$s2, rho), kern_prior(j, t, V))
				Z[k, col + 1:R] = drop(chol(cg$P) %*% (fin$V[j, , t] - cg$mean)); col = col + R
			}
		}
	}
	# 72 coordinates: the largest of 72 standard normals rarely passes 4.5
	expect_lt(max(abs(colMeans(Z) * sqrt(N))), 4.5)
	expect_lt(max(abs((apply(Z, 2, var) - 1) / sqrt(2 / N))), 4.5)
})

test_that("t kernel weights the transitions by the local precisions", {
	st = kern_state(6, 2, 3, seed = 14)
	lu = matrix(1, 6, 3); lv = matrix(1, 6, 3)
	lu[1, 2] = 3
	rho = 0.4
	cg = kern_cond(kern_data_u(1, 1, st$U, st$V, st$E, kern_pars$s2, rho), kern_prior(1, 1, st$U, lam = lu))
	kern_check_first(
		function() rUV_dynamic_t_fc(
			st$U,
			st$V,
			st$E,
			kern_pars$rho_uv,
			kern_pars$sigma_uv,
			kern_pars$s2,
			kern_pars$nu,
			lu,
			lv,
			shrink = FALSE,
			rho = rho
		)$U[1, , 1],
		cg,
		4000,
		3e6
	)
})

test_that("snap kernel cuts the forward link after a reset and weighs the indicator exactly", {
	rho = 0.4
	st = kern_state(6, 2, 3, seed = 15)
	for (reset in c(0, 1)) {
		du = matrix(0, 6, 3); du[1, 2] = reset
		cg = kern_cond(
			kern_data_u(1, 1, st$U, st$V, st$E, kern_pars$s2, rho),
			kern_prior(1, 1, st$U, fwd = reset == 0)
		)
		kern_check_first(
			function() rUV_dynamic_snap_fc(
				st$U,
				st$V,
				st$E,
				kern_pars$rho_uv,
				kern_pars$sigma_uv,
				kern_pars$s2,
				kern_pars$kappa,
				kern_pars$pi_snap,
				du,
				matrix(0, 6, 3),
				shrink = FALSE,
				rho = rho
			)$U[1, , 1],
			cg,
			3000,
			4e6 + 1e5 * reset
		)
	}
	# the indicator of actor 1 at period 2 with the next period drifting: the
	# exact odds carry the density of x_3 with x_2 integrated out under each branch
	st = kern_state(6, 2, 3, seed = 16)
	rho = 0.7; N = 3000
	ru = kern_pars$rho_uv; s2u = kern_pars$sigma_uv^2; k2 = kern_pars$kappa^2
	d_k = numeric(N); p_ex = numeric(N)
	for (k in seq_len(N)) {
		set.seed(5e6 + k)
		fin = rUV_dynamic_snap_fc(
			st$U,
			st$V,
			st$E,
			ru,
			kern_pars$sigma_uv,
			kern_pars$s2,
			kern_pars$kappa,
			kern_pars$pi_snap,
			shrink = FALSE,
			rho = rho
		)
		U = st$U; V = st$V
		U[, , 1] = fin$U[, , 1]; V[, , 1] = fin$V[, , 1]
		dp = kern_data_u(1, 2, U, V, st$E, kern_pars$s2, rho)
		xp = fin$U[1, , 1]; xn = st$U[1, , 3]
		log_m = function(Pp, hp, cst) {
			P = dp$P + Pp; h = dp$h + hp
			cst - 0.5 * as.numeric(determinant(P)$modulus) + 0.5 * sum(h * solve(P, h))
		}
		fwd_c = -log(2 * pi * s2u) - 0.5 * sum(xn^2) / s2u
		l0 = log_m(
			diag((1 + ru^2) / s2u, 2),
			ru * (xp + xn) / s2u,
			-log(2 * pi * s2u) - 0.5 * sum((ru * xp)^2) / s2u + fwd_c
		)
		l1 = log_m(diag(1 / k2 + ru^2 / s2u, 2), ru * xn / s2u, -log(2 * pi * k2) + fwd_c)
		p_ex[k] = plogis(qlogis(kern_pars$pi_snap) + l1 - l0)
		d_k[k] = fin$delta_u[1, 2]
	}
	z = (sum(d_k) - sum(p_ex)) / sqrt(sum(p_ex * (1 - p_ex)))
	expect_lt(abs(z), 4.5)
})

test_that("symmetric ar1 kernel is a gibbs sweep of the product u_t l u_t'", {
	skip_on_cran()
	# every row of every period standardised by the exact conditional of the
	# symmetric model given the current coordinates of the other actors (rows
	# already drawn in the sweep at their new values), off-diagonal cells only
	N = 1200; n = 6; TT = 3; R = 2; L = c(1.3, -0.7)
	st = kern_state(n, R, TT, seed = 17, sym = TRUE)
	for (t in 1:TT) st$V[, , t] = st$U[, , t] %*% diag(L)
	Z = matrix(0, N, n * TT * R)
	for (k in seq_len(N)) {
		set.seed(6e6 + k)
		# rho is not used by a symmetric fit
		fin = rUV_dynamic_fc(
			st$U,
			st$V,
			st$E,
			kern_pars$rho_uv,
			kern_pars$sigma_uv,
			kern_pars$s2,
			shrink = FALSE,
			symmetric = TRUE,
			rho = 0.5
		)
		U = st$U
		col = 0
		for (t in seq_len(TT)) {
			if (t > 1) U[, , t - 1] = fin$U[, , t - 1]
			for (i in seq_len(n)) {
				if (i > 1) U[i - 1, , t] = fin$U[i - 1, , t]
				Vo = (U[, , t] %*% diag(L))[-i, , drop = FALSE]
				cg = kern_cond(
					list(P = crossprod(Vo) / kern_pars$s2, h = drop(st$E[i, -i, t] %*% Vo) / kern_pars$s2),
					kern_prior(i, t, U)
				)
				Z[k, col + 1:R] = drop(chol(cg$P) %*% (fin$U[i, , t] - cg$mean)); col = col + R
			}
			U[n, , t] = fin$U[n, , t]
		}
		# the kernel hands back the coordinates in both slots
		if (k == 1) expect_identical(fin$V, fin$U)
	}
	expect_lt(max(abs(colMeans(Z) * sqrt(N))), 4.5)
	expect_lt(max(abs((apply(Z, 2, var) - 1) / sqrt(2 / N))), 4.5)
})

test_that("the snap and t kernels refuse a symmetric network", {
	st = kern_state(6, 2, 3, seed = 18, sym = TRUE)
	expect_error(
		rUV_dynamic_snap_fc(st$U, st$V, st$E, 0.7, 0.5, 1.3, 2, 0.2, symmetric = TRUE),
		class = "lame_symmetric_kind"
	)
	expect_error(
		rUV_dynamic_t_fc(st$U, st$V, st$E, 0.7, 0.5, 1.3, 4, symmetric = TRUE),
		class = "lame_symmetric_kind"
	)
})
