# binary dynamic moves that read the predictor, checked with an rw1 intercept path in the predictor
skip_on_cran()

ez_pairs = function(n, rho) {
	up = which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	e1 = rnorm(nrow(up)); e2 = rho * e1 + sqrt(1 - rho^2) * rnorm(nrow(up))
	z = matrix(0, n, n); z[up] = e1; z[up[, 2:1, drop = FALSE]] = e2
	z
}
ez_icpt = function(nt, v0 = 1.4, sb = .5) cumsum(c(rnorm(1, 0, sqrt(v0)), rnorm(nt - 1L, 0, sb)))
ez_z = function(d) colMeans(d) / (apply(d, 2, sd) / sqrt(nrow(d)))
# weighted before/after difference of the statistics, z scores
ez_wz = function(D, W) {
	wn = W / sum(W); dm = colSums(D * wn); se = sqrt(colSums(wn^2 * sweep(D, 2, dm)^2))
	ifelse(se > 0, dm / se, 0)
}

test_that("the collapsed rho update keeps the prior of rho with an rw1 intercept path in the predictor", {
	set.seed(61021)
	d = NULL
	for (k in 1:8000) {
		rho = lame:::.lame_rho_prior_draw()
		path = ez_icpt(2L)
		eta = array(rnorm(32, 0, .8), c(4, 4, 2))
		for (t in 1:2) eta[, , t] = eta[, , t] + path[t]
		Z = eta + array(c(ez_pairs(4L, rho), ez_pairs(4L, rho)), c(4, 4, 2))
		Y = array(1 * (Z > 0), c(4, 4, 2)); for (t in 1:2) diag(Y[, , t]) = NA
		cx = lame:::.lame_binary_dynamic_collapsed_rho_context(Y)
		r = lame:::.lame_binary_dynamic_collapsed_rho_mh(rho, eta, cx)$rho
		d = rbind(d, c(r, r^2, abs(r), r^3) - c(rho, rho^2, abs(rho), rho^3))
	}
	expect_lt(max(abs(ez_z(d))), 4.6)
})

# directed binary panel of three actors and two periods, with a fixed rw1 intercept path and additive effects
ez_panel = function(seed) {
	set.seed(seed)
	n = 3L; nt = 2L
	path = ez_icpt(nt)
	base = array(0, c(n, n, nt))
	a = rnorm(n, 0, .4); b = rnorm(n, 0, .4)
	for (t in 1:nt) base[, , t] = path[t] + outer(a, b, '+')
	Y = array(rbinom(n * n * nt, 1, .45), c(n, n, nt)); for (t in 1:nt) diag(Y[, , t]) = NA
	obs = which(is.finite(Y)); sY = 2 * Y[obs] - 1
	list(n = n, nt = nt, base = base, Y = Y, obs = obs, sY = sY, path = path)
}

test_that("the non-centred factor scale keeps the posterior of the positions and their scale with an rw1 intercept path", {
	p = ez_panel(7101); n = p$n; nt = p$nt
	rho = .8; shape = 2; scale = .3
	stat = function(U, V, s) c(s, s^2, sum(U^2), sum(V^2), U[1, 1] * V[2, 1], U[3, 2] * V[1, 2])
	D = NULL; W = NULL
	for (k in 1:30000) {
		s = sqrt(1 / rgamma(1, shape, rate = scale))
		ar = function() {
			x = matrix(0, n, nt)
			x[, 1] = rnorm(n, 0, s / sqrt(1 - rho^2))
			x[, 2] = rho * x[, 1] + rnorm(n, 0, s)
			x
		}
		U = ar(); V = ar()
		uv = array(0, c(n, n, nt)); for (t in 1:nt) uv[, , t] = outer(U[, t], V[, t])
		eta = p$base + uv
		w = exp(sum(pnorm(p$sY * eta[p$obs], log.p = TRUE)))
		q = lame:::.lame_binary_dynamic_uv_scale_step(
			eta[p$obs],
			uv[p$obs],
			0,
			p$sY,
			s,
			shape,
			scale,
			tries = 3L
		)
		D = rbind(D, stat(q$factor * U, q$factor * V, q$sigma) - stat(U, V, s)); W = c(W, w)
	}
	expect_gt(sum(W)^2 / sum(W^2), 1000)
	expect_lt(max(abs(ez_wz(D, W))), 4.6)
})

test_that("the non-centred AR step of the additive deviations keeps its posterior with an rw1 intercept path", {
	p = ez_panel(7102); n = p$n; nt = p$nt
	m = .7; sdr = .2; shape = 2; scale = .2
	stat = function(x, r, s) c(r, r^2, s, s^2, sum(x^2), x[1, 1] * x[1, 2], x[4, 2])
	K = 250000L; D = matrix(0, K, 7); W = numeric(K)
	for (k in seq_len(K)) {
		repeat { r = rnorm(1, m, sdr); if (abs(r) < .99) break }
		s = sqrt(1 / rgamma(1, shape, rate = scale))
		x = matrix(0, 2 * n, nt)
		x[, 1] = rnorm(2 * n, 0, s / sqrt(1 - r^2))
		x[, 2] = r * x[, 1] + rnorm(2 * n, 0, s)
		like = function(paths) {
			eta = p$base
			for (t in 1:nt) {
				eta[, , t] = eta[, , t] + outer(paths[seq_len(n), t], paths[n + seq_len(n), t], '+')
			}
			list(loglik = sum(pnorm(p$sY * eta[p$obs], log.p = TRUE)), eta = eta)
		}
		w = exp(like(x)$loglik)
		z = lame:::.lame_binary_dynamic_nc_ar_step(x, r, s, like, m, sdr, shape, scale, tries = 3L)
		D[k, ] = stat(z$paths, z$rho, z$sigma) - stat(x, r, s); W[k] = w
	}
	expect_gt(sum(W)^2 / sum(W^2), 500)
	expect_lt(max(abs(ez_wz(D, W))), 4.6)
})

test_that("the marginal actor-factor slice keeps the posterior with an rw1 intercept path in the predictor", {
	bvn = .lame_binary_dynamic_bvn_log
	set.seed(921743); n = 4L; nt = 2L; rank = 2L; Sab = matrix(c(.3, .08, .08, .25), 2)
	for (side in 1:2) {
		rho = if (side == 1) .45 else -.3
		U = array(runif(n * rank * nt, -.2, .2), c(n, rank, nt))
		V = array(runif(n * rank * nt, -.2, .2), c(n, rank, nt))
		levels = matrix(rnorm(n * 2, sd = .1), n, 2); levels = sweep(levels, 2, colMeans(levels), '-')
		A = B = matrix(rnorm(n * nt, sd = .15), n, nt)
		A = sweep(A, 2, colMeans(A), '-') + levels[, 1]; B = sweep(B, 2, colMeans(B), '-') + levels[, 2]
		path = ez_icpt(nt)
		base = array(rnorm(n * n * nt, sd = .1), c(n, n, nt))
		for (t in 1:nt) base[, , t] = base[, , t] + path[t]
		predictor = function(U, V, A, B) {
			z = base
			for (tt in 1:nt) {
				uv = matrix(U[, , tt], n, rank) %*% t(matrix(V[, , tt], n, rank))
				z[, , tt] = z[, , tt] + outer(A[, tt], B[, tt], '+') + uv
			}
			z
		}
		eta = predictor(U, V, A, B)
		Y = array(NA_real_, c(n, n, nt))
		Y[2, 1, 1] = 1; Y[1, 2, 1] = 0; Y[4, 3, 2] = 1; Y[2, 4, 2] = 0
		transpose = function(z) aperm(z, c(2, 1, 3))
		pos = if (side == 1) U else V; opp = if (side == 1) V else U; mar = if (side == 1) A else B
		lev = levels[, side]; other = levels[, 3 - side]
		yy = if (side == 1) Y else transpose(Y); ez = if (side == 1) eta else transpose(eta)
		cx = .lame_binary_actor_slice_context(rank, nt, Sab, side, .3, .25, .65, .22)
		ss = .lame_binary_actor_slice_system(pos, opp, mar, lev, other, yy, ez, c(2L, 4L), rho, cx)
		loglik = function(x) .lame_binary_actor_slice_loglik(
			ss$offset1 + drop(ss$D1 %*% x),
			ss$offset2 + drop(ss$D2 %*% x),
			ss,
			bvn
		)
		N = 6000L; delta = matrix(0, N, 7); weights = numeric(N)
		summ = function(x) c(
			x[1],
			mean(x[cx$ix$dev]),
			mean(x[cx$ix$first]),
			mean(x[cx$ix$second]),
			mean(x[cx$ix$first]^2),
			mean(x[cx$ix$second]^2),
			x[1]^2
		)
		for (k in 1:N) {
			xx = ss$prior_mean + drop(cx$root %*% rnorm(cx$p))
			ss$current = xx
			weights[k] = exp(loglik(xx))
			z = .lame_binary_actor_slice_step(ss, max_abs = Inf, bvn = bvn)
			delta[k, ] = (summ(z$x) - summ(xx)) * weights[k]
		}
		zs = colMeans(delta) / (apply(delta, 2, sd) / sqrt(N))
		expect_gt(sum(weights)^2 / sum(weights^2), 1000)
		expect_lt(max(abs(zs)), 4.5)
	}
})
