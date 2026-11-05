test_that("the bipartite joint trace includes the static interaction information", {
	set.seed(80019);nr = 3L;nc = 2L;nt = 3L;ra = 1L;rb = 2L;u = array(
		rnorm(nr * ra * nt),
		c(nr, ra, nt)
	);v = array(
		rnorm(nc * rb * nt),
		c(nc, rb, nt)
	);g = matrix(rnorm(ra * rb), ra, rb);lambda = .8;rho = .7;level = lame:::.lda_level_prec
	block = nr * ra + nc * rb;ng = ra * rb;ig = block * nt + seq_len(ng)
	theta = c(
		unlist(
			lapply(1:nt, function(t) c(as.vector(t(matrix(u[, , t], nr, ra))), as.vector(t(matrix(v[, , t], nc, rb)))))
		),
		c(g)
	)
	prediction = function(theta) {gg = matrix(theta[ig], ra, rb);unlist(
		lapply(
			1:nt,
			function(t) {idx = (t - 1) * block;uu = matrix(
				theta[idx + seq_len(nr * ra)],
				nr,
				ra,
				byrow = TRUE
			);vv = matrix(
				theta[idx + nr * ra + seq_len(nc * rb)],
				nc,
				rb,
				byrow = TRUE
			);uu %*% gg %*% t(vv)}
		)
	)}
	h = 1e-5;J = vapply(
		seq_along(theta),
		function(k) {e = numeric(length(theta));e[k] = h;(prediction(theta + e) - prediction(theta - e)) / (2 * h)},
		numeric(nr * nc * nt)
	);K = crossprod(J)
	D = matrix(0, nt - 1, nt);for (t in 2:nt) {D[t - 1, t] = 1;D[t - 1, t - 1] = -rho};temporal = lambda * crossprod(D) + diag(level + 1e-6, nt)
	P = matrix(
		0,
		length(theta),
		length(theta)
	);for (k in 1:block) {ix = k + block * (0:(nt - 1));P[ix, ix] = temporal};P[ig, ig] = diag(ng)
	expected = sum(diag(solve(K + P, K)))
	observed = lame:::.lda_bip_joint_trace(
		u,
		v,
		g,
		rep(list(matrix(0, nr, nc)), nt),
		array(TRUE, c(nr, nc, nt)),
		array(1, c(nr, nc, nt)),
		lambda,
		rho
	)
	expect_equal(observed, expected, tolerance = 1e-7)
})

test_that("a stalled dynamic Gaussian fit receives a certified joint refinement", {
	skip_on_cran()
	# the sweep limit is set below what the fit needs, so the refinement finishes what the block sweeps leave
	set.seed(1574001L);nr = 30L;nc = 20L;nt = 8L;ranks = c(2L, 2L)
	paths = function(n, r, sd, rho) {x = array(0, c(n, r, nt));x[, , 1] = rnorm(
		n * r,
		0,
		sd / sqrt(1 - rho^2)
	);for (t in 2:nt) x[, , t] = rho * x[, , t - 1] + rnorm(n * r, 0, sd);x}
	u = paths(nr, 2, .4, .7);v = paths(nc, 2, .4, .7);g = diag(2)
	a = matrix(rep(rnorm(nr, 0, .4), nt), nr) + matrix(paths(nr, 1, .15, .7), nr, nt);b = matrix(rep(rnorm(nc, 0, .4), nt), nc) + matrix(paths(nc, 1, .15, .7), nc, nt)
	Y = X = vector('list', nt)
	for (t in 1:nt) {xx = matrix(rnorm(nr * nc), nr);yy = .2 + .5 * xx + outer(a[, t], b[, t], '+') + u[, , t] %*% g %*% t(v[, , t]) + matrix(rnorm(nr * nc), nr);dimnames(yy) = list(
		paste0('r', 1:nr),
		paste0('c', 1:nc)
	);Y[[t]] = yy;X[[t]] = array(
		xx,
		c(nr, nc, 1),
		dimnames = c(dimnames(yy), list('x'))
	)}
	f = lame(
		Y,
		Xdyad = X,
		R = 2,
		R_row = 2,
		R_col = 2,
		family = 'normal',
		mode = 'bipartite',
		method = 'als',
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		dynamic_beta = FALSE,
		als_max_iter = 100L,
		als_tol = 1e-6,
		verbose = FALSE,
		bias_correct = 'none'
	)
	expect_true(f$converged);expect_true(f$diagnostics$stationarity$pass)
	refinement = f$diagnostics$joint_refinement
	expect_true(refinement$accepted);expect_lte(refinement$objective_after, refinement$objective_before)
	expect_lt(unname(f$VC['ve']), 1.1);expect_gt(unname(f$VC['ve']), .9)
	expect_equal(tail(f$objective_history, 1), refinement$objective_after, tolerance = 1e-10)
})
