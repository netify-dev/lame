# the ridge of a row or column covariate with the additive effects of its margin (R/covariate_ridge_move.R)

# stationarity on a 6-actor panel: draw a state from the prior, apply one visit, and every statistic keeps its mean
# mode: 'dyn' (dynamic effects, absent actor-periods zero), 'static' or 'ame' (no centring); returns the largest z score
cr_exact = function(mode, ndraw = 8000L, seed = 21L, mutant = 'none') {
	set.seed(seed)
	n = 6L; N = if (mode == 'ame') 1L else 3L; center = mode != 'ame'; dyn = mode == 'dyn'
	rho_ab = .5; sigma_ab = .7; rho = .4; s2 = 1.3; p = 4L
	Sab = matrix(c(1, .4, .4, 1.5), 2)
	P0 = crossprod(matrix(rnorm(p * p), p)) / 3 + diag(.5, p)
	xr = matrix(rnorm(n * N), n, N); xc = matrix(rnorm(n * N), n, N); xd = array(rnorm(n * n * N), c(n, n, N))
	if (!dyn) { xr = matrix(xr[, 1], n, N); xc = matrix(xc[, 1], n, N) }
	if (dyn) { xr = xr + matrix(
		c(0, 2, -1.5),
		n,
		N,
		byrow = TRUE
	); xc = xc + matrix(c(1, 0, -1), n, N, byrow = TRUE); xr[3:4, 2] = 0; xc[5, 3] = 0 }
	Xlist = lapply(
		seq_len(N),
		function(t) { X = array(0, c(n, n, p)); X[, , 1] = 1; X[, , 2] = xr[, t]; X[, , 3] = matrix(xc[, t], n, n, byrow = TRUE); X[, , 4] = xd[, , t]; X }
	)
	M = array(runif(n * n * N) < .6, c(n, n, N)); for (t in seq_len(N)) diag(M[, , t]) = FALSE
	Qar = lame:::.lame_cov_ridge_ar_precision(N, rho_ab, sigma_ab)
	Lq = t(chol(Sab)); Lar = t(chol(solve(Qar))); Lb = t(chol(solve(P0)))
	ctr = function(m) if (center) m - rep(colMeans(m), each = nrow(m)) else m
	draw = function() {
		beta = drop(Lb %*% rnorm(p))
		lev = t(Lq %*% matrix(rnorm(2 * n), 2)); al = lev[, 1]; bl = lev[, 2]
		if (center) { al = al - mean(al); bl = bl - mean(bl) }
		if (dyn) { am = al + ctr(
			t(Lar %*% matrix(rnorm(N * n), N))
		); bm = bl + ctr(t(Lar %*% matrix(rnorm(N * n), N))) } else { am = NULL; bm = NULL }
		at = if (dyn) am else matrix(al, n, N); bt = if (dyn) bm else matrix(bl, n, N)
		EZ = Z = array(0, c(n, n, N))
		for (t in seq_len(N)) {
			E = matrix(0, n, n)
			for (i in 1:(n - 1L)) for (j in (i + 1L):n) { u = rnorm(2); E[i, j] = sqrt(s2) * u[1]; E[j, i] = sqrt(s2) * (rho * u[1] + sqrt(1 - rho^2) * u[2]) }
			diag(E) = rnorm(n)
			eta = Reduce(
				`+`,
				lapply(seq_len(p), function(k) beta[k] * Xlist[[t]][, , k])
			) + outer(at[, t], bt[, t], '+')
			EZ[, , t] = eta; Z[, , t] = eta + E
		}
		list(beta = beta, al = al, bl = bl, am = am, bm = bm, EZ = EZ, Z = Z)
	}
	step = function(s) {
		o = lame:::.lame_cov_ridge_visit(
			Xlist,
			s$beta,
			P0,
			s$al,
			s$bl,
			s$am,
			s$bm,
			if (mutant == 'sab') 4 * diag(2) else Sab,
			rho_ab,
			sigma_ab,
			TRUE,
			TRUE,
			center = center,
			E = s$Z - s$EZ,
			M = M,
			rho_t = rep(rho, N),
			s2 = s2,
			skip = if (center) 1L else integer(0)
		)
		for (t in seq_len(N)) {
			s$EZ[, , t] = s$EZ[, , t] + o$shift[t]
			zt = s$Z[, , t]
			if (mutant != 'noz') zt[!M[, , t]] = zt[!M[, , t]] + o$shift[t]
			s$Z[, , t] = zt
		}
		s$beta = o$beta; s$al = o$a_lev; s$bl = o$b_lev; s$am = o$a_mat; s$bm = o$b_mat
		s
	}
	free1 = which(!M[, , 1] & row(M[, , 1]) != col(M[, , 1]))[1]
	stat = function(s) c(
		s$beta,
		s$al[1:3],
		s$bl[1:3],
		if (dyn) c(s$am[1:2, 2], s$bm[2:3, 3]),
		s$beta[2] * s$al[1],
		s$beta[3]^2,
		(s$Z - s$EZ)[1, 2, 1],
		(s$Z - s$EZ)[2, 1, 1]^2,
		s$Z[2, 3, 1],
		s$EZ[1, 2, 1],
		sum(s$beta^2),
		sum(s$al^2),
		sum(s$bl^2),
		if (dyn) sum((s$am - s$al)^2),
		(s$Z - s$EZ)[free1]^2,
		(s$Z - s$EZ)[1, 1, N]^2
	)
	S0 = S1 = NULL
	for (r in seq_len(ndraw)) { s = draw(); S0 = rbind(S0, stat(s)); S1 = rbind(S1, stat(step(s))) }
	d = S1 - S0
	sdv = apply(d, 2L, sd)
	z = colMeans(d)[sdv > 0] / (sdv[sdv > 0] / sqrt(ndraw))
	max(abs(z))
}

test_that('one visit of the covariate ridge leaves the prior of the state, the latent array included, unchanged (dynamic effects, static effects, ame)', {
	skip_on_cran()
	expect_lt(cr_exact('dyn'), 4.2)
	expect_lt(cr_exact('static'), 4.2)
	expect_lt(cr_exact('ame', 4000L), 4.2)
})

test_that('a wrong level covariance, a missing data term, unmoved free latent cells and a wrong persistence are caught', {
	skip_on_cran()
	expect_gt(cr_exact('dyn', 6000L, mutant = 'sab'), 5)
	expect_gt(cr_exact('dyn', 6000L, mutant = 'noz'), 5)
	testthat::local_mocked_bindings(.lame_cov_ridge_data = function(E, M, rho, s2) c(q = 0, c = 0), .package = 'lame')
	expect_gt(cr_exact('dyn', 6000L), 5)
})

test_that('the covariate types and the AR precision are what the visit assumes', {
	n = 5L
	Xl = lapply(
		1:2,
		function(t) { X = array(1, c(n, n, 4)); X[, , 2] = rnorm(n); X[, , 3] = matrix(rnorm(n), n, n, byrow = TRUE); X[, , 4] = rnorm(n * n); X }
	)
	ty = lame:::.lame_cov_ridge_types(Xl, 1L)
	expect_equal(ty[[2]]$type, 'row'); expect_equal(ty[[3]]$type, 'col'); expect_null(ty[[4]]); expect_null(ty[[1]])
	Q = lame:::.lame_cov_ridge_ar_precision(4L, .5, .8)
	S = .64 / (1 - .25) * .5^abs(outer(1:4, 1:4, '-'))
	expect_equal(Q %*% (S / .64 * .64), diag(4), tolerance = 1e-8)
	expect_equal(lame:::.lame_cov_ridge_ar_precision(1L, .5, .8), matrix((1 - .25) / .64))
})

cr_panel = function(fam, n = 12L, N = 3L, vary = TRUE, seed = 4L) {
	set.seed(seed)
	nm = letters[seq_len(n)]
	mk = function() { Y = switch(
		fam,
		normal = matrix(rnorm(n * n), n),
		poisson = matrix(rpois(n * n, 1), n),
		binary = matrix(rbinom(n * n, 1, .3), n),
		ordinal = matrix(sample(1:3, n * n, TRUE), n)
	); diag(Y) = NA; dimnames(Y) = list(nm, nm); Y }
	xr = matrix(rnorm(n), n, 1, dimnames = list(nm, 'xr')); xc = matrix(rnorm(n), n, 1, dimnames = list(nm, 'xc'))
	rr = lapply(seq_len(N), function(t) if (vary && t > 1L) matrix(rnorm(n), n, 1, dimnames = list(nm, 'xr')) else xr)
	cc = lapply(seq_len(N), function(t) if (vary && t > 1L) matrix(rnorm(n), n, 1, dimnames = list(nm, 'xc')) else xc)
	nmp = paste0('p', seq_len(N))
	list(Y = setNames(lapply(seq_len(N), function(t) mk()), nmp), Xrow = setNames(rr, nmp), Xcol = setNames(cc, nmp))
}

test_that('lame() and ame() run the visit for every family, with static and dynamic effects, and the switch omits it', {
	skip_on_cran()
	for (fam in c('binary', 'normal', 'poisson', 'ordinal')) for (dab in c(FALSE, TRUE)) {
		d = cr_panel(fam, vary = dab)
		f = suppressWarnings(
			suppressMessages(
				lame(
					d$Y,
					Xrow = d$Xrow,
					Xcol = d$Xcol,
					family = fam,
					R = 0,
					dynamic_ab = dab,
					burn = 10,
					nscan = 20,
					odens = 1,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					seed = 1
				)
			)
		)
		expect_equal(unname(f$cov_ridge_move['visits']), 30)
		expect_gt(f$cov_ridge_move[['ridge_moved']], 0)
	}
	d = cr_panel('ordinal', vary = TRUE)
	f0 = suppressWarnings(
		suppressMessages(
			lame(
				d$Y,
				Xrow = d$Xrow,
				Xcol = d$Xcol,
				family = 'ordinal',
				R = 0,
				dynamic_ab = TRUE,
				burn = 10,
				nscan = 20,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(cov_ridge_move = FALSE)
			)
		)
	)
	expect_null(f0$cov_ridge_move)
	d = cr_panel('binary', N = 1L)
	fa = suppressWarnings(
		suppressMessages(
			ame(
				d$Y[[1]],
				Xrow = d$Xrow[[1]],
				Xcol = d$Xcol[[1]],
				family = 'binary',
				R = 0,
				burn = 10,
				nscan = 20,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_equal(unname(fa$cov_ridge_move['visits']), 30)
})
