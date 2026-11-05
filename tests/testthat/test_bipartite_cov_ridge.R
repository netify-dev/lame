# the ridge of a nodal covariate in a bipartite ame() fit (R/bipartite_cov_ridge.R): a row covariate trades with the
# row additive effects at no cost to any predictor, and the move draws that trade exactly

ba_ns = function() asNamespace('lame')

test_that('the ridge of a nodal covariate leaves the joint law of the coefficient and the effects unchanged, and moves the coefficient by more than the conditional draws do', {
	ns = ba_ns()
	set.seed(61)
	nr = 7L; nc = 5L; gs = 1.5; va = 0.6; vb = 0.8
	x_r = rnorm(nr); x_c = rnorm(nc)
	X = array(0, c(nr, nc, 3)); X[, , 1] = 1; X[, , 2] = matrix(x_r, nr, nc); X[, , 3] = matrix(x_c, nr, nc, byrow = TRUE)
	nd = if (identical(Sys.getenv('NOT_CRAN'), 'true')) 9000L else 600L; D = matrix(0, nd, 11)
	stat = function(be, a, b) c(be, be^2, sum(a^2), sum(b^2), sum(a * x_r), sum(b * x_c), be[2] * be[3])
	for (d in 1:nd) {
		be = rnorm(3, 0, sqrt(gs)); a = rnorm(nr, 0, sqrt(va)); b = rnorm(nc, 0, sqrt(vb))
		o = ns$.lame_bip_cov_ridge(X, be, a, b, gs, va, vb, TRUE, TRUE)
		# the predictor of every cell is unchanged
		if (d <= 40) expect_equal(
			outer(o$a, o$b, '+') + o$beta[2] * X[, , 2] + o$beta[3] * X[, , 3],
			outer(a, b, '+') + be[2] * X[, , 2] + be[3] * X[, , 3]
		)
		D[d, ] = stat(o$beta, o$a, o$b) - stat(be, a, b)
	}
	mu = colMeans(D); z = mu / (apply(D, 2, sd) / sqrt(nd))
	expect_lt(max(abs(z[apply(D, 2, sd) > 1e-12])), 4.5)
	# the intercept slice (constant over the whole matrix) is left alone
	expect_identical(
		ns$.lame_bip_cov_ridge(X[, , 1, drop = FALSE], 0.3, rnorm(nr), rnorm(nc), gs, va, vb, TRUE, TRUE)$beta,
		0.3
	)
})

ba_panel = function(family, nt = 1L, seed = 3L, nA = 14L, nB = 10L) {
	set.seed(seed)
	u = matrix(rnorm(nA * 2), nA, 2); v = matrix(rnorm(nB * 2), nB, 2)
	nms_r = sprintf('r%02d', seq_len(nA)); nms_c = sprintf('c%02d', seq_len(nB))
	Y = X = vector('list', nt)
	for (t in seq_len(nt)) {
		x1 = matrix(rnorm(nA * nB), nA, nB)
		z = -0.5 + 0.6 * x1 + outer(
			rnorm(nA, 0, .5),
			rnorm(nB, 0, .5),
			'+'
		) + 0.9 * u %*% diag(c(1, .6)) %*% t(v) + matrix(rnorm(nA * nB), nA, nB)
		y = switch(
			family,
			normal = z,
			binary = 1 * (z > 0),
			ordinal = matrix(findInterval(z, c(-0.3, 0.8, 1.6)), nA, nB),
			cbin = { y = 0 * z; for (i in 1:nA) y[i, order(z[i, ], decreasing = TRUE)[1:3]] = 1; y },
			frn = { y = 0 * z; for (i in 1:nA) y[i, order(z[i, ], decreasing = TRUE)[1:3]] = 3:1; y }
		)
		dimnames(y) = list(nms_r, nms_c); Y[[t]] = y
		X[[t]] = array(x1, c(nA, nB, 1), dimnames = list(nms_r, nms_c, 'x1'))
	}
	names(Y) = names(X) = as.character(2000 + seq_len(nt))
	list(Y = Y, X = X)
}

ba_fit = function(family, lame_fit, prior = list(), nscan = 40, burn = 20, seed = 5) {
	d = ba_panel(family, nt = if (lame_fit) 2L else 1L)
	extra = if (family %in% c('cbin', 'frn')) list(odmax = 3L) else list()
	if (lame_fit) suppressMessages(suppressWarnings(do.call(
		lame::lame,
		c(
			list(
				d$Y,
				Xdyad = d$X,
				family = family,
				mode = 'bipartite',
				R = 2,
				burn = burn,
				nscan = nscan,
				odens = 2,
				seed = seed,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				prior = prior
			),
			extra
		)
	)))
	else suppressMessages(suppressWarnings(do.call(
		lame::ame,
		c(
			list(
				d$Y[[1]],
				Xdyad = d$X[[1]],
				family = family,
				mode = 'bipartite',
				R = 2,
				burn = burn,
				nscan = nscan,
				odens = 2,
				seed = seed,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				prior = prior
			),
			extra
		)
	)))
}


test_that('with cov_ridge_move = FALSE the draws differ from the defaults, and the switch is checked', {
	skip_on_cran()
	d = ba_panel('normal', seed = 5L)
	xr = matrix(rnorm(14), ncol = 1, dimnames = list(rownames(d$Y[[1]]), 'xr'))
	fit = function(prior) suppressMessages(suppressWarnings(lame::ame(
		d$Y[[1]],
		Xdyad = d$X[[1]],
		Xrow = xr,
		family = 'normal',
		mode = 'bipartite',
		R = 2,
		burn = 20,
		nscan = 40,
		odens = 2,
		seed = 5,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior
	)))
	off = list(cov_ridge_move = FALSE)
	expect_identical(fit(off)$BETA, fit(off)$BETA)
	expect_false(identical(fit(off)$BETA, fit(list())$BETA))
	expect_error(fit(list(cov_ridge_move = 'x')), 'cov_ridge_move')
})

test_that('on a panel with a row covariate the ridge move gives the coefficient and the row effects far more effective draws', {
	skip_on_cran()
	set.seed(21)
	nA = 30L; nB = 20L
	rn = sprintf('r%02d', 1:nA); cn = sprintf('c%02d', 1:nB)
	xr = rnorm(nA); a = 0.8 * xr + rnorm(nA, 0, 0.3)
	y = matrix(
		-0.3 + outer(a, rep(1, nB)) + outer(rep(0, nA), rnorm(nB, 0, .4), '+') + matrix(rnorm(nA * nB), nA, nB),
		nA,
		nB,
		dimnames = list(rn, cn)
	)
	xr = matrix(xr, ncol = 1, dimnames = list(rn, 'xr'))
	ess = function(prior) {
		f = suppressMessages(
			suppressWarnings(
				lame::ame(
					y,
					Xrow = xr,
					family = 'normal',
					mode = 'bipartite',
					R = 0,
					burn = 200,
					nscan = 4000,
					odens = 1,
					seed = 3,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE,
					prior = prior
				)
			)
		)
		coda::effectiveSize(f$BETA[, grep('xr', colnames(f$BETA))[1]])
	}
	expect_gt(ess(list()), 2 * ess(list(cov_ridge_move = FALSE)))
})
