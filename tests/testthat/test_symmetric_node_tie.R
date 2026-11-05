# the row and column coefficients of a nodal covariate in a symmetric binary or ordinal fit are tied after every draw

stn_sym = function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m }

# small symmetric panel with a time-varying nodal covariate, a dyadic covariate and dynamic additive effects
stn_panel = function(fam = 'binary', n = 30L, Tn = 4L, seed = 1L) {
	set.seed(500 + seed)
	rn = sprintf('a%02d', seq_len(n))
	lev = rnorm(n, 0, .7)
	dev = matrix(0, n, Tn)
	dev[, 1] = rnorm(n, 0, .4)
	for (t in 2:Tn) dev[, t] = .7 * dev[, t - 1] + rnorm(n, 0, .4 * sqrt(1 - .49))
	xs = rnorm(n)
	zz = lapply(seq_len(Tn), function(t) {
		xr = sqrt(.5) * xs + sqrt(.5) * rnorm(n)
		xd = stn_sym(matrix(rnorm(n * n), n, n))
		at = lev + dev[, t]
		list(
			xr = xr,
			xd = xd,
			z = .6 * xd + .35 * outer(xr, xr, '+') + outer(at, at, '+') + stn_sym(matrix(rnorm(n * n), n, n))
		)
	})
	b0 = -unname(quantile(unlist(lapply(zz, function(o) o$z[upper.tri(o$z)])), if (fam == 'poisson') .9 else .92))
	Y = X = Xr = vector('list', Tn)
	for (t in seq_len(Tn)) {
		zt = zz[[t]]$z + b0
		y = if (fam == 'binary') 1 * (zt > 0) else matrix(rpois(n * n, exp(zt - 1)), n, n)
		y = stn_sym(y)
		diag(y) = NA
		dimnames(y) = list(rn, rn)
		Y[[t]] = y
		X[[t]] = array(zz[[t]]$xd, c(n, n, 1), dimnames = list(rn, rn, 'xd'))
		Xr[[t]] = matrix(zz[[t]]$xr, n, 1, dimnames = list(rn, 'xn'))
	}
	names(Y) = names(X) = names(Xr) = paste0('p', seq_len(Tn))
	list(Y = Y, Xdyad = X, Xrow = Xr)
}
stn_fit = function(d, fam, prior = list(), burn = 10L, nscan = 40L, R = 0L, ...) {
	suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$Xdyad,
		Xrow = d$Xrow,
		symmetric = TRUE,
		family = fam,
		R = R,
		burn = burn,
		nscan = nscan,
		odens = 1L,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior,
		...
	)))
}

test_that('the tie sets the two coefficients of a nodal covariate to their mean and leaves the others', {
	b = c(-1, .2, .8, .5, 1.5)
	expect_equal(.lame_sym_tie(b, 2L, 3L), c(-1, .5, .5, .5, 1.5))
	expect_equal(.lame_sym_tie(b, integer(0), integer(0)), b)
	expect_equal(.lame_sym_tie(c(.2, .8, 1), 1L, 2L), c(.5, .5, 1))
	# two nodal covariates: row slices first, then column slices, as in the design
	expect_equal(.lame_sym_tie(c(0, .1, .3, .9, .5, 2), 2:3, 4:5), c(0, .5, .4, .5, .4, 2))
	# both slices dynamic: paths are tied period by period
	path = cbind(c(1, 2, 3), c(3, 2, 1))
	st = .lame_sym_tie_state(path, c(.5, 1.5), c(2L, 3L), c(1L, 4L), 2L, 3L)
	expect_equal(st$path, cbind(c(2, 2, 2), c(2, 2, 2)))
	expect_equal(st$stat, c(.5, 1.5))
	# both slices static: static entries are tied and a dynamic path is left alone
	st = .lame_sym_tie_state(matrix(1:3, 3, 1), c(.2, .8, 1), 1L, 2:4, 2L, 3L)
	expect_equal(st$stat, c(.5, .5, 1))
	expect_equal(st$path, matrix(1:3, 3, 1))
})


test_that('the exact symmetric draw cannot separate the two slices of a nodal covariate, and the tied draw has the law of their mean', {
	skip_on_cran()
	set.seed(3)
	n = 8L; Tn = 2L; p = 4L
	xr = rnorm(n)
	X = lapply(seq_len(Tn), function(t) {
		a = array(0, c(n, n, p))
		a[, , 1] = 1
		a[, , 2] = matrix(xr, n, n)
		a[, , 3] = matrix(xr, n, n, byrow = TRUE)
		a[, , 4] = stn_sym(matrix(rnorm(n * n), n, n))
		a
	})
	Q = crossprod(matrix(rnorm(16), 4)) / 20 + diag(.05, 4)
	cx = .lame_symmetric_gaussian_context(X, Q, actor = FALSE, center = FALSE)
	# the symmetrised design has the same column for the row and column slices
	expect_equal(cx$D[, 2], cx$D[, 3])
	Z = array(unlist(lapply(seq_len(Tn), function(t) stn_sym(matrix(rnorm(n * n), n, n)))), c(n, n, Tn))
	# posterior of the coefficients given Z, from its definition
	y = unlist(lapply(seq_len(Tn), function(t) (Z[, , t][cx$index] + Z[, , t][cx$mirror]) / 2))
	Lam = crossprod(cx$D) + Q
	Sig = solve(Lam)
	m = drop(Sig %*% crossprod(cx$D, y))
	Tm = diag(p); Tm[2:3, 2:3] = .5
	mu = drop(Tm %*% m)
	V = Tm %*% Sig %*% t(Tm)
	draw = function(tie) replicate(20000, { b = .lame_symmetric_gaussian_draw(Z, 1, cx)$beta; tie(b) })
	one = draw(function(b) .lame_sym_tie(b, 2L, 3L))
	expect_equal(one[2, ], one[3, ])
	emp = rowMeans(one)
	se = sqrt(diag(V) / 20000 * 3)
	expect_lt(max(abs(emp - mu) / se), 4.5)
	expect_lt(max(abs(cov(t(one)) - V) / sqrt(outer(diag(V), diag(V)))), .03)
	# a tie that keeps one of the two rather than their mean has a different law
	bad = draw(function(b) { b[3] = b[2]; b })
	expect_gt(max(abs(rowMeans(bad) - mu) / se), 6)
})


test_that('the tie applies to the binary and ordinal families with a nodal covariate, and to nothing else', {
	skip_on_cran()
	d = stn_panel('binary', n = 12L, Tn = 3L, seed = 2L)
	set.seed(2)
	o = d
	o$Y = lapply(
		d$Y,
		function(y) { v = y[upper.tri(y)] + rnorm(sum(upper.tri(y)), 0, .3); cu = quantile(v, c(.8, .9, .96)); m = matrix(NA_real_, nrow(y), ncol(y)); m[upper.tri(m)] = findInterval(v, cu); m[lower.tri(m)] = t(m)[lower.tri(m)]; dimnames(m) = dimnames(y); m }
	)
	pd = stn_panel('poisson', n = 12L, Tn = 3L, seed = 2L)
	calls = 0L
	local_mocked_bindings(
		.lame_sym_tie = function(beta, rb, cb) { calls <<- calls + 1L; beta },
		.lame_sym_tie_state = function(
			path,
			stat,
			dyn_idx,
			stat_idx,
			rb,
			cb
		) { calls <<- calls + 1L; list(path = path, stat = stat) },
		.package = 'lame'
	)
	count = function(...) { calls <<- 0L; f = stn_fit(...); expect_true(all(is.finite(f$BETA))); calls }
	# binary and ordinal fits with dynamic additive effects, positions or coefficients
	expect_gt(count(d, 'binary', dynamic_ab = TRUE), 0)
	expect_gt(count(o, 'ordinal', dynamic_ab = TRUE), 0)
	expect_gt(count(d, 'binary', dynamic_uv = TRUE, R = 1L), 0)
	expect_gt(count(d, 'binary', dynamic_beta = c('intercept', 'dyad')), 0)
	expect_gt(count(d, 'binary', dynamic_beta = TRUE), 0)
	# the switch, the normal and Poisson families, and a fit without a nodal covariate
	expect_equal(count(d, 'binary', dynamic_ab = TRUE, prior = list(symmetric_node_tie = FALSE)), 0)
	expect_equal(count(pd, 'poisson', dynamic_ab = TRUE), 0)
	dn = d
	dn$Y = lapply(
		d$Y,
		function(y) { m = y + matrix(
			rnorm(length(y), 0, .5),
			nrow(y)
		); m = (m + t(m)) / 2; diag(m) = NA; dimnames(m) = dimnames(y); m }
	)
	expect_equal(count(dn, 'normal', dynamic_ab = TRUE), 0)
	d$Xrow = NULL
	expect_equal(count(d, 'binary', dynamic_ab = TRUE), 0)
})

test_that('the switch is validated in lame()', {
	d = stn_panel('binary', n = 10L, Tn = 2L, seed = 3L)
	expect_error(stn_fit(d, 'binary', prior = list(symmetric_node_tie = 'no')), 'symmetric_node_tie')
	f = stn_fit(d, 'binary', prior = list(symmetric_node_tie = FALSE))
	expect_true(all(is.finite(f$BETA)))
})

test_that('dynamic additive effects of a symmetric binary or ordinal fit with a nodal covariate stay at their scale', {
	skip_on_cran()
	d = stn_panel('binary', n = 30L, Tn = 4L, seed = 1L)
	# innovation sd of the paths in the data is 0.28
	on = stn_fit(d, 'binary', burn = 300L, nscan = 1700L, dynamic_ab = TRUE)
	expect_lt(mean(tail(on$sigma_ab, 100)), 1)
	expect_lt(max(abs(on$BETA[, 'xn_node'])), 2)
	# without the tie the coefficient walks along its ridge with the effect paths and the scale runs away;
	# the chain policy refuses that arm at the first flag, so the fit is taken from the condition
	off = tryCatch(
		stn_fit(
			d,
			'binary',
			prior = list(symmetric_node_tie = FALSE, chain_extensions = 0),
			burn = 300L,
			nscan = 1700L,
			dynamic_ab = TRUE
		),
		lame_chain_unsettled = function(e) e$fit
	)
	expect_gt(mean(tail(off$sigma_ab, 100)), 3)
	# the same for the ordinal family
	set.seed(2)
	o = d
	o$Y = lapply(
		d$Y,
		function(y) { v = y[upper.tri(y)] + rnorm(sum(upper.tri(y)), 0, .3); cu = quantile(v, c(.8, .9, .96)); m = matrix(NA_real_, nrow(y), ncol(y)); m[upper.tri(m)] = findInterval(v, cu); m[lower.tri(m)] = t(m)[lower.tri(m)]; dimnames(m) = dimnames(y); m }
	)
	on = stn_fit(o, 'ordinal', burn = 300L, nscan = 1200L, dynamic_ab = TRUE)
	expect_lt(mean(tail(on$sigma_ab, 100)), 1.5)
})
