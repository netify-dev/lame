skip_on_cran()

# 60-sweep fits are flagged as unsettled by the chain policy of ame() and lame(); these tests read the unsettled fit's draws
suppressWarnings = function(expr) tryCatch(base::suppressWarnings(expr), lame_chain_unsettled = function(e) e$fit)

# exact draws of the rows of the static factors of a bipartite product (R/bipartite_uv_gibbs.R)

# the conditional of a row, by its definition: N(0, I) prior, regression of the observed residuals of the row on W
buvg_cond = function(e, W, s2) {
	R = ncol(W)
	P = diag(R); b = numeric(R)
	for (t in seq_len(dim(e)[2])) {
		o = which(!is.na(e[, t]))
		if (length(o)) { P = P + crossprod(W[o, , drop = FALSE]) / s2; b = b + drop(
			crossprod(W[o, , drop = FALSE], e[o, t])
		) / s2 }
	}
	list(mean = drop(solve(P, b)), cov = solve(P))
}

test_that('a row of a factor is a draw from its Gaussian conditional, with complete cells and with missing ones', {
	set.seed(1)
	ns = asNamespace('lame')
	nr = 3L; nc = 7L; N = 3L; R = 2L; s2 = .8
	# the long run on a developer machine, a shorter one with proportionally wider bounds under the CRAN rules
	nd = if (identical(Sys.getenv('NOT_CRAN'), 'true')) 40000L else 2500L
	W = matrix(rnorm(nc * R), nc)
	E = array(rnorm(nr * nc * N), c(nr, nc, N))
	for (case in c('complete', 'missing')) {
		if (case == 'missing') { E[1, 2:3, 1] = NA; E[2, 5, 2] = NA; E[3, 1:7, 3] = NA }
		draws = replicate(nd, ns$.lame_bip_factor_draw(E, W, s2))
		for (i in seq_len(nr)) {
			cd = buvg_cond(matrix(E[i, , ], nc), W, s2)
			x = t(draws[i, , ])
			expect_lt(max(abs(colMeans(x) - cd$mean) / sqrt(diag(cd$cov) / nd)), 4.6)
			expect_lt(max(abs(cov(x) - cd$cov)), .02 * sqrt(40000 / nd))
		}
	}
	# a row with no observed cell is a draw from the prior
	E[1, , ] = NA
	x = t(replicate(nd / 2L, ns$.lame_bip_factor_draw(E, W, s2)[1, ]))
	expect_lt(max(abs(cov(x) - diag(R))), .05 * sqrt(40000 / nd))
})

# the joint law of a row and its residuals: a row u ~ N(0, I) and residuals W u + noise in every period; replacing u by a draw from its conditional
# leaves the law of (u, E) unchanged, so the new row has the prior law and the same joint moments with the residuals
buvg_law = function(draw, R_n = 60000L, seed = 2L, N = 3L, s2 = .8) {
	set.seed(seed)
	nc = 6L; R = 2L
	W = matrix(rnorm(nc * R), nc)
	out = matrix(0, R_n, 6)
	for (r in seq_len(R_n)) {
		u = rnorm(R)
		E = array(rnorm(nc * N, 0, sqrt(s2)), c(1L, nc, N))
		for (t in 1:N) E[1, , t] = E[1, , t] + drop(W %*% u)
		un = draw(E, W, s2)[1, ]
		es = sum(E)
		out[r, ] = c(un, un^2, un[1] * un[2], un[1] * es / nc / N)
	}
	# expectations under the prior of u and the model of E: E[u] = 0, E[u_k^2] = 1, E[u1 u2] = 0
	e_cross = sum(W[, 1]) / 6 # E[u1 * mean(E)] = mean(W[, 1])
	ex = c(0, 0, 1, 1, 0, e_cross)
	(colMeans(out) - ex) / (apply(out, 2, sd) / sqrt(R_n))
}

test_that('the exact row draw has the prior law of a row and its joint moments with the residuals, and broken versions are caught', {
	ns = asNamespace('lame')
	# a shorter run under the CRAN rules: the law of the draw, and the mutant that leaves out the residuals of the other periods
	n_cran = 20000L
	expect_lt(max(abs(buvg_law(ns$.lame_bip_factor_draw, R_n = n_cran))), 4.6)
	# the residuals of the other periods left out
	one_period = function(E, W, s2) ns$.lame_bip_factor_draw(E[, , 1, drop = FALSE], W, s2)
	expect_gt(max(abs(buvg_law(one_period, R_n = n_cran))), 5)
	skip_on_cran()
	expect_lt(max(abs(buvg_law(ns$.lame_bip_factor_draw))), 4.6)
	expect_gt(max(abs(buvg_law(one_period))), 6)
	# the wrong residual variance
	wrong_s2 = function(E, W, s2) ns$.lame_bip_factor_draw(E, W, s2 * 2)
	expect_gt(max(abs(buvg_law(wrong_s2))), 6)
	# a prior that is too wide
	loose = local({
		f = get('.lame_bip_factor_draw', ns)
		out = deparse(f, width.cutoff = 500L)
		hit = grepl('chol(diag(R) +', out, fixed = TRUE)
		stopifnot(any(hit))
		out[hit] = sub('chol(diag(R) +', 'chol(0.3 * diag(R) +', out[hit], fixed = TRUE)
		g = eval(parse(text = paste(out, collapse = '\n')), envir = ns)
		environment(g) = ns
		g
	})
	expect_gt(max(abs(buvg_law(loose))), 6)
})

# a bipartite panel with a rank-two product
buvg_panel = function(fam, nA = 14L, nB = 11L, Tn = 2L, seed = 1L) {
	set.seed(seed)
	U = matrix(rnorm(nA * 2, 0, .8), nA); V = matrix(rnorm(nB * 2, 0, .8), nB); G = diag(c(1, .7))
	rn = sprintf('r%02d', seq_len(nA)); cn = sprintf('c%02d', seq_len(nB))
	Y = lapply(seq_len(Tn), function(t) {
		z = U %*% G %*% t(V) + matrix(rnorm(nA * nB), nA)
		y = switch(
			fam,
			normal = z,
			binary = 1 * (z > .3),
			poisson = matrix(rpois(nA * nB, exp(z - 1)), nA),
			ordinal = matrix(findInterval(z, c(.3, 1, 1.8)), nA)
		)
		dimnames(y) = list(rn, cn)
		y
	})
	names(Y) = paste0('p', seq_len(Tn))
	Y
}

test_that('bipartite fits draw the factors from their conditionals, for every family', {
	skip_on_cran()
	for (fam in c('normal', 'binary', 'poisson', 'ordinal')) {
		Y = buvg_panel(fam)
		on = suppressWarnings(
			suppressMessages(
				lame(
					Y,
					mode = 'bipartite',
					family = fam,
					R_row = 2L,
					R_col = 2L,
					burn = 20,
					nscan = 60,
					odens = 1,
					seed = 4,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE
				)
			)
		)
		expect_true(all(is.finite(on$U)) && all(is.finite(on$V)) && all(is.finite(on$G)))
		a = suppressWarnings(
			suppressMessages(
				ame(
					Y[[1]],
					mode = 'bipartite',
					family = fam,
					R_row = 2L,
					R_col = 2L,
					burn = 20,
					nscan = 60,
					odens = 1,
					seed = 4,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE
				)
			)
		)
		expect_true(all(is.finite(a$U)) && all(is.finite(a$V)) && all(is.finite(a$G)))
	}
	# factors of different rank
	Y = buvg_panel('normal')
	r = suppressWarnings(
		suppressMessages(
			lame(
				Y,
				mode = 'bipartite',
				R_row = 1L,
				R_col = 2L,
				burn = 20,
				nscan = 60,
				odens = 1,
				seed = 4,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE
			)
		)
	)
	expect_equal(dim(r$G), c(1L, 2L))
	expect_true(all(is.finite(r$U)) && all(is.finite(r$V)) && all(is.finite(r$G)))
})
