# block draw of the persistence and the stationary variance of AR(1) coordinate paths

uvrw_paths = function(seed, nA = 20L, nB = 25L, N = 6L, rho = .96, sig = .5) {
	set.seed(seed)
	mk = function(n) { a = array(0, c(n, 1L, N)); a[, , 1] = rnorm(n, 0, sig / sqrt(1 - rho^2)); for (t in 2:N) a[, , t] = rho * a[, , t - 1] + rnorm(n, 0, sig); a }
	list(U = mk(nA), V = mk(nB))
}
# joint conditional of (rho, w) given the paths by numerical integration: the marginal of rho on a grid
# (w integrated out as a truncated inverse gamma) and E[w | rho]
uvrw_grid = function(p, cap, pm = .9, ps = .1, shape = 2, scale = 1, smax = 10, bound = .99, lu = NULL, lv = NULL) {
	N = dim(p$U)[3L]
	S1 = n1 = Stt = Sxy = Syy = ntr = 0
	for (k in 1:2) {
		cu = list(p$U, p$V)[[k]]; lam = list(lu, lv)[[k]]
		S1 = S1 + sum(cu[, , 1]^2); n1 = n1 + prod(dim(cu)[1:2])
		# the transition of actor i into period t, with its weight (one without weights)
		for (t in 2:N) for (i in seq_len(dim(cu)[1])) {
			w = if (is.null(lam)) 1 else lam[i, t]
			if (!is.finite(w) || w <= 0) next
			x = cu[i, , t]; y = cu[i, , t - 1]
			Stt = Stt + w * sum(x^2); Sxy = Sxy + w * sum(x * y); Syy = Syy + w * sum(y^2); ntr = ntr + length(x)
		}
	}
	rg = seq(-bound + 1e-3, bound - 1e-3, length.out = 1979)
	shape_w = (n1 + ntr) / 2 + shape
	one = lapply(rg, function(r) {
		om = 1 - r^2
		rate = (S1 + (Stt - 2 * r * Sxy + r^2 * Syy + 2 * scale) / om) / 2
		wmax = min(cap, smax^2 / om)
		lint = pgamma(
			1 / wmax,
			shape = shape_w,
			rate = rate,
			lower.tail = FALSE,
			log.p = TRUE
		) + lgamma(shape_w) - shape_w * log(rate)
		# E[w | rho] for w = 1 / g, g ~ Gamma(shape_w, rate) on g > 1 / wmax: E[1 / g] = rate / (shape_w - 1) * P(Gamma(shape_w - 1) > 1 / wmax) / P(Gamma(shape_w) > 1 / wmax)
		ew = rate / (shape_w - 1) * exp(
			pgamma(
				1 / wmax,
				shape = shape_w - 1,
				rate = rate,
				lower.tail = FALSE,
				log.p = TRUE
			) - pgamma(1 / wmax, shape = shape_w, rate = rate, lower.tail = FALSE, log.p = TRUE)
		)
		c(-(ntr / 2 + shape) * log(om) + dnorm(r, pm, ps, log = TRUE) + lint, ew)
	})
	lw = sapply(one, `[`, 1); ew = sapply(one, `[`, 2)
	wt = exp(lw - max(lw)); wt = wt / sum(wt)
	list(rho_mean = sum(wt * rg), rho_sd = sqrt(sum(wt * rg^2) - sum(wt * rg)^2), w_mean = sum(wt * ew))
}
uvrw_run = function(p, cap, n, seed, pm = .9, ps = .1, smax = 10, start = c(.9, .3), lu = NULL, lv = NULL) {
	set.seed(seed)
	x = matrix(NA_real_, n, 2L); s = start
	for (i in seq_len(n)) { o = .lame_uv_rho_w_step(
		p$U,
		p$V,
		s[1],
		s[2],
		pm,
		ps,
		2,
		1,
		sigma_max = smax,
		stat_var_cap = cap,
		lambda_u = lu,
		lambda_v = lv
	); s = c(o$rho, o$sigma); x[i, ] = s }
	cbind(rho = x[, 1], sigma = x[, 2], w = x[, 2]^2 / (1 - x[, 1]^2))
}
uvrw_ess = function(x) { x = x - mean(x); n = length(x); a = acf(
	x,
	lag.max = min(n - 1L, 400L),
	plot = FALSE
)$acf[-1]; k = which(a < .05)[1]; if (is.na(k)) k = length(a); n / (1 + 2 * sum(a[seq_len(k - 1L)])) }

test_that('the block draw leaves the joint conditional of the persistence and the stationary variance as it is', {
	skip_on_cran()
	for (case in list(list(cap = 4, seed = 1L), list(cap = 1, seed = 2L), list(cap = Inf, seed = 3L))) {
		p = uvrw_paths(case$seed)
		g = uvrw_grid(p, case$cap)
		x = uvrw_run(p, case$cap, 12000L, 10L + case$seed)[-(1:500), ]
		# the mean of rho and the mean of w against the density, with the Monte Carlo error of the chain
		e_rho = uvrw_ess(x[, 'rho']); e_w = uvrw_ess(x[, 'w'])
		expect_lt(abs(mean(x[, 'rho']) - g$rho_mean) / (sd(x[, 'rho']) / sqrt(e_rho)), 4)
		expect_lt(abs(sd(x[, 'rho']) / g$rho_sd - 1), .12)
		expect_lt(abs(mean(x[, 'w']) - g$w_mean) / (sd(x[, 'w']) / sqrt(e_w)), 4)
		# the draws stay in the support
		expect_true(all(abs(x[, 'rho']) < .99))
		expect_true(all(x[, 'sigma'] <= 10) && all(x[, 'w'] <= case$cap * (1 + 1e-9)))
	}
})

test_that('the block draw with the precision weights of the Student-t kind leaves the joint conditional as it is', {
	skip_on_cran()
	for (case in list(list(cap = 4, seed = 11L), list(cap = 1, seed = 12L), list(cap = Inf, seed = 13L))) {
		p = uvrw_paths(case$seed)
		set.seed(100L + case$seed)
		lu = matrix(rgamma(dim(p$U)[1L] * dim(p$U)[3L], 3, 3), dim(p$U)[1L], dim(p$U)[3L])
		lv = matrix(rgamma(dim(p$V)[1L] * dim(p$V)[3L], 3, 3), dim(p$V)[1L], dim(p$V)[3L])
		# a transition without a weight is left out
		lu[2, 3] = 0; lv[4, 5] = NA
		g = uvrw_grid(p, case$cap, lu = lu, lv = lv)
		g0 = uvrw_grid(p, case$cap)
		expect_gt(abs(g$rho_mean - g0$rho_mean) + abs(g$w_mean - g0$w_mean), 1e-4)
		x = uvrw_run(p, case$cap, 12000L, 20L + case$seed, lu = lu, lv = lv)[-(1:500), ]
		e_rho = uvrw_ess(x[, 'rho']); e_w = uvrw_ess(x[, 'w'])
		expect_lt(abs(mean(x[, 'rho']) - g$rho_mean) / (sd(x[, 'rho']) / sqrt(e_rho)), 4)
		expect_lt(abs(sd(x[, 'rho']) / g$rho_sd - 1), .12)
		expect_lt(abs(mean(x[, 'w']) - g$w_mean) / (sd(x[, 'w']) / sqrt(e_w)), 4)
		expect_true(all(abs(x[, 'rho']) < .99))
	}
	# weights of one give the unweighted draw, draw for draw
	p = uvrw_paths(4L)
	ones = matrix(1, dim(p$U)[1L], dim(p$U)[3L]); onesv = matrix(1, dim(p$V)[1L], dim(p$V)[3L])
	expect_equal(uvrw_run(p, 4, 50L, 3L), uvrw_run(p, 4, 50L, 3L, lu = ones, lv = onesv), tolerance = 1e-6)
})

test_that('the block draw moves along the cap', {
	skip_on_cran()
	p = uvrw_paths(2L)
	g = uvrw_grid(p, 1)
	new = uvrw_run(p, 1, 6000L, 5L)[-(1:200), ]
	# the persistence explores a conditional of sd about .009
	expect_gt(uvrw_ess(new[, 'rho']), 500)
	expect_lt(abs(mean(new[, 'rho']) - g$rho_mean), 4 * g$rho_sd / sqrt(uvrw_ess(new[, 'rho'])))
})

test_that('the block draw leaves a state it cannot update, a path of one period and a symmetric network as it should', {
	p = uvrw_paths(4L)
	one = list(U = p$U[, , 1, drop = FALSE], V = p$V[, , 1, drop = FALSE])
	o = .lame_uv_rho_w_step(one$U, one$V, .9, .3, .9, .1, 2, 1)
	expect_equal(unlist(o), c(rho = .9, sigma = .3))
	bad = p$U; bad[1, 1, 2] = NA
	o = .lame_uv_rho_w_step(bad, p$V, .9, .3, .9, .1, 2, 1)
	expect_equal(unlist(o), c(rho = .9, sigma = .3))
	# a symmetric network has one set of coordinates: V is ignored
	set.seed(6)
	s1 = .lame_uv_rho_w_step(p$U, NULL, .9, .3, .9, .1, 2, 1, symmetric = TRUE)
	set.seed(6)
	s2 = .lame_uv_rho_w_step(p$U, p$V * 1000, .9, .3, .9, .1, 2, 1, symmetric = TRUE)
	expect_identical(s1, s2)
	# a start outside the support (rho beyond the bound) is pulled in and the sigma is drawn within the cap
	o = .lame_uv_rho_w_step(p$U, p$V, .999, .3, .9, .1, 2, 1, stat_var_cap = 1)
	expect_lte(abs(o$rho), .99); expect_lte(o$sigma^2 / (1 - o$rho^2), 1 + 1e-9)
	# and a persistence at the bound moves
	set.seed(7)
	moved = vapply(1:20, function(i) .lame_uv_rho_w_step(p$U, p$V, .99, .3, .9, .1, 2, 1)$rho, numeric(1))
	expect_gt(length(unique(moved)), 15L)
})

test_that('a fit draws the persistence and the stationary variance as a block', {
	skip_on_cran()
	set.seed(9)
	nA = 14L; nB = 12L; N = 4L
	Y = lapply(
		seq_len(N),
		function(t) { y = matrix(rpois(nA * nB, .6), nA, nB); y[1:4, ] = 0; dimnames(y) = list(
			sprintf('r%02d', 1:nA),
			sprintf('c%02d', 1:nB)
		); y }
	)
	fit = function(prior) suppressWarnings(suppressMessages(lame(
		Y,
		family = 'poisson',
		mode = 'bipartite',
		R_row = 1L,
		R_col = 1L,
		dynamic_uv = TRUE,
		burn = 30,
		nscan = 100,
		odens = 2,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = prior
	)))
	on = fit(list())
	expect_true(all(is.finite(on$rho_uv)) && all(is.finite(on$sigma_uv)))
	expect_true(all(on$sigma_uv^2 / (1 - on$rho_uv^2) <= 4 * (1 + 1e-9)))
})

test_that('a directed fit of the Student-t kind draws the persistence and the stationary variance as a block too', {
	skip_on_cran()
	set.seed(10)
	n = 16L; N = 4L
	Y = lapply(
		seq_len(N),
		function(t) { y = matrix(rpois(n * n, .5), n, n); diag(y) = NA; y[1:3, ] = 0; y[, 1:3] = 0; dimnames(y) = list(
			sprintf('a%02d', 1:n),
			sprintf('a%02d', 1:n)
		); y }
	)
	fit = function(prior) suppressWarnings(suppressMessages(lame(
		Y,
		family = 'poisson',
		R = 1L,
		dynamic_uv = TRUE,
		dynamic_uv_kind = 't',
		burn = 30,
		nscan = 100,
		odens = 2,
		seed = 3,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		prior = c(list(als_start = FALSE), prior)
	)))
	on = fit(list())
	expect_true(all(is.finite(on$rho_uv)) && all(is.finite(on$sigma_uv)))
	expect_true(all(on$sigma_uv^2 / (1 - on$rho_uv^2) <= 4 * (1 + 1e-9)))
})
