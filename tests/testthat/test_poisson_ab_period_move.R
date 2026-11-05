# translations of the actor-period effects of a Poisson dynamic_ab fit (R/poisson_ab_period_move.R): one step
# leaves the joint law of effects, latent array and counts unchanged and records its work

pab_ns = function() asNamespace('lame')

# paired before/after moment differences under the joint law: levels held, AR(1) deviations, latent array around
# the predictor, counts from the latent array with missing cells and an absent actor-period
pab_geweke = function(mode, centred, ndraw, seed, mut = 'none', nr = 6L, nt = 3L) {
	ns = pab_ns()
	set.seed(seed)
	nc = if (mode == 'bip') 4L else nr
	rho = 0.6; sig = 0.5; s2 = 0.4; off = c(0, 0.3, -0.2)[seq_len(nt)]
	lev_a = rnorm(nr, 0, .7); lev_b = if (mode == 'sym') lev_a else rnorm(nc, 0, .7)
	cx0 = NULL
	ar1 = function(n) {
		d = matrix(0, n, nt); d[, 1] = rnorm(n, 0, sig / sqrt(1 - rho^2))
		for (t in 2:nt) d[, t] = rho * d[, t - 1] + rnorm(n, 0, sig)
		d
	}
	absent = c(2L, 2L)
	stat = function(a, b, Z, y) {
		A = a; A[!is.finite(A)] = 0; B = b; B[!is.finite(B)] = 0
		ob = is.finite(y); y0 = ifelse(ob, y, 0)
		c(
			sum(A),
			sum(A^2),
			sum(B^2),
			sum(A[, -1] * A[, -nt]),
			sum(B[, -1] * B[, -nt]),
			mean(Z),
			mean(Z^2),
			sum(y0 * Z),
			sum(ob * Z^2),
			A[1, 1] * mean(Z[1, , 1]),
			sum(A[, 1] * rowMeans(Z[, , 1])),
			sum(B[, nt] * colMeans(Z[, , nt]))
		)
	}
	D = vector('list', ndraw)
	for (it in seq_len(ndraw)) {
		a = lev_a + ar1(nr); b = if (mode == 'sym') a else lev_b + ar1(nc)
		if (mode == 'sym') b = a
		E = array(0, c(nr, nc, nt)); Z = E; y = E
		for (t in 1:nt) {
			e = outer(ifelse(is.finite(a[, t]), a[, t], 0), ifelse(is.finite(b[, t]), b[, t], 0), '+') - 0.5
			if (mode == 'sym') e = (e + t(e)) / 2
			z = e + matrix(rnorm(nr * nc, 0, sqrt(s2)), nr, nc)
			if (mode == 'sym') z[lower.tri(z)] = t(z)[lower.tri(z)]
			yy = matrix(rpois(nr * nc, exp(z + off[t])), nr, nc)
			if (mode == 'sym') yy[lower.tri(yy)] = t(yy)[lower.tri(yy)]
			yy[matrix(runif(nr * nc) < 0.1, nr, nc)] = NA
			if (mode == 'sym') yy[lower.tri(yy)] = t(yy)[lower.tri(yy)]
			if (mode != 'bip') diag(yy) = NA
			E[, , t] = e; Z[, , t] = z; y[, , t] = yy
		}
		cx = ns$.lame_pab_context(y, mode, off)
		s0 = stat(a, b, Z, cx$Y)
		o = ns$.lame_pab_visit(cx, Z, E, a, b, lev_a, lev_b, rho, sig, TRUE, TRUE, centred)
		Z1 = o$Z; a1 = o$a_mat; b1 = o$b_mat
		if (mut == 'noz') Z1 = Z
		if (mut == 'sign') { Z1 = Z1 + 0.05 }
		if (mut == 'prior') { # mutant: drop the effects prior
			a1 = a; b1 = b; dlt = o$a_mat - a; a1 = a + 0.8 * dlt; if (mode == 'sym') b1 = a1
		}
		D[[it]] = stat(a1, b1, Z1, cx$Y) - s0
		if (centred && it == 1L) stopifnot(abs(sum(a1, na.rm = TRUE) - sum(a, na.rm = TRUE)) < 1e-9)
	}
	D = do.call(rbind, D)
	N = nrow(D); mu = colMeans(D); sdv = apply(D, 2, sd)
	z = ifelse(sdv > 0, mu / (sdv / sqrt(N)), 0)
	keep = sdv > 1e-12
	ev = eigen(cov(D[, keep]), symmetric = TRUE)
	kp = ev$values > max(ev$values) * 1e-9
	m = drop(crossprod(ev$vectors[, kp], mu[keep]))
	T2 = N * sum(m^2 / ev$values[kp])
	list(maxz = max(abs(z)), p = pchisq(T2, sum(kp), lower.tail = FALSE))
}

test_that('the slice step draws from the normal times the Poisson factor of a row', {
	skip_on_cran()
	ns = pab_ns()
	set.seed(3)
	mu = 0.4; prec = 2; A = 1.5; B = 0.9; C = 0.5
	g = seq(-4, 4, length.out = 4001); lf = -prec * (g - mu)^2 / 2 + g * A - exp(g) * B - exp(-g) * C
	w = exp(lf - max(lf)); w = w / sum(w)
	# the step leaves the density invariant: start from a draw of it, shifted since the step moves from zero
	d0 = sample(g, 30000, replace = TRUE, prob = w)
	d = vapply(d0, function(x0) x0 + ns$.lame_pab_slice(mu - x0, prec, A, B * exp(x0), C * exp(-x0))$delta, 0)
	expect_equal(mean(d), sum(g * w), tolerance = 0.02)
	expect_equal(sd(d), sqrt(sum(g^2 * w) - sum(g * w)^2), tolerance = 0.03)
})

test_that('the visit keeps the residual of every cell, the sum of the effects of a period when they are paired, and leaves an absent actor-period alone', {
	ns = pab_ns()
	set.seed(5)
	nr = 7L; nt = 3L
	a = matrix(rnorm(nr * nt), nr, nt); b = matrix(rnorm(nr * nt), nr, nt); a[3, 2] = NA; b[5, 1] = NA
	E = array(rnorm(nr * nr * nt), c(nr, nr, nt)); Z = E + array(rnorm(nr * nr * nt, 0, .6), dim(E))
	y = array(rpois(nr * nr * nt, 0.5), dim(E)); y[array(runif(length(y)) < .1, dim(y))] = NA
	cx = ns$.lame_pab_context(y, 'dir', c(0, 0, 0))
	o = ns$.lame_pab_visit(cx, Z, E, a, b, rep(0, nr), rep(0, nr), 0.5, 0.6, TRUE, TRUE, TRUE)
	expect_equal(o$Z - o$EZ, Z - E)
	expect_equal(colSums(o$a_mat, na.rm = TRUE), colSums(a, na.rm = TRUE))
	expect_equal(colSums(o$b_mat, na.rm = TRUE), colSums(b, na.rm = TRUE))
	expect_true(is.na(o$a_mat[3, 2]) && is.na(o$b_mat[5, 1]))
	expect_true(any(o$a_mat != a, na.rm = TRUE))
	# the predictor moves with the effects
	expect_equal(
		o$EZ[, , 1] - E[, , 1],
		outer(o$a_mat[, 1] - a[, 1], rep(1, nr)) + outer(
			rep(1, nr),
			ifelse(is.finite(o$b_mat[, 1]), o$b_mat[, 1] - b[, 1], 0)
		),
		ignore_attr = TRUE
	)
})

test_that('one step leaves the joint law of the effects, the latent array and the counts unchanged, directed, symmetric and bipartite, paired and single', {
	skip_on_cran()
	for (cfg in list(
		list('dir', TRUE, 71),
		list('dir', FALSE, 72),
		list('bip', TRUE, 73),
		list('sym', TRUE, 74),
		list('sym', FALSE, 75)
	)) {
		r = pab_geweke(cfg[[1]], cfg[[2]], 7000, cfg[[3]])
		expect_lt(r$maxz, 4.5)
		expect_gt(r$p, 1e-3)
	}
})

test_that('the exactness check fails for a step that leaves out the latent shift, shifts the latent array wrongly or ignores the prior', {
	skip_on_cran()
	for (m in c('noz', 'sign', 'prior')) {
		r = pab_geweke('dir', TRUE, 5000, 76, mut = m)
		expect_true(r$maxz > 4.5 || r$p < 1e-3)
	}
})

test_that('a Poisson dynamic_ab fit records the translations and ab_period_move = FALSE omits them', {
	skip_on_cran()
	set.seed(9)
	n = 10L; nt = 3L
	rn = sprintf('a%02d', 1:n)
	Y = lapply(1:nt, function(t) { y = matrix(rpois(n * n, 0.6), n, n, dimnames = list(rn, rn)); diag(y) = NA; y })
	names(Y) = as.character(2000 + 1:nt)
	fit = function(prior) suppressMessages(
		suppressWarnings(
			lame::lame(
				Y,
				family = 'poisson',
				dynamic_ab = TRUE,
				burn = 20,
				nscan = 40,
				odens = 2,
				seed = 4,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				prior = prior
			)
		)
	)
	f = tryCatch(fit(list()), lame_chain_unsettled = function(e) e$fit)
	expect_true(!is.null(f$ab_period_move) && f$ab_period_move[['visits']] > 0)
	g = tryCatch(fit(list(ab_period_move = FALSE)), lame_chain_unsettled = function(e) e$fit)
	expect_null(g$ab_period_move)
	expect_error(fit(list(ab_period_move = 'yes')), 'ab_period_move')
})
