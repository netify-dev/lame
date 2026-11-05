.collapsed_translation_fixture = function() {
	set.seed(492761)
	n = 5L; rank = 2L; nt = 3L; H = qr.Q(qr(contr.helmert(n)))
	ru = .72; su = .5; ra = .5; sa = .3
	iV = matrix(c(1.4, .2, .2, .8), 2); Sab = matrix(c(.7, .35, .35, 1.1), 2)
	covariance = function(rho, sigma) sigma^2 / (1 - rho^2) * rho^abs(outer(1:nt, 1:nt, '-'))
	Cu = covariance(ru, su); Ka = covariance(ra, sa)
	draw = function(bound = Inf) {
		repeat {
			U = array(matrix(rnorm(n * rank * nt), n * rank, nt) %*% chol(Cu), c(n, rank, nt))
			V = array(matrix(rnorm(n * rank * nt), n * rank, nt) %*% chol(Cu), c(n, rank, nt))
			if (max(abs(c(U, V))) < bound) break
		}
		lv = H %*% (matrix(rnorm((n - 1) * 2), n - 1, 2) %*% chol(Sab))
		list(
			U = U,
			V = V,
			margin = lv[, 2] + H %*% (matrix(rnorm((n - 1) * nt), n - 1, nt) %*% chol(Ka)),
			level = lv[, 2],
			other_level = lv[, 1],
			beta = drop(rnorm(2) %*% chol(solve(iV)))
		)
	}
	step = function(x, bound = Inf) {
		y = .lame_binary_dynamic_translation_collapsed_step(
			x$U,
			x$V,
			x$margin,
			x$level,
			x$other_level,
			x$beta,
			iV,
			ru,
			su,
			ra,
			sa,
			Sab,
			2L,
			bound
		)
		y$V = x$V; y$other_level = x$other_level; y
	}
	predictor = function(x) sapply(
		seq_len(nt),
		function(t) c(
			tcrossprod(x$U[, , t], x$V[, , t]) + matrix(x$margin[, t], n, n, byrow = TRUE) + x$beta[1]
		)
	)
	density = function(x) {
		u = matrix(x$U, ncol = nt); d = x$margin - Sab[2, 1] / Sab[1, 1] * x$other_level
		-.5 * (sum((u %*% solve(Cu)) * u) +
				sum((d %*% solve(Ka + Sab[2, 2] - Sab[2, 1]^2 / Sab[1, 1])) * d) + drop(x$beta %*% iV %*% x$beta))
	}
	statistic = function(x) c(
		x$beta[1],
		x$beta[1]^2,
		mean(x$U),
		mean(x$U^2),
		x$margin[1, 1],
		x$margin[1, 1]^2,
		mean(x$margin^2),
		mean(x$U) * x$beta[1],
		x$level[1],
		mean(x$level^2),
		mean(x$level * x$other_level),
		mean((x$margin - x$level)^2)
	)
	list(
		draw = draw,
		step = step,
		predictor = predictor,
		density = density,
		statistic = statistic,
		iV = iV,
		ru = ru,
		su = su,
		ra = ra,
		sa = sa,
		Sab = Sab
	)
}

test_that('collapsed actor translation agrees with a dense marginal Gaussian target', {
	f = .collapsed_translation_fixture(); x = f$draw()
	z = .lame_binary_dynamic_translation_collapsed_system(
		x$U,
		x$V,
		x$margin,
		x$other_level,
		x$beta,
		f$iV,
		f$ru,
		f$su,
		f$ra,
		f$sa,
		f$Sab,
		2L
	)
	q = rnorm(ncol(z$basis))
	y = .lame_binary_dynamic_translation_apply(x$U, x$V, x$margin, x$beta, z$basis, q)
	y$V = x$V; y$other_level = x$other_level
	expect_equal(f$predictor(x), f$predictor(y), tolerance = 1e-10)
	expect_equal(
		f$density(y) - f$density(x),
		-sum(q * z$score) - drop(q %*% z$precision %*% q) / 2,
		tolerance = 1e-10
	)
})

test_that('collapsed actor translation preserves the full prior including coordinate bounds', {
	skip_on_cran()
	f = .collapsed_translation_fixture()
	for (bound in c(Inf, 1.6)) {
		N = 2000L
		d = replicate(N, {
			x = f$draw(bound); y = f$step(x, bound)
			c(
				f$statistic(y) - f$statistic(x),
				y$bound_rejected,
				max(abs(colSums(y$margin))),
				abs(sum(y$level)),
				max(abs(y$U))
			)
		})
		expect_lt(max(d[14:15, ]), 1e-10)
		if (is.finite(bound)) {
			expect_gt(sum(d[13, ]), 0); expect_lte(max(d[16, ]), bound)
		}
		delta = d[1:12, , drop = FALSE]
		mcse = rowMeans(delta) / (apply(delta, 1, sd) / sqrt(N))
		expect_true(all(is.finite(mcse))); expect_lt(max(abs(mcse)), 6)
	}
})
