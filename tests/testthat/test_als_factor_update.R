# penalized factor updates and feasible reweighting steps

test_that("compiled factor sweeps solve the weighted ridge conditionals", {
	set.seed(307)
	for (rank in c(1L, 3L)) {
		nr = 8L; nc = 11L
		u = matrix(rnorm(nr * rank), nr, rank)
		v = matrix(rnorm(nc * rank), nc, rank)
		z = matrix(rnorm(nr * nc), nr, nc)
		w = matrix(runif(nr * nc), nr, nc)
		w[w < 0.2] = 0
		w[2, ] = 0
		su = seq_len(rank) / 2
		sv = seq_len(rank) / 3
		u_before = u; v_before = v
		rng_before = .Random.seed
		f = lame:::als_factor_update_cpp(z, w, u, v, su, sv, 1, 2 * rank + 2, 1L, 1e-6)
		trace_sum = function(d, weights, scale) {
			Reduce(`+`, lapply(seq_len(ncol(weights)), function(k) {
				diag(chol2inv(chol(crossprod(d, weights[, k] * d) + diag(1 / scale, rank))))
			}))
		}
		su_new = (1 + colSums(u^2) + trace_sum(v, t(w), su)) / (nr + 2 * rank + 2)
		sv_new = (1 + colSums(v^2) + trace_sum(u, w, sv)) / (nc + 2 * rank + 2)
		# qr on the augmented design, independent of the normal-equation solver
		ridged = function(d, y, weight, scale) {
			as.vector(qr.solve(
				rbind(sqrt(weight) * d, diag(1 / sqrt(scale), rank)),
				c(sqrt(weight) * y, rep(0, rank))
			))
		}
		u_qr = do.call(rbind, lapply(seq_len(nr), function(i) ridged(v, z[i, ], w[i, ], su_new)))
		v_qr = do.call(rbind, lapply(seq_len(nc), function(j) ridged(u_qr, z[, j], w[, j], sv_new)))
		expect_equal(f$U, u_qr, tolerance = 1e-11)
		expect_equal(f$V, v_qr, tolerance = 1e-11)
		expect_equal(as.vector(f$sig2_u), su_new, tolerance = 1e-11)
		expect_equal(as.vector(f$sig2_v), sv_new, tolerance = 1e-11)
		expect_equal(f$Omat, tcrossprod(u_qr, v_qr), tolerance = 1e-11)
		expect_equal(f$U[2, ], rep(0, rank))
		expect_identical(u, u_before)
		expect_identical(v, v_before)
		expect_identical(.Random.seed, rng_before)
	}
})

test_that("damped reweighting steps keep the factor product they evaluate", {
	set.seed(311)
	state = function(symmetric) {
		u = matrix(rnorm(18), 9, 2)
		v = if (symmetric) u else matrix(rnorm(14), 7, 2)
		l = if (symmetric) c(2, -1) else NULL
		if (symmetric) u = qr.Q(qr(u))
		if (symmetric) v = u
		list(
			mu = 1,
			beta = 0.5,
			a = rep(0, 9),
			b = rep(0, nrow(v)),
			U = u,
			V = v,
			L = l,
			Omat = if (symmetric) u %*% (l * t(u)) else tcrossprod(u, v),
			sig2_u = c(1, 2),
			sig2_v = c(2, 1)
		)
	}
	for (symmetric in c(FALSE, TRUE)) {
		f0 = state(symmetric); f1 = state(symmetric)
		for (step in c(0.5, 0.25, 0.0625)) {
			f = lame:::.ae_irls_blend(f0, f1, step, symmetric)
			product = if (symmetric) f$U %*% (f$L * t(f$U)) else tcrossprod(f$U, f$V)
			expect_equal(f$Omat, product, tolerance = 1e-12)
			expect_lte(sum(svd(f$Omat, nu = 0, nv = 0)$d > 1e-9), 2L)
			if (symmetric) expect_equal(crossprod(f$U), diag(2), tolerance = 1e-12)
		}
	}
})
