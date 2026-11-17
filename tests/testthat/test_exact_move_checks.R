# exact checks of the sampler moves and the chain statistics, small enough to run on CRAN

test_that('the snap mirror of the rotation is a reflection and its own inverse on the orbit', {
	mirror = lame:::.lame_uvrd_mirror
	set.seed(31)
	err = replicate(500, {
		B = matrix(rnorm(4), 2)
		E = mirror(B)
		# orthogonal with determinant -1, and the state the jump reaches, E'B, jumps back by E' = E^-1
		c(max(abs(crossprod(E) - diag(2))), abs(det(E) + 1), max(abs(mirror(crossprod(E, B)) - t(E))))
	})
	expect_lt(max(err[1:2, ]), 1e-12)
	expect_lt(max(err[3, ]), 1e-10)
	# a mutated mirror without the rotation is not self-inverse
	broken = function(B) {
		t2 = atan2(B[2, 1] + B[1, 2], B[1, 1] - B[2, 2])
		matrix(c(cos(t2), sin(t2), sin(t2), -cos(t2)), 2, 2)
	}
	B = matrix(c(1, 0.3, -0.8, 0.4), 2)
	E = broken(B)
	expect_gt(max(abs(broken(crossprod(E, B)) - t(E))), 0.1)
})

test_that('the rotation of dynamic paths keeps every product to rounding and leaves a degenerate state alone', {
	f = lame:::.lame_uv_rotation_dynamic_draw
	set.seed(3)
	U = array(rnorm(60), c(5, 3, 4))
	V = array(rnorm(60), c(5, 3, 4))
	o = f(U, V, 0.8, 0.6)
	expect_gt(o$moved, 0L)
	for (t in 1:4) {
		expect_lt(max(abs(tcrossprod(o$U[, , t], o$V[, , t]) - tcrossprod(U[, , t], V[, , t]))), 1e-12)
	}
	same = function(o, U, V) {
		expect_identical(o$U, U)
		expect_identical(o$V, V)
		expect_identical(o$moved, 0L)
	}
	U1 = U[, 1, , drop = FALSE]
	V1 = V[, 1, , drop = FALSE]
	same(f(U1, V1, 0.8, 0.6), U1, V1)
	same(f(U, V[, 1:2, ], 0.8, 0.6), U, V[, 1:2, ])
	same(f(U, V, 1, 0.6), U, V)
	same(f(U, V, 0.8, 0), U, V)
	same(f(U, V, 0.8, 0.6, w_u = matrix(-1, 5, 4)), U, V)
	U_one = U[, , 1, drop = FALSE]
	V_one = V[, , 1, drop = FALSE]
	same(f(U_one, V_one, 0.8, 0.6), U_one, V_one)
	Un = U
	Un[1] = NA
	same(f(Un, V, 0.8, 0.6), Un, V)
	# a state outside the bound stays
	same(f(U, V, 0.8, 0.6, max_abs = 1e-3), U, V)
})

test_that('the split of dynamic paths leaves a degenerate state alone', {
	f = lame:::.lame_uv_split_dynamic_draw
	set.seed(1)
	U = array(rnorm(12), c(3, 2, 2))
	V = array(rnorm(8), c(4, 1, 2))
	same = function(o) {
		expect_identical(o$U, U)
		expect_identical(o$V, V)
		expect_identical(o$moved, 0L)
	}
	Un = U
	Un[1] = NA
	o = f(Un, V, 0.8, 0.6)
	expect_identical(o$U, Un)
	expect_identical(o$moved, 0L)
	same(f(U, V, 1, 0.6))
	same(f(U, V, 0.8, 0))
	same(f(U, V, 0.8, NA))
	same(f(U, V, 0.8, 0.6, w_u = matrix(-1, 3, 2)))
	# a zero side carries no split
	U0 = U * 0
	o = f(U0, V, 0.8, 0.6)
	expect_identical(o$U, U0)
	expect_identical(o$V, V)
	# a state outside the bound stays
	same(f(U, V, 0.8, 0.6, max_abs = 1e-3))
	expect_identical(f(U, V[, , 1, drop = FALSE], 0.8, 0.6)$moved, 0L)
})

test_that('the R-hat parts of fixed draws have the split layout and are rank invariant', {
	parts = lame:::.lame_rhat_parts
	set.seed(7)
	x = matrix(rnorm(400), 200, 2)
	p = parts(x)
	expect_named(p, c('rhat', 'ess', 'ess_chain', 'n', 'chains'))
	expect_identical(unname(p[c('n', 'chains')]), c(100, 4))
	expect_identical(unname(p['rhat']), lame:::.lame_rhat(x))
	expect_lte(p['ess_chain'], p['ess'])
	expect_lte(p['ess'], 100)
	expect_gt(p['ess'], 50)
	# ranks only: a rescaling and a reordering of the chains change nothing
	expect_identical(parts(4 * x), p)
	expect_identical(parts(x[, 2:1]), p)
	# persistent chains hold fewer effective draws
	y = x
	for (i in 2:200) y[i, ] = 0.9 * y[i - 1, ] + sqrt(1 - 0.81) * x[i, ]
	expect_lt(parts(y)['ess'], 0.5 * p['ess'])
	# a constant series, one chain or too few draws has no R-hat
	expect_true(all(is.na(parts(matrix(1, 50, 2)))))
	expect_true(all(is.na(parts(x[, 1, drop = FALSE]))))
	expect_true(all(is.na(parts(x[1:6, ]))))
	skip_if_not_installed('posterior')
	expect_equal(unname(p['rhat']), as.numeric(posterior::rhat(x)))
})
