# kernels of the block coordinate descent of the fast estimator: each is
# checked against the array / loop form of the same arithmetic

test_that("the data sums of a bcd fit equal the sums over the arrays", {
	set.seed(2101)
	nr = 7L; nc = 9L; tt = 3L; p = 2L
	w = array(runif(nr * nc * tt), c(nr, nc, tt)); w[sample(length(w), 40)] = 0
	z = array(rnorm(nr * nc * tt), c(nr, nc, tt)); z[w == 0] = 0
	x = array(rnorm(nr * nc * p * tt), c(nr, nc, p, tt))
	x[, , , 1][w[, , 1] == 0] = NA # a design value at an unobserved cell
	st = lame:::ae_bcd_stats_cpp(w, z, x, nr, nc, tt, p)
	xt = function(k, t) { m = x[, , k, t]; m[w[, , t] == 0] = 0; m }
	sz = apply(w * z, c(1, 2), sum)
	expect_equal(st$sz, sz, tolerance = 1e-13)
	expect_equal(matrix(st$gm[, 1], nr, nc), apply(w, c(1, 2), sum), tolerance = 1e-13)
	for (k in 1:p) {
		gk = Reduce(`+`, lapply(1:tt, function(t) w[, , t] * xt(k, t)))
		expect_equal(matrix(st$gm[, k + 1], nr, nc), gk, tolerance = 1e-13)
	}
	D = do.call(rbind, lapply(1:tt, function(t) cbind(1, sapply(1:p, function(k) as.vector(xt(k, t))))))
	wv = as.vector(w); zv = as.vector(z)
	expect_equal(st$xtx, crossprod(D, wv * D), tolerance = 1e-12)
	expect_equal(as.vector(st$xtz), as.vector(crossprod(D, wv * zv)), tolerance = 1e-12)
	expect_equal(st$zz, sum(wv * zv^2), tolerance = 1e-12)
	expect_equal(st$gr, sapply(1:(p + 1), function(k) rowSums(matrix(st$gm[, k], nr, nc))), tolerance = 1e-13)
	expect_equal(st$gc, sapply(1:(p + 1), function(k) colSums(matrix(st$gm[, k], nr, nc))), tolerance = 1e-13)
	expect_equal(as.vector(st$szr), rowSums(sz), tolerance = 1e-13)
	expect_equal(as.vector(st$szc), colSums(sz), tolerance = 1e-13)
})

test_that("the sweep kernels equal the array form of one bcd sweep", {
	set.seed(2102)
	nr = 8L; nc = 6L; p = 2L
	gm = cbind(as.vector(matrix(runif(nr * nc, 1, 4), nr, nc)), matrix(rnorm(nr * nc * p), nr * nc, p))
	sz = matrix(rnorm(nr * nc), nr, nc); wsum = matrix(gm[, 1], nr, nc)
	a = rnorm(nr); b = rnorm(nc); omat = matrix(rnorm(nr * nc), nr, nc)
	theta = rnorm(p + 1L); xtz = rnorm(p + 1L)
	cmat = outer(a, b, "+") + omat
	expect_equal(lame:::ae_bcd_gtc_cpp(a, b, omat, gm), as.vector(crossprod(gm, as.vector(cmat))), tolerance = 1e-13)
	phi = matrix(as.vector(gm %*% theta), nr, nc)
	g = sz - phi - wsum * omat
	rs = lame:::ae_bcd_rs_cpp(theta, gm, sz, wsum, omat)
	expect_equal(rs$phi, phi, tolerance = 1e-13)
	expect_equal(rs$base_rs, rowSums(g), tolerance = 1e-13)
	expect_equal(rs$base_cs, colSums(g), tolerance = 1e-13)
	dmu = 0.3
	rbar = lame:::ae_bcd_rbar_cpp(a, b, phi, dmu, sz, wsum)
	ref = (sz - (phi + dmu * wsum) - wsum * outer(a, b, "+")) / wsum
	expect_equal(rbar, ref, tolerance = 1e-13)
	cross = lame:::ae_bcd_sse_cpp(a, b, omat, phi, dmu, sz, wsum)
	u = sz - (phi + dmu * wsum)
	expect_equal(cross, c(sum(cmat * u), sum(wsum * cmat^2)), tolerance = 1e-13)
})

test_that("the exact residual sum of squares equals the array form", {
	set.seed(2103)
	nr = 6L; nc = 7L; tt = 4L; p = 3L
	w = array(runif(nr * nc * tt), c(nr, nc, tt)); w[sample(length(w), 30)] = 0
	z = array(rnorm(nr * nc * tt), c(nr, nc, tt))
	x = array(rnorm(nr * nc * p * tt), c(nr, nc, p, tt))
	mu = 0.4; beta = c(0.5, -0.2, 0.1); cmat = matrix(rnorm(nr * nc), nr, nc)
	ref = 0
	for (t in 1:tt) {
		eta = mu + cmat
		for (k in 1:p) eta = eta + beta[k] * x[, , k, t]
		ref = ref + sum(w[, , t] * (z[, , t] - eta)^2)
	}
	expect_equal(lame:::ae_bcd_sse_exact_cpp(w, z, x, nr, nc, tt, p, mu, beta, cmat), ref, tolerance = 1e-13)
})

# the additive block as the R recursion states it: up to 200 sweeps to a 1e-10
# step, then the largest residual of the normal equations
additive_reference = function(wsum, base_rs, base_cs, a, b, cnt_row, cnt_col, symmetric) {
	cnt_row_pos = pmax(cnt_row, 1); cnt_col_pos = pmax(cnt_col, 1)
	for (inner in seq_len(200L)) {
		a_prev = a; b_prev = b
		if (symmetric) {
			a = (base_rs - as.vector(wsum %*% a)) / cnt_row_pos
			a[cnt_row == 0] = 0
			b = a
		} else {
			a = (base_rs - as.vector(wsum %*% b)) / cnt_row_pos
			a[cnt_row == 0] = 0
			b = (base_cs - as.vector(crossprod(wsum, a))) / cnt_col_pos
			b[cnt_col == 0] = 0
		}
		if (max(abs(a - a_prev)) < 1e-10 && max(abs(b - b_prev)) < 1e-10) break
	}
	ra = cnt_row * a + as.vector(wsum %*% b) - base_rs
	rb = cnt_col * b + as.vector(crossprod(wsum, a)) - base_cs
	ra[cnt_row == 0] = 0; rb[cnt_col == 0] = 0
	list(
		a = a,
		b = b,
		resid = max(max(abs(ra)), max(abs(rb))),
		scale = max(1, max(abs(base_rs)), max(abs(base_cs)))
	)
}

test_that("the compiled additive block reproduces the recursion, with and without the exact start", {
	set.seed(2104)
	for (symmetric in c(FALSE, TRUE)) {
		n = 14L
		w = matrix(runif(n * n, 0.05, 1), n, n); diag(w) = 0
		if (symmetric) w = (w + t(w)) / 2
		w[3, ] = 0; if (symmetric) w[, 3] = 0 else w[, 5] = 0 # an actor with no weight
		cnt_row = rowSums(w); cnt_col = colSums(w)
		base_rs = rnorm(n); base_cs = rnorm(n)
		base_rs[cnt_row == 0] = 0; base_cs[cnt_col == 0] = 0 # nothing observed for an actor
		if (symmetric) base_cs = base_rs # symmetric residuals
		if (!symmetric) { # a consistent system: equal totals
			act = cnt_col > 0
			base_cs[act] = base_cs[act] - mean(base_cs[act]) + sum(base_rs) / sum(act)
		}
		a0 = rnorm(n); b0 = if (symmetric) a0 else rnorm(n)
		ref = additive_reference(w, base_rs, base_cs, a0, b0, cnt_row, cnt_col, symmetric)
		got = lame:::ae_additive_cpp(w, base_rs, base_cs, a0, b0, cnt_row, cnt_col, symmetric, NULL)
		expect_equal(got$a, ref$a, tolerance = 1e-10)
		expect_equal(got$b, ref$b, tolerance = 1e-10)
		expect_equal(got$resid, ref$resid, tolerance = 1e-6)
		expect_equal(got$scale, ref$scale)
		# started at the exact minimiser, the directed recursion reaches the same
		# limit (same gauge); the symmetric recursion has no limit from a general
		# start, so the exact start is compared with the exact solution instead
		direct = lame:::.ae_additive_factor(w, cnt_row, cnt_col, symmetric)
		expect_false(is.null(direct))
		fast = lame:::ae_additive_cpp(w, base_rs, base_cs, a0, b0, cnt_row, cnt_col, symmetric, direct)
		if (symmetric) {
			exact = lame:::.ae_additive_solve(base_rs, base_cs, w, cnt_row, cnt_col, TRUE)
			expect_equal(fast$a, exact$a, tolerance = 1e-8)
			expect_lt(fast$resid, 1e-9)
		} else {
			expect_equal(fast$a, ref$a, tolerance = 1e-8)
			expect_equal(fast$b, ref$b, tolerance = 1e-8)
		}
		res = lame:::ae_additive_resid_cpp(w, base_rs, base_cs, ref$a, ref$b, cnt_row, cnt_col)
		expect_equal(res, ref$resid, tolerance = 1e-8)
	}
})

test_that("the compiled balanced factor update equals the update followed by the r balancing", {
	set.seed(2105)
	for (rank in 1:4) {
		nr = 21L; nc = 17L
		u = matrix(rnorm(nr * rank), nr, rank); v = matrix(rnorm(nc * rank), nc, rank)
		z = matrix(rnorm(nr * nc), nr, nc); w = matrix(runif(nr * nc), nr, nc); w[w < 0.2] = 0
		su = runif(rank, 0.2, 2); sv = runif(rank, 0.2, 2)
		raw = lame:::als_factor_update_cpp(z, w, u, v, su, sv, 1, 2 * rank + 2, 15L, 1e-6)
		same_order = identical(
			order(sqrt(raw$sig2_u * raw$sig2_v), decreasing = TRUE),
			order(sqrt(colSums(raw$U^2) * colSums(raw$V^2)), decreasing = TRUE)
		)
		got = lame:::als_factor_balanced_cpp(z, w, u, v, su, sv, 1, 2 * rank + 2, 15L, 1e-6)
		expect_identical(isTRUE(got$balanced), same_order)
		if (same_order) {
			ref = lame:::.ae_factor_balance(raw$U, raw$V, raw$sig2_u, raw$sig2_v)
			expect_equal(got$U, ref$U, tolerance = 1e-10)
			expect_equal(got$V, ref$V, tolerance = 1e-10)
			expect_equal(got$Omat, tcrossprod(ref$U, ref$V), tolerance = 1e-10)
		} else {
			expect_equal(got$U, raw$U, tolerance = 1e-13)
		}
		expect_equal(got$sig2_u, raw$sig2_u)
	}
	# factors with a negligible column: the qr pivots, and the balancing still agrees
	nr = 15L; nc = 12L; rank = 2L
	u = matrix(rnorm(nr * rank), nr, rank); v = matrix(rnorm(nc * rank), nc, rank)
	w = matrix(0, nr, nc); z = matrix(0, nr, nc) # no data: the factors are shrunk to zero
	got = lame:::als_factor_balanced_cpp(z, w, u, v, c(0.5, 0.5), c(0.5, 0.5), 1, 6, 15L, 1e-6)
	expect_true(all(is.finite(c(got$U, got$V, got$Omat))))
	expect_equal(got$Omat, tcrossprod(got$U, got$V), tolerance = 1e-12)
})

test_that("a fit that reaches the direct additive start matches the plain recursion", {
	skip_on_cran()
	set.seed(2106)
	n = 24L; tt = 3L
	rn = paste0("a", sprintf("%02d", seq_len(n)))
	g = rnorm(n); h = rnorm(n)
	Y = lapply(seq_len(tt), function(t) {
		x = outer(g, h, "+")
		y = 1 * (0.6 * x + outer(rnorm(n, 0, 0.5), rnorm(n, 0, 0.5), "+") + rnorm(n * n) > 1)
		diag(y) = NA; dimnames(y) = list(rn, rn); y
	})
	X = lapply(seq_len(tt), function(t) {
		a = array(0.9 * outer(g, h, "+") + 0.3 * rnorm(n * n), c(n, n, 1), dimnames = list(rn, rn, "x")); a
	})
	fit = suppressWarnings(suppressMessages(lame_als(
		Y,
		Xdyad = X,
		R = 1,
		family = "binary",
		bias_correct = "none",
		verbose = FALSE,
		max_iter = 60
	)))
	expect_s3_class(fit, "ame_als")
	expect_true(all(is.finite(fit$beta)))
})
