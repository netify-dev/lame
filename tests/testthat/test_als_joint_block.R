skip_on_cran()

# the structured joint information of the final standard errors (R/ame_als_joint_block.R) against
# the dense matrix of .ae_eb_extend() and its inverse, on small tall, wide and directed shapes

kn = function(name) get(name, envir = asNamespace("lame"))

# a joint information at a random state: two periods, one dyadic covariate, node covariates on both
# sides; `vb` adds the variational row covariances and prior blocks of .ae_mf_finish()
jb_case = function(kind, nr, nc, R = 2, seed = 1, vb = FALSE, prior = "scalar") {
	set.seed(seed)
	Tt = 2L; p = 1L
	if (kind == "directed") nc = nr
	W = array(runif(nr * nc * Tt, 0.05, 0.3), c(nr, nc, Tt))
	if (kind == "directed") for (t in 1:Tt) diag(W[, , t]) = 0
	W[sample(length(W), round(0.1 * length(W)))] = 0
	X = array(rnorm(nr * nc * Tt), c(nr, nc, p, Tt))
	Wr = cbind(rnorm(nr)); Wc = cbind(rnorm(nc))
	ab = list(kind = kind, nr = nr, nc = nc, q = nr + nc, nu = 1, Wr = Wr, Wc = Wc, eta0 = 3)
	st = lame:::ae_bcd_stats_cpp(W, W, X, nr, nc, Tt, p, 1L)
	ctx = kn(".ae_eb_context")(ab, st, rowSums(W, dims = 2L), p)
	sab = if (kind == "directed") matrix(c(0.5, 0.15, 0.15, 0.4), 2) else diag(c(0.5, 0.4))
	fit = list(sab = sab, U = matrix(rnorm(nr * R, 0, 0.8), nr), V = matrix(rnorm(nc * R, 0, 0.8), nc),
		sig2_u = rep(0.7, R), sig2_v = rep(0.9, R))
	v = NULL
	if (vb) {
		rc = function(n) array(vapply(1:n, function(i) { A = matrix(rnorm(R * R, 0, 0.2), R); crossprod(A) + diag(0.05, R) }, matrix(0, R, R)), c(R, R, n))
		pr = if (prior == "joint") {
			A = matrix(rnorm(4 * R * R, 0, 0.3), 2 * R)
			list(kind = "joint", Lam = crossprod(A) + diag(1, 2 * R))
		} else list(kind = "scalar", tau = 0.8)
		v = list(S = list(Om = rc(nc), Sig = rc(nr), prior = pr))
	}
	list(ctx = ctx, ab = ab, fit = fit, R = R, vb = v)
}

dense_of = function(cs) {
	sys = kn(".ae_jb_system")(cs$ctx, cs$ab, cs$fit, cs$R, cs$vb)
	list(sys = sys, H = kn(".ae_jb_dense")(sys))
}

shapes = list(
	tall = list(kind = "bipartite", nr = 30, nc = 6),
	wide = list(kind = "bipartite", nr = 5, nc = 24),
	directed = list(kind = "directed", nr = 12, nc = 12)
)

test_that("the block form, its pivots and its log determinant are those of the dense matrix", {
	for (nm in names(shapes)) for (R in c(0L, 1L, 3L)) for (vb in c(FALSE, TRUE)) {
		sh = shapes[[nm]]
		cs = jb_case(sh$kind, sh$nr, sh$nc, R = R, seed = 3 + R, vb = vb && R > 0L, prior = if (sh$kind == "directed") "joint" else "scalar")
		dh = dense_of(cs)
		J = dh$sys$J
		lab = paste(nm, R, vb)
		# the blocks in the stacked order are the dense matrix
		Hb = matrix(0, J$d, J$d)
		ln = J$l_nat
		for (i in seq_len(J$n1)) Hb[ln[i, ], ln[i, ]] = J$D[, , i]
		rows = as.vector(t(ln))
		Hb[rows, J$r_nat] = J$Bt
		Hb[J$r_nat, rows] = t(J$Bt)
		Hb[J$r_nat, J$r_nat] = J$Hrr
		expect_equal(Hb, dh$H, tolerance = 1e-13, label = lab)
		f = kn(".ae_jb_factor")(J)
		expect_true(f$ok)
		p0 = kn(".ae_eb_h0_pivots")(cs$ctx, dh$sys$Lam)
		expect_equal(kn(".ae_jb_pivots")(J, f, p0), diag(chol(dh$H))^2, tolerance = 1e-10, label = lab)
		iz = setdiff(seq_len(J$d), c(cs$ctx$i_th, cs$ctx$i_g))
		expect_equal(kn(".ae_jb_logdet_z")(dh$sys), 2 * sum(log(diag(chol(dh$H[iz, iz])))), tolerance = 1e-12, label = lab)
	}
})

test_that("entries, solves and submatrices of the block inverse are those of the dense inverse", {
	for (nm in names(shapes)) {
		sh = shapes[[nm]]
		cs = jb_case(sh$kind, sh$nr, sh$nc, R = 2L, seed = 11, vb = nm == "tall")
		dh = dense_of(cs)
		S = kn(".ae_jb_inverse")(dh$sys)
		expect_identical(S$form, "block")
		Sd = solve(dh$H)
		d = nrow(Sd)
		expect_equal(kn(".ae_ji_full")(S), Sd, tolerance = 1e-10, label = nm)
		i = sample(d, 200, TRUE); j = sample(d, 200, TRUE)
		expect_equal(kn(".ae_ji_ent")(S, i, j), Sd[cbind(i, j)], tolerance = 1e-10, label = nm)
		x = matrix(rnorm(d * 3), d)
		expect_equal(kn(".ae_ji_mult")(S, x), Sd %*% x, tolerance = 1e-10, label = nm)
		expect_equal(kn(".ae_ji_mult")(S, x[, 1]), as.vector(Sd %*% x[, 1]), tolerance = 1e-10, label = nm)
		ii = sample(d, 15); jj = sample(d, 9)
		expect_equal(kn(".ae_ji_sub")(S, ii, jj), Sd[ii, jj], tolerance = 1e-10, label = nm)
		# given the coefficients: the inverse of the effect and factor block
		iz = setdiff(seq_len(d), c(cs$ctx$i_th, cs$ctx$i_g))
		Sz = kn(".ae_jb_inverse_z")(dh$sys, S)
		expect_equal(kn(".ae_ji_sub")(Sz, iz, iz), solve(dh$H[iz, iz]), tolerance = 1e-10, label = nm)
		# the alpha covariance blocks, the factor covariance parts and the sab information
		ctx = cs$ctx
		P = kn(".ae_eb_css")(ctx, Sd)
		bl = kn(".ae_eb_css_blocks")(S, ctx)
		ia = seq_len(ctx$nr); ib = ctx$nr + seq_len(ctx$nc)
		full = function(x) if (is.matrix(x)) x else diag(x$dg, length(x$dg)) + tcrossprod(x$ZS, x$Z)
		expect_equal(full(bl$aa), P[ia, ia], tolerance = 1e-10, label = nm)
		expect_equal(bl$ab, P[ia, ib], tolerance = 1e-10, label = nm)
		expect_equal(full(bl$bb), P[ib, ib], tolerance = 1e-10, label = nm)
		Sdo = kn(".ae_ji_dense")(Sd)
		expect_equal(kn(".ae_eb_sab_info")(S, ctx, cs$ab, cs$fit$sab), kn(".ae_eb_sab_info")(Sdo, ctx, cs$ab, cs$fit$sab), tolerance = 1e-9, label = nm)
		m0 = ctx$m
		pt = kn(".ae_eb_factor_cov_parts")(S, cs$fit$U, cs$fit$V, rep(1, 2), FALSE, m0)
		jf = (m0 + 1):d
		Sc = kn(".ae_eb_factor_cov")(Sd[jf, jf], kn(".ae_eb_gauge_dirs")(cs$fit$U, cs$fit$V, rep(1, 2), FALSE), c(rep(ctx$nr, 2), rep(ctx$nc, 2)))
		at_u = function(r) (r - 1) * ctx$nr + ia
		at_v = function(r) 2 * ctx$nr + (r - 1) * ctx$nc + seq_len(ctx$nc)
		for (r in 1:2) for (s in 1:2) {
			expect_equal(pt$uu[[r]][[s]], diag(Sc[at_u(r), at_u(s)]), tolerance = 1e-9, label = nm)
			expect_equal(pt$vv[[r]][[s]], diag(Sc[at_v(r), at_v(s)]), tolerance = 1e-9, label = nm)
			expect_equal(pt$uv[[r]][[s]], Sc[at_u(r), at_v(s)], tolerance = 1e-9, label = nm)
		}
	}
})

test_that("a matrix whose pivots fail the test takes the route of the dense routine", {
	# a factor column of zeros with a tight prior: the pivots fail the first test of .ae_eb_inverse() and
	# pass that of the equilibrated matrix, whose inverse is the one of the block factor
	cs = jb_case("bipartite", 20, 5, R = 2L, seed = 5)
	cs$fit$U[, 1] = 0
	cs$fit$V[, 1] = 0
	cs$fit$sig2_u = c(1e-9, 0.7); cs$fit$sig2_v = c(1e9, 0.9)
	dh = dense_of(cs)
	pv = diag(chol(dh$H))^2
	expect_false(min(pv) > 1e-12 * max(pv))
	f = kn(".ae_jb_factor")(dh$sys$J)
	expect_equal(kn(".ae_jb_pivots")(dh$sys$J, f, kn(".ae_eb_h0_pivots")(cs$ctx, dh$sys$Lam)), pv, tolerance = 1e-8)
	S = kn(".ae_jb_inverse")(dh$sys)
	expect_identical(S$form, "block")
	Sd = kn(".ae_eb_inverse")(dh$H)
	expect_equal(kn(".ae_ji_full")(S), Sd, tolerance = 1e-8)
	# a collinear node covariate: singular, so the dense pseudo-inverse with its mark
	cs = jb_case("bipartite", 20, 5, R = 1L, seed = 6)
	cs$ab$Wr = cbind(cs$ab$Wr, 2 * cs$ab$Wr)
	cs$ctx = kn(".ae_eb_context")(cs$ab, list(xtx = cs$ctx$Hth, gr = cs$ctx$C[1:20, , drop = FALSE], gc = cs$ctx$C[21:25, , drop = FALSE], gm = cs$ctx$gm), cs$ctx$wsum, 1L)
	dh = dense_of(cs)
	S = kn(".ae_jb_inverse")(dh$sys)
	expect_identical(S$form, "dense")
	Sd = kn(".ae_eb_inverse")(dh$H)
	expect_false(is.null(attr(Sd, "ae_singular")))
	expect_identical(attr(S, "ae_singular"), attr(Sd, "ae_singular"))
	expect_equal(S$S, Sd)
	# an undirected shape has no block form
	ctx = cs$ctx; ctx$kind = "symmetric"
	expect_false(kn(".ae_jb_blockform")(ctx))
})

test_that("the gauge traces of the effect noise variance are those of the dense gauge inverse", {
	set.seed(7)
	for (dims in list(c(40, 9), c(6, 30), c(15, 15))) {
		W = matrix(runif(prod(dims), 0, 2), dims[1])
		W[sample(length(W), 20)] = 0
		rs = rowSums(W); cs = colSums(W)
		K = rbind(cbind(diag(rs), W), cbind(t(W), diag(cs)))
		Si = kn(".ae_gauge_inverse")(K, c(rep(1, dims[1]), rep(-1, dims[2])))
		tr = kn(".ae_gauge_traces")(rs, cs, W)
		ia = seq_len(dims[1]); ib = dims[1] + seq_len(dims[2])
		expect_equal(c(tr$a$trace, tr$a$sum, tr$b$trace, tr$b$sum),
			c(sum(diag(Si[ia, ia])), sum(Si[ia, ia]), sum(diag(Si[ib, ib])), sum(Si[ib, ib])), tolerance = 1e-10)
	}
	# a disconnected design fails the pivot test and keeps the pseudo-inverse
	W = matrix(0, 8, 6); W[1:4, 1:3] = 1; W[5:8, 4:6] = 1
	expect_null(kn(".ae_gauge_traces")(rowSums(W), colSums(W), W))
})

test_that("the log determinant of a blend of three or more actor covariances is that of the summed matrix", {
	for (dims in list(c(40, 5), c(4, 36))) {
		Cs = lapply(1:4, function(s) {
			cs = jb_case("bipartite", dims[1], dims[2], R = 0L, seed = 20 + s)
			kn(".ae_ak_inverse")(kn(".ae_ak_build")(cs$ctx, kn(".ae_eb_precision")(cs$fit$sab, cs$ab)))
		})
		# the sum of two weighted inverses a blend keeps when its blended precision fails the factor test; blends
		# from a sum add terms
		wt = function(C, w) { t = C$terms[[1]]; t$w = w; t }
		B = Cs[[1]]
		B$dg = B$dg + 0.5 * (Cs[[2]]$dg - B$dg)
		B$cross = B$cross + 0.5 * (Cs[[2]]$cross - B$cross)
		B$terms = list(wt(Cs[[1]], 0.5), wt(Cs[[2]], 0.5))
		M = kn(".ae_acov_full")(Cs[[1]])
		M = M + 0.5 * (kn(".ae_acov_full")(Cs[[2]]) - M)
		for (s in 3:4) {
			B = kn(".ae_acov_mix")(B, Cs[[s]], 0.3 + 0.1 * s)
			M = M + (0.3 + 0.1 * s) * (kn(".ae_acov_full")(Cs[[s]]) - M)
			expect_equal(kn(".ae_acov_logdet_sum")(B$terms), 2 * sum(log(diag(chol(M)))), tolerance = 1e-10)
		}
		expect_length(B$terms, 4L)
		expect_equal(kn(".ae_acov_logdet")(B), 2 * sum(log(diag(chol(M)))), tolerance = 1e-10)
	}
})
