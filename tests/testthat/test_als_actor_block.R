skip_on_cran()

library(lame)

# the structured actor block of the empirical-Bayes ALS estimator (R/ame_als_actor_block.R) against the
# dense algebra it replaces, on small random cases of every shape: tall and wide bipartite, directed,
# undirected, and directed with the cross weights of a dyadic correlation

kn = function(f) get(f, envir = asNamespace("lame"))

ab_case = function(kind, nr, nc, pgn = 2L, seed = 1L, cross = FALSE) {
	set.seed(seed)
	if (kind != "bipartite") nc = nr
	Y = array(rbinom(nr * nc, 1, .3), c(nr, nc, 1)); if (kind != "bipartite") diag(Y[, , 1]) = NA
	p = 2L; X = array(rnorm(nr * nc * p), c(nr, nc, p, 1))
	prep = list(
		Yarr = Y,
		bip = kind == "bipartite",
		nr = nr,
		nc = nc,
		Tt = 1L,
		p_dyad = p,
		X = X,
		W_row = if (pgn) matrix(rnorm(nr * pgn), nr) else NULL,
		W_col = if (pgn && kind == "bipartite") matrix(rnorm(nc * pgn), nc) else NULL
	)
	ab = kn(".ae_eb_prior")(prep, kind == "symmetric")
	W = array(runif(nr * nc, .02, .3), c(nr, nc, 1)); Z = array(rnorm(nr * nc), c(nr, nc, 1))
	if (kind != "bipartite") { W[cbind(1:nr, 1:nr, 1)] = 0; Z[cbind(1:nr, 1:nr, 1)] = 0 }
	if (kind == "symmetric") { W[, , 1] = (W[, , 1] + t(W[, , 1])) / 2 }
	st = kn("ae_bcd_stats_cpp")(W, Z, X, nr, nc, 1L, p, 1L)
	ctx = kn(".ae_eb_context")(ab, st, rowSums(W, dims = 2L), p)
	if (cross) {
		cw = array(-runif(nr * nr, 0, .05), c(nr, nr, 1)); cw[, , 1] = (cw[, , 1] + t(cw[, , 1])) / 2; diag(cw[, , 1]) = 0
		ctx = kn(".ae_pair_add_cross")(ctx, kn(".ae_pair_cross")(list(nr = nr, Tt = 1L, p = p, X = X), cw))
	}
	sab = if (kind == "directed") matrix(c(.9, -.3, -.3, .6), 2) else diag(c(.8, .5))
	list(ab = ab, ctx = ctx, Lam = kn(".ae_eb_precision")(sab, ab), st = st)
}

shapes = list(
	tall = list("bipartite", 17L, 6L),
	wide = list("bipartite", 5L, 13L),
	directed = list("directed", 9L, 9L),
	symmetric = list("symmetric", 9L, 9L)
)

test_that("the precision acts as its dense matrix", {
	for (sh in shapes) {
		cs = ab_case(sh[[1]], sh[[2]], sh[[3]])
		Ld = .ae_prec_dense(cs$Lam)
		x = rnorm(nrow(Ld)); Xm = matrix(rnorm(3 * nrow(Ld)), ncol = 3)
		expect_equal(kn(".ae_prec_mult")(cs$Lam, x), as.vector(Ld %*% x), tolerance = 1e-13)
		expect_equal(kn(".ae_prec_mult")(cs$Lam, Xm), Ld %*% Xm, tolerance = 1e-13)
		expect_equal(kn(".ae_prec_logdet")(cs$Lam), kn(".ae_mf_logdet")(Ld), tolerance = 1e-12)
	}
})

test_that("the factor of the actor block gives the dense solve, inverse, log determinant and pivots", {
	for (nm in c("tall", "wide", "directed")) {
		sh = shapes[[nm]]
		cs = ab_case(sh[[1]], sh[[2]], sh[[3]], seed = 3L)
		K = kn(".ae_ak_build")(cs$ctx, cs$Lam)
		expect_identical(K$form, "block")
		Kd = kn(".ae_ak_dense")(K)
		# the dense matrix is the actor block of the dense information
		Hd = kn(".ae_eb_dense_H")(cs$ctx, cs$Lam)
		expect_identical(Kd, Hd[cs$ctx$i_s, cs$ctx$i_s])
		f = kn(".ae_ak_factor")(K)
		Ki = chol2inv(chol(Kd))
		B = matrix(rnorm(nrow(Kd) * 2), ncol = 2)
		expect_equal(kn(".ae_ak_solve")(f, B), Ki %*% B, tolerance = 1e-12, label = nm)
		expect_equal(f$dg, diag(Ki), tolerance = 1e-12)
		expect_equal(f$cross, Ki[seq_len(K$nr), K$nr + seq_len(K$nc)], tolerance = 1e-12)
		expect_equal(kn(".ae_ak_full")(f), Ki, tolerance = 1e-12)
		expect_equal(f$ld, 2 * sum(log(diag(chol(Kd)))), tolerance = 1e-12)
		expect_equal(kn(".ae_ak_logdet")(K), f$ld)
		expect_equal(kn(".ae_ak_pivots")(f), diag(chol(Kd))^2, tolerance = 1e-12, label = nm)
	}
	# the rank-r recurrence of the pivots of D + Z' P Z
	set.seed(4); D = runif(12, 2, 3); Zt = matrix(rnorm(36), 3); P = -crossprod(matrix(rnorm(9), 3)) / 200
	expect_equal(kn("ae_dlr_pivots_cpp")(D, Zt, P), diag(chol(diag(D) + t(Zt) %*% P %*% Zt))^2, tolerance = 1e-12)
	# it stops at the first pivot that is not positive
	piv = kn("ae_dlr_pivots_cpp")(c(1, 1, 1), matrix(c(0, 2, 0), 1), matrix(-1))
	expect_true(piv[2] <= 0 && is.na(piv[3]))
})

test_that("the inverse of the penalised information is the dense inverse for every shape", {
	cases = c(shapes, list(directed_cross = list("directed", 9L, 9L, TRUE)))
	for (nm in names(cases)) for (pgn in c(0L, 2L)) for (schur in c(TRUE, FALSE)) {
		sh = cases[[nm]]
		cs = ab_case(sh[[1]], sh[[2]], sh[[3]], pgn = pgn, seed = 5L, cross = length(sh) > 3L)
		ctx = cs$ctx
		Hd = kn(".ae_eb_dense_H")(ctx, cs$Lam)
		Si = solve(Hd)
		inv = kn(".ae_eb_system_inverse")(ctx, cs$Lam, schur = schur)
		rhs = rnorm(ctx$m)
		lab = paste(nm, pgn, schur)
		expect_equal(inv$solve(rhs), as.vector(Si %*% rhs), tolerance = 1e-10, label = lab)
		Css = kn(".ae_eb_css")(ctx, Si)
		expect_equal(inv$css$dg, diag(Css), tolerance = 1e-10, label = lab)
		if (sh[[1]] == "directed") {
			expect_equal(inv$css$pair, Css[cbind(1:ctx$nr, ctx$nr + 1:ctx$nr)], tolerance = 1e-10, label = lab)
		}
		if (schur) {
			Sc = solve(Hd[ctx$i_s, ctx$i_s])
			expect_equal(inv$Sc$dg, diag(Sc), tolerance = 1e-10)
			expect_equal(kn(".ae_acov_full")(inv$Sc), Sc, tolerance = 1e-10)
			expect_equal(kn(".ae_acov_logdet")(inv$Sc), -2 * sum(log(diag(chol(Hd[ctx$i_s, ctx$i_s])))), tolerance = 1e-10)
			x = rnorm(ctx$q)
			expect_equal(kn(".ae_acov_mult")(inv$Sc, x), as.vector(Sc %*% x), tolerance = 1e-10)
		}
	}
})

test_that("the pivots of the whole information are those of its dense Cholesky factor", {
	for (nm in c("tall", "wide", "directed")) for (pgn in c(0L, 2L)) {
		sh = shapes[[nm]]
		cs = ab_case(sh[[1]], sh[[2]], sh[[3]], pgn = pgn, seed = 7L)
		ctx = cs$ctx
		f = kn(".ae_ak_factor")(kn(".ae_ak_build")(ctx, cs$Lam))
		Hd = kn(".ae_eb_dense_H")(ctx, cs$Lam)
		r = c(ctx$i_th, ctx$i_g)
		B = Hd[r, ctx$i_s, drop = FALSE]
		Sg = Hd[r, r] - B %*% solve(Hd[ctx$i_s, ctx$i_s], t(B))
		expect_equal(kn(".ae_eb_system_pivots")(f, B, ctx$Hth, Sg), diag(chol(Hd))^2, tolerance = 1e-10, label = paste(nm, pgn))
	}
})

test_that("a blend of block-form actor covariances blends their precisions and stays one term", {
	for (nm in c("tall", "wide", "directed")) {
		sh = shapes[[nm]]
		a1 = ab_case(sh[[1]], sh[[2]], sh[[3]], seed = 8L); a2 = ab_case(sh[[1]], sh[[2]], sh[[3]], seed = 9L)
		K1 = kn(".ae_ak_build")(a1$ctx, a1$Lam); K2 = kn(".ae_ak_build")(a2$ctx, a1$Lam); K3 = kn(".ae_ak_build")(a1$ctx, a2$Lam)
		C1 = kn(".ae_ak_inverse")(K1); C2 = kn(".ae_ak_inverse")(K2); C3 = kn(".ae_ak_inverse")(K3)
		D1 = kn(".ae_ak_dense")(K1); D2 = kn(".ae_ak_dense")(K2); D3 = kn(".ae_ak_dense")(K3)
		ld = function(M) 2 * sum(log(diag(chol(M))))
		# a full step, a damped step and a damped step from a blend: (1 - rho) K_S + rho K_T each time
		B1 = kn(".ae_acov_mix")(C1, C2, 1)
		expect_equal(kn(".ae_acov_full")(B1), solve(D2), tolerance = 1e-10)
		B2 = kn(".ae_acov_mix")(C1, C2, 0.3)
		P2 = 0.7 * D1 + 0.3 * D2
		B3 = kn(".ae_acov_mix")(B2, C3, 0.5)
		P3 = 0.5 * P2 + 0.5 * D3
		B4 = kn(".ae_acov_mix")(B3, C2, 0.4)
		P4 = 0.6 * P3 + 0.4 * D2
		for (bl in list(list(B2, P2), list(B3, P3), list(B4, P4))) {
			M = solve(bl[[2]])
			expect_length(bl[[1]]$terms, 1L)
			expect_equal(bl[[1]]$terms[[1]]$w, 1)
			expect_equal(kn(".ae_acov_logdet")(bl[[1]]), -ld(bl[[2]]), tolerance = 1e-10, label = nm)
			expect_equal(kn(".ae_acov_full")(bl[[1]]), M, tolerance = 1e-10)
			expect_equal(bl[[1]]$dg, diag(M), tolerance = 1e-10)
			expect_equal(kn(".ae_mf_vact")(bl[[1]], a1$ab), kn(".ae_mf_vact")(kn(".ae_acov_dense")(M, a1$ab), a1$ab), tolerance = 1e-10)
			expect_equal(kn(".ae_acov_trace")(a2$Lam, bl[[1]]), sum(.ae_prec_dense(a2$Lam) * M), tolerance = 1e-10)
		}
	}
})

test_that("a blend of dense or summed actor covariances keeps the trace and the log determinant of the blended matrix", {
	ld = function(M) 2 * sum(log(diag(chol(M))))
	check = function(bl, M, ab, Lam, nm) {
		expect_equal(kn(".ae_acov_logdet")(bl), ld(M), tolerance = 1e-10, label = nm)
		expect_equal(kn(".ae_acov_full")(bl), M, tolerance = 1e-12)
		expect_equal(bl$dg, diag(M), tolerance = 1e-12)
		expect_equal(kn(".ae_mf_vact")(bl, ab), kn(".ae_mf_vact")(kn(".ae_acov_dense")(M, ab), ab), tolerance = 1e-12)
		expect_equal(kn(".ae_acov_trace")(Lam, bl), sum(.ae_prec_dense(Lam) * M), tolerance = 1e-12)
	}
	for (nm in c("tall", "wide", "directed", "symmetric")) {
		sh = shapes[[nm]]
		a1 = ab_case(sh[[1]], sh[[2]], sh[[3]], seed = 8L); a2 = ab_case(sh[[1]], sh[[2]], sh[[3]], seed = 9L)
		C1 = kn(".ae_ak_inverse")(kn(".ae_ak_build")(a1$ctx, a1$Lam))
		C2 = kn(".ae_ak_inverse")(kn(".ae_ak_build")(a2$ctx, a1$Lam))
		C3 = kn(".ae_ak_inverse")(kn(".ae_ak_build")(a1$ctx, a2$Lam))
		F1 = kn(".ae_acov_full")(C1); F2 = kn(".ae_acov_full")(C2); F3 = kn(".ae_acov_full")(C3)
		M2 = F1 + 0.3 * (F2 - F1)
		# a dense covariance mixes as a matrix; a block-form start is a sum of two weighted inverses (the state a
		# blend takes when the blended precision fails its factor test), and blends from it add terms
		B2 = if (is.null(C1$full)) {
			wt = function(C, w) { t = C$terms[[1]]; t$w = w; t }
			list(kind = C1$kind, nr = C1$nr, nc = C1$nc, dg = C1$dg + 0.3 * (C2$dg - C1$dg),
				cross = if (!is.null(C1$cross)) C1$cross + 0.3 * (C2$cross - C1$cross), terms = list(wt(C1, 0.7), wt(C2, 0.3)))
		} else {
			kn(".ae_acov_mix")(C1, C2, 0.3)
		}
		B3 = kn(".ae_acov_mix")(B2, C3, 0.5)
		M3 = M2 + 0.5 * (F3 - M2)
		B4 = kn(".ae_acov_mix")(B3, C2, 0.4)
		M4 = M3 + 0.4 * (F2 - M3)
		if (is.null(C1$full)) expect_length(B4$terms, 4L)
		for (bl in list(list(B2, M2), list(B3, M3), list(B4, M4))) check(bl[[1]], bl[[2]], a1$ab, a2$Lam, nm)
	}
})

test_that("a singular information takes the dense routine and its pseudo-inverse", {
	cs = ab_case("bipartite", 15L, 6L, pgn = 2L, seed = 10L)
	ctx = cs$ctx
	# a row covariate collinear with the other makes gamma unidentified
	ctx$M[, 2] = 2 * ctx$M[, 1]
	for (schur in c(TRUE, FALSE)) {
		inv = kn(".ae_eb_system_inverse")(ctx, cs$Lam, schur = schur)
		Hd = kn(".ae_eb_dense_H")(ctx, cs$Lam)
		S = if (schur) kn(".ae_eb_inverse_pair")(Hd, ctx$i_s)$S else kn(".ae_eb_inverse")(Hd)
		expect_type(attr(kn(".ae_eb_inverse")(Hd), "ae_singular"), "character")
		rhs = rnorm(ctx$m)
		expect_equal(inv$solve(rhs), as.vector(S %*% rhs), tolerance = 1e-12)
		expect_equal(inv$css$dg, diag(kn(".ae_eb_css")(ctx, S)), tolerance = 1e-12)
	}
})
