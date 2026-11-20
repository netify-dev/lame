# the actor block of the empirical-Bayes ALS estimator in structured form. the stacked effects
# s = (a, b) of a bipartite fit have a diagonal prior precision, those of a directed fit 2 x 2 blocks
# pairing a_i with b_i, and the information of both has a diagonal row block and a diagonal column
# block: only the row-by-column block E is dense. the penalised information K = [diag(da), E; E',
# diag(db)] is factored by eliminating the larger diagonal block, which leaves a Schur complement of
# the size of the smaller side; the coefficients (mu, beta, gamma) are eliminated after the effects.
# every result is the exact algebra of the dense inverse it replaces, and a factor that fails the
# test of the dense routine (.ae_eb_inverse(), .ae_eb_inverse_pair(), on the pivots of the same
# matrix in the same order) hands the dense matrix to that routine, so its fallbacks and notes are
# unchanged. an undirected fit, or one whose dyadic correlation couples the effects of the two actors
# of a dyad (.ae_pair_add_cross()), has a dense K and uses the dense routines

# prior precision of the stacked effects: the diagonal `d`, and for a directed fit the entries o12
# (row a_i, column b_i) and o21 of the 2 x 2 blocks
.ae_eb_precision <- function(sab, ab) {
	if (identical(ab$kind, "symmetric")) return(list(d = rep(1 / sab[1L, 1L], ab$q)))
	nr <- ab$nr; nc <- ab$nc
	if (identical(ab$kind, "bipartite")) return(list(d = c(rep(1 / sab[1L, 1L], nr), rep(1 / sab[2L, 2L], nc))))
	P <- solve(sab)
	list(d = c(rep(P[1L, 1L], nr), rep(P[2L, 2L], nr)), nr = nr, o12 = P[1L, 2L], o21 = P[2L, 1L])
}

# the precision times a vector or a matrix with q rows
.ae_prec_mult <- function(Lam, x) {
	if (is.null(Lam$o12)) return(Lam$d * x)
	ia <- seq_len(Lam$nr); ib <- Lam$nr + ia
	if (!is.matrix(x)) return(c(Lam$d[ia] * x[ia] + Lam$o12 * x[ib], Lam$o21 * x[ia] + Lam$d[ib] * x[ib]))
	out <- x
	out[ia, ] <- Lam$d[ia] * x[ia, , drop = FALSE] + Lam$o12 * x[ib, , drop = FALSE]
	out[ib, ] <- Lam$o21 * x[ia, , drop = FALSE] + Lam$d[ib] * x[ib, , drop = FALSE]
	out
}

# log determinant of the precision (the upper triangle, as chol() reads it)
.ae_prec_logdet <- function(Lam) {
	if (is.null(Lam$o12)) return(sum(log(Lam$d)))
	ia <- seq_len(Lam$nr)
	sum(log(Lam$d[ia])) + sum(log(Lam$d[Lam$nr + ia] - Lam$o12^2 / Lam$d[ia]))
}

# what the actor block is built from: the dyad sums of the weights, the shape of the fit and the
# cross weights of a dyadic correlation (Csum and its row sums d)
.ae_ak_ctx <- function(wsum, ab, cross = NULL) {
	list(wsum = wsum, kind = ab$kind, nu = ab$nu, nr = ab$nr, nc = ab$nc, cross = cross)
}

# the information of the actor effects plus the precision `Lam` (none: the information alone):
# form "block" keeps (da, db, W) and the precision's cross entries, form "dense" the matrix
.ae_ak_build <- function(ctx, Lam = NULL) {
	nr <- ctx$nr; nc <- ctx$nc
	cr <- ctx$cross
	if (identical(ctx$kind, "symmetric") || !is.null(cr)) {
		K <- .ae_eb_info_s(ctx$wsum, ctx)
		if (!is.null(cr)) {
			ia <- seq_len(nr); ib <- nr + ia
			K[ia, ia] <- K[ia, ia] + cr$Csum
			K[ib, ib] <- K[ib, ib] + cr$Csum
			K[cbind(ia, ib)] <- K[cbind(ia, ib)] + cr$d
			K[cbind(ib, ia)] <- K[cbind(ib, ia)] + cr$d
		}
		if (!is.null(Lam)) {
			dd <- cbind(seq_len(nrow(K)), seq_len(nrow(K)))
			K[dd] <- K[dd] + Lam$d
			if (!is.null(Lam$o12)) {
				ia <- seq_len(nr); ib <- nr + ia
				K[cbind(ia, ib)] <- K[cbind(ia, ib)] + Lam$o12
				K[cbind(ib, ia)] <- K[cbind(ib, ia)] + Lam$o21
			}
		}
		return(list(form = "dense", kind = ctx$kind, nr = nr, nc = nc, K = K))
	}
	rs <- rowSums(ctx$wsum); cs <- colSums(ctx$wsum)
	K <- list(form = "block", kind = ctx$kind, nr = nr, nc = nc, W = ctx$wsum, da = rs, db = cs)
	if (!is.null(Lam)) {
		K$da <- rs + Lam$d[seq_len(nr)]
		K$db <- cs + Lam$d[nr + seq_len(nc)]
		K$o12 <- Lam$o12
		K$o21 <- Lam$o21
	}
	K
}

# the matrix itself
.ae_ak_dense <- function(K) {
	if (identical(K$form, "dense")) return(K$K)
	nr <- K$nr; nc <- K$nc
	out <- rbind(cbind(diag(K$da, nr), K$W), cbind(t(K$W), diag(K$db, nc)))
	if (!is.null(K$o12)) {
		ia <- seq_len(nr)
		out[cbind(ia, nr + ia)] <- out[cbind(ia, nr + ia)] + K$o12
		out[cbind(nr + ia, ia)] <- out[cbind(nr + ia, ia)] + K$o21
	}
	out
}

# factor of a block-form K: side 1 (the larger of the row and column blocks, the rows when equal) is
# eliminated, S2 = diag(d2) - E12' diag(1 / d1) E12 is the Schur complement of side 2 with Cholesky
# factor L2, G = diag(1 / d1) E12; with `inverse`, the diagonal `dg` of K^-1 and its
# row-by-column block `cross`. `ok` is FALSE when S2 is not positive definite
.ae_ak_factor <- function(K, inverse = TRUE) {
	nr <- K$nr; nc <- K$nc
	E <- K$W
	if (!is.null(K$o12)) {
		dd <- cbind(seq_len(nr), seq_len(nr))
		E[dd] <- E[dd] + K$o12
	}
	big_a <- nr >= nc
	if (big_a) {
		d1 <- K$da; d2 <- K$db; E12 <- E
	} else {
		d1 <- K$db; d2 <- K$da; E12 <- t(E)
	}
	G <- E12 / d1
	S2 <- diag(d2, length(d2)) - crossprod(E12, G)
	L2 <- tryCatch(chol(S2), error = function(e) NULL)
	f <- list(form = "block", kind = K$kind, nr = nr, nc = nc, K = K, E = E, big_a = big_a, ok = !is.null(L2))
	if (!f$ok) return(f)
	f$d1 <- d1; f$G <- G; f$S2 <- S2; f$L2 <- L2
	f$ld <- sum(log(d1)) + 2 * sum(log(diag(L2)))
	if (inverse) {
		S2i <- chol2inv(L2)
		GS <- G %*% S2i
		dg1 <- 1 / d1 + rowSums(GS * G)
		f$S2i <- S2i
		f$dg <- if (big_a) c(dg1, diag(S2i)) else c(diag(S2i), dg1)
		f$cross <- if (big_a) -GS else -t(GS)
	}
	f
}

# the pivots (squared diagonal of the Cholesky factor) of K in its own order, rows first: those of
# the row block are da; those of the column block are the pivots of its Schur complement, the
# factor's S2 when the rows are eliminated first and otherwise a rank-nr recurrence (ae_dlr_pivots_cpp)
.ae_ak_pivots <- function(f) {
	if (f$big_a) return(c(f$K$da, diag(f$L2)^2))
	c(f$K$da, ae_dlr_pivots_cpp(f$K$db, f$E, diag(-1 / f$K$da, f$nr)))
}

# the test of .ae_eb_inverse() on the pivots of a Cholesky factorisation: all positive and the
# smallest above 1e-12 of the largest
.ae_ak_well <- function(p) all(is.finite(p)) && all(p > 0) && min(p) > 1e-12 * max(p)

# K^-1 B for a vector or a matrix with q rows, in the stacked order (rows, then columns)
.ae_ak_solve <- function(f, B) {
	vec <- !is.matrix(B)
	B <- as.matrix(B)
	ia <- seq_len(f$nr); ib <- f$nr + seq_len(f$nc)
	i1 <- if (f$big_a) ia else ib
	i2 <- if (f$big_a) ib else ia
	B1 <- B[i1, , drop = FALSE]
	x2 <- f$S2i %*% (B[i2, , drop = FALSE] - crossprod(f$G, B1))
	out <- B
	out[i1, ] <- B1 / f$d1 - f$G %*% x2
	out[i2, ] <- x2
	if (vec) as.vector(out) else out
}

# K^-1 as a q x q matrix
.ae_ak_full <- function(f) {
	n1 <- length(f$d1)
	GS <- f$G %*% f$S2i
	A11 <- tcrossprod(GS, f$G)
	A11[cbind(seq_len(n1), seq_len(n1))] <- A11[cbind(seq_len(n1), seq_len(n1))] + 1 / f$d1
	out <- rbind(cbind(A11, -GS), cbind(-t(GS), f$S2i))
	if (f$big_a) return(out)
	# side 1 holds the columns: back to rows first
	i <- c(n1 + seq_len(f$nr), seq_len(n1))
	out[i, i]
}

# log determinant of K, which fails where chol() of the matrix fails
.ae_ak_logdet <- function(K) {
	if (identical(K$form, "block")) {
		f <- .ae_ak_factor(K, inverse = FALSE)
		if (f$ok) return(f$ld)
	}
	2 * sum(log(diag(chol(.ae_ak_dense(K)))))
}

# K^-1 as an actor covariance (.ae_acov()), by .ae_eb_inverse() when K is dense or its factor fails
# that routine's test
.ae_ak_inverse <- function(K) {
	if (identical(K$form, "block")) {
		f <- .ae_ak_factor(K)
		if (f$ok && .ae_ak_well(.ae_ak_pivots(f))) return(.ae_acov(f))
	}
	.ae_acov_dense(.ae_eb_inverse(.ae_ak_dense(K)), K)
}

# the covariance of the actor effects in the variational state and the Laplace correction. it keeps
# what its users read: the diagonal `dg`, the row-by-column block `cross` (the variance of a cell's
# additive part, .ae_mf_vact()) and either the matrix `full` (dense K, a pseudo-inverse) or `terms`,
# a weighted sum of inverses of block-form K (one term, or a sum where .ae_acov_mix() cannot blend
# the precisions)
.ae_acov <- function(f) {
	list(kind = f$kind, nr = f$nr, nc = f$nc, dg = f$dg, cross = f$cross, terms = list(list(w = 1, f = f)))
}
.ae_acov_dense <- function(M, x) {
	nr <- x$nr
	list(
		kind = x$kind,
		nr = nr,
		nc = x$nc,
		dg = diag(M),
		cross = if (!identical(x$kind, "symmetric")) M[seq_len(nr), nr + seq_len(x$nc), drop = FALSE],
		full = M
	)
}

# the state a fraction `rho` of the way from A to B. two single block-form covariances blend their
# precisions, (1 - rho) K_A + rho K_B, so the state stays one block-form term however many damped
# steps follow; dense covariances, sums, and a blended K whose factor fails .ae_ak_well() mix the
# covariances (the diagonal and the cross block as the matrices do)
.ae_acov_mix <- function(A, B, rho) {
	single <- function(x) is.null(x$full) && length(x$terms) == 1L && x$terms[[1L]]$w == 1
	if (single(A) && single(B)) {
		ka <- A$terms[[1L]]$f$K; kb <- B$terms[[1L]]$f$K
		km <- ka
		for (nm in c("da", "db", "W", "o12", "o21")) {
			if (!is.null(ka[[nm]])) km[[nm]] <- (1 - rho) * ka[[nm]] + rho * kb[[nm]]
		}
		f <- .ae_ak_factor(km)
		if (isTRUE(f$ok) && .ae_ak_well(.ae_ak_pivots(f))) return(.ae_acov(f))
	}
	mix <- function(x, y) x + rho * (y - x)
	out <- list(kind = A$kind, nr = A$nr, nc = A$nc, dg = mix(A$dg, B$dg))
	if (!is.null(A$cross)) out$cross <- mix(A$cross, B$cross)
	if (is.null(A$full) && is.null(B$full)) {
		wt <- function(tm, c) lapply(tm, function(t) {
			t$w <- c * t$w
			t
		})
		out$terms <- Filter(function(t) t$w != 0, c(wt(A$terms, 1 - rho), wt(B$terms, rho)))
	} else {
		out$full <- mix(.ae_acov_full(A), .ae_acov_full(B))
	}
	out
}

.ae_acov_full <- function(A) {
	if (!is.null(A$full)) return(A$full)
	out <- 0
	for (t in A$terms) out <- out + t$w * .ae_ak_full(t$f)
	out
}

# the covariance times a vector
.ae_acov_mult <- function(A, x) {
	if (!is.null(A$full)) return(as.vector(A$full %*% x))
	out <- 0
	for (t in A$terms) out <- out + t$w * .ae_ak_solve(t$f, x)
	out
}

# log determinant: one term w K^-1 has q log w - log |K|, two terms
# w1 K1^-1 + w2 K2^-1 = K1^-1 (w1 K2 + w2 K1) K2^-1, whose middle factor is a block-form K again;
# more terms take .ae_acov_logdet_sum(), and the terms are summed as matrices where a factor fails
.ae_acov_logdet <- function(A) {
	if (!is.null(A$full)) return(.ae_mf_logdet(A$full))
	tm <- A$terms
	if (length(tm) == 1L) return(length(A$dg) * log(tm[[1L]]$w) - tm[[1L]]$f$ld)
	if (length(tm) == 2L) {
		w1 <- tm[[1L]]$w; w2 <- tm[[2L]]$w
		k1 <- tm[[1L]]$f$K; k2 <- tm[[2L]]$f$K
		km <- k1
		for (nm in c("da", "db", "W", "o12", "o21")) {
			if (!is.null(k1[[nm]])) km[[nm]] <- w1 * k2[[nm]] + w2 * k1[[nm]]
		}
		fm <- .ae_ak_factor(km, inverse = FALSE)
		if (fm$ok) return(fm$ld - tm[[1L]]$f$ld - tm[[2L]]$f$ld)
	}
	if (length(tm) > 2L) {
		ld <- .ae_acov_logdet_sum(tm)
		if (is.finite(ld)) return(ld)
	}
	.ae_mf_logdet(.ae_acov_full(A))
}

# log |M| of M = sum_t w_t K_t^-1 for J terms. side 1 of every K_t^-1 is diag(1 / d1_t) + G_t S2_t^-1 G_t',
# so M_11 = diag(Dg) + Z N Z' with Z = [G_1 ... G_J] and N = blockdiag(w_t S2_t^-1), M_12 = -Z N E and
# M_22 = E' N E, E the J stacked identities of side 2. with C = N^-1 + Z' Dg^-1 Z, log |M_11| = sum log Dg +
# log |N| + log |C| and, by the Woodbury identity N - N Z' M_11^-1 Z N = C^-1, the Schur complement of side
# 2 is E' C^-1 E. NA where a factor fails or where side 1 is not the larger dimension of Z
.ae_acov_logdet_sum <- function(tm) {
	n1 <- length(tm[[1L]]$f$d1)
	n2 <- ncol(tm[[1L]]$f$G)
	J <- length(tm)
	if (J * n2 >= n1) return(NA_real_)
	Dg <- 0
	Ni <- matrix(0, J * n2, J * n2)
	ld_n <- 0
	for (j in seq_len(J)) {
		f <- tm[[j]]$f; w <- tm[[j]]$w
		Dg <- Dg + w / f$d1
		ix <- (j - 1L) * n2 + seq_len(n2)
		Ni[ix, ix] <- f$S2 / w
		ld_n <- ld_n + n2 * log(w) - 2 * sum(log(diag(f$L2)))
	}
	Z <- do.call(cbind, lapply(tm, function(t) t$f$G))
	C <- Ni + crossprod(Z, Z / Dg)
	Lc <- tryCatch(chol((C + t(C)) / 2), error = function(e) NULL)
	if (is.null(Lc) || !all(Dg > 0)) return(NA_real_)
	Y <- backsolve(Lc, do.call(rbind, rep(list(diag(n2)), J)), transpose = TRUE)
	Ls <- tryCatch(chol(crossprod(Y)), error = function(e) NULL)
	if (is.null(Ls)) return(NA_real_)
	sum(log(Dg)) + ld_n + 2 * sum(log(diag(Lc))) + 2 * sum(log(diag(Ls)))
}

# tr(Lam Sc) for the precision `Lam` of .ae_eb_precision()
.ae_acov_trace <- function(Lam, A) {
	out <- sum(Lam$d * A$dg)
	if (is.null(Lam$o12)) return(out)
	ia <- seq_len(Lam$nr); ib <- Lam$nr + ia
	if (!is.null(A$full)) {
		return(out + Lam$o12 * sum(A$full[cbind(ia, ib)]) + Lam$o21 * sum(A$full[cbind(ib, ia)]))
	}
	out + (Lam$o12 + Lam$o21) * sum(diag(A$cross))
}

# the penalised information of (mu, beta | actor effects | gamma) as a matrix (`Lam` NULL: the
# likelihood part alone), for the joint information of the final standard errors and the dense routines
.ae_eb_dense_H <- function(ctx, Lam = NULL) {
	H <- matrix(0, ctx$m, ctx$m)
	H[ctx$i_th, ctx$i_th] <- ctx$Hth
	H[ctx$i_s, ctx$i_th] <- ctx$C
	H[ctx$i_th, ctx$i_s] <- t(ctx$C)
	H[ctx$i_s, ctx$i_s] <- .ae_ak_dense(.ae_ak_build(ctx, Lam))
	if (!is.null(Lam) && ctx$pg > 0L) {
		LM <- .ae_prec_mult(Lam, ctx$M)
		H[ctx$i_s, ctx$i_g] <- -LM
		H[ctx$i_g, ctx$i_s] <- t(-LM)
		H[ctx$i_g, ctx$i_g] <- crossprod(ctx$M, LM)
	}
	H
}

# the diagonal of the covariance of alpha = s - M gamma and, for a directed fit, its (a_i, b_i)
# entries (.ae_eb_update_sab()), from a dense inverse `S` of the information
.ae_eb_css_parts <- function(ctx, S) {
	Css <- .ae_eb_css(ctx, S)
	nr <- ctx$nr
	list(dg = diag(Css), pair = if (identical(ctx$kind, "directed")) Css[cbind(seq_len(nr), nr + seq_len(nr))])
}

# the pivots of the whole information in its own order, as chol() of .ae_eb_dense_H() gives them:
# (mu, beta), then the rows and the columns given them (K - Y Y' with Y = C Lt^-1, Lt the factor of
# their block), then gamma given everything (from `Sg`, the Schur complement of (mu, beta, gamma)
# given the effects). `f` is the factor of K, `B` the rows of (mu, beta, gamma) against the effects;
# NA where a factorisation fails
.ae_eb_system_pivots <- function(f, B, Hth, Sg) {
	K <- f$K
	nr <- f$nr; nc <- f$nc
	pa <- nrow(Hth); pg <- nrow(Sg) - pa
	th <- seq_len(pa); g <- pa + seq_len(pg)
	ia <- seq_len(nr); ib <- nr + seq_len(nc)
	Lt <- tryCatch(chol(Hth), error = function(e) NULL)
	if (is.null(Lt)) return(NA_real_)
	Yt <- backsolve(Lt, B[th, , drop = FALSE], transpose = TRUE)
	Yta <- Yt[, ia, drop = FALSE]
	p_a <- ae_dlr_pivots_cpp(K$da, Yta, -diag(pa))
	if (!.ae_ak_well(p_a)) return(p_a)
	# the columns given (mu, beta) and the rows: S2 - Z N^-1 Z' with N = I - Ya' Da^-1 Ya and
	# Z = Yb - E' Da^-1 Ya
	Gy <- t(Yta) / K$da
	Ni <- tryCatch(chol2inv(chol(diag(pa) - Yta %*% Gy)), error = function(e) NULL)
	if (is.null(Ni)) return(NA_real_)
	Z <- t(Yt[, ib, drop = FALSE]) - crossprod(f$E, Gy)
	if (f$big_a) {
		Lb <- tryCatch(chol(f$S2 - Z %*% Ni %*% t(Z)), error = function(e) NULL)
		if (is.null(Lb)) return(NA_real_)
		p_b <- diag(Lb)^2
	} else {
		P0 <- matrix(0, nr + pa, nr + pa)
		P0[cbind(ia, ia)] <- -1 / K$da
		P0[nr + th, nr + th] <- -Ni
		p_b <- ae_dlr_pivots_cpp(K$db, rbind(f$E, t(Z)), P0)
	}
	p_g <- numeric(0)
	if (pg > 0L) {
		Lg <- tryCatch(
			chol(Sg[g, g, drop = FALSE] - Sg[g, th, drop = FALSE] %*% solve(Sg[th, th, drop = FALSE], Sg[th, g, drop = FALSE])),
			error = function(e) NULL
		)
		if (is.null(Lg)) return(NA_real_)
		p_g <- diag(Lg)^2
	}
	c(diag(Lt)^2, p_a, p_b, p_g)
}

# the rows B of (mu, beta, gamma) of the information against the effects, BD = (K^-1 B')' and the Schur
# complement Sg of the coefficients given the effects, from the factor `f` of K
.ae_eb_system_parts <- function(ctx, Lam, f) {
	pa <- ctx$pa; pg <- ctx$pg
	th <- seq_len(pa); g <- pa + seq_len(pg)
	LM <- if (pg > 0L) .ae_prec_mult(Lam, ctx$M)
	B <- rbind(t(ctx$C), if (pg > 0L) -t(LM))
	Hrr <- matrix(0, pa + pg, pa + pg)
	Hrr[th, th] <- ctx$Hth
	if (pg > 0L) Hrr[g, g] <- crossprod(ctx$M, LM)
	BD <- t(.ae_ak_solve(f, t(B)))
	list(B = B, BD = BD, Sg = Hrr - tcrossprod(BD, B))
}

# the inverse of the penalised information of (mu, beta | actor effects | gamma) at the precision
# `Lam`, as what its users need: `solve(rhs)`, the solution of the normal equations; `css`, the
# parts of the covariance of alpha (.ae_eb_css_parts()); and with `schur`, `Sc`, the covariance of
# the effects given the coefficients (.ae_acov()). with `schur` the test is that of
# .ae_eb_inverse_pair() (the effect block, then the coefficients given the effects), without it that of
# .ae_eb_inverse() on the whole matrix in its own order (mu and beta, rows, columns, gamma)
.ae_eb_system_inverse <- function(ctx, Lam, schur = TRUE) {
	s <- ctx$i_s
	r <- c(ctx$i_th, ctx$i_g)
	dense <- function() {
		H <- .ae_eb_dense_H(ctx, Lam)
		if (schur) {
			inv <- .ae_eb_inverse_pair(H, s)
			S <- inv$S
		} else {
			S <- .ae_eb_inverse(H)
		}
		list(
			solve = function(rhs) as.vector(S %*% rhs),
			css = .ae_eb_css_parts(ctx, S),
			Sc = if (schur) .ae_acov_dense(inv$Sc, ctx)
		)
	}
	K <- .ae_ak_build(ctx, Lam)
	if (identical(K$form, "dense")) return(dense())
	f <- .ae_ak_factor(K)
	if (!f$ok) return(dense())
	pK <- .ae_ak_pivots(f)
	if (schur && !.ae_ak_well(pK)) return(dense())
	pa <- ctx$pa; pg <- ctx$pg
	th <- seq_len(pa); g <- pa + seq_len(pg)
	nr <- ctx$nr; nc <- ctx$nc
	ia <- seq_len(nr); ib <- nr + seq_len(nc)
	sp <- .ae_eb_system_parts(ctx, Lam, f)
	B <- sp$B; BD <- sp$BD; Sg <- sp$Sg
	chS <- tryCatch(chol((Sg + t(Sg)) / 2), error = function(e) NULL)
	if (is.null(chS)) return(dense())
	if (schur) {
		if (!.ae_ak_well(c(pK, diag(chS)^2))) return(dense())
	} else if (!.ae_ak_well(.ae_eb_system_pivots(f, B, ctx$Hth, Sg))) {
		return(dense())
	}
	Srr <- chol2inv(chS)
	# the covariance of alpha: K^-1 + Gm Srr Gm', Gm the effects' loadings on (mu, beta, gamma)
	Gm <- t(BD)
	if (pg > 0L) Gm[, g] <- Gm[, g] + ctx$M
	GS <- Gm %*% Srr
	list(
		solve = function(rhs) {
			xr <- as.vector(Srr %*% (rhs[r] - BD %*% rhs[s]))
			x <- numeric(length(rhs))
			x[r] <- xr
			x[s] <- .ae_ak_solve(f, rhs[s]) - as.vector(crossprod(BD, xr))
			x
		},
		css = list(
			dg = f$dg + rowSums(GS * Gm),
			pair = if (identical(ctx$kind, "directed")) diag(f$cross) + rowSums(GS[ia, , drop = FALSE] * Gm[ib, , drop = FALSE])
		),
		Sc = if (schur) .ae_acov(f)
	)
}
