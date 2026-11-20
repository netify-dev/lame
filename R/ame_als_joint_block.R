# the joint information of (mu, beta | effects | gamma | factors) of the final standard errors of the
# empirical-Bayes ALS estimator in structured form. on the longer side of a bipartite or directed fit
# (the rows when the sides are equal) every actor's effect and its R factor entries couple only to
# each other and to the rest: the coefficients, the other side's effects and factors and gamma. so the
# information is block diagonal there, with blocks D_i of size R + 1 bordered by B, and the rest has the
# Schur complement Sch = Hrr - B D^-1 B' of size n_short (R + 1) + pa + pg. its inverse is kept as
# D^-1, Gt = D^-1 B', Srr = Sch^-1 and X = Gt Srr, from which the users read entries, solves and
# submatrices (.ae_ji_ent(), .ae_ji_mult(), .ae_ji_sub()) without the matrix of size m + R q.
# .ae_eb_inverse() of the dense matrix is used where its Cholesky pivots, in the matrix's own order
# (.ae_jb_pivots()), fail its test, and for an undirected fit or one with a dyadic correlation, whose
# effects and factors couple across actors (.ae_eb_extend())

# an actor block exists without an undirected effect or a dyadic correlation
.ae_jb_blockform <- function(ctx) !identical(ctx$kind, "symmetric") && is.null(ctx$cross)

# the information at the fitted state `fit` with R factor columns (0: the coefficients and effects
# alone), as .ae_eb_extend() forms it, and its block form when the shape has one
.ae_jb_system <- function(ctx, ab, fit, R = 0L, vb = NULL) {
	Lam <- .ae_eb_precision(fit$sab, ab)
	sys <- list(ctx = ctx, ab = ab, fit = fit, R = R, vb = vb, Lam = Lam)
	if (.ae_jb_blockform(ctx)) sys$J <- .ae_jb_info(ctx, Lam, fit, R, vb)
	sys
}

# the dense matrix (stacked order: mu, beta, effects, gamma, factors)
.ae_jb_dense <- function(sys) {
	H0 <- .ae_eb_dense_H(sys$ctx, sys$Lam)
	if (sys$R > 0L) .ae_eb_extend(sys$ctx, sys$ab, H0, sys$fit, sys$vb) else H0
}

# the blocks: D (k x k x n1, k = R + 1: the effect, then the factor entries), the border Bt (n1 k x ns,
# the k rows of an actor together) and the rest Hrr (ns x ns), the rest in its stacked order: mu and
# beta, the short side's effects, gamma, its factors. each entry is the upper-triangle entry of the
# dense matrix, which is what chol() reads. `l_nat` and `r_nat` place the blocks in the stacked vector
.ae_jb_info <- function(ctx, Lam, fit = NULL, R = 0L, vb = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; pa <- ctx$pa; pg <- ctx$pg; q <- ctx$q; m0 <- ctx$m
	W <- ctx$wsum
	K <- .ae_ak_build(ctx, Lam)
	E <- W
	if (!is.null(K$o12)) {
		dd <- cbind(seq_len(nr), seq_len(nr))
		E[dd] <- E[dd] + K$o12
	}
	LM <- if (pg > 0L) .ae_prec_mult(Lam, ctx$M)
	Gk <- lapply(seq_len(pa), function(k) matrix(ctx$gm[, k], nr, nc))
	ia <- seq_len(nr); ib <- nr + seq_len(nc)
	rows <- nr >= nc
	pb <- if (R > 0L && !is.null(vb)) .ae_mf_prior_blocks(vb$S$prior, R)
	joint <- R > 0L && !is.null(vb) && identical(vb$S$prior$kind, "joint")
	sub_rows <- function(A, i) if (is.null(A)) NULL else A[i, , drop = FALSE]
	if (rows) {
		n1 <- nr; n2 <- nc
		v <- list(W = W, E = E, d1 = K$da, d2 = K$db, C1 = ctx$C[ia, , drop = FALSE], C2 = ctx$C[ib, , drop = FALSE],
			LM1 = sub_rows(LM, ia), LM2 = sub_rows(LM, ib), G = Gk)
		if (R > 0L) {
			v <- c(v, list(F1 = fit$U, F2 = fit[["V"]], Cs = vb$S$Om, Cl = vb$S$Sig, p11 = pb$uu, p22 = pb$vv,
				p12 = if (joint) pb$uv, s1 = fit$sig2_u, s2 = fit$sig2_v))
		}
	} else {
		n1 <- nc; n2 <- nr
		v <- list(W = t(W), E = t(E), d1 = K$db, d2 = K$da, C1 = ctx$C[ib, , drop = FALSE], C2 = ctx$C[ia, , drop = FALSE],
			LM1 = sub_rows(LM, ib), LM2 = sub_rows(LM, ia), G = lapply(Gk, t))
		if (R > 0L) {
			v <- c(v, list(F1 = fit[["V"]], F2 = fit$U, Cs = vb$S$Sig, Cl = vb$S$Om, p11 = pb$vv, p22 = pb$uu,
				p12 = if (joint) t(pb$uv), s1 = fit$sig2_v, s2 = fit$sig2_u))
		}
	}
	k <- R + 1L
	ns <- pa + n2 + pg + R * n2
	r_th <- seq_len(pa); r_e <- pa + seq_len(n2); r_g <- pa + n2 + seq_len(pg)
	r_f <- function(s) pa + n2 + pg + (s - 1L) * n2 + seq_len(n2)
	rl <- function(sl) (seq_len(n1) - 1L) * k + 1L + sl
	Bt <- matrix(0, n1 * k, ns)
	Bt[rl(0L), r_th] <- v$C1
	Bt[rl(0L), r_e] <- v$E
	D <- array(0, c(k, k, n1))
	D[1L, 1L, ] <- v$d1
	Hrr <- matrix(0, ns, ns)
	Hrr[r_th, r_th] <- ctx$Hth
	Hrr[r_th, r_e] <- t(v$C2)
	Hrr[cbind(r_e, r_e)] <- v$d2
	if (pg > 0L) {
		Bt[rl(0L), r_g] <- -v$LM1
		Hrr[r_e, r_g] <- -v$LM2
		Hrr[r_g, r_g] <- crossprod(ctx$M, LM)
	}
	for (r in seq_len(R)) {
		f1 <- v$F1[, r]; f2 <- v$F2[, r]
		D[1L, 1L + r, ] <- as.vector(v$W %*% f2)
		for (kk in seq_len(pa)) {
			Bt[rl(r), r_th[kk]] <- as.vector(v$G[[kk]] %*% f2)
			Hrr[r_th[kk], r_f(r)] <- as.vector(crossprod(v$G[[kk]], f1))
		}
		Bt[rl(r), r_e] <- v$W * rep(f2, each = n1)
		Bt[rl(0L), r_f(r)] <- v$W * f1
		Hrr[cbind(r_e, r_f(r))] <- as.vector(crossprod(v$W, f1))
		for (s in seq_len(R)) {
			# factor r of a long-side actor against factor s of a short-side one pairs the other entries
			# of the two ranks: w_ij F1_is F2_jr
			Bt[rl(r), r_f(s)] <- v$W * outer(v$F1[, s], f2)
			if (!is.null(v$p12)) {
				ix <- cbind(rl(r), r_f(s))
				Bt[ix] <- Bt[ix] + v$p12[r, s]
			}
			if (s < r) next
			d1 <- as.vector(v$W %*% (f2 * v$F2[, s]))
			d2 <- as.vector(crossprod(v$W, f1 * v$F1[, s]))
			if (!is.null(vb)) {
				d1 <- d1 + as.vector(v$W %*% v$Cs[r, s, ]) + v$p11[r, s]
				d2 <- d2 + as.vector(crossprod(v$W, v$Cl[r, s, ])) + v$p22[r, s]
			} else if (r == s) {
				d1 <- d1 + 1 / v$s1[r]
				d2 <- d2 + 1 / v$s2[r]
			}
			D[1L + r, 1L + s, ] <- d1
			D[1L + s, 1L + r, ] <- d1
			Hrr[cbind(r_f(r), r_f(s))] <- d2
		}
	}
	for (r in seq_len(R)) D[1L + r, 1L, ] <- D[1L, 1L + r, ]
	Hrr[lower.tri(Hrr)] <- t(Hrr)[lower.tri(Hrr)]
	if (rows) {
		l_nat <- cbind(pa + seq_len(nr), if (R > 0L) m0 + outer(seq_len(nr), (seq_len(R) - 1L) * nr, "+"))
		nat_e <- pa + nr + seq_len(nc)
		nat_f <- if (R > 0L) m0 + R * nr + seq_len(R * nc)
	} else {
		l_nat <- cbind(pa + nr + seq_len(nc), if (R > 0L) m0 + R * nr + outer(seq_len(nc), (seq_len(R) - 1L) * nc, "+"))
		nat_e <- pa + seq_len(nr)
		nat_f <- if (R > 0L) m0 + seq_len(R * nr)
	}
	nat_g <- pa + q + seq_len(pg)
	M <- ctx$M
	list(
		rows = rows,
		n1 = n1,
		n2 = n2,
		k = k,
		R = R,
		m0 = m0,
		d = m0 + R * q,
		ns = ns,
		D = D,
		Bt = Bt,
		Hrr = Hrr,
		l_nat = l_nat,
		r_nat = c(seq_len(pa), nat_e, nat_g, nat_f),
		nat_e = nat_e,
		nat_g = nat_g,
		M1 = if (pg > 0L) M[if (rows) ia else ib, , drop = FALSE],
		M2 = if (pg > 0L) M[if (rows) ib else ia, , drop = FALSE]
	)
}

# the system with only the rest variables `keep` (positions in the rest)
.ae_jb_keep <- function(J, keep) {
	J$Bt <- J$Bt[, keep, drop = FALSE]
	J$Hrr <- J$Hrr[keep, keep, drop = FALSE]
	J$r_nat <- J$r_nat[keep]
	J$ns <- length(keep)
	J
}

# factor of the blocks: each D_i by Cholesky (`ld1` the sum of their log determinants), the Schur
# complement Sch of the rest and its Cholesky factor `Lr`; with `inverse`, Srr = Sch^-1, X = Gt Srr
# and the diagonal blocks Sll of the inverse on the long side. `ok` is FALSE where a factor fails
.ae_jb_factor <- function(J, inverse = TRUE) {
	bl <- ae_jb_blocks_cpp(J$D, J$Bt)
	if (bl$bad > 0L) return(list(ok = FALSE))
	Sch <- J$Hrr - crossprod(J$Bt, bl$Gt)
	Sch <- (Sch + t(Sch)) / 2
	Lr <- tryCatch(chol(Sch), error = function(e) NULL)
	f <- list(ok = !is.null(Lr), ld1 = bl$ld, Dinv = bl$Dinv, Gt = bl$Gt, Sch = Sch, Lr = Lr)
	if (!f$ok || !inverse) return(f)
	k <- J$k
	f$Srr <- chol2inv(Lr)
	f$X <- f$Gt %*% f$Srr
	Sll <- f$Dinv
	at <- function(a) seq.int(a, by = k, length.out = J$n1)
	for (a in seq_len(k)) for (b in a:k) {
		x <- Sll[a, b, ] + rowSums(f$X[at(a), , drop = FALSE] * f$Gt[at(b), , drop = FALSE])
		Sll[a, b, ] <- x
		Sll[b, a, ] <- x
	}
	f$Sll <- Sll
	f
}

# the pivots (squared diagonal of the Cholesky factor) of the information in its own order, as chol() of
# .ae_jb_dense() gives them: those of the coefficients and effects (`p0`, .ae_eb_h0_pivots()), then the
# factor entries. the rest variables that precede the long side's factors (mu, beta, the short side's
# effects, gamma and, when the columns are the long side, the row factors) have the Schur complement Sc0
# given the long side's effects, whose own pivots are those of its factor; the long side's factors take
# the recurrence of ae_jb_pivots_cpp() on Sc0; the factors after them are pivots of Sch. NA where a
# factorisation fails
.ae_jb_pivots <- function(J, f, p0) {
	if (J$R == 0L) return(p0)
	lead <- which(J$r_nat < min(J$l_nat[, -1L]))
	B0 <- J$Bt[seq.int(1L, by = J$k, length.out = J$n1), lead, drop = FALSE]
	Sc0 <- J$Hrr[lead, lead, drop = FALSE] - crossprod(B0, B0 / J$D[1L, 1L, ])
	Lc <- tryCatch(chol((Sc0 + t(Sc0)) / 2), error = function(e) NULL)
	if (is.null(Lc)) return(NA_real_)
	p_lead <- diag(Lc)^2
	p_long <- ae_jb_pivots_cpp(J$D, J$Bt[, lead, drop = FALSE], chol2inv(Lc), J$R)
	after <- setdiff(seq_len(J$ns), lead)
	c(p0, p_lead[J$r_nat[lead] > J$m0], p_long, diag(f$Lr)[after]^2)
}

# the inverse of the information of `sys` (.ae_jb_system()) as an object for .ae_ji_*(). .ae_eb_inverse()
# takes the Cholesky inverse when the pivots pass its test, else that of the equilibrated matrix (unit
# diagonal), whose pivots are the pivots over the diagonal: either is the inverse the block factor
# gives. a matrix that fails both goes to that routine
.ae_jb_inverse <- function(sys) {
	J <- sys$J
	if (!is.null(J)) {
		f <- .ae_jb_factor(J)
		if (f$ok) {
			p <- .ae_jb_pivots(J, f, .ae_eb_h0_pivots(sys$ctx, sys$Lam))
			if (.ae_ak_well(p)) return(.ae_ji_block(J, f))
			h <- numeric(J$d)
			h[J$l_nat] <- t(apply(J$D, 3L, diag))
			h[J$r_nat] <- diag(J$Hrr)
			if (length(p) == J$d && all(is.finite(h)) && all(h > 0) && .ae_ak_well(p / h)) return(.ae_ji_block(J, f))
		}
	}
	.ae_ji_dense(.ae_eb_inverse(.ae_jb_dense(sys)))
}

# the inverse of the information of the effects and factors alone (mu, beta and gamma left out), the
# covariance given the coefficients: from the block form when its factor is positive definite, else
# from the dense inverse `S` of the whole information, as the Schur complement of the coefficients
.ae_jb_inverse_z <- function(sys, S) {
	J <- sys$J
	ctx <- sys$ctx
	if (!is.null(J) && identical(S$form, "block")) {
		Jz <- .ae_jb_keep(J, which(!(J$r_nat %in% c(ctx$i_th, ctx$i_g))))
		f <- .ae_jb_factor(Jz)
		if (f$ok) return(.ae_ji_block(Jz, f))
	}
	Sf <- .ae_ji_full(S)
	iz <- setdiff(seq_len(nrow(Sf)), c(ctx$i_th, ctx$i_g))
	it <- c(ctx$i_th, ctx$i_g)
	Sz <- Sf[iz, iz, drop = FALSE] -
		Sf[iz, it, drop = FALSE] %*% solve(Sf[it, it, drop = FALSE], Sf[it, iz, drop = FALSE])
	.ae_ji_dense(Sz, iz, nrow(Sf))
}

# log determinant of the information of the effects and factors alone, which fails where chol() of the
# dense matrix fails
.ae_jb_logdet_z <- function(sys) {
	J <- sys$J
	ctx <- sys$ctx
	if (!is.null(J)) {
		f <- .ae_jb_factor(.ae_jb_keep(J, which(!(J$r_nat %in% c(ctx$i_th, ctx$i_g)))), inverse = FALSE)
		if (f$ok) return(f$ld1 + 2 * sum(log(diag(f$Lr))))
	}
	H <- .ae_jb_dense(sys)
	iz <- setdiff(seq_len(nrow(H)), c(ctx$i_th, ctx$i_g))
	2 * sum(log(diag(chol(H[iz, iz, drop = FALSE]))))
}

# the pivots of the information of the coefficients and effects in its own order (.ae_eb_system_pivots())
.ae_eb_h0_pivots <- function(ctx, Lam) {
	f <- .ae_ak_factor(.ae_ak_build(ctx, Lam))
	if (!f$ok) return(NA_real_)
	sp <- .ae_eb_system_parts(ctx, Lam, f)
	.ae_eb_system_pivots(f, sp$B, ctx$Hth, sp$Sg)
}

# the inverse as an object: `form` "block" (from .ae_jb_factor()) or "dense" (a matrix `S` whose rows are
# the stacked positions `map` of a vector of length `d`); `d` the size of the stacked vector
.ae_ji_block <- function(J, f) {
	d <- J$d; k <- J$k
	ln <- as.vector(t(J$l_nat))
	lrow <- ridx <- act <- slot <- integer(d)
	lrow[ln] <- seq_along(ln)
	act[ln] <- rep(seq_len(J$n1), each = k)
	slot[ln] <- rep(seq_len(k), J$n1)
	ridx[J$r_nat] <- seq_len(J$ns)
	list(
		form = "block",
		d = d,
		k = k,
		n1 = J$n1,
		ln = ln,
		rn = J$r_nat,
		lrow = lrow,
		ridx = ridx,
		act = act,
		slot = slot,
		Dinv = f$Dinv,
		Gt = f$Gt,
		X = f$X,
		Srr = f$Srr,
		Sll = f$Sll,
		info = list(rows = J$rows, n2 = J$n2, r_e = match(J$nat_e, J$r_nat), r_g = match(J$nat_g, J$r_nat), M1 = J$M1, M2 = J$M2)
	)
}
.ae_ji_dense <- function(S, map = NULL, d = nrow(S)) {
	pos <- integer(d)
	pos[if (is.null(map)) seq_len(d) else map] <- seq_len(nrow(S))
	out <- list(form = "dense", d = d, S = S, pos = pos)
	attr(out, "ae_singular") <- attr(S, "ae_singular")
	out
}

# the matrix itself (a dense inverse; tests and the fallback of .ae_jb_inverse_z())
.ae_ji_full <- function(S) {
	if (identical(S$form, "dense")) return(S$S)
	.ae_ji_sub(S, seq_len(S$d), seq_len(S$d))
}

# the entries S[i[k], j[k]] (stacked positions)
.ae_ji_ent <- function(S, i, j) {
	if (identical(S$form, "dense")) return(S$S[cbind(S$pos[i], S$pos[j])])
	li <- S$lrow[i]; lj <- S$lrow[j]; ri <- S$ridx[i]; rj <- S$ridx[j]
	out <- numeric(length(i))
	w <- ri > 0L & rj > 0L
	out[w] <- S$Srr[cbind(ri[w], rj[w])]
	w <- li > 0L & rj > 0L
	out[w] <- -S$X[cbind(li[w], rj[w])]
	w <- ri > 0L & lj > 0L
	out[w] <- -S$X[cbind(lj[w], ri[w])]
	w <- which(li > 0L & lj > 0L)
	if (length(w)) {
		ai <- S$act[i[w]]; aj <- S$act[j[w]]
		same <- ai == aj
		ws <- w[same]
		out[ws] <- S$Sll[cbind(S$slot[i[ws]], S$slot[j[ws]], ai[same])]
		wd <- w[!same]
		if (length(wd)) out[wd] <- rowSums(S$X[li[wd], , drop = FALSE] * S$Gt[lj[wd], , drop = FALSE])
	}
	out
}

# S %*% x for a vector or a matrix with d rows
.ae_ji_mult <- function(S, x) {
	vec <- !is.matrix(x)
	x <- as.matrix(x)
	out <- matrix(0, S$d, ncol(x))
	if (identical(S$form, "dense")) {
		at <- which(S$pos > 0L)
		xs <- matrix(0, nrow(S$S), ncol(x))
		xs[S$pos[at], ] <- x[at, , drop = FALSE]
		out[at, ] <- (S$S %*% xs)[S$pos[at], , drop = FALSE]
		return(if (vec) as.vector(out) else out)
	}
	xl <- x[S$ln, , drop = FALSE]
	yr <- S$Srr %*% (x[S$rn, , drop = FALSE] - crossprod(S$Gt, xl))
	yl <- -(S$Gt %*% yr)
	k <- S$k
	for (a in seq_len(k)) {
		ra <- seq.int(a, by = k, length.out = S$n1)
		for (b in seq_len(k)) {
			rb <- seq.int(b, by = k, length.out = S$n1)
			yl[ra, ] <- yl[ra, ] + S$Dinv[a, b, ] * xl[rb, , drop = FALSE]
		}
	}
	out[S$ln, ] <- yl
	out[S$rn, ] <- yr
	if (vec) as.vector(out) else out
}

# the submatrix S[i, j]; a long-side by long-side part is formed as a product, so it is kept small
.ae_ji_sub <- function(S, i, j) {
	if (identical(S$form, "dense")) return(S$S[S$pos[i], S$pos[j], drop = FALSE])
	out <- matrix(0, length(i), length(j))
	li <- S$lrow[i]; lj <- S$lrow[j]; ri <- S$ridx[i]; rj <- S$ridx[j]
	il <- which(li > 0L); ir <- which(ri > 0L); jl <- which(lj > 0L); jr <- which(rj > 0L)
	if (length(ir) && length(jr)) out[ir, jr] <- S$Srr[ri[ir], rj[jr], drop = FALSE]
	if (length(il) && length(jr)) out[il, jr] <- -S$X[li[il], rj[jr], drop = FALSE]
	if (length(ir) && length(jl)) out[ir, jl] <- -t(S$X[lj[jl], ri[ir], drop = FALSE])
	if (length(il) && length(jl)) {
		blk <- S$X[li[il], , drop = FALSE] %*% t(S$Gt[lj[jl], , drop = FALSE])
		same <- which(outer(S$act[i[il]], S$act[j[jl]], "=="), arr.ind = TRUE)
		if (length(same)) {
			blk[same] <- S$Sll[cbind(S$slot[i[il]][same[, 1L]], S$slot[j[jl]][same[, 2L]], S$act[i[il]][same[, 1L]])]
		}
		out[il, jl] <- blk
	}
	out
}
