# static actor effects under sum-to-zero constraints, with the intercept held fixed
.lda_update_static_ab_centered <- function(Z, W, Xb, Oarr, symmetric = FALSE) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	nt <- dim(Z)[3L]
	observed <- is.finite(Z) & is.finite(W) & W > 0
	weight <- W
	weight[!observed] <- 0
	base <- array(0, dim(Z))
	for (tt in seq_len(nt)) {
		x <- Z[, , tt] - Xb[, , tt] - .lda_o_slice(Oarr, tt)
		x[!observed[, , tt]] <- 0
		base[, , tt] <- weight[, , tt] * x
	}
	wsum <- apply(weight, c(1, 2), sum)
	rhs_a <- apply(base, 1L, sum)
	rhs_b <- apply(base, 2L, sum)
	contrast <- function(n) rbind(diag(n - 1L), rep(-1, n - 1L))
	qa <- contrast(nr)
	if (symmetric) {
		H <- diag(rowSums(wsum) + colSums(wsum)) + wsum + t(wsum)
		a <- as.vector(
			qa %*% .ae_safe_solve(
				crossprod(qa, H %*% qa),
				as.vector(crossprod(qa, rhs_a + rhs_b))
			)
		)
		b <- a
	} else {
		qb <- contrast(nc)
		A <- crossprod(qa, rowSums(wsum) * qa)
		B <- crossprod(qb, colSums(wsum) * qb)
		C <- crossprod(qa, wsum %*% qb)
		H <- rbind(cbind(A, C), cbind(t(C), B))
		theta <- .ae_safe_solve(H, c(crossprod(qa, rhs_a), crossprod(qb, rhs_b)))
		a <- as.vector(qa %*% theta[seq_len(nr - 1L)])
		b <- as.vector(qb %*% theta[nr - 1L + seq_len(nc - 1L)])
	}
	list(a = matrix(a, nr, nt), b = matrix(b, nc, nt), level = 0)
}
