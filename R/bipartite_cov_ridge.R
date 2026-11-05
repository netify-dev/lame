# joint move of a nodal covariate's coefficient and the matching additive effects in a bipartite fit.
# shifting beta by -delta and the row (or column) effects by delta x leaves the likelihood unchanged,
# so delta is drawn exactly from its normal conditional under the N(0, gs) and N(0, va or vb) priors
.lame_bip_cov_ridge <- function(X, beta, a, b, gs, va, vb, rvar, cvar) {
	p <- if (length(dim(X)) == 3L) dim(X)[3L] else 0L
	for (k in seq_len(p)) {
		Xk <- X[, , k]
		if (!all(is.finite(Xk))) next
		is_row <- rvar && all(Xk == Xk[, 1L])
		is_col <- cvar && all(t(Xk) == Xk[1L, ])
		if (is_row && is_col) next
		if (is_row) {
			r <- Xk[, 1L]
			prec <- 1 / gs + sum(r^2) / va
			delta <- stats::rnorm(1L, (beta[k] / gs - sum(a * r) / va) / prec, 1 / sqrt(prec))
			beta[k] <- beta[k] - delta
			a <- a + delta * r
		} else if (is_col) {
			r <- Xk[1L, ]
			prec <- 1 / gs + sum(r^2) / vb
			delta <- stats::rnorm(1L, (beta[k] / gs - sum(b * r) / vb) / prec, 1 / sqrt(prec))
			beta[k] <- beta[k] - delta
			b <- b + delta * r
		}
	}
	list(beta = beta, a = a, b = b)
}
