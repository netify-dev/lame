# exact Gibbs draw of one factor of a static bipartite product U G V'. given the
# other factor and G, each row is a Gaussian regression coefficient with
# precision I + W'W / s2 summed over observed cells, and rows are independent.
# E is the residual array of the drawn margin (transposed for V), W its design
# (V G' for U, U G for V) and s2 the residual variance; returns rows x R
.lame_bip_factor_draw <- function(E, W, s2) {
	nr <- dim(E)[1L]
	R <- ncol(W)
	N <- dim(E)[3L]
	if (!anyNA(E)) {
		ch <- chol(diag(R) + N * crossprod(W) / s2)
		rinv <- backsolve(ch, diag(R))
		K <- rowSums(E, dims = 2L) %*% W / s2
		return(K %*% tcrossprod(rinv) + matrix(stats::rnorm(nr * R), nr) %*% t(rinv))
	}
	out <- matrix(0, nr, R)
	for (i in seq_len(nr)) {
		P <- diag(R)
		b <- numeric(R)
		for (t in seq_len(N)) {
			o <- which(!is.na(E[i, , t]))
			if (length(o)) {
				Wo <- W[o, , drop = FALSE]
				P <- P + crossprod(Wo) / s2
				b <- b + drop(crossprod(Wo, E[i, o, t])) / s2
			}
		}
		ch <- chol(P)
		out[i, ] <- backsolve(ch, backsolve(ch, b, transpose = TRUE) + stats::rnorm(R))
	}
	out
}
