# exact draw of the scale split of the static factors of a directed network.
# U V' is unchanged by U -> c U, V -> V / c, so only the prior identifies c.
# with Suv integrated out the rows have the matrix-t density
#     p(U, V) proportional to |Psi + W'W|^{-(n + df) / 2},   W = [U, V],
# and c is drawn along that orbit (Liu and Wu 1999) by a slice step on log c
# (Neal 2003). while the covariance is held at n s I (s = Psi / df) the rows
# are N(0, n s I) instead.

# `df` and `Psi` (scalar or 2R x 2R) define the inverse-Wishart prior; `hier`
# is FALSE while the covariance is held at n s I. returns U, V and c
.lame_uv_split_draw <- function(U, V, df, Psi, hier = TRUE, w = 0.5) {
	still <- list(U = U, V = V, c = 1)
	R <- ncol(U)
	n <- nrow(U)
	if (R < 1L || n < 2L || !all(is.finite(U)) || !all(is.finite(V))) return(still)
	if (hier) {
		S <- crossprod(cbind(U, V))
		Psi <- if (length(Psi) == 1L) Psi * diag(2L * R) else Psi
		if (!is.matrix(Psi) || !all(dim(Psi) == 2L * R) || !is.finite(df) || df <= 0) return(still)
		lg <- function(u) {
			d <- rep(c(exp(u), exp(-u)), each = R)
			ch <- tryCatch(chol(Psi + S * outer(d, d)), error = function(e) NULL)
			if (is.null(ch)) return(-Inf)
			-(n + df) * sum(log(diag(ch)))
		}
	} else {
		# the scale of the held covariance, as rUV_rep_fc() holds it
		s <- (if (length(Psi) == 1L) Psi else mean(diag(as.matrix(Psi)))) / df
		if (!is.finite(s) || s <= 0) return(still)
		a <- sum(U^2)
		b <- sum(V^2)
		lg <- function(u) -(exp(2 * u) * a + exp(-2 * u) * b) / (2 * n * s)
	}
	k <- exp(.lame_bso_slice(lg, w)$u)
	if (!is.finite(k) || k == 1) return(still)
	list(U = k * U, V = V / k, c = k)
}
