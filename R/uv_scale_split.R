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

# the rest of the gauge: U -> U A, V -> V A^-T leaves U V' unchanged for every invertible A, and the split above
# moves only A = c I. the column scales A = I + (c - 1) e_k e_k' and the shears A = I + h e_k e_l' are
# one-parameter groups with unit Jacobian and Haar measure d log c and dh, so a slice step on each draws along
# its orbit under the same prior (Liu and Sabatti 2000); together they reach every A of positive determinant.
# W'W moves to B'W'WB with B = diag(A, A^-T), so a step costs a few 2R x 2R Cholesky factors

# arguments as for .lame_uv_split_draw; returns U, V and the number of steps that moved
.lame_uv_gauge_draw <- function(U, V, df, Psi, hier = TRUE, w = 0.5) {
	still <- list(U = U, V = V, moved = 0L)
	R <- ncol(U)
	n <- nrow(U)
	if (R < 2L || n < 2L || !all(is.finite(U)) || !all(is.finite(V))) return(still)
	if (hier) {
		Psi <- if (length(Psi) == 1L) Psi * diag(2L * R) else Psi
		if (!is.matrix(Psi) || !all(dim(Psi) == 2L * R) || !is.finite(df) || df <= 0) return(still)
		# Psi + S is positive definite, so the log determinant needs no guard
		lp <- function(S) -0.5 * (n + df) * determinant(Psi + S)$modulus[[1L]]
	} else {
		s <- (if (length(Psi) == 1L) Psi else mean(diag(as.matrix(Psi)))) / df
		if (!is.finite(s) || s <= 0) return(still)
		lp <- function(S) -sum(diag(S)) / (2 * n * s)
	}
	S <- crossprod(cbind(U, V))
	moved <- 0L
	# column scales: U_k -> c U_k, V_k -> V_k / c
	for (k in seq_len(R)) {
		pos <- c(k, R + k)
		lg <- function(u) {
			d <- rep(1, 2L * R)
			d[pos] <- c(exp(u), exp(-u))
			lp(S * outer(d, d))
		}
		u <- .lame_bso_slice(lg, w)$u
		if (!is.finite(u) || u == 0) next
		d <- rep(1, 2L * R)
		d[pos] <- c(exp(u), exp(-u))
		U[, k] <- U[, k] * d[k]
		V[, k] <- V[, k] * d[R + k]
		S <- S * outer(d, d)
		moved <- moved + 1L
	}
	# shears: U_l -> U_l + h U_k, V_k -> V_k - h V_l; U_k and V_l stay, so the width in h is fixed along the orbit
	for (k in seq_len(R)) for (l in seq_len(R)[-k]) {
		if (!(S[k, k] > 0) || !(S[R + l, R + l] > 0)) next
		shear <- function(h) {
			B <- diag(2L * R)
			B[k, l] <- h
			B[R + l, R + k] <- -h
			B
		}
		lg_shear <- function(h) { B <- shear(h); lp(crossprod(B, S %*% B)) }
		h <- .lame_bso_slice(lg_shear, w * sqrt(S[R + l, R + l] / S[k, k]))$u
		if (!is.finite(h) || h == 0) next
		B <- shear(h)
		U[, l] <- U[, l] + h * U[, k]
		V[, k] <- V[, k] - h * V[, l]
		S <- crossprod(B, S %*% B)
		moved <- moved + 1L
	}
	list(U = U, V = V, moved = moved)
}

# the scale split of the static factors of a symmetric network, V = U L with L diagonal: U_r -> c U_r,
# L_r -> L_r / c^2 (so V_r -> V_r / c) leaves U L U' unchanged. with the precision of the column integrated
# out (Gamma(1, 1/2), as rUV_sym_fc() draws it) U_r has density (1 + ||U_r||^2)^(-(n + 2) / 2), N(0, n I)
# while it is held at 1 / n; L_r ~ N(0, l_var), and the map has Jacobian c^(n - 2). each log c is drawn by
# a slice step along its orbit

# `l_var` the prior variance of an eigen-weight, `hier` FALSE while the precision is held; returns U, V
# and the number of columns that moved
.lame_uv_sym_split_draw <- function(U, V, l_var = 1, hier = TRUE, w = 0.5) {
	still <- list(U = U, V = V, moved = 0L)
	R <- ncol(U)
	n <- nrow(U)
	if (R < 1L || n < 2L || !all(is.finite(U)) || !all(is.finite(V))) return(still)
	if (!is.finite(l_var) || l_var <= 0) return(still)
	moved <- 0L
	for (r in seq_len(R)) {
		a <- sum(U[, r]^2)
		if (!(a > 0)) next
		l2 <- sum(V[, r]^2) / a
		prior_u <- if (hier) function(u) -0.5 * (n + 2) * log1p(exp(2 * u) * a) else
			function(u) -exp(2 * u) * a / (2 * n)
		lg <- function(u) prior_u(u) - exp(-4 * u) * l2 / (2 * l_var) + (n - 2) * u
		u <- .lame_bso_slice(lg, w)$u
		if (!is.finite(u) || u == 0) next
		U[, r] <- U[, r] * exp(u)
		V[, r] <- V[, r] * exp(-u)
		moved <- moved + 1L
	}
	list(U = U, V = V, moved = moved)
}
