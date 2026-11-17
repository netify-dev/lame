# the gauge of dynamic position paths beyond the scale split and the rotation: U_t -> U_t A,
# V_t -> V_t A^-T for every period together leaves every U_t V_t' unchanged (bipartite paths
# take G -> A^-1 G instead, symmetric ones L -> A^-1 L A^-T). the path prior sees A only
# through M = sum_it w_it e_it' e_it, e_it = x_it - k_it rho x_i,t-1 (the first period weighted
# by 1 - rho^2; w, k as in R/uv_split_dynamic.R), as Q = tr(A' M_U A) + tr(A^-1 M_V A^-T). the
# column scales A = I + (c - 1) e_k e_k' and the shears A = I + h e_k e_l' are one-parameter groups
# (Haar measure d log c and dh, Jacobian c^(N (n_U - n_V)) and 1), so each is drawn along its orbit
# (Liu and Sabatti 2000): log c by a slice step, h exactly as a Gaussian. every step keeps the
# coordinates inside max_abs

# the R x R moment M of one side, paths X (actors x R x N)
.lame_uvg_moment <- function(X, rho, w = NULL, k = NULL) {
	n <- dim(X)[1L]
	R <- dim(X)[2L]
	N <- dim(X)[3L]
	if (is.null(w)) w <- matrix(1, n, N)
	if (is.null(k)) k <- matrix(1, n, N)
	x1 <- matrix(X[, , 1L], n, R)
	M <- (1 - rho * rho) * crossprod(x1, w[, 1L] * x1)
	for (t in seq_len(N)[-1L]) {
		e <- matrix(X[, , t], n, R) - rho * k[, t] * matrix(X[, , t - 1L], n, R)
		M <- M + crossprod(e, w[, t] * e)
	}
	M
}

# bounds on log c keeping |c x| <= max_abs for the column `x` scaled up and `y` scaled down
.lame_uvg_bounds <- function(x, y, max_abs) {
	mx <- max(abs(x))
	my <- max(abs(y))
	c(if (my > 0) log(my / max_abs) else -Inf, if (mx > 0) log(max_abs / mx) else Inf)
}

# directed unipartite paths `U_cube`, `V_cube` (n x R x N), arguments as for
# .lame_uv_split_dynamic_draw; returns the paths and the number of steps that moved
.lame_uv_gauge_dynamic_draw <- function(
	U_cube,
	V_cube,
	rho,
	sigma,
	w_u = NULL,
	w_v = NULL,
	k_u = NULL,
	k_v = NULL,
	max_abs = Inf,
	width = 0.25
) {
	still <- list(U = U_cube, V = V_cube, moved = 0L)
	dU <- dim(U_cube)
	dV <- dim(V_cube)
	if (length(dU) != 3L || length(dV) != 3L || dU[2L] != dV[2L] || dU[3L] != dV[3L] || dU[2L] < 2L) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(V_cube))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	R <- dU[2L]
	s2 <- sigma * sigma
	lam <- dU[3L] * (dU[1L] - dV[1L])
	Mu <- .lame_uvg_moment(U_cube, rho, w_u, k_u)
	Mv <- .lame_uvg_moment(V_cube, rho, w_v, k_v)
	moved <- 0L
	# column scales: U_k -> c U_k, V_k -> V_k / c in every period
	for (k in seq_len(R)) {
		a <- Mu[k, k]
		b <- Mv[k, k]
		if (!(a > 0) || !(b > 0)) next
		lm <- .lame_uvg_bounds(U_cube[, k, ], V_cube[, k, ], max_abs)
		if (lm[1L] > 0 || lm[2L] < 0) next
		lg <- function(u) -(exp(2 * u) * a + exp(-2 * u) * b) / (2 * s2) + lam * u
		u <- .lame_nc_slice(lg, 0, lm[1L], lm[2L], w = width)$x
		if (!is.finite(u) || u == 0) next
		U_cube[, k, ] <- U_cube[, k, ] * exp(u)
		V_cube[, k, ] <- V_cube[, k, ] * exp(-u)
		Mu[k, ] <- Mu[k, ] * exp(u); Mu[, k] <- Mu[, k] * exp(u)
		Mv[k, ] <- Mv[k, ] * exp(-u); Mv[, k] <- Mv[, k] * exp(-u)
		moved <- moved + 1L
	}
	# shears: U_l -> U_l + h U_k, V_k -> V_k - h V_l in every period; Q is quadratic in h
	for (k in seq_len(R)) for (l in seq_len(R)[-k]) {
		prec <- Mu[k, k] + Mv[l, l]
		if (!(prec > 0)) next
		h <- stats::rnorm(1L, -(Mu[k, l] - Mv[k, l]) / prec, sigma / sqrt(prec))
		u_l <- U_cube[, l, ] + h * U_cube[, k, ]
		v_k <- V_cube[, k, ] - h * V_cube[, l, ]
		# a draw outside the bound keeps the state, exact for the bounded prior
		if (max(abs(u_l), abs(v_k)) > max_abs) next
		U_cube[, l, ] <- u_l
		V_cube[, k, ] <- v_k
		A <- diag(R)
		A[k, l] <- h
		B <- diag(R)
		B[l, k] <- -h
		Mu <- crossprod(A, Mu %*% A)
		Mv <- crossprod(B, Mv %*% B)
		moved <- moved + 1L
	}
	list(U = U_cube, V = V_cube, moved = moved)
}

# bipartite paths `U_cube` (nA x RA x N), `V_cube` (nB x RB x N) and the static G (RA x RB) with an
# N(0, g_var) prior cut at g_max: on each side the columns against the rows (columns) of G,
# U_t -> U_t A with G -> A^-1 G, V_t -> V_t C with G -> G C^-T. a column scale has Jacobian
# c^(nA N - RB), a shear h is Gaussian; `w_u`, `w_v`, `k_u`, `k_v` as for .lame_uv_gauge_dynamic_draw;
# returns the paths, G and the number of steps that moved
.lame_uv_gauge_bip_dynamic_draw <- function(
	U_cube,
	V_cube,
	G,
	rho,
	sigma,
	g_var = 1,
	g_max = Inf,
	max_abs = Inf,
	width = 0.25,
	w_u = NULL,
	w_v = NULL,
	k_u = NULL,
	k_v = NULL
) {
	still <- list(U = U_cube, V = V_cube, G = G, moved = 0L)
	dU <- dim(U_cube)
	dV <- dim(V_cube)
	if (length(dU) != 3L || length(dV) != 3L || dU[3L] != dV[3L]) return(still)
	if (!is.matrix(G) || !identical(dim(G), c(dU[2L], dV[2L]))) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(V_cube)) || !all(is.finite(G))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	if (!is.finite(g_var) || g_var <= 0) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	if (is.null(g_max) || !is.finite(g_max) || g_max <= 0) g_max <- Inf
	s2 <- sigma * sigma
	N <- dU[3L]
	# one side: paths X (actors x R x N) against the rows of Gm (R x R_other)
	side <- function(X, Gm, w, k) {
		R <- dim(X)[2L]
		M <- .lame_uvg_moment(X, rho, w, k)
		lam <- N * dim(X)[1L] - ncol(Gm)
		hits <- 0L
		for (k in seq_len(R)) {
			a <- M[k, k]
			b <- sum(Gm[k, ]^2)
			if (!(a > 0) || !(b > 0)) next
			# |c x| <= max_abs on the paths, |G_k / c| <= g_max on the interaction
			lm <- c(log(max(abs(Gm[k, ])) / g_max), log(max_abs / max(abs(X[, k, ]))))
			if (lm[1L] > 0 || lm[2L] < 0) next
			lg <- function(u) -exp(2 * u) * a / (2 * s2) - exp(-2 * u) * b / (2 * g_var) + lam * u
			u <- .lame_nc_slice(lg, 0, lm[1L], lm[2L], w = width)$x
			if (!is.finite(u) || u == 0) next
			X[, k, ] <- X[, k, ] * exp(u)
			Gm[k, ] <- Gm[k, ] * exp(-u)
			M[k, ] <- M[k, ] * exp(u); M[, k] <- M[, k] * exp(u)
			hits <- hits + 1L
		}
		for (k in seq_len(R)) for (l in seq_len(R)[-k]) {
			prec <- M[k, k] / s2 + sum(Gm[l, ]^2) / g_var
			if (!(prec > 0)) next
			m <- -(M[k, l] / s2 - sum(Gm[k, ] * Gm[l, ]) / g_var) / prec
			h <- stats::rnorm(1L, m, 1 / sqrt(prec))
			x_l <- X[, l, ] + h * X[, k, ]
			g_k <- Gm[k, ] - h * Gm[l, ]
			if (max(abs(x_l)) > max_abs || max(abs(g_k)) > g_max) next
			X[, l, ] <- x_l
			Gm[k, ] <- g_k
			A <- diag(R)
			A[k, l] <- h
			M <- crossprod(A, M %*% A)
			hits <- hits + 1L
		}
		list(X = X, G = Gm, hits = hits)
	}
	su <- side(U_cube, G, w_u, k_u)
	sv <- side(V_cube, t(su$G), w_v, k_v)
	list(U = su$X, V = sv$X, G = t(sv$G), moved = su$hits + sv$hits)
}

# symmetric paths `U_cube` (n x R x N) with weights `Ld` (prior N(0, l_var)), V_t = U_t diag(Ld):
# U_r -> c U_r, L_r -> L_r / c^2 in every period, Jacobian c^(n N - 2); `w`, `k` the transition weights;
# returns U, Ld and the number of columns that moved
.lame_uv_sym_split_dynamic_draw <- function(U_cube, Ld, rho, sigma, l_var = 1, max_abs = Inf, width = 0.25, w = NULL, k = NULL) {
	still <- list(U = U_cube, Ld = Ld, moved = 0L)
	dU <- dim(U_cube)
	if (length(dU) != 3L || length(Ld) != dU[2L] || dU[2L] < 1L) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(Ld))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	if (!is.finite(l_var) || l_var <= 0) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	s2 <- sigma * sigma
	lam <- dU[1L] * dU[3L] - 2
	M <- .lame_uvg_moment(U_cube, rho, w, k)
	moved <- 0L
	for (r in seq_len(dU[2L])) {
		a <- M[r, r]
		if (!(a > 0)) next
		l2 <- Ld[r]^2
		hi <- log(max_abs / max(abs(U_cube[, r, ])))
		if (hi < 0) next
		lg <- function(u) -exp(2 * u) * a / (2 * s2) - exp(-4 * u) * l2 / (2 * l_var) + lam * u
		u <- .lame_nc_slice(lg, 0, -Inf, hi, w = width)$x
		if (!is.finite(u) || u == 0) next
		U_cube[, r, ] <- U_cube[, r, ] * exp(u)
		Ld[r] <- Ld[r] * exp(-2 * u)
		moved <- moved + 1L
	}
	list(U = U_cube, Ld = Ld, moved = moved)
}

# the rest of the gauge of symmetric paths: U_t L U_t' is unchanged by U_t -> U_t A whenever A L A' = L. in the
# scaled coordinates a_r = |L_r|^(1/2) u_r a pair of columns turns by a hyperbolic rotation when the weights have
# opposite signs (a a' - b b' is kept) and by a circular one when they share a sign, so A = D H D^-1 with
# D = diag(|L_r|^(1/2), |L_q|^(1/2)) and det A = 1. the paths move together (x_s -> x_s A for every period s) or
# from a period t on (s = t, ..., N); the unit-variance AR(1) prior then sees A only through the transition into t
# and the innovations after it,
#     Q(A) = tr(A' M_t A) - 2 rho tr(A' C_t),   M_t = sum_i x_it' x_it + sum_{s > t} e_is' e_is,   C_t = sum_i x_it' x_i,t-1
# (t = 1: M_1 the moment of .lame_uvg_moment and C_1 = 0), and the angle, Haar measure on the one-parameter group,
# has density exp(-Q / (2 (1 - rho^2))). each angle is drawn by a slice step from the current one, the circular
# one over its period, and a step that would take a coordinate past max_abs keeps the state

# the 2 x 2 map of the pair with weights `l` at angle `x`
.lame_uv_sym_pair_map <- function(x, l) {
	d <- sqrt(abs(l))
	H <- if (l[1L] * l[2L] < 0) {
		matrix(c(cosh(x), sinh(x), sinh(x), cosh(x)), 2L, 2L)
	} else {
		matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2L, 2L)
	}
	H * outer(d, 1 / d)
}

# symmetric paths `U_cube` (n x R x N) with weights `Ld`; returns the paths and the number of steps that moved
.lame_uv_sym_pair_dynamic_draw <- function(U_cube, Ld, rho, max_abs = Inf, width = 0.25) {
	still <- list(U = U_cube, moved = 0L)
	dU <- dim(U_cube)
	if (length(dU) != 3L || dU[2L] < 2L || length(Ld) != dU[2L]) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(Ld))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	n <- dU[1L]
	N <- dU[3L]
	s2 <- 1 - rho * rho
	moved <- 0L
	pairs <- utils::combn(dU[2L], 2L)
	for (j in seq_len(ncol(pairs))) {
		pq <- pairs[, j]
		l <- Ld[pq]
		if (any(l == 0)) next
		hyper <- l[1L] * l[2L] < 0
		for (t in seq_len(N)) {
			later <- t:N
			X <- U_cube[, pq, later, drop = FALSE]
			if (t == 1L) {
				M <- .lame_uvg_moment(X, rho)
				C <- matrix(0, 2L, 2L)
			} else {
				prev <- matrix(U_cube[, pq, t - 1L], n, 2L)
				xt <- matrix(X[, , 1L], n, 2L)
				M <- crossprod(xt)
				for (s in seq_along(later)[-1L]) {
					e <- matrix(X[, , s], n, 2L) - rho * matrix(X[, , s - 1L], n, 2L)
					M <- M + crossprod(e)
				}
				C <- rho * crossprod(xt, prev)
			}
			if (!all(is.finite(M)) || !all(is.finite(C))) next
			Xm <- matrix(aperm(X, c(1L, 3L, 2L)), ncol = 2L)
			top <- max(abs(Xm))
			lg <- function(x) {
				A <- .lame_uv_sym_pair_map(x, l)
				# rows of Xm A stay within 2 top max|A|, so only a wide map needs the coordinates
				if (2 * top * max(abs(A)) > max_abs && max(abs(Xm %*% A)) > max_abs) return(-Inf)
				-(sum(A * (M %*% A)) - 2 * sum(A * C)) / (2 * s2)
			}
			x <- if (hyper) .lame_nc_slice(lg, 0, -Inf, Inf, w = width)$x else
				.lame_nc_slice(lg, 0, -Inf, Inf, w = 2 * pi, m = 1L)$x
			if (!is.finite(x) || x == 0) next
			A <- .lame_uv_sym_pair_map(x, l)
			for (s in later) U_cube[, pq, s] <- matrix(U_cube[, pq, s], n, 2L) %*% A
			moved <- moved + 1L
		}
	}
	list(U = U_cube, moved = moved)
}
