# exact draw of the scale split of dynamic latent paths: U_t -> c U_t, V_t -> V_t / c leaves
# U_t G V_t', for all periods together or one period alone, so only the path prior sees c.
# the paths have density exp(-(Q_U + Q_V) / (2 sigma^2)) with
#     Q_X = sum_i [ w_i1 (1 - rho^2) ||x_i1||^2 + sum_{t > 1} w_it ||x_it - k_it rho x_i,t-1||^2 ]
# (w = k = 1 for ar1, w the scale-mixture weights for t, k = 0 at a snap reset), and u = log c is
# drawn by a slice step from u = 0, cut where a coordinate would leave max_abs

# sums over the actors of one side: sq[, t] = ||x_it||^2, cr[, t] = <x_it, x_i,t-1>
.lame_uvsd_sums <- function(X) {
	N <- dim(X)[3L]
	sq <- apply(X * X, c(1L, 3L), sum)
	cr <- matrix(0, dim(X)[1L], N)
	if (N > 1L) cr[, -1L] <- apply(X[, , -1L, drop = FALSE] * X[, , -N, drop = FALSE], c(1L, 3L), sum)
	list(sq = matrix(sq, dim(X)[1L], N), cr = cr)
}

# weighted squared norm `a` and cross product `b` of period t of one side
.lame_uvsd_period <- function(sm, w, k, rho, p1, t) {
	N <- ncol(sm$sq)
	own <- if (t == 1L) p1 * w[, 1L] else w[, t]
	a <- sum(own * sm$sq[, t])
	b <- 0
	if (t > 1L) b <- b + rho * sum(w[, t] * k[, t] * sm$cr[, t])
	if (t < N) {
		a <- a + rho * rho * sum(w[, t + 1L] * k[, t + 1L] * k[, t + 1L] * sm$sq[, t])
		b <- b + rho * sum(w[, t + 1L] * k[, t + 1L] * sm$cr[, t + 1L])
	}
	c(a = a, b = b)
}

# the quadratic form Q of one side
.lame_uvsd_q <- function(sm, w, k, rho, p1) {
	N <- ncol(sm$sq)
	q <- p1 * sum(w[, 1L] * sm$sq[, 1L])
	for (t in seq_len(N)[-1L]) {
		q <- q + sum(w[, t] * (sm$sq[, t] - 2 * rho * k[, t] * sm$cr[, t] +
			rho * rho * k[, t] * k[, t] * sm$sq[, t - 1L]))
	}
	q
}

# paths `U_cube` (nA x RA x N) and `V_cube` (nB x RB x N), precision weights `w_u`, `w_v` and
# persistence indicators `k_u`, `k_v` (actors x N); returns the paths and the number of steps that moved
.lame_uv_split_dynamic_draw <- function(
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
	if (length(dU) != 3L || length(dV) != 3L || dU[3L] != dV[3L]) return(still)
	if (any(dU == 0L) || any(dV == 0L)) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(V_cube))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	N <- dU[3L]
	one <- function(x, n) if (is.null(x)) matrix(1, n, N) else matrix(as.numeric(x), n, N)
	w_u <- one(w_u, dU[1L]); w_v <- one(w_v, dV[1L])
	k_u <- one(k_u, dU[1L]); k_v <- one(k_v, dV[1L])
	if (!all(is.finite(c(w_u, w_v, k_u, k_v))) || any(c(w_u, w_v) < 0)) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	s2 <- sigma * sigma
	p1 <- 1 - rho * rho
	lam <- dU[1L] * dU[2L] - dV[1L] * dV[2L]
	moved <- 0L
	# bounds on u keeping |c x| <= max_abs on both sides
	lims <- function(mu, mv) {
		c(if (mv > 0) log(mv / max_abs) else -Inf, if (mu > 0) log(max_abs / mu) else Inf)
	}
	su <- .lame_uvsd_sums(U_cube)
	sv <- .lame_uvsd_sums(V_cube)
	# joint step over every period
	qu <- .lame_uvsd_q(su, w_u, k_u, rho, p1)
	qv <- .lame_uvsd_q(sv, w_v, k_v, rho, p1)
	if (is.finite(qu) && is.finite(qv) && qu > 0 && qv > 0) {
		lg <- function(u) -(exp(2 * u) * qu + exp(-2 * u) * qv) / (2 * s2) + N * lam * u
		lm <- lims(max(abs(U_cube)), max(abs(V_cube)))
		if (lm[1L] <= 0 && lm[2L] >= 0) {
			u <- .lame_nc_slice(lg, 0, lm[1L], lm[2L], w = width)$x
			if (is.finite(u) && u != 0) {
				U_cube <- U_cube * exp(u); V_cube <- V_cube * exp(-u)
				su <- lapply(su, function(m) m * exp(2 * u)); sv <- lapply(sv, function(m) m * exp(-2 * u))
				moved <- moved + 1L
			}
		}
	}
	# one period at a time, its neighbours held
	for (t in seq_len(N)) {
		pu <- .lame_uvsd_period(su, w_u, k_u, rho, p1, t)
		pv <- .lame_uvsd_period(sv, w_v, k_v, rho, p1, t)
		if (!all(is.finite(c(pu, pv))) || pu[["a"]] <= 0 || pv[["a"]] <= 0) next
		au <- pu[["a"]]; bu <- pu[["b"]]; av <- pv[["a"]]; bv <- pv[["b"]]
		lg <- function(u) {
			-(au * exp(2 * u) - 2 * bu * exp(u) + av * exp(-2 * u) - 2 * bv * exp(-u)) / (2 * s2) + lam * u
		}
		lm <- lims(max(abs(U_cube[, , t])), max(abs(V_cube[, , t])))
		if (lm[1L] > 0 || lm[2L] < 0) next
		u <- .lame_nc_slice(lg, 0, lm[1L], lm[2L], w = width)$x
		if (!is.finite(u) || u == 0) next
		cu <- exp(u); cv <- exp(-u)
		U_cube[, , t] <- U_cube[, , t] * cu
		V_cube[, , t] <- V_cube[, , t] * cv
		# rescale the sums that touch period t
		su$sq[, t] <- su$sq[, t] * cu * cu; sv$sq[, t] <- sv$sq[, t] * cv * cv
		su$cr[, t] <- su$cr[, t] * cu; sv$cr[, t] <- sv$cr[, t] * cv
		if (t < N) { su$cr[, t + 1L] <- su$cr[, t + 1L] * cu; sv$cr[, t + 1L] <- sv$cr[, t + 1L] * cv }
		moved <- moved + 1L
	}
	list(U = U_cube, V = V_cube, moved = moved)
}

# precision weights w and persistence indicators k of the split and rotation moves: the
# scale-mixture weights for t (the first period's included), k = 0 and w = snap_w at a snap
# reset, NULL (all ones) for ar1
.lame_uv_dyn_kind_weights <- function(
	kind,
	lambda_u = NULL,
	lambda_v = NULL,
	delta_u = NULL,
	delta_v = NULL,
	snap_w = 1
) {
	out <- list(w_u = NULL, w_v = NULL, k_u = NULL, k_v = NULL)
	if (identical(kind, "t")) {
		out$w_u <- lambda_u; out$w_v <- lambda_v
	} else if (identical(kind, "snap")) {
		k_u <- 1 * (delta_u < 0.5); k_v <- 1 * (delta_v < 0.5)
		k_u[, 1L] <- 1; k_v[, 1L] <- 1
		out$k_u <- k_u; out$k_v <- k_v
		out$w_u <- ifelse(k_u > 0.5, 1, snap_w); out$w_v <- ifelse(k_v > 0.5, 1, snap_w)
	}
	out
}
