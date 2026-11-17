# exact draw of the orientation of directed dynamic latent paths (?lame, prior$uv_rotation_move).
# x_is -> x_is Q for s = t, ..., N and every actor, Q orthogonal, leaves each U_s V_s' and changes
# the path prior (w, k as in R/uv_split_dynamic.R) only through the transition into t, so Q has the
# matrix Fisher density exp(tr(Q' M_t)) on O(R) with
#     M_t = (rho / sigma^2) sum over sides and actors w_it k_it x_it' x_i,t-1
# and is drawn one coordinate plane at a time from its O(2) law (.lame_uv_o2), keeping the state
# when a tail coordinate would pass max_abs. under kind snap each plane then tries its mirror
# image F P' by a metropolis step on f_t, the transition into t with its resets summed out,
#     m(y | x) = (1 - pi) N(y; rho x, sigma^2 I) + pi N(y; 0, kappa^2 I),
# and the resets of that transition are drawn again

# M_t of one side: the cross products of period t with the weighted period before
.lame_uvrd_m <- function(X, w, k, t) {
	n <- dim(X)[1L]
	R <- dim(X)[2L]
	crossprod(matrix(X[, , t], n, R), (w[, t] * k[, t]) * matrix(X[, , t - 1L], n, R))
}

# log densities of the drift (a) and the reset (b) of the rows y given rows x
.lame_uvrd_ab <- function(y, x, rho, s2, k2, lp, lq) {
	list(
		a = lq - rowSums((y - rho * x)^2) / (2 * s2) - 0.5 * ncol(y) * log(s2),
		b = lp - rowSums(y * y) / (2 * k2) - 0.5 * ncol(y) * log(k2)
	)
}

# sum over the rows of log m(y | x), the transition with its reset summed out (lp, lq = log pi, log(1 - pi))
.lame_uvrd_mix <- function(y, x, rho, s2, k2, lp, lq) {
	ab <- .lame_uvrd_ab(y, x, rho, s2, k2, lp, lq)
	m <- pmax(ab$a, ab$b)
	sum(m + log(exp(ab$a - m) + exp(ab$b - m)))
}

# log f_t of one side, the transition into period t holding `xt`
.lame_uvrd_snap_ll <- function(X, xt, t, rho, s2, k2, lp, lq) {
	.lame_uvrd_mix(xt, matrix(X[, , t - 1L], dim(X)[1L], dim(X)[2L]), rho, s2, k2, lp, lq)
}

# a draw of the resets of transition s of one side, given its paths
.lame_uvrd_reset <- function(X, s, rho, s2, k2, lp, lq) {
	n <- dim(X)[1L]
	R <- dim(X)[2L]
	ab <- .lame_uvrd_ab(matrix(X[, , s], n, R), matrix(X[, , s - 1L], n, R), rho, s2, k2, lp, lq)
	1 * (stats::runif(n) < stats::plogis(ab$b - ab$a))
}

# step at t: x_s -> x_s E for s = t, ..., N, plane by plane, from the matrix Fisher law, or with
# `sn` (the reset law; `st` then also holds the resets du, dv) the mirror jump and the resets of t
.lame_uvrd_step <- function(st, t, rho, s2, w_u, w_v, k_u, k_v, pairs, max_abs, sn = NULL) {
	N <- dim(st$U)[3L]
	nA <- dim(st$U)[1L]
	nB <- dim(st$V)[1L]
	R <- dim(st$U)[2L]
	M <- rho / s2 * (.lame_uvrd_m(st$U, w_u, k_u, t) + .lame_uvrd_m(st$V, w_v, k_v, t))
	if (!all(is.finite(M))) return(st)
	ut <- matrix(st$U[, , t], nA, R)
	vt <- matrix(st$V[, , t], nB, R)
	later <- seq_len(N)[-seq_len(t)]
	if (!is.null(sn)) ll <- .lame_uvrd_snap_ll(st$U, ut, t, rho, s2, sn$k2, sn$lp, sn$lq) +
		.lame_uvrd_snap_ll(st$V, vt, t, rho, s2, sn$k2, sn$lp, sn$lq)
	for (j in seq_len(ncol(pairs))) {
		pq <- pairs[, j]
		B <- M[pq, pq, drop = FALSE]
		E2 <- if (is.null(sn)) .lame_uv_o2(B) else .lame_uvrd_mirror(B)
		if (identical(E2, diag(2))) next
		u_new <- ut[, pq, drop = FALSE] %*% E2
		v_new <- vt[, pq, drop = FALSE] %*% E2
		top <- max(abs(u_new), abs(v_new))
		if (is.finite(max_abs)) for (s in later) {
			top <- max(
				top,
				abs(matrix(st$U[, pq, s], nA, 2L) %*% E2),
				abs(matrix(st$V[, pq, s], nB, 2L) %*% E2)
			)
		}
		if (top > max_abs) next
		if (!is.null(sn)) {
			ut_p <- ut; vt_p <- vt
			ut_p[, pq] <- u_new; vt_p[, pq] <- v_new
			ll_p <- .lame_uvrd_snap_ll(st$U, ut_p, t, rho, s2, sn$k2, sn$lp, sn$lq) +
				.lame_uvrd_snap_ll(st$V, vt_p, t, rho, s2, sn$k2, sn$lp, sn$lq)
			if (!is.finite(ll_p) || log(stats::runif(1)) >= ll_p - ll) next
			ll <- ll_p
		}
		ut[, pq] <- u_new
		vt[, pq] <- v_new
		# turn the rest of the tail with period t
		for (s in later) {
			st$U[, pq, s] <- matrix(st$U[, pq, s], nA, 2L) %*% E2
			st$V[, pq, s] <- matrix(st$V[, pq, s], nB, 2L) %*% E2
		}
		# rotate the rows of M with the coordinates of period t
		M[pq, ] <- crossprod(E2, M[pq, , drop = FALSE])
		st$hit <- TRUE
	}
	st$U[, , t] <- ut
	st$V[, , t] <- vt
	if (!is.null(sn)) {
		# draw the resets of transition t given the paths
		st$du[, t] <- .lame_uvrd_reset(st$U, t, rho, s2, sn$k2, sn$lp, sn$lq)
		st$dv[, t] <- .lame_uvrd_reset(st$V, t, rho, s2, sn$k2, sn$lp, sn$lq)
	}
	st
}

# the mirror jump F P' of a 2 x 2 block: P the rotation and F the reflection that
# maximise tr(E' B), with the conventions of .lame_uv_o2
.lame_uvrd_mirror <- function(B) {
	t1 <- atan2(B[2L, 1L] - B[1L, 2L], B[1L, 1L] + B[2L, 2L])
	t2 <- atan2(B[2L, 1L] + B[1L, 2L], B[1L, 1L] - B[2L, 2L])
	P <- matrix(c(cos(t1), sin(t1), -sin(t1), cos(t1)), 2L, 2L)
	refl <- matrix(c(cos(t2), sin(t2), sin(t2), -cos(t2)), 2L, 2L)
	tcrossprod(refl, P)
}

# arguments as for .lame_uv_split_dynamic_draw; kind snap also passes
# `snap = list(pi, kappa, delta_u, delta_v)`. returns the paths, the steps that moved and the resets
.lame_uv_rotation_dynamic_draw <- function(
	U_cube,
	V_cube,
	rho,
	sigma,
	w_u = NULL,
	w_v = NULL,
	k_u = NULL,
	k_v = NULL,
	max_abs = Inf,
	snap = NULL
) {
	still <- list(U = U_cube, V = V_cube, moved = 0L)
	if (!is.null(snap)) still[c("delta_u", "delta_v")] <- list(snap$delta_u, snap$delta_v)
	dU <- dim(U_cube)
	dV <- dim(V_cube)
	if (length(dU) != 3L || length(dV) != 3L || dU[3L] != dV[3L] || dU[2L] != dV[2L]) return(still)
	if (dU[2L] < 2L || any(dU == 0L) || any(dV == 0L)) return(still)
	if (!all(is.finite(U_cube)) || !all(is.finite(V_cube))) return(still)
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma) || sigma <= 0) return(still)
	N <- dU[3L]
	R <- dU[2L]
	if (N < 2L || rho == 0) return(still)
	one <- function(x, n) if (is.null(x)) matrix(1, n, N) else matrix(as.numeric(x), n, N)
	w_u <- one(w_u, dU[1L]); w_v <- one(w_v, dV[1L])
	k_u <- one(k_u, dU[1L]); k_v <- one(k_v, dV[1L])
	if (!all(is.finite(c(w_u, w_v, k_u, k_v))) || any(c(w_u, w_v) < 0)) return(still)
	if (is.null(max_abs) || !is.finite(max_abs) || max_abs <= 0) max_abs <- Inf
	s2 <- sigma * sigma
	st <- list(U = U_cube, V = V_cube)
	sn <- NULL
	if (!is.null(snap)) {
		pi_s <- as.numeric(snap$pi)
		k2 <- as.numeric(snap$kappa)^2
		if (length(pi_s) != 1L || !is.finite(pi_s) || pi_s < 0 || pi_s > 1) return(still)
		if (length(k2) != 1L || !is.finite(k2) || k2 <= 0) return(still)
		if (!identical(dim(snap$delta_u), dU[c(1L, 3L)]) || !identical(dim(snap$delta_v), dV[c(1L, 3L)])) return(still)
		sn <- list(k2 = k2, lp = log(pi_s), lq = log1p(-pi_s))
		st$du <- matrix(as.numeric(snap$delta_u > 0.5), dU[1L], N)
		st$dv <- matrix(as.numeric(snap$delta_v > 0.5), dV[1L], N)
		st$du[, 1L] <- 0; st$dv[, 1L] <- 0
		# the mirror jump reads the block with every transition a drift
		one_u <- matrix(1, dU[1L], N)
		one_v <- matrix(1, dV[1L], N)
		snap_w <- s2 / k2
	}
	pairs <- utils::combn(R, 2L)
	moved <- 0L
	for (t in seq_len(N)[-1L]) {
		st$hit <- FALSE
		if (is.null(sn)) {
			st <- .lame_uvrd_step(st, t, rho, s2, w_u, w_v, k_u, k_v, pairs, max_abs)
		} else {
			# step given the resets, then the mirror jump with the resets of transition t summed out
			k_u <- 1 - st$du; k_v <- 1 - st$dv
			w_u <- ifelse(k_u > 0.5, 1, snap_w); w_v <- ifelse(k_v > 0.5, 1, snap_w)
			st <- .lame_uvrd_step(st, t, rho, s2, w_u, w_v, k_u, k_v, pairs, max_abs)
			st <- .lame_uvrd_step(st, t, rho, s2, one_u, one_v, one_u, one_v, pairs, max_abs, sn)
		}
		moved <- moved + st$hit
	}
	out <- list(U = st$U, V = st$V, moved = as.integer(moved))
	if (!is.null(sn)) out[c("delta_u", "delta_v")] <- list(st$du, st$dv)
	out
}
