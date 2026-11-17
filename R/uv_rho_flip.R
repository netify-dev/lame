# exact mirror move of the persistence of dynamic latent paths (?lame, prior$uv_rho_flip_move).
# F: rho -> -rho, x_t -> (-1)^(t - 1) x_t for every actor of both sides is an involution with unit
# jacobian. it keeps every U_t V_t', U_t G V_t' and U_t L U_t', so the likelihood, the latent array
# and the dyadic correlation; each transition w_t ||x_t - k_t rho x_t-1||^2 (k = 0 at a snap reset,
# w the t weights), since x_t and rho x_t-1 change sign together; the stationary first period and the
# reset scale, which see rho only through rho^2; the resets, the weights and the box max_abs. only the
# normal prior of rho, truncated to a symmetric interval, changes, so the metropolis ratio is
#     pi(-rho) / pi(rho) = exp(-2 rho m / s^2)
# the uniform is recycled from the persistence draw just before it, so the move draws no random
# numbers: the slice height over the density at the new rho, or the metropolis uniform rescaled within
# its accept or reject interval, is uniform and independent of the new state. a draw that used no
# uniform skips the move; the chance of that depends on the state only through |rho| and the
# proposal scale, which F keeps, so the skip leaves the mirror ratio exact.
# the data see |rho| only, so lame() stores every draw in the gauge rho >= 0: a state at negative
# rho is saved as its mirror F(state), which every reported quantity takes as the same draw, and the
# sampler goes on from the state itself

# a metropolis uniform rescaled within its accept or reject interval, a fresh uniform independent of
# the decision and so of the new state (`la` the log acceptance ratio)
.lame_recycle_u <- function(u, la, accepted) {
	a <- exp(min(la, 0))
	if (accepted) u / a else (u - a) / (1 - a)
}

# `U_cube`, `V_cube` the paths (actors x rank x N; V_cube may be NULL), `rho` the persistence,
# `prior_mean`, `prior_sd` its normal prior, `u` a uniform independent of the state (NULL or NA
# skips the move); returns the state, whether the move was tried and whether it moved
.lame_uv_rho_flip <- function(U_cube, V_cube, rho, prior_mean, prior_sd, u) {
	out <- list(U = U_cube, V = V_cube, rho = rho, visited = FALSE, moved = FALSE)
	N <- dim(U_cube)[3L]
	if (length(u) != 1L || !is.finite(u) || length(rho) != 1L || !is.finite(rho) ||
			is.null(N) || N < 2L) return(out)
	out$visited <- TRUE
	if (log(u) >= -2 * rho * prior_mean / prior_sd^2) return(out)
	out$U <- .lame_uv_mirror(U_cube)
	if (!is.null(V_cube)) out[["V"]] <- .lame_uv_mirror(V_cube)
	out$rho <- -rho
	out$moved <- TRUE
	out
}

# the paths of the mirror state, x_t -> (-1)^(t - 1) x_t (`X` actors x rank x N)
.lame_uv_mirror <- function(X) sweep(X, 3L, rep_len(c(1, -1), dim(X)[3L]), "*")

# period signs of each coordinate path of a symmetric fit (?lame, prior$uv_sym_sign_move). with L
# diagonal, x_tr -> -x_tr for one coordinate r at one period t keeps U_t L U_t', so the likelihood,
# and the unit-variance AR(1) prior of the paths sees the signs only through the adjacent products
# k_tr = rho / (1 - rho^2) <x_t-1,r, x_tr>: the relative signs d_tr = c_t-1,r c_tr are independent
# with P(d_tr = 1) = plogis(2 k_tr). the draw is exact and leaves every coordinate bound.
# `U_cube` the paths (actors x rank x N), `rho` the persistence, `keep` (actors x N) 1 where a period
# drifts from the last and 0 at a snap reset, whose position is symmetric about zero; returns the paths
# and the number of coordinate periods whose sign changed
.lame_uv_sym_sign_draw <- function(U_cube, rho, keep = NULL) {
	d <- dim(U_cube)
	N <- d[3L]
	out <- list(U = U_cube, moved = 0L)
	if (is.null(N) || N < 2L || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 ||
			!all(is.finite(U_cube))) return(out)
	w <- rho / (1 - rho * rho)
	R <- d[2L]
	for (r in seq_len(R)) {
		x <- matrix(U_cube[, r, ], d[1L], N)
		k <- if (is.null(keep)) w * colSums(x[, -N, drop = FALSE] * x[, -1L, drop = FALSE])
		else w * colSums(keep[, -1L, drop = FALSE] * x[, -N, drop = FALSE] * x[, -1L, drop = FALSE])
		signs <- cumprod(c(1, ifelse(stats::runif(N - 1L) < stats::plogis(2 * k), 1, -1)))
		if (all(signs == 1)) next
		out$moved <- out$moved + sum(signs < 0)
		out$U[, r, ] <- sweep(x, 2L, signs, "*")
	}
	out
}
