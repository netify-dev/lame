# draw the actor effects alone (no coefficients), keeping structural zeros
.lame_static_actor_draw <- function(
	residual,
	Sab,
	rho,
	s2,
	symmetric = FALSE,
	rvar = TRUE,
	cvar = TRUE,
	center = FALSE
) {
	n <- dim(residual)[1L]; nt <- dim(residual)[3L]
	S <- .lame_sab_gls(Sab, symmetric)
	if (!isTRUE(rvar)) S[1L, ] <- S[, 1L] <- 0
	if (!isTRUE(cvar)) S[2L, ] <- S[, 2L] <- 0
	.lame_replicated_gaussian_draw(
		residual,
		array(numeric(0), c(n * n, 0L, nt)),
		S,
		as.numeric(rho),
		s2,
		100,
		isTRUE(center)
	)
}

# joint draw of the static levels for a symmetric dynamic actor model, deviations held fixed
# each dyad involves two levels, so the levels cannot be drawn one row at a time
.lame_symmetric_level_draw <- function(
	residual,
	deviation,
	variance,
	s2,
	center = FALSE
) {
	n <- nrow(deviation)
	nt <- ncol(deviation)
	Q <- diag(1 / variance, n)
	h <- numeric(n)
	for (t in seq_len(nt)) {
		y <- residual[, , t] - outer(deviation[, t], deviation[, t], "+")
		keep <- upper.tri(y) & is.finite(y)
		weight <- (keep + t(keep)) / s2
		y[!keep] <- 0
		Q <- Q + weight
		diag(Q) <- diag(Q) + rowSums(weight)
		h <- h + rowSums(y + t(y)) / s2
	}
	L <- chol(Q)
	mu <- backsolve(L, forwardsolve(t(L), h))
	x <- drop(mu + backsolve(L, rnorm(n)))
	if (isTRUE(center)) {
		u <- backsolve(L, forwardsolve(t(L), rep(1, n)))
		x <- x - u * sum(x) / sum(u)
	}
	x
}
