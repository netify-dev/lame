# Gaussian sufficient statistics for actor-specific slopes; reciprocity
# couples the slopes on the two directions of each unordered dyad
.lame_per_actor_suffstats <- function(
	residual,
	covariate,
	kind,
	s2,
	rho = 0,
	bipartite = FALSE
) {
	if (bipartite || rho == 0) {
		z <- .per_actor_suffstats_period(residual, covariate, kind, s2)
		z$precision <- NULL
		return(z)
	}
	n <- nrow(residual)
	if (ncol(residual) != n || !is.finite(rho) || abs(rho) >= 1 || s2 <= 0)
		stop("Invalid actor-specific dyadic Gaussian likelihood.")
	P <- matrix(0, n, n); score <- numeric(n)
	for (i in seq_len(n)) {
		if (is.finite(residual[i, i]) && is.finite(covariate[i, i])) {
			P[i, i] <- P[i, i] + covariate[i, i]^2 / (s2 * (1 + rho))
			score[i] <- score[i] + covariate[i, i] * residual[i, i] / (s2 * (1 + rho))
		}
		if (i == n) next
		for (j in (i + 1L):n) {
			x <- c(covariate[i, j], covariate[j, i])
			r <- c(residual[i, j], residual[j, i])
			actors <- if (identical(kind, "row")) c(i, j) else c(j, i)
			ok <- is.finite(x) & is.finite(r)
			if (all(ok)) {
				W <- matrix(c(1, -rho, -rho, 1), 2) / (s2 * (1 - rho^2))
				P[actors, actors] <- P[actors, actors] + W * tcrossprod(x)
				score[actors] <- score[actors] + x * drop(W %*% r)
			} else if (any(ok)) {
				k <- actors[ok]; P[k, k] <- P[k, k] + x[ok]^2 / s2
				score[k] <- score[k] + x[ok] * r[ok] / s2
			}
		}
	}
	list(H = diag(P), h = score, precision = P)
}

.lame_per_actor_offset <- function(theta, covariates, kind) {
	out <- array(0, c(dim(covariates[[1L]]), length(covariates)))
	for (t in seq_along(covariates)) {
		out[, , t] <- sweep(
			covariates[[t]],
			if (identical(kind, "row")) 1L else 2L,
			theta[, t],
			`*`
		)
	}
	out
}

# base predictor once all blocks have moved, before the actor-specific paths are drawn
.lame_per_actor_base_predictor <- function(
	Xlist,
	beta_path,
	a_path,
	b_path,
	U,
	V,
	G = NULL,
	rank_active = TRUE
) {
	nt <- length(Xlist); na <- nrow(a_path); nb <- nrow(b_path)
	out <- array(0, c(na, nb, nt))
	for (t in seq_len(nt)) {
		X <- Xlist[[t]]; p <- dim(X)[3L]
		eta <- outer(a_path[, t], b_path[, t], `+`)
		if (p > 0L) eta <- eta + matrix(matrix(X, na * nb, p) %*% beta_path[t, ], na, nb)
		if (rank_active) {
			u <- if (length(dim(U)) == 3L) matrix(U[, , t], na, dim(U)[2L]) else U
			v <- if (length(dim(V)) == 3L) matrix(V[, , t], nb, dim(V)[2L]) else V
			eta <- eta + if (is.null(G)) tcrossprod(u, v) else u %*% G %*% t(v)
		}
		out[, , t] <- eta
	}
	out
}
