# centred per-actor dynamic-beta sweep: each actor's slope path under an AR(1)
# prior, moved by pairwise contrasts that keep the per-period sum to zero

# per-period sufficient statistics for one actor side: precision H and
# cross-product h per actor
.per_actor_suffstats_period <- function(R_t, X_t, kind, s2) {
	mask <- is.finite(R_t) & is.finite(X_t)
	if (kind == "row") {
		Hs <- rowSums((X_t^2) * mask, na.rm = TRUE) / s2
		hs <- rowSums(X_t * R_t * mask, na.rm = TRUE) / s2
	} else {
		Hs <- colSums((X_t^2) * mask, na.rm = TRUE) / s2
		hs <- colSums(X_t * R_t * mask, na.rm = TRUE) / s2
	}
	list(H = Hs, h = hs)
}

# AR(1) precision matrix of a stationary path of length T_per
.ar1_path_precision <- function(T_per, rho, sigma2) {
	if (T_per < 1L) return(matrix(0, 0L, 0L))
	if (!is.finite(rho) || abs(rho) >= 1 || !is.finite(sigma2) || sigma2 <= 0)
		stop("Invalid actor-specific AR parameters.")
	K <- matrix(0, T_per, T_per)
	K[1L, 1L] <- (1 - rho^2) / sigma2
	if (T_per > 1L) for (t in seq_len(T_per - 1L)) {
		K[t, t] <- K[t, t] + rho^2 / sigma2
		K[t + 1L, t + 1L] <- K[t + 1L, t + 1L] + 1 / sigma2
		K[t, t + 1L] <- K[t + 1L, t] <- -rho / sigma2
	}
	K
}

# draw the zero-sum contrast of one actor pair and apply the shift
.sample_pair_contrast <- function(
	theta,
	i,
	ell,
	H_mat,
	h_mat,
	Q_AR1,
	precision_list = NULL
) {
	T_per <- ncol(theta)
	if (is.null(precision_list)) {
		H_data <- H_mat[i, ] + H_mat[ell, ]
		h_data <- (h_mat[i, ] - H_mat[i, ] * theta[i, ]) -
			(h_mat[ell, ] - H_mat[ell, ] * theta[ell, ])
	} else {
		H_data <- vapply(precision_list, function(P) P[i, i] + P[ell, ell] - 2 * P[i, ell], numeric(1))
		h_data <- vapply(seq_len(T_per), function(t) {
			P <- precision_list[[t]]
			h_mat[i, t] - h_mat[ell, t] - sum((P[i, ] - P[ell, ]) * theta[, t])
		}, numeric(1))
	}
	Q_delta <- 2 * Q_AR1 + diag(H_data, T_per)
	h_delta <- h_data - as.numeric(Q_AR1 %*% theta[i, ]) +
		as.numeric(Q_AR1 %*% theta[ell, ])
	L <- chol(Q_delta)
	m <- backsolve(L, forwardsolve(t(L), h_delta))
	as.numeric(m + backsolve(L, stats::rnorm(T_per)))
}

# one randomised sweep of zero-sum pairwise updates over the actors
.sweep_per_actor <- function(theta, H_mat, h_mat, Q_AR1, precision_list = NULL) {
	n_actors <- nrow(theta)
	if (n_actors < 2L) return(theta)
	idx <- sample.int(n_actors)
	n_pairs <- floor(n_actors / 2L)
	for (p in seq_len(n_pairs)) {
		i <- idx[2L * p - 1L]
		ell <- idx[2L * p]
		delta <- .sample_pair_contrast(theta, i, ell, H_mat, h_mat, Q_AR1, precision_list)
		theta[i, ] <- theta[i, ] + delta
		theta[ell, ] <- theta[ell, ] - delta
	}
	theta
}

# inverse-gamma draw of the per-actor AR(1) innovation variance
.sample_sigma_actor <- function(theta, rho_actor, prior_shape = 2, prior_scale = 1) {
	dimension <- max(0L, nrow(theta) - 1L)
	T_per <- ncol(theta)
	Q_unit <- .ar1_path_precision(T_per, rho_actor, 1)
	quad <- sum(theta * (theta %*% Q_unit))
	variance <- 1 / stats::rgamma(
		1L,
		shape = prior_shape + dimension * T_per / 2,
		rate = prior_scale + .5 * quad
	)
	if (!is.finite(variance) || variance <= 0) stop("Invalid actor-specific innovation variance draw.")
	variance
}

# slice sampler for the persistence on the whole prior support (-0.99, 0.99),
# shrinking toward the current point; exact for non-log-concave conditionals
.sample_rho_actor <- function(
	theta,
	sigma2,
	rho_current,
	prior_a = 2,
	prior_b = 2
) {
	dimension <- max(0L, nrow(theta) - 1L)
	log_density <- function(rho) {
		Q <- .ar1_path_precision(ncol(theta), rho, sigma2)
		.5 * dimension * log1p(-rho^2) - .5 * sum(theta * (theta %*% Q)) +
			stats::dbeta((rho + 1) / 2, prior_a, prior_b, log = TRUE)
	}
	if (!is.finite(rho_current) || abs(rho_current) >= .99)
		stop("Actor-specific persistence is outside its prior support.")
	level <- log_density(rho_current) + log(stats::runif(1L))
	lower <- -.99; upper <- .99
	for (attempt in seq_len(10000L)) {
		proposal <- stats::runif(1L, lower, upper)
		if (log_density(proposal) >= level) return(proposal)
		if (proposal < rho_current) lower <- proposal else upper <- proposal
	}
	stop("Actor-specific persistence slice update failed to resolve its conditional.")
}
