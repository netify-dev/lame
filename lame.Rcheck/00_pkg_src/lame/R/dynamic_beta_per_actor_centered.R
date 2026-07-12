# centered per-actor dynamic-beta sweep. samples each subject's length-t
# slope path under an ar(1) prior using pairwise contrasts that preserve the
# per-period sum-to-zero constraint exactly.

# per-period observation sufficient statistics for one actor block: precision
# h and cross-product h per actor (row or column side).
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

# ar(1) path precision matrix for a length-t_per latent path.
.ar1_path_precision <- function(T_per, rho, sigma2) {
	if (T_per < 1L) return(matrix(0, 0L, 0L))
	if (T_per == 1L) return(matrix(1 / sigma2, 1L, 1L))
	rho <- max(min(rho, 0.99999), -0.99999)
	sigma2 <- max(sigma2, 1e-12)
	K <- matrix(0, T_per, T_per)
	for (t in seq_len(T_per - 1L)) {
		K[t, t] <- K[t, t] + rho^2 / sigma2
		K[t + 1, t + 1] <- K[t + 1, t + 1] + 1 / sigma2
		K[t, t + 1] <- K[t, t + 1] - rho / sigma2
		K[t + 1, t] <- K[t + 1, t] - rho / sigma2
	}
	K[1, 1] <- K[1, 1] + (1 - rho^2) / sigma2
	K
}

# sample the zero-sum contrast for one actor pair and return the shift.
.sample_pair_contrast <- function(theta, i, ell, H_mat, h_mat, Q_AR1) {
	T_per <- ncol(theta)
	H_data <- H_mat[i, ] + H_mat[ell, ]
	h_data <- (h_mat[i, ] - H_mat[i, ] * theta[i, ]) -
		(h_mat[ell, ] - H_mat[ell, ] * theta[ell, ])
	Q_delta <- 2 * Q_AR1 + diag(H_data, T_per)
	h_delta <- h_data - as.numeric(Q_AR1 %*% theta[i, ]) +
		as.numeric(Q_AR1 %*% theta[ell, ])
	L <- tryCatch(chol(Q_delta), error = function(e) NULL)
	if (is.null(L)) {
		L <- chol(Q_delta + diag(1e-06 * max(diag(Q_delta), 1), T_per))
	}
	m <- backsolve(L, forwardsolve(t(L), h_delta))
	z <- stats::rnorm(T_per)
	delta <- as.numeric(m + backsolve(L, z))
	delta
}

# one randomised sweep of zero-sum pairwise updates over the actors.
.sweep_per_actor <- function(theta, H_mat, h_mat, Q_AR1) {
	n_actors <- nrow(theta)
	if (n_actors < 2L) return(theta)
	idx <- sample.int(n_actors)
	n_pairs <- floor(n_actors / 2L)
	for (p in seq_len(n_pairs)) {
		i <- idx[2L * p - 1L]
		ell <- idx[2L * p]
		delta <- .sample_pair_contrast(theta, i, ell, H_mat, h_mat, Q_AR1)
		theta[i, ] <- theta[i, ] + delta
		theta[ell, ] <- theta[ell, ] - delta
	}
	theta
}

# inverse-gamma draw of the per-actor ar(1) innovation variance.
.sample_sigma_actor <- function(theta, rho_actor, prior_shape = 2, prior_scale = 1) {
	n_actors <- nrow(theta)
	T_per <- ncol(theta)
	rho_actor <- max(min(rho_actor, 0.99999), -0.99999)
	if (T_per < 2L) {
		quad <- sum(theta^2)
		shape_post <- prior_shape + n_actors / 2
		scale_post <- prior_scale + 0.5 * quad
	} else {
		init_quad <- (1 - rho_actor^2) * sum(theta[, 1L]^2)
		diffs <- theta[, -1L, drop = FALSE] - rho_actor * theta[, -T_per, drop = FALSE]
		trans_quad <- sum(diffs^2)
		shape_post <- prior_shape + n_actors * T_per / 2
		scale_post <- prior_scale + 0.5 * (init_quad + trans_quad)
	}
	1 / stats::rgamma(1L, shape = shape_post, rate = scale_post)
}

# random-walk metropolis update of the per-actor ar(1) persistence.
.sample_rho_actor <- function(theta, sigma2, rho_current, prior_a = 2, prior_b = 2,
                              sd_prop = 0.05) {
	T_per <- ncol(theta)
	if (T_per < 2L) return(rho_current)
	rho_prop <- rho_current + stats::rnorm(1L, 0, sd_prop)
	if (rho_prop > 0.99) rho_prop <- 0.99 - (rho_prop - 0.99)
	if (rho_prop < -0.99) rho_prop <- -0.99 - (rho_prop + 0.99)
	rho_prop <- max(min(rho_prop, 0.99), -0.99)
	loglik <- function(rho) {
		diffs <- theta[, -1L, drop = FALSE] - rho * theta[, -T_per, drop = FALSE]
		-0.5 * sum(diffs^2) / sigma2
	}
	logprior <- function(rho) {
		stats::dbeta((rho + 1) / 2, prior_a, prior_b, log = TRUE)
	}
	log_ratio <- loglik(rho_prop) + logprior(rho_prop) -
		loglik(rho_current) - logprior(rho_current)
	if (log(stats::runif(1L)) < log_ratio) rho_prop else rho_current
}
