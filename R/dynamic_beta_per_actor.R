# actor-specific AR(1) paths keep their full temporal covariance, since the
# sum-to-zero constraint in every period binds jointly across time

# univariate FFBS draw of one actor's path with its full smoothed covariance
.actor_ffbs_path <- function(H, h, rho_actor, sigma_actor2) {
	T_per <- length(H)
	if (!T_per) return(list(theta = numeric(), V = numeric(), covariance = matrix(0, 0, 0)))
	if (length(h) != T_per || any(!is.finite(c(H, h))) || any(H < 0))
		stop("Invalid actor-specific Gaussian sufficient statistics.")
	.ar1_path_precision(T_per, rho_actor, sigma_actor2)
	rho <- rho_actor; innovation <- sigma_actor2
	m_filt <- P_filt <- numeric(T_per)
	m_prior <- 0; P_prior <- innovation / (1 - rho^2)
	for (t in seq_len(T_per)) {
		denominator <- 1 + H[t] * P_prior
		P_filt[t] <- P_prior / denominator
		m_filt[t] <- (m_prior + P_prior * h[t]) / denominator
		m_prior <- rho * m_filt[t]
		P_prior <- rho^2 * P_filt[t] + innovation
	}
	theta <- numeric(T_per); covariance <- matrix(0, T_per, T_per)
	theta[T_per] <- m_filt[T_per] + sqrt(P_filt[T_per]) * stats::rnorm(1L)
	covariance[T_per, T_per] <- P_filt[T_per]
	if (T_per > 1L) for (t in (T_per - 1L):1L) {
		prediction <- rho^2 * P_filt[t] + innovation
		gain <- rho * P_filt[t] / prediction
		conditional_variance <- P_filt[t] * innovation / prediction
		theta[t] <- m_filt[t] + gain * (theta[t + 1L] - rho * m_filt[t]) +
			sqrt(conditional_variance) * stats::rnorm(1L)
		future <- (t + 1L):T_per
		covariance[t, future] <- gain * covariance[t + 1L, future]
		covariance[future, t] <- covariance[t, future]
		covariance[t, t] <- conditional_variance + gain^2 * covariance[t + 1L, t + 1L]
	}
	list(theta = theta, V = diag(covariance), covariance = covariance)
}

# condition actor paths on a zero sum in every period; V is a matrix of
# per-period variances or a list of each actor's temporal covariance
.exact_center_per_actor <- function(theta_star, V) {
	n <- nrow(theta_star); nt <- ncol(theta_star)
	covariance <- if (is.matrix(V)) lapply(seq_len(n), function(i) diag(V[i, ], nt)) else V
	if (length(covariance) != n) stop("Actor covariance layout does not match its paths.")
	if (n < 2L) return(matrix(0, n, nt))
	L <- chol(Reduce(`+`, covariance))
	adjustment <- backsolve(L, forwardsolve(t(L), colSums(theta_star)))
	out <- theta_star
	for (i in seq_len(n)) out[i, ] <- theta_star[i, ] - drop(covariance[[i]] %*% adjustment)
	out
}

.sweep_per_actor_exact <- function(
	H_mat,
	h_mat,
	rho_actor,
	sigma_actor2,
	precision_list = NULL
) {
	n <- nrow(H_mat); nt <- ncol(H_mat)
	if (n < 2L || nt == 0L) return(matrix(0, n, nt))
	if (!is.null(precision_list)) {
		# reciprocity couples actors, so draw jointly on the centered subspace
		Q <- kronecker(diag(n), .ar1_path_precision(nt, rho_actor, sigma_actor2))
		for (t in seq_len(nt)) {
			idx <- t + (seq_len(n) - 1L) * nt
			Q[idx, idx] <- Q[idx, idx] + precision_list[[t]]
		}
		basis <- kronecker(qr.Q(qr(stats::contr.helmert(n))), diag(nt))
		L <- chol(crossprod(basis, Q %*% basis))
		score <- drop(crossprod(basis, c(t(h_mat))))
		mean <- backsolve(L, forwardsolve(t(L), score))
		return(t(matrix(drop(basis %*% (mean + backsolve(L, stats::rnorm(length(mean))))), nt, n)))
	}
	theta <- matrix(0, n, nt); covariance <- vector("list", n)
	for (i in seq_len(n)) {
		z <- .actor_ffbs_path(H_mat[i, ], h_mat[i, ], rho_actor, sigma_actor2)
		theta[i, ] <- z$theta; covariance[[i]] <- z$covariance
	}
	.exact_center_per_actor(theta, covariance)
}
