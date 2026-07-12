# per-actor dynamic beta with exact gaussian conditioning for the
# sum-to-zero centering constraint. the path:
#
#   (1) for each actor i, run a univariate carter-kohn ffbs on the
#       length-t path theta_{i, 1:t} given period-wise (h_i(t), h_i(t))
#       and ar(1) hyperparameters (rho_actor, sigma_actor^2).
#       this is the unconstrained per-actor posterior conditional on
#       hyperparameters and the other actors (in the lame model, actor
#       slopes are conditionally independent given hypers).
#
#   (2) project the unconstrained sample theta* onto the sum-to-zero
#       manifold at each period using the variance-weighted lagrangian:
#         theta_i(t) <- theta*_i(t) - v_i(t) / sum_j v_j(t) * sum_j theta*_j(t)
#       where v_i(t) is the *marginal* posterior variance of theta_i(t)
#       under the ffbs posterior. this is the exact conditioning under
#       the block-diagonal-across-actors covariance assumption, which
#       holds in the lame per-actor model.
#
# derivation: if (theta*_1, ..., theta*_n) ~ n(mu, sigma) with
# sigma block-diagonal across actors (sigma_i is the per-actor t x t
# carter-kohn covariance), the conditional given sum_i theta_i = 0 is
#   theta_i | constraint = theta*_i - sigma_i 1 (1' sigma 1)^{-1} 1' theta*
# where the projection 1' theta* is the period-wise sum. per period t:
#   theta_i(t) = theta*_i(t) - v_i(t) * (sum_j theta*_j(t)) / sum_j v_j(t)
# this is exact (not heuristic) under block-diagonality.

#' Univariate Carter-Kohn FFBS for one actor's length-T slope path
#'
#' Given period-wise sufficient statistics (H[t] = sum_j X^2_ij / s2,
#' h[t] = sum_j X_ij R_ij / s2) and AR(1) hyperparameters, returns one
#' joint draw from the posterior N(m, P) where the prior is AR(1) and
#' the observation is the Gaussian likelihood implied by (H, h).
#'
#' @param H length-T vector of period-wise observation precisions
#' @param h length-T vector of period-wise observation cross-products
#' @param rho_actor AR(1) coefficient
#' @param sigma_actor2 AR(1) innovation variance
#' @return list with \code{theta} (length-T draw) and \code{V} (length-T
#'   marginal posterior variance per period -- needed for the exact
#'   centering projection).
#' @keywords internal
.actor_ffbs_path <- function(H, h, rho_actor, sigma_actor2) {
	T_per <- length(H)
	if (T_per < 1L) return(list(theta = numeric(0), V = numeric(0)))
	# clamp rho away from boundaries for numerical stability
	rho <- max(min(rho_actor, 0.99999), -0.99999)
	s2 <- max(sigma_actor2, 1e-12)
	# stationary prior var at t = 1
	v_stat <- s2 / max(1 - rho^2, 1e-6)
	# forward filter
	m_filt <- numeric(T_per)
	P_filt <- numeric(T_per)
	m_prior <- 0
	P_prior <- v_stat
	for (t in seq_len(T_per)) {
		# information form: p_post^{-1} = p_prior^{-1} + h[t]
		P_t <- 1 / (1 / P_prior + H[t])
		m_t <- P_t * (m_prior / P_prior + h[t])
		m_filt[t] <- m_t
		P_filt[t] <- P_t
		if (t < T_per) {
			m_prior <- rho * m_t
			P_prior <- (rho^2) * P_t + s2
		}
	}
	# backward sample (carter-kohn)
	theta <- numeric(T_per)
	V_marg <- numeric(T_per)  # marginal posterior variance per period
	# draw theta[t]
	theta[T_per] <- m_filt[T_per] + sqrt(P_filt[T_per]) * stats::rnorm(1)
	V_marg[T_per] <- P_filt[T_per]
	if (T_per >= 2L) for (t in (T_per - 1L):1L) {
		# smoother
		P_pred_next <- (rho^2) * P_filt[t] + s2
		J <- rho * P_filt[t] / P_pred_next
		m_smooth <- m_filt[t] + J * (theta[t + 1L] - rho * m_filt[t])
		P_smooth <- P_filt[t] - J * (rho * P_filt[t])
		P_smooth <- max(P_smooth, 1e-12)
		theta[t] <- m_smooth + sqrt(P_smooth) * stats::rnorm(1)
		V_marg[t] <- P_smooth
	}
	list(theta = theta, V = V_marg)
}

#' Exact variance-weighted sum-to-zero projection for the per-actor block
#'
#' Given an unconstrained sample \code{theta_star} (n_actors x T) and
#' per-actor per-period posterior variances \code{V} (same shape),
#' projects to the sum-to-zero manifold per period using the
#' variance-weighted Lagrangian projection (exact under block-diagonal
#' covariance across actors).
#'
#' @keywords internal
.exact_center_per_actor <- function(theta_star, V) {
	stopifnot(dim(theta_star) == dim(V))
	T_per <- ncol(theta_star)
	out <- theta_star
	for (t in seq_len(T_per)) {
		s_t <- sum(theta_star[, t])
		V_sum <- sum(V[, t])
		if (V_sum < 1e-12) next  # degenerate (no variance signal); skip
		shrink <- s_t / V_sum
		out[, t] <- theta_star[, t] - V[, t] * shrink
	}
	out
}

#' Per-actor sweep: unconstrained per-actor FFBS + exact projection
#'
#' One full sweep over actors: each actor's path is sampled from its
#' conditional AR(1) posterior via .actor_ffbs_path(), and the joint
#' sample is then projected to satisfy sum_i theta_i(t) = 0 per period
#' using the variance-weighted Lagrangian.
#'
#' @param H_mat n_actors x T precision sufficient statistics
#' @param h_mat n_actors x T cross-product sufficient statistics
#' @param rho_actor AR(1) coefficient
#' @param sigma_actor2 AR(1) innovation variance
#' @return n_actors x T matrix of centered draws.
#' @keywords internal
.sweep_per_actor_exact <- function(H_mat, h_mat, rho_actor, sigma_actor2) {
	n_actors <- nrow(H_mat)
	T_per <- ncol(H_mat)
	if (n_actors == 0L || T_per == 0L) return(matrix(0, n_actors, T_per))
	theta_star <- matrix(0, n_actors, T_per)
	V <- matrix(0, n_actors, T_per)
	for (i in seq_len(n_actors)) {
		out_i <- .actor_ffbs_path(H_mat[i, ], h_mat[i, ],
		                          rho_actor, sigma_actor2)
		theta_star[i, ] <- out_i$theta
		V[i, ] <- out_i$V
	}
	# exact projection to sum-to-zero manifold
	.exact_center_per_actor(theta_star, V)
}
