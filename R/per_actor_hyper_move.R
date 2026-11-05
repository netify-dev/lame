# interweaving step for the AR(1) hyperparameters (sigma, rho) of the
# actor-specific slopes of a dynamic_beta_per_actor fit. with the standardised
# innovations w held fixed,
#     theta_i1 = sigma / sqrt(1 - rho^2) w_i1,   theta_it = rho theta_i,t-1 + sigma w_it,
# a symmetric walk on (log sigma, logit rho) rebuilds the paths and accepts with
# the Metropolis ratio of prior times Gaussian likelihood, which crosses the
# funnel of the centred sampler. rebuilt paths keep their zero sum over actors,
# and the walk scale comes from a fixed ladder so the step keeps no state. the
# exact conditional draws of sigma^2 and rho follow at every sweep.

# the paths of fixed standardised innovations under (sigma, rho): `w` is the actors x
# periods matrix of innovations, `sigma` the standard deviation of an innovation
.lame_pa_rebuild <- function(w, rho, sigma) {
	th <- w
	th[, 1L] <- sigma / sqrt(1 - rho^2) * w[, 1L]
	if (ncol(w) > 1L) for (t in 2:ncol(w)) th[, t] <- rho * th[, t - 1L] + sigma * w[, t]
	th
}

# the standardised innovations of the paths `theta` under (sigma, rho)
.lame_pa_innovations <- function(theta, rho, sigma) {
	w <- theta
	w[, 1L] <- theta[, 1L] * sqrt(1 - rho^2) / sigma
	if (ncol(theta) > 1L) for (t in 2:ncol(theta)) w[, t] <- (theta[, t] - rho * theta[, t - 1L]) / sigma
	w
}

# log likelihood of the paths under the Gaussian statistics (up to a constant):
# `H_mat`, `h_mat` the actors x periods precision diagonal and score, `P_list`
# NULL for independent cells or the period precision matrices (reciprocity)
.lame_pa_loglik <- function(theta, H_mat, h_mat, P_list = NULL) {
	if (is.null(P_list)) return(sum(h_mat * theta) - 0.5 * sum(H_mat * theta^2))
	out <- 0
	for (t in seq_len(ncol(theta))) {
		th <- theta[, t]
		out <- out + sum(h_mat[, t] * th) - 0.5 * sum(th * (P_list[[t]] %*% th))
	}
	out
}

# one sweep of the step. `theta` the actors x periods paths (zero sum over the
# actors in every period), `sigma2` the innovation variance, `rho` the
# persistence, `shape`/`scale` the inverse-gamma prior of sigma2 and `a`/`b` the
# shapes of the scaled beta prior of (rho + 1) / 2 on (-0.99, 0.99)
.lame_pa_nc_step <- function(
	theta,
	H_mat,
	h_mat,
	P_list,
	rho,
	sigma2,
	shape,
	scale,
	a,
	b,
	tries = 6L,
	scales = c(0.03, 0.1, 0.3, 0.8)
) {
	out <- list(theta = theta, sigma2 = sigma2, rho = rho, attempts = 0L, accepts = 0L)
	nt <- ncol(theta)
	if (nrow(theta) < 2L || nt < 1L) return(out)
	lo <- -0.99
	rng <- 1.98
	log_prior <- function(sig, r) {
		u <- (r - lo) / rng
		if (!is.finite(u) || u <= 0 || u >= 1) return(-Inf)
		-2 * shape * log(sig) - scale / sig^2 + stats::dbeta((r + 1) / 2, a, b, log = TRUE) + log(u) + log1p(-u)
	}
	sig0 <- sqrt(sigma2)
	rho0 <- rho
	w <- .lame_pa_innovations(theta, rho0, sig0)
	ll_cur <- .lame_pa_loglik(theta, H_mat, h_mat, P_list)
	lp_cur <- log_prior(sig0, rho0)
	if (!is.finite(ll_cur) || !is.finite(lp_cur)) return(out)
	th_cur <- theta
	for (k in seq_len(tries)) {
		sc <- scales[sample.int(length(scales), 1L)]
		mv <- sample.int(3L, 1L)
		sig1 <- if (mv != 2L) sig0 * exp(sc * stats::rnorm(1L)) else sig0
		rho1 <- rho0
		if (mv >= 2L) {
			th <- stats::qlogis((rho0 - lo) / rng) + 3 * sc * stats::rnorm(1L)
			rho1 <- lo + rng * stats::plogis(th)
		}
		lp_new <- log_prior(sig1, rho1)
		if (!is.finite(lp_new)) next
		th_new <- .lame_pa_rebuild(w, rho1, sig1)
		out$attempts <- out$attempts + 1L
		if (!all(is.finite(th_new))) next
		ll_new <- .lame_pa_loglik(th_new, H_mat, h_mat, P_list)
		if (is.finite(ll_new) && log(stats::runif(1L)) < ll_new - ll_cur + lp_new - lp_cur) {
			th_cur <- th_new
			ll_cur <- ll_new
			lp_cur <- lp_new
			sig0 <- sig1
			rho0 <- rho1
			out$accepts <- out$accepts + 1L
		}
	}
	out$theta <- th_cur
	out$sigma2 <- sig0^2
	out$rho <- rho0
	out
}
