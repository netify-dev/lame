# hierarchical pool across blocks for the dynamic_beta hyperparameters.
# sigma: sigma_b^2 ~ IG(alpha, beta_share), beta_share with a shared gamma
# hyperprior and a closed-form gamma conditional.
# rho: a Gaussian hierarchy on the persistence mapped from its interval to the
# real line, with proposals on that scale.

#' @noRd
# draw the shared beta_share hyperparameter of sigma_b^2
.pool_update_sigma_scale <- function(
	sigma_beta_by_group,
	alpha = 2.0,
	hyper_a = 2.0,
	hyper_b = 0.5
) {
	B <- length(sigma_beta_by_group)
	if (B == 0L) return(hyper_b)
	if (any(!is.finite(sigma_beta_by_group)) || any(sigma_beta_by_group <= 0))
		stop("Coefficient innovation standard deviations must be finite and positive.")
	inv_var <- 1 / sigma_beta_by_group^2
	if (any(!is.finite(inv_var))) stop("Coefficient innovation precision is outside numerical range.")
	stats::rgamma(
		1L,
		shape = hyper_a + B * alpha,
		rate = hyper_b + sum(inv_var)
	)
}

# Gibbs draws of (mu, tau) of the Gaussian hierarchy on
# theta = logit((rho - lower) / (upper - lower)), mu given the current tau
.pool_update_rho_hyper <- function(
	rho_beta_by_group,
	prior_mu_mean = 0,
	prior_mu_sd = 10,
	prior_tau_a = 2,
	prior_tau_b = .5,
	current_tau = 1,
	rho_lower = 0,
	rho_upper = .999
) {
	theta <- stats::qlogis((rho_beta_by_group - rho_lower) / (rho_upper - rho_lower))
	if (any(!is.finite(theta))) stop("Pooled persistence must lie strictly within its support.")
	B <- length(theta)
	precision <- 1 / prior_mu_sd^2 + B / current_tau^2
	mean <- (prior_mu_mean / prior_mu_sd^2 + sum(theta) / current_tau^2) / precision
	mu <- rnorm(1, mean, 1 / sqrt(precision))
	tau <- sqrt(
		1 / rgamma(
			1,
			prior_tau_a + B / 2,
			prior_tau_b + sum((theta - mu)^2) / 2
		)
	)
	list(mu = mu, tau = tau)
}

# random-walk Metropolis on theta per block, under the Gaussian hierarchy and
# the gap-aware transition density
.pool_sample_rho <- function(
	beta_path,
	group_id,
	n_groups,
	Lambda,
	sigma_by_coef,
	rho_current,
	mu,
	tau,
	rho_lower,
	rho_upper,
	time_positions = seq_len(nrow(beta_path))
) {
	gaps <- diff(time_positions)
	theta <- stats::qlogis((rho_current - rho_lower) / (rho_upper - rho_lower))
	for (g in seq_len(n_groups)) {
		idx <- which(group_id == g)
		x <- beta_path[, idx, drop = FALSE]
		log_target <- function(th) {
			rho <- rho_lower + (rho_upper - rho_lower) * stats::plogis(th)
			if (rho < 0 && any(abs(gaps - round(gaps)) > 1e-10)) return(-Inf)
			z <- .lame_beta_gap_transition(rho, gaps)
			err <- x[-1L, , drop = FALSE] - z$f * x[-nrow(x), , drop = FALSE]
			ss <- colSums(err^2 / z$q) / (diag(Lambda)[idx] * sigma_by_coef[idx]^2)
			dnorm(th, mu, tau, log = TRUE) - .5 * (sum(ss) + length(idx) * sum(log(z$q)))
		}
		for (step in seq_len(4L)) {
			proposal <- theta[g] + rnorm(1, sd = .5)
			lr <- log_target(proposal) - log_target(theta[g])
			if (is.finite(lr) && log(runif(1)) < lr) theta[g] <- proposal
		}
	}
	rho_lower + (rho_upper - rho_lower) * stats::plogis(theta)
}
