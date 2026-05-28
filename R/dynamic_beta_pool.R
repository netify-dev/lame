# hierarchical pool across blocks for the dynamic_beta hyperparameters.
# sigma case: sigma_b^2 ~ ig(alpha, beta_share) with beta_share drawn
# from a gamma hyperprior shared across blocks. closed-form gamma
# conditional on sum_b 1/sigma_b^2.
# rho case: logit(rho_b) ~ n(mu_rho, tau_rho^2). normal density on
# the logit scale; no jacobian correction.

#' @noRd
# update the shared beta_share hyperparameter for sigma_b^2.
.pool_update_sigma_scale <- function(sigma_beta_by_group,
                                       alpha = 2.0,
                                       hyper_a = 2.0,
                                       hyper_b = 0.5) {
	B <- length(sigma_beta_by_group)
	if (B == 0L) return(hyper_b)
	# inverse-variance contributions across blocks
	inv_var <- 1 / pmax(sigma_beta_by_group^2, 1e-12)
	stats::rgamma(1L,
		shape = hyper_a + B * alpha,
		rate  = hyper_b + sum(inv_var))
}

#' @noRd
# update (mu_rho, tau_rho^2) governing per-block rho_b on the logit scale.
.pool_update_rho_hyper <- function(rho_beta_by_group,
                                     prior_mu_mean = 0,
                                     prior_mu_sd   = 10,
                                     prior_tau_a   = 2.0,
                                     prior_tau_b   = 0.5) {
	# transform: theta_b = logit((rho + 1) / 2) so theta_b in R
	r <- pmin(pmax(rho_beta_by_group, -0.999), 0.999)
	theta <- log((r + 1) / (1 - r) / 2 + 1 / (1 - r))  # logit((r+1)/2)
	# simpler & equivalent: theta_b = 0.5 * log((1+r)/(1-r))
	theta <- 0.5 * log((1 + r) / (1 - r))
	B <- length(theta)
	# mu | tau, theta : normal
	prior_prec_mu <- 1 / prior_mu_sd^2
	# tau | mu, theta : inverse-gamma (need an initial tau to use here; use sample variance)
	tau_sq <- max(stats::var(theta), 1e-6)
	post_prec_mu <- prior_prec_mu + B / tau_sq
	post_mean_mu <- (prior_prec_mu * prior_mu_mean + sum(theta) / tau_sq) / post_prec_mu
	mu_new <- stats::rnorm(1, post_mean_mu, 1 / sqrt(post_prec_mu))
	# tau update: 1/tau^2 | mu, theta ~ Gamma(prior_tau_a + B/2,
	#                                          prior_tau_b + 0.5 * sum((theta - mu)^2))
	inv_tau_sq <- stats::rgamma(1L,
		shape = prior_tau_a + B / 2,
		rate  = prior_tau_b + 0.5 * sum((theta - mu_new)^2))
	tau_new <- 1 / sqrt(inv_tau_sq)
	list(mu = mu_new, tau = tau_new)
}
