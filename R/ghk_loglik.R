# closed-form pointwise log-likelihoods for the observed_ghk scoring route (rank estimators are in
# R/ghk_loglik_rank.R)

#' @noRd
# pointwise ordinal log-likelihood; cutpoints are probit quantiles of the empirical cdf of y,
# matching the ordinal data-augmentation step
.ordinal_pointwise_loglik <- function(y_obs, ez_obs) {
	# recode y as integers 0..k
	y_int <- as.integer(round(y_obs - min(y_obs, na.rm = TRUE)))
	K_max <- max(y_int, na.rm = TRUE)
	# cutpoints on the probit scale from the empirical cdf, with gamma_0 = -Inf and gamma_{K+1} = Inf
	tbl <- as.numeric(table(factor(y_int, levels = seq_len(K_max + 1L) - 1L)))
	cdf <- cumsum(tbl) / sum(tbl)
	cuts <- c(-Inf, stats::qnorm(cdf[seq_len(K_max)]), Inf)
	# category k: log(Phi(gamma_{k+1} - ez) - Phi(gamma_k - ez))
	up <- cuts[y_int + 2L] # upper cutpoint for the observed category
	lo <- cuts[y_int + 1L]
	p <- stats::pnorm(up - ez_obs) - stats::pnorm(lo - ez_obs)
	p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	log(p)
}

#' @noRd
# closed-form pointwise log-likelihood for the families the rank dispatcher
# (.pointwise_loglik_observed_ghk_rank) delegates: normal, binary and Poisson-lognormal
.pointwise_loglik_observed_ghk <- function(
	y_obs,
	ez_obs,
	z_obs,
	obs_idx,
	family,
	s2,
	dyad_rho = 0,
	n_mc = 64L
) {
	switch(
		family,
		normal = stats::dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE),
		binary = .binary_pointwise_loglik(y_obs, ez_obs),
		# overdispersed Poisson marginal: z ~ N(eta, s2) integrated out by adaptive Gauss-Hermite, the
		# same estimand as the default log_lik
		poisson = .lame_pois_marginal_loglik(round(y_obs), ez_obs, s2),
		cli::cli_abort("{.fn .pointwise_loglik_observed_ghk} scores the normal, binary and Poisson families only.")
	)
}
