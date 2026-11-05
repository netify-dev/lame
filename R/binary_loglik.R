# stable probit pointwise log likelihood; taking the observed tail avoids clipping and 0 * -Inf
.binary_pointwise_loglik <- function(y_obs, ez_obs) {
	stats::pnorm(ifelse(y_obs > 0, ez_obs, -ez_obs), log.p = TRUE)
}
