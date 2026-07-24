#' Metropolis-Hastings update for the within-dyad correlation
#'
#' Draws a new value of the dyadic correlation rho from its full conditional
#' using a random-walk-style Metropolis-Hastings step. Each dyad's pair of
#' directed residuals is modelled as bivariate normal with unit variances
#' (after scaling by the dyadic variance) and correlation rho; the proposal is
#' a normal truncated to (-1, 1) and the acceptance ratio combines the
#' bivariate-normal log-likelihood, an optional arc-sine prior, and the
#' truncated-proposal normaliser.
#'
#' @param Z n x n normal relational matrix.
#' @param rho current value of rho.
#' @param s2 current value of the dyadic variance.
#' @param offset matrix matching \code{Z}; \code{Z - offset} is treated as
#'   dyadic noise (so it should absorb any additive/multiplicative effects).
#' @param asp logical; use an arc-sine prior (\code{TRUE}, the default when
#'   \code{NULL}) or a uniform prior (\code{FALSE}).
#'
#' @return a new (scalar) value of rho.
#' @keywords internal
#' @author lame authors
#' @export rrho_mh
rrho_mh <-
function(Z, rho, s2 = 1, offset = 0, asp = NULL) {

	if (is.null(asp)) { asp <- TRUE }

	## --- standardized directed residual pairs, one per dyad -------------
	resid <- (Z - offset) / sqrt(s2)
	ut <- upper.tri(resid)
	dir_out <- resid[ut]			# e_ij for i < j
	dir_in  <- t(resid)[ut]			# e_ji for i < j
	n_dyad <- length(dir_out)

	## sufficient statistics for the bivariate-normal likelihood
	cross_sum <- sum(dir_out * dir_in)			# sum e_ij e_ji
	sq_sum    <- sum(dir_out^2) + sum(dir_in^2)	# sum e_ij^2 + e_ji^2

	## --- truncated-normal proposal on (-1, 1) ---------------------------
	prop_sd <- 2 * (1 - cor(dir_out, dir_in)^2) / sqrt(n_dyad)
	lo <- pnorm((-1 - rho) / prop_sd)
	hi <- pnorm(( 1 - rho) / prop_sd)
	rho_new <- rho + prop_sd * qnorm(runif(1, lo, hi))

	## proposals outside the stable interior are rejected outright, keeping
	## the chain within the truncated support (-0.995, 0.995)
	if (!is.finite(rho_new) || abs(rho_new) > 0.995) { return(rho) }

	## joint bivariate-normal log-likelihood (unit variances, corr r)
	loglik <- function(r) {
		omr2 <- 1 - r^2
		-0.5 * (n_dyad * log(omr2) + (sq_sum - 2 * r * cross_sum) / omr2)
	}

	## arc-sine prior log-density (up to an additive constant)
	logprior <- function(r) {
		if (asp) { -0.5 * log(1 - r^2) } else { 0 }
	}

	## log normaliser of the (-1, 1)-truncated normal proposal centred at r;
	## required so q(rho | rho_new) / q(rho_new | rho) enters the MH ratio
	logqnorm <- function(r) {
		log(pnorm((1 - r) / prop_sd) - pnorm((-1 - r) / prop_sd))
	}

	log_accept <- (loglik(rho_new)   - loglik(rho)) +
				  (logprior(rho_new) - logprior(rho)) +
				  (logqnorm(rho)     - logqnorm(rho_new))

	if (log_accept > log(runif(1))) { rho <- rho_new }

	rho
}
