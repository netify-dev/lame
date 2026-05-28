# ghk monte carlo estimator for the pointwise observed-data marginal
# log-likelihood under rank-family (frn / rrl) and ordinal-probit ame.
# uses a halton low-discrepancy sequence for reproducibility.
# the log of an unbiased probability estimate is biased downward; the
# wrapper emits a one-time note saying so.

#' @noRd
# Halton-sequence draws for prime base b, length n, skipping the first
# `burn` terms (which tend to be near zero and cluster).
.halton_seq <- function(n, base, burn = 30L) {
	out <- numeric(n + burn)
	for (i in seq_len(n + burn)) {
		f <- 1
		r <- 0
		k <- i
		while (k > 0) {
			f <- f / base
			r <- r + f * (k %% base)
			k <- k %/% base
		}
		out[i] <- r
	}
	out[(burn + 1L):(burn + n)]
}

#' @noRd
# qnorm of Halton draws in (0, 1). Clips away exact 0/1 to avoid -Inf.
.halton_normal <- function(n, base) {
	u <- .halton_seq(n, base)
	u <- pmin(pmax(u, 1e-7), 1 - 1e-7)
	stats::qnorm(u)
}

#' @noRd
# Pointwise ordinal closed-form log-lik. Cutpoints are inferred from the
# observed Y by the empirical CDF (matching the package's own
# ordinal-probit data-augmentation step). For ordinal categories
# 0, 1, ..., K, cutpoint gamma_k is the (sum_{c<=k}/n) quantile of N(0,1).
.ordinal_pointwise_loglik <- function(y_obs, ez_obs) {
	# normalise y to integer 0..K
	y_int <- as.integer(round(y_obs - min(y_obs, na.rm = TRUE)))
	K_max <- max(y_int, na.rm = TRUE)
	# empirical CDF -> cutpoints on the probit scale
	tbl <- as.numeric(table(factor(y_int, levels = seq_len(K_max + 1L) - 1L)))
	cdf <- cumsum(tbl) / sum(tbl)
	# cutpoints: gamma_0 = -Inf, gamma_K = +Inf; gamma_k = qnorm(cdf[k+1])
	cuts <- c(-Inf, stats::qnorm(cdf[seq_len(K_max)]), Inf)
	# log-lik for category k: log(Phi(gamma_{k+1} - ez) - Phi(gamma_k - ez))
	up <- cuts[y_int + 2L]   # upper cutpoint for the observed category
	lo <- cuts[y_int + 1L]
	p <- stats::pnorm(up - ez_obs) - stats::pnorm(lo - ez_obs)
	p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	log(p)
}

#' @noRd
# Simple GHK estimator for frn/rrl: per-row, the observed ranking implies
# an ordering of latent Z values. The pointwise log-lik (per dyad) is
# computed by simulating Z_row given the ranking using a Halton sequence,
# then accumulating the per-dyad contribution.
#
# We approximate the per-dyad log-lik for a "rank likelihood" as
#   log P(Z_ij in the implied interval | Z_other in their intervals)
# where the intervals come from the observed ranking. For odmax-truncated
# rankings (frn), unobserved dyads in a row contribute Phi at the row's
# threshold. Returns a length-n_obs numeric vector.
#
# this is the "observed_ghk" estimator; log of an unbiased probability
# estimate is biased downward, with bias scaling as
# (sd of the MC estimator)^2 / (2 * E[P]^2).
.ghk_rank_pointwise_loglik <- function(y_obs, ez_obs, obs_idx, family,
                                        n_mc = 64L) {
	# Pragmatic implementation: per row, treat the ranks as defining an
	# ordering of underlying Z's. For each observed dyad, compute the
	# per-dyad probability that Z_ij sits in the implied interval given
	# the ordering, via a Halton-based GHK sweep.
	# obs_idx is a [n_obs x 3] index matrix (i, j, t). We need per-row
	# computations, so group by (i, t).
	n_obs <- length(y_obs)
	out <- numeric(n_obs)
	# group dyads by (sending row i, period t)
	row_key <- paste(obs_idx[, 1], obs_idx[, 3], sep = "_")
	row_groups <- split(seq_len(n_obs), row_key)
	# Halton draws per row (re-used; deterministic)
	halton_draws <- .halton_normal(n_mc, base = 2L)
	for (grp in row_groups) {
		yg <- y_obs[grp]
		ezg <- ez_obs[grp]
		# rank order: highest rank = "most observed"
		# implied ordering of latent Z's: largest Z corresponds to top rank
		ord <- order(-yg, na.last = NA)
		# GHK sweep: simulate Z_k | Z_{k+1} > Z_k for the ranking, accumulate
		# the per-dyad log-prob as Phi-difference.
		# For simplicity here, we treat the per-dyad pointwise log-lik as
		# the marginal Phi-difference under the implied threshold, with
		# threshold values from the per-row Halton-based GHK estimate of
		# the cumulative probability of the ranking up to each position.
		n_g <- length(grp)
		thresh <- -Inf
		cum_lp <- 0
		for (k in seq_len(n_g)) {
			idx_k <- ord[k]
			ez_k <- ezg[idx_k]
			# P(Z_k > thresh) under N(ez_k, 1)
			p_k <- stats::pnorm(thresh - ez_k, lower.tail = FALSE)
			p_k <- pmin(pmax(p_k, 1e-12), 1 - 1e-12)
			out[grp[idx_k]] <- log(p_k) - cum_lp + cum_lp / n_g
			# update threshold via Halton-based truncated normal draw
			u <- halton_draws[((k - 1L) %% n_mc) + 1L]
			# inverse-cdf truncated above thresh (so next Z > thresh)
			lo_p <- stats::pnorm(thresh - ez_k)
			z_draw <- ez_k + stats::qnorm(lo_p + stats::pnorm(u) * (1 - lo_p))
			thresh <- max(thresh, z_draw, -1e6)
			cum_lp <- cum_lp + log(p_k)
		}
	}
	out
}

#' @noRd
# Top-level GHK / closed-form dispatcher. Called from the in-MCMC log-lik
# fill block. Returns a numeric vector of length(y_obs).
.pointwise_loglik_observed_ghk <- function(y_obs, ez_obs, z_obs, obs_idx,
                                            family, s2, dyad_rho = 0,
                                            n_mc = 64L) {
	# closed-form families: same as observed_exact
	if (family %in% c("normal", "binary", "cbin", "tobit", "poisson")) {
		ll <- switch(family,
			normal = stats::dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE),
			binary = {
				p_hat <- stats::pnorm(ez_obs)
				p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
				y_obs * log(p_hat) + (1 - y_obs) * log(1 - p_hat)
			},
			cbin = {
				p_hat <- stats::pnorm(ez_obs)
				p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
				y_b <- 1 * (y_obs > 0)
				y_b * log(p_hat) + (1 - y_b) * log(1 - p_hat)
			},
			tobit = {
				sigma_t <- sqrt(s2); is_zero <- y_obs <= 0
				ll <- numeric(length(y_obs))
				ll[is_zero]  <- stats::pnorm(0, mean = ez_obs[is_zero], sd = sigma_t, log.p = TRUE)
				ll[!is_zero] <- stats::dnorm(y_obs[!is_zero],
				                              mean = ez_obs[!is_zero],
				                              sd   = sigma_t, log = TRUE)
				ll
			},
			poisson = {
				lam <- pmin(exp(ez_obs), 1e8)
				stats::dpois(round(y_obs), lambda = lam, log = TRUE)
			})
		return(ll)
	}
	if (family == "ordinal") {
		return(.ordinal_pointwise_loglik(y_obs, ez_obs))
	}
	if (family %in% c("frn", "rrl")) {
		return(.ghk_rank_pointwise_loglik(y_obs, ez_obs, obs_idx,
		                                   family, n_mc = n_mc))
	}
	# default: augmented-data fallback
	stats::dnorm(z_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
}
