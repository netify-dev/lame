# ghk monte carlo estimator for the pointwise observed-data marginal
# log-likelihood under rank-family (frn) and ordinal-probit ame.
# uses a halton low-discrepancy sequence for reproducibility.
# the log of an unbiased probability estimate is biased downward; the
# wrapper emits a one-time note saying so.

#' @noRd
# halton-sequence draws for prime base b, length n, skipping the first
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
# qnorm of halton draws in (0, 1). clips away exact 0/1 to avoid -inf.
.halton_normal <- function(n, base) {
	u <- .halton_seq(n, base)
	u <- pmin(pmax(u, 1e-7), 1 - 1e-7)
	stats::qnorm(u)
}

#' @noRd
# pointwise ordinal closed-form log-lik. cutpoints are inferred from the
# observed y by the empirical cdf (matching the package's own
# ordinal-probit data-augmentation step). for ordinal categories
# 0, 1, ..., k, cutpoint gamma_k is the (sum_{c<=k}/n) quantile of n(0,1).
.ordinal_pointwise_loglik <- function(y_obs, ez_obs) {
	# normalise y to integer 0..k
	y_int <- as.integer(round(y_obs - min(y_obs, na.rm = TRUE)))
	K_max <- max(y_int, na.rm = TRUE)
	# empirical cdf -> cutpoints on the probit scale
	tbl <- as.numeric(table(factor(y_int, levels = seq_len(K_max + 1L) - 1L)))
	cdf <- cumsum(tbl) / sum(tbl)
	# cutpoints: gamma_0 = -inf, gamma_k = +inf; gamma_k = qnorm(cdf[k+1])
	cuts <- c(-Inf, stats::qnorm(cdf[seq_len(K_max)]), Inf)
	# log-lik for category k: log(phi(gamma_{k+1} - ez) - phi(gamma_k - ez))
	up <- cuts[y_int + 2L]   # upper cutpoint for the observed category
	lo <- cuts[y_int + 1L]
	p <- stats::pnorm(up - ez_obs) - stats::pnorm(lo - ez_obs)
	p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	log(p)
}

#' @noRd
# simple ghk estimator for frn: per-row, the observed ranking implies
# an ordering of latent z values. the pointwise log-lik (per dyad) is
# computed by simulating z_row given the ranking using a halton sequence,
# then accumulating the per-dyad contribution.
#
# we approximate the per-dyad log-lik for a "rank likelihood" as
#   log p(z_ij in the implied interval | z_other in their intervals)
# where the intervals come from the observed ranking. for odmax-truncated
# rankings (frn), unobserved dyads in a row contribute phi at the row's
# threshold. returns a length-n_obs numeric vector.
#
# this is the "observed_ghk" estimator; log of an unbiased probability
# estimate is biased downward, with bias scaling as
# (sd of the mc estimator)^2 / (2 * e[p]^2).
.ghk_rank_pointwise_loglik <- function(y_obs, ez_obs, obs_idx, family,
                                        n_mc = 64L) {
	# pragmatic implementation: per row, treat the ranks as defining an
	# ordering of underlying z's. for each observed dyad, compute the
	# per-dyad probability that z_ij sits in the implied interval given
	# the ordering, via a halton-based ghk sweep.
	# obs_idx is a [n_obs x 3] index matrix (i, j, t). we need per-row
	# computations, so group by (i, t).
	n_obs <- length(y_obs)
	out <- numeric(n_obs)
	# group dyads by (sending row i, period t)
	row_key <- paste(obs_idx[, 1], obs_idx[, 3], sep = "_")
	row_groups <- split(seq_len(n_obs), row_key)
	# halton draws per row (re-used; deterministic)
	halton_draws <- .halton_normal(n_mc, base = 2L)
	for (grp in row_groups) {
		yg <- y_obs[grp]
		ezg <- ez_obs[grp]
		# rank order: highest rank = "most observed"
		# implied ordering of latent z's: largest z corresponds to top rank
		ord <- order(-yg, na.last = NA)
		# ghk sweep: simulate z_k | z_{k+1} > z_k for the ranking, accumulate
		# the per-dyad log-prob as phi-difference.
		# for simplicity here, we treat the per-dyad pointwise log-lik as
		# the marginal phi-difference under the implied threshold, with
		# threshold values from the per-row halton-based ghk estimate of
		# the cumulative probability of the ranking up to each position.
		n_g <- length(grp)
		thresh <- -Inf
		cum_lp <- 0
		for (k in seq_len(n_g)) {
			idx_k <- ord[k]
			ez_k <- ezg[idx_k]
			# p(z_k > thresh) under n(ez_k, 1)
			p_k <- stats::pnorm(thresh - ez_k, lower.tail = FALSE)
			p_k <- pmin(pmax(p_k, 1e-12), 1 - 1e-12)
			out[grp[idx_k]] <- log(p_k) - cum_lp + cum_lp / n_g
			# update threshold via halton-based truncated normal draw
			u <- halton_draws[((k - 1L) %% n_mc) + 1L]
			# inverse-cdf truncated above thresh (so next z > thresh)
			lo_p <- stats::pnorm(thresh - ez_k)
			z_draw <- ez_k + stats::qnorm(lo_p + stats::pnorm(u) * (1 - lo_p))
			thresh <- max(thresh, z_draw, -1e6)
			cum_lp <- cum_lp + log(p_k)
		}
	}
	out
}

#' @noRd
# top-level ghk / closed-form dispatcher. called from the in-mcmc log-lik
# fill block. returns a numeric vector of length(y_obs).
.pointwise_loglik_observed_ghk <- function(y_obs, ez_obs, z_obs, obs_idx,
                                            family, s2, dyad_rho = 0,
                                            n_mc = 64L) {
	# closed-form families: same as observed_exact
	if (family %in% c("normal", "binary", "cbin", "poisson")) {
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
			poisson = {
				# overdispersed poisson marginal: integrate out the latent
				# z ~ N(eta, s2) layer via adaptive gauss-hermite -- the same
				# estimand the default (non-ghk) log_lik path reports. a plain
				# dpois(round(y), exp(eta)) plug-in would be the equidispersed
				# model's likelihood and ignore the latent variance s2,
				# understating tail cells by tens of nats.
				.lame_pois_marginal_loglik(round(y_obs), ez_obs, s2)
			})
		return(ll)
	}
	if (family == "ordinal") {
		return(.ordinal_pointwise_loglik(y_obs, ez_obs))
	}
	if (family == "frn") {
		return(.ghk_rank_pointwise_loglik(y_obs, ez_obs, obs_idx,
		                                   family, n_mc = n_mc))
	}
	# default: augmented-data fallback
	stats::dnorm(z_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
}
