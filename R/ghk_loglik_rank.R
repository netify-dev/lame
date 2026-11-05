# GHK, Halton and Gauss-Hermite pointwise log likelihoods for the rank families.
#
# frn rows use a Geweke-Hajivassiliou-Keane simulator on the implied ordering:
# sequential conditioning with truncated normal draws, Halton sequences
# randomised by an Owen-style digital shift stored on the fit, and a cap of 15
# dimensions above which a row uses adjacent pairwise comparisons. a row's
# estimate is split evenly across its observed dyads, keeping the row total.
# cbin rows integrate the implicit row intercept by Gauss-Hermite quadrature.
# every Monte Carlo component is deterministic given fit$log_lik_meta$seed.

#' @noRd
# coprime prime bases for a multidimensional Halton sequence: integration
# dimension k of the GHK estimator uses base k (enough for the dimension cap
# plus the bank offset ladder)
.ghk_prime_bases <- c(
	2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L,
	37L, 41L, 43L, 47L, 53L, 59L, 61L, 67L, 71L, 73L,
	79L, 83L, 89L, 97L
)

# randomised Halton (generalised van der Corput) sequence in base b with an
# Owen-style digital shift (Faure 1992)
.halton_seq_randomised <- function(n, base, burn = 30L, shift = 0) {
	out <- numeric(n + burn)
	for (i in seq_len(n + burn)) {
		f <- 1; r <- 0; k <- i
		while (k > 0) {
			f <- f / base
			# digital shift of each digit, mod base
			r <- r + f * ((k %% base) + shift) %% base
			k <- k %/% base
		}
		out[i] <- r %% 1
	}
	out[(burn + 1L):(burn + n)]
}

#' @noRd
# row-level GHK estimate of log P(z_{perm[1]} > ... > z_{perm[d]}) for the frn
# log likelihood, given the expected predictors ez and the observed ordering
# perm (perm[1] = top rank), averaged over n_mc Halton draws
.ghk_row_loglik <- function(
	ez,
	perm,
	n_mc = 64L,
	halton_shift = 0,
	halton_base = 2L
) {
	D <- length(ez)
	if (D == 0L) return(0)
	if (D == 1L) return(0) # single-dyad row is trivially in order
	if (D > 15L) {
		# above the dimension cap, sum the log probabilities of the adjacent
		# comparisons (biased but tractable)
		ez_ord <- ez[perm]
		lp <- 0
		for (k in seq_len(D - 1L)) {
			# z_k - z_{k+1} ~ N(ez_k - ez_{k+1}, 2)
			diff_mu <- ez_ord[k] - ez_ord[k + 1L]
			lp <- lp + stats::pnorm(diff_mu / sqrt(2), log.p = TRUE)
		}
		return(lp)
	}
	# GHK over the implied ordering, with Halton-normal draws for the d - 1
	# integration dimensions. each dimension uses its own coprime base, since a
	# single base gives a correlated lattice that biases the rank probability
	# upward; halton_base offsets the ladder so bank members stay coprime
	pr_ladder <- .ghk_prime_bases
	off <- match(halton_base, pr_ladder, nomatch = 1L) - 1L
	hn <- matrix(0, D - 1L, n_mc)
	for (d in seq_len(D - 1L)) {
		base_d <- pr_ladder[((off + d - 1L) %% length(pr_ladder)) + 1L]
		hn[d, ] <- .halton_seq_randomised(n_mc, base_d, shift = halton_shift)
	}
	# clip away exact 0/1
	hn <- pmin(pmax(hn, 1e-7), 1 - 1e-7)
	ez_ord <- ez[perm]
	lp_mc <- numeric(n_mc)
	for (m in seq_len(n_mc)) {
		thresh <- -Inf
		lp <- 0
		# sweep ranks from k = d down to 1 (perm[1] is the top rank), thresh the
		# running maximum of the lower ranks: add log(1 - Phi(thresh - ez_perm[k])),
		# then draw z_perm[k] from TN(ez_perm[k], 1; thresh, Inf)
		for (k_idx in D:1L) {
			ez_k <- ez_ord[k_idx]
			if (is.infinite(thresh)) {
				p_k <- 1 # unconstrained probability
			} else {
				p_k <- stats::pnorm(thresh - ez_k, lower.tail = FALSE)
				p_k <- pmin(pmax(p_k, 1e-15), 1 - 1e-15)
				lp <- lp + log(p_k)
			}
			# draw z_perm[k] from the truncated normal
			if (k_idx > 1L) {
				lo_cdf <- if (is.infinite(thresh)) 0 else stats::pnorm(thresh - ez_k)
				# Halton draw for dimension (k_idx - 1, m)
				u_draw <- hn[((D - k_idx)) %% (D - 1L) + 1L, m]
				z_draw <- ez_k + stats::qnorm(lo_cdf + u_draw * (1 - lo_cdf))
				if (is.finite(z_draw)) thresh <- z_draw
			}
		}
		lp_mc[m] <- lp
	}
	# log-mean-exp over the Monte Carlo sample
	mx <- max(lp_mc)
	if (!is.finite(mx)) return(-1e10)
	mx + log(mean(exp(lp_mc - mx)))
}

#' @noRd
# per-row GHK with per-dyad attribution for frn. obs_idx is [n_obs x 3] (i, j, t);
# y_obs holds the observed ranks. each (i, t) row's log likelihood is split
# evenly across its observed dyads
.ghk_rank_pointwise <- function(
	y_obs,
	ez_obs,
	obs_idx,
	n_mc = 64L,
	halton_shift = 0,
	n_bank = 8L,
	n_mc_fixed = FALSE
) {
	n_obs <- length(y_obs)
	out <- numeric(n_obs)
	if (n_obs == 0L) return(out)
	# bank of Halton sequences with distinct coprime bases; pooling their row
	# probability estimates by log-mean-exp cuts the Monte Carlo variance, and the
	# bank stays a deterministic function of the fit seed
	bank_bases <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L)
	bank_bases <- bank_bases[seq_len(min(n_bank, length(bank_bases)))]
	# group by (row i, period t)
	row_key <- paste(obs_idx[, 1], obs_idx[, 3], sep = "_")
	row_groups <- split(seq_len(n_obs), row_key)
	for (grp in row_groups) {
		yg <- y_obs[grp]
		ezg <- ez_obs[grp]
		# rank permutation, highest rank first: higher rank means larger z; ties ordered by ez
		perm <- order(-yg, -ezg, na.last = NA)
		# drop NA-rank dyads (frn unobserved past odmax)
		valid <- which(!is.na(yg[perm]))
		perm_valid <- perm[valid]
		if (length(perm_valid) <= 1L) {
			out[grp] <- 0
			next
		}
		# the log of a Monte Carlo probability has a downward (Jensen) bias that grows
		# with row length, so the sample count scales as 8 * 2^D (at least n_mc, at
		# most 2048); n_mc_fixed = TRUE (prior$ghk_n_mc in lame()) uses n_mc exactly
		D_row <- length(perm_valid)
		n_mc_row <- if (isTRUE(n_mc_fixed)) {
			as.integer(n_mc)
		} else {
			as.integer(min(
				2048L,
				max(
					as.integer(n_mc),
					8L * 2L^D_row
				)
			))
		}
		lp_bank <- vapply(
			bank_bases,
			function(hb)
				.ghk_row_loglik(
					ezg[perm_valid],
					seq_along(perm_valid),
					n_mc = n_mc_row,
					halton_shift = halton_shift,
					halton_base = hb
				),
			numeric(1)
		)
		mx <- max(lp_bank)
		row_lp <- if (is.finite(mx)) mx + log(mean(exp(lp_bank - mx))) else -1e10
		# split evenly across the row's observed dyads (keeps the row total, gives a
		# per-dyad value for loo)
		out[grp[perm_valid]] <- row_lp / length(perm_valid)
	}
	out
}

#' @noRd
# Gauss-Hermite alternative for cbin: a 1-d integral over the implicit row
# intercept, with statmod::gauss.quad.prob nodes when available
.cbin_gauss_hermite_loglik <- function(yg, ezg, n_nodes = 10L) {
	if (requireNamespace("statmod", quietly = TRUE)) {
		gh <- statmod::gauss.quad.prob(n_nodes, dist = "normal")
		nodes <- gh$nodes; wts <- gh$weights
	} else {
		# fallback: evenly spaced quantiles of N(0, 1)
		u <- seq(0.025, 0.975, length.out = n_nodes)
		nodes <- stats::qnorm(u)
		wts <- rep(1 / n_nodes, n_nodes)
	}
	# integrate over an additive row intercept a; the row-conditional log
	# likelihood is sum_j [y_j log Phi(ez_j + a) + (1 - y_j) log Phi(-(ez_j + a))]
	lp_nodes <- numeric(n_nodes)
	for (k in seq_len(n_nodes)) {
		a <- nodes[k]
		p_j <- stats::pnorm(ezg + a)
		p_j <- pmin(pmax(p_j, 1e-12), 1 - 1e-12)
		lp_nodes[k] <- sum(yg * log(p_j) + (1 - yg) * log(1 - p_j))
	}
	mx <- max(lp_nodes)
	mx + log(sum(wts * exp(lp_nodes - mx)))
}

#' @noRd
# GHK dispatcher for the rank families; closed-form families delegate to
# .pointwise_loglik_observed_ghk
.pointwise_loglik_observed_ghk_rank <- function(
	y_obs,
	ez_obs,
	z_obs,
	obs_idx,
	family,
	s2,
	dyad_rho = 0,
	n_mc = 64L,
	halton_shift = 0,
	n_mc_fixed = FALSE
) {
	if (family %in% c("normal", "binary", "poisson")) {
		return(.pointwise_loglik_observed_ghk(
			y_obs,
			ez_obs,
			z_obs,
			obs_idx,
			family,
			s2,
			dyad_rho,
			n_mc
		))
	}
	if (family == "ordinal") {
		return(.ordinal_pointwise_loglik(y_obs, ez_obs))
	}
	# cbin: Gauss-Hermite over the row intercept
	if (family == "cbin") {
		row_key <- paste(obs_idx[, 1], obs_idx[, 3], sep = "_")
		row_groups <- split(seq_len(length(y_obs)), row_key)
		out <- numeric(length(y_obs))
		for (grp in row_groups) {
			yg <- y_obs[grp]; ezg <- ez_obs[grp]
			row_lp <- .cbin_gauss_hermite_loglik(yg, ezg, n_nodes = 10L)
			out[grp] <- row_lp / length(grp)
		}
		return(out)
	}
	# frn: GHK on the implied row ordering
	if (family == "frn") {
		return(.ghk_rank_pointwise(
			y_obs,
			ez_obs,
			obs_idx,
			n_mc = n_mc,
			halton_shift = halton_shift,
			n_mc_fixed = n_mc_fixed
		))
	}
	# any other family: normal density of the latent value
	stats::dnorm(z_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
}
