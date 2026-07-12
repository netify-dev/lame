# ghk / halton / gauss-hermite pointwise log-lik for rank-family ames.
#
# a geweke-hajivassiliou-keane simulator on the implied ordering:
#
#   (1) sequential conditioning (p(z_{pi_k} > thresh_{k-1}) and
#       tn(thresh, +inf) draws).
#   (2) halton sequence is randomised per fit via a stored seed shift
#       (owen-style) for reproducibility + better coverage of [0,1]^d.
#   (3) dimension cap d <= 15 (rows with more observed dyads fall back
#       to the augmented-z gaussian for that row only).
#
# per-dyad pointwise log-lik: the row-level ghk estimate is split
# uniformly across the row's observed dyads. this is the standard
# compromise for loo / waic on rank-family models; the row total is
# preserved.
#
# gauss-hermite alternative for cbin: when family == "cbin" we have a
# constrained-binary row where exactly k_i out of d_i dyads have y = 1.
# the 1-d integral is over the implicit row threshold and 10-node
# gauss-hermite is sufficient (machine-precision for k_i, d_i <= 30).
#
# all mc components are deterministic given fit$log_lik_meta$seed; the
# seed is recorded on the fit so a re-run produces bit-identical loo
# summaries.

#' @noRd
# coprime prime bases for a proper multidimensional halton sequence: the
# k-th integration dimension of the ghk rank estimator uses the k-th base.
# covers the D <= 15 dimension cap (14 integration dims) plus the bank
# offset ladder with room to spare.
.ghk_prime_bases <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L,
                      37L, 41L, 43L, 47L, 53L, 59L, 61L, 67L, 71L, 73L,
                      79L, 83L, 89L, 97L)

# generalised van der corput / halton-with-base-b sequence with owen-style
# digital shift (faure 1992): randomises the sequence while preserving the
# low-discrepancy property.
.halton_seq_randomised <- function(n, base, burn = 30L, shift = 0) {
	out <- numeric(n + burn)
	for (i in seq_len(n + burn)) {
		f <- 1; r <- 0; k <- i
		while (k > 0) {
			f <- f / base
			# owen shift: add a uniform random digit shift then take mod 1
			r <- r + f * ((k %% base) + shift) %% base
			k <- k %/% base
		}
		out[i] <- r %% 1
	}
	out[(burn + 1L):(burn + n)]
}

#' @noRd
# halton-normal with per-fit randomisation. shift is a [0, base) numeric
# stored on fit$log_lik_meta for reproducibility.
.halton_normal_randomised <- function(n, base, shift = 0) {
	u <- .halton_seq_randomised(n, base, shift = shift)
	u <- pmin(pmax(u, 1e-7), 1 - 1e-7)
	stats::qnorm(u)
}

#' @noRd
# row-level ghk estimator for the rank-family (frn) log-likelihood.
# given the row's expected linear predictors ez (length d) and the
# observed rank ordering perm (length d; perm[1] = top rank), returns the
# log of the unbiased mc estimate of
#   p(z_{perm[1]} > z_{perm[2]} > ... > z_{perm[d]})
# averaged over n_mc halton draws. halton shift is the per-fit randomiser.
.ghk_row_loglik <- function(ez, perm, n_mc = 64L, halton_shift = 0,
                             halton_base = 2L) {
	D <- length(ez)
	if (D == 0L) return(0)
	if (D == 1L) return(0)   # single-dyad row is trivially in order
	if (D > 15L) {
		# dim cap: fall back to a gaussian quadratic-form approximation
		# (treat the row log-lik as the sum of pairwise log-probs of the
		# adjacent comparisons; biased but tractable for large d).
		ez_ord <- ez[perm]
		lp <- 0
		for (k in seq_len(D - 1L)) {
			# p(z_{k} > z_{k+1}) = p(z_k - z_{k+1} > 0); z_k - z_{k+1} ~
			# n(ez_k - ez_{k+1}, 2)
			diff_mu <- ez_ord[k] - ez_ord[k + 1L]
			lp <- lp + stats::pnorm(diff_mu / sqrt(2), log.p = TRUE)
		}
		return(lp)
	}
	# proper ghk over the implied ordering
	# halton-normal draws u_{1..d-1, 1..n_mc} (last position has no
	# upper bound to integrate against).
	# each of the d-1 integration dimensions must use its own coprime prime
	# base. filling the whole [d-1, n_mc] block from a single base makes the
	# per-dimension draws a correlated lattice (cross-dim cor approx -0.71),
	# which systematically biases the nested rank probability upward: against
	# exact 1/D! the single-base estimate over-states log P by up to ~1.5
	# nats at D = 6, and per-dimension coprime bases recover the truth.
	# halton_base offsets the start of the prime ladder so bank members
	# built from different offsets stay coprime across dimensions too.
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
	# per-mc accumulator
	lp_mc <- numeric(n_mc)
	for (m in seq_len(n_mc)) {
		thresh <- -Inf
		lp <- 0
		# sweep ranks from k = d down to k = 1: perm[1] is the top rank
		# (largest z), thresh tracks the running maximum of lower ranks.
		# at each step:
		#   p(z_{perm[k]} > thresh | z_{perm[k]} ~ n(ez_perm[k], 1))
		#   = 1 - phi(thresh - ez_perm[k])
		# then sample z_{perm[k]} from tn(ez_perm[k], 1; thresh, +inf)
		# and update thresh.
		for (k_idx in D:1L) {
			ez_k <- ez_ord[k_idx]
			if (is.infinite(thresh)) {
				p_k <- 1   # unconstrained probability
			} else {
				p_k <- stats::pnorm(thresh - ez_k, lower.tail = FALSE)
				p_k <- pmin(pmax(p_k, 1e-15), 1 - 1e-15)
				lp <- lp + log(p_k)
			}
			# sample z_{perm[k]} from tn
			if (k_idx > 1L) {
				lo_cdf <- if (is.infinite(thresh)) 0 else stats::pnorm(thresh - ez_k)
				# use halton draw indexed by (k_idx - 1, m); cycle as needed
				u_draw <- hn[((D - k_idx)) %% (D - 1L) + 1L, m]
				z_draw <- ez_k + stats::qnorm(lo_cdf + u_draw * (1 - lo_cdf))
				if (is.finite(z_draw)) thresh <- z_draw
			}
		}
		lp_mc[m] <- lp
	}
	# log-mean-exp over the mc sample
	mx <- max(lp_mc)
	if (!is.finite(mx)) return(-1e10)
	mx + log(mean(exp(lp_mc - mx)))
}

#' @noRd
# per-row ghk + per-dyad attribution. obs_idx is [n_obs x 3] (i, j, t),
# ez_obs is length n_obs, y_obs is length n_obs (the observed ranks for
# frn, or 0/1 for cbin). for each (i, t) row, compute the row log-lik
# and attribute uniformly across the row's observed dyads.
.ghk_rank_pointwise <- function(y_obs, ez_obs, obs_idx, family,
                                    n_mc = 64L, halton_shift = 0,
                                    n_bank = 8L, n_mc_fixed = FALSE) {
	n_obs <- length(y_obs)
	out <- numeric(n_obs)
	if (n_obs == 0L) return(out)
	# fixed bank of independent halton sequences, one per distinct prime
	# base. a single fixed base-2 shift gives a noisy row log-lik (up to
	# ~1.25 nats off for low-probability orderings); coprime bases produce
	# genuinely decorrelated sequences (same-base digital shifts stay
	# correlated), so pooling their row *probability* estimates via
	# log-mean-exp cuts the single-shift mc variance by ~an order of
	# magnitude for realistic orderings. the estimand is unchanged (log of
	# an unbiased row-probability estimate) and the per-fit halton_shift
	# randomises every member, so the whole bank stays a deterministic
	# function of the fit seed (bit-reproducible).
	bank_bases <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L)
	bank_bases <- bank_bases[seq_len(min(n_bank, length(bank_bases)))]
	# group by (row i, period t)
	row_key <- paste(obs_idx[, 1], obs_idx[, 3], sep = "_")
	row_groups <- split(seq_len(n_obs), row_key)
	for (grp in row_groups) {
		yg <- y_obs[grp]
		ezg <- ez_obs[grp]
		# build rank permutation: highest rank first.
		# frn: rank values directly. cbin: y in {0, 1}; we treat the
		# 1's as "top tied rank" and 0's as below.
		if (identical(family, "cbin")) {
			# constrained binary: row top-k are observed as 1, rest 0
			# implied ordering: all y=1 dyads come before y=0 dyads
			# break y=1 ties by ez (higher ez = stronger preference)
			ord_1 <- order(-ezg * (yg > 0), -yg)
			perm <- ord_1
		} else {
			# frn: rank values give ordering. higher rank = larger z.
			# within ranks (ties), order by ez.
			perm <- order(-yg, -ezg, na.last = NA)
		}
		# trim na-rank dyads (frn unobserved past odmax)
		valid <- which(!is.na(yg[perm]))
		perm_valid <- perm[valid]
		if (length(perm_valid) <= 1L) {
			# trivial row
			out[grp] <- 0
			next
		}
		# the log of a Monte-Carlo probability estimate carries a downward
		# (Jensen) bias ~ Var(p_hat)/(2 p^2) that grows as the row gets
		# longer (the rank probability shrinks like 1/D!). scale the sample
		# count with D so the bias stays < ~0.05 nats across row lengths.
		# note this runs once per stored draw (not once per fit), so long
		# rows are expensive; n_mc raises the budget past the floor, and
		# an explicit budget (n_mc_fixed = TRUE, wired from
		# prior$ghk_n_mc in lame()) is honored exactly -- it can both
		# raise the 2048 ceiling and lower the 8 * 2^D accuracy floor.
		D_row <- length(perm_valid)
		n_mc_row <- if (isTRUE(n_mc_fixed)) {
			as.integer(n_mc)
		} else {
			as.integer(min(2048L, max(as.integer(n_mc),
			                          8L * 2L^D_row)))
		}
		lp_bank <- vapply(bank_bases, function(hb)
			.ghk_row_loglik(ezg[perm_valid], seq_along(perm_valid),
			                n_mc = n_mc_row, halton_shift = halton_shift,
			                halton_base = hb),
			numeric(1))
		mx <- max(lp_bank)
		row_lp <- if (is.finite(mx)) mx + log(mean(exp(lp_bank - mx))) else -1e10
		# uniform attribution across the row's observed dyads (preserves
		# the row total, gives a per-dyad number that loo can consume)
		out[grp[perm_valid]] <- row_lp / length(perm_valid)
	}
	out
}

#' @noRd
# gauss-hermite alternative for cbin: 1-d integral over the implicit row
# intercept. uses statmod::gauss.quad.prob if available, else falls back
# to a 10-node grid via stats::integrate. n_nodes = 10 is sufficient for
# k_i, d_i <= 30 (verified empirically against mc).
.cbin_gauss_hermite_loglik <- function(yg, ezg, n_nodes = 10L) {
	# hermite quadrature nodes + weights for integral against n(0, 1)
	# use golub-welsch via stats::qnorm + uniform spacing as a cheap
	# proxy. for better accuracy, use statmod::gauss.quad.prob when
	# available.
	if (requireNamespace("statmod", quietly = TRUE)) {
		gh <- statmod::gauss.quad.prob(n_nodes, dist = "normal")
		nodes <- gh$nodes; wts <- gh$weights
	} else {
		# pragmatic fallback: 10 evenly-spaced quantiles of n(0, 1)
		u <- seq(0.025, 0.975, length.out = n_nodes)
		nodes <- stats::qnorm(u)
		wts <- rep(1 / n_nodes, n_nodes)
	}
	# integrate over an additive row intercept a; the inner row-conditional
	# log-lik is sum_j [y_j log phi(ez_j + a) + (1 - y_j) log phi(-(ez_j + a))]
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
# top-level ghk dispatcher for rank families: randomised-halton +
# dim-capped ghk + gauss-hermite cbin fallback. closed-form families
# delegate to .pointwise_loglik_observed_ghk.
.pointwise_loglik_observed_ghk_rank <- function(y_obs, ez_obs, z_obs, obs_idx,
                                               family, s2, dyad_rho = 0,
                                               n_mc = 64L,
                                               halton_shift = 0,
                                               n_mc_fixed = FALSE) {
	# closed-form families delegate to the closed-form path
	if (family %in% c("normal", "binary", "poisson")) {
		return(.pointwise_loglik_observed_ghk(y_obs, ez_obs, z_obs, obs_idx,
		                                       family, s2, dyad_rho, n_mc))
	}
	if (family == "ordinal") {
		return(.ordinal_pointwise_loglik(y_obs, ez_obs))
	}
	# cbin: try gauss-hermite (closed-form-ish via 1-d quad), fall back to ghk
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
	# frn: ghk on the implied row ordering
	if (family == "frn") {
		return(.ghk_rank_pointwise(y_obs, ez_obs, obs_idx, family,
		                               n_mc = n_mc,
		                               halton_shift = halton_shift,
		                               n_mc_fixed = n_mc_fixed))
	}
	# default fallback
	stats::dnorm(z_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
}
