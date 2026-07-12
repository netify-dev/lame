# cowles (1996) z-marginalised metropolis-hastings update for explicit
# ordinal cutpoints, an alternative to the
# data-induced cutpoint convention. default behaviour unchanged
# (`ordinal_cutpoints = "data_induced"`); the explicit path is opt-in
# via `ordinal_cutpoints = "explicit"` and stores a posterior on
# `fit$alpha` with `alpha_1 = 0` as the identification anchor.
#
# reference: cowles, m. k. (1996). accelerating monte carlo markov chain
# convergence for cumulative-link generalized linear models. statistics
# and computing, 6, 101--111.

#' Numerically stable log(Phi(hi) - Phi(lo))
#'
#' Computes log of the standard-normal CDF difference. Naive
#' \code{log(pnorm(hi) - pnorm(lo))} returns \code{-Inf} whenever both
#' arguments are above ~7 or below ~-7 (catastrophic cancellation).
#' Uses \code{pnorm(., log.p = TRUE)} and a stable log-subtraction.
#'
#' @param hi numeric, upper bound
#' @param lo numeric, lower bound (must be <= hi elementwise)
#' @return numeric, log(Phi(hi) - Phi(lo))
#' @keywords internal
log_phi_diff <- function(hi, lo) {
	# stable in both tails. for points in the lower tail we use the
	# lower-tail log cdf and compute log(exp(a) - exp(b)) where a > b.
	# for the upper tail we flip to the upper-tail log cdf (1 - phi) and
	# compute log(exp(b') - exp(a')) where b' > a' because the survival
	# function is monotone-decreasing. we dispatch elementwise off the
	# midpoint sign so an array with both regimes is handled correctly.
	out <- rep(-Inf, length(hi))
	mid <- (hi + lo) / 2
	lower <- which(is.finite(hi) & is.finite(lo) & mid <= 0)
	upper <- which(is.finite(hi) & is.finite(lo) & mid > 0)
	# infinite-bound case: lo = -inf or hi = +inf collapses to a single cdf
	inf_lo <- which(!is.finite(lo) & is.finite(hi))
	inf_hi <- which(is.finite(lo) & !is.finite(hi))
	if (length(lower)) {
		a <- stats::pnorm(hi[lower], log.p = TRUE)
		b <- stats::pnorm(lo[lower], log.p = TRUE)
		out[lower] <- a + log1p(-exp(b - a))
	}
	if (length(upper)) {
		# survival-function form: p(lo < x <= hi) = s(lo) - s(hi),
		# s(.) = pnorm(., lower.tail = false). both s values are
		# representable as log-tail probs without cancellation.
		ap <- stats::pnorm(hi[upper], lower.tail = FALSE, log.p = TRUE)
		bp <- stats::pnorm(lo[upper], lower.tail = FALSE, log.p = TRUE)
		# bp > ap because s is decreasing and lo < hi
		out[upper] <- bp + log1p(-exp(ap - bp))
	}
	if (length(inf_lo)) {
		out[inf_lo] <- stats::pnorm(hi[inf_lo], log.p = TRUE)
	}
	if (length(inf_hi)) {
		out[inf_hi] <- stats::pnorm(lo[inf_hi], lower.tail = FALSE, log.p = TRUE)
	}
	out[!is.finite(out)] <- -Inf
	out
}

# alpha <-> delta reparameterisation
.alpha_to_delta <- function(alpha) {
	# alpha is length (k-1) with alpha[1] = 0; delta is length (k-2)
	if (length(alpha) <= 1L) return(numeric(0))
	log(diff(alpha))
}

.delta_to_alpha <- function(delta) {
	# returns length (k-1), anchored at alpha[1] = 0
	c(0, cumsum(exp(delta)))
}

# initial alpha from empirical cdf (probit inverse), anchored at alpha[1] = 0
.init_alpha_from_data <- function(Y) {
	obs <- Y[is.finite(Y)]
	lvls <- sort(unique(obs))
	K <- length(lvls)
	if (K < 2L) {
		cli::cli_abort(c(
			"ordinal_cutpoints = 'explicit' requires at least 2 distinct observed categories.",
			"i" = "Found {K}; check {.arg Y} or use a different family."))
	}
	# guard against empty interior class (cutpoint posterior improper)
	counts <- tabulate(match(obs, lvls))
	if (any(counts < 5L) && K >= 3L) {
		cli::cli_warn(c(
			"ordinal_cutpoints = 'explicit': category counts {.val {counts}} have at least one cell below 5.",
			"i" = "Cutpoint posterior may be weakly identified. Consider collapsing categories."))
	}
	# probit inverse of cumulative proportions, shifted so alpha[1] = 0
	props <- cumsum(counts[-K]) / length(obs)
	props <- pmin(pmax(props, 1e-3), 1 - 1e-3)
	alpha_raw <- stats::qnorm(props)
	alpha_raw - alpha_raw[1L]
}

# z-marginalised ordinal log-likelihood at cutpoints alpha:
# sum_{ij} log{pnorm(alpha_y - eta) - pnorm(alpha_y_minus_1 - eta)}
# under z ~ N(eta, 1) the marginal category probability is
# pnorm(alpha_k - eta) - pnorm(alpha_{k-1} - eta) for both the symmetric
# and asymmetric models (variance 1). the symmetric path counts each
# dyad once via the upper-triangle mask in sample_alpha_cowles, so no
# variance rescaling is applied.
# `y_int` is the integer-recoded y in 1..k (na for missing); `alpha_full`
# is c(-inf, alpha, inf) of length k+1.
.ord_marginal_loglik <- function(alpha_full, eta_obs, Y_int_obs, symmetric) {
	scale <- 1
	hi <- alpha_full[Y_int_obs + 1L]
	lo <- alpha_full[Y_int_obs]
	log_p <- log_phi_diff(scale * (hi - eta_obs), scale * (lo - eta_obs))
	sum(log_p[is.finite(log_p)])
}

#' Cowles MH update for explicit ordinal cutpoints (Z-marginalised)
#'
#' One Metropolis-Hastings sweep on \code{alpha} with Z integrated out
#' of the conditional. The proposal is a symmetric random walk in
#' \code{delta = log(diff(alpha))} space; the acceptance ratio uses
#' the marginal ordinal-probit likelihood. Mixing of \code{alpha}
#' under this update is dramatically faster than the data-induced
#' Gibbs convention (Cowles 1996 reports ESS gains of an order of
#' magnitude on standard ordinal-probit benchmarks).
#'
#' @param alpha numeric, current cutpoints with \code{alpha[1] = 0},
#'   length \code{K-1}.
#' @param Y_int integer matrix (or array) of ordinal categories coded
#'   1..K, with \code{NA} for missing cells.
#' @param EZ linear predictor at the dyad level, same shape as
#'   \code{Y_int}.
#' @param tau_prop numeric, RW proposal SD in delta space (will be
#'   adapted in burn-in by the caller).
#' @param symmetric logical, whether the network is symmetric (restricts
#'   the likelihood to upper-triangle cells so each dyad is counted once;
#'   variance 1, no precision rescaling).
#' @return list with elements \code{alpha}, \code{delta}, \code{accept}
#'   (logical).
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
sample_alpha_cowles <- function(alpha, Y_int, EZ, tau_prop, symmetric = FALSE) {
	K_minus_1 <- length(alpha)
	if (K_minus_1 < 2L) {
		# k = 2 (binary): no free cutpoint to update
		return(list(alpha = alpha, delta = numeric(0), accept = FALSE))
	}

	# mask: observed cells; for symmetric, restrict to upper triangle so the
	# information per cutpoint is not double-counted
	obs_mask <- !is.na(Y_int)
	if (isTRUE(symmetric)) {
		ut <- upper.tri(Y_int)
		if (length(dim(Y_int)) == 2L) obs_mask <- obs_mask & ut
	}
	if (!any(obs_mask)) {
		return(list(alpha = alpha, delta = .alpha_to_delta(alpha), accept = FALSE))
	}
	Y_obs <- Y_int[obs_mask]
	EZ_obs <- EZ[obs_mask]

	delta <- .alpha_to_delta(alpha)
	delta_prop <- delta + tau_prop * stats::rnorm(length(delta))
	alpha_prop <- .delta_to_alpha(delta_prop)

	K <- K_minus_1 + 1L  # number of categories
	af_curr <- c(-Inf, alpha,      Inf)
	af_prop <- c(-Inf, alpha_prop, Inf)

	ll_curr <- .ord_marginal_loglik(af_curr, EZ_obs, Y_obs, symmetric)
	ll_prop <- .ord_marginal_loglik(af_prop, EZ_obs, Y_obs, symmetric)

	# jacobian of delta -> alpha is product of exp(delta_k), so log |j*|/|j|
	# = sum(delta* - delta). flat prior on ordered alpha; the proposal
	# itself is symmetric in delta.
	log_jac <- sum(delta_prop) - sum(delta)
	log_a <- (ll_prop - ll_curr) + log_jac
	# numerical sanity: nan means support violated
	if (!is.finite(log_a)) {
		return(list(alpha = alpha, delta = delta, accept = FALSE))
	}
	if (log(stats::runif(1)) < log_a) {
		list(alpha = alpha_prop, delta = delta_prop, accept = TRUE)
	} else {
		list(alpha = alpha, delta = delta, accept = FALSE)
	}
}

#' Sample Z given explicit cutpoints alpha (asymmetric ordinal)
#'
#' Rectangular / square truncated-normal sample on observed cells.
#' This is the explicit-cutpoint version of \code{\link{rZ_ord_fc}},
#' parameterised by the explicit \code{alpha} vector rather than
#' deriving cutpoints from Z order statistics.
#'
#' @param Z current latent matrix
#' @param EZ expected value of Z (linear predictor)
#' @param Y_int integer-recoded ordinal Y (1..K)
#' @param alpha length-(K-1) cutpoints with \code{alpha[1] = 0}
#' @return updated Z
#' @keywords internal
rZ_ord_explicit_fc <- function(Z, EZ, Y_int, alpha) {
	af <- c(-Inf, alpha, Inf)
	obs <- !is.na(Y_int)
	if (!any(obs)) return(Z)
	# truncated normal per cell
	lo <- af[Y_int[obs]]
	hi <- af[Y_int[obs] + 1L]
	ez_obs <- EZ[obs]
	u <- stats::runif(sum(obs))
	# inverse-cdf on the truncation interval
	plo <- stats::pnorm(lo - ez_obs)
	phi <- stats::pnorm(hi - ez_obs)
	# guard against degenerate (plo == phi) by widening minimally
	deg <- (phi - plo) < 1e-12
	if (any(deg)) {
		# midpoint fallback at the cell mean
		mid <- (lo[deg] + hi[deg]) / 2
		mid[!is.finite(mid)] <- ez_obs[deg][!is.finite(mid)]
		Z[obs][deg] <- mid
	}
	if (any(!deg)) {
		z_new <- ez_obs[!deg] + stats::qnorm(
			plo[!deg] + u[!deg] * (phi[!deg] - plo[!deg]))
		Z[obs][!deg] <- z_new
	}
	Z
}

#' Symmetric ordinal Z sample given explicit cutpoints (R >= 0 safe)
#'
#' Mirror image of \code{\link{rZ_ord_sym_fc}} (variance 1 on the upper
#' triangle, mirrored to the lower) parameterised by an explicit
#' \code{alpha} vector instead of data-induced cutpoints.
#'
#' @param Z current symmetric latent matrix
#' @param EZ expected value of Z
#' @param Y_int integer ordinal Y (symmetric)
#' @param alpha length-(K-1) cutpoints with \code{alpha[1] = 0}
#' @return updated symmetric Z (diagonal sampled from unconstrained prior)
#' @keywords internal
rZ_ord_sym_explicit_fc <- function(Z, EZ, Y_int, alpha) {
	n <- nrow(Y_int)
	ut <- upper.tri(Z)
	af <- c(-Inf, alpha, Inf)
	sd_sym <- 1
	# upper triangle observed cells
	obs <- ut & !is.na(Y_int)
	if (any(obs)) {
		Y_obs <- Y_int[obs]
		ez_avg <- (EZ[obs] + t(EZ)[obs]) / 2
		lo <- af[Y_obs]
		hi <- af[Y_obs + 1L]
		plo <- stats::pnorm((lo - ez_avg) / sd_sym)
		phi <- stats::pnorm((hi - ez_avg) / sd_sym)
		u <- stats::runif(length(ez_avg))
		# stable truncation interval sample
		z_new <- ez_avg + sd_sym * stats::qnorm(plo + u * (phi - plo))
		# guard degenerate
		bad <- !is.finite(z_new)
		if (any(bad)) z_new[bad] <- (lo[bad] + hi[bad]) / 2
		Z[obs] <- z_new
		# mirror to lower triangle
		mirror_idx <- which(obs, arr.ind = TRUE)[, c(2L, 1L), drop = FALSE]
		Z[mirror_idx] <- z_new
	}
	# diagonal: unconstrained prior draw (matches existing convention)
	diag(Z) <- stats::rnorm(n, diag(EZ), 1)
	Z
}
