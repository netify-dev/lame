# Z-marginalised Metropolis-Hastings update of explicit ordinal cutpoints
# (ordinal_cutpoints = "explicit"), stored in fit$alpha with alpha_1 = 0 as the anchor;
# Cowles (1996), Statistics and Computing 6, 101-111

#' Numerically stable log(Phi(hi) - Phi(lo))
#'
#' Computes log of the standard-normal CDF difference. Naive
#' \code{log(pnorm(hi) - pnorm(lo))} returns \code{-Inf} whenever both
#' arguments are above ~7 or below ~-7 (catastrophic cancellation).
#' Uses \code{pnorm(., log.p = TRUE)} and a stable log-subtraction.
#'
#' @param hi numeric, upper bound
#' @param lo numeric, lower bound (must be <= hi elementwise)
#' @param width optional original interval widths, before centering the bounds.
#'   Retains legitimate short intervals whose centered endpoints round together.
#' @return numeric, log(Phi(hi) - Phi(lo))
#' @keywords internal
log_phi_diff <- function(hi, lo, width = NULL) {
	if (length(hi) != length(lo) || anyNA(hi) || anyNA(lo))
		stop("Normal interval bounds must have matching lengths and no missing values.")
	if (is.null(width)) width <- hi - lo
	if (length(width) != length(hi) || anyNA(width))
		stop("Normal interval widths must match their bounds.")
	out <- rep(-Inf, length(hi))
	mid <- hi / 2 + lo / 2
	lower <- which(is.finite(hi) & is.finite(lo) & hi > lo & mid <= 0)
	upper <- which(is.finite(hi) & is.finite(lo) & hi > lo & mid > 0)
	if (length(lower)) {
		a <- stats::pnorm(hi[lower], log.p = TRUE)
		b <- stats::pnorm(lo[lower], log.p = TRUE)
		out[lower] <- a + log(-expm1(b - a))
	}
	if (length(upper)) {
		a <- stats::pnorm(hi[upper], lower.tail = FALSE, log.p = TRUE)
		b <- stats::pnorm(lo[upper], lower.tail = FALSE, log.p = TRUE)
		out[upper] <- b + log(-expm1(a - b))
	}
	inf_lo <- which(lo == -Inf & is.finite(hi))
	inf_hi <- which(is.finite(lo) & hi == Inf)
	out[inf_lo] <- stats::pnorm(hi[inf_lo], log.p = TRUE)
	out[inf_hi] <- stats::pnorm(lo[inf_hi], lower.tail = FALSE, log.p = TRUE)
	out[lo == -Inf & hi == Inf] <- 0
	# a short interval whose shifted endpoints round together is integrated about its midpoint:
	# P = width * phi(mid) * integral[-1/2,1/2] exp(-mid*width*u - (width*u)^2/2) du,
	# restricted to a smooth, bounded integrand
	short <- which(
		is.finite(mid) & is.finite(width) & width > 0 &
			width <= 1e-3 & abs(mid * width) <= 1e-2
	)
	for (k in short) {
		w <- width[k]; m <- mid[k]
		v <- stats::integrate(
			function(u) exp(-m * w * u - (w * u)^2 / 2),
			-.5,
			.5,
			rel.tol = 1e-12,
			abs.tol = 0,
			subdivisions = 100L
		)$value
		out[k] <- stats::dnorm(m, log = TRUE) + log(w) + log(v)
	}
	unresolved <- is.finite(hi) & is.finite(lo) & is.finite(width) &
		width > 0 & !is.finite(out)
	if (any(unresolved))
		stop("Positive normal interval probability cannot be resolved numerically.")
	out
}

# alpha <-> delta reparameterisation
.alpha_to_delta <- function(alpha) {
	# alpha has length K - 1 with alpha[1] = 0; delta has length K - 2
	if (length(alpha) <= 1L) return(numeric(0))
	log(diff(alpha))
}

.delta_to_alpha <- function(delta) {
	# length K - 1, anchored at alpha[1] = 0
	c(0, cumsum(exp(delta)))
}

# initial alpha from the empirical CDF (probit inverse), anchored at alpha[1] = 0
.init_alpha_from_data <- function(Y) {
	obs <- Y[is.finite(Y)]
	lvls <- sort(unique(obs))
	K <- length(lvls)
	if (K < 2L) {
		cli::cli_abort(c(
			"ordinal_cutpoints = 'explicit' requires at least 2 distinct observed categories.",
			"i" = "Found {K}; check {.arg Y} or use a different family."
		))
	}
	# warn about sparse categories, whose cutpoints are weakly identified
	counts <- tabulate(match(obs, lvls))
	if (any(counts < 5L) && K >= 3L) {
		cli::cli_warn(c(
			"ordinal_cutpoints = 'explicit': category counts {.val {counts}} have at least one cell below 5.",
			"i" = "A category with fewer than 5 observations leaves its cutpoint weakly identified; collapsing categories removes it."
		))
	}
	# probit inverse of cumulative proportions, shifted so alpha[1] = 0
	props <- cumsum(counts[-K]) / length(obs)
	props <- pmin(pmax(props, 1e-3), 1 - 1e-3)
	alpha_raw <- stats::qnorm(props)
	alpha_raw - alpha_raw[1L]
}

# Z-marginalised ordinal log-likelihood at cutpoints alpha,
# sum log{pnorm(alpha_y - eta) - pnorm(alpha_{y-1} - eta)}; `alpha_full` is
# c(-Inf, alpha, Inf), and the symmetric caller already counts each dyad once
.ord_marginal_loglik <- function(alpha_full, eta_obs, Y_int_obs) {
	hi <- alpha_full[Y_int_obs + 1L]
	lo <- alpha_full[Y_int_obs]
	log_p <- log_phi_diff(
		hi - eta_obs,
		lo - eta_obs,
		width = hi - lo
	)
	sum(log_p)
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
		# K = 2 (binary): no free cutpoint to update
		return(list(alpha = alpha, delta = numeric(0), accept = FALSE))
	}

	# mask: observed cells; for symmetric, restrict to upper triangle so the
	# information per cutpoint is not double-counted
	obs_mask <- !is.na(Y_int)
	if (isTRUE(symmetric)) {
		ut <- upper.tri(matrix(FALSE, dim(Y_int)[1L], dim(Y_int)[2L]))
		if (length(dim(Y_int)) == 3L) ut <- array(ut, dim(Y_int))
		obs_mask <- obs_mask & ut
	}
	if (!any(obs_mask)) {
		return(list(alpha = alpha, delta = .alpha_to_delta(alpha), accept = FALSE))
	}
	Y_obs <- Y_int[obs_mask]
	EZ_obs <- if (isTRUE(symmetric)) {
		EZ_mirror <- if (length(dim(EZ)) == 3L) aperm(EZ, c(2L, 1L, 3L)) else t(EZ)
		(EZ[obs_mask] + EZ_mirror[obs_mask]) / 2
	} else EZ[obs_mask]

	delta <- .alpha_to_delta(alpha)
	delta_prop <- delta + tau_prop * stats::rnorm(length(delta))
	alpha_prop <- .delta_to_alpha(delta_prop)

	af_curr <- c(-Inf, alpha, Inf)
	af_prop <- c(-Inf, alpha_prop, Inf)

	ll_curr <- .ord_marginal_loglik(af_curr, EZ_obs, Y_obs)
	ll_prop <- .ord_marginal_loglik(af_prop, EZ_obs, Y_obs)

	# log Jacobian of delta -> alpha; the prior on ordered alpha is flat and the proposal
	# is symmetric in delta
	log_jac <- sum(delta_prop) - sum(delta)
	log_a <- (ll_prop - ll_curr) + log_jac
	# a non-finite ratio means the support is violated
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
#' @param rho directed within-dyad correlation; zero for independent cells.
#' @return updated Z
#' @keywords internal
rZ_ord_explicit_fc <- function(Z, EZ, Y_int, alpha, rho = 0) {
	if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
		stop("Explicit ordinal correlation must be finite and inside (-1, 1).")
	if (rho != 0) return(.rZ_ord_explicit_correlated_fc(Z, EZ, Y_int, alpha, rho))
	af <- c(-Inf, alpha, Inf)
	obs <- !is.na(Y_int)
	if (any(!is.finite(EZ)) || any(!is.finite(alpha)) ||
			any(diff(alpha) <= 0) || any(Y_int[obs] != as.integer(Y_int[obs])) ||
			any(Y_int[obs] < 1L | Y_int[obs] > length(alpha) + 1L)) {
		stop("Invalid explicit ordinal latent-state inputs.")
	}
	if (any(obs)) {
		lo <- af[Y_int[obs]]; hi <- af[Y_int[obs] + 1L]
		draw <- rtnorm_interval_logp(EZ[obs], lo, hi)
		if (any(!is.finite(draw)) || any(draw <= lo | draw >= hi)) {
			stop("Explicit ordinal truncated-normal draw is numerically unresolved.")
		}
		Z[obs] <- draw
	}
	# missing cells, including nuisance diagonals, have no rank constraint
	if (any(!obs)) Z[!obs] <- stats::rnorm(sum(!obs), EZ[!obs])
	if (any(!is.finite(Z))) stop("Explicit ordinal latent draw is outside numerical range.")
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
	obs_input <- !is.na(Y_int)
	if (any(!is.finite(EZ)) || any(!is.finite(alpha)) ||
			any(diff(alpha) <= 0) || any(Y_int[obs_input] != as.integer(Y_int[obs_input])) ||
			any(Y_int[obs_input] < 1L | Y_int[obs_input] > length(alpha) + 1L)) {
		stop("Invalid explicit ordinal latent-state inputs.")
	}
	ut <- upper.tri(Z)
	af <- c(-Inf, alpha, Inf)
	sd_sym <- 1
	# observed upper-triangle cells
	obs <- ut & !is.na(Y_int)
	if (any(obs)) {
		Y_obs <- Y_int[obs]
		ez_avg <- (EZ[obs] + t(EZ)[obs]) / 2
		lo <- af[Y_obs]
		hi <- af[Y_obs + 1L]
		z_new <- sd_sym * rtnorm_interval_logp(ez_avg / sd_sym, lo / sd_sym, hi / sd_sym)
		if (any(!is.finite(z_new)) || any(z_new <= lo | z_new >= hi)) {
			stop("Symmetric explicit ordinal truncated-normal draw is numerically unresolved.")
		}
		Z[obs] <- z_new
		# mirror to lower triangle
		mirror_idx <- which(obs, arr.ind = TRUE)[, c(2L, 1L), drop = FALSE]
		Z[mirror_idx] <- z_new
	}
	# missing unordered dyads carry no category constraint
	missing <- ut & is.na(Y_int)
	if (any(missing)) {
		Z[missing] <- stats::rnorm(sum(missing), (EZ[missing] + t(EZ)[missing]) / 2)
		mirror <- which(missing, arr.ind = TRUE)[, c(2L, 1L), drop = FALSE]
		Z[mirror] <- Z[missing]
	}
	# diagonal: unconstrained prior draw
	diag(Z) <- stats::rnorm(n, diag(EZ), 1)
	if (any(!is.finite(Z))) stop("Symmetric explicit ordinal latent draw is outside numerical range.")
	Z
}
