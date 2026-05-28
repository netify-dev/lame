# symmetric ordinal Z sampler. Z_ij = Z_ji throughout; the upper triangle
# is sampled from a truncated normal with the symmetric-doubled precision,
# then mirrored to the lower triangle.
#
# Math. Under the symmetric ordinal model Y_ij = Y_ji = k with
# Z_ij = Z_ji constrained equal across the directed dyad, the joint
# conditional kernel (before truncation) is
#   phi(z - eta_ij) * phi(z - eta_ji) propto exp(-(z - mu_eff)^2 / (2 * 0.5))
# where mu_eff = (eta_ij + eta_ji) / 2 and the variance is 0.5 (precision
# doubled relative to the directed conditional). Ordinal truncation
# restricts z to (alpha_{k-1}, alpha_k].
#
# note on the empirical SE ratio: a naive "symmetric posterior SD is
# sqrt(0.5) ~ 0.71 of the asymmetric posterior SD" prediction does not
# hold cleanly in practice (the measured ratio is closer to ~0.85), and
# the asymmetric posterior on symmetric data is also shifted in mean, so
# the two posteriors are not simple rescalings of each other. both
# samplers are valid; the symmetric path is more
# efficient than the asymmetric one on truly symmetric data, but the
# precise efficiency gain is data- and prior-dependent.
#
# the cutpoint convention matches the rest of the package's ordinal path:
# cutpoints are induced by the order statistics of Z within neighbouring
# rank classes, not separately parameterised. Cowles MH on explicit alpha
# is a future addition.

#' Sample Z under symmetric ordinal data (square symmetric matrix)
#'
#' Symmetric-Z drop-in for \code{\link{rZ_ord_fc}}: each unordered dyad
#' {i, j} carries a single latent value (Z_ij = Z_ji) drawn from a
#' truncated normal with mean equal to the average linear predictor
#' (EZ_ij + EZ_ji) / 2 and variance 1/2 (precision 2). The lower
#' triangle is mirrored after the upper-triangle draws so the full
#' matrix is exactly symmetric at every sweep.
#'
#' @param Z current symmetric latent matrix (n x n, Z = t(Z)).
#' @param EZ expected value of Z (regression + additive + multiplicative
#'   contributions); can be asymmetric, the helper averages the two
#'   directed contributions.
#' @param Y observed symmetric ordinal matrix (n x n, Y = t(Y));
#'   diagonal is ignored.
#' @return updated symmetric n x n latent matrix with \code{NA} on the
#'   diagonal.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_ord_sym_fc <- function(Z, EZ, Y) {
	n <- nrow(Y)
	uY <- sort(unique(c(Y)))
	W <- matrix(match(Y, uY), n, n)
	W[is.na(W)] <- -1L

	ut <- upper.tri(Z)
	sd_sym <- sqrt(0.5)

	for (w in sample(c(-1L, seq_along(uY)))) {
		lb <- suppressWarnings(max(Z[!is.na(W) & W == w - 1L], na.rm = TRUE))
		ub <- suppressWarnings(min(Z[!is.na(W) & W == w + 1L], na.rm = TRUE))
		up <- ut & W == w
		if (!any(up)) next
		# symmetric conditional mean: average of (i,j) and (j,i) linear
		# predictors. precision 2 -> sd sqrt(1/2).
		ez_avg <- (EZ[up] + t(EZ)[up]) / 2
		z_new <- ez_avg + sd_sym * stats::qnorm(stats::runif(sum(up),
		                                                       stats::pnorm((lb - ez_avg) / sd_sym),
		                                                       stats::pnorm((ub - ez_avg) / sd_sym)))
		Z[up] <- z_new
		# mirror to (j,i) lower-triangle cells. up has the (i,j) indices in
		# column-major order; the mirror is at (j,i), so swap arr.ind columns
		mirror_idx <- which(up, arr.ind = TRUE)[, c(2L, 1L), drop = FALSE]
		Z[mirror_idx] <- z_new
	}
	# diagonal: sample from the unconstrained prior (matches the convention
	# in rZ_ord_fc / rZ_bin_fc). Setting `NA` here propagates into the
	# symmetrised residual `E <- 0.5*(E + t(E))` consumed by rUV_sym_fc_cpp,
	# which then errors at R >= 1 with "Mat::min(): object has no elements"
	# because every column of the residual has at least one NA.
	diag(Z) <- stats::rnorm(nrow(Z), diag(EZ), 1)
	Z
}
