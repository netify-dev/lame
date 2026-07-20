# symmetric ordinal z sampler. z_ij = z_ji throughout; the upper triangle
# is sampled from a truncated normal (variance 1), then mirrored to the
# lower triangle.
#
# the symmetric ordinal model carries one latent value per
# unordered dyad {i, j}: z_ij = z_ji ~ N(mu_ij, 1) truncated to the
# cutpoint interval (alpha_{k-1}, alpha_k], with mu_ij the symmetric
# linear predictor. under the symmetric fit eta_ij = eta_ji, so the
# conditional mean is mu_eff = (eta_ij + eta_ji) / 2 and the variance
# is 1, the model's conditional. this matches the variance-1 convention
# used by the asymmetric sampler, by rbeta_ab_fc's rho->1 GLS (which
# assumes s2 = 1), and by every other family.
#
# the cutpoint convention matches the rest of the package's ordinal path:
# cutpoints are induced by the order statistics of z within neighbouring
# rank classes, not separately parameterised.

#' Sample Z under symmetric ordinal data (square symmetric matrix)
#'
#' Symmetric-Z analogue of \code{\link{rZ_ord_fc}}: each unordered actor pair
#' carries one latent value (\code{Z[i, j] = Z[j, i]}) drawn from a
#' truncated normal with mean equal to the average linear predictor
#' (EZ_ij + EZ_ji) / 2 and variance 1 (the model's conditional). The
#' lower triangle is mirrored after the upper-triangle draws so the full
#' matrix is exactly symmetric at every sweep.
#'
#' @param Z current symmetric latent matrix (n x n, Z = t(Z)).
#' @param EZ expected value of Z (regression + additive + multiplicative
#'   contributions); can be asymmetric, the helper averages the two
#'   directed contributions.
#' @param Y observed symmetric ordinal matrix (n x n, Y = t(Y));
#'   diagonal is ignored.
#' @return updated symmetric n x n latent matrix; the diagonal is refreshed
#'   with unconstrained normal draws around \code{diag(EZ)} so downstream
#'   samplers never see missing values.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_ord_sym_fc <- function(Z, EZ, Y) {
	n <- nrow(Y)
	uY <- sort(unique(c(Y)))
	W <- matrix(match(Y, uY), n, n)
	W[is.na(W)] <- -1L

	ut <- upper.tri(Z)
	sd_sym <- 1

	for (w in sample(c(-1L, seq_along(uY)))) {
		lb <- suppressWarnings(max(Z[!is.na(W) & W == w - 1L], na.rm = TRUE))
		ub <- suppressWarnings(min(Z[!is.na(W) & W == w + 1L], na.rm = TRUE))
		up <- ut & W == w
		if (!any(up)) next
		# symmetric conditional mean: average of (i,j) and (j,i) linear
		# predictors. variance 1 -> sd 1 (the model's conditional).
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
	# in rz_ord_fc / rz_bin_fc). setting `na` here propagates into the
	# symmetrised residual `e <- 0.5*(e + t(e))` consumed by ruv_sym_fc_cpp,
	# which then errors at r >= 1 with "mat::min(): object has no elements"
	# because every column of the residual has at least one na.
	diag(Z) <- stats::rnorm(nrow(Z), diag(EZ), 1)
	Z
}
