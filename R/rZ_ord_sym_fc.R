# symmetric ordinal Z sampler: one latent value per unordered dyad, drawn on
# the upper triangle from a unit-variance truncated normal around the averaged
# linear predictor and mirrored to the lower triangle; cutpoints are implied by
# the neighbouring rank classes

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
#' @param rho dyadic correlation carried by the fit; it only sets the
#'   variance of the diagonal draw, \code{1 + rho}, the convention of the
#'   coefficient kernel.
#' @return updated symmetric n x n latent matrix; the diagonal is refreshed
#'   with unconstrained normal draws around \code{diag(EZ)} so downstream
#'   samplers never see missing values.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_ord_sym_fc <- function(Z, EZ, Y, rho = 0) {
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
		# conditional mean averages the (i,j) and (j,i) linear predictors
		ez_avg <- (EZ[up] + t(EZ)[up]) / 2
		# log-tail truncated draw that errors on an empty interval
		z_new <- .lame_checked_truncated_gaussian(ez_avg, sd_sym, lb, ub)
		Z[up] <- z_new
		# mirror to the (j,i) cells
		mirror_idx <- which(up, arr.ind = TRUE)[, c(2L, 1L), drop = FALSE]
		Z[mirror_idx] <- z_new
	}
	# unconstrained diagonal draw; an NA diagonal would break the symmetric UV update
	diag(Z) <- stats::rnorm(nrow(Z), diag(EZ), sqrt(1 + rho))
	Z
}
