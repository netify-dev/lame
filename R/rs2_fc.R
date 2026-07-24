#' Gibbs update for dyadic variance
#'
#' Draw the dyadic residual variance \code{s2} from its inverse-gamma full
#' conditional. The residual matrix \code{E = Z - offset} is decorrelated
#' within reciprocal dyad pairs and along the diagonal, and the pooled sum of
#' squares (together with an optional scaled-inverse-chi-square prior) yields
#' the posterior draw.
#'
#' @param Z n x n normal relational matrix
#' @param rho current within-dyad correlation
#' @param offset matrix conformable with \code{Z}; \code{Z - offset} is treated
#'   as pure dyadic noise (so it should carry the additive and multiplicative
#'   mean structure)
#' @param nu0 prior degrees of freedom (defaults to 1)
#' @param s20 prior scale for s2; when \code{NULL} the empirical mean of the
#'   decorrelated squared residuals is used so the prior sits at the data scale
#' @return a single draw of s2
#' @keywords internal
#' @author lame authors
#' @export rs2_fc
rs2_fc <-
function(Z, rho, offset=0, nu0=NULL, s20=NULL) {

	E <- Z - offset
	if (is.null(nu0)) { nu0 <- 1 }

	# The exchangeable 2x2 dyad covariance [[1, rho], [rho, 1]] is diagonalised
	# by the orthonormal sum / difference contrasts: the sum eigenvector
	# (1, 1)/sqrt(2) has eigenvalue 1 + rho and the difference eigenvector
	# (1, -1)/sqrt(2) has eigenvalue 1 - rho. Rescaling each contrast by the
	# square root of its eigenvalue whitens the pair, so the accumulated squares
	# equal the Mahalanobis quadratic form E[ij,ji] %*% R^{-1} %*% E[ij,ji].
	above <- upper.tri(E)
	e_up  <- E[above]
	e_lo  <- t(E)[above]

	w_sum  <- (e_up + e_lo) / sqrt(2 * (1 + rho))
	w_diff <- (e_up - e_lo) / sqrt(2 * (1 - rho))

	# Self-loops (the diagonal) carry variance proportional to 1 + rho.
	w_diag <- diag(E) / sqrt(1 + rho)

	whitened <- c(w_sum, w_diff, w_diag)
	n_terms  <- length(whitened)
	ss       <- sum(whitened^2)

	if (is.null(s20)) {
		ok  <- is.finite(whitened)
		bar <- sum(whitened[ok]^2) / max(sum(ok), 1)
		s20 <- if (is.finite(bar) && bar > 0) bar else 1
	}

	shape <- (n_terms + nu0) / 2
	rate  <- (ss + nu0 * s20) / 2

	1 / rgamma(1, shape, rate)
}
