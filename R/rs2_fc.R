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
#' @param s20 fixed positive prior scale for s2 (default 1). Model-fitting
#'   callers compute a data-scaled value once before sampling; do not
#'   recompute this scale from current residuals within a chain.
#' @return a single draw of s2
#' @keywords internal
#' @author lame authors
#' @export rs2_fc
rs2_fc <-
	function(Z, rho, offset = 0, nu0 = NULL, s20 = NULL) {

		E <- Z - offset
		if (is.null(nu0)) { nu0 <- 1 }
		if (!is.numeric(nu0) || length(nu0) != 1L || !is.finite(nu0) || nu0 < 0) {
			cli::cli_abort("{.arg nu0} must be a finite non-negative scalar.")
		}
		s20 <- .lame_validate_s20(if (is.null(s20)) 1 else s20)

		# the sum and difference contrasts diagonalise the dyad correlation (eigenvalues
		# 1 + rho and 1 - rho), so their rescaled squares give the Mahalanobis form
		above <- upper.tri(E)
		e_up <- E[above]
		e_lo <- t(E)[above]

		w_sum <- (e_up + e_lo) / sqrt(2 * (1 + rho))
		w_diff <- (e_up - e_lo) / sqrt(2 * (1 - rho))

		# the diagonal carries variance proportional to 1 + rho
		w_diag <- diag(E) / sqrt(1 + rho)

		whitened <- c(w_sum, w_diff, w_diag)
		n_terms <- length(whitened)
		ss <- sum(whitened^2)

		shape <- (n_terms + nu0) / 2
		rate <- (ss + nu0 * s20) / 2

		1 / rgamma(1, shape, rate)
	}
