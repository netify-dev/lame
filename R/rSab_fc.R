#' Gibbs update for additive effects covariance
#'
#' Gibbs sampling for the covariance matrix of additive row and column
#' effects in the AME model. This function implements the inverse-Wishart
#' posterior update for the covariance matrix Sab.
#'
#' @usage rSab_fc(a, b, Sab0=NULL, eta0=NULL, rvar=TRUE, cvar=TRUE, symmetric=FALSE,
#'   center=FALSE, tiny=1e-8)
#' @param a vector of row random effects (additive sender effects)
#' @param b vector of column random effects (additive receiver effects)
#' @param Sab0 prior scale matrix. \strong{Convention (Hoff/amen):}
#'   \code{Sab0} is passed \emph{raw}; the inverse-Wishart scale used
#'   internally is \eqn{\eta_0 \cdot Sab0}, so the prior mean of Sab is
#'   \eqn{\eta_0 Sab0/(\eta_0-3) \approx Sab0} and the prior mode is
#'   \eqn{\eta_0 Sab0/(\eta_0+3)}. Every call site in the package (ame()
#'   unipartite and bipartite, lame() unipartite and bipartite, and the
#'   raSab_bin/cbin/frn joint updates) passes \code{Sab0} raw under this
#'   same convention, so a given \code{prior$Sab0} means the same prior
#'   everywhere. Default is diag(2) (weakly informative).
#' @param eta0 prior degrees of freedom for the prior distribution.
#'   Default is 4, the smallest integer giving a finite prior mean for a 2x2 covariance.
#' @param rvar logical: should row variance be updated? (default TRUE)
#' @param cvar logical: should column variance be updated? (default TRUE)
#' @param symmetric logical: is this a symmetric network? (default FALSE)
#' @param center Logical: are the supplied actor effects constrained to sum
#'   to zero? If so, remove one degree of freedom from the covariance update.
#'   Defaults to FALSE for unconstrained effects.
#' @param tiny the variance a margin that is switched off (\code{rvar} or \code{cvar} FALSE) keeps: 1e-8 in
#'   the unit of the response (default 1e-8, the value for a link-fixed latent scale).
#' @return Updated covariance matrix Sab (2x2 matrix with variances on diagonal
#'   and covariance off-diagonal)
#' @details
#' The function implements different update strategies:
#' \itemize{
#'   \item Full update: When both rvar and cvar are TRUE, updates the full 2x2
#'     covariance matrix using an inverse-Wishart distribution
#'   \item Row variance only: When only rvar is TRUE, updates only Sab\[1,1\]
#'   \item Column variance only: When only cvar is TRUE, updates only Sab\[2,2\]
#'   \item Symmetric case: When symmetric is TRUE, draws a single shared
#'     variance for the row and column effects (inverse-gamma update) with
#'     zero covariance
#' }
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
#' @export rSab_fc
rSab_fc <- function(
	a,
	b,
	Sab0 = NULL,
	eta0 = NULL,
	rvar = TRUE,
	cvar = TRUE,
	symmetric = FALSE,
	center = FALSE,
	tiny = 1e-8
) {
	if (is.null(Sab0)) { Sab0 <- diag(2) }
	if (is.null(eta0)) { eta0 <- 4 }

	# actors unobserved in every period have effects pinned at 0; leave them out of
	# the inverse-Wishart degrees of freedom so va/vb are not biased down
	informative <- (a != 0) | (b != 0)
	n <- sum(informative)
	if (n == 0L) n <- length(a) # degenerate start: use the full count
	a_inf <- a[informative]; b_inf <- b[informative]
	if (isTRUE(center)) n <- max(n - 1L, 0L)
	Sab <- matrix(0, 2, 2)

	# full covariance update
	if (rvar & cvar & !symmetric) {
		Sab <- solve(rwish(solve(eta0 * Sab0 + crossprod(cbind(a_inf, b_inf))), eta0 + n))
	}

	# row variance only
	if (rvar & !cvar & !symmetric) {
		sum_a2 <- sum(a_inf^2)
		if (!is.finite(sum_a2) || sum_a2 < 0) {
			stop("Nonfinite additive scatter cannot define an inverse-gamma conditional.")
		}
		rate_param <- (eta0 * Sab0[1, 1] + sum_a2) / 2
		if (!is.finite(rate_param) || rate_param <= 0) {
			stop("Inverse-gamma rate must be finite and positive.")
		}
		Sab[1, 1] <- 1 / rgamma(1, (eta0 + n) / 2, rate_param)
		Sab[2, 2] <- tiny
		Sab[1, 2] <- Sab[2, 1] <- 0
	}

	# column variance only
	if (!rvar & cvar & !symmetric) {
		Sab[1, 1] <- tiny
		sum_b2 <- sum(b_inf^2)
		if (!is.finite(sum_b2) || sum_b2 < 0) {
			stop("Nonfinite additive scatter cannot define an inverse-gamma conditional.")
		}
		rate_param <- (eta0 * Sab0[2, 2] + sum_b2) / 2
		if (!is.finite(rate_param) || rate_param <= 0) {
			stop("Inverse-gamma rate must be finite and positive.")
		}
		Sab[2, 2] <- 1 / rgamma(1, (eta0 + n) / 2, rate_param)
		Sab[1, 2] <- Sab[2, 1] <- 0
	}

	# symmetric case
	if (symmetric) {
		sum_a2 <- sum(a_inf^2)
		if (!is.finite(sum_a2) || sum_a2 < 0) {
			stop("Nonfinite additive scatter cannot define an inverse-gamma conditional.")
		}
		rate_param <- (eta0 * Sab0[1, 1] + sum_a2) / 2
		if (!is.finite(rate_param) || rate_param <= 0) {
			stop("Inverse-gamma rate must be finite and positive.")
		}
		var_ab <- 1 / rgamma(1, (eta0 + n) / 2, rate_param)
		Sab[1, 1] <- Sab[2, 2] <- var_ab
		Sab[1, 2] <- Sab[2, 1] <- 0
	}

	Sab <- (Sab + t(Sab)) / 2

	if (any(!is.finite(Sab)) ||
			((isTRUE(rvar) || isTRUE(symmetric)) && Sab[1, 1] <= 0) ||
			((isTRUE(cvar) || isTRUE(symmetric)) && Sab[2, 2] <= 0)) {
		stop("Unable to represent the additive covariance conditional draw.")
	}
	return(Sab)
}
