#' Gibbs update of regression coefficients and additive effects for
#' replicated relational data
#'
#' Draws jointly from the full conditional distribution of the regression
#' coefficients \code{beta} together with the additive sender/receiver
#' effects \code{a} and \code{b}, pooling the information contained in the
#' replicate slices of a longitudinal design that share a common
#' \code{beta}, \code{a} and \code{b}. The actor effects are eliminated
#' block by block over an orthonormal Helmert basis, the coefficients are
#' drawn from their marginal Gaussian conditional, and the actor effects are
#' then drawn given the coefficients. The coefficient prior is the g-prior
#' with \code{g = n^2}.
#'
#' @usage rbeta_ab_rep_fc(Z.T,Sab,rho,X.T,s2=1,center=FALSE)
#' @param Z.T n x n x T array of latent relations with multiplicative
#' effects already removed; the third margin indexes the replicates
#' @param Sab Positive semidefinite actor covariance, including zero or singular
#'   directions for disabled or tied actor effects. Every positive scale is retained.
#' @param rho within-dyad (reciprocity) correlation
#' @param X.T n x n x p x T design array
#' @param s2 dyadic variance
#' @param center Condition both actor margins to sum to zero.
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author lame authors
#' @keywords internal
#' @export rbeta_ab_rep_fc
rbeta_ab_rep_fc <- function(Z.T, Sab, rho, X.T, s2 = 1, center = FALSE) {
	n <- dim(Z.T)[1L]; nt <- dim(Z.T)[3L]; p <- dim(X.T)[3L]
	if (is.null(nt) || !identical(dim(Z.T), c(n, n, nt)) ||
			!identical(dim(X.T), c(n, n, p, nt)))
		cli::cli_abort("Z.T must be n x n x T and X.T must be n x n x p x T.")
	if (!is.matrix(Sab) || !identical(dim(Sab), c(2L, 2L)) || any(!is.finite(Sab)))
		cli::cli_abort("Sab must be a finite 2 x 2 covariance matrix.")
	mX <- array(X.T, c(n * n, p, nt))
	.lame_replicated_gaussian_draw(
		Z.T,
		mX,
		Sab,
		as.numeric(rho),
		s2,
		n * n,
		isTRUE(center)
	)
}
