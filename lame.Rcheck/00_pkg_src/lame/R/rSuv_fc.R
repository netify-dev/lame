#' Gibbs update for multiplicative effects covariance
#'
#' Draws the covariance matrix of the stacked multiplicative row/column
#' effects \eqn{[U, V]} from its full conditional inverse-Wishart
#' distribution in the AME model.
#'
#' @usage rSuv_fc(U, V, Suv0=NULL, kappa0=NULL)
#' @param U matrix of multiplicative row effects (n x R).
#' @param V matrix of multiplicative column effects (n x R).
#' @param Suv0 prior scale matrix (2R x 2R). Defaults to the identity,
#'   a weakly informative choice.
#' @param kappa0 prior degrees of freedom. Defaults to 2 + 2R, the
#'   smallest value giving a proper prior for a 2R x 2R covariance.
#' @return The sampled 2R x 2R covariance matrix for \eqn{[U, V]}: the
#'   leading R x R block is the covariance of U, the trailing R x R block
#'   is the covariance of V, and the off-diagonal blocks are the U-V
#'   cross-covariances.
#' @details
#' Stacking the effects columnwise as \eqn{W = [U, V]}, the conjugate
#' inverse-Wishart update combines the prior scale \code{kappa0 * Suv0}
#' with the residual cross-product \code{crossprod(W)} and adds the n
#' observed rows to the degrees of freedom. A draw from the inverse
#' Wishart is obtained by drawing from the Wishart with the inverted
#' scale matrix (via \code{rwish}) and inverting the result.
#' @author lame authors
#' @keywords internal
#' @export rSuv_fc
rSuv_fc <- function(U, V, Suv0=NULL, kappa0=NULL)  {
	n <- nrow(U)
	R <- ncol(U)
	dim_uv <- 2 * R

	if (is.null(Suv0)) { Suv0 <- diag(dim_uv) }
	if (is.null(kappa0)) { kappa0 <- 2 + dim_uv }

	stacked <- cbind(U, V)
	cross <- crossprod(stacked)

	post_scale <- kappa0 * Suv0 + cross
	post_df <- kappa0 + n

	wishart_draw <- rwish(solve(post_scale), post_df)
	Suv <- solve(wishart_draw)

	Suv
}
