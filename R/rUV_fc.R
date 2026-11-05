#' Gibbs sampling of U and V
#'
#' A Gibbs sampler for updating the multiplicative effect matrices U and V
#'
#' @param Z n X n normal relational matrix
#' @param U current value of U
#' @param V current value of V
#' @param Suv covariance of (U V)
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param offset a matrix of the same dimension as Z. It is assumed that
#' Z-offset is equal to the multiplicative effects plus dyadic noise, so the
#' offset should contain any additive effects (such as \code{Xbeta(X,beta+
#' outer(a,b,"+")  }  )
#'
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @keywords internal
#' @author lame authors
#' @export rUV_fc
rUV_fc <- function(Z, U, V, Suv, rho, s2 = 1, offset = 0) {
	resid <- Z - offset
	E <- array(resid, c(nrow(resid), ncol(resid), 1L))
	.lame_factor_validate_coordinates(E, U, V, rho, s2)
	n_actor <- nrow(U); dim_uv <- ncol(U)
	if (dim_uv == 0L) return(list(U = U, V = V))
	if (!identical(dim(Suv), c(2L * dim_uv, 2L * dim_uv)))
		cli::cli_abort("Factor prior covariance dimensions must match U and V.")
	precision_uv <- .lame_factor_prior_precision(Suv)
	for (r in sample(seq_len(dim_uv))) {
		other <- setdiff(seq_len(dim_uv), r)
		partial_resid <- resid - U[, other, drop = FALSE] %*% t(V[, other, drop = FALSE])
		signal <- .lame_factor_common_contrast(array(partial_resid, c(n_actor, n_actor, 1L)), rho, s2)
		free <- setdiff(seq_len(2L * dim_uv), r)
		prior_precision <- precision_uv[r, r]
		reg_coef <- -precision_uv[r, free] / prior_precision
		cond_mean <- drop(cbind(U[, other, drop = FALSE], V) %*% reg_coef)
		U[, r] <- .lame_factor_gaussian_draw(V[, r], signal, cond_mean, prior_precision)
		pos <- dim_uv + r; free <- setdiff(seq_len(2L * dim_uv), pos)
		prior_precision <- precision_uv[pos, pos]
		reg_coef <- -precision_uv[pos, free] / prior_precision
		cond_mean <- drop(cbind(U, V[, other, drop = FALSE]) %*% reg_coef)
		V[, r] <- .lame_factor_gaussian_draw(
			U[, r],
			signal,
			cond_mean,
			prior_precision,
			transpose = TRUE
		)
	}
	list(U = U, V = V)
}
