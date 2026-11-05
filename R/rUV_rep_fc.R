#' Gibbs sampling of U and V from replicated relational data
#'
#' Draws the multiplicative-effect matrices U and V from their full
#' conditional distributions, pooling information across the replicate
#' slices of a residual array. Each factor column is updated in turn from
#' its Gaussian full conditional, optionally under a hierarchical
#' shrinkage prior on the stacked (U, V) rows.
#'
#' @usage rUV_rep_fc(E.T,U,V,rho,s2=1,shrink=TRUE,prior_df=NULL,prior_psi=NULL)
#' @param E.T Array of square residual relational matrices (additive
#'   effects and covariates removed). The third margin indexes replicates.
#' @param U current value of U
#' @param V current value of V
#' @param rho Dyadic correlation: a scalar or one value per replicate.
#' @param s2 dyadic variance
#' @param shrink adaptively shrink the factors with a hierarchical prior
#' @param prior_df degrees of freedom of the inverse-Wishart prior of the
#'   covariance of the stacked (U, V) rows (default \code{2R + 2}; the
#'   samplers pass the data-scaled default of \code{ame()}, see
#'   \code{prior$kappa0} there)
#' @param prior_psi the inverse-Wishart scale matrix \code{prior_psi * I}
#'   (default \code{prior_df}, the unit scale \code{Suv0 = I})
#'
#' @details
#' Given its partner column, a factor column has a Gaussian conditional
#' whose precision has one eigenvalue along the partner and a repeated
#' eigenvalue perpendicular to it; the draw uses this closed form, with the
#' symmetric and antisymmetric parts of the residual weighted by the
#' reciprocal variances \code{s2 (1 + rho)} and \code{s2 (1 - rho)}.
#'
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @keywords internal
#' @author lame authors
#' @export rUV_rep_fc
rUV_rep_fc <- function(E.T, U, V, rho, s2 = 1, shrink = TRUE, prior_df = NULL, prior_psi = NULL) {
	.lame_factor_validate_coordinates(E.T, U, V, rho, s2)
	grid_dim <- dim(E.T); n <- grid_dim[1L]; n_rep <- grid_dim[3L]; R <- ncol(U)
	if (R == 0L) return(list(U = U, V = V))
	# inverse-Wishart prior on the factor covariance (R/uv_prior.R), scale df * I by default
	if (is.null(prior_df)) prior_df <- .lame_uv_prior_df(R)
	if (is.null(prior_psi)) prior_psi <- prior_df
	# covariance of the stacked (U, V) rows: an inverse-Wishart draw under shrinkage, else fixed and diffuse
	precision_uv <- NULL
	if (shrink) {
		# iW(df, psi I) on a 2R-dimensional covariance must have df > 2R-1
		if (R >= 3L || sqrt(sum(U^2) + sum(V^2)) > 1e4) {
			# draw the Wishart precision from a QR factor so psi I + W'W is never inverted,
			# which stays stable when a factor direction is very large
			q <- qr(rbind(sqrt(prior_psi) * diag(2L * R), cbind(U, V)), LAPACK = TRUE)
			root <- backsolve(
				qr.R(q),
				matrix(
					rnorm((n + prior_df) * 2L * R),
					nrow = 2L * R
				)
			)
			precision_uv <- matrix(0, 2L * R, 2L * R)
			precision_uv[q$pivot, q$pivot] <- tcrossprod(root)
		} else {
			cov_uv <- rSuv_fc(U, V, Suv0 = diag(2 * R) * prior_psi / prior_df, kappa0 = prior_df)
		}
	} else {
		# diffuse covariance for the first half of burn-in: n times the prior scale
		cov_uv <- diag(n * prior_psi / prior_df, nrow = 2 * R)
	}


	if (is.null(precision_uv)) precision_uv <- .lame_factor_prior_precision(cov_uv)
	if (any(!is.finite(precision_uv)) || any(diag(precision_uv) <= 0))
		cli::cli_abort("Factor conditional prior precision must be finite and positive.")
	prior_pieces <- function(pos) {
		others <- setdiff(seq_len(2L * R), pos)
		list(
			reg = -precision_uv[pos, others] / precision_uv[pos, pos],
			precision = precision_uv[pos, pos]
		)
	}
	for (r in sample.int(R)) {
		others <- tcrossprod(U[, -r, drop = FALSE], V[, -r, drop = FALSE])
		partial <- E.T
		for (tt in seq_len(n_rep)) partial[, , tt] <- E.T[, , tt] - others
		signal <- .lame_factor_common_contrast(partial, rho, s2)
		pu <- prior_pieces(r)
		mu_prior <- drop(cbind(U[, -r, drop = FALSE], V) %*% pu$reg)
		U[, r] <- .lame_factor_gaussian_draw(V[, r], signal, mu_prior, pu$precision)
		pv <- prior_pieces(R + r)
		mv_prior <- drop(cbind(U, V[, -r, drop = FALSE]) %*% pv$reg)
		V[, r] <- .lame_factor_gaussian_draw(
			U[, r],
			signal,
			mv_prior,
			pv$precision,
			transpose = TRUE
		)
	}
	list(U = U, V = V)
}
