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
rUV_fc <-
function(Z, U, V, Suv, rho, s2 = 1, offset = 0) {

	resid <- Z - offset
	n_actor <- nrow(resid)
	dim_uv <- ncol(U)

	# Symmetric inverse square-root of the 2x2 dyadic noise covariance
	# Se = s2 * [[1,rho],[rho,1]].  Its eigenpairs are (1,1)/sqrt(2) with
	# eigenvalue s2*(1+rho) and (1,-1)/sqrt(2) with eigenvalue s2*(1-rho),
	# so Se^{-1/2} = [[wg,wd],[wd,wg]] with the closed form below.
	inv_root_plus <- 1 / sqrt(s2 * (1 + rho))
	inv_root_minus <- 1 / sqrt(s2 * (1 - rho))
	wg <- 0.5 * (inv_root_plus + inv_root_minus)
	wd <- 0.5 * (inv_root_plus - inv_root_minus)
	diag_scale <- wg^2 + wd^2
	cross_scale <- 2 * wg * wd

	# Draw a single actor-column of a factor matrix from its Gaussian full
	# conditional N(mu, Lambda^{-1}), with
	#   Lambda = diag_scale * ||partner||^2 * I + cross_scale * partner partner'
	#            + (1/prior_var) * I
	#   info   = whitened_signal + prior_mean / prior_var
	# using an explicit precision matrix and a Cholesky-based sample.
	draw_factor_column <- function(partner, whitened_signal, prior_mean, prior_var) {
		partner_ss <- sum(partner^2)
		precision <- diag(diag_scale * partner_ss + 1 / prior_var, n_actor)
		precision <- precision + cross_scale * tcrossprod(partner)
		info <- whitened_signal + prior_mean / prior_var
		root <- chol(precision)
		post_mean <- backsolve(root, backsolve(root, info, transpose = TRUE))
		post_mean + backsolve(root, rnorm(n_actor))
	}

	for (r in sample(seq_len(dim_uv))) {

		# Residual attributable to the r-th latent dimension, and its
		# precision-whitened form Es = diag_scale*Er + cross_scale*Er'.
		other <- setdiff(seq_len(dim_uv), r)
		partial_resid <- resid - U[, other, drop = FALSE] %*% t(V[, other, drop = FALSE])
		whitened_resid <- diag_scale * partial_resid + cross_scale * t(partial_resid)

		## --- update column r of U, conditioning on the current V ---
		u_index <- r
		free <- setdiff(seq_len(2 * dim_uv), u_index)
		reg_coef <- as.numeric(Suv[u_index, free] %*% solve(Suv[free, free]))
		cond_var <- as.numeric(Suv[u_index, u_index] - reg_coef %*% Suv[free, u_index])
		cond_mean <- cbind(U[, other, drop = FALSE], V) %*% reg_coef
		v_col <- V[, r]
		U[, r] <- draw_factor_column(
			partner = v_col,
			whitened_signal = whitened_resid %*% v_col,
			prior_mean = cond_mean,
			prior_var = cond_var)

		## --- update column r of V, conditioning on the refreshed U ---
		v_index <- dim_uv + r
		free <- setdiff(seq_len(2 * dim_uv), v_index)
		reg_coef <- as.numeric(Suv[v_index, free] %*% solve(Suv[free, free]))
		cond_var <- as.numeric(Suv[v_index, v_index] - reg_coef %*% Suv[free, v_index])
		cond_mean <- cbind(U, V[, other, drop = FALSE]) %*% reg_coef
		u_col <- U[, r]
		V[, r] <- draw_factor_column(
			partner = u_col,
			whitened_signal = crossprod(whitened_resid, u_col),
			prior_mean = cond_mean,
			prior_var = cond_var)
	}

	list(U = U, V = V)
}