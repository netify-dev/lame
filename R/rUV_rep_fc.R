#' Gibbs sampling of U and V from replicated relational data
#'
#' Draws the multiplicative-effect matrices U and V from their full
#' conditional distributions, pooling information across the replicate
#' slices of a residual array. Each factor column is updated in turn from
#' its Gaussian full conditional, optionally under a hierarchical
#' shrinkage prior on the stacked (U, V) rows.
#'
#' @usage rUV_rep_fc(E.T,U,V,rho,s2=1,shrink=TRUE)
#' @param E.T Array of square residual relational matrices (additive
#'   effects and covariates removed). The third margin indexes replicates.
#' @param U current value of U
#' @param V current value of V
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param shrink adaptively shrink the factors with a hierarchical prior
#'
#' @return \item{U}{a new value of U} \item{V}{a new value of V}
#' @keywords internal
#' @author lame authors
#' @export rUV_rep_fc
rUV_rep_fc <-
	function(E.T,U,V,rho,s2=1,shrink=TRUE) {

		grid_dim <- dim(E.T)
		n <- grid_dim[1]
		n_rep <- grid_dim[3]
		R <- ncol(U)

		## covariance of the stacked (U,V) rows: either resampled from its
		## inverse-Wishart full conditional (shrinkage) or fixed diffuse
		if(shrink) {
			cov_uv <- rSuv_fc(U, V, Suv0=diag(2*R)/(R+2), kappa0=R+2)
		} else {
			cov_uv <- diag(n, nrow=2*R)
		}

		## symmetric inverse square root of the 2x2 dyadic covariance
		## Se = s2 * matrix(c(1,rho,rho,1),2); the eigenbasis is (1,1)/sqrt2
		## and (1,-1)/sqrt2, giving a closed form for its whitening scalars
		inv_root_sum <- 1/sqrt(s2*(1+rho))
		inv_root_dif <- 1/sqrt(s2*(1-rho))
		w_self <- (inv_root_sum + inv_root_dif)/2
		w_swap <- (inv_root_sum - inv_root_dif)/2
		coef_iso <- w_self^2 + w_swap^2
		coef_mix <- 2*w_self*w_swap

		## residual summed over replicate slices
		E_sum <- rowSums(E.T, dims=2)

		## sample an n-vector from N(prec^{-1} rhs, prec^{-1}) using the
		## Cholesky factor of the precision matrix
		gaussian_from_precision <- function(prec, rhs) {
			upper <- chol(prec)
			post_mean <- backsolve(upper, backsolve(upper, rhs, transpose=TRUE))
			post_mean + backsolve(upper, rnorm(n))
		}

		## conditional prior mean-coefficients and variance for row-block
		## `pos` given the remaining 2R-1 factor coordinates
		prior_pieces <- function(pos) {
			others <- setdiff(seq_len(2*R), pos)
			reg <- as.numeric(cov_uv[pos, others] %*% solve(cov_uv[others, others]))
			cond_var <- as.numeric(cov_uv[pos, pos] - crossprod(reg, cov_uv[others, pos]))
			list(reg=reg, cond_var=cond_var)
		}

		for(r in sample.int(R)) {

			## residual with the r-th rank-one term stripped out, then
			## whitened and symmetrised (pooled across replicates)
			resid_sum <- E_sum - n_rep*tcrossprod(U[,-r,drop=FALSE], V[,-r,drop=FALSE])
			white_sum <- coef_iso*resid_sum + coef_mix*t(resid_sum)

			## ---- U[,r] | everything else ----
			vr <- V[,r]
			sq_v <- max(sum(vr^2), 1e-6)
			pu <- prior_pieces(r)
			mu_prior <- cbind(U[,-r,drop=FALSE], V) %*% pu$reg
			prec_u <- diag(n_rep*coef_iso*sq_v + 1/pu$cond_var, n) +
				(n_rep*coef_mix)*tcrossprod(vr)
			rhs_u <- white_sum %*% vr + mu_prior/pu$cond_var
			U[,r] <- gaussian_from_precision(prec_u, rhs_u)

			## ---- V[,r] | everything else ----
			ur <- U[,r]
			sq_u <- max(sum(ur^2), 1e-6)
			pv <- prior_pieces(R + r)
			mv_prior <- cbind(U, V[,-r,drop=FALSE]) %*% pv$reg
			prec_v <- diag(n_rep*coef_iso*sq_u + 1/pv$cond_var, n) +
				(n_rep*coef_mix)*tcrossprod(ur)
			rhs_v <- crossprod(white_sum, ur) + mv_prior/pv$cond_var
			V[,r] <- gaussian_from_precision(prec_v, rhs_v)
		}

		list(U=U, V=V)
	}
