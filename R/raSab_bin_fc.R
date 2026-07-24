#' Sample additive row effects and their covariance for the binary family
#'
#' Joint Gibbs update of the social-relations additive row effects \code{a}
#' and the additive-effect covariance \code{Sab} under a probit (binary)
#' likelihood.  Given the latent normal scores \code{Z} and the observed 0/1
#' matrix \code{Y}, each row effect is redrawn from its Gaussian full
#' conditional (conditioned on the column effects \code{b} through \code{Sab})
#' truncated to the interval that keeps every latent score sign-compatible with
#' the observed responses.  The covariance \code{Sab} is then refreshed from its
#' inverse-Wishart full conditional.
#'
#' @param Z square matrix of current latent normal scores
#' @param Y square binary relational matrix (NA entries allowed)
#' @param a current value of the row effects
#' @param b current value of the column effects
#' @param Sab current 2x2 covariance of the additive effects
#' @param Sab0 prior inverse-scale matrix (defaults to the 2x2 identity)
#' @param eta0 prior degrees of freedom (defaults to 4)
#' @param SS number of inner Gibbs sweeps
#' @return list with the updated \code{Z}, \code{a} and \code{Sab}
#' @keywords internal
#' @author lame authors
#' @export raSab_bin_fc
raSab_bin_fc <-
function(Z,Y,a,b,Sab,Sab0=NULL,eta0=NULL,SS=round(sqrt(nrow(Z)))) {
	if(is.null(Sab0)){ Sab0 <- diag(2) }
	if(is.null(eta0)){ eta0 <- 4 }

	n <- nrow(Z)

	# single inverse-Wishart draw with scale matrix psi and df degrees of
	# freedom, built from a Bartlett-decomposed Wishart on the inverse scale
	draw_inv_wishart <- function(psi, df) {
		p <- nrow(psi)
		root <- t(chol(solve(psi)))          # lower factor: solve(psi) = root %*% t(root)
		bart <- matrix(0, p, p)
		diag(bart) <- sqrt(rchisq(p, df - seq_len(p) + 1))
		if(p > 1L){ bart[lower.tri(bart)] <- rnorm(p * (p - 1L) / 2L) }
		wfac <- root %*% bart                # wfac %*% t(wfac) ~ Wishart(solve(psi), df)
		solve(tcrossprod(wfac))              # its inverse is the IW(psi, df) draw
	}

	# strip the current row effect from each latent score; the remainder is
	# held fixed while a is resampled, so Z_ij = carry_ij + a_i throughout
	carry <- Z - matrix(a, n, n)

	# probit sign constraints on Z_ij = carry_ij + a_i give, per row i,
	#   Y_ij == 1  ->  a_i > -carry_ij   (a lower limit)
	#   Y_ij == 0  ->  a_i < -carry_ij   (an upper limit)
	# missing cells contribute no limit
	cut <- -carry
	low_src <- cut ; low_src[is.na(Y) | Y != 1] <- -Inf
	high_src <- cut ; high_src[is.na(Y) | Y != 0] <- Inf
	a_low <- apply(low_src, 1, max)
	a_high <- apply(high_src, 1, min)
	a_low[is.na(a_low)] <- -Inf
	a_high[is.na(a_high)] <- Inf

	for(ss in seq_len(SS)) {
		# a | b is Gaussian: regress a on b through the 2x2 covariance Sab
		slope <- Sab[1, 2] / Sab[2, 2]
		cond_mean <- slope * b
		cond_sd <- sqrt(Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2])

		# inverse-CDF sampling of the truncated normal on (a_low, a_high)
		p_low <- pnorm((a_low - cond_mean) / cond_sd)
		p_high <- pnorm((a_high - cond_mean) / cond_sd)
		a <- cond_mean + cond_sd * qnorm(runif(n, p_low, p_high))

		# Sab ~ inverse-Wishart(eta0 * Sab0 + [a b]'[a b], eta0 + n)
		ab <- cbind(a, b)
		Sab <- draw_inv_wishart(eta0 * Sab0 + crossprod(ab), eta0 + n)
	}

	list(Z = carry + matrix(a, n, n), a = a, Sab = Sab)
}
