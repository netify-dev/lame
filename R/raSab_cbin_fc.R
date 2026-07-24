#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#'
#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#'
#'
#' @usage raSab_cbin_fc(Z, Y, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,SS =
#' round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each individual
#' @param odobs observed outdegree
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @keywords internal
#' @author lame authors
#' @export raSab_cbin_fc
raSab_cbin_fc <-
function(Z, Y, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL, SS=round(sqrt(nrow(Z)))) {

	# --- prior defaults for the inverse-Wishart on Sab ---
	if(is.null(Sab0)) { Sab0 <- diag(2) }
	if(is.null(eta0)) { eta0 <- 4 }

	n <- nrow(Z)

	# Strip the current additive row effect so that, for a proposed row
	# effect a_i, the reconstructed latent value is resid[i,j] + a_i.
	resid <- Z - matrix(a, nrow=n, ncol=n)

	# For a proposed a_i, entry (i,j) has latent value a_i + resid[i,j];
	# a positive value encodes Y_ij = 1 and a negative value Y_ij = 0.
	# Hence each observed edge state translates into a one-sided bound on
	# a_i, namely  a_i > -resid[i,j]  (Y=1)  or  a_i < -resid[i,j]  (Y=0).
	thr <- -resid

	# Assemble the tightest interval [lo_i, hi_i] each row effect must obey.
	# Only realized ones push the floor up; only realized zeros pull the
	# ceiling down. Missing cells and cells of the wrong sign contribute a
	# vacuous (+/-Inf) bound.
	one_cell  <- !is.na(Y) & (Y == 1)
	zero_cell <- !is.na(Y) & (Y == 0)

	lo_mat <- matrix(-Inf, nrow=n, ncol=n)
	hi_mat <- matrix( Inf, nrow=n, ncol=n)
	lo_mat[one_cell]  <- thr[one_cell]
	hi_mat[zero_cell] <- thr[zero_cell]

	lo <- apply(lo_mat, 1, max)
	hi <- apply(hi_mat, 1, min)

	# Censored-binary correction: a row that has already spent its full
	# nomination budget (odobs == odmax) may have suppressed genuine ones,
	# so its zeros carry no information and the ceiling is released.
	hi[odobs == odmax] <- Inf

	# Gibbs sweeps alternating a | (b, Sab) and Sab | (a, b).
	for(iter in seq_len(SS)) {

		# Conditional of a given b under the bivariate normal (a,b) ~ N(0,Sab):
		#   E[a|b] = b * Sab_ab / Sab_bb ,  Var[a|b] = Sab_aa - Sab_ab^2 / Sab_bb.
		cond_mean <- b * (Sab[1, 2] / Sab[2, 2])
		cond_sd   <- sqrt(Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2])

		# Draw each a_i from that normal truncated to [lo_i, hi_i] via the
		# inverse-CDF (probability-integral) method.
		z_lo <- pnorm((lo - cond_mean) / cond_sd)
		z_hi <- pnorm((hi - cond_mean) / cond_sd)
		u    <- runif(n, z_lo, z_hi)
		a    <- cond_mean + cond_sd * qnorm(u)

		# Inverse-Wishart update for Sab. With posterior scale
		# S_post = eta0*Sab0 + [a b]'[a b] and df = eta0 + n, sample a
		# Wishart precision from the inverted scale and invert it back.
		ab_mat <- cbind(a, b)
		S_post <- eta0 * Sab0 + t(ab_mat) %*% ab_mat
		Sab    <- solve(rwish(solve(S_post), eta0 + n))
	}

	# Rebuild Z with the refreshed row effect folded back in.
	list(Z = resid + matrix(a, nrow=n, ncol=n), a = a, Sab = Sab)
}