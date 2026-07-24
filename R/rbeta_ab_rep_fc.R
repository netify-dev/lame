#' Gibbs update of regression coefficients and additive effects for
#' replicated relational data
#'
#' Draws jointly from the full conditional distribution of the regression
#' coefficients \code{beta} together with the additive sender/receiver
#' effects \code{a} and \code{b}, pooling the information contained in the
#' replicate slices of a longitudinal design that share a common
#' \code{beta}, \code{a} and \code{b}. The draw is obtained by forming the
#' complete joint Gaussian posterior precision matrix over
#' \code{(beta, a, b)} from sufficient statistics accumulated across slices,
#' after whitening the within-dyad correlation, and sampling from it.
#'
#' @usage rbeta_ab_rep_fc(Z.T,Sab,rho,X.T,s2=1)
#' @param Z.T n x n x T array of latent relations with multiplicative
#' effects already removed; the third margin indexes the replicates
#' @param Sab 2 x 2 covariance matrix of the additive row/column effects
#' @param rho within-dyad (reciprocity) correlation
#' @param X.T n x n x p x T design array
#' @param s2 dyadic variance
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author lame authors
#' @keywords internal
#' @export rbeta_ab_rep_fc
rbeta_ab_rep_fc <-
	function(Z.T,Sab,rho,X.T,s2=1)  {
		n_rep <- dim(X.T)[4]
		p <- dim(X.T)[3]
		n <- dim(Z.T)[1]

		# --- whitening of the within-dyad covariance -------------------
		# Each ordered pair (i,j) carries error covariance s2*[[1,rho],
		# [rho,1]] with its mirror (j,i). Left-multiplying the stacked
		# pair by W = (s2*[[1,rho],[rho,1]])^{-1/2} makes every directed
		# observation unit-variance and mutually uncorrelated. W is
		# symmetric with a common diagonal (w_d) and off-diagonal (w_o).
		dyad_cov <- s2 * matrix(c(1, rho, rho, 1), 2, 2)
		W <- mhalf(solve(dyad_cov))
		w_d <- W[1, 1]
		w_o <- W[1, 2]

		# --- prior precision of a single actor's (a_i,b_i) pair ---------
		Sab <- (Sab + t(Sab)) / 2
		if (any(!is.finite(Sab))) { Sab <- diag(2) }
		lam_min <- min(eigen(Sab, symmetric = TRUE, only.values = TRUE)$values)
		if (!is.finite(lam_min) || lam_min < 1e-10) {
			Sab <- Sab + diag(1e-6, 2)
		}
		Pab <- solve(Sab)

		# --- sufficient statistics accumulated over the replicates -----
		# In the whitened frame directed observation (i,j) has mean
		#   x_ij' beta + w_d*a_i + w_o*a_j + w_d*b_j + w_o*b_i .
		Pbb  <- matrix(0, p, p)   # beta cross-product
		g_b  <- numeric(p)        # beta linear term
		Dr   <- matrix(0, n, p)   # whitened row-summed design, actor x cov
		Dc   <- matrix(0, n, p)   # whitened col-summed design
		yr   <- numeric(n)        # whitened row sums of the data
		yc   <- numeric(n)        # whitened col sums of the data

		for (r in seq_len(n_rep)) {
			Z  <- Z.T[, , r]
			Zw <- w_d * Z + w_o * t(Z)
			yr <- yr + rowSums(Zw)
			yc <- yc + colSums(Zw)

			if (p > 0) {
				Xr <- array(X.T[, , , r], dim = c(n, n, p))
				# flatten each covariate slice and its transpose
				Fx  <- apply(Xr, 3, as.vector)                       # (n*n) x p
				Fxt <- apply(aperm(Xr, c(2, 1, 3)), 3, as.vector)    # transposed
				Dw  <- w_d * Fx + w_o * Fxt                          # whitened design
				Pbb <- Pbb + crossprod(Dw)
				# faint g-prior style ridge on beta (matches reference scaling)
				Pbb <- Pbb + crossprod(Fx) / (n * n) / n_rep
				g_b <- g_b + crossprod(Dw, as.vector(Zw))

				rowSum <- apply(Xr, c(1, 3), sum)   # sum over receivers, n x p
				colSum <- apply(Xr, c(2, 3), sum)   # sum over senders,   n x p
				Dr <- Dr + (w_d * rowSum + w_o * colSum)
				Dc <- Dc + (w_d * colSum + w_o * rowSum)
			}
		}

		# --- precision blocks for the additive effects -----------------
		# Summing the whitened loading vectors over all ordered pairs and
		# replicates gives closed forms in terms of I and the all-ones J.
		q_ii <- (w_d^2 + w_o^2) * n * n_rep   # diagonal weight
		q_j  <- 2 * w_d * w_o * n_rep         # J weight for like blocks
		I_n  <- diag(n)
		J_n  <- matrix(1, n, n)

		Paa <- q_ii * I_n + q_j * J_n + Pab[1, 1] * I_n
		Pbb_ab <- q_ii * I_n + q_j * J_n + Pab[2, 2] * I_n
		Pab_x <- (w_d^2 + w_o^2) * n_rep * J_n +
			2 * w_d * w_o * n * n_rep * I_n + Pab[1, 2] * I_n

		# cross precision between beta and the additive effects
		Pb_a <- w_d * t(Dr) + w_o * t(Dc)   # p x n
		Pb_b <- w_d * t(Dc) + w_o * t(Dr)   # p x n

		# linear terms for a and b
		g_a <- w_d * yr + w_o * yc
		g_bb <- w_d * yc + w_o * yr

		# --- assemble the joint posterior and draw ---------------------
		if (p > 0) {
			Prec <- rbind(
				cbind(Pbb,       Pb_a,   Pb_b),
				cbind(t(Pb_a),   Paa,    Pab_x),
				cbind(t(Pb_b),   t(Pab_x), Pbb_ab)
			)
			g <- c(g_b, g_a, g_bb)
			off <- p
		} else {
			Prec <- rbind(
				cbind(Paa,    Pab_x),
				cbind(t(Pab_x), Pbb_ab)
			)
			g <- c(g_a, g_bb)
			off <- 0
		}

		post_cov  <- solve(Prec)
		post_mean <- post_cov %*% g
		draw <- c(rmvnorm(1, c(post_mean), post_cov))

		beta <- if (p > 0) draw[seq_len(p)] else numeric(0)
		a <- draw[off + seq_len(n)]
		b <- draw[off + n + seq_len(n)]

		list(beta = beta, a = a, b = b)
	}
