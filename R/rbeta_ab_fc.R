#' Joint Gibbs update of regression and additive effects (single relation)
#'
#' Draws jointly from the full conditional of the regression coefficients
#' \code{beta} together with the additive sender/receiver effects
#' \code{a}, \code{b} for a single social-relations regression matrix. The
#' within-dyad reciprocity covariance is whitened so every directed cell
#' becomes a unit-variance Gaussian observation, the additive-effect prior
#' \code{N(0,Sab)} is written through a square-root factor so a
#' rank-deficient \code{Sab} is handled without inverting it, and the
#' complete joint Gaussian precision over \code{(beta, factors)} is
#' assembled from cell-level sufficient statistics and sampled in one draw.
#'
#' @param Z n x n relational matrix (multiplicative effects already removed)
#' @param Sab 2 x 2 covariance of the additive row/column effects
#' @param rho within-dyad (reciprocity) correlation
#' @param X n x n x p covariate design array
#' @param s2 dyadic variance
#' @param offset matrix subtracted from \code{Z} before sampling
#' @param iV0 prior precision for \code{beta}; a g-prior is used when NULL
#' @param m0 prior mean for \code{beta}; zero when NULL
#' @param g g-prior variance scale used when \code{iV0} is NULL
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @keywords internal
#' @author lame authors
#' @export rbeta_ab_fc
rbeta_ab_fc <-
	function(Z,Sab,rho,X=NULL,s2=1,offset=0,iV0=NULL,m0=NULL,g=length(Z)) {

		# response is measured relative to the supplied offset
		Z <- Z - offset
		n <- nrow(Z)

		# a bare additive-effects design is used when no covariates are given
		if (is.null(X)) { X <- design_array(intercept = FALSE, n = n) }

		# cell-level sufficient statistics must be attached to the design
		if (is.null(attributes(X)$XX)) {
			X <- precomputeX(X)
			cli::cli_warn(paste0("Summary statistics of X are not precomputed. ",
					"Run X<-precomputeX(X) to speed up calculations."))
		}

		# statistics summed over all n^2 directed cells (produced by precomputeX)
		mX  <- attributes(X)$mX     # (n^2) x p   stacked design
		mXt <- attributes(X)$mXt    # (n^2) x p   stacked transposed design
		XX  <- attributes(X)$XX     # sum_ij x_ij x_ij'
		XXt <- attributes(X)$XXt    # sum_ij x_ij x_ji'
		Xr  <- attributes(X)$Xr     # n x p       out-sum of the design per actor
		Xc  <- attributes(X)$Xc     # n x p       in-sum of the design per actor
		p   <- dim(X)[3]

		# ---- prior precision for beta ---------------------------------
		# When the caller gives none, use a Zellner g-prior. A vanishing
		# ridge (a fixed fraction of the g-prior, so it shares the same
		# 1/scale(Y)^2 units) conditions the precision without ever
		# dominating it; an absolute ridge would gain a spurious scale
		# factor and crush the slopes on small-scale responses.
		if (p > 0 && is.null(iV0)) {
			iV0 <- .beta_prior_precision(XX, mX, g)
		}
		if (is.null(m0)) { m0 <- rep(0, p) }

		# ---- whiten the within-dyad reciprocity covariance ------------
		# Directed cell (i,j) is scaled so that, summed over every ordered
		# cell, the objective is the GLS quadratic form for the dyadic
		# covariance s2*[[1,rho],[rho,1]]. W is symmetric with diagonal wd
		# and off-diagonal wo.
		dyad_cov <- s2 * matrix(c(1, rho, rho, 1), 2, 2)
		W  <- mhalf(solve(dyad_cov))
		wd <- W[1, 1]
		wo <- W[1, 2]

		# whitened response and its per-actor out/in margins
		Zw <- wd * Z + wo * t(Z)
		Zw[is.na(Zw)] <- 0
		z_out <- rowSums(Zw)
		z_in  <- colSums(Zw)

		# whitened design sufficient statistics
		XXw <- (wd^2 + wo^2) * XX + 2 * wd * wo * XXt   # beta cross-product
		mXw <- wd * mX + wo * mXt                       # stacked whitened design
		Dout <- wd * Xr + wo * Xc                       # whitened out-sum design
		Din  <- wd * Xc + wo * Xr                       # whitened in-sum design

		# ---- square-root factor of the additive-effect prior ----------
		# Write (a_i,b_i) = Gab %*% f_i with f_i ~ N(0,I_k) and
		# Gab Gab' = Sab, where k = rank(Sab). Sampling the factors f
		# stays well defined even when Sab is singular.
		Sab <- (Sab + t(Sab)) / 2
		ev  <- eigen(Sab, symmetric = TRUE)
		k   <- sum(zapsmall(ev$values) > 0)
		if (k == 0) {
			cli::cli_warn(paste("k=0, no positive eigenvalues in Sab. Eigenvalues:",
						  paste(round(ev$values, 6), collapse = ", ")))
		}

		beta <- if (p > 0) NULL else numeric(0)
		ab   <- matrix(0, n, 2)

		if (k > 0) {
			Gab <- ev$vectors[, 1:k, drop = FALSE] %*%
				diag(sqrt(ev$values[1:k]), nrow = k)

			# loadings of factor f_i when actor i is the sender (send) and
			# when it is the receiver (recv) of a directed cell
			send <- wd * Gab[1, ] + wo * Gab[2, ]     # length k
			recv <- wo * Gab[1, ] + wd * Gab[2, ]     # length k

			# factor-factor precision: identical k x k block on the diagonal
			# plus a common cross block on every actor pair, closed forms
			# obtained by summing the loadings over all directed cells.
			cross <- outer(send, recv)
			Moff  <- cross + t(cross)                              # every block
			Mdiag <- n * (outer(send, send) + outer(recv, recv)) + diag(k)
			Qff <- kronecker(diag(n), Mdiag) +
				kronecker(matrix(1, n, n), Moff)

			# factor linear term: actor i contributes z_out[i]*send +
			# z_in[i]*recv, laid out actor-major
			lf <- as.vector(t(outer(z_out, send) + outer(z_in, recv)))

			if (p > 0) {
				# beta precision / linear term (whitened GLS + prior)
				lb  <- crossprod(mXw, as.vector(Zw)) + iV0 %*% m0
				Qbb <- XXw + iV0

				# beta-factor cross precision, actor-major columns in f
				Qbf <- matrix(0, p, n * k)
				for (i in seq_len(n)) {
					cols <- ((i - 1) * k + 1):(i * k)
					Qbf[, cols] <- outer(Dout[i, ], send) +
						outer(Din[i, ], recv)
				}

				Prec <- rbind(cbind(Qbb, Qbf),
							  cbind(t(Qbf), Qff))
				lin  <- c(lb, lf)
			} else {
				Prec <- Qff
				lin  <- lf
			}

			post_cov  <- solve(Prec)
			post_mean <- post_cov %*% lin
			draw <- c(rmvnorm(1, post_mean, post_cov))

			if (p > 0) {
				beta  <- draw[seq_len(p)]
				fdraw <- draw[-seq_len(p)]
			} else {
				fdraw <- draw
			}

			# recover (a_i,b_i) = Gab %*% f_i for every actor
			Fmat <- matrix(fdraw, nrow = n, ncol = k, byrow = TRUE)
			ab   <- Fmat %*% t(Gab)

		} else if (p > 0) {
			# no additive effects: whitened GLS regression on beta alone
			lb  <- crossprod(mXw, as.vector(Zw)) + iV0 %*% m0
			Qbb <- XXw + iV0
			post_cov  <- solve(Qbb)
			post_mean <- post_cov %*% lb
			beta <- c(rmvnorm(1, post_mean, post_cov))
		}

		list(beta = beta, a = ab[, 1], b = ab[, 2])
	}

####

# default prior precision for beta: a zellner g-prior with a vanishing ridge
# (a fixed fraction of the g-prior, so it shares the same 1/scale(Y)^2 units)
# that conditions the precision without ever dominating it; an absolute ridge
# would gain a spurious scale factor and crush the slopes on small-scale
# responses. the intercept is left almost unshrunk: its prior variance is
# widened to ~100x the per-cell g-prior value. XX[1,1] counts the design
# cells carrying the all-ones column, so g/XX[1,1] is that per-cell variance
# and 100x puts the intercept prior sd near ten data sd's under
# near-complete data. shared by rbeta_ab_fc() and the poisson level move in
# ame_unipartite() so both see the same prior.
.beta_prior_precision <- function(XX, mX, g) {
	iV0 <- (XX + 1e-6 * diag(diag(XX), nrow = nrow(XX))) / g
	if (all(mX[, 1] == 1)) {
		V0 <- solve(iV0)
		V0[1, 1] <- V0[1, 1] + 100 * g / XX[1, 1]
		iV0 <- solve(V0)
	}
	iV0
}
