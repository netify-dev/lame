#' Joint Gibbs update of regression and additive effects (single relation)
#'
#' Draws jointly from the full conditional of the regression coefficients
#' \code{beta} together with the additive sender/receiver effects
#' \code{a}, \code{b} for a single social-relations regression matrix. The
#' within-dyad reciprocity covariance is whitened through the symmetric and
#' antisymmetric parts of the relation, the additive-effect prior
#' \code{N(0,Sab)} is written through a square-root factor so a
#' rank-deficient \code{Sab} is handled without inverting it, and the actor
#' effects are eliminated block by block over an orthonormal Helmert basis.
#' Regression coefficients are drawn from their marginal Gaussian
#' conditional, followed by the conditionally independent actor coordinates.
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
rbeta_ab_fc <- function(
	Z,
	Sab,
	rho,
	X = NULL,
	s2 = 1,
	offset = 0,
	iV0 = NULL,
	m0 = NULL,
	g = length(Z)
) {
	Z <- Z - offset
	if (!is.matrix(Z) || nrow(Z) != ncol(Z) || nrow(Z) < 1L)
		cli::cli_abort("Gaussian conditional response must be a nonempty square matrix.")
	n <- nrow(Z)
	if (!is.matrix(Sab) || !identical(dim(Sab), c(2L, 2L)) || any(!is.finite(Sab)))
		cli::cli_abort("Sab must be a finite two by two covariance matrix.")
	if (is.null(X)) X <- design_array(intercept = FALSE, n = n)
	if (is.null(attributes(X)$XX)) {
		X <- precomputeX(X)
		cli::cli_warn(paste0(
			"Summary statistics of X are not precomputed. ",
			"Run X<-precomputeX(X) to speed up calculations."
		))
	}
	mX <- attributes(X)$mX
	mXt <- attributes(X)$mXt
	XX <- attributes(X)$XX
	p <- dim(X)[3L]
	if (length(dim(X)) != 3L || !identical(dim(X)[1:2], c(n, n)) ||
			!identical(dim(mX), c(n * n, p)) || !identical(dim(mXt), c(n * n, p)))
		cli::cli_abort("Gaussian conditional design dimensions must match the response.")
	if (p > 0L && is.null(iV0)) iV0 <- .beta_prior_precision(XX, mX, g)
	if (p == 0L) iV0 <- matrix(0, 0L, 0L)
	if (is.null(m0)) m0 <- rep(0, p)
	if (length(m0) != p || any(!is.finite(m0)) ||
			(p > 0L && (!is.matrix(iV0) || !identical(dim(iV0), c(p, p)))))
		cli::cli_abort("Gaussian coefficient prior dimensions and mean must match the design.")
	if (length(s2) != 1L || !is.finite(s2) || s2 <= 0 ||
			length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
		cli::cli_abort("Residual variance must be positive and rho must lie strictly between -1 and 1.")
	if (any(is.infinite(Z)) || any(!is.finite(mX)) || any(!is.finite(mXt)))
		cli::cli_abort("Gaussian conditional residuals and design must be finite apart from missing residuals.")
	# missing dyads and their reciprocals enter as zero on the complete grid
	Z[is.na(Z) | is.na(t(Z))] <- 0
	out <- replicated_gaussian_draw_cpp(
		.lame_replicated_design(mX, n, p, 1L),
		array(Z, c(n, n, 1L)),
		Sab,
		rho,
		s2,
		iV0,
		m0,
		FALSE
	)
	list(beta = as.numeric(out$beta), a = as.numeric(out$a), b = as.numeric(out$b))
}

# default beta prior precision: a Zellner g-prior plus a ridge proportional to diag(XX),
# so it keeps the 1/scale(Y)^2 units; the intercept prior variance is widened by
# 100 per-cell g-prior variances. shared by rbeta_ab_fc() and ame_unipartite()
.beta_prior_precision <- function(XX, mX, g) {
	iV0 <- (XX + 1e-6 * diag(diag(XX), nrow = nrow(XX))) / g
	if (all(mX[, 1] == 1)) {
		V0 <- solve(iV0)
		V0[1, 1] <- V0[1, 1] + 100 * g / XX[1, 1]
		iV0 <- solve(V0)
	}
	iV0
}
