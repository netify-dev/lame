# joint Gaussian full conditional of the coefficients and additive effects for replicated
# unipartite data; the design-only pieces are cached while the same design is passed
.lame_gaussian_design_cache <- new.env(parent = emptyenv())

.lame_replicated_design <- function(mX, n, p, nt) {
	if (!p) return(replicated_gaussian_design_cpp(array(0, c(n * n, 0L, nt)), n))
	cache <- .lame_gaussian_design_cache
	if (!is.null(cache$design) && identical(cache$key, mX)) return(cache$design)
	design <- replicated_gaussian_design_cpp(array(as.numeric(mX), c(n * n, p, nt)), n)
	cache$key <- mX
	cache$design <- design
	design
}

# default coefficient prior: the replicated g-prior N(0, (sum_t X_t'X_t / (T g))^-1)
.lame_replicated_prior <- function(mX, XX, g, p, nt) {
	if (!p) return(matrix(0, 0L, 0L))
	if (is.null(XX)) {
		XX <- array(0, c(p, p, nt))
		for (tt in seq_len(nt)) XX[, , tt] <- crossprod(matrix(mX[, , tt], ncol = p))
	}
	if (!identical(dim(XX), as.integer(c(p, p, nt))) || any(!is.finite(XX)))
		stop("Invalid replicated coefficient prior statistics.")
	P <- apply(XX, c(1L, 2L), sum) / (nt * g)
	P / 2 + t(P) / 2
}

.lame_replicated_gaussian_args <- function(ZT, mX, Sab, rho, s2, g, XX) {
	n <- dim(ZT)[1L]; nt <- dim(ZT)[3L]; p <- dim(mX)[2L]
	if (is.null(nt) || dim(ZT)[2L] != n || !identical(dim(mX), c(n * n, p, nt)))
		stop("Invalid replicated Gaussian response or design dimensions.")
	if (length(g) != 1L || !is.finite(g) || g <= 0)
		stop("Coefficient prior scale g must be positive and finite.")
	list(
		design = .lame_replicated_design(mX, n, p, nt),
		P0 = .lame_replicated_prior(mX, XX, g, p, nt),
		m0 = numeric(p),
		rho = as.numeric(rho)
	)
}

.lame_replicated_gaussian_draw <- function(
	ZT,
	mX,
	Sab,
	rho,
	s2,
	g,
	center = FALSE,
	XX = NULL
) {
	x <- .lame_replicated_gaussian_args(ZT, mX, Sab, rho, s2, g, XX)
	ans <- replicated_gaussian_draw_cpp(
		x$design,
		ZT,
		Sab,
		x$rho,
		s2,
		x$P0,
		x$m0,
		isTRUE(center)
	)
	lapply(ans, as.numeric)
}

# exact mean and covariance of (beta, a, b) under the same conditional
.lame_replicated_gaussian_moments <- function(
	ZT,
	mX,
	Sab,
	rho,
	s2,
	g,
	center = FALSE,
	XX = NULL
) {
	x <- .lame_replicated_gaussian_args(ZT, mX, Sab, rho, s2, g, XX)
	replicated_gaussian_moments_cpp(
		x$design,
		ZT,
		Sab,
		x$rho,
		s2,
		x$P0,
		x$m0,
		isTRUE(center)
	)
}
