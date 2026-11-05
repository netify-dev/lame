# exact draw along the ridge of a row or column covariate in a directed
# unipartite fit, for every family.
#
# moving beta by delta and the row effects by -delta x_it leaves every predictor
# unchanged, so only the priors see the move. under lame()'s centring the
# direction also carries the intercept:
#   beta        += delta (e_k - m e_1)              m the mean of the period means of x
#   a_it        -= delta (x_it - xbar_t)             effect of actor i in period t
#   level_i     -= delta (xbar_i - mean_i xbar_i)    the part constant over time
#   deviation   -= the rest
# the predictor of period t then moves by delta r_t, r_t = xbar_t - m, and where
# r_t is not zero the data-constrained cells keep their latent value. every term
# is Gaussian in delta, so the conditional is drawn in one step. static
# coefficients only; a static additive effect needs a time-constant covariate.
# the bipartite ame() case is in R/bipartite_cov_ridge.R.

# AR(1) precision matrix of a stationary path of length N
.lame_cov_ridge_ar_precision <- function(N, rho, sigma) {
	if (N == 1L) return(matrix((1 - rho^2) / sigma^2, 1L, 1L))
	Q <- matrix(0, N, N)
	d <- rep(1 + rho^2, N)
	d[c(1L, N)] <- 1
	diag(Q) <- d
	Q[cbind(seq_len(N - 1L), 2:N)] <- -rho
	Q[cbind(2:N, seq_len(N - 1L))] <- -rho
	Q / sigma^2
}

# row or column covariates of the design: for each coefficient k not in `skip`,
# the type ("row", "col" or NULL) and the n x N matrix of its values
.lame_cov_ridge_types <- function(Xlist, skip = integer(0)) {
	p <- dim(Xlist[[1L]])[3L]
	n <- dim(Xlist[[1L]])[1L]
	N <- length(Xlist)
	out <- vector("list", p)
	for (k in setdiff(seq_len(p), skip)) {
		xr <- matrix(0, n, N)
		xc <- matrix(0, n, N)
		is_row <- is_col <- TRUE
		for (t in seq_len(N)) {
			Xk <- Xlist[[t]][, , k]
			if (!all(is.finite(Xk))) { is_row <- is_col <- FALSE; break }
			if (is_row && max(abs(Xk - Xk[, 1L])) > 1e-10) is_row <- FALSE
			if (is_col && max(abs(t(Xk) - Xk[1L, ])) > 1e-10) is_col <- FALSE
			xr[, t] <- Xk[, 1L]
			xc[, t] <- Xk[1L, ]
			if (!is_row && !is_col) break
		}
		if (is_row && is_col) next
		if (is_row) out[[k]] <- list(type = "row", x = xr)
		else if (is_col) out[[k]] <- list(type = "col", x = xc)
	}
	out
}

# data term of one period: quadratic and linear coefficients of the Gaussian
# dyad density in r for the data-constrained cells `M`, the other cells moving
# with the predictor. `E` the residuals Z - EZ
.lame_cov_ridge_data <- function(E, M, rho, s2) {
	iu <- which(upper.tri(E))
	# il: the transposed partner of each cell of iu, in the same order
	nn <- nrow(E)
	ij <- arrayInd(iu, c(nn, nn))
	il <- ij[, 2L] + (ij[, 1L] - 1L) * nn
	m1 <- M[iu]; m2 <- M[il]
	e1 <- E[iu]; e2 <- E[il]
	e1[!is.finite(e1)] <- 0; e2[!is.finite(e2)] <- 0
	den <- s2 * (1 - rho^2)
	c(q = sum(m1 + m2 - 2 * rho * m1 * m2) / den, c = sum(m1 * (e1 - rho * e2) + m2 * (e2 - rho * e1)) / den)
}

# one visit of every row and column covariate of the design
#   Xlist   the design, one n x n x p array per period
#   beta    the coefficients (intercept first when `center`)
#   P0      the coefficient prior precision
#   a_lev, b_lev, a_mat, b_mat   levels and totals of the dynamic effects (a_mat /
#           b_mat NULL for static effects, a_lev / b_lev then the effects)
#   Sab     covariance of the pair of levels; rho_ab, sigma_ab of the deviations
#   rvar, cvar   which margins exist
#   center  lame() (centred effects, intercept takes the constant) or ame()
#   E, M, rho_t, s2   latent residuals, data-constrained cells (NULL: none),
#           per-period dyadic correlation and variance, for the data term
# returns the updated state, the per-period predictor shift and the counters
.lame_cov_ridge_visit <- function(
	Xlist,
	beta,
	P0,
	a_lev,
	b_lev,
	a_mat,
	b_mat,
	Sab,
	rho_ab,
	sigma_ab,
	rvar,
	cvar,
	center = TRUE,
	E = NULL,
	M = NULL,
	rho_t = NULL,
	s2 = 1,
	skip = integer(0)
) {
	n <- dim(Xlist[[1L]])[1L]
	N <- length(Xlist)
	dyn <- !is.null(a_mat)
	stats <- c(ridge_calls = 0, ridge_moved = 0)
	shift_t <- numeric(N)
	typ <- .lame_cov_ridge_types(Xlist, skip)
	Pab <- tryCatch(solve(Sab), error = function(e) NULL)
	if (is.null(Pab)) return(
		list(beta = beta, a_lev = a_lev, b_lev = b_lev, a_mat = a_mat, b_mat = b_mat, shift = shift_t, stats = stats)
	)
	Qar <- if (dyn) .lame_cov_ridge_ar_precision(N, rho_ab, sigma_ab) else NULL
	ks <- which(!vapply(typ, is.null, logical(1L)))
	if (!length(ks)) return(
		list(beta = beta, a_lev = a_lev, b_lev = b_lev, a_mat = a_mat, b_mat = b_mat, shift = shift_t, stats = stats)
	)
	for (k in ks[sample.int(length(ks))]) {
		ty <- typ[[k]]
		row <- identical(ty$type, "row")
		if (row && !isTRUE(rvar)) next
		if (!row && !isTRUE(cvar)) next
		x <- ty$x
		if (!dyn && any(abs(x - x[, 1L]) > 1e-10)) next
		xt <- colMeans(x)
		m <- if (center) mean(xt) else 0
		r <- if (center) xt - m else numeric(N)
		xbar <- rowMeans(x)
		L <- if (center) xbar - mean(xbar) else xbar
		dB <- numeric(length(beta))
		dB[k] <- 1
		if (center) dB[1L] <- -m
		lev <- if (row) a_lev else b_lev
		oth <- if (row) b_lev else a_lev
		both <- isTRUE(rvar) && isTRUE(cvar)
		p_ll <- if (!both) 1 / (if (row) Sab[1L, 1L] else Sab[2L, 2L]) else if (row) Pab[1L, 1L] else Pab[2L, 2L]
		p_lo <- if (both) Pab[1L, 2L] else 0
		# log density along the line: -(Q delta^2 - 2 C delta) / 2
		Q <- drop(crossprod(dB, P0 %*% dB)) + p_ll * sum(L^2)
		C <- -drop(crossprod(dB, P0 %*% beta)) + p_ll * sum(lev * L) + p_lo * sum(oth * L)
		if (dyn) {
			tot <- if (row) a_mat else b_mat
			dev <- tot - lev
			S <- if (center) x - matrix(xt, n, N, byrow = TRUE) else x
			D <- S - L
			Q <- Q + sum(D * (D %*% Qar))
			C <- C + sum(dev * (D %*% Qar))
		}
		if (center && any(abs(r) > 1e-10)) {
			if (is.null(E)) next
			for (t in seq_len(N)) {
				if (abs(r[t]) <= 1e-10) next
				dd <- .lame_cov_ridge_data(
					E[, , t],
					if (is.null(M)) array(FALSE, dim(E)[1:2]) else M[, , t],
					rho_t[t],
					s2
				)
				Q <- Q + r[t]^2 * dd[["q"]]
				C <- C + r[t] * dd[["c"]]
			}
		}
		if (!is.finite(Q) || Q <= 0 || !is.finite(C)) next
		delta <- stats::rnorm(1L, C / Q, 1 / sqrt(Q))
		stats["ridge_calls"] <- stats["ridge_calls"] + 1
		beta <- beta + delta * dB
		lev <- lev - delta * L
		if (row) a_lev <- lev else b_lev <- lev
		if (dyn) {
			if (row) a_mat <- a_mat - delta * S else b_mat <- b_mat - delta * S
		}
		shift_t <- shift_t + delta * r
		stats["ridge_moved"] <- stats["ridge_moved"] + 1
		# constrained cells keep their latent value, so their residual moves
		if (!is.null(E) && !is.null(M)) for (t in seq_len(N)) if (abs(r[t]) > 1e-10) E[, , t][M[, , t]] <- E[, , t][M[, , t]] - delta * r[t]
	}
	list(beta = beta, a_lev = a_lev, b_lev = b_lev, a_mat = a_mat, b_mat = b_mat, shift = shift_t, stats = stats)
}
