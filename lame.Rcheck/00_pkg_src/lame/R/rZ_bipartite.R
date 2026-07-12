# bipartite (rectangular) versions of the family z-sampling helpers.
#
# the unipartite-asymmetric samplers in rz_<family>_fc.r hardcode the
# square n x n assumption: they split z into upper/lower triangles,
# couple each (i,j) with its reciprocal (j,i) through `rho` and a
# transposed term `t(z)`, and sample the diagonal separately. none of
# that applies bipartite: y is na x nb rectangular, there is no
# "reciprocal cell" (no t(z) term), there is no diagonal, and rho is
# always 0 (the bipartite path drops dyadic correlation).
#
# each rectangular sampler below mirrors its unipartite sibling's
# constraint structure but drops the triangle / diagonal / rho machinery.
# conditional draws are independent truncated normals once the row-level
# rank-cone constraints (for cbin/frn) have been resolved.
#
# the cutpoint convention for `ordinal` matches the existing unipartite
# helper: cutpoints are *data-induced* (the boundaries are the order
# statistics of z within adjacent rank classes), not separately
# parameterised. this keeps the bipartite path consistent with the
# square ordinal sampler the rest of the codebase assumes.

#' Exact truncated-normal draw on (lb, ub) via the log.p-scale inverse CDF
#'
#' Draws \code{z ~ N(ez, 1)} truncated to \code{(lb, ub)} for each cell,
#' working entirely on the log-probability scale in the tail that
#' contains the truncation interval. This stays exact for arbitrarily
#' extreme \code{ez} / bounds (log-scale \code{pnorm}/\code{qnorm} are
#' accurate to \code{|ez|} of several hundred), so no probability clamp
#' is needed. Infinite bounds reduce to the
#' one-sided case automatically.
#'
#' @param ez vector of conditional means.
#' @param lb lower bounds (recycled to \code{length(ez)}).
#' @param ub upper bounds (recycled to \code{length(ez)}).
#' @return vector of truncated-normal draws, same length as \code{ez}.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rtnorm_interval_logp <- function(ez, lb, ub) {
	m <- length(ez)
	lb <- rep_len(lb, m)
	ub <- rep_len(ub, m)
	z <- numeric(m)
	# guard the open unit interval (r's runif never returns 0 or 1, but
	# the clamp costs nothing and removes the theoretical log(0) case)
	u <- pmin(pmax(stats::runif(m), .Machine$double.xmin), 1 - 1e-16)
	# non-finite ez falls through to the upper-tail branch and yields
	# nan instead of erroring in the branch test
	lo <- !is.na(ub - ez) & ((ub - ez) <= 0)
	if (any(lo)) {
		# interval in the lower tail: work with lower-tail log masses
		la <- stats::pnorm(lb[lo] - ez[lo], log.p = TRUE)
		lc <- stats::pnorm(ub[lo] - ez[lo], log.p = TRUE)
		r <- exp(la - lc)
		z[lo] <- ez[lo] + stats::qnorm(lc + log(u[lo] + (1 - u[lo]) * r),
			log.p = TRUE)
	}
	if (any(!lo)) {
		# interval reaching into the upper tail: mirrored survival form
		lqa <- stats::pnorm(lb[!lo] - ez[!lo], lower.tail = FALSE,
			log.p = TRUE)
		lqb <- stats::pnorm(ub[!lo] - ez[!lo], lower.tail = FALSE,
			log.p = TRUE)
		r <- exp(lqb - lqa)
		z[!lo] <- ez[!lo] + stats::qnorm(lqa + log((1 - u[!lo]) + u[!lo] * r),
			lower.tail = FALSE, log.p = TRUE)
	}
	z
}

#' Sample Z under bipartite Poisson (rectangular MH step)
#'
#' Rectangular analogue of \code{\link{rZ_pois_fc}}. Each \code{(i,j)}
#' cell is updated independently with a Metropolis-Hastings step on
#' the Poisson log-link; there is no upper/lower-triangle coupling and
#' no diagonal because the bipartite Y is \code{nA x nB} with disjoint
#' row and column actor sets.
#'
#' @param Z current rectangular latent matrix (nA x nB).
#' @param EZ expected value (regression + additive + multiplicative
#'   effects), nA x nB.
#' @param s2 dyadic variance (proposal scale on the latent log-mean).
#' @param Y observed count matrix (nA x nB); \code{NA} entries are
#'   resampled from the prior.
#' @param log_exposure optional per-period log-exposure offset (scalar
#'   or matrix conformable with Y). Default \code{0}.
#' @return updated \code{nA x nB} latent matrix.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_pois_bip_fc <- function(Z, EZ, s2, Y, log_exposure = 0) {
	nA <- nrow(Z); nB <- ncol(Z)
	# propose
	Zp <- Z + matrix(stats::rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
	# log-prior under n(ez, s2)
	lp_curr <- -0.5 / s2 * (Z  - EZ)^2
	lp_prop <- -0.5 / s2 * (Zp - EZ)^2
	# log-likelihood: poisson(exp(z + log_exposure)); guarded against
	# overflow at exp(z) when z is large
	obs <- !is.na(Y)
	mu_curr <- exp(pmin(Z  + log_exposure, 30))
	mu_prop <- exp(pmin(Zp + log_exposure, 30))
	ll_curr <- matrix(0, nA, nB)
	ll_prop <- matrix(0, nA, nB)
	ll_curr[obs] <- stats::dpois(Y[obs], mu_curr[obs], log = TRUE)
	ll_prop[obs] <- stats::dpois(Y[obs], mu_prop[obs], log = TRUE)
	# accept ratio
	lr <- (lp_prop + ll_prop) - (lp_curr + ll_curr)
	lh <- log(matrix(stats::runif(nA * nB), nA, nB))
	accept <- lr > lh
	Z[accept] <- Zp[accept]
	# resample na cells from the prior (parity with rz_pois_fc_cpp behaviour
	# for missing cells; uninformative draws keep ez-conditional bookkeeping
	# valid downstream)
	if (any(!obs)) {
		Z[!obs] <- stats::rnorm(sum(!obs), EZ[!obs], sqrt(s2))
	}
	Z
}


#' Sample Z under bipartite ordinal data (rectangular)
#'
#' Rectangular analogue of \code{\link{rZ_ord_fc}}. Cutpoints are
#' data-induced (the boundary between rank \code{w} and rank \code{w+1}
#' is the maximum Z in rank \code{w} and minimum Z in rank \code{w+1}),
#' matching the existing unipartite convention. With rho = 0 in
#' bipartite, every cell at rank \code{w} is drawn independently from
#' a truncated normal centered at \code{EZ}.
#'
#' @param Z current rectangular latent matrix (nA x nB).
#' @param EZ expected value (nA x nB).
#' @param Y observed ordinal matrix (nA x nB).
#' @return updated nA x nB latent matrix.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_ord_bip_fc <- function(Z, EZ, Y) {
	uY <- sort(unique(c(Y)))
	nA <- nrow(Y); nB <- ncol(Y)
	W <- matrix(match(Y, uY), nA, nB)
	W[is.na(W)] <- -1L

	for (w in sample(c(-1L, seq_along(uY)))) {
		lb <- suppressWarnings(max(Z[!is.na(W) & W == w - 1L], na.rm = TRUE))
		ub <- suppressWarnings(min(Z[!is.na(W) & W == w + 1L], na.rm = TRUE))
		up <- W == w
		if (!any(up)) next
		ez <- EZ[up]
		# tail-stable truncated draw: the plain pnorm/qnorm form emits
		# +/-inf once a bound saturates the normal cdf
		Z[up] <- rtnorm_interval_logp(ez, lb, ub)
	}
	Z
}


#' Sample Z under bipartite censored binary nominations (rectangular)
#'
#' Rectangular analogue of \code{\link{rZ_cbin_fc}}. Row-wise
#' constraints: \code{Y_ij = 1 => Z_ij > 0}; \code{Y_ij = 0 & odobs_i <
#' odmax_i => Z_ij < 0}; nominated alters dominate non-nominated alters
#' within the row. With rho = 0 there is no reciprocity term.
#'
#' @param Z current rectangular latent matrix (nA x nB).
#' @param EZ expected value (nA x nB).
#' @param Y observed 0/1 matrix (nA x nB); \code{NA} entries resampled
#'   from the unconstrained prior.
#' @param odmax row-wise maximum nominations (length nA or scalar).
#' @param odobs row-wise observed outdegree (length nA).
#' @return updated nA x nB latent matrix.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_cbin_bip_fc <- function(Z, EZ, Y, odmax, odobs) {
	nA <- nrow(Z); nB <- ncol(Z)
	for (y in sample(0:1)) {
		if (y == 1L) {
			ub_mat <- matrix(Inf, nA, nB)
			# the lower bound for a y=1 cell is the max z among y=0 cells in the
			# same row (so nominees dominate non-nominees), with floor 0. rows
			# with no y=0 cells yield max(numeric(0)) = -inf, then the pmax(., 0)
			# floor makes lb = 0 anyway, which is what we want
			z_zeros <- Z - (Y != 0) * (Inf ^ (Y != 0))
			lbm_row <- suppressWarnings(apply(z_zeros, 1, max, na.rm = TRUE))
			lb_mat <- matrix(pmax(lbm_row, 0), nA, nB)
		} else {
			lb_mat <- matrix(-Inf, nA, nB)
			z_ones <- Z + (Y != 1) * (Inf ^ (Y != 1))
			ubm_row <- suppressWarnings(apply(z_ones, 1, min, na.rm = TRUE))
			ub_mat <- matrix(ubm_row, nA, nB)
			# saturated rows: non-nominees cap below the min nominee z;
			# unsaturated rows: non-nominees are < 0. defensively guard
			# against na from rows where odmax or odobs is missing.
			unsat <- odobs < odmax
			unsat[is.na(unsat)] <- FALSE
			if (any(unsat)) {
				ub_mat[unsat, ] <- pmin(ub_mat[unsat, , drop = FALSE], 0)
			}
		}
		up <- Y == y
		if (!any(up, na.rm = TRUE)) next
		up[is.na(up)] <- FALSE
		ez <- EZ[up]
		lb <- lb_mat[up]; ub <- ub_mat[up]
		# guard against degenerate (lb >= ub) caused by inf/-inf edge cases
		ok <- is.finite(ez) & (ub > lb)
		if (any(ok)) {
			z_curr <- Z[up]
			# exact truncated-normal draw on (lb, ub) via the log.p-scale
			# inverse cdf: tail-stable for bounds arbitrarily far from ez.
			z_curr[ok] <- rtnorm_interval_logp(ez[ok], lb[ok], ub[ok])
			Z[up] <- z_curr
		}
	}
	# na cells: unconstrained draw from the prior
	if (any(is.na(Y))) {
		miss <- is.na(Y)
		Z[miss] <- stats::rnorm(sum(miss), EZ[miss], 1)
	}
	Z
}


#' Sample Z under bipartite fixed-rank nominations (rectangular)
#'
#' Rectangular analogue of \code{\link{rZ_frn_fc}}. Row-wise rank
#' constraints: ranked nominees are positive and ordered by rank;
#' non-nominees fall below all ranked cells in the row (and below 0
#' when the row is unsaturated). With rho = 0 the unipartite triangle
#' coupling drops out.
#'
#' \strong{Monotonicity convention (important for simulation).} The
#' sampler enforces \emph{higher Y -> higher Z}: when \code{Y[i, j] >
#' Y[i, k]}, the corresponding latent values satisfy \code{Z[i, j] >
#' Z[i, k]}. When simulating data for recovery tests, build Y via
#' \code{rank(Z[i, ])} (low rank = small Y), not \code{rank(-Z[i, ])}.
#' the inverted convention identifies the negative of the true beta.
#'
#' @param Z current rectangular latent matrix (nA x nB).
#' @param EZ expected value (nA x nB).
#' @param Y observed rank matrix (nA x nB); 0 = unranked,
#'   1..ncol(YL) = rank order.
#' @param YL per-row ranked-receiver list (nA x max_rank).
#' @param odmax row-wise max nominations.
#' @param odobs row-wise observed outdegree.
#' @return updated nA x nB latent matrix.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @keywords internal
rZ_frn_bip_fc <- function(Z, EZ, Y, YL, odmax, odobs) {
	nA <- nrow(Z); nB <- ncol(Z)
	Y2 <- Y; Y2[is.na(Y2)] <- -1L
	max_rank <- ncol(YL)

	for (y in sample(c((-1L):max_rank))) {
		if (y < 2L) {
			if (y <= 0L) {
				lbm <- rep(-Inf, nA)
			} else {
				# y == 1: lower bound for the top-ranked cell is the max z among
				# non-nominees in the row, floored at 0. rows with no non-nominees
				# yield max(numeric(0)) = -inf which the pmax(., 0) lifts to 0.
				lbm <- pmax(0, suppressWarnings(apply(Z - (Y2 != 0) * (Inf ^ (Y2 != 0)), 1, max,
				                                       na.rm = TRUE)))
			}
		} else {
			lbm <- Z[cbind(seq_len(nA), YL[, y - 1L])]
		}
		if (y == -1L) {
			ubm <- rep(Inf, nA)
		} else if (y < max_rank && y >= 0L) {
			ubm <- Z[cbind(seq_len(nA), YL[, y + 1L])]
		} else {
			ubm <- rep(Inf, nA)
		}
		if (y == 0L) {
			unsat <- odobs < odmax
			unsat[is.na(unsat)] <- FALSE
			ubm[unsat] <- 0
		}
		ubm[is.na(ubm)] <- Inf
		lbm[is.na(lbm)] <- -Inf

		up <- Y2 == y
		if (!any(up)) next
		# row index of each true cell, then map to its row-level bounds
		row_idx <- row(Z)[up]
		lb <- lbm[row_idx]; ub <- ubm[row_idx]
		ez <- EZ[up]
		ok <- is.finite(ez) & (ub > lb)
		if (any(ok)) {
			z_curr <- Z[up]
			# exact truncated-normal draw on (lb, ub) via the log.p-scale
			# inverse cdf: tail-stable for bounds arbitrarily far from ez.
			z_curr[ok] <- rtnorm_interval_logp(ez[ok], lb[ok], ub[ok])
			Z[up] <- z_curr
		}
	}
	Z
}
