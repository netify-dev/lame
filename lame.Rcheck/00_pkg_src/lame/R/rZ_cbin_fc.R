#' Draw latent Z for censored-binary nomination data
#'
#' Gibbs refresh of the latent normal matrix \code{Z} for the
#' censored-binary (\code{cbin}) family, where each actor emits a capped
#' set of binary nominations. Within every sending row the update honours
#' three order constraints: a nominated tie has a positive latent value,
#' every nominated tie dominates every non-nominated tie, and a
#' non-nominated tie in a row that has not spent its full nomination
#' budget (\code{odobs < odmax}) is held below zero. Reciprocal cells are
#' coupled through the dyadic correlation \code{rho}; the self-loop
#' diagonal and any missing cells are redrawn from their unconstrained
#' conditionals.
#'
#' @param Z current latent matrix (square).
#' @param EZ conditional mean of \code{Z}.
#' @param rho dyadic (within-dyad) correlation.
#' @param Y observed nomination matrix, entries in \code{{0, 1}} with
#'   \code{NA} for unobserved cells.
#' @param odmax nomination cap, a scalar or length-\code{nrow(Z)} vector.
#' @param odobs observed out-degree per actor.
#' @return the updated latent matrix, same shape and dimnames as \code{Z}.
#' @keywords internal
#' @author lame authors
rZ_cbin_fc <- function(Z, EZ, rho, Y, odmax, odobs) {
	n <- nrow(Z)
	cond_sd <- sqrt(1 - rho^2)

	# work with a copy of Y in which unobserved cells carry a sentinel
	# class, so that class masks never index with NA and every missing
	# cell is naturally excluded from the 0- and 1-classes.
	Yc <- Y
	Yc[is.na(Yc)] <- -1L

	# row / column indices of every matrix cell, used both to split the
	# off-diagonal into its two triangles and to look up per-row bounds.
	ri <- row(Z)
	ci <- col(Z)
	tri_hi <- which(ri < ci)   # (i, j), i < j
	tri_lo <- which(ri > ci)   # (i, j), i > j

	# reciprocal linear position of a set of cells: cell (i, j) -> (j, i).
	mirror <- function(pos) ci[pos] + (ri[pos] - 1L) * n

	# rows that may still gain nominations: their non-nominated ties are
	# capped at zero rather than at the weakest nominee.
	has_room <- rep(FALSE, n)
	has_room[which(odobs < odmax)] <- TRUE

	# refresh a block of cells with a N(ez, cond_sd^2) draw truncated to
	# (lo, hi). the standardized log.p-scale inverse cdf stays finite even
	# when the interval sits deep in a tail; any non-finite result retains
	# the previous latent value.
	refresh <- function(pos, lo, hi) {
		if (!length(pos)) return(invisible(NULL))
		mp <- mirror(pos)
		ez <- EZ[pos] + rho * (Z[mp] - EZ[mp])
		val <- cond_sd * rtnorm_interval_logp(ez / cond_sd,
			lo / cond_sd, hi / cond_sd)
		bad <- !is.finite(val)
		if (any(bad)) val[bad] <- Z[pos][bad]
		Z[pos] <<- val
	}

	# sweep the two nomination classes in random order; the bounds for one
	# class read the current Z of the other, so a fresh snapshot is taken
	# at the top of each class pass.
	for (cls in sample(c(0L, 1L))) {
		if (cls == 1L) {
			# nominee floor: strictly positive and above the strongest
			# non-nominee in the row.
			zt <- Z
			zt[Yc != 0L] <- -Inf
			lo_row <- pmax(0, apply(zt, 1L, max))
			hi_row <- rep(Inf, n)
		} else {
			# non-nominee ceiling: below the weakest nominee in the row,
			# and below zero for rows with unused budget.
			zt <- Z
			zt[Yc != 1L] <- Inf
			hi_row <- apply(zt, 1L, min)
			hi_row[has_room] <- 0
			lo_row <- rep(-Inf, n)
		}

		# upper triangle first, then lower, so the second block sees the
		# just-updated reciprocal values (a sequential Gibbs pass).
		hit_hi <- tri_hi[Yc[tri_hi] == cls]
		refresh(hit_hi, lo_row[ri[hit_hi]], hi_row[ri[hit_hi]])

		hit_lo <- tri_lo[Yc[tri_lo] == cls]
		refresh(hit_lo, lo_row[ri[hit_lo]], hi_row[ri[hit_lo]])
	}

	# self-loops carry the marginal dyad variance; genuinely missing cells
	# are drawn from the unit-variance prior (this second write governs any
	# NA diagonal).
	diag(Z) <- stats::rnorm(n, diag(EZ), sqrt(1 + rho))
	gap <- is.na(Y)
	if (any(gap)) Z[gap] <- stats::rnorm(sum(gap), EZ[gap], 1)

	Z
}