#' Full-conditional latent draw for the ordinal family
#'
#' Refreshes the latent sociomatrix \code{Z} for an ordinal outcome. The
#' observed ordinal categories in \code{Y} impose only rank information on
#' \code{Z}: every latent value tied to category \eqn{c} must sit above all
#' latent values of the next-lower category and below all latent values of the
#' next-higher category. Within that rank-implied window each entry is redrawn
#' from its dyadic full conditional, a truncated normal whose mean borrows the
#' transpose partner through the within-dyad correlation \code{rho} and whose
#' variance is \eqn{1-\rho^2}. Missing dyads carry no rank constraint and are
#' drawn from the untruncated conditional. The diagonal is nuisance and is
#' refreshed from an independent normal.
#'
#' @param Z square latent matrix, current state.
#' @param EZ square matrix of latent means.
#' @param rho scalar within-dyad correlation.
#' @param Y square matrix of ordinal outcomes (may contain \code{NA}).
#' @return the updated square latent matrix \code{Z}.
#' @keywords internal
#' @author lame authors
rZ_ord_fc <- function(Z, EZ, rho, Y) {
	n <- nrow(Z)
	# ordered distinct categories, and a 1..K category index per cell (NA kept)
	ranks <- sort(unique(c(Y)))
	n_cat <- length(ranks)
	cat_of <- matrix(match(Y, ranks), n, n)

	cond_sd <- sqrt(1 - rho^2)
	halves <- list(upper.tri(Z), lower.tri(Z))

	# visit the missing group (0) and each category in a random order so the
	# rank windows are read from a freshly stirred configuration each sweep
	for (g in sample(c(0L, seq_len(n_cat)))) {
		if (g == 0L) {
			pick <- is.na(cat_of)
			floor_z <- -Inf
			ceil_z <- Inf
		} else {
			pick <- !is.na(cat_of) & cat_of == g
			under <- !is.na(cat_of) & cat_of == g - 1L
			over <- !is.na(cat_of) & cat_of == g + 1L
			floor_z <- if (any(under)) max(Z[under]) else -Inf
			ceil_z <- if (any(over)) min(Z[over]) else Inf
		}
		# refresh the two triangles separately: updating one triangle changes
		# the transpose-based conditional mean seen by the other
		for (mask in halves) {
			cells <- pick & mask
			if (!any(cells)) next
			cmean <- EZ[cells] + rho * (t(Z)[cells] - t(EZ)[cells])
			# work on the unit-variance scale, then rescale by cond_sd
			draw <- cond_sd * rtnorm_interval_logp(cmean / cond_sd,
				floor_z / cond_sd, ceil_z / cond_sd)
			# defend against a non-finite draw from an extreme window by
			# holding the current value in place
			stuck <- !is.finite(draw)
			if (any(stuck)) draw[stuck] <- Z[cells][stuck]
			Z[cells] <- draw
		}
	}
	diag(Z) <- rnorm(n, diag(EZ), sqrt(1 + rho))
	Z
}