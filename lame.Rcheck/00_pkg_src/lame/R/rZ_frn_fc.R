#' Draw the latent Z matrix for fixed-rank-nomination data
#'
#' Gibbs update of the latent sociomatrix \code{Z} for a fixed-rank-nomination
#' (FRN) outcome. Each row's latent values are drawn from truncated normals
#' whose bounds encode three rank constraints: nominated ties outrank
#' non-nominations, higher observed ranks map to larger latents, and a
#' non-nomination made with spare out-degree capacity must stay negative.
#'
#' The constraints realised are (1) \code{Y[i,j] > Y[i,k]} implies
#' \code{Z[i,j] > Z[i,k]}, (2) \code{Y[i,j] > 0} implies \code{Z[i,j] > 0},
#' and (3) \code{Y[i,j] == 0} with \code{odobs[i] < odmax[i]} implies
#' \code{Z[i,j] < 0}.
#'
#' @param Z current latent sociomatrix (square).
#' @param EZ conditional mean matrix for \code{Z}.
#' @param rho within-dyad correlation.
#' @param Y square matrix of observed ranked nominations (0 = no tie, NA on
#'   the diagonal / missing).
#' @param YL matrix whose \code{r}-th column gives, per row, the column index
#'   of the individual holding rank \code{r} (ascending preference).
#' @param odmax scalar or per-row maximum number of nominations allowed.
#' @param odobs per-row observed out-degree.
#' @return the updated square latent matrix \code{Z}.
#' @keywords internal
#' @author lame authors
#' @export rZ_frn_fc
rZ_frn_fc <- function(Z, EZ, rho, Y, YL, odmax, odobs) {

	n <- nrow(Z)
	cond_sd <- sqrt(1 - rho^2)
	row_of <- row(Z)                 # entry [i, j] carries the value i
	upper <- upper.tri(Z)
	lower <- lower.tri(Z)
	max_rank <- ncol(YL)
	rows <- seq_len(n)

	# recode the diagonal / missing cells to the sentinel rank -1 so they can
	# be swept through unconstrained alongside the genuine rank strata
	Yr <- Y
	Yr[is.na(Yr)] <- -1

	# per-row latent value of whoever currently holds rank `r`
	held_at <- function(r) Z[cbind(rows, YL[, r])]

	# lower/upper truncation limits (one value per row) for a rank stratum `y`
	limits_for <- function(y) {
		if (y <= 0) {
			lo <- rep.int(-Inf, n)
		} else if (y == 1) {
			# the least-preferred nominee must sit above 0 and above every
			# non-nominated latent in its row
			pool <- Z
			pool[Yr != 0] <- -Inf
			lo <- pmax(0, apply(pool, 1, max, na.rm = TRUE))
		} else {
			lo <- held_at(y - 1)
		}

		if (y == -1 || y == max_rank) {
			hi <- rep.int(Inf, n)
		} else {
			hi <- held_at(y + 1)
			if (y == 0) {
				# unused nomination slots pin the ceiling at zero
				hi[odobs < odmax] <- 0
			}
		}

		lo[is.na(lo)] <- -Inf
		hi[is.na(hi)] <- Inf
		list(lo = lo, hi = hi)
	}

	# redraw every entry of one triangle that belongs to stratum `y`
	sweep_triangle <- function(mask, y, lo, hi) {
		cells <- mask & (Yr == y)
		if (!any(cells)) {
			return(invisible(NULL))
		}
		r <- row_of[cells]
		mu <- EZ[cells] + rho * (t(Z)[cells] - t(EZ)[cells])
		vals <- cond_sd * rtnorm_interval_logp(mu / cond_sd,
			lo[r] / cond_sd, hi[r] / cond_sd)
		spoiled <- !is.finite(vals)
		if (any(spoiled)) {
			vals[spoiled] <- Z[cells][spoiled]
		}
		Z[cells] <<- vals
		invisible(NULL)
	}

	# visit strata and triangles in random order so the Gibbs sweep does not
	# privilege any fixed traversal
	for (y in sample((-1):max_rank)) {
		lim <- limits_for(y)
		for (side in sample(c("U", "L"))) {
			sweep_triangle(if (side == "U") upper else lower,
				y, lim$lo, lim$hi)
		}
	}

	diag(Z) <- rnorm(n, diag(EZ), sqrt(1 + rho))
	Z
}