#' Full-conditional draw of the dyadic variance for replicated relational data
#'
#' Samples the dyadic variance \code{s2} from its Gibbs full conditional when
#' the residual cube stacks several independent replicate networks that share a
#' common within-dyad correlation \code{rho}. Each ordered dyad pair is whitened
#' by the symmetric inverse square root of the 2x2 correlation matrix and the
#' pooled whitened sum of squares drives an inverse-gamma draw.
#'
#' @param E.T numeric array of dimension n x n x N holding the residual
#'   matrices along the third margin (one square slice per replicate).
#' @param rho current value of the within-dyad correlation.
#' @return a single numeric draw of \code{s2}.
#' @keywords internal
#' @author lame authors
#' @export rs2_rep_fc
rs2_rep_fc <-
	function(E.T, rho) {
		n_rep <- dim(E.T)[3]

		# whitening map: symmetric inverse square root of the dyad correlation
		# matrix, so a whitened (e_ij, e_ji) pair is unit-scaled and
		# uncorrelated under the model
		corr <- matrix(c(1, rho, rho, 1), 2, 2)
		whiten <- mhalf(solve(corr))

		ssq <- 0				# pooled whitened sum of squares
		count <- 0L				# pooled count of whitened scalars

		for (r in seq_len(n_rep)) {
			slice <- E.T[, , r]
			above <- upper.tri(slice)
			# stack the two ordered residuals of each dyad, then whiten
			dyad_pairs <- cbind(slice[above], t(slice)[above])
			whitened <- dyad_pairs %*% whiten
			ssq <- ssq + sum(whitened * whitened)
			count <- count + length(whitened)
		}

		# a single pseudo-observation placed at the empirical mean square keeps
		# the draw invariant to the overall scale of the residuals
		mean_sq <- if (count > 0) ssq / count else 1
		if (!is.finite(mean_sq) || mean_sq <= 0) mean_sq <- 1

		shape <- (count + 1) / 2
		rate <- (ssq + mean_sq) / 2
		1 / rgamma(1, shape, rate)
	}
