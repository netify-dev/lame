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
#' @param s20 fixed positive prior scale for the residual variance (default 1).
#' @return a single numeric draw of \code{s2}.
#' @keywords internal
#' @author lame authors
#' @export rs2_rep_fc
rs2_rep_fc <-
	function(E.T, rho, s20 = 1) {
		s20 <- .lame_validate_s20(s20)
		n_rep <- dim(E.T)[3]

		# whiten each dyad pair by the inverse square root of its correlation matrix
		corr <- matrix(c(1, rho, rho, 1), 2, 2)
		whiten <- mhalf(solve(corr))

		ssq <- 0 # pooled whitened sum of squares
		count <- 0L # pooled count of whitened scalars

		for (r in seq_len(n_rep)) {
			slice <- E.T[, , r]
			above <- upper.tri(slice)
			# stack the two ordered residuals of each dyad, then whiten
			dyad_pairs <- cbind(slice[above], t(slice)[above])
			whitened <- dyad_pairs %*% whiten
			ssq <- ssq + sum(whitened * whitened)
			count <- count + length(whitened)
		}

		# inverse-gamma draw with one prior pseudo-observation at scale s20
		shape <- (count + 1) / 2
		rate <- (ssq + s20) / 2
		1 / rgamma(1, shape, rate)
	}
