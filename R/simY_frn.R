#' Simulate a relational matrix under a fixed rank nomination scheme
#'
#' Simulate a sociomatrix of fixed rank nominations from a social relations
#' model. A latent matrix \code{Z} is drawn from the SRM, and each sender's
#' outgoing ties are converted into ranked nominations: at most \code{odmax}
#' partners are nominated, only positive latent affinities qualify, and the
#' retained nominations are numbered from weakest (1) to strongest.
#'
#' @usage simY_frn(EZ, rho, odmax, YO = NULL, bipartite = FALSE)
#' @param EZ a square matrix giving the expected value of the latent Z matrix
#' @param rho dyadic correlation
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each node
#' @param YO a square matrix identifying where missing values should be
#' maintained
#' @param bipartite whether row and column actors are disjoint; square
#'   diagonal cells remain eligible for nomination.
#' @return a square matrix, where higher values represent stronger
#' relationships
#' @keywords internal
#' @author lame authors
#' @export simY_frn
simY_frn <-
	function(EZ, rho, odmax, YO = NULL, bipartite = FALSE) {
		n_send <- nrow(EZ)
		n_recv <- ncol(EZ)

		# per-sender nomination budget
		cap <- if (length(odmax) == 1) rep(odmax, n_send) else odmax

		# latent affinities under the social relations model
		z_lat <- simZ(EZ, rho, bipartite = bipartite)

		# self-ties are never eligible (only meaningful when square)
		square <- !isTRUE(bipartite) && (n_send == n_recv)
		if (square) diag(z_lat) <- -Inf

		# structurally missing cells are removed from contention
		miss <- if (!is.null(YO)) is.na(YO) else NULL
		if (!is.null(miss)) z_lat[miss] <- -Inf

		if (.frn_rankable(z_lat)) {
			nominations <- .frn_top_ranks(z_lat, cap)
		} else {
			# missing or tied affinities: row by row, ties sharing their average rank
			nominations <- matrix(0, n_send, n_recv)
			for (i in seq_len(n_send)) {
				affinity <- z_lat[i, ]

				# ascending ranks: the largest affinity earns the top rank
				asc_rank <- rank(affinity)
				eligible <- (asc_rank > n_recv - cap[i]) & (affinity > 0)

				if (any(eligible)) {
					nominations[i, eligible] <- rank(affinity[eligible])
				}
			}
		}

		# restore missingness markers in the returned matrix
		if (square) diag(nominations) <- NA
		if (!is.null(miss)) nominations[miss] <- NA

		nominations
	}

# whether .frn_top_ranks() reproduces the row-by-row ranking: no missing
# affinities and no ties among the positive ones
.frn_rankable <- function(z) {
	!anyNA(z) && !anyDuplicated(z[z > 0])
}

# ranked nominations of every row at once: the positive affinities among each
# sender's top cap[i] are numbered 1..k by increasing strength, all else is 0
.frn_top_ranks <- function(z, cap) {
	n_send <- nrow(z)
	n_recv <- ncol(z)
	rows <- rep.int(seq_len(n_send), n_recv)
	# the place of each cell in its row from the strongest down
	ord <- order(rows, -z)
	place <- integer(length(z))
	place[ord] <- seq_along(ord) - (rows[ord] - 1L) * n_recv
	eligible <- (place <= cap[rows]) & (z > 0)
	n_kept <- tabulate(rows[eligible], n_send)
	ranks <- matrix(0, n_send, n_recv)
	ranks[eligible] <- n_kept[rows[eligible]] - place[eligible] + 1
	ranks
}
