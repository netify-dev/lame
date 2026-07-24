#' Simulate a relational matrix under a fixed rank nomination scheme
#'
#' Simulate a sociomatrix of fixed rank nominations from a social relations
#' model. A latent matrix \code{Z} is drawn from the SRM, and each sender's
#' outgoing ties are converted into ranked nominations: at most \code{odmax}
#' partners are nominated, only positive latent affinities qualify, and the
#' retained nominations are numbered from weakest (1) to strongest.
#'
#' @usage simY_frn(EZ, rho, odmax, YO)
#' @param EZ a square matrix giving the expected value of the latent Z matrix
#' @param rho dyadic correlation
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each node
#' @param YO a square matrix identifying where missing values should be
#' maintained
#' @return a square matrix, where higher values represent stronger
#' relationships
#' @keywords internal
#' @author lame authors
#' @export simY_frn
simY_frn <-
function(EZ,rho,odmax,YO=NULL) {
	n_send <- nrow(EZ)
	n_recv <- ncol(EZ)

	# per-sender nomination budget
	cap <- if(length(odmax)==1) rep(odmax,n_send) else odmax

	# latent affinities under the social relations model
	z_lat <- simZ(EZ,rho)

	# self-ties are never eligible (only meaningful when square)
	square <- (n_send==n_recv)
	if(square) diag(z_lat) <- -Inf

	# structurally missing cells are removed from contention
	miss <- if(!is.null(YO)) is.na(YO) else NULL
	if(!is.null(miss)) z_lat[miss] <- -Inf

	nominations <- matrix(0,n_send,n_recv)

	for(i in seq_len(n_send)) {
		affinity <- z_lat[i,]

		# ascending ranks: the largest affinity earns the top rank
		asc_rank <- rank(affinity)

		# eligible = among this sender's top-`cap[i]` targets AND strictly positive
		eligible <- (asc_rank > n_recv - cap[i]) & (affinity > 0)

		if(any(eligible)) {
			# renumber the retained targets 1..k by increasing strength
			nominations[i,eligible] <- rank(affinity[eligible])
		}
	}

	# restore missingness markers in the returned matrix
	if(square) diag(nominations) <- NA
	if(!is.null(miss)) nominations[miss] <- NA

	nominations
}
