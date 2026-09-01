#' Edgelist to sociomatrix
#'
#' Construct a sociomatrix from an edgelist. Each row of \code{el} gives the
#' two endpoints of an edge and, optionally, a weight in a third column
#' (unweighted edges default to 1).
#'
#' @usage el2sm(el,directed=TRUE,nadiag=all(el[,1]!=el[,2]))
#' @param el a matrix in which each row contains the indices of an edge and possibly the weight for the edge
#' @param directed if FALSE, then a relation is placed in both entry ij and ji of the sociomatrix, for each edge ij (or ji)
#' @param nadiag put NAs on the diagonal
#' @return a sociomatrix
#' @author lame authors
#' @examples
#'
#' Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
#' E<-sm2el(Y)
#' el2sm(E) - Y
#'
#' @export el2sm
el2sm <- function(el, directed = TRUE, nadiag = all(el[, 1] != el[, 2])) {
	el <- as.matrix(el)

	# edge weights come from a third column when present; otherwise every
	# edge counts once. coerce explicitly so a character edgelist does not
	# quietly yield a character-valued sociomatrix.
	weights <- if (ncol(el) >= 3L) {
		w <- suppressWarnings(as.numeric(el[, 3L]))
		if (anyNA(w)) {
			cli::cli_abort("el2sm(): the weight column (column 3) is not numeric.")
		}
		w
	} else {
		rep.int(1, nrow(el))
	}

	# node universe from the endpoint columns only. integer-labelled nodes
	# span 1..max so isolates keep their slot; otherwise use the sorted set
	# of distinct labels.
	from <- el[, 1L]
	to   <- el[, 2L]
	endpoints <- c(from, to)
	integer_labels <- is.numeric(endpoints) && all(endpoints == round(endpoints))
	nodes <- if (integer_labels) seq_len(max(endpoints)) else sort(unique(endpoints))
	n <- length(nodes)

	# map each endpoint to its row/column position and deposit the weights
	ri <- match(from, nodes)
	ci <- match(to,   nodes)
	sm <- matrix(0, n, n, dimnames = list(nodes, nodes))
	sm[cbind(ri, ci)] <- weights

	if (!directed) sm <- sm + t(sm)
	if (nadiag) diag(sm) <- NA
	sm
}
