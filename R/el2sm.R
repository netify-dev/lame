#' Edgelist to sociomatrix 
#' 
#' Construction of a sociomatrix from an edgelist
#' 
#' @usage el2sm(el,directed=TRUE,nadiag=all(el[,1]!=el[,2]))
#' @param el a matrix in which each row contains the indices of an edge and possibly the weight for the edge
#' @param directed if FALSE, then a relation is placed in both entry ij and ji of the sociomatrix, for each edge ij (or ji)
#' @param nadiag put NAs on the diagonal 
#' @return a sociomatrix 
#' @author Peter Hoff
#' @examples
#' 
#' Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
#' E<-sm2el(Y) 
#' el2sm(E) - Y 
#' 
#' @export el2sm
el2sm<-function(el,directed=TRUE,nadiag=all(el[,1]!=el[,2])) {
	el <- as.matrix(el)            # accept a data.frame edgelist too

	# edge weights: column 3 if present, else 1. Coerce to numeric so a
	# character edgelist (e.g. from igraph::as_edgelist) does not silently
	# produce a character sociomatrix.
	w <- rep(1, nrow(el))
	if(ncol(el) > 2){
		w <- suppressWarnings(as.numeric(el[,3]))
		if(any(is.na(w))){
			stop("el2sm(): the weight column (column 3) is not numeric.",
			     call. = FALSE)
		}
	}

	# node set: from the first two (id) columns ONLY -- never the weight column
	ids <- el[,1:2]
	numeric_ids <- is.numeric(ids) && all(round(ids) == ids)
	nodes <- if(numeric_ids) 1:max(ids) else sort(unique(c(ids)))

	elij <- cbind( match(el[,1],nodes), match(el[,2],nodes) )

	n  <- length(nodes)
	sm <- matrix(0,n,n)            # construct sociomatrix
	sm[elij[,1:2]] <- w            # fill in
	if(nadiag) { diag(sm) <- NA  } # set diagonal to NA
	if(!directed){ sm<-sm+t(sm) }
	dimnames(sm)[[1]]<-dimnames(sm)[[2]]<-nodes
	sm
}
