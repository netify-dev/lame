#' Sociomatrix to edgelist
#'
#' Construction of an edgelist from a sociomatrix
#'
#' @usage sm2el(sm,directed=TRUE)
#' @param sm a sociomatrix with possibly valued relations
#' @param directed if TRUE, only use the upper triangular part of the matrix to enumerate edges
#' @return an edgelist
#' @author lame authors
#' @examples
#'
#' Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
#' E<-sm2el(Y)
#' el2sm(E) - Y
#'
#' @keywords internal
#' @export sm2el
sm2el<-function(sm,directed=TRUE) {
	mat<-sm
	if(!directed){
		mat[lower.tri(mat)]<-0
	}
	# linear positions of the (non-missing) nonzero relations, column-major
	nz<-which(mat!=0)
	nr<-nrow(mat)
	# recover the (row,col) coordinate of each linear position
	rows<-((nz-1L)%%nr)+1L
	cols<-((nz-1L)%/%nr)+1L
	weights<-mat[nz]
	out<-cbind(row=rows,col=cols)
	# only carry weights when the relation is genuinely valued
	if(stats::var(weights)>0){
		out<-cbind(out,w=weights)
	}
	# emit sender-ordered; ties keep column-major (i.e. receiver) order
	out[order(rows),]
}
