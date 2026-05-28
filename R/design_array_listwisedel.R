#' Computes the design socioarray of covariate values
#'
#' Computes the design socioarray of covariate values for an AME fit
#'
#'
#' @usage design_array_listwisedel(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n)
#' @param Xrow an n x pr matrix of row covariates
#' @param Xcol an n x pc matrix of column covariates
#' @param Xdyad an n x n x pd array of dyadic covariates
#' @param intercept logical
#' @param n number of rows/columns
#' @return an n x n x (pr+pc+pd+intercept) 3-way array
#' @author Peter Hoff
#' @export design_array_listwisedel
design_array_listwisedel<-function(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n) {

	####
	# covariate array
	pr<-length(Xrow)/n
	pc<-length(Xcol)/n
	pd<-length(Xdyad)/n^2
	X<-array(dim=c(n,n,pr+pc+pd))
	dnX<-NULL
	####

	####
	# row covariates
	if(pr>0) {
		Xrow<-as.matrix(Xrow)
		Xrowa<-array(dim=c(n,n,pr))
		for( j in 1:pr ){ Xrowa[,,j]<-matrix( Xrow[,j], n,n) }
		X[,,1:pr]<- Xrowa
		row_nms <- colnames(Xrow)
		if (is.null(row_nms)) row_nms <- paste0("Xrow", seq_len(pr))
		dnX <- c(dnX, .lame_apply_suffix(row_nms, "row"))
	}
	####

	####
	# column covariates
	if(pc>0) {
		Xcol<-as.matrix(Xcol)
		Xcola<-array(dim=c(n,n,pc))
		for( j in 1:pc ){ Xcola[,,j]<-t(matrix( Xcol[,j], n,n)) }
		X[,,pr+1:pc]<- Xcola
		col_nms <- colnames(Xcol)
		if (is.null(col_nms)) col_nms <- paste0("Xcol", seq_len(pc))
		dnX <- c(dnX, .lame_apply_suffix(col_nms, "col"))
	}
	####

	####
	# dyadic covariates
	if(pd>0)                                                {
		dyad_names <- dimnames(Xdyad)[[3]]
		if(pd==1){ Xdyad<-array(Xdyad,dim=c(n,n,pd)) }
		X[,,pr+pc+1:pd]<-Xdyad
		if(is.null(dyad_names)) dyad_names <- paste0("Xdyad", 1:pd)
		dnX <- c(dnX, .lame_apply_suffix(dyad_names, "dyad"))
	}
	####

	####
	# add intercept
	if(!any(apply(X,3,function(x){var(c(x),na.rm=TRUE)})==0) & intercept ) {
		X1<-array(dim=c(0,0,1)+dim(X))
		X1[,,1]<-1 ; X1[,,-1]<-X
		X<-X1
		dnX<-c("intercept",dnX)
	}
	####

	####
	# set variable names
	if(dim(X)[[3]]>1) { dimnames(X)[[3]]<- dnX }
	if(dim(X)[[3]]==1){ dimnames(X)[[3]]<- list(dnX) }
	####

	X
}
