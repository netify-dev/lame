#' Computes the design socioarray of covariate values
#'
#' Computes the design socioarray of covariate values for an AME fit
#'
#'
#' @usage design_array(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n,warn=TRUE)
#' @param Xrow an n x pr matrix of row covariates
#' @param Xcol an n x pc matrix of column covariates
#' @param Xdyad an n x n x pd array of dyadic covariates
#' @param intercept logical
#' @param n number of rows/columns
#' @param warn logical; warn when missing covariate values are zero-filled
#'   (default \code{TRUE}). Set \code{FALSE} when the caller has already
#'   propagated the missingness into the response.
#' @return an n x n x (pr+pc+pd+intercept) 3-way array
#' @author Peter Hoff
#' @export design_array
design_array<-function(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n,warn=TRUE) {

	####
	# coerce data.frames up front -- length(df) returns ncol(df), not nrow*ncol,
	# which would give a fractional pr/pc and a subscript-out-of-bounds crash
	# in the array allocation below. (the later as.matrix() calls happen after
	# pr/pc are already wrong.)
	if (!is.null(Xrow) && is.data.frame(Xrow)) Xrow <- as.matrix(Xrow)
	if (!is.null(Xcol) && is.data.frame(Xcol)) Xcol <- as.matrix(Xcol)

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
		row_names <- colnames(Xrow)
		if(is.null(row_names)) {
			row_names <- paste0("Xrow", 1:pr)
		}
		dnX<-c(dnX, .lame_apply_suffix(row_names, "row"))
	}
	####

	####
	# column covariates
	if(pc>0) {
		Xcol<-as.matrix(Xcol)
		Xcola<-array(dim=c(n,n,pc))
		for( j in 1:pc ){ Xcola[,,j]<-t(matrix( Xcol[,j], n,n)) }
		X[,,pr+1:pc]<- Xcola
		col_names <- colnames(Xcol)
		if(is.null(col_names)) {
			col_names <- paste0("Xcol", 1:pc)
		}
		dnX<-c(dnX, .lame_apply_suffix(col_names, "col"))
	}
	####

	####
	# dyadic covariates
	if(pd>0)                                                {
		# capture the covariate name before the single-slice rebuild
		dyad_names <- dimnames(Xdyad)[[3]]
		if(pd==1){ Xdyad<-array(Xdyad,dim=c(n,n,pd)) }
		X[,,pr+pc+1:pd]<-Xdyad
		if(is.null(dyad_names)) {
			dyad_names <- paste0("Xdyad", 1:pd)
		}
		dnX<-c(dnX, .lame_apply_suffix(dyad_names, "dyad"))
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

	####
	# replace missing values with zeros
	if( warn && sum(is.na(X)) > sum( is.na(apply(X,3,diag)) ) ) {
		if(getOption("lame.warn.na", TRUE)) {
			warning("Replacing NAs in design matrix with zeros", call. = FALSE)
		}
	}
	X[is.na(X)]<-0
	####

	X
}

#' Propagate covariate missingness into the response
#'
#' For an AME fit, a dyad whose covariate value is missing cannot contribute a
#' covariate-coefficient observation. Rather than silently imputing the missing
#' covariate as 0 (which biases the coefficient when the covariate is not
#' mean-centred), the dyad is treated as an unobserved tie and handled by data
#' augmentation. This helper sets the affected \code{Y} cells to \code{NA};
#' \code{\link{lame}} already does this internally, and this brings
#' \code{\link{ame}} into line.
#'
#' @param Y an n x n response matrix.
#' @param Xrow,Xcol,Xdyad row, column and dyadic covariates (or \code{NULL}).
#' @return \code{Y} with covariate-missing cells set to \code{NA}.
#' @keywords internal
.ame_propagate_cov_na <- function(Y, Xrow = NULL, Xcol = NULL, Xdyad = NULL) {
	n <- nrow(Y)
	if (!is.null(Xrow)) {
		Xr  <- as.matrix(Xrow)
		bad <- which(!stats::complete.cases(Xr))
		if (length(bad)) Y[bad, ] <- NA
	}
	if (!is.null(Xcol)) {
		Xc  <- as.matrix(Xcol)
		bad <- which(!stats::complete.cases(Xc))
		if (length(bad)) Y[, bad] <- NA
	}
	if (!is.null(Xdyad)) {
		Xd <- if (length(dim(Xdyad)) == 3L) Xdyad else
			array(Xdyad, dim = c(n, n, 1L))
		miss <- apply(Xd, c(1, 2), function(z) any(is.na(z)))
		Y[miss] <- NA
	}
	Y
}
