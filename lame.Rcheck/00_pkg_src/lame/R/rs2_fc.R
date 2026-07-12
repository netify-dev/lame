#' Gibbs update for dyadic variance
#' 
#' Gibbs update for dyadic variance
#' 
#' 
#' @usage rs2_fc(Z, rho,offset=0,nu0=NULL,s20=NULL)
#' @param Z n X n normal relational matrix
#' @param rho current value of rho 
#' @param offset matrix of the same dimension as Z. It is assumed that Z-offset
#' is equal to dyadic noise, so the offset should contain any additive and 
#' multiplicative effects (such as \code{Xbeta(X,beta+ U\%*\%t(V) +  
#' outer(a,b,"+")  } )
#' @param nu0 prior degrees of freedom
#' @param s20 prior estimate of s2; when \code{NULL} (default) the
#'   empirical mean-square residual is used, so the prior carries
#'   \code{nu0} pseudo-observations at the scale of the data rather than
#'   an absolute guess of 1 (which corrupts s2 -- and through it rho --
#'   whenever Y is not on unit scale)
#' @return a new value of s2
#' @author Peter Hoff
#' @export rs2_fc
rs2_fc <-
function(Z,rho,offset=0,nu0=NULL,s20=NULL) {
		E<-Z-offset
		if(is.null(nu0)){ nu0<-1 }

		H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
		EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H
		ED<-diag(E)/sqrt(1+rho)
		if(is.null(s20)){
			ms <- (sum(EM^2,na.rm=TRUE)+sum(ED^2,na.rm=TRUE)) /
				max(sum(is.finite(EM))+sum(is.finite(ED)), 1)
			s20 <- if (is.finite(ms) && ms > 0) ms else 1
		}
		1/rgamma(1, (length(EM)+length(ED)+nu0)/2 , (sum(EM^2)+sum(ED^2)+nu0*s20)/2 )
}