#' Metropolis update for dyadic correlation
#' 
#' Metropolis update for dyadic correlation
#' 
#' 
#' @usage rrho_mh(Z, rho, s2 = 1,offset=0, asp=NULL)
#' @param Z n X n normal relational matrix
#' @param rho current value of rho
#' @param s2 current value of s2
#' @param offset matrix of the same dimension as Z. It is assumed that 
#' Z-offset is equal to dyadic noise, so the offset should contain any 
#' additive and multiplicative effects (such as 
#' \code{Xbeta(X,beta+ U\%*\%t(V) +  outer(a,b,"+")  }   ) 
#' @param asp use arc sine prior (TRUE) or uniform prior (FALSE) 
#' 
#' @return a new value of rho
#' @author Peter Hoff
#' @export rrho_mh
rrho_mh <-
function(Z,rho,s2=1,offset=0,asp=NULL) {
	if(is.null(asp)){ asp<-TRUE } 

	E<-Z-offset
	EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] )/sqrt(s2)
	emcp<-sum(EM[,1]*EM[,2])
	emss<-sum(EM^2)

	m<- nrow(EM)
	sr<- 2*(1-cor(EM)[1,2]^2)/sqrt(m)

	rho1<-rho+sr*qnorm( runif(1,pnorm( (-1-rho)/sr),pnorm( (1-rho)/sr)))

	# reject-and-keep: proposals outside (-.995, .995) are auto-rejected so
	# the chain stays in the truncated support. Combined with the
	# (-1,1)-normalized proposal factor Z below this is exact for pi
	# truncated to (-.995, .995).
	if(!is.finite(rho1) || abs(rho1) > .995) { return(rho) }

	# truncated-normal proposal normalizer Z(r) = Phi((1-r)/sr) - Phi((-1-r)/sr).
	# The proposal is a normal truncated to (-1,1), so the Metropolis-Hastings
	# ratio must include q(rho|rho1)/q(rho1|rho) = Z(rho)/Z(rho1); omitting it
	# targets pi(rho)*Z(rho), biasing the chain away from |rho| -> 1.
	lZ <- function(r) log(pnorm((1-r)/sr) - pnorm((-1-r)/sr))

	lhr<-(-.5*(m*log(1-rho1^2)+(emss-2*rho1*emcp)/(1-rho1^2)))-
			 (-.5*(m*log(1-rho^2 )+(emss-2*rho*emcp )/(1-rho^2 )))   +
			 asp*( (-.5*log(1-rho1^2)) - (-.5*log(1-rho^2)) )        +
			 ( lZ(rho) - lZ(rho1) )

	if(log(runif(1))<lhr) { rho<-rho1 }
	rho
}