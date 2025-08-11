#' Griddy Gibbs update for dyadic correlation
#'
#' Simulation of dyadic correlation from its approximate full conditional 
#' distribution using griddy Gibbs sampling
#'
#' @param Z n X n normal relational matrix
#' @param Sab covariance of additive effects
#' @param s2 residual variance
#' @param offset matrix of the same dimension as Z. It is assumed that 
#' Z-offset follows an SRM distribution, so the offset should contain any 
#' regression terms and multiplicative effects (such as 
#' \code{Xbeta(X,beta+ U\%*\%t(V) }   )
#' @param ngp the number of points for an unevenly-spaced 
#' grid on which to approximate the full conditional distribution 
#' @param asp use arc sine prior (TRUE) or uniform prior (FALSE) 
#' 
#' @return a value of rho
#' @author Peter Hoff
#' @export rrho_fc
rrho_fc<-function(Z,Sab,s2=1,offset=0,ngp=100,asp=NULL)
{ 
  if(is.null(asp)){ asp<-TRUE } 

  E<-Z-offset

  ## first obtain rough estimate of rho and its sd
  m<-mean(E,na.rm=TRUE)
  a<-apply(E,1,mean,na.rm=TRUE) - m
  b<-apply(E,2,mean,na.rm=TRUE) - m
  R<-E-(m+outer(a,b,"+"))
  RM<-cbind(R[upper.tri(R)],t(R)[upper.tri(R)] )/sqrt(s2)
  emcp<-sum(RM[,1]*RM[,2])
  emss<-sum(RM^2)
  nd<- nrow(RM)
  rho0<-cor(RM[,1],RM[,2])
  srho0<- 2*(1-rho0^2)/sqrt(nd)

  ## grid concentrated near likely regions of high probability 
  rhos<-qnorm( (1:ngp)/(ngp+1), rho0,2*srho0 )
  rhos<-rhos[ -1<rhos & rhos<1 ]

  ## griddy Gibbs
  ll<-llsrmRho(E,Sab,rhos,s2)
  prho<-exp( ll-max(ll) - asp*.5*log(1-rhos^2) )
  Frho<-c(0,cumsum(prho)/sum(prho),1)
  rhos<-c(0,rhos,1)
  f<-runif(1)
  k<-max(which(Frho<f))
  rho<-rhos[k] + (rhos[k+1]-rhos[k])*(f-Frho[k])/(Frho[k+1]-Frho[k])
  min(abs(rho),.995)*sign(rho)
}