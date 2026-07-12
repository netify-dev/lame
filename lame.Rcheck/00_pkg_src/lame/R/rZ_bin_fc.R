#' Simulate Z based on a probit model
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and a binary relational matrix Y
#' 
#' 
#' @usage rZ_bin_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y square binary relational matrix
#' @return a square matrix , the new value of Z
#' @author Peter Hoff
#' @export rZ_bin_fc
rZ_bin_fc <-
function(Z,EZ,rho,Y) {

	####
	# truncated normal gibbs update by triangle. one-sided truncated draws
	# use the log.p-scale inverse cdf (matching draw_tn_cell in
	# src/rZ_bin_fused.cpp): z = ez + sz * qnorm(log(u) + log-tail-mass,
	# log.p=TRUE) with pnorm/qnorm on the log scale in the tail that holds
	# the truncation interval. This is exact for arbitrary |ez| (the naive
	# qnorm(runif(pnorm(lo), pnorm(hi))) collapses to a 0-width interval once
	# |ez| exceeds ~8.2, which freezes constraint-violating cells at their
	# initial, possibly sign-wrong, value for the whole chain).
	sz<-sqrt(1-rho^2)
	ut<-upper.tri(EZ)
	lt<-lower.tri(EZ)

	Y[is.na(Y)]<- -1
	for(y in c((-1):1)) {
		lb<-c(-Inf,-Inf,0)[y+2] ; ub<-c(Inf,0,Inf)[y+2]

		for(tri in 1:2) {
			if(tri==1){ up<-ut & Y==y }
			if(tri==2){ up<-lt & Y==y }
			nup<-sum(up)
			if(nup==0) next

			ez<- EZ[up] + rho*( t(Z)[up]  - t(EZ)[up] )
			u <- runif(nup)
			if(lb == -Inf && ub == Inf) {
				# unconstrained (missing cells)
				zup <- ez + sz*qnorm(u)
			} else if(ub == Inf) {
				# [lb, Inf): upper-tail inverse cdf on the log scale
				logS <- pnorm((lb-ez)/sz, lower.tail=FALSE, log.p=TRUE)
				zup <- ez + sz*qnorm(log(u) + logS, lower.tail=FALSE, log.p=TRUE)
			} else {
				# (-Inf, ub]: lower-tail inverse cdf on the log scale
				logF <- pnorm((ub-ez)/sz, lower.tail=TRUE, log.p=TRUE)
				zup <- ez + sz*qnorm(log(u) + logF, lower.tail=TRUE, log.p=TRUE)
			}
			# backstop: a non-finite draw needs u at an exact endpoint,
			# which R's rng never returns; keep the current value if so
			zerr<-which(!is.finite(zup))
			if(length(zerr)>0){ zup[zerr]<-(Z[up])[zerr] }
			Z[up]<-zup
		}
	}
	####

	####
	# joint proposal for consistent dyads
	c<-(sqrt(1+rho) + sqrt(1-rho))/2
	d<-(sqrt(1+rho) - sqrt(1-rho))/2
	E<-matrix(rnorm(nrow(Y)^2),nrow(Y),nrow(Y))
	ZP<-EZ + c*E + d*t(E)
	A<-( (Y== -1) | ( sign(ZP) == sign(Y-.5)) ) ; diag(A)<-TRUE
	A<-A & t(A)
	Z[A]<-ZP[A]
	####

	diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))

	Z
}