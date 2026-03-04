#' Gibbs sampling of additive row and column effects and regression coefficient
#' with independent replicate relational data
#'
#' Simulates from the joint full conditional distribution of (a,b,beta),
#' assuming same additive row and column effects and regression coefficient
#' across replicates.
#'
#'
#' @usage rbeta_ab_rep_fc(Z.T,Sab,rho,X.T,s2=1)
#' @param Z.T n x n x T array, with the third dimension for replicates.
#' Each slice of the array is a (latent) normal relational matrix, with
#' multiplicative effects subtracted out
#' @param Sab row and column covariance
#' @param rho dyadic correlation
#' @param X.T n x n x p x T covariate array
#' @param s2 dyadic variance
#' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
#' \item{b}{additive column effects}
#' @author Peter Hoff, Yanjun He
#' @export rbeta_ab_rep_fc
rbeta_ab_rep_fc <-
	function(Z.T,Sab,rho,X.T,s2=1)  {
		####
		# decorrelation setup
		N<-dim(X.T)[4]
		p<-dim(X.T)[3]
		Se<-matrix(c(1,rho,rho,1),2,2)*s2
		iSe2<-mhalf(solve(Se))
		td<-iSe2[1,1] ; to<-iSe2[1,2]

		if(any(!is.finite(Sab))) {
			Sab <- diag(c(1, 1))
		}
		Sab <- (Sab + t(Sab))/2

		Sabs<-iSe2%*%Sab%*%iSe2
		Sabs <- (Sabs + t(Sabs))/2

		# ridge for numerical stability
		min_eig_check <- min(eigen(Sabs, symmetric=TRUE, only.values=TRUE)$values)
		if(!is.finite(min_eig_check) || min_eig_check < 1e-10) {
			Sabs <- Sabs + diag(1e-6, nrow(Sabs))
		}

		tmp<-eigen(Sabs, symmetric=TRUE)
		k<-sum(zapsmall(tmp$val)>0 )
		####

		####
		# accumulate sufficient statistics across replicates
		lb<-Qb<-Zr.T<-Zc.T<-Xr.T<-Xc.T<-0

		for (t in 1:N){
			Z<-Z.T[,,t]
			X<-array(X.T[,,,t],dim=dim(X.T)[1:3])
			p<-dim(X)[3]
			Xr<-apply(X,c(1,3),sum)
			Xc<-apply(X,c(2,3),sum)
			mX<- apply(X,3,c)
			mXt<-apply(aperm(X,c(2,1,3)),3,c)
			XX<-t(mX)%*%mX
			XXt<-t(mX)%*%mXt

			mXs<-td*mX+to*mXt
			XXs<-(to^2+td^2)*XX + 2*to*td*XXt
			Zs<-td*Z+to*t(Z)
			zr<-rowSums(Zs) ; zc<-colSums(Zs) ; zs<-sum(zc) ; n<-length(zr)

			if(p>0) {
				lb<-lb+crossprod(mXs,c(Zs))
				Qb<-Qb+XXs + (XX/nrow(mXs))/N
			}

			Xsr<-td*Xr + to*Xc
			Xsc<-td*Xc + to*Xr

			Zr.T<-Zr.T+zr
			Zc.T<-Zc.T+zc
			Xr.T<-Xr.T+Xsr
			Xc.T<-Xc.T+Xsc
		}
		####

		####
		# row and column reduction
		ab<-matrix(0,nrow(Z),2)
		if(k>0) {
			n<-nrow(Z.T[,,1])
			G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
			K<-matrix(c(0,1,1,0),2,2)

			A<-N*n*t(G)%*%G + diag(k)
			B<-N*t(G)%*%K%*%G
			iA0<-solve(A)
			C0<- -solve(A+ n*B)%*%B%*%iA0

			iA<-G%*%iA0%*%t(G)
			C<-G%*%C0%*%t(G)

			H<-iA%x%diag(n)+C%x%matrix(1,nrow=n,ncol=n)
			Hrr<-H[1:n,1:n]
			Hrc<-H[1:n,(n+1):(2*n)]
			Hcr<-H[(n+1):(2*n),1:n]
			Hcc<-H[(n+1):(2*n),(n+1):(2*n)]
			Qb<-Qb-t(Xr.T)%*%Hrr%*%Xr.T-t(Xc.T)%*%Hcr%*%Xr.T-t(Xr.T)%*%Hrc%*%Xc.T-t(Xc.T)%*%Hcc%*%Xc.T
			lb<-lb-t(Xr.T)%*%Hrr%*%Zr.T-t(Xc.T)%*%Hcr%*%Zr.T-t(Xr.T)%*%Hrc%*%Zc.T-t(Xc.T)%*%Hcc%*%Zc.T
		}
		####

		####
		# sample beta
		if(p>0)  {
			V.b<-solve(Qb)
			M.b<-V.b%*%lb
			beta<-c(rmvnorm(1,M.b,V.b))
		}
		if(p==0){ beta<-numeric(0) }
		####

		####
		# sample a, b
		Rr.T<-Zr.T-Xr.T%*%matrix(beta,ncol=1)
		Rc.T<-Zc.T-Xc.T%*%matrix(beta,ncol=1)

		m<-t(t(crossprod(rbind(c(Rr.T),c(Rc.T)),t(iA0%*%t(G)))) + rowSums(sum(Rr.T)*C0%*%t(G)) )
		hiA0<-mhalf(iA0)
		e<-matrix(rnorm(n*k),n,k)
		w<-m+ t( t(e%*%hiA0) - c(((hiA0-mhalf(iA0+n*C0))/n)%*% colSums(e) ) )
		ab.vec<- t(w%*%t(G)%*%solve(iSe2))
		####

		list(beta=beta,a=ab.vec[1,],b=ab.vec[2,] )
	}
