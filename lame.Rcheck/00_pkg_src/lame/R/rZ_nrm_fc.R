#' Simulate missing values in a normal AME model
#' 
#' Simulates missing values of a sociomatrix under a normal AME model
#' 
#' 
#' @usage rZ_nrm_fc(Z, EZ, rho,s2, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param Y square relational matrix
#' @return a square matrix, equal to Y at non-missing values
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export rZ_nrm_fc
rZ_nrm_fc<-function(Z,EZ,rho,s2,Y) {
	ZS<-simY_nrm(EZ,rho,s2)
	# diagonal only applies to square (non-bipartite) cases
	if(nrow(Y) == ncol(Y)) {
		diag(ZS)<-rnorm(nrow(Y),diag(EZ),sqrt(s2*(1+rho)))
	}
	A<-is.na(Y)
	Z[A]<-ZS[A]
	# mixed pairs (z_ij missing but its reciprocal y_ji observed): the
	# bulk draw above pairs ZS[i,j] with fresh noise, so it is marginal
	# in the reciprocal; the model full conditional is
	# N(EZ_ij + rho*(y_ji - EZ_ji), s2*(1 - rho^2)). redraw those cells
	# exactly. both-missing pairs (including the diagonal) are already
	# handled correctly by the correlated joint draw in simY_nrm.
	if(nrow(Y) == ncol(Y) && abs(rho) > 0) {
		mixed<-A & !t(A)
		if(any(mixed)) {
			idx<-which(mixed, arr.ind = TRUE)
			jdx<-idx[, c(2, 1), drop = FALSE]
			Z[idx]<-rnorm(nrow(idx), EZ[idx] + rho*(Y[jdx] - EZ[jdx]),
				sqrt(s2*(1 - rho^2)))
		}
	}
	Z
}
