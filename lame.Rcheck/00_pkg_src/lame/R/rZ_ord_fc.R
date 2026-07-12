#' Simulate Z given the partial ranks
#' 
#' Simulates a random latent matrix Z given its expectation, dyadic correlation
#' and partial rank information provided by W
#' 
#' 
#' @usage rZ_ord_fc(Z, EZ, rho, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param Y matrix of ordinal data
#' @return a square matrix, the new value of Z
#' @author Peter Hoff
#' @export rZ_ord_fc
rZ_ord_fc<-
	function (Z, EZ, rho,Y) {
		uY<-sort(unique(c(Y)))
		W<-matrix(match(Y,uY),nrow(Y),nrow(Y))
		W[is.na(W)] <- -1
		
		sz <- sqrt(1 - rho^2)
		ut <- upper.tri(Z)
		lt <- lower.tri(Z)
		for (w in sample(c(-1, 1:length(uY)))) {
			lb <- suppressWarnings(max(Z[!is.na(W) & W == w - 1],
																 na.rm = TRUE))
			ub <- suppressWarnings(min(Z[!is.na(W) & W == w + 1],
																 na.rm = TRUE))
			up <- ut & W == w
			ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
			zup <- sz * rtnorm_interval_logp(ez / sz, lb / sz, ub / sz)
			# tail-stable two-sided truncated draw (rZ_bipartite.R sibling);
			# a naive qnorm(runif(pnorm(a), pnorm(b))) collapses to NaN/Inf
			# once the (lb, ub) interval sits > ~8 sd from ez. backstop
			# keeps the current value for any non-finite draw.
			zerr <- which(!is.finite(zup))
			if (length(zerr)) zup[zerr] <- (Z[up])[zerr]
			Z[up] <- zup
			up <- lt & W == w
			ez <- EZ[up] + rho * (t(Z)[up] - t(EZ)[up])
			zup <- sz * rtnorm_interval_logp(ez / sz, lb / sz, ub / sz)
			# tail-stable two-sided truncated draw (rZ_bipartite.R sibling);
			# a naive qnorm(runif(pnorm(a), pnorm(b))) collapses to NaN/Inf
			# once the (lb, ub) interval sits > ~8 sd from ez. backstop
			# keeps the current value for any non-finite draw.
			zerr <- which(!is.finite(zup))
			if (length(zerr)) zup[zerr] <- (Z[up])[zerr]
			Z[up] <- zup
		}
		diag(Z)<-rnorm(nrow(Z),diag(EZ),sqrt(1+rho))
		Z
	}
