// Natural-order Cholesky pivots of a diagonal matrix plus a low-rank term, for the structured actor
// block of the empirical-bayes ALS estimator (R/ame_als_actor_block.R): the row and column blocks of the
// penalised information of the actor effects are diagonal, so the pivots of a block that follows them
// come from a rank-r recurrence instead of a dense factorisation.

#include <RcppArmadillo.h>

// pivots d_k (the squared diagonal of the Cholesky factor) of D + Z' P Z, with D = diag(`D`), the
// columns of `Zt` (r x n) the rows of Z' and P symmetric r x r:
//   d_k = D_k + z_k' P_k z_k,  P_{k+1} = P_k - (P_k z_k)(P_k z_k)' / d_k.
// the recurrence stops at the first pivot that is not positive (where a Cholesky factorisation fails);
// the later entries are NA
// [[Rcpp::export]]
Rcpp::NumericVector ae_dlr_pivots_cpp(const arma::vec& D, const arma::mat& Zt, arma::mat P) {
	const arma::uword n = D.n_elem;
	Rcpp::NumericVector d(n, NA_REAL);
	arma::vec pz(P.n_rows);
	for (arma::uword k = 0; k < n; ++k) {
		pz = P * Zt.col(k);
		const double dk = D[k] + arma::dot(Zt.col(k), pz);
		d[k] = dk;
		if (!(dk > 0.0)) break;
		P -= (pz * pz.t()) / dk;
	}
	return d;
}
