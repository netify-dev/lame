// Per-actor blocks of the joint information of the effects and factors of the empirical-bayes ALS
// estimator (R/ame_als_joint_block.R): on the longer side of a bipartite or directed fit every actor's
// effect and its R factor entries couple only to each other and to the other side, so the information
// is block diagonal there with blocks of size R + 1, bordered by the rest.

#include <RcppArmadillo.h>

// Cholesky factor of every block D_i (k x k x n, the upper triangle read) and, for the border rows
// `Bt` (n k x ns, the k rows of actor i together), Gt_i = D_i^-1 Bt_i. `bad` is the first actor
// (1-based) whose block is not positive definite, 0 when all are; `ld` the sum of log |D_i|
// [[Rcpp::export]]
Rcpp::List ae_jb_blocks_cpp(const arma::cube& D, const arma::mat& Bt) {
	const arma::uword k = D.n_rows, n = D.n_slices;
	arma::cube Di(k, k, n);
	arma::mat Gt(Bt.n_rows, Bt.n_cols);
	double ld = 0.0;
	int bad = 0;
	arma::mat L(k, k), Li(k, k);
	for (arma::uword i = 0; i < n; ++i) {
		arma::mat Ds = arma::symmatu(D.slice(i));
		bool ok = arma::chol(L, Ds, "upper");
		if (ok) for (arma::uword a = 0; a < k; ++a) if (!(L(a, a) > 0.0)) ok = false;
		if (!ok) {
			bad = static_cast<int>(i) + 1;
			break;
		}
		for (arma::uword a = 0; a < k; ++a) ld += 2.0 * std::log(L(a, a));
		Li = arma::inv(arma::trimatu(L));
		Di.slice(i) = arma::symmatu(Li * Li.t());
		Gt.rows(i * k, i * k + k - 1) = Di.slice(i) * Bt.rows(i * k, i * k + k - 1);
	}
	return Rcpp::List::create(
		Rcpp::Named("bad") = bad,
		Rcpp::Named("ld") = ld,
		Rcpp::Named("Dinv") = Di,
		Rcpp::Named("Gt") = Gt
	);
}

// natural-order Cholesky pivots of the factor entries of the long side, which follow its effects and
// the border variables `Bt` couples them to (n k x m, the columns of those variables only): u_ri for
// r = 1..R and then i = 1..n. `Si` is the inverse of the Schur complement of those variables given the
// effects of the long side. Adding u_ri to the block of actor i lowers that Schur complement by f f' / s,
// with s = h - d' D^-1 d and f = e - B' D^-1 d (d, e its couplings to the block and to the border); its
// pivot is s - f' Si f, and Si takes the rank-one update of the Sherman-Morrison formula. the
// recurrence stops at the first pivot that is not positive; the later entries are NA
// [[Rcpp::export]]
Rcpp::NumericVector ae_jb_pivots_cpp(const arma::cube& D, const arma::mat& Bt, arma::mat Si, int R) {
	const arma::uword k = D.n_rows, n = D.n_slices;
	Rcpp::NumericVector out(n * static_cast<arma::uword>(R), NA_REAL);
	for (arma::uword r = 1; r <= static_cast<arma::uword>(R); ++r) {
		for (arma::uword i = 0; i < n; ++i) {
			const arma::mat& Dsl = D.slice(i);
			arma::mat Dc = arma::symmatu(Dsl.submat(0, 0, r - 1, r - 1));
			arma::vec dv = Dsl.submat(0, r, r - 1, r);
			arma::vec sol;
			if (!arma::solve(sol, Dc, dv, arma::solve_opts::no_approx)) return out;
			const double s = Dsl(r, r) - arma::dot(dv, sol);
			arma::vec f = Bt.row(i * k + r).t() - Bt.rows(i * k, i * k + r - 1).t() * sol;
			arma::vec g = Si * f;
			const double piv = s - arma::dot(f, g);
			out[(r - 1) * n + i] = piv;
			if (!(piv > 0.0)) return out;
			Si += (g * g.t()) / piv;
		}
	}
	return out;
}
