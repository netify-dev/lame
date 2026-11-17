// [[Rcpp::depends(RcppArmadillo)]]
// gibbs sweeps of the snap kind of dynamic_uv for bipartite and symmetric panels. each position u_it is
// drawn with its reset indicator delta_it from their joint conditional given the residual array, the other
// positions and the next period's indicator, as rUV_dynamic_snap_fc_cpp does for directed panels: the
// indicator by the ratio of the two Gaussian marginals (drift N(rho u_i,t-1, sigma^2 I), reset
// N(0, kappa^2 I)), then the position from the chosen conditional. a bipartite row regresses on V_t G' and
// a column on U_t G; a symmetric row regresses on the off-diagonal cells of its row with partner
// weights l * u_j (the self cell is quadratic in u_i and left out). state arrays are copied, so a refused
// sweep leaves the caller's positions as they were
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// shared with src/dynamic_latent_factors.cpp
double snap_forward_log_ratio(const arma::vec& x_prev, const arma::vec& x_next, double rho_uv,
	double sigma2, double kappa2, int R);
int snap_draw_row(arma::rowvec& update_row, const arma::mat& VtV, const arma::vec& ei_data,
	const arma::vec& back_mean, double s2_sigma2_inv, double s2_kappa_inv, double fwd_diag,
	const arma::vec& fwd_ei, double logit_pi, double s2, int R, bool has_back, double c_stat,
	double fwd_log_ratio);

namespace {

struct SnapPar {
	double s2, sigma2, kappa2, s2_sigma2_inv, s2_rho_sigma2, s2_rho2_sigma2, s2_kappa_inv, s2_stat_init, logit_pi, rho;
};

SnapPar snap_par(double rho_uv, double sigma_uv, double s2, double kappa, double pi_snap) {
	if (!std::isfinite(rho_uv) || std::abs(rho_uv) >= 1.0 || !std::isfinite(sigma_uv) || sigma_uv <= 0.0 ||
		!std::isfinite(s2) || s2 <= 0.0 || !std::isfinite(kappa) || kappa <= 0.0)
		Rcpp::stop("Invalid snap position state or scale.");
	if (!std::isfinite(pi_snap) || pi_snap < 0.0 || pi_snap > 1.0)
		Rcpp::stop("Snap probability must lie between zero and one.");
	const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
	SnapPar p;
	p.s2 = s2;
	p.sigma2 = sigma_uv * sigma_uv;
	p.kappa2 = kappa * kappa;
	p.s2_sigma2_inv = s2 * sigma2_inv;
	p.s2_rho_sigma2 = s2 * rho_uv * sigma2_inv;
	p.s2_rho2_sigma2 = s2 * rho_uv * rho_uv * sigma2_inv;
	p.s2_kappa_inv = s2 / (kappa * kappa);
	p.s2_stat_init = s2 * (1.0 - rho_uv * rho_uv) * sigma2_inv;
	p.logit_pi = std::log(pi_snap) - std::log1p(-pi_snap);
	p.rho = rho_uv;
	return p;
}

// draw row i of period t of the paths X (actors x R x T) with its reset, given the data precision P and the
// data numerator ei; `dcur` holds the current resets (read at t + 1), `dnew` takes the draw at t
int snap_row(arma::cube& X, int i, int t, const arma::mat& P, const arma::vec& ei, const SnapPar& p,
	const arma::mat& dcur, arma::mat& dnew, bool bounded, double max_abs, int& rejections) {
	const int R = X.n_cols;
	const int T = X.n_slices;
	double fwd_diag = 0.0;
	arma::vec fwd_ei(R, fill::zeros);
	const bool fwd_on = (t < T - 1) && (dcur(i, t + 1) < 0.5);
	if (fwd_on) {
		fwd_diag = p.s2_rho2_sigma2;
		fwd_ei = p.s2_rho_sigma2 * X.slice(t + 1).row(i).t();
	}
	arma::vec back_mean(R, fill::zeros);
	if (t > 0) back_mean = p.s2_rho_sigma2 * X.slice(t - 1).row(i).t();
	double fwd_lr = 0.0;
	if (t > 0 && fwd_on) {
		fwd_lr = snap_forward_log_ratio(X.slice(t - 1).row(i).t(), X.slice(t + 1).row(i).t(), p.rho,
			p.sigma2, p.kappa2, R);
	}
	const arma::rowvec old = X.slice(t).row(i);
	arma::rowvec new_row = old;
	int d = snap_draw_row(new_row, P, ei, back_mean, p.s2_sigma2_inv, p.s2_kappa_inv, fwd_diag, fwd_ei,
		p.logit_pi, p.s2, R, t > 0, p.s2_stat_init, fwd_lr);
	if (!new_row.is_finite()) Rcpp::stop("Non-finite snap position draw.");
	if (bounded && arma::abs(new_row).max() > max_abs) {
		new_row = old;
		d = static_cast<int>(dcur(i, t));
		++rejections;
	}
	X.slice(t).row(i) = new_row;
	dnew(i, t) = d;
	return d;
}

}

//' Bipartite snap-shift dynamic UV Gibbs update
//'
//' @param U_state 3D array (nA x RA x T)
//' @param V_state 3D array (nB x RB x T)
//' @param E 3D array of residuals (nA x nB x T)
//' @param G Interaction matrix (RA x RB)
//' @param rho_uv AR(1) persistence of the drift
//' @param sigma_uv Innovation standard deviation of the drift
//' @param s2 Dyadic variance
//' @param kappa Standard deviation of a reset position
//' @param pi_snap Reset probability
//' @param delta_u_current Current row resets (nA x T)
//' @param delta_v_current Current column resets (nB x T)
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite values disable it
//' @return List with U, V, delta_u, delta_v and bound_rejections
// [[Rcpp::export]]
List rUV_dynamic_snap_bip_fc_cpp(const arma::cube& U_state, const arma::cube& V_state,
	const arma::cube& E, const arma::mat& G, double rho_uv, double sigma_uv, double s2, double kappa,
	double pi_snap, const arma::mat& delta_u_current, const arma::mat& delta_v_current,
	double max_abs = -1.0) {
	arma::cube U = U_state;
	arma::cube V = V_state;
	const int nA = U.n_rows, nB = V.n_rows, T = U.n_slices;
	if (!U.is_finite() || !V.is_finite() || !G.is_finite() || G.n_rows != U.n_cols || G.n_cols != V.n_cols ||
		(int)V.n_slices != T || (int)delta_u_current.n_rows != nA || (int)delta_v_current.n_rows != nB ||
		(int)delta_u_current.n_cols != T || (int)delta_v_current.n_cols != T)
		Rcpp::stop("Invalid bipartite snap position state.");
	const SnapPar p = snap_par(rho_uv, sigma_uv, s2, kappa, pi_snap);
	const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
	int rejections = 0;
	arma::mat du(nA, T, fill::zeros), dv(nB, T, fill::zeros);
	for (int t = 0; t < T; t++) {
		arma::mat E_t = E.slice(t);
		E_t.replace(datum::nan, 0.0);
		// rows on W_t = V_t G'
		arma::mat W = V.slice(t) * G.t();
		arma::mat WtW = W.t() * W;
		arma::mat EW = E_t * W;
		for (int i = 0; i < nA; i++) {
			snap_row(U, i, t, WtW, EW.row(i).t(), p, delta_u_current, du, bounded, max_abs, rejections);
		}
		// columns on Q_t = U_t G
		arma::mat Q = U.slice(t) * G;
		arma::mat QtQ = Q.t() * Q;
		arma::mat EQ = E_t.t() * Q;
		for (int j = 0; j < nB; j++) {
			snap_row(V, j, t, QtQ, EQ.row(j).t(), p, delta_v_current, dv, bounded, max_abs, rejections);
		}
	}
	return List::create(Named("U") = U, Named("V") = V, Named("delta_u") = du, Named("delta_v") = dv,
		Named("bound_rejections") = rejections);
}

//' Symmetric snap-shift dynamic UV Gibbs update
//'
//' @param U_state 3D array of unit-variance paths (n x R x T)
//' @param Ld Weights of the product U_t diag(Ld) U_t'
//' @param E 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) persistence of the drift
//' @param sigma_uv Innovation standard deviation of the drift
//' @param s2 Dyadic variance
//' @param kappa Standard deviation of a reset position
//' @param pi_snap Reset probability
//' @param delta_current Current resets (n x T)
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite values disable it
//' @return List with U, V (equal to U), delta_u, delta_v (equal to delta_u) and bound_rejections
// [[Rcpp::export]]
List rUV_dynamic_snap_sym_fc_cpp(const arma::cube& U_state, const arma::vec& Ld, const arma::cube& E,
	double rho_uv, double sigma_uv, double s2, double kappa, double pi_snap, const arma::mat& delta_current,
	double max_abs = -1.0) {
	arma::cube U = U_state;
	const int n = U.n_rows, R = U.n_cols, T = U.n_slices;
	if (!U.is_finite() || !Ld.is_finite() || (int)Ld.n_elem != R || (int)delta_current.n_rows != n ||
		(int)delta_current.n_cols != T)
		Rcpp::stop("Invalid symmetric snap position state.");
	const SnapPar p = snap_par(rho_uv, sigma_uv, s2, kappa, pi_snap);
	const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
	int rejections = 0;
	arma::mat du(n, T, fill::zeros);
	const arma::rowvec lw = Ld.t();
	for (int t = 0; t < T; t++) {
		const arma::mat& E_t = E.slice(t);
		arma::mat Vs = U.slice(t).each_row() % lw;
		arma::mat VtVs = Vs.t() * Vs;
		arma::mat Es = E_t * Vs;
		for (int i = 0; i < n; i++) {
			const arma::rowvec vi = Vs.row(i);
			arma::vec ei = Es.row(i).t() - E_t(i, i) * vi.t();
			arma::mat P = VtVs - vi.t() * vi;
			snap_row(U, i, t, P, ei, p, delta_current, du, bounded, max_abs, rejections);
			// carry the regression sums to the next row
			const arma::rowvec vnew = U.slice(t).row(i) % lw;
			VtVs += vnew.t() * vnew - vi.t() * vi;
			Es += E_t.col(i) * (vnew - vi);
			Vs.row(i) = vnew;
		}
	}
	return List::create(Named("U") = U, Named("V") = U, Named("delta_u") = du, Named("delta_v") = du,
		Named("bound_rejections") = rejections);
}
