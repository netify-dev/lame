// The position sweep of the dynamic empirical-Bayes estimator for binary panels
// (R/lame_dynamic_als_eb_pos.R): the update of the mean-field Gaussian posterior of
// the moving latent positions of a directed or bipartite network, one side given the
// other. lde_pos_sweep_cpp is the arithmetic of .lde_pos_sweep_ref() (the R version,
// kept as the reference of the tests): per actor, the precision of the stacked path is
// the ar(1) prior plus, period by period, the weighted second moments of the other
// side's positions, and the mean, the diagonal and the lag-one blocks of the
// covariance come from one cholesky factor of it.
// Column-major storage throughout: means n x R x Tt, covariance blocks R x R x n x Tt
// and lag-one blocks R x R x n x (Tt - 1).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

// the precision, for unit innovation variance, of a stationary ar(1) path of length Tt
arma::mat ar_precision(int Tt, double rho) {
	arma::mat Q(Tt, Tt, arma::fill::zeros);
	if (Tt == 1) {
		Q(0, 0) = 1.0 - rho * rho;
		return Q;
	}
	for (int t = 0; t < Tt; ++t) Q(t, t) = (t == 0 || t == Tt - 1) ? 1.0 : 1.0 + rho * rho;
	for (int t = 0; t + 1 < Tt; ++t) {
		Q(t, t + 1) = -rho;
		Q(t + 1, t) = -rho;
	}
	return Q;
}

// the upper cholesky factor of S; on failure the matrix is symmetrised and a ridge
// (1e-10 of the larger of its diagonal and one, times 100 at every retry) is added, as the
// information blocks of the R code are: they are positive definite by construction and a
// failure is rounding in a near-degenerate block
bool chol_ridge(arma::mat& U, const arma::mat& S) {
	if (arma::chol(U, S, "upper")) return true;
	arma::mat Ss = 0.5 * (S + S.t());
	double dmax = 1.0;
	for (arma::uword i = 0; i < Ss.n_rows; ++i) dmax = std::max(dmax, std::fabs(Ss(i, i)));
	double ridge = 1e-10 * dmax;
	for (int k = 0; k < 12; ++k) {
		arma::mat Sr = Ss;
		Sr.diag() += ridge;
		if (arma::chol(U, Sr, "upper")) return true;
		ridge *= 100.0;
	}
	return false;
}

// the same for a scale-mixture transition: the step into period t (t >= 2) has the weight w[t]; the
// innovation variance is s2 / w[t]. d_1 = (1 - rho^2) + rho^2 w_2, d_t = w_t + rho^2 w_t+1, d_Tt = w_Tt,
// q_t,t+1 = -rho w_t+1 (w is the actor's row of the weights, indexed by period)
arma::mat ar_precision_w(int Tt, double rho, const std::vector<double>& w) {
	arma::mat Q(Tt, Tt, arma::fill::zeros);
	if (Tt == 1) {
		Q(0, 0) = 1.0 - rho * rho;
		return Q;
	}
	Q(0, 0) = (1.0 - rho * rho) + rho * rho * w[1];
	for (int t = 1; t + 1 < Tt; ++t) Q(t, t) = w[t] + rho * rho * w[t + 1];
	Q(Tt - 1, Tt - 1) = w[Tt - 1];
	for (int t = 0; t + 1 < Tt; ++t) {
		Q(t, t + 1) = -rho * w[t + 1];
		Q(t + 1, t) = -rho * w[t + 1];
	}
	return Q;
}

struct Side {
	NumericVector mean, sig, cl;
	double ldp, trp, ldp0;
	arma::mat P0;
};

// q of one side's positions given the other's. The cell (r, o, t) of the weight and
// residual arrays sits at r * sr + o * so + t * st. w, z, eta, obs are the cell arrays of
// the whole panel (a cell counts where observed and finite)
Side side_update(const double* w, const double* z, const double* eta, const int* obs,
                 R_xlen_t sr, R_xlen_t so, R_xlen_t st, int n, int m, int R, int Tt,
                 const double* vh, const double* om, double s2, double rho, const double* tw) {
	const int d = R * Tt;
	const int npair = R * (R + 1) / 2;
	Side out;
	out.mean = NumericVector((R_xlen_t) n * R * Tt);
	out.sig = NumericVector((R_xlen_t) R * R * n * Tt);
	out.cl = NumericVector((R_xlen_t) R * R * n * std::max(Tt - 1, 1));
	arma::mat P0 = arma::kron(ar_precision(Tt, rho) / s2, arma::eye(R, R));
	out.P0 = P0;
	arma::mat U0;
	if (!chol_ridge(U0, P0)) stop("The prior precision of a position path is not positive definite.");
	out.ldp0 = n * 2.0 * arma::sum(arma::log(U0.diag()));
	out.ldp = 0.0;
	out.trp = 0.0;
	// the weighted second moments of the other side, by actor, pair and period, and the
	// weighted targets
	std::vector<double> lam((size_t) npair * n * Tt, 0.0), hh((size_t) R * n * Tt, 0.0);
	std::vector<double> ev(npair), v(R);
	for (int t = 0; t < Tt; ++t) {
		for (int o = 0; o < m; ++o) {
			for (int k = 0; k < R; ++k) v[k] = vh[o + (R_xlen_t) m * (k + (R_xlen_t) R * t)];
			int q = 0;
			for (int k = 0; k < R; ++k) {
				for (int l = k; l < R; ++l, ++q) {
					ev[q] = v[k] * v[l] + om[k + (R_xlen_t) R * (l + (R_xlen_t) R * (o + (R_xlen_t) m * t))];
				}
			}
			for (int r = 0; r < n; ++r) {
				const R_xlen_t idx = (R_xlen_t) r * sr + (R_xlen_t) o * so + (R_xlen_t) t * st;
				if (!obs[idx]) continue;
				const double wc = w[idx];
				if (!std::isfinite(wc) || wc == 0.0) continue;
				double res = z[idx] - eta[idx];
				if (!std::isfinite(res)) res = 0.0;
				double* lr = &lam[((size_t) r + (size_t) n * t) * npair];
				for (int q = 0; q < npair; ++q) lr[q] += wc * ev[q];
				double* hr = &hh[((size_t) r + (size_t) n * t) * R];
				const double wr = wc * res;
				for (int k = 0; k < R; ++k) hr[k] += wr * v[k];
			}
		}
	}
	arma::mat Pi(d, d), Ui, Linv, Sg, P0i;
	arma::vec h(d), mu;
	std::vector<double> wi(Tt, 1.0);
	for (int i = 0; i < n; ++i) {
		// this actor's prior precision: the common one, or the one of its transition weights
		if (tw != nullptr) {
			for (int t = 0; t < Tt; ++t) wi[t] = tw[i + (R_xlen_t) n * t];
			P0i = arma::kron(ar_precision_w(Tt, rho, wi) / s2, arma::eye(R, R));
		}
		const arma::mat& P0a = (tw != nullptr) ? P0i : P0;
		Pi = P0a;
		for (int t = 0; t < Tt; ++t) {
			const double* lr = &lam[((size_t) i + (size_t) n * t) * npair];
			int q = 0;
			for (int k = 0; k < R; ++k) {
				for (int l = k; l < R; ++l, ++q) {
					Pi(t * R + k, t * R + l) += lr[q];
					if (l != k) Pi(t * R + l, t * R + k) += lr[q];
				}
			}
			const double* hr = &hh[((size_t) i + (size_t) n * t) * R];
			for (int k = 0; k < R; ++k) h(t * R + k) = hr[k];
		}
		if (!chol_ridge(Ui, Pi)) stop("A dynamic empirical-Bayes information block is not positive definite.");
		out.ldp += 2.0 * arma::sum(arma::log(Ui.diag()));
		Linv = arma::inv(arma::trimatu(Ui));
		Sg = Linv * Linv.t();
		out.trp += arma::accu(P0a % Sg);
		mu = Sg * h;
		for (int t = 0; t < Tt; ++t) {
			for (int k = 0; k < R; ++k) out.mean[i + (R_xlen_t) n * (k + (R_xlen_t) R * t)] = mu(t * R + k);
			for (int l = 0; l < R; ++l) {
				for (int k = 0; k < R; ++k) {
					out.sig[k + (R_xlen_t) R * (l + (R_xlen_t) R * (i + (R_xlen_t) n * t))] = Sg(t * R + k, t * R + l);
					if (t + 1 < Tt) {
						out.cl[k + (R_xlen_t) R * (l + (R_xlen_t) R * (i + (R_xlen_t) n * t))] = Sg(t * R + k, (t + 1) * R + l);
					}
				}
			}
		}
	}
	return out;
}

List side_list(const Side& s, int n, int R, int Tt) {
	NumericVector mean = clone(s.mean), sig = clone(s.sig), cl = clone(s.cl);
	mean.attr("dim") = IntegerVector::create(n, R, Tt);
	sig.attr("dim") = IntegerVector::create(R, R, n, Tt);
	cl.attr("dim") = IntegerVector::create(R, R, n, std::max(Tt - 1, 1));
	return List::create(Named("mean") = mean, Named("sig") = sig, Named("cl") = cl,
	                    Named("ldp") = s.ldp, Named("trp") = s.trp, Named("ldp0") = s.ldp0,
	                    Named("P0") = s.P0);
}

}  // namespace

// One sweep over both sides of a directed or bipartite network: q(U) given q(V), then q(V)
// given the fresh q(U). w, z, eta_add, obs are the na x nb x Tt arrays of the working
// weights, working responses, additive predictor and observation indicator; the incoming
// positions are given by their means and covariance blocks.
// [[Rcpp::export(rng = false)]]
List lde_pos_sweep_cpp(const NumericVector& w, const NumericVector& z, const NumericVector& eta_add,
                       const IntegerVector& obs, const NumericVector& umean, const NumericVector& usig,
                       const NumericVector& vmean, const NumericVector& vsig,
                       int na, int nb, int R, int Tt, double s2, double rho,
                       const NumericVector& lam_u, const NumericVector& lam_v) {
	const R_xlen_t ncell = (R_xlen_t) na * nb * Tt;
	if (w.size() != ncell || z.size() != ncell || eta_add.size() != ncell || obs.size() != ncell) {
		stop("The cell arrays do not match the dimensions of the panel.");
	}
	if (vmean.size() != (R_xlen_t) nb * R * Tt || vsig.size() != (R_xlen_t) R * R * nb * Tt) {
		stop("The receiver positions do not match the dimensions of the panel.");
	}
	if (umean.size() != (R_xlen_t) na * R * Tt || usig.size() != (R_xlen_t) R * R * na * Tt) {
		stop("The sender positions do not match the dimensions of the panel.");
	}
	// the transition weights of a scale-mixture transition (n x Tt), none for the Gaussian one
	if (lam_u.size() != 0 && lam_u.size() != (R_xlen_t) na * Tt) stop("The sender transition weights do not match the panel.");
	if (lam_v.size() != 0 && lam_v.size() != (R_xlen_t) nb * Tt) stop("The receiver transition weights do not match the panel.");
	// senders: rows i, others j; cell (i, j, t) at i + na j + na nb t
	Side su = side_update(w.begin(), z.begin(), eta_add.begin(), obs.begin(), 1, na, (R_xlen_t) na * nb,
	                      na, nb, R, Tt, vmean.begin(), vsig.begin(), s2, rho,
	                      lam_u.size() == 0 ? nullptr : lam_u.begin());
	// receivers: rows j, others i
	Side sv = side_update(w.begin(), z.begin(), eta_add.begin(), obs.begin(), na, 1, (R_xlen_t) na * nb,
	                      nb, na, R, Tt, su.mean.begin(), su.sig.begin(), s2, rho,
	                      lam_v.size() == 0 ? nullptr : lam_v.begin());
	return List::create(Named("qu") = side_list(su, na, R, Tt), Named("qv") = side_list(sv, nb, R, Tt));
}
