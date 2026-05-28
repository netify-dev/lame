// Dynamic regression coefficients for LAME.
//
// This file ships five Rcpp-exported helpers that implement an AR(1)
// state-space model on (a subset of) the regression coefficients beta:
//
//   beta_b,t = rho_b * beta_b,t-1 + N(0, sigma_b^2 * Lambda_b)
//
// where b indexes a coefficient BLOCK (intercept / dyad / row / col) and
// t indexes time periods. The samplers are:
//
//   1) sample_beta_dynamic_cpp   -- joint forward-filter / backward-sample
//                                   draw of the entire (T x p_dyn) path,
//                                   conditional on the static block, the
//                                   per-period sufficient statistics, and
//                                   the AR(1) parameters.
//   2) sample_beta_static_cpp    -- conjugate Gaussian update of the static
//                                   block, conditional on the dynamic path.
//   3) sample_rho_beta_cpp       -- per-block truncated-Normal full
//                                   conditional for the AR(1) rho parameter.
//   4) sample_sigma_beta_cpp     -- per-block inverse-Gamma full conditional
//                                   for the AR(1) innovation variance.
//   5) get_EZ_dynamic_beta_cpp   -- rebuild EZ when beta is time-varying,
//                                   per-period beta_t feeds into a unipartite
//                                   or bipartite EZ slice. (Used for both
//                                   inside-MCMC EZ refresh and at
//                                   fit-assembly time.)
//
// Compiles with the rest of src/. Byte-identical default: none of these
// functions runs when the user does not pass dynamic_beta=TRUE.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// --- shared utilities --------------------------------------------------------

// Symmetric positive-definite inverse with progressive jitter. Returns true
// on success and writes the inverse into `out`; returns false on total
// failure (so the caller can bump a failure counter and use a fallback).
// The jitter ladder lets us tolerate the same near-singular Omega situations
// that the existing rbeta_ab_rep_fc_cpp handles via its "ridge" line.
static bool safe_sympd_inverse(const arma::mat& M, arma::mat& out,
                               int* n_jitter_used = NULL) {
	arma::mat A = 0.5 * (M + M.t());
	double dmean = std::max(1e-12, std::abs(arma::trace(A)) /
	                                  std::max<arma::uword>(1, A.n_rows));
	double jitters[5] = {0.0, 1e-10, 1e-7, 1e-4, 1e-2};
	for (int k = 0; k < 5; ++k) {
		arma::mat A2 = A;
		if (jitters[k] > 0.0) A2.diag() += jitters[k] * dmean;
		bool ok = arma::inv_sympd(out, A2);
		if (ok && out.is_finite()) {
			if (n_jitter_used != NULL && k > 0) *n_jitter_used = k;
			return true;
		}
	}
	// final fallback: general inv
	bool ok = arma::inv(out, A + arma::eye<arma::mat>(A.n_rows, A.n_cols) *
	                         (1e-2 * dmean));
	if (n_jitter_used != NULL) *n_jitter_used = 5;
	return ok && out.is_finite();
}

// Draw N(mu, Sigma) using a Cholesky of Sigma. On Cholesky failure, fall
// back to symmetric eigen and clamp negative eigenvalues -- matches the
// degraded behaviour of the existing rmvnorm_cpp without ever returning
// NaNs from a singular Sigma. Increments *n_chol_fail when the cholesky
// fallback fires (we expose this through tryErrorChecks).
static arma::vec rmvnorm_chol(const arma::vec& mu, const arma::mat& Sigma,
                              int* n_chol_fail = NULL) {
	int p = mu.n_elem;
	arma::vec eps(p, fill::randn);
	arma::mat S = 0.5 * (Sigma + Sigma.t());
	arma::mat L;
	bool ok = arma::chol(L, S, "lower");
	if (!ok) {
		arma::vec d; arma::mat U2;
		bool eok = arma::eig_sym(d, U2, S);
		if (!eok) {
			if (n_chol_fail != NULL) (*n_chol_fail)++;
			return mu; // give up; return the mean
		}
		for (int k = 0; k < (int)d.n_elem; ++k) {
			if (d(k) < 0.0 || !std::isfinite(d(k))) d(k) = 0.0;
		}
		L = U2 * arma::diagmat(arma::sqrt(d));
		if (n_chol_fail != NULL) (*n_chol_fail)++;
	}
	return mu + L * eps;
}

// Build the per-period dyadic "weight precision" coefficient. For
// non-dyad-correlated families and bipartite/symmetric cases this is just
// 1/s2. For unipartite directed with dyadic correlation, the dyad-level
// likelihood reduces to a 2x2 Gaussian on (Z_ij, Z_ji); the implied
// per-coefficient sufficient statistics use the inverse-2x2 covariance,
// which the caller passes pre-decomposed as (c_diag, c_off):
//   c_diag = 1 / (s2 * (1 - rho^2))
//   c_off  = -rho / (s2 * (1 - rho^2))
// The flag use_dyad_rho tells the sampler which branch to take.

// Per-period sufficient statistics for the dynamic block:
//   Omega_t = X_D,t' W_t X_D,t   (p_dyn x p_dyn)
//   h_t     = X_D,t' W_t (Z_t - offset_t - X_S,t beta_S)   (p_dyn)
// where W_t encodes the dyadic-correlation precision. `bipartite` and
// `symmetric` toggle whether the design is folded across i!=j (unipartite
// directed) or over half-triangle (symmetric); for bipartite there's no
// pairing.
static void accumulate_Omega_h(const arma::mat& Xdyn_t,
                               const arma::mat& Xstat_t,
                               const arma::mat& Z_t,
                               const arma::mat& offset_t,
                               const arma::vec& beta_static,
                               double s2, double dyad_rho,
                               bool bipartite, bool symmetric,
                               bool use_dyad_rho,
                               arma::mat& Omega_out, arma::vec& h_out) {
	const int n_obs = Xdyn_t.n_rows;
	const int p_dyn = Xdyn_t.n_cols;
	// y_D,t = Z_t - offset_t - X_S,t beta_S, as a long vector aligned with
	// the rows of Xdyn_t (which is also a long vector reshape).
	arma::vec y = arma::vectorise(Z_t) - arma::vectorise(offset_t);
	if (Xstat_t.n_cols > 0 && beta_static.n_elem == Xstat_t.n_cols) {
		y -= Xstat_t * beta_static;
	}
	// drop NA observations (treated as missing)
	arma::uvec keep = arma::find_finite(y);
	arma::vec y_obs = y.elem(keep);
	arma::mat X_obs = Xdyn_t.rows(keep);
	if (!use_dyad_rho || bipartite || symmetric) {
		// W = (1/s2) * I on the observed dyads
		double w = 1.0 / s2;
		Omega_out = w * X_obs.t() * X_obs;
		h_out     = w * X_obs.t() * y_obs;
	} else {
		// directed unipartite with dyadic correlation: we need to pair
		// (i,j) with (j,i). The caller passed Xdyn_t / Xstat_t already
		// vectorised in n*n order with row-major storage (i.e. arma's
		// arma::vectorise of a [n x n] matrix), and the (i,j) entry sits
		// at index k = (j-1)*n + i. So the pair partner of index k is the
		// index for (j,i), which is (i-1)*n + j. We accumulate the 2x2
		// inverse-covariance Q on each ORDERED pair (i,j) with i < j once.
		const int n = Z_t.n_rows;
		if ((int)Z_t.n_cols != n) {
			// fallback to isotropic if we somehow get a non-square unipartite
			double w = 1.0 / s2;
			Omega_out = w * X_obs.t() * X_obs;
			h_out     = w * X_obs.t() * y_obs;
			return;
		}
		double det_w = s2 * s2 * (1.0 - dyad_rho * dyad_rho);
		// guard near-singular det
		if (!std::isfinite(det_w) || det_w < 1e-12) det_w = 1e-12;
		double c_diag = s2 / det_w;           // = 1 / (s2 (1-rho^2))
		double c_off  = -dyad_rho * s2 / det_w; // = -rho/(s2(1-rho^2))
		// build vectorised offset
		arma::mat E_t = Z_t - offset_t;
		if (Xstat_t.n_cols > 0 && beta_static.n_elem == Xstat_t.n_cols) {
			arma::mat Xs_beta_mat = arma::reshape(Xstat_t * beta_static, n, n);
			E_t -= Xs_beta_mat;
		}
		Omega_out.zeros(p_dyn, p_dyn);
		h_out.zeros(p_dyn);
		// iterate unordered pairs only so each dyad is processed exactly
		// once. otherwise a half-observed dyad (y_ij obs, y_ji NA) would
		// be counted twice because both (i, j) and (j, i) hit the NA
		// branch and re-add the observed half.
		double w_iso = 1.0 / s2;
		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				// per-cell row of the long-format X (column-major: idx = j*n + i)
				int idx_ij = j * n + i;
				int idx_ji = i * n + j;
				double y_ij = E_t(i, j);
				double y_ji = E_t(j, i);
				bool obs_ij = std::isfinite(y_ij);
				bool obs_ji = std::isfinite(y_ji);
				if (obs_ij && obs_ji) {
					// full 2x2 dyad-correlated contribution
					arma::rowvec x_ij = Xdyn_t.row(idx_ij);
					arma::rowvec x_ji = Xdyn_t.row(idx_ji);
					Omega_out += c_diag * (x_ij.t() * x_ij + x_ji.t() * x_ji);
					Omega_out += c_off  * (x_ij.t() * x_ji + x_ji.t() * x_ij);
					h_out     += c_diag * (x_ij.t() * y_ij + x_ji.t() * y_ji);
					h_out     += c_off  * (x_ij.t() * y_ji + x_ji.t() * y_ij);
				} else if (obs_ij) {
					// only y_ij observed: isotropic on that half-dyad
					arma::rowvec x_ij = Xdyn_t.row(idx_ij);
					Omega_out += w_iso * x_ij.t() * x_ij;
					h_out     += w_iso * x_ij.t() * y_ij;
				} else if (obs_ji) {
					// only y_ji observed: isotropic on that half-dyad
					arma::rowvec x_ji = Xdyn_t.row(idx_ji);
					Omega_out += w_iso * x_ji.t() * x_ji;
					h_out     += w_iso * x_ji.t() * y_ji;
				}
				// both missing: no contribution
			}
		}
	}
}

// --- 1) FFBS draw of the dynamic-block beta path -----------------------------

//' Sample the dynamic-block beta path via FFBS
//'
//' Forward-filter / backward-sample the AR(1) state-space model for the
//' dynamic-block beta coefficients. Returns the joint draw of beta_dyn at
//' every time period.
//'
//' @param Xdyn_list T-length list of (n*n) x p_dyn long-format design
//'   matrices for the dynamic block (column-major reshape per period).
//' @param Xstat_list T-length list of (n*n) x p_static long-format design
//'   matrices for the static block.
//' @param Z_list T-length list of (n x n) latent matrices.
//' @param offset_list T-length list of (n x n) offset matrices (a_i + b_j +
//'   U_i'V_j contributions; everything that's not in X*beta).
//' @param beta_static Length p_static current static beta vector.
//' @param rho_by_coef Length p_dyn vector of AR(1) rho values for each
//'   dynamic coefficient. (Per-block but expanded per-column for vectorised
//'   indexing.)
//' @param sigma_by_coef Length p_dyn vector of AR(1) innovation standard
//'   deviations for each dynamic coefficient.
//' @param Lambda Block-diagonal innovation scale matrix (p_dyn x p_dyn).
//'   Combined with sigma^2 to give Q = sigma^2 * Lambda.
//' @param beta0_mean Length p_dyn prior mean for beta at t=0.
//' @param beta0_cov  p_dyn x p_dyn prior covariance for beta at t=0.
//' @param s2 Dyadic variance.
//' @param dyad_rho Dyadic correlation (ignored when use_dyad_rho=FALSE).
//' @param bipartite Whether the network is bipartite.
//' @param symmetric Whether the network is symmetric.
//' @param use_dyad_rho Whether to use the dyad-corr branch (TRUE only for
//'   unipartite, asymmetric, with a non-zero rho).
//' @return List with: path -- a (T x p_dyn) matrix of beta draws (one row
//'   per period); chol_fail -- integer count of Cholesky failures.
// [[Rcpp::export]]
List sample_beta_dynamic_cpp(List Xdyn_list, List Xstat_list,
                             List Z_list, List offset_list,
                             arma::vec beta_static,
                             arma::vec rho_by_coef,
                             arma::vec sigma_by_coef,
                             arma::mat Lambda,
                             arma::vec beta0_mean,
                             arma::mat beta0_cov,
                             double s2, double dyad_rho,
                             bool bipartite, bool symmetric,
                             bool use_dyad_rho) {
	const int T = Xdyn_list.size();
	if (T == 0) {
		return List::create(Named("path") = arma::mat(0, 0),
		                    Named("chol_fail") = 0);
	}
	arma::mat X0 = as<arma::mat>(Xdyn_list[0]);
	const int p_dyn = X0.n_cols;
	if (p_dyn == 0) {
		return List::create(Named("path") = arma::mat(T, 0),
		                    Named("chol_fail") = 0);
	}
	// AR(1) state evolution: F = diag(rho_by_coef), Q = sigma^2 * Lambda
	// (sigma_by_coef is one entry per coef expanded from per-block sigmas)
	arma::mat F = arma::diagmat(rho_by_coef);
	arma::mat Q = arma::diagmat(sigma_by_coef % sigma_by_coef);
	Q = Q * Lambda;
	Q = 0.5 * (Q + Q.t());

	// storage for forward filter
	std::vector<arma::vec> m_filt(T);
	std::vector<arma::mat> C_filt(T);
	std::vector<arma::vec> m_pred(T);
	std::vector<arma::mat> C_pred(T);

	int chol_fail = 0;

	// when ρ is close to 1 the stationary AR(1) variance σ²Λ/(1-ρ²)
	// dominates and the supplied beta0_cov (default diag(10)) would be
	// tighter than the implied stationary prior, anchoring the t=0
	// state. inflate beta0_cov to max(beta0_cov, stationary variance)
	// element-wise on the diagonal so the prior is at least as diffuse as the
	// AR(1) stationary distribution implies.
	arma::mat beta0_cov_eff = beta0_cov;
	{
		arma::mat F2 = F * F;  // diag(rho^2)
		arma::vec diag_F2 = F2.diag();
		// stationary diagonal variance = diag(Q) / (1 - rho^2), clamped at floor
		// to avoid divide-by-near-zero when rho ~ 0.999.
		arma::vec diag_Q = Q.diag();
		arma::vec stat_var(p_dyn);
		for (int k = 0; k < p_dyn; ++k) {
			double denom = std::max(1.0 - diag_F2(k), 1e-3);
			stat_var(k) = diag_Q(k) / denom;
		}
		// inflate t=0 prior variance only where stationary > supplied
		for (int k = 0; k < p_dyn; ++k) {
			if (stat_var(k) > beta0_cov_eff(k, k)) {
				beta0_cov_eff(k, k) = stat_var(k);
			}
		}
	}

	for (int t = 0; t < T; ++t) {
		arma::vec m_prev;
		arma::mat C_prev;
		if (t == 0) {
			m_prev = beta0_mean;
			C_prev = beta0_cov_eff;
		} else {
			m_prev = m_filt[t-1];
			C_prev = C_filt[t-1];
		}
		// predict: m_t|t-1 = F m_prev, C_t|t-1 = F C_prev F' + Q
		arma::vec m_p = F * m_prev;
		arma::mat C_p = F * C_prev * F.t() + Q;
		C_p = 0.5 * (C_p + C_p.t());
		m_pred[t] = m_p;
		C_pred[t] = C_p;
		// accumulate Omega_t and h_t
		arma::mat Omega_t(p_dyn, p_dyn, fill::zeros);
		arma::vec h_t(p_dyn, fill::zeros);
		arma::mat Xdyn_t = as<arma::mat>(Xdyn_list[t]);
		arma::mat Xstat_t;
		if (Xstat_list.size() > t) {
			SEXP sx = Xstat_list[t];
			if (!Rf_isNull(sx)) Xstat_t = as<arma::mat>(sx);
		}
		arma::mat Z_t = as<arma::mat>(Z_list[t]);
		arma::mat offset_t = as<arma::mat>(offset_list[t]);
		accumulate_Omega_h(Xdyn_t, Xstat_t, Z_t, offset_t, beta_static,
		                   s2, dyad_rho, bipartite, symmetric, use_dyad_rho,
		                   Omega_t, h_t);
		// posterior precision K = inv(C_pred) + Omega, then C_filt = inv(K)
		arma::mat C_p_inv;
		if (!safe_sympd_inverse(C_p, C_p_inv)) {
			chol_fail++;
			C_p_inv = arma::diagmat(arma::ones<arma::vec>(p_dyn) * 1e-6);
		}
		arma::mat K = C_p_inv + Omega_t;
		arma::mat C_f;
		if (!safe_sympd_inverse(K, C_f)) {
			chol_fail++;
			C_f = C_p; // fallback to predictive
		}
		arma::vec m_f = C_f * (C_p_inv * m_p + h_t);
		C_filt[t] = C_f;
		m_filt[t] = m_f;
	}
	// backward sample
	arma::mat path(T, p_dyn);
	// draw beta_T from N(m_filt[T-1], C_filt[T-1])
	arma::vec beta_t = rmvnorm_chol(m_filt[T-1], C_filt[T-1], &chol_fail);
	path.row(T-1) = beta_t.t();
	for (int t = T - 2; t >= 0; --t) {
		// smoothing: J_t = C_filt[t] F' inv(C_pred[t+1])
		arma::mat C_p_next_inv;
		if (!safe_sympd_inverse(C_pred[t+1], C_p_next_inv)) {
			chol_fail++;
			C_p_next_inv = arma::eye<arma::mat>(p_dyn, p_dyn) * 1e-6;
		}
		arma::mat J_t = C_filt[t] * F.t() * C_p_next_inv;
		arma::vec m_smooth = m_filt[t] + J_t * (beta_t - m_pred[t+1]);
		arma::mat C_smooth = C_filt[t] - J_t * F * C_filt[t];
		C_smooth = 0.5 * (C_smooth + C_smooth.t());
		beta_t = rmvnorm_chol(m_smooth, C_smooth, &chol_fail);
		path.row(t) = beta_t.t();
	}
	return List::create(Named("path") = path,
	                    Named("chol_fail") = chol_fail);
}

// --- 2) static-block conjugate update ----------------------------------------

//' Sample the static-block beta conditional on the dynamic path
//'
//' Conjugate Gaussian update for the static coefficient block, treating the
//' dynamic path as known. Uses a flat-ish ridge-style prior `prior_prec` on
//' the static beta (typically diag(1/g) to match the existing g-prior path).
//'
//' @param Xdyn_list T-length list of long-format dynamic design matrices.
//' @param Xstat_list T-length list of long-format static design matrices.
//' @param Z_list T-length list of latent (n x n) matrices.
//' @param offset_list T-length list of (n x n) offsets.
//' @param beta_dyn_path (T x p_dyn) dynamic path matrix.
//' @param prior_mean Length p_static prior mean.
//' @param prior_prec p_static x p_static prior precision.
//' @param s2 Dyadic variance.
//' @param dyad_rho Dyadic correlation.
//' @param bipartite TRUE/FALSE.
//' @param symmetric TRUE/FALSE.
//' @param use_dyad_rho TRUE/FALSE.
//' @return List with: beta -- length p_static; chol_fail -- integer.
// [[Rcpp::export]]
List sample_beta_static_cpp(List Xdyn_list, List Xstat_list,
                            List Z_list, List offset_list,
                            arma::mat beta_dyn_path,
                            arma::vec prior_mean,
                            arma::mat prior_prec,
                            double s2, double dyad_rho,
                            bool bipartite, bool symmetric,
                            bool use_dyad_rho) {
	const int T = Z_list.size();
	if (T == 0 || Xstat_list.size() == 0) {
		return List::create(Named("beta") = arma::vec(arma::uword(0)),
		                    Named("chol_fail") = 0);
	}
	arma::mat Xs0 = as<arma::mat>(Xstat_list[0]);
	const int p_static = Xs0.n_cols;
	if (p_static == 0) {
		return List::create(Named("beta") = arma::vec(arma::uword(0)),
		                    Named("chol_fail") = 0);
	}
	arma::mat Omega_post = prior_prec;
	arma::vec h_post     = prior_prec * prior_mean;
	int chol_fail = 0;
	for (int t = 0; t < T; ++t) {
		arma::mat Xs_t = as<arma::mat>(Xstat_list[t]);
		arma::mat Z_t  = as<arma::mat>(Z_list[t]);
		arma::mat off_t = as<arma::mat>(offset_list[t]);
		// adjust offset for the dynamic-block contribution at period t
		arma::mat Xdyn_t = as<arma::mat>(Xdyn_list[t]);
		if (Xdyn_t.n_cols > 0) {
			arma::vec b_dyn_t = beta_dyn_path.row(t).t();
			int n = Z_t.n_rows; int m = Z_t.n_cols;
			arma::vec contrib = Xdyn_t * b_dyn_t;
			off_t = off_t + arma::reshape(contrib, n, m);
		}
		arma::mat Omega_t(p_static, p_static, fill::zeros);
		arma::vec h_t(p_static, fill::zeros);
		// Build the "static treated as the only design" view: pass an
		// empty Xdyn-effective and zero beta-static
		arma::mat Xempty(Xs_t.n_rows, 0);
		arma::vec b_empty(arma::uword(0));
		accumulate_Omega_h(Xs_t, Xempty, Z_t, off_t, b_empty,
		                   s2, dyad_rho, bipartite, symmetric, use_dyad_rho,
		                   Omega_t, h_t);
		Omega_post += Omega_t;
		h_post     += h_t;
	}
	arma::mat Vpost;
	if (!safe_sympd_inverse(Omega_post, Vpost)) {
		chol_fail++;
		Vpost = arma::eye<arma::mat>(p_static, p_static);
	}
	arma::vec mpost = Vpost * h_post;
	arma::vec beta_new = rmvnorm_chol(mpost, Vpost, &chol_fail);
	return List::create(Named("beta") = beta_new,
	                    Named("chol_fail") = chol_fail);
}

// --- 3) rho_beta per-block update --------------------------------------------

//' Sample the AR(1) rho for each dynamic block
//'
//' Truncated-Normal full conditional, one rho per block.
//'
//' @param beta_path (T x p_dyn) matrix of beta draws.
//' @param group_id Length p_dyn integer vector (1-based) of block IDs.
//' @param n_groups Number of distinct block IDs.
//' @param Lambda_inv (p_dyn x p_dyn) inverse of the (full) Lambda scale.
//' @param sigma_by_coef Length p_dyn vector of per-coef sigma (block-shared).
//' @param rho_current Length n_groups vector of current rho values.
//' @param rho_prior_mean Length n_groups vector of prior means.
//' @param rho_prior_sd   Length n_groups vector of prior SDs.
//' @param rho_lower Lower truncation bound (typically 0).
//' @param rho_upper Upper truncation bound (typically 0.999).
//' @return Length n_groups vector of new rho values.
// [[Rcpp::export]]
arma::vec sample_rho_beta_cpp(arma::mat beta_path,
                              arma::ivec group_id,
                              int n_groups,
                              arma::mat Lambda_inv,
                              arma::vec sigma_by_coef,
                              arma::vec rho_current,
                              arma::vec rho_prior_mean,
                              arma::vec rho_prior_sd,
                              double rho_lower = 0.0,
                              double rho_upper = 0.999) {
	const int T = beta_path.n_rows;
	const int p = beta_path.n_cols;
	arma::vec rho_new = rho_current;
	if (T < 2 || p == 0) return rho_new;
	for (int g = 1; g <= n_groups; ++g) {
		arma::uvec cols = arma::find(group_id == g);
		if (cols.n_elem == 0) continue;
		arma::mat path_g = beta_path.cols(cols);
		// per-block Lambda_inv slice
		arma::mat Linv_g = Lambda_inv.submat(cols, cols);
		// share a sigma_g across this block; pull the first coef's sigma
		double sigma_g = sigma_by_coef(cols(0));
		double s2_g = sigma_g * sigma_g;
		// sufficient stats: A_g = sum_t beta_t-1' Linv_g beta_t-1
		//                   B_g = sum_t beta_t-1' Linv_g beta_t
		double A_g = 0.0, B_g = 0.0;
		for (int t = 1; t < T; ++t) {
			arma::vec bp = path_g.row(t-1).t();
			arma::vec bc = path_g.row(t).t();
			A_g += as_scalar(bp.t() * Linv_g * bp);
			B_g += as_scalar(bp.t() * Linv_g * bc);
		}
		double prior_v = rho_prior_sd(g - 1) * rho_prior_sd(g - 1);
		double inv_prior_v = 1.0 / prior_v;
		double post_v = 1.0 / (A_g / s2_g + inv_prior_v);
		double post_m = post_v * (B_g / s2_g + rho_prior_mean(g - 1) * inv_prior_v);
		// rejection-sample truncated normal. when 100 rejections happen
		// in a row (data-implied rho outside [rho_lower, rho_upper]) the
		// fallback (i) clamps post_m to the bounds (posterior mode within
		// the truncation, not the prior mean) and (ii) draws a final
		// candidate via inverse-CDF truncated-normal so the result is
		// still a draw from the truncated posterior rather than a point.
		double prop = post_m;
		double sd_post = std::sqrt(std::max(post_v, 1e-12));
		bool got = false;
		for (int k = 0; k < 100; ++k) {
			double cand = R::rnorm(post_m, sd_post);
			if (cand > rho_lower && cand < rho_upper) {
				prop = cand;
				got = true;
				break;
			}
		}
		if (!got) {
			// inverse-CDF truncated normal: Phi^{-1}(u * (F(hi) - F(lo)) + F(lo))
			// with mean post_m, sd sd_post, on the [rho_lower, rho_upper] window.
			double z_lo = (rho_lower - post_m) / sd_post;
			double z_hi = (rho_upper - post_m) / sd_post;
			double F_lo = R::pnorm(z_lo, 0.0, 1.0, 1, 0);
			double F_hi = R::pnorm(z_hi, 0.0, 1.0, 1, 0);
			double u = R::runif(0.0, 1.0) * (F_hi - F_lo) + F_lo;
			// guard against numerical Phi^-1(0) or Phi^-1(1)
			u = std::max(1e-12, std::min(1.0 - 1e-12, u));
			double z = R::qnorm(u, 0.0, 1.0, 1, 0);
			prop = post_m + sd_post * z;
			// final safety clamp -- inverse-CDF should land in (lo, hi) by
			// construction but be defensive against numerical edge cases
			prop = std::max(rho_lower + 1e-12, std::min(rho_upper - 1e-12, prop));
		}
		rho_new(g - 1) = prop;
	}
	return rho_new;
}

// --- 4) sigma_beta per-block update ------------------------------------------

//' Sample the AR(1) innovation sigma for each dynamic block
//'
//' Inverse-Gamma full conditional, one sigma per block.
//'
//' @param beta_path (T x p_dyn) matrix.
//' @param group_id Length p_dyn integer vector (1-based) of block IDs.
//' @param n_groups Number of distinct block IDs.
//' @param Lambda_inv (p_dyn x p_dyn) inverse of the (full) Lambda scale.
//' @param rho_by_group Length n_groups vector of current rho values.
//' @param prior_shape Length n_groups vector of IG shape parameters.
//' @param prior_scale Length n_groups vector of IG scale parameters.
//' @return Length n_groups vector of new sigma values.
// [[Rcpp::export]]
arma::vec sample_sigma_beta_cpp(arma::mat beta_path,
                                arma::ivec group_id,
                                int n_groups,
                                arma::mat Lambda_inv,
                                arma::vec rho_by_group,
                                arma::vec prior_shape,
                                arma::vec prior_scale) {
	const int T = beta_path.n_rows;
	const int p = beta_path.n_cols;
	arma::vec sigma_new(n_groups);
	for (int g = 0; g < n_groups; ++g) sigma_new(g) = 0.25;
	if (T < 2 || p == 0) return sigma_new;
	for (int g = 1; g <= n_groups; ++g) {
		arma::uvec cols = arma::find(group_id == g);
		if (cols.n_elem == 0) continue;
		arma::mat path_g = beta_path.cols(cols);
		arma::mat Linv_g = Lambda_inv.submat(cols, cols);
		double rho_g = rho_by_group(g - 1);
		double S_g = 0.0;
		int q_g = cols.n_elem;
		for (int t = 1; t < T; ++t) {
			arma::vec d = path_g.row(t).t() - rho_g * path_g.row(t-1).t();
			S_g += as_scalar(d.t() * Linv_g * d);
		}
		double a_post = prior_shape(g - 1) + 0.5 * q_g * (T - 1);
		double b_post = prior_scale(g - 1) + 0.5 * S_g;
		// IG draw: variance = b / Gamma(a, 1)
		double gv = R::rgamma(a_post, 1.0);
		double var_new = (gv > 0) ? (b_post / gv) : (b_post);
		sigma_new(g - 1) = std::sqrt(std::max(var_new, 1e-12));
	}
	return sigma_new;
}

// --- 5) EZ rebuild with per-period beta --------------------------------------

//' Compute EZ when beta is time-varying
//'
//' Rebuild the (n x n x T) or (nA x nB x T) EZ cube using a per-period beta
//' vector. Mirrors get_EZ_cpp / get_EZ_bip_cpp but accepts a (T x p) beta
//' matrix where each row is the beta for that period.
//'
//' @param Xlist T-length list of design arrays (each n x n x p or nA x nB x p).
//' @param beta_full_path (T x p) matrix of per-period beta. For coefficients
//'   that are NOT in the dynamic block, the rows are identical (the static
//'   beta replicated across all periods).
//' @param a_mat (n_a x T) row effects (or repmat-ed static effects).
//' @param b_mat (n_b x T) column effects.
//' @param U_cube (n_u x R x T) row latent positions (or replicated-static).
//' @param V_cube (n_v x R x T) column latent positions.
//' @param G (R x R or RA x RB) interaction matrix (bipartite); identity for
//'   unipartite.
//' @param bipartite TRUE/FALSE.
//' @param symmetric TRUE/FALSE.
//' @return n_a x n_b x T cube of EZ values.
// [[Rcpp::export]]
arma::cube get_EZ_dynamic_beta_cpp(List Xlist,
                                   arma::mat beta_full_path,
                                   arma::mat a_mat,
                                   arma::mat b_mat,
                                   arma::cube U_cube,
                                   arma::cube V_cube,
                                   arma::mat G,
                                   bool bipartite,
                                   bool symmetric) {
	const int T = Xlist.size();
	if (T == 0) return arma::cube(0, 0, 0);
	arma::cube X0 = as<arma::cube>(Xlist[0]);
	const int n_a = X0.n_rows;
	const int n_b = X0.n_cols;
	const int p   = X0.n_slices;
	arma::cube EZ(n_a, n_b, T, fill::zeros);
	for (int t = 0; t < T; ++t) {
		arma::cube Xt = as<arma::cube>(Xlist[t]);
		arma::mat Et(n_a, n_b, fill::zeros);
		// X*beta_t
		if (p > 0 && beta_full_path.n_cols >= (arma::uword)p) {
			arma::rowvec b_row = beta_full_path.row(t);
			for (int k = 0; k < p; ++k) {
				double bk = b_row(k);
				if (bk != 0.0) Et += bk * Xt.slice(k);
			}
		}
		// additive effects
		if ((int)a_mat.n_rows == n_a && (int)a_mat.n_cols >= t + 1) {
			Et.each_col() += a_mat.col(t);
		}
		if ((int)b_mat.n_rows == n_b && (int)b_mat.n_cols >= t + 1) {
			Et.each_row() += b_mat.col(t).t();
		}
		// latent factors
		if (U_cube.n_rows == (arma::uword)n_a && U_cube.n_slices >= (arma::uword)(t+1)
		    && V_cube.n_rows == (arma::uword)n_b && V_cube.n_slices >= (arma::uword)(t+1)) {
			arma::mat Ut = U_cube.slice(t);
			arma::mat Vt = V_cube.slice(t);
			if (bipartite) {
				if (G.n_rows == Ut.n_cols && G.n_cols == Vt.n_cols
				    && Ut.n_cols > 0 && Vt.n_cols > 0) {
					Et += Ut * G * Vt.t();
				}
			} else {
				if (Ut.n_cols > 0 && Ut.n_cols == Vt.n_cols) {
					if (symmetric) {
						// L is folded into U for the symmetric path; the caller
						// supplies it that way (V = U for symmetric).
						Et += Ut * Vt.t();
					} else {
						Et += Ut * Vt.t();
					}
				}
			}
		}
		EZ.slice(t) = Et;
	}
	return EZ;
}
