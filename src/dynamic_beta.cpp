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
#include "ar_truncated_normal.h"
using namespace arma;
using namespace Rcpp;

// --- shared utilities --------------------------------------------------------

// a precision whose diagonally scaled inverse reaches 0.01 / eps leaves under
// two digits to resolve it, and whether its inversion succeeds depends on rounding
static bool resolved_sympd_inverse(const arma::mat& A, const arma::mat& inverse) {
  return inverse.is_finite() &&
    arma::max(inverse.diag() % A.diag()) * 100.0 * arma::datum::eps < 1.0;
}

// Preserve the supplied Gaussian precision. A diagonally scaled solve can
// improve its conditioning; adding a ridge would instead change the target.
static bool safe_sympd_inverse(const arma::mat& M, arma::mat& out) {
  arma::mat A = 0.5 * (M + M.t());
  if (!A.is_finite()) return false;
  if (arma::inv_sympd(out, A) && out.is_finite()) return resolved_sympd_inverse(A, out);
  arma::vec diagonal = A.diag();
  if (arma::any(diagonal <= 0.0)) return false;
  arma::vec scale = 1.0 / arma::sqrt(diagonal);
  arma::mat normalized = A;
  normalized.each_col() %= scale;
  normalized.each_row() %= scale.t();
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, normalized) || !inverse.is_finite()) return false;
  inverse.each_col() %= scale;
  inverse.each_row() %= scale.t();
  out = inverse;
  return resolved_sympd_inverse(A, out);
}

// A Gaussian draw must retain the supplied covariance. Numerical failure
// is explicit rather than replacing negative eigenvalues or returning a mean.
static arma::vec rmvnorm_chol(const arma::vec& mu, const arma::mat& Sigma,
                              int* n_chol_fail = NULL) {
	arma::vec eps(mu.n_elem, fill::randn);
	arma::mat S = 0.5 * (Sigma + Sigma.t());
	arma::mat L;
	if (!mu.is_finite() || !S.is_finite() || !arma::chol(L, S, "lower"))
		Rcpp::stop("Unable to resolve coefficient Gaussian covariance.");
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
	const int p_dyn = Xdyn_t.n_cols;
	// y_D,t = Z_t - offset_t - X_S,t beta_S, as a long vector aligned with
	// the rows of Xdyn_t (which is also a long vector reshape).
	arma::vec y = arma::vectorise(Z_t) - arma::vectorise(offset_t);
	if (Xstat_t.n_cols > 0 && beta_static.n_elem == Xstat_t.n_cols) {
		y -= Xstat_t * beta_static;
	}
	// drop NA observations (treated as missing)
	arma::uvec keep = arma::find_finite(y);
	// self-ties are outside the unipartite likelihood, and z_ii is a
	// finite nuisance draw, so the diagonal is not an observation here
	if (!bipartite && (int) Z_t.n_cols == (int) Z_t.n_rows) {
		const arma::uword n = Z_t.n_rows;
		arma::uvec offd(keep.n_elem);
		arma::uword m0 = 0;
		for (arma::uword q = 0; q < keep.n_elem; ++q) {
			if ((keep[q] % n) != (keep[q] / n)) offd[m0++] = keep[q];
		}
		keep = offd.head(m0);
	}
	// symmetric networks have one latent per UNORDERED dyad: z_ij = z_ji.
	// y is vectorise(Z_t) in column-major order, so entry (i,j) sits at
	// index (j-1)*n + i. Keeping both (i,j) and (j,i) would count each
	// dyad twice, doubling Omega/h and shrinking the coefficient posterior
	// sd by sqrt(2). Restrict to i < j so each dyad contributes once.
	if (symmetric && !bipartite && (int) Z_t.n_cols == (int) Z_t.n_rows) {
		const arma::uword n = Z_t.n_rows;
		arma::uvec upper(keep.n_elem);
		arma::uword m = 0;
		for (arma::uword q = 0; q < keep.n_elem; ++q) {
			arma::uword k = keep[q];
			arma::uword i = k % n, j = k / n;
			if (i < j) upper[m++] = k;
		}
		keep = upper.head(m);
	}
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
			Rcpp::stop("Reciprocal coefficient updates require square outcomes.");
		}
		double det_w = s2 * s2 * (1.0 - dyad_rho * dyad_rho);
		if (!std::isfinite(s2) || s2 <= 0.0 || !std::isfinite(dyad_rho) || std::abs(dyad_rho) >= 1.0)
			Rcpp::stop("Invalid reciprocal coefficient variance.");
		// Preserve the existing arithmetic when its determinant is representable.
		double c_diag = (std::isfinite(det_w) && det_w > 0.0) ? s2 / det_w : 1.0 / s2 / (1.0 - dyad_rho * dyad_rho);
		double c_off = (std::isfinite(det_w) && det_w > 0.0) ? -dyad_rho * s2 / det_w : -dyad_rho * c_diag;
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

// --- 0) per-period sufficient statistics of the dynamic block ---------------

//' Per-period sufficient statistics of the dynamic coefficient block
//'
//' The statistics the forward filter of \code{sample_beta_dynamic_cpp} reads,
//' for the interweaving step of the path hyperparameters: given the latent
//' array, the offsets and the static block, the log density of the dynamic
//' path is, up to a constant, the sum over periods of
//' h_t' beta_t - beta_t' Omega_t beta_t / 2.
//'
//' @inheritParams sample_beta_dynamic_cpp
//' @return List with: Omega -- a (p_dyn x p_dyn x T) array; h -- a (T x p_dyn)
//'   matrix.
//' @noRd
// [[Rcpp::export]]
List beta_dynamic_stats_cpp(List Xdyn_list, List Xstat_list,
                            List Z_list, List offset_list,
                            arma::vec beta_static,
                            double s2, double dyad_rho,
                            bool bipartite, bool symmetric,
                            bool use_dyad_rho) {
	const int T = Xdyn_list.size();
	if (T == 0) {
		return List::create(Named("Omega") = arma::cube(0, 0, 0),
		                    Named("h") = arma::mat(0, 0));
	}
	const int p_dyn = as<arma::mat>(Xdyn_list[0]).n_cols;
	arma::cube Omega(p_dyn, p_dyn, T, fill::zeros);
	arma::mat h(T, p_dyn, fill::zeros);
	for (int t = 0; t < T; ++t) {
		arma::mat Xdyn_t = as<arma::mat>(Xdyn_list[t]);
		arma::mat Xstat_t;
		if (Xstat_list.size() > t) {
			SEXP sx = Xstat_list[t];
			if (!Rf_isNull(sx)) Xstat_t = as<arma::mat>(sx);
		}
		arma::mat Z_t = as<arma::mat>(Z_list[t]);
		arma::mat offset_t = as<arma::mat>(offset_list[t]);
		arma::mat Omega_t(p_dyn, p_dyn, fill::zeros);
		arma::vec h_t(p_dyn, fill::zeros);
		accumulate_Omega_h(Xdyn_t, Xstat_t, Z_t, offset_t, beta_static,
		                   s2, dyad_rho, bipartite, symmetric, use_dyad_rho,
		                   Omega_t, h_t);
		Omega.slice(t) = Omega_t;
		h.row(t) = h_t.t();
	}
	return List::create(Named("Omega") = Omega, Named("h") = h);
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
//' @param beta0_mean Length p_dyn prior mean for the first state beta_1
//'   (the prior is placed directly on beta_1, no predict step at t = 1).
//' @param beta0_cov  p_dyn x p_dyn prior covariance for beta_1. Must be a
//'   fixed matrix that does not depend on the current (rho, sigma) draw so
//'   the transition-only hyperparameter conditionals remain exact.
//' @param s2 Dyadic variance.
//' @param dyad_rho Dyadic correlation (ignored when use_dyad_rho=FALSE).
//' @param bipartite Whether the network is bipartite.
//' @param symmetric Whether the network is symmetric.
//' @param use_dyad_rho Whether to use the dyad-corr branch (TRUE only for
//'   unipartite, asymmetric, with a non-zero rho).
//' @param Omega_in,h_in Optional per-period sufficient statistics, as
//'   returned by \code{beta_dynamic_stats_cpp} for the same arguments (a
//'   p_dyn x p_dyn x T array and a T x p_dyn matrix); the filter then reads
//'   them instead of accumulating them again. The default (NULL) accumulates.
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
                             bool use_dyad_rho,
                             Rcpp::Nullable<arma::cube> Omega_in = R_NilValue,
                             Rcpp::Nullable<arma::mat> h_in = R_NilValue) {
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

	// the per-period statistics, when the caller computed them already
	const bool have_stats = Omega_in.isNotNull() && h_in.isNotNull();
	arma::cube Omega_given;
	arma::mat h_given;
	if (have_stats) {
		Omega_given = Rcpp::as<arma::cube>(Omega_in);
		h_given = Rcpp::as<arma::mat>(h_in);
		if ((int) Omega_given.n_slices != T || (int) Omega_given.n_rows != p_dyn ||
		    (int) Omega_given.n_cols != p_dyn || (int) h_given.n_rows != T ||
		    (int) h_given.n_cols != p_dyn) {
			Rcpp::stop("Sufficient statistics do not match the coefficient design.");
		}
	}

	// initial-condition coherence: the prior is placed DIRECTLY on beta_1
	// (m_pred[0] = beta0_mean, C_pred[0] = beta0_cov, no F/Q predict step),
	// and beta0_cov must be a fixed, (rho, sigma)-independent matrix -- the
	// caller freezes any stationary inflation once at initialization from
	// the prior means. This makes p(beta_1) parameter-free so the
	// transition-only rho/sigma conditionals (sample_rho_beta_cpp /
	// sample_sigma_beta_cpp) are exact full conditionals. Propagating a
	// (rho, sigma)-dependent beta0_cov_eff through F/Q would make the
	// hyperparameter updates incoherent with the state draw.

	for (int t = 0; t < T; ++t) {
		arma::vec m_p;
		arma::mat C_p;
		if (t == 0) {
			// fixed prior directly on the first state
			m_p = beta0_mean;
			C_p = beta0_cov;
		} else {
			// predict: m_t|t-1 = F m_prev, C_t|t-1 = F C_prev F' + Q
			m_p = F * m_filt[t-1];
			C_p = F * C_filt[t-1] * F.t() + Q;
		}
		C_p = 0.5 * (C_p + C_p.t());
		m_pred[t] = m_p;
		C_pred[t] = C_p;
		// accumulate Omega_t and h_t (or read the statistics the caller passed)
		arma::mat Omega_t(p_dyn, p_dyn, fill::zeros);
		arma::vec h_t(p_dyn, fill::zeros);
		if (have_stats) {
			Omega_t = Omega_given.slice(t);
			h_t = h_given.row(t).t();
		} else {
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
		}
		// posterior precision K = inv(C_pred) + Omega, then C_filt = inv(K)
		arma::mat C_p_inv;
		if (!safe_sympd_inverse(C_p, C_p_inv)) {
			Rcpp::stop("Unable to resolve coefficient predictive covariance.");
		}
		arma::mat K = C_p_inv + Omega_t;
		arma::mat C_f;
		if (!safe_sympd_inverse(K, C_f)) {
			Rcpp::stop("Unable to resolve coefficient posterior precision.");
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
			Rcpp::stop("Unable to resolve coefficient smoothing covariance.");
		}
		arma::mat J_t = C_filt[t] * F.t() * C_p_next_inv;
		arma::vec m_smooth = m_filt[t] + J_t * (beta_t - m_pred[t+1]);
		arma::mat C_smooth = C_filt[t] - J_t * F * C_filt[t];
		// If subtraction loses a small innovation variance, use the equivalent
		// sum-of-squares covariance; it does not subtract nearly equal matrices.
		if (arma::any(C_smooth.diag() < 1e-6 * C_filt[t].diag())) {
			arma::mat residual = arma::eye<arma::mat>(p_dyn, p_dyn) - J_t * F;
			C_smooth = residual * C_filt[t] * residual.t() + J_t * Q * J_t.t();
		}
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
		Rcpp::stop("Unable to resolve static coefficient posterior precision.");
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
		// Retain ordinary rejection draws and their RNG sequence. If the
		// truncation makes rejection inefficient, use an exact log-tail
		// inverse CDF; neither probabilities nor draws are clamped.
		double prop = post_m;
		double sd_post = std::sqrt(post_v);
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
			prop = lame_ar_truncnorm_draw(post_m, sd_post, rho_lower, rho_upper);
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
		if (!(gv > 0.0) || !std::isfinite(gv))
			Rcpp::stop("Invalid dynamic coefficient innovation precision draw.");
		double var_new = b_post / gv;
		if (!(var_new > 0.0) || !std::isfinite(var_new))
			Rcpp::stop("Dynamic coefficient innovation variance is outside numerical range.");
		sigma_new(g - 1) = std::sqrt(var_new);
	}
	return sigma_new;
}
