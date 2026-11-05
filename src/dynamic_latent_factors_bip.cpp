// Bipartite dynamic latent factor updates
// Replaces the nested R loops in lame.R for bipartite dynamic_uv
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

static inline bool vec_is_finite(const arma::vec& x) {
  return x.is_finite();
}

// Fast 2x2 symmetric positive definite inverse (direct formula)
// For R=2, this avoids LAPACK overhead entirely
static inline arma::mat inv_sympd_small(const arma::mat& A) {
  if (!A.is_finite() || arma::any(A.diag() <= 0.0)) Rcpp::stop("Invalid bipartite Gaussian matrix.");
  if(A.n_rows == 1) {
    arma::mat result(1, 1);
    result(0, 0) = 1.0 / A(0, 0);
    return result;
  }
  if(A.n_rows == 2) {
    double a = A(0,0), b = A(0,1), d = A(1,1);
    double det = a * d - b * b;
    if(!std::isfinite(det) || det <= 0.0 || std::abs(det) < 1e-15) {
      return arma::inv_sympd(A);  // fallback
    }
    arma::mat result(2, 2);
    double inv_det = 1.0 / det;
    result(0,0) = d * inv_det;
    result(0,1) = -b * inv_det;
    result(1,0) = -b * inv_det;
    result(1,1) = a * inv_det;
    return result;
  }
  // General case
  return arma::inv_sympd(A);
}

// Fast lower Cholesky for small matrices
static inline arma::mat chol_lower_small(const arma::mat& A) {
  if (!A.is_finite() || arma::any(A.diag() <= 0.0)) Rcpp::stop("Invalid bipartite Gaussian matrix.");
  if(A.n_rows == 1) {
    arma::mat L(1, 1);
    L(0, 0) = sqrt(A(0, 0));
    return L;
  }
  if(A.n_rows == 2) {
    arma::mat L(2, 2, fill::zeros);
    L(0,0) = sqrt(A(0,0));
    if(std::isfinite(L(0,0)) && L(0,0) > 0.0) {
      L(1,0) = A(1,0) / L(0,0);
      double tmp = A(1,1) - L(1,0) * L(1,0);
      if (!std::isfinite(tmp) || tmp <= 0.0) Rcpp::stop("Unable to resolve bipartite Gaussian covariance.");
      L(1,1) = sqrt(tmp);
    } else {
      Rcpp::stop("Invalid bipartite Gaussian covariance.");
    }
    return L;
  }
  arma::mat L;
  if (!arma::chol(L, A, "lower")) Rcpp::stop("Unable to resolve bipartite Gaussian covariance.");
  return L;
}

//' Bipartite dynamic UV Gibbs update
//'
//' Replaces the nested R loops for bipartite dynamic_uv in lame.R.
//' Uses direct 2x2 inverse formula for common R=2 case.
//'
//' @param U_state 3D array (nA x RA x T)
//' @param V_state 3D array (nB x RB x T)
//' @param E 3D array of residuals (nA x nB x T)
//' @param G Interaction matrix (RA x RB)
//' @param rho_uv AR(1) persistence parameter
//' @param sigma_uv Innovation standard deviation
//' @param s2 Dyadic variance
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite
//'   values disable it. Each actor conditional respects the bound.
//' @return List with updated U_cube, V_cube
// [[Rcpp::export]]
List rUV_dynamic_bip_fc_cpp(const arma::cube& U_state, const arma::cube& V_state,
                             const arma::cube& E, const arma::mat& G,
                             double rho_uv, double sigma_uv, double s2, double max_abs = -1.0) {
  // rcpparmadillo wraps the caller's array instead of copying it, so the
  // kernel draws into copies of its own: a sweep the caller refuses has to
  // leave the positions it was given exactly as they were
  arma::cube U_cube = U_state;
  arma::cube V_cube = V_state;
  const int nA = U_cube.n_rows;
  const int RA = U_cube.n_cols;
  const int nB = V_cube.n_rows;
  const int RB = V_cube.n_cols;
  const int T = U_cube.n_slices;
  const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
  int bound_rejections = 0;

  if(!U_cube.is_finite() || !V_cube.is_finite() || !G.is_finite() || !std::isfinite(rho_uv) || std::abs(rho_uv) >= 1.0 ||
     !std::isfinite(sigma_uv) || sigma_uv <= 0 ||
     !std::isfinite(s2) || s2 <= 0) {
    Rcpp::stop("Invalid bipartite Gaussian state or scale.");
  }

  // ar(1) prior terms need s2 scaling to match the s2-parameterized
  // posterior: prec = W'W + s2*(ar terms), var = s2 * inv(prec)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_s2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_s2 = s2 * rho_uv * rho_uv * sigma2_inv;
  // stationary AR(1) initial condition at t = 0 (s2-scaled prior precision)
  const double s2_stat_init = s2 * (1.0 - rho_uv * rho_uv) * sigma2_inv;

  for(int t = 0; t < T; t++) {
    // Get slices
    arma::mat E_t = E.slice(t);
    // Replace NAs with 0
    E_t.replace(datum::nan, 0.0);

    arma::mat& V_t = V_cube.slice(t);

    // U update: W_t = V_t * G^T (nB x RA)
    arma::mat W_t = V_t * G.t();
    arma::mat WtW = W_t.t() * W_t;  // RA x RA

    // Batch compute: E_all_U = E_t * W_t gives nA x RA
    arma::mat E_all_U = E_t * W_t;

    for(int i = 0; i < nA; i++) {
      arma::vec ei = E_all_U.row(i).t();
      arma::mat prec = WtW;

      if(t > 0) {
        prec.diag() += s2_sigma2_inv;
        ei += s2_rho_s2 * U_cube.slice(t-1).row(i).t();
      } else {
        // stationary AR(1) prior at t = 0 (the model's marginal init)
        prec.diag() += s2_stat_init;
      }
      if(t < T-1) {
        prec.diag() += s2_rho2_s2;
        ei += s2_rho_s2 * U_cube.slice(t+1).row(i).t();
      }
      // The AR prior supplies the precision; an extra ridge changes its law.

      try {
        arma::mat iprec = inv_sympd_small(prec);
        arma::vec mu_i = iprec * ei;
        arma::mat ch = chol_lower_small(s2 * iprec);
        arma::vec z(RA, fill::randn);
        arma::vec draw = mu_i + ch * z;
        if(vec_is_finite(draw)) {
          if (!bounded || arma::abs(draw).max() <= max_abs) {
            U_cube.slice(t).row(i) = draw.t();
          } else {
            ++bound_rejections;
          }
        } else {
          Rcpp::stop("Non-finite bipartite U Gaussian draw.");
        }
      } catch(...) {
        Rcpp::stop("Unable to resolve bipartite U Gaussian conditional.");
      }
    }

    // V update: Q_t = U_t * G (nA x RB)
    arma::mat& U_t = U_cube.slice(t);
    arma::mat Q_t = U_t * G;
    arma::mat QtQ = Q_t.t() * Q_t;  // RB x RB

    // Batch compute: E_all_V = E_t^T * Q_t gives nB x RB
    arma::mat E_all_V = E_t.t() * Q_t;

    for(int j = 0; j < nB; j++) {
      arma::vec ej = E_all_V.row(j).t();
      arma::mat prec = QtQ;

      if(t > 0) {
        prec.diag() += s2_sigma2_inv;
        ej += s2_rho_s2 * V_cube.slice(t-1).row(j).t();
      } else {
        prec.diag() += s2_stat_init;
      }
      if(t < T-1) {
        prec.diag() += s2_rho2_s2;
        ej += s2_rho_s2 * V_cube.slice(t+1).row(j).t();
      }
      // The AR prior supplies the precision; an extra ridge changes its law.

      try {
        arma::mat iprec = inv_sympd_small(prec);
        arma::vec mu_j = iprec * ej;
        arma::mat ch = chol_lower_small(s2 * iprec);
        arma::vec z(RB, fill::randn);
        arma::vec draw = mu_j + ch * z;
        if(vec_is_finite(draw)) {
          if (!bounded || arma::abs(draw).max() <= max_abs) {
            V_cube.slice(t).row(j) = draw.t();
          } else {
            ++bound_rejections;
          }
        } else {
          Rcpp::stop("Non-finite bipartite V Gaussian draw.");
        }
      } catch(...) {
        Rcpp::stop("Unable to resolve bipartite V Gaussian conditional.");
      }
    }
  }

  // no in-kernel U/V column scale balancing: a state-dependent scale map
  // does not preserve the AR(1) prior's law, so any rebalancing belongs
  // on stored draws only.

  if (bounded) return List::create(Named("U") = U_cube, Named("V") = V_cube,
    Named("bound_rejections") = bound_rejections);
  return List::create(Named("U") = U_cube, Named("V") = V_cube);
}
