// Bipartite dynamic latent factor updates
// Replaces the nested R loops in lame.R for bipartite dynamic_uv
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Fast 2x2 symmetric positive definite inverse (direct formula)
// For R=2, this avoids LAPACK overhead entirely
static inline arma::mat inv_sympd_small(const arma::mat& A) {
  if(A.n_rows == 1) {
    arma::mat result(1, 1);
    result(0, 0) = 1.0 / A(0, 0);
    return result;
  }
  if(A.n_rows == 2) {
    double a = A(0,0), b = A(0,1), d = A(1,1);
    double det = a * d - b * b;
    if(std::abs(det) < 1e-15) {
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
  if(A.n_rows == 1) {
    arma::mat L(1, 1);
    L(0, 0) = sqrt(A(0, 0));
    return L;
  }
  if(A.n_rows == 2) {
    arma::mat L(2, 2, fill::zeros);
    L(0,0) = sqrt(A(0,0));
    if(L(0,0) > 1e-15) {
      L(1,0) = A(1,0) / L(0,0);
      double tmp = A(1,1) - L(1,0) * L(1,0);
      L(1,1) = sqrt(std::max(tmp, 1e-15));
    } else {
      L(1,1) = sqrt(A(1,1));
    }
    return L;
  }
  arma::mat L;
  arma::chol(L, A, "lower");
  return L;
}

//' Bipartite dynamic UV Gibbs update
//'
//' Replaces the nested R loops for bipartite dynamic_uv in lame.R.
//' Uses direct 2x2 inverse formula for common R=2 case.
//'
//' @param U_cube 3D array (nA x RA x T)
//' @param V_cube 3D array (nB x RB x T)
//' @param E 3D array of residuals (nA x nB x T)
//' @param G Interaction matrix (RA x RB)
//' @param rho_uv AR(1) persistence parameter
//' @param sigma_uv Innovation standard deviation
//' @param s2 Dyadic variance
//' @return List with updated U_cube, V_cube
// [[Rcpp::export]]
List rUV_dynamic_bip_fc_cpp(arma::cube U_cube, arma::cube V_cube,
                             const arma::cube& E, const arma::mat& G,
                             double rho_uv, double sigma_uv, double s2) {
  const int nA = U_cube.n_rows;
  const int RA = U_cube.n_cols;
  const int nB = V_cube.n_rows;
  const int RB = V_cube.n_cols;
  const int T = U_cube.n_slices;

  // ar(1) prior terms need s2 scaling to match the s2-parameterized
  // posterior: prec = W'W + s2*(ar terms), var = s2 * inv(prec)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_s2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_s2 = s2 * rho_uv * rho_uv * sigma2_inv;

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
      }
      if(t < T-1) {
        prec.diag() += s2_rho2_s2;
        ei += s2_rho_s2 * U_cube.slice(t+1).row(i).t();
      }
      prec.diag() += 1.0;

      arma::mat iprec = inv_sympd_small(prec);
      arma::vec mu_i = iprec * ei;
      arma::mat ch = chol_lower_small(s2 * iprec);
      arma::vec z(RA, fill::randn);
      U_cube.slice(t).row(i) = (mu_i + ch * z).t();
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
      }
      if(t < T-1) {
        prec.diag() += s2_rho2_s2;
        ej += s2_rho_s2 * V_cube.slice(t+1).row(j).t();
      }
      prec.diag() += 1.0;

      arma::mat iprec = inv_sympd_small(prec);
      arma::vec mu_j = iprec * ej;
      arma::mat ch = chol_lower_small(s2 * iprec);
      arma::vec z(RB, fill::randn);
      V_cube.slice(t).row(j) = (mu_j + ch * z).t();
    }
  }

  return List::create(Named("U") = U_cube, Named("V") = V_cube);
}
