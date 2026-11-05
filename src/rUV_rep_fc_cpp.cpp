//Includes/namespaces
#include <RcppArmadillo.h>
#include "factor_gaussian_stable.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// helper: build index vector excluding one element
static inline arma::uvec indices_excluding(int total, int exclude) {
  arma::uvec idx(total - 1);
  int k = 0;
  for(int i = 0; i < total; i++) {
    if(i != exclude) {
      idx(k++) = i;
    }
  }
  return idx;
}

// [[Rcpp::export]]
arma::mat rwish_cpp(const arma::mat& S0, int nu) {
  const int n = S0.n_rows;
  if (n == 0 || S0.n_cols != S0.n_rows || !S0.is_finite()) {
    Rcpp::stop("Wishart scale must be a finite nonempty square matrix.");
  }
  if (nu <= 0) {
    Rcpp::stop("Wishart degrees of freedom must be positive.");
  }
  // retain every positive scale direction; a ridge changes the prior.
  arma::mat S0_sym = 0.5 * (S0 + S0.t());
  S0_sym = 0.5 * (S0_sym + S0_sym.t());
  arma::mat sS0;
  if (!arma::chol(sS0, S0_sym, "lower")) {
    Rcpp::stop("Wishart scale must be positive definite.");
  }

  arma::mat Z = arma::randn(nu, n) * sS0.t();

  return Z.t() * Z;
}

// Column-wise Gibbs update of the multiplicative factors. Each factor column
// has an exact Gaussian conditional given its partner column (see
// factor_gaussian_stable.h); iSe2 and maxmargin are unused and retained
// for the call signature.
// [[Rcpp::export]]
List rUV_rep_fc_cpp(
    const arma::cube& ET, const arma::mat& U_state, const arma::mat& V_state,
    double rho, double s2, const arma::mat& iSe2,
    double maxmargin, bool shrink, const NumericVector& rLoopIDs,
    double prior_df, double prior_psi) {
  (void)iSe2; (void)maxmargin;
  // the inverse-Wishart prior of the factor covariance: prior_df degrees of freedom and the scale matrix
  // prior_psi * I (R/uv_prior.R)
  if (!std::isfinite(prior_df) || prior_df < 1 || !std::isfinite(prior_psi) || prior_psi <= 0)
    Rcpp::stop("The factor covariance prior needs positive degrees of freedom and scale.");
  const int df0 = static_cast<int>(std::lround(prior_df));
  if (!ET.n_rows || ET.n_rows != ET.n_cols || !ET.n_slices ||
      ET.n_rows != U_state.n_rows || V_state.n_rows != U_state.n_rows ||
      V_state.n_cols != U_state.n_cols || !ET.is_finite() || !U_state.is_finite() || !V_state.is_finite() ||
      !std::isfinite(rho) || std::abs(rho) >= 1 || !std::isfinite(s2) || s2 <= 0)
    Rcpp::stop("Factor conditionals require finite complete residuals, matching coordinates, valid rho and variance.");
  arma::mat U_work = U_state, V_work = V_state;
  const int R = U_work.n_cols, n = U_work.n_rows;
  for (R_xlen_t i = 0; i < rLoopIDs.size(); ++i)
    if (!std::isfinite(rLoopIDs[i]) || rLoopIDs[i] != std::floor(rLoopIDs[i]) || rLoopIDs[i] < 0 || rLoopIDs[i] >= R)
      Rcpp::stop("Factor column schedule must contain valid zero-based indices.");
  if (R == 0) return List::create(Named("U") = U_work, Named("V") = V_work);
  arma::mat UV = join_rows(U_work, V_work);
  arma::mat Suv;
  arma::mat Puv;
  // switch algebra before large factors lose the identity in I + W'W.
  // both routes have the same inverse-Wishart conditional law.
  const bool use_precision = shrink &&
    (R >= 3 || arma::norm(UV, "fro") > 1e4);

  if (use_precision) {
    // psi I + W'W = A'A, evaluated by QR so a large factor direction
    // cannot erase the prior contribution through cancellation.
    arma::mat Q, upper;
    arma::mat augmented = join_cols(std::sqrt(prior_psi) * eye<mat>(2 * R, 2 * R), UV);
    if (!arma::qr_econ(Q, upper, augmented) || !upper.is_finite())
      Rcpp::stop("Factor covariance augmentation QR failed.");
    arma::mat normals = arma::randn(2 * R, n + df0);
    arma::mat root;
    if (!arma::solve(root, arma::trimatu(upper), normals) || !root.is_finite())
      Rcpp::stop("Factor covariance augmentation root solve failed.");
    Puv = root * root.t();
  } else if (shrink) {
    // inverse-Wishart draw of the factor covariance.
    arma::mat UV_prod = UV.t() * UV;
    arma::mat inv_mat = inv_sympd(prior_psi * eye<mat>(R*2, R*2) + UV_prod);
    Suv = inv_sympd(rwish_cpp(inv_mat, n + df0));
  } else {
    // the diffuse covariance of the first half of burn-in, n times the scale of the prior (prior_psi / prior_df: 1 for a
    // unit response, the unit of the factor rows for the normal family)
    Suv = n * (prior_psi / prior_df) * eye<mat>(R*2, R*2);
  }


  if (!use_precision) Puv = lame_factor::covariance_precision(Suv);
  if (!Puv.is_finite() || arma::any(Puv.diag() <= 0))
    Rcpp::stop("Factor conditional prior precision must be finite and positive.");
  const arma::mat E_mean = arma::mean(ET, 2);
  for (R_xlen_t i = 0; i < rLoopIDs.size(); ++i) {
    const int r = rLoopIDs[i], rv = R + r;
    arma::uvec other = indices_excluding(R, r);
    arma::mat Usmall = U_work.cols(other), Vsmall = V_work.cols(other);
    arma::mat other_product = Usmall * Vsmall.t();
    lame_factor::Signal signal = lame_factor::pooled_signal(E_mean, other_product, ET.n_slices, rho, s2);
    arma::uvec free = indices_excluding(2 * R, r), one(1); one(0) = r;
    const double prior_precision = Puv(r, r);
    arma::vec reg = -Puv(free, one) / prior_precision;
    arma::vec prior_mean = join_rows(Usmall, V_work) * reg;
    U_work.col(r) = lame_factor::draw(V_work.col(r), signal, prior_mean, prior_precision, false);
    free = indices_excluding(2 * R, rv); one(0) = rv;
    const double precision_v = Puv(rv, rv);
    reg = -Puv(free, one) / precision_v;
    prior_mean = join_rows(U_work, Vsmall) * reg;
    V_work.col(r) = lame_factor::draw(U_work.col(r), signal, prior_mean, precision_v, true);
  }
  return List::create(Named("U") = U_work, Named("V") = V_work);
}
