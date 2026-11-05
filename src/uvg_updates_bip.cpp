// UV and G updates for bipartite networks
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
using namespace arma;
using namespace Rcpp;

// Preserve the supplied Gaussian precision, without a numerical ridge.
static inline mat inv_spd(const mat& A) {
  mat As = 0.5 * (A + A.t());  // Ensure symmetry
  if (!As.is_finite()) Rcpp::stop("Nonfinite bipartite Gaussian precision.");
  if (arma::any(As.diag() <= 0.0)) Rcpp::stop("Invalid bipartite Gaussian precision.");
  arma::vec scale = 1.0 / arma::sqrt(As.diag());
  arma::mat normalized = As;
  normalized.each_col() %= scale;
  normalized.each_row() %= scale.t();
  // the precision is prior plus data, so positive definite; refuse only
  // when its scaled condition leaves no digit to resolve the inverse
  const double condition = arma::rcond(normalized);
  if (!std::isfinite(condition) || condition <= std::numeric_limits<double>::epsilon())
    Rcpp::stop("Bipartite Gaussian precision is numerically unresolved.");
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, As) || !inverse.is_finite())
    Rcpp::stop("Unable to resolve bipartite Gaussian precision.");
  return inverse;
}

// Sample G (static) given U, V and residual cube R
// Prior: vec(G) ~ N(0, (lambdaG * I)^-1)
// [[Rcpp::export]]
arma::mat sample_G_bip_cpp(const arma::cube& R,
                           const arma::cube& U,
                           const arma::cube& V,
                           const double lambdaG,
                           const arma::vec& s2) {
  const unsigned int T = R.n_slices;
  const unsigned int RA = U.n_cols;
  const unsigned int RB = V.n_cols;

  mat Prec(RA * RB, RA * RB, fill::zeros);
  vec rhs(RA * RB, fill::zeros);

  for (unsigned int t = 0; t < T; ++t) {
    double w = (s2.n_elem == 1) ? 1.0/s2(0) : 1.0/s2(t);
    mat Ut = U.slice(t);  // nA x RA
    mat Vt = V.slice(t);  // nB x RB
    mat At = Ut.t() * Ut;  // RA x RA
    mat Bt = Vt.t() * Vt;  // RB x RB

    // Kronecker product for vec(G) precision
    Prec += kron(Bt, At) * w;

    // Sufficient statistic
    mat St = Ut.t() * R.slice(t) * Vt;  // RA x RB
    rhs += vectorise(St) * w;
  }

  // Add prior precision
  Prec += lambdaG * eye(RA * RB, RA * RB);

  // Sample from posterior
  mat Sigma = inv_spd(Prec);
  vec mean = Sigma * rhs;
  vec draw = mean + chol(Sigma, "lower") * randn<vec>(RA * RB);

  if (!draw.is_finite()) Rcpp::stop("Nonfinite bipartite G Gaussian draw.");

  // Reshape to matrix form
  return reshape(draw, RA, RB);
}
