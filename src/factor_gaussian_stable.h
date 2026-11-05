#ifndef LAME_FACTOR_GAUSSIAN_STABLE_H
#define LAME_FACTOR_GAUSSIAN_STABLE_H
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
namespace lame_factor {
struct Signal {
  arma::mat symmetric;
  arma::mat antisymmetric;
  double plus;
  double minus;
};
struct System {
  arma::vec mean;
  arma::vec reflector;
  arma::vec inverse_roots;
};
inline arma::vec reflect(const arma::vec& x, const arma::vec& h) {
  if (h.n_elem == 0) return x;
  return x - 2.0 * h * arma::dot(h, x);
}
// Precision of the stacked factor prior covariance. The covariance is
// symmetrized and equilibrated before its Cholesky factor; a covariance
// that is not numerically positive definite is inverted through its
// eigendecomposition with a relative eigenvalue floor.
inline arma::mat covariance_precision(const arma::mat& covariance) {
  const arma::uword d = covariance.n_rows;
  if (!d || covariance.n_cols != d || !covariance.is_finite() || arma::any(covariance.diag() <= 0))
    Rcpp::stop("Factor prior covariance must be finite with positive diagonal.");
  arma::vec scales = arma::sqrt(covariance.diag());
  arma::mat outer = scales * scales.t();
  arma::mat balanced = covariance / outer;
  balanced = 0.5 * (balanced + balanced.t());
  arma::mat precision, root;
  if (arma::chol(root, balanced)) {
    arma::mat inverse = arma::solve(arma::trimatu(root), arma::eye<arma::mat>(d, d));
    precision = inverse * inverse.t();
  } else {
    arma::vec ev; arma::mat V;
    if (!arma::eig_sym(ev, V, balanced)) Rcpp::stop("Factor prior covariance eigendecomposition failed.");
    const double floor = std::max(ev.max(), 1.0) * std::sqrt(std::numeric_limits<double>::epsilon());
    ev.elem(arma::find(ev < floor)).fill(floor);
    precision = V * arma::diagmat(1 / ev) * V.t();
  }
  precision = 0.5 * (precision + precision.t()) / outer;
  if (!precision.is_finite()) Rcpp::stop("Factor prior precision is not finite.");
  return precision;
}
// Same signal from the period average of the residuals, which carries all
// the information when the reciprocal weights are equal across periods.
inline Signal pooled_signal(const arma::mat& E_mean, const arma::mat& other, arma::uword nt,
                            double rho, double s2) {
  Signal signal;
  signal.plus = nt * ((1.0 / s2) / (1.0 + rho));
  signal.minus = nt * ((1.0 / s2) / (1.0 - rho));
  if (!std::isfinite(signal.plus) || !std::isfinite(signal.minus) || signal.plus <= 0 || signal.minus <= 0)
    Rcpp::stop("Factor Gaussian likelihood weights cannot be represented.");
  arma::mat partial = E_mean - other;
  arma::mat partial_t = partial.t();
  signal.symmetric = 0.5 * (partial + partial_t);
  signal.antisymmetric = 0.5 * (partial - partial_t);
  return signal;
}
inline System conditional_system(const arma::vec& partner, const Signal& signal,
                                 const arma::vec& prior_mean, double prior_precision,
                                 bool transpose = false) {
  const arma::uword n = partner.n_elem;
  if (!n || prior_mean.n_elem != n || !partner.is_finite() || !prior_mean.is_finite() ||
      !std::isfinite(prior_precision) || prior_precision <= 0 ||
      signal.symmetric.n_rows != n || signal.symmetric.n_cols != n ||
      signal.antisymmetric.n_rows != n || signal.antisymmetric.n_cols != n ||
      !signal.symmetric.is_finite() || !signal.antisymmetric.is_finite())
    Rcpp::stop("Factor Gaussian signal, partner and conditional prior must be valid.");
  System system;
  const double scale = arma::abs(partner).max();
  if (scale == 0) {
    system.mean = prior_mean; system.reflector.reset();
    system.inverse_roots = arma::vec(n, arma::fill::value(1.0 / std::sqrt(prior_precision)));
    return system;
  }
  arma::vec scaled = partner / scale;
  const double norm_scaled = std::sqrt(arma::dot(scaled, scaled));
  arma::vec direction = scaled / norm_scaled;
  system.reflector = direction;
  system.reflector(0) += direction(0) >= 0 ? 1.0 : -1.0;
  system.reflector /= std::sqrt(arma::dot(system.reflector, system.reflector));
  const double plus = signal.plus, minus = signal.minus, self = plus / 2.0 + minus / 2.0;
  const double wp = (std::sqrt(plus) * scale) * norm_scaled;
  const double wd = (std::sqrt(self) * scale) * norm_scaled;
  const double parallel = prior_precision + wp * wp, perpendicular = prior_precision + wd * wd;
  if (!std::isfinite(parallel) || !std::isfinite(perpendicular) || parallel <= 0 || perpendicular <= 0)
    Rcpp::stop("Factor Gaussian eigenprecision cannot be represented.");
  arma::vec prior = reflect(prior_mean, system.reflector);
  arma::vec symmetric = reflect(signal.symmetric * direction, system.reflector);
  arma::vec antisymmetric = reflect(signal.antisymmetric * direction, system.reflector);
  if (transpose) antisymmetric = -antisymmetric;
  const auto coefficient = [scale, norm_scaled](double weight, double precision) {
    return (((std::sqrt(weight) * scale) * norm_scaled) / std::sqrt(precision)) *
      (std::sqrt(weight) / std::sqrt(precision));
  };
  arma::vec coordinates(n);
  coordinates(0) = coefficient(plus, parallel) * symmetric(0) +
    (prior_precision / parallel) * prior(0);
  for (arma::uword i = 1; i < n; ++i)
    coordinates(i) = coefficient(plus, perpendicular) * symmetric(i) +
      coefficient(minus, perpendicular) * antisymmetric(i) + (prior_precision / perpendicular) * prior(i);
  system.mean = reflect(coordinates, system.reflector);
  system.inverse_roots = arma::vec(n, arma::fill::value(1.0 / std::sqrt(perpendicular)));
  system.inverse_roots(0) = 1.0 / std::sqrt(parallel);
  if (!system.mean.is_finite() || !system.inverse_roots.is_finite() || arma::any(system.inverse_roots <= 0))
    Rcpp::stop("Factor Gaussian moments cannot be represented.");
  return system;
}
inline arma::vec draw(const arma::vec& partner, const Signal& signal,
                      const arma::vec& prior_mean, double prior_precision, bool transpose = false) {
  System system = conditional_system(partner, signal, prior_mean, prior_precision, transpose);
  arma::vec e = arma::randn<arma::vec>(partner.n_elem);
  arma::vec noise = reflect(system.inverse_roots % e, system.reflector);
  arma::vec value = system.mean + noise;
  if (!value.is_finite()) Rcpp::stop("Factor Gaussian draw is not finite.");
  return value;
}
}
#endif
