#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// evaluate the Poisson log mass without clipping or cancellation near its mode.
// the log mass at the mode, log p(y | y), is R::dpois at an integer y and
// y log(y) - y - lgamma(y + 1) otherwise: R::dpois is -Inf at a non-integer
// y, so every proposal was rejected and Z never moved (the R entry points
// refuse non-integer counts; see .lame_check_poisson_counts). a negative y
// has no mass.
static double poisson_log_mass(double y, double log_rate) {
  if (std::isnan(log_rate)) return NA_REAL;
  if (y < 0) return R_NegInf;
  if (!std::isfinite(log_rate)) {
    return log_rate < 0 && y == 0 ? 0.0 : R_NegInf;
  }
  if (y == 0) return -std::exp(log_rate);
  double delta = log_rate - std::log(y);
  double at_mode = (y == std::floor(y)) ? R::dpois(y, y, true) :
    y * std::log(y) - y - std::lgamma(y + 1.0);
  return at_mode - y * (std::expm1(delta) - delta);
}

// Helper function to compute log density for GBME models
// [[Rcpp::export]]
arma::mat ldZgbme_cpp(const arma::mat& Z, const arma::mat& Y,
                      const arma::mat& EZ, double rho, double s2) {
  // From log p(z|ez,rho,s2)
  double c = 0.5 * (1.0/std::sqrt(1+rho) + 1.0/std::sqrt(1-rho)) / std::sqrt(s2);
  double d = 0.5 * (1.0/std::sqrt(1+rho) - 1.0/std::sqrt(1-rho)) / std::sqrt(s2);

  arma::mat E = Z - EZ;
  arma::mat lpZ = -0.5 * arma::square(c * E + d * E.t());
  lpZ = lpZ + lpZ.t();
  lpZ.diag() = lpZ.diag() / 2.0;

  // From log p(y|z) for Poisson model
  arma::mat llZ(Y.n_rows, Y.n_cols);
  for(unsigned int i = 0; i < Y.n_rows; i++) {
    for(unsigned int j = 0; j < Y.n_cols; j++) {
      if(!std::isnan(Y(i,j))) {
        llZ(i,j) = poisson_log_mass(Y(i,j), Z(i,j));
      } else {
        llZ(i,j) = 0;
      }
    }
  }

  // a missing cell contributes nothing to either member of its dyad, so
  // the mask runs before the two cells are combined
  arma::uvec na_indices = arma::find_nonfinite(Y);
  llZ.elem(na_indices).zeros();
  llZ = llZ + llZ.t();

  return lpZ + llZ;
}

// [[Rcpp::export]]
arma::mat rZ_pois_fc_cpp(const arma::mat& Z, const arma::mat& EZ,
                         double rho, double s2, const arma::mat& Y) {
  // this kernel assumes a square (unipartite) network; the rectangular
  // bipartite panel is served by rZ_pois_bip_fc() in R/rZ_bipartite.R
  if(Z.n_rows != Z.n_cols) {
    Rcpp::stop("rZ_pois_fc_cpp requires square matrices (unipartite). "
               "Got %d x %d.", Z.n_rows, Z.n_cols);
  }
  int n = Y.n_rows;

  // Copy Z for updating
  arma::mat Z_new = Z;

  // Propose candidate Z
  arma::mat Zp = Z + arma::randn(n, n) * std::sqrt(s2);

  // Compute log densities for current and proposed values
  arma::mat lr_current = ldZgbme_cpp(Z, Y, EZ, rho, s2);
  arma::mat lr_proposed = ldZgbme_cpp(Zp, Y, EZ, rho, s2);

  // Compute acceptance ratio
  arma::mat lr = lr_proposed - lr_current;

  // Generate symmetric matrix of log uniform random variables
  arma::mat lh = arma::log(arma::randu(n, n));
  lh = arma::trimatl(lh) + arma::trimatl(lh, -1).t();

  // Update dyads where lr > lh
  arma::umat accept = (lr > lh);
  Z_new.elem(arma::find(accept)) = Zp.elem(arma::find(accept));

  return Z_new;
}

// -- poisson with period exposure ---------------------------------------
// Z is the unexposed latent log-rate; exposure enters only the Poisson
// observation likelihood via lambda = exposure * exp(Z). the unscaled
// rZ_pois_fc_cpp / ldZgbme_cpp / simY_pois are used when period_exposure
// is NULL or all-ones; the exposure-aware variants below are dispatched
// from R when any exposure differs from 1.

// log-density helper that adds the log-exposure offset to the Poisson
// data term. Mirrors ldZgbme_cpp exactly when log_exposure == 0.
// [[Rcpp::export]]
arma::mat ldZgbme_pois_exposure_cpp(const arma::mat& Z, const arma::mat& Y,
                                     const arma::mat& EZ, double rho,
                                     double s2, double log_exposure) {
  double c = 0.5 * (1.0/std::sqrt(1+rho) + 1.0/std::sqrt(1-rho)) / std::sqrt(s2);
  double d = 0.5 * (1.0/std::sqrt(1+rho) - 1.0/std::sqrt(1-rho)) / std::sqrt(s2);

  arma::mat E = Z - EZ;
  arma::mat lpZ = -0.5 * arma::square(c * E + d * E.t());
  lpZ = lpZ + lpZ.t();
  lpZ.diag() = lpZ.diag() / 2.0;

  // Poisson log-likelihood with exposure offset: lambda = exp(log_exp + Z)
  // -> ll = y * (log_exp + Z) - exp(log_exp + Z) - lgamma(y+1)
  // The y * log_exp term drops out of the MH ratio because it does not
  // depend on Z (the proposal moves Z, not log_exp).
  arma::mat llZ(Y.n_rows, Y.n_cols);
  for(unsigned int i = 0; i < Y.n_rows; i++) {
    for(unsigned int j = 0; j < Y.n_cols; j++) {
      if(!std::isnan(Y(i,j))) {
        llZ(i,j) = poisson_log_mass(Y(i,j), log_exposure + Z(i,j));
      } else {
        llZ(i,j) = 0;
      }
    }
  }
  // a missing cell contributes nothing to either member of its dyad, so
  // the mask runs before the two cells are combined
  arma::uvec na_indices = arma::find_nonfinite(Y);
  llZ.elem(na_indices).zeros();
  llZ = llZ + llZ.t();

  return lpZ + llZ;
}

// MH update for Z with period exposure threaded into the data term.
// [[Rcpp::export]]
arma::mat rZ_pois_fc_exposure_cpp(const arma::mat& Z, const arma::mat& EZ,
                                   double rho, double s2,
                                   const arma::mat& Y, double log_exposure) {
  int n = Y.n_rows;
  arma::mat Z_new = Z;
  arma::mat Zp = Z + arma::randn(n, n) * std::sqrt(s2);
  arma::mat lr_current  = ldZgbme_pois_exposure_cpp(Z,  Y, EZ, rho, s2, log_exposure);
  arma::mat lr_proposed = ldZgbme_pois_exposure_cpp(Zp, Y, EZ, rho, s2, log_exposure);
  arma::mat lr = lr_proposed - lr_current;
  arma::mat lh = arma::log(arma::randu(n, n));
  lh = arma::trimatl(lh) + arma::trimatl(lh, -1).t();
  arma::umat accept = (lr > lh);
  Z_new.elem(arma::find(accept)) = Zp.elem(arma::find(accept));
  return Z_new;
}

// Posterior-predictive simulate variant honoring period exposure.
// [[Rcpp::export]]
arma::mat simY_pois_exposure(const arma::mat& EZ, double exposure, bool bipartite = false) {
  if (!std::isfinite(exposure) || exposure < 0) {
    Rcpp::stop("Poisson exposure must be finite and nonnegative.");
  }
  int n_a = EZ.n_rows;
  int n_b = EZ.n_cols;
  bool square = (n_a == n_b);
  arma::mat Y(n_a, n_b);
  for (int i = 0; i < n_a; i++) {
    for (int j = 0; j < n_b; j++) {
      if ((square && !bipartite && i == j) || std::isnan(EZ(i,j))) {
        Y(i,j) = NA_REAL;
      } else {
        double lambda = exposure == 0 ? 0.0 : std::exp(EZ(i,j) + std::log(exposure));
        if (!std::isfinite(lambda)) {
          Rcpp::stop("The simulated Poisson rate overflowed; no clipped count draw was substituted.");
        }
        Y(i,j) = R::rpois(lambda);
      }
    }
  }
  return Y;
}

// draw from the unexposed log rate; square bipartite diagonals are observations.
// [[Rcpp::export]]
arma::mat simY_pois(const arma::mat& EZ, bool bipartite = false) {
  return simY_pois_exposure(EZ, 1.0, bipartite);
}
