#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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
        // Poisson log likelihood: y*log(lambda) - lambda - log(y!)
        double lambda = std::exp(Z(i,j));
        if(lambda > 0) {
          llZ(i,j) = Y(i,j) * Z(i,j) - lambda - std::lgamma(Y(i,j) + 1.0);
        } else {
          llZ(i,j) = -1e10; // Very small log likelihood for numerical stability
        }
      } else {
        llZ(i,j) = 0;
      }
    }
  }
  
  llZ = llZ + llZ.t();
  
  // Set NA positions to 0
  arma::uvec na_indices = arma::find_nonfinite(Y);
  llZ.elem(na_indices).zeros();
  
  return lpZ + llZ;
}

// [[Rcpp::export]]
arma::mat rZ_pois_fc_cpp(const arma::mat& Z, const arma::mat& EZ, 
                         double rho, double s2, const arma::mat& Y) {
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

// Simplified version for normal family (for efficiency)
// [[Rcpp::export]]
arma::mat ldZgbme_nrm_cpp(const arma::mat& Z, const arma::mat& Y, 
                          const arma::mat& EZ, double rho, double s2) {
  // From log p(z|ez,rho,s2)
  double c = 0.5 * (1.0/std::sqrt(1+rho) + 1.0/std::sqrt(1-rho)) / std::sqrt(s2);
  double d = 0.5 * (1.0/std::sqrt(1+rho) - 1.0/std::sqrt(1-rho)) / std::sqrt(s2);
  
  arma::mat E = Z - EZ;
  arma::mat lpZ = -0.5 * arma::square(c * E + d * E.t());
  lpZ = lpZ + lpZ.t();
  lpZ.diag() = lpZ.diag() / 2.0;
  
  // For normal model, log p(y|z) is just normal density
  arma::mat llZ = -0.5 * arma::square(Y - Z) / s2;
  llZ = llZ + llZ.t();
  
  // Set NA positions to 0
  arma::uvec na_indices = arma::find_nonfinite(Y);
  llZ.elem(na_indices).zeros();
  
  return lpZ + llZ;
}

// Function to simulate Y from Poisson model given Z
// [[Rcpp::export]]
arma::mat simY_pois(const arma::mat& EZ) {
  int n = EZ.n_rows;
  arma::mat Y(n, n);
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      if(i != j) {
        double lambda = std::exp(EZ(i,j));
        // Cap lambda to avoid numerical issues
        lambda = std::min(lambda, 1e6);
        Y(i,j) = R::rpois(lambda);
      } else {
        Y(i,j) = NA_REAL;
      }
    }
  }
  
  return Y;
}