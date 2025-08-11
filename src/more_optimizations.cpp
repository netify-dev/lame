#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat rmvnorm_cpp(int n, const arma::vec& mu, const arma::mat& Sigma) {
  int p = mu.n_elem;
  arma::mat Sigma_chol;
  
  // Handle non-finite or non-positive definite cases
  if(Sigma.has_nan() || Sigma.has_inf()) {
    Sigma_chol = arma::eye(p, p);
  } else {
    bool chol_success = false;
    
    // Try Cholesky decomposition
    chol_success = arma::chol(Sigma_chol, Sigma, "upper");
    
    if(!chol_success) {
      // Try eigenvalue decomposition
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, Sigma);
      
      // Fix negative eigenvalues
      eigval.elem(arma::find(eigval < 1e-10)).fill(1e-10);
      
      // Reconstruct Sigma_chol
      Sigma_chol = eigvec * arma::diagmat(arma::sqrt(eigval));
    }
  }
  
  // Generate random normal matrix
  arma::mat E = arma::randn(n, p);
  
  // Transform and add mean
  arma::mat result = E * Sigma_chol.t();
  result.each_row() += mu.t();
  
  return result;
}

// [[Rcpp::export]]
arma::mat simZ_cpp(const arma::mat& EZ, double rho, double s2 = 1.0) {
  double c = (std::sqrt(1 + rho) + std::sqrt(1 - rho)) / 2.0;
  double d = (std::sqrt(1 + rho) - std::sqrt(1 - rho)) / 2.0;
  
  int n = EZ.n_rows;
  arma::mat EC = arma::randn(n, n);
  EC = std::sqrt(s2) * (c * EC + d * EC.t());
  
  return EZ + EC;
}

// [[Rcpp::export]]
arma::mat simY_nrm_cpp(const arma::mat& EY, double rho, double s2) {
  arma::mat YS = simZ_cpp(EY, rho, s2);
  YS.diag().fill(NA_REAL);
  return YS;
}

// [[Rcpp::export]]
arma::mat rZ_nrm_fc_cpp(const arma::mat& Z, const arma::mat& EZ, 
                        double rho, double s2, const arma::mat& Y) {
  arma::mat ZS = simY_nrm_cpp(EZ, rho, s2);
  
  // Sample diagonal elements
  int n = Y.n_rows;
  arma::vec diag_vals = arma::randn(n) * std::sqrt(s2 * (1 + rho)) + EZ.diag();
  ZS.diag() = diag_vals;
  
  // Copy Z and update missing values
  arma::mat Z_new = Z;
  arma::uvec na_indices = arma::find_nonfinite(Y);
  Z_new.elem(na_indices) = ZS.elem(na_indices);
  
  return Z_new;
}

// [[Rcpp::export]]
arma::mat mhalf_cpp(const arma::mat& M) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, M);
  
  // Handle negative eigenvalues
  eigval.elem(arma::find(eigval < 1e-10)).fill(1e-10);
  
  return eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t();
}

// [[Rcpp::export]]
arma::cube design_array_cpp(const arma::mat& Xrow, const arma::mat& Xcol, 
                           const arma::cube& Xdyad, bool intercept, 
                           bool n, bool symmetric) {
  int nrow = Xrow.n_rows;
  int ncol = Xcol.n_rows;
  
  // Calculate dimensions
  int pr = Xrow.n_cols;
  int pc = Xcol.n_cols;
  int pd = Xdyad.n_slices;
  
  int p = intercept + pr + pc + pd;
  if(symmetric) {
    p = intercept + pr + pd;
  }
  
  // Initialize design array
  arma::cube X(nrow, ncol, p, arma::fill::zeros);
  
  int idx = 0;
  
  // Add intercept
  if(intercept) {
    X.slice(idx).fill(1.0);
    idx++;
  }
  
  // Add row effects
  for(int k = 0; k < pr; k++) {
    for(int i = 0; i < nrow; i++) {
      X.slice(idx).row(i).fill(Xrow(i, k));
    }
    idx++;
  }
  
  // Add column effects (if not symmetric)
  if(!symmetric) {
    for(int k = 0; k < pc; k++) {
      for(int j = 0; j < ncol; j++) {
        X.slice(idx).col(j).fill(Xcol(j, k));
      }
      idx++;
    }
  }
  
  // Add dyadic effects
  for(int k = 0; k < pd; k++) {
    X.slice(idx) = Xdyad.slice(k);
    idx++;
  }
  
  return X;
}