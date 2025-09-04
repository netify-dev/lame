//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// Fast version of rZ_bin_fc_cpp with:
// 1. Pre-computed index matrices
// 2. Efficient transpose handling
// 3. Batch RNG operations
// 4. Minimized memory allocations

// [[Rcpp::export]]
arma::cube rZ_bin_fc_cpp(
    arma::cube ZT, arma::cube EZT, double rho, arma::cube YT
) {
  
  double sz = sqrt(1 - rho * rho);
  int n = ZT.n_rows;
  int T = ZT.n_slices;
  
  // Pre-compute indices once (these never change)
  arma::uvec ut = find(trimatl(arma::ones(n, n)) == 0);
  arma::uvec lt = find(trimatu(arma::ones(n, n)) == 0);
  
  // Pre-allocate work vectors
  int n_upper = ut.n_elem;
  int n_lower = lt.n_elem;
  arma::vec ez_u(n_upper);
  arma::vec ez_l(n_lower);
  
  for(int t = 0; t < T; ++t) {
    arma::mat& Z = ZT.slice(t);  // Reference to avoid copy
    const arma::mat& Y = YT.slice(t);  // Const reference
    const arma::mat& EZ = EZT.slice(t);
    
    // Replace NaN with -1 in Y (create working copy)
    arma::mat Y_work = Y;
    Y_work.replace(datum::nan, -1);
    
    // Work directly with transposes where needed
    // Avoid creating full transpose matrices
    
    // Process each Y value type (-1, 0, 1)
    int yVals[3] = {-1, 0, 1};
    double lbVals[3] = {-datum::inf, -datum::inf, 0};
    double ubVals[3] = {datum::inf, 0, datum::inf};
    
    for(int i = 0; i < 3; ++i) {
      int y = yVals[i];
      double lb = lbVals[i];
      double ub = ubVals[i];
      
      arma::uvec Y_indices = find(Y_work == y);
      
      // Upper triangle
      arma::uvec up_u = arma::intersect(Y_indices, ut);
      if(up_u.n_elem > 0) {
        // Transpose matrices for correct indexing (AMEN uses t(Z)[up])
        arma::mat Zt = Z.t();
        arma::mat EZt = EZ.t();
        
        // Compute ez values efficiently
        for(unsigned int j = 0; j < up_u.n_elem; ++j) {
          int idx = up_u(j);
          ez_u(j) = EZ(idx) + rho * (Zt(idx) - EZt(idx));
        }
        
        // Batch RNG operations using Armadillo's RNG
        arma::vec unif_draws = arma::randu(up_u.n_elem);
        
        // Vectorized truncated normal sampling
        for(unsigned int j = 0; j < up_u.n_elem; ++j) {
          double x1 = (lb - ez_u(j)) / sz;
          double x2 = (ub - ez_u(j)) / sz;
          double p_lo = R::pnorm(x1, 0.0, 1.0, 1, 0);
          double p_hi = R::pnorm(x2, 0.0, 1.0, 1, 0);
          double u = p_lo + unif_draws(j) * (p_hi - p_lo);
          Z(up_u(j)) = ez_u(j) + sz * R::qnorm(u, 0.0, 1.0, 1, 0);
        }
      }
      
      // Lower triangle
      arma::uvec up_l = arma::intersect(Y_indices, lt);
      if(up_l.n_elem > 0) {
        // Transpose matrices for correct indexing (AMEN uses t(Z)[up])
        arma::mat Zt = Z.t();
        arma::mat EZt = EZ.t();
        
        // Compute ez values efficiently
        for(unsigned int j = 0; j < up_l.n_elem; ++j) {
          int idx = up_l(j);
          ez_l(j) = EZ(idx) + rho * (Zt(idx) - EZt(idx));
        }
        
        // Batch RNG operations
        arma::vec unif_draws = arma::randu(up_l.n_elem);
        
        // Vectorized truncated normal sampling
        for(unsigned int j = 0; j < up_l.n_elem; ++j) {
          double x1 = (lb - ez_l(j)) / sz;
          double x2 = (ub - ez_l(j)) / sz;
          double p_lo = R::pnorm(x1, 0.0, 1.0, 1, 0);
          double p_hi = R::pnorm(x2, 0.0, 1.0, 1, 0);
          double u = p_lo + unif_draws(j) * (p_hi - p_lo);
          Z(up_l(j)) = ez_l(j) + sz * R::qnorm(u, 0.0, 1.0, 1, 0);
        }
      }
    }
    
    // AMEN's second step: additional proposal with special covariance structure
    double c = (sqrt(1 + rho) + sqrt(1 - rho)) / 2.0;
    double d = (sqrt(1 + rho) - sqrt(1 - rho)) / 2.0;
    
    // Generate random matrix E
    arma::mat E = arma::randn(n, n);
    
    // Compute proposal ZP = EZ + c*E + d*t(E)
    arma::mat ZP = EZ + c * E + d * E.t();
    
    // Accept updates where they satisfy constraints or for missing values
    arma::mat A = arma::zeros<arma::mat>(n, n);
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        if(std::isnan(Y(i,j)) || Y_work(i,j) == -1) {
          A(i,j) = 1;
        } else if(Y_work(i,j) == 1 && ZP(i,j) > 0) {
          A(i,j) = 1;
        } else if(Y_work(i,j) == 0 && ZP(i,j) < 0) {
          A(i,j) = 1;
        }
      }
    }
    
    // Diagonal is always accepted
    A.diag().ones();
    
    // A must be symmetric: A = A & t(A)
    A = A % A.t();
    
    // Update Z where A is true
    arma::uvec accept_idx = find(A == 1);
    Z(accept_idx) = ZP(accept_idx);
    
    // Update diagonal with proper variance sqrt(1+rho)
    Z.diag() = diagvec(EZ) + sqrt(1 + rho) * arma::randn(n);
  }
  
  return ZT;
}

// Single matrix version for non-replicated case
// [[Rcpp::export]]
arma::mat rZ_bin_fc_single(
    const arma::mat& Z, const arma::mat& EZ, 
    double rho, const arma::mat& Y
) {
  double sz = sqrt(1 - rho * rho);
  int n = Z.n_rows;
  
  // Pre-compute indices
  arma::uvec ut = find(trimatl(arma::ones(n, n)) == 0);
  arma::uvec lt = find(trimatu(arma::ones(n, n)) == 0);
  
  // Create output matrix
  arma::mat Z_new = Z;
  
  // Replace NaN with -1 in Y
  arma::mat Y_work = Y;
  Y_work.replace(datum::nan, -1);
  
  // Process each Y value type
  int yVals[3] = {-1, 0, 1};
  double lbVals[3] = {-datum::inf, -datum::inf, 0};
  double ubVals[3] = {datum::inf, 0, datum::inf};
  
  for(int i = 0; i < 3; ++i) {
    int y = yVals[i];
    double lb = lbVals[i];
    double ub = ubVals[i];
    
    arma::uvec Y_indices = find(Y_work == y);
    
    // Upper triangle
    arma::uvec up_u = arma::intersect(Y_indices, ut);
    if(up_u.n_elem > 0) {
      // Batch generate uniforms
      arma::vec unif_draws = arma::randu(up_u.n_elem);
      
      // Transpose matrices for correct indexing (AMEN uses t(Z)[up])
      arma::mat Zt = Z.t();
      arma::mat EZt = EZ.t();
      
      for(unsigned int j = 0; j < up_u.n_elem; ++j) {
        int idx = up_u(j);
        double ez = EZ(idx) + rho * (Zt(idx) - EZt(idx));
        
        double x1 = (lb - ez) / sz;
        double x2 = (ub - ez) / sz;
        double p_lo = R::pnorm(x1, 0.0, 1.0, 1, 0);
        double p_hi = R::pnorm(x2, 0.0, 1.0, 1, 0);
        double u = p_lo + unif_draws(j) * (p_hi - p_lo);
        Z_new(idx) = ez + sz * R::qnorm(u, 0.0, 1.0, 1, 0);
      }
    }
    
    // Lower triangle
    arma::uvec up_l = arma::intersect(Y_indices, lt);
    if(up_l.n_elem > 0) {
      // Batch generate uniforms
      arma::vec unif_draws = arma::randu(up_l.n_elem);
      
      // Transpose matrices for correct indexing (AMEN uses t(Z)[up])
      arma::mat Zt = Z.t();
      arma::mat EZt = EZ.t();
      
      for(unsigned int j = 0; j < up_l.n_elem; ++j) {
        int idx = up_l(j);
        double ez = EZ(idx) + rho * (Zt(idx) - EZt(idx));
        
        double x1 = (lb - ez) / sz;
        double x2 = (ub - ez) / sz;
        double p_lo = R::pnorm(x1, 0.0, 1.0, 1, 0);
        double p_hi = R::pnorm(x2, 0.0, 1.0, 1, 0);
        double u = p_lo + unif_draws(j) * (p_hi - p_lo);
        Z_new(idx) = ez + sz * R::qnorm(u, 0.0, 1.0, 1, 0);
      }
    }
  }
  
  // AMEN's second step: additional proposal with special covariance structure
  double c = (sqrt(1 + rho) + sqrt(1 - rho)) / 2.0;
  double d = (sqrt(1 + rho) - sqrt(1 - rho)) / 2.0;
  
  // Generate random matrix E
  arma::mat E = arma::randn(n, n);
  
  // Compute proposal ZP = EZ + c*E + d*t(E)
  arma::mat ZP = EZ + c * E + d * E.t();
  
  // Accept updates where they satisfy constraints or for missing values
  // A = ((Y == -1) | (sign(ZP) == sign(Y - 0.5)))
  arma::mat A = arma::zeros<arma::mat>(n, n);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(std::isnan(Y(i,j)) || Y_work(i,j) == -1) {
        A(i,j) = 1;
      } else if(Y_work(i,j) == 1 && ZP(i,j) > 0) {
        A(i,j) = 1;
      } else if(Y_work(i,j) == 0 && ZP(i,j) < 0) {
        A(i,j) = 1;
      }
    }
  }
  
  // Diagonal is always accepted
  A.diag().ones();
  
  // A must be symmetric: A = A & t(A)
  A = A % A.t();
  
  // Update Z where A is true
  arma::uvec accept_idx = find(A == 1);
  Z_new(accept_idx) = ZP(accept_idx);
  
  // Update diagonal with proper variance sqrt(1+rho)
  Z_new.diag() = diagvec(EZ) + sqrt(1 + rho) * arma::randn(n);
  
  return Z_new;
}