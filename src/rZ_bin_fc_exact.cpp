// Exact implementation of AMEN's rZ_bin_fc
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::export]]
arma::mat rZ_bin_fc_exact(
    arma::mat Z, const arma::mat& EZ, 
    double rho, arma::mat Y
) {
  int n = Z.n_rows;
  double sz = sqrt(1 - rho * rho);
  
  // Replace NaN with -1 in Y
  Y.replace(datum::nan, -1);
  
  // Create transpose matrices once
  arma::mat Zt = Z.t();
  arma::mat EZt = EZ.t();
  
  // Process each Y value (-1, 0, 1)
  std::vector<int> yvals = {-1, 0, 1};
  std::vector<double> lb_vals = {-datum::inf, -datum::inf, 0};
  std::vector<double> ub_vals = {datum::inf, 0, datum::inf};
  
  for(int yi = 0; yi < 3; ++yi) {
    int y = yvals[yi];
    double lb = lb_vals[yi];
    double ub = ub_vals[yi];
    
    // Process upper and lower triangles
    for(int tri = 1; tri <= 2; ++tri) {
      arma::uvec up;
      
      if(tri == 1) {
        // Upper triangle: i < j
        up = find((trimatu(Y, 1) == y));
      } else {
        // Lower triangle: i > j  
        up = find((trimatl(Y, -1) == y));
      }
      
      if(up.n_elem > 0) {
        // Compute ez = EZ[up] + rho*(t(Z)[up] - t(EZ)[up])
        arma::vec ez = EZ.elem(up) + rho * (Zt.elem(up) - EZt.elem(up));
        
        // Truncated normal sampling
        arma::vec u_rand = randu(up.n_elem);
        arma::vec zup(up.n_elem);
        
        for(unsigned int j = 0; j < up.n_elem; ++j) {
          double plo = R::pnorm((lb - ez(j))/sz, 0.0, 1.0, 1, 0);
          double phi = R::pnorm((ub - ez(j))/sz, 0.0, 1.0, 1, 0);
          double u = plo + u_rand(j) * (phi - plo);
          zup(j) = ez(j) + sz * R::qnorm(u, 0.0, 1.0, 1, 0);
          
          // Handle infinities
          if(!std::isfinite(zup(j))) {
            zup(j) = Z(up(j));
          }
        }
        
        // Update Z
        Z.elem(up) = zup;
      }
    }
  }
  
  // Second step: proposal with special covariance structure  
  double c = (sqrt(1 + rho) + sqrt(1 - rho)) / 2.0;
  double d = (sqrt(1 + rho) - sqrt(1 - rho)) / 2.0;
  arma::mat E = randn(n, n);
  arma::mat Et = E.t();
  arma::mat ZP = EZ + c * E + d * Et;
  
  // Create acceptance matrix
  // A = ((Y == -1) | (sign(ZP) == sign(Y - 0.5)))
  arma::umat A = zeros<umat>(n, n);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(Y(i,j) == -1) {
        A(i,j) = 1;
      } else if(Y(i,j) == 1 && ZP(i,j) > 0) {
        A(i,j) = 1;  
      } else if(Y(i,j) == 0 && ZP(i,j) < 0) {
        A(i,j) = 1;
      }
    }
  }
  A.diag().ones(); // Diagonal always accepted
  
  // Make A symmetric: A = A & t(A)
  arma::umat At = A.t();
  A = A % At;  // Element-wise AND
  
  // Update Z where A is true
  arma::uvec accept_idx = find(A == 1);
  Z.elem(accept_idx) = ZP.elem(accept_idx);
  
  // Update diagonal
  arma::vec diag_mean = diagvec(EZ);
  Z.diag() = diag_mean + sqrt(1 + rho) * randn(n);
  
  return Z;
}