// Fixed implementation that exactly matches AMEN's rZ_bin_fc
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::export]]
arma::mat rZ_bin_fc_fixed(
    arma::mat Z, const arma::mat& EZ, 
    double rho, arma::mat Y
) {
  int n = Z.n_rows;
  double sz = sqrt(1 - rho * rho);
  
  // Replace NaN with -1 in Y (modifies Y)
  Y.replace(datum::nan, -1);
  
  // Get upper and lower triangle indices
  arma::umat ut = trimatl(ones<umat>(n, n), -1) == 0;  // upper.tri
  arma::umat lt = trimatu(ones<umat>(n, n), 1) == 0;   // lower.tri
  
  // Process each Y value (-1, 0, 1)
  int yVals[3] = {-1, 0, 1};
  double lbVals[3] = {-datum::inf, -datum::inf, 0};
  double ubVals[3] = {datum::inf, 0, datum::inf};
  
  for(int yi = 0; yi < 3; ++yi) {
    double y = yVals[yi];
    double lb = lbVals[yi];
    double ub = ubVals[yi];
    
    // Process upper triangle
    arma::uvec up = find((ut == 1) && (Y == y));
    if(up.n_elem > 0) {
      // AMEN uses t(Z)[up] and t(EZ)[up] for the cross terms
      arma::mat Zt = Z.t();
      arma::mat EZt = EZ.t();
      
      arma::vec ez = EZ(up) + rho * (Zt(up) - EZt(up));
      
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
      Z(up) = zup;
    }
    
    // Process lower triangle
    up = find((lt == 1) && (Y == y));
    if(up.n_elem > 0) {
      // AMEN uses t(Z)[up] and t(EZ)[up] for the cross terms
      arma::mat Zt = Z.t();
      arma::mat EZt = EZ.t();
      
      arma::vec ez = EZ(up) + rho * (Zt(up) - EZt(up));
      
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
      Z(up) = zup;
    }
  }
  
  // Second step: proposal with special covariance structure
  double c = (sqrt(1 + rho) + sqrt(1 - rho)) / 2.0;
  double d = (sqrt(1 + rho) - sqrt(1 - rho)) / 2.0;
  arma::mat E = randn(n, n);
  arma::mat ZP = EZ + c * E + d * E.t();
  
  // Create acceptance matrix
  // A = ((Y == -1) | (sign(ZP) == sign(Y - 0.5)))
  arma::umat A = (Y == -1) || (sign(ZP) == sign(Y - 0.5));
  A.diag().ones(); // Diagonal always accepted
  
  // Make A symmetric: A = A & t(A)
  A = A && A.t();
  
  // Update Z where A is true
  Z.elem(find(A == 1)) = ZP.elem(find(A == 1));
  
  // Update diagonal
  Z.diag() = diagvec(EZ) + sqrt(1 + rho) * randn(n);
  
  return Z;
}