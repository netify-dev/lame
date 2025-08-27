//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// Fast Wishart sampler with numerical stability - handles all edge cases
// [[Rcpp::export]]
arma::mat rwish_cpp(const arma::mat& S0, int nu) {
  int n = S0.n_rows;
  
  // Handle non-finite values by replacing with reasonable defaults
  arma::mat S0_clean = S0;
  if (!S0.is_finite()) {
    // Replace non-finite values with identity matrix structure
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(!std::isfinite(S0_clean(i,j))) {
          S0_clean(i,j) = (i == j) ? 1.0 : 0.0;
        }
      }
    }
  }
  
  // Ensure S0 is symmetric
  arma::mat S0_sym = 0.5 * (S0_clean + S0_clean.t());
  
  // Check for near-zero or negative eigenvalues and fix them
  arma::vec eigval;
  arma::mat eigvec;
  bool eig_success = arma::eig_sym(eigval, eigvec, S0_sym);
  
  if (!eig_success) {
    // If eigen decomposition fails completely, use identity matrix
    S0_sym = arma::eye(n, n);
    arma::eig_sym(eigval, eigvec, S0_sym);
  }
  
  // Fix eigenvalues to ensure positive definiteness
  double min_eig = eigval.min();
  double max_eig = eigval.max();
  
  // Handle extreme eigenvalue ratios (near-singular matrices)
  if (max_eig > 0 && min_eig/max_eig < 1e-10) {
    // Regularize by setting minimum eigenvalue
    double regularization = max_eig * 1e-8;
    for(int i = 0; i < eigval.n_elem; i++) {
      if(eigval(i) < regularization) {
        eigval(i) = regularization;
      }
    }
    // Reconstruct matrix with regularized eigenvalues
    S0_sym = eigvec * arma::diagmat(eigval) * eigvec.t();
  } else if (min_eig <= 0) {
    // Add ridge to all eigenvalues
    double ridge = std::max(1e-6, 1e-6 - min_eig);
    eigval += ridge;
    S0_sym = eigvec * arma::diagmat(eigval) * eigvec.t();
  }
  
  // Ensure symmetry after reconstruction
  S0_sym = 0.5 * (S0_sym + S0_sym.t());
  
  // Generate Wishart sample using Cholesky or eigenvalue decomposition
  arma::mat sS0;
  bool chol_success = arma::chol(sS0, S0_sym, "lower");
  
  if (!chol_success) {
    // Use eigenvalue decomposition for generation
    // This should rarely happen after the fixes above
    arma::eig_sym(eigval, eigvec, S0_sym);
    eigval = arma::max(eigval, arma::ones(n) * 1e-10);
    arma::mat sqrt_eigval = arma::diagmat(arma::sqrt(eigval));
    sS0 = eigvec * sqrt_eigval;
  }
  
  // Generate random matrix
  arma::mat Z = arma::randn(nu, n) * sS0.t();
  
  return Z.t() * Z;
}

// Fast UV update with efficient memory usage
// [[Rcpp::export]]
List rUV_rep_fc_cpp(
    const arma::cube& ET, arma::mat U, arma::mat V, 
    double rho, double s2, const arma::mat& iSe2, 
    double maxmargin, bool shrink, const NumericVector& rLoopIDs) {
  
  int Time = ET.n_slices;
  int R = U.n_cols;
  int n = U.n_rows;
  
  // Pre-allocate work matrices
  arma::mat UV = join_rows(U, V);
  arma::mat Suv;
  
  // Compute Suv once
  if(shrink) {
    arma::mat UV_prod = UV.t() * UV;
    arma::mat inv_mat = inv_sympd(eye<mat>(R*2, R*2) + UV_prod);
    Suv = inv_sympd(rwish_cpp(inv_mat, n + R + 2));
  } else {
    Suv = eye<mat>(R*2, R*2);
  }
  
  double g = iSe2(0,0);
  double d = iSe2(0,1);
  double g2d2 = g*g + d*d;
  
  // Pre-allocate Est matrix
  arma::mat Est(n * n, Time);
  
  // Pre-compute some matrices that don't change in the loop
  arma::mat Usmall(n, R-1);
  arma::mat Vsmall(n, R-1);
  
  for(int i = 0; i < rLoopIDs.size(); i++) {
    int r = rLoopIDs[i];
    
    // Build Usmall and Vsmall efficiently
    if(r == 0) {
      Usmall = U.cols(1, R-1);
      Vsmall = V.cols(1, R-1);
    } else if(r == R-1) {
      Usmall = U.cols(0, R-2);
      Vsmall = V.cols(0, R-2);
    } else {
      Usmall = join_rows(U.cols(0, r-1), U.cols(r+1, R-1));
      Vsmall = join_rows(V.cols(0, r-1), V.cols(r+1, R-1));
    }
    
    arma::mat UVmr = Usmall * Vsmall.t();
    
    // Vectorized computation of Est
    for(int t = 0; t < Time; t++) {
      arma::mat ert = ET.slice(t) - UVmr;
      Est.col(t) = vectorise(g2d2 * ert + 2 * g * d * ert.t());
    }
    
    // Update U
    arma::vec vr = V.col(r);
    
    // Extract relevant parts of Suv efficiently
    arma::rowvec Suvsmall = Suv.row(r);
    Suvsmall.shed_col(r);
    
    arma::mat Suvsmall2 = Suv;
    Suvsmall2.shed_col(r);
    Suvsmall2.shed_row(r);
    
    arma::vec b0 = Suvsmall * inv_sympd(Suvsmall2);
    double v0 = Suv(r,r) - dot(b0, Suv.col(r).subvec(0, Suv.n_rows-2));
    
    arma::vec m0 = join_rows(Usmall, V) * b0;
    
    double sumvr2 = dot(vr, vr);
    double ssv = (sumvr2 >= maxmargin) ? sumvr2 : maxmargin;
    
    double a = Time * g2d2 * ssv + 1/v0;
    double c = -2 * Time * g * d / (a*a + a * 2 * Time * g * d * ssv);
    
    arma::vec Esv = Est * vr;
    arma::vec m1 = Esv/a + c * vr * dot(Esv + m0/v0, vr) + m0/(a*v0);
    
    double ah = sqrt(1/a);
    double bh = (sqrt(1/a + ssv*c) - sqrt(1/a)) / ssv;
    
    arma::vec e = randn(n);
    U.col(r) = m1 + ah * e + bh * vr * dot(vr, e);
    
    // Update V (similar process)
    arma::vec ur = U.col(r);
    int rv = R + r;
    
    arma::rowvec Suvsmall_v = Suv.row(rv);
    Suvsmall_v.shed_col(rv);
    
    arma::mat Suvsmall2_v = Suv;
    Suvsmall2_v.shed_col(rv);
    Suvsmall2_v.shed_row(rv);
    
    arma::vec b0_v = Suvsmall_v * inv_sympd(Suvsmall2_v);
    double v0_v = Suv(rv,rv) - dot(b0_v, Suv.col(rv).subvec(0, Suv.n_rows-2));
    
    arma::vec m0_v = join_rows(U, Vsmall) * b0_v;
    
    double sumur2 = dot(ur, ur);
    double ssu = (sumur2 >= maxmargin) ? sumur2 : maxmargin;
    
    double a_v = Time * g2d2 * ssu + 1/v0_v;
    double c_v = -2 * Time * g * d / (a_v*a_v + a_v * 2 * Time * g * d * ssu);
    
    arma::vec tEsu = Est.t() * ur;
    arma::vec m1_v = tEsu/a_v + c_v * ur * dot(tEsu + m0_v/v0_v, ur) + m0_v/(a_v*v0_v);
    
    double ah_v = sqrt(1/a_v);
    double bh_v = (sqrt(1/a_v + ssu*c_v) - sqrt(1/a_v)) / ssu;
    
    arma::vec e_v = randn(n);
    V.col(r) = m1_v + ah_v * e_v + bh_v * ur * dot(ur, e_v);
  }
  
  return List::create(
    Named("U") = U,
    Named("V") = V
  );
}