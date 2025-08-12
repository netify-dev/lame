// Bipartite network EZ computation with rectangular G matrix
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Compute expected value for bipartite networks
// base: nA x nB x T with mu + X*beta already added
// a: nA x T (row/sender effects)
// b: nB x T (column/receiver effects)
// U: nA x RA x T (row latent factors)
// V: nB x RB x T (column latent factors)
// G: RA x RB (interaction matrix, possibly rectangular)
//
// [[Rcpp::export]]
arma::cube get_EZ_bip_cpp(const arma::cube& base,
                          const arma::mat& a,
                          const arma::mat& b,
                          const arma::cube& U,
                          const arma::cube& V,
                          const arma::mat& G) {
  const unsigned int nA = base.n_rows;
  const unsigned int nB = base.n_cols;
  const unsigned int T = base.n_slices;
  const unsigned int RA = U.n_cols;
  const unsigned int RB = V.n_cols;
  
  // Validate dimensions
  if(a.n_rows != nA || a.n_cols != T) {
    stop("Dimension mismatch: a must be nA x T");
  }
  if(b.n_rows != nB || b.n_cols != T) {
    stop("Dimension mismatch: b must be nB x T");
  }
  if(U.n_rows != nA || U.n_slices != T) {
    stop("Dimension mismatch: U must be nA x RA x T");
  }
  if(V.n_rows != nB || V.n_slices != T) {
    stop("Dimension mismatch: V must be nB x RB x T");
  }
  if(G.n_rows != RA || G.n_cols != RB) {
    stop("Dimension mismatch: G must be RA x RB");
  }
  
  arma::cube EZ = base; // copy
  
  for (unsigned int t = 0; t < T; ++t) {
    arma::mat Et = EZ.slice(t);
    
    // Add row effects (a_i to each column)
    Et.each_col() += a.col(t);
    
    // Add column effects (b_j to each row)  
    Et.each_row() += b.col(t).t();
    
    // Add bilinear term U * G * V^T
    if(RA > 0 && RB > 0) {
      const arma::mat Ut = U.slice(t);  // nA x RA
      const arma::mat Vt = V.slice(t);  // nB x RB
      Et += Ut * G * Vt.t();            // nA x nB
    }
    
    EZ.slice(t) = Et;
  }
  
  return EZ;
}

// Helper to build additive effects outer product for bipartite
// a: nA x T, b: nB x T
// Returns nA x nB x T array
// [[Rcpp::export]]
arma::cube outer_ab_bip_cpp(const arma::mat& a, const arma::mat& b) {
  const unsigned int nA = a.n_rows;
  const unsigned int nB = b.n_rows;
  const unsigned int T = a.n_cols;
  
  arma::cube out(nA, nB, T);
  
  for(unsigned int t = 0; t < T; ++t) {
    arma::mat slice_t(nA, nB);
    slice_t.each_col() = a.col(t);      // Broadcast a across columns
    slice_t.each_row() += b.col(t).t(); // Add b across rows
    out.slice(t) = slice_t;
  }
  
  return out;
}