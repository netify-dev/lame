// Optimized UV sampling functions
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Optimized Wishart sampler - avoid unnecessary copies
// [[Rcpp::export]]
arma::mat rwish_opt_cpp(const arma::mat& S0, int nu) {
  int p = S0.n_rows;
  // Use in-place operations where possible
  mat Z(nu, p, fill::randn);
  mat cholS0 = chol(S0);
  // In-place multiplication
  Z *= cholS0;
  return Z.t() * Z;
}

// [[Rcpp::export]]
List rUV_rep_opt_cpp(
    const arma::cube& ET, arma::mat U, arma::mat V, 
    double rho, double s2, const arma::mat& iSe2, 
    double maxmargin, bool shrink, const arma::uvec& rLoopIDs) {
  
  const int Time = ET.n_slices;
  const int R = U.n_cols;
  const int n = U.n_rows;
  const int n2 = n * n;
  
  // Pre-allocate all matrices we'll reuse
  mat UV(n, 2*R);
  mat Suv(2*R, 2*R);
  mat UVmr(n, n);
  mat Est(n2, Time);
  mat ert(n, n);
  vec Esvvec(n2);
  mat Esv(n, 1);
  mat tEsu(n, 1);
  
  // Combine U and V efficiently
  UV.cols(0, R-1) = U;
  UV.cols(R, 2*R-1) = V;
  
  // Compute shrinkage matrix once
  if(shrink) {
    mat UVtUV = UV.t() * UV;
    UVtUV.diag() += 1.0;
    Suv = inv_sympd(rwish_opt_cpp(inv_sympd(UVtUV), n+R+2));
  } else {
    Suv.eye();
  }
  
  // Pre-compute constants
  const double g = iSe2(0,0);
  const double d = iSe2(0,1);
  const double g2d2 = g*g + d*d;
  const double gd2 = 2.0*g*d;
  
  // Pre-allocate working vectors
  vec e(n);
  vec vr(n);
  vec ur(n);
  
  // Main loop - process each latent dimension
  for(unsigned int i = 0; i < rLoopIDs.n_elem; i++) {
    const int r = rLoopIDs[i];
    
    // Compute UVmr = U_{-r} * V_{-r}^T efficiently
    UVmr = U * V.t();
    // Subtract contribution of r-th component
    UVmr -= U.col(r) * V.col(r).t();
    
    // Compute Est matrix for all time points
    for(int t = 0; t < Time; t++) {
      ert = ET.slice(t) - UVmr;
      // Vectorize and compute in one step
      Est.col(t) = vectorise(g2d2 * ert + gd2 * ert.t());
    }
    
    // === Update U ===
    vr = V.col(r);
    const double sumvr2 = std::max(dot(vr, vr), maxmargin);
    
    // Extract relevant parts of Suv efficiently
    const int idx_u = r;
    vec Suv_row = Suv.row(idx_u).t();
    mat Suv_sub = Suv;
    Suv_sub.shed_row(idx_u);
    Suv_sub.shed_col(idx_u);
    
    vec Suv_col = Suv.col(idx_u);
    Suv_col.shed_row(idx_u);
    
    vec Suv_row_sub = Suv_row;
    Suv_row_sub.shed_row(idx_u);
    
    // Compute b0 and v0 using solve for efficiency
    vec b0 = solve(Suv_sub.t(), Suv_row_sub);
    double v0 = Suv(idx_u, idx_u) - dot(b0, Suv_col);
    
    // Compute m0 efficiently
    UV.shed_col(idx_u);
    vec m0 = UV * b0;
    UV.insert_cols(idx_u, U.col(r)); // Restore for next iteration
    
    // Compute update parameters
    const double a = Time * g2d2 * sumvr2 + 1.0/v0;
    const double c = -gd2 / (a*a + a * gd2 * Time * sumvr2);
    
    // Sum across time efficiently
    Esvvec = sum(Est, 1);
    Esv = reshape(Esvvec, n, n) * vr;
    
    // Compute mean
    const double vr_dot = dot(vr, Esv + m0/v0);
    vec m1 = Esv/a + c * vr * vr_dot + m0/(a*v0);
    
    // Sample new U
    const double ah = 1.0/sqrt(a);
    const double bh = (sqrt(1.0/a + sumvr2*c) - ah) / sumvr2;
    e.randn();
    U.col(r) = m1 + ah * e + bh * vr * dot(vr, e);
    
    // === Update V ===
    ur = U.col(r);
    const double sumur2 = std::max(dot(ur, ur), maxmargin);
    
    // Extract relevant parts for V update
    const int idx_v = R + r;
    Suv_row = Suv.row(idx_v).t();
    Suv_sub = Suv;
    Suv_sub.shed_row(idx_v);
    Suv_sub.shed_col(idx_v);
    
    Suv_col = Suv.col(idx_v);
    Suv_col.shed_row(idx_v);
    
    Suv_row_sub = Suv_row;
    Suv_row_sub.shed_row(idx_v);
    
    // Compute b0_v and v0_v
    vec b0_v = solve(Suv_sub.t(), Suv_row_sub);
    double v0_v = Suv(idx_v, idx_v) - dot(b0_v, Suv_col);
    
    // Compute m0_v
    UV.shed_col(idx_v);
    vec m0_v = UV * b0_v;
    UV.insert_cols(idx_v, V.col(r)); // Restore
    
    // Compute update parameters for V
    const double a_v = Time * g2d2 * sumur2 + 1.0/v0_v;
    const double c_v = -gd2 / (a_v*a_v + a_v * gd2 * Time * sumur2);
    
    // Compute tEsu efficiently
    tEsu = reshape(Esvvec, n, n).t() * ur;
    
    // Compute mean for V
    const double ur_dot = dot(ur, tEsu + m0_v/v0_v);
    vec m1_v = tEsu/a_v + c_v * ur * ur_dot + m0_v/(a_v*v0_v);
    
    // Sample new V
    const double ah_v = 1.0/sqrt(a_v);
    const double bh_v = (sqrt(1.0/a_v + sumur2*c_v) - ah_v) / sumur2;
    e.randn();
    V.col(r) = m1_v + ah_v * e + bh_v * ur * dot(ur, e);
  }
  
  return List::create(Named("U") = U, Named("V") = V);
}

// Optimized symmetric UV sampler
// [[Rcpp::export]]
List rUV_sym_opt_cpp(
    const arma::mat& E, arma::mat U, arma::mat V, 
    double s2, bool shrink, const arma::uvec& uLoopIDs) {
  
  const int R = U.n_cols;
  const int n = U.n_rows;
  
  // Initialize L efficiently
  vec L_diag(R);
  for(int r = 0; r < R; r++) {
    L_diag(r) = (U(0,r) != 0) ? V(0,r) / U(0,r) : 0.0;
  }
  mat L = diagmat(L_diag);
  
  // Pre-allocate matrices
  mat ivDiagMat(R, R);
  vec ivU(R);
  mat eui(n, R);
  mat iQ(R, R);
  vec l(R);
  vec randNorm(R);
  
  // Compute shrinkage
  if(shrink) {
    const double shape = (2.0 + n) / 2.0;
    vec scale = 0.5 * (1.0 + sum(square(U), 0).t());
    for(int r = 0; r < R; r++) {
      ivU(r) = R::rgamma(shape, 1.0/scale(r));
    }
    ivDiagMat = diagmat(ivU);
  } else {
    ivDiagMat.eye();
    ivDiagMat *= n;
  }
  
  // Pre-compute U^T * U once
  mat UtU = U.t() * U;
  
  // Main loop for updating U
  for(unsigned int s = 0; s < uLoopIDs.n_elem; s++) {
    const int i = uLoopIDs(s);
    
    // Compute eui = U * E[i,:]^T efficiently
    const vec ei = E.row(i).t();
    for(int r = 0; r < R; r++) {
      eui.col(r) = U.col(r) * ei(i);
      for(int j = 0; j < n; j++) {
        if(j != i) {
          eui(j,r) = U(j,r) * ei(j);
        }
      }
    }
    
    // Compute update parameters
    rowvec euisum = sum(eui, 0);
    l = L * (euisum.t() - U.row(i).t() * E(i,i)) / s2;
    
    // Update UtU temporarily
    mat UtU_temp = UtU - U.row(i).t() * U.row(i);
    
    // Compute inverse using Cholesky for symmetric positive definite
    mat Q = ivDiagMat + (L * UtU_temp * L) / s2;
    iQ = inv_sympd(Q);
    
    // Sample new U[i,:]
    // Ensure iQ is symmetric before Cholesky
    iQ = 0.5 * (iQ + iQ.t());
    
    // Add numerical safeguard for positive definiteness
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, iQ);
    double min_eigval = eigval.min();
    if (min_eigval < 1e-10) {
      iQ += (1e-10 - min_eigval + 1e-6) * eye(R, R);
    }
    
    randNorm.randn();
    mat cholIQ;
    bool chol_success = chol(cholIQ, iQ, "lower");
    if (!chol_success) {
      // Fallback: use diagonal approximation
      cholIQ = diagmat(sqrt(diagvec(iQ)));
    }
    U.row(i) = (iQ * l + cholIQ * randNorm).t();
    
    // Update UtU for next iteration
    UtU = UtU_temp + U.row(i).t() * U.row(i);
  }
  
  // Pre-compute lower triangular mask once
  mat tmponesmat = trimatl(ones<mat>(n, n));
  uvec tmpindex = find(tmponesmat == 0);
  
  // Update L diagonal elements
  mat Usmall(n, R-1);
  mat Lsmall(R-1, R-1);
  mat Er(n, n);
  mat uut(n, n);
  
  for(int r = 0; r < R; r++) {
    // Create Usmall and Lsmall without column/row r
    int col_idx = 0;
    for(int c = 0; c < R; c++) {
      if(c != r) {
        Usmall.col(col_idx) = U.col(c);
        int row_idx = 0;
        for(int rr = 0; rr < R; rr++) {
          if(rr != r) {
            Lsmall(row_idx, col_idx) = L(rr, c);
            row_idx++;
          }
        }
        col_idx++;
      }
    }
    
    // Compute residual
    Er = E - Usmall * Lsmall * Usmall.t();
    
    // Compute outer product once
    uut = U.col(r) * U.col(r).t();
    
    // Compute likelihood and precision
    double l_val = accu(Er.elem(tmpindex) % uut.elem(tmpindex)) / s2;
    double iq = 1.0 / (1.0 + accu(square(uut.elem(tmpindex))) / s2);
    
    // Sample new L[r,r]
    L(r,r) = R::rnorm(iq * l_val, sqrt(iq));
  }
  
  // Compute V = U * L
  V = U * L;
  
  return List::create(Named("U") = U, Named("V") = V);
}