// UV and G updates for bipartite networks
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Utility: symmetric inverse with ridge for numerical stability
static inline mat inv_spd(const mat& A) {
  mat As = 0.5 * (A + A.t());  // Ensure symmetry
  As.diag() += 1e-8;            // Small ridge for stability
  return inv_sympd(As);
}

// Sample U rows given V, G and residual cube R (nA x nB x T)
// lambdaU: RA-length vector of prior precisions (>=0)
// s2: T-length vector of residual variances (per-slice); pass scalar as length-1
// [[Rcpp::export]]
arma::cube sample_U_bip_cpp(const arma::cube& R,
                            const arma::cube& V,
                            const arma::mat& G,
                            const arma::vec& lambdaU,
                            const arma::vec& s2) {
  const unsigned int nA = R.n_rows;
  const unsigned int nB = R.n_cols;
  const unsigned int T = R.n_slices;
  const unsigned int RB = V.n_cols;
  const unsigned int RA = G.n_rows;
  
  arma::cube U(nA, RA, T, fill::zeros);
  
  // Precompute per-t crossprod of M_t = V_t * G.t() : (nB x RA)
  std::vector<mat> MtXMt(T);
  for (unsigned int t = 0; t < T; ++t) {
    mat Mt = V.slice(t) * G.t();  // nB x RA
    double var_t = (s2.n_elem == 1) ? s2(0) : s2(t);
    MtXMt[t] = Mt.t() * Mt / var_t;  // RA x RA
  }
  
  mat LambdaU_diag = diagmat(lambdaU);  // RA x RA
  
  // Sample each row of U independently
  for (unsigned int i = 0; i < nA; ++i) {
    mat Prec(RA, RA, fill::zeros);
    vec rhs(RA, fill::zeros);
    
    for (unsigned int t = 0; t < T; ++t) {
      double w = (s2.n_elem == 1) ? 1.0/s2(0) : 1.0/s2(t);
      mat Mt = V.slice(t) * G.t();  // nB x RA
      Prec += MtXMt[t];
      rhs += Mt.t() * R.slice(t).row(i).t() * w;
    }
    
    Prec += LambdaU_diag;
    mat Sigma = inv_spd(Prec);
    vec mean = Sigma * rhs;
    
    // Draw from multivariate normal
    vec draw = mean + chol(Sigma, "lower") * randn<vec>(RA);
    
    // Copy to all time slices (static U) or per-t if dynamic
    for (unsigned int t = 0; t < T; ++t) {
      for(unsigned int r = 0; r < RA; ++r) {
        U(i, r, t) = draw(r);
      }
    }
  }
  
  return U;
}

// Sample V rows given U, G and residual cube R
// [[Rcpp::export]]
arma::cube sample_V_bip_cpp(const arma::cube& R,
                            const arma::cube& U,
                            const arma::mat& G,
                            const arma::vec& lambdaV,
                            const arma::vec& s2) {
  const unsigned int nA = R.n_rows;
  const unsigned int nB = R.n_cols;
  const unsigned int T = R.n_slices;
  const unsigned int RA = U.n_cols;
  const unsigned int RB = G.n_cols;
  
  arma::cube V(nB, RB, T, fill::zeros);
  
  // Precompute per-t crossprod of N_t = U_t * G : (nA x RB)
  std::vector<mat> NtXNt(T);
  for (unsigned int t = 0; t < T; ++t) {
    mat Nt = U.slice(t) * G;  // nA x RB
    double var_t = (s2.n_elem == 1) ? s2(0) : s2(t);
    NtXNt[t] = Nt.t() * Nt / var_t;  // RB x RB
  }
  
  mat LambdaV_diag = diagmat(lambdaV);  // RB x RB
  
  // Sample each row of V independently
  for (unsigned int j = 0; j < nB; ++j) {
    mat Prec(RB, RB, fill::zeros);
    vec rhs(RB, fill::zeros);
    
    for (unsigned int t = 0; t < T; ++t) {
      double w = (s2.n_elem == 1) ? 1.0/s2(0) : 1.0/s2(t);
      mat Nt = U.slice(t) * G;  // nA x RB
      Prec += NtXNt[t];
      rhs += Nt.t() * R.slice(t).col(j) * w;
    }
    
    Prec += LambdaV_diag;
    mat Sigma = inv_spd(Prec);
    vec mean = Sigma * rhs;
    
    // Draw from multivariate normal
    vec draw = mean + chol(Sigma, "lower") * randn<vec>(RB);
    
    // Copy to all time slices
    for (unsigned int t = 0; t < T; ++t) {
      for(unsigned int r = 0; r < RB; ++r) {
        V(j, r, t) = draw(r);
      }
    }
  }
  
  return V;
}

// Sample G (static) given U, V and residual cube R
// Prior: vec(G) ~ N(0, (lambdaG * I)^-1)
// [[Rcpp::export]]
arma::mat sample_G_bip_cpp(const arma::cube& R,
                           const arma::cube& U,
                           const arma::cube& V,
                           const double lambdaG,
                           const arma::vec& s2) {
  const unsigned int T = R.n_slices;
  const unsigned int RA = U.n_cols;
  const unsigned int RB = V.n_cols;
  
  mat Prec(RA * RB, RA * RB, fill::zeros);
  vec rhs(RA * RB, fill::zeros);
  
  for (unsigned int t = 0; t < T; ++t) {
    double w = (s2.n_elem == 1) ? 1.0/s2(0) : 1.0/s2(t);
    mat Ut = U.slice(t);  // nA x RA
    mat Vt = V.slice(t);  // nB x RB
    mat At = Ut.t() * Ut;  // RA x RA
    mat Bt = Vt.t() * Vt;  // RB x RB
    
    // Kronecker product for vec(G) precision
    Prec += kron(Bt, At) * w;
    
    // Sufficient statistic
    mat St = Ut.t() * R.slice(t) * Vt;  // RA x RB
    rhs += vectorise(St) * w;
  }
  
  // Add prior precision
  Prec += lambdaG * eye(RA * RB, RA * RB);
  
  // Sample from posterior
  mat Sigma = inv_spd(Prec);
  vec mean = Sigma * rhs;
  vec draw = mean + chol(Sigma, "lower") * randn<vec>(RA * RB);
  
  // Reshape to matrix form
  return reshape(draw, RA, RB);
}

// Canonical orientation via SVD of G
// Ensures identifiability by setting G to rectangular diagonal
// [[Rcpp::export]]
List canon_orient_bip_cpp(const arma::cube& U,
                         const arma::cube& V,
                         const arma::mat& G) {
  // SVD of G
  mat Uc, Vc;
  vec s;
  svd(Uc, s, Vc, G);
  
  // Create rectangular diagonal matrix
  mat S(G.n_rows, G.n_cols, fill::zeros);
  unsigned int min_dim = std::min(G.n_rows, G.n_cols);
  for(unsigned int i = 0; i < min_dim; ++i) {
    S(i, i) = s(i);
  }
  
  // Transform U and V
  const unsigned int T = U.n_slices;
  arma::cube U_new = U;
  arma::cube V_new = V;
  
  for(unsigned int t = 0; t < T; ++t) {
    U_new.slice(t) = U.slice(t) * Uc;
    V_new.slice(t) = V.slice(t) * Vc;
  }
  
  return List::create(
    Named("U") = U_new,
    Named("V") = V_new,
    Named("G") = S
  );
}