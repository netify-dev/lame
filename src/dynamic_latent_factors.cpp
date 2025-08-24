// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Update dynamic latent positions using AR(1) process
//' 
//' @param U_current Current 3D array of U positions (n x R x T)
//' @param V_current Current 3D array of V positions (n x R x T) 
//' @param ET 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) autoregressive parameter for latent positions
//' @param sigma_uv Innovation standard deviation for latent positions
//' @param s2 Dyadic variance
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @return List with updated U and V arrays
// [[Rcpp::export]]
List rUV_dynamic_fc_cpp(arma::cube U_current, arma::cube V_current,
                        const arma::cube& ET, double rho_uv, 
                        double sigma_uv, double s2, bool shrink,
                        bool symmetric) {
  
  const int n = U_current.n_rows;
  const int R = U_current.n_cols; 
  const int T = U_current.n_slices;
  
  // Pre-compute constants
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double rho_sigma2 = rho_uv * sigma2_inv;
  const double rho2_sigma2 = rho_uv * rho_uv * sigma2_inv;
  const double s2_inv = 1.0 / s2;
  
  // Work with references to avoid copies
  arma::cube& U_new = U_current;
  arma::cube& V_new = V_current;
  
  // Pre-allocate working matrices
  arma::mat VtV(R, R);
  arma::mat UtU(R, R);
  arma::mat iVtV(R, R);
  arma::mat iUtU(R, R);
  arma::vec ei(R);  // Changed from n to R - stores latent dimension effects
  arma::vec ej(R);  // Changed from n to R - stores latent dimension effects
  arma::vec mu(R);
  arma::vec u_new(R);
  arma::vec v_new(R);
  arma::mat cholDecomp(R, R);
  arma::vec randVec(R);
  
  // Update latent positions for each time point
  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);
    
    // Pre-compute V^T * V for this time slice
    VtV = V_t.t() * V_t;
    
    // Gibbs update for U given V at time t
    for(int i = 0; i < n; i++) {
      // Compute the contribution from the data: V' * E[i,:]
      ei = V_t.t() * E_t.row(i).t();
      
      // Copy VtV for modification
      arma::mat VtV_mod = VtV;
      
      // Add AR(1) prior contributions
      if(t > 0) {
        VtV_mod.diag() += sigma2_inv;
        ei += rho_sigma2 * U_current.slice(t-1).row(i).t();
      }
      if(t < T-1) {
        VtV_mod.diag() += rho2_sigma2;
        ei += rho_sigma2 * U_current.slice(t+1).row(i).t();
      }
      
      // Add s2 regularization
      VtV_mod.diag() += s2_inv;
      
      // Compute posterior parameters using Cholesky
      iVtV = inv_sympd(VtV_mod);
      mu = iVtV * ei;
      
      // Sample from posterior using Cholesky decomposition
      cholDecomp = chol(s2 * iVtV, "lower");
      randVec.randn();
      U_t.row(i) = (mu + cholDecomp * randVec).t();
    }
    
    // Gibbs update for V given U at time t (if not symmetric)
    if(!symmetric) {
      // Pre-compute U^T * U for this time slice
      UtU = U_t.t() * U_t;
      
      for(int j = 0; j < n; j++) {
        // Compute the contribution from the data: U' * E[:,j]
        ej = U_t.t() * E_t.col(j);
        
        // Copy UtU for modification
        arma::mat UtU_mod = UtU;
        
        // Add AR(1) prior contributions
        if(t > 0) {
          UtU_mod.diag() += sigma2_inv;
          ej += rho_sigma2 * V_current.slice(t-1).row(j).t();
        }
        if(t < T-1) {
          UtU_mod.diag() += rho2_sigma2;
          ej += rho_sigma2 * V_current.slice(t+1).row(j).t();
        }
        
        // Add s2 regularization
        UtU_mod.diag() += s2_inv;
        
        // Compute posterior parameters using Cholesky
        iUtU = inv_sympd(UtU_mod);
        mu = iUtU * ej;
        
        // Sample from posterior
        cholDecomp = chol(s2 * iUtU, "lower");
        randVec.randn();
        V_t.row(j) = (mu + cholDecomp * randVec).t();
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
    }
  }
  
  // Apply shrinkage if requested
  if(shrink) {
    arma::vec s_vals;
    arma::mat U_svd, V_svd;
    
    for(int t = 0; t < T; t++) {
      arma::mat& U_t = U_new.slice(t);
      arma::mat& V_t = V_new.slice(t);
      
      // SVD for shrinkage - use economical SVD
      svd_econ(U_svd, s_vals, V_svd, U_t * V_t.t());
      
      // Reconstruct with shrinkage
      const int rank = std::min(R, (int)s_vals.n_elem);
      if(rank > 0) {
        arma::mat S_sqrt = diagmat(sqrt(s_vals.head(rank)));
        U_t = U_svd.head_cols(rank) * S_sqrt;
        V_t = V_svd.head_cols(rank) * S_sqrt;
      }
    }
  }
  
  return List::create(Named("U") = U_new, Named("V") = V_new);
}

//' Initialize dynamic latent positions with AR(1) structure
//' 
//' @param n Number of actors
//' @param R Latent dimension
//' @param T Number of time points
//' @param rho_uv AR(1) parameter
//' @param sigma_uv Innovation standard deviation
//' @return 3D array of latent positions (n x R x T)
// [[Rcpp::export]]
arma::cube init_dynamic_positions(int n, int R, int T, 
                                  double rho_uv, double sigma_uv) {
  arma::cube positions(n, R, T);
  
  // Initialize first time point from N(0,1)
  positions.slice(0) = randn(n, R);
  
  // Generate subsequent time points using AR(1)
  for(int t = 1; t < T; t++) {
    for(int i = 0; i < n; i++) {
      for(int r = 0; r < R; r++) {
        double innov = R::rnorm(0.0, sigma_uv);
        positions(i, r, t) = rho_uv * positions(i, r, t-1) + innov;
      }
    }
  }
  
  return positions;
}

//' Sample AR(1) parameter for dynamic latent factors
//' 
//' @param U_cube 3D array of U positions (n x R x T)
//' @param V_cube 3D array of V positions (n x R x T)
//' @param sigma_uv Innovation standard deviation
//' @param rho_current Current value of rho
//' @param symmetric Whether network is symmetric
//' @return Updated rho value
// [[Rcpp::export]]
double sample_rho_uv(const arma::cube& U_cube, const arma::cube& V_cube,
                     double sigma_uv, double rho_current, bool symmetric) {
  
  const int n = U_cube.n_rows;
  const int R = U_cube.n_cols;
  const int T = U_cube.n_slices;
  
  // Calculate sufficient statistics using vectorized operations
  double sum_yt_yt1 = 0.0;
  double sum_yt1_yt1 = 0.0;
  
  // Process each time slice
  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);
    
    // Use element-wise operations for efficiency
    sum_yt_yt1 += accu(U_t % U_t1);
    sum_yt1_yt1 += accu(U_t1 % U_t1);
    
    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);
      
      sum_yt_yt1 += accu(V_t % V_t1);
      sum_yt1_yt1 += accu(V_t1 % V_t1);
    }
  }
  
  // Posterior parameters (conjugate normal)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double var_post = 1.0 / (sum_yt1_yt1 * sigma2_inv + 1.0);
  const double mean_post = var_post * sum_yt_yt1 * sigma2_inv;
  
  // Sample and truncate to ensure stationarity
  double rho_new = R::rnorm(mean_post, sqrt(var_post));
  return std::max(-0.99, std::min(0.99, rho_new));
}

//' Sample innovation variance for dynamic latent factors
//' 
//' @param U_cube 3D array of U positions (n x R x T)
//' @param V_cube 3D array of V positions (n x R x T)
//' @param rho_uv AR(1) parameter
//' @param symmetric Whether network is symmetric
//' @return Updated sigma_uv value
// [[Rcpp::export]]
double sample_sigma_uv(const arma::cube& U_cube, const arma::cube& V_cube,
                       double rho_uv, bool symmetric) {
  
  const int n = U_cube.n_rows;
  const int R = U_cube.n_cols;
  const int T = U_cube.n_slices;
  
  // Calculate sum of squared innovations using vectorized operations
  double ss = 0.0;
  
  // Process each time slice
  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);
    
    // Compute innovations and sum squares
    arma::mat innov_U = U_t - rho_uv * U_t1;
    ss += accu(innov_U % innov_U);
    
    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);
      
      arma::mat innov_V = V_t - rho_uv * V_t1;
      ss += accu(innov_V % innov_V);
    }
  }
  
  // Count total number of parameters
  const int count = n * R * (T - 1) * (symmetric ? 1 : 2);
  
  // Sample from inverse gamma posterior
  const double shape = count / 2.0 + 1.0;
  const double scale = ss / 2.0 + 1.0;
  const double sigma2_new = 1.0 / R::rgamma(shape, 1.0/scale);
  
  return sqrt(sigma2_new);
}