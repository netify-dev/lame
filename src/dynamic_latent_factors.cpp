// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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

  // ar(1) prior terms need s2 scaling to match the s2-parameterized
  // posterior: VtV_mod = V'V + s2*(ar terms), var = s2 * inv(VtV_mod)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_sigma2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2 * rho_uv * rho_uv * sigma2_inv;

  arma::cube& U_new = U_current;
  arma::cube& V_new = V_current;

  arma::mat VtV(R, R);
  arma::mat UtU(R, R);
  arma::mat cholDecomp(R, R);

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    VtV = V_t.t() * V_t;

    // Phase 1C: Batch compute all ei = V_t^T * E_t^T as a single BLAS call
    // E_all_U = E_t * V_t gives n x R, each row i = E_t.row(i) * V_t
    arma::mat E_all_U = E_t * V_t;  // n x R (one BLAS gemm instead of n gemv calls)

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();  // Extract row (O(R) not O(nR))
      arma::mat VtV_mod = VtV;

      if(t > 0) {
        VtV_mod.diag() += s2_sigma2_inv;
        ei += s2_rho_sigma2 * U_current.slice(t-1).row(i).t();
      }
      if(t < T-1) {
        VtV_mod.diag() += s2_rho2_sigma2;
        ei += s2_rho_sigma2 * U_current.slice(t+1).row(i).t();
      }

      // regularization prior: N(0, s2*I) => scaled precision = 1
      VtV_mod.diag() += 1.0;

      arma::mat iVtV = inv_sympd(VtV_mod);
      arma::vec mu = iVtV * ei;

      cholDecomp = chol(s2 * iVtV, "lower");
      arma::vec randVec(R, fill::randn);
      U_t.row(i) = (mu + cholDecomp * randVec).t();
    }

    if(!symmetric) {
      UtU = U_t.t() * U_t;

      // Phase 1C: Batch compute all ej = U_t^T * E_t as single BLAS call
      arma::mat E_all_V = E_t.t() * U_t;  // n x R

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();  // Extract row
        arma::mat UtU_mod = UtU;

        if(t > 0) {
          UtU_mod.diag() += s2_sigma2_inv;
          ej += s2_rho_sigma2 * V_current.slice(t-1).row(j).t();
        }
        if(t < T-1) {
          UtU_mod.diag() += s2_rho2_sigma2;
          ej += s2_rho_sigma2 * V_current.slice(t+1).row(j).t();
        }

        // regularization prior: N(0, s2*I) => scaled precision = 1
        UtU_mod.diag() += 1.0;

        arma::mat iUtU = inv_sympd(UtU_mod);
        arma::vec mu = iUtU * ej;

        cholDecomp = chol(s2 * iUtU, "lower");
        arma::vec randVec(R, fill::randn);
        V_t.row(j) = (mu + cholDecomp * randVec).t();
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
    }
  }

  if(shrink) {
    arma::vec s_vals;
    arma::mat U_svd, V_svd;

    // normalize using the time-averaged product to preserve temporal alignment
    arma::mat avg_prod(n, n, fill::zeros);
    for(int t = 0; t < T; t++) {
      avg_prod += U_new.slice(t) * V_new.slice(t).t();
    }
    avg_prod /= T;

    svd_econ(U_svd, s_vals, V_svd, avg_prod);
    const int rank = std::min(R, (int)s_vals.n_elem);

    if(rank > 0) {
      arma::mat S_sqrt = diagmat(sqrt(s_vals.head(rank)));
      // compute rotation from average to align all time slices consistently
      arma::mat U_ref = U_svd.head_cols(rank) * S_sqrt;
      arma::mat V_ref = V_svd.head_cols(rank) * S_sqrt;

      for(int t = 0; t < T; t++) {
        arma::mat& U_t = U_new.slice(t);
        arma::mat& V_t = V_new.slice(t);

        // procrustes rotation of U_t toward U_ref
        arma::mat M = U_ref.t() * U_t;
        arma::mat Mu, Mv;
        arma::vec Ms;
        svd_econ(Mu, Ms, Mv, M);
        arma::mat rot = Mv * Mu.t();

        U_t = U_t * rot;
        V_t = V_t * rot;
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

  positions.slice(0) = randn(n, R);

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

  double sum_yt_yt1 = 0.0;
  double sum_yt1_yt1 = 0.0;

  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);

    sum_yt_yt1 += accu(U_t % U_t1);
    sum_yt1_yt1 += accu(U_t1 % U_t1);

    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);

      sum_yt_yt1 += accu(V_t % V_t1);
      sum_yt1_yt1 += accu(V_t1 % V_t1);
    }
  }

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double var_post = 1.0 / (sum_yt1_yt1 * sigma2_inv + 1.0);
  const double mean_post = var_post * sum_yt_yt1 * sigma2_inv;

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

  double ss = 0.0;

  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);

    arma::mat innov_U = U_t - rho_uv * U_t1;
    ss += accu(innov_U % innov_U);

    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);

      arma::mat innov_V = V_t - rho_uv * V_t1;
      ss += accu(innov_V % innov_V);
    }
  }

  const int count = n * R * (T - 1) * (symmetric ? 1 : 2);

  const double shape = count / 2.0 + 1.0;
  const double scale = ss / 2.0 + 1.0;
  const double sigma2_new = 1.0 / R::rgamma(shape, 1.0/scale);

  return sqrt(sigma2_new);
}
