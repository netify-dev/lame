// Batch Z sampling across time periods to eliminate per-t R-to-C++ overhead
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;
using namespace Rcpp;

//' Batch normal Z sampling across all time periods
//'
//' Replaces the per-time-period R loop with a single C++ call that loops
//' internally, reducing R-to-C++ transition overhead.
//'
//' @param Z 3D array of latent values (n x n x T)
//' @param EZ 3D array of expected values (n x n x T)
//' @param rho Dyadic correlation parameter
//' @param s2 Dyadic variance
//' @param Y 3D array of observed values (n x n x T)
//' @return List with updated Z and E_nrm (residuals)
// [[Rcpp::export]]
List rZ_nrm_batch_cpp(arma::cube Z, const arma::cube& EZ,
                       double rho, double s2, const arma::cube& Y) {
  const int n = Z.n_rows;
  const int m = Z.n_cols;
  const int T = Z.n_slices;
  const double sz = sqrt(1.0 - rho * rho) * sqrt(s2);
  const double sqrts2 = sqrt(s2);

  arma::cube E_nrm(n, m, T);

  for(int t = 0; t < T; t++) {
    arma::mat& Z_t = Z.slice(t);
    const arma::mat& EZ_t = EZ.slice(t);
    const arma::mat& Y_t = Y.slice(t);

    // For square matrices with rho != 0, use dyadic correlation structure
    if(n == m && std::abs(rho) > 1e-10) {
      arma::mat Zt = Z_t.t();
      arma::mat EZt = EZ_t.t();

      for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
          if(i != j) {
            double ez_ij = EZ_t(i,j) + rho * (Zt(i,j) - EZt(i,j));
            if(std::isnan(Y_t(i,j))) {
              Z_t(i,j) = R::rnorm(ez_ij, sz);
            } else {
              Z_t(i,j) = Y_t(i,j);
            }
          }
        }
      }
      // Diagonal
      for(int i = 0; i < n; i++) {
        Z_t(i,i) = std::isnan(Y_t(i,i)) ?
          R::rnorm(EZ_t(i,i), sqrt((1.0 + rho) * s2)) : Y_t(i,i);
      }
    } else {
      // Bipartite or rho=0: simple normal sampling
      for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
          if(std::isnan(Y_t(i,j))) {
            Z_t(i,j) = R::rnorm(EZ_t(i,j), sqrts2);
          } else {
            Z_t(i,j) = Y_t(i,j);
          }
        }
      }
    }

    E_nrm.slice(t) = Z_t - EZ_t;
  }

  return List::create(Named("Z") = Z, Named("E_nrm") = E_nrm);
}

//' Batch binary Z sampling across all time periods (bipartite, rho=0)
//'
//' Vectorized probit update for bipartite binary networks without dyadic correlation.
//'
//' @param Z 3D array of latent values (nA x nB x T)
//' @param EZ 3D array of expected values (nA x nB x T)
//' @param Y 3D array of observed values (nA x nB x T)
//' @return Updated Z array
// [[Rcpp::export]]
arma::cube rZ_bin_bip_batch_cpp(arma::cube Z, const arma::cube& EZ,
                                 const arma::cube& Y) {
  const int nA = Z.n_rows;
  const int nB = Z.n_cols;
  const int T = Z.n_slices;

  for(int t = 0; t < T; t++) {
    arma::mat& Z_t = Z.slice(t);
    const arma::mat& EZ_t = EZ.slice(t);
    const arma::mat& Y_t = Y.slice(t);

    for(int i = 0; i < nA; i++) {
      for(int j = 0; j < nB; j++) {
        double ez = EZ_t(i,j);
        if(std::isnan(Y_t(i,j))) {
          Z_t(i,j) = R::rnorm(ez, 1.0);
        } else if(Y_t(i,j) > 0.5) {
          // Y=1: truncated normal from below at 0
          double p_lo = R::pnorm(-ez, 0.0, 1.0, 1, 0);
          double u = p_lo + R::runif(0.0, 1.0) * (1.0 - p_lo);
          Z_t(i,j) = ez + R::qnorm(u, 0.0, 1.0, 1, 0);
        } else {
          // Y=0: truncated normal from above at 0
          double p_hi = R::pnorm(-ez, 0.0, 1.0, 1, 0);
          double u = R::runif(0.0, 1.0) * p_hi;
          Z_t(i,j) = ez + R::qnorm(u, 0.0, 1.0, 1, 0);
        }
      }
    }
  }

  return Z;
}

//' Batch tobit Z sampling across all time periods (bipartite, rho=0)
//'
//' @param Z 3D array of latent values (nA x nB x T)
//' @param EZ 3D array of expected values (nA x nB x T)
//' @param s2 Dyadic variance
//' @param Y 3D array of observed values (nA x nB x T)
//' @return Updated Z array
// [[Rcpp::export]]
arma::cube rZ_tob_bip_batch_cpp(arma::cube Z, const arma::cube& EZ,
                                 double s2, const arma::cube& Y) {
  const int nA = Z.n_rows;
  const int nB = Z.n_cols;
  const int T = Z.n_slices;
  const double sqrts2 = sqrt(s2);

  for(int t = 0; t < T; t++) {
    arma::mat& Z_t = Z.slice(t);
    const arma::mat& EZ_t = EZ.slice(t);
    const arma::mat& Y_t = Y.slice(t);

    for(int i = 0; i < nA; i++) {
      for(int j = 0; j < nB; j++) {
        double ez = EZ_t(i,j);
        if(std::isnan(Y_t(i,j))) {
          Z_t(i,j) = R::rnorm(ez, sqrts2);
        } else if(Y_t(i,j) > 0) {
          Z_t(i,j) = Y_t(i,j);
        } else {
          // Y=0: truncated normal from above at 0
          double z_prop = R::rnorm(ez, sqrts2);
          Z_t(i,j) = std::min(z_prop, 0.0);
        }
      }
    }
  }

  return Z;
}
