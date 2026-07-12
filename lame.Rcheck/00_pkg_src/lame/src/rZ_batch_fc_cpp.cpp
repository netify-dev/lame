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

    // For square matrices with rho != 0, use dyadic correlation structure.
    // dyad-pair block update: observed cells are fixed at y; a missing
    // cell whose reciprocal is observed is drawn from the exact
    // conditional N(EZ_ij + rho*(y_ji - EZ_ji), s2*(1 - rho^2)); pairs
    // with both cells missing get a joint bivariate draw with
    // within-pair correlation rho. (the earlier snapshot-based scheme
    // updated both cells of a both-missing pair synchronously from the
    // stale transpose, driving their stationary within-pair correlation
    // to 0 instead of rho.)
    if(n == m && std::abs(rho) > 1e-10) {
      const double cc = (std::sqrt(1.0 + rho) + std::sqrt(1.0 - rho)) / 2.0;
      const double dd = (std::sqrt(1.0 + rho) - std::sqrt(1.0 - rho)) / 2.0;

      for(int j = 1; j < m; j++) {
        for(int i = 0; i < j; i++) {
          const bool mis_ij = std::isnan(Y_t(i,j));
          const bool mis_ji = std::isnan(Y_t(j,i));
          if(!mis_ij) Z_t(i,j) = Y_t(i,j);
          if(!mis_ji) Z_t(j,i) = Y_t(j,i);
          if(mis_ij && mis_ji) {
            const double x1 = norm_rand();
            const double x2 = norm_rand();
            Z_t(i,j) = EZ_t(i,j) + sqrts2 * (cc * x1 + dd * x2);
            Z_t(j,i) = EZ_t(j,i) + sqrts2 * (cc * x2 + dd * x1);
          } else if(mis_ij) {
            Z_t(i,j) = R::rnorm(EZ_t(i,j) + rho * (Y_t(j,i) - EZ_t(j,i)), sz);
          } else if(mis_ji) {
            Z_t(j,i) = R::rnorm(EZ_t(j,i) + rho * (Y_t(i,j) - EZ_t(i,j)), sz);
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
        if(!std::isfinite(ez)) { Z_t(i,j) = 0.0; continue; }
        if(std::isnan(Y_t(i,j))) {
          Z_t(i,j) = R::rnorm(ez, 1.0);
        } else if(Y_t(i,j) > 0.5) {
          // Y=1: truncated normal on [0, inf). exact log.p-scale
          // inverse cdf: tail-stable for arbitrary |ez| with no
          // probability clamp (the earlier 1e-12 clamp produced
          // sign-violating draws once ez > ~7.03 and capped the
          // likelihood's restoring force)
          double logS = R::pnorm(-ez, 0.0, 1.0, 0, 1);
          Z_t(i,j) = ez + R::qnorm(std::log(unif_rand()) + logS,
                                   0.0, 1.0, 0, 1);
        } else {
          // Y=0: truncated normal on (-inf, 0]
          double logF = R::pnorm(-ez, 0.0, 1.0, 1, 1);
          Z_t(i,j) = ez + R::qnorm(std::log(unif_rand()) + logF,
                                   0.0, 1.0, 1, 1);
        }
        // backstop for a non-finite draw (requires unif_rand() to hit
        // an exact endpoint, which r's rng never does): fall back to
        // the truncation boundary for observed cells (sign-consistent
        // for either y), the conditional mean for missing cells
        if(!std::isfinite(Z_t(i,j))) {
          Z_t(i,j) = std::isnan(Y_t(i,j)) ? ez : 0.0;
        }
      }
    }
  }

  return Z;
}
