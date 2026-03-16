// Bipartite beta/a/b Gibbs update in C++
// Replaces the nested R loops in lame.R for bipartite covariate updates
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Compute X'X and X'y for bipartite covariate regression
//'
//' Replaces the O(T x p^2 x n^2) nested R loop for bipartite XtX and Xty
//' computation with a single C++ call.
//'
//' @param Xlist List of T arrays, each nA x nB x p
//' @param resid 3D array of residuals nA x nB x T
//' @param p Number of covariates
//' @return List with XtX (p x p) and Xty (p vector)
// [[Rcpp::export]]
List compute_XtX_Xty_bip_cpp(Rcpp::List Xlist, const arma::cube& resid, int p) {
  const int T = Xlist.size();

  arma::mat XtX(p, p, fill::zeros);
  arma::vec Xty(p, fill::zeros);

  for(int t = 0; t < T; t++) {
    arma::cube Xt = Xlist[t];
    const arma::mat& resid_t = resid.slice(t);

    for(int k1 = 0; k1 < p; k1++) {
      const arma::mat& Xk1 = Xt.slice(k1);

      // X'y: sum of element-wise product X[,,k1] * resid
      // Only count non-NA elements
      double xty_val = 0.0;
      for(unsigned int i = 0; i < Xk1.n_rows; i++) {
        for(unsigned int j = 0; j < Xk1.n_cols; j++) {
          if(std::isfinite(resid_t(i,j)) && std::isfinite(Xk1(i,j))) {
            xty_val += Xk1(i,j) * resid_t(i,j);
          }
        }
      }
      Xty(k1) += xty_val;

      // X'X: sum of element-wise product X[,,k1] * X[,,k2]
      for(int k2 = k1; k2 < p; k2++) {
        const arma::mat& Xk2 = Xt.slice(k2);
        double val = 0.0;
        for(unsigned int i = 0; i < Xk1.n_rows; i++) {
          for(unsigned int j = 0; j < Xk1.n_cols; j++) {
            if(std::isfinite(Xk1(i,j)) && std::isfinite(Xk2(i,j))) {
              val += Xk1(i,j) * Xk2(i,j);
            }
          }
        }
        XtX(k1, k2) += val;
        if(k1 != k2) XtX(k2, k1) += val;
      }
    }
  }

  return List::create(Named("XtX") = XtX, Named("Xty") = Xty);
}

//' Compute Xbeta product for bipartite networks
//'
//' Computes sum_k beta_k * X_k for a single time slice
//'
//' @param X 3D array (nA x nB x p) of covariates for one time period
//' @param beta Coefficient vector of length p
//' @return nA x nB matrix
// [[Rcpp::export]]
arma::mat Xbeta_bip_cpp(const arma::cube& X, const arma::vec& beta) {
  const int nA = X.n_rows;
  const int nB = X.n_cols;
  const int p = beta.n_elem;

  arma::mat result(nA, nB, fill::zeros);
  for(int k = 0; k < p; k++) {
    if(beta(k) != 0.0) {
      result += beta(k) * X.slice(k);
    }
  }
  return result;
}

//' Full bipartite Gibbs update for beta, a, b
//'
//' Single C++ function replacing the bipartite beta/a/b update block in lame.R
//'
//' @param Z 3D array (nA x nB x T)
//' @param Xlist List of T arrays (nA x nB x p)
//' @param UV_eff nA x nB matrix (U*G*V' or U*V')
//' @param a_current Current row effects (length nA)
//' @param b_current Current column effects (length nB)
//' @param s2 Dyadic variance
//' @param g_prior G-prior parameter
//' @param va Row effect variance (diagonal element of Sab)
//' @param vb Column effect variance (diagonal element of Sab)
//' @param rvar Whether to update row effects
//' @param cvar Whether to update column effects
//' @return List with beta, a, b
// [[Rcpp::export]]
List rbeta_ab_bip_gibbs_cpp(const arma::cube& Z,
                             Rcpp::List Xlist,
                             const arma::mat& UV_eff,
                             const arma::vec& a_current,
                             const arma::vec& b_current,
                             double s2, double g_prior,
                             double va, double vb,
                             bool rvar, bool cvar) {
  const int nA = Z.n_rows;
  const int nB = Z.n_cols;
  const int T = Z.n_slices;
  const int p = (Xlist.size() > 0) ? ((arma::cube)Xlist[0]).n_slices : 0;

  arma::vec beta_new(p, fill::zeros);
  arma::vec a_new(nA, fill::zeros);
  arma::vec b_new(nB, fill::zeros);

  // Step 1: Compute residuals = Z - UV_eff - a - b
  arma::cube resid_b(nA, nB, T);
  for(int t = 0; t < T; t++) {
    arma::mat r = Z.slice(t) - UV_eff;
    if(rvar) r.each_col() -= a_current;
    if(cvar) r.each_row() -= b_current.t();
    resid_b.slice(t) = r;
  }

  // Step 2: Update beta
  if(p > 0) {
    List xtx_xty = compute_XtX_Xty_bip_cpp(Xlist, resid_b, p);
    arma::mat XtX = xtx_xty["XtX"];
    arma::vec Xty = xtx_xty["Xty"];

    arma::mat V_post_inv = XtX / s2 + arma::eye(p, p) / (g_prior * s2);
    arma::mat V_post;
    bool inv_ok = arma::inv_sympd(V_post, V_post_inv);
    if(!inv_ok) {
      V_post = arma::diagmat(arma::vec(p, fill::value(s2)));
    }
    arma::vec m_post = V_post * (Xty / s2);

    // Sample from posterior
    V_post = 0.5 * (V_post + V_post.t());
    arma::mat chol_V;
    bool chol_ok = arma::chol(chol_V, V_post, "lower");
    if(chol_ok) {
      arma::vec z(p, fill::randn);
      beta_new = m_post + chol_V * z;
    } else {
      beta_new = m_post;
    }
  }

  // Step 3: Compute residuals for a update (subtract beta*X and b)
  if(rvar) {
    for(int i = 0; i < nA; i++) {
      double sum_resid = 0.0;
      int n_obs = 0;
      for(int t = 0; t < T; t++) {
        for(int j = 0; j < nB; j++) {
          double r = Z(i, j, t) - UV_eff(i, j);
          if(cvar) r -= b_current(j);
          if(p > 0) {
            arma::cube Xt = Xlist[t];
            for(int k = 0; k < p; k++) {
              r -= beta_new(k) * Xt(i, j, k);
            }
          }
          if(std::isfinite(r)) {
            sum_resid += r;
            n_obs++;
          }
        }
      }
      double prec_a = n_obs / s2 + 1.0 / va;
      double mean_a = (sum_resid / s2) / prec_a;
      a_new(i) = R::rnorm(mean_a, sqrt(1.0 / prec_a));
    }
  }

  // Step 4: Compute residuals for b update (subtract beta*X and a_new)
  if(cvar) {
    for(int j = 0; j < nB; j++) {
      double sum_resid = 0.0;
      int n_obs = 0;
      for(int t = 0; t < T; t++) {
        for(int i = 0; i < nA; i++) {
          double r = Z(i, j, t) - UV_eff(i, j);
          if(rvar) r -= a_new(i);
          if(p > 0) {
            arma::cube Xt = Xlist[t];
            for(int k = 0; k < p; k++) {
              r -= beta_new(k) * Xt(i, j, k);
            }
          }
          if(std::isfinite(r)) {
            sum_resid += r;
            n_obs++;
          }
        }
      }
      double prec_b = n_obs / s2 + 1.0 / vb;
      double mean_b = (sum_resid / s2) / prec_b;
      b_new(j) = R::rnorm(mean_b, sqrt(1.0 / prec_b));
    }
  }

  return List::create(Named("beta") = beta_new,
                      Named("a") = a_new,
                      Named("b") = b_new);
}
