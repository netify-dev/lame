// Bipartite beta/a/b Gibbs update in C++
// Replaces the nested R loops in lame.R for bipartite covariate updates
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
using namespace arma;
using namespace Rcpp;

// Compute X'X and X'y for bipartite covariate regression
//
// Replaces the O(T x p^2 x n^2) nested R loop for bipartite XtX and Xty
// computation with a single C++ call.
//
// @param Xlist List of T arrays, each nA x nB x p
// @param resid 3D array of residuals nA x nB x T
// @param p Number of covariates
// @return List with XtX (p x p) and Xty (p vector)
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
//' @param center Condition each actor margin to sum to zero.
//' @return List with beta, a, b
// [[Rcpp::export]]
List rbeta_ab_bip_gibbs_cpp(const arma::cube& Z,
                             Rcpp::List Xlist,
                             const arma::mat& UV_eff,
                             const arma::vec& a_current,
                             const arma::vec& b_current,
                             double s2, double g_prior,
                             double va, double vb,
                             bool rvar, bool cvar, bool center = false) {
  const int nA = Z.n_rows;
  const int nB = Z.n_cols;
  const int T = Z.n_slices;
  const int p = (Xlist.size() > 0) ? ((arma::cube)Xlist[0]).n_slices : 0;
  if (!std::isfinite(s2) || s2 <= 0.0 ||
      (p > 0 && (!std::isfinite(g_prior) || g_prior <= 0.0)) ||
      (rvar && (!std::isfinite(va) || va <= 0.0)) ||
      (cvar && (!std::isfinite(vb) || vb <= 0.0)))
    Rcpp::stop("Bipartite conditional variances and coefficient scale must be finite and positive.");

  arma::vec beta_new(p, fill::zeros);
  arma::vec a_new(nA, fill::zeros);
  arma::vec b_new(nB, fill::zeros);
  arma::vec variance_a(nA, fill::zeros), variance_b(nB, fill::zeros);

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
    // Forming a Gram matrix can lose the prior even if Cholesky succeeds.
    // Below sqrt(machine epsilon), use the augmented root before squaring
    // the design conditioning. This is a numerical route, not a prior floor.
    const double qr_threshold = std::sqrt(std::numeric_limits<double>::epsilon());
    const double reciprocal_condition = V_post_inv.is_finite() ? arma::rcond(V_post_inv) : 0.0;
    bool inv_ok = std::isfinite(reciprocal_condition) && reciprocal_condition > qr_threshold &&
      arma::inv_sympd(V_post, V_post_inv);
    arma::vec m_post;
    arma::mat chol_V;
    bool ordinary_ok = false;
    if (inv_ok && V_post.is_finite()) {
      m_post = V_post * (Xty / s2);
      V_post = 0.5 * (V_post + V_post.t());
      ordinary_ok = m_post.is_finite() && arma::chol(chol_V, V_post, "lower");
    }
    if (ordinary_ok) {
      // Preserve the original successful draw and its RNG sequence.
      arma::vec z(p, fill::randn);
      beta_new = m_post + chol_V * z;
    } else {
      // The positive prior can be lost when X'X is formed at a large scale.
      // QR of [X; I/sqrt(g)] retains its information before squaring the
      // condition number. Keep common sqrt(s2) outside the factorization:
      // scaling huge nearly equal X entries first loses their differences.
      // This is the SAME Gaussian conditional, including its original prior.
      const arma::uword cells = static_cast<arma::uword>(nA) * nB * T;
      arma::mat augmented(cells + p, p, fill::zeros);
      arma::vec response(cells + p, fill::zeros);
      const double prior_scale = 1.0 / std::sqrt(g_prior);
      if (!std::isfinite(prior_scale) || prior_scale <= 0.0)
        Rcpp::stop("Unable to represent the bipartite coefficient prior precision.");
      arma::uword row = 0;
      for (int t = 0; t < T; ++t) {
        arma::cube Xt = Xlist[t];
        for (int j = 0; j < nB; ++j) for (int i = 0; i < nA; ++i, ++row) {
          if (std::isfinite(resid_b(i,j,t))) response(row) = resid_b(i,j,t);
          for (int k = 0; k < p; ++k)
            if (std::isfinite(Xt(i,j,k))) augmented(row,k) = Xt(i,j,k);
        }
      }
      for (int k = 0; k < p; ++k) augmented(cells+k,k) = prior_scale;
      arma::mat q_factor, precision_root;
      if (!augmented.is_finite() || !response.is_finite() ||
          !arma::qr_econ(q_factor, precision_root, augmented) ||
          !precision_root.is_finite() || arma::any(arma::abs(precision_root.diag()) == 0.0))
        Rcpp::stop("Unable to resolve the bipartite coefficient Gaussian precision.");
      const double root_condition = arma::rcond(precision_root);
      // Beyond this range, likelihood roundoff can overwhelm prior rows.
      // Fail explicitly rather than silently drawing a different Gaussian.
      if (!std::isfinite(root_condition) ||
          root_condition <= 64.0 * std::numeric_limits<double>::epsilon())
        Rcpp::stop("Augmented bipartite Gaussian system is numerically unresolved.");
      // Expose small difference directions before QR. The SAME reversible
      // right map is applied to likelihood and prior rows, so no Gaussian
      // prior information is substituted. Accumulation precision is checked
      // below; platforms with double-sized long double refuse unresolved
      // systems instead of returning inaccurate finite conditional draws.
      arma::mat left_basis, right_basis;
      arma::vec singular;
      const double overall_scale = arma::abs(augmented).max();
      if (!(overall_scale > 0.0) || !std::isfinite(overall_scale) ||
          !arma::svd_econ(left_basis, singular, right_basis, augmented / overall_scale) ||
          !right_basis.is_finite())
        Rcpp::stop("Unable to resolve the augmented Gaussian column basis.");
      arma::mat transformed(augmented.n_rows, p, fill::zeros);
      arma::mat error_bound(augmented.n_rows, p, fill::zeros);
      const long double product_epsilon = std::numeric_limits<long double>::epsilon();
      // One fused product-plus-sum rounding per term gives gamma_p.
      const long double product_gamma = (static_cast<long double>(p) * product_epsilon) /
        (1.0L - static_cast<long double>(p) * product_epsilon);
      for (arma::uword i = 0; i < augmented.n_rows; ++i) for (int j = 0; j < p; ++j) {
        long double value = 0.0L, magnitude = 0.0L;
        for (int k = 0; k < p; ++k) {
          const long double term = static_cast<long double>(augmented(i,k)) *
            static_cast<long double>(right_basis(k,j));
          value = std::fma(static_cast<long double>(augmented(i,k)),
            static_cast<long double>(right_basis(k,j)), value);
          magnitude += std::abs(term);
        }
        transformed(i,j) = static_cast<double>(value);
        error_bound(i,j) = static_cast<double>(product_gamma * magnitude +
          std::numeric_limits<double>::epsilon() * std::abs(value));
      }
      arma::mat coordinate_map = right_basis;
      for (int j = 0; j < p; ++j) {
        const double column_scale = arma::abs(transformed.col(j)).max();
        if (!(column_scale > 0.0) || !std::isfinite(column_scale))
          Rcpp::stop("Unable to scale the augmented Gaussian column basis.");
        transformed.col(j) /= column_scale;
        error_bound.col(j) /= column_scale;
        error_bound.col(j) += 2.0 * std::numeric_limits<double>::epsilon() * arma::abs(transformed.col(j));
        coordinate_map.col(j) /= column_scale;
      }
      if (!transformed.is_finite() || !coordinate_map.is_finite() ||
          !arma::qr_econ(q_factor, precision_root, transformed) ||
          !arma::svd(singular, precision_root) || !singular.is_finite())
        Rcpp::stop("Unable to resolve the preconditioned Gaussian system.");
      const double perturbation = arma::norm(error_bound, "fro") / singular.min();
      if (!std::isfinite(perturbation) || perturbation >= 1.0 ||
          2.0 * perturbation / ((1.0 - perturbation) * (1.0 - perturbation)) > 1e-5)
        Rcpp::stop("Augmented bipartite Gaussian system is numerically unresolved.");
      arma::mat upper = arma::trimatu(precision_root);
      arma::mat covariance_root;
      if (!arma::solve(m_post, arma::trimatu(upper), q_factor.t() * response, arma::solve_opts::fast + arma::solve_opts::no_approx) ||
          !arma::solve(covariance_root, arma::trimatu(upper), arma::eye(p,p), arma::solve_opts::fast + arma::solve_opts::no_approx) ||
          !m_post.is_finite() || !covariance_root.is_finite())
        Rcpp::stop("Unable to resolve the bipartite coefficient Gaussian draw.");
      m_post = coordinate_map * m_post;
      covariance_root = std::sqrt(s2) * coordinate_map * covariance_root;
      if (!m_post.is_finite() || !covariance_root.is_finite())
        Rcpp::stop("Nonfinite mapped bipartite Gaussian conditional.");
      arma::vec z(p, fill::randn);
      beta_new = m_post + covariance_root * z;
    }
    if (!beta_new.is_finite()) Rcpp::stop("Nonfinite bipartite coefficient Gaussian draw.");
  }

  // Step 3: Compute residuals for a update (subtract beta*X and b)
  // the per-period covariate cubes are converted from the r list once,
  // not once per actor and period inside the loops below
  std::vector<arma::cube> Xcache;
  if(p > 0) {
    Xcache.resize(T);
    for(int t = 0; t < T; t++) Xcache[t] = Rcpp::as<arma::cube>(Xlist[t]);
  }
  if(rvar) {
    for(int i = 0; i < nA; i++) {
      double sum_resid = 0.0;
      int n_obs = 0;
      for(int t = 0; t < T; t++) {
        for(int j = 0; j < nB; j++) {
          double r = Z(i, j, t) - UV_eff(i, j);
          if(cvar) r -= b_current(j);
          if(p > 0) {
            const arma::cube& Xt = Xcache[t];
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
      if (!std::isfinite(prec_a) || prec_a <= 0.0 || !std::isfinite(mean_a))
        Rcpp::stop("Unable to resolve the bipartite row Gaussian conditional.");
      a_new(i) = R::rnorm(mean_a, sqrt(1.0 / prec_a));
      if (!std::isfinite(a_new(i))) Rcpp::stop("Nonfinite bipartite row Gaussian draw.");
      variance_a(i) = 1.0 / prec_a;
    }
  }

  if(center && rvar && accu(variance_a) > 0)
    a_new -= variance_a * (accu(a_new) / accu(variance_a));

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
            const arma::cube& Xt = Xcache[t];
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
      if (!std::isfinite(prec_b) || prec_b <= 0.0 || !std::isfinite(mean_b))
        Rcpp::stop("Unable to resolve the bipartite column Gaussian conditional.");
      b_new(j) = R::rnorm(mean_b, sqrt(1.0 / prec_b));
      if (!std::isfinite(b_new(j))) Rcpp::stop("Nonfinite bipartite column Gaussian draw.");
      variance_b(j) = 1.0 / prec_b;
    }
  }

  if(center && cvar && accu(variance_b) > 0)
    b_new -= variance_b * (accu(b_new) / accu(variance_b));

  return List::create(Named("beta") = beta_new,
                      Named("a") = a_new,
                      Named("b") = b_new);
}
