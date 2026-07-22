// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Sample dynamic additive effects with AR(1) evolution
//'
//' Gibbs update of the per-period additive effects under the AR(1) state
//' model with stationary initial condition. For each actor and period the
//' full conditional combines the AR(1) bridge prior (stationary init at
//' t = 1) with the dyadic residual likelihood: resid_ij = a_i + b_j + e_ij,
//' e_ij ~ N(0, s2). The reciprocal effect (b_j for the a-step, the freshly
//' updated a_i for the b-step) is subtracted from each residual, missing
//' residuals are skipped, and for symmetric networks each dyad contributes
//' exactly once (resid(i,j) = a_i + a_j + e_ij).
//'
//' @param a_current Current 2D array of row effects (n x T)
//' @param b_current Current 2D array of column effects (n x T)
//' @param Z_array 3D array of latent positions (n x n x T)
//' @param EZ_array 3D array of expected values without additive effects (n x n x T)
//' @param rho_ab AR(1) parameter for additive effects
//' @param sigma_ab Innovation standard deviation
//' @param s2 Dyadic residual variance
//' @param symmetric Whether the network is symmetric
//' @return List with updated a and b arrays
// [[Rcpp::export]]
List sample_dynamic_ab_cpp(arma::mat a_current, arma::mat b_current,
                           const arma::cube& Z_array, const arma::cube& EZ_array,
                           double rho_ab, double sigma_ab,
                           double s2, bool symmetric) {

  const int n = a_current.n_rows;
  const int T = a_current.n_cols;

  arma::mat a_new = a_current;
  arma::mat b_new = b_current;

  const double var_innov = std::max(sigma_ab * sigma_ab, 1e-8);
  const double s2_safe = std::max(s2, 1e-8);
  const double rho2 = std::min(rho_ab * rho_ab, 0.9801);

  // t-outer sweep with in-place updates so each conditional sees the
  // freshest values of the other effects (proper sequential-scan Gibbs)
  for(int t = 0; t < T; t++) {

    arma::mat resid_t = Z_array.slice(t) - EZ_array.slice(t);

    // a-step (symmetric folds b into a)
    for(int i = 0; i < n; i++) {

      // ar(1) bridge prior with stationary initial condition
      double prior_prec, prior_mean;
      if(T == 1) {
        prior_prec = (1.0 - rho2) / var_innov;
        prior_mean = 0.0;
      } else if(t == 0) {
        prior_prec = 1.0 / var_innov;
        prior_mean = rho_ab * a_new(i, 1);
      } else if(t == T - 1) {
        prior_prec = 1.0 / var_innov;
        prior_mean = rho_ab * a_new(i, T - 2);
      } else {
        prior_prec = (1.0 + rho2) / var_innov;
        prior_mean = rho_ab * (a_new(i, t - 1) + a_new(i, t + 1)) / (1.0 + rho2);
      }

      double dsum = 0.0;
      int nobs = 0;
      for(int j = 0; j < n; j++) {
        if(j == i) continue;
        double r = resid_t(i, j);
        if(!std::isfinite(r)) continue;
        // subtract the reciprocal effect: b_j (directed) or a_j (symmetric,
        // where each dyad appears exactly once in row i)
        double other = symmetric ? a_new(j, t) : b_new(j, t);
        if(!std::isfinite(other)) other = 0.0;
        dsum += r - other;
        nobs++;
      }

      double post_prec = prior_prec + nobs / s2_safe;
      double post_mean = (prior_prec * prior_mean + dsum / s2_safe) / post_prec;

      if(std::isfinite(post_prec) && post_prec > 0.0 && std::isfinite(post_mean)) {
        a_new(i, t) = R::rnorm(post_mean, 1.0 / sqrt(post_prec));
      } else {
        a_new(i, t) = std::isfinite(prior_mean) ? prior_mean : 0.0;
      }

      if(symmetric) b_new(i, t) = a_new(i, t);
    }

    // b-step: uses the freshly updated a values for all rows
    if(!symmetric) {
      for(int j = 0; j < n; j++) {

        double prior_prec, prior_mean;
        if(T == 1) {
          prior_prec = (1.0 - rho2) / var_innov;
          prior_mean = 0.0;
        } else if(t == 0) {
          prior_prec = 1.0 / var_innov;
          prior_mean = rho_ab * b_new(j, 1);
        } else if(t == T - 1) {
          prior_prec = 1.0 / var_innov;
          prior_mean = rho_ab * b_new(j, T - 2);
        } else {
          prior_prec = (1.0 + rho2) / var_innov;
          prior_mean = rho_ab * (b_new(j, t - 1) + b_new(j, t + 1)) / (1.0 + rho2);
        }

        double dsum = 0.0;
        int nobs = 0;
        for(int i = 0; i < n; i++) {
          if(i == j) continue;
          double r = resid_t(i, j);
          if(!std::isfinite(r)) continue;
          double ai = a_new(i, t);
          if(!std::isfinite(ai)) ai = 0.0;
          dsum += r - ai;
          nobs++;
        }

        double post_prec = prior_prec + nobs / s2_safe;
        double post_mean = (prior_prec * prior_mean + dsum / s2_safe) / post_prec;

        if(std::isfinite(post_prec) && post_prec > 0.0 && std::isfinite(post_mean)) {
          b_new(j, t) = R::rnorm(post_mean, 1.0 / sqrt(post_prec));
        } else {
          b_new(j, t) = std::isfinite(prior_mean) ? prior_mean : 0.0;
        }
      }
    }
  }

  return List::create(
    Named("a") = a_new,
    Named("b") = b_new
  );
}

//' Sample AR(1) parameter for dynamic additive effects
//'
//' @param a_mat Matrix of row effects (n x T)
//' @param b_mat Matrix of column effects (n x T)
//' @param sigma_ab Innovation standard deviation
//' @param rho_current Current value of rho
//' @param symmetric Whether the network is symmetric
//' @param prior_mean Prior mean for rho. Used only when \code{prior_sd > 0}.
//' @param prior_sd Prior SD for rho. \code{prior_sd < 0} (the default)
//'   selects a Jeffreys-like flat prior; a positive value switches to a
//'   truncated Normal(prior_mean, prior_sd^2) prior.
//' @return Updated rho value
// [[Rcpp::export]]
double sample_rho_ab_cpp(const arma::mat& a_mat, const arma::mat& b_mat,
                         double sigma_ab, double rho_current, bool symmetric,
                         double prior_mean = 0.0, double prior_sd = -1.0) {
  // prior_sd < 0 (the default) selects the historical Jeffreys-like prior
  // -0.5 * log(1 - rho^2); a positive prior_sd switches to a truncated
  // Normal(prior_mean, prior_sd^2) prior on rho, as documented in
  // ?lame for prior$rho_ab_mean / prior$rho_ab_sd.

  const int n = a_mat.n_rows;
  const int T = a_mat.n_cols;

  if(T <= 1) return rho_current;

  double sum_prod = 0.0;
  double sum_sq_lag = 0.0;
  double sum_sq_curr = 0.0;
  int count = 0;

  // stationary-init sufficient statistics (t = 1 states); the state sampler
  // uses the stationary AR(1) initial condition, so its (1 - rho^2) factor
  // must appear in this conditional too for a coherent joint
  double q1 = 0.0;
  int p1 = 0;

  for(int i = 0; i < n; i++) {
    q1 += a_mat(i,0) * a_mat(i,0);
    p1++;
    if(!symmetric) {
      q1 += b_mat(i,0) * b_mat(i,0);
      p1++;
    }
    for(int t = 1; t < T; t++) {
      sum_prod += a_mat(i,t) * a_mat(i,t-1);
      sum_sq_lag += a_mat(i,t-1) * a_mat(i,t-1);
      sum_sq_curr += a_mat(i,t) * a_mat(i,t);
      count++;

      if(!symmetric) {
        sum_prod += b_mat(i,t) * b_mat(i,t-1);
        sum_sq_lag += b_mat(i,t-1) * b_mat(i,t-1);
        sum_sq_curr += b_mat(i,t) * b_mat(i,t);
        count++;
      }
    }
  }

  double proposal_sd = 0.1;
  double rho_proposal = R::rnorm(rho_current, proposal_sd);

  if(std::abs(rho_proposal) >= 1.0) {
    return rho_current;
  }

  const double sigma2 = sigma_ab * sigma_ab;

  double ll_current = -0.5 * count * log(2.0 * M_PI * sigma2);
  ll_current -= 0.5 * (sum_sq_curr - 2.0 * rho_current * sum_prod +
                       rho_current * rho_current * sum_sq_lag) / sigma2;
  // stationary init factor: (1 - rho^2)^{p1/2} exp(-(1 - rho^2) q1 / (2 sigma^2))
  ll_current += 0.5 * p1 * log(1.0 - rho_current * rho_current);
  ll_current -= 0.5 * (1.0 - rho_current * rho_current) * q1 / sigma2;

  double ll_proposal = -0.5 * count * log(2.0 * M_PI * sigma2);
  ll_proposal -= 0.5 * (sum_sq_curr - 2.0 * rho_proposal * sum_prod +
                        rho_proposal * rho_proposal * sum_sq_lag) / sigma2;
  ll_proposal += 0.5 * p1 * log(1.0 - rho_proposal * rho_proposal);
  ll_proposal -= 0.5 * (1.0 - rho_proposal * rho_proposal) * q1 / sigma2;

  double log_prior_current, log_prior_proposal;
  if (prior_sd > 0.0) {
    // truncated Normal(prior_mean, prior_sd^2); truncation constant cancels
    const double inv2v = 0.5 / (prior_sd * prior_sd);
    log_prior_current  = -inv2v * (rho_current  - prior_mean) * (rho_current  - prior_mean);
    log_prior_proposal = -inv2v * (rho_proposal - prior_mean) * (rho_proposal - prior_mean);
  } else {
    // Jeffreys-like flat prior on rho
    log_prior_current  = -0.5 * log(1.0 - rho_current * rho_current);
    log_prior_proposal = -0.5 * log(1.0 - rho_proposal * rho_proposal);
  }

  double log_ratio = (ll_proposal + log_prior_proposal) - (ll_current + log_prior_current);

  if(log(R::runif(0,1)) < log_ratio) {
    return rho_proposal;
  } else {
    return rho_current;
  }
}

//' Sample innovation variance for dynamic additive effects
//'
//' @param a_mat Matrix of row effects (n x T)
//' @param b_mat Matrix of column effects (n x T)
//' @param rho_ab AR(1) parameter
//' @param symmetric Whether the network is symmetric
//' @param prior_shape Shape parameter for inverse gamma prior
//' @param prior_scale Scale parameter for inverse gamma prior
//' @return Updated sigma_ab value
// [[Rcpp::export]]
double sample_sigma_ab_cpp(const arma::mat& a_mat, const arma::mat& b_mat,
                           double rho_ab, bool symmetric,
                           double prior_shape = 2.0, double prior_scale = 1.0) {

  const int n = a_mat.n_rows;
  const int T = a_mat.n_cols;

  if(T <= 1) return 0.1;

  double ss_innov = 0.0;
  int count = 0;

  // stationary-init contribution (t = 1 states): each contributes
  // (sigma^2)^{-1/2} exp(-(1 - rho^2) state^2 / (2 sigma^2)) for coherence
  // with the state sampler's stationary initial condition
  const double one_m_rho2 = 1.0 - std::min(rho_ab * rho_ab, 0.9801);
  for(int i = 0; i < n; i++) {
    ss_innov += one_m_rho2 * a_mat(i,0) * a_mat(i,0);
    count++;
    if(!symmetric) {
      ss_innov += one_m_rho2 * b_mat(i,0) * b_mat(i,0);
      count++;
    }
  }

  for(int i = 0; i < n; i++) {
    for(int t = 1; t < T; t++) {
      double innov_a = a_mat(i,t) - rho_ab * a_mat(i,t-1);
      ss_innov += innov_a * innov_a;
      count++;

      if(!symmetric) {
        double innov_b = b_mat(i,t) - rho_ab * b_mat(i,t-1);
        ss_innov += innov_b * innov_b;
        count++;
      }
    }
  }

  double post_shape = prior_shape + 0.5 * count;
  double post_scale = prior_scale + 0.5 * ss_innov;

  double variance = post_scale / R::rgamma(post_shape, 1.0);

  return sqrt(variance);
}

//' Initialize dynamic additive effects with AR(1) structure
//'
//' @param n Number of actors
//' @param Tn Number of time points
//' @param rho_ab AR(1) parameter
//' @param sigma_ab Innovation standard deviation
//' @param mean_a Mean for row effects
//' @param mean_b Mean for column effects
//' @return List with initialized a and b matrices
// [[Rcpp::export]]
List init_dynamic_ab_cpp(int n, int Tn, double rho_ab, double sigma_ab,
                         double mean_a = 0.0, double mean_b = 0.0) {

  arma::mat a_mat(n, Tn);
  arma::mat b_mat(n, Tn);

  double var_stat = sigma_ab * sigma_ab / (1.0 - rho_ab * rho_ab);
  double sd_stat = sqrt(var_stat);

  for(int i = 0; i < n; i++) {
    a_mat(i, 0) = mean_a + R::rnorm(0, sd_stat);
    b_mat(i, 0) = mean_b + R::rnorm(0, sd_stat);

    for(int t = 1; t < Tn; t++) {
      a_mat(i, t) = rho_ab * a_mat(i, t-1) + R::rnorm(0, sigma_ab);
      b_mat(i, t) = rho_ab * b_mat(i, t-1) + R::rnorm(0, sigma_ab);
    }
  }

  return List::create(
    Named("a") = a_mat,
    Named("b") = b_mat
  );
}
