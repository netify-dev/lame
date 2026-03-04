// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

//' Sample dynamic additive effects with AR(1) evolution
//'
//' @param a_current Current 2D array of row effects (n x T)
//' @param b_current Current 2D array of column effects (n x T)
//' @param Z_array 3D array of latent positions (n x n x T)
//' @param EZ_array 3D array of expected values without additive effects (n x n x T)
//' @param rho_ab AR(1) parameter for additive effects
//' @param sigma_ab Innovation standard deviation
//' @param Sab Covariance matrix for a and b (2x2)
//' @param symmetric Whether the network is symmetric
//' @return List with updated a and b arrays
// [[Rcpp::export]]
List sample_dynamic_ab_cpp(arma::mat a_current, arma::mat b_current,
                           const arma::cube& Z_array, const arma::cube& EZ_array,
                           double rho_ab, double sigma_ab,
                           const arma::mat& Sab, bool symmetric) {

  const int n = a_current.n_rows;
  const int T = a_current.n_cols;

  arma::mat a_new = a_current;
  arma::mat b_new = b_current;

  double var_innov = sigma_ab * sigma_ab;
  double var_stationary = var_innov / (1.0 - rho_ab * rho_ab);

  arma::mat Sab_safe = Sab;
  double min_eigenval = arma::min(arma::eig_sym(Sab_safe));
  if (min_eigenval <= 0) {
    Sab_safe += arma::eye(2, 2) * (1e-6 - min_eigenval);
  }

  arma::mat L_Sab;
  bool chol_success = arma::chol(L_Sab, Sab_safe, "lower");
  if (!chol_success) {
    L_Sab = arma::diagmat(arma::sqrt(arma::diagvec(Sab_safe)));
  }

  for(int i = 0; i < n; i++) {

    arma::vec a_i = a_current.row(i).t();
    arma::vec b_i = b_current.row(i).t();

    for(int t = 0; t < T; t++) {

      arma::mat resid_t = Z_array.slice(t) - EZ_array.slice(t);

      double prior_mean_a = 0.0;
      double prior_mean_b = 0.0;
      double prior_prec = 1.0 / var_stationary;

      if(t > 0) {
        prior_mean_a = rho_ab * a_i(t-1);
        prior_mean_b = rho_ab * b_i(t-1);
        prior_prec = 1.0 / var_innov;
      }

      if(t < T-1) {
        prior_mean_a = 0.5 * (prior_mean_a + rho_ab * a_i(t+1));
        prior_mean_b = 0.5 * (prior_mean_b + rho_ab * b_i(t+1));
        prior_prec *= 2.0;
      }

      arma::vec row_sum = sum(resid_t, 1);
      arma::vec col_sum = sum(resid_t, 0).t();

      row_sum(i) -= n * a_i(t);
      col_sum(i) -= n * b_i(t);

      double data_prec_a = n / Sab(0,0);
      double post_prec_a = prior_prec + data_prec_a;
      double post_mean_a = (prior_prec * prior_mean_a + data_prec_a * row_sum(i)/n) / post_prec_a;
      double post_sd_a = 1.0 / sqrt(post_prec_a);

      a_i(t) = R::rnorm(post_mean_a, post_sd_a);

      if(!symmetric) {
        double data_prec_b = n / Sab(1,1);
        double post_prec_b = prior_prec + data_prec_b;
        double post_mean_b = (prior_prec * prior_mean_b + data_prec_b * col_sum(i)/n) / post_prec_b;
        double post_sd_b = 1.0 / sqrt(post_prec_b);

        b_i(t) = R::rnorm(post_mean_b, post_sd_b);
      } else {
        b_i(t) = a_i(t);
      }
    }

    a_new.row(i) = a_i.t();
    b_new.row(i) = b_i.t();
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
//' @return Updated rho value
// [[Rcpp::export]]
double sample_rho_ab_cpp(const arma::mat& a_mat, const arma::mat& b_mat,
                         double sigma_ab, double rho_current, bool symmetric) {

  const int n = a_mat.n_rows;
  const int T = a_mat.n_cols;

  if(T <= 1) return rho_current;

  double sum_prod = 0.0;
  double sum_sq_lag = 0.0;
  double sum_sq_curr = 0.0;
  int count = 0;

  for(int i = 0; i < n; i++) {
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

  double rho_mle = sum_prod / sum_sq_lag;
  rho_mle = std::max(-0.99, std::min(0.99, rho_mle));

  double proposal_sd = 0.1;
  double rho_proposal = R::rnorm(rho_current, proposal_sd);

  if(std::abs(rho_proposal) >= 1.0) {
    return rho_current;
  }

  double ll_current = -0.5 * count * log(2.0 * M_PI * sigma_ab * sigma_ab);
  ll_current -= 0.5 * (sum_sq_curr - 2.0 * rho_current * sum_prod +
                       rho_current * rho_current * sum_sq_lag) / (sigma_ab * sigma_ab);

  double ll_proposal = -0.5 * count * log(2.0 * M_PI * sigma_ab * sigma_ab);
  ll_proposal -= 0.5 * (sum_sq_curr - 2.0 * rho_proposal * sum_prod +
                        rho_proposal * rho_proposal * sum_sq_lag) / (sigma_ab * sigma_ab);

  double log_prior_current = -0.5 * log(1.0 - rho_current * rho_current);
  double log_prior_proposal = -0.5 * log(1.0 - rho_proposal * rho_proposal);

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
//' @param T Number of time points
//' @param rho_ab AR(1) parameter
//' @param sigma_ab Innovation standard deviation
//' @param mean_a Mean for row effects
//' @param mean_b Mean for column effects
//' @return List with initialized a and b matrices
// [[Rcpp::export]]
List init_dynamic_ab_cpp(int n, int T, double rho_ab, double sigma_ab,
                         double mean_a = 0.0, double mean_b = 0.0) {

  arma::mat a_mat(n, T);
  arma::mat b_mat(n, T);

  double var_stat = sigma_ab * sigma_ab / (1.0 - rho_ab * rho_ab);
  double sd_stat = sqrt(var_stat);

  for(int i = 0; i < n; i++) {
    a_mat(i, 0) = mean_a + R::rnorm(0, sd_stat);
    b_mat(i, 0) = mean_b + R::rnorm(0, sd_stat);

    for(int t = 1; t < T; t++) {
      a_mat(i, t) = rho_ab * a_mat(i, t-1) + R::rnorm(0, sigma_ab);
      b_mat(i, t) = rho_ab * b_mat(i, t-1) + R::rnorm(0, sigma_ab);
    }
  }

  return List::create(
    Named("a") = a_mat,
    Named("b") = b_mat
  );
}
