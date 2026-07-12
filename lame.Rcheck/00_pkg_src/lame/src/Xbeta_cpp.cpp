#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat Xbeta_cpp(const arma::cube& X, const arma::vec& beta) {
  int n = X.n_rows;
  int m = X.n_cols;
  int p = beta.n_elem;
  if (p != static_cast<int>(X.n_slices)) {
    Rcpp::stop("length(beta) must match dim(X)[3].");
  }
  
  arma::mat XB(n, m, fill::zeros);
  
  for(int k = 0; k < p; k++) {
    XB += beta(k) * X.slice(k);
  }
  
  return XB;
}

// [[Rcpp::export]]
arma::mat outer_cpp(const arma::vec& a, const arma::vec& b) {
  return a * b.t();
}

// [[Rcpp::export]]
arma::mat gof_stats_cpp(const arma::mat& Y) {
  arma::vec stats(5, fill::zeros);

  // Self-ties are outside the unipartite likelihood. Treat the diagonal as
  // missing even if a caller supplied finite placeholder values.
  arma::mat Y_clean = Y;
  Y_clean.diag().fill(arma::datum::nan);

  int n = Y_clean.n_rows;
  arma::vec row_means(n);
  int n_valid_rows = 0;
  for(int i = 0; i < n; i++) {
    arma::rowvec row_i = Y_clean.row(i);
    arma::uvec finite_idx = arma::find_finite(row_i);
    if(finite_idx.n_elem > 0) {
      row_means(n_valid_rows) = arma::mean(row_i(finite_idx));
      n_valid_rows++;
    }
  }
  if(n_valid_rows > 1) {
    stats(0) = arma::stddev(row_means.head(n_valid_rows));
  } else {
    stats(0) = 0.0;
  }

  int m = Y_clean.n_cols;
  arma::vec col_means(m);
  int n_valid_cols = 0;
  for(int j = 0; j < m; j++) {
    arma::vec col_j = Y_clean.col(j);
    arma::uvec finite_idx = arma::find_finite(col_j);
    if(finite_idx.n_elem > 0) {
      col_means(n_valid_cols) = arma::mean(col_j(finite_idx));
      n_valid_cols++;
    }
  }
  if(n_valid_cols > 1) {
    stats(1) = arma::stddev(col_means.head(n_valid_cols));
  } else {
    stats(1) = 0.0;
  }
  
  // Dyadic dependence (correlation between Y and Y')
  arma::vec y_vec = arma::vectorise(Y_clean);
  arma::vec yt_vec = arma::vectorise(Y_clean.t());
  
  arma::uvec valid = arma::intersect(arma::find_finite(y_vec), arma::find_finite(yt_vec));
  if(valid.n_elem > 1) {
    arma::vec y_valid = y_vec(valid);
    arma::vec yt_valid = yt_vec(valid);
    double y_sd_dyad = arma::stddev(y_valid);
    double yt_sd_dyad = arma::stddev(yt_valid);
    if(std::isfinite(y_sd_dyad) && std::isfinite(yt_sd_dyad) &&
       y_sd_dyad > 0.0 && yt_sd_dyad > 0.0) {
      double rho = arma::as_scalar(arma::cor(y_valid, yt_valid));
      stats(2) = std::isfinite(rho) ? rho : 0.0;
    }
  }

  arma::uvec idx_f = arma::find_finite(y_vec);
  double y_mean = idx_f.n_elem ? arma::mean(y_vec(idx_f)) : 0.0;
  double y_sd = idx_f.n_elem > 1 ? arma::stddev(y_vec(idx_f)) : 1.0;
  
  arma::mat E = Y_clean - y_mean;
  E.elem(arma::find_nonfinite(E)).zeros();
  
  arma::mat D(Y_clean.n_rows, Y_clean.n_cols, fill::ones);
  D.elem(arma::find_nonfinite(Y_clean)).zeros();
  
  // Cycle dependence
  arma::mat E3 = E * E * E;
  arma::mat D3 = D * D * D;
  double cycle_denom = arma::trace(D3) * std::pow(y_sd, 3);
  if(std::isfinite(cycle_denom) && cycle_denom > 0.0) {
    double cycle = arma::trace(E3) / cycle_denom;
    stats(3) = std::isfinite(cycle) ? cycle : 0.0;
  }

  // Transitive dependence
  arma::mat EtE = E * E.t() * E;
  arma::mat DtD = D * D.t() * D;
  double trans_denom = arma::trace(DtD) * std::pow(y_sd, 3);
  if(std::isfinite(trans_denom) && trans_denom > 0.0) {
    double trans = arma::trace(EtE) / trans_denom;
    stats(4) = std::isfinite(trans) ? trans : 0.0;
  }

  return stats.t();
}
