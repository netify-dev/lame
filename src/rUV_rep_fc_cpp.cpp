//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// Helper: build index vector excluding one element (Phase 1D optimization)
static inline arma::uvec indices_excluding(int total, int exclude) {
  arma::uvec idx(total - 1);
  int k = 0;
  for(int i = 0; i < total; i++) {
    if(i != exclude) {
      idx(k++) = i;
    }
  }
  return idx;
}

// [[Rcpp::export]]
arma::mat rwish_cpp(const arma::mat& S0, int nu) {
  int n = S0.n_rows;

  arma::mat S0_clean = S0;
  if (!S0.is_finite()) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(!std::isfinite(S0_clean(i,j))) {
          S0_clean(i,j) = (i == j) ? 1.0 : 0.0;
        }
      }
    }
  }

  arma::mat S0_sym = 0.5 * (S0_clean + S0_clean.t());

  arma::vec eigval;
  arma::mat eigvec;
  bool eig_success = arma::eig_sym(eigval, eigvec, S0_sym);

  if (!eig_success) {
    S0_sym = arma::eye(n, n);
    arma::eig_sym(eigval, eigvec, S0_sym);
  }

  double min_eig = eigval.min();
  double max_eig = eigval.max();

  if (max_eig > 0 && min_eig/max_eig < 1e-10) {
    double regularization = max_eig * 1e-8;
    for(int i = 0; i < (int)eigval.n_elem; i++) {
      if(eigval(i) < regularization) {
        eigval(i) = regularization;
      }
    }
    S0_sym = eigvec * arma::diagmat(eigval) * eigvec.t();
  } else if (min_eig <= 0) {
    double ridge = std::max(1e-6, 1e-6 - min_eig);
    eigval += ridge;
    S0_sym = eigvec * arma::diagmat(eigval) * eigvec.t();
  }

  S0_sym = 0.5 * (S0_sym + S0_sym.t());

  arma::mat sS0;
  bool chol_success = arma::chol(sS0, S0_sym, "lower");

  if (!chol_success) {
    arma::eig_sym(eigval, eigvec, S0_sym);
    eigval = arma::max(eigval, arma::ones(n) * 1e-10);
    arma::mat sqrt_eigval = arma::diagmat(arma::sqrt(eigval));
    sS0 = eigvec * sqrt_eigval;
  }

  arma::mat Z = arma::randn(nu, n) * sS0.t();

  return Z.t() * Z;
}

// Fast UV update with efficient memory usage
// [[Rcpp::export]]
List rUV_rep_fc_cpp(
    const arma::cube& ET, arma::mat U, arma::mat V,
    double rho, double s2, const arma::mat& iSe2,
    double maxmargin, bool shrink, const NumericVector& rLoopIDs) {

  int Time = ET.n_slices;
  int R = U.n_cols;
  int n = U.n_rows;

  // Guard: this function assumes square (unipartite) networks
  if(ET.n_rows != ET.n_cols) {
    Rcpp::stop("rUV_rep_fc_cpp requires square residual matrices (unipartite). "
               "Got %d x %d. Use the R fallback for bipartite networks.",
               ET.n_rows, ET.n_cols);
  }

  arma::mat UV = join_rows(U, V);
  arma::mat Suv;

  if(shrink) {
    arma::mat UV_prod = UV.t() * UV;
    arma::mat inv_mat = inv_sympd(eye<mat>(R*2, R*2) + UV_prod);
    Suv = inv_sympd(rwish_cpp(inv_mat, n + R + 2));
  } else {
    Suv = eye<mat>(R*2, R*2);
  }

  double g = iSe2(0,0);
  double d = iSe2(0,1);
  double g2d2 = g*g + d*d;

  arma::mat Est(n * n, Time);
  arma::mat Usmall(n, R-1);
  arma::mat Vsmall(n, R-1);

  // Phase 1D: Pre-compute index vectors for submat extraction
  int dim2R = 2 * R;

  for(int i = 0; i < rLoopIDs.size(); i++) {
    int r = rLoopIDs[i];

    if(r == 0) {
      Usmall = U.cols(1, R-1);
      Vsmall = V.cols(1, R-1);
    } else if(r == R-1) {
      Usmall = U.cols(0, R-2);
      Vsmall = V.cols(0, R-2);
    } else {
      Usmall = join_rows(U.cols(0, r-1), U.cols(r+1, R-1));
      Vsmall = join_rows(V.cols(0, r-1), V.cols(r+1, R-1));
    }

    arma::mat UVmr = Usmall * Vsmall.t();

    for(int t = 0; t < Time; t++) {
      arma::mat ert = ET.slice(t) - UVmr;
      Est.col(t) = vectorise(g2d2 * ert + 2 * g * d * ert.t());
    }

    arma::vec vr = V.col(r);

    // Phase 1D: Use submat with index vectors instead of shed_col/shed_row
    arma::uvec idx_excl_r = indices_excluding(dim2R, r);
    arma::uvec idx_r_only(1); idx_r_only(0) = r;

    arma::rowvec Suvsmall = Suv(idx_r_only, idx_excl_r);
    arma::mat Suvsmall2 = Suv(idx_excl_r, idx_excl_r);

    arma::vec b0 = (Suvsmall * inv_sympd(Suvsmall2)).t();
    arma::vec Suv_col_r = Suv(idx_excl_r, idx_r_only);
    double v0 = Suv(r,r) - dot(b0, Suv_col_r);

    arma::vec m0 = join_rows(Usmall, V) * b0;

    double sumvr2 = dot(vr, vr);
    double ssv = (sumvr2 >= maxmargin) ? sumvr2 : maxmargin;

    double a = Time * g2d2 * ssv + 1/v0;
    double c = -2 * Time * g * d / (a*a + a * 2 * Time * g * d * ssv);

    arma::vec Esv_vec = sum(Est, 1);
    arma::mat Esv_mat = arma::reshape(Esv_vec, n, n);
    arma::vec Esv = Esv_mat * vr;
    arma::vec m1 = Esv/a + c * vr * dot(Esv + m0/v0, vr) + m0/(a*v0);

    double ah = sqrt(1/a);
    double bh = (sqrt(1/a + ssv*c) - sqrt(1/a)) / ssv;

    arma::vec e = randn(n);
    U.col(r) = m1 + ah * e + bh * vr * dot(vr, e);

    arma::vec ur = U.col(r);
    int rv = R + r;

    // Phase 1D: Use submat for V update too
    arma::uvec idx_excl_rv = indices_excluding(dim2R, rv);
    arma::uvec idx_rv_only(1); idx_rv_only(0) = rv;

    arma::rowvec Suvsmall_v = Suv(idx_rv_only, idx_excl_rv);
    arma::mat Suvsmall2_v = Suv(idx_excl_rv, idx_excl_rv);

    arma::vec b0_v = (Suvsmall_v * inv_sympd(Suvsmall2_v)).t();
    arma::vec Suv_col_rv = Suv(idx_excl_rv, idx_rv_only);
    double v0_v = Suv(rv,rv) - dot(b0_v, Suv_col_rv);

    arma::vec m0_v = join_rows(U, Vsmall) * b0_v;

    double sumur2 = dot(ur, ur);
    double ssu = (sumur2 >= maxmargin) ? sumur2 : maxmargin;

    double a_v = Time * g2d2 * ssu + 1/v0_v;
    double c_v = -2 * Time * g * d / (a_v*a_v + a_v * 2 * Time * g * d * ssu);

    arma::vec tEsu = Esv_mat.t() * ur;
    arma::vec m1_v = tEsu/a_v + c_v * ur * dot(tEsu + m0_v/v0_v, ur) + m0_v/(a_v*v0_v);

    double ah_v = sqrt(1/a_v);
    double bh_v = (sqrt(1/a_v + ssu*c_v) - sqrt(1/a_v)) / ssu;

    arma::vec e_v = randn(n);
    V.col(r) = m1_v + ah_v * e_v + bh_v * ur * dot(ur, e_v);
  }

  return List::create(
    Named("U") = U,
    Named("V") = V
  );
}
