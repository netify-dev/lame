// Additional optimized functions for lame package
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Optimized precomputeX - much faster than R version for large arrays
// [[Rcpp::export]]
List precomputeX_cpp(const arma::cube& X) {
  const int n = X.n_rows;
  const int p = X.n_slices;
  const int n2 = n * n;
  
  // Pre-allocate outputs
  arma::mat Xr(n, p);     // row sums
  arma::mat Xc(n, p);     // col sums
  arma::mat mX(n2, p);    // matricized version
  arma::mat mXt(n2, p);   // dyad-transposed matricized version
  
  // Compute sums and matricized versions in one pass
  for(int k = 0; k < p; k++) {
    const arma::mat& Xk = X.slice(k);
    
    // Row and column sums
    Xr.col(k) = sum(Xk, 1);
    Xc.col(k) = sum(Xk, 0).t();
    
    // Matricized versions
    mX.col(k) = vectorise(Xk);
    mXt.col(k) = vectorise(Xk.t());
  }
  
  // Compute cross products
  arma::mat XX = mX.t() * mX;
  arma::mat XXt = mX.t() * mXt;
  
  return List::create(
    Named("Xr") = Xr,
    Named("Xc") = Xc,
    Named("mX") = mX,
    Named("mXt") = mXt,
    Named("XX") = XX,
    Named("XXt") = XXt
  );
}

// Optimized llsrmRho - vectorized computation over rho grid
// [[Rcpp::export]]
arma::vec llsrmRho_cpp(const arma::mat& Y, const arma::mat& Sab,
                        const arma::vec& rhos, double s2) {
  const int n = Y.n_rows;
  const int nr = rhos.n_elem;
  
  // Pre-compute Y statistics once
  const double sumY = accu(Y);
  const double sumY2 = accu(Y % Y);
  const double sumYtY = accu(Y % Y.t());
  const arma::vec rsy = sum(Y, 1);
  const arma::vec csy = sum(Y, 0).t();
  const double rr = dot(rsy, rsy);
  const double rc = dot(rsy, csy);
  const double cc = dot(csy, csy);
  
  arma::vec ll(nr);
  
  // Process each rho value
  for(int r = 0; r < nr; r++) {
    const double rho = rhos(r);
    const double sqrt_one_plus = sqrt(1.0 + rho);
    const double sqrt_one_minus = sqrt(1.0 - rho);
    const double sqrt_s2 = sqrt(s2);
    
    // Compute c and d
    const double c = 0.5 * (1.0/sqrt_one_plus + 1.0/sqrt_one_minus) / sqrt_s2;
    const double d = 0.5 * (1.0/sqrt_one_plus - 1.0/sqrt_one_minus) / sqrt_s2;
    
    // Build ihSe matrix
    arma::mat ihSe(2, 2);
    ihSe(0,0) = ihSe(1,1) = c;
    ihSe(0,1) = ihSe(1,0) = d;
    
    // Compute Sabt = ihSe * Sab * ihSe
    arma::mat Sabt = ihSe * Sab * ihSe;
    
    // Compute iSabt (inverse of 2x2 matrix)
    const double det = Sabt(0,0) * Sabt(1,1) - Sabt(0,1) * Sabt(1,0);
    arma::mat iSabt(2, 2);
    iSabt(0,0) = Sabt(1,1) / det;
    iSabt(1,1) = Sabt(0,0) / det;
    iSabt(0,1) = iSabt(1,0) = -Sabt(0,1) / det;
    
    // Compute G = inv(iSabt + n*I)
    arma::mat G = inv_sympd(iSabt + n * eye(2, 2));
    
    // Compute pH = -inv(iSabt + n*ones)
    arma::mat ones22(2, 2, arma::fill::ones);
    arma::mat pH = -inv_sympd(iSabt + n * ones22);
    
    // Compute H matrix elements
    const double sumH = pH(0,0) * G(1,0) + pH(0,1) * G(0,0) +
                       pH(0,0) * G(1,1) + pH(0,1) * G(0,1) +
                       pH(1,0) * G(1,0) + pH(1,1) * G(0,0) +
                       pH(1,0) * G(1,1) + pH(1,1) * G(0,1);
    
    // Compute log determinant of variance
    double ldV = n * n * log(s2) +
                 0.5 * n * (n-1) * log(1.0 + rho) +
                 0.5 * n * (n-1) * log(1.0 - rho);
    
    // Add determinant terms
    ldV += (n-1) * log(arma::det(eye(2,2) + n * Sabt));
    ldV += log(arma::det(eye(2,2) + n * Sabt * ones22));
    
    // Compute quadratic form
    const double sumZ = (c + d) * sumY;
    const double sumZ2 = (c*c + d*d) * sumY2 + 2*c*d * sumYtY;
    
    const double B = G(0,0) * (c*c * rr + d*d * cc + 2*c*d * rc) +
                     2 * G(0,1) * ((c*c + d*d) * rc + c*d * (rr + cc)) +
                     G(1,1) * (c*c * cc + d*d * rr + 2*c*d * rc);
    
    const double yPy = sumZ2 - B - sumZ * sumZ * sumH;
    
    ll(r) = -0.5 * (ldV + yPy);
  }
  
  return ll;
}

// Optimized rbeta_ab_fc - core computation in C++
// [[Rcpp::export]]
List rbeta_ab_fc_cpp(const arma::mat& Z, const arma::mat& Sab, double rho,
                     const arma::cube& X, double s2, const arma::mat& offset,
                     const arma::mat& iV0, const arma::vec& m0, double g) {
  
  const int n = Z.n_rows;
  const int p = X.n_slices;
  
  // Compute Z - offset
  arma::mat Z_off = Z - offset;
  
  // Precompute X statistics if not already done
  arma::mat Xr(n, p), Xc(n, p), mX(n*n, p), mXt(n*n, p), XX, XXt;
  
  for(int k = 0; k < p; k++) {
    const arma::mat& Xk = X.slice(k);
    Xr.col(k) = sum(Xk, 1);
    Xc.col(k) = sum(Xk, 0).t();
    mX.col(k) = vectorise(Xk);
    mXt.col(k) = vectorise(Xk.t());
  }
  XX = mX.t() * mX;
  XXt = mX.t() * mXt;
  
  // Decorrelation
  arma::mat Se(2, 2);
  Se(0,0) = Se(1,1) = 1.0;
  Se(0,1) = Se(1,0) = rho;
  Se *= s2;
  
  // Compute matrix square root of inverse
  arma::mat iSe_inv = inv_sympd(Se);
  arma::vec eigval_se;
  arma::mat eigvec_se;
  eig_sym(eigval_se, eigvec_se, iSe_inv);
  arma::mat iSe2 = eigvec_se * diagmat(sqrt(eigval_se)) * eigvec_se.t();
  const double td = iSe2(0,0);
  const double to = iSe2(0,1);
  
  arma::mat Sabs = iSe2 * Sab * iSe2;
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, Sabs);
  
  // Find positive eigenvalues
  int k = 0;
  for(int i = 0; i < eigval.n_elem; i++) {
    if(eigval(i) > 1e-10) k++;
  }
  
  // Transformed quantities
  arma::mat mXs = td * mX + to * mXt;
  arma::mat XXs = (to*to + td*td) * XX + 2*to*td * XXt;
  arma::mat Zs = td * Z_off + to * Z_off.t();
  arma::vec zr = sum(Zs, 1);
  arma::vec zc = sum(Zs, 0).t();
  double zs = accu(zc);
  
  // Initialize outputs
  arma::vec beta(p);
  arma::vec a(n, fill::zeros);
  arma::vec b(n, fill::zeros);
  
  if(p > 0) {
    // Dyadic and prior contributions
    arma::vec lb = mXs.t() * vectorise(Zs) + iV0 * m0;
    arma::mat Qb = XXs + iV0;
    
    // Row and column reduction if k > 0
    if(k > 0) {
      arma::mat G = eigvec.cols(eigvec.n_cols - k, eigvec.n_cols - 1);
      G.each_col() %= sqrt(eigval.tail(k));
      
      arma::mat A = n * G.t() * G + eye(k, k);
      arma::mat K(2, 2);
      K(0,0) = K(1,1) = 0;
      K(0,1) = K(1,0) = 1;
      arma::mat B = G.t() * K * G;
      
      arma::mat iA0 = inv_sympd(A);
      arma::mat C0 = -inv_sympd(A + n*B) * B * iA0;
      
      arma::mat iA = G * iA0 * G.t();
      arma::mat C = G * C0 * G.t();
      
      // Update lb and Qb
      arma::mat Xsr = td * Xr + to * Xc;
      arma::mat Xsc = td * Xc + to * Xr;
      arma::vec Xss = sum(Xsc, 0).t();
      
      lb -= iA(0,0) * Xsr.t() * zr + iA(1,1) * Xsc.t() * zc +
            iA(0,1) * (Xsr.t() * zc + Xsc.t() * zr) +
            accu(C) * Xss * zs;
      
      arma::mat tmp = Xsr.t() * Xsc;
      Qb -= iA(0,0) * Xsr.t() * Xsr + iA(1,1) * Xsc.t() * Xsc +
            iA(1,0) * (tmp + tmp.t()) + accu(C) * Xss * Xss.t();
    }
    
    // Sample beta
    arma::mat V = inv_sympd(Qb);
    arma::vec m = V * lb;
    beta = mvnrnd(m, V);
    
    // Sample a, b if k > 0
    if(k > 0) {
      // Compute residuals
      arma::mat E = Zs;
      for(int i = 0; i < p; i++) {
        E -= beta(i) * (td * X.slice(i) + to * X.slice(i).t());
      }
      
      arma::vec er = sum(E, 1);
      arma::vec ec = sum(E, 0).t();
      double es = accu(ec);
      
      arma::mat G = eigvec.cols(eigvec.n_cols - k, eigvec.n_cols - 1);
      G.each_col() %= sqrt(eigval.tail(k));
      
      arma::mat A = n * G.t() * G + eye(k, k);
      arma::mat iA0 = inv_sympd(A);
      // Ensure iA0 is symmetric before Cholesky
      iA0 = 0.5 * (iA0 + iA0.t());
      
      arma::mat m_ab = iA0 * G.t() * join_cols(er, ec);
      arma::mat hiA0 = chol(iA0, "lower");
      
      arma::mat e(n, k, fill::randn);
      arma::mat w = e * hiA0.t();
      
      arma::mat ab = (w + repmat(m_ab.t(), n, 1)) * G.t() * inv(iSe2);
      a = ab.col(0);
      b = ab.col(1);
    }
  }
  
  return List::create(
    Named("beta") = beta,
    Named("a") = a,
    Named("b") = b
  );
}

// Optimized ldZgbme - for GBME models
// [[Rcpp::export]]
arma::mat ldZgbme_opt_cpp(const arma::mat& Z, const arma::mat& Y,
                          const arma::mat& EZ, double rho, double s2) {
  // Constants
  const double sqrt_one_plus = sqrt(1.0 + rho);
  const double sqrt_one_minus = sqrt(1.0 - rho);
  const double sqrt_s2 = sqrt(s2);
  
  const double c = 0.5 * (1.0/sqrt_one_plus + 1.0/sqrt_one_minus) / sqrt_s2;
  const double d = 0.5 * (1.0/sqrt_one_plus - 1.0/sqrt_one_minus) / sqrt_s2;
  
  // Compute E = Z - EZ
  arma::mat E = Z - EZ;
  
  // Compute lpZ efficiently
  arma::mat lpZ = -0.5 * square(c * E + d * E.t());
  lpZ += lpZ.t();
  lpZ.diag() *= 0.5;
  
  // For normal likelihood (most common case)
  // Can be extended for other distributions
  arma::mat llZ = -0.5 * square(Y - Z) / s2;
  llZ += llZ.t();
  
  // Handle missing values
  arma::uvec na_idx = find_nonfinite(Y);
  llZ.elem(na_idx).zeros();
  
  return lpZ + llZ;
}

// Optimized array to list conversion
// [[Rcpp::export]]
List array_to_list_cpp(const arma::cube& arr, 
                       const List& actorByYr,
                       const CharacterVector& pdLabs,
                       const IntegerMatrix& actorIndices) {
  const int T = arr.n_slices;
  List result(T);
  
  for(int t = 0; t < T; t++) {
    CharacterVector actors_t = actorByYr[t];
    const int n_t = actors_t.size();
    
    // Extract submatrix for actors at time t
    arma::mat mat_t(n_t, n_t);
    
    // Use precomputed indices
    IntegerVector idx = actorIndices(_, t);
    
    // Extract submatrix
    for(int i = 0; i < n_t; i++) {
      for(int j = 0; j < n_t; j++) {
        mat_t(i,j) = arr(idx[i] - 1, idx[j] - 1, t);  // R uses 1-based indexing
      }
    }
    
    // Set dimnames
    NumericMatrix mat_wrap = wrap(mat_t);
    rownames(mat_wrap) = actors_t;
    colnames(mat_wrap) = actors_t;
    
    result[t] = mat_wrap;
  }
  
  result.names() = pdLabs;
  return result;
}

// Optimized rrho_fc using efficient grid search
// [[Rcpp::export]]
double rrho_fc_cpp(const arma::mat& Z, const arma::mat& Sab, double s2,
                   const arma::mat& offset, int ngp, bool asp) {
  
  arma::mat E = Z - offset;
  
  // Rough estimate of rho and its sd
  double m = mean(vectorise(E));
  arma::vec a = mean(E, 1) - m;
  arma::vec b = mean(E, 0).t() - m;
  
  arma::mat R = E - (m + a * ones(1, E.n_cols) + ones(E.n_rows, 1) * b.t());
  
  // Extract upper triangular elements
  arma::uvec upper_idx = find(trimatu(ones<mat>(E.n_rows, E.n_cols), 1));
  arma::vec R_upper = R.elem(upper_idx);
  arma::mat Rt = R.t();
  arma::vec Rt_upper = Rt.elem(upper_idx);
  
  arma::mat RM = join_rows(R_upper, Rt_upper) / sqrt(s2);
  
  double emcp = dot(RM.col(0), RM.col(1));
  double emss = accu(square(RM));
  int nd = RM.n_rows;
  double rho0 = dot(RM.col(0), RM.col(1)) / sqrt(dot(RM.col(0), RM.col(0)) * dot(RM.col(1), RM.col(1)));
  double srho0 = 2.0 * (1.0 - rho0*rho0) / sqrt(nd);
  
  // Grid concentrated near likely regions
  arma::vec rhos(ngp);
  for(int i = 0; i < ngp; i++) {
    double q = (i + 1.0) / (ngp + 1.0);
    rhos(i) = R::qnorm(q, rho0, 2*srho0, true, false);
  }
  
  // Keep only valid rhos
  rhos = rhos.elem(find(rhos > -1 && rhos < 1));
  
  // Compute log likelihood on grid
  arma::vec ll = llsrmRho_cpp(E, Sab, rhos, s2);
  
  // Apply prior if requested
  arma::vec prho = exp(ll - max(ll));
  if(asp) {
    prho %= pow(1.0 - square(rhos), -0.5);
  }
  
  // Sample from discrete distribution
  arma::vec Frho = cumsum(prho) / sum(prho);
  double f = R::runif(0, 1);
  
  int k = 0;
  for(int i = 0; i < Frho.n_elem; i++) {
    if(Frho(i) > f) {
      k = i;
      break;
    }
  }
  
  // Linear interpolation
  double rho;
  if(k == 0) {
    rho = rhos(0);
  } else if(k == rhos.n_elem - 1) {
    rho = rhos(k);
  } else {
    double w = (f - Frho(k-1)) / (Frho(k) - Frho(k-1));
    rho = rhos(k-1) + w * (rhos(k) - rhos(k-1));
  }
  
  return std::min(std::abs(rho), 0.995) * ((rho >= 0) ? 1 : -1);
}