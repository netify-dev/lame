//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// Linear combinations of submatrices of an array
// Computes a matrix of expected values based on an array X of predictors,
// vector beta of regression coefficients, the outer product of the a and b
// random effects, and the third order effects U and V.

 // [[Rcpp::export]]

 arma::cube get_EZ_cpp(
     Rcpp::List Xlist, arma::vec beta, arma::mat ab,
     arma::mat U, arma::mat V
 ) {

   int N = Xlist.size();
   int n = ab.n_rows;
   int p = beta.size();

   arma::cube EZ(n, n, N);

   // Pre-compute UV product ONCE (Phase 1B optimization)
   arma::mat UVt = U * V.t();

   // Pre-compute ab + UVt (constant across time)
   arma::mat ab_uv = ab + UVt;

   // Allocate Xbeta ONCE outside loop (Phase 1B optimization)
   arma::mat Xbeta(n, n);

   for(int t = 0; t < N; ++t) {

     if(p > 0 && Xlist.size() > 0) {
       Xbeta.zeros();  // Reset without reallocation
       arma::cube Xt = Xlist[t];
       int num_covs = std::min(p, (int)Xt.n_slices);

       for(int i = 0; i < num_covs; ++i) {
         if(beta(i) != 0.0) {
           Xbeta += beta(i) * Xt.slice(i);  // In-place accumulation
         }
       }
       EZ.slice(t) = Xbeta + ab_uv;
     } else {
       EZ.slice(t) = ab_uv;
     }
   }

   return(EZ);

 }
