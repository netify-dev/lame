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
   
   arma::cube EZ = arma::zeros(n,n,N);
   
   for(int t=0 ; t<N ; ++t){
     
     arma::cube Xt = Xlist[t];
     arma::mat Xbeta = arma::zeros(n,n);
     
     for(int i=0 ; i<p ; ++i){
       Xbeta = Xbeta + beta(i) * Xt.slice(i);
     }
     
     EZ.slice(t) = Xbeta + ab + U*V.t();
   }	
   
   return( EZ );
   
 }