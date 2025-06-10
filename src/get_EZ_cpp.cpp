//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Linear combinations of submatrices of an array
 //' 
 //' Computes a matrix of expected values based on an array X of predictors,
 //' vector beta of regression coefficients, the outer product of the a and b
 //' random effects, and the third order effects U and V.
 //' 
 //' @usage Xbeta_cpp(X, beta)
 //' @param X a N length list filled with n by n by p arrays
 //' @param beta a p by 1 vector
 //' @param ab a n by n matrix calculated by taking the outer 
 //' product of the a and b random effect vectors
 //' @param U an n by k matrix of multiplicative row effects
 //' @param V an n by k matrix of multiplicative column effects
 //' @return An n by n by N array
 //' @author Peter Hoff, Shahryar Minhas
 //' @export get_EZ_cpp
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