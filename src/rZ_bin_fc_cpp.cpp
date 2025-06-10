//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Simulate Z based on a probit model
 //' 
 //' Simulates a random latent matrix Z given its expectation, dyadic correlation
 //' and a binary relational matrix Y
 //' 
 //' 
 //' @usage rZ_bin_fc_cpp(Z, EZ, rho, Y)
 //' @param Z a square matrix, the current value of Z
 //' @param EZ expected value of Z
 //' @param rho dyadic correlation
 //' @param Y square binary relational matrix
 //' @return a square matrix , the new value of Z
 //' @author Peter Hoff, Shahryar Minhas
 //' @export rZ_bin_fc_cpp
 // [[Rcpp::export]]
 
 arma::cube rZ_bin_fc_cpp(
     arma::cube ZT, arma::cube EZT, double rho, arma::cube YT
 ) {
   
   for(int t=0 ; t<ZT.n_slices; ++t){
     
     arma::mat Z = ZT.slice(t);
     arma::mat Y = YT.slice(t);
     arma::mat EZ = EZT.slice(t);
     
     double sz = pow( ( 1 - pow(rho,2) ), .5);
     
     arma::uvec ut = find(trimatl(arma::ones(EZ.n_rows, EZ.n_cols))==0);
     IntegerVector utIndex = as<IntegerVector>(wrap( ut ));	
     arma::uvec lt = find(trimatu(arma::ones(EZ.n_rows, EZ.n_cols))==0);
     IntegerVector ltIndex = as<IntegerVector>(wrap( lt ));	
     
     Y.replace(datum::nan, -1);	
     
     arma::mat Zt = Z.t(); arma::mat EZt = EZ.t();
     
     IntegerVector yVals = IntegerVector::create(-1,0,1);
     NumericVector lbVals = NumericVector::create(-1*datum::inf,-1*datum::inf,0);
     NumericVector ubVals = NumericVector::create(datum::inf,0,datum::inf);
     for(int i=0; i<yVals.size(); ++i){
       int y=yVals[i]; double lb=lbVals[i]; double ub=ubVals[i];
       IntegerVector Yyindex = as<IntegerVector>(wrap( find(Y==y) ));
       
       // upper triangle
       arma::uvec up_u = as<arma::uvec>(wrap( intersect(Yyindex, utIndex) ));
       arma::vec ez_u = EZ.elem(up_u) + rho * ( Zt.elem(up_u) - EZt.elem(up_u) );
       
       NumericVector x1_u; x1_u = (lb-ez_u)/sz; NumericVector x2_u; x2_u = (ub-ez_u)/sz;
       NumericVector rUnifLo_u = pnorm(x1_u,0.0,1.0); NumericVector rUnifHi_u = pnorm(x2_u,0.0,1.0);
       NumericVector rUnifDraw_u(up_u.n_elem);
       for(int ii=0;ii<up_u.n_elem;++ii){ rUnifDraw_u[ii]=runif(1,rUnifLo_u[ii],rUnifHi_u[ii])[0]; }
       NumericVector qnormDraw_u = qnorm(rUnifDraw_u,0.0,1.0); arma::vec szQnorm_u = sz * qnormDraw_u;
       Z.elem(up_u) = ez_u + szQnorm_u;
       
       // lower triangle
       arma::uvec up_l = as<arma::uvec>(wrap( intersect(Yyindex, ltIndex) ));
       arma::vec ez_l = EZ.elem(up_l) + rho * ( Zt.elem(up_l) - EZt.elem(up_l) );
       
       NumericVector x1_l; x1_l = (lb-ez_l)/sz; NumericVector x2_l; x2_l = (ub-ez_l)/sz;
       NumericVector rUnifLo_l = pnorm(x1_l,0.0,1.0); NumericVector rUnifHi_l = pnorm(x2_l,0.0,1.0);
       NumericVector rUnifDraw_l(up_l.n_elem);
       for(int ii=0;ii<up_l.n_elem;++ii){ rUnifDraw_l[ii]=runif(1,rUnifLo_l[ii],rUnifHi_l[ii])[0]; }
       NumericVector qnormDraw_l = qnorm(rUnifDraw_l,0.0,1.0); arma::vec szQnorm_l = sz * qnormDraw_l;
       Z.elem(up_l) = ez_l + szQnorm_l;
     }
     
     arma::vec ezDiag = diagvec(EZ);
     arma::vec zDiag(Z.n_rows); for(int i=0; i<Z.n_rows; ++i){ zDiag[i]=rnorm(1,ezDiag[i],1)[0]; }
     Z.diag() = zDiag;
     
     ZT.slice(t) = Z;
   }
   
   return( ZT );
 }
 
 // add in fns for other distributions in the future
 // this fn is too slow to use in ame_rep right now