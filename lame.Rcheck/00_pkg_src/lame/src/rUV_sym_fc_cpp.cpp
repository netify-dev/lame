//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

arma::mat matMultVec(arma::mat x, arma::vec y){
  int nRows = x.n_rows;
  int nCols = x.n_cols;
  
  arma::mat xy = arma::zeros(nRows,nCols);
  if(nRows <= nCols){
    for(int r = 0 ; r<nRows ; r++){
      xy.row(r) = x.row(r) % y.t();
    }		
  }
  if(nCols < nRows){
    for(int c = 0 ; c<nCols ; c++){
      xy.col(c) = x.col(c) % y;
    }		
  }		
  
  return(xy);
}

//' Gibbs sampling of U and V
 //' 
 //' A Gibbs sampler for updating the multiplicative effect matrices U and V
 //' in the symmetric case. In this case \code{U\%*\%t(V)} is symmetric, so
 //' this is parameterized as \code{V=U\%*\%L} where \code{L} is the 
 //' diagonal matrix of eigenvalues of \code{U\%*\%t(V)}. 
 //' 
 //' @usage rUV_sym_fc_cpp(E, U, V, s2 = 1, shrink=TRUE)
 //' @param E square residual relational matrix
 //' @param U current value of U
 //' @param V current value of V
 //' @param s2 dyadic variance
 //' @param shrink adaptively shrink the factors with a hierarchical prior
 //' @param uLoopIDs vector of rows to iteratively loop through
 //' when updating U
 //' @return \item{U}{a new value of U} \item{V}{a new value of V}
 //' @author Peter Hoff, Shahryar Minhas
 //' @examples
 //' 
 //' U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
 //' E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
 //' rUV_sym_fc_cpp(E,U0,V0,s2=1,shrink=TRUE,uLoopIDs=rep( sample(1:nrow(E)),4)-1 )
 //' 
 //' @export rUV_sym_fc_cpp
 // [[Rcpp::export]]
 
 List rUV_sym_fc_cpp(
     arma::mat E, arma::mat U, arma::mat V, 
     double s2, bool shrink, NumericVector uLoopIDs) {
   
   int R = U.n_cols; int n = U.n_rows;
   arma::mat L = diagmat( V.row(0)/U.row(0)  );
   L.replace(datum::nan, 0);
   
   arma::vec ivU = arma::zeros(R);
   arma::mat ivDiagMat;
   arma::rowvec scale;
   double shape;
   
   if(shrink){
     shape = (2+n)/2;
     scale = (1+sum(pow(U,2),0))/2;
     for(int r=0 ; r<R ; r++){
       ivU[r] = rgamma(1, shape, 1/scale[r] )[0];
     }
     ivDiagMat = diagmat(ivU);
   }
   
   double ugh = n; arma::vec tmp(R) ; tmp.fill( 1/ugh );
   if(!shrink){
     ivDiagMat = diagmat(tmp);
   }
   
   for(int s = 0; s<uLoopIDs.size() ; ++s){
     int i = uLoopIDs[s];
     arma::vec erow = E.row(i).t();
     arma::mat eui = matMultVec(U, erow);
     
     arma::rowvec euicolsum = sum(eui,0);
     
     arma::mat l = L * trans((euicolsum - U.row(i) * E(i,i))/s2);
     arma::mat iQ = inv( ivDiagMat + L * ( (U.t() * U) - (U.row(i).t() * U.row(i)) ) * L/s2 );
     arma::vec randNormDraw = rnorm(R);
     U.row(i) = trans(iQ * l + trans(chol(iQ)) * randNormDraw);
   }
   
   arma::mat tmponesmat = trimatl(arma::ones(E.n_rows, E.n_cols));
   arma::uvec tmpindex = find(tmponesmat==0);	
   
   for(int r = 0 ; r<R ; r++){
     arma::mat Usmall = U; Usmall.shed_col(r);
     arma::mat Lsmall = L; Lsmall.shed_col(r); Lsmall.shed_row(r);
     arma::mat Er = E - Usmall * Lsmall * Usmall.t();
     arma::mat lMat = Er % ( U.col(r) * U.col(r).t() );
     double l = accu(lMat.elem( tmpindex ))/s2;
     arma::mat uut2 = pow(U.col(r) * U.col(r).t(), 2);
     double iq = 1/(1+accu(uut2.elem(tmpindex))/s2);
     L(r,r) = rnorm(1, iq*l, pow(iq,.5) )[0];
   }
   
   return(
     Rcpp::List::create(
       Rcpp::Named("U")=U,
       Rcpp::Named("V")=U*L
     )		
   );	
 }
