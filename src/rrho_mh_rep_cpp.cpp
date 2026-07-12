//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// Metropolis update for dyadic correlation with independent replicate data
// Metropolis update for dyadic correlation with independent replicate data. 
 //' 
 //' 
 //' @usage rrho_mh_rep_cpp(E.T, rho, s2 = 1)
 //' @param E.T Array of square residual relational matrix series. The third
 //' dimension of the array is for different replicates. Each slice of the array
 //' according to the third dimension is a square residual relational matrix. 
 //' @param rho current value of rho
 //' @param s2 current value of s2
 //' @return a new value of rho
 //' @author Shahryar Minhas
 //' @keywords internal

 // [[Rcpp::export]]
 
 double rrho_mh_rep_cpp(
     arma::cube ET, double rho, double s2
 ) {
   
   int N = ET.n_slices;
   arma::mat tmp = trimatl(arma::ones(ET.n_rows, ET.n_cols));
   arma::uvec tmpindex = find(tmp==0);
   arma::mat EM(0, 2);
   
   for(int t=0 ; t<N ; ++t){
     arma::mat E = ET.slice(t); arma::mat Et = E.t();
     arma::vec eUpperTri = E.elem( tmpindex );
     arma::vec etUpperTri = Et.elem( tmpindex );
     EM = join_cols(EM, join_rows(eUpperTri, etUpperTri) / pow(s2, .5));
   }
   
   double emcp = accu( EM.col(0) % EM.col(1) );
   double emss = accu( pow(EM, 2) );
   
   int m = EM.n_rows;
   arma::mat emCor = arma::cor(EM);
   double sr = 2*( 1 - pow(emCor(0,1), 2) )/pow(m, .5);
   NumericVector x1; x1 = (-1-rho)/sr; NumericVector x2; x2 = (1-rho)/sr;
   
   double runiflo = Rcpp::pnorm(x1,0.0,1.0,1,0)[0];
   double runifhi = Rcpp::pnorm(x2,0.0,1.0,1,0)[0];
   
   if (runifhi <= runiflo || !R_finite(runiflo) || !R_finite(runifhi)) {
     return rho;
   }
   
   double runifdraw = R::runif(runiflo, runifhi);
   double qnormdraw = R::qnorm(runifdraw, 0.0, 1.0, 1, 0);
   double rho1 = rho + sr*qnormdraw;
   // truncated-normal proposal normalizer Z(r) = Phi((1-r)/sr) - Phi((-1-r)/sr).
   // The proposal is a normal truncated to (-1,1), so the Metropolis-Hastings
   // ratio must include q(rho|rho1)/q(rho1|rho) = Z(rho)/Z(rho1); omitting it
   // targets pi(rho)*Z(rho), biasing the chain away from |rho| -> 1.
   double lZ0 = std::log(R::pnorm((1-rho)/sr,0.0,1.0,1,0)  - R::pnorm((-1-rho)/sr,0.0,1.0,1,0));
   double lZ1 = std::log(R::pnorm((1-rho1)/sr,0.0,1.0,1,0) - R::pnorm((-1-rho1)/sr,0.0,1.0,1,0));
   double lhr = ( -.5*(m*log(1-pow(rho1,2))+(emss-2*rho1*emcp)/(1-pow(rho1,2))) ) -
     (-.5*(m*log(1-pow(rho,2) )+(emss-2*rho*emcp )/(1-pow(rho,2) ))) +
     ( (-.5*log(1-pow(rho1,2))) - (-.5*log(1-pow(rho,2))) )         +
     ( lZ0 - lZ1 );
   double rhoNew;
   if( log(runif(1,0,1)[0]) < lhr ){ rhoNew = rho1; } else { rhoNew = rho; }
   return( rhoNew );
 }