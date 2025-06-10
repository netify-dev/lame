//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

//' Simulation from a Wishart distribution
 //' 
 //' Simulates a random Wishart-distributed matrix
 //' 
 //' 
 //' @usage rwish_cpp(S0, nu = dim(S0)[1] + 1)
 //' @param S0 a positive definite matrix
 //' @param nu a positive integer
 //' @return a positive definite matrix
 //' @author Peter Hoff, Shahryar Minhas
 //' @examples
 //' 
 //' ## The expectation is S0*nu
 //' 
 //' S0<-rwish(diag(3)) 
 //' SS<-matrix(0,3,3) 
 //' for(s in 1:1000) { SS<-SS+rwish(S0,5) }
 //' SS/s 
 //' S0*5
 //' 
 //' @export rwish_cpp
 // [[Rcpp::export]]
 
 arma::mat rwish_cpp(
     arma::mat S0, int nu
 ){
   arma::mat sS0 = chol(S0);
   arma::mat normMat = rnorm(nu * S0.n_rows) ; normMat.reshape(nu, S0.n_rows);
   arma::mat Z = normMat * sS0;
   arma::mat wishPull = Z.t() * Z;
   return(wishPull);
 }
 
 //' Gibbs sampling of U and V
 //' 
 //' A Gibbs sampler for updating the multiplicative effect matrices U and V,
 //' assuming they are the same across replicates. 
 //' 
 //' 
 //' @usage rUV_rep_fc_cpp(ET,U,V,rho,s2=1,maxmargin=1e-6,shrink=TRUE,rLoopIDs=sample(1:R)-1)
 //' @param ET Array of square residual relational matrix series with additive
 //' effects and covariates subtracted out. The third dimension of the array is
 //' for different replicates. Each slice of the array according to the third
 //' dimension is a square residual relational matrix. 
 //' @param U current value of U
 //' @param V current value of V
 //' @param rho dyadic correlation
 //' @param s2 dyadic variance
 //' @param iSe2 rho dyadic variance combo
 //' @param maxmargin small float
 //' @param shrink adaptively shrink the factors with a hierarchical prior
 //' @param rLoopIDs shuffled vector of latent space dimension sequence
 //' @return \item{U}{a new value of U} \item{V}{a new value of V}
 //' @author Peter Hoff, Yanjun He, Shahryar Minhas
 //' @export rUV_rep_fc_cpp
 // [[Rcpp::export]]
 
 List rUV_rep_fc_cpp(
     arma::cube ET, arma::mat U, arma::mat V, 
     double rho, double s2, arma::mat iSe2, 
     double maxmargin, bool shrink, NumericVector rLoopIDs){
   
   int Time = ET.n_slices; int R = U.n_cols; int n = U.n_rows;
   arma::mat UV = join_rows(U, V);
   
   arma::mat Suv; 
   if(shrink){
     Suv = inv( rwish_cpp(
       inv( eye<mat>(R*2,R*2) + UV.t()*UV ), n+R+2) );
   }
   if(!shrink){
     Suv = eye<mat>(R*2,R*2);
   }
   
   double g = iSe2(0,0); double d = iSe2(0,1);
   double g2d2 = pow(g,2) + pow(d,2);
   
   for(int i = 0; i<rLoopIDs.size(); i++){
     int r = rLoopIDs[i];
     arma::mat Usmall = U; Usmall.shed_col(r);
     arma::mat Vsmall = V; Vsmall.shed_col(r);
     arma::mat UVmr =  Usmall * Vsmall.t();
     arma::mat Est = arma::zeros(n * n, Time);
     for(int t=0 ; t<Time ; t++){
       arma::mat ert = ET.slice(t) - UVmr;
       Est.col(t) = vectorise(g2d2*ert + 2*g*d*ert.t());
     }
     
     // update u
     arma::mat vr = V.col(r);
     arma::mat Suvsmall = Suv.row(r); Suvsmall.shed_col(r);
     arma::mat Suvsmall2 = Suv; Suvsmall2.shed_col(r); Suvsmall2.shed_row(r);
     arma::mat b0 = Suvsmall * inv(Suvsmall2);
     arma::mat Suvsmall3 = Suv.col(r); Suvsmall3.shed_row(r) ;	
     arma::vec v0 = vectorise( Suv(r,r) - b0*Suvsmall3 );
     arma::mat m0 = join_rows(Usmall, V) * b0.t();
     double sumvr2 = accu(pow(vr,2)); double ssv;
     if( sumvr2 >= maxmargin ){ ssv=sumvr2; } else { ssv=maxmargin; }
     arma::vec a = Time*g2d2*ssv+1/v0 ; arma::vec c = -2*Time*g*d/(pow(a,2)+a*2*Time*g*d*ssv);
     arma::vec Esvvec = sum(Est, 1);
     int nEsv = pow(Esvvec.size(), .5);
     arma::mat Esvx; Esvx.insert_cols(0, Esvvec); Esvx.reshape(nEsv,nEsv);
     arma::mat Esv = Esvx * vr;
     arma::mat m1 = Esv/a[0] + c[0]*vr*accu( (Esv+m0/v0[0]) % vr) + m0/(a[0]*v0[0]); 
     arma::vec ah=pow(1/a, .5); arma::vec bh=(pow(1/a+ssv*c,.5)-pow(1/a,.5))/ssv;
     arma::vec e = rnorm(ET.n_rows);
     U.col(r) = m1 + ah[0]*e + bh[0] * vr * accu(vr % e);
     
     // update v
     arma::mat ur = U.col(r);
     int rv = R + r;
     arma::mat Suvsmall_v = Suv.row(rv); Suvsmall_v.shed_col(rv);
     arma::mat Suvsmall2_v = Suv; Suvsmall2_v.shed_col(rv); Suvsmall2_v.shed_row(rv);
     arma::mat b0_v = Suvsmall_v * inv(Suvsmall2_v);
     arma::mat Suvsmall3_v = Suv.col(rv); Suvsmall3_v.shed_row(rv);
     arma::vec v0_v = vectorise( Suv(rv,rv) - b0_v*Suvsmall3_v );
     arma::mat m0_v = join_rows(U, Vsmall) * b0_v.t();
     double sumur2 = accu(pow(ur,2)); double ssu;
     if( sumur2 >= maxmargin ){ ssu=sumur2; } else { ssu=maxmargin; }
     arma::vec a_v = Time*g2d2*ssu+1/v0_v ; arma::vec c_v = -2*Time*g*d/(pow(a_v,2)+a_v*2*Time*g*d*ssu);
     arma::mat tEsu = Esvx.t() * ur;
     arma::mat m1_v = tEsu/a_v[0] + c_v[0]*ur*accu( (tEsu+m0_v/v0_v[0]) % ur) + m0_v/(a_v[0]*v0_v[0]); 
     arma::vec ah_v=pow(1/a_v, .5); arma::vec bh_v=(pow(1/a_v+ssu*c_v,.5)-pow(1/a_v,.5))/ssu;
     arma::vec e_v = rnorm(ET.n_rows);
     V.col(r) = m1_v + ah_v[0]*e_v + bh_v[0] * ur * accu(ur % e_v);
   }
   
   return(
     Rcpp::List::create(
       Rcpp::Named("U")=U,
       Rcpp::Named("V")=V
     )
   );
   
 }