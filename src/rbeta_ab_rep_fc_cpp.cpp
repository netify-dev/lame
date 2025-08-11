//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

arma::mat mhalf_cpp(
    arma::mat M
) {
  
  arma::vec eigVal;
  arma::mat eigVec;
  arma::eig_sym(eigVal, eigVec, M);
  
  arma::mat eigValDiagMat = pow(diagmat(eigVal), .5);
  arma::mat tmp = eigVec * eigValDiagMat * eigVec.t();
  
  return(tmp);
}

arma::vec rmvnorm_cpp(
    arma::vec mu, arma::mat Sigma
) {
  int n = 1;
  arma::mat E = rnorm( n * mu.size() ) ; E.reshape(n, mu.size());
  arma::vec tmp = ( E * chol(Sigma) ).t();
  for(int i=0 ; i < tmp.n_rows ; i++){
    tmp.row(i) = tmp.row(i) + mu[i];
  }
  
  return( tmp );
}

// Gibbs sampling of additive row and column effects and regression coefficient
// with independent replicate relational data 
 //' Simulates from the joint full conditional distribution of (a,b,beta),
 //' assuming same additive row and column effects and regression coefficient
 //' across replicates. 
 //' 
 //' 
 //' @param ZT n x n x T array, with the third dimension for replicates.
 //' Each slice of the array is a (latent) normal relational matrix, with
 //' multiplicative effects subtracted out
 //' @param Xr n x p x T row covariate array generated within ame_repL fn
 //' @param Xc n x p x T column covariate array generated within ame_repL fn
 //' @param mX n^2 x p x T design array generated within ame_repL fn
 //' @param mXt n^2 x p x T transposed design array generated within ame_repL fn
 //' @param XX p x p x T regression sum of squares array generated within ame_repl fn
 //' @param XXt p x p x T transposed regression sum of squares array generated 
 //' within ame_repl fn
 //' @param iSe2 variance matrix
 //' @param Sabs row and column covariance
 //' @param k dimensions for row and column random effects
 //' @param G eigenvalue calcs from Sab
 //' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
 //' \item{b}{additive column effects}
 //' @author Peter Hoff, Yanjun He, Shahryar Minhas
 //' @keywords internal

 // [[Rcpp::export]]
 
 List rbeta_ab_rep_fc_cpp(
     arma::cube ZT, arma::cube Xr, arma::cube Xc, 
     arma::cube mX, arma::cube mXt, 
     arma::cube XX, arma::cube XXt,
     arma::mat iSe2, arma::mat Sabs, int k, arma::mat G
 ) {
   
   double td = iSe2(0,0);
   double to = iSe2(0,1);	
   
   const int N = ZT.n_slices;  
   const int n = ZT.n_rows;
   const int p = Xr.n_cols;
   
   arma::vec lb = arma::zeros(p);
   arma::mat Qb = arma::zeros(p,p);
   arma::vec ZrT = arma::zeros(n);
   arma::vec ZcT = arma::zeros(n);
   arma::mat XrT = arma::zeros(n,p);
   arma::mat XcT = arma::zeros(n,p);
   
   for(int t=0 ; t < N ; ++t){
     arma::mat Z = ZT.slice(t);
     arma::mat mXs = td * mX.slice(t) + to * mXt.slice(t);
     arma::mat XXs = (pow(to,2)+pow(td,2))*XX.slice(t) + 2*to*td*XXt.slice(t);
     arma::mat Zs = td*Z + to*trans(Z);
     arma::vec zr = sum(Zs,1); arma::vec zc=trans(sum(Zs,0));
     
     if(p>0){
       lb = lb + (trans(mXs) * vectorise(Zs));
       Qb = Qb + XXs + (XX.slice(t)/mXs.n_rows)/ZT.n_slices;
     }
     
     arma::mat Xsr = td*Xr.slice(t) + to*Xc.slice(t);
     arma::mat Xsc = td*Xc.slice(t) + to*Xr.slice(t);
     
     ZrT = ZrT+zr;
     ZcT = ZcT+zc;
     XrT = XrT+Xsr;
     XcT = XcT+Xsc;
   }
   
   // dyadic and prior contributions
   
   // row and column reduction
   arma::mat ab = arma::zeros(n,2);
   
   arma::mat K; arma::mat idmatk; arma::mat A; arma::mat B;
   arma::mat iA0; arma::mat C0; arma::mat iA; arma::mat C;
   arma::mat idmatn; arma::mat onesmatn; arma::mat H;
   arma::mat Hrr; arma::mat Hrc; arma::mat Hcr; arma::mat Hcc;
   if(k > 0){
     K = arma::mat("0 1; 1 0");
     idmatk = eye<mat>(k,k);
     
     A = N*n*G.t()*G + idmatk;
     B = N*G.t()*K*G;
     iA0 = inv(A);
     C0 = -inv(A + n*B)*B*iA0;
     
     iA = G * iA0 * G.t();
     C = G * C0 * G.t();
     
     idmatn = eye<mat>(n,n);
     onesmatn = ones<mat>(n,n);
     H = arma::kron(iA, idmatn) + arma::kron(C,onesmatn);
     Hrr = H.submat(0,0,n-1,n-1);
     Hrc = H.submat(0,n,n-1,2*n-1);
     Hcr = H.submat(n,0,2*n-1,n-1);
     Hcc = H.submat(n,n,2*n-1,2*n-1);
     Qb = Qb-XrT.t()*Hrr*XrT-XcT.t()*Hcr*XrT-XrT.t()*Hrc*XcT-XcT.t()*Hcc*XcT;
     lb = lb-XrT.t()*Hrr*ZrT-XcT.t()*Hcr*ZrT-XrT.t()*Hrc*ZcT-XcT.t()*Hcc*ZcT;
   }
   
   arma::mat Vb; arma::mat Mb; arma::vec beta;
   if(p > 0){
     Vb = inv(Qb); Mb = Vb * lb;
     beta = rmvnorm_cpp(Mb,Vb);
   }
   if(p==0){ beta = arma::zeros(1); }
   
   arma::mat RrT = ZrT - XrT * beta; 
   arma::mat RcT = ZcT - XcT * beta; 
   
   arma::mat RTcrossiA0G = join_rows(RrT, RcT) * (iA0 * G.t()).t();
   arma::vec RrTC0G = C0 * G.t() * accu(RrT);
   arma::mat m = arma::zeros(n,RTcrossiA0G.n_cols);
   for( int r=0 ; r < RTcrossiA0G.n_cols ; r++ ) {
     m.col(r) = RTcrossiA0G.col(r) + RrTC0G[r];
   }
   
   arma::mat hiA0 = mhalf_cpp(iA0);
   arma::mat e = rnorm( n*k ); e.reshape(n,k);
   arma::vec colE = sum(e, 0).t();
   arma::mat ehiA0 = (e * hiA0).t();
   arma::mat iA0nCo = mhalf_cpp(iA0+n*C0);
   arma::mat hiA0nCo = (hiA0-iA0nCo)/n;
   arma::vec ugh = hiA0nCo * colE;
   arma::mat w = arma::zeros(n,RTcrossiA0G.n_cols);
   for( int r=0 ; r < RTcrossiA0G.n_cols ; r++ ) {
     w.col(r) = m.col(r) + (ehiA0.row(r) - ugh[r]).t();
   }
   
   arma::mat abVec = w * G.t() * inv(iSe2);
   arma::vec a = abVec.col(0);
   arma::vec b = abVec.col(1);
   
   return(
     Rcpp::List::create(
       Rcpp::Named("beta")=beta,
       Rcpp::Named("a")=a,
       Rcpp::Named("b")=b
     )
   );
   
 }
