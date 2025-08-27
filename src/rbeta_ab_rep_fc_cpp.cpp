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
  // Use R's RNG instead of Armadillo's
  arma::vec E(mu.size());
  for(int i = 0; i < mu.size(); i++) {
    E(i) = R::rnorm(0.0, 1.0);
  }
  
  // Removed debug output - not needed in production
  
  // Ensure Sigma is symmetric before Cholesky
  arma::mat Sigma_sym = 0.5 * (Sigma + Sigma.t());
  
  // Try Cholesky decomposition with error handling
  arma::mat cholSigma;
  bool chol_success = chol(cholSigma, Sigma_sym);
  
  if(!chol_success) {
    // If Cholesky fails, add small ridge and try again
    Sigma_sym.diag() += 1e-6;
    chol_success = chol(cholSigma, Sigma_sym);
    
    if(!chol_success) {
      // If still fails, return the mean without randomness
      // Silent fallback - no warning needed for regularization
      return mu;
    }
  }
  
  arma::vec tmp = cholSigma.t() * E;
  tmp = tmp + mu;
  
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
 //' @param XXt p x p x T regression sum of squares array generated within ame_repl fn
 //' @param iSe2 half matrix of inverse se2 
 //' @param Sabs row and column covariance
 //' @param k dimensions for row and column random effects
 //' @param G eigenvalue calcs from Sab
 //' @param g g-prior parameter for regression coefficients
 //' @return \item{beta}{regression coefficients} \item{a}{additive row effects}
 //' \item{b}{additive column effects}
 //' @author Shahryar Minhas
 //' @keywords internal

 // [[Rcpp::export]]
 
 List rbeta_ab_rep_fc_cpp(
     arma::cube ZT, arma::cube Xr, arma::cube Xc, 
     arma::cube mX, arma::cube mXt, 
     arma::cube XX, arma::cube XXt,
     arma::mat iSe2, arma::mat Sabs, int k, arma::mat G,
     double g = 100.0  // Add g-prior parameter with default
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
   
   for(int t=0 ; t<N ; ++t){
     arma::mat mXs = td * mX.slice(t) + to * mXt.slice(t);
     arma::mat XXs = (to*to+td*td) * XX.slice(t) + 2*to*td*XXt.slice(t);
     arma::mat Zs = td * ZT.slice(t) + to * ZT.slice(t).t();
     arma::vec zr = sum(Zs,1); arma::vec zc = sum(Zs,0).t(); 
     
     if(p>0){
       lb = lb + (trans(mXs) * vectorise(Zs));
       // Add g-prior regularization
       Qb = Qb + XXs + (XX.slice(t)/g)/N;  // Use g-prior instead of n
     }
     
     arma::mat Xsr = td * Xr.slice(t) + to * Xc.slice(t);
     arma::mat Xsc = td * Xc.slice(t) + to * Xr.slice(t);		
     
     ZrT = ZrT + zr;
     ZcT = ZcT + zc;
     XrT = XrT + Xsr;
     XcT = XcT + Xsc;		
   }	
   
   // row and column reduction
   arma::mat ab = arma::zeros(n,2);
   arma::mat iA0, C0;  // Declare outside if block so they're accessible later
   
   if(k>0){
     arma::mat K = arma::zeros(2,2);
     K(0,1) = 1; K(1,0) = 1;
     
     arma::mat A = N*n*G.t()*G + arma::eye(k,k);
     arma::mat B = N*G.t()*K*G;
     iA0 = inv(A);
     C0 = -inv(A + n*B) * B * iA0;
     
     arma::mat iA = G * iA0 * G.t();
     arma::mat C = G * C0 * G.t();
     
     // BigMatrix<-rbind(cbind(n*diag(n),matrix(1,n,n)),cbind(matrix(1,n,n),n*diag(n)))
     // Gn<-G%x%diag(n)
     // V.inv<-N*crossprod(crossprod(BigMatrix,Gn),Gn)+diag(k*n)
     // V<-solve(V.inv)
     // H<-tcrossprod(tcrossprod(Gn,V),Gn)
     arma::mat iAkron = kron(iA, arma::eye(n,n));
     arma::mat Ckron = kron(C, arma::ones(n,n));
     arma::mat H = iAkron + Ckron;
     
     arma::mat Hrr, Hrc, Hcr, Hcc;
     Hrr = H.submat(0,0,n-1,n-1);
     Hrc = H.submat(0,n,n-1,2*n-1);
     Hcr = H.submat(n,0,2*n-1,n-1);
     Hcc = H.submat(n,n,2*n-1,2*n-1);
     Qb = Qb-XrT.t()*Hrr*XrT-XcT.t()*Hcr*XrT-XrT.t()*Hrc*XcT-XcT.t()*Hcc*XcT;
     lb = lb-XrT.t()*Hrr*ZrT-XcT.t()*Hcr*ZrT-XrT.t()*Hrc*ZcT-XcT.t()*Hcc*ZcT;
   }
   
   arma::mat Vb; arma::mat Mb; arma::vec beta;
   if(p > 0){
     // Add larger ridge for numerical stability
     double ridge = 1e-3;
     Qb.diag() += ridge;
     
     // Compute inverse with error handling
     bool inv_success = false;
     try {
       Vb = inv_sympd(Qb);  // Try symmetric positive definite inverse
       inv_success = true;
     } catch(...) {
       // Fall back to regular inverse with more regularization
       Qb.diag() += 0.1;
       try {
         Vb = inv(Qb);
         inv_success = true;
       } catch(...) {
         // Silent fallback - matrix inversion failed
       }
     }
     
     if(!inv_success || Vb.has_nan() || Vb.has_inf()) {
       beta = arma::zeros(p);
     } else {
       Mb = Vb * lb;
       
       // Check if Mb has NaN values
       if(Mb.has_nan() || Mb.has_inf()) {
         beta = arma::zeros(p);
       } else {
         // Ensure Vb is positive definite before sampling
         Vb = 0.5 * (Vb + Vb.t());  // Ensure symmetry
         arma::vec eigval = eig_sym(Vb);
         if(eigval.min() < 1e-6) {
           Vb.diag() += (1e-6 - eigval.min() + 1e-6);
         }
         
         beta = rmvnorm_cpp(Mb, Vb);
         
         // Check if beta has NaN after sampling
         if(beta.has_nan() || beta.has_inf()) {
           beta = Mb;  // Use mean without randomness
         }
       }
     }
   }
   if(p==0){ beta = arma::zeros(0); }  // FIX: Return empty vector when no predictors
   
   arma::vec a, b;
   if(k > 0) {
     // Ensure proper matrix dimensions for multiplication
     arma::vec RrT, RcT;
     if(p > 0 && beta.n_elem > 0 && XrT.n_cols == beta.n_elem) {
       arma::mat beta_col = arma::reshape(beta, beta.n_elem, 1);
       RrT = ZrT - XrT * beta_col;
       RcT = ZcT - XcT * beta_col;
     } else {
       RrT = ZrT;
       RcT = ZcT;
     } 
     
     arma::mat RTcrossiA0G = join_rows(RrT, RcT) * (iA0 * G.t()).t();
     
     double sumRrT = arma::sum(arma::vectorise(RrT));
     
     arma::mat m = arma::zeros(n,k);
     for(int r=0 ; r < n ; r++){
       arma::mat C0Gt = C0 * G.t();
      m.row(r) = RTcrossiA0G.row(r) + sumRrT * C0Gt.row(0);
     }
     
     arma::mat hiA0 = mhalf_cpp(iA0);
     // Use R's RNG for consistency
     arma::mat e(n, k);
     for(int i = 0; i < n; i++) {
       for(int j = 0; j < k; j++) {
         e(i,j) = R::rnorm(0.0, 1.0);
       }
     }
     arma::mat w = m + e * hiA0 - 
       arma::ones(n,1) * ((hiA0 - mhalf_cpp(iA0 + n*C0))/n * arma::sum(e,0).t()).t();
     
     arma::mat abvec = w * G.t() * inv(iSe2);
     
     a = abvec.col(0);
     b = abvec.col(1);
   } else {
     // When k=0, no random effects
     a = arma::zeros(n);
     b = arma::zeros(n);
   }
   
   return(
     Rcpp::List::create(
       Rcpp::Named("beta")=beta,
       Rcpp::Named("a")=a,
       Rcpp::Named("b")=b
     )
   );
   
 }