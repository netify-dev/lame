//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp; 

// [[Rcpp::export]]
void test_rbeta_inputs(
    arma::cube ZT, arma::cube Xr, arma::cube Xc, 
    arma::cube mX, arma::cube mXt, 
    arma::cube XX, arma::cube XXt,
    arma::mat iSe2, arma::mat Sabs, int k, arma::mat G,
    double g = 100.0
) {
  const int N = ZT.n_slices;  
  const int n = ZT.n_rows;
  const int p = Xr.n_cols;
  
  Rcpp::Rcout << "Dimensions: n=" << n << " p=" << p << " N=" << N << " k=" << k << std::endl;
  Rcpp::Rcout << "ZT: " << ZT.n_rows << "x" << ZT.n_cols << "x" << ZT.n_slices << std::endl;
  Rcpp::Rcout << "Xr: " << Xr.n_rows << "x" << Xr.n_cols << "x" << Xr.n_slices << std::endl;
  
  // Check for NaN in inputs
  if(ZT.has_nan()) Rcpp::Rcout << "WARNING: ZT has NaN" << std::endl;
  if(Xr.has_nan()) Rcpp::Rcout << "WARNING: Xr has NaN" << std::endl;
  if(Xc.has_nan()) Rcpp::Rcout << "WARNING: Xc has NaN" << std::endl;
  if(XX.has_nan()) Rcpp::Rcout << "WARNING: XX has NaN" << std::endl;
  
  // Test lb and Qb construction
  arma::vec lb = arma::zeros(p);
  arma::mat Qb = arma::zeros(p,p);
  
  double td = iSe2(0,0);
  double to = iSe2(0,1);
  
  for(int t=0; t<N; ++t){
    arma::mat mXs = td * mX.slice(t) + to * mXt.slice(t);
    arma::mat XXs = (to*to+td*td) * XX.slice(t) + 2*to*td*XXt.slice(t);
    arma::mat Zs = td * ZT.slice(t) + to * ZT.slice(t).t();
    
    // Check for NaN in intermediate values
    if(Zs.has_nan()) {
      Rcpp::Rcout << "  Zs has NaN at t=" << t << std::endl;
    }
    if(mXs.has_nan()) {
      Rcpp::Rcout << "  mXs has NaN at t=" << t << std::endl;
    }
    
    if(p>0){
      arma::vec lb_add = trans(mXs) * vectorise(Zs);
      if(lb_add.has_nan()) {
        Rcpp::Rcout << "  lb_add has NaN at t=" << t << std::endl;
        Rcpp::Rcout << "    vectorise(Zs) has NaN: " << vectorise(Zs).has_nan() << std::endl;
        Rcpp::Rcout << "    First non-nan Zs: " << Zs(0,0) << std::endl;
      }
      lb = lb + lb_add;
      Qb = Qb + XXs + (XX.slice(t)/g)/N;
    }
  }
  
  Rcpp::Rcout << "After loop:" << std::endl;
  Rcpp::Rcout << "  lb has NaN: " << lb.has_nan() << std::endl;
  Rcpp::Rcout << "  Qb has NaN: " << Qb.has_nan() << std::endl;
  if(!lb.has_nan()) {
    Rcpp::Rcout << "  lb first 3: " << lb(0) << " " << lb(1) << " " << lb(2) << std::endl;
  }
  if(!Qb.has_nan()) {
    Rcpp::Rcout << "  Qb(0,0): " << Qb(0,0) << std::endl;
  }
}