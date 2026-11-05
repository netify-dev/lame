//Includes/namespaces
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

arma::mat mhalf_cpp2(
   arma::mat M
) {

  arma::vec eigVal;
  arma::mat eigVec;
  arma::eig_sym(eigVal, eigVec, M);

  arma::mat eigValDiagMat = pow(diagmat(eigVal), .5);
  arma::mat tmp = eigVec * eigValDiagMat * eigVec.t();

  return(tmp);
}

// Gibbs update for dyadic variance with independent replicate relational data
// Gibbs update for dyadic variance with independent replicate relational data
 //'
 //'
 //' @usage rs2_rep_fc_cpp(E.T, rhoMat=solve(matrix(c(1,rho,rho,1),2,2)) )
 //' @param E.T Array of square residual relational matrix series. The third
 //' dimension of the array is for different replicates. Each slice of the array
 //' according to the third dimension is a square residual relational matrix
 //' @param rhoMat inverted off-diagnoal matrix of rho
 //' @return a new value of s2
 //' @author Shahryar Minhas
 //' @keywords internal

 // [[Rcpp::export]]

 double rs2_rep_fc_cpp(
     arma::cube ET, arma::mat rhoMat, double s20 = 1.0,
     double beta_df = 0.0, double beta_ss = 0.0
 ) {

   // Guard: this function assumes square (unipartite) residual matrices
   if(ET.n_rows != ET.n_cols) {
     Rcpp::stop("rs2_rep_fc_cpp requires square residual matrices (unipartite). "
                "Got %d x %d.", ET.n_rows, ET.n_cols);
   }

   if (!(s20 > 0) || !std::isfinite(s20)) {
     Rcpp::stop("s20 must be a finite positive scalar.");
   }
   if (!(beta_df >= 0) || !std::isfinite(beta_df) ||
       !(beta_ss >= 0) || !std::isfinite(beta_ss)) {
     Rcpp::stop("Coefficient-prior degrees of freedom and sum of squares must be finite and non-negative.");
   }
   int N = ET.n_slices;
   arma::mat H = mhalf_cpp2(rhoMat);
   arma::mat tmp = trimatl(arma::ones(ET.n_rows, ET.n_cols));
   arma::uvec tmpindex = find(tmp==0);
   arma::mat EM(0, 2);
   for(int t=0 ; t<N ; ++t){
     arma::mat E = ET.slice(t); arma::mat Et = E.t();
     arma::vec eUpperTri = E.elem( tmpindex );
     arma::vec etUpperTri = Et.elem( tmpindex );
     EM = join_cols(EM, join_rows(eUpperTri, etUpperTri) * H);
   }

   // The model caller fixes this positive prior scale before sampling.
   // Recomputing it from current residuals makes conditionals incompatible.
   double ss = accu(pow(EM, 2));
   // force floating-point arithmetic: EM.n_elem is unsigned, so
   // (n_elem + 1)/2 truncates and loses the half pseudo-observation
   double shape = ((double) EM.n_elem + 1.0 + beta_df)/2.0;
   double scale = (ss + s20 + beta_ss)/2;
   double s2 = 1/rgamma(1, shape, 1/scale)[0];
   return( s2 );
 }
