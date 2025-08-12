// Four-cycle counting for bipartite networks
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Count 4-cycles per slice in a (nA x nB x T) cube.
// A 4-cycle in bipartite is: i->j->k->l->i where i,k are in A and j,l are in B
// Counts via row-pair common neighbors: sum_{i<k} C(c_{ik}, 2) where c_{ik} = |N(i) ∩ N(k)|
// Treats any nonzero as an edge. Missing values should be pre-imputed to 0.
// [[Rcpp::export]]
arma::vec count_four_cycles_bip_cpp(const arma::cube& Y) {
  const unsigned int T = Y.n_slices;
  arma::vec out(T, fill::zeros);
  
  for (unsigned int t = 0; t < T; ++t) {
    // Binarize the network slice
    mat Yt = Y.slice(t);
    Yt.transform([](double v) { return (std::abs(v) > 0.0) ? 1.0 : 0.0; });
    
    // Row-pair common-neighbor counts: C = Y * Y^T
    mat C = Yt * Yt.t();  // nA x nA
    C.diag().zeros();     // Remove self-loops (i==k)
    
    // Extract strictly upper triangular (i < k)
    mat Cu = trimatu(C, 1);
    
    // Count 4-cycles: sum_{i<k} c_{ik} * (c_{ik} - 1) / 2
    mat term = (Cu % (Cu - 1.0)) * 0.5;
    out(t) = accu(term);
  }
  
  return out;
}

// Compute bipartite degree distributions
// Returns a list with row_degrees and col_degrees for each time slice
// [[Rcpp::export]]
List compute_degrees_bip_cpp(const arma::cube& Y) {
  const unsigned int nA = Y.n_rows;
  const unsigned int nB = Y.n_cols;
  const unsigned int T = Y.n_slices;
  
  arma::mat row_degrees(nA, T);
  arma::mat col_degrees(nB, T);
  
  for(unsigned int t = 0; t < T; ++t) {
    mat Yt = Y.slice(t);
    // Binarize
    Yt.transform([](double v) { return (std::abs(v) > 0.0) ? 1.0 : 0.0; });
    
    // Row degrees (out-degrees from A nodes)
    row_degrees.col(t) = sum(Yt, 1);
    
    // Column degrees (in-degrees to B nodes)
    col_degrees.col(t) = sum(Yt, 0).t();
  }
  
  return List::create(
    Named("row_degrees") = row_degrees,
    Named("col_degrees") = col_degrees
  );
}

// Compute bipartite clustering coefficient
// For bipartite networks, this is based on closed 4-cycles
// [[Rcpp::export]]
arma::vec clustering_coef_bip_cpp(const arma::cube& Y) {
  const unsigned int T = Y.n_slices;
  arma::vec out(T);
  
  for(unsigned int t = 0; t < T; ++t) {
    mat Yt = Y.slice(t);
    Yt.transform([](double v) { return (std::abs(v) > 0.0) ? 1.0 : 0.0; });
    
    // Compute common neighbors matrix
    mat C = Yt * Yt.t();  // nA x nA
    C.diag().zeros();
    
    // Count pairs with at least 2 common neighbors (potential for 4-cycles)
    double num_pairs = 0;
    double num_cycles = 0;
    
    for(unsigned int i = 0; i < C.n_rows; ++i) {
      for(unsigned int j = i+1; j < C.n_cols; ++j) {
        if(C(i,j) >= 2) {
          num_pairs += 1;
          num_cycles += C(i,j) * (C(i,j) - 1) / 2;
        }
      }
    }
    
    out(t) = (num_pairs > 0) ? num_cycles / num_pairs : 0;
  }
  
  return out;
}