// Kernel of the scale moves of R/ordinal_scale_move.R: the sums of the products of the latent array Z, its predictor E and the part M of
// the predictor that the scale does not carry, in the forms the metric of the covariance of the latent array is made of. Everything
// else the moves need (the residual sum of squares, the mean residual of a period, the quadratic and linear terms of the scale of a
// period) is arithmetic on these sums, so one pass over the cells serves the whole visit.
#include <Rcpp.h>
using namespace Rcpp;

// Z, E, M arrays nr x nc x nt in column-major order (M of length 0 stands for zeros); mode 1 directed (square), 2 symmetric (square),
// 3 bipartite. Row t of the result holds, for the pairs (Z,Z), (Z,E), (Z,M), (E,E), (E,M), (M,M) in turn, three sums, then the sums
// of Z, E and M (21 values):
//   directed   P = sum_{i<j} (x_ij + x_ji)(y_ij + y_ji), Q = sum_{i<j} (x_ij - x_ji)(y_ij - y_ji), D = sum_i x_ii y_ii; sums over all cells
//   symmetric  U = sum_{i<j} xs_ij ys_ij, 0, 0, with xs_ij = (x_ij + x_ji) / 2: the conditional mean of the latent value of a pair is the mean of
//              the two directed predictors, which differ when a row or column covariate is in the design; sums over the pairs
//   bipartite  S = sum x_ij y_ij, 0, 0; sums over all cells
// [[Rcpp::export]]
NumericMatrix osm_gram_cpp(const NumericVector& Z, const NumericVector& E, const NumericVector& M, const IntegerVector& dims, int mode) {
	const R_xlen_t nr = dims[0], nc = dims[1], nt = dims[2];
	const R_xlen_t cells = nr * nc;
	if (Z.size() != cells * nt || E.size() != cells * nt) Rcpp::stop("osm_gram_cpp: the latent array and the predictor must have the stated dimensions");
	const bool has_m = M.size() > 0;
	if (has_m && M.size() != cells * nt) Rcpp::stop("osm_gram_cpp: the unscaled part must have the stated dimensions");
	if (mode < 1 || mode > 3) Rcpp::stop("osm_gram_cpp: mode must be 1, 2 or 3");
	if (mode != 3 && nr != nc) Rcpp::stop("osm_gram_cpp: a square network needs a square array");
	NumericMatrix out(nt, 21);
	for (R_xlen_t t = 0; t < nt; ++t) {
		const double* z = &Z[0] + t * cells;
		const double* e = &E[0] + t * cells;
		const double* m = has_m ? &M[0] + t * cells : (const double*) 0;
		double s[6][3] = {{0}};
		double tot[3] = {0, 0, 0};
		if (mode == 3) {
			for (R_xlen_t k = 0; k < cells; ++k) {
				const double a = z[k], b = e[k], c = has_m ? m[k] : 0.0;
				s[0][0] += a * a; s[1][0] += a * b; s[2][0] += a * c; s[3][0] += b * b; s[4][0] += b * c; s[5][0] += c * c;
				tot[0] += a; tot[1] += b; tot[2] += c;
			}
		} else {
			for (R_xlen_t j = 0; j < nc; ++j) {
				for (R_xlen_t i = 0; i < j; ++i) {
					const R_xlen_t u = i + j * nr, l = j + i * nr;
					const double a1 = z[u], a2 = z[l], b1 = e[u], b2 = e[l], c1 = has_m ? m[u] : 0.0, c2 = has_m ? m[l] : 0.0;
					if (mode == 1) {
						const double ap = a1 + a2, bp = b1 + b2, cp = c1 + c2, aq = a1 - a2, bq = b1 - b2, cq = c1 - c2;
						s[0][0] += ap * ap; s[1][0] += ap * bp; s[2][0] += ap * cp; s[3][0] += bp * bp; s[4][0] += bp * cp; s[5][0] += cp * cp;
						s[0][1] += aq * aq; s[1][1] += aq * bq; s[2][1] += aq * cq; s[3][1] += bq * bq; s[4][1] += bq * cq; s[5][1] += cq * cq;
						tot[0] += ap; tot[1] += bp; tot[2] += cp;
					} else {
						// one latent value per unordered pair, around the mean of the two directed predictors
						const double as = 0.5 * (a1 + a2), bs = 0.5 * (b1 + b2), cs = 0.5 * (c1 + c2);
						s[0][0] += as * as; s[1][0] += as * bs; s[2][0] += as * cs; s[3][0] += bs * bs; s[4][0] += bs * cs; s[5][0] += cs * cs;
						tot[0] += as; tot[1] += bs; tot[2] += cs;
					}
				}
			}
			if (mode == 1) {
				for (R_xlen_t i = 0; i < nr; ++i) {
					const R_xlen_t d = i + i * nr;
					const double a = z[d], b = e[d], c = has_m ? m[d] : 0.0;
					s[0][2] += a * a; s[1][2] += a * b; s[2][2] += a * c; s[3][2] += b * b; s[4][2] += b * c; s[5][2] += c * c;
					tot[0] += a; tot[1] += b; tot[2] += c;
				}
			}
		}
		for (int p = 0; p < 6; ++p) for (int f = 0; f < 3; ++f) out(t, 3 * p + f) = s[p][f];
		for (int k = 0; k < 3; ++k) out(t, 18 + k) = tot[k];
	}
	return out;
}
