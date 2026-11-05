// kernels of the ordinal position move (R/sparse_position_move.R): whether a shifted line of latent values is a
// possible latent line for the rank likelihood or the cutpoints, and the extremes of its cells by category
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ZM, WM: the cells of the line by periods (cells x N; WM the integer category, NA for a cell without one); omax,
// omin: the highest and lowest value of the other cells of each category in each period (K x N, -Inf / Inf where the
// category has no other cell); af: the cutpoints c(-Inf, alpha, Inf), empty for the rank likelihood
// [[Rcpp::export]]
bool spm_ord_ok_cpp(NumericMatrix ZM, IntegerMatrix WM, NumericMatrix omax, NumericMatrix omin, NumericVector af) {
	const int n = ZM.nrow(), N = ZM.ncol(), K = omax.nrow();
	if (af.size() > 0) {
		for (int t = 0; t < N; ++t) for (int a = 0; a < n; ++a) {
			int w = WM(a, t);
			if (w == NA_INTEGER) continue;
			double z = ZM(a, t);
			if (!(z > af[w - 1] && z < af[w])) return false;
		}
		return true;
	}
	if (K < 2) return true;
	const double inf = std::numeric_limits<double>::infinity();
	std::vector<double> mx(K), mn(K);
	for (int t = 0; t < N; ++t) {
		for (int k = 0; k < K; ++k) { mx[k] = omax(k, t); mn[k] = omin(k, t); }
		for (int a = 0; a < n; ++a) {
			int w = WM(a, t);
			if (w == NA_INTEGER) continue;
			double z = ZM(a, t);
			if (z > mx[w - 1]) mx[w - 1] = z;
			if (z < mn[w - 1]) mn[w - 1] = z;
		}
		double run = mx[0];
		for (int k = 1; k < K; ++k) {
			if (run >= mn[k]) return false;
			if (mx[k] > run) run = mx[k];
		}
		(void) inf;
	}
	return true;
}

// the extremes of the cells of a line by category and period, with their positions in the period's matrix; i the
// actor (1-based), row whether the line is the actor's row, nr the number of rows
// [[Rcpp::export]]
List spm_ord_line_cpp(NumericMatrix ZM, IntegerMatrix WM, int K, int i, bool row, int nr) {
	const int n = ZM.nrow(), N = ZM.ncol();
	NumericMatrix mmax(K, N), mmin(K, N);
	IntegerMatrix lm(K, N), lmn(K, N);
	std::fill(mmax.begin(), mmax.end(), NA_REAL);
	std::fill(mmin.begin(), mmin.end(), NA_REAL);
	std::fill(lm.begin(), lm.end(), NA_INTEGER);
	std::fill(lmn.begin(), lmn.end(), NA_INTEGER);
	for (int t = 0; t < N; ++t) for (int a = 0; a < n; ++a) {
		int w = WM(a, t);
		if (w == NA_INTEGER) continue;
		double z = ZM(a, t);
		int loc = row ? i + a * nr : (a + 1) + (i - 1) * nr;
		int k = w - 1;
		if (ISNAN(mmax(k, t)) || z > mmax(k, t)) { mmax(k, t) = z; lm(k, t) = loc; }
		if (ISNAN(mmin(k, t)) || z < mmin(k, t)) { mmin(k, t) = z; lmn(k, t) = loc; }
	}
	return List::create(Named("mmax") = mmax, Named("mmin") = mmin, Named("lm") = lm, Named("lmn") = lmn);
}
