// The row check of .lame_nom_feasible (R/nomination_feasible.R): which rows of a latent slice break the constraints of the nomination
// samplers. Z latent slice, Y the nominations (NA missing), YL the column of each rank (frn; a 0 x 0 matrix for cbin), room the rows with
// unused nomination slots, bip whether the slice is bipartite (a unipartite diagonal holds no cell). A row with hidden cells has its room
// from them as well, so the room rule is not applied to it (see the R code).
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector nom_bad_rows_cpp(NumericMatrix Z, NumericMatrix Y, IntegerMatrix YL, LogicalVector room, bool bip) {
	const int nr = Z.nrow(), nc = Z.ncol();
	const bool frn = YL.nrow() > 0;
	LogicalVector bad(nr);
	std::vector<int> colobs(nc, 0), rowobs(nr, 0);
	for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i) {
		if (!bip && i == j) continue;
		if (!std::isnan(Y(i, j))) { ++colobs[j]; ++rowobs[i]; }
	}
	for (int i = 0; i < nr; ++i) {
		double mn = R_PosInf, mx = R_NegInf;
		bool hid = false;
		for (int j = 0; j < nc; ++j) {
			const double y = Y(i, j);
			if (std::isnan(y)) {
				if (!(!bip && i == j) && colobs[j] > 0 && rowobs[i] > 0) hid = true;
				continue;
			}
			if (y > 0) { if (Z(i, j) < mn) mn = Z(i, j); }
			else if (y == 0) { if (Z(i, j) > mx) mx = Z(i, j); }
		}
		if (!std::isfinite(mn)) continue;
		const bool room_chk = room[i] == TRUE && !hid;
		bool b = mn <= 0 || (std::isfinite(mx) && mx >= mn) || (room_chk && std::isfinite(mx) && mx >= 0);
		if (!b && frn) {
			double prev = 0;
			bool have = false;
			for (int q = 0; q < YL.ncol(); ++q) {
				const int c = YL(i, q);
				if (c == NA_INTEGER) continue;
				const double z = Z(i, c - 1);
				if (have && z - prev <= 0) { b = true; break; }
				prev = z;
				have = true;
			}
		}
		bad[i] = b;
	}
	return bad;
}
