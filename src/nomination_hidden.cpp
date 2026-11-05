// hidden cells of the ranked-nomination rows (cbin, and frn read as its
// nominations): the implied labels and the per-cell conditional draws of
// R/rZ_nomination_hidden.R, compiled. the arithmetic and the order of the
// uniform draws follow the R version it replaces, so a fixed seed gives the
// same draws: sums that R forms in long double (rowSums, cumsum) are formed
// in long double here too.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
using namespace Rcpp;

// the nomination labels the current latent values imply for the hidden
// cells: the highest positive hidden cells of each row, as many as the cap
// leaves after the row's observed nominations
// [[Rcpp::export]]
NumericMatrix cbin_implied_labels_cpp(NumericMatrix Z, NumericMatrix Y,
		NumericVector cap, LogicalMatrix hid) {
	int nr = Y.nrow(), nc = Y.ncol();
	NumericMatrix Ya = clone(Y);
	for (int i = 0; i < nr; i++) {
		int n1 = 0, nh = 0;
		for (int j = 0; j < nc; j++) {
			double y = Y(i, j);
			if (!ISNAN(y) && y > 0) n1++;
			if (hid(i, j)) nh++;
		}
		if (nh == 0) continue;
		double r = std::max(cap[i] - n1, 0.0);
		std::vector<int> pos;
		for (int j = 0; j < nc; j++) {
			if (hid(i, j)) {
				Ya(i, j) = 0.0;
				if (Z(i, j) > 0) pos.push_back(j);
			}
		}
		if (r > 0 && !pos.empty()) {
			// the top min(r, #positive) positive hidden cells, by value
			int take = (int) std::min((double) pos.size(), r);
			std::vector<bool> used(pos.size(), false);
			for (int t = 0; t < take; t++) {
				int best = -1;
				for (size_t q = 0; q < pos.size(); q++) {
					if (used[q]) continue;
					if (best < 0 || Z(i, pos[q]) > Z(i, pos[best])) best = (int) q;
				}
				used[best] = true;
				Ya(i, pos[best]) = 1.0;
			}
		}
	}
	return Ya;
}

// log of the standard normal mass on (lo, hi), as .lame_log_norm_mass
static double log_norm_mass(double lo, double hi) {
	if (!(hi > lo)) return R_NegInf;
	if (lo >= 0) {
		double a = R::pnorm(lo, 0.0, 1.0, 0, 1);
		double b = R::pnorm(hi, 0.0, 1.0, 0, 1);
		return a + std::log1p(-std::exp(b - a));
	}
	double a = R::pnorm(hi, 0.0, 1.0, 1, 1);
	double b = R::pnorm(lo, 0.0, 1.0, 1, 1);
	return a + std::log1p(-std::exp(b - a));
}

// one standard normal draw on (lb, ub) given its uniform, as
// rtnorm_interval_logp with ez = 0
static double tn_draw(double lb, double ub, double u) {
	u = std::min(std::max(u, std::numeric_limits<double>::min()), 1 - 1e-16);
	if (!ISNAN(ub) && ub <= 0) {
		double la = R::pnorm(lb, 0.0, 1.0, 1, 1);
		double lc = R::pnorm(ub, 0.0, 1.0, 1, 1);
		double r = std::exp(la - lc);
		return 0.0 + R::qnorm(lc + std::log(u + (1 - u) * r), 0.0, 1.0, 1, 1);
	}
	double lqa = R::pnorm(lb, 0.0, 1.0, 0, 1);
	double lqb = R::pnorm(ub, 0.0, 1.0, 0, 1);
	double r = std::exp(lqb - lqa);
	return 0.0 + R::qnorm(lqa + std::log((1 - u) + u * r), 0.0, 1.0, 0, 1);
}

// each hidden cell redrawn from its conditional given the rest of its row
// and the observed cells (see .lame_cbin_hidden_moves in
// R/rZ_nomination_hidden.R for the condition). hl holds each row's hidden
// columns (1-based) in the random order drawn in R; one hidden cell per row
// per round, the upper triangle before the lower one when sq
// [[Rcpp::export]]
NumericMatrix cbin_hidden_moves_cpp(NumericMatrix Z, NumericMatrix EZ,
		double rho, bool sq, NumericVector cap, NumericMatrix Y,
		LogicalMatrix hid, List hl) {
	int nr = Z.nrow(), nc = Z.ncol();
	NumericMatrix Zn = clone(Z);
	double cond_sd = sq ? std::sqrt(1 - rho * rho) : 1.0;
	std::vector<double> m1(nr, R_PosInf), M0(nr, R_NegInf), r(nr);
	std::vector<int> n1(nr, 0), P0(nr, 0);
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < nc; j++) {
			if (sq && i == j) continue;
			double y = Y(i, j);
			if (ISNAN(y)) continue;
			if (y > 0) {
				n1[i]++;
				if (Zn(i, j) < m1[i]) m1[i] = Zn(i, j);
			} else if (y == 0) {
				if (Zn(i, j) > M0[i]) M0[i] = Zn(i, j);
				if (Zn(i, j) > 0) P0[i]++;
			}
		}
		r[i] = cap[i] - n1[i];
	}
	std::vector<IntegerVector> H(nr);
	int n_round = 0;
	for (int i = 0; i < nr; i++) {
		H[i] = as<IntegerVector>(hl[i]);
		if (H[i].size() > n_round) n_round = H[i].size();
	}
	if (n_round == 0) return Zn;
	int nblock = sq ? 2 : 1;
	std::vector<int> bi, bj;
	std::vector<double> bmu, bx, blo, bhi;
	std::vector<double> blm;
	for (int k = 0; k < n_round; k++) {
		for (int blk = 0; blk < nblock; blk++) {
			bi.clear(); bj.clear(); bmu.clear(); bx.clear(); blm.clear(); blo.clear(); bhi.clear();
			for (int i = 0; i < nr; i++) {
				if (H[i].size() <= k) continue;
				int j = H[i][k] - 1;
				if (sq) {
					if (blk == 0 && !(i < j)) continue;
					if (blk == 1 && !(i > j)) continue;
				}
				double x_old = Zn(i, j);
				double M = M0[i], m = m1[i];
				// the other hidden cells of the row above 0, M0 and m1
				int a0 = 0, aM = 0, am = 0;
				for (int hq = 0; hq < H[i].size(); hq++) {
					double zc = Zn(i, H[i][hq] - 1);
					if (zc > 0) a0++;
					if (zc > M) aM++;
					if (zc > m) am++;
				}
				if (x_old > 0) a0--;
				if (x_old > M) aM--;
				if (x_old > m) am--;
				double Mp = M > 0 ? M : 0.0;
				if (ISNAN(M)) Mp = 0.0;
				double lo[4] = {R_NegInf, 0.0, Mp, m};
				double hi[4] = {0.0, Mp, m, R_PosInf};
				double mu = EZ(i, j);
				if (sq && rho != 0) mu = mu + rho * (Zn(j, i) - EZ(j, i));
				double lm[4];
				for (int g = 0; g < 4; g++) {
					bool gt0 = g >= 1;
					bool gtM = g >= 2 && M > 0;
					bool gtm = g == 3;
					double K = n1[i] + a0 + P0[i] + (gt0 ? 1 : 0);
					bool ok_le = K <= cap[i] && P0[i] == 0;
					bool ok_gt = K > cap[i] && (am + (gtm ? 1 : 0)) <= r[i] &&
						(M <= 0 || (aM + (gtM ? 1 : 0)) >= r[i]);
					bool feas = (ok_le || ok_gt) && hi[g] > lo[g];
					double v = log_norm_mass((lo[g] - mu) / cond_sd, (hi[g] - mu) / cond_sd);
					lm[g] = feas ? v : R_NegInf;
				}
				// as R's max: a NaN weight makes the maximum NaN and the cell is
				// left where it is
				bool has_nan = false;
				double mx = lm[0];
				for (int g = 0; g < 4; g++) {
					if (ISNAN(lm[g])) has_nan = true;
					if (lm[g] > mx) mx = lm[g];
				}
				if (has_nan || !R_FINITE(mx)) continue;
				bi.push_back(i); bj.push_back(j); bmu.push_back(mu); bx.push_back(x_old);
				for (int g = 0; g < 4; g++) {
					blm.push_back(lm[g] - mx);
					blo.push_back(lo[g]); bhi.push_back(hi[g]);
				}
			}
			int ng = bi.size();
			if (ng == 0) continue;
			// region: u against the cumulated normalised weights
			std::vector<int> pick(ng);
			std::vector<double> w(4 * ng);
			for (int q = 0; q < ng; q++) {
				long double s = 0.0L;
				for (int g = 0; g < 4; g++) { w[4 * q + g] = std::exp(blm[4 * q + g]); s += w[4 * q + g]; }
				double sd = (double) s;
				for (int g = 0; g < 4; g++) w[4 * q + g] = w[4 * q + g] / sd;
			}
			for (int q = 0; q < ng; q++) {
				double u = unif_rand();
				long double cs = 0.0L;
				int cnt = 0;
				for (int g = 0; g < 4; g++) {
					cs += w[4 * q + g];
					if (u > (double) cs) cnt++;
				}
				int gp = 1 + cnt;
				if (gp > 4) gp = 4;
				pick[q] = gp - 1;
			}
			std::vector<double> uu(ng);
			for (int q = 0; q < ng; q++) uu[q] = unif_rand();
			for (int q = 0; q < ng; q++) {
				int g = pick[q];
				double lb = blo[4 * q + g], ub = bhi[4 * q + g];
				double z = tn_draw((lb - bmu[q]) / cond_sd, (ub - bmu[q]) / cond_sd, uu[q]);
				double val = bmu[q] + cond_sd * z;
				if (!R_FINITE(val)) val = bx[q];
				Zn(bi[q], bj[q]) = val;
			}
		}
	}
	return Zn;
}
