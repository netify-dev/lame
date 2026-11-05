// receiver-level move of the censored-binary (cbin) and fixed-rank (frn)
// nomination samplers, compiled: the sweep over the columns of
// .lame_nom_col_shift_r in R/nomination_actor_move.R, with the same visiting
// order, the same uniforms and the same arithmetic. a column of every period
// translates by one amount together with its receiver level: the amount is
// drawn from the prior of the shifted level cut to the interval the order
// constraints of the cells of that column leave (see R/nomination_actor_move.R
// for the derivation). the classes of the cells and the hidden-cell labels
// are prepared in R; the sweep keeps the row aggregates its bounds read (the
// strongest non-nominee and the weakest nominee of each row) current as
// columns move.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
using namespace Rcpp;

// one draw of N(ez, 1) on (lb, ub) from the uniform u, as rtnorm_interval_logp
static double nom_tn_draw(double ez, double lb, double ub, double u) {
	u = std::min(std::max(u, std::numeric_limits<double>::min()), 1 - 1e-16);
	if (!ISNAN(ub - ez) && (ub - ez) <= 0) {
		double la = R::pnorm(lb - ez, 0.0, 1.0, 1, 1);
		double lc = R::pnorm(ub - ez, 0.0, 1.0, 1, 1);
		double r = std::exp(la - lc);
		return ez + R::qnorm(lc + std::log(u + (1 - u) * r), 0.0, 1.0, 1, 1);
	}
	double lqa = R::pnorm(lb - ez, 0.0, 1.0, 0, 1);
	double lqb = R::pnorm(ub - ez, 0.0, 1.0, 0, 1);
	double r = std::exp(lqb - lqa);
	return ez + R::qnorm(lqa + std::log((1 - u) + u * r), 0.0, 1.0, 0, 1);
}

// kind 1 cbin, 2 frn without hidden cells, 3 frn with hidden cells
static inline bool nom_is_nominee(int kind, int cl) {
	return kind == 1 ? cl == 1 : (kind == 2 ? cl >= 1 : (cl >= 1 || cl == -2));
}

// [[Rcpp::export]]
List nom_col_shift_cpp(NumericVector Z, IntegerVector cls, LogicalMatrix flag,
		IntegerVector map, IntegerVector kind, IntegerVector mr,
		IntegerVector ord, NumericVector u, NumericVector b, NumericVector m,
		double s2, bool center, NumericVector beta, NumericVector iV0row) {
	IntegerVector dz = Z.attr("dim");
	const int nr = dz[0], nc = dz[1], T = dz[2];
	IntegerVector dm = map.attr("dim");
	const int mw = dm[1];
	NumericVector Zc = clone(Z);
	NumericVector bc = clone(b);
	NumericVector betac = clone(beta);
	const double inf = R_PosInf;
	const size_t snr = (size_t) nr, snc = (size_t) nc;
	auto zi = [&](int i, int j, int t) -> double& {
		return Zc[(size_t) i + snr * ((size_t) j + snc * (size_t) t)];
	};
	auto ci = [&](int i, int j, int t) -> int {
		return cls[(size_t) i + snr * ((size_t) j + snc * (size_t) t)];
	};
	// the column (0-based) holding rank / order position q (1-based) of row i, -1 if none
	auto mp = [&](int i, int q, int t) -> int {
		if (q < 1 || q > mw) return -1;
		int v = map[(size_t) i + snr * ((size_t) (q - 1) + (size_t) mw * (size_t) t)];
		return (v == NA_INTEGER || v < 1) ? -1 : v - 1;
	};
	// row aggregates: the strongest non-nominee, and the strongest of -z over the nominees
	std::vector<double> M0(snr * T, R_NegInf), Nq(snr * T, R_NegInf);
	auto rowmax0 = [&](int i, int t) {
		double mx = R_NegInf;
		for (int c = 0; c < nc; c++)
			if (ci(i, c, t) == 0 && zi(i, c, t) > mx) mx = zi(i, c, t);
		return mx;
	};
	auto rowmaxn = [&](int i, int t) {
		double mx = R_NegInf;
		const int kd = kind[t];
		for (int c = 0; c < nc; c++)
			if (nom_is_nominee(kd, ci(i, c, t)) && -zi(i, c, t) > mx) mx = -zi(i, c, t);
		return mx;
	};
	for (int t = 0; t < T; t++)
		for (int i = 0; i < nr; i++) {
			M0[(size_t) i + snr * t] = rowmax0(i, t);
			Nq[(size_t) i + snr * t] = rowmaxn(i, t);
		}
	double calls = 0, skips = 0;
	for (int k = 0; k < ord.size(); k++) {
		const int j = ord[k] - 1;
		double dlo = -inf, dhi = inf;
		for (int t = 0; t < T; t++) {
			const int kd = kind[t];
			for (int i = 0; i < nr; i++) {
				const int cl = ci(i, j, t);
				if (cl == -1) continue;
				const double M0i = M0[(size_t) i + snr * t], Nqi = Nq[(size_t) i + snr * t];
				double lo = -inf, hi = inf;
				if (kd == 1) {
					if (cl == 1) lo = std::max(0.0, M0i);
					else if (cl == 0) hi = flag(i, t) ? 0.0 : -Nqi;
				} else if (kd == 2) {
					if (cl == 0) {
						if (flag(i, t)) hi = 0.0;
						else {
							int c1 = mp(i, 1, t);
							hi = (c1 < 0) ? inf : zi(i, c1, t);
						}
					} else if (cl >= 1) {
						if (cl == 1) lo = std::max(0.0, M0i);
						else {
							int cp = mp(i, cl - 1, t);
							lo = (cp < 0) ? -inf : zi(i, cp, t);
						}
						if (cl < mr[t]) {
							int cn = mp(i, cl + 1, t);
							hi = (cn < 0) ? inf : zi(i, cn, t);
						}
					}
				} else {
					if (cl == 0) hi = flag(i, t) ? -Nqi : 0.0;
					else if (cl == -2) lo = std::max(0.0, M0i);
					else if (cl >= 1) {
						if (cl == 1) lo = std::max(0.0, M0i);
						else {
							int cp = mp(i, cl - 1, t);
							lo = (cp < 0) ? -inf : zi(i, cp, t);
						}
						if (cl < mr[t]) {
							int cn = mp(i, cl + 1, t);
							hi = (cn < 0) ? inf : zi(i, cn, t);
						}
					}
				}
				const double zj = zi(i, j, t);
				if (std::isfinite(lo)) dlo = std::max(dlo, lo - zj);
				if (std::isfinite(hi)) dhi = std::min(dhi, hi - zj);
			}
		}
		if (dlo > 0 || dhi < 0) {
			skips += 1;
			continue;
		}
		double delta;
		if (center) {
			double mb = 0;
			for (int c = 0; c < nc; c++) mb += bc[c] - m[c];
			mb /= nc;
			double ib = 0;
			for (int q = 0; q < betac.size(); q++) ib += iV0row[q] * betac[q];
			const double a2 = (nc - 1.0) / nc / s2 + iV0row[0] / ((double) nc * nc);
			const double a1 = ((bc[j] - m[j]) - mb) / s2 + ib / nc;
			const double sd_d = 1.0 / std::sqrt(a2);
			delta = sd_d * nom_tn_draw((-a1 / a2) / sd_d, dlo / sd_d, dhi / sd_d, u[k]);
		} else {
			const double sd_d = std::sqrt(s2);
			const double val = sd_d * nom_tn_draw(m[j] / sd_d, (bc[j] + dlo) / sd_d, (bc[j] + dhi) / sd_d, u[k]);
			delta = val - bc[j];
		}
		if (!std::isfinite(delta)) {
			skips += 1;
			continue;
		}
		calls += 1;
		for (int t = 0; t < T; t++) {
			const int kd = kind[t];
			for (int i = 0; i < nr; i++) {
				const int cl = ci(i, j, t);
				const double old = zi(i, j, t), nw = old + delta;
				zi(i, j, t) = nw;
				const size_t ai = (size_t) i + snr * t;
				if (cl == 0) {
					if (nw >= M0[ai]) M0[ai] = nw;
					else if (old >= M0[ai]) M0[ai] = rowmax0(i, t);
				}
				if (nom_is_nominee(kd, cl)) {
					if (-nw >= Nq[ai]) Nq[ai] = -nw;
					else if (-old >= Nq[ai]) Nq[ai] = rowmaxn(i, t);
				}
			}
		}
		if (center) {
			for (int c = 0; c < nc; c++) bc[c] -= delta / nc;
			bc[j] += delta;
			betac[0] += delta / nc;
		} else {
			bc[j] += delta;
		}
	}
	return List::create(Named("Z") = Zc, Named("b") = bc, Named("beta") = betac,
		Named("calls") = calls, Named("skips") = skips);
}
