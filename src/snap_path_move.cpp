// compiled visit of the snap path move of a dynamic_uv fit (R/dynamic_position_move.R,
// .lame_dpm_snap_visit). every sender (row) path and then every receiver (column) path, each side in a random
// order (a symmetric panel has one side),
// takes five exact steps given the rest of the state:
//   1. an elliptical slice step of the path with its resets summed out, from a Gaussian reference prior
//      (the AR(1) path, or the path that resets at one period drawn at random);
//   2. Metropolis steps on the scale of the path from a period drawn at random on;
//   3. the resets given the path;
//   4. an elliptical slice step of the path given its resets;
//   5. the resets given the path again.
// the likelihood factor of a line is a probit factor per observed cell (binary: the latent array integrated
// out given the shared dyad factor), the Poisson factor of its cells (the latent array moving with the
// predictor, so no residual changes) or the Gaussian density of its latent cells given their reciprocal
// partners (normal, ordinal, cbin and frn, the latent array held)
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

// the cells of one actor's line over the periods
struct Line {
	int rk = 0;
	std::vector<int> t;
	std::vector<double> b, s, y, e, w, v;
};

// hyperparameters of the snap prior of a path
struct Snap {
	double rho, s2, k2, lp, lq;
};

// log densities of the drift (a) and the reset (b) of transition t of the path x (rk by nt, by column)
inline void snap_ab(const std::vector<double>& x, int rk, int t, const Snap& p, double& a, double& b) {
	double sa = 0, sb = 0;
	for (int r = 0; r < rk; ++r) {
		const double y = x[r + rk * t];
		const double d = y - p.rho * x[r + rk * (t - 1)];
		sa += d * d;
		sb += y * y;
	}
	a = p.lq - sa / (2 * p.s2) - 0.5 * rk * std::log(p.s2);
	b = p.lp - sb / (2 * p.k2) - 0.5 * rk * std::log(p.k2);
}

inline double lse(double a, double b) {
	const double m = a > b ? a : b;
	return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// log of the reset-summed transitions over those of the reference prior (reset at t_ref when t_ref > 0)
double snap_ratio(const std::vector<double>& x, int rk, int nt, const Snap& p, int t_ref) {
	double out = 0;
	for (int t = 1; t < nt; ++t) {
		double a, b;
		snap_ab(x, rk, t, p, a, b);
		out += lse(a, b) - (t == t_ref ? b - p.lp : a - p.lq);
	}
	return out;
}

// log of the snap prior of a path with its resets summed out, up to a constant
double snap_log_prior(const std::vector<double>& x, int rk, int nt, const Snap& p) {
	const double v1 = p.s2 / (1 - p.rho * p.rho);
	double out = -0.5 * rk * std::log(v1);
	for (int r = 0; r < rk; ++r) out -= x[r] * x[r] / (2 * v1);
	for (int t = 1; t < nt; ++t) {
		double a, b;
		snap_ab(x, rk, t, p, a, b);
		out += lse(a, b);
	}
	return out;
}

// lower Cholesky root (nt by nt, by column) of the covariance of a coordinate path given the resets dl
bool snap_root(int nt, const Snap& p, const std::vector<int>& dl, std::vector<double>& L) {
	std::vector<double> v(nt), cf(nt, 0), C(nt * nt, 0);
	v[0] = p.s2 / (1 - p.rho * p.rho);
	for (int t = 1; t < nt; ++t) {
		cf[t] = dl[t] ? 0 : p.rho;
		v[t] = dl[t] ? p.k2 : p.rho * p.rho * v[t - 1] + p.s2;
	}
	for (int t = 0; t < nt; ++t) {
		C[t + nt * t] = v[t];
		double f = 1;
		for (int s = t - 1; s >= 0; --s) {
			f *= cf[s + 1];
			C[s + nt * t] = C[t + nt * s] = f * v[s];
		}
	}
	L.assign(nt * nt, 0);
	for (int j = 0; j < nt; ++j) {
		double d = C[j + nt * j];
		for (int k = 0; k < j; ++k) d -= L[j + nt * k] * L[j + nt * k];
		if (!(d > 0)) return false;
		L[j + nt * j] = std::sqrt(d);
		for (int i = j + 1; i < nt; ++i) {
			double s = C[i + nt * j];
			for (int k = 0; k < j; ++k) s -= L[i + nt * k] * L[j + nt * k];
			L[i + nt * j] = s / L[j + nt * j];
		}
	}
	return true;
}

// log likelihood factor of a line when its path moves by d (rk by nt); `kind` 0 probit, 1 Poisson, 2 Gaussian
double line_ll(const Line& ln, const std::vector<double>& d, int kind, double sd_cell) {
	const int rk = ln.rk;
	double out = 0;
	for (size_t c = 0; c < ln.t.size(); ++c) {
		double dl = 0;
		for (int r = 0; r < rk; ++r) dl += ln.w[c * rk + r] * d[r + rk * ln.t[c]];
		if (kind == 1) out += ln.y[c] * dl - ln.e[c] * std::expm1(dl);
		else if (kind == 2) out -= (ln.b[c] - dl) * (ln.b[c] - dl) / (2 * ln.v[c]);
		else out += R::pnorm(ln.s[c] * (ln.b[c] + dl) / sd_cell, 0, 1, 1, 1);
	}
	return out;
}

// the resets of a path given the path (first period ineligible)
void snap_resets(const std::vector<double>& x, int rk, int nt, const Snap& p, std::vector<int>& dl) {
	dl.assign(nt, 0);
	for (int t = 1; t < nt; ++t) {
		double a, b;
		snap_ab(x, rk, t, p, a, b);
		dl[t] = unif_rand() < 1 / (1 + std::exp(a - b)) ? 1 : 0;
	}
}

}

// one visit. EZ, Z: latent mean and array (nr x nc x N; Z moves with the predictor under the Poisson factor); off:
// shared dyad factor of a probit cell (directed only); mask: +1 event, -1 non-event, 0 no datum (probit) or 1, 0
// (Poisson, Gaussian); Y, offset: counts and period offsets of the Poisson factor; U, V: paths (nr x RA x N,
// nc x RB x N); G: interaction matrix of a bipartite product U_t G V_t'; Ld: weights of a symmetric product
// U_t diag(Ld) U_t', whose V is not read; delta_u, delta_v: resets (actors x N); sd_cell: probit cell sd given
// the dyad factor; s2, rho_t: latent variance and per-period dyadic correlation of the Gaussian factor (directed
// only); max_abs: coordinate rail (non-finite: none); kind: 0 probit, 1 Poisson, 2 Gaussian; mode: 0 directed,
// 1 bipartite, 2 symmetric. a symmetric line holds the cells (i, j), j != i, of row i, each pair once; its probit
// factor reads the pair mean of the two triangles, and the move shifts both triangles and the diagonal cell
static List snap_visit(
	NumericVector EZ,
	NumericVector Z,
	NumericVector off,
	NumericVector mask,
	NumericVector Y,
	NumericVector offset,
	NumericVector U,
	NumericVector V,
	NumericMatrix G,
	NumericVector Ld,
	double rho,
	double sigma,
	double kappa,
	double pi,
	NumericMatrix delta_u,
	NumericMatrix delta_v,
	double sd_cell,
	double max_abs,
	double s2,
	NumericVector rho_t,
	int kind,
	int mode
) {
	const bool poisson = kind == 1;
	const bool bip = mode == 1, sym = mode == 2;
	IntegerVector dz = EZ.attr("dim");
	IntegerVector du_dim = U.attr("dim");
	const int nr = dz[0], nc = dz[1], nt = dz[2];
	const int rku = du_dim[1];
	int rkv = rku;
	if (!sym) {
		IntegerVector dv_dim = V.attr("dim");
		rkv = dv_dim[1];
	}
	NumericVector ez = clone(EZ), zz = clone(Z), uu = clone(U), vv = clone(V);
	NumericMatrix du = clone(delta_u), dv = clone(delta_v);
	double calls = 0, evals = 0, stuck = 0, scales = 0, scales_acc = 0;
	const Snap p = {rho, sigma * sigma, kappa * kappa, std::log(pi), std::log1p(-pi)};
	const bool rail = std::isfinite(max_abs) && max_abs > 0;
	auto cix = [nr, nc](int i, int j, int t) { return (size_t)i + (size_t)nr * j + (size_t)nr * nc * t; };
	auto pix_u = [nr, rku](int i, int r, int t) { return (size_t)i + (size_t)nr * r + (size_t)nr * rku * t; };
	auto pix_v = [nc, rkv](int i, int r, int t) { return (size_t)i + (size_t)nc * r + (size_t)nc * rkv * t; };
	int rk = rku, dim = rku * nt;
	std::vector<double> cur, x, d, nu, z, prop, L, wv;
	std::vector<int> dl, dl_ref(nt);
	auto in_rail = [&](const std::vector<double>& q) {
		if (!rail) return true;
		for (int k = 0; k < dim; ++k) if (std::fabs(q[k]) > max_abs) return false;
		return true;
	};
	for (int side = 0; side < (sym ? 1 : 2); ++side) {
		const bool row = side == 0;
		NumericVector& own = row ? uu : vv;
		const int n_own = row ? nr : nc, n_line = row ? nc : nr;
		rk = row ? rku : rkv;
		dim = rk * nt;
		cur.assign(dim, 0); x.assign(dim, 0); d.assign(dim, 0); nu.assign(dim, 0); z.assign(dim, 0); prop.assign(dim, 0);
		wv.assign(rk, 0);
		auto own_ix = [&](int i, int r, int t) { return row ? pix_u(i, r, t) : pix_v(i, r, t); };
		// partner weights of cell (i, j) of actor i's line at period t, into wv
		auto weights = [&](int i, int j, int t) {
			if (sym) {
				for (int r = 0; r < rk; ++r) wv[r] = Ld[r] * uu[pix_u(j, r, t)];
			} else if (bip && row) {
				for (int r = 0; r < rk; ++r) {
					double s = 0;
					for (int q = 0; q < rkv; ++q) s += G(r, q) * vv[pix_v(j, q, t)];
					wv[r] = s;
				}
			} else if (bip) {
				for (int r = 0; r < rk; ++r) {
					double s = 0;
					for (int q = 0; q < rku; ++q) s += uu[pix_u(j, q, t)] * G(q, r);
					wv[r] = s;
				}
			} else {
				const NumericVector& partner = row ? vv : uu;
				for (int r = 0; r < rk; ++r) wv[r] = partner[pix_u(j, r, t)];
			}
		};
		// a random order of the actors of the side
		std::vector<int> ord(n_own);
		for (int k = 0; k < n_own; ++k) ord[k] = k;
		for (int k = n_own - 1; k > 0; --k) {
			const int m = (int)(unif_rand() * (k + 1));
			std::swap(ord[k], ord[m < k + 1 ? m : k]);
		}
		for (int oi = 0; oi < n_own; ++oi) {
			const int i = ord[oi];
			Line ln;
			ln.rk = rk;
			for (int t = 0; t < nt; ++t) for (int j = 0; j < n_line; ++j) {
				const size_t c = row ? cix(i, j, t) : cix(j, i, t);
				if ((!bip && j == i) || mask[c] == 0) continue;
				ln.t.push_back(t);
				weights(i, j, t);
				for (int r = 0; r < rk; ++r) ln.w.push_back(wv[r]);
				if (poisson) {
					ln.y.push_back(Y[c]);
					ln.e.push_back(std::exp(zz[c] + offset[t]));
				} else if (kind == 2) {
					// residual of the cell given its reciprocal partner, whose mean does not move with this path
					double rh = 0, pq = 0;
					if (mode == 0) {
						const size_t q = row ? cix(j, i, t) : cix(i, j, t);
						rh = std::isfinite(zz[q]) ? rho_t[t] : 0;
						pq = rh == 0 ? 0 : zz[q] - ez[q];
					}
					ln.b.push_back(zz[c] - ez[c] - rh * pq);
					ln.v.push_back(s2 * (1 - rh * rh));
				} else {
					if (mode == 0) ln.b.push_back(ez[c] + off[c]);
					else if (sym) ln.b.push_back(0.5 * (ez[c] + ez[cix(j, i, t)]));
					else ln.b.push_back(ez[c]);
					ln.s.push_back(mask[c]);
				}
			}
			for (int t = 0; t < nt; ++t) for (int r = 0; r < rk; ++r) cur[r + rk * t] = own[own_ix(i, r, t)];
			x = cur;
			auto ll_at = [&](const std::vector<double>& q) {
				for (int k = 0; k < dim; ++k) d[k] = q[k] - cur[k];
				return line_ll(ln, d, kind, sd_cell);
			};
			// one elliptical slice step from the reference root L with log likelihood f
			auto ess = [&](auto f) {
				for (int k = 0; k < dim; ++k) z[k] = norm_rand();
				for (int t = 0; t < nt; ++t) for (int r = 0; r < rk; ++r) {
					double s = 0;
					for (int q = 0; q <= t; ++q) s += L[t + nt * q] * z[r + rk * q];
					nu[r + rk * t] = s;
				}
				const double thr = f(x) + std::log(unif_rand());
				double th = 2 * M_PI * unif_rand(), lo = th - 2 * M_PI, hi = th;
				for (int ev = 1; ev <= 100; ++ev) {
					for (int k = 0; k < dim; ++k) prop[k] = x[k] * std::cos(th) + nu[k] * std::sin(th);
					const double lp_ = f(prop);
					evals += 1;
					if (std::isfinite(lp_) && lp_ > thr) { x = prop; return; }
					if (th < 0) lo = th; else hi = th;
					th = lo + (hi - lo) * unif_rand();
				}
				stuck += 1;
			};
			calls += 1;
			// slice step with the resets summed out, from a reference that resets at t_ref (0: none)
			const int t_ref = (int)(unif_rand() * nt);
			for (int t = 0; t < nt; ++t) dl_ref[t] = (t == t_ref && t > 0) ? 1 : 0;
			if (snap_root(nt, p, dl_ref, L)) {
				ess([&](const std::vector<double>& q) {
					if (!in_rail(q)) return NEG_INF;
					return ll_at(q) + snap_ratio(q, rk, nt, p, t_ref);
				});
			}
			// scale of the path from t0 on, with the Jacobian c^(rk (nt - t0))
			const int t0 = (int)(unif_rand() * nt);
			auto lt = [&](const std::vector<double>& q) {
				if (!in_rail(q)) return NEG_INF;
				return ll_at(q) + snap_log_prior(q, rk, nt, p);
			};
			double lt_x = lt(x);
			const double taus[3] = {0.1, 0.4, 1.0};
			for (double tau : taus) {
				const double eps = tau * norm_rand();
				prop = x;
				for (int k = rk * t0; k < dim; ++k) prop[k] = std::exp(eps) * x[k];
				const double lt_p = lt(prop);
				const double la = lt_p - lt_x + rk * (nt - t0) * eps;
				scales += 1;
				if (std::isfinite(la) && std::log(unif_rand()) < la) {
					x = prop;
					lt_x = lt_p;
					scales_acc += 1;
				}
			}
			// the resets, the path given them, the resets again
			snap_resets(x, rk, nt, p, dl);
			if (snap_root(nt, p, dl, L)) {
				ess([&](const std::vector<double>& q) {
					if (!in_rail(q)) return NEG_INF;
					return ll_at(q);
				});
				snap_resets(x, rk, nt, p, dl);
			}
			NumericMatrix& dd = row ? du : dv;
			for (int t = 0; t < nt; ++t) dd(i, t) = dl[t];
			// move the paths and the line's cells (every cell, the diagonal included)
			for (int k = 0; k < dim; ++k) d[k] = x[k] - cur[k];
			for (int t = 0; t < nt; ++t) {
				for (int r = 0; r < rk; ++r) own[own_ix(i, r, t)] = x[r + rk * t];
				for (int j = 0; j < n_line; ++j) {
					double dl_c = 0;
					if (sym && j == i) {
						// the diagonal cell is quadratic in the path
						for (int r = 0; r < rk; ++r) {
							const double a = x[r + rk * t], b = cur[r + rk * t];
							dl_c += Ld[r] * (a * a - b * b);
						}
					} else {
						weights(i, j, t);
						for (int r = 0; r < rk; ++r) dl_c += wv[r] * d[r + rk * t];
					}
					const size_t c = row ? cix(i, j, t) : cix(j, i, t);
					ez[c] += dl_c;
					if (poisson) zz[c] += dl_c;
					if (sym && j != i) {
						const size_t q = cix(j, i, t);
						ez[q] += dl_c;
						if (poisson) zz[q] += dl_c;
					}
				}
			}
		}
	}
	return List::create(
		Named("U") = uu,
		Named("V") = vv,
		Named("EZ") = ez,
		Named("Z") = zz,
		Named("delta_u") = du,
		Named("delta_v") = dv,
		Named("stats") = NumericVector::create(
			Named("dpm_calls") = calls,
			Named("dpm_evals") = evals,
			Named("dpm_stuck") = stuck,
			Named("snap_scales") = scales,
			Named("snap_scales_accepted") = scales_acc
		)
	);
}

// the visit of a directed panel (arguments as for snap_visit)
// [[Rcpp::export]]
List snap_path_visit_cpp(
	NumericVector EZ,
	NumericVector Z,
	NumericVector off,
	NumericVector mask,
	NumericVector Y,
	NumericVector offset,
	NumericVector U,
	NumericVector V,
	double rho,
	double sigma,
	double kappa,
	double pi,
	NumericMatrix delta_u,
	NumericMatrix delta_v,
	double sd_cell,
	double max_abs,
	double s2,
	NumericVector rho_t,
	int kind
) {
	return snap_visit(
		EZ, Z, off, mask, Y, offset, U, V, NumericMatrix(0, 0), NumericVector(0), rho, sigma, kappa, pi,
		delta_u, delta_v, sd_cell, max_abs, s2, rho_t, kind, 0
	);
}

// the visit of any panel: `mode` 0 directed, 1 bipartite (with G), 2 symmetric (with Ld; V not read)
// [[Rcpp::export]]
List snap_path_visit_mode_cpp(
	NumericVector EZ,
	NumericVector Z,
	NumericVector off,
	NumericVector mask,
	NumericVector Y,
	NumericVector offset,
	NumericVector U,
	NumericVector V,
	NumericMatrix G,
	NumericVector Ld,
	double rho,
	double sigma,
	double kappa,
	double pi,
	NumericMatrix delta_u,
	NumericMatrix delta_v,
	double sd_cell,
	double max_abs,
	double s2,
	NumericVector rho_t,
	int kind,
	int mode
) {
	if (mode < 0 || mode > 2) Rcpp::stop("Unknown panel mode of the snap path move.");
	if (mode == 1 && (G.nrow() != INTEGER(U.attr("dim"))[1] || G.ncol() != INTEGER(V.attr("dim"))[1]))
		Rcpp::stop("The interaction matrix does not match the bipartite paths.");
	if (mode == 2 && Ld.size() != INTEGER(U.attr("dim"))[1])
		Rcpp::stop("The symmetric weights do not match the paths.");
	return snap_visit(
		EZ, Z, off, mask, Y, offset, U, V, G, Ld, rho, sigma, kappa, pi,
		delta_u, delta_v, sd_cell, max_abs, s2, rho_t, kind, mode
	);
}
