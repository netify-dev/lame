// compiled visits of the position moves of the cbin and frn samplers: the static positions (R/sparse_position_move.R) and the
// coordinate paths of a dynamic_uv fit (R/nomination_path_move.R). A visit takes the senders (or receivers) in a given order;
// for each actor it draws the position (path) by an elliptical slice step from its prior with a log likelihood of zero inside the
// region where every shifted latent row (column) is still a row the data allow and minus infinity outside it. The arithmetic and the
// order of the random numbers are those of the R versions (.lame_ess_step and the closures of the visits), so that the compiled
// and the R engines give the same draws from the same seed. The classes of the cells (cls) are the ones of .lame_nom_col_prep for a
// panel without hidden cells: -1 free, 0 non-nominee, 1 nominee (cbin), the rank (frn).
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

static const double INF = std::numeric_limits<double>::infinity();

// the panel of one call: latent array Z (nr x nc x N), the classes, the kind (1 cbin, 2 frn), the rows with room, the rank holders
struct Panel {
	int nr, nc, N, kind, maxrank;
	double* Z;
	const int* cls;
	const int* room;
	const int* YL;
	inline size_t ix(int i, int j, int t) const { return (size_t)i + (size_t)j * nr + (size_t)t * nr * nc; }
};

// whether the row i of period t, with the values z, is a latent row the data allow (the condition of .lame_spm_row_ok)
static bool row_ok(const Panel& P, int i, int t, const double* z, std::vector<std::pair<int, double>>& chain) {
	double m0 = -INF, nmin = INF;
	bool any_nom = false;
	chain.clear();
	for (int j = 0; j < P.nc; ++j) {
		const int c = P.cls[P.ix(i, j, t)];
		if (c == 0) { if (z[j] > m0) m0 = z[j]; }
		else if (P.kind == 1 ? c == 1 : c >= 1) {
			any_nom = true;
			if (z[j] < nmin) nmin = z[j];
			if (P.kind == 2) chain.push_back(std::make_pair(c, z[j]));
		}
	}
	if (any_nom) {
		if (!(nmin > std::max(0.0, m0))) return false;
		if (P.kind == 2 && chain.size() > 1) {
			std::sort(chain.begin(), chain.end());
			for (size_t k = 1; k < chain.size(); ++k) if (chain[k].second - chain[k - 1].second <= 0) return false;
		}
	}
	if (P.room[i + (size_t)t * P.nr] && m0 >= 0) return false;
	return true;
}

// the interval (lo, hi) that the cell (i, j) of period t may take given the other cells of its row (.lame_nom_col_limits)
static void cell_limits(const Panel& P, int i, int j, int t, double& lo, double& hi) {
	lo = -INF; hi = INF;
	const int c = P.cls[P.ix(i, j, t)];
	if (c == -1) return;
	double m0 = -INF, nmin = INF;
	for (int q = 0; q < P.nc; ++q) {
		const int cq = P.cls[P.ix(i, q, t)];
		const double zq = P.Z[P.ix(i, q, t)];
		if (cq == 0) { if (zq > m0) m0 = zq; }
		else if (P.kind == 1 ? cq == 1 : cq >= 1) { if (zq < nmin) nmin = zq; }
	}
	const bool room = P.room[i + (size_t)t * P.nr] != 0;
	if (P.kind == 1) {
		if (c == 1) lo = std::max(0.0, m0);
		else if (c == 0) hi = room ? 0.0 : nmin;
	} else {
		auto holder = [&](int y) -> double {  // value of the cell of rank y in the row, NaN when absent
			const int col = P.YL[i + (size_t)(y - 1) * P.nr + (size_t)t * P.nr * P.maxrank];
			return col > 0 ? P.Z[P.ix(i, col - 1, t)] : NAN;
		};
		if (c == 0) {
			double h = holder(1);
			if (std::isnan(h)) h = INF;
			if (room) h = 0.0;
			hi = h;
		} else {
			if (c == 1) lo = std::max(0.0, m0);
			else { double v = holder(c - 1); lo = std::isnan(v) ? -INF : v; }
			if (c < P.maxrank) { double v = holder(c + 1); hi = std::isnan(v) ? INF : v; }
		}
	}
}

// one elliptical slice step with a log likelihood that is 0 inside the region (feasible(x)) and -Inf outside; the random numbers in the
// order of .lame_ess_step: the normals, the threshold uniform, the angle, then one uniform per shrink
template <class F>
static void ess_region(std::vector<double>& x, const std::vector<double>& mu, const std::vector<double>& nu, F feasible, int& evals, bool& stuck) {
	const int q = (int)x.size();
	const double thr = std::log(unif_rand());
	double th = 2 * M_PI * unif_rand();
	double lo = th - 2 * M_PI, hi = th;
	std::vector<double> xp(q);
	evals = 0; stuck = false;
	for (;;) {
		const double ct = std::cos(th), st = std::sin(th);
		for (int k = 0; k < q; ++k) xp[k] = mu[k] + (x[k] - mu[k]) * ct + nu[k] * st;
		evals++;
		if (feasible(xp) && 0.0 > thr) { x = xp; return; }
		if (th < 0) lo = th; else hi = th;
		if (evals >= 100) { stuck = true; return; }
		th = lo + (hi - lo) * unif_rand();
	}
}

// ---- static positions ------------------------------------------------------------------------------------------------------------

// one phase of the visit of the static positions. X (n x R) the positions being drawn (senders when rows = TRUE), W (n x R) the other
// set, reg (R x R) the regression of the prior mean on the partner's position, root (R x R) the lower Cholesky factor of the prior
// covariance, order (1-based) the visiting order. Z is updated in place (a copy is returned).
// [[Rcpp::export]]
List spm_nom_phase_cpp(NumericVector Z, IntegerVector cls, int kind, IntegerVector room, IntegerVector YL, int maxrank, int nr, int nc,
		NumericMatrix X, NumericMatrix W, NumericMatrix reg, NumericMatrix root, IntegerVector order, bool rows) {
	const int N = Z.size() / ((size_t)nr * nc);
	const int R = X.ncol();
	NumericVector Zc = clone(Z);
	NumericMatrix Xc = clone(X);
	Panel P{nr, nc, N, kind, maxrank, Zc.begin(), cls.begin(), room.begin(), YL.size() ? YL.begin() : nullptr};
	const int n = rows ? nr : nc;
	std::vector<double> cur(R), mu(R), nu(R), e(R), d(n), zrow(nc);
	std::vector<std::pair<int, double>> chain;
	double calls = 0, evals_tot = 0, stuck_tot = 0;
	for (int s = 0; s < order.size(); ++s) {
		const int a = order[s] - 1;
		for (int k = 0; k < R; ++k) cur[k] = Xc(a, k);
		for (int k = 0; k < R; ++k) { double m = 0; for (int q = 0; q < R; ++q) m += reg(k, q) * W(a, q); mu[k] = m; }
		for (int k = 0; k < R; ++k) e[k] = R::norm_rand();
		for (int k = 0; k < R; ++k) { double m = 0; for (int q = 0; q < R; ++q) m += root(k, q) * e[q]; nu[k] = m; }
		// the shift of the cells of the line for a candidate x: D_b = (x - cur)' W_b
		auto shifts = [&](const std::vector<double>& x) {
			for (int b = 0; b < n; ++b) { double m = 0; for (int k = 0; k < R; ++k) m += W(b, k) * (x[k] - cur[k]); d[b] = m; }
		};
		bool ok_state = true;
		std::vector<double> x = cur;
		int evals; bool stuck;
		if (rows) {
			auto feasible = [&](const std::vector<double>& xp) {
				shifts(xp);
				for (int t = 0; t < N; ++t) {
					for (int j = 0; j < nc; ++j) zrow[j] = P.Z[P.ix(a, j, t)] + d[j];
					if (!row_ok(P, a, t, zrow.data(), chain)) return false;
				}
				return true;
			};
			ess_region(x, mu, nu, feasible, evals, stuck);
		} else {
			auto feasible = [&](const std::vector<double>& xp) {
				shifts(xp);
				for (int t = 0; t < N; ++t) for (int i = 0; i < nr; ++i) {
					double lo, hi;
					cell_limits(P, i, a, t, lo, hi);
					const double zn = P.Z[P.ix(i, a, t)] + d[i];
					if (zn <= lo || zn >= hi) return false;
				}
				return true;
			};
			ess_region(x, mu, nu, feasible, evals, stuck);
		}
		(void) ok_state;
		calls += 1; evals_tot += evals; stuck_tot += stuck ? 1 : 0;
		bool moved = false;
		for (int k = 0; k < R; ++k) if (x[k] != cur[k]) moved = true;
		if (!moved) continue;
		shifts(x);
		for (int t = 0; t < N; ++t) for (int b = 0; b < n; ++b) {
			if (rows) P.Z[P.ix(a, b, t)] += d[b]; else P.Z[P.ix(b, a, t)] += d[b];
		}
		for (int k = 0; k < R; ++k) Xc(a, k) = x[k];
	}
	return List::create(Named("Z") = Zc, Named("X") = Xc, Named("calls") = calls, Named("evals") = evals_tot, Named("stuck") = stuck_tot);
}

// ---- coordinate paths ------------------------------------------------------------------------------------------------------------

// one phase of the visit of the coordinate paths. Ucube (nr x RA x N), Vcube (nc x RB x N), G (RA x RB; empty for a directed network, whose
// product is U V'), root (N x N lower Cholesky factor of the AR(1) covariance of a path), order (1-based), rows whether the senders are
// drawn, max_abs the bound of the coordinates
// [[Rcpp::export]]
List npm_nom_phase_cpp(NumericVector Z, IntegerVector cls, int kind, IntegerVector room, IntegerVector YL, int maxrank, int nr, int nc,
		NumericVector Ucube, NumericVector Vcube, int RA, int RB, NumericMatrix G, NumericMatrix root, IntegerVector order, bool rows, double max_abs) {
	const int N = Z.size() / ((size_t)nr * nc);
	const bool bip = G.nrow() > 0;
	NumericVector Zc = clone(Z), Uc = clone(Ucube), Vc = clone(Vcube);
	Panel P{nr, nc, N, kind, maxrank, Zc.begin(), cls.begin(), room.begin(), YL.size() ? YL.begin() : nullptr};
	const bool rail = std::isfinite(max_abs) && max_abs > 0;
	// the factor the moving side sees: rows: W_t[j, a] = sum_b V[j, b, t] G[a, b] (bipartite) or V[j, a, t]; columns: W_t[i, b] = sum_a U[i, a, t] G[a, b] or U[i, b, t]
	const int Rm = rows ? RA : RB;           // coordinates of the moving side
	const int np = rows ? nc : nr;           // partners (cells of the line)
	const int nm = rows ? nr : nc;           // moving actors
	NumericVector& Mv = rows ? Uc : Vc;
	const NumericVector& Ot = rows ? Vc : Uc;
	const int Ro = rows ? RB : RA;
	const int nOt = rows ? nc : nr;
	std::vector<double> W((size_t)np * Rm * N);
	auto build_W = [&]() {
		for (int t = 0; t < N; ++t) for (int a = 0; a < Rm; ++a) for (int b = 0; b < np; ++b) {
			double m;
			if (!bip) m = Ot[(size_t)b + (size_t)a * nOt + (size_t)t * nOt * Ro];
			else {
				m = 0;
				for (int c = 0; c < Ro; ++c) m += Ot[(size_t)b + (size_t)c * nOt + (size_t)t * nOt * Ro] * (rows ? G(a, c) : G(c, a));
			}
			W[(size_t)b + (size_t)a * np + (size_t)t * np * Rm] = m;
		}
	};
	build_W();
	const int q = Rm * N;
	std::vector<double> cur(q), x(q), nu(q), e(q), zrow(nc), d((size_t)np * N);
	std::vector<std::pair<int, double>> chain;
	double calls = 0, evals_tot = 0, stuck_tot = 0;
	for (int s = 0; s < order.size(); ++s) {
		const int a = order[s] - 1;
		// path of the actor as an Rm x N matrix, column-major (a_r + Rm * t)
		for (int t = 0; t < N; ++t) for (int r = 0; r < Rm; ++r) cur[r + (size_t)Rm * t] = Mv[(size_t)a + (size_t)r * nm + (size_t)t * nm * Rm];
		// the prior draw: normals as matrix(rnorm(q), Rm), then nu[r, t] = sum_u root[t, u] e[r, u]
		for (int k = 0; k < q; ++k) e[k] = R::norm_rand();
		for (int t = 0; t < N; ++t) for (int r = 0; r < Rm; ++r) { double m = 0; for (int u = 0; u <= t; ++u) m += root(t, u) * e[r + (size_t)Rm * u]; nu[r + (size_t)Rm * t] = m; }
		// the shifts d[b, t] = W_t[b, ] (x_t - cur_t)
		auto shifts = [&](const std::vector<double>& xx) {
			for (int t = 0; t < N; ++t) for (int b = 0; b < np; ++b) {
				double m = 0;
				for (int r = 0; r < Rm; ++r) m += W[(size_t)b + (size_t)r * np + (size_t)t * np * Rm] * (xx[r + (size_t)Rm * t] - cur[r + (size_t)Rm * t]);
				d[(size_t)b + (size_t)np * t] = m;
			}
		};
		auto feasible = [&](const std::vector<double>& xp) {
			if (rail) { double mx = 0; for (int k = 0; k < q; ++k) mx = std::max(mx, std::fabs(xp[k])); if (mx > max_abs) return false; }
			shifts(xp);
			if (rows) {
				for (int t = 0; t < N; ++t) {
					for (int j = 0; j < nc; ++j) zrow[j] = P.Z[P.ix(a, j, t)] + d[(size_t)j + (size_t)np * t];
					if (!row_ok(P, a, t, zrow.data(), chain)) return false;
				}
			} else {
				for (int t = 0; t < N; ++t) for (int i = 0; i < nr; ++i) {
					double lo, hi;
					cell_limits(P, i, a, t, lo, hi);
					const double zn = P.Z[P.ix(i, a, t)] + d[(size_t)i + (size_t)np * t];
					if (zn <= lo || zn >= hi) return false;
				}
			}
			return true;
		};
		// the ESS has prior mean 0 here: xp = cur cos + nu sin
		std::vector<double> mu0(q, 0.0);
		x = cur;
		int evals; bool stuck;
		ess_region(x, mu0, nu, feasible, evals, stuck);
		calls += 1; evals_tot += evals; stuck_tot += stuck ? 1 : 0;
		bool moved = false;
		for (int k = 0; k < q; ++k) if (x[k] != cur[k]) moved = true;
		if (!moved) continue;
		shifts(x);
		for (int t = 0; t < N; ++t) for (int b = 0; b < np; ++b) {
			if (rows) P.Z[P.ix(a, b, t)] += d[(size_t)b + (size_t)np * t]; else P.Z[P.ix(b, a, t)] += d[(size_t)b + (size_t)np * t];
		}
		for (int t = 0; t < N; ++t) for (int r = 0; r < Rm; ++r) Mv[(size_t)a + (size_t)r * nm + (size_t)t * nm * Rm] = x[r + (size_t)Rm * t];
		// the partner factors of the other phase see the new path; the line factor W of this phase does not change
	}
	return List::create(Named("Z") = Zc, Named("U") = Uc, Named("V") = Vc, Named("calls") = calls, Named("evals") = evals_tot, Named("stuck") = stuck_tot);
}
