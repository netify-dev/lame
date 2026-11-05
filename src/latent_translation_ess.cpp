// Elliptical slice steps of the Poisson translation moves of R/latent_translation_move.R, with the log-likelihood closure of .lame_ess_step
// replaced by its formula. The random numbers are drawn in the order of .lame_ess_step (rnorm for the auxiliary draw, runif for the slice
// height, runif for the angle, one runif per rejected angle) and the sums are accumulated in long double as sum() does.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// the level of one actor: prior N(mu, sd^2) for the shift from 0, log factor ys dl - s expm1(dl)
// [[Rcpp::export]]
List ess_pois_scalar_cpp(double ys, double s, double mu, double sd) {
	const double twopi = 2.0 * M_PI;
	const double nu = sd * norm_rand();
	const double xc = 0.0 - mu;
	const double thr = 0.0 + std::log(unif_rand());
	double theta = R::runif(0.0, twopi);
	double lo = theta - twopi, hi = theta;
	int evals = 0;
	for (;;) {
		const double xp = mu + xc * std::cos(theta) + nu * std::sin(theta);
		const double llp = ys * xp - s * std::expm1(xp);
		++evals;
		if (std::isfinite(llp) && llp > thr) return List::create(_["x"] = xp, _["evals"] = evals, _["stuck"] = false);
		if (theta < 0) lo = theta; else hi = theta;
		if (evals >= 100) return List::create(_["x"] = 0.0, _["evals"] = evals, _["stuck"] = true);
		theta = R::runif(lo, hi);
	}
}

// the position of one actor: state cur (R), prior N(mu, root root'), opp (n x R) the matrix that carries a change of the position to the
// shifts of the line's cells, w the weights, ys and s the counts and exp(z + offset) sums of the cells of the line, `skip` the cell (1-based,
// 0 for none) whose shift is held at zero (a symmetric network's diagonal)
// [[Rcpp::export]]
List ess_pois_uv_cpp(NumericVector cur, NumericVector mu, NumericMatrix root, NumericMatrix opp, NumericVector w, NumericVector ys, NumericVector s,
		int skip) {
	const int R = cur.size();
	const int n = opp.nrow();
	const double twopi = 2.0 * M_PI;
	bool trivial = true;
	for (int j = 0; j < n; ++j) if (ys[j] > 0 || s[j] > 0) { trivial = false; break; }
	std::vector<double> z(R), nu(R), xc(R), xp(R), t(R), d(n);
	for (int a = 0; a < R; ++a) z[a] = norm_rand();
	for (int a = 0; a < R; ++a) { double acc = 0; for (int b = 0; b < R; ++b) acc += root(a, b) * z[b]; nu[a] = acc; }
	for (int a = 0; a < R; ++a) xc[a] = cur[a] - mu[a];
	auto ll = [&](const std::vector<double>& x) {
		if (trivial) return 0.0;
		for (int a = 0; a < R; ++a) t[a] = w[a] * (x[a] - cur[a]);
		for (int j = 0; j < n; ++j) { double acc = 0; for (int a = 0; a < R; ++a) acc += opp(j, a) * t[a]; d[j] = acc; }
		if (skip > 0) d[skip - 1] = 0.0;
		long double acc = 0;
		for (int j = 0; j < n; ++j) acc += ys[j] * d[j] - s[j] * std::expm1(d[j]);
		return (double)acc;
	};
	const double thr = 0.0 + std::log(unif_rand());
	double theta = R::runif(0.0, twopi);
	double lo = theta - twopi, hi = theta;
	int evals = 0;
	for (;;) {
		for (int a = 0; a < R; ++a) xp[a] = mu[a] + xc[a] * std::cos(theta) + nu[a] * std::sin(theta);
		const double llp = ll(xp);
		++evals;
		if (std::isfinite(llp) && llp > thr) return List::create(_["x"] = NumericVector(xp.begin(), xp.end()), _["evals"] = evals, _["stuck"] = false);
		if (theta < 0) lo = theta; else hi = theta;
		if (evals >= 100) return List::create(_["x"] = clone(cur), _["evals"] = evals, _["stuck"] = true);
		theta = R::runif(lo, hi);
	}
}

// ---- whole phases of the Poisson translation visits (non-centred, static effects: ame() and the static lame() fits) ----

struct Ess1 { double x; int evals; bool stuck; };

static Ess1 ess1(double ys, double s, double mu, double sd) {
	const double twopi = 2.0 * M_PI;
	const double nu = sd * norm_rand();
	const double xc = 0.0 - mu;
	const double thr = 0.0 + std::log(unif_rand());
	double theta = R::runif(0.0, twopi);
	double lo = theta - twopi, hi = theta;
	int evals = 0;
	for (;;) {
		const double xp = mu + xc * std::cos(theta) + nu * std::sin(theta);
		const double llp = ys * xp - s * std::expm1(xp);
		++evals;
		if (std::isfinite(llp) && llp > thr) return {xp, evals, false};
		if (theta < 0) lo = theta; else hi = theta;
		if (evals >= 100) return {0.0, evals, true};
		theta = R::runif(lo, hi);
	}
}

// One margin of .lame_ltm_visit for a Poisson fit with static, uncentred levels. S the sums of exp(z + offset) over each actor's line,
// ys the counts of the line, E (symmetric only: the nr x nr x N array of exp(z + offset) on the observed cells, updated), x the levels of
// the margin and o those of the other margin (directed), mode 0 directed / 1 symmetric / 2 bipartite, iS the 2 x 2 precision of a directed
// pair of levels, v the variance of a symmetric or bipartite level, perm the actors in order. Returns the new levels, the shifts of the
// cells' predictor (dsh) and the counters (calls, evals, stuck, skipped)
// [[Rcpp::export]]
List ltm_level_phase_cpp(NumericVector S, NumericVector ys, NumericVector E, IntegerVector dimE, NumericVector x, NumericVector o, NumericMatrix iS,
		double v, int mode, bool row, IntegerVector perm) {
	NumericVector S2 = clone(S), x2 = clone(x), E2 = clone(E);
	const int n = x.size();
	NumericVector dsh(n);
	double calls = 0, evals = 0, stuck = 0, skipped = 0;
	const int N = mode == 1 ? dimE[2] : 0;
	const size_t nn = mode == 1 ? (size_t)dimE[0] : 0;
	for (int q = 0; q < perm.size(); ++q) {
		const int i = perm[q] - 1;
		double a1, a2;
		if (mode == 0) {
			const double li = x2[i], bi = o[i];
			if (row) { a1 = iS(0, 0) * li + iS(0, 1) * bi; a2 = iS(0, 0); }
			else { a1 = iS(1, 1) * li + iS(0, 1) * bi; a2 = iS(1, 1); }
		} else { a1 = x2[i] / v; a2 = 1.0 / v; }
		if (!(std::isfinite(a1) && std::isfinite(a2)) || a2 <= 0) { skipped += 1; continue; }
		Ess1 r = ess1(ys[i], S2[i], -a1 / a2, 1.0 / std::sqrt(a2));
		calls += 1; evals += r.evals; stuck += r.stuck ? 1 : 0;
		const double d = r.x;
		if (d == 0) continue;
		const double ed = std::exp(d);
		if (mode == 1) {
			// the pairs {i, j} are shared with actor j: its sum gains the change of the pair
			std::vector<double> rs(nn);
			for (size_t j = 0; j < nn; ++j) {
				long double a = 0;
				for (int t = 0; t < N; ++t) a += E2[i + j * nn + (size_t)t * nn * nn];
				rs[j] = (double)a;
			}
			for (size_t j = 0; j < nn; ++j) S2[j] = S2[j] + (ed - 1) * rs[j];
			std::vector<double> ei((size_t)nn * N);
			for (int t = 0; t < N; ++t) for (size_t j = 0; j < nn; ++j) ei[j + (size_t)t * nn] = E2[i + j * nn + (size_t)t * nn * nn];
			for (int t = 0; t < N; ++t) for (size_t j = 0; j < nn; ++j) E2[i + j * nn + (size_t)t * nn * nn] = ei[j + (size_t)t * nn] * ed;
			for (int t = 0; t < N; ++t) for (size_t j = 0; j < nn; ++j) E2[j + i * nn + (size_t)t * nn * nn] = ei[j + (size_t)t * nn] * ed;
		}
		S2[i] = S2[i] * ed;
		dsh[i] = d;
		x2[i] = x2[i] + d;
	}
	return List::create(_["x"] = x2, _["dsh"] = dsh, _["stats"] = NumericVector::create(calls, evals, stuck, skipped));
}

// One margin of .lame_ltm_uv_visit (Poisson). mode 0 directed, 1 symmetric (V = U diag(w)), 2 bipartite. mov the positions that move
// (U for the rows, V for the columns), oth the other positions (directed: the prior mean of a position is reg %*% oth[i, ]), opp the matrix
// that carries a change of a position to the shifts of the cells of its line (not used for a symmetric network, where it is mov itself),
// reg and root the conditional prior of the position (reg is empty for a zero mean), w the weights, sj the sums of exp(z + offset) over the
// periods, ysj the counts over the periods, row whether the line is a row. Returns the moved positions, the updated sums, the shifts of the cells
// and the counters (calls, evals, stuck)
// [[Rcpp::export]]
List ltm_uv_phase_cpp(NumericMatrix mov, NumericMatrix oth, NumericMatrix opp, NumericMatrix reg, NumericMatrix root, NumericVector w,
		NumericMatrix sj, NumericMatrix ysj, int mode, bool row, IntegerVector perm) {
	NumericMatrix M = clone(mov), SJ = clone(sj);
	const int nr = sj.nrow(), nc = sj.ncol();
	const int R = M.ncol();
	NumericMatrix Dsh(nr, nc);
	const bool sym = mode == 1;
	const bool use_reg = reg.nrow() > 0;
	double calls = 0, evals = 0, stuck = 0;
	const double twopi = 2.0 * M_PI;
	std::vector<double> cur(R), mu(R), z(R), nu(R), xc(R), xp(R), t(R), dx(R);
	for (int q = 0; q < perm.size(); ++q) {
		const int i = perm[q] - 1;
		const NumericMatrix& O = sym ? M : opp;
		const int nl = O.nrow();
		for (int a = 0; a < R; ++a) cur[a] = M(i, a);
		for (int a = 0; a < R; ++a) {
			if (use_reg) { double acc = 0; for (int b = 0; b < R; ++b) acc += reg(a, b) * oth(i, b); mu[a] = acc; } else mu[a] = 0.0;
		}
		std::vector<double> ys(nl), s(nl), d(nl);
		for (int j = 0; j < nl; ++j) {
			if (row || sym) { ys[j] = ysj(i, j); s[j] = SJ(i, j); } else { ys[j] = ysj(j, i); s[j] = SJ(j, i); }
		}
		bool trivial = true;
		for (int j = 0; j < nl; ++j) if (ys[j] > 0 || s[j] > 0) { trivial = false; break; }
		// the draw of .lame_ess_step
		for (int a = 0; a < R; ++a) z[a] = norm_rand();
		for (int a = 0; a < R; ++a) { double acc = 0; for (int b = 0; b < R; ++b) acc += root(a, b) * z[b]; nu[a] = acc; }
		for (int a = 0; a < R; ++a) xc[a] = cur[a] - mu[a];
		auto ll = [&](const std::vector<double>& x) {
			if (trivial) return 0.0;
			for (int a = 0; a < R; ++a) t[a] = w[a] * (x[a] - cur[a]);
			for (int j = 0; j < nl; ++j) { double acc = 0; for (int a = 0; a < R; ++a) acc += O(j, a) * t[a]; d[j] = acc; }
			if (sym) d[i] = 0.0;
			long double acc = 0;
			for (int j = 0; j < nl; ++j) acc += ys[j] * d[j] - s[j] * std::expm1(d[j]);
			return (double)acc;
		};
		const double thr = 0.0 + std::log(unif_rand());
		double theta = R::runif(0.0, twopi);
		double lo = theta - twopi, hi = theta;
		int ev = 0;
		bool moved = false;
		for (;;) {
			for (int a = 0; a < R; ++a) xp[a] = mu[a] + xc[a] * std::cos(theta) + nu[a] * std::sin(theta);
			const double llp = ll(xp);
			++ev;
			if (std::isfinite(llp) && llp > thr) { moved = true; break; }
			if (theta < 0) lo = theta; else hi = theta;
			if (ev >= 100) { stuck += 1; break; }
			theta = R::runif(lo, hi);
		}
		calls += 1; evals += ev;
		if (!moved) continue;
		bool any = false;
		for (int a = 0; a < R; ++a) { dx[a] = xp[a] - cur[a]; if (dx[a] != 0) any = true; }
		if (!any) continue;
		for (int a = 0; a < R; ++a) t[a] = w[a] * dx[a];
		for (int j = 0; j < nl; ++j) { double acc = 0; for (int a = 0; a < R; ++a) acc += O(j, a) * t[a]; d[j] = acc; }
		if (sym) {
			long double dgl = 0;
			for (int a = 0; a < R; ++a) dgl += w[a] * (xp[a] * xp[a] - cur[a] * cur[a]);
			const double dg = (double)dgl;
			d[i] = 0.0;
			for (int j = 0; j < nl; ++j) { SJ(i, j) = SJ(i, j) * std::exp(d[j]); }
			for (int j = 0; j < nl; ++j) { SJ(j, i) = SJ(j, i) * std::exp(d[j]); }
			for (int j = 0; j < nl; ++j) Dsh(i, j) = Dsh(i, j) + d[j];
			for (int j = 0; j < nl; ++j) Dsh(j, i) = Dsh(j, i) + d[j];
			Dsh(i, i) = Dsh(i, i) + dg;
		} else if (row) {
			for (int j = 0; j < nl; ++j) { SJ(i, j) = SJ(i, j) * std::exp(d[j]); Dsh(i, j) = d[j]; }
		} else {
			for (int j = 0; j < nl; ++j) { SJ(j, i) = SJ(j, i) * std::exp(d[j]); Dsh(j, i) = d[j]; }
		}
		for (int a = 0; a < R; ++a) M(i, a) = xp[a];
	}
	return List::create(_["M"] = M, _["sj"] = SJ, _["Dsh"] = Dsh, _["stats"] = NumericVector::create(calls, evals, stuck));
}

// the shift of a coefficient (or of one period of its path) in a Poisson fit: prior N(mu, sd^2) for the shift from 0, log factor
// d ysx - sum(ez expm1(d x)) over the cells (the sum in double precision as ppm_shift_sum_cpp forms it)
// [[Rcpp::export]]
List ess_pois_cells_cpp(const NumericVector& ez, const NumericVector& x, double ysx, double mu, double sd) {
	const double twopi = 2.0 * M_PI;
	const R_xlen_t n = ez.size();
	if (x.size() != n) Rcpp::stop("ess_pois_cells_cpp: the cells and the design differ in length");
	const double nu = sd * norm_rand();
	const double xc = 0.0 - mu;
	const double thr = 0.0 + std::log(unif_rand());
	double theta = R::runif(0.0, twopi);
	double lo = theta - twopi, hi = theta;
	int evals = 0;
	for (;;) {
		const double xp = mu + xc * std::cos(theta) + nu * std::sin(theta);
		double sm = 0.0;
		for (R_xlen_t i = 0; i < n; ++i) sm += ez[i] * std::expm1(xp * x[i]);
		const double llp = xp * ysx - sm;
		++evals;
		if (std::isfinite(llp) && llp > thr) return List::create(_["x"] = xp, _["evals"] = evals, _["stuck"] = false);
		if (theta < 0) lo = theta; else hi = theta;
		if (evals >= 100) return List::create(_["x"] = 0.0, _["evals"] = evals, _["stuck"] = true);
		theta = R::runif(lo, hi);
	}
}
